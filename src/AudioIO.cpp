/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "AudioIO.h"
#include "scl.h"

#define SAFE_FREE(ptr) if(ptr){ free(ptr); ptr = 0; }

namespace gam{


AudioIOData::AudioIOData(void * user) :
	user(user),
	mStream(0),
	mFramesPerBuffer(0), mFramesPerSecond(0),
	mBufI(0), mBufO(0), mBufA(0), mNumI(0), mNumO(0), mNumA(0)
{
}

AudioIOData::~AudioIOData(){
	SAFE_FREE(mBufI);
	SAFE_FREE(mBufO);
	SAFE_FREE(mBufA);
}


AudioDevice::AudioDevice(int deviceNum)
: mID(-1), mImpl(0){
	setImpl(deviceNum);
}

AudioDevice::~AudioDevice(){}

const char * AudioDevice::name() const { return ((const PaDeviceInfo*)mImpl)->name; }
int AudioDevice::maxInputChannels() const { return ((const PaDeviceInfo*)mImpl)->maxInputChannels; }
int AudioDevice::maxOutputChannels() const { return ((const PaDeviceInfo*)mImpl)->maxOutputChannels; }
double AudioDevice::defaultSampleRate() const { return ((const PaDeviceInfo*)mImpl)->defaultSampleRate; }
void AudioDevice::setImpl(int deviceNum){ initDevices(); mImpl = Pa_GetDeviceInfo(deviceNum); mID=deviceNum; }
AudioDevice AudioDevice::defaultInput(){ return AudioDevice(Pa_GetDefaultInputDevice()); }
AudioDevice AudioDevice::defaultOutput(){ return AudioDevice(Pa_GetDefaultOutputDevice()); }

struct InitSingleton{
	InitSingleton(): mCleanUp(paNoError == Pa_Initialize()){}
	~InitSingleton(){ if(mCleanUp){ Pa_Terminate(); } }
	bool mCleanUp;
};

void AudioDevice::initDevices(){
	static InitSingleton dummy;
}

int AudioDevice::numDevices(){ initDevices(); return Pa_GetDeviceCount(); }

void AudioDevice::print() const{

	//if(deviceNum == paNoDevice){ printf("No device\n"); return; }

	//const AudioDevice dev(deviceNum);
	if(!valid()){ printf("Invalid device\n"); return; }

	printf("[%2d] %s, ", id(), name());
	
	int chans = maxInputChannels();
	if(chans > 0) printf("%2i in, ", chans);
	chans = maxOutputChannels();
	if(chans > 0) printf("%2i out, ", chans);

	printf("%.0f Hz\n", defaultSampleRate());
	
//	PaSampleFormat sampleFormats = info->nativeSampleFormats;
//	
//	printf("[ ");
//	if(0 != sampleFormats & paFloat32)		printf("f32 ");
//	if(0 != sampleFormats & paInt32)		printf("i32 ");
//	if(0 != sampleFormats & paInt24)		printf("i24 ");
//	if(0 != sampleFormats & paInt16)		printf("i16 ");
//	if(0 != sampleFormats & paInt8)			printf("i8 ");
//	if(0 != sampleFormats & paUInt8)		printf("ui8 ");
//	printf("], ");
	
//	if(info->numSampleRates != -1){
//		printf("[");
//		for(int i=0; i<info->numSampleRates; i++){
//			printf("%f ", info->sampleRates[i]);
//		}
//		printf("] Hz");
//	}
//	else{
//		printf("[%.0f <-> %.0f] Hz", info->sampleRates[0], info->sampleRates[1]);
//	}
//	printf("\n");
}

void AudioDevice::printAll(){
	for(int i=0; i<numDevices(); i++){
		printf("[%2d] ", i);
		AudioDevice dev(i);
		dev.print();
		//print(i);
	}
}




void (* AudioIO::callback)(AudioIOData &) = 0;

AudioIO::AudioIO(
	uint32_t framesPerBuf, double framesPerSec, void (* callbackA)(AudioIOData &), void * user,
	int outChansA, int inChansA )
:	AudioIOData(user),
	mIsOpen(false), mIsRunning(false), mInResizeDeferred(false), mOutResizeDeferred(false),
	mKillNANs(true),
	mOutDevice(AudioDevice::defaultOutput()), mInDevice(AudioDevice::defaultInput())
{
	callback = callbackA;
	
	init();
	this->setFramesPerBuffer(framesPerBuf);
	chans(inChansA, false);
	chans(outChansA, true);
	this->setFPS(framesPerSec);
}

		
AudioIO::~AudioIO(){
	close();
}


void AudioIO::init(){

	// Choose default devices for now...
	inDevice(defaultInDevice());
	outDevice(defaultOutDevice());
	
	// Setup input stream parameters
	const PaDeviceInfo * dInfo = Pa_GetDeviceInfo(mInParams.device);
	if(dInfo) mInParams.suggestedLatency = dInfo->defaultLowInputLatency; // for RT
	mInParams.sampleFormat = paFloat32;// | paNonInterleaved;
	//mInParams.sampleFormat = paInt16;
	mInParams.hostApiSpecificStreamInfo = NULL;

	// Setup output stream parameters
	dInfo = Pa_GetDeviceInfo(mOutParams.device);
	if(dInfo) mOutParams.suggestedLatency = dInfo->defaultLowOutputLatency; // for RT
	mOutParams.sampleFormat = paFloat32;// | paNonInterleaved;
	mOutParams.hostApiSpecificStreamInfo = NULL;

	setInDeviceChans(0);
	setOutDeviceChans(0);
}


void AudioIO::auxChans(uint32_t num){
	if(mem::resize(mBufA, mNumA * mFramesPerBuffer, num * mFramesPerBuffer)){
		mNumA = num;
	}
}


void AudioIO::chans(int num, bool forOutput){
	
	PaStreamParameters * params = forOutput ? &mOutParams : &mInParams;
	
	if(num == 0){
		params->device = paNoDevice;
		params->channelCount = 0;
		return;
	}

	const PaDeviceInfo * info = Pa_GetDeviceInfo(params->device);
	if(0 == info) return;	// this particular device is not open, so return

	// compute number of channels to give PortAudio
	int maxChans = 
		(int)(forOutput ? info->maxOutputChannels : info->maxInputChannels);
	
	// -1 means open all channels
	if(num == -1){
		num = maxChans;
	}
	
	int currentNum = chans(forOutput);
	
	if(num != currentNum){

		params->channelCount = scl::min(num, maxChans);

		forOutput ? mNumO = num : mNumI = num;
		
		deferBufferResize(forOutput);
	}
}


bool AudioIO::close(){
	mErrNum = paNoError;
	
	if(mIsOpen) mErrNum = Pa_CloseStream(mStream);	
	
	if(paNoError == mErrNum){
		mIsOpen = false;
		mIsRunning = false;
	}
	
	return paNoError == mErrNum;
}


void AudioIO::deferBufferResize(bool forOutput){
	if(forOutput)	mOutResizeDeferred = true;
	else			mInResizeDeferred = true;
}


bool AudioIO::open(){

	mErrNum = paNoError;

	if(!(mIsOpen || mIsRunning)){

		resizeBuffer(false);
		resizeBuffer(true);
		
		PaStreamParameters * inParams = &mInParams;
		PaStreamParameters * outParams = &mOutParams;
		
		// Must pass in 0s for input- or output-only streams.
		if(paNoDevice ==  inParams->device) inParams  = 0;
		if(paNoDevice == outParams->device) outParams = 0;

		mErrNum = Pa_OpenStream(
			&mStream,			// PortAudioStream **
			inParams,			// PaStreamParameters * in
			outParams,			// PaStreamParameters * out
			mFramesPerSecond,	// frames/sec (double)
			mFramesPerBuffer,	// frames/buffer (unsigned long)
            paNoFlag,			// paNoFlag, paClipOff, paDitherOff
			paCallback,			// static callback function (PaStreamCallback *)
			this
		);

		mIsOpen = paNoError == mErrNum;
	}
	//printf("AudioIO::open()\n"); printError();
	return paNoError == mErrNum;
}


int AudioIO::paCallback(const void *input,
						void *output,
						unsigned long frameCount,
						const PaStreamCallbackTimeInfo* timeInfo,
						PaStreamCallbackFlags statusFlags,
						void * userData )
{

	AudioIO& io = *(AudioIO *)userData;

	const float * paI = (const float *)input;
	float * paO = (float *)output;

	bool deinterleave = true;

	if(deinterleave){
		mem::deinterleave((float *)io.in(0),  paI, io.numFrames(), io.inDeviceChans() );
		mem::deinterleave(io.out(0), paO, io.numFrames(), io.outDeviceChans());
	}

	io();	// call callback

	// kill pesky nans so we don't hurt anyone's ears
	if(io.killNANs()){
		for(uint32_t i=0; i<io.numFrames()*io.outDeviceChans(); ++i){
			float& s = io.out(0)[i];
			if(isnan(s)) s = 0.f;
		}
	}

	if(deinterleave){
		mem::interleave(paO, io.out(0), io.numFrames(), io.outDeviceChans());
	}

	return 0;
}


void AudioIO::reopen(){
	if(mIsRunning)  { close(); start(); }
	else if(mIsOpen){ close(); open(); }
}

void AudioIO::resizeBuffer(bool forOutput){

	float ** buffer = forOutput ? &mBufO : &mBufI;
	uint32_t * chans   = forOutput ? &mNumO : &mNumI;
	bool * deferred = forOutput ? &mOutResizeDeferred : &mInResizeDeferred;

	if(*deferred){
		if(*chans > 0){
			float * ptr = (float *)realloc(*buffer, *chans * mFramesPerBuffer * sizeof(float));
			if(ptr){
				*buffer = ptr;
				*deferred = false;
			}
			else{
				*chans = 0;
				*buffer = 0;
			}
		}
		else{
			SAFE_FREE(*buffer);
			*deferred = false;
		}
	}
}


void AudioIO::setFPS(double v){	//printf("AudioIO::fps(%f)\n", v);
	if(AudioIOData::fps() != v){
                
		if(!supportsFPS(v)) v = mOutDevice.defaultSampleRate();

		mFramesPerSecond = v;
		reopen();
	}
}


void AudioIO::setFramesPerBuffer(uint32_t num){
	if(numFrames() != num){
		mFramesPerBuffer = num;
		auxChans(AudioIOData::auxChans());
		reopen();
	}
}


bool AudioIO::start(){
	mErrNum = paNoError;
	
	if(!mIsOpen) open();
	if(mIsOpen && !mIsRunning)	mErrNum = Pa_StartStream(mStream);
	if(paNoError == mErrNum)	mIsRunning = true;
	
	return paNoError == mErrNum;
}


bool AudioIO::stop(){
	mErrNum = paNoError;
	
	if(mIsRunning)				mErrNum = Pa_StopStream(mStream);
	if(paNoError == mErrNum)	mIsRunning = false;
	
	return paNoError == mErrNum;
}


bool AudioIO::supportsFPS(double fps){

	PaStreamParameters * pi = AudioIOData::inDeviceChans() == 0 ? 0 : &mInParams;
	PaStreamParameters * po = AudioIOData::outDeviceChans() == 0 ? 0 : &mOutParams;	
	mErrNum = Pa_IsFormatSupported(pi, po, fps);
	
	if(error()){ printf("AudioIO error: "); printError(); }
	
	return paFormatIsSupported == mErrNum;
}

//void AudioIO::virtualChans(uint32_t num, bool forOutput){
//
//	uint32_t * currNum = forOutput ? &mVOChans : &mVIChans;
//
//	if(num != *currNum){
//		*currNum = num;
//		deferBufferResize(forOutput);
//	}
//}

//void assignBufferAccessors(){
//	for(int i=0; i<mVIChans; ++i) mAIBufs[i + mDIChans] = mVIBufs + i * mFramesPerBuffer;
//	for(int i=0; i<mVOChans; ++i) mAOBufs[i + mDOChans] = mVOBufs + i * mFramesPerBuffer;
//}


void AudioIO::print(){
	if(mInDevice.id() == mOutDevice.id()){
		printf("I/O Device:  "); mInDevice.print();
	}
	else{
		printf("In Device:   "); mInDevice.print();
		printf("Out Device:  "); mOutDevice.print();
	}

		printf("In Chans:    %d (%dD + %dV)\n", inChans(), inDeviceChans(), inChans() - inDeviceChans());
		printf("Out Chans:   %d (%dD + %dV)\n", outChans(), outDeviceChans(), outChans() - outDeviceChans());

	const PaStreamInfo * sInfo = Pa_GetStreamInfo(mStream);
	if(sInfo){
		printf("In Latency:  %.0f ms\nOut Latency: %0.f ms\nSample Rate: %0.f Hz\n",
			sInfo->inputLatency * 1000., sInfo->outputLatency * 1000., sInfo->sampleRate);
	}
	printf("Frames/Buf:  %d\n", mFramesPerBuffer);
}


void AudioIO::printError(){
	printf("%s \n", errorText(mErrNum));
}

} // gam::
