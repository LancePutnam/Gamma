/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "AudioIO.h"
#include "scl.h"

#define ULONG unsigned long
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


void (* AudioIO::callback)(AudioIOData &) = 0;

AudioIO::AudioIO(
	ULONG framesPerBuf, double framesPerSec, void (* callbackA)(AudioIOData &), void * user,
	int outChansA, int inChansA )
:	AudioIOData(user),
	mIsOpen(false), mIsRunning(false), mInResizeDeferred(false), mOutResizeDeferred(false),
	mKillNANs(true)
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
	Pa_Terminate();
}


void AudioIO::init(){	
	initDevices();

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


int AudioIO::initDevices(){ return Pa_Initialize(); }


void AudioIO::auxChans(ULONG num){
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
		for(int i=0; i<io.numFrames()*io.outDeviceChans(); ++i){
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
	ULONG * chans   = forOutput ? &mNumO : &mNumI;
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
                
		if(!supportsFPS(v)) v = Pa_GetDeviceInfo(mOutDevice)->defaultSampleRate;

		mFramesPerSecond = v;
		reopen();
	}
}


void AudioIO::setFramesPerBuffer(ULONG num){
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

//void AudioIO::virtualChans(ULONG num, bool forOutput){
//
//	ULONG * currNum = forOutput ? &mVOChans : &mVIChans;
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
	PaDeviceIndex dI = mInParams.device;
	PaDeviceIndex dO = mOutParams.device;
	if(dI == dO){
		printf("I/O Device:  ");	printDevice(dI);
	}
	else{
		printf("In Device:   ");	printDevice(dI);
		printf("Out Device:  ");	printDevice(dO);
	}

		printf("In Chans:    %lu (%luD + %luV)\n", inChans(), inDeviceChans(), inChans() - inDeviceChans());
		printf("Out Chans:   %lu (%luD + %luV)\n", outChans(), outDeviceChans(), outChans() - outDeviceChans());

	const PaStreamInfo * sInfo = Pa_GetStreamInfo(mStream);
	if(sInfo){
		printf("In Latency:  %.0f ms\nOut Latency: %0.f ms\nSample Rate: %0.f Hz\n",
			sInfo->inputLatency * 1000., sInfo->outputLatency * 1000., sInfo->sampleRate);
	}
	printf("Frames/Buf:  %lu\n", mFramesPerBuffer);
}

void AudioIO::printDevice(int deviceNo){

	if(deviceNo == paNoDevice){ printf("No device\n"); return; }

	const PaDeviceInfo * info = Pa_GetDeviceInfo(deviceNo);
	
	if(0 == info){ printf("Invalid device\n"); return; }

	printf("[%2d] %s, ", deviceNo, info->name);
	
	int chans = info->maxInputChannels;
	if(chans > 0) printf("%2i in, ", chans);
	chans = info->maxOutputChannels;
	if(chans > 0) printf("%2i out, ", chans);
	
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

	printf("%.0f Hz\n", info->defaultSampleRate);
}

void AudioIO::printDevices(){
	for(int i=0; i<numDevices(); i++){
		printf("[%2d] ", i);
		printDevice(i);
	}
}

void AudioIO::printError(){
	printf("%s \n", errorText(mErrNum));
}

} // end namespace gam

#undef ULONG
