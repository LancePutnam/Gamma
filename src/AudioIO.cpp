/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>		/* memset() */
#include <strings.h>
#include <math.h>

#include "portaudio.h"
#include "Gamma/AudioIO.h"

namespace gam{

/*
static void err(const char * msg, const char * src, bool exits){
	fprintf(stderr, "%s%serror: %s\n", src, src[0]?" ":"", msg);
	if(exits) exit(EXIT_FAILURE);
}
*/

static void warn(const char * msg, const char * src){
	fprintf(stderr, "%s%swarning: %s\n", src, src[0]?" ":"", msg);
}

template <class T>
static void deleteBuf(T *& buf){ delete[] buf; buf=0; }

static inline int min(int x, int y){ return x<y?x:y; }

template <class T>
static int resize(T *& buf, int n){
	deleteBuf(buf);
	buf = new T[n];
	return n;
}

template <class T>
static inline void zero(T * buf, int n){ bzero(buf, n*sizeof(T)); }

template <class T>
static void deinterleave(T * dst, const T * src, int numFrames, int numChannels){
	int numSamples = numFrames * numChannels;
	for(int c=0; c < numChannels; c++){
		for(int i=c; i < numSamples; i+=numChannels){
			*dst++ = src[i];
		}
	}
}

template <class T>
static void interleave(T * dst, const T * src, int numFrames, int numChannels){
	int numSamples = numFrames * numChannels;
	for(int c=0; c < numChannels; c++){
		for(int i=c; i < numSamples; i+=numChannels){
			dst[i] = *src++;
		}
	}
}


//==============================================================================
AudioDevice::AudioDevice(int deviceNum)
:	mID(-1), mImpl(0)
{
	setImpl(deviceNum);
}

AudioDevice::AudioDevice(const std::string& nameKeyword, bool input, bool output)
:	mID(-1), mImpl(0)
{
	for(int i=0; i<numDevices(); ++i){
		AudioDevice d(i);
		std::string n = d.name();
		if(	((input & d.hasInput()) || (output & d.hasOutput())) &&
			n.find(nameKeyword) != std::string::npos
		){
			setImpl(i);
			break;
		}
	}
}

AudioDevice::~AudioDevice(){}

const char * AudioDevice::name() const { return ((const PaDeviceInfo*)mImpl)->name; }
int AudioDevice::channelsInMax() const { return ((const PaDeviceInfo*)mImpl)->maxInputChannels; }
int AudioDevice::channelsOutMax() const { return ((const PaDeviceInfo*)mImpl)->maxOutputChannels; }
double AudioDevice::defaultSampleRate() const { return ((const PaDeviceInfo*)mImpl)->defaultSampleRate; }
bool AudioDevice::hasInput() const { return channelsInMax()>0; }
bool AudioDevice::hasOutput() const { return channelsOutMax()>0; }
void AudioDevice::setImpl(int deviceNum){ initDevices(); mImpl = Pa_GetDeviceInfo(deviceNum); mID=deviceNum; }
AudioDevice AudioDevice::defaultInput(){ initDevices(); return AudioDevice(Pa_GetDefaultInputDevice()); }
AudioDevice AudioDevice::defaultOutput(){ initDevices(); return AudioDevice(Pa_GetDefaultOutputDevice()); }

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
	
	int chans = channelsInMax();
	if(chans > 0) printf("%2i in, ", chans);
	chans = channelsOutMax();
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


//==============================================================================
struct AudioIOData::Impl{
	Impl(): mStream(0), mErrNum(0), mIsOpen(false), mIsRunning(false){}

	bool error() const { return mErrNum != paNoError; }

	void printError() const { printf("%s \n", Pa_GetErrorText(mErrNum)); }

	bool supportsFPS(double fps) const {
		const PaStreamParameters * pi = mInParams.channelCount  == 0 ? 0 : &mInParams;
		const PaStreamParameters * po = mOutParams.channelCount == 0 ? 0 : &mOutParams;	
		mErrNum = Pa_IsFormatSupported(pi, po, fps);
		if(error()){ printf("AudioIO error: "); printError(); }
		return paFormatIsSupported == mErrNum;
	}

	void inDevice(PaDeviceIndex index){ mInParams.device = index; }
	void outDevice(PaDeviceIndex index){ mOutParams.device = index; }
	void setInDeviceChans(int num){ mInParams.channelCount = num; }
	void setOutDeviceChans(int num){mOutParams.channelCount = num; }

	bool close(){
		mErrNum = paNoError;
		if(mIsOpen) mErrNum = Pa_CloseStream(mStream);
		if(paNoError == mErrNum){
			mIsOpen = false;
			mIsRunning = false;
		}
		return paNoError == mErrNum;
	}

	bool stop(){
		mErrNum = paNoError;
		if(mIsRunning)				mErrNum = Pa_StopStream(mStream);
		if(paNoError == mErrNum)	mIsRunning = false;
		return paNoError == mErrNum;
	}

	PaStreamParameters mInParams, mOutParams;	// Input and output stream parameters
	PaStream * mStream;					// i/o stream
	mutable PaError mErrNum;			// Most recent error number
	bool mIsOpen;						// An audio device is open
	bool mIsRunning;					// An audio stream is running
};

AudioIOData::AudioIOData(void * userData)
:	mImpl(new Impl),
	mUser(userData),
	mFramesPerBuffer(0), mFramesPerSecond(0),
	mBufI(0), mBufO(0), mBufB(0), mBufT(0), mNumI(0), mNumO(0), mNumB(0)
{}

AudioIOData::~AudioIOData(){
	deleteBuf(mBufI);
	deleteBuf(mBufO);
	deleteBuf(mBufB);
	deleteBuf(mBufT);
}

void AudioIOData::zeroBus(){ zero(mBufB, framesPerBuffer() * mNumB); }
void AudioIOData::zeroOut(){ zero(mBufO, channelsOut() * framesPerBuffer()); }

int AudioIOData::channelsIn () const { return mNumI; }
int AudioIOData::channelsOut() const { return mNumO; }
int AudioIOData::channelsBus() const { return mNumB; }
int AudioIOData::channelsInDevice() const { return (int)mImpl->mInParams.channelCount; }
int AudioIOData::channelsOutDevice() const { return (int)mImpl->mOutParams.channelCount; }

double AudioIOData::framesPerSecond() const { return mFramesPerSecond; }
double AudioIOData::time() const { return Pa_GetStreamTime(mImpl->mStream); }
double AudioIOData::time(int frame) const { return (double)frame / framesPerSecond() + time(); }
int AudioIOData::framesPerBuffer() const { return mFramesPerBuffer; }
double AudioIOData::secondsPerBuffer() const { return (double)framesPerBuffer() / framesPerSecond(); }


//==============================================================================
static int paCallback(	const void *input,
						void *output,
						unsigned long frameCount,
						const PaStreamCallbackTimeInfo* timeInfo,
						PaStreamCallbackFlags statusFlags,
						void *userData );

AudioIO::AudioIO(
	int framesPerBuf, double framesPerSec, void (* callbackA)(AudioIOData &), void * userData,
	int outChansA, int inChansA )
:	AudioIOData(userData),
	callback(callbackA),
	mInDevice(AudioDevice::defaultInput()), mOutDevice(AudioDevice::defaultOutput()),
	mInResizeDeferred(false), mOutResizeDeferred(false),
	mZeroNANs(true), mClipOut(true), mAutoZeroOut(true)
{
	init();
	this->framesPerBuffer(framesPerBuf);
	channels(inChansA, false);
	channels(outChansA, true);
	this->framesPerSecond(framesPerSec);
}

		
AudioIO::~AudioIO(){
	close();
}


void AudioIO::init(){

	// Choose default devices for now...
	deviceIn(AudioDevice::defaultInput());
	deviceOut(AudioDevice::defaultOutput());
	
//	inDevice(defaultInDevice());
//	outDevice(defaultOutDevice());
//	
//	// Setup input stream parameters
//	const PaDeviceInfo * dInfo = Pa_GetDeviceInfo(mInParams.device);	
//	if(dInfo) mInParams.suggestedLatency = dInfo->defaultLowInputLatency; // for RT
//	mInParams.sampleFormat = paFloat32;// | paNonInterleaved;
//	//mInParams.sampleFormat = paInt16;
//	mInParams.hostApiSpecificStreamInfo = NULL;
//
//	// Setup output stream parameters
//	dInfo = Pa_GetDeviceInfo(mOutParams.device);
//	if(dInfo) mOutParams.suggestedLatency = dInfo->defaultLowOutputLatency; // for RT
//	mOutParams.sampleFormat = paFloat32;// | paNonInterleaved;
//	mOutParams.hostApiSpecificStreamInfo = NULL;

	mImpl->setInDeviceChans(0);
	mImpl->setOutDeviceChans(0);
}

void AudioIO::deviceIn(const AudioDevice& v){

	if(v.valid() && v.hasInput()){
		mImpl->inDevice(v.id());
		const PaDeviceInfo * dInfo = Pa_GetDeviceInfo(mImpl->mInParams.device);	
		if(dInfo) mImpl->mInParams.suggestedLatency = dInfo->defaultLowInputLatency; // for RT
		mImpl->mInParams.sampleFormat = paFloat32;// | paNonInterleaved;
		//mInParams.sampleFormat = paInt16;
		mImpl->mInParams.hostApiSpecificStreamInfo = NULL;
	}
	else{
		warn("attempt to set input device to a device without inputs", "io::AudioIO");
	}
}

void AudioIO::deviceOut(const AudioDevice& v){
	if(v.valid() && v.hasOutput()){
		mImpl->outDevice(v.id());
		const PaDeviceInfo * dInfo = Pa_GetDeviceInfo(mImpl->mOutParams.device);
		if(dInfo) mImpl->mOutParams.suggestedLatency = dInfo->defaultLowOutputLatency; // for RT
		mImpl->mOutParams.sampleFormat = paFloat32;// | paNonInterleaved;
		mImpl->mOutParams.hostApiSpecificStreamInfo = NULL;
	}
	else{
		warn("attempt to set output device to a device without outputs", "io::AudioIO");
	}
}


void AudioIO::channelsBus(int num){
	mNumB = resize(mBufB, num * mFramesPerBuffer);
}


void AudioIO::channels(int num, bool forOutput){
	
	PaStreamParameters * params = forOutput ? &mImpl->mOutParams : &mImpl->mInParams;
	
	if(num == 0){
		//params->device = paNoDevice;
		params->channelCount = 0;
		return;
	}

	const PaDeviceInfo * info = Pa_GetDeviceInfo(params->device);
	if(0 == info){
		if(forOutput)	warn("attempt to set number of channels on invalid output device", "io::AudioIO");
		else			warn("attempt to set number of channels on invalid input device", "io::AudioIO");
		return;	// this particular device is not open, so return
	}


	// compute number of channels to give PortAudio
	int maxChans = 
		(int)(forOutput ? info->maxOutputChannels : info->maxInputChannels);
	
	// -1 means open all channels
	if(num == -1){
		num = maxChans;
	}
	
	int currentNum = channels(forOutput);
	
	if(num != currentNum){

		params->channelCount = min(num, maxChans);

		forOutput ? mNumO = num : mNumI = num;
		
		deferBufferResize(forOutput);
	}
}


bool AudioIO::close(){ return mImpl->close(); }


void AudioIO::deferBufferResize(bool forOutput){
	if(forOutput)	mOutResizeDeferred = true;
	else			mInResizeDeferred = true;
}


bool AudioIO::open(){
	Impl& i = *mImpl;

	i.mErrNum = paNoError;

	if(!(i.mIsOpen || i.mIsRunning)){

		resizeBuffer(false);
		resizeBuffer(true);

		resize(mBufT, mFramesPerBuffer);
		
		PaStreamParameters * inParams = &i.mInParams;
		PaStreamParameters * outParams = &i.mOutParams;
		
		// Must pass in 0s for input- or output-only streams.
		// Stream will not be opened if no device or channel count is zero
		if((paNoDevice ==  inParams->device) || (0 ==  inParams->channelCount)) inParams  = 0;
		if((paNoDevice == outParams->device) || (0 == outParams->channelCount)) outParams = 0;

		i.mErrNum = Pa_OpenStream(
			&i.mStream,			// PortAudioStream **
			inParams,			// PaStreamParameters * in
			outParams,			// PaStreamParameters * out
			mFramesPerSecond,	// frames/sec (double)
			mFramesPerBuffer,	// frames/buffer (unsigned long)
            paNoFlag,			// paNoFlag, paClipOff, paDitherOff
			paCallback,			// static callback function (PaStreamCallback *)
			this
		);

		i.mIsOpen = paNoError == i.mErrNum;
	}
	//printf("AudioIO::open()\n"); printError();
	return paNoError == i.mErrNum;
}


int paCallback(
	const void *input,
	void *output,
	unsigned long frameCount,
	const PaStreamCallbackTimeInfo* timeInfo,
	PaStreamCallbackFlags statusFlags,
	void * userData
){
	AudioIO& io = *(AudioIO *)userData;

	const float * paI = (const float *)input;
	float * paO = (float *)output;

	bool bDeinterleave = true;

	if(bDeinterleave){
		deinterleave(const_cast<float *>(&io.in(0,0)),  paI, io.framesPerBuffer(), io.channelsInDevice() );
		//deinterleave(&io.out(0,0), paO, io.framesPerBuffer(), io.channelsOutDevice());
	}
	
	if(io.autoZeroOut()) io.zeroOut();

	io();	// call callback

	// kill pesky nans so we don't hurt anyone's ears
	if(io.zeroNANs()){
		for(int i=0; i<io.framesPerBuffer()*io.channelsOutDevice(); ++i){
			float& s = (&io.out(0,0))[i];
			if(isnan(s)) s = 0.f;
		}
	}
	
	if(io.clipOut()){
		for(int i=0; i<io.framesPerBuffer()*io.channelsOutDevice(); ++i){
			float& s = (&io.out(0,0))[i];
			if		(s<-1.f) s =-1.f;
			else if	(s> 1.f) s = 1.f;
		}		
	}

	if(bDeinterleave){
		interleave(paO, &io.out(0,0), io.framesPerBuffer(), io.channelsOutDevice());
	}

	return 0;
}


void AudioIO::reopen(){
	if(mImpl->mIsRunning)  { close(); start(); }
	else if(mImpl->mIsOpen){ close(); open(); }
}

void AudioIO::resizeBuffer(bool forOutput){

	float *& buffer = forOutput ? mBufO : mBufI;
	int& chans      = forOutput ? mNumO : mNumI;
	bool& deferred  = forOutput ? mOutResizeDeferred : mInResizeDeferred;

	if(deferred){
		if(chans > 0){			
			int n = resize(buffer, chans * mFramesPerBuffer);
			if(n){	deferred = false; }
			else{	chans = 0; }
		}
		else{
			deleteBuf(buffer);
			deferred = false;
		}
	}
}


void AudioIO::framesPerSecond(double v){	//printf("AudioIO::fps(%f)\n", v);
	if(AudioIOData::framesPerSecond() != v){
                
		if(!supportsFPS(v)) v = mOutDevice.defaultSampleRate();

		mFramesPerSecond = v;
		reopen();
	}
}


void AudioIO::framesPerBuffer(int n){
	if(framesPerBuffer() != n){
		mFramesPerBuffer = n;
		channelsBus(AudioIOData::channelsBus());
		reopen();
	}
}


bool AudioIO::start(){
	Impl& i = *mImpl;
	i.mErrNum = paNoError;
	if(!i.mIsOpen) open();
	if(i.mIsOpen && !i.mIsRunning)	i.mErrNum = Pa_StartStream(i.mStream);
	if(paNoError == i.mErrNum)	mImpl->mIsRunning = true;
	return paNoError == i.mErrNum;
}

bool AudioIO::stop(){ return mImpl->stop(); }

bool AudioIO::supportsFPS(double fps) const { return mImpl->supportsFPS(fps); }

void AudioIO::print(){
	if(mInDevice.id() == mOutDevice.id()){
		printf("I/O Device:  "); mInDevice.print();
	}
	else{
		printf("Device In:   "); mInDevice.print();
		printf("Device Out:  "); mOutDevice.print();
	}

		printf("Chans In:    %d (%dD + %dV)\n", channelsIn(), channelsInDevice(), channelsIn() - channelsInDevice());
		printf("Chans Out:   %d (%dD + %dV)\n", channelsOut(), channelsOutDevice(), channelsOut() - channelsOutDevice());

	const PaStreamInfo * sInfo = Pa_GetStreamInfo(mImpl->mStream);
	if(sInfo){
		printf("In Latency:  %.0f ms\nOut Latency: %0.f ms\nSample Rate: %0.f Hz\n",
			sInfo->inputLatency * 1000., sInfo->outputLatency * 1000., sInfo->sampleRate);
	}
	printf("Frames/Buf:  %d\n", mFramesPerBuffer);
}


void AudioIO::operator()(){ frame(0); if(callback) callback(*this); }

int AudioIO::channels(bool forOutput) const { return forOutput ? channelsOut() : channelsIn(); }
double AudioIO::cpu() const { return Pa_GetStreamCpuLoad(mImpl->mStream); }
bool AudioIO::zeroNANs() const { return mZeroNANs; }


} // gam::
