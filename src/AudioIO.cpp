/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring> // memset

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
static void deleteBuf(T *& buf){ delete[] buf; buf=nullptr; }

template <class T>
static int resize(T *& buf, int n){
	deleteBuf(buf);
	buf = new T[n];
	return n;
}

template <class T>
static inline void zero(T * buf, int n){ std::memset(buf, 0, n*sizeof(T)); }

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


void AudioDevice::initDevices(){
	struct InitSingleton{
		InitSingleton(): mCleanUp(paNoError == Pa_Initialize()){}
		~InitSingleton(){ if(mCleanUp){ Pa_Terminate(); } }
		bool mCleanUp;
	};
	static InitSingleton dummy;
}

AudioDevice::AudioDevice(int deviceNum){
	setImpl(deviceNum);
}

AudioDevice::AudioDevice(const std::string& nameKeyword, StreamMode stream){
	for(int i=0; i<numDevices(); ++i){
		AudioDevice d(i);
		bool bi = (stream &  INPUT) && d.hasInput();
		bool bo = (stream & OUTPUT) && d.hasOutput();
		std::string n = d.name();

		if(	(bi || bo) && n.find(nameKeyword) != std::string::npos){
			setImpl(i);
			break;
		}
	}
}

AudioDevice::~AudioDevice(){}

bool AudioDevice::valid() const { return 0!=mImpl; }
int AudioDevice::id() const { return mID; }
const char * AudioDevice::name() const { return mName; }
int AudioDevice::channelsInMax() const { return mChanIMax; }
int AudioDevice::channelsOutMax() const { return mChanOMax; }
double AudioDevice::defaultSampleRate() const { return mDefSampleRate; }
bool AudioDevice::hasInput() const { return channelsInMax()>0; }
bool AudioDevice::hasOutput() const { return channelsOutMax()>0; }

void AudioDevice::setImpl(int deviceNum){
	initDevices();
	auto * info = Pa_GetDeviceInfo(deviceNum);
	mName = info->name;
	mChanIMax = info->maxInputChannels;
	mChanOMax = info->maxOutputChannels;
	mDefSampleRate = info->defaultSampleRate;
	mImpl = info;
	mID = deviceNum;
}

/*static*/ AudioDevice AudioDevice::defaultInput(){
	initDevices();
	return AudioDevice(Pa_GetDefaultInputDevice());
}

/*static*/ AudioDevice AudioDevice::defaultOutput(){
	initDevices();
	return AudioDevice(Pa_GetDefaultOutputDevice());
}

/*static*/ int AudioDevice::numDevices(){
	initDevices();
	return Pa_GetDeviceCount();
}

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
	/*
	PaSampleFormat sampleFormats = info->nativeSampleFormats;
	
	printf("[ ");
	if(0 != sampleFormats & paFloat32)		printf("f32 ");
	if(0 != sampleFormats & paInt32)		printf("i32 ");
	if(0 != sampleFormats & paInt24)		printf("i24 ");
	if(0 != sampleFormats & paInt16)		printf("i16 ");
	if(0 != sampleFormats & paInt8)			printf("i8 ");
	if(0 != sampleFormats & paUInt8)		printf("ui8 ");
	printf("], ");
	
	if(info->numSampleRates != -1){
		printf("[");
		for(int i=0; i<info->numSampleRates; i++){
			printf("%f ", info->sampleRates[i]);
		}
		printf("] Hz");
	}
	else{
		printf("[%.0f <-> %.0f] Hz", info->sampleRates[0], info->sampleRates[1]);
	}
	printf("\n");
	*/
}

/*static*/ void AudioDevice::printAll(){
	for(int i=0; i<numDevices(); i++){
		printf("[%2d] ", i);
		AudioDevice dev(i);
		dev.print();
		//print(i);
	}
}


//==============================================================================
class AudioIOData::Impl{
public:
	bool error() const { return mErrNum != paNoError; }

	void printError(const char * text = "") const {
		if(error()){
			fprintf(stderr, "%s: %s\n", text, Pa_GetErrorText(mErrNum));
		}
	}

	bool supportsFPS(double fps) const {
		const PaStreamParameters * pi = mInParams.channelCount  == 0 ? 0 : &mInParams;
		const PaStreamParameters * po = mOutParams.channelCount == 0 ? 0 : &mOutParams;	
		mErrNum = Pa_IsFormatSupported(pi, po, fps);
		printError("AudioIO::Impl::supportsFPS");
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
	PaStream * mStream = nullptr;		// i/o stream
	mutable PaError mErrNum = 0;		// Most recent error number
	bool mIsOpen = false;				// An audio device is open
	bool mIsRunning = false;			// An audio stream is running
};

AudioIOData::AudioIOData(void * userData)
:	mImpl(new Impl),
	mUser(userData)
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
AudioIO::AudioIO(
	int framesPerBuf, double framesPerSec, void (* callbackA)(AudioIOData &), void * userData,
	int outChansA, int inChansA)
:	AudioIOData(userData),
	callback(callbackA),
	mInDevice(AudioDevice::defaultInput()), mOutDevice(AudioDevice::defaultOutput())
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
	/*
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
	*/
	mImpl->setInDeviceChans(0);
	mImpl->setOutDeviceChans(0);
}

AudioIO& AudioIO::append(AudioCallback& v){
	mAudioCallbacks.push_back(&v);
	return *this;
}

AudioIO& AudioIO::prepend(AudioCallback& v){
	mAudioCallbacks.insert(mAudioCallbacks.begin(), &v);
	return *this;
}

AudioIO& AudioIO::remove(AudioCallback& v){
	// the proper way to do it:
	mAudioCallbacks.erase(std::remove(mAudioCallbacks.begin(), mAudioCallbacks.end(), &v), mAudioCallbacks.end());
	return *this;
}

void AudioIO::deviceIn(const AudioDevice& v){

	if(v.valid() && v.hasInput()){
		//printf("deviceIn: %s, %d\n", v.name(), v.id());
		mInDevice = v;
		mImpl->inDevice(v.id());
		const PaDeviceInfo * dInfo = Pa_GetDeviceInfo(mImpl->mInParams.device);	
		if(dInfo) mImpl->mInParams.suggestedLatency = dInfo->defaultLowInputLatency; // for RT
		mImpl->mInParams.sampleFormat = paFloat32;// | paNonInterleaved;
		//mInParams.sampleFormat = paInt16;
		mImpl->mInParams.hostApiSpecificStreamInfo = NULL;
	}
	else{
		warn("attempt to set input device to a device without inputs", "AudioIO");
	}
}

void AudioIO::deviceOut(const AudioDevice& v){
	if(v.valid() && v.hasOutput()){
		mOutDevice = v;
		mImpl->outDevice(v.id());
		const PaDeviceInfo * dInfo = Pa_GetDeviceInfo(mImpl->mOutParams.device);
		if(dInfo) mImpl->mOutParams.suggestedLatency = dInfo->defaultLowOutputLatency; // for RT
		mImpl->mOutParams.sampleFormat = paFloat32;// | paNonInterleaved;
		mImpl->mOutParams.hostApiSpecificStreamInfo = NULL;
	}
	else{
		warn("attempt to set output device to a device without outputs", "AudioIO");
	}
}

void AudioIO::device(const AudioDevice& v){
	deviceIn(v); deviceOut(v);
}


void AudioIO::channelsBus(int num){

	if(mImpl->mIsOpen){
		warn("the number of channels cannnot be set with the stream open", "AudioIO");
		return;
	}

	resize(mBufB, num * mFramesPerBuffer);
	mNumB = num;
}


void AudioIO::channels(int num, bool forOutput){

	if(mImpl->mIsOpen){
		warn("the number of channels cannnot be set with the stream open", "AudioIO");
		return;
	}

	PaStreamParameters * params = forOutput ? &mImpl->mOutParams : &mImpl->mInParams;
	
	if(num == 0){
		//params->device = paNoDevice;
		params->channelCount = 0;
		return;
	}

	const PaDeviceInfo * info = Pa_GetDeviceInfo(params->device);
	if(0 == info){
		if(forOutput)	warn("attempt to set number of channels on invalid output device", "AudioIO");
		else			warn("attempt to set number of channels on invalid input device", "AudioIO");
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
		params->channelCount = std::min(num, maxChans);
		forOutput ? mNumO = num : mNumI = num;
		resizeBuffer(forOutput);
	}
}


bool AudioIO::close(){ return mImpl->close(); }


static int paCallback(
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


	io.processAudio();	// call callback


	// apply smoothly-ramped gain to all output channels
	if(io.usingGain()){
	
		float dgain = (io.mGain-io.mGainPrev) / io.framesPerBuffer();
	
		for(int j=0; j<io.channelsOutDevice(); ++j){
			float * out = io.outBuffer(j);
			float gain = io.mGainPrev;
			
			for(int i=0; i<io.framesPerBuffer(); ++i){
				out[i] *= gain;
				gain += dgain;
			}
		}
		
		io.mGainPrev = io.mGain;
	}

	// kill pesky nans so we don't hurt anyone's ears
	if(io.zeroNANs()){
		for(int i=0; i<io.framesPerBuffer()*io.channelsOutDevice(); ++i){
			float& s = (&io.out(0,0))[i];
			//if(isnan(s)) s = 0.f;
			if(s != s) s = 0.f; // portable isnan; only nans do not equal themselves
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

bool AudioIO::open(){
	Impl& i = *mImpl;

	i.mErrNum = paNoError;

	if(!(i.mIsOpen || i.mIsRunning)){
		
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

	i.printError("Error in AudioIO::open()");
	return paNoError == i.mErrNum;
}

void AudioIO::reopen(){
	if(mImpl->mIsRunning)  { close(); start(); }
	else if(mImpl->mIsOpen){ close(); open(); }
}

void AudioIO::resizeBuffer(bool forOutput){
	float *& buffer = forOutput ? mBufO : mBufI;
	int& chans      = forOutput ? mNumO : mNumI;

	if(chans > 0){			
		int n = resize(buffer, chans * mFramesPerBuffer);
		if(0 == n) chans = 0;
	}
	else{
		deleteBuf(buffer);
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
	if(mImpl->mIsOpen){
		warn("the number of frames/buffer cannnot be set with the stream open", "AudioIO");
		return;
	}

	if(framesPerBuffer() != n){
		mFramesPerBuffer = n;
		resizeBuffer(true);
		resizeBuffer(false);
		channelsBus(AudioIOData::channelsBus());
		resize(mBufT, mFramesPerBuffer);
	}
}


bool AudioIO::start(){
	Impl& i = *mImpl;
	i.mErrNum = paNoError;
	if(!i.mIsOpen) open();
	if(i.mIsOpen && !i.mIsRunning)	i.mErrNum = Pa_StartStream(i.mStream);
	if(paNoError == i.mErrNum)	mImpl->mIsRunning = true;
	i.printError("Error in AudioIO::start()");
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


//void AudioIO::processAudio(){ frame(0); if(callback) callback(*this); }
void AudioIO::processAudio(){ 
	frame(0); 
	if(callback) callback(*this); 
	
	std::vector<AudioCallback *>::iterator iter = mAudioCallbacks.begin(); 
	while(iter != mAudioCallbacks.end()){
		frame(0); 
		(*iter++)->onAudio(*this);
	}
}

int AudioIO::channels(bool forOutput) const { return forOutput ? channelsOut() : channelsIn(); }
double AudioIO::cpu() const { return Pa_GetStreamCpuLoad(mImpl->mStream); }
bool AudioIO::zeroNANs() const { return mZeroNANs; }

} // gam::
