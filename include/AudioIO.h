#ifndef GAMMA_AUDIOIO_H_INC
#define GAMMA_AUDIOIO_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include <string.h>		/* memset() */
#include <portaudio.h>
#include "mem.h"

#define ULONG unsigned long

namespace gam{


/// Audio data to be sent to callback
class AudioIOData {
public:
	AudioIOData(void * user);
	~AudioIOData();

	void * user;						///< User specified data
	
	float *       aux(ULONG channel);	///< Returns an aux channel buffer.
	const float * in (ULONG channel);	///< Returns an in channel buffer.
	float *       out(ULONG channel);	///< Returns an out channel buffer.
	
	ULONG inChans() const;				///< Returns effective number of input channels
	ULONG outChans() const;				///< Returns effective number of output channels
	ULONG auxChans() const;				///< Returns number of aux channels

	ULONG inDeviceChans() const;		///< Returns number of channels opened on input device
	ULONG outDeviceChans() const;		///< Returns number of channels opened on output device
	ULONG numFrames() const;			///< Returns frames/buffer of audio I/O stream
	double secondsPerBuffer() const;	///< Returns seconds/buffer of audio I/O stream
	double fps() const;					///< Returns frames/second of audio I/O streams
	double spf() const;					///< Returns seconds/frame of audio I/O streams
	double time() const;				///< Returns current stream time in seconds
	double time(ULONG frame) const;		///< Returns current stream time in seconds of frame
	void zeroAux();						///< Zeros all the aux buffers
	void zeroOut();						///< Zeros all the internal output buffers
	
protected:
	PaStreamParameters mInParams, mOutParams;	// Input and output stream parameters.
	PaStream * mStream;
	ULONG mFramesPerBuffer;
	double mFramesPerSecond;
	
	float *mBufI, *mBufO, *mBufA;	// input, output, and aux buffers
	ULONG mNumI, mNumO, mNumA;		// input, output, and aux channels
};



class AudioDevice{
public:
	AudioDevice(int deviceNum);
	~AudioDevice();

	bool valid() const { return 0 != mImpl; }
	int id() const { return mID; }
	const char * name() const;
	int maxInputChannels() const;
	int maxOutputChannels() const;
	double defaultSampleRate() const;
	
	void print() const;	/// Prints info about specific i/o device to stdout.

	static AudioDevice defaultInput();
	static AudioDevice defaultOutput();
	static int numDevices();			///< Returns number of audio i/o devices available.
	static void printAll();				///< Prints info about all available i/o devices to stdout.

private:
	void setImpl(int deviceNum);
	static void initDevices();
	int mID;
	const void * mImpl;
};


/// Audio callback type
typedef void (*audioCallback)(AudioIOData& io);


/// Audio input/output streaming.

/// This is a C++ wrapper around the PortAudio v1.9 library.
/// 
class AudioIO : public AudioIOData {
public:

	/// Creates AudioIO using default I/O devices.
	///
	/// @param[in]	framesPerBuf	Number of sample frames to process per callback
	/// @param[in]	framesPerSec	Frame rate.  Unsupported values will use default rate of device.
	/// @param[in]	callback		Audio processing callback
	/// @param[in]	userData		Pointer to user data accessible within callback
	/// @param[in]	outChans		Number of output channels to open
	/// @param[in]	inChans			Number of input channels to open
	/// If the number of input or output channels is greater than the device
	/// supports, virtual buffers will be created.
	AudioIO(
		ULONG framesPerBuf=64, double framesPerSec=44100.0,
		void (* callback)(AudioIOData &) = 0, void * userData = 0,
		int outChans = 2, int inChans = 0 );

	~AudioIO();
		
	static audioCallback callback;	///< User specified callback function.
	
	void operator()();	// Calls callback manually

	bool open();	///< Opens audio device.
	bool close();	///< Closes audio device. Will stop active IO.
	bool start();	///< Starts the audio IO.  Will open audio device if necessary.
	bool stop();	///< Stops the audio IO.
	
	/// Sets number of effective channels on input or output device depending on 'forOutput' flag.
	
	/// An effective channel is either a real device channel or virtual channel 
	/// depending on how many channels the device supports. Passing in -1 for
	/// the number of channels opens all available channels.
	void chans(int num, bool forOutput);
	void setFPS(double fps);
	void setFramesPerBuffer(ULONG numFrames);
	void auxChans(ULONG num);

	ULONG chans(bool forOutput) const;
	double cpu() const;					///< Returns current CPU usage of audio thread
	bool killNANs() const;
	bool supportsFPS(double fps);

	void print();				///< Prints info about current i/o devices to stdout.
	void printError();			///< Prints info about current error status to stdout.

	static const char * errorText(int errNum);		// Returns error string.
	
private:
	PaError mErrNum;							// Most recent error number.
	//PaDeviceIndex mInDevice, mOutDevice;		// Input and output device ids.
	AudioDevice mInDevice, mOutDevice;

	bool mIsOpen;			// An audio device is open
	bool mIsRunning;		// An audio stream is running
	bool mInResizeDeferred, mOutResizeDeferred;
	bool mKillNANs;			// whether to zero NANs

	void init();		// Initializes PortAudio and member variables.
	
	static int paCallback(	const void *input,
							void *output,
							unsigned long frameCount,
							const PaStreamCallbackTimeInfo* timeInfo,
							PaStreamCallbackFlags statusFlags,
							void *userData );

	PaDeviceIndex defaultInDevice();
	PaDeviceIndex defaultOutDevice();
	
	bool error() const;
	void inDevice(PaDeviceIndex index);		// directly set input device
	void outDevice(PaDeviceIndex index);	// directly set output device
	void setInDeviceChans(ULONG num);			// directly set # device input channels
	void setOutDeviceChans(ULONG num);			// directly set # device output channels
	//void virtualChans(ULONG num, bool forOutput);
	
	void deferBufferResize(bool forOutput);
	void resizeBuffer(bool forOutput);

	void reopen();		// reopen stream (restarts stream if needed)

};


// Implementation_______________________________________________________________

inline void AudioIOData::zeroAux(){ mem::zero(mBufA, numFrames() * mNumA); }
inline void AudioIOData::zeroOut(){ mem::zero(mBufO, outChans() * numFrames()); }

inline float *       AudioIOData::aux(ULONG num){ return mBufA + num * numFrames(); }
inline const float * AudioIOData::in (ULONG chn){ return mBufI + chn * numFrames(); }
inline float *       AudioIOData::out(ULONG chn){ return mBufO + chn * numFrames(); }

inline ULONG AudioIOData:: inChans() const { return mNumI; }
inline ULONG AudioIOData::outChans() const { return mNumO; }
inline ULONG AudioIOData::auxChans() const { return mNumA; }
inline ULONG AudioIOData::inDeviceChans() const { return (ULONG)mInParams.channelCount; }
inline ULONG AudioIOData::outDeviceChans() const { return (ULONG)mOutParams.channelCount; }

inline double AudioIOData::fps() const { return mFramesPerSecond; }
inline double AudioIOData::spf() const { return 1. / fps(); }
inline double AudioIOData::time() const { return Pa_GetStreamTime(mStream); }
inline double AudioIOData::time(ULONG frame) const { return (double)frame * spf() + time(); }
inline ULONG AudioIOData::numFrames() const { return mFramesPerBuffer; }
inline double AudioIOData::secondsPerBuffer() const { return (double)numFrames() * spf(); }

inline void AudioIO::operator()(){ if(callback) callback(*this); }

inline ULONG AudioIO::chans(bool forOutput) const { return forOutput ? outChans() : inChans(); }
inline double AudioIO::cpu() const { return Pa_GetStreamCpuLoad(mStream); }
inline bool AudioIO::killNANs() const { return mKillNANs; }
inline const char * AudioIO::errorText(int errNum){ return Pa_GetErrorText(errNum); }

inline bool AudioIO::error() const { return mErrNum != paNoError; }
inline void AudioIO::inDevice(PaDeviceIndex index){ mInParams.device = index; }
inline void AudioIO::outDevice(PaDeviceIndex index){ mOutParams.device = index; }
inline void AudioIO::setInDeviceChans(ULONG num){ mInParams.channelCount = num; }
inline void AudioIO::setOutDeviceChans(ULONG num){ mOutParams.channelCount = num; }

inline PaDeviceIndex AudioIO::defaultInDevice(){ return Pa_GetDefaultInputDevice(); }
inline PaDeviceIndex AudioIO::defaultOutDevice(){ return Pa_GetDefaultOutputDevice(); }

/* Sample callback
void audioCB(AudioIO * io){
	unsigned long numFrames = io->framesPerBuffer();
	float * output		= io->outBuffer(0);
	float * input		= io->inBuffer(0);
	
	MyApplication * app = (MyApplication *)io->userData;
	
	float * out = output;
	float * in = input;
	
	for(unsigned long i=numFrames; i>0; --i){

	}
}
*/

} // end namespace gam

#undef ULONG

#endif

