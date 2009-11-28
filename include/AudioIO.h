#ifndef GAMMA_AUDIOIO_H_INC
#define GAMMA_AUDIOIO_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */


/*	This is a simple example demonstrating how to set up a callback
	and process input and output buffers.
	
	struct MyStuff{};
	
	void audioCB(AudioIOData& io){
		float * out1 = io.out(0);
		float * out2 = io.out(1);
		const float * in1 = io.in(0);
		const float * in2 = io.in(1);
		
		MyStuff& stuff = *(MyStuff *)io.user;

		for(unsigned i=0; i<io.framesPerBuffer(); ++i){

			float inSample1 = in1[i];
			float inSample2 = in2[i];

			out1[i] = -inSample1;
			out2[i] = -inSample2;
		}
	}
	
	int main(){
		MyStuff stuff;
		
		AudioIO audioIO(128, 44100, audioCB, &stuff, 2,2);
		audioIO.start();
	}
*/


#include <string.h>		/* memset() */
#include <portaudio.h>
#include "mem.h"


namespace gam{


/// Audio data to be sent to callback
class AudioIOData {
public:
	AudioIOData(void * user);
	virtual ~AudioIOData();

	void * user;							///< User specified data
	
	float *       aux(uint32_t channel);	///< Returns an aux channel buffer
	const float * in (uint32_t channel);	///< Returns an in channel buffer
	float *       out(uint32_t channel);	///< Returns an out channel buffer
	float *		  temp();					///< Returns single channel temporary buffer
	
	uint32_t inChans() const;				///< Returns effective number of input channels
	uint32_t outChans() const;				///< Returns effective number of output channels
	uint32_t auxChans() const;				///< Returns number of aux channels

	uint32_t inDeviceChans() const;			///< Returns number of channels opened on input device
	uint32_t outDeviceChans() const;		///< Returns number of channels opened on output device
	uint32_t framesPerBuffer() const;		///< Returns frames/buffer of audio I/O stream
	double framesPerSecond() const;			///< Returns frames/second of audio I/O streams
	double secondsPerBuffer() const;		///< Returns seconds/buffer of audio I/O stream
	double time() const;					///< Returns current stream time in seconds
	double time(uint32_t frame) const;		///< Returns current stream time in seconds of frame
	void zeroAux();							///< Zeros all the aux buffers
	void zeroOut();							///< Zeros all the internal output buffers
	
protected:
	PaStreamParameters mInParams, mOutParams;	// Input and output stream parameters.
	PaStream * mStream;
	uint32_t mFramesPerBuffer;
	double mFramesPerSecond;
	
	float *mBufI, *mBufO, *mBufA;		// input, output, and aux buffers
	float * mBufT;						// temporary one channel buffer
	uint32_t mNumI, mNumO, mNumA;		// input, output, and aux channels
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
	using AudioIOData::inChans;
	using AudioIOData::outChans;
	using AudioIOData::framesPerBuffer;
	using AudioIOData::framesPerSecond;

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
		uint32_t framesPerBuf=64, double framesPerSec=44100.0,
		void (* callback)(AudioIOData &) = 0, void * userData = 0,
		int outChans = 2, int inChans = 0 );

	virtual ~AudioIO();
		
	static audioCallback callback;	///< User specified callback function.
	
	void operator()();				///< Calls callback manually

	bool open();					///< Opens audio device.
	bool close();					///< Closes audio device. Will stop active IO.
	bool start();					///< Starts the audio IO.  Will open audio device if necessary.
	bool stop();					///< Stops the audio IO.

	/// Sets number of effective channels on input or output device depending on 'forOutput' flag.
	
	/// An effective channel is either a real device channel or virtual channel 
	/// depending on how many channels the device supports. Passing in -1 for
	/// the number of channels opens all available channels.
	void chans(int num, bool forOutput);
	void clipOut(bool v){ mClipOut=v; }			///< Set whether to clip output between -1 and 1
	void framesPerSecond(double v);				///< Set number of frames per second
	void framesPerBuffer(uint32_t n);			///< Set number of frames per processing buffer
	void inChans(uint32_t n){ chans(n,false); }	///< Set number of input channels
	void outChans(uint32_t n){ chans(n,true); }	///< Set number of output channels
	void auxChans(uint32_t num);				///< Set number of auxiliary channels
	void zeroNANs(bool v){ mZeroNANs=v; }		///< Set whether to zero NANs in output buffer

	uint32_t chans(bool forOutput) const;
	bool clipOut() const { return mClipOut; }	///< Returns clipOut setting
	double cpu() const;							///< Returns current CPU usage of audio thread
	bool supportsFPS(double fps);				///< Return true if fps supported, otherwise false
	bool zeroNANs() const;						///< Returns zeroNANs setting

	void print();								///< Prints info about current i/o devices to stdout.
	void printError();							///< Prints info about current error status to stdout.

	static const char * errorText(int errNum);		// Returns error string.
	
private:
	PaError mErrNum;							// Most recent error number.
	//PaDeviceIndex mInDevice, mOutDevice;		// Input and output device ids.
	AudioDevice mInDevice, mOutDevice;

	bool mIsOpen;			// An audio device is open
	bool mIsRunning;		// An audio stream is running
	bool mInResizeDeferred, mOutResizeDeferred;
	bool mZeroNANs;			// whether to zero NANs
	bool mClipOut;			// whether to clip output between -1 and 1

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
	void setInDeviceChans(uint32_t num);			// directly set # device input channels
	void setOutDeviceChans(uint32_t num);			// directly set # device output channels
	//void virtualChans(uint32_t num, bool forOutput);
	
	void deferBufferResize(bool forOutput);
	void resizeBuffer(bool forOutput);

	void reopen();		// reopen stream (restarts stream if needed)

};


// Implementation_______________________________________________________________

inline void AudioIOData::zeroAux(){ mem::zero(mBufA, framesPerBuffer() * mNumA); }
inline void AudioIOData::zeroOut(){ mem::zero(mBufO, outChans() * framesPerBuffer()); }

inline float *       AudioIOData::aux(uint32_t num){ return mBufA + num * framesPerBuffer(); }
inline const float * AudioIOData::in (uint32_t chn){ return mBufI + chn * framesPerBuffer(); }
inline float *       AudioIOData::out(uint32_t chn){ return mBufO + chn * framesPerBuffer(); }
inline float *       AudioIOData::temp(){ return mBufT; }

inline uint32_t AudioIOData:: inChans() const { return mNumI; }
inline uint32_t AudioIOData::outChans() const { return mNumO; }
inline uint32_t AudioIOData::auxChans() const { return mNumA; }
inline uint32_t AudioIOData::inDeviceChans() const { return (uint32_t)mInParams.channelCount; }
inline uint32_t AudioIOData::outDeviceChans() const { return (uint32_t)mOutParams.channelCount; }

inline double AudioIOData::framesPerSecond() const { return mFramesPerSecond; }
inline double AudioIOData::time() const { return Pa_GetStreamTime(mStream); }
inline double AudioIOData::time(uint32_t frame) const { return (double)frame / framesPerSecond() + time(); }
inline uint32_t AudioIOData::framesPerBuffer() const { return mFramesPerBuffer; }
inline double AudioIOData::secondsPerBuffer() const { return (double)framesPerBuffer() / framesPerSecond(); }

inline void AudioIO::operator()(){ if(callback) callback(*this); }

inline uint32_t AudioIO::chans(bool forOutput) const { return forOutput ? outChans() : inChans(); }
inline double AudioIO::cpu() const { return Pa_GetStreamCpuLoad(mStream); }
inline bool AudioIO::zeroNANs() const { return mZeroNANs; }
inline const char * AudioIO::errorText(int errNum){ return Pa_GetErrorText(errNum); }

inline bool AudioIO::error() const { return mErrNum != paNoError; }
inline void AudioIO::inDevice(PaDeviceIndex index){ mInParams.device = index; }
inline void AudioIO::outDevice(PaDeviceIndex index){ mOutParams.device = index; }
inline void AudioIO::setInDeviceChans(uint32_t num){ mInParams.channelCount = num; }
inline void AudioIO::setOutDeviceChans(uint32_t num){ mOutParams.channelCount = num; }

inline PaDeviceIndex AudioIO::defaultInDevice(){ return Pa_GetDefaultInputDevice(); }
inline PaDeviceIndex AudioIO::defaultOutDevice(){ return Pa_GetDefaultOutputDevice(); }


} // gam::

#endif
