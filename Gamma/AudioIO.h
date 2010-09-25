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
		
		MyStuff& stuff = *(MyStuff *)io.user();

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

#include <string>

namespace gam{

/// Audio data to be sent to callback
class AudioIOData {
public:

	/// Audio frame iterator
	class Iterator{
	public:
		Iterator(const AudioIOData& io, int frameStart=0)
		:	mFrame(frameStart-1), mIO(io)
		{}
		
		/// Get current frame number
		int frame() const { return mFrame; }

		const AudioIOData& io() const { return mIO; }

		/// Get current input sample on specified channel
		const float& in (int chan) const { return mIO.in (chan)[frame()]; }
		
		/// Get current output sample on specified channel
		float&       out(int chan) const { return mIO.out(chan)[frame()]; }

		/// Get current bus sample on specified channel
		float&       bus(int chan) const { return mIO.bus(chan)[frame()]; }
		
		/// Add value to current output sample on specified channel
		void out(float v, int chan) const { out(chan)+=v; }
		
		/// Add value to current output sample on specified channels
		void out(float v, int ch1, int ch2) const { out(v, ch1); out(v,ch2); }
	
		/// Iterate frame counter, returning true while more frames
		bool operator()() const { return (++mFrame)<mIO.framesPerBuffer(); }
	
	protected:
		mutable int mFrame;
		const AudioIOData& mIO;
	};

	/// Constructor
	AudioIOData(void * user);

	virtual ~AudioIOData();
	
	void * user() const{ return mUser; } ///< Get pointer to user data
	
	float *      bus(int chan) const;	///< Get a bus channel buffer
	const float* in (int chan) const;	///< Get an in channel buffer
	float *      out(int chan) const;	///< Get an out channel buffer
	float *		 temp() const;			///< Get a single channel temporary buffer
	
	int channelsIn () const;			///< Get effective number of input channels
	int channelsOut() const;			///< Get effective number of output channels
	int channelsBus() const;			///< Get number of allocated bus channels

	int channelsInDevice() const;		///< Get number of channels opened on input device
	int channelsOutDevice() const;		///< Get number of channels opened on output device
	int framesPerBuffer() const;		///< Get frames/buffer of audio I/O stream
	double framesPerSecond() const;		///< Get frames/second of audio I/O streams
	double secondsPerBuffer() const;	///< Get seconds/buffer of audio I/O stream
	double time() const;				///< Get current stream time in seconds
	double time(int frame) const;		///< Get current stream time in seconds of frame
	void zeroBus();						///< Zeros all the bus buffers
	void zeroOut();						///< Zeros all the internal output buffers

protected:
	class Impl; Impl * mImpl;
	void * mUser;					// User specified data
	int mFramesPerBuffer;
	double mFramesPerSecond;
	float *mBufI, *mBufO, *mBufB;	// input, output, and aux buffers
	float * mBufT;					// temporary one channel buffer
	int mNumI, mNumO, mNumB;		// input, output, and aux channels
};


/// Audio device abstraction
class AudioDevice{
public:

	/// @param[in] deviceNum	Device enumeration number
	AudioDevice(int deviceNum);
	
	/// @param[in] nameKeyword	Keyword to search for in device name
	/// @param[in] input		Whether to search input devices
	/// @param[in] ouput		Whether to search output devices
	AudioDevice(const std::string& nameKeyword, bool input=true, bool output=true);

	~AudioDevice();

	bool valid() const { return 0!=mImpl; }	///< Returns whether device is valid
	int id() const { return mID; }			///< Get device unique ID
	const char * name() const;				///< Get device name
	int channelsInMax() const;				///< Get maximum number of input channels supported
	int channelsOutMax() const;				///< Get maximum number of output channels supported
	double defaultSampleRate() const;		///< Get default sample rate
	
	bool hasInput() const;					///< Returns whether device has input
	bool hasOutput() const;					///< Returns whether device has output
	
	void print() const;						///< Prints info about specific i/o device to stdout

	static AudioDevice defaultInput();		///< Get system's default input device
	static AudioDevice defaultOutput();		///< Get system's default output device
	static int numDevices();				///< Returns number of audio i/o devices available
	static void printAll();					///< Prints info about all available i/o devices to stdout

private:
	void setImpl(int deviceNum);
	static void initDevices();
	int mID;
	const void * mImpl;
};


/// Audio callback type
typedef void (*audioCallback)(AudioIOData& io);


/// Audio input/output streaming.

/// This is a wrapper around the PortAudio v1.9 library.
/// 
class AudioIO : public AudioIOData {
public:

	/// Creates AudioIO using default I/O devices.
	///
	/// @param[in] framesPerBuf		Number of sample frames to process per callback
	/// @param[in] framesPerSec		Frame rate.  Unsupported values will use default rate of device.
	/// @param[in] callback			Audio processing callback
	/// @param[in] userData			Pointer to user data accessible within callback
	/// @param[in] outChans			Number of output channels to open
	/// @param[in] inChans			Number of input channels to open
	/// If the number of input or output channels is greater than the device
	/// supports, virtual buffers will be created.
	AudioIO(
		int framesPerBuf=64, double framesPerSec=44100.0,
		void (* callback)(AudioIOData &) = 0, void * userData = 0,
		int outChans = 2, int inChans = 0 );

	virtual ~AudioIO();

	using AudioIOData::channelsIn;
	using AudioIOData::channelsOut;
	using AudioIOData::framesPerBuffer;
	using AudioIOData::framesPerSecond;

	audioCallback callback;						///< User specified callback function.

	bool autoZeroOut() const { return mAutoZeroOut; }
	int channels(bool forOutput) const;
	bool clipOut() const { return mClipOut; }	///< Returns clipOut setting
	double cpu() const;							///< Returns current CPU usage of audio thread
	bool supportsFPS(double fps) const;			///< Return true if fps supported, otherwise false
	bool zeroNANs() const;						///< Returns whether to zero NANs in output buffer going to DAC
	
	void operator()();							///< Call callback manually
	bool open();								///< Opens audio device.
	bool close();								///< Closes audio device. Will stop active IO.
	bool start();								///< Starts the audio IO.  Will open audio device if necessary.
	bool stop();								///< Stops the audio IO.

	void autoZeroOut(bool v){ mAutoZeroOut=v; }

	/// Sets number of effective channels on input or output device depending on 'forOutput' flag.
	
	/// An effective channel is either a real device channel or virtual channel 
	/// depending on how many channels the device supports. Passing in -1 for
	/// the number of channels opens all available channels.
	void channels(int num, bool forOutput);

	void channelsIn(int n){channels(n,false);}	///< Set number of input channels
	void channelsOut(int n){channels(n,true);}	///< Set number of output channels
	void channelsBus(int num);					///< Set number of bus channels
	void clipOut(bool v){ mClipOut=v; }			///< Set whether to clip output between -1 and 1
	void deviceIn(const AudioDevice& v);		///< Set input device
	void deviceOut(const AudioDevice& v);		///< Set output device
	void framesPerSecond(double v);				///< Set number of frames per second
	void framesPerBuffer(int n);				///< Set number of frames per processing buffer
	void zeroNANs(bool v){ mZeroNANs=v; }		///< Set whether to zero NANs in output buffer going to DAC

	void print();								///< Prints info about current i/o devices to stdout.

	static const char * errorText(int errNum);		// Returns error string.

private:
	AudioDevice mInDevice, mOutDevice;

	bool mInResizeDeferred, mOutResizeDeferred;
	bool mZeroNANs;			// whether to zero NANs
	bool mClipOut;			// whether to clip output between -1 and 1
	bool mAutoZeroOut;		// whether to automatically zero output buffers each block

	void init();		// Initializes PortAudio and member variables.
	void deferBufferResize(bool forOutput);
	void reopen();		// reopen stream (restarts stream if needed)
	void resizeBuffer(bool forOutput);
};

} // gam::

#endif
