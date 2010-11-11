#ifndef GAMMA_RECORDER_H_INC
#define GAMMA_RECORDER_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include <stdlib.h>
#include <vector>

namespace gam{

/// Sound recorder
class Recorder {
public:

	/// @param[in] chans	number of channels
	/// @param[in] frames	number of (multi-)channel frames
	Recorder(int channels=1, int frames=8192);


	/// Get number of recording channels
	int channels() const { return mChans; }
	
	/// Get number of multi-channel recording frames
	int frames() const { return size()/channels(); }

	/// Get total number of samples (frames x channels) in buffer
	int size() const { return mRing.size(); }

	/// Write sample into ring buffer without advancing write tap
	void overwrite(float v, int chan){
		mRing[mIW+chan] = v;
	}

	/// Write sample into ring buffer and advance write tap
	void write(float v, int chan=0){
		overwrite(v,chan);
		if((mIW+=channels()) >= (int)mRing.size()) mIW=0;		
	}

	/// Write samples into ring buffer and advance write tap
	/// Call this from the audio thread
	void write(float v1, float v2, int chan=0){
		overwrite(v1,chan);
		overwrite(v2,chan+1);
		if((mIW+=channels()) >= (int)mRing.size()) mIW=0;		
	}
	
	/// Empty buffer of most recent samples written from audio thread.

	/// Returns number of frames copied to buffer. If the number of 
	/// frames returned is 0, then no samples were read and 'buf'
	/// is unmodified.
	/// This should be called from a lower priority thread.
	int read(float *& buf);

	/// Resize buffers
	
	/// @param[in] chans	number of channels
	/// @param[in] frames	number of (multi-)channel frames
	void resize(int chans, int frames);

protected:
	int mChans;	// no. of interleaved channels
	int mIW;	// index of next sample to write
	int mIR;	// start index for reading
	std::vector<float> mRing;
	std::vector<float> mRead;
};

} // gam::
#endif
