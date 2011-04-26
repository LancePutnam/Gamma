#ifndef GAMMA_PLAYER_H_INC
#define GAMMA_PLAYER_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include "Gamma/Containers.h"
#include "Gamma/ipl.h"
#include "Gamma/scl.h"
#include "Gamma/SoundFile.h"
#include "Gamma/Strategy.h"
#include "Gamma/Sync.h"

namespace gam{


/// Sample buffer player

///	The number of frames in the sample should not exceed 2^32.  This equates
///	to 27 hours at 44.1 kHz.
template <class T=gam::real, template<class> class Tipol=ipl::Trunc, class Ttap=tap::Clip>
class Player: public Synced, public Array<T>{
public:
	using Array<T>::size; using Array<T>::elems;

	Player()
	:	Array<T>(defaultArray()), mSampleRate(1), mChans(1),
		mPos(0), mInc(0), mRate(0), mMin(0), mMax(1)
	{}


	/// @param[in] src		Another Player to read data from
	/// @param[in] rate		Playback rate
	explicit Player(const Player<T>& src, double rate=1)
	:	Array<T>(src), mSampleRate(src.sampleRate()), mChans(src.channels()),
		mPos(0), mInc(1), mRate(rate), mMin(0), mMax(src.size())
	{ initSynced(); }


	/// @param[in] src		Sample array to reference
	/// @param[in] smpRate	Sample rate of samples
	/// @param[in] rate		Playback rate
	Player(const Array<T>& src, double smpRate, double rate=1)
	:	Array<T>(src), mSampleRate(smpRate), mChans(1),
		mPos(0), mInc(1), mRate(rate), mMin(0), mMax(src.size())
	{
		initSynced();
		sampleRate(smpRate);
	}


	/// @param[in] path		Path to sound file
	/// @param[in] rate		Playback rate
	template<class Char>
	explicit Player(const Char * path, double rate=1);


	/// Increment read tap
	void advance(){
		mPos = mTap(pos(), mInc, max(), min()); // update read tap, in frames
	}

	/// Returns sample at current position on specified channel and increments phase
	T operator()(int channel=0){ T r = read(channel); advance(); return r; }

	/// Returns sample at current position on specified channel (without incrementing phase)
	T read(int channel) const {
		uint32_t posi = uint32_t(pos());
		int Nframes= frames();
		int offset = channel*Nframes;
		return mIpol(*this, posi+offset, pos()-posi, offset+Nframes-1, offset);
	}

	/// Set sample buffer
	
	/// @param[in] src		Sample buffer (if multichannel, must be deinterleaved)
	/// @param[in] smpRate	Sample rate of samples
	/// @param[in] channels	Number of channels in sample buffer
	void buffer(const Array<T>& src, double smpRate, int channels);

	void free();							///< Free sample buffer (if owner)
	void max(double v);						///< Set interval max, in frames
	void min(double v);						///< Set interval min, in frames
	void pos(double v);						///< Set current read position, in frames
	void phase(double v);					///< Set current read position [0, 1)
	void rate(double v);					///< Set playback rate scalar
	void range(double phs, double period);	///< Set interval start phase and period
	void reset();							///< Reset playback head

	double max() const { return mMax; }		///< Get interval max
	double min() const { return mMin; }		///< Get interval min
	double period() const;					///< Get total period of sample data
	double pos() const { return mPos; }		///< Get current read position, in frames
	double posInInterval(double frac) const;///< Get position from fraction within interval
	double rate() const { return mRate; }	///< Get playback rate
	double sampleRate() const { return mSampleRate; } ///< Get sample rate of sample buffer

	int channels() const { return mChans; }	///< Get number of channels

	virtual void onResync(double r){ sampleRate(mSampleRate); }

protected:
	static Array<T>& defaultArray(){
		static T v = T(0);
		static Array<T> a(&v, 1);
		return a;
	}

	Tipol<T> mIpol;
	Ttap mTap;

	double mPos, mInc;			// real index position and increment
	double mSampleRate;			// sample rate of array data
	int mChans;					// number of channels
	double mRate, mMin, mMax;
	
	void sampleRate(double v){
		mSampleRate = v;
		rate(mRate);
	}

	int frames() const { return size()/channels(); }
};

#define PRE template <class T, template<class> class Ti, class Tt>
#define CLS Player<T,Ti,Tt>

PRE
template<class Char>
CLS::Player(const Char * path, double rate)
:	Array<T>(), mPos(0), mInc(1), mChans(1), mRate(rate), mMin(0), mMax(1)
{
	SoundFile sf(path);
	
	if(sf.openRead()){
		Array<T>::resize(sf.samples());
		sf.readAllD(elems());
		sampleRate(sf.frameRate());
		mChans = sf.channels();
		mMax = frames();
		sf.close();
	}
	else{
		source(defaultArray());
	}
}

PRE void CLS::buffer(const Array<T>& src, double smpRate, int channels){
	source(src);
	sampleRate(smpRate);	// sets mSampleRate, mRate, and mInc
	mChans = channels;
	mMin = 0;
	mMax = frames();
	mPos = 0;
}

PRE inline void CLS::pos(double v){	mPos = v; }
PRE inline void CLS::phase(double v){ pos(v * frames()); }
PRE inline void CLS::min(double v){	mMin = scl::clip<double>(v, frames()); }	
PRE inline void CLS::max(double v){ mMax = scl::clip<double>(v, frames()); }

PRE void CLS::free(){ this->freeElements(); }
PRE inline void CLS::rate(double v){ mRate = v; mInc = v * scaleSPU(); }
PRE inline void CLS::range(double posn, double period){
	phase(posn);
	min(pos());
	max(pos() + period * spu());	
}

PRE inline void CLS::reset(){
	pos(rate()<0 ? max() : min());
	mTap.reset();
}

PRE inline double CLS::period() const { return frames() * ups(); }
PRE inline double CLS::posInInterval(double frac) const { return min() + (max() - min()) * frac; }

#undef PRE
#undef CLS

} // gam::

#endif

