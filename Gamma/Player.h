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


/// Sample player

///	The number of frames in the sample should not exceed 2^32.  This equates
///	to 27 hours at 44.1 kHz.

/*
			Max sample duration
phasor		44.1					96						192
32-bit		1623.19 m / 27.05 h		745.65 m / 12.43 h		372.83 m / 6.21 h
64-bit		6.97e12 m / 1.16e11 h	3.2e12 m / 5.33e10 h	1.6e12 m / 2.67e10 h
*/

template <class T=gam::real, template<class> class Tipol=ipl::Trunc, class Ttap=tap::Clip>
class Player: public Synced, public Array<T>{
public:
	using Array<T>::size; using Array<T>::elems;

	/// @param[in] path			Path to sound file
	/// @param[in] rate			Playback rate scalar
	Player(const char * path, double rate=1);

	/// @param[in] src			Another sampler to read data from
	/// @param[in] rate			Playback rate scalar	
	Player(const Player<T>& src, double rate=1)
	:	Array<T>(src), mSampleRate(src.sampleRate()), mChans(src.channels()),
		mPos(0), mInc(1), mRate(rate), mMin(0), mMax(src.size())
	{ initSynced(); }
	
	Player(const Array<T>& src, double smpRate, double rate)
	:	Array<T>(src), mSampleRate(smpRate), mChans(1),
		mPos(0), mInc(1), mRate(rate), mMin(0), mMax(src.size())
	{
		initSynced();
		sampleRate(smpRate);
	}

	/// Increment read tap
	void advance(){
		mPos = mTap(pos(), mInc, max(), min()); // update read tap, in frames
	}

	/// Returns sample at current position on specified channel and increments phase
	T operator()(int channel=0){ T r = read(channel); advance(); return r; }

	/// Returns sample at current position on specified channel (without incrementing phase)
	T read(int channel) const {
		uint32_t posi = (uint32_t)pos();
		int Nframes= frames();
		int offset = channel*Nframes;
		return mIpol(*this, posi+offset, pos()-posi, offset+Nframes-1, offset);
	}

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
	static T sDummyElement;
	static Array<T> sDummyArray;

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

PRE T CLS::sDummyElement = (T)0;
PRE Array<T> CLS::sDummyArray(&sDummyElement, 1);

PRE CLS::Player(const char * path, double rate)
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
		sDummyArray.source(&sDummyElement, 1);	// make sure this gets set
		source(sDummyArray);
	}
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

