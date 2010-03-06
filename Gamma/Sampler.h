#ifndef GAMMA_SAMPLER_H_INC
#define GAMMA_SAMPLER_H_INC

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

template <class T=gam::real, class Tipol=ipl::Trunc, class Ttap=tap::Wrap>
class Player: public Synced, public Array<T>{
public:
	using Array<T>::size; using Array<T>::elems;

	/// @param[in] path			Path to sound file
	/// @param[in] rate			Playback rate scalar
	Player(const char * path, double rate=1);

	/// @param[in] src			Another sampler to read data from
	/// @param[in] rate			Playback rate scalar	
	Player(const Player<T>& src, double rate=1)
	:	Array<T>(src), mSampleRate(src.sampleRate()), 
		mPos(0), mInc(1), mRate(rate), mMin(0), mMax(src.size())
	{ initSynced(); }
	
	Player(const Array<T>& src, double smpRate, double rate)
	:	Array<T>(src), mSampleRate(smpRate), 
		mPos(0), mInc(1), mRate(rate), mMin(0), mMax(src.size())
	{
		initSynced();
		sampleRate(smpRate);
	}

	/// Generate next sample
	T operator()(){
		double p = mPos;
		mPos = mTap(mPos + mInc, max(), min());
		uint32_t i = (uint32_t)p;
		return mIpol(*this, i, p-i, size()-1);
	}

	void free();							///< Free sample buffer (if owner)
	void max(double v);						///< Set max range point (sample)
	void min(double v);						///< Set min range point (sample)
	void pos(double v);						///< Set current read position (sample)
	void phase(double v);					///< Set current read position [0, 1)
	void rate(double v);					///< Set playback rate scalar
	void range(double phs, double period);	///< Set range start phase and period
	
	double max() const;						///< Get max range point
	double min() const;						///< Get min range point
	double period() const;					///< Get total period of sample
	double pos() const;						///< Get current read position (samples)
	double posInRange(double frac) const;	///< Get position from a fraction inside set range
	double rate() const;					///< Get playback rate
	double sampleRate() const;				///< Get sample rate of sample buffer

	virtual void onResync(double r){ sampleRate(mSampleRate); }

protected:
	static T sDummyElement;
	static Array<T> sDummyArray;

	Tipol mIpol;
	Ttap mTap;

	double mPos, mInc;
	double mSampleRate;
	double mRate, mMin, mMax;
	
	void sampleRate(double v){
		mSampleRate = v;
		rate(mRate);
	}

};

#define PRE template <class T, class Tipol, class Ttap>
#define CLS Player<T, Tipol, Ttap>

PRE T CLS::sDummyElement = (T)0;
PRE Array<T> CLS::sDummyArray(&sDummyElement, 1);

PRE CLS::Player(const char * path, double rate)
:	Array<T>(), mPos(0), mInc(1), mRate(rate), mMin(0), mMax(1)
{
	
	SoundFile sf(path);
	
	if(sf.openRead()){
		Array<T>::resize(sf.samples());
		sf.readAllD(elems());
		sampleRate(sf.frameRate());
		mMax = size();
		sf.close();
	}
	else{
		sDummyArray.source(&sDummyElement, 1);	// make sure this gets set
		source(sDummyArray);
	}
}

#define FSIZE (double)size()
PRE inline void CLS::pos(double value){	mPos = value; }
PRE inline void CLS::phase(double value){ pos(value * FSIZE); }
PRE inline void CLS::min(double value){	mMin = scl::clip(value, FSIZE); }	
PRE inline void CLS::max(double value){ mMax = scl::clip(value, FSIZE); }
#undef FSIZE

PRE void CLS::free(){ this->freeElements(); }
PRE inline void CLS::rate(double v){ mRate = v; mInc = v * scaleSPU(); }
PRE inline void CLS::range(double posn, double period){
	phase(posn);
	min(pos());
	max(pos() + period * spu());	
}

PRE inline double CLS::max() const { return mMax; }
PRE inline double CLS::min() const { return mMin; }
PRE inline double CLS::period() const { return size() * ups(); }
PRE inline double CLS::pos() const { return mPos; }
PRE inline double CLS::posInRange(double frac) const { return min() + (max() - min()) * frac; }
PRE inline double CLS::rate() const { return mRate; }
PRE double CLS::sampleRate() const { return mSampleRate; }

#undef PRE
#undef CLS


} // gam::

#endif

