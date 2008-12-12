#ifndef GAMMA_FIL_H_INC
#define GAMMA_FIL_H_INC

/*	Synz (synthesize) - Signal processing library
	See COPYRIGHT file for authors and license information */

#include "scl.h"

namespace gam{
namespace fil{



/// One element delay.
template<class T = gam::real> 
class Delay1{
public:
	/// @param[in] iprev	Initial previous value.
	Delay1(T iprev = 0) : prev(iprev){}

	T prev;				///< The previous input sample.

	/// Returns next delayed sample.
	T operator()(T v){ T p=prev; prev=v; return p; }
};



/// Two element delay.
template <class T=gam::real>
class Delay2{
public:
	Delay2(const T& v2=(T)0, const T& v1=(T)0, const T& v0=(T)0) : v0(v0), v1(v1), v2(v2){}
	T operator()(const T& v){ v2=v1; v1=v0; v0=v; return v0; }
	T operator()(){ return v0; }		///< Get current value
	T d1f(){ return v0-v1; }			///< Get first forward difference
	T d1b(){ return v1-v2; }			///< Get first backward difference
	T d1c(){ return (v0-v2)*0.5; }		///< Get first center difference
	T d2(){ return d1f() - d1b(); }		///< Get second difference

private: T v0, v1, v2;
};


// Direct domain lowpass filter
template <class T=gam::real>
class Lowpass{
public:
	Lowpass(const T& cen=T(0), const T& wid=T(0.1), const T& peak=T(0))
	:	mCenter(cen), mPeak(peak)
	{ width(wid); }

	T operator()(const T& x) const {
		const T c = x - mCenter;
		return num(c)/den(c)*0.5+0.5;
	}
	
//	T den(const T& x) const { return mWidth+scl::abs(x); }
//	T num(const T& x) const { return mPeak*mWidth-x; }

	T den(const T& x) const { return mWidth+x*x; }
	T num(const T& x) const { return mPeak*mWidth-x*scl::abs(x); }

	void width(const T& v){
		static const T eps = 0.000001;	// width must be greater than zero
		mWidth = scl::abs(v);
		if(mWidth < eps) mWidth = eps;
	}

	void set(const T& cen, const T& wid, const T& peak){
		mCenter=cen; width(wid); mPeak=peak;
	}

protected:
	T mCenter, mWidth, mPeak;
};


// Direct domain inverse comb filter
template <class T=gam::real>
class CombI{
public:

	CombI(const T& frq=T(1), const T& resn=T(0.9), const T& phs=T(0))
	:	mFreq(frq), mPhase(phs)
	{ res(resn); }

	T operator()(const T& x) const { return num(x)/den(x); }
	T den(const T& x) const { return T(1) + mRes*(mRes - T(2)*cos(M_2PI*(mFreq*x + mPhase))); }
	T num(const T& x) const { return mNum; }

	void res(const T& v){
		static const T eps = 0.999999;	// width must be greater than zero
		mRes = scl::abs(v);
		if(mRes>eps) mRes=eps;
		mNum = T(1) + mRes*(mRes - T(2));
	}
	
	void set(const T& frq, const T& resn, const T& phs){
		mFreq=frq; res(resn); mPhase=phs;
	}
	
protected:
	T mFreq, mRes, mPhase;
	T mNum;
};

} // fil::
} // gam::

#endif
