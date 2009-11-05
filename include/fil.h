#ifndef GAMMA_FIL_H_INC
#define GAMMA_FIL_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include "scl.h"

namespace gam{
namespace fil{



/// One element delay
template<class T = gam::real> 
class Delay1{
public:
	/// @param[in] iprev	Initial previous value.
	Delay1(T iprev = 0): prev(iprev){}

	T prev;				///< The previous input sample.

	/// Returns next delayed sample.
	T operator()(T v){ T p=prev; prev=v; return p; }
};



/// Two element delay
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



/// Returns every nth input n times.
template <class T>
struct Hold{

	Hold(uint32_t n=1): mCount(0), mPeriod(n){}

	T operator()(const T& v){
		if(mCount==0) mVal=v;
		if((++mCount) >= mPeriod) mCount=0;
		return mVal;
	}

	void period(uint32_t v){ mPeriod=v; }
	void reset(){ mCount=0; }

private:
	uint32_t mCount;
	uint32_t mPeriod;
	T mVal;
};



/// Transfer function of an arbitrary difference equation.
class TransferFunc {
public:

	/// @param[in] gain		overall filter gain
	TransferFunc(double gain=1): mGain(gain){}

	/// Add feedforward sample delay
	TransferFunc& addX(double c, double d){ mX.push_back(DelayUnit(c,d)); return *this; }

	/// Add feedback sample delay
	TransferFunc& addY(double c, double d){ mY.push_back(DelayUnit(c,d)); return *this; }
	
	/// Clear sample delays
	TransferFunc& clear(){ mX.clear(); mY.clear(); return *this; }
	
	/// Set overall gain factor
	TransferFunc& gain(double v){ mGain=v; return *this; }
	
	/// Returns frequency response at unit frequency [-0.5, 0.5]
	Complexd operator()(double f){
		Complexd X(0,0), Y(1,0);
		f *= M_2PI;
		for(uint32_t i=0; i<mX.size(); ++i) X += mX[i].response(f);
		for(uint32_t i=0; i<mY.size(); ++i) Y -= mY[i].response(f);
		return X/Y * mGain; // H(z) = X(z)/Y(z)
	}

protected:

	struct DelayUnit{
		DelayUnit(double c_, double d_): c(c_), d(d_){}
		Complexd response(double f){ return Complexd::Polar(c, f*d); }
		double c, d;
	};

	std::vector<DelayUnit> mX, mY;
	double mGain;
};



/// Integrator
template <class T>
class Integrator{

	Integrator(const T& v=T(0)): o1(v){}

	T operator()(const T& i0){ return o1+=i0; }

protected:
	T o1;
};



/// Complex resonator

/// The complex resonator consists of two complex numbers- one is an absolute
/// phase and amplitude and the other is a relative (differential) frequency 
/// and decay/growth factor.
template <class T=double>
class Resonator : public Complex<T>{
public:
	typedef Complex<T> C;
	using C::operator();
	using C::operator=;

	/// @param[in] frq	unit frequency
	/// @param[in] amp	amplitude
	/// @param[in] phs	unit phase
	/// @param[in] dec	unit decay factor
	Resonator(const T& frq=T(0), const T& amp=T(1), const T& phs=T(0), const T& dec=T(1)){
		set(frq, amp, phs);
	}

	/// Generate next value
	C operator()(){ return (*this)*=mFreq; }

	/// Filter input
	C operator()(const C& v){ return (*this) = (*this)*mFreq + v; }
	C operator()(const T& v){ return (*this) = (*this)*mFreq + v; }

	/// Set 60 dB decay interval
	//void decay(const Tv& v){ width(T(2.198806796637603 /* -ln(0.001)/pi */)/v); }
	
	/// Set unit bandwidth
	//void width(const Tv& v){ mDecay=::exp(-M_PI*v); freq(freq()); }
	
	/// Set amount to decay after N iterations
	void decay(const T& target, const T& N=T(1)){
		factor(freq(), (1==N) ? target : ::pow(target, 1./N));
	}

	/// Set recursive multiplication factor (frequency and decay/growth factor)
	void factor(const T& frq, const T& dec=T(1)){
		mFreq.fromPolar(dec, frq*M_2PI);
	}

	/// Set unit frequency
	void freq(const T& v){ factor(v, decay()); }

	/// Set unit frequency, amplitude, and unit phase
	
	/// The phase state will be rewound 1 iteration so the first function call
	/// will return a complex number at the desired phase.
	void set(const T& frq, const T& amp, const T& phs, const T& dec=T(1)){
		this->fromPolar(amp, phs*M_2PI);
		factor(frq, dec);
		(*this) = behind();
	}

	void set(const T& frq, const Complex<T>& phs){
		(*this)(phs.r, phs.i); freq(frq);
	}

	/// Get value one iteration ahead of current state
	C ahead() const { return (*this)*mFreq; }
	
	/// Get value one iteration behind current state
	C behind() const { return (*this)/mFreq; }

	/// Get unit decay
	T decay() const { return mFreq.mag(); }
	
	/// Get unit frequency
	T freq() const { return mFreq.phase()*M_1_2PI; }

protected:
	C mFreq;
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
