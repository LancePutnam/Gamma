#ifndef GAMMA_FIL_H_INC
#define GAMMA_FIL_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include <vector>
#include "Gamma/scl.h"
#include "Gamma/Types.h"

namespace gam{

/// Filter function objects
namespace fil{


/// Fixed N-element delay
template <uint32_t N, class T>
class FixedDelay{
public:
	#define DO for(uint32_t i=0; i<N; ++i)

	/// @param[in] v	Initial value of elements
	FixedDelay(const T& v=T()){ DO mVals[i]=v; }

	/// Set nth delayed element
	T& operator[](uint32_t i){ return mVals[i]; }
	
	/// Get nth delayed element
	const T& operator[](uint32_t i) const { return mVals[i]; }

	/// Input element and return Nth delayed element
	T operator()(const T& v) const {
		const T r = mVals[N-1];
		for(uint32_t i=N-1; i>0; --i) mVals[i] = mVals[i-1];
		mVals[0]=v;
		return r;
	}

	/// Returns delay size number of elements
	uint32_t size() const { return N; }

	#undef DO

protected:
	mutable T mVals[N];
};



/// One element delay
template<class T=gam::real> 
class Delay1 : public FixedDelay<1,T>{
public:

	/// @param[in] v	Initial value of elements
	Delay1(const T& v=T()): FixedDelay<1,T>(v){}
};


/// Two element delay
template<class T=gam::real> 
class Delay2 : public FixedDelay<2,T>{
public:

	/// @param[in] v	Initial value of elements
	Delay2(const T& v=T()): FixedDelay<2,T>(v){}
	
	/// @param[in] v2	Initial value of 2nd delayed element
	/// @param[in] v1	Initial value of 1st delayed element
	Delay2(const T& v2, const T& v1){ (*this)[1]=v2; (*this)[0]=v1; }
};



/// Returns every nth input n times.
template <class T=gam::real>
struct Hold{

	/// @param[in] n	Period (in elements) of hold
	Hold(uint32_t n=1): mCount(0), mPeriod(n){}

	/// Filter input
	T operator()(const T& v) const {
		if(mCount==0) mVal=v;
		if((++mCount) >= mPeriod) mCount=0;
		return mVal;
	}

	/// Set period (in elements) of hold
	void period(uint32_t v){ mPeriod=v; }
	
	/// Reset hold counter
	void reset(){ mCount=0; }

private:
	mutable uint32_t mCount;
	uint32_t mPeriod;
	mutable T mVal;
};



//template <int N, class Tv, class Tc>
//class FixedSeries{
//public:
//
//protected:
//	
//};
//
//
//template <int Ni, int No, class Tv, class Tp>
//class LCCD{
//public:
//
//private:
//	Tv mi[Ni];
//	Tp ma[Ni];
//	Tv mo[No];
//	Tp mb[No];
//};



/// Integrates input with previous output
template <class Tv=double, class Tp=double>
class Integrator{
public:

	/// @param[in] v	Initial value
	Integrator(const Tp& leakCoef=Tp(1), const Tv& v=Tv(0)){
		mo[0]=v;
		leak(leakCoef);
	}

	/// Filter input value
	Tv operator()(const Tv& i0) const { return mo[0]=mo[0]*mb[0] + i0; }
	
	Integrator& leak(const Tp& v){ mb[0]=v; return *this; }
	Integrator& zero(){ mo[0]=Tv(0); return *this; }

protected:
	mutable Tv mo[1];
	Tp mb[1];
};


/// One-pole filter
template <class Tv=double, class Tp=double>
class OnePole{
public:

	OnePole(const Tv& prev=Tv(0), const Tp& smooth=Tp(1), const Tp& gain=Tp(1))
	{
		mo[1] = prev;
		set(smooth, gain);
	}

	OnePole& set(const Tp& smooth, const Tp& gain=Tp(1)){
		ma[0] = gain;
		mb[0] = smooth;
		return *this;
	}

	Tv operator()(const Tv& v) const {
		return mo[0] = mo[0]*mb[0] + v*ma[0];
	}

protected:
	mutable Tv mo[1];
	Tp ma[1], mb[1];
};




/// Complex resonator

/// A complex resonator consists of two complex numbers- one is an absolute
/// phase and amplitude and the other is a relative (differential) frequency 
/// and decay/grow factor.
template <class T=double>
class Reson : public Complex<T>{
public:
	typedef Complex<T> C;
	using C::operator();
	using C::operator=;

	/// @param[in] frq	unit frequency
	/// @param[in] amp	amplitude
	/// @param[in] phs	unit phase
	/// @param[in] dec	unit decay/grow factor
	Reson(const T& frq=T(0), const T& amp=T(1), const T& phs=T(0), const T& dec=T(1)){
		set(frq, amp, phs);
	}

	/// Advance one iteration and return value
	C operator()(){ return (*this)*=mFactor; }

	/// Filter input
	C operator()(const C& v){ return (*this) = (*this)*mFactor + v; }
	C operator()(const T& v){ return (*this) = (*this)*mFactor + v; }
	
	/// Recede one iteration and return value
	C recede(){ return (*this)/=mFactor; }

	/// Set amplitude
	void amp(const T& v){ (*this).fromPolar(v, this->arg()); }
	
	/// Set amplitude decay/grow factor after N iterations
	void decay(const T& target, const T& N=T(1)){
		// NOTE: this handles negative decays, thought better to leave this to frequency component
		//factor(freq(), (1==N) ? target : (::pow(scl::abs(target), 1./N)*scl::sgn(target)));
		factor(freq(), (1==N) ? target : (::pow(scl::abs(target), 1./N)));
	}

	/// Set recursive multiplication factor (frequency and decay/growth factor)
	void factor(const T& frq, const T& dec=T(1)){
		mFactor.fromPolar(dec, frq*M_2PI);
	}
	
	void factor(const Complex<T>& v){ mFactor=v; }

	/// Set unit frequency
	void freq(const T& v){ factor(v, decay()); }

	/// Set unit frequency, amplitude, unit phase, and decay/grow factor
	
	/// The phase state will be rewound 1 iteration so the first function call
	/// will return a complex number at the desired phase.
	void set(const T& frq, const T& amp, const T& phs, const T& dec=T(1)){
		this->fromPolar(amp, phs*M_2PI);
		factor(frq, dec);
		recede();
	}

	void set(const T& frq, const Complex<T>& phs){
		(*this)(phs.r, phs.i); freq(frq);
	}

	/// Get value one iteration ahead of current state
	C ahead() const { return (*this)*mFactor; }
	
	/// Get value one iteration behind current state
	C behind() const { return (*this)/mFactor; }

	/// Get unit decay
	T decay() const { return mFactor.mag(); }
	
	/// Get unit frequency
	T freq() const { return mFactor.phase()*M_1_2PI; }

	const C& factor() const { return mFactor; }
	C& factor(){ return mFactor; }

protected:
	C mFactor;
	
	// Set 60 dB decay interval
	//void decay(const Tv& v){ width(T(2.198806796637603 /* -ln(0.001)/pi */)/v); }
	
	// Set unit bandwidth
	//void width(const Tv& v){ mDecay=::exp(-M_PI*v); freq(freq()); }
};

typedef Reson<float>	Resonf;
typedef Reson<double>	Resond;



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
	gam::Complex<double> operator()(double f){
		gam::Complex<double> X(0,0), Y(1,0);
		f *= M_2PI;
		for(uint32_t i=0; i<mX.size(); ++i) X += mX[i].response(f);
		for(uint32_t i=0; i<mY.size(); ++i) Y -= mY[i].response(f);
		return X/Y * mGain; // H(z) = Y(z)/X(z)
	}

protected:

	struct DelayUnit{
		/// param[in] c		weighting coefficient
		/// param[in] d		delay in samples
		DelayUnit(double c_, double d_): c(c_), d(d_){}
		gam::Complex<double> response(double f){ return gam::Polar<double>(c, f*d); }
		double c, d;
	};

	std::vector<DelayUnit> mX, mY;
	double mGain;
};




// Direct domain lowpass filter
//template <class T=gam::real>
//class Lowpass{
//public:
//	Lowpass(const T& cen=T(0), const T& wid=T(0.1), const T& peak=T(0))
//	:	mCenter(cen), mPeak(peak)
//	{ width(wid); }
//
//	T operator()(const T& x) const {
//		const T c = x - mCenter;
//		return num(c)/den(c)*0.5+0.5;
//	}
//	
////	T den(const T& x) const { return mWidth+scl::abs(x); }
////	T num(const T& x) const { return mPeak*mWidth-x; }
//
//	T den(const T& x) const { return mWidth+x*x; }
//	T num(const T& x) const { return mPeak*mWidth-x*scl::abs(x); }
//
//	void width(const T& v){
//		static const T eps = 0.000001;	// width must be greater than zero
//		mWidth = scl::abs(v);
//		if(mWidth < eps) mWidth = eps;
//	}
//
//	void set(const T& cen, const T& wid, const T& peak){
//		mCenter=cen; width(wid); mPeak=peak;
//	}
//
//protected:
//	T mCenter, mWidth, mPeak;
//};
//
//
//// Direct domain inverse comb filter
//template <class T=gam::real>
//class CombI{
//public:
//
//	CombI(const T& frq=T(1), const T& resn=T(0.9), const T& phs=T(0))
//	:	mFreq(frq), mPhase(phs)
//	{ res(resn); }
//
//	T operator()(const T& x) const { return num(x)/den(x); }
//	T den(const T& x) const { return T(1) + mRes*(mRes - T(2)*cos(M_2PI*(mFreq*x + mPhase))); }
//	T num(const T& x) const { return mNum; }
//
//	void res(const T& v){
//		static const T eps = 0.999999;	// width must be greater than zero
//		mRes = scl::abs(v);
//		if(mRes>eps) mRes=eps;
//		mNum = T(1) + mRes*(mRes - T(2));
//	}
//	
//	void set(const T& frq, const T& resn, const T& phs){
//		mFreq=frq; res(resn); mPhase=phs;
//	}
//	
//protected:
//	T mFreq, mRes, mPhase;
//	T mNum;
//};

} // fil::
} // gam::

#endif
