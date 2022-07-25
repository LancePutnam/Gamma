#ifndef GAMMA_ANALYSIS_H_INC
#define GAMMA_ANALYSIS_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include <cmath> // abs, fpclassify
#include "Gamma/Filter.h"

namespace gam{

///\defgroup Analysis


/// Envelope follower

/// This object produces an estimate of the amplitude envelope of a signal
/// by feeding a full-wave rectification of the signal through a low-pass filter.
///\ingroup Filter Envelope Analysis
template <class Tv=real, class Tp=real, class Td=GAM_DEFAULT_DOMAIN>
class EnvFollow : public Td {
public:

	/// \param[in] freq		Cutoff frequency of smoothing filter
	EnvFollow(Tp freq=10)
	:	lpf(freq){}


	/// Filter next sample
	Tv operator()(Tv i0){
		using std::abs;
		return lpf(abs(i0));
	}

	/// Returns current amplitude estimate
	Tv value() const { return lpf.last(); }

	/// Set lag length of filter
	EnvFollow& lag(Tp v){ lpf.lag(v); return *this; }

	/// Checks if current estimate is less than a threshold
	bool done(Tv eps=0.001) const { return value() < eps; }

	OnePole<Tv,Tp,Td> lpf;	///< Low-pass filter
};



/// Silence detector

/// This returns true if the magnitude of the input signal remains less than
/// some threshold over a specified number of samples.
/// \ingroup Analysis
class SilenceDetect{
public:
	SilenceDetect(unsigned count = 1000)
	:	mCount(count)
	{}


	/// Set number of samples required to trigger silence

	/// This is the number of contiguous samples that must be below the
	/// threshold magnitude in order to trigger a silence detection.
	SilenceDetect& count(unsigned v){
		mCount=v; return *this;
	}

	/// Reset the silence counter
	void reset(){ mNumSilent = 0; }

	/// Detect silence in input signal
	/// \param[in] input		The input signal
	/// \param[in] threshold	Magnitude below which a signal is considered silent
	/// \returns true if silence was detected, otherwise false
	template <typename T>
	bool operator()(const T& input, const T& threshold=T(0.001)){
		using std::abs;
		if(abs(input) < threshold){
			++mNumSilent;
			return done();
		}
		reset();
		return false;
	}

	/// Returns true if silence is being detected
	bool done() const { return mNumSilent >= mCount; }

private:
	unsigned mNumSilent=0, mCount;
};



/// Compares signal magnitude to a threshold

/// This filter compares the input magnitude to a threshold and returns 1 if 
/// it's greater than the threshold and 0 otherwise. The output is sent through
/// a one-pole low-pass filter.
///\ingroup Analysis
template <class T=gam::real>
class Threshold{
public:
	/// \param[in] thresh	Comparing threshold
	/// \param[in] freq		Cutoff frequency of output smoother
	Threshold(T thresh, T freq=10):lpf(freq), thresh(thresh){}
	
	T operator()(T in, T hi, T lo){
		using std::abs;
		return lpf(abs(in) > thresh ? hi : lo);
	}
	
	/// Returns 0 if less than threshold, 1 otherwise
	T operator()(T in){ return (*this)(in, T(1), T(0)); }
	
	/// Returns 1 if less than threshold, 0 otherwise
	T inv(T in){ return (*this)(in, T(0), T(1)); }
	
	OnePole<T> lpf;	///< Output smoother
	T thresh;		///< Threshold value
};



/// Zero-crossing detector

/// This object determines when a zero crossing occurs. It can distinguish 
/// between both positive (rising) and negative (falling) zero crossings.
///\ingroup Analysis
template<class Tv=gam::real>
class ZeroCross{
public:

	ZeroCross(Tv prev = Tv(0)): mPrev(prev){}

	/// Detect zero crossing

	/// \returns 
	///		 0 if no zero crossing,
	///		-1 if a negative (falling) zero crossing, or	
	///		 1 if a positive (rising) zero crossing.
	int operator()(Tv input){
		int pzc = int((input > Tv(0)) && (mPrev <= Tv(0)));
		int nzc =-int((input < Tv(0)) && (mPrev >= Tv(0)));
		mPrev = input;
		return pzc + nzc;
	}

private:
	Tv mPrev;
};



/// Periodic counter

/// Triggers periodically after reaching a specified maximum count. Useful for
/// determining when a window of samples has passed.
class PCounter{
public:

	PCounter(unsigned period=256): mPeriod(period){}

	/// Set period
	PCounter& period(unsigned n){ mPeriod=n; return *this; }
	unsigned period() const { return mPeriod; }

	/// Whether counter just finished one period
	bool cycled() const { return 0==mCount; }

	/// Get the current count
	unsigned count() const { return mCount; }

	/// Reset count to zero
	PCounter& reset(){ mCount=0; return *this; }

	/// Returns true when end of period is reached, otherwise false
	bool operator()(){
		++mCount;
		if(mCount == mPeriod){
			mCount = 0;
			return true;
		}
		return false;
	}

private:
	unsigned mPeriod;
	unsigned mCount = 0;
};



/// Maximum absolute value over window

/// This produces an accurate estimate of a signal's amplitude envelope by 
/// measuring the maximum absolute value over a given window size.
/// Unlike estimators that use low-pass filters, this estimator is less affected 
/// by the frequency content of the input signal, namely high frequencies.
template <class Tv=gam::real>
class MaxAbs : public PCounter {
public:

	MaxAbs(int winSize=256): PCounter(winSize){}

	const Tv& operator()(const Tv& in){
		using std::abs;
		auto inAbs = abs(in);
		if(inAbs > mMaxCalc) mMaxCalc = inAbs;
		if(PCounter::operator()()){
			mInc = (mMaxCalc - mMaxPrev) / period();
			mVal = mMaxPrev - mInc;
			mMaxPrev = mMaxCalc;
			mMaxCalc = Tv(0);
		}
		mVal += mInc;
		return mVal;
	}

	/// Current (non-smoothed) value
	const Tv& value() const { return mMaxPrev; }

	/// Current linearly smoothed value
	const Tv& valueL() const { return mVal; }

private:
	Tv mMaxCalc = Tv(0);
	Tv mMaxPrev = Tv(0);
	Tv mVal = Tv(0);
	Tv mInc = Tv(0);
};



/// Computes the zero-crossing rate of an input signal

/// The zero-crossing rate (ZCR) is a measure proportional to how many times a 
/// signal crosses zero over a given number of samples. The result very closely 
/// matches the spectral centroid of the input signal. High ZCR correlates to
/// noisy, bright and high-pitched signals whereas low ZCR correlates to dark
/// and low-frequency signals. If applied to sampled audio, low amplitude
/// sections should be rejected as the true signal is likely buried in noise.
///\ingroup Analysis
template <class Tv=gam::real>
class ZeroCrossRate : public PCounter {
public:

	/// \param[in] winSize		size of analysis window
	ZeroCrossRate(int winSize=256)
	:	PCounter(winSize)
	{}

	/// Get the current zero-crossing rate, in [0, 1]
	float value() const { return mRate; }

	/// Input next sample and return current zero-crossing rate
	Tv operator()(Tv input){
		if(0 != mDetector(input)){
			++mCrosses;
		}
		if(PCounter::operator()()){
			mRate = float(mCrosses) / period();
			mCrosses = 0;
		}
		return mRate;
	}

private:
	ZeroCross<Tv> mDetector;
	float mRate = 0.f;
	unsigned mCrosses = 0;
};



/// Analyzes a signal for potential errors or bad states
class Inspector {
public:
	/// Error status
	enum {
		ERR_DC			= 1<<0,	///< Signal has DC
		ERR_CLIP		= 1<<1,	///< Signal out of [-1,1] range
		ERR_DENORMAL	= 1<<2,	///< Signal in denormal (subnormal) range
		ERR_INF			= 1<<3,	///< Signal has an infinity
		ERR_NAN			= 1<<4	///< Signal has a NaN (not-a-number)
	};

	/// Input next sample; \returns error status
	unsigned operator()(float v){
		auto fpstate = std::fpclassify(v);
		switch(fpstate){
		case FP_INFINITE:  mErr |= ERR_INF; break;
		case FP_NAN:       mErr |= ERR_NAN; break;
		case FP_SUBNORMAL: mErr |= ERR_DENORMAL; break;
		default:;
		}

		mRunSum += v;
		if(mDCCount()){
			mDC = mRunSum / mDCCount.period();
			mRunSum = 0.;
			auto absDC = std::abs(mDC);
			if(absDC > mMaxDC) mMaxDC = absDC;
			if(mMaxDC > 0.005) mErr |= ERR_DC;
		}

		auto absv = std::abs(v);
		if(absv > mPeak) mPeak = absv;
		if(mPeak > 1.) mErr |= ERR_CLIP;

		return mErr;
	}
	
	unsigned status() const { return mErr; }	///< Get error status
	operator bool() const { return mErr; }
	float peak() const { return mPeak; }		///< Get running peak value
	float DC() const { return mDC; }			///< Get current DC value
	float maxDC() const { return mMaxDC; }		///< Get running max DC value

	/// Reset analysis and error states
	Inspector& reset(){
		mErr = 0;
		mDC = mMaxDC = mRunSum = 0.;
		mDCCount.reset();
		mPeak = 0.;
		return *this;
	}

	void print() const {
		printf("peak:%5.3f DC:% 7.4f ", mPeak, mDC);

		if(mErr){
			printf("[");
			if(hasErr(ERR_DC      )) printf("DC:% 7.4f ", mMaxDC);
			if(hasErr(ERR_CLIP    )) printf("clip:%5.3f ", mPeak);
			if(hasErr(ERR_DENORMAL)) printf("denormal ");
			if(hasErr(ERR_INF     )) printf("inf ");
			if(hasErr(ERR_NAN     )) printf("NaN ");
			printf("\b]");
		}
		printf("\n");
	}

private:
	unsigned mErr = 0;
	float mPeak = 0.;
	float mDC = 0., mMaxDC = 0., mRunSum = 0.;
	PCounter mDCCount{4096};
	bool hasErr(unsigned e) const { return mErr&e; }
};

} // gam::
#endif
