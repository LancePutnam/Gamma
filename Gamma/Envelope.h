#ifndef GAMMA_ENVELOPE_H_INC
#define GAMMA_ENVELOPE_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include <math.h>
#include <float.h>

#include "Gamma/gen.h"
#include "Gamma/ipl.h"
#include "Gamma/scl.h"
#include "Gamma/Strategy.h"
#include "Gamma/Sync.h"

namespace gam{

/// Exponential curve with variable curvature

/// This curve will return values in the interval [0, end] starting from 0 and
/// ending on 'end' over its length in samples.  The last point is exclusive, so
/// it takes length + 1 samples to reach 'end' inclusively.  For iterations 
/// exceeding the specified length, the values returned will be unbounded.
template <class T=gam::real>
class Curve{
public:
	Curve();
	
	/// @param[in] length	length of curve in samples
	/// @param[in] curve	curvature; pos. approaches slowly, neg. approaches rapidly, 0 approaches linearly
	/// @param[in] end		end value
	Curve(T length, T curve, T end = T(1));

	/// Returns whether curve has gone past end value
	bool done() const;

	/// Get end value
	T end() const { return mEnd; }

	/// Get current value
	T value() const;

	T operator()();					///< Generates next value
	void reset();					///< Reset envelope

	/// Set length and curvature
	
	/// @param[in] length	length of curve in samples
	/// @param[in] curve	curvature; pos. approaches slowly, neg. approaches rapidly, 0 approaches linearly
	/// @param[in] end		end value
	void set(T length, T curve, T end = T(1));
	
	void value(const T& v);

protected:
	T mEnd;
	T mMul, mA, mB;
	
	T eps() const;
};



template <int N, class T=gam::real>
class CurveEnv{
public:

	CurveEnv(): mRelease(N)
	{ reset(); }

	int size() const { return N; }
	int position() const { return mPos; }
	int releasePoint() const { return mRelease; }
	int stage() const { return mStage; }
	
	T value() const {
		return mValues[mStage] + mCurve.value();
	}

	bool done() const { return mStage == size(); }
	bool released() const { return mRelease < 0; }
	bool sustained() const { return (mStage == mRelease) && !released(); }

	T operator()(){
		//begin:
		if(sustained()){
			return mValues[mStage];
		}
		else if(mPos < mLen){
			++mPos;
			return mValues[mStage] + mCurve();
		}
//		else if(mStage < size()-1){
		else if(mStage < size()){
			++mStage;

			if(!done()){
				mPos = 0;
				mLen = mLengths[mStage];
				mCurve.set(mLen, mCurves[mStage], mValues[mStage+1]-mValues[mStage]);
				//goto begin;
				(*this)();
			}
		}
//		return mValues[mStage+1];
		return mValues[mStage];
	}
	
	void release(){
		mRelease=-scl::abs(mRelease);

		// begin release portion immediately starting at current level
//		T curVal = value();
//		mStage = -mRelease;
//		mPos = 0;
//		mLen = mLengths[mStage];
//		mCurve.set(mLen, mCurves[mStage], 2*mValues[mStage+1]-mValues[mStage]-curVal);
//		mCurve.value(curVal-mValues[mStage]);
	}

	/// Sets the point at which the envelope holds its value until released
	void releasePoint(int v){ mRelease = v; }

	/// Reset envelope to starting point
	void reset(){
		mPos = 0xFFFFFFFF;
		mLen = 0;
		mStage = -1;
		mRelease = scl::abs(mRelease);
	}

	/// Set length and curvature of a segment
	CurveEnv& segment(int i, T length, T curve){
		mLengths[i]=length;
		mCurves [i]=curve;
		return *this;
	}

	/// Set length and curvature of many segments
	template <class V>
	CurveEnv& segments(const V* lengths, const V* curves, int len, int begin=0){
		int max = size() - begin;
		int n = len < max ? len : max;
		for(int i=0; i<n; ++i){
			segment(i+begin, lengths[i], curves[i]);
		}
		return *this;
	}
	
	CurveEnv& segments(T la, T ca, T lb, T cb){ T l[]={la,lb}; T c[]={ca,cb}; return segments(l,c,2); }
	CurveEnv& segments(T la, T ca, T lb, T cb, T lc, T cc){ T l[]={la,lb,lc}; T c[]={ca,cb,cc}; return segments(l,c,3); }
	CurveEnv& segments(T la, T ca, T lb, T cb, T lc, T cc, T ld, T cd){ T l[]={la,lb,lc,ld}; T c[]={ca,cb,cc,cd}; return segments(l,c,4); }

	CurveEnv& point(int i, T val){ mValues[i]=val; return *this; }

	template <class V>
	CurveEnv& points(const V* vals, int len){
		int n = len <= size() ? len : size()+1;
		for(int i=0; i<n; ++i) point(i, vals[i]);
		return *this;
	}
	
	CurveEnv& points(T a, T b){ T v[]={a,b}; return points(v,2); }
	CurveEnv& points(T a, T b, T c){ T v[]={a,b,c}; return points(v,3); }
	CurveEnv& points(T a, T b, T c, T d){ T v[]={a,b,c,d}; return points(v,4); }
	CurveEnv& points(T a, T b, T c, T d, T e){ T v[]={a,b,c,d,e}; return points(v,5); }

protected:
	Curve<T> mCurve;
	T mLengths[N];			// segment lengths
	T mCurves[N];			// segment curvatures
	T mValues[N+1];			// break point values
	
	uint32_t mPos, mLen;	// position in and length of current segment, in samples
	int mStage;				// the current curve segment
	int mRelease;			// index of point before release portion
};



/// Attack-decay envelope
template <class T=gam::real, class Ts=Synced>
class AD : public Ts{
public:

	/// @param[in] lengthA		Attack length
	/// @param[in] lengthD		Decay length
	/// @param[in] curveA		Attack curvature
	/// @param[in] curveD		Decay curvature
	/// @param[in] amp			Amplitude
	AD(T lengthA = 0.01, T lengthD = 2, T curveA =-4, T curveD = 4, T amp = 1);
	
	bool done() const;			///< Returns whether value is below threshold
	T amp() const;				///< Get amplitude (maximum value)
	T value() const;			///< Get current value
	
	T operator()();				///< Generates next sample
	
	void amp(T v);				///< Set maximum amplitude
	void attack(T units);		///< Set attack units
	void curve(T valA, T valD);	///< Set attack/decay curvatures
	void decay(T units);		///< Set decay units
	void length(T unitsA, T unitsD);	///< Set attack/decay units
	void set(T lengthA, T lengthD, T curveA, T curveD, T amp); ///< Set all envelope parameters
	void reset();				///< Reset envelope
	
	virtual void onResync(double r);

protected:
	Curve<T> mFncA, mFncD;
	T mLenA, mLenD;
	T mCrvA, mCrvD;
	T mAmp;
	uint32_t mStage, mSmpsA, mCntA;
};




/// Exponentially decaying  curve.
template <class T=gam::real, class Ts=Synced>
class Decay : public Ts{
public:

	/// @param[in] decay	Number of units until initial value decays -60 dB
	/// @param[in] val		Intial value
	Decay(T decay=T(1), T val=T(1));

	T decay() const;		///< Returns -60 dB decay length.
	bool done(T thresh=T(0.001)) const; ///< Returns whether value is below threshold
	T value() const;		///< Returns current value.

	T operator()();			///< Generates next sample.
	
	void decay(T val);		///< Set number of units for curve to decay -60 dB.
	void reset();			///< Set current value to 1.
	void value(T val);		///< Set current value.

	virtual void onResync(double r);
	
protected:
	T mVal, mMul, mDcy;
};


/// Binary gate controlled by threshold comparison
template <class T=gam::real, class Ts=Synced>
class Gate : public Ts{
public:

	/// @param[in] closingDelay		units to wait before closing while under threshold
	/// @param[in] threshold		threshold below which gate closes
	Gate(double closingDelay=0, double threshold=0.001)
	:	mDelay(closingDelay), mRemain(closingDelay), mThresh(threshold), mClosed(0)
	{}

	/// Check whether gate is closed
	bool done() const { return mClosed; }

	/// Filter value
	T operator()(const T& v){		
		if(gam::norm(v) < mThresh){
			mRemain -= Synced::ups();
			if(mRemain <= 0) return close();
		}
		else{
			mRemain = mDelay;
		}
		return open();
	}
	
	/// Set closing delay
	Gate& delay(double v){ mDelay=mRemain=v; return *this; }

protected:
	double mDelay, mRemain;
	double mThresh;
	int mClosed;
	T close(){ mClosed=1; return T(0); }
	T  open(){ mClosed=0; return T(1); }
};


/// Interpolation envelope segment
template <class Tv=gam::real, template <class> class Si=iplSeq::Linear, class Tp=gam::real, class Ts=Synced>
class Seg : public Ts{
public:

	/// @param[in] len		Length of segment in domain units
	/// @param[in] start	Start value
	/// @param[in] end		End value
	Seg(Tp len=0.5, Tv start=1, Tv end=0, Tp phase=0):
		mFreq((Tp)1/len), mAcc(0, phase), mIpl(start)
	{
		mIpl.push(end);
		Ts::initSynced();
	}
	
	/// Returns whether envelope is done
	bool done() const { return mAcc.val >= Tp(1); }
	
	/// Generate next value
	Tv operator()(){
		Tp f = mAcc.val;
		if(done()) return mIpl.val();
		mAcc();
		return mIpl(scl::min(f, Tp(1)));
	}
	
	/// Set new end value.  Start value is set to current value.
	void operator= (Tv v){
		mIpl.val(mIpl(scl::min(mAcc.val, Tp(1))));
		mIpl.push(v);
		reset();
	}

	/// Generates a new end point from a generator when the segment end is reached
	
	/// This is useful for creating pitched noise from a random number generator.
	///
	template <class G>
	Tv operator()(G& g){
		Tp f = mAcc.val;
		Tv v;
		if(f >= Tp(1)){	v = mIpl.val(); (*this) = g(); }
		else{			v = val(); }
		mAcc();
		return v;
	}
	
	/// Set frequency of envelope
	void freq(Tp v){ mFreq = v; mAcc.add = v * Ts::ups(); }
	
	/// Set length in domain units.
	void period(Tp v){ freq((Tp)1/v); }

	/// Set phase along segment
	void phase(Tp v){ mAcc = v; }

	/// Reset envelope
	void reset(){ phase((Tp)0); }
	
	Tv val() const { return mIpl(scl::min(mAcc.val, Tp(1))); }
	
	Si<Tv>& ipol(){ return mIpl; }
	
	virtual void onResync(double r){ freq(mFreq); }
	
protected:
	Tp mFreq;
	gen::RAdd<Tp> mAcc;
	Si<Tv> mIpl;
};




/// Exponential envelope segment for smoothing out value changes.
template <class T=gam::real, class Ts=Synced>
class SegExp : public Ts{
public:

	/// @param[in] len		Length of segment in domain units
	/// @param[in] crv		Curvature of segment
	/// @param[in] start	Start value
	/// @param[in] end		End value
	SegExp(T len, T crv=-3, T start=1, T end=0):
		mLen(len), mCrv(crv), mVal1(start), mVal0(end)
	{
		Ts::initSynced();
	}
	
	/// Returns whether envelope is done
	bool done() const { return mFnc.value() >= T(1); }
	
	/// Generate next value
	T operator()(){
		if(done()) return mVal0;
		return ipl::linear(scl::min(mFnc(), T(1)), mVal1, mVal0);
	}
	
	/// Set new end value.  Start value is set to current value.
	void operator= (T v){
		mVal1 = ipl::linear(scl::min(mFnc.value(), T(1)), mVal1, mVal0);
		mVal0 = v;
		mFnc.reset();
	}
	
	/// Set curvature.  Negative gives faster change, positive gives slower change.
	void curve(T v){ set(mLen, v); }
	
	/// Set length in domain units.
	void period(T v){ set(v, mCrv); }

	void reset(){ mFnc.reset(); }

	/// Set length and curvature
	void set(T len, T crv){
		mLen = len; mCrv = crv;
		mFnc.set(len * Ts::spu(), crv);
	}
	
	virtual void onResync(double r){ set(mLen, mCrv); }
	
protected:
	T mLen, mCrv, mVal1, mVal0;
	Curve<T> mFnc;
};





// Implementation_______________________________________________________________

#define TEM template <class T>

//---- Curve
TEM Curve<T>::Curve(): mEnd(1), mMul(1), mA(0), mB(0){}

TEM Curve<T>::Curve(T length, T curve, T amp){
	set(length, curve, amp);
}

TEM inline bool Curve<T>::done() const { return scl::abs(mA - mB*mMul) >= scl::abs(end()); }

TEM inline T Curve<T>::value() const { return mA - mB; }

// dividing by mMul goes back one step
TEM inline void Curve<T>::reset(){ mB = mA / mMul; }

TEM inline T Curve<T>					::eps() const { return T(0.00001  ); }
template<> inline double Curve<double>	::eps() const { return   0.00000001; }


// hack to get proper max floating point value
namespace{
	template<class T> inline T maxReal(){ return DBL_MAX; }
	template<> inline float maxReal<float>(){ return FLT_MAX; }
}

TEM void Curve<T>::set(T len, T crv, T end){
	static const T EPS = eps();

	if(len == T(0)){ // if length is 0, return end value immediately
		mEnd = end;
		mMul = maxReal<T>();
		mA = end;
		mB = 0;
		return;
	}

	// Avoid discontinuity when curve = 0 (a line)
	if(crv < EPS && crv > -EPS){
		crv = crv < T(0) ? -EPS : EPS;
	}
	
	T crvOverLen = crv / len;
	
	if(crvOverLen < EPS && crvOverLen > -EPS){
		crvOverLen = crvOverLen < T(0) ? -EPS : EPS;
		crv = crvOverLen * len;
	}

	mEnd = end;
	mMul = ::exp(crvOverLen);
	mA = mEnd / (T(1) - ::exp(crv));
	reset();

//	level = start;
//
//	grow = exp(curve / dur);
//	
//	a1 = (end - start) / (1 - exp(curve));
//	a2 = start + a1;
//	b1 = a1;
}

TEM inline void Curve<T>::value(const T& v){ mB = mA-v; }

TEM inline T Curve<T>::operator()(){
	mB *= mMul;
	return value();
}
#undef TEM



#define TM1 template <class T, class Ts>
#define TM2 T,Ts

//---- AD

TM1 AD<TM2>::AD(T lenA, T lenD, T crvA, T crvD, T amp) :
	mLenA(lenA), mLenD(lenD), mCrvA(crvA), mCrvD(crvD), mAmp(amp),
	mStage(0), mCntA(0)
{
	Ts::initSynced();
}


TM1 inline bool AD<TM2>::done() const { return 2 == mStage; }

TM1 inline T AD<TM2>::amp() const { return mAmp; }

TM1 inline T AD<TM2>::value() const {
	switch(mStage){
		case 0:	return scl::min(mFncA.value(), amp());
		case 1: return scl::max(mFncD.value(), T(0));
		default: return T(0);
	}
}

TM1 inline T AD<TM2>::operator()(){
	switch(mStage){
	case 0:		
		if(mCntA++ < mSmpsA)	return scl::min(mFncA(), amp());
		else					mStage = 1;
		
	case 1:{
		T v = amp() - mFncD();
		if(v > T(0))			return v;
		else					mStage = 2;
	}
	default: return T(0);
	}
}

TM1 void AD<TM2>::amp(T v){
	mAmp = v;
	attack(mLenA);
	decay(mLenD);	
}

TM1 void AD<TM2>::curve(T valA, T valD){
	mCrvA = valA; mCrvD = -valD;
	attack(mLenA);
	decay(mLenD);
}

TM1 void AD<TM2>::attack(T v){
	mLenA = v;
	T smps = v * Ts::spu();
	mSmpsA = (uint32_t)smps;
	mFncA.set(smps, mCrvA, amp());
}

TM1 void AD<TM2>::decay(T v){
	mLenD = v;
	mFncD.set(v * Ts::spu(), mCrvD, amp());
}

TM1 void AD<TM2>::length(T a, T d){ attack(a); decay(d); }

TM1 void AD<TM2>::set(T lenA, T lenD, T curveA, T curveD, T amp){
	mAmp = amp;
	mCrvA = curveA; mCrvD = -curveD;
	attack(lenA);
	decay(lenD);	
}

TM1 inline void AD<TM2>::reset(){
	mStage = 0;
	mCntA = 0;
	mFncA.reset();
	mFncD.reset();
}

TM1 void AD<TM2>::onResync(double r){
	//printf("AD: onSyncChange()\n");
	curve(mCrvA, mCrvD);
}


//---- Decay

TM1 Decay<TM2>::Decay(T decay_, T val)
:	mVal(val)
{
	Ts::initSynced();
	decay(decay_);
}

TM1 inline T Decay<TM2>::operator()(){ T o = mVal; mVal *= mMul; return o; }

TM1 inline void Decay<TM2>::decay(T v){
	mDcy = v;
	mMul = scl::t60(v * Ts::spu());
}

TM1 inline void Decay<TM2>::reset(){ mVal = 1; }
TM1 inline void Decay<TM2>::value(T v){ mVal = v; }

TM1 inline T Decay<TM2>::decay() const { return mDcy; }
TM1 inline bool Decay<TM2>::done(T thr) const { return mVal < thr; }
TM1 inline T Decay<TM2>::value() const { return mVal; }

TM1 void Decay<TM2>::onResync(double r){ decay(mDcy); }

#undef TM1
#undef TM2



/*
// deprecated in favor of Seg
template <class Tv=gam::real, class Tp=gam::real, class Ts=Synced>
class LineSeg : public Ts{
public:

	/// @param[in] len		Length of segment in domain units
	/// @param[in] start	Start value
	/// @param[in] end		End value
	LineSeg(Tp len=0.5, Tv start=1, Tv end=0):
		mFreq((Tp)1/len), mVal1(start), mVal0(end), mAcc(0,0)
	{
		Ts::initSynced();
	}
	
	/// Generate next value
	Tv operator()(){
		Tp v = mAcc.val;
		if(v >= (Tp)1) return mVal0;
		mAcc();
		return ipl::linear(scl::min(v, (Tp)1), mVal1, mVal0);
	}
	
	/// Set new end value.  Start value is set to current value.
	void operator= (Tv v){
		mVal1 = ipl::linear(scl::min(mAcc.val, (Tp)1), mVal1, mVal0);
		mVal0 = v;
		reset();
	}
	
	/// Set frequency of envelope
	void freq(Tp v){ mFreq = v; mAcc.add = v * Ts::ups(); }
	
	/// Set length in domain units.
	void length(Tp v){ freq((Tp)1/v); }

	/// Reset envelope
	void reset(){ mAcc = (Tp)0; }
	
	Tv& end(){ return mVal0; }
	
	virtual void onResync(double r){ freq(mFreq); }
	
protected:
	Tv mVal1, mVal0;
	Tp mFreq;
	gen::RAdd<Tp> mAcc;
};
*/


} // gam::
#endif
