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

/// Variable curve functor.

/// This curve will return samples in the interval [0, 1] starting from 0 and
/// ending on 1 over its length in samples.  The last point is exclusive, so
/// it takes length + 1 samples to reach 1 inclusively.  For iterations 
/// exceeding the specified length, the values returned will be unbounded.
template <class T=gam::real>
class Curve{
public:
	Curve();
	
	/// @param[in] length	Length of curve in samples
	/// @param[in] curve	Curvature value
	Curve(T length, T curve);

	void reset();					///< Reset envelope
	void set(T length, T curve);	///< Set length and curvature
	
	T operator()();	///< Generates next value
	T value();		///< Returns current value

protected:
	T mValue;
	T mMul, a1, b1;
};



/// Attack-decay envelope
template <class T=gam::real, class Ts=Synced>
class AD : public Ts{
public:

	/// @param[in] unitsA	Attack length
	/// @param[in] unitsD	Decay length
	/// @param[in] crvA		Attack curvature
	/// @param[in] crvD		Decay curvature
	AD(T unitsA, T unitsD, T crvA = 4, T crvD = -4);
	
	bool done() const;			///< Returns whether value is below threshold
	
	T operator()();				///< Generates next sample
	
	void attack(T units);		///< Set attack units
	void curve(T valA, T valD);	///< Set attack/decay curvatures
	void decay(T units);		///< Set decay units
	void reset();				///< Reset envelope
	void set(T unitsA, T unitsD);	///< Set attack/decay units
	
	virtual void onResync(double r);

protected:
	Curve<T> mFncA, mFncD;
	T mUnitsA, mUnitsD;
	T mCrvA, mCrvD;
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
	bool done() const { printf("%g\n", mRemain); return mClosed; }

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
	
	/// Generate next value
	Tv operator()(){
		Tp f = mAcc.val;
		if(f >= (Tp)1) return mIpl.val();
		mAcc();
		return mIpl(scl::min(f, (Tp)1));
	}
	
	/// Set new end value.  Start value is set to current value.
	void operator= (Tv v){
		mIpl.val(mIpl(scl::min(mAcc.val, (Tp)1)));
		mIpl.push(v);
		reset();
	}

	/// Generates a new end point from a generator when the segment end is reached
	
	/// This is useful for creating pitched noise
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
	
	/// Generate next value
	T operator()(){
		if(mFnc.value() >= (T)1) return mVal0;
		return ipl::linear(scl::min(mFnc(), (T)1), mVal1, mVal0);
	}
	
	/// Set new end value.  Start value is set to current value.
	void operator= (T v){
		mVal1 = ipl::linear(scl::min(mFnc.value(), (T)1), mVal1, mVal0);
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



// Exponentially attacking and decaying curve.
//TEM class AttDec : public Synced{
//public:
//	/// @param[in] dec60	Number of units until value decays -60db
//	/// @param[in] dec60	Number of units until value decays -60db
//	AttDec(double att60, double dec60);
//	
//	T next();				///< Returns next curve sample.
//	T operator () ();		///< Returns next curve sample.
//	
//	void decay(T value);	///< Set decay multiplication factor.
//	void value();
//	
//	/// Set number of units for curve to decay -60 dB.
//	void decay60(double value);
//	
//	T attack() const;
//	T decay() const;		///< Returns decay multiplication factor.
//	T value() const;		///< Returns current value.
//	
//	virtual void onSyncChange();
//	
//protected:
//	T mValA, mValD, mAtt, mAtt60, mDcy, mDcy60;
//};
//
//TEM inline void Decay<T>::attack60(double value){
//	mAtt60 = value;
//	mAtt = (T)scl::t60(value * spu());
//}
//
//TEM inline void Decay<T>::decay60(double value){
//	mDcy60 = value;
//	mDcy = (T)scl::t60(value * spu());
//}
//
//TEM inline void Decay<T>::next(){
//	if(mValA > (T)0.001){
//		T o = (T)1 - mValA;
//		mValA *= mAtt;
//		return o;
//	}
//	else{
//		T o = mValD;
//		mValD *= mDcy;
//		return o;
//	}
//}
//
//TEM inline void Decay<T>::next(){
//	if(mValA > (T)0.001){
//		T o = (T)1 - mValA;
//		mValA *= mAtt;
//		return o;
//	}
//	else{
//		T o = mValD;
//		mValD *= mDcy;
//		return o;
//	}
//}
//o1 = scl::dot2(o1, mStored, co1, ci0);



// Implementation_______________________________________________________________

#define TEM template <class T>

//---- Curve
TEM Curve<T>::Curve(): mValue(0), mMul(0), a1(0), b1(0){}

TEM Curve<T>::Curve(T length, T curve){
	set(length, curve);
	mValue = (T)0;
}

TEM inline void Curve<T>::reset(){ mValue = (T)0; b1 = a1; }

#define EPS (T)0.00001
TEM void Curve<T>::set(T length, T curve){

	// Avoid discontinuity when curve = 0 (a line)
	if(curve < EPS && curve > -EPS){
		curve = curve < (T)0 ? -EPS : EPS;
	}
	
	T curve_length = curve / length;
	
	if(curve_length < EPS && curve_length > -EPS){
		curve_length = curve_length < (T)0 ? -EPS : EPS;
		curve = curve_length * length;
	}

	mMul = ::exp(curve_length);
	a1 = (T)1 / ((T)1 - ::exp(curve));
	b1 = a1;

//	level = start;
//
//	grow = exp(curve / dur);
//	
//	a1 = (end - start) / (1 - exp(curve));
//	a2 = start + a1;
//	b1 = a1;
}

#undef EPS

TEM inline T Curve<T>::operator()(){
	mValue = a1 - b1;
	b1 *= mMul;
	return mValue;
}

TEM inline T Curve<T>::value(){ return mValue; }

#undef TEM



#define TM1 template <class T, class Ts>
#define TM2 T,Ts

//---- AD

TM1 AD<TM2>::AD(T unitsA, T unitsD, T crvA, T crvD) :
	mUnitsA(unitsA), mUnitsD(unitsD), mCrvA(crvA), mCrvD(crvD), mStage(0), mCntA(0)
{
	Ts::initSynced();
}

TM1 void AD<TM2>::attack(T units){
	mUnitsA = units;
	T smps = units * Ts::spu();
	mSmpsA = (uint32_t)smps;
	mFncA.set(smps, mCrvA);
}

TM1 void AD<TM2>::curve(T valA, T valD){
	mCrvA = valA; mCrvD = -valD;
	attack(mUnitsA);
	decay(mUnitsD);
}

TM1 void AD<TM2>::decay(T units){
	mUnitsD = units;
	mFncD.set(units * Ts::spu(), mCrvD);
}

TM1 inline void AD<TM2>::reset(){
	mStage = 0;
	mCntA = 0;
	mFncA.reset();
	mFncD.reset();
}

TM1 void AD<TM2>::set(T unitsA, T unitsD){ attack(unitsA); decay(unitsD); }

TM1 inline bool AD<TM2>::done() const { return mStage == 2; }

TM1 inline T AD<TM2>::operator()(){
	switch(mStage){
	case 0:		
		if(mCntA++ < mSmpsA)	return scl::min(mFncA(), (T)1);
		else					mStage = 1;
		
	case 1:{
		T v = (T)1 - mFncD();
		if(v > (T)0)			return v;
		else					mStage = 2;
	}
	default: return (T)0;
	}
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
