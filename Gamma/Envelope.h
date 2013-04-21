#ifndef GAMMA_ENVELOPE_H_INC
#define GAMMA_ENVELOPE_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

/// \defgroup Envelopes
/// Everything in Gamma having to do with envelopes.


#include <math.h>
#include <float.h>

#include "Gamma/Types.h"
#include "Gamma/gen.h"
#include "Gamma/ipl.h"
#include "Gamma/scl.h"
#include "Gamma/Strategy.h"
#include "Gamma/Sync.h"

namespace gam{

/// Exponential curve with variable curvature

/// This curve will return values in the interval [start, end] starting from 0
/// and ending on 'end' over its length in samples.  The last point is
/// exclusive, so it takes length + 1 samples to reach 'end' inclusively.
/// For iterations exceeding the specified length, the values returned will be
/// unbounded. \n\n
/// Given any two points, as long as they don't have the same value, there are
/// infinitely many exponential segments starting at the first and ending at the
/// second. One of these is a straight line between the two. Hence the "variable
/// curvature" of the Curve object. \n\n
/// Curve touches both its start and end points while Decay asymptotically
/// approaches zero.
///
/// \tparam Tv	value (sample) type
/// \tparam Tp	parameter type
/// \ingroup Envelopes
/// \sa Decay
template <class Tv=real, class Tp=real>
class Curve{
public:
	Curve();
	
	/// \param[in] length	length of curve in samples
	/// \param[in] curve	curvature, c, where 
	///						c > 0 approaches slowly (accelerates),
	///						c < 0 approaches rapidly (decelerates), and
	///						c = 0 approaches linearly
	/// \param[in] start	start value
	/// \param[in] end		end value
	Curve(Tp length, Tp curve, Tv start=Tv(1), Tv end=Tv(0));

	bool done() const;				///< Returns whether curve has gone past end value
	Tv end() const { return mEnd; }	///< Get end value
	Tv value() const;				///< Get current value

	Tv operator()();				///< Generates next value
	Curve& reset(Tv start=Tv(0));	///< Reset envelope
	Curve& value(const Tv& v);		///< Set value

	/// Set length and curvature
	
	/// \param[in] length	length of curve in samples
	/// \param[in] curve	curvature; pos. approaches slowly, neg. approaches rapidly, 0 approaches linearly
	/// \param[in] start	start value
	/// \param[in] end		end value
	Curve& set(Tp length, Tp curve, Tv start=Tv(0), Tv end=Tv(1));

protected:
	Tv mEnd, mA, mB;
	Tp mMul;
};



/// Envelope with a fixed number of exponential segments and a sustain point

/// The envelope consists of N exponential curve 'segments' and N+1 break-point
/// levels. The curvature and length of each segment and the break-point levels
/// can be controlled independently.
/// This class can be used to construct many specialized envelopes such as an AD 
/// (Attack Decay), an ADSR (Attack Decay Sustain Release), and an ADSHR (Attack
/// Decay Sustain Hold Release). The number of envelope segments is fixed to
/// ensure better memory locality.
///
/// \tparam N	number of segments
/// \tparam Tv	value (sample) type
/// \tparam Tp	parameter type
/// \ingroup Envelopes 
template <int N, class Tv=real, class Tp=real, class Ts=Synced>
class Env : public Ts{
public:

	Env()
	:	mSustain(N), mLoop(0){
		for(int i=0; i<N; ++i){
			mLengths[i]= 1e-8;
			mCurves[i] =-4;
			mLevels[i] = 1e-8;
		}	mLevels[N] = 1e-8;
		reset();
	}

	Env(Tp lvl1, Tp len1, Tp lvl2)
	:	mSustain(N), mLoop(0){
		levels(lvl1,lvl2);
		lengths()[0] = len1;
		curve(-4);
		reset();
	}

	Env(Tp lvl1, Tp len1, Tp lvl2, Tp len2, Tp lvl3)
	:	mSustain(N), mLoop(0){
		levels(lvl1,lvl2,lvl3);
		lengths(len1,len2);
		curve(-4);
		reset();
	}

	Env(Tp lvl1, Tp len1, Tp lvl2, Tp len2, Tp lvl3, Tp len3, Tp lvl4)
	:	mSustain(N), mLoop(0){
		levels(lvl1,lvl2,lvl3,lvl4);
		lengths(len1,len2,len3);
		curve(-4);
		reset();
	}

	/// Get the number of segments
	int size() const { return N; }
	
	/// Get the position, in samples, within the current segment
	int position() const { return mPos; }
	
	/// Get the sustain break-point
	int sustainPoint() const { return mSustain; }
	
	/// Get the envelope's current segment
	int stage() const { return mStage; }
	
	/// Get the current envelope value
	Tv value() const { return mCurve.value(); }

	/// Returns whether the envelope is done
	bool done() const { return mStage == size(); }
	
	/// Returns whether the envelope is released
	bool released() const { return mSustain < 0; }
	
	/// Returns whether the envelope is currently sustained
	bool sustained() const { return (mStage == mSustain); }


	/// Generate next value
	Tv operator()(){

		// Sustain stage:
		if(sustained()){
			return mLevels[mStage];
		}

		// Interpolating along segment:
		else if(mPos < mLen){
			++mPos;
			return mCurve();
		}

		// Just went past end of current segment and there are more left:
		else if(mStage < size()){
			++mStage;
			if(mLoop && done()) mStage=0;
			if(!done()){
				mPos = 0;
				setLen(mStage);
				mCurve.set(mLen, mCurves[mStage], mLevels[mStage], mLevels[mStage+1]);

				// Immediately return start level of new stage
				return (*this)();
			}
		}

		// Envelope is done:
		return mLevels[mStage];
	}	

	/// Release the envelope
	void release(){
		mSustain = -scl::abs(mSustain);

		// begin release portion immediately starting at current level
		Tv curVal = value();
		mStage = -mSustain;
		if(!done()){
			mPos = 0;
			setLen(mStage);
			mCurve.set(mLen, mCurves[mStage], curVal, mLevels[mStage+1]);
		}
	}


	/// Set whether envelope loops
	Env& loop(bool v){ mLoop=v; return *this; }

	/// Sets the point at which the envelope holds its value until released
	Env& sustainPoint(int v){ mSustain=v; return *this; }

	/// Disable sustain
	Env& sustainDisable(){ return sustainPoint(N); }

	/// Reset envelope to starting point
	void reset(){
		// this forces a stage increment upon first iteration
		mPos = 0xFFFFFFFF;
		mLen = 0;
		mStage = -1;
		mSustain = scl::abs(mSustain);
	}

	
	/// Get segment lengths array
	Tp * lengths(){ return mLengths; }

	/// Set break-point values
	template <class V>
	Env& lengths(const V* vals, int len){
		int n = len <= size() ? len : size();
		for(int i=0; i<n; ++i) lengths()[i] = vals[i];
		return *this;
	}
	
	/// Set first two segment lengths
	Env& lengths(Tp a, Tp b){ Tp v[]={a,b}; return lengths(v,2); }

	/// Set first three segment lengths
	Env& lengths(Tp a, Tp b, Tp c){ Tp v[]={a,b,c}; return lengths(v,3); }

	/// Set first four segment lengths
	Env& lengths(Tp a, Tp b, Tp c, Tp d){ Tp v[]={a,b,c,d}; return lengths(v,4); }

	/// Set first five segment lengths
	Env& lengths(Tp a, Tp b, Tp c, Tp d, Tp e){ Tp v[]={a,b,c,d,e}; return lengths(v,5); }

	/// Get total length of all envelope segments
	Tp totalLength() const {
		Tp sum=Tp(0);
		for(int i=0;i<size();++i) sum += mLengths[i];
		return sum;
	}
	
	/// Set total length of envelope by adjusting one segment length
	
	/// \param[in] length		desired length
	/// \param[in] modSegment	segment whose length is modified to match desired length
	Env& totalLength(Tp length, int modSegment){
		mLengths[modSegment] = Tp(0);
		mLengths[modSegment] = length - totalLength();
		return *this;
	}

	/// Set total length of envelope by scaling all segment lengths
	Env& totalLength(Tp length){
		Tp mul = length / totalLength();
		for(int i=0; i<size(); ++i){
			lengths()[i] *= mul;
		}
		return *this;
	}


	/// Get segment curvature array
	Tp * curves(){ return mCurves; }
	
	/// Set curvature of all segments
	Env& curve(Tp v){
		for(int i=0; i<N; ++i) curves()[i]=v;
		return *this;
	}


	/// Set length and curvature of a segment
	Env& segment(int i, Tp len, Tp crv){
		mLengths[i]=len;
		mCurves [i]=crv;
		return *this;
	}

	/// Set length and curvature of many segments
	template <class V>
	Env& segments(const V* lens, const V* crvs, int len, int begin=0){
		int max = size() - begin;
		int n = len < max ? len : max;
		for(int i=0; i<n; ++i){
			segment(i+begin, lens[i], crvs[i]);
		}
		return *this;
	}
	
	/// Set length and curvature of first two segments
	Env& segments(Tp la, Tp ca, Tp lb, Tp cb){
		Tp l[]={la,lb}; Tp c[]={ca,cb}; return segments(l,c,2); }

	/// Set length and curvature of first three segments
	Env& segments(Tp la, Tp ca, Tp lb, Tp cb, Tp lc, Tp cc){
		Tp l[]={la,lb,lc}; Tp c[]={ca,cb,cc}; return segments(l,c,3); }
	
	/// Set length and curvature of first four segments
	Env& segments(Tp la, Tp ca, Tp lb, Tp cb, Tp lc, Tp cc, Tp ld, Tp cd){
		Tp l[]={la,lb,lc,ld}; Tp c[]={ca,cb,cc,cd}; return segments(l,c,4); }


	/// Get break-point levels array
	Tv * levels(){ return mLevels; }

	/// Set break-point values
	template <class V>
	Env& levels(const V* vals, int len){
		int n = len <= size() ? len : size()+1;
		for(int i=0; i<n; ++i) levels()[i] = vals[i];
		return *this;
	}
	
	/// Set first two break-point levels
	Env& levels(Tv a, Tv b){ Tv v[]={a,b}; return levels(v,2); }

	/// Set first three break-point levels
	Env& levels(Tv a, Tv b, Tv c){ Tv v[]={a,b,c}; return levels(v,3); }

	/// Set first four break-point levels
	Env& levels(Tv a, Tv b, Tv c, Tv d){ Tv v[]={a,b,c,d}; return levels(v,4); }

	/// Set first five break-point levels
	Env& levels(Tv a, Tv b, Tv c, Tv d, Tv e){ Tv v[]={a,b,c,d,e}; return levels(v,5); }


	Env& maxLevel(Tv v){
		using namespace gam::scl;
		Tv mx(0);
		for(int i=0; i<N+1; ++i) mx = max(abs(mLevels[i]), mx);
		v = v/mx;
		for(int i=0; i<N+1; ++i) mLevels[i] *= v;
		return *this;
	}

protected:
	Curve<Tv,Tp> mCurve;
	Tp mLengths[N];			// segment lengths, in samples
	Tp mCurves[N];			// segment curvatures
	Tv mLevels[N+1];		// break-point levels

	uint32_t mPos, mLen;	// position in and length of current segment, in samples
	int mStage;				// the current curve segment
	int mSustain;			// index of sustain point
	int mLoop;

	void setLen(int i){ mLen=mLengths[i]*Ts::spu(); }
};



/// AD (Attack, Decay) envelope

/// \tparam Tv	value (sample) type
/// \tparam Tp	parameter type
/// \tparam Ts	sync type
/// \ingroup Envelopes 
template <class Tv=real, class Tp=real, class Ts=Synced>
class AD : public Env<2,Tv,Tp,Ts>{
public:
	using Env<2,Tv,Tp,Ts>::release;

	/// \param[in] att	Attack length
	/// \param[in] dec	Decay length
	/// \param[in] amp	Amplitude
	/// \param[in] crv	Curvature of all segments
	AD(Tp att =Tp(0.01), Tp dec =Tp(0.1), Tv amp = Tv(1), Tp crv =Tp(-4))
	{
		attack(att).decay(dec);
		this->levels(0,amp,0);
		this->curve(crv);
	}

	/// Set attack length
	AD& attack(Tp len){ return setLen(0,len); }

	/// Set decay length
	AD& decay(Tp len){ return setLen(1,len); }

	/// Set amplitude
	AD& amp(Tv v){ this->levels()[1]=v; return *this; }
	
protected:
	AD& setLen(int i, Tp v){
		this->lengths()[i] = v; return *this;
	}
};


/// ADSR (Attack, Decay, Sustain, Release) envelope

/// This is a three segment envelope that rises to one and then falls back
/// to zero. The attack is the rise length, the decay is the length until
/// hitting the sustain level, and the release is the length until hitting
/// zero again.
/// This envelope is most useful when the duration of the envelope is not known 
/// in advance. However, it can be easily converted into a fixed length ADR by
/// calling the sustainDisable() method. The decay segment then acts as a pseudo
/// steady state portion.
///
/// \tparam Tv	value (sample) type
/// \tparam Tp	parameter type
/// \tparam Ts	sync type
/// \ingroup Envelopes 
template <class Tv=real, class Tp=real, class Ts=Synced>
class ADSR : public Env<3,Tv,Tp,Ts>{
public:
	using Env<3,Tv,Tp,Ts>::release;

	/// \param[in] att	Attack length
	/// \param[in] dec	Decay length
	/// \param[in] sus	Sustain level (as factor of amplitude)
	/// \param[in] rel	Release length
	/// \param[in] amp	Amplitude
	/// \param[in] crv	Curvature of all segments
	ADSR(
		Tp att =Tp(0.01), Tp dec =Tp(0.1), Tv sus =Tv(0.7), Tp rel =Tp(1.),
		Tv amp =Tv( 1),
		Tp crv =Tp(-4)
	)
	{
		this->sustainPoint(2);
		this->levels(0,amp,sus*amp,0);
		attack(att).decay(dec).release(rel);
		this->curve(crv);
	}

	/// Set attack length
	ADSR& attack(Tp len){ return setLen(0,len); }

	/// Set decay length
	ADSR& decay(Tp len){ return setLen(1,len); }

	/// Set sustain level
	ADSR& sustain(Tv val){
		this->levels()[2] = val * this->levels()[1];
		return *this;
	}

	/// Set release length
	ADSR& release(Tp len){ return setLen(2,len); }
	
	/// Set amplitude
	ADSR& amp(Tv v){ return this->maxLevel(v); }
	
protected:
	ADSR& setLen(int i, Tp v){
		this->lengths()[i] = v; return *this;
	}
};



/// Exponentially decaying curve

/// This envelope exponentially decays towards zero starting from an initial
/// value. Because zero is never reached, the decay length determines when the
/// envelope is -60 dB down from its initial value. This envelope is one of the 
/// most computationally efficient envelopes requiring only a single multiply
/// per iteration. \n\n Compare to Curve which touches exactly both start and end points.
/// \ingroup Envelopes
/// \sa Curve
template <class T=real, class Ts=Synced>
class Decay : public Ts{
public:

	/// \param[in] decay	Number of units until initial value decays -60 dB
	/// \param[in] val		Intial value
	Decay(T decay=T(1), T val=T(1));

	T decay() const;		///< Returns -60 dB decay length
	bool done(T thresh=T(0.001)) const; ///< Returns whether value is below threshold
	T value() const;		///< Returns current value

	T operator()();			///< Generate next sample
	
	void decay(T val);		///< Set number of units for curve to decay -60 dB
	void reset();			///< Set current value to 1
	void value(T val);		///< Set current value
	
protected:
	T mVal, mMul, mDcy;

	virtual void onResync(double r);
};



/// Binary gate controlled by threshold comparison. The gate closes if the (norm of the) input value remains below threshold for at least closingDelay units of time (typically seconds).
    
/// \ingroup Envelopes
template <class T=real, class Ts=Synced>
class Gate : public Ts{
public:

	/// \param[in] closingDelay		units to wait before closing while under threshold
	/// \param[in] threshold		threshold below which gate closes
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
    
/// \ingroup Envelopes
template <
	class Tv=real,
	template <class> class Si=iplSeq::Linear,
	class Tp=real,
	class Ts=Synced
>
class Seg : public Ts{
public:

	/// \param[in] len		Length of segment in domain units
	/// \param[in] start	Start value
	/// \param[in] end		End value
	/// \param[in] phase	Start phase along segment, in [0,1)
	Seg(Tp len=0.5, Tv start=1, Tv end=0, Tp phase=0):
		mFreq((Tp)1/len), mAcc(0, phase), mIpl(start)
	{
		mIpl.push(end);
		Ts::initSynced();
	}
	
	/// Returns whether envelope is done
	bool done() const { return mAcc.val >= Tp(1); }

	/// Get current value
	Tv val() const { return mIpl(scl::min(mAcc.val, Tp(1))); }


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
	
	/// Set length, in domain units
	void length(Tp v){ freq((Tp)1/v); }

	/// Set length, in domain units
	void period(Tp v){ length(v); }

	/// Set phase along segment
	void phase(Tp v){ mAcc = v; }

	/// Reset envelope
	void reset(){ phase((Tp)0); }

	
	Si<Tv>& ipol(){ return mIpl; }
	
protected:
	Tp mFreq;
	gen::RAdd<Tp> mAcc;
	Si<Tv> mIpl;
	
	virtual void onResync(double r){ freq(mFreq); }
};




/// Exponential envelope segment for smoothing out value changes.
    
/// \ingroup Envelopes
template <class T=gam::real, class Ts=Synced>
class SegExp : public Ts{
public:

	/// \param[in] len		Length of segment in domain units
	/// \param[in] crv		Curvature of segment
	/// \param[in] start	Start value
	/// \param[in] end		End value
	SegExp(T len, T crv=-3, T start=1, T end=0):
		mLen(len), mCrv(crv), mVal1(start), mVal0(end)
	{
		Ts::initSynced();
	}
	
	/// Returns whether envelope is done
	bool done() const { return mCurve.value() >= T(1); }
	
	/// Generate next value
	T operator()(){
		if(done()) return mVal0;
		return ipl::linear(scl::min(mCurve(), T(1)), mVal1, mVal0);
	}
	
	/// Set new end value.  Start value is set to current value.
	void operator= (T v){
		mVal1 = ipl::linear(scl::min(mCurve.value(), T(1)), mVal1, mVal0);
		mVal0 = v;
		mCurve.reset();
	}
	
	/// Set curvature.  Negative gives faster change, positive gives slower change.
	void curve(T v){ set(mLen, v); }
	
	/// Set length in domain units.
	void period(T v){ set(v, mCrv); }

	void reset(){ mCurve.reset(); }

	/// Set length and curvature
	void set(T len, T crv){
		mLen = len; mCrv = crv;
		mCurve.set(len * Ts::spu(), crv);
	}
	
	virtual void onResync(double r){ set(mLen, mCrv); }
	
protected:
	T mLen, mCrv, mVal1, mVal0;
	Curve<T,T> mCurve;
};





// Implementation_______________________________________________________________


//---- Curve
template <class Tv,class Tp>
Curve<Tv,Tp>::Curve(): mEnd(Tv(1)), mA(Tv(0)), mB(Tv(0)), mMul(Tp(1))
{}

template <class Tv,class Tp>
Curve<Tv,Tp>::Curve(Tp length, Tp curve, Tv start, Tv end){
	set(length, curve, start, end);
}

template <class Tv,class Tp>
inline bool Curve<Tv,Tp>::done() const {
	Tv dv = mB - mB*mMul; // linear apx of derivative
	if(dv > Tv(0))	return value() >= end();
	else			return value() <= end();
}

template <class Tv,class Tp>
inline Tv Curve<Tv,Tp>::value() const { return mA - mB; }

// dividing by mMul goes back one step
template <class Tv,class Tp>
inline Curve<Tv,Tp>& Curve<Tv,Tp>::reset(Tv start){ mB = (mA-start) / mMul; return *this; }

// hack to get proper max floating point value
namespace{
	template<class T> inline T	eps(){ return T(0.00001  ); }
	template<> inline double	eps(){ return   0.00000001; }
	template<class T> inline T	maxReal(){ return DBL_MAX; }
	template<> inline float		maxReal<float>(){ return FLT_MAX; }
}

template <class Tv,class Tp>
Curve<Tv,Tp>& Curve<Tv,Tp>::set(Tp len, Tp crv, Tv start, Tv end){
	static const Tp EPS = eps<Tp>();

	if(len == Tp(0)){ // if length is 0, return end value immediately
		mEnd = end;
		mMul = maxReal<Tp>();
		mA = end;
		mB = Tv(0);
		return *this;
	}

	// Avoid discontinuity when curve = 0 (a line)
	if(crv < EPS && crv > -EPS){
		crv = crv < Tp(0) ? -EPS : EPS;
	}
	
	Tp crvOverLen = crv / len;
	
	if(crvOverLen < EPS && crvOverLen > -EPS){
		crvOverLen = crvOverLen < Tp(0) ? -EPS : EPS;
		crv = crvOverLen * len;
	}

	/*
	This algorithm uses an exponential curve in [0,1] to linearly interpolate
	between the start and end values.
	
	    1 - e^(cx)
	y = ---------- * (end-start) + start
	     1 - e^c

	     delta             delta
	  = ------- + start - ------- e^(cx) 
	    1 - e^c           1 - e^c
	*/

	mEnd = end;
	mMul = ::exp(crvOverLen);
	mA = (end-start) / (Tp(1) - ::exp(crv));
	mB = mA / mMul;
	mA+= start;
	return *this;
}

template <class Tv,class Tp>
inline Curve<Tv,Tp>& Curve<Tv,Tp>::value(const Tv& v){ mB = mA-v; return *this; }

template <class Tv,class Tp>
inline Tv Curve<Tv,Tp>::operator()(){
	mB *= mMul;
	return value();
}

    

//---- Decay

template <class T, class Ts>
Decay<T,Ts>::Decay(T decay_, T val)
:	mVal(val)
{
	Ts::initSynced();
	decay(decay_);
}

template <class T, class Ts>
inline T Decay<T,Ts>::operator()(){ T o = mVal; mVal *= mMul; return o; }

template <class T, class Ts>
inline void Decay<T,Ts>::decay(T v){
	mDcy = v;
	mMul = scl::t60(v * Ts::spu());
}

template <class T, class Ts>
inline void Decay<T,Ts>::reset(){ mVal = 1; }
    
template <class T, class Ts>
inline void Decay<T,Ts>::value(T v){ mVal = v; }

template <class T, class Ts>
inline T Decay<T,Ts>::decay() const { return mDcy; }
    
template <class T, class Ts>
inline bool Decay<T,Ts>::done(T thr) const { return mVal < thr; }
    
template <class T, class Ts>
inline T Decay<T,Ts>::value() const { return mVal; }

template <class T, class Ts>
void Decay<T,Ts>::onResync(double r){ decay(mDcy); }


} // gam::
#endif
