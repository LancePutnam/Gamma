#ifndef GAMMA_ENVELOPE_H_INC
#define GAMMA_ENVELOPE_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include <cfloat> /* DBL_MAX, FLT_MAX */
#include "Gamma/gen.h"
#include "Gamma/ipl.h"
#include "Gamma/scl.h"
#include "Gamma/Domain.h"
#include "Gamma/Strategy.h"

namespace gam{

/// Non-periodic, slowly varying modulation sources

/// \defgroup Envelope


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
/// \ingroup Envelope
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
/// \ingroup Envelope 
template <int N, class Tv=real, class Tp=real, class Td=GAM_DEFAULT_DOMAIN>
class Env : public Td{
public:

	Env();

	Env(Tp lvl1, Tp len1, Tp lvl2);

	Env(Tp lvl1, Tp len1, Tp lvl2, Tp len2, Tp lvl3);

	Env(Tp lvl1, Tp len1, Tp lvl2, Tp len2, Tp lvl3, Tp len3, Tp lvl4);


	/// Get the number of segments
	int size() const { return N; }
	
	/// Get the position, in samples, within the current segment
	unsigned position() const { return mPos; }
	
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
	Tv operator()();

	/// Release the envelope
	void release();


	/// Set whether envelope loops

	/// Note that when this is activated the last segment moves from the second
	/// to last level to the first level. Thus, the very last level is ignored.
	/// This is done to avoid clicking due to mismatched levels when the
	/// envelope wraps around.
	Env& loop(bool v){ mLoop=v; return *this; }

	/// Sets the point at which the envelope holds its value until released
	Env& sustainPoint(int v){ mSustain=v; return *this; }

	/// Disable sustain
	Env& sustainDisable(){ return sustainPoint(N); }

	/// Reset envelope to starting point
	void reset();

	/// Reset envelope to starting point setting first level to current level
	void resetSoft();

	/// Jump to end of envelope
	void finish();

	
	/// Get segment lengths array
	Tp * lengths(){ return mLengths; }
	const Tp * lengths() const { return mLengths; }

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

	template <unsigned i>
	Env& length(Tp v){
		static_assert(i<N, "Invalid segment index");
		lengths()[i] = v;
		return *this;
	}

	/// Get total length of all envelope segments
	Tp totalLength() const;
	
	/// Set total length of envelope by adjusting one segment length
	
	/// \param[in] length		desired length
	/// \param[in] modSegment	segment whose length is modified to match desired length
	Env& totalLength(Tp length, int modSegment);

	/// Set total length of envelope by scaling all segment lengths
	Env& totalLength(Tp length);


	/// Get segment curvature array
	Tp * curves(){ return mCurves; }
	const Tp * curves() const { return mCurves; }
	
	/// Set curvature of all segments
	Env& curve(Tp v){
		for(int i=0; i<N; ++i) curves()[i]=v;
		return *this;
	}

	/// Set segment curve amounts
	template <class V>
	Env& curves(const V* vals, int len){
		int n = len < size() ? len : size();
		for(int i=0; i<n; ++i) mCurves[i] = vals[i];
		return *this;
	}

	/// Set curvature of first two segments
	Env& curves(Tp ca, Tp cb){
		Tp c[]={ca,cb}; return curves(c,2); }

	/// Set curvature of first three segments
	Env& curves(Tp ca, Tp cb, Tp cc){
		Tp c[]={ca,cb,cc}; return curves(c,3); }

	/// Set curvature of first four segments
	Env& curves(Tp ca, Tp cb, Tp cc, Tp cd){
		Tp c[]={ca,cb,cc,cd}; return curves(c,4); }

	template <unsigned i>
	Env& curve(Tp v){
		static_assert(i<N, "Invalid curve index");
		curves()[i] = v;
		return *this;
	}


	/// Set length and curvature of a segment
	Env& segment(int i, Tp len, Tp crv){
		mLengths[i]=len;
		mCurves [i]=crv;
		return *this;
	}

	template <unsigned i>
	Env& segment(Tp len, Tp crv){
		static_assert(i<N, "Invalid segment index");
		return segment(i, len,crv);
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
	const Tv * levels() const { return mLevels; }

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

	template <unsigned i>
	Env& level(Tv v){
		static_assert(i<(N+1), "Invalid level index");
		levels()[i] = v;
		return *this;
	}

	/// Set maximum level
	Env& maxLevel(Tv v);

protected:
	Curve<Tv,Tp> mCurve;
	Tp mLengths[N];		// segment lengths, in samples
	Tp mCurves[N];		// segment curvatures
	Tv mLevels[N+1];	// break-point levels

	unsigned mPos, mLen;// position in and length of current segment, in samples
	int mStage;			// the current curve segment
	int mSustain=N;		// index of sustain point
	int mLoop=0;

	void setLen(int i){ mLen=unsigned(mLengths[i]*Td::spu()); }
};



/// AD (Attack, Decay) envelope

/// \tparam Tv	value (sample) type
/// \tparam Tp	parameter type
/// \tparam Td	domain observer type
/// \ingroup Envelope 
template <class Tv=real, class Tp=real, class Td=GAM_DEFAULT_DOMAIN>
class AD : public Env<2,Tv,Tp,Td>{
public:
	using Env<2,Tv,Tp,Td>::release;

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


	/// Get attack length
	Tp attack() const { return this->lengths()[0]; }

	/// Get decay length
	Tp decay() const { return this->lengths()[1]; }

	/// Get amplitude
	Tv amp() const { return this->levels()[1]; }

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
/// \tparam Td	domain observer type
/// \ingroup Envelope 
template <class Tv=real, class Tp=real, class Td=GAM_DEFAULT_DOMAIN>
class ADSR : public Env<3,Tv,Tp,Td>{
public:
	using Env<3,Tv,Tp,Td>::release;

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
	ADSR& amp(Tv v){ this->maxLevel(v); return *this; }


	/// Get attack length
	Tp attack() const { return this->lengths()[0]; }

	/// Get decay length
	Tp decay() const { return this->lengths()[1]; }

	/// Get sustain level
	Tv sustain() const { return this->levels()[2] / this->levels()[1]; }

	/// Get release length
	Tp release() const { return this->lengths()[2]; }

	/// Get amplitude
	Tv amp() const { return this->levels()[1]; }

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
/// per iteration. That said, after a certain number iterations, denormalized
/// floats may be generated which can dramatically decrease performance. To
/// avoid denormals get the next sample via: decay.done() ? 0.f : decay().
/// \n\n Compare to Curve which touches exactly both start and end points.
/// \ingroup Envelope
/// \sa Curve
template <class T=real, class Td=GAM_DEFAULT_DOMAIN>
class Decay : public Td{
public:

	/// \param[in] decay	Number of units until initial value decays -60 dB
	/// \param[in] val		Intial value
	Decay(T decay=T(1), T val=T(1));

	T decay() const;		///< Returns -60 dB decay length
	bool done(T thresh=T(0.001)) const; ///< Returns whether value is below threshold
	T value() const;		///< Returns current value

	T operator()();			///< Generate next sample
	
	void decay(T v);		///< Set number of units for curve to decay -60 dB

	void value(T v);		///< Set current value

	void reset(T amp=T(1));	///< Reset envelope and assign amplitude
	void finish(T amp=T(0.001)); ///< Jump to end of envelope

	void onDomainChange(double r);

protected:
	T mVal, mMul, mDcy;
};



/// Binary gate controlled by threshold comparison

/// The gate returns 1 if the input magnitude is above a specified threshold,
/// otherwise it returns 0. A closing delay can also be specified to determine
/// the window of time the input must be below the threshold before the gate
/// closes. This is equivalent to what is known in electronics as a comparator
/// with hysteresis.
///
/// \ingroup Envelope
template <class T=real, class Td=GAM_DEFAULT_DOMAIN>
class Gate : public Td{
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
			mRemain -= Td::ups();
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

/// \ingroup Envelope Interpolation
template <
	class Tv=real,
	template <class> class Si=iplSeq::Linear,
	class Tp=real,
	class Td=GAM_DEFAULT_DOMAIN
>
class Seg : public Td{
public:

	/// \param[in] len		Length of segment in domain units
	/// \param[in] start	Start value
	/// \param[in] end		End value
	/// \param[in] phase	Start phase along segment, in [0,1)
	Seg(Tp len=0.5, Tv start=1, Tv end=0, Tp phase=0):
		mFreq((Tp)1/len), mAcc(0, phase), mIpl(start)
	{
		mIpl.push(end);
		onDomainChange(1);
	}


	/// Returns whether envelope is done
	bool done() const { return mAcc.val >= Tp(1); }

	/// Get current value
	Tv val() const { return mIpl(scl::min(mAcc.val, Tp(1))); }

	
	/// Set new end value.  Start value is set to current value.
	void operator= (Tv v){
		mIpl.val(mIpl(scl::min(mAcc.val, Tp(1))));
		mIpl.push(v);
		reset();
	}

	/// Generate next value
	Tv operator()(){
		if(done()) return mIpl.val();
		Tp f = mAcc.val;
		mAcc();
		return mIpl(f);
	}

	/// Generates a new end point from a generator when the segment end is reached
	
	/// This can be used to upsample and interpolate a lower-rate signal,
	/// e.g., creating pitched noise from a random number generator.
	template <class Gen>
	Tv operator()(Gen& g){
		if(done()){
			mIpl.push(g());
			mAcc.val = mAcc.val - Tp(1); // wrap phase
		}
		Tp f = mAcc.val;
		mAcc();
		return mIpl(f);
	}
	
	/// Set frequency of envelope
	void freq(Tp v){ mFreq = v; mAcc.add = v * Tp(Td::ups()); }
	
	/// Set length, in domain units
	void length(Tp v){ freq(Tp(1)/v); }

	/// Set length, in domain units
	void period(Tp v){ length(v); }

	/// Set phase along segment
	void phase(Tp v){ mAcc = v; }

	/// Reset envelope
	void reset(){ phase(Tp(0)); }


	Si<Tv>& ipol(){ return mIpl; }

	void onDomainChange(double /*r*/){ freq(mFreq); }

protected:
	Tp mFreq;
	gen::RAdd<Tp> mAcc;
	Si<Tv> mIpl;
};



/// Exponential envelope segment for smoothing out value changes.

/// \ingroup Envelope Interpolation
template <class T=gam::real, class Td=GAM_DEFAULT_DOMAIN>
class SegExp : public Td{
public:

	/// \param[in] len		Length of segment in domain units
	/// \param[in] crv		Curvature of segment
	/// \param[in] start	Start value
	/// \param[in] end		End value
	SegExp(T len, T crv=-3, T start=1, T end=0):
		mLen(len), mCrv(crv), mVal1(start), mVal0(end)
	{
		onDomainChange(1);
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
		mCurve.set(len * Td::spu(), crv);
	}
	
	void onDomainChange(double r){ set(mLen, mCrv); }
	
protected:
	T mLen, mCrv, mVal1, mVal0;
	Curve<T,T> mCurve;
};



// Implementation_______________________________________________________________

template <class Tv,class Tp>
Curve<Tv,Tp>::Curve()
:	mEnd(Tv(1)), mA(Tv(0)), mB(Tv(0)), mMul(Tp(1))
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

// dividing by mMul goes back one step
template <class Tv,class Tp>
Curve<Tv,Tp>& Curve<Tv,Tp>::reset(Tv start){
	mB = (mA-start) / mMul;
	return *this;
}

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
	mMul = exp(crvOverLen);
	mA = (end-start) / (Tp(1) - exp(crv));
	mB = mA / mMul;
	mA+= start;
	return *this;
}

template <class Tv,class Tp>
inline Tv Curve<Tv,Tp>::value() const {
	return mA - mB;
}

template <class Tv,class Tp>
inline Curve<Tv,Tp>& Curve<Tv,Tp>::value(const Tv& v){
	mB = mA-v;
	return *this;
}

template <class Tv,class Tp>
inline Tv Curve<Tv,Tp>::operator()(){
	mB *= mMul;
	return value();
}



template <int N,class Tv,class Tp,class Td>
Env<N,Tv,Tp,Td>::Env(){
	for(int i=0; i<N; ++i){
		mLengths[i]= Tp(1e-8);
		mCurves[i] = Tp(-4);
		mLevels[i] = Tv();
	}
	mLevels[N] = Tv();
	reset();
}

template <int N,class Tv,class Tp,class Td>
Env<N,Tv,Tp,Td>::Env(Tp lvl1, Tp len1, Tp lvl2){
	levels(lvl1,lvl2);
	lengths()[0] = len1;
	curve(-4);
	reset();
}

template <int N,class Tv,class Tp,class Td>
Env<N,Tv,Tp,Td>::Env(Tp lvl1, Tp len1, Tp lvl2, Tp len2, Tp lvl3){
	levels(lvl1,lvl2,lvl3);
	lengths(len1,len2);
	curve(-4);
	reset();
}

template <int N,class Tv,class Tp,class Td>
Env<N,Tv,Tp,Td>::Env(Tp lvl1, Tp len1, Tp lvl2, Tp len2, Tp lvl3, Tp len3, Tp lvl4){
	levels(lvl1,lvl2,lvl3,lvl4);
	lengths(len1,len2,len3);
	curve(-4);
	reset();
}

template <int N,class Tv,class Tp,class Td>
inline Tv Env<N,Tv,Tp,Td>::operator()(){

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
			int nextStage = mStage+1;
			// If looping, ensure we wrap back around to first level
			if(mLoop && (nextStage==size())) nextStage = 0;
			mCurve.set(Tp(mLen), mCurves[mStage], mLevels[mStage], mLevels[nextStage]);

			// Immediately return start level of new stage
			return (*this)();
		}
	}

	// Envelope is done:
	return mLevels[mStage];
}	

template <int N,class Tv,class Tp,class Td>
void Env<N,Tv,Tp,Td>::release(){

	if(released()) return;

	// begin release portion immediately starting at current level
	Tv curVal = value();

	/* pre-sustain release?
	auto oldStage = mStage;
	if(mStage < mSustain-1){
		//if(mSustain >= 1) mStage = mSustain-1;
		// go back until non-zero length segment is found...
		if(mSustain >= 1){
			mStage = mSustain-1;
			while(mStage > 0){
				if(mLengths[mStage] > 0.) break;
				--mStage;
			}
		}
		//printf("stage set to %d\n", mStage);
	}//*/

	mSustain = -mSustain; // neg. sustain means released

	mStage = -mSustain; // old version

	if(!done() /*&& oldStage != mStage*/){
		mPos = 0;
		setLen(mStage);
		mCurve.set(Tp(mLen), mCurves[mStage], curVal, mLevels[mStage+1]);
	}
}

template <int N,class Tv,class Tp,class Td>
void Env<N,Tv,Tp,Td>::reset(){
	// this forces a stage increment upon first iteration
	mPos = 0;//0xFFFFFFFF;
	mLen = 0;
	mStage = -1;
	mSustain = scl::abs(mSustain);
}

template <int N,class Tv,class Tp,class Td>
void Env<N,Tv,Tp,Td>::resetSoft(){
	Tv curVal = value();
	mPos = 0;
	mStage = 0;
	setLen(mStage);
	mCurve.set(Tp(mLen), mCurves[mStage], curVal, mLevels[mStage+1]);
	mSustain = scl::abs(mSustain);
}

template <int N,class Tv,class Tp,class Td>
void Env<N,Tv,Tp,Td>::finish(){
	mStage = size();
	mPos = mLen;
	mSustain = scl::abs(mSustain);
}

template <int N,class Tv,class Tp,class Td>
Tp Env<N,Tv,Tp,Td>::totalLength() const {
	Tp sum=Tp(0);
	for(int i=0;i<size();++i) sum += mLengths[i];
	return sum;
}

template <int N,class Tv,class Tp,class Td>
Env<N,Tv,Tp,Td>& Env<N,Tv,Tp,Td>::totalLength(Tp length, int modSegment){
	mLengths[modSegment] = Tp(0);
	mLengths[modSegment] = length - totalLength();
	return *this;
}

template <int N,class Tv,class Tp,class Td>
Env<N,Tv,Tp,Td>& Env<N,Tv,Tp,Td>::totalLength(Tp length){
	Tp mul = length / totalLength();
	for(int i=0; i<size(); ++i){
		lengths()[i] *= mul;
	}
	return *this;
}

template <int N,class Tv,class Tp,class Td>
Env<N,Tv,Tp,Td>& Env<N,Tv,Tp,Td>::maxLevel(Tv v){
	Tv mx(0);
	for(int i=0; i<N+1; ++i) mx = scl::max(scl::abs(mLevels[i]), mx);
	v = v/mx;
	for(int i=0; i<N+1; ++i) mLevels[i] *= v;
	return *this;
}



template <class T, class Td>
Decay<T,Td>::Decay(T decay_, T val)
:	mVal(val)
{
	onDomainChange(1);
	decay(decay_);
}

template <class T, class Td>
inline T Decay<T,Td>::operator()(){
	T o = mVal;
	mVal *= mMul;
	return o;
}

template <class T, class Td>
void Decay<T,Td>::decay(T v){
	mDcy = v;
	mMul = scl::t60(v * Td::spu());
}

template <class T, class Td>
void Decay<T,Td>::reset(T amp){ mVal = amp; }

template <class T, class Td>
void Decay<T,Td>::finish(T amp){ mVal = amp; }

template <class T, class Td>
void Decay<T,Td>::value(T v){ mVal = v; }

template <class T, class Td>
T Decay<T,Td>::decay() const { return mDcy; }

template <class T, class Td>
inline bool Decay<T,Td>::done(T thr) const { return mVal < thr; }

template <class T, class Td>
T Decay<T,Td>::value() const { return mVal; }

template <class T, class Td>
void Decay<T,Td>::onDomainChange(double /*r*/){ decay(mDcy); }


} // gam::
#endif
