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
	Curve(Tp length, Tp curve, Tv start=Tv(0), Tv end=Tv(1));

	bool done() const;					///< Returns whether curve has gone past end value
	Tv start() const { return mStart; }	///< Get start value
	Tv end() const { return mEnd; }		///< Get end value
	Tv value() const;					///< Get current value

	Tv operator()();					///< Generates next value
	Curve& reset();						///< Reset envelope
	Curve& value(const Tv& v);			///< Set value

	/// Set length and curvature
	
	/// \param[in] length	length of curve in samples
	/// \param[in] curve	curvature; pos. approaches slowly, neg. approaches rapidly, 0 approaches linearly
	/// \param[in] start	start value
	/// \param[in] end		end value
	Curve& set(Tp length, Tp curve, Tv start=Tv(0), Tv end=Tv(1));

	/// Set start and end levels
	Curve& levels(Tv start, Tv end);

protected:
	Tv mStart, mEnd, mA, mB;
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

	/// Set envelope as constant value with all lengths equal to zero
	Env(Tv lvl = Tv());

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

	template <int i>
	Env& sustainPoint(){
		static_assert(0<=i && i<=N, "Invalid sustain index");
		return sustainPoint(i);
	}

	/// Disable sustain
	Env& sustainDisable(){ return sustainPoint<N>(); }

	/// Reset envelope to starting point
	void reset();

	/// Reset envelope to starting point setting first level to current level
	void resetSoft();

	/// Jump to end of envelope
	void finish();

	
	/// Get segment lengths array
	Tp * lengths(){ return mLengths; }
	const Tp * lengths() const { return mLengths; }

	/// Get length of a segment (compile-time checked)
	template <unsigned i>
	const Tp& length() const { return at<i>(mLengths); }

	/// Set length of a segment (compile-time checked)
	template <unsigned i>
	Env& length(Tp v){ at<i>(mLengths) = v; return *this; }

	/// Set lengths of segments
	template <class V>
	Env& lengths(const V* vals, int len){
		int n = len <= size() ? len : size();
		for(int i=0; i<n; ++i) lengths()[i] = vals[i];
		return *this;
	}

	/// Set lengths of first two segments
	Env& lengths(Tp a, Tp b){ Tp v[]={a,b}; return lengths(v,2); }

	/// Set lengths of first three segments
	Env& lengths(Tp a, Tp b, Tp c){ Tp v[]={a,b,c}; return lengths(v,3); }

	/// Set lengths of first four segments
	Env& lengths(Tp a, Tp b, Tp c, Tp d){ Tp v[]={a,b,c,d}; return lengths(v,4); }

	/// Set lengths of first five segments
	Env& lengths(Tp a, Tp b, Tp c, Tp d, Tp e){ Tp v[]={a,b,c,d,e}; return lengths(v,5); }


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

	/// Get curvature of a segment (compile-time checked)
	template <unsigned i>
	const Tp& curve() const { return at<i>(mCurves); }

	/// Set curvature of a segment (compile-time checked)
	template <unsigned i>
	Env& curve(Tp v){ at<i>(mCurves) = v; return *this; }

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


	/// Set length and curvature of a segment
	Env& segment(int i, Tp len, Tp crv){
		mLengths[i]=len;
		mCurves [i]=crv;
		return *this;
	}

	/// Set length and curvature of a segment (compile-time checked)
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

	/// Get level at index (compile-time checked)
	template <unsigned i>
	const Tv& level() const { return at<i>(mLevels); }

	/// Set level at index (compile-time checked)
	template <unsigned i>
	Env& level(Tv v){ at<i>(mLevels) = v; return *this; }

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

	/// Set maximum level
	Env& maxLevel(Tv v);

	/// Configure as hold-attack-decay envelope
	Env& setHAD(Tp h, Tp a, Tp d, Tv amp = Tv(1)){
		static_assert(N==3, "Requires three segments");
		levels(0,0,amp,0);
		lengths(h,a,d);
		return *this;
	}

	/// Push new segment onto end and pop the first segment

	/// This allows for a kind of "streaming" envelope where new segments are
	/// added to the end while segments at the start get ejected. Adding a new
	/// segment will in most cases not interrupt the envelope. If currently in 
	/// the first (popped) segment, then a soft reset is performed. The number 
	/// of segments that can be pushed without popping the current segment is
	/// equal to the envelope stage.
	/// \param[in] len		Length of new segment
	/// \param[in] lvl		End level of new segment
	/// \param[in] crv		Curvature of new segment
	Env& pushSegment(Tp len, Tv lvl, Tv crv);

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

	static Tp curveDefault(){ return {-4}; }

	template <unsigned Idx, class T, int M>
	static T& at(T (&arr)[M]){
		static_assert(Idx < M, "Index out of bounds");
		return arr[Idx];
	}
	template <unsigned Idx, class T, int M>
	static const T& at(const T (&arr)[M]){
		return at<Idx>(const_cast<T(&)[M]>(arr));
	}
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
	AD& attack(Tp v){ this->template length<0>(v); return *this; }

	/// Set decay length
	AD& decay(Tp v){ this->template length<1>(v); return *this; }

	/// Set amplitude
	AD& amp(Tv v){ this->template level<1>(v); return *this; }


	/// Get attack length
	Tp attack() const { return this->template length<0>(); }

	/// Get decay length
	Tp decay() const { return this->template length<1>(); }

	/// Get amplitude
	Tv amp() const { return this->template level<1>(); }

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
		this->template sustainPoint<2>();
		this->levels(0,amp,sus*amp,0);
		attack(att).decay(dec).release(rel);
		this->curve(crv);
	}

	/// Set attack length
	ADSR& attack(Tp v){ this->template length<0>(v); return *this; }

	/// Set decay length
	ADSR& decay(Tp v){ this->template length<1>(v); return *this; }

	/// Set sustain level
	ADSR& sustain(Tv v){ this->template level<2>(v * amp()); return *this; }

	/// Set release length
	ADSR& release(Tp v){ this->template length<2>(v); return *this; }
	
	/// Set amplitude
	ADSR& amp(Tv v){ this->maxLevel(v); return *this; }


	/// Get attack length
	Tp attack() const { return this->template length<0>(); }

	/// Get decay length
	Tp decay() const { return this->template length<1>(); }

	/// Get sustain level
	Tv sustain() const { return this->template level<2>() / amp(); }

	/// Get release length
	Tp release() const { return this->template length<2>(); }

	/// Get amplitude
	Tv amp() const { return this->template level<1>(); }
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

	/// \param[in] decay_	Number of units until initial value decays -60 dB
	/// \param[in] val		Intial value
	Decay(T decay_=T(1), T val=T(1))
	:	mVal(val)
	{
		onDomainChange(1);
		decay(decay_);
	}

	/// Set number of units for curve to decay -60 dB
	Decay& decay(T v){
		mDcy = v;
		mMul = scl::t60(v * Td::spu());
		return *this;
	}

	/// Set current value
	Decay& value(T v){ mVal = v; return *this; }

	/// Reset envelope and assign amplitude
	Decay& reset(T amp=T(1)){ return value(amp); }

	/// Jump to end of envelope
	Decay& finish(T amp=T(0.001)){ return value(amp); }


	/// Generate next sample
	T operator()(){
		T o = mVal;
		mVal *= mMul;
		return o;
	}

	/// Returns -60 dB decay length
	T decay() const { return mDcy; }

	/// Returns whether value is below threshold
	bool done(T thresh=T(0.001)) const { return mVal < thresh; }

	/// Returns current value
	T value() const { return mVal; }


	void onDomainChange(double /*r*/){ decay(mDcy); }

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


	/// Set closing delay
	Gate& delay(double v){ mDelay=mRemain=v; return *this; }

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


	/// Check whether gate is closed
	bool done() const { return mClosed; }

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
	Seg(Tp len=0.5, Tv start=1, Tv end=0, Tp phase=0){
		length(len);
		levels(start, end);
		this->phase(phase);
	}


	/// Set frequency of envelope
	Seg& freq(Tp v){
		mFreq = scl::abs(v);
		mAcc.add = v * Tp(Td::ups());
		return *this;
	}

	/// Set length, in domain units
	Seg& length(Tp v){
		return v > Tp(Td::ups()) ? freq(Tp(1)/v) : freq(Td::spu());
	}

	/// Set length, in domain units
	Seg& period(Tp v){ return length(v); }

	/// Set phase along segment
	Seg& phase(Tp v){ mAcc = v; return *this; }

	/// Reset envelope
	Seg& reset(){ return phase(Tp(0)); }

	/// Set start and end values and reset
	Seg& levels(Tv start, Tv end){
		mIpl.val(start);
		mIpl.push(end);
		return reset();
	}

	/// Set maximum level
	Seg& maxLevel(Tv v){
		auto maxAbs = mIpl.maxAbs();
		if(maxAbs == Tv(0)) mIpl.set(v);
		else{
			auto mul = v/maxAbs;
			for(auto& val : mIpl.vals) val *= mul;
		}
		return *this;
	}

	/// Get maximum level
	Tv maxLevel() const { return mIpl.maxAbs(); }

	/// Set new end value (start value set to current value)
	Seg& operator= (Tv v){
		return levels(mIpl(scl::min(mAcc.val, Tp(1))), v);
	}

	/// Generate next value
	Tv operator()(){
		if(done()) return mIpl.val();
		return next();
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
		return next();
	}


	/// Returns whether envelope is done
	bool done() const { return mAcc.val >= Tp(1); }

	/// Get current value
	Tv value() const { return mIpl(scl::min(mAcc.val, Tp(1))); }

	Si<Tv>& ipol(){ return mIpl; }

	void onDomainChange(double /*r*/){ freq(mFreq); }

protected:
	Tp mFreq;
	gen::RAdd<Tp> mAcc;
	Si<Tv> mIpl;

	Tv next(){
		auto frac = mAcc.val;
		mAcc();
		return mIpl(frac);
	}
};



/// Exponential envelope segment for smoothing out value changes.

/// \ingroup Envelope Interpolation
template <
	class Tv=gam::real,
	class Tp=gam::real,
	class Td=GAM_DEFAULT_DOMAIN
>
class SegExp : public Td{
public:

	/// \param[in] len		Length of segment in domain units
	/// \param[in] crv		Curvature of segment
	/// \param[in] start	Start value
	/// \param[in] end		End value
	SegExp(Tp len, Tp crv=Tp(-3), Tv start=Tv(1), Tv end=Tv(0))
	:	mLen(len), mCrv(crv), mStart(start), mEnd(end)
	{
		onDomainChange(1);
	}


	/// Set curvature (negative for fast approach, positive for slow approach)
	SegExp& curve(Tp v){ return set(mLen, v); }

	/// Get curvature
	Tp curve() const { return mCrv; }
	/// Get curvature with negative sign (decelerating approach)
	Tp curveDecel() const { return -curveAccel(); }
	/// Get curvature with positive sign (accelerating approach)
	Tp curveAccel() const { return scl::abs(mCrv); }

	/// Set length in domain units
	SegExp& length(Tp v){ return set(v, mCrv); }
	SegExp& period(Tp v){ return length(v); }

	/// Set length and curvature

	/// It is more efficient to set both parameters simultaneously rather than
	/// via their individual setters. This call will reset the internal
	/// interpolator.
	SegExp& set(Tp len, Tp crv){
		mLen = len; mCrv = crv;
		mCurve.set(len * Td::spu(), crv, Tv(0), Tv(1));
		return *this;
	}

	/// Set start and end values
	SegExp& levels(Tv start, Tv end){
		mStart = start;
		mEnd = end;
		return *this;
	}

	/// Reset envelope
	SegExp& reset(){ mCurve.reset(); return *this; }

	/// Set new end value (start value set to current value)

	/// If the curve or length parameters are to be changed to accompany this
	/// call, then they should be called \e after this call since they will 
	/// reset the internal interpolator.
	SegExp& operator= (Tv v){
		mStart = value();
		mEnd = v;
		mCurve.reset();
		return *this;
	}

	/// Generate next value
	Tv operator()(){
		if(done()) return mEnd;
		return tween(mCurve());
	}

	/// Get current value
	Tv value() const { return tween(mCurve.value()); }

	/// Returns whether envelope is done
	bool done() const { return mCurve.value() >= Tp(1); }

	void onDomainChange(double r){ set(mLen, mCrv); }
	
protected:
	Tv mStart, mEnd;
	Tp mLen, mCrv;
	Curve<Tp,Tp> mCurve; // specialized version with levels (0,1)

	Tv tween(Tp f) const { return ipl::linear(scl::min(f, Tp(1)), mStart, mEnd); }
};



// Implementation_______________________________________________________________

template <class Tv,class Tp>
Curve<Tv,Tp>::Curve()
:	mStart(Tv(0)), mEnd(Tv(1)), mA(Tv(0)), mB(Tv(0)), mMul(Tp(1))
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
Curve<Tv,Tp>& Curve<Tv,Tp>::reset(){
	mB = (mA-mStart) / mMul;
	return *this;
}

// hack to get proper max floating point value
namespace detail{
	template<class T> inline T	eps(){ return T(0.00001  ); }
	template<> inline double	eps(){ return   0.00000001; }
	template<class T> inline T	maxReal(){ return DBL_MAX; }
	template<> inline float		maxReal<float>(){ return FLT_MAX; }
}

template <class Tv,class Tp>
Curve<Tv,Tp>& Curve<Tv,Tp>::set(Tp len, Tp crv, Tv start, Tv end){
	static const Tp EPS = detail::eps<Tp>();

	mStart = start;
	mEnd = end;

	if(len == Tp(0)){ // if length is 0, return end value immediately
		mMul = detail::maxReal<Tp>();
		mA = mEnd;
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

	mMul = exp(crvOverLen);
	mA = (mEnd-mStart) / (Tp(1)-exp(crv)) + mStart;
	return reset();
}

template <class Tv,class Tp>
Curve<Tv,Tp>& Curve<Tv,Tp>::levels(Tv start, Tv end){
	auto c = (mA-mStart) / (mEnd-mStart); // 1 / (1-e^c)
	mStart = start;
	mEnd = end;
	mA = (mEnd-mStart)*c + mStart;
	return reset();
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
Env<N,Tv,Tp,Td>::Env(Tv lvl){
	for(int i=0; i<N; ++i){
		mLengths[i]= Tp(1e-8);
		mCurves[i] = curveDefault();
		mLevels[i] = lvl;
	}
	at<N>(mLevels) = lvl;
	reset();
}

template <int N,class Tv,class Tp,class Td>
Env<N,Tv,Tp,Td>::Env(Tp lvl1, Tp len1, Tp lvl2){
	levels(lvl1,lvl2);
	at<0>(mLengths) = len1;
	curve(curveDefault());
	reset();
}

template <int N,class Tv,class Tp,class Td>
Env<N,Tv,Tp,Td>::Env(Tp lvl1, Tp len1, Tp lvl2, Tp len2, Tp lvl3){
	levels(lvl1,lvl2,lvl3);
	lengths(len1,len2);
	curve(curveDefault());
	reset();
}

template <int N,class Tv,class Tp,class Td>
Env<N,Tv,Tp,Td>::Env(Tp lvl1, Tp len1, Tp lvl2, Tp len2, Tp lvl3, Tp len3, Tp lvl4){
	levels(lvl1,lvl2,lvl3,lvl4);
	lengths(len1,len2,len3);
	curve(curveDefault());
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
		mLengths[i] *= mul;
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

template <int N,class Tv,class Tp,class Td>
Env<N,Tv,Tp,Td>& Env<N,Tv,Tp,Td>::pushSegment(Tp len, Tv lvl, Tv crv){
	for(int i=1; i<N; ++i){
		mLengths[i-1] = mLengths[i];
		mCurves[i-1] = mCurves[i];
	}
	for(int i=1; i<N+1; ++i) mLevels[i-1] = mLevels[i];
	length<N-1>(len);
	curve<N-1>(crv);
	level<N>(lvl);

	if(!mStage) resetSoft();
	else{ done() ? mStage-=2 : --mStage; }

	return *this;
}

} // gam::
#endif
