#ifndef GAMMA_DELAY_H_INC
#define GAMMA_DELAY_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include "Gamma/arr.h"
#include "Gamma/ipl.h"
#include "Gamma/mem.h"
#include "Gamma/scl.h"
#include "Gamma/tbl.h"

#include "Gamma/Containers.h"
#include "Gamma/Strategy.h"
#include "Gamma/Sync.h"
#include "Gamma/Types.h"

namespace gam{


/// Some filter types.
namespace Filter{
	enum { LP, HP, BP, BPC, BR, AP };
	
//	/// Returns filter type as human readable string.
//	const char * string(int type){
//		switch(type){
//		#define CS(type, string) case Filter::type: return #string;
//		CS(LP, Low Pass) CS(HP, High Pass) CS(BP, Band Pass) CS(BPC, Band Pass)
//		CS(BR, Band Reject) CS(AP, All Pass)
//		default: return "Unknown";
//		#undef CS
//		}
//	}
};



/// First-order all-pass filter.

/// This filter is essentially a "phase" filter in that its only effect is to 
/// shift frequencies from 0 to Nyquist from 0 to -180 degress, respectively.
/// The cutoff frequency is the frequency where the phase is shifted by -90 
/// degrees.  A lowpass or highpass filter can be constructed by adding or 
/// subtracting, respectively, the output of the filter to or from the input.
/// When the center frequency is fs/4, the filter acts as a unit delay.
/// When the center frequency is 0, the filter acts as a no-op.
/// When the center frequency is fs/2, the filter acts as an inverter.
template<class Tv=gam::real, class Tp=gam::real, class Ts=Synced>
class AllPass1 : public Ts {
public:
	///
	/// @param[in]	frq		Initial cutoff frequency.
	AllPass1(Tp frq=1000);

	void freq (Tp v);	///< Set cutoff frequency
	void freqF(Tp v);	///< Faster, but slightly less accurate than freq()	
	void zero(){ d1=Tv(0); }
	
	Tv operator()(Tv input);	///< Filters sample
	
	Tv high(Tv input);		///< High-pass filters sample
	Tv low (Tv input);		///< Low-pass filters sample
	
	Tp freq();				///< Get current cutoff frequency
	
	virtual void onResync(double r);
	
protected:
	Tv d1;		// once delayed value
	Tp c;		// feed coefficient
	Tp mFreq;
};


/// 2-pole/2-zero IIR filter.

/// The biquadratic filter contains 2 zeroes and 2 poles in its transfer
/// function. The zeroes provide a better response near the DC and Nyquist
/// frequencies than an all-pole filter would. Second-order IIRs have a 12 
/// dB/octave cutoff slope and are normally cascaded (run in series) to obtain
/// a sharper response. This particular implementation utilizes the Butterworth
/// design.
template <class Tv=gam::real, class Tp=gam::real, class Ts=Synced>
class Biquad : public Ts{
public:

	/// @param[in]	frq		Initial center frequency in Hz.
	/// @param[in]	res		Initial resonance.
	/// @param[in]	frq		Initial type of filter.
	Biquad(Tp frq = Tp(1000), Tp res = Tp(1), int filterType = Filter::LP);

	/// Set input (a's) and output (b's) coefficients directly
	void coef(Tp a0, Tp a1, Tp a2, Tp b0, Tp b1, Tp b2);

	void freq(Tp v);					///< Set center frequency. 
	void res(Tp v);						///< Set resonance.
	void set(Tp frq, Tp res);			///< Set filter center frequency and resonance.
	void set(Tp frq, Tp res, int type);	///< Set all filter params.
	void type(int type);				///< Set type of filter (see Filter::type)
	void zero();						///< Zero internal delays.

	Tv operator()(Tv i0);				///< Return next filter output.
	Tv nextBP(Tv i0);					///< Optimized for band-pass (BP, BPC) types.
	
	Tp freq() const { return mFreq; }	///< Get center frequency
	Tp res() const { return mRes; }		///< Get resonance
	int type() const { return mType; }	///< Get filter type
	
	virtual void onResync(double r);

protected:
	Tp mA0, mA1, mA2, mB0, mB1, mB2;	// ffd and fbk coefficients
	Tv d1, d2;		// inner sample delays
	Tp mFreq, mRes;	// center frequency, resonance
	int mType;
	Tp mReal, mImag;	// real, imag components of center frequency
	Tp mAlpha;
	Tp mFrqToRad;
};



/// DC frequency blocker
template <class Tv=gam::real, class Tp=gam::real, class Ts=Synced>
class BlockDC : public Ts{
public:

	/// @param[in] width	Bandwidth of pole
	BlockDC(Tp width=35): d1(0), mWidth(width){ Ts::initSynced(); }

	/// Filter sample
	Tv operator()(Tv i0){		
		i0 += d1*mB1; Tv o0 = i0-d1; d1 = i0; return o0;
	}

	/// Set bandwidth of pole
	void width(Tp v){
		mWidth = v; mB1 = scl::poleRadius(v, Ts::ups());
	}

	void zero(){ d1=0; }

	virtual void onResync(double r){ width(mWidth); }

protected:
	Tv d1; Tp mWidth, mB1;
};




/// Nyquist frequency blocker
template <class Tv=gam::real, class Tp=gam::real, class Ts=Synced>
class BlockNyq : public BlockDC<Tv,Tp,Ts>{
public:

	/// @param[in] width	Bandwidth of pole
	BlockNyq(Tp width=35): Base(width){}

	/// Filter sample
	Tv operator()(Tv i0){		
		i0 += d1*mB1; Tv o0 = i0-d1; d1 =-i0; return o0;
	}

protected:
	typedef BlockDC<Tv,Tp,Ts> Base;
	using Base::d1; using Base::mB1;
};




/// Variable length delay-line.
template <class Tv=gam::real, template<class> class Si=ipl::Linear, class Ts=Synced>
class Delay : public ArrayPow2<Tv>, Ts{
public:

	Delay();

	/// @param[in]	delay		Initial delay length
	/// The size of the delay buffer will be the smallest possible power of two.
	Delay(float delay);

	/// @param[in]	maxDelay	Maximum delay length
	/// @param[in]	delay		Initial delay length
	/// The size of the delay buffer will be the smallest possible power of two.
	Delay(float maxDelay, float delay);
	
	virtual ~Delay(){}

	void delay(float v);						///< Set delay length
	void delayUnit(float u);					///< Set delay as (0, 1) of buffer size
	void freq(float v);							///< Set natural frequency (1/delay())
	void ipolType(int v){ mIpol.type(v); }		///< Set interpolation type
	void maxDelay(float v);						///< Set maximum delay length
	void zero();								///< Sets all elements to zero

	Tv operator()(const Tv& v);					///< Returns next filtered value
	Tv operator()() const;						///< Reads delayed element from buffer
	Tv read(float ago);							///< Returns element 'ago' units ago
	void write(const Tv& v);					///< Writes element into buffer. Tap is post-incremented.
	void writePre(const Tv& v);					///< Writes element into buffer. Tap is pre-incremented.
	
	float delay() const;						///< Returns current delay length
	uint32_t delayIndex(uint32_t delay) const;	///< Returns index of delayed element
	float delayUnit() const;					///< Returns unit delay (to max delay)
	uint32_t indexBack() const;					///< Returns index of backmost element
	float maxDelay() const;						///< Returns maximum delay length units

	virtual void onResize();
	virtual void onResync(double r);

	void print();

protected:
	Si<Tv> mIpol;

	float mMaxDelay;				// maximum delay length
	float mDelayFactor;				// multiplication factor when setting delay
	float mDelayLength;				// current delay length
	uint32_t mPhase;				// write tap
	uint32_t mPhaseInc;				// phase increment
	uint32_t mDelay;				// read tap as delay from write tap

	void incPhase();				// increment phase
	void refreshDelayFactor();
	uint32_t delayFToI(float v);	// convert f.p. delay to fixed-point
};


// Delay partial template definitions
#define DELAY_DEF(l, Ti)\
template <class Tv=gam::real, class Ts=Synced>\
struct Delay##l : public Delay<Tv,Ti,Ts>{\
	typedef Delay<Tv,Ti,Ts> Base;\
	Delay##l(): Base(){}\
	Delay##l(float delay): Base(delay){}\
	Delay##l(float maxDelay, float delay): Base(maxDelay, delay){}\
};

DELAY_DEF(R, ipl::Round)
DELAY_DEF(T, ipl::Trunc)
DELAY_DEF(C, ipl::Cubic)
DELAY_DEF(A, ipl::AllPass)

#undef DELAY_DEF




/// Variable delay-line with multiple read taps.
template <class Tv=gam::real, template <class> class Si=ipl::Linear, class Ts=Synced>
class Delays : public Delay<Tv,Si,Ts> {
public:

	/// @param[in]	delay		Initial delay length. The size of the delay line will be the smallest possible power of two.
	/// @param[in]	numTaps		Number of reader taps.
	Delays(float delay, uint32_t numTaps);
	virtual ~Delays();

	void delay(float length, uint32_t tap);		///< Set delay length.
	//void delayNorm(float norm, uint32_t tap);		///< Set normalized delay as (0, 1) of buffer size.

	Tv read(uint32_t tap) const;		///< Read sample from tap.

	void taps(uint32_t numTaps);		///< Set number of read taps.
	uint32_t taps() const;				///< Get number of read taps.		
	
protected:
	uint32_t mTaps;		// number of delay read taps
	uint32_t * mDelays;	// read taps as delays from write tap
};



/// Variable length delay-line with feedback and/or feedforward.

/// The general comb filter transfer function provides N evenly spaced poles
/// and/or zeroes around the unit circle.  Feedback and feedforward produce
/// evenly spaced resonances and notches, respectively, in the frequency
/// response. Positive feeds result in even harmonics and negative feeds give
/// odd harmonics. If the feedback and feedforward amounts are inverses of each
/// other, an Nth order all-pass filter results. Comb filters are stable as
/// long as |feedback| < 1.
template <class Tv=gam::real, template <class> class Si=ipl::Linear, class Tf=gam::real, class Ts=Synced>
class Comb : public Delay<Tv,Si,Ts> {

private:
	typedef Delay<Tv,Si,Ts> Base;

public:
	using Base::operator();

	Comb();

	/// @param[in]	delay		Initial delay length. The size of the delay line will be the smallest possible power of two.
	/// @param[in]	ffd			Initial feedforward amount.
	/// @param[in]	fbk			Initial feedback amount.
	Comb(float delay, const Tf& ffd, const Tf& fbk);
	
	/// @param[in]	maxDelay	Maximum delay length. The size of the delay line will be the smallest possible power of two.
	/// @param[in]	delay		Initial delay length.
	/// @param[in]	ffd			Initial feedforward amount.
	/// @param[in]	fbk			Initial feedback amount.
	Comb(float maxDelay, float delay, const Tf& ffd, const Tf& fbk);
	
	virtual ~Comb(){}
	
	/// Set number of units for response to decay to end value.
	
	/// The sign of the decay length determines the sign of the feedback coefficient.
	/// The default end value of 0.001 (-60 dB) is the reverberation time of 
	/// the filter.  Setting the decay amount effects only the feedback value.
	/// The decay must be updated whenever the delay length of the filter changes.
	void decay(float units, float end = 0.001f);
	void fbk(const Tf& v);					///< Set feedback amount (-1, 1).
	void fbkAllPass(const Tf& v);			///< Set feedback amount (-1, 1) with feedforward set to opposite.
	void ffd(const Tf& v);					///< Set feedforward amount [-1, 1].
	void feeds(const Tf& fwd, const Tf& bwd){ ffd(fwd); fbk(bwd); }

	void set(float delay, const Tf& ffd, const Tf& fbk); ///< Set several parameters.

	Tv operator()(const Tv& i0);					///< Returns filtered value.
	Tv operator()(const Tv& i0, const Tv& oN);		///< Circulate filter with ffd & fbk
	Tv circulateFbk(const Tv& i0, const Tv& oN);	///< Circulate filter with fbk only	

	/// Filters sample (feedback only).
	Tv nextFbk(const Tv& i0);
	
	float norm() const;				///< Returns unity gain scale factor.
	float normFbk() const;			///< Returns unity gain scale factor due to feedback.
	float normFfd() const;			///< Returns unity gain scale factor due to feedforward.
	Tf ffd() const;				///< Returns feedforward amount.
	Tf fbk() const;				///< Returns feedback amount.

protected:
	Tf mFFD, mFBK;
};



// Base class for 2-pole or 2-zero filter
template <class Tv=gam::real, class Tp=gam::real, class Ts=Synced>
class Filter2 : public Ts{
public:
	Filter2(): d2(0), d1(0){}

	void freq(Tp v){ freqRef(v); }

	void width(Tp v){
		mWidth = v;
		mRad = scl::poleRadius(mWidth, Ts::ups());
		cd2 = -mRad * mRad;
		computeCoef();
	}

	/// Zero delay elements
	void zero(){ d2=d1=(Tv)0; }
	
	virtual void onResync(double r){ freq(mFreq); width(mWidth); }

protected:

	Filter2(Tp frq, Tp wid): mFreq(frq), mWidth(wid){ zero(); }

	void freqRef(Tp& v){
		mFreq = v;		
		v = scl::clip<Tp>(v * Ts::ups(), 0.5);
		mCos = scl::cosP3<Tp>(v);
		//mCos = scl::cosT8<T>(v * M_2PI);
		computeCoef();
	}
	
	void computeCoef(){	cd1 = (Tp)2 * mRad * mCos; }
	void delay(Tv v){ d2 = d1; d1 = v; }
	
	Tp mFreq, mWidth;
	Tp mB0, cd1, cd2;
	Tp mCos, mRad;
	Tv d2, d1;
};



/// Second-order all-pass filter.
template <class Tv=gam::real, class Tp=gam::real, class Ts=Synced>
class AllPass2 : public Filter2<Tv,Tp,Ts>{
public:

	/// @param[in] frq	Center frequency
	/// @param[in] wid	Bandwidth
	AllPass2(Tp frq, Tp wid=100): Base(frq, wid){ Ts::initSynced(); }

	/// Filter sample
	Tv operator()(Tv i0){
		i0 += d1*cd1 + d2*cd2;
		Tv o0 = i0*-cd2 - d1*cd1 + d2;
		delay(i0); return o0;
	}

protected:
	#define USE(var) using Base::var;
	typedef Filter2<Tv,Tp,Ts> Base; USE(d2) USE(d1) USE(mB0) USE(cd1) USE(cd2)
	#undef USE
};




/// Two-zero notch
template <class Tv=gam::real, class Tp=gam::real, class Ts=Synced>
class Notch : public Filter2<Tv,Tp,Ts>{
public:
	
	/// @param[in] frq	Center frequency
	/// @param[in] wid	Bandwidth
	Notch(Tp frq=1000, Tp wid=100): Base(frq, wid){ Ts::initSynced(); }

	/// Set center frequency
	void freq(Tp v){ Base::freq(v); computeGain(); }

	/// Set bandwidth
	void width(Tp v){ Base::width(v); computeGain(); }

	/// Filter sample
	Tv operator()(Tv i0){
		i0 *= mB0;
		Tv o0 = i0 - d1*cd1 - d2*cd2; delay(i0); return o0;
	}

	virtual void onResync(double r){ freq(mFreq); width(mWidth); }

protected:

	#define USE(var) using Base::var;
	typedef Filter2<Tv,Tp> Base;
	USE(mFreq) USE(mWidth) USE(d2) USE(d1) USE(mB0) USE(cd1) USE(cd2) USE(mRad)
	#undef USE

	// compute constant gain factor
	void computeGain(){ mB0 = Tp(1) / (Tp(1) + scl::abs(cd1) - cd2); }
};




/// Two-pole resonator
template <class Tv=gam::real, class Tp=gam::real, class Ts=Synced>
class Reson : public Filter2<Tv,Tp,Ts>{
public:

	/// @param[in] frq	Center frequency
	/// @param[in] wid	Bandwidth	
	Reson(Tp frq=440, Tp wid=100): Base(frq, wid){ Ts::initSynced(); }

	/// Set center frequency
	void freq(Tp v){
		Base::freqRef(v);
		mSin = scl::cosP3<Tp>(scl::foldOnce<Tp>(v - Tp(0.25), Tp(0.5)));
		computeGain();
	}

	/// Set bandwidth
	void width(Tp v){ Base::width(v); computeGain(); }

	void set(Tp frq, Tp wid){ Base::width(wid); freq(frq); }

	/// Filter sample
	Tv operator()(Tv i0){
		i0 *= mB0;
		i0 += d1*cd1 + d2*cd2; delay(i0); return i0; 
	}

	void onResync(double r){ freq(mFreq); width(mWidth); }

protected:
	#define USE(var) using Base::var;
	typedef Filter2<Tv,Tp,Ts> Base;
	USE(mFreq) USE(mWidth) USE(d2) USE(d1) USE(mB0) USE(cd1) USE(cd2) USE(mRad)
	#undef USE

	Tp mSin;
	
	// compute constant gain factor
	void computeGain(){ mB0 = (Tp(1) - mRad*mRad) * mSin; }
};



/// Hilbert transformer, converts real signal into complex
template <class Tv=gam::real, class Tp=gam::real, class Ts=Synced>
class Hilbert : public Ts{
public:
	#define SR 44100.
	Hilbert()
	:	cf0(   5.4135/SR), cf1(  41.118 /SR), cf2(  167.3595/SR),	/* for SR=44100 */
		cf3( 671.3715/SR), cf4(2694.363 /SR), cf5(11976.867 /SR),
		sf0(  18.786 /SR), sf1(  83.5065/SR), sf2(  335.1345/SR), 
		sf3(1344.4065/SR), sf4(5471.871 /SR), sf5(41551.671 /SR)
	{
		//Ts::initSynced();
	}
	#undef SR

	/// Convert input from real to complex
	Complex<Tv> operator()(const Tv& i){
		return Complex<Tv>(
			cf0(cf1(cf2(cf3(cf4(cf5(i)))))),
			-sf0(sf1(sf2(sf3(sf4(sf5(i))))))
		);
	}

//	virtual void onResync(double r){
//		//#define DO(o) o.freq(o.freq()*r); printf("%g\n", o.freq());
//		#define DO(o) printf("%g\n", o.freq());
//		DO(cf0) DO(cf1) DO(cf2) DO(cf3) DO(cf4) DO(cf5)
//		DO(sf0) DO(sf1) DO(sf2) DO(sf3) DO(sf4) DO(sf5)
//		#undef DO
//	}

protected:
	AllPass1<Tv, Tp, Synced1> cf0,cf1,cf2,cf3,cf4,cf5, sf0,sf1,sf2,sf3,sf4,sf5;
};



template <class Tv=gam::real>
class MovingAvg : public DelayN<Tv>{
public:

	MovingAvg(uint32_t size): Base(size), mSum(0), mRSize(0){ onResize(); }
	
	MovingAvg& operator=(const Tv& v){ DelayN<Tv>::operator=(v); return *this; }
	
	Tv operator()(const Tv& i0){
		return (mSum += i0 - Base::operator()(i0)) * mRSize;
	}
	
	void zero(){ assign(Tv(0)); }
	
	virtual void onResize(){ mRSize = 1./Base::size(); }

protected:
	typedef DelayN<Tv> Base;
	
	Tv mSum;
	double mRSize;
};




/// One-pole smoothing filter.
template<class Tv=gam::real, class Tp=gam::real, class Ts=Synced>
class OnePole : public Ts{ 
public:
	OnePole();

	/// @param[in]	freq	Smoothing frequency
	/// @param[in]	stored	Initial stored value
	OnePole(Tp freq, const Tv& stored=0);

	void operator  = (const Tv& val);	///< Stores input value for operator().
	void operator *= (const Tv& val);	///< Multiplies stored value by value.
	void freq(Tp val);					///< Set -3 dB bandwidth of pole.
	void smooth(Tp val);				///< Set smoothing coefficient directly.
	void zero(){ o1=0; }				///< Zero internal delay

	const Tv& operator()();				///< Returns filtered output using stored value (see store()).
	const Tv& operator()(const Tv& input);		///< Returns filtered output from input value.
	const Tv& last() const;				///< Returns last output.
	const Tv& stored() const;			///< Returns stored value.
	Tv& stored();						///< Returns stored value.
	bool zeroing(Tv eps=0.0001) const;	///< Returns whether the filter is outputting zeros.
	
	virtual void onResync(double r);

protected:
	Tp mFreq;
	Tp mA0, mB1;
	Tv mStored;
	Tv o1;
};




// Implementation_______________________________________________________________

#define TEM template <class Tv>
#define T_VPS template <class Tv, class Tp, class Ts>

//---- AllPass1

T_VPS AllPass1<Tv,Tp,Ts>::AllPass1(Tp frq)
:	d1(Tv(0)), mFreq(frq)
{	Ts::initSynced(); }

T_VPS inline void AllPass1<Tv,Tp,Ts>::freq(Tp v){
	mFreq = v;
	c = tan( scl::freqToRad(v, Ts::ups()) - M_PI_4); // valid ang in [-pi/4, pi/4]
}

T_VPS inline void AllPass1<Tv,Tp,Ts>::freqF(Tp v){
	mFreq = v;
	c = scl::freqToRad(v, Ts::ups());
	c = Tp(1.27323954474) * c - Tp(1);		// 4/pi * c - 1, linear apx
	c = c * (Tp(0.76) + Tp(0.24) * c * c);
}

T_VPS inline Tv AllPass1<Tv,Tp,Ts>::operator()(Tv i0){ 
	i0 -= d1 * c;			// o0 = c * (i0 - c * d1) + d1
	Tv o0 = i0 * c + d1;
	d1 = i0;
	return o0;
}

T_VPS inline Tv AllPass1<Tv,Tp,Ts>::high(Tv i0){ return (i0 - operator()(i0)) * Tv(0.5); }
T_VPS inline Tv AllPass1<Tv,Tp,Ts>::low (Tv i0){ return (i0 + operator()(i0)) * Tv(0.5); }

T_VPS inline Tp AllPass1<Tv,Tp,Ts>::freq(){ return mFreq; }
T_VPS void AllPass1<Tv,Tp,Ts>::onResync(double r){ freq(mFreq); }




//---- Biquad
T_VPS Biquad<Tv,Tp,Ts>::Biquad(Tp frq, Tp res, int type)
:	d1(0), d2(0)
{
	Ts::initSynced();
	set(frq, res, type);
}

T_VPS void Biquad<Tv,Tp,Ts>::onResync(double r){
	mFrqToRad = M_2PI * Ts::ups();
	freq(mFreq);
}

T_VPS void Biquad<Tv,Tp,Ts>::set(Tp freqA, Tp resA, int typeA){
	mRes = resA;
	mType = typeA;
	freq(freqA);
}

T_VPS void Biquad<Tv,Tp,Ts>::zero(){ d1=d2=(Tv)0; }

T_VPS void Biquad<Tv,Tp,Ts>::coef(Tp a0, Tp a1, Tp a2, Tp b0, Tp b1, Tp b2){
	mA0=a0; mA1=a1; mA2=a2; mB0=b0; mB1=b1; mB2=b2;
}

T_VPS inline void Biquad<Tv,Tp,Ts>::freq(Tp v){
//	float w = freqA * mFrqToRad;	// radial frequency [0, pi)
//	mReal = cos(w);
//	mImag = sin(w);

	mFreq = v;
	float phs = scl::clip(mFreq * mFrqToRad, 3.11f);
	mReal = scl::cosT8(phs);
	mImag = scl::sinT7(phs);
	res(mRes);
}

T_VPS inline void Biquad<Tv,Tp,Ts>::res(Tp v){
	mRes = v;
	mAlpha = mImag / mRes;
	mB0 = Tp(1) / (Tp(1) + mAlpha);	// reciprocal of mB0
	mB1 = Tp(-2) * mReal * mB0;
	mB2 = (Tp(1) - mAlpha) * mB0;
	
	type(mType);
}

T_VPS inline void Biquad<Tv,Tp,Ts>::set(Tp frq, Tp res){ set(frq, res, mType); }

/*
void setAllpass(float fr, float R){
  R= 1.f/R;
  double cs = -2*cos((2*PI*fr)/sr); 
  mA1 = cs*R;
  mA2 = R*R;
  mB1 = cs/R;
  mB2 = 1/m_a2;
  mA0 = 1;
}
*/

T_VPS inline void Biquad<Tv,Tp,Ts>::type(int typeA){
	mType = typeA;
	
	switch(mType){
	case Filter::LP:
		mA1 = (Tp(1) - mReal) * mB0;
		mA0 = mA1 * Tp(0.5);
		mA2 = mA0;
		break;
	case Filter::HP:
		mA1 = -(Tp(1) + mReal) * mB0;
		mA0 = -mA1 * Tp(0.5);
		mA2 = mA0;
		break;
	case Filter::BP:
		mA0 = mImag * Tp(0.5) * mB0;
		mA1 = Tp(0);
		mA2 = -mA0;
		break;
	case Filter::BPC:
        mA0 = mAlpha * mB0;
        mA1 = Tp(0);
        mA2 = -mA0;
		break;
	case Filter::BR:
        mA0 = mB0;	// 1.f * a0
        mA1 = mB1;
        mA2 = mB0;	// 1.f * a0
		break;
	case Filter::AP:
		mA0 = mB2;
		mA1 = mB1;
		mA2 = mB0;
		break;
	default:;
	};
}

T_VPS inline Tv Biquad<Tv,Tp,Ts>::operator()(Tv i0){
	// Direct form II
	i0 = i0 - d1*mB1 - d2*mB2;
	Tv o0 = i0*mA0 + d1*mA1 + d2*mA2;
	d2 = d1; d1 = i0;
	return o0;
}

T_VPS inline Tv Biquad<Tv,Tp,Ts>::nextBP(Tv i0){
	i0 = i0 - d1*mB1 - d2*mB2;	
	Tv o0 = (i0 - d2)*mA0;
	d2 = d1; d1 = i0;
	return o0;
}



//---- Delay
#define TM1 template <class Tv, template <class> class Ti, class Ts>
#define TM2 Tv,Ti,Ts

#define DELAY_INIT mPhase(0)

TM1 Delay<TM2>::Delay()
:	ArrayPow2<Tv>(), DELAY_INIT
{	Ts::initSynced(); }

TM1 Delay<TM2>::Delay(float maxDelay, float delay)
:	ArrayPow2<Tv>(), DELAY_INIT
{
	Ts::initSynced();
	this->maxDelay(maxDelay);
	this->delay(delay);
}

TM1 Delay<TM2>::Delay(float delay)
:	ArrayPow2<Tv>(), DELAY_INIT
{	//printf("Delay::Delay(float)\n");
	Ts::initSynced();
	this->maxDelay(delay);
	this->delay(delay);
}

#undef DELAY_INIT

TM1 void Delay<TM2>::maxDelay(float length){ //printf("Delay::maxDelay(%f)\n", length);
	mMaxDelay = length;
	if(Ts::sync() && Ts::sync()->hasBeenSet()){
		//printf("Delay::maxDelay(): resize to %d\n", (uint32_t)(mMaxDelay * spu()));
		this->resize((uint32_t)(mMaxDelay * Ts::spu())); // calls onResize() -> onResync(double r)
	}
}

TM1 void Delay<TM2>::zero(){ 
	for(uint32_t i=0; i<ArrayPow2<Tv>::size(); ++i) ArrayPow2<Tv>::elems()[i] = 0;
}


TM1 inline Tv Delay<TM2>::operator()() const{ return mIpol(*this, mPhase - mDelay); }

TM1 inline Tv Delay<TM2>::operator()(const Tv& i0){
//	writePre(i0);
//	return (*this)();
	Tv o0 = (*this)();	// read delayed element
	write(i0);			// write input element
	return o0;
}

TM1 inline uint32_t Delay<TM2>::delayFToI(float v){
	return castIntRound((v * mDelayFactor) * 4294967296.);
	//return scl::unitToUInt(v * mDelayFactor);
}

TM1 inline void Delay<TM2>::incPhase(){ mPhase += mPhaseInc; }

TM1 void Delay<TM2>::onResize(){ //printf("Delay::onResize %d elements\n", this->size());
	mPhaseInc = this->oneIndex();
	for(uint32_t i=0; i<this->size(); ++i) (*this)[i] = Tv(0);
	onResync(1);
}

TM1 void Delay<TM2>::onResync(double r){ //printf("Delay::onSyncChange\n");
	maxDelay(mMaxDelay);
	refreshDelayFactor();
	delay(mDelayLength);
}

TM1 inline Tv Delay<TM2>::read(float ago){ return mIpol(*this, mPhase - delayFToI(ago)); }

TM1 void Delay<TM2>::refreshDelayFactor(){ mDelayFactor = 1. / maxDelay(); }

TM1 inline void Delay<TM2>::write(const Tv& v){
	mem::put(this->elems(), this->fracBits(), mPhase, v);
	incPhase();
}

TM1 inline void Delay<TM2>::writePre(const Tv& v){
	incPhase();
	mem::put(this->elems(), this->fracBits(), mPhase, v);
}

TM1 inline void Delay<TM2>::delay(float v){
	mDelayLength = v;
	mDelay = delayFToI(v);
}

TM1 inline float Delay<TM2>::delay() const { return mDelayLength; }

TM1 inline float Delay<TM2>::delayUnit() const {
	return uintToUnit<float>(mDelay);
}

TM1 inline void Delay<TM2>::delayUnit(float n){ mDelay = unitToUInt(n); }

TM1 inline uint32_t Delay<TM2>::delayIndex(uint32_t delay) const {
	return this->index(mPhase - (delay << this->fracBits()));
}

TM1 inline void Delay<TM2>::freq(float v){ delay(1.f/v); }

TM1 inline uint32_t Delay<TM2>::indexBack() const {
	return this->index(mPhase + this->oneIndex());
}

TM1 inline float Delay<TM2>::maxDelay() const { return this->size() * Ts::ups(); }

TM1 void Delay<TM2>::print(){
	printf("SPU:       %f\n", Ts::spu());
	printf("Buffer:    %p\n", this->elems());
	printf("BufBits:   %d\n", this->log2Size());
	printf("FracBits:  %d\n", this->fracBits());
	printf("Phase:     %d\n", mPhase);
	printf("PhaseInc:  %d\n", mPhaseInc);
	printf("Delay:     %d, %f\n", mDelay, mDelayLength);
	printf("Max Delay: %f\n",  mMaxDelay);
	printf("DlyFactor: %f\n",  mDelayFactor);
}



//---- Delays
TM1 Delays<TM2>::Delays(float delay, uint32_t numTaps)
:	Delay<TM2>(delay), mTaps(0), mDelays(0)
{	taps(numTaps); }

TM1 Delays<TM2>::~Delays(){
	mem::free(mDelays);
}

TM1 void Delays<TM2>::taps(uint32_t numTaps){
	if(mem::resize(mDelays, mTaps, numTaps)) mTaps = numTaps;
}

TM1 inline void Delays<TM2>::delay(float v, uint32_t tap){
	mDelays[tap] = this->delayFToI(v);
}

TM1 inline Tv Delays<TM2>::read(uint32_t tap) const {
	return mIpol(*this, this->mPhase - mDelays[tap]);
}

TM1 inline uint32_t Delays<TM2>::taps() const { return mTaps; }

#undef TM1
#undef TM2


//---- Comb
#define TM1 template<class Tv, template<class> class Si, class Tf, class Ts>
#define TM2 Tv,Si,Tf,Ts
TM1 Comb<TM2>::Comb(): Delay<Tv,Si,Ts>(), mFFD(0), mFBK(0){}

TM1 Comb<TM2>::Comb(float delay, const Tf& ffd, const Tf& fbk)
:	Delay<Tv,Si,Ts>(delay), mFFD(ffd), mFBK(fbk)
{}

TM1 Comb<TM2>::Comb(float delayMax, float delay, const Tf& ffd, const Tf& fbk)
:	Delay<Tv,Si,Ts>(delayMax, delay), mFFD(ffd), mFBK(fbk)
{}

TM1 inline Tv Comb<TM2>::operator()(const Tv& i0){ return (*this)(i0, (*this)()); }

TM1 inline Tv Comb<TM2>::operator()(const Tv& i0, const Tv& oN){
	Tv t = i0 + oN * mFBK;
	this->write(t);
	return oN + t * mFFD;
}

TM1 inline Tv Comb<TM2>::circulateFbk(const Tv& i0, const Tv& oN){
	this->write(i0 + oN * mFBK);
	return oN;
}

TM1 inline Tv Comb<TM2>::nextFbk(const Tv& i0){ return circulateFbk(i0, (*this)()); }

TM1 inline void Comb<TM2>::decay(float units, float end){
	mFBK = ::pow(end, this->delay() / scl::abs(units));
	if(units < 0.f) mFBK = -mFBK;
}

TM1 inline void Comb<TM2>::fbk(const Tf& v){ mFBK=v; }
TM1 inline void Comb<TM2>::fbkAllPass(const Tf& v){ fbk(v); ffd(-v); }
TM1 inline void Comb<TM2>::ffd(const Tf& v){ mFFD=v; }
TM1 inline void Comb<TM2>::set(float d, const Tf& ff, const Tf& fb){ this->delay(d); ffd(ff); fbk(fb); }

TM1 inline Tf Comb<TM2>::fbk() const { return mFBK; }
TM1 inline Tf Comb<TM2>::ffd() const { return mFFD; }
TM1 inline float Comb<TM2>::norm() const { return (1.f - scl::abs(fbk()))/(1.f + scl::abs(ffd())); }
TM1 inline float Comb<TM2>::normFbk() const { return 1.f - scl::abs(fbk()); }
TM1 inline float Comb<TM2>::normFfd() const { return 1.f/(1.f + scl::abs(ffd())); }

#undef TM1
#undef TM2



//---- OnePole
T_VPS OnePole<Tv,Tp,Ts>::OnePole()
:	mFreq(10), mStored(Tv(0)), o1(Tv(0))
{	Ts::initSynced(); }

T_VPS OnePole<Tv,Tp,Ts>::OnePole(Tp frq, const Tv& stored)
:	mFreq(frq), mStored(stored), o1(stored)
{	Ts::initSynced(); }

T_VPS void OnePole<Tv,Tp,Ts>::onResync(double r){ freq(mFreq); }

T_VPS inline void OnePole<Tv,Tp,Ts>::operator  = (const Tv& v){ mStored  = v; }
T_VPS inline void OnePole<Tv,Tp,Ts>::operator *= (const Tv& v){ mStored *= v; }

// f = ln(mB1) / -M_2PI * SR  ( @ 44.1 f = ln(c01) * -7018.733)
// @ 44.1 : 0.9     <=> 739.5
// @ 44.1 : 0.99    <=>  70.54
// @ 44.1 : 0.999   <=>   7.022
// @ 44.1 : 0.9999  <=>   0.702
// @ 44.1 : 0.99999 <=>   0.0702
T_VPS inline void OnePole<Tv,Tp,Ts>::freq(Tp v){
	mFreq = v;	
	v = scl::max(v, Tp(0));	// ensure positive freq
	
	// freq is half the bandwidth of a pole at 0
	smooth( scl::poleRadius(Tp(2) * v, Ts::ups()) );
	//printf("%f, %f, %f\n", spu(), mB1, v);
}

T_VPS inline void OnePole<Tv,Tp,Ts>::smooth(Tp v){ mB1=v; mA0=Tp(1) - scl::abs(v); }

T_VPS inline const Tv& OnePole<Tv,Tp,Ts>::operator()(){ return (*this)(mStored); }
T_VPS inline const Tv& OnePole<Tv,Tp,Ts>::operator()(const Tv& i0){ o1 = o1*mB1 + i0*mA0; return o1; }

T_VPS inline const Tv& OnePole<Tv,Tp,Ts>::last() const { return o1; }
T_VPS inline const Tv& OnePole<Tv,Tp,Ts>::stored() const { return mStored; }
T_VPS inline Tv& OnePole<Tv,Tp,Ts>::stored(){ return mStored; }
T_VPS inline bool OnePole<Tv,Tp,Ts>::zeroing(Tv eps) const { return scl::abs(o1) < eps && mStored == Tv(0); }


#undef TEM
#undef T_VPS

} // gam::
#endif
