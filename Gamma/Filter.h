#ifndef GAMMA_FILTER_H_INC
#define GAMMA_FILTER_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include "Gamma/ipl.h"
#include "Gamma/scl.h"
#include "Gamma/Containers.h"
#include "Gamma/Domain.h"
#include "Gamma/Types.h"

namespace gam{

/// Signal transformers

/// \defgroup Filter
/// A linear filter boosts or attenuates the amplitudes of frequency components.
/// In contrast, a non-linear filter typically produces new frequency components.
/// Both types of filters also typically alter the phases of frequencies.
/// Most of the filters in this group are linear filters.


/// Filter types
enum FilterType{
	LOW_PASS,			/**< Low-pass */
	HIGH_PASS,			/**< High-pass */
	BAND_PASS,			/**< Band-pass */
	RESONANT,			/**< Resonant band-pass */
	BAND_REJECT,		/**< Band-reject */
	ALL_PASS,			/**< All-pass */
	PEAKING,			/**< Peaking */
	LOW_SHELF,			/**< Low-shelf */
	HIGH_SHELF,			/**< High-shelf */
	SMOOTHING,			/**< Smoothing */
	THRU_PASS			/**< Thru-pass (no filter) */
};



/// Returns pole radius given a unit bandwidth
template <class T>
inline T poleRadius(T bw){ return exp(-M_PI * bw); }

/// Returns pole radius given a bandwidth and sampling interval
template <class T>
inline T poleRadius(T bw, double ups){ return poleRadius(bw * ups); }
//return (T)1 - (M_2PI * bw * ups); // linear apx for fn < ~0.02

/// Convert domain frequency to radians clipped to interval [0, pi).
template <class T>
inline T freqToRad(T freq, double ups){ return scl::clip(freq * ups, 0.499) * M_PI; }

namespace{
	#ifndef GAM_FILTER_PRECISION
	#define GAM_FILTER_PRECISION 1
	#endif

	#if GAM_FILTER_PRECISION >= 2
	inline double getFreqToAng(double ups){ return M_2PI*ups; }
	template <class T> inline T getRe(T ang){ return std::cos(ang); }
	template <class T> inline T getIm(T ang){ return std::sin(ang); }
	template <class T> inline T clipAng(T ang){ return scl::clip<T>(ang, 3.14); }
	#elif GAM_FILTER_PRECISION == 1
	// sinFast domain [-2,2] -> [-pi,pi]
	inline double getFreqToAng(double ups){ return 4.*ups; }
	template <class T> inline T getRe(T ang){ return scl::sinFast(T(1)-ang); }
	template <class T> inline T getIm(T ang){ return scl::sinFast02(ang); }
	template <class T> inline T clipAng(T ang){ return scl::clip<T>(ang, 2.); }
	#endif

	template <class T>
	inline void getPole(T& re, T& im, T ang){
		ang = clipAng(ang);
		re = getRe(ang);
		im = getIm(ang);
	}

	template <class T>
	inline T getPoleRe(T ang){
		return getRe(clipAng(ang));
	}
};


/// First-order all-pass filter

/// This filter has the transfer function H(z) = (a + z^-1) / (1 + a z^-1).
/// It is essentially a "phase" filter in that its only effect is to 
/// shift frequencies from 0 to Nyquist from 0 to -180 degress, respectively.
/// The cutoff frequency is the frequency where the phase is shifted by -90 
/// degrees.  A lowpass or highpass filter can be constructed by adding or 
/// subtracting, respectively, the output of the filter to or from the input.
/// When the center frequency is fs/4, the filter acts as a unit delay.
/// When the center frequency is 0, the filter acts as a no-op.
/// When the center frequency is fs/2, the filter acts as an inverter.
/// \tparam Tv	Value (sample) type
/// \tparam Tp	Parameter type
/// \tparam Td	Domain type
/// \ingroup Filter
template<class Tv=gam::real, class Tp=gam::real, class Td=GAM_DEFAULT_DOMAIN>
class AllPass1 : public Td {
public:
	///
	/// \param[in]	frq		Center frequency
	AllPass1(Tp frq = Tp(1000));

	void freq (Tp v);		///< Set cutoff frequency
	void freqF(Tp v);		///< Faster, but slightly less accurate than freq()	
	void zero(){ d1=Tv(0); }
	
	Tv operator()(Tv in);	///< Filters sample
	
	Tv high(Tv in);			///< High-pass filters sample
	Tv low (Tv in);			///< Low-pass filters sample
	
	Tp freq();				///< Get current cutoff frequency
	
	void onDomainChange(double r);
	
protected:
	Tv d1;		// once delayed value
	Tp c;		// feed coefficient
	Tp mFreq;
};


/// 2-pole/2-zero IIR filter

/// The biquadratic filter contains 2 zeros and 2 poles in its transfer
/// function. The zeros provide a better response near the DC and Nyquist
/// frequencies than an all-pole filter would. Second-order IIRs have a 12 
/// dB/octave cutoff slope and are normally cascaded (run in series) to obtain
/// a sharper response. This particular implementation utilizes the RBJ 
/// Sallen-Key/Butterworth design from:
///
/// http://www.musicdsp.org/files/Audio-EQ-Cookbook.txt
///
/// \tparam Tv	Value (sample) type
/// \tparam Tp	Parameter type
/// \tparam Td	Domain type
/// \ingroup Filter
template <class Tv=gam::real, class Tp=gam::real, class Td=GAM_DEFAULT_DOMAIN>
class Biquad : public Td{
public:

	/// \param[in]	frq		Center frequency
	/// \param[in]	res		Resonance (Q) amount in (0, inf)
	/// \param[in]	type	Type of filter
	/// \param[in]	lvl		Amplitude level for PEAKING, LOW_SHELF, HIGH_SHELF
	Biquad(
		Tp frq = Tp(1000),
		Tp res = Tp(0.707),
		FilterType type = LOW_PASS,
		Tp lvl = Tp(1)
	);


	/// Get array of 3 feedforward coefficients
	const Tp * a() const { return mA; }
	Tp * a(){ return mA; }
	
	/// Get array of 3 feedback coefficients (first element, b0, is not used)
	const Tp * b() const { return mB; }
	Tp * b(){ return mB; }

	/// Set feedforward (a) and feedback (b) coefficients; b0 assumed to be 1
	void coef(Tp a0, Tp a1, Tp a2, Tp b1, Tp b2);


	void freq(Tp v);					///< Set center frequency
	void res(Tp v);						///< Set resonance (Q); gain at center frequency
	void level(Tp v);					///< Set level (PEAKING, LOW_SHELF, HIGH_SHELF types only)
	void set(Tp frq, Tp res);			///< Set filter center frequency and resonance
	void set(Tp frq, Tp res, FilterType type);	///< Set all filter params
	void type(FilterType type);			///< Set type of filter
	void zero();						///< Zero internal delays

	Tv operator()(Tv in);				///< Filter next sample
	Tv nextBP(Tv in);					///< Optimized for band-pass types
	
	Tp freq() const;					///< Get center frequency
	Tp res() const;						///< Get resonance (Q)
	Tp level() const;					///< Get level
	FilterType type() const;			///< Get filter type
	
	void onDomainChange(double r);

protected:
	Tp mA[3];			// feedforward coefs
	Tp mB[3];			// feedback coefs (first element used to scale coefs)
	Tv d1, d2;			// inner sample delays
	Tp mFreq, mResRecip;// center frequency, 1/resonance
	Tp mLevel;			// amplitude level (for peaking)
	FilterType mType;
	Tp mReal, mImag;	// real, imag components of center frequency
	Tp mAlpha, mBeta;
	Tp mFreqToAng;

	void levelNoResUpdate(Tp v);
	void resRecip(Tp v);
};



/// DC frequency blocker

/// \tparam Tv	Value (sample) type
/// \tparam Tp	Parameter type
/// \tparam Td	Domain type
/// \ingroup Filter
template <class Tv=gam::real, class Tp=gam::real, class Td=GAM_DEFAULT_DOMAIN>
class BlockDC : public Td{
public:

	/// \param[in] bwidth	Bandwidth of pole
	BlockDC(Tp bwidth = Tp(35))
	:	d1(0), mWidth(bwidth)
	{
		onDomainChange(1);
	}

	/// Filter sample
	Tv operator()(Tv in){
		Tv t = in + d1*mB1;
		Tv o0 = t - d1;
		d1 = t;
		return o0;
	}

	/// Set bandwidth of pole
	void width(Tp v){
		mWidth = v;
		mB1 = poleRadius(v, Td::ups());
	}

	void zero(){ d1=0; }

	void onDomainChange(double /*r*/){ width(mWidth); }

protected:
	Tv d1;
	Tp mWidth, mB1;
};



/// Nyquist frequency blocker

/// \tparam Tv	Value (sample) type
/// \tparam Tp	Parameter type
/// \tparam Td	Domain type
/// \ingroup Filter  
template <class Tv=gam::real, class Tp=gam::real, class Td=GAM_DEFAULT_DOMAIN>
class BlockNyq : public BlockDC<Tv,Tp,Td>{
public:

	/// \param[in] bwidth	Bandwidth of pole
	BlockNyq(Tp bwidth = Tp(35)){
		width(bwidth);
	}

	/// Filter sample
	Tv operator()(Tv in){
		Tv t = in + d1*mB1;
		Tv o0 = t + d1;
		d1 = t;
		return o0;
	}

	/// Set bandwidth of pole
	void width(Tp v){
		Base::width(v);
		mB1 = -mB1;
	}

protected:
	typedef BlockDC<Tv,Tp,Td> Base;
	using Base::d1; using Base::mB1;
};



/// Abstract base class for 2-pole or 2-zero filter

/// \tparam Tv	Value (sample) type
/// \tparam Tp	Parameter type
/// \tparam Td	Domain type
/// \ingroup Filter
template <class Tv=gam::real, class Tp=gam::real, class Td=GAM_DEFAULT_DOMAIN>
class Filter2 : public Td{
public:

	/// Set frequency
	void freq(Tp v){
		mFreq = v;		
		/*v = scl::clip<Tp>(v * Td::ups(), 0.5);
		//mReal = scl::cosP3<Tp>(v);
		mReal = scl::cosT8<Tp>(v * M_2PI);*/
		mReal = getPoleRe(v * mFreqToAng);
		computeCoef1();
	}

	/// Set bandwidth
	void width(Tp v){
		mWidth = v;
		mRad = poleRadius(mWidth, Td::ups());
		mC[2] = -mRad * mRad;
		computeCoef1();
	}

	/// Zero delay elements
	void zero(){ d2=d1=Tv(0); }
	
	void onDomainChange(double r){
		mFreqToAng = getFreqToAng(Td::ups());
		onCacheVars();
	}

protected:

	Filter2(Tp frq, Tp wid)
	:	mFreq(frq), mWidth(wid)
	{
		zero();
		onDomainChange(1);
	}

	virtual void onCacheVars(){ freq(mFreq); width(mWidth); }

	void computeCoef1(){ mC[1] = Tp(2) * mRad * mReal; }
	void delay(Tv v){ d2=d1; d1=v; }
	Tp& gain(){ return mC[0]; }
	
	Tp mFreq, mWidth;
	Tp mFreqToAng;
	Tp mC[3];			// coefficients
	Tp mReal, mRad;
	Tv d2, d1;			// 2- and 1-sample delays
};


#define INHERIT_FILTER2 \
typedef Filter2<Tv,Tp,Td> Base;\
using Base::onDomainChange;\
using Base::computeCoef1;\
using Base::mFreq;\
using Base::mFreqToAng;\
using Base::mWidth;\
using Base::gain;\
using Base::mC;\
using Base::mRad;\
using Base::mReal;\
using Base::d2;\
using Base::d1;


/// Second-order all-pass filter

/// This all-pass filter shifts phases from 0 to 2 pi from 0 to Nyquist. The
/// center frequency controls where the phase is shifted by pi.
/// \tparam Tv	Value (sample) type
/// \tparam Tp	Parameter type
/// \tparam Td	Domain type
/// \ingroup Filter
template <class Tv=gam::real, class Tp=gam::real, class Td=GAM_DEFAULT_DOMAIN>
class AllPass2 : public Filter2<Tv,Tp,Td>{
public:

	/// \param[in] frq	Center frequency
	/// \param[in] wid	Bandwidth
	AllPass2(Tp frq = Tp(1000), Tp wid = Tp(100))
	:	Base(frq, wid){}

	/// Filter sample
	Tv operator()(Tv in){
		Tv  t = in + d1*mC[1] + d2*mC[2];
		Tv o0 = d2 - d1*mC[1] -  t*mC[2];
		this->delay(t);
		return o0;
	}

protected:
	INHERIT_FILTER2;
};



/// Two-zero notch

/// \tparam Tv	Value (sample) type
/// \tparam Tp	Parameter type
/// \tparam Td	Domain type
/// \ingroup Filter
template <class Tv=gam::real, class Tp=gam::real, class Td=GAM_DEFAULT_DOMAIN>
class Notch : public Filter2<Tv,Tp,Td>{
public:
	
	/// \param[in] frq	Center frequency
	/// \param[in] wid	Bandwidth
	Notch(Tp frq = Tp(1000), Tp wid = Tp(100))
	:	Base(frq, wid){}

	/// Set center frequency
	void freq(Tp v){ Base::freq(v); computeGain(); }

	/// Set bandwidth
	void width(Tp v){ Base::width(v); computeGain(); }

	/// Filter sample
	Tv operator()(Tv in){
		Tv t = in * gain();
		Tv o0 = t - d1*mC[1] - d2*mC[2];
		this->delay(t);
		return o0;
	}

protected:
	INHERIT_FILTER2;

	virtual void onCacheVars(){ freq(mFreq); width(mWidth); }

	// compute constant gain factor
	void computeGain(){ gain() = Tp(1) / (Tp(1) + scl::abs(mC[1]) - mC[2]); }
};



/// Two-pole resonator

/// \tparam Tv	Value (sample) type
/// \tparam Tp	Parameter type
/// \tparam Td	Domain type
/// \ingroup Filter
template <class Tv=gam::real, class Tp=gam::real, class Td=GAM_DEFAULT_DOMAIN>
class Reson : public Filter2<Tv,Tp,Td>{
public:

	/// \param[in] frq	Center frequency
	/// \param[in] wid	Bandwidth	
	Reson(Tp frq = Tp(1000), Tp wid = Tp(100))
	:	Base(frq, wid){}

	/// Set center frequency
	void freq(Tp v){
		getPole(mReal, mImag, v * mFreqToAng);
		computeCoef1();
		computeGain();
	}

	/// Set bandwidth
	void width(Tp v){ Base::width(v); computeGain(); }

	void set(Tp frq, Tp wid){ Base::width(wid); freq(frq); }

	/// Filter sample
	Tv operator()(Tv in){
		Tv t = in * gain() + d1*mC[1] + d2*mC[2];
		this->delay(t);
		return t;
	}

protected:
	INHERIT_FILTER2;
	Tp mImag;

	virtual void onCacheVars(){ freq(mFreq); width(mWidth); }

	// compute constant gain factor
	void computeGain(){ gain() = (Tp(1) - mRad*mRad) * mImag; }
};



/// Hilbert transform filter

/// This filter converts a real signal into a complex (analytic) signal.
/// Corresponding harmonics of the real and imaginary parts of the output are
/// 90 degrees out of phase.
/// The analytic signal can be used for frequency shifting or determining the
/// amplitude envelope of the input signal.
/// Technically speaking, a Hilbert transform converts a real signal into its
/// harmonic conjugate. The input and output of the Hilbert transform, comprise
/// the real and imaginary components of a complex (analytic) signal.
///
/// \tparam Tv	Value (sample) type
/// \tparam Tp	Parameter type
/// \ingroup Filter
template <class Tv=gam::real, class Tp=gam::real>
class Hilbert {
public:
	#define SR 44100.
	Hilbert()
	:	cf0(   5.4135/SR), cf1(  41.118 /SR), cf2(  167.3595/SR),	/* for SR=44100 */
		cf3( 671.3715/SR), cf4(2694.363 /SR), cf5(11976.867 /SR),
		sf0(  18.786 /SR), sf1(  83.5065/SR), sf2(  335.1345/SR), 
		sf3(1344.4065/SR), sf4(5471.871 /SR), sf5(41551.671 /SR)
	{}
	#undef SR

	/// Convert input from real to complex
	Complex<Tv> operator()(Tv in){
		return Complex<Tv>(
			 cf0(cf1(cf2(cf3(cf4(cf5(in)))))),
			-sf0(sf1(sf2(sf3(sf4(sf5(in))))))
		);
	}
	
	void zero()
	{
		cf0.zero(); cf1.zero(); cf2.zero(); cf3.zero(); cf4.zero(); cf5.zero();
		sf0.zero(); sf1.zero(); sf2.zero(); sf3.zero(); sf4.zero(); sf5.zero();
	}

protected:
	AllPass1<Tv, Tp, Domain1> cf0,cf1,cf2,cf3,cf4,cf5, sf0,sf1,sf2,sf3,sf4,sf5;
};



/// Leaky integrator

/// \tparam Tv	Value (sample) type
/// \tparam Tp	Parameter type
/// \ingroup Filter
template <class Tv=double, class Tp=double>
class Integrator{
public:

	/// \param[in] leakCoef		Leak coefficient, in [0,1)
	/// \param[in] initVal		Initial value
	Integrator(Tp leakCoef = Tp(1), Tv initVal = Tv(0)){
		mo[0] = initVal;
		leak(leakCoef);
	}

	/// Filter input value
	Tv operator()(Tv in) const {
		return mo[0] = mo[0]*mb[0] + in;
	}
	
	Integrator& leak(Tp v){ mb[0]=v; return *this; }
	Integrator& zero(){ mo[0]=Tv(0); return *this; }

protected:
	mutable Tv mo[1];
	Tp mb[1];
};



/// Differencing filter

/// Returns difference between successive samples.
/// The frequency response has a zero at DC.
/// \ingroup Filter
template <class T=gam::real>
class Differencer{
public:
	T operator()(T in){
		T res = in-mPrev;
		mPrev = in;
		return res;
	}
private:
	T mPrev = T(0);
};



/// Moving average filter

/// A moving average (MA) filter is a low-pass FIR filter. The filter kernel is
/// a rectangle function with magnitude equal to the reciprocal of the kernel
/// size, N. The cutoff frequency is controlled by the kernel size and is 
/// approximately fs/N where fs is the sampling frequency. Due to the symmetry 
/// of the window, the moving average filter can be implemented efficiently
/// using a single delay line with O(1) processing time complexity.
/// \ingroup Filter
template <class Tv=gam::real>
class MovingAvg : public DelayN<Tv>{
public:

	MovingAvg()
	:	mSum(0), mRSize(0)
	{}

	/// \param[in] size		Kernel size, greater than 1
	explicit MovingAvg(unsigned size)
	:	Base(size), mSum(0), mRSize(0)
	{	onResize(); }

	MovingAvg& operator=(const Tv& v){ DelayN<Tv>::operator=(v); return *this; }
	
	Tv operator()(Tv in){
		mSum += in - Base::operator()(in);
		return mSum * mRSize;
	}

	virtual void onResize(){
		mRSize = 1./Base::size();
		//printf("mRSize = %g\n", mRSize);
	}

protected:
	typedef DelayN<Tv> Base;
	
	Tv mSum;		// the moving sum
	double mRSize;	// reciprocal of the size
};



/// One-pole filter

/// This filter uses a single pole at either DC or Nyquist to create a low-pass
/// or high-pass response, respectively.
///
/// \tparam Tv	Value (sample) type
/// \tparam Tp	Parameter type
/// \tparam Td	Domain type
/// \ingroup Filter
template<class Tv=gam::real, class Tp=gam::real, class Td=GAM_DEFAULT_DOMAIN>
class OnePole : public Td{ 
public:

	/// \param[in]	freq	Cutoff frequency
	/// \param[in]	type	Filter type (LOW_PASS or HIGH_PASS)
	/// \param[in]	stored	Initial stored value
	OnePole(Tp freq = Tp(1000), FilterType type = LOW_PASS, const Tv& stored = Tv(0));

	/// \param[in]	freq	Cutoff frequency
	/// \param[in]	stored	Initial stored value
	OnePole(Tp freq, const Tv& stored);


	const Tp& freq() const { return mFreq; }	///< Get cutoff frequency

	void type(FilterType type);			///< Set type of filter (gam::LOW_PASS or gam::HIGH_PASS)
	void freq(Tp val);					///< Set cutoff frequency (-3 dB bandwidth of pole)

	/// Set lag length of low-pass response

	/// \param[in] length	Length of lag with zero for no lag
	/// \param[in] thresh	Value to which a downward unit step reaches after
	///						the lag length. Must be greater than 0.
	void lag(Tp length, Tp thresh=Tp(0.001));

	void smooth(Tp val);				///< Set smoothing coefficient directly
	void zero(){ mO1=0; }				///< Zero internal delay
	void reset(Tv v = Tv(0)){ mO1=mStored=v; }

	const Tv& operator()();				///< Returns filtered output using stored value
	const Tv& operator()(Tv in);		///< Returns filtered output from input value

	void operator  = (Tv val);			///< Stores input value for operator()
	void operator *= (Tv val);			///< Multiplies stored value by value

	const Tv& last() const;				///< Returns last output
	const Tv& stored() const;			///< Returns stored value
	Tv& stored();						///< Returns stored value

	bool zeroing(Tv eps=0.0001) const;	///< Returns whether the filter is outputting zeros
	
	void onDomainChange(double r);

protected:
	FilterType mType;
	Tp mFreq, mA0, mB1;
	Tp mFreqToAng;
	Tv mStored, mO1;
};




// Implementation_______________________________________________________________

//---- AllPass1
template <class Tv, class Tp, class Td>
AllPass1<Tv,Tp,Td>::AllPass1(Tp frq)
:	d1(Tv(0)), mFreq(frq)
{
	onDomainChange(1);
}

template <class Tv, class Tp, class Td>
inline void AllPass1<Tv,Tp,Td>::freq(Tp v){
	mFreq = v;
	c = tan( freqToRad(v, Td::ups()) - M_PI_4); // valid ang in [-pi/4, pi/4]
}

template <class Tv, class Tp, class Td>
inline void AllPass1<Tv,Tp,Td>::freqF(Tp v){
	mFreq = v;
	c = freqToRad(v, Td::ups());
	c = Tp(1.27323954474) * c - Tp(1);		// 4/pi * c - 1, linear apx
	c = c * (Tp(0.76) + Tp(0.24) * c * c);
}

template <class Tv, class Tp, class Td> 
inline Tv AllPass1<Tv,Tp,Td>::operator()(Tv i0){
/*
	// Direct Form I
	//o0 = i0*c + i1 - o1*c;
	o0 = (i0 - o1)*c + i1;
	i1 = i0;
	o1 = o0;
*/
	i0 -= d1 * c;			// o0 = c * (i0 - c * d1) + d1
	Tv o0 = i0 * c + d1;	//    = c * i0 + (1 - c*c) * d1
	d1 = i0;
	return o0;
}

template <class Tv, class Tp, class Td>
inline Tv AllPass1<Tv,Tp,Td>::high(Tv i0){ return (i0 - operator()(i0)) * Tv(0.5); }

template <class Tv, class Tp, class Td>
inline Tv AllPass1<Tv,Tp,Td>::low (Tv i0){ return (i0 + operator()(i0)) * Tv(0.5); }

template <class Tv, class Tp, class Td>
inline Tp AllPass1<Tv,Tp,Td>::freq(){ return mFreq; }

template <class Tv, class Tp, class Td>
void AllPass1<Tv,Tp,Td>::onDomainChange(double /*r*/){ freq(mFreq); }


//---- Biquad
template <class Tv, class Tp, class Td>
Biquad<Tv,Tp,Td>::Biquad(Tp frq, Tp res, FilterType type, Tp lvl)
:	d1(0), d2(0), mType(type)
{
	mFreqToAng = getFreqToAng(Td::ups());
	levelNoResUpdate(lvl);
	set(frq, res);
}

template <class Tv, class Tp, class Td>
void Biquad<Tv,Tp,Td>::onDomainChange(double /*r*/){
	mFreqToAng = getFreqToAng(Td::ups());
	freq(mFreq);
}

template <class Tv, class Tp, class Td>
Tp Biquad<Tv,Tp,Td>::freq() const { return mFreq; }

template <class Tv, class Tp, class Td>
Tp Biquad<Tv,Tp,Td>::res() const { return Tp(0.5)/mResRecip; }

template <class Tv, class Tp, class Td>
Tp Biquad<Tv,Tp,Td>::level() const { return mLevel; }

template <class Tv, class Tp, class Td>
FilterType Biquad<Tv,Tp,Td>::type() const { return mType; }

template <class Tv, class Tp, class Td>
void Biquad<Tv,Tp,Td>::set(Tp freq_, Tp res_, FilterType type_){
	mType = type_;
	mResRecip = Tp(0.5)/res_;
	freq(freq_); // updates coefs based on mResRecip
}

template <class Tv, class Tp, class Td>
inline void Biquad<Tv,Tp,Td>::set(Tp frq, Tp res){ set(frq, res, mType); }

template <class Tv, class Tp, class Td>
void Biquad<Tv,Tp,Td>::zero(){ d1=d2=Tv(0); }

template <class Tv, class Tp, class Td>
void Biquad<Tv,Tp,Td>::coef(Tp a0, Tp a1, Tp a2, Tp b1, Tp b2){
	mA[0]=a0; mA[1]=a1; mA[2]=a2; mB[1]=b1; mB[2]=b2;
}

template <class Tv, class Tp, class Td>
inline void Biquad<Tv,Tp,Td>::freq(Tp v){
	mFreq = v;
	getPole(mReal, mImag, mFreq * mFreqToAng);
	resRecip(mResRecip);
}

template <class Tv, class Tp, class Td>
inline void Biquad<Tv,Tp,Td>::levelNoResUpdate(Tp v){
	mLevel=v;
	switch(mType){
	case PEAKING:
		mBeta = Tp(1)/(std::max(mLevel, Tp(0.0001)));
		break;
	case LOW_SHELF:
	case HIGH_SHELF:
		mBeta = Tp(2)*std::pow(mLevel, Tp(0.25));
		break;
	default:
		mBeta = Tp(1);
	}
}

template <class Tv, class Tp, class Td>
inline void Biquad<Tv,Tp,Td>::level(Tp v){
	levelNoResUpdate(v);
	resRecip(mResRecip);
}

/*
// Bandwidth, in octaves
template <class Tv, class Tp, class Td>
inline void Biquad<Tv,Tp,Td>::width(Tp v){
  alpha = sin(w0)/(2*Q)                                       (case: Q)
		= sin(w0)*sinh( ln(2)/2 * BW * w0/sin(w0) )           (case: BW)
        = sin(w0)/2 * sqrt( (A + 1/A)*(1/S - 1) + 2 )         (case: S)

        FYI: The relationship between bandwidth and Q is
             1/Q = 2*sinh(ln(2)/2*BW*w0/sin(w0))     (digital filter w BLT)
        or   1/Q = 2*sinh(ln(2)/2*BW)             (analog filter prototype)

	Tp oneOverQ = 2*sinh(0.34657359028 * v * mFreq*mFrqToRad / mImag);
	qRecip(oneOverQ);
}
*/

template <class Tv, class Tp, class Td>
inline void Biquad<Tv,Tp,Td>::resRecip(Tp v){
	mResRecip = v; // 1 / (2 res)
	mAlpha = mImag * mResRecip * mBeta;

	switch(mType){
	case LOW_SHELF: case HIGH_SHELF: case THRU_PASS: break; // coefs computed in type()
	default:
		// Note: b_0 is assumed to be equal to 1 in the difference equation.
		// For this reason, we divide all other coefficients by b_0.
		mB[0] = Tp(1) / (Tp(1) + mAlpha);	// 1/b_0
		mB[1] = Tp(-2) * mReal * mB[0];
		mB[2] = (Tp(1) - mAlpha) * mB[0];
	}

	type(mType);
}

template <class Tv, class Tp, class Td>
inline void Biquad<Tv,Tp,Td>::res(Tp v){
	resRecip(Tp(0.5)/v);
}

template <class Tv, class Tp, class Td>
inline void Biquad<Tv,Tp,Td>::type(FilterType typeA){
	mType = typeA;
	
	switch(mType){
	case LOW_PASS:
		mA[1] = (Tp(1) - mReal) * mB[0];
		mA[0] = mA[1] * Tp(0.5);
		mA[2] = mA[0];
		break;
	case HIGH_PASS: // low-pass with odd k a_k and freq flipped
		mA[1] = (Tp(-1) - mReal) * mB[0];
		mA[0] = mA[1] * Tp(-0.5);
		mA[2] = mA[0];
		break;
	case RESONANT:
		mA[0] = mImag * Tp(0.5) * mB[0];
		mA[1] = Tp(0);
		mA[2] =-mA[0];
		break;
	case BAND_PASS:
		mA[0] = mAlpha * mB[0];
		mA[1] = Tp(0);
		mA[2] =-mA[0];
		break;
	case BAND_REJECT:
		mA[0] = mB[0];
		mA[1] = mB[1];
		mA[2] = mB[0];
		break;
	case ALL_PASS:
		mA[0] = mB[2];
		mA[1] = mB[1];
		mA[2] = Tp(1);
		break;
	case PEAKING:{
		Tp alpha_A_b0 = mAlpha * mLevel * mB[0];
		mA[0] = mB[0] + alpha_A_b0;
		mA[1] = mB[1];
		mA[2] = mB[0] - alpha_A_b0;
		}
		break;
	case LOW_SHELF:{
		Tp A = mBeta*mBeta*Tp(0.25); // sqrt(mLevel)
		Tp Ap1 = A + Tp(1), Am1 = A - Tp(1);
		mB[0] =    Tp(1)/(Ap1 + Am1*mReal + mAlpha); // 1/b_0
		mB[1] =   Tp(-2)*(Am1 + Ap1*mReal         ) * mB[0];
		mB[2] =          (Ap1 + Am1*mReal - mAlpha) * mB[0];
		mA[0] =        A*(Ap1 - Am1*mReal + mAlpha) * mB[0];
		mA[1] = Tp( 2)*A*(Am1 - Ap1*mReal         ) * mB[0];
		mA[2] =        A*(Ap1 - Am1*mReal - mAlpha) * mB[0];
		}
		break;
	case HIGH_SHELF:{
		Tp A = mBeta*mBeta*Tp(0.25); // sqrt(mLevel)
		Tp Ap1 = A + Tp(1), Am1 = A - Tp(1);
		mB[0] =   Tp(1)/(Ap1 - Am1*mReal + mAlpha); // 1/b_0
		mB[1] =   Tp(2)*(Am1 - Ap1*mReal         ) * mB[0];
		mB[2] =         (Ap1 - Am1*mReal - mAlpha) * mB[0];
		mA[0] =       A*(Ap1 + Am1*mReal + mAlpha) * mB[0];
		mA[1] =Tp(-2)*A*(Am1 + Ap1*mReal         ) * mB[0];
		mA[2] =       A*(Ap1 + Am1*mReal - mAlpha) * mB[0];
		}
		break;
	case THRU_PASS:
		mB[0] = Tp(1);
		mB[1] = mB[2] = Tp(0);
		mA[0] = Tp(1); // = mLevel for free gain?
		mA[1] = mA[2] = Tp(0);
		break;
	default:;
	};
}

template <class Tv, class Tp, class Td>
inline Tv Biquad<Tv,Tp,Td>::operator()(Tv i0){
	// Direct form II
	i0 = i0 - d1*mB[1] - d2*mB[2];
	Tv o0 = i0*mA[0] + d1*mA[1] + d2*mA[2];
	d2 = d1; d1 = i0;
	return o0;
}

template <class Tv, class Tp, class Td>
inline Tv Biquad<Tv,Tp,Td>::nextBP(Tv i0){
	i0 = i0 - d1*mB[1] - d2*mB[2];	
	Tv o0 = (i0 - d2)*mA[0];
	d2 = d1; d1 = i0;
	return o0;
}


//---- OnePole
template <class Tv, class Tp, class Td>
OnePole<Tv,Tp,Td>::OnePole(Tp frq, FilterType type, const Tv& stored)
:	mType(type), mFreq(frq), mStored(stored), mO1(stored)
{
	onDomainChange(1);
}

template <class Tv, class Tp, class Td>
OnePole<Tv,Tp,Td>::OnePole(Tp frq, const Tv& stored)
:	OnePole(frq, LOW_PASS, stored)
{}

template <class Tv, class Tp, class Td>
void OnePole<Tv,Tp,Td>::onDomainChange(double /*r*/){
	mFreqToAng = getFreqToAng(Td::ups());
	freq(mFreq);
}

template <class Tv, class Tp, class Td>
inline void OnePole<Tv,Tp,Td>::type(FilterType v){
	mType = v;
	freq(mFreq);
}

template <class Tv, class Tp, class Td>
inline void OnePole<Tv,Tp,Td>::freq(Tp v){
	mFreq = v;

	// |H(w)| = |a0| / sqrt(1 + b1^2 + 2 b1 cos w)

	switch(mType){
	default:{ // low-pass
		// cutoff based on pole at DC (inaccurate with large bandwidth)
		//mB1 = poleRadius(Tp(2) * v * Td::ups());
		// b1 found by setting |H(w)| = 0.707
		Tp re = getPoleRe(v * mFreqToAng);
		Tp p1 = re - Tp(2);
		mB1 = -(p1 + sqrt(p1*p1 - Tp(1)));
		mA0 = Tp(1) - mB1;
		}
		break;
	case HIGH_PASS:{
		// cutoff based on pole at Nyquist (inaccurate with large bandwidth)
		//mB1 = -poleRadius(Tp(1) - Tp(2) * v * Td::ups());
		// b1 found by setting |H(w)| = 0.707
		Tp re = getPoleRe(v * mFreqToAng);
		Tp p1 = -Tp(2) - re; // -re flips cutoff
		mB1 = (p1 + sqrt(p1*p1 - Tp(1)));
		mA0 = Tp(1) + mB1;
		//printf("%g\n", mB1);
		}
		break;
	case SMOOTHING:
		lag(Tp(1)/mFreq);
		break;
	case THRU_PASS:
		mB1 = Tp(0);
		mA0 = Tp(1);
		break;
	}

	//printf("In OnePole::freq: mB1 = %g, mA0 = %g\n", mB1, mA0);
}

template <class Tv, class Tp, class Td>
inline void OnePole<Tv,Tp,Td>::lag(Tp length, Tp thresh){
	mType = SMOOTHING;
	if(length > Tp(0)){
		mFreq = Tp(1)/length;
		// Since b^(length f_s) = thresh,
		// 	b = thresh ^ (1 / (length f_s))
		mB1 = pow(thresh, Tp(Td::ups()*mFreq));
	} else {
		mFreq = Td::spu()*0.5;
		mB1 = Tp(0);
	}
	mA0 = Tp(1) - mB1;
}

template <class Tv, class Tp, class Td>
inline void OnePole<Tv,Tp,Td>::smooth(Tp v){ mB1=v; mA0=Tp(1) - scl::abs(v); }

template <class Tv, class Tp, class Td>
inline const Tv& OnePole<Tv,Tp,Td>::operator()(){ return (*this)(mStored); }

template <class Tv, class Tp, class Td>
inline const Tv& OnePole<Tv,Tp,Td>::operator()(Tv in){
	mO1 = mO1*mB1 + in*mA0;
	return mO1;
}

template <class Tv, class Tp, class Td>
inline void OnePole<Tv,Tp,Td>::operator  = (Tv v){ mStored  = v; }

template <class Tv, class Tp, class Td>
inline void OnePole<Tv,Tp,Td>::operator *= (Tv v){ mStored *= v; }

template <class Tv, class Tp, class Td>
inline const Tv& OnePole<Tv,Tp,Td>::last() const { return mO1; }

template <class Tv, class Tp, class Td>
inline const Tv& OnePole<Tv,Tp,Td>::stored() const { return mStored; }

template <class Tv, class Tp, class Td>
inline Tv& OnePole<Tv,Tp,Td>::stored(){ return mStored; }

template <class Tv, class Tp, class Td>
inline bool OnePole<Tv,Tp,Td>::zeroing(Tv eps) const {
	return scl::abs(mO1) < eps && mStored == Tv(0);
}
} // gam::
#endif
