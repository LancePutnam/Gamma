#ifndef GAMMA_OSCILLATOR_H_INC
#define GAMMA_OSCILLATOR_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include "Gamma/gen.h"
#include "Gamma/scl.h"
#include "Gamma/tbl.h"
#include "Gamma/Strategy.h"
#include "Gamma/Sync.h"
#include "Gamma/Types.h"

namespace gam{


/// Fixed-point phase accumulator

/// This is a linear phase accumulator that uses integer (fixed-point) 
/// arithmetic. The advantage of using fixed-point versus floating-point is that
/// the phase is wrapped automatically when the integer overflows. As long as
/// we used unsigned integers, this wrapping behavior is well-defined- all 
/// results of addition are taken modulo the maximum size of the integer.
/// \tparam Stap	Read tap strategy (tap::Clip, tap::Fold, tap::Rep, or tap::Wrap)
/// \tparam Ts		Synced type
template <class Stap=tap::Wrap, class Ts=Synced>
class Accum : public Ts {
public:

	/// @param[in] frq		Frequency
	/// @param[in] phs		Phase in [0, 1)
	Accum(float frq=0, float phs=0);


	void freq(float v);				///< Set frequency
	void phase(float v);			///< Set phase from [0, 1) of one period
	void phaseMax();				///< Set phase to maximum value
	void phaseAdd(float v);			///< Add value to phase [0, 1)
	void period(float v);			///< Set period length
	void reset(){ mPhaseI=0; mTap.reset(); }	///< Reset phase accumulator
	Stap& tap(){ return mTap; }

	/// Returns true if tap is done
	bool done() const { return mTap.done(phaseI()); }

	float freq() const;				///< Get frequency
	uint32_t freqI() const;			///< Get fixed-point frequency
	float freqUnit() const;			///< Get frequency in [0, 1)
	float phase() const;			///< Get phase in [0, 1)
	uint32_t phaseI() const;		///< Get fixed-point phase

	/// Returns 0x80000000 on phase wrap, 0 otherwise
	
	/// The return value can be used as a bool.  It's an integer because it
	/// saves a conditional check converting to a bool.
	uint32_t cycle();
	uint32_t operator()();			///< Alias of cycle()

	uint32_t nextPhase();			///< Increment phase and return pre-incremented phase
	uint32_t nextPhase(float freqOffset);
	uint32_t nextPhasePost();		///< Increment phase and return post-incremented phase
	uint32_t cycles();				///< Get 1 to 0 transitions of all accumulator bits
	uint32_t once();

	virtual void onResync(double r);

protected:
	float mFreq;		// Current frequency
	uint32_t mPhaseI;	// Current fixed-point phase in [0, 2^32)
	uint32_t mFreqI;	// Current fixed-point frequency
	Stap mTap;
	
	uint32_t mapFI(float v) const;	// convert unit floating-point to fixed-point integer
	float mapIF(uint32_t v) const;	// convert fixed-point integer to unit floating-point
	uint32_t mapFreq(float v) const;
};

#define ACCUM_INHERIT\
	using Accum<Stap,Ts>::phaseI;\
	using Accum<Stap,Ts>::freqI;\
	using Accum<Stap,Ts>::nextPhase;\
	using Accum<Stap,Ts>::nextPhasePost;



/// Linear sweep in interval [0,1)
template <class Stap=tap::Wrap, class Ts=Synced>
class Sweep : public Accum<Stap, Ts> {
public:
	/// @param[in] frq		Frequency
	/// @param[in] phs		Phase in [0,1)
	Sweep(float frq=440, float phs=0): Base(frq, phs){}

	float operator()(){ Base::cycle(); return Base::phase(); }

private: typedef Accum<Stap,Ts> Base;
};



/// Floating-point phase accumulator with output in [-pi, pi).
template <class Tv=gam::real, class Ts=Synced>
class AccumPhase : public Ts{
public:
	/// @param[in]	frq		Frequency
	/// @param[in]	phs		Phase in [0, 1)
	AccumPhase(Tv frq=440, Tv phs=0);

	
	/// Generate next sample. Stored phase is post-incremented.
	Tv nextPhase();
	
	/// Generate next sample with a frequency offset
	Tv nextPhase(Tv frqOffset);

	void freq(Tv v);		///< Set frequency
	void period(Tv v);		///< Set period length
	void phase(Tv v);		///< Set phase from [0, 1) of one period
	void phaseAdd(Tv v);	///< Add value to unit phase
	
	Tv freq();				///< Get frequency
	Tv period();			///< Get period
	Tv phase();				///< Get normalized phase in [0, 1)
	
	virtual void onResync(double r);
	void print(const char * append = "\n", FILE * fp = stdout);
	
protected:
	Tv mFreq, mPhase;
	Tv m2PiUPS;
	void recache();
	Tv mapFreq(Tv v) const;
	Tv mapPhase(Tv v) const;
	Tv nextPhaseUsing(Tv frq);
};



/// Tabulated function oscillator

/// This generator produces a periodic signal by reading values from a table.
/// Its advantage over other types of waveform generators is that it can
/// produce arbitrary periodic waveforms with a fixed computational cost. Its
/// main weakness is lack of parametric control over the waveform timbre.
/// This generator is named after the generator of the same name in the MUSIC
/// series of compilers. [Mathews, M. (1969). The Technology of Computer Music. 
/// The M.I.T. Press, Boston.]
/// \tparam Tv		table element type
/// \tparam Sipol	interpolation strategy
/// \tparam Stap	table reading strategy
template <class Tv=gam::real, template<class> class Sipol=ipl::Linear, class Stap=tap::Wrap, class Ts=Synced>
class Osc : public Accum<Stap,Ts>, public ArrayPow2<Tv>{
public:

	/// Constructor that allocates an internal table

	/// @param[in]	frq			Frequency
	/// @param[in]	phs			Phase in [0, 1)
	/// @param[in]	size		Size of table (actual number is power of 2 ceiling)
	Osc(float frq=440, float phs=0, uint32_t size=512)
	:	Base(frq, phs), ArrayPow2<Tv>(size)
	{}

	/// Constructor that references an external table

	/// @param[in]	frq			Frequency
	/// @param[in]	phs			Phase in [0, 1)
	/// @param[in]	src			A table to use as a reference
	Osc(float frq, float phs, ArrayPow2<Tv>& src)
	:	Base(frq, phs), ArrayPow2<Tv>(src.elems(), src.size())
	{}


	/// Generate next sample
	Tv operator()(){
		Tv r = val(); this->nextPhasePost(); return r;
	}

	/// Get current value
	Tv val() const { return mIpol(*this, phaseI()); }
	
	/// Add sine to table
	
	/// @param[in] cycles	number of cycles
	/// @param[in] amp		amplitude
	/// @param[in] phs		unit phase, [0, 1)
	Osc& addSine(double cycles, double amp=1, double phs=0){
		double f = cycles/this->size();
		for(unsigned i=0; i<this->size(); ++i){
			double p = (f*i + phs)*M_2PI;
			(*this)[i] += sin(p) * amp;
		}
		return *this;
	}

	/// Zero table elements
	void zero(){ for(unsigned i=0; i<this->size(); ++i) (*this)[i] = Tv(0); }

	ArrayPow2<Tv>& table(){ return *this; }

//	using ArrayPow2<Tv>::elems; using ArrayPow2<Tv>::size;
protected:
	Sipol<Tv> mIpol;
private:
	ACCUM_INHERIT
	typedef Accum<Stap,Ts> Base;
};



/// Complex sinusoid oscillator

/// This oscillator outputs a sine and cosine wave simultaneously.  Efficiency 
/// wise, it's comparable to a non-interpolating table oscillator, but gives a 
/// quadrature (90 degree phase shifted) wave for free.  The sinusoids
/// are computed by multiplying two complex numbers to rotate a phasor.
/// This requires only four multiplications and two additions per iteration.
/// The main disadvantage of this oscillator is that it is expensive to change
/// its frequency.  This is implemented from Mathews, M., Smith, J. 2003.
/// "Methods for synthesizing very high Q parametrically well behaved two pole 
/// filters."
template<class Tv=gam::real, class Ts=Synced>
class CSine : public Ts{
public:

	typedef Complex<Tv> complex;

	/// @param[in] frq		Frequency
	/// @param[in] amp		Amplitude
	/// @param[in] dcy		Decay time (negative means no decay)
	/// @param[in] phs		Phase in [0, 1)
	CSine(Tv frq=440, Tv amp=1, Tv dcy=-1, Tv phs=0);


	complex val;			///< Current complex output
	
	void amp(Tv val);		///< Set amplitude
	void decay(Tv val);		///< Set number of units to decay -60 dB. Negative = no decay.
	void freq(Tv val);		///< Set frequency
	void freq(const complex& val){ mInc=val; }
	void phase(Tv radians);	///< Set phase
	void reset();			///< Resets amplitude and sets phase to 0
	void set(Tv frq, Tv phs, Tv amp, Tv dcy);

	complex operator()();	///< Iterate and return current complex output.
	Tv freq();				///< Return frequency.

	virtual void onResync(double r);

protected:
	Tv mAmp, mFreq, mDcy60;
	Tv mDcy;				// phasor amp
	complex mInc;		// rotation increment
};



/// Computed sine wave oscillator.

/// This oscillator uses a 7th order Taylor series approximation to compute sine 
/// values. Computation time is about as much as a linearly-interpolating table 
/// lookup. In most cases, Taylor series are also more spectrally pure than 
/// table lookup methods since the distortion arises as harmonics.
template<class Tv=gam::real, class Ts=Synced>
class Sine : public AccumPhase<Tv,Ts> {
public:
	/// @param[in]	frq		Frequency
	/// @param[in]	phs		Phase in [0, 1)
	Sine(Tv frq=440, Tv phs=0) : AccumPhase<Tv,Ts>(frq, phs){}
	
	/// Generate next sample with a frequency offset
	Tv operator()(Tv frqOffset = Tv(0)){
		return scl::sinP9(this->nextPhase(frqOffset) * M_1_PI);
	}
};



/// Sine oscillator based on an efficient recursion equation.

/// This oscillator is based on a recursion equation requiring only one
/// multiply and add per sample computation. While calculation is fast, frequency
/// and phase updates are rather expensive and 64-bit precision is required 
/// to prevent growing or decaying in amplitude over time.  This generator is 
/// ideal in situations where a stationary sinusoid is all that is required, 
/// e.g. a grain or modulator.
template <class Tv=double, class Ts=Synced>
class SineR : public gen::RSin<Tv>, Ts{
public:

	/// @param[in]	frq		Frequency
	/// @param[in]	amp		Amplitude
	/// @param[in]	phs		Phase in [0, 1)
	SineR(Tv frq=440, Tv amp=1, Tv phs=0){ set(frq, amp, phs); }

	/// Get frequency
	Tv freq() const { return Base::freq() * Ts::spu(); }

	/// Set amplitude and phase
	void ampPhase(Tv a=1, Tv p=0){ set(freq(), a, p); }

	void freq(const Tv& v){ Base::freq(v*Ts::ups()); }

	/// Set all control parameters
	void set(Tv frq, Tv amp, Tv phs=0){ Base::set(frq*Ts::ups(), phs, amp); }


	virtual void onResync(double ratio){ Base::freq(Base::freq()/ratio); }

private:
	typedef gen::RSin<Tv> Base;
};



/// Multiple SineRs

/// For efficiency reasons, this object does not keep its frequencies synchronized
/// with the sample rate. If the sample rate changes, each oscillator must be
/// manually re-set.
template <class Tv=double, class Ts=Synced>
class SineRs : public Array<SineR<Tv, Synced1> >, Ts{
public:

	/// @param[in]	num		Number of resonators
	SineRs(uint32_t num): Base(num){ Ts::initSynced(); }

	/// Generate next sum of all oscillators
	Tv operator()(){
		Tv r= Tv(0);
		for(uint32_t j=0; j<this->size(); ++j) r+=(*this)[j]();
		return r;
	}

	/// Get last output of oscillator i
	Tv last(uint32_t i) const { return (*this)[i].val; }

	/// Set all control parameters of oscillator i
	void set(uint32_t i, Tv frq, Tv amp=1, Tv phs=0){
		(*this)[i].set(frq*Ts::ups(), amp, phs); }

private:
	typedef Array<SineR<Tv, Synced1> > Base;
};



/// Damped sine oscillator based on an efficient recursion equation.

/// This oscillator is similar to SineR, however, it has an extra multiply
/// in its sample generation to allow the oscillator to decay.
template <class Tv=double, class Ts=Synced>
class SineD : public gen::RSin2<Tv>, Ts{
public:

	/// @param[in]	frq		Frequency
	/// @param[in]	amp		Amplitude
	/// @param[in]	dcy		T60 decay length
	/// @param[in]	phs		Phase in [0, 1)
	SineD(Tv frq=440, Tv amp=1, Tv dcy=-1, Tv phs=0){ set(frq, amp, dcy, phs); }

	/// Get frequency
	Tv freq() const { return Base::freq() * Ts::spu(); }

	/// Set amplitude and phase
	void ampPhase(Tv a=1, Tv p=0){ set(freq(), a, Base::decay(), p); }
	
	/// Set all control parameters
	void set(Tv frq, Tv amp, Tv dcy, Tv phs=0){
		Base::set(frq*Ts::ups(), phs, dcy > Tv(0) ? (Tv)scl::radius60(dcy, Ts::ups()) : Tv(1), amp);
	}

	virtual void onResync(double ratio){
		Base::freq(Base::freq()/ratio);
		//printf("%g\n", Base::decay());
		//printf("%g %g %g\n", Base::decay(), ratio, ::pow(Base::decay(), 1./ratio));
		// double radius60(double dcy, double ups){ return ::exp(M_LN001/dcy * ups); }
		//		same as (0.001)^(ups/dcy)
		Base::decay(::pow(Base::decay(), 1./ratio));
	}

private:
	typedef gen::RSin2<Tv> Base;
};



/// Multiple SineDs

/// For efficiency reasons, this object does not keep its frequencies synchronized
/// with the sample rate. If the sample rate changes, each oscillator must be
/// manually re-set.
template <class Tv=double, class Ts=Synced>
class SineDs : public Array<SineD<Tv, Synced1> >, Ts{
public:

	/// @param[in]	num		Number of resonators
	SineDs(uint32_t num): Base(num){
		Ts::initSynced(); 
		for(uint32_t i=0; i<num; ++i) set(i, 0,0,0);
	}

	/// Generate next sum of all oscillators
	Tv operator()(){
		Tv r=Tv(0);
		for(uint32_t j=0; j<this->size(); ++j) r+=(*this)[j]();
		return r;
	}

	/// Get last output of oscillator i
	Tv last(uint32_t i) const { return (*this)[i].val; }

	/// Set all control parameters of oscillator i
	void set(uint32_t i, Tv frq, Tv amp, Tv dcy, Tv phs=0){
		(*this)[i].set(frq*Ts::ups(), amp, dcy*Ts::spu(), phs); }

private:
	typedef Array<SineD<Tv, Synced1> > Base;
};



/// Lookup table sine oscillator

/// This oscillator looks up values in a table containing the sine function
/// in [0, pi/2]. Doing a non-interpolating table lookup is very fast and
/// stable compared to other methods. The downsides are that the waveform is
/// generally not as spectrally pure and additional memory needs to be allocated
/// to store the lookup table (although it's relatively small and only allocated
/// once).
template <class Stap=tap::Wrap, class Ts=Synced>
class TableSine : public Accum<Stap,Ts> {
public:

	/// @param[in]	frq		Frequency
	/// @param[in]	phase	Phase in [0, 1)
	TableSine(float frq=440, float phase=0);

	float operator()(float freqOffset=0);	///< Return next linearly-interpolated sample
	float nextN(float freqOffset=0);		///< Return next non-interpolated sample
	float nextL(float freqOffset=0);		///< Return next linearly-interpolated sample

	/// Resize global sine table
	
	/// @param[in] bits		set effective table size to be 2^bits
	///
	/// This sets the effective table size with only one quarter the amount of
	/// memory actually being allocated. For example, if bits=10, the effective 
	/// table size is 2^10 = 1024, but the amount of allocated memory is only
	/// 1024/4 = 256. This call is not thread safe.
	static void resize(uint32_t bits);

protected:
//	static ArrayPow2<float> cTable; // can't use because need 2**N+1 table
	static float * cTable;		// Reference to my sample table. Must be 1<<bits.
	static uint32_t cTblBits;
	static uint32_t cFracBits;	// # of bits in fractional part of accumulator
	static uint32_t cOneIndex;
private:
	typedef Accum<Stap,Ts> Base;
	ACCUM_INHERIT
};




/// Low-frequency oscillator (non band-limited).

/// This object generates various waveform types by mapping the output of a 
/// an accumulator through mathematical functions.
template <class Stap=tap::Wrap, class Ts=Synced>
class LFO : public Accum<Stap,Ts>{
public:

	LFO();
	
	/// @param[in] frq		Frequency
	/// @param[in] phase	Phase in [0, 1)
	/// @param[in] mod		Modifier amount in [0, 1)
	LFO(float frq, float phase=0, float mod=0.5);

	uint32_t modi;			///< Modifier parameter

	/// Set frequency, phase and modifier amount
	LFO& set(float f, float p, float m);

	LFO& mod(double n);	///< Sets modifier parameter of waveform from unit value


	float cos();		///< Cosine based on 3rd order polynomial
	float down();		///< Downward ramp (1 to -1)
	float even3();		///< Even harmonic sine-like wave (3rd order)
	float even5();		///< Even harmonic sine-like wave (5th order)
	float imp();		///< Impulse (occurs at beginning of cycle)
	float line2();		///< 2-segment line. mod changes wave from down to tri to up
	float para();		///< Parabolic wave (triangle wave with all harmonics)
	float pulse();		///< Pulse (up + down). 'mod' controls pulse width
	float sinPara();	///< Sine-like wave constructed from parabolas (odd harmonics)
	float stair();		///< Stair (square + square). 'mod' controls pulse width
	float sqr();		///< Square (-1 to 1)
	float tri();		///< Triangle (starts at 1 goes down to -1 then up to 1)
	float up();			///< Upward ramp
	float up2();		///< Dual upward ramp (up + up). 'mod' controls pulse width.

	float cosU();		///< Unipolar cosine based on 3rd order polynomial
	float downU();		///< Unipolar downward ramp
	float hann();		///< Hann window
	float line2U();		///< Unipolar line2
	float paraU();		///< Unipolar parabolic wave
	float pulseU();		///< Unipolar pulse
	float stairU();		///< Unipolar stair
	float sqrU();		///< Unipolar square
	float triU();		///< Unipolar triangle (starts at 1 going down then up)
	float upU();		///< Unipolar upward ramp
	float up2U();		///< Unipolar upward ramp2

	float patU();
	float patU(uint32_t mul);
	float sineT9();
	float sineP9();

	bool seq();			// Returns 'mod' as sequence of triggers

	ACCUM_INHERIT
private:
	typedef Accum<Stap,Ts> Base;
};



/// Base class for producing band-limited waveforms

/// This produces a finite sum of equi-amplitude cosine waves that approach the
/// shape of a periodic impulse train. 
/// Due to numerical issues, this generator should not be used for producing 
/// very low frequency modulation signals. For that purpose, it is better to use
/// the LFO class.
template<class Tv=gam::real, class Ts=Synced>
class Buzz : public AccumPhase<Tv,Ts> {
public:

	/// @param[in]	frq			Frequency
	/// @param[in]	phase		Phase in [0, 1)
	/// @param[in]	harmonics	Number of harmonics
	Buzz(Tv frq=440, Tv phase=0, Tv harmonics=8);
	virtual ~Buzz(){}

	void antialias();			///< Adjust number of harmonics to prevent aliasing
	void harmonics(Tv num);		///< Set number of harmonics
	void harmonicsMax();		///< Set number of harmonics to fill Nyquist range

	Tv operator()();			///< Returns next sample of all harmonic impulse
	Tv odd();					///< Returns next sample of odd harmonic impulse
	Tv saw(Tv intg=0.997);		///< Returns next sample of saw waveform
	Tv square(Tv intg=0.997);	///< Returns next sample of square waveform
	
	Tv maxHarmonics();			///< Get number of harmonics below Nyquist based on current settings

	virtual void onResync(double r);

protected:
	Tv mAmp;			// amplitude normalization factor
	Tv mN;				// # harmonics
	Tv mNDesired;		// desired number of harmonics
	Tv mNFrac;		
	Tv mSPU_2;			// cached locals
	Tv mPrev;			// previous output for integration
	void setAmp();
private: typedef AccumPhase<Tv,Ts> Base;
};



/// Band-limited impulse train

/// This produces a Fourier representation of an impulse train where the number of
/// harmonics is adjusted automatically to prevent aliasing.
/// Due to numerical issues, this generator should not be used for producing 
/// very low frequency modulation signals. For that purpose, it is better to use
/// the LFO class.
template <class Tv=gam::real, class Ts=Synced>
struct Impulse : public Buzz<Tv,Ts>{

private: typedef Buzz<Tv,Ts> Base;

public:
	using Base::freq;

	/// @param[in] frq		Frequency
	/// @param[in] phs		Phase, in [0, 1)
	Impulse(Tv frq=440, Tv phs=0): Base(frq, phs){ onResync(1); }

	/// Set frequency
	void freq(Tv v){ Base::freq(v); Base::harmonicsMax(); }

	virtual void onResync(double r){ Base::onResync(r); freq(AccumPhase<Tv, Ts>::freq()); }
};



/// Band-limited saw wave

/// This produces a Fourier representation of a saw wave where the number of
/// harmonics is adjusted automatically to prevent aliasing.
/// Due to numerical issues, this generator should not be used for producing 
/// very low frequency modulation signals. For that purpose, it is better to use
/// the LFO class.
template <class Tv=gam::real, class Ts=Synced>
struct Saw : public Impulse<Tv,Ts> {

	/// @param[in] frq		Frequency
	/// @param[in] phs		Phase, in [0, 1)
	Saw(Tv frq=440, Tv phs=0): Impulse<Tv, Ts>(frq, phs){}

	/// Generate next sample
	
	/// @param[in] itg		Integration amount
	///
	Tv operator()(Tv intg=0.997){ return Impulse<Tv,Ts>::saw(intg); }
};



/// Band-limited square wave

/// This produces a Fourier representation of a square wave where the number of
/// harmonics is adjusted automatically to prevent aliasing.
/// Due to numerical issues, this generator should not be used for producing 
/// very low frequency modulation signals. For that purpose, it is better to use
/// the LFO class.
template <class Tv=gam::real, class Ts=Synced>
struct Square : public Impulse<Tv,Ts> {

	/// @param[in] frq		Frequency
	/// @param[in] phs		Phase, in [0, 1)
	Square(Tv frq=440, Tv phs=0) : Impulse<Tv,Ts>(frq, phs){}

	/// Generate next sample
	
	/// @param[in] itg		Integration amount
	///
	Tv operator()(Tv intg=0.997){ return Impulse<Tv,Ts>::square(intg); }
};



/// Discrete summation formula (DSF) oscillator

/// This produces a finite set of harmonics whose amplitudes follow a geometric
/// series. The amplitude of harmonic i is ar^i where 'ar' is called the 
/// amplitude ratio. The frequency of harmonic i is (i * fr + 1) where 'fr' is
/// called the frequency ratio. Harmonics run from i=0 (the fundamental) to
/// the maximum specified harmonic.
template<class Tv=gam::real, class Ts=Synced>
class DSF : public AccumPhase<Tv,Ts> {
public:

	/// @param[in]	frq			Frequency
	/// @param[in]	freqRatio	Frequency ratio of partials
	/// @param[in]	ampRatio	Amplitude ratio of partials
	/// @param[in]	harmonics	Number of harmonics
	DSF(Tv frq=440, Tv freqRatio=1, Tv ampRatio=0.5, Tv harmonics=8);
	
	Tv operator()();			///< Generate next sample
	
	void ampRatio(Tv v);		///< Set amplitude ratio of partials
	void antialias();			///< Adjust harmonics so partials do not alias
	void freq(Tv v);			///< Set frequency of fundamental
	void freqRatio(Tv v);		///< Set frequency ratio of partials
	void harmonics(Tv v);		///< Set number of harmonics
	void harmonicsMax();		///< Set number of harmonics to fill Nyquist range

	Tv ampRatio();				///< Get amplitude ratio
	Tv freqRatio();				///< Get frequency ratio
	Tv harmonics();				///< Get current number of harmonics
	Tv maxHarmonics();			///< Get maximum number of harmonics for current settings
	
	virtual void onResync(double r);

protected:
	typedef AccumPhase<Tv,Ts> Base;

	Tv mN, mNDesired;		// actual and desired # harmonics
	Tv mFreqRatio;			// frequency ratio
	Tv mA;					// Partial amplitude ratio
	Tv mBeta, mBetaInc;		// "detune" accumulator
	Tv mAPow, mASqP1;		// cached vars
	
	void updateAPow();
	void updateBetaInc();
};



// Simple band-limited impulse generator

// This uses a fast, simplified formula for generating a band-limited impulse,
// but only operates at integer divisions of the Nyquist frequency.
class ImpulseFast : public Synced {
public:
	ImpulseFast(): mPhase(0), mOffset(0){ freq(0); }


	/// Set frequency
	void freq(double v){
		double samples = spu() / v;
		
		uint32_t period = (uint32_t)(samples);
		//period &= 0xfffffffe;		// force period to be even
									// odd periods introduce DC

		mPeriod = (double)period;
		
		if(scl::even(period))	mOffset = 0.;
		//else					mOffset = 0.5f / mPeriod;
	}


	/// Generate next sample
	float operator()(){
		float v = 0.f;

		if(mPhase >= mPeriod){
			mPhase -= mPeriod;
			v = 1.f;
		}
		else if(scl::even((uint32_t)mPhase)){
			v = -1.f/(mPeriod * 0.5f - 1.f);
		}
		
		++mPhase;
		return v + mOffset;
	}
	
protected:
	double mPhase;		// phase in samples
	double mPeriod;		// period in samples;
	float mOffset;		// DC compensation
};



// Implementation_______________________________________________________________

#define TEM template<class Tv, class Ts>
#define TEMS template<class Ts>
#define TEMTS template<class St, class Ts>

//---- Accum
#define TACCUM	Accum<St,Ts>

TEMTS TACCUM::Accum(float f, float p): mFreq(f){
	Ts::initSynced();
	(p >= 1.f) ? phaseMax() : this->phase(p);
	//printf("%u\n", mPhaseI);
}

TEMTS inline uint32_t TACCUM::mapFI(float v) const {
	//return scl::unitToUInt(v);
	//return (uint32_t)(v * 4294967296.);
	return castIntRound(v * 4294967296.);
}

TEMTS inline float TACCUM::mapIF(uint32_t v) const {
	//return v/4294967296.;
	return uintToUnit<float>(v);
}

TEMTS inline uint32_t TACCUM::mapFreq(float v) const {
	return mapFI(v * Ts::ups());
}

TEMTS void TACCUM::onResync(double r){ //printf("Accum: onSyncChange (%p)\n", this);
	freq(mFreq);
}

TEMTS inline void TACCUM::freq(float v){
	mFreq  = v;
	mFreqI = mapFreq(v);
}

TEMTS inline void TACCUM::period(float v){ freq(1.f/v); }
TEMTS inline void TACCUM::phase(float v){ mPhaseI = mapFI(v); }
TEMTS inline void TACCUM::phaseAdd(float v){ mTap(mPhaseI, mapFI(v)); }
TEMTS inline void TACCUM::phaseMax(){ mPhaseI = 0xffffffff; }

TEMTS inline float TACCUM::freq() const { return mFreq; }
TEMTS inline uint32_t TACCUM::freqI() const { return mFreqI; }
TEMTS inline float TACCUM::freqUnit() const { return mapIF(mFreqI); }
TEMTS inline float TACCUM::phase() const { return mapIF(mPhaseI); }
TEMTS inline uint32_t TACCUM::phaseI() const { return mPhaseI; }

TEMTS inline uint32_t TACCUM::nextPhase(float frqOffset){
	uint32_t r=mPhaseI;
	mTap(mPhaseI, mFreqI + mapFreq(frqOffset));
	return r;
}

TEMTS inline uint32_t TACCUM::nextPhase(){ uint32_t r=mPhaseI; nextPhasePost(); return r; }
TEMTS inline uint32_t TACCUM::nextPhasePost(){ return mTap(mPhaseI, mFreqI); }



TEMTS inline uint32_t TACCUM::operator()(){ return cycle(); }

TEMTS inline uint32_t TACCUM::cycle(){ return cycles() & 0x80000000; }

//TEMTS inline uint32_t TACCUM::cycle(uint32_t mask){
//	return cycles() & mask;
//}

TEMTS inline uint32_t TACCUM::cycles(){
	uint32_t prev = phaseI();
	nextPhasePost();	
	return ~phaseI() & prev;
}

TEMTS inline uint32_t TACCUM::once(){
	uint32_t prev = phaseI();
	uint32_t c = cycle();
	if(c) mPhaseI = prev;
	return c;
}

#undef TACCUM



//---- AccumPhase

TEM AccumPhase<Tv, Ts>::AccumPhase(Tv f, Tv p)
:	mFreq(f), m2PiUPS(1)
{
	Ts::initSynced();
	this->phase(p);
}

TEM inline Tv AccumPhase<Tv, Ts>::mapFreq(Tv v) const { return v*m2PiUPS; }
TEM inline Tv AccumPhase<Tv, Ts>::mapPhase(Tv v) const { return v*Tv(M_2PI); }

TEM inline Tv AccumPhase<Tv, Ts>::nextPhaseUsing(Tv frq){
	mPhase = scl::wrapPhase(mPhase); // guarantees that result is in [-pi, pi)
	Tv r = mPhase;
	mPhase += frq;
	return r;
}

TEM inline Tv AccumPhase<Tv, Ts>::nextPhase(){
	return nextPhaseUsing(mFreq);
}

TEM inline Tv AccumPhase<Tv, Ts>::nextPhase(Tv frqMod){
	return nextPhaseUsing(mFreq + mapFreq(frqMod));
}

TEM inline void AccumPhase<Tv, Ts>::freq(Tv v){ mFreq = mapFreq(v); }
TEM inline void AccumPhase<Tv, Ts>::period(Tv v){ freq(Tv(1)/v); }
TEM inline void AccumPhase<Tv, Ts>::phase(Tv v){ mPhase = mapPhase(v); }
TEM inline void AccumPhase<Tv, Ts>::phaseAdd(Tv v){ mPhase += mapPhase(v); }

TEM inline Tv AccumPhase<Tv, Ts>::freq(){ return mFreq/m2PiUPS; } //mFreq; }
TEM inline Tv AccumPhase<Tv, Ts>::period(){ return Tv(1) / freq(); }
TEM inline Tv AccumPhase<Tv, Ts>::phase(){ return mPhase * Tv(M_1_2PI); }

TEM void AccumPhase<Tv, Ts>::onResync(double r){ Tv f=freq(); recache(); freq(f); }
TEM void AccumPhase<Tv, Ts>::recache(){ m2PiUPS = Tv(Ts::ups() * M_2PI); }

TEM void AccumPhase<Tv, Ts>::print(const char * append, FILE * fp){
//	fprintf(fp, "%f %f %f%s", freq(), phase(), mFreqI, append);
	fprintf(fp, "%f %f %f%s", freq(), phase(), mFreq, append);
}


//---- CSine

TEM CSine<Tv, Ts>::CSine(Tv f, Tv a, Tv dcy60, Tv p)
	: val(a, 0), mAmp(a), mFreq(f), mDcy60(dcy60)
{
	Ts::initSynced();
	this->phase(p);
}

TEM void CSine<Tv, Ts>::amp(Tv v){
	if(scl::abs(mAmp) > Tv(0.000001)){
		Tv factor = v / mAmp;
		val *= factor;
	} else {
		val(v, Tv(0));
	}
	mAmp = v;
}

TEM void CSine<Tv, Ts>::decay(Tv v){
	mDcy60 = v;
	mDcy = v > Tv(0) ? Tv(scl::t60(v * Ts::spu())) : Tv(1);
	freq(mFreq);
}

TEM void CSine<Tv, Ts>::freq(Tv v){
	mFreq = v;
	Tv phaseInc = v * Ts::ups() * Tv(M_2PI);
	mInc.fromPolar(mDcy, phaseInc);
	//printf("%f %f %f %f\n", phaseInc, mDcy, c1, s1);
}

TEM void CSine<Tv, Ts>::phase(Tv v){
	// set phase without changing current magnitude
	val.fromPolar(val.norm(), v*Tv(M_2PI));
}

TEM void CSine<Tv, Ts>::reset(){ val(mAmp, Tv(0)); }

TEM void CSine<Tv, Ts>::set(Tv frq, Tv phase, Tv amp, Tv dcy60){
	mFreq = frq;
	decay(dcy60);
	this->amp(amp);
	this->phase(phase);
}

TEM inline Complex<Tv> CSine<Tv, Ts>::operator()(){
	complex c = val;
	val *= mInc;
	return c;
}

TEM inline Tv CSine<Tv, Ts>::freq(){ return mFreq; }

TEM void CSine<Tv, Ts>::onResync(double r){
	decay(mDcy60); // this sets frequency as well
}


//---- TableSine

#define TTABLESINE TableSine<St,Ts>
TEMTS uint32_t TTABLESINE::cTblBits  = 0;	
TEMTS uint32_t TTABLESINE::cFracBits = 0;
TEMTS uint32_t TTABLESINE::cOneIndex = 0;
TEMTS float * TTABLESINE::cTable = 0;

TEMTS TTABLESINE::TableSine(float f, float p): Base(f, p){
	// initialize global table ONCE
	if(0 == cTable){ resize(11); }
}

TEMTS void TTABLESINE::resize(uint32_t bits){
	if(bits != cTblBits){
		if(cTable) delete[] cTable;

		cTblBits = bits;
		cFracBits = 32UL - cTblBits;
		cOneIndex = 1<<cFracBits;
		uint32_t size = 1<<(cTblBits-2);
		cTable = new float[size + 1];
		tbl::sinusoid(cTable, size, 0, 0.25);
		cTable[size] = 1;		
	}
}

TEMTS inline float TTABLESINE::operator()(float df){ return nextL(df); }

TEMTS inline float TTABLESINE::nextN(float df){	
	return tbl::atQ(cTable, cFracBits, nextPhase(df));
}

TEMTS inline float TTABLESINE::nextL(float df){
	uint32_t P = nextPhase(df);

	return ipl::linear(
		gam::fraction(cTblBits, P),
		tbl::atQ(cTable, cFracBits, P),
		tbl::atQ(cTable, cFracBits, P + cOneIndex)
	);
}

#undef TTABLESINE



//---- LFO
#define TLFO LFO<St,Ts>
TEMTS TLFO::LFO(): Base(){ mod(0.5); }
TEMTS TLFO::LFO(float f, float p, float m): Base(f, p){ mod(m); }

TEMTS inline TLFO& TLFO::set(float f, float p, float m){ this->freq(f); this->phase(p); return mod(m); }
TEMTS inline TLFO& TLFO::mod(double v){ modi = castIntRound(v*4294967296.); return *this; }

TEMTS inline float TLFO::line2(){
	using namespace gam::scl;
	
//	// Starts at 1
//	float r1 = rampDown(phaseI());
//	float r2 = rampDown(phaseI() + modi);

	// Starts at -1 (better for creating attack/decay like envelopes)
	uint32_t m = scl::clip<uint32_t>(modi, 0xffefffff, 512); // avoid div by zero
	float r1 = rampDown(phaseI() - m);
	float r2 = rampDown(phaseI());
	float p  = rampUpU(m);

	float r = (r1*r1 - r2*r2)/(4.f*p*(1.f - p));
	nextPhasePost();
	return r;
}

TEMTS inline float TLFO::line2U(){ return line2()*0.5f+0.5f; }

#define DEF(name, exp) TEMTS inline float TLFO::name{ float r = exp; return r; }
//DEF(cos(),		tri(); r *= 0.5f * r*r - 1.5f)
//DEF(cos(),		up(); r=scl::abs(r*r) )//r = -1.f - scl::pow2(2.f*r)*(scl::abs(r)-1.5f) )
DEF(cos(),		up(); r = -1.f - r*r*(4.f*scl::abs(r)-6.f) )
DEF(down(),		scl::rampDown(nextPhase()))
DEF(even3(),	up(); static const float c=-1.50f*sqrtf(3.f); r *= (1.f-r*r)*c;)
DEF(even5(),	up(); static const float c=-1.25f*::powf(5.f,0.25f); r *= (1.f-scl::pow4(r))*c;)
DEF(imp(),		scl::pulseU(nextPhase(), this->freqI()) )
DEF(para(),		paraU()*1.5f - 0.5f)
DEF(pulse(),	scl::pulse(nextPhase(), modi))
DEF(sinPara(),	scl::sinPara(nextPhase()))
DEF(stair(),	scl::stair(nextPhase(), modi))
DEF(sqr(),		scl::square(nextPhase()))
DEF(tri(),		scl::triangle(nextPhase()))
DEF(up(),		scl::rampUp(nextPhase()))
DEF(up2(),		scl::rampUp2(nextPhase(), modi))
DEF(cosU(),		tri(); r = scl::mapSinSU(r))
DEF(downU(),	scl::rampDownU(nextPhase()))
DEF(hann(),		tri(); r = r * (0.25f * r*r - 0.75f) + 0.5f)
DEF(paraU(),	up(); r*=r;)
DEF(pulseU(),	scl::pulseU(nextPhase(), modi))
DEF(sqrU(),		scl::squareU(nextPhase()))
DEF(stairU(),	scl::stairU(nextPhase(), modi))
DEF(triU(),		scl::triangleU(nextPhase()))
DEF(upU(),		scl::rampUpU(nextPhase()))
DEF(up2U(),		scl::rampUp2U(nextPhase(), modi))
DEF(patU(),		scl::rampUpU(nextPhase() & modi))

DEF(patU(uint32_t mul), scl::rampUpU((nextPhase() & modi) * mul))

DEF(sineT9(),	up(); r = scl::sinT9(r * M_PI))
DEF(sineP9(),	up(); r = scl::sinP9(r))

#undef DEF

TEMTS inline bool TLFO::seq(){
	uint32_t prev = nextPhase();
	if( (phaseI() ^ prev) & 0xf8000000 ){
		if( (modi >> (phaseI()>>27)) & 0x1 ) return true;
	}
	return false;
}

#undef TLFO


//---- Buzz

TEM Buzz<Tv,Ts>::Buzz(Tv f, Tv p, Tv harmonics)
:	Base(f, p), mAmp(0), mPrev(Tv(0))
{
	onResync(1);
	this->harmonics(harmonics);
}

TEM inline void Buzz<Tv,Ts>::harmonics(Tv v){
	mN = mNDesired = scl::floor(v);
	setAmp();
	mNFrac = v - mN;
}

TEM inline void Buzz<Tv,Ts>::harmonicsMax(){ harmonics(maxHarmonics()); }

TEM inline void Buzz<Tv,Ts>::antialias(){
	float maxN = scl::floor(maxHarmonics());
	mN = mNDesired > maxN ? maxN : mNDesired;
	setAmp();
}

TEM inline Tv Buzz<Tv,Ts>::maxHarmonics(){ return mSPU_2 / this->freq(); }

TEM inline void Buzz<Tv,Ts>::setAmp(){
	// Normally, the amplitude is 1/(2N), but we will linearly interpolate
	// based on fractional harmonics to avoid sudden changes in amplitude to
	// the lower harmonics which is very noticeable.
	mAmp = (mN != Tv(0)) ? (Tv(0.5) / (mN+mNFrac)) : 0;
	//mAmp = (mN != Tv(0)) ? (Tv(0.5) / (mN)) : 0;
}

#define EPS 0.000001
TEM inline Tv Buzz<Tv, Ts>::operator()(){
	/*        1   / sin((N+0.5)x)    \
	   f(x) = -- |  ------------- - 1 |
	          2N  \   sin(0.5x)      /
	*/
	Tv theta = this->nextPhase();
	Tv result;
	Tv denom = scl::sinT9(theta * Tv(0.5));

	// denominator goes to zero when theta is an integer multiple of 2 pi
	if(scl::abs(denom) < Tv(EPS)){
		result = Tv(2) * mN * mAmp;
		//printf("Impulse::(): oops\n");
	}
	else{
		Tv nphase = scl::wrapPhase(theta * (mN + Tv(0.5)));
		//result = (scl::sinT7(nphase) / denom - Tv(1)) * mAmp;
		result = ((scl::sinT7(nphase) - denom) / denom) * mAmp;
	}

	//Tv fund = ((denom*denom)-0.5)*-4*mAmp;
	//result -= fund; // drop first harmonic

	// Mitigate pops when number of harmonics is changed by a small amount
	//result -= 2.f * mAmp * cos(theta * mN) * (1.f - scl::pow2(mNFrac));

	return result;
}

TEM inline Tv Buzz<Tv,Ts>::odd(){
	/*        1   / sin(2N x) \
	   f(x) = -- |  ---------  |
	          2N  \   sin(x)  /
	*/
	Tv theta = this->nextPhase();
	Tv result;

	Tv n2 = scl::roundAway(mN*0.5) * 2;	
	Tv n2frac = ((mN + mNFrac) - (n2-1)); // fraction, in [0,2), btw odd harmonics

	// cos is more precise near zero-crossings
	Tv denom = scl::cosT8(scl::wrapPhaseOnce(theta - M_PI_2));
	//Tv denom = scl::sinT9(theta);
	if( scl::abs(denom) < Tv(EPS) ){
		if( theta > M_PI ) theta -= M_2PI;
		Tv A = n2 / (n2 + n2frac);
		result = (theta > -M_PI_2 && theta < M_PI_2) ? A : -A;
		//printf("Impulse::odd(): oops\n");
	}
	else result = scl::sinT7(scl::wrapPhase(n2 * theta)) / (denom * (n2 + n2frac));
	
	return result;
}
#undef EPS

TEM inline Tv Buzz<Tv,Ts>::saw(Tv i){ return mPrev=(*this)()*0.125 + i*mPrev; }
TEM inline Tv Buzz<Tv,Ts>::square(Tv i){ return mPrev=odd()*0.125 + i*mPrev; }

TEM void Buzz<Tv,Ts>::onResync(double r){
	Base::onResync(r);
	mSPU_2 = Tv(Synced::spu() * 0.5);
}



//---- DSF

TEM DSF<Tv,Ts>::DSF(Tv frq, Tv freqRatioA, Tv ampRatioA, Tv harmonicsA)
	: Base(frq)
{
	freq(frq);
	freqRatio(freqRatioA);
	harmonics(harmonicsA);
	ampRatio(ampRatioA);

	mBeta = 0.f;
}

TEM inline void DSF<Tv,Ts>::freq(Tv v){
	Base::freq(v);
	updateBetaInc();	
}

TEM inline void DSF<Tv,Ts>::freqRatio(Tv v){
	mFreqRatio = v;
	updateBetaInc();
}

TEM inline void DSF<Tv,Ts>::ampRatio(Tv v){
	if(v != mA){
		// if near 1, nudge away
		static const Tv epslt = 0.9997;
		static const Tv epsgt = 1/epslt;
				if(v >= 1 && v < epsgt) v = epsgt;
		else	if(v <= 1 && v > epslt) v = epslt;
		mA = v;
		mASqP1 = mA * mA + 1.f;
		updateAPow();
	}
}

TEM inline void DSF<Tv,Ts>::harmonics(Tv v){
	if(v != mN){
		mN = mNDesired = v;
		updateAPow();
	}
}

TEM inline void DSF<Tv,Ts>::harmonicsMax(){ harmonics(maxHarmonics()); }

TEM inline void DSF<Tv,Ts>::antialias(){
	Tv maxN = maxHarmonics();
	if(mNDesired > maxN)	mN = maxN;
	else					mN = mNDesired;
	updateAPow();
}

TEM inline Tv DSF<Tv,Ts>::ampRatio(){ return mA; }
TEM inline Tv DSF<Tv,Ts>::freqRatio(){ return mFreqRatio; }
TEM inline Tv DSF<Tv,Ts>::harmonics(){ return mN; }

TEM inline Tv DSF<Tv,Ts>::maxHarmonics(){
	return scl::floor((Tv(this->spu()) * Tv(0.5)/this->freq() - Tv(1))/freqRatio() + Tv(1));
}

TEM inline void DSF<Tv,Ts>::updateAPow(){ mAPow = ::pow(mA, mN); }
TEM inline void DSF<Tv,Ts>::updateBetaInc(){ mBetaInc = this->mFreq * mFreqRatio; }

// Generalized DSF formula:
// sum{k=0, N}( a^k sin(T + k B) )
//		=  sin(T) - a sin(T - B) - a^(N+1) [sin(T + (N+1) B) - a sin(T + N B))]
//			/ 1 + a^2 - 2a cos(B)

#define SIN scl::sinT7
#define COS scl::cosT8
//#define SIN sin
//#define COS cos
TEM inline Tv DSF<Tv,Ts>::operator()(){
	Tv theta = Base::nextPhase();
	mBeta = scl::wrapPhase(mBeta);

	Tv phs2 = scl::wrapPhaseOnce(theta - mBeta);
	Tv phs3 = scl::wrapPhase(theta + mN * mBeta);
	Tv phs4 = scl::wrapPhaseOnce(phs3 - mBeta);

	Tv result = SIN(theta) - mA * SIN(phs2) - mAPow * (SIN(phs3) - mA * SIN(phs4));	
	result /= mASqP1 - Tv(2) * mA * COS(mBeta);
	//result /= mA * (mA - Tv(2) * COS(mBeta)) + Tv(1);

	mBeta += mBetaInc;
	return result;
}
#undef SIN
#undef COS

TEM void DSF<Tv,Ts>::onResync(double r){
	Base::onResync(r);
	freq(Base::freq());
	harmonics(mNDesired);
}

#undef TEM
#undef TEMS
#undef TEMTS

} // gam::
#endif
