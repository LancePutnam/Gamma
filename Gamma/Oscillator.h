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


/// Fixed-point phase accumulator.
template <class Stap=tap::Wrap, class Ts=Synced>
class Accum : public Ts {
public:
	Accum();

	/// @param[in] frq		Initial frequency.
	/// @param[in] phs		Initial phase [0, 1).
	Accum(float frq, float phs=0);
	virtual ~Accum(){}

	void freq(float v);				///< Set frequency.
	void phase(float u);			///< Set phase from [0, 1) of one period.
	void phaseMax();				///< Set phase to maximum value.
	void phaseAdd(float u);			///< Add value to phase [0, 1).		
	void period(float v);			///< Set period length.
	void reset(){ mTap.reset(); }	///< Reset phase accumulator

	bool done(){ return mTap.done(phaseI()); }
//	Stap& tap(){ return mTap; }

	float freq() const;				///< Returns frequency.
	float phase() const;			///< Returns current unit phase
	uint32_t phaseI() const;		///< Returns current fixed-point phase value
	float phaseInc() const;			///< Returns unit phase increment in [0, 1).
	uint32_t phaseIncI() const;		///< Returns current fixed-point phase increment value

	uint32_t operator()();			///< Alias of cycle().

	/// Returns 0x80000000 on phase wrap, 0 otherwise.
	
	/// The return value can be used as a bool.  It's an integer because it
	/// saves a conditional check converting to a bool.
	uint32_t cycle();

	uint32_t cycles();		///< Returns 1 to 0 transitions of all accumulator bits
	uint32_t once();
	uint32_t incPhase();	///< Increment phase and return post-incremented phase.

	virtual void onResync(double r);

protected:
	float mFreq;			// Current frequency.
	uint32_t mPhase;		// Current phase btw [0, 2^32)
	uint32_t mPhaseInc;
	Stap mTap;
	
	uint32_t phaseFI(float v) const;	// convert unit floating-point to fixed-point integer
	float phaseIF(uint32_t v) const;	// convert fixed-point integer to unit floating-point
};

#define ACCUM_INHERIT\
	using Accum<Stap,Ts>::phaseI;\
	using Accum<Stap,Ts>::phaseIncI;\
	using Accum<Stap,Ts>::incPhase;


/// Linear sweep in interval [0,1)
template <class Stap=tap::Wrap, class Ts=Synced>
class Sweep : public Accum<Stap, Ts> {
public:
	/// @param[in] frq		Initial frequency
	/// @param[in] phs		Initial unit phase in [0,1)
	Sweep(float frq=440, float phs=0): Base(frq, phs){}

	float operator()(){ Base::cycle(); return Base::phase(); }

private: typedef Accum<Stap,Ts> Base;
};


/// Floating point phase accumulator in [-pi, pi).
template <class Tv=gam::real, class Ts=Synced>
class AccumPhase : public Ts{
public:
	/// @param[in]	frq		Initial frequency
	/// @param[in]	phs		Initial unit phase in [0, 1)
	AccumPhase(Tv frq=440, Tv phs=0);
	virtual ~AccumPhase(){}
	
	Tv nextPhase();			///< Return next phase. Stored phase is post-incremented.

	void freq(Tv v);		///< Set frequency.
	void period(Tv v);		///< Set period length.
	void phase(Tv v);		///< Set phase from [0, 1) of one period.
	void phaseAdd(Tv v);	///< Add value to unit phase.
	
	Tv freq();				///< Return frequency.
	Tv period();			///< Return period.
	Tv phase();				///< Return phase [0, 1).
	
	virtual void onResync(double r);
	
	void print(const char * append = "\n", FILE * fp = stdout);
	
protected:
	Tv mFreq, mPhase, mPhaseInc;
	Tv m2PiUPS;
	void recache();
};



/// Tabulated function oscillator

/// Tv is the table element type, Sipol is an interpolation strategy, and
/// Stap is a table reading strategy.
/// Mathews, M. (1969). The Technology of Computer Music. The M.I.T. Press, Boston.
template <class Tv=gam::real, template<class> class Sipol=ipl::Linear, class Stap=tap::Wrap, class Ts=Synced>
class Osc : public Accum<Stap,Ts>, public ArrayPow2<Tv>{
public:

	/// Constructor that alocates an internal table

	/// @param[in]	frq			Initial frequency
	/// @param[in]	phs			Initial unit phase in [0, 1)
	/// @param[in]	size		Number of table elements (actual number is power of 2 ceiling)
	Osc(float frq, float phs=0, uint32_t size=512)
	:	Base(frq, phs), ArrayPow2<Tv>(size)
	{}

	/// Constructor that references an external table

	/// @param[in]	frq			Initial frequency
	/// @param[in]	phs			Initial unit phase in [0, 1)
	/// @param[in]	src			A table to use as a reference
	Osc(float frq, float phs, const ArrayPow2<Tv> & src)
	:	Base(frq, phs), ArrayPow2<Tv>(src.elems(), src.size())
	{}
	
	virtual ~Osc(){}

	/// Generate next sample
	Tv operator()(){
		Tv o0 = val(); mTap(this->mPhase, phaseIncI()); return o0;
	}
	
	void zero(){ for(unsigned i=0; i<this->size(); ++i) (*this)[i] = 0; }

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

//	using ArrayPow2<Tv>::elems; using ArrayPow2<Tv>::size;
protected:
	Sipol<Tv> mIpol;
private:
	ACCUM_INHERIT
	typedef Accum<Stap,Ts> Base;
};



/// Sine-cosine quadrature wave oscillator.

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
class Quadra : public Ts{
public:

	typedef Complex<Tv> complex;

	/// @param[in] frq		Initial frequency
	/// @param[in] amp		Initial amplitude
	/// @param[in] dcy		Initial decay time. Negative means no decay.
	/// @param[in] phs		Initial unit phase in [0, 1)
	Quadra(Tv frq=440, Tv amp=1, Tv dcy=-1, Tv phs=0);
	//virtual ~Quadra(){}

	complex val;			///< Current complex output
	
	void amp(Tv val);		///< Set amplitude
	void decay(Tv val);		///< Set number of units to decay -60 dB. Negative = no decay.
	void freq(Tv val);		///< Set frequency
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
class Sine : public AccumPhase<Tv, Ts> {
public:
	/// @param[in]	frq		Initial frequency
	/// @param[in]	phs		Initial unit phase in [0, 1)
	Sine(Tv frq=440, Tv phs=0) : AccumPhase<Tv, Ts>(frq, phs){}

	/// Return next sample.
	Tv operator()(){ return scl::sinP9(AccumPhase<Tv, Ts>::nextPhase() * M_1_PI); }
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
	/// @param[in]	phs		Unit phase in [0, 1)
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
	Tv operator()(){ Tv r=(Tv)0; for(uint32_t j=0; j<this->size(); ++j) r+=(*this)[j](); return r; }
	
	/// Generate next samples adding into a buffer
	template <class V>
	void add(V * dst, uint32_t n){
		for(uint32_t j=0; j<this->size(); ++j){
			for(uint32_t i=0; i<n; ++i)	dst[i] += (*this)[j]();
		}
	}

	/// Returns ith oscillator's last value
	Tv last(uint32_t i){ return (*this)[i].val; }

	/// Set all control parameters
	void set(uint32_t i, Tv frq, Tv amp, Tv phs=0){ (*this)[i].set(frq*Ts::ups(), amp, phs); }

private:
	typedef Array<SineR<Tv, Synced1> > Base;
};



/// Damped sine oscillator based on an efficient recursion equation.

/// This oscillator is similar to SineR, however, it has an extra multiply
/// in its sample generation to allow the oscillator to decay.
template <class Tv=double, class Ts=Synced>
class SineD : public gen::RSin2<Tv>, Ts{
public:

	/// @param[in]	frq		Initial frequency
	/// @param[in]	amp		Initial amplitude
	/// @param[in]	dcy		Initial T60 decay length
	/// @param[in]	phs		Initial unit phase in [0, 1)
	SineD(Tv frq=440, Tv amp=1, Tv dcy=-1, Tv phs=0){ set(frq, amp, dcy, phs); }

	/// Returns frequency
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
	Tv operator()(){ Tv r=(Tv)0; for(uint32_t j=0; j<this->size(); ++j) r+=(*this)[j](); return r; }
	
	/// Generate next samples adding into a buffer
	template <class V>
	void add(V * dst, uint32_t n){
		for(uint32_t j=0; j<this->size(); ++j){
			for(uint32_t i=0; i<n; ++i)	dst[i] += (*this)(j);
		}
	}

	/// Returns ith oscillator's last value
	Tv last(uint32_t i){ return (*this)[i].val; }

	/// Set all control parameters
	void set(uint32_t i, Tv frq, Tv amp, Tv dcy, Tv phs=0){ (*this)[i].set(frq*Ts::ups(), amp, dcy*Ts::spu(), phs); }

private:
	typedef Array<SineD<Tv, Synced1> > Base;
};



/// Lookup table sine oscillator.

/// This oscillator looks up values in a table containing the sine function
/// in [0, pi/2].  Doing a non-interpolating table lookup is very fast and
/// stable compared to other methods.  The downsides are that the waveform is
/// generally not as spectrally pure and additional memory needs to be allocated
/// to store the lookup table (although it's relatively small and only allocated
/// once).
template <class Stap=tap::Wrap, class Ts=Synced>
class TableSine : public Accum<Stap,Ts> {
public:

	/// @param[in]	frq		Initial frequency
	/// @param[in]	phase	Initial unit phase in [0, 1)
	TableSine(float frq=440, float phase=0);

	float operator()();
	float nextN();				///< Return next non-interpolated sample
	float nextL();				///< Return next linearly-interpolated sample
	
protected:
//	static ArrayPow2<float> mTable; // can't use because need 2**N+1 table
	static float * mTable;		// Reference to my sample table. Must be 1<<bits.
	static uint32_t mTblBits;
	static uint32_t mFracBits;	// # of bits in fractional part of phasor
	static uint32_t mOneIndex;
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
	
	/// @param[in] frq		Initial frequency
	/// @param[in] phase	Initial phase in [0, 1)
	/// @param[in] mod		Initial modifier amount in [0, 1)
	LFO(float frq, float phase=0, float mod=0.5);

	uint32_t modi;			///< Modifier parameter

	void operator()(float f, float p, float m);
	void mod(double n);	///< Sets modifier parameter of waveform from unit value

	float cos();		///< Cosine based on 3rd order polynomial
	float down();		///< Downward ramp (1 to -1)
	float even3();		///< Even harmonic sine-like wave (3rd order)
	float even5();		///< Even harmonic sine-like wave (5th order)
	float imp();		///< Impulse (occurs at beginning of cycle)
	float line2();		///< 2-segment line. mod changes wave from down to tri to up
	float pulse();		///< Pulse (up + down). 'mod' controls pulse width
	float stair();		///< Stair (square + square). 'mod' controls pulse width
	float sqr();		///< Square (-1 to 1)
	float tri();		///< Triangle (starts at 1 going down to -1 then up to 1)
	float up();			///< Upward ramp
	float up2();		///< Dual upward ramp (up + up). 'mod' controls pulse width.

	float cosU();		///< Unipolar cosine based on 3rd order polynomial
	float downU();		///< Unipolar downward ramp
	float hann();		///< Hann window
	float line2U();		///< Unipolar line2
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

private:
	ACCUM_INHERIT
	typedef Accum<Stap,Ts> Base;
};



/// Variable harmonic impulse wave
template<class Tv=gam::real, class Ts=Synced>
class Buzz : public AccumPhase<Tv,Ts> {
public:

	/// @param[in]	frq			Initial frequency
	/// @param[in]	phase		Initial unit phase in [0, 1)
	/// @param[in]	harmonics	Initial number of harmonics
	Buzz(Tv frq=440, Tv phase=0, Tv harmonics=8);
	virtual ~Buzz(){}

	void antialias();			///< Adjust number of harmonics to prevent aliasing
	void harmonics(Tv num);		///< Set number of harmonics
	void harmonicsMax();		///< Set number of harmonics to fill Nyquist range

	Tv operator()();			///< Returns next sample of all harmonic impulse
	Tv odd();					///< Returns next sample of odd harmonic impulse
	Tv saw(Tv intg=0.993);		///< Returns next sample of saw waveform
	Tv square(Tv intg=0.993);	///< Returns next sample of square waveform
	
	Tv maxHarmonics();			///< Returns number of harmonics below Nyquist based on current settings

	virtual void onResync(double r);

protected:
	Tv mAmp;			// amplitude normalization factor
	Tv mN;				// # harmonics
	Tv mNDesired;		// desired number of harmonics
	Tv mNFrac;		
	Tv mSPU_2;			// cached locals
	Tv mPrev;			// previous output for integration
	void recache();
	void setAmp();
private: typedef AccumPhase<Tv,Ts> Base;
};



/// Band-limited impulse wave
template <class Tv=gam::real, class Ts=Synced>
struct Impulse : public Buzz<Tv,Ts>{
public:
	Impulse(Tv frq=440, Tv phase=0): Base(frq, phase){ onResync(1); }

	/// Set frequency
	void freq(Tv v){ Base::freq(v); Base::harmonicsMax(); }

	virtual void onResync(double r){ Base::recache(); freq(AccumPhase<Tv, Ts>::freq()); }

private: typedef Buzz<Tv,Ts> Base; using Base::freq;
};



/// Band-limited saw wave
template <class Tv=gam::real, class Ts=Synced>
struct Saw : public Impulse<Tv,Ts> {
	Saw(Tv frq=440, Tv phase=0): Impulse<Tv, Ts>(frq, phase){}
	Tv operator()(Tv intg=0.993){ return Impulse<Tv,Ts>::saw(intg); }
};



/// Band-limited square wave
template <class Tv=gam::real, class Ts=Synced>
struct Square : public Impulse<Tv,Ts> {
	Square(Tv frq=440, Tv phase=0) : Impulse<Tv,Ts>(frq, phase){}
	Tv operator()(Tv intg=0.993){ return Impulse<Tv,Ts>::square(intg); }
};

// Variable harmonic saw wave.

// This generator integrates a band-limited impulse using a leaky integrater
// to suppress DC build-up.  The integration amount can be adjusted to go
// between a pure impulse to a saw.
//class Saw : public Impulse{
//public:
//	/// @param[in]	frq		Initial frequency in Hz.
//	/// @param[in]	harmonics	Initial number of harmonics.
//	Saw(float frq=440.f, float harmonics=8.f, float integration=0.993);
//	virtual ~Saw(){}
//
//	float next();					///< Returns next sample.
//	
//	/// Set integration amount
//	
//	/// This amount, in [0, 1), determines how much integration of the impulse
//	/// wave occurs and how long it takes for DC to decay.
//	/// Higher values approach a more ideal saw wave, but result
//	/// in more DC offset.	Lower values suppress the DC but pass more of the 
//	/// impulse through resulting in a brighter sound.  A reasonable comprise
//	/// is ~0.993 at a 44.1k sampling rate.
//	void integration(float normal);
//	
//protected:
//	float mInt;
//};



/// Discrete summation formula (DSF) oscillator
template<class Tv=gam::real, class Ts=Synced>
class DSF : public AccumPhase<Tv,Ts> {
public:

	/// @param[in]	frq		Initial frequency in Hz
	/// @param[in]	freqRatio	Initial frequency ratio of partials
	/// @param[in]	ampRatio	Initial amplitude ratio of partials
	/// @param[in]	harmonics	Initial number of harmonics
	DSF(Tv frq=440, Tv freqRatio=1, Tv ampRatio=0.5, Tv harmonics=8);
	virtual ~DSF(){}
	
	Tv operator()();			///< Returns next sample
	
	void ampRatio(Tv val);		///< Set amplitude ratio of partials
	void antialias();			///< Adjust harmonics so partials do not alias
	void freq(Tv val);			///< Set frequency of fundamental
	void freqRatio(Tv val);		///< Set frequency ratio of partials
	void harmonics(Tv val);		///< Set number of harmonics
	void harmonicsMax();		///< Set number of harmonics to fill Nyquist range

	Tv ampRatio();				///< Returns amplitude ratio
	Tv freqRatio();				///< Returns frequency ratio
	Tv harmonics();				///< Returns current number of harmonics
	Tv maxHarmonics();			///< Returns maximum number of harmonics for current settings
	
	virtual void onResync(double r);

protected:
	typedef AccumPhase<Tv,Ts> Base;

	Tv mN, mNDesired;		// actual and desired # harmonics
	Tv mFreqRatio;			// Frequency ratio
	Tv mA;					// Partial amplitude ratio
	Tv mBeta, mBetaInc;		// "detune" phasor
	Tv mAPow, mASqP1;		// cached vars
	
	void updateAPow();
	void updateBetaInc();
};



// Simple band-limited impulse generator.

// This uses a fast, simplified formula for generating a band-limited impulse,
// but only operates at integer divisions of the Nyquist frequency.
class ImpulseFast : public Synced {
public:
	/// Constructor.
	ImpulseFast(): mPhase(0), mOffset(0){ freq(0); }


	/// Set frequency of oscillation.
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

TEMTS TACCUM::Accum(): mFreq(0){
	Ts::initSynced();
	this->phase(0);
}

TEMTS TACCUM::Accum(float freq, float phase): mFreq(freq){
	Ts::initSynced();
	(phase >= 1.f) ? phaseMax() : this->phase(phase);
}

TEMTS inline uint32_t TACCUM::phaseFI(float v) const {
	//return scl::unitToUInt(v);
	//return (uint32_t)(v * 4294967296.);
	return castIntRound(v * 4294967296.);
}

TEMTS inline float TACCUM::phaseIF(uint32_t v) const {
	return uintToUnit<float>(v);
}

TEMTS void TACCUM::onResync(double r){ //printf("Accum: onSyncChange\n");
	freq(mFreq);
}

TEMTS inline void TACCUM::freq(float v){
	mFreq = v;
	mPhaseInc = phaseFI(v * Ts::ups());
}

TEMTS inline void TACCUM::period(float value){ freq(1.f / value); }
TEMTS inline void TACCUM::phase(float v){ mPhase = phaseFI(v); }
TEMTS inline void TACCUM::phaseAdd(float v){ mTap(mPhase, phaseFI(v)); }
TEMTS inline void TACCUM::phaseMax(){ mPhase = 0xffffffff; }

TEMTS inline float TACCUM::freq() const { return mFreq; }
TEMTS inline float TACCUM::phase() const { return phaseIF(phaseI()); }
TEMTS inline uint32_t TACCUM::phaseI() const { return mPhase; }
TEMTS inline float TACCUM::phaseInc() const { return phaseIF(phaseIncI()); }
TEMTS inline uint32_t TACCUM::phaseIncI() const { return mPhaseInc; }
//TEMTS inline uint32_t TACCUM::incPhase(){ return mPhase += phaseIncI(); }
TEMTS inline uint32_t TACCUM::incPhase(){ return mTap(mPhase, phaseIncI()); }

TEMTS inline uint32_t TACCUM::operator()(){ return cycle(); }

TEMTS inline uint32_t TACCUM::cycle(){ return cycles() & 0x80000000; }

//TEMTS inline uint32_t TACCUM::cycle(uint32_t mask){
//	return cycles() & mask;
//}

TEMTS inline uint32_t TACCUM::cycles(){
	uint32_t prev = phaseI();
	incPhase();	
	return ~phaseI() & prev;
}

TEMTS inline uint32_t TACCUM::once(){
	uint32_t prev = phaseI();
	uint32_t c = cycle();
	if(c) mPhase = prev;
	return c;
}

#undef TACCUM



//---- AccumPhase

TEM AccumPhase<Tv, Ts>::AccumPhase(Tv frq, Tv phase)
:	mFreq(frq)
{
	Ts::initSynced();
	this->phase(phase);
}

TEM inline Tv AccumPhase<Tv, Ts>::nextPhase(){
	Tv r = mPhase;
	mPhase = scl::wrapPhase(mPhase + mPhaseInc);
	return r;
}

TEM inline void AccumPhase<Tv, Ts>::freq(Tv v){
	mFreq = v;
	mPhaseInc = v * m2PiUPS;
}

TEM inline void AccumPhase<Tv, Ts>::period(Tv v){ freq(Tv(1)/v); }
TEM inline void AccumPhase<Tv, Ts>::phase(Tv u){ mPhase = u * Tv(M_2PI); }
TEM inline void AccumPhase<Tv, Ts>::phaseAdd(Tv u){ mPhase += u * Tv(M_2PI); }

TEM inline Tv AccumPhase<Tv, Ts>::freq(){ return mFreq; }
TEM inline Tv AccumPhase<Tv, Ts>::period(){ return Tv(1) / mFreq; }
TEM inline Tv AccumPhase<Tv, Ts>::phase(){ return mPhase * Tv(M_1_2PI); }

TEM void AccumPhase<Tv, Ts>::onResync(double r){ recache(); freq(mFreq); }
TEM void AccumPhase<Tv, Ts>::recache(){ m2PiUPS = Tv(Ts::ups() * M_2PI); }

TEM void AccumPhase<Tv, Ts>::print(const char * append, FILE * fp){
	fprintf(fp, "%f %f %f%s", freq(), phase(), mPhaseInc, append);
}


//---- Quadra

TEM Quadra<Tv, Ts>::Quadra(Tv frq, Tv amp, Tv dcy60, Tv phase)
	: val(amp, 0), mAmp(amp), mFreq(frq), mDcy60(dcy60)
{
	Ts::initSynced();
	this->phase(phase);
}

TEM void Quadra<Tv, Ts>::amp(Tv v){
	if(scl::abs(mAmp) > (Tv)0.000001){
		Tv factor = v / mAmp;
		val *= factor;
	} else {
		val(v, (Tv)0);
	}
	mAmp = v;
}

TEM void Quadra<Tv, Ts>::decay(Tv v){
	mDcy60 = v;
	mDcy = v > (Tv)0 ? (Tv)scl::t60(v * Ts::spu()) : (Tv)1;
	freq(mFreq);
}

TEM void Quadra<Tv, Ts>::freq(Tv v){
	mFreq = v;
	Tv phaseInc = v * Ts::ups() * (Tv)M_2PI;
	mInc.fromPolar(mDcy, phaseInc);
	//printf("%f %f %f %f\n", phaseInc, mDcy, c1, s1);
}

TEM void Quadra<Tv, Ts>::phase(Tv v){
	// set phase without changing current magnitude
	val.fromPolar(val.norm(), v*(Tv)M_2PI);
}

TEM void Quadra<Tv, Ts>::reset(){ val(mAmp, (Tv)0); }

TEM void Quadra<Tv, Ts>::set(Tv frq, Tv phase, Tv amp, Tv dcy60){
	mFreq = frq;
	decay(dcy60);
	this->amp(amp);
	this->phase(phase);
}

TEM inline Complex<Tv> Quadra<Tv, Ts>::operator()(){
	complex c = val;
	val *= mInc;
	return c;
}

TEM inline Tv Quadra<Tv, Ts>::freq(){ return mFreq; }

TEM void Quadra<Tv, Ts>::onResync(double r){
	decay(mDcy60); // this sets frequency as well
}



//---- TableSine

#define TTABLESINE TableSine<St,Ts>
TEMTS uint32_t TTABLESINE::mTblBits  = 9UL;	// actual table memory is a quarter of this
TEMTS uint32_t TTABLESINE::mFracBits = 32UL - TableSine::mTblBits;
TEMTS uint32_t TTABLESINE::mOneIndex = 0;
TEMTS float * TTABLESINE::mTable = 0;

TEMTS TTABLESINE::TableSine(float freq, float phase): Base(freq, phase){
	if(0 == mTable){
		mOneIndex = 1<<mFracBits;
		uint32_t size = 1<<(mTblBits-2);
		mTable = new float[size + 1];	// does this need to be deleted?
		tbl::sinusoid(mTable, size, 0, 0.25);
		mTable[size] = 1;
	}
}

TEMTS inline float TTABLESINE::operator()(){ return nextL(); }

TEMTS inline float TTABLESINE::nextN(){
	float output = tbl::atQ(mTable, mFracBits, phaseI());
	incPhase();
	return output;
}

TEMTS inline float TTABLESINE::nextL(){
	float output = ipl::linear(
		gam::fraction(mTblBits, phaseI()),
		tbl::atQ(mTable, mFracBits, phaseI()),
		tbl::atQ(mTable, mFracBits, phaseI() + mOneIndex)
	);
	incPhase();
	return output;
}

#undef TTABLESINE



//---- LFO
#define TLFO LFO<St,Ts>
TEMTS TLFO::LFO(): Base(){ mod(0.5); }
TEMTS TLFO::LFO(float f, float p, float m): Base(f, p){ mod(m); }

TEMTS inline void TLFO::operator()(float f, float p, float m){ this->freq(f); this->phase(p); mod(m); }
TEMTS inline void TLFO::mod(double n){ modi = castIntRound(n * 4294967296.); }

#define DEF(name, exp) TEMTS inline float TLFO::name{ float r = exp; incPhase(); return r; }


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
	incPhase();
	return r;
}

TEMTS inline float TLFO::line2U(){
	return line2()*0.5f+0.5f;
}

DEF(down(),		scl::rampDown(phaseI()))
DEF(even3(),	scl::rampUp(phaseI()); static float c=-1.50*sqrt(3.); r *= (1-r*r)*c;)
DEF(even5(),	scl::rampUp(phaseI()); static float c=-1.25*::pow(5.,0.25); r *= (1-scl::pow4(r))*c;)
DEF(pulse(),	scl::pulse(phaseI(), modi))
DEF(cos(),		scl::triangle(phaseI()); r *= 0.5f * r*r - 1.5f)
DEF(stair(),	scl::stair(phaseI(), modi))
DEF(sqr(),		scl::square(phaseI()))
DEF(tri(),		scl::triangle(phaseI()))
DEF(up(),		scl::rampUp(phaseI()))
DEF(up2(),		scl::rampUp2(phaseI(), modi))
DEF(cosU(),		scl::triangle(phaseI()); r = r * (0.25f * r*r - 0.75f) + 0.5f)
DEF(downU(),	scl::rampDownU(phaseI()))
DEF(hann(),		scl::triangle(phaseI()); r = scl::warpSinSU(r))
DEF(pulseU(),	scl::pulseU(phaseI(), modi))
DEF(sqrU(),		scl::squareU(phaseI()))
DEF(stairU(),	scl::stairU(phaseI(), modi))
DEF(triU(),		scl::triangleU(phaseI()))
DEF(upU(),		scl::rampUpU(phaseI()))
DEF(up2U(),		scl::rampUp2U(phaseI(), modi))
DEF(patU(),		scl::rampUpU(phaseI() & modi))

DEF(patU(uint32_t mul), scl::rampUpU((phaseI() & modi) * mul))

DEF(sineT9(),	scl::rampUp(phaseI()); r = scl::sinT9(r * M_PI))
DEF(sineP9(),	scl::rampUp(phaseI()); r = scl::sinP9(r))

#undef DEF

//TEMS inline float TLFO::imp(){ return this->cycle() ? 1.f : 0.f; }
TEMTS inline float TLFO::imp(){ 
	float r = phaseI() < this->phaseIncI() ? 1.f : 0.f;
	incPhase();
	return r; 
}

TEMTS inline bool TLFO::seq(){
	uint32_t prev = phaseI();
	incPhase();
	if( (phaseI() ^ prev) & 0xf8000000 ){
		if( (modi >> (phaseI()>>27)) & 0x1 ) return true;
	}
	return false;
}

#undef TLFO


//---- Buzz

TEM Buzz<Tv,Ts>::Buzz(Tv frq, Tv phase, Tv harmonics)
:	Base(frq, phase), mAmp(0), mPrev((Tv)0)
{
	onResync(1);
	this->harmonics(harmonics);
}

TEM inline void Buzz<Tv,Ts>::harmonics(Tv value){
	mN = mNDesired = scl::floor(value);
	setAmp();
	mNFrac = value - mN;
}

TEM inline void Buzz<Tv,Ts>::harmonicsMax(){ harmonics(maxHarmonics()); }

TEM inline void Buzz<Tv,Ts>::antialias(){
	float maxN = scl::floor(maxHarmonics());
	mN = mNDesired > maxN ? maxN : mNDesired;
	setAmp();
}

TEM inline Tv Buzz<Tv,Ts>::maxHarmonics(){ return mSPU_2 / this->mFreq; }

TEM inline void Buzz<Tv,Ts>::setAmp(){ mAmp = (mN != Tv(0)) ? (Tv(0.5) / mN) : 0; }

#define EPS 0.000001
TEM inline Tv Buzz<Tv, Ts>::operator()(){
	Tv theta = this->nextPhase();

	Tv result;
	Tv denom = scl::sinT9(theta * (Tv)0.5);
	if( scl::abs(denom) < (Tv)EPS ){ result = Tv(1); /*printf("Impulse::(): oops\n");*/ }
	else{
		Tv nphase = scl::wrapPhase(theta * (mN + Tv(0.5)));
		//result = (scl::sinT7(nphase) / denom - Tv(1)) * Tv(0.5);
		
		//result = ((scl::sinT7(nphase) - denom) / (denom*mN)) * Tv(0.5);
		result = ((scl::sinT7(nphase) - denom) / denom) * mAmp;
		
//		(n/d - 1)*0.5/N
//		n/d - d/d
//		(n-d)/d*0.5/N
	}
	
	//result -= 2.f * cos(theta * mN) * (1.f - mNFrac);	// removes clicks!
	
	return result;
}

TEM inline Tv Buzz<Tv,Ts>::odd(){
	Tv theta = this->nextPhase();
	
	Tv result;
	Tv n = scl::ceil(mN, Tv(2), Tv(0.5));	// get next highest even for formula
										// wave has odd harmonics 1,3, ..., n-1 with peak amp n

	Tv denom = scl::cosT8(scl::wrapPhaseOnce(theta - M_PI_2));	// this is more precise near zero-crossings
	//Tv denom = scl::sinT9(theta);
	if( scl::abs(denom) < Tv(EPS) ){
		if( theta > M_PI ) theta -= M_2PI;
		result = (theta > -M_PI_2 && theta < M_PI_2) ? 1 : -1;
		//printf("Impulse::odd(): oops\n");
	}
	else result = scl::sinT7(scl::wrapPhase(n * theta)) / (denom * n);
	
	return result;
}
#undef EPS

TEM inline Tv Buzz<Tv,Ts>::saw(Tv i){ return mPrev=(*this)()*0.125 + i*mPrev; }
TEM inline Tv Buzz<Tv,Ts>::square(Tv i){ return mPrev=odd()*0.125 + i*mPrev; }

TEM void Buzz<Tv,Ts>::onResync(double r){
	Base::recache();
	Base::freq(Base::freq());
}

TEM void Buzz<Tv,Ts>::recache(){
	Base::recache();
	mSPU_2 = (Tv)(Synced::spu() * 0.5);
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

TEM inline void DSF<Tv,Ts>::freq(Tv value){
	Base::freq(value);
	updateBetaInc();	
}

TEM inline void DSF<Tv,Ts>::freqRatio(Tv value){
	mFreqRatio = value;
	updateBetaInc();
}

TEM inline void DSF<Tv,Ts>::ampRatio(Tv value){
	if(value != mA){
		mA = value;
		mASqP1 = mA * mA + 1.f;
		updateAPow();
	}
}

TEM inline void DSF<Tv,Ts>::harmonics(Tv value){
	if(value != mN){
		mN = mNDesired = value;
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
	return scl::floor((Tv(this->spu()) * Tv(0.5)/this->mFreq - Tv(1))/mFreqRatio + Tv(1));
}

// Generalized DSF formula:
// sum{k=0, N}( a^k sin(theta + k beta) )
//		=  sin(theta) - a sin(theta - beta) - a^(N+1) (sin(theta + (N+1) beta) - a sin(theta + N beta)))
//			/ 1 + a^2 - 2a cos(beta)
//
// The denominator is the frequency response of a resonator w/ resonance a.

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
	result /= mASqP1 - (Tv)2 * mA * COS(mBeta);

	mBeta += mBetaInc;
	
	return result;
}
#undef SIN
#undef COS

TEM inline void DSF<Tv,Ts>::updateAPow(){ mAPow = ::pow(mA, mN); }
TEM inline void DSF<Tv,Ts>::updateBetaInc(){ mBetaInc = this->mPhaseInc * mFreqRatio; }

TEM void DSF<Tv,Ts>::onResync(double r){
	Base::onResync(r);
	freq(Base::mFreq);
	harmonics(mNDesired);
}


// This object stores a accumulator increment value.
template <class Tv=gam::real, class Ts=Synced>
class Increment : public Ts{
public:
	Increment(double frq): mInc(0){ Ts::initSynced(); freq(frq); }
	
	void freq(double v){ mInc = (Tv)(v * Ts::ups()); }
	Tv& operator()(Tv& v){ return v = scl::wrap(v + mInc); }
	virtual void onResync(double r){ mInc/=r; }
	
protected:
	Tv mInc;
};



#undef TEM
#undef TEMS
#undef TEMTS

} // gam::
#endif
