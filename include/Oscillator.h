#ifndef GAMMA_OSCILLATOR_H_INC
#define GAMMA_OSCILLATOR_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include "arr.h"
#include "gen.h"
#include "scl.h"
#include "tbl.h"
#include "Strategy.h"
#include "Sync.h"
#include "Types.h"

namespace gam{


/// Fixed-point phase accumulator.
template <class Ts=Synced>
class Accum : public Ts {
public:
	Accum();

	/// @param[in] frq		Initial frequency.
	/// @param[in] phs		Initial phase [0, 1).
	Accum(float frq, float phs=0);
	virtual ~Accum(){}

	void freq(float value);			///< Set frequency.
	void phase(float normal);		///< Set phase from [0, 1) of one period.
	void phaseMax();				///< Set phase to maximum value.
	void phaseAdd(float normal);	///< Add value to phase [0, 1).		
	void period(float value);		///< Set period length.

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
	
	uint32_t phaseFI(float v) const;	// convert unit floating-point to fixed-point integer
	float phaseIF(uint32_t v) const;	// convert fixed-point integer to unit floating-point
};

#define ACCUM_INHERIT\
	using Accum<Ts>::phaseI;\
	using Accum<Ts>::phaseIncI;\
	using Accum<Ts>::incPhase;



/// Floating point phase accumulator in [-pi, pi).
template <class Tv=gam::real, class Ts=Synced>
class AccumPhase : public Ts{
public:
	/// @param[in]	frq		Initial frequency.
	/// @param[in]	phs		Initial phase [0, 1).
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



/// General-purpose table oscillator

/// Tv is the table element type, Sipol is an interpolation strategy, and
/// Stap is a table reading strategy.
template <class Tv=gam::real, class Sipol=ipl::Linear, class Stap=tap::Wrap, class Ts=Synced>
class Osc : public Accum<Ts>, public ArrayPow2<Tv>{
public:

	ACCUM_INHERIT
	using ArrayPow2<Tv>::elems; using ArrayPow2<Tv>::size;

	/// Constructor that alocates an internal table

	/// @param[in]	frq			Initial frequency
	/// @param[in]	phs			Initial normalized phase [0, 1)
	/// @param[in]	size		Number of table elements (actual number is power of 2 ceiling)
	Osc(float frq, float phs=0, uint32_t size=512)
	:	Accum<Ts>(frq, phs), ArrayPow2<Tv>(size)
	{}

	/// Constructor that references an external table

	/// @param[in]	frq			Initial frequency
	/// @param[in]	phs			Initial normalized phase [0, 1)
	/// @param[in]	src			A table to use as a reference
	Osc(float frq, float phs, const ArrayPow2<Tv> & src)
	:	Accum<Ts>(frq, phs), ArrayPow2<Tv>(src.elems(), src.size())
	{}
	
	virtual ~Osc(){}

	/// Generate next sample
	Tv operator()(){
		Tv o0 = val(); mTap(this->mPhase, phaseIncI()); return o0;
	}
	
	void zero(){ for(int i=0; i<this->size(); ++i) (*this)[i] = 0; }

	Tv val() const { return mIpol(*this, phaseI()); }
	
	bool done(){ return tap().done(phaseI()); }
	
	Stap& tap(){ return mTap; }
	
	// Add [harmonic, amp, phs] triple to transfer function
	template <class T1, class T2, class T3>
	Osc& operator<<(const Tup3<T1, T2, T3> & t){
		Tv * e = elems();
		arr::add(e, gen::Sin<Tv>(t.v1 * M_2PI / (Tv)size(), t.v3 * M_2PI, t.v2), size());
		return *this;
	}

protected:
	Sipol mIpol;
	Stap mTap;
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

	/// @param[in] frq		Initial frequency.
	/// @param[in] amp		Initial amplitude.
	/// @param[in] dcy		Initial decay time. Negative means no decay.
	/// @param[in] phs		Initial phase in [0, 1).
	Quadra(Tv frq=440, Tv amp=1, Tv dcy=-1, Tv phs=0);
	//virtual ~Quadra(){}

	complex val;			///< Current complex output.
	
	void amp(Tv val);		///< Set amplitude of sinusoids.
	void decay(Tv val);		///< Set number of units to decay -60 dB. Negative = no decay.
	void freq(Tv val);		///< Set frequency of sinusoids.
	void phase(Tv radians);	///< Set phase of sinusoids.
	void reset();			///< Resets amplitude and sets phase to 0.
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
	/// @param[in]	frq		Initial frequency.
	/// @param[in]	phs		Initial normalized phase [0, 1).
	Sine(Tv frq=440, Tv phs=0) : AccumPhase<Tv, Ts>(frq, phs){}

	/// Return next sample.
	Tv operator()(){ return scl::sinP9(AccumPhase<Tv, Ts>::nextPhase() * M_1_PI); }
};



/// Sine oscillator based on an efficient recursion equation.

/// This oscillator is based on a recursion equation requiring only one
/// multiply and add per sample computation.  The downsides are that frequency
/// and phase updates are extremely expensive and 64-bit precision is required 
/// to prevent growing or decaying in amplitude over time.  This generator is 
/// ideal in situations where a stationary sinusoid is all that is required, 
/// e.g. a grain or modulator.
template <class Tv=double, class Ts=Synced>
class SineR : public gen::RSin<Tv>, Ts{
public:
	typedef gen::RSin<Tv> super;

	SineR(Tv frq=440, Tv amp=1, Tv phs=0){ Ts::initSynced(); set(frq, amp, phs); }

	/// Returns frequency
	Tv freq() const { return super::freq() * Ts::spu(); }

	/// Set amplitude and phase
	void ampPhase(Tv a=1, Tv p=0){ set(freq(), a, p); }

	/// Set all control parameters
	void set(Tv frq, Tv amp, Tv phs=0){ gen::RSin<Tv>::set(frq*Ts::ups(), phs, amp); }
	
	/// Generate next samples adding into a buffer
	template <class V>
	void add(V * dst, uint32_t n){ for(uint32_t i=0; i<n; ++i) dst[i] += (*this)(); }

	// This might not be a good idea because we don't know amplitude...
	//virtual void onResync(double ratio){ set(gen::RSin<Tv>::freq()/ratio, 1, 0); }
};



/// Multiple SineRs

/// For efficiency reasons, this object does not keep its frequencies synchronized
/// with the sample rate. If the sample rate changes, each oscillator must be
/// manually re-set.
template <class Tv=double, class Ts=Synced>
class SineRs : public Array<SineR<Tv, Synced1> >, Ts{
public:

	typedef Array<SineR<Tv, Synced1> > super;

	/// @param[in]	num		Number of resonators
	SineRs(uint32_t num): super(num){ Ts::initSynced(); }

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
};



/// Damped sine oscillator based on an efficient recursion equation.

/// This oscillator is similar to SineR, however, it has an extra multiply
/// in its sample generation to allow the oscillator to decay.
template <class Tv=double, class Ts=Synced>
class SineD : public gen::RSin2<Tv>, Ts{
public:
	typedef gen::RSin2<Tv> super;

	/// @param[in]	frq		Initial frequency.
	/// @param[in]	amp		Initial amplitude.
	/// @param[in]	dcy		Initial T60 decay length.
	/// @param[in]	phs		Initial unit phase [0, 1).
	SineD(Tv frq=440, Tv amp=1, Tv dcy=-1, Tv phs=0){ Ts::initSynced(); set(frq, amp, dcy, phs); }

	/// Returns frequency
	Tv freq() const { return super::freq() * Ts::spu(); }

	/// Set amplitude and phase
	void ampPhase(Tv a=1, Tv p=0){ set(freq(), a, super::decay(), p); }
	
	/// Set all control parameters
	void set(Tv frq, Tv amp, Tv dcy, Tv phs=0){
		super::set(frq*Ts::ups(), phs, dcy > (Tv)0 ? (Tv)scl::radius60(dcy, Ts::ups()) : (Tv)1, amp);
	}
	
	/// Generate next samples adding into a buffer
	template <class V>
	void add(V * dst, uint32_t n){ for(uint32_t i=0; i<n; ++i) dst[i] += (*this)(); }
};



/// Multiple SineDs

/// For efficiency reasons, this object does not keep its frequencies synchronized
/// with the sample rate. If the sample rate changes, each oscillator must be
/// manually re-set.
template <class Tv=double, class Ts=Synced>
class SineDs : public Array<SineD<Tv, Synced1> >, Ts{
public:

	typedef Array<SineD<Tv, Synced1> > super;

	/// @param[in]	num		Number of resonators
	SineDs(uint32_t num): super(num){
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
};





/// Oscillator that reads values from a power-of-two table. Deprecated: use Osc.
template<class Tv=gam::real, class Ts=Synced>
class TableOsc : public Accum<Ts>{
public:
	ACCUM_INHERIT

	/// @param[in]	table		Reference to array the oscillator will read from
	/// @param[in]	log2Size	Number of table elements. The number of elements is 2^log2Size.
	/// @param[in]	frq			Initial frequency
	/// @param[in]	phase		Initial normalized phase [0, 1)
	TableOsc(Tv * table, uint32_t size, float frq=440, float phase=0);
	virtual ~TableOsc(){}
	
	void table(Tv * src, uint32_t size);	///< Set my table reference.
	
	Tv nextN();				///< Return next non-interpolated sample.
	Tv nextL();				///< Return next linearly-interpolated sample.
	Tv onceL();
	
	//virtual void onResync(double r);
	
protected:
	Tv * mTable;			// Reference to my sample table. Must be 1<<bits.
	uint32_t mTblBits;
	uint32_t mFracBits;	// # of bits in fractional part of phasor
	uint32_t mOneIndex;
};



/// Lookup table sine oscillator.

/// This oscillator looks up values in a table containing the sine function
/// in [0, pi/2].  Doing a non-interpolating table lookup is very fast and
/// stable compared to other methods.  The downsides are that the waveform is
/// generally not as spectrally pure and additional memory needs to be allocated
/// to store the lookup table (although it's relatively small and only allocated
/// once).
template <class Ts=Synced>
class TableSine : public Accum<Ts> {
public:
	ACCUM_INHERIT

	/// @param[in]	frq		Initial frequency
	/// @param[in]	phase	Initial normalized phase [0, 1)
	TableSine(float frq=440, float phase=0);
	//virtual ~TableSine(){}

	float operator()();
	float nextN();				///< Return next non-interpolated sample
	float nextL();				///< Return next linearly-interpolated sample
	
protected:
//	static ArrayPow2<float> mTable; // can't use because need 2**N+1 table
	static float * mTable;		// Reference to my sample table. Must be 1<<bits.
	static uint32_t mTblBits;
	static uint32_t mFracBits;	// # of bits in fractional part of phasor
	static uint32_t mOneIndex;
};



/// Multiple table oscillator for making band-limited waveforms.

/// Each table is mapped to an octave in the frequency spectrum.  The highest
/// table maps to normalized frequencies [1/4, 1/2) and the lowest to [0, 1/m) 
/// where m = 2 ^ (numTables + 1).
template<class Tv=gam::real, class Ts=Synced>
class MultiTableOsc : public TableOsc<Tv, Ts> {
public:

	/// @param[in]	table		Reference to matrix the oscillator will read from
	/// @param[in]	log2Size	Log base-2 of the number of elements in each waveform table
	/// @param[in]	numTables	Number of waveform tables
	/// @param[in]	frq		Initial frequency
	/// @param[in]	phase		Initial normalized phase [0, 1)
	MultiTableOsc(Tv * table, uint32_t log2Size, uint32_t numTables, float frq=220, float phase=0);
	virtual ~MultiTableOsc(){}
	
	void freq(float value);		///< Set frequency of oscillation.
	void freqLL(float value);	///< Set frequency of oscillation when using nextLL().

	/// Return next linearly-interpolated sample from linearly-interpolated table.
	
	/// Use freqLL() to update frequency.
	///
	Tv nextLL();				
	
	void table(Tv * src, uint32_t log2Size, uint32_t numTables);	///< Set my table reference.
	
	float freq(){ return this->mFreq; }
	
	virtual void onResync(double r);
	
protected:
	uint32_t mNumTables;
	Tv * mMultiTable;
	Tv * mTableLo;
	float mFracTbl;
};



/// Linear function oscillator (non band-limited).

/// This object generates various waveform types by modifying
/// the output of a single upward ramping accumulator.
template <class Ts=Synced>
class LFO : public Accum<Ts>{
public:
	ACCUM_INHERIT

	LFO();
	
	/// @param[in] frq		Initial frequency
	/// @param[in] phase	Initial normalized phase [0, 1).
	/// @param[in] mod		Initial normalized modifier amount
	LFO(float frq, float phase=0, float mod=0.5);

	uint32_t modi;			///< Modifier parameter

	void operator()(float f, float p, float m);
	void mod(double n);	///< Sets modifier parameter of waveform from normal

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

};



/// Variable harmonic impulse using the Dirichlet kernel.
template<class Tv=gam::real, class Ts=Synced>
class Impulse : public AccumPhase<Tv, Ts> {
public:

	/// @param[in]	frq		Initial frequency
	/// @param[in]	phase		Initial normalized phase [0, 1)
	/// @param[in]	harmonics	Initial number of harmonics
	Impulse(Tv frq=440, Tv phase=0, Tv harmonics=8);
	virtual ~Impulse(){}

	void antialias();			///< Adjust number of harmonics to prevent aliasing
	void harmonics(Tv num);		///< Set number of harmonics
	void harmonicsMax();		///< Set number of harmonics to fill Nyquist range

	Tv operator()();				///< Returns next sample of even harmonic impulse
	Tv odd();					///< Returns next sample of odd harmonic impulse
	Tv saw(Tv intg = 0.993);		///< Returns next sample of saw waveform
	Tv square(Tv intg = 0.993);	///< Returns next sample of square waveform
	
	Tv maxHarmonics();			///< Returns number of harmonics below Nyquist based on current settings

	virtual void onResync(double r);

protected:
	Tv mN;			// # harmonics
	Tv mNDesired;	// desired number of harmonics
	Tv mNFrac;		
	Tv mSPU_2;			// cached locals
	Tv mPrev;			// previous output for integration
	void recache();
};



// Super-class for band-limited oscillators
template <class Tv, class Ts=Synced>
struct OscBL : public Impulse<Tv, Ts>{

	typedef Impulse<Tv, Ts> super; using super::freq;

	OscBL(Tv frq=440, Tv phase=0): Impulse<Tv, Ts>(frq, phase){ onResync(1); }

	void freq(Tv v){ super::freq(v); super::harmonicsMax(); }

	virtual void onResync(double r){ super::recache(); freq(AccumPhase<Tv, Ts>::freq()); }
};

template <class Tv=gam::real, class Ts=Synced>
struct Saw : public OscBL<Tv, Ts> {
	Saw(Tv frq=440, Tv phase=0): OscBL<Tv, Ts>(frq, phase){}
	Tv operator()(Tv intg = 0.993){ return Impulse<Tv, Ts>::saw(intg); }
};

template <class Tv=gam::real, class Ts=Synced>
struct Square : public OscBL<Tv, Ts> {
	Square(Tv frq=440, Tv phase=0) : OscBL<Tv, Ts>(frq, phase){}
	Tv operator()(Tv intg = 0.993){ return Impulse<Tv, Ts>::square(intg); }
};



/// Variable harmonic saw wave.

/// This generator integrates a band-limited impulse using a leaky integrater
/// to suppress DC build-up.  The integration amount can be adjusted to go
/// between a pure impulse to a saw.
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



/// Discrete summation formula (DSF) oscillator.
template<class Tv=gam::real, class Ts=Synced>
class DSF : public AccumPhase<Tv, Ts> {
public:

	/// @param[in]	frq		Initial frequency in Hz
	/// @param[in]	freqRatio	Initial frequency ratio of partials
	/// @param[in]	ampRatio	Initial amplitude ratio of partials
	/// @param[in]	harmonics	Initial number of harmonics
	DSF(Tv frq=440, Tv freqRatio=1, Tv ampRatio=0.5, Tv harmonics=8);
	virtual ~DSF(){}
	
	Tv operator()();				///< Returns next sample
	
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
class ImpulseLP : public Synced {
public:
	/// Constructor.
	ImpulseLP(): mPhase(0), mOffset(0){ freq(0); }


	/// Set frequency of oscillation.
	void freq(float v){
		float samples = spu() / v;
		
		uint32_t period = (uint32_t)(samples);
		//period &= 0xfffffffe;		// force period to be even
									// odd periods introduce DC

		mPeriod = (float)period;
		mPeriod = floor(samples);
		
		if(scl::even(period))	mOffset = 0.f;
		//else					mOffset = 0.5f / mPeriod;
	}


	/// Generate next sample
	float operator()(){
		float result = 0.f;

		if(mPhase >= mPeriod){
			mPhase -= mPeriod;
			result = 1.f;
		}
		else if(scl::even((uint32_t)mPhase)){
			result = -1.f/(mPeriod * 0.5f - 1.f);
		}
		
		mPhase++;
		return result + mOffset;
	}
	
protected:
	float mPhase;		// phase in samples
	float mPeriod;		// period in samples;
	float mOffset;		// DC compensation
};







// Implementation_______________________________________________________________

#define TEM template <class Tv, class Ts>
#define TEMS template <class Ts>

//---- Accum

TEMS Accum<Ts>::Accum(): mFreq(0){
	Ts::initSynced();
	this->phase(0);
}

TEMS Accum<Ts>::Accum(float freq, float phase): mFreq(freq){
	Ts::initSynced();
	(phase >= 1.f) ? phaseMax() : this->phase(phase);
}

TEMS inline uint32_t Accum<Ts>::phaseFI(float v) const {
	//return scl::normalToUInt(v);
	//return (uint32_t)(v * 4294967296.);
	return scl::castIntRound(v * 4294967296.);
}

TEMS inline float Accum<Ts>::phaseIF(uint32_t v) const {
	return scl::uintToNormal<float>(v);
}

TEMS void Accum<Ts>::onResync(double r){ //printf("Accum: onSyncChange\n");
	freq(mFreq);
}

TEMS inline void Accum<Ts>::freq(float v){
	mFreq = v;
	mPhaseInc = phaseFI(v * Ts::ups());
}

TEMS inline void Accum<Ts>::period(float value){ freq(1.f / value); }
TEMS inline void Accum<Ts>::phase(float v){ mPhase = phaseFI(v); }
TEMS inline void Accum<Ts>::phaseAdd(float v){ mPhase += phaseFI(v); }
TEMS inline void Accum<Ts>::phaseMax(){ mPhase = 0xffffffff; }

TEMS inline float Accum<Ts>::freq() const { return mFreq; }
TEMS inline float Accum<Ts>::phase() const { return phaseIF(phaseI()); }
TEMS inline uint32_t Accum<Ts>::phaseI() const { return mPhase; }
TEMS inline float Accum<Ts>::phaseInc() const { return phaseIF(phaseIncI()); }
TEMS inline uint32_t Accum<Ts>::phaseIncI() const { return mPhaseInc; }
TEMS inline uint32_t Accum<Ts>::incPhase(){ return mPhase += phaseIncI(); }

TEMS inline uint32_t Accum<Ts>::operator()(){ return cycle(); }

TEMS inline uint32_t Accum<Ts>::cycle(){ return cycles() & 0x80000000; }

//inline uint32_t Accum::cycle(uint32_t mask){
//	return cycles() & mask;
//}

TEMS inline uint32_t Accum<Ts>::cycles(){
	uint32_t prev = phaseI();
	incPhase();	
	return ~phaseI() & prev;
}

TEMS inline uint32_t Accum<Ts>::once(){
	uint32_t prev = phaseI();
	uint32_t c = cycle();
	if(c) mPhase = prev;
	return c;
}



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

TEM inline void AccumPhase<Tv, Ts>::freq(Tv value){
	mFreq = value;
	mPhaseInc = value * m2PiUPS;
}

TEM inline void AccumPhase<Tv, Ts>::period(Tv value){ freq((Tv)1 / value); }
TEM inline void AccumPhase<Tv, Ts>::phase(Tv normal){ mPhase = normal * (Tv)M_2PI; }
TEM inline void AccumPhase<Tv, Ts>::phaseAdd(Tv normal){ mPhase += normal * (Tv)M_2PI; }

TEM inline Tv AccumPhase<Tv, Ts>::freq(){ return mFreq; }
TEM inline Tv AccumPhase<Tv, Ts>::period(){ return (Tv)1 / mFreq; }
TEM inline Tv AccumPhase<Tv, Ts>::phase(){ return mPhase * (Tv)M_1_2PI; }

TEM void AccumPhase<Tv, Ts>::onResync(double r){ recache(); freq(mFreq); }
TEM void AccumPhase<Tv, Ts>::recache(){ m2PiUPS = (Tv)(Ts::ups() * M_2PI); }

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
	val.fromPolar(sqrt(val.dot()), v*(Tv)M_2PI);
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



//---- TableOsc

TEM TableOsc<Tv, Ts>::TableOsc(Tv * tableA, uint32_t log2Size, float frq, float phase)
	: Accum<Ts>(frq, phase)
{
	table(tableA, log2Size);
}

TEM void TableOsc<Tv, Ts>::table(Tv * src, uint32_t log2Size){
	mTable = src;
	mTblBits = log2Size;
	mFracBits = 32U - log2Size;		// for full-wave lookup
	mOneIndex = 1<<mFracBits;
}

TEM inline Tv TableOsc<Tv, Ts>::nextN(){
	Tv output = mem::at(mTable, mFracBits, phaseI());
	incPhase();
	return output;
}

TEM inline Tv TableOsc<Tv, Ts>::nextL(){
	Tv output = ipl::linear(
		tbl::fraction(mTblBits, phaseI()),
		mem::at(mTable, mFracBits, phaseI()),
		mem::at(mTable, mFracBits, phaseI() + mOneIndex)
	);
	incPhase();
	return output;
}

TEM inline Tv TableOsc<Tv, Ts>::onceL(){
	Tv output = ipl::linear(
		tbl::fraction(mTblBits, phaseI()),
		mem::at(mTable, mFracBits, phaseI()),
		mem::at(mTable, mFracBits, phaseI() + mOneIndex)
	);
	
	if(this->cycle()) this->phaseMax();
	
	return output;
}




//---- TableSine

TEMS uint32_t TableSine<Ts>::mTblBits  = 9UL;	// actual table memory is a quarter of this
TEMS uint32_t TableSine<Ts>::mFracBits = 32UL - TableSine::mTblBits;
TEMS uint32_t TableSine<Ts>::mOneIndex = 0;
TEMS float * TableSine<Ts>::mTable = 0;

TEMS TableSine<Ts>::TableSine(float freq, float phase) : Accum<Ts>(freq, phase){
	if(0 == mTable){
		mOneIndex = 1<<mFracBits;
		uint32_t size = 1<<(mTblBits-2);
		mTable = new float[size + 1];	// does this need to be deleted?
		tbl::sinusoid(mTable, size, 0, 0.25);
		mTable[size] = 1;
	}
}

TEMS inline float TableSine<Ts>::operator()(){ return nextL(); }

TEMS inline float TableSine<Ts>::nextN(){
	float output = tbl::atQ(mTable, mFracBits, phaseI());
	incPhase();
	return output;
}

TEMS inline float TableSine<Ts>::nextL(){
	float output = ipl::linear(
		tbl::fraction(mTblBits, phaseI()),
		tbl::atQ(mTable, mFracBits, phaseI()),
		tbl::atQ(mTable, mFracBits, phaseI() + mOneIndex)
	);
	incPhase();
	return output;
}



//---- MultiTableOsc

TEM MultiTableOsc<Tv, Ts>::MultiTableOsc(Tv * table, uint32_t log2Size, uint32_t numTables, float frq, float phase)
	: TableOsc<Tv, Ts>(table, log2Size, frq, phase)
{	
	this->table(table, log2Size, numTables);	
	mFracTbl = 0.f;
	onResync(1);
}

TEM void MultiTableOsc<Tv, Ts>::table(Tv * src, uint32_t log2Size, uint32_t numTables){
	TableOsc<Tv, Ts>::table(src, log2Size);
	mMultiTable = src;
	mTableLo = src;
	mNumTables = numTables;
}

#define MTO MultiTableOsc<Tv, Ts>
TEM void MultiTableOsc<Tv, Ts>::onResync(double r){
	//MTO::mUPS = MTO::ups();
	freqLL(MTO::mFreq);
}

TEM inline void MultiTableOsc<Tv, Ts>::freq(float value){
	MTO::mFreq = value;
	float ratio = value * MTO::mUPS;	// 0 to 1 of sample rate
	if(ratio >= 0.5f || ratio <= -0.5f) ratio = 0.49999998f;
	MTO::mPhaseInc = scl::normalToUInt(ratio);
	uint32_t tableNum = scl::floatExponent(ratio) - 126 + mNumTables;
	if(tableNum >= mNumTables) tableNum = 0;
	MTO::mTable = mMultiTable + tableNum * (1<<MTO::mTblBits);
}

// kinda getting big to inline???
TEM inline void MultiTableOsc<Tv, Ts>::freqLL(float value){
	MTO::mFreq = value;
	float ratio = value * MTO::ups();	// 0 to 1 of sample rate
	if(ratio >= 0.5f || ratio <= -0.5f) ratio = 0.49999998f;
	MTO::mPhaseInc = scl::normalToUInt(ratio);
	uint32_t tableNum = scl::floatExponent(ratio) - 126 + mNumTables;
	if(tableNum >= mNumTables){
		MTO::mTable = mMultiTable;
		mTableLo = MTO::mTable;
		mFracTbl = 0.f;
	}
	else{
		MTO::mTable = mMultiTable + tableNum * (1<<MTO::mTblBits);
		
		tableNum++;
		if(tableNum >= mNumTables) tableNum = mNumTables - 1;
		mTableLo = mMultiTable + tableNum * (1<<MTO::mTblBits);;
		mFracTbl = scl::floatMantissa(ratio);
		mFracTbl *= mFracTbl;	// warp fraction for faster fade-in
		mFracTbl *= mFracTbl;
	}
}

TEM inline Tv MultiTableOsc<Tv, Ts>::nextLL(){
	float fracSmp = tbl::fraction(MTO::mTblBits, MTO::mPhase);
	Tv output = ipl::linear(
		fracSmp,
		mem::at(MTO::mTable, MTO::mFracBits, MTO::mPhase),
		mem::at(MTO::mTable, MTO::mFracBits, MTO::mPhase + MTO::mOneIndex)
	);
	Tv outputLo = ipl::linear(
		fracSmp,
		mem::at(mTableLo, MTO::mFracBits, MTO::mPhase),
		mem::at(mTableLo, MTO::mFracBits, MTO::mPhase + MTO::mOneIndex)
	);
	output = ipl::linear(mFracTbl, output, outputLo);
	MTO::incPhase();
	return output;
}
#undef MTO



//---- LFO

TEMS LFO<Ts>::LFO(): Accum<Ts>(){ mod(0.5); }
TEMS LFO<Ts>::LFO(float f, float p, float m): Accum<Ts>(f, p){ mod(m); }

TEMS inline void LFO<Ts>::operator()(float f, float p, float m){ this->freq(f); this->phase(p); mod(m); }
TEMS inline void LFO<Ts>::mod(double n){ modi = scl::castIntRound(n * 4294967296.); }

#define DEF(name, exp) TEMS inline float LFO<Ts>::name{ float r = exp; incPhase(); return r; }


TEMS inline float LFO<Ts>::line2(){
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

TEMS inline float LFO<Ts>::line2U(){
	return line2()*0.5f+0.5f;
}



DEF(down(),		scl::rampDown(phaseI()))
DEF(even3(),	scl::rampUp(phaseI()); static float c=-1.50*sqrt(3.); r *= (1-r*r)*c;)
DEF(even5(),	scl::rampUp(phaseI()); static float c=-1.25*pow(5.,0.25); r *= (1-scl::pow4(r))*c;)
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

//TEMS inline float LFO<Ts>::imp(){ return this->cycle() ? 1.f : 0.f; }
TEMS inline float LFO<Ts>::imp(){ 
	float r = phaseI() < this->phaseIncI() ? 1.f : 0.f;
	incPhase();
	return r; 
}

TEMS inline bool LFO<Ts>::seq(){
	uint32_t prev = phaseI();
	incPhase();
	if( (phaseI() ^ prev) & 0xf8000000 ){
		if( (modi >> (phaseI()>>27)) & 0x1 ) return true;
	}
	return false;
}



//---- Impulse

TEM Impulse<Tv, Ts>::Impulse(Tv frq, Tv phase, Tv harmonics)
:	AccumPhase<Tv, Ts>(frq, phase), mPrev((Tv)0)
{
	onResync(1);
	this->harmonics(harmonics);
}

TEM inline void Impulse<Tv, Ts>::harmonics(Tv value){
	mN = mNDesired = scl::floor(value);
	mNFrac = value - mN;
}

TEM inline void Impulse<Tv, Ts>::harmonicsMax(){ harmonics(maxHarmonics()); }

TEM inline void Impulse<Tv, Ts>::antialias(){
	float maxN = scl::floor(maxHarmonics());
	mN = mNDesired > maxN ? maxN : mNDesired;
}

TEM inline Tv Impulse<Tv, Ts>::maxHarmonics(){ return mSPU_2 / this->mFreq; }

#define EPS 0.000001
TEM inline Tv Impulse<Tv, Ts>::operator()(){
	Tv theta = this->nextPhase();

	Tv result;
	Tv denom = scl::sinT9(theta * (Tv)0.5);
	if( scl::abs(denom) < (Tv)EPS ){ result = (Tv)1; /*printf("Impulse::(): oops\n");*/ }
	else{
		Tv nphase = scl::wrapPhase(theta * (mN + (Tv)0.5));
		//result = (scl::sinT7(nphase) / denom - (Tv)1) * (Tv)0.5;
		
		result = ((scl::sinT7(nphase) - denom) / (denom*mN)) * (Tv)0.5;
		
//		(n/d - 1)*0.5/N
//		n/d - d/d
//		(n-d)/d*0.5/N
	}
	
	//result -= 2.f * cos(theta * mN) * (1.f - mNFrac);	// removes clicks!
	
	return result;
}

TEM inline Tv Impulse<Tv, Ts>::odd(){
	Tv theta = this->nextPhase();
	
	Tv result;
	Tv n = scl::ceil(mN, (Tv)2, (Tv)0.5);	// get next highest even for formula
										// wave has odd harmonics 1,3, ..., n-1 with peak amp n

	Tv denom = scl::cosT8(scl::wrapPhaseOnce(theta - M_PI_2));	// this is more precise near zero-crossings
	//Tv denom = scl::sinT9(theta);
	if( scl::abs(denom) < (Tv)EPS ){
		if( theta > M_PI ) theta -= M_2PI;
		result = (theta > -M_PI_2 && theta < M_PI_2) ? 1 : -1;
		//printf("Impulse::odd(): oops\n");
	}
	else result = scl::sinT7(scl::wrapPhase(n * theta)) / (denom * n);
	
	return result;
}
#undef EPS

TEM inline Tv Impulse<Tv, Ts>::saw(Tv i){ return mPrev = (*this)() + i * mPrev; }
TEM inline Tv Impulse<Tv, Ts>::square(Tv i){ return mPrev = odd() + i * mPrev; }

TEM void Impulse<Tv, Ts>::onResync(double r){
	recache();
	freq((*this).freq());
}

TEM void Impulse<Tv, Ts>::recache(){
	AccumPhase<Tv, Ts>::recache();
	mSPU_2 = (Tv)(Synced::spu() * 0.5);
}



//---- DSF

TEM DSF<Tv, Ts>::DSF(Tv frq, Tv freqRatioA, Tv ampRatioA, Tv harmonicsA)
	: AccumPhase<Tv, Ts>(frq)
{
	freq(frq);
	freqRatio(freqRatioA);
	harmonics(harmonicsA);
	ampRatio(ampRatioA);

	mBeta = 0.f;
}

TEM inline void DSF<Tv, Ts>::freq(Tv value){
	AccumPhase<Tv, Ts>::freq(value);
	updateBetaInc();	
}

TEM inline void DSF<Tv, Ts>::freqRatio(Tv value){
	mFreqRatio = value;
	updateBetaInc();
}

TEM inline void DSF<Tv, Ts>::ampRatio(Tv value){
	if(value != mA){
		mA = value;
		mASqP1 = mA * mA + 1.f;
		updateAPow();
	}
}

TEM inline void DSF<Tv, Ts>::harmonics(Tv value){
	if(value != mN){
		mN = mNDesired = value;
		updateAPow();
	}
}

TEM inline void DSF<Tv, Ts>::harmonicsMax(){ harmonics(maxHarmonics()); }

TEM inline void DSF<Tv, Ts>::antialias(){
	Tv maxN = maxHarmonics();
	if(mNDesired > maxN)	mN = maxN;
	else					mN = mNDesired;
	updateAPow();
}

TEM inline Tv DSF<Tv, Ts>::ampRatio(){ return mA; }
TEM inline Tv DSF<Tv, Ts>::freqRatio(){ return mFreqRatio; }
TEM inline Tv DSF<Tv, Ts>::harmonics(){ return mN; }

TEM inline Tv DSF<Tv, Ts>::maxHarmonics(){
	return scl::floor(((Tv)this->spu() * (Tv)0.5 / this->mFreq - (Tv)1) / mFreqRatio + (Tv)1);
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
TEM inline Tv DSF<Tv, Ts>::operator()(){
	Tv theta = this->nextPhase();
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

TEM inline void DSF<Tv, Ts>::updateAPow(){ mAPow = pow(mA, mN); }
TEM inline void DSF<Tv, Ts>::updateBetaInc(){ mBetaInc = this->mPhaseInc * mFreqRatio; }

TEM void DSF<Tv, Ts>::onResync(double r){
	AccumPhase<Tv, Ts>::onResync(r);
	freq(this->mFreq);
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




/*
template <class T=float, class Ts=Synced>
class Sine : public Ts{
public:
	Sine(T f, T p=0){ set(f,p); Ts::initSynced(); }

	void freq (T v){ p1 = v*M_2PI; }
	void phase(T v){ p0 = v*M_2PI; }
	void set(T f, T p){ freq(f); phase(p); }

	float operator()(){
		p0 = scl::wrap<T>(p0, M_2PI);
		T r = sin(p0);
		p0 += p1;
		return r;
	}
	
	// change in sampling rate (ratio = new/old)
	virtual void onResync(double ratio){ p1 /= ratio; }

protected:
	T p0, p1;	// 0 and 1 derivative of phase
};



template <class T=double, class Ts=Synced>
class SineAM : public Ts{
public:
	SineAM(T f1, T f2, T p1=0, T p2=0): osc1(f1, p1), osc2(f2, p2){ Ts::initSynced(); }

	T operator()(){ return osc1() * osc2(); }

	void freq(T v1, T v2){
		osc1.freq(v1*Ts::ups());
		osc2.freq(v2*Ts::ups());
	}

	virtual void onResync(double ratio){
		osc1.onResync(ratio);
		osc2.onResync(ratio);
	}

protected:
	Sine<T, Synced1> osc1, osc2;
};



// Free resonator
template <class T=double, class Ts=Synced>
class SineRes : public gen::RSin<T>, Ts{
public:

	SineRes(T frq=440, T amp=1, T phs=0){ Ts::initSynced(); set(frq, amp, phs); }

	/// Set all control parameters
	void set(T frq, T amp, T phs=0){ gen::RSin<T>::set(frq*Ts::ups(), phs, amp); }
	
	/// Generate next samples adding into a buffer
	template <class V>
	void add(V * dst, uint32_t n){ for(uint32_t i=0; i<n; ++i) dst[i] += (*this)(); }

	// This might not be a good idea because we don't know amplitude...
	// It also doesn't even work right!
	//virtual void onResync(double ratio){ set(gen::RSin<T>::freq()/ratio, 1, 0); }

};


// Multiple SineRes (has-a)
// This object contains an array of SineRes's.
// Good:	can synchronize SineRes easily
// Bad:		must rewrite all SineRes parameter setting methods

template <class T=double, class Ts=Synced>
class SineResN : public Ts{
public:

	SineResN(uint32_t n){ mSineRes.resize(n); }

	/// Set all control parameters of an oscillator
	void set(uint32_t i, T frq, T amp, T phs=0){ mSineRes[i].set(frq*Ts::ups(), amp, phs); }

protected:
	Array<SineRes<T, Synced1> > mSineRes;
};


// Multiple SineRes (is-a)
// This object is an array of SineRes's.
//
// Good:	do not have to rewrite parameter setting methods
// Bad:		how do SineRes know about the sampling rate?
//
// Here we have to explicitly rescale frequency parameters, i.e.:
//		sineResN[i].set(f*Sync::master().spu(), 1, 0);
//
// The frequency value is rescaled to the interval [0, 0.5].

template <class T=double, class Ts=Synced>
class SineResN : public Array<SineRes<T, Synced1> >, Ts{
public:
	typedef Array<SineRes<T, Synced1> > super;

	SineResN(uint32_t n): super(n){}
};



// An older, more complicated version...
template <class T=double>
class SineRes : public Synced{
public:

	/// @param[in]	frq		Initial frequency.
	/// @param[in]	amp		Initial amplitude.
	/// @param[in]	phs		Initial unit phase [0, 1).
	SineRes(T frq=440, T amp=1, T phs=0): mFreq(frq), mAmp(amp){ initSynced(); phase(phs); }
	
	/// Set unit phase
	void phase(T v){		
		v = v*M_2PI - mAngle;
		
		// compute previous two phases
		rsin.val2 = sin(v - mAngle) * mAmp;
		rsin.val  = sin(v         ) * mAmp;
	}
	
	/// Set all control parameters
	void set(T frq, T amp, T phs=0){
		freq(frq);
		mAmp = amp;
		phase(phs);
	}
	
	/// Generate next sample
	T operator()(){	return rsin(); }
	
	/// Generate next samples adding into a buffer
	template <class V>
	void add(V * dst, uint32_t n){ for(uint32_t i=0; i<n; ++i) dst[i] += (*this)(); }
	
	virtual void onResync(double r){ freq(mFreq); }

protected:
	T mFreq, mAngle, mAmp;
	gen::RSin<T> mRSin;	// recursive sinusoid generator


	/// Set frequency
	void freq(T v){
		mFreq = v;
		mAngle = M_2PI * v * ups();
		rsin.mul = (T)2 * cos(mAngle);	// [-2, 2]
		
		//recomputePrev();	// recompute previous values so amplitude doesn't change
	}

	// recompute o2 assuming o1 was the last output
	void recomputePrev(){
		if(mAmp == (T)0) return;
		T& v1 = rsin.val;
		T& v2 = rsin.val2;
		
		// crux of the problem is we need to get the current phase
		float phs = asin(v1 / mAmp);	// normalize o1 and get its angle
		
		// NOT WORKING: this still has degenerate cases
		if(v1 < v2) phs = M_PI - phs;	// flip angle if real < 0
		v2 = sin(phs - mAngle) * mAmp;
	}
};

*/



#undef TEM
#undef TEMS

} // gam::
#endif

