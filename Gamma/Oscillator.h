#ifndef GAMMA_OSCILLATOR_H_INC
#define GAMMA_OSCILLATOR_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

///\defgroup Oscillators

#include "Gamma/gen.h"
#include "Gamma/scl.h"
#include "Gamma/tbl.h"
#include "Gamma/Strategy.h"
#include "Gamma/Domain.h"
#include "Gamma/Types.h"

namespace gam{


/// Fixed-point phase accumulator

/// This is a linear phase accumulator that uses integer (fixed-point) 
/// arithmetic. The advantage of using fixed-point versus floating-point is that
/// the phase is wrapped automatically when the integer overflows. As long as
/// we used unsigned integers, this wrapping behavior is well-defined--all 
/// results of addition are taken modulo the maximum size of the integer.
/// The frequency resolution is SR/2^32. For SR=44100 Hz, this equates to
/// a frequency resolution of ~10^-5 Hz which has a period of ~27 hours.
///
/// \tparam Sp	Phase increment strategy (e.g., phsInc::Loop, phsInc::Oneshot)
/// \tparam Td	Domain type
/// \ingroup Oscillators     
template <class Sp = phsInc::Loop, class Td = DomainObserver>
class Accum : public Td {
public:

	/// \param[in] frq		Frequency
	/// \param[in] phs		Phase in [0, 1)
	Accum(float frq=0, float phs=0);


	void freq(float v);				///< Set frequency
	void freqI(uint32_t v);			///< Set fixed-point frequency
	void phase(float v);			///< Set phase from [0, 1) of one period
	void phaseMax();				///< Set phase to maximum value
	void phaseAdd(float v);			///< Add value to phase [0, 1)
	void period(float v);			///< Set period length
	void reset(){ mPhaseI=0; mSp.reset(); }	///< Reset phase accumulator
	Sp& phsInc(){ return mSp; }		///< Get phase increment strategy

	/// Returns true if tap is done
	bool done() const { return mSp.done(phaseI()); }

	float freq() const;				///< Get frequency
	uint32_t freqI() const;			///< Get fixed-point frequency
	float freqUnit() const;			///< Get frequency in [0, 1)
	float period() const;			///< Get period
	float phase() const;			///< Get phase in [0, 1)
	uint32_t phaseI() const;		///< Get fixed-point phase

	/// Returns 0x80000000 on phase wrap, 0 otherwise
	
	/// The return value can be used as a bool.  It's an integer because it
	/// saves a conditional check converting to a bool.
	uint32_t cycle();
	uint32_t operator()();			///< Alias of cycle()

	uint32_t nextPhase();			///< Increment phase and return updated phase
	uint32_t nextPhase(float freqOffset);
	uint32_t cycles();				///< Get 1 to 0 transitions of all accumulator bits
	uint32_t once();

	/// Returns sequence of 32 triggers based on a pattern of bits

	/// The 5 MSBs of the phase are used to determine how much to right
	/// bit-shift the pattern value to check its bit.
	/// This means we get a pattern of 32 triggers over one oscillation period.
	/// If the oscillation frequency is positive then the pattern bits are
	/// scanned from the LSB (bit 0) to the MSB (bit 31). If the frequency is
	/// negative, then the bits are scanned (perhaps more intuitively) from the 
	/// MSB to the LSB.
	/// The following shows all 4-bit patterns made from hex values:
	///\verbatim
	///		hex pattern		hex	pattern		hex pattern		hex pattern
	///		0	. . . .		4	. / . .		8	/ . . .		c	/ / . .
	///		1	. . . /		5	. / . /		9	/ . . /		d	/ / . /
	///		2	. . / .		6	. / / .		a	/ . / .		e	/ / / .
	///		3	. . / /		7	. / / /		b	/ . / /		f	/ / / /			\endverbatim
	bool seq(uint32_t pattern);

	virtual void onDomainChange(double r);

//protected:
private:
	float mFreq;		// Current frequency
	uint32_t mPhaseI;	// Current fixed-point phase in [0, 2^32)
	uint32_t mFreqI;	// Current fixed-point frequency
	Sp mSp;

	uint32_t mapFI(float v) const;	// convert unit floating-point to fixed-point integer
	double mapIF(uint32_t v) const;	// convert fixed-point integer to unit floating-point
	uint32_t mapFreq(float v) const;
};

#define ACCUM_INHERIT\
	using Accum<Sp,Td>::phaseI;\
	using Accum<Sp,Td>::freqI;\
	using Accum<Sp,Td>::nextPhase;


/// Linear sweep in interval [0,1)

/// \tparam Sp	Phase increment strategy (e.g., phsInc::Loop, phsInc::Oneshot)
/// \tparam Td	Domain type
///\ingroup Oscillators 
template <class Sp = phsInc::Loop, class Td = DomainObserver>
class Sweep : public Accum<Sp, Td> {
public:
	/// \param[in] frq		Frequency
	/// \param[in] phs		Phase in [0,1)
	Sweep(float frq=440, float phs=0): Base(frq, phs){}

	float operator()(){ Base::cycle(); return Base::phase(); }

private: typedef Accum<Sp,Td> Base;
};

    
/// Floating-point phase accumulator with output in [-A, A)

/// \tparam Tv	Value (sample) type
/// \tparam Td	Domain type
///\ingroup Oscillators     
template <class Tv = gam::real, class Td = DomainObserver>
class AccumPhase : public Td{
public:
	/// \param[in] frq		Frequency
	/// \param[in] phs		Phase in [0, 1)
	/// \param[in] amp		Amplitude
	AccumPhase(Tv frq=440, Tv phs=0, Tv amp=M_PI);

	
	/// Generate next sample. Stored phase is post-incremented.
	Tv nextPhase();
	
	/// Generate next sample with a frequency offset
	Tv nextPhase(Tv frqOffset);

	void freq(Tv v);		///< Set frequency
	void period(Tv v);		///< Set period length
	void phase(Tv v);		///< Set phase from [0, 1) of one period
	void phaseAdd(Tv v);	///< Add value to unit phase
	void amp(Tv v);			///< Set amplitude
	
	Tv freq() const;		///< Get frequency
	Tv period() const;		///< Get period
	Tv phase() const ;		///< Get normalized phase in [0, 1)
	Tv amp() const;			///< Get amplitude
	
	virtual void onDomainChange(double r);
	void print(FILE * fp = stdout, const char * append = "\n");
	
protected:
	Tv mInc, mPhase, mAmp;
	Tv mFreqToInc;
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
///
/// \tparam Tv	Table element type
/// \tparam Si	Interpolation strategy
/// \tparam Sp	Phase increment strategy
/// \tparam Td	Domain type
/// \ingroup Oscillators
/// \ingroup Envelopes
/// \sa Other ways to synthesize sine waves: TableSine, CSine, LFO, Sine, SineR
/// \sa Functions for building waveforms in tables with additive synthesis: 
///		addSine, addSines, addSinesPow, addWave
template<
	class Tv = gam::real,
	template<class> class Si = ipl::Linear,
	class Sp = phsInc::Loop,
	class Td = DomainObserver
>
class Osc : public Accum<Sp,Td>, public ArrayPow2<Tv>{
public:

	/// Constructor that allocates an internal table

	/// \param[in]	frq			Frequency
	/// \param[in]	phs			Phase in [0, 1)
	/// \param[in]	size		Size of table (actual number is power of 2 ceiling)
	Osc(float frq=440, float phs=0, uint32_t size=512)
	:	Base(frq, phs), ArrayPow2<Tv>(size, Tv())
	{}

	/// Constructor that references an external table

	/// \param[in]	frq			Frequency
	/// \param[in]	phs			Phase in [0, 1)
	/// \param[in]	src			A table to use as a reference
	Osc(float frq, float phs, ArrayPow2<Tv>& src)
	:	Base(frq, phs), ArrayPow2<Tv>(src.elems(), src.size())
	{}


	/// Generate next sample
	Tv operator()(){
		this->nextPhase(); return val();
	}

	/// Get current value
	Tv val() const { return mIpol(*this, phaseI()); }
	
	/// Add sine to table
	
	/// \param[in] cycles	number of cycles
	/// \param[in] amp		amplitude
	/// \param[in] phs		unit phase, [0, 1)
	Osc& addSine(double cycles, double amp=1, double phs=0){
		double f = cycles/this->size();
		for(unsigned i=0; i<this->size(); ++i){
			double p = (f*i + phs)*M_2PI;
			(*this)[i] += sin(p) * amp;
		}
		return *this;
	}

	/// Zero table elements
	void zero(){ this->assign(Tv(0)); }

	/// Get reference to table
	ArrayPow2<Tv>& table(){ return *this; }

//	using ArrayPow2<Tv>::elems; using ArrayPow2<Tv>::size;
protected:
	Si<Tv> mIpol;
private:
	ACCUM_INHERIT
	typedef Accum<Sp,Td> Base;
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
///
/// \tparam Tv	Value (sample) type
/// \tparam Td	Domain type
/// \ingroup Oscillators 
///  \sa Osc, TableSine, CSine, LFO, Sine, SineR
template<class Tv = gam::real, class Td = DomainObserver>
class CSine : public Td{
public:

	typedef Complex<Tv> complex;

	/// \param[in] frq		Frequency
	/// \param[in] amp		Amplitude
	/// \param[in] dcy		-60 dB decay length, in units (or negative for no decay)
	/// \param[in] phs		Phase in [0, 1)
	CSine(Tv frq=440, Tv amp=1, Tv dcy=-1, Tv phs=0);


	complex val;				///< Current value

	complex operator()();		///< Generate next sample

	void amp(Tv val);			///< Set amplitude
	void decay(Tv val);			///< Set -60 dB decay length, in units (or negative for no decay)
	void freq(Tv val);			///< Set frequency
	void freq(const complex& val){ mInc=val; }
	void phase(Tv radians);		///< Set phase
	void reset();				///< Reset amplitude and set phase to 0
	void set(Tv frq, Tv phs, Tv amp, Tv dcy);

	Tv amp() const {return mAmp;}		///< Get amplitude
	Tv decay() const {return mDcy60;}	///< Get decay length
	Tv freq() const {return mFreq;}		///< Get frequency

	virtual void onDomainChange(double r);

protected:
	Tv mAmp, mFreq, mDcy60;
	Tv mDcy;			// phasor amp
	complex mInc;		// rotation increment
};



/// Computed sine wave oscillator.

/// This oscillator uses a polynomial approximation to compute sine values. 
/// Computation time is about as much as a linearly-interpolating table lookup.
/// In addition, polynomial approximations are often more spectrally pure than 
/// table lookup methods since the distortion arises as harmonics.
///
/// \tparam Tv	Value (sample) type
/// \tparam Td	Domain type
/// \ingroup Oscillators 
/// \sa Osc, TableSine, CSine, LFO, SineR
template<class Tv = gam::real, class Td = DomainObserver>
class Sine : public AccumPhase<Tv,Td> {
public:
	/// \param[in]	frq		Frequency
	/// \param[in]	phs		Phase in [0, 1)
	Sine(Tv frq=440, Tv phs=0) : AccumPhase<Tv,Td>(frq, phs, Tv(1)){}

	/// Generate next sample with a frequency offset
	Tv operator()(Tv frqOffset = Tv(0)){
		//return ::sin(this->nextPhase(frqOffset)/*M_PI*/); // 260%
		//return scl::sinT7(this->nextPhase(frqOffset)/*M_PI*/); // 140%
		return scl::sinP9(this->nextPhase(frqOffset)); // 110%
		//return scl::sinP7(this->nextPhase(frqOffset)); // 100%
		//return scl::sinFast(this->nextPhase(frqOffset)/*M_PI*/); // 95%
	}
};



/// Sine oscillator based on an efficient recursion equation.

/// This oscillator is based on a recursion equation requiring only one
/// multiply and add per sample computation. While computation time and quality
/// are near ideal, parameter updates are expensive and 64-bit precision is
/// required to prevent growing or decaying in amplitude over time.  This
/// generator is ideal in situations where a stationary sinusoid is all that is
/// required, e.g. a grain or modulator.
///
/// \tparam Tv	Value (sample) type
/// \tparam Td	Domain type
/// \ingroup Oscillators 
/// \sa Osc, TableSine, CSine, LFO, Sine, SineRs (Synthesizes multiple sines)
template <class Tv = double, class Td = DomainObserver>
class SineR : public gen::RSin<Tv>, Td{
public:

	/// \param[in]	frq		Frequency
	/// \param[in]	amp		Amplitude
	/// \param[in]	phs		Phase in [0, 1)
	SineR(Tv frq=440, Tv amp=1, Tv phs=0){ set(frq, amp, phs); }

	/// Get frequency
	Tv freq() const { return Base::freq() * Td::spu(); }

	/// Set amplitude and phase
	void ampPhase(Tv a=1, Tv p=0){ set(freq(), a, p); }

	/// Set frequency
	void freq(Tv v){ Base::freq(v*Td::ups()); }

	/// Set all control parameters
	void set(Tv frq, Tv amp, Tv phs=0){ Base::set(frq*Td::ups(), phs, amp); }


	virtual void onDomainChange(double ratio){ Base::freq(Base::freq()/ratio); }

private:
	typedef gen::RSin<Tv> Base;
};



/// Multiple SineRs

/// This is a dynamically-sized bank of SineR objects.
///
/// \tparam Tv	Value (sample) type
/// \tparam Td	Domain type
/// \ingroup Oscillators 
/// \sa SineR (synthesizes a single sine)
template <class Tv = double, class Td = DomainObserver>
class SineRs : public Array<SineR<Tv, Domain1> >, Td{
public:

	SineRs(){}

	/// \param[in]	num		Number of resonators
	SineRs(unsigned num): Base(num){ Td::refreshDomain(); }

	/// Generate next sum of all oscillators
	Tv operator()(){
		Tv r = Tv(0);
		for(unsigned i=0; i<this->size(); ++i) r += (*this)[i]();
		return r;
	}

	/// Get last output of oscillator i
	Tv last(unsigned i) const { return (*this)[i].val; }

	/// Set all control parameters of oscillator i
	void set(unsigned i, Tv frq, Tv amp=1, Tv phs=0){
		(*this)[i].set(frq*Td::ups(), amp, phs);
	}


	virtual void onDomainChange(double ratio){
		for(unsigned i=0; i<this->size(); ++i){
			(*this)[i].onDomainChange(ratio);
		}
	}
private:
	typedef Array<SineR<Tv, Domain1> > Base;
};



/// Damped sine oscillator based on an efficient recursion equation.

/// This oscillator is similar to SineR, however, it has an extra multiply
/// in its sample generation to allow the oscillator to decay.
///
/// \tparam Tv	Value (sample) type
/// \tparam Td	Domain type
/// \ingroup Oscillators 
/// \sa SineR, SineDs
template <class Tv = double, class Td = DomainObserver>
class SineD : public gen::RSin2<Tv>, Td{
public:

	/// \param[in]	frq		Frequency
	/// \param[in]	amp		Amplitude
	/// \param[in]	dcy		T60 decay length (negative == no decay)
	/// \param[in]	phs		Phase in [0, 1)
	SineD(Tv frq=440, Tv amp=1, Tv dcy=-1, Tv phs=0)
	{ set(frq, amp, dcy, phs); }


	/// Get frequency
	Tv freq() const { return Base::freq() * Td::spu(); }

	/// Set amplitude and phase
	void ampPhase(Tv a=1, Tv p=0){ set(freq(), a, Base::decay(), p); }
	
	/// Set all control parameters

	/// \param[in]	frq		Frequency
	/// \param[in]	amp		Amplitude
	/// \param[in]	dcy		T60 decay length (negative == no decay)
	/// \param[in]	phs		Phase in [0, 1)
	void set(Tv frq, Tv amp, Tv dcy, Tv phs=0){
		Base::set(
			frq * Td::ups(),
			phs,
			dcy > Tv(0) ? Tv(scl::radius60(dcy, Td::ups())) : Tv(1),
			amp
		);
	}

	virtual void onDomainChange(double ratio){
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

/// This is a dynamically-sized bank of SineD objects.
///
/// \tparam Tv	Value (sample) type
/// \tparam Td	Domain type
/// \ingroup Oscillators 
/// \sa SineD
template <class Tv = double, class Td = DomainObserver>
class SineDs : public Array<SineD<Tv, Domain1> >, Td{
public:

	SineDs(){}

	/// \param[in]	num		Number of resonators
	SineDs(unsigned num): Base(num){
		Td::refreshDomain(); 
		for(unsigned i=0; i<num; ++i) set(i, 0,0,0);
	}


	/// Generate next sum of all oscillators
	Tv operator()(){
		Tv r=Tv(0);
		for(unsigned j=0; j<this->size(); ++j) r+=(*this)[j]();
		return r;
	}

	/// Get last output of oscillator i
	Tv last(unsigned i) const { return (*this)[i].val; }

	/// Set all control parameters of oscillator i

	/// \param[in] i		index of oscillator
	/// \param[in] frq		Frequency
	/// \param[in] amp		Amplitude
	/// \param[in] dcy		T60 decay length (negative == no decay)
	/// \param[in] phs		Phase in [0, 1)
	void set(unsigned i, Tv frq, Tv amp, Tv dcy, Tv phs=0){
		(*this)[i].set(frq*Td::ups(), amp, dcy*Td::spu(), phs);
	}


	virtual void onDomainChange(double ratio){
		for(unsigned i=0; i<this->size(); ++i){
			(*this)[i].onDomainChange(ratio);
		}
	}
private:
	typedef Array<SineD<Tv, Domain1> > Base;
};



/// Low-frequency oscillator

/// This object generates various waveform types by mapping the output of a 
/// an accumulator through mathematical functions. The resulting waveforms are
/// non-band-limited.
///
/// \tparam Sp	Phase increment strategy (e.g., phsInc::Loop, phsInc::Oneshot)
/// \tparam Td	Domain type
/// \ingroup Oscillators 
/// \sa Osc, TableSine, CSine, Sine, SineR
template <class Sp = phsInc::Loop, class Td = DomainObserver>
class LFO : public Accum<Sp,Td>{
public:

	LFO();
	
	/// \param[in] frq		Frequency
	/// \param[in] phase	Phase in [0, 1)
	/// \param[in] mod		Modifier amount in [0, 1)
	LFO(float frq, float phase=0, float mod=0.5);


	/// Set frequency, phase and modifier amount
	LFO& set(float f, float p, float m);

	LFO& mod(double n);		///< Set modifier from unit value
	LFO& modI(uint32_t v);	///< Set modifier from integer

	/// Get modifier value
	uint32_t modI() const { return mMod; }
	double mod() const { return mMod / 4294967296.; }

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

	ACCUM_INHERIT
private:
	typedef Accum<Sp,Td> Base;
	uint32_t mMod;			// Modifier parameter
};



/// Sum of cosine waves

/// This produces a finite sum of harmonic, equi-amplitude cosine waves that 
/// approach the shape of a periodic impulse train. The algorithm uses a 
/// closed-form solution to the sum thus its computational complexity is always
/// O(1) (constant) regardless of the number of harmonics.
///
/// \tparam Tv	Value (sample) type
/// \tparam Td	Domain type
/// \ingroup Oscillators 
template<class Tv = gam::real, class Td = DomainObserver>
class Buzz : public AccumPhase<Tv,Td> {
public:

	/// \param[in]	frq			Frequency
	/// \param[in]	phase		Phase in [0, 1)
	/// \param[in]	harmonics	Number of harmonics
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

	virtual void onDomainChange(double r);

protected:
	Tv mAmp;			// amplitude normalization factor
	Tv mN;				// # harmonics
	Tv mNDesired;		// desired number of harmonics
	Tv mNFrac;		
	Tv mSPU_2;			// cached locals
	Tv mPrev;			// previous output for integration
	void setAmp();
private: typedef AccumPhase<Tv,Td> Base;
};



/// Band-limited impulse train

/// This produces a Fourier representation of an impulse train where the number
/// of harmonics is adjusted automatically to prevent aliasing.
///
/// \tparam Tv	Value (sample) type
/// \tparam Td	Domain type
/// \ingroup Oscillators
/// \sa Buzz
template <class Tv = gam::real, class Td = DomainObserver>
struct Impulse : public Buzz<Tv,Td>{
public:

	/// \param[in] frq		Frequency
	/// \param[in] phs		Phase, in [0, 1)
	Impulse(Tv frq=440, Tv phs=0): Base(frq, phs){ onDomainChange(1); }

	/// Set frequency
	void freq(Tv v){ Base::freq(v); Base::harmonicsMax(); }

	virtual void onDomainChange(double r){
		Base::onDomainChange(r); freq(AccumPhase<Tv,Td>::freq()); }

private: typedef Buzz<Tv,Td> Base;
};



/// Band-limited saw wave

/// This produces a Fourier representation of a saw wave where the number of
/// harmonics is adjusted automatically to prevent aliasing.
/// Due to numerical issues, this generator should not be used for producing 
/// very low frequency modulation signals. For that purpose, it is better to use
/// the LFO class.
///
/// \tparam Tv	Value (sample) type
/// \tparam Td	Domain type
/// \ingroup Oscillators 
template <class Tv = gam::real, class Td = DomainObserver>
struct Saw : public Impulse<Tv,Td> {

	/// \param[in] frq		Frequency
	/// \param[in] phs		Phase, in [0, 1)
	Saw(Tv frq=440, Tv phs=0): Impulse<Tv, Td>(frq, phs){}

	/// Generate next sample
	
	/// \param[in] itg		Leaky integration factor
	///
	Tv operator()(Tv itg=0.997){ return Impulse<Tv,Td>::saw(itg); }
};



/// Band-limited square wave

/// This produces a Fourier representation of a square wave where the number of
/// harmonics is adjusted automatically to prevent aliasing.
/// Due to numerical issues, this generator should not be used for producing 
/// very low frequency modulation signals. For that purpose, it is better to use
/// the LFO class.
///
/// \tparam Tv	Value (sample) type
/// \tparam Td	Domain type
/// \ingroup Oscillators 
template <class Tv = gam::real, class Td = DomainObserver>
struct Square : public Impulse<Tv,Td> {

	/// \param[in] frq		Frequency
	/// \param[in] phs		Phase, in [0, 1)
	Square(Tv frq=440, Tv phs=0) : Impulse<Tv,Td>(frq, phs){}

	/// Generate next sample
	
	/// \param[in] itg		Leaky integration factor
	///
	Tv operator()(Tv itg=0.997){ return Impulse<Tv,Td>::square(itg); }
};



/// Discrete summation formula (DSF) oscillator

/// This produces a finite sum of harmonics whose amplitudes follow a geometric
/// series. The amplitude of harmonic i is ar^i where 'ar' is called the 
/// amplitude ratio. The frequency of harmonic i is (i * fr + 1) where 'fr' is
/// called the frequency ratio. Harmonics run from i=0 (the fundamental) to
/// the maximum specified harmonic.
///
/// \tparam Tv	Value (sample) type
/// \tparam Td	Domain type
/// \ingroup Oscillators 
template<class Tv = gam::real, class Td = DomainObserver>
class DSF : public AccumPhase<Tv,Td> {
public:

	/// \param[in]	frq			Frequency
	/// \param[in]	freqRatio	Frequency ratio of partials
	/// \param[in]	ampRatio	Amplitude ratio of partials
	/// \param[in]	harmonics	Number of harmonics
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
	
	virtual void onDomainChange(double r);

protected:
	typedef AccumPhase<Tv,Td> Base;

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
/// \ingroup Oscillators 
class ImpulseFast : public DomainObserver {
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

//---- Accum
    
template<class St, class Td> Accum<St,Td>::Accum(float f, float p)
:	mFreq(f), mFreqI(0)
{
	Td::refreshDomain();
	phase(p);
	//(p >= 1.f) ? phaseMax() : this->phase(p);
}

// 32-bit float is good enough here since [0.f, 1.f) uses 29 bits.
template<class St, class Td> inline uint32_t Accum<St,Td>::mapFI(float v) const {
	//return scl::unitToUInt(v);
	//return (uint32_t)(v * 4294967296.);
	return castIntRound(v * 4294967296.);
}

template<class St, class Td> inline double Accum<St,Td>::mapIF(uint32_t v) const {
	return v/4294967296.;
	//return uintToUnit<float>(v); // not enough precision
}

template<class St, class Td> inline uint32_t Accum<St,Td>::mapFreq(float v) const {
	return mapFI(v * Td::ups());
}

template<class St, class Td> void Accum<St,Td>::onDomainChange(double r){ //printf("Accum: onDomainChange (%p)\n", this);
	uint32_t fprev = mFreqI;
	freq(mFreq);
	
	// ensure phase will be correct value upon next increment
	mPhaseI = mPhaseI + fprev - mFreqI;
}

template<class St, class Td> inline void Accum<St,Td>::freq(float v){
	mFreq = v;
	mFreqI= mapFreq(v);
}

template<class St, class Td> inline void Accum<St,Td>::freqI(uint32_t v){
	mFreqI= v;
	mFreq = mapIF(v) * Td::spu();
}

template<class St, class Td> inline void Accum<St,Td>::period(float v){ freq(1.f/v); }
template<class St, class Td> inline void Accum<St,Td>::phase(float v){ mPhaseI = mapFI(v) - mFreqI; }
template<class St, class Td> inline void Accum<St,Td>::phaseAdd(float v){ mSp(mPhaseI, mapFI(v)); }
template<class St, class Td> inline void Accum<St,Td>::phaseMax(){ mPhaseI = 0xffffffff; }

template<class St, class Td> inline float Accum<St,Td>::freq() const { return mFreq; }
template<class St, class Td> inline uint32_t Accum<St,Td>::freqI() const { return mFreqI; }
template<class St, class Td> inline float Accum<St,Td>::freqUnit() const { return mapIF(mFreqI); }
template<class St, class Td> inline float Accum<St,Td>::period() const { return 1.f/freq(); }
template<class St, class Td> inline float Accum<St,Td>::phase() const { return mapIF(mPhaseI); }
template<class St, class Td> inline uint32_t Accum<St,Td>::phaseI() const { return mPhaseI; }

template<class St, class Td> inline uint32_t Accum<St,Td>::nextPhase(float frqOffset){
	return mSp(mPhaseI, mFreqI + mapFreq(frqOffset));
}

template<class St, class Td> inline uint32_t Accum<St,Td>::nextPhase(){ return mSp(mPhaseI, mFreqI); }

template<class St, class Td> inline uint32_t Accum<St,Td>::operator()(){ return cycle(); }

template<class St, class Td> inline uint32_t Accum<St,Td>::cycle(){ return cycles() & 0x80000000; }

//template<class St, class Td> inline uint32_t Accum<St,Td>::cycle(uint32_t mask){
//	return cycles() & mask;
//}

template<class St, class Td> inline uint32_t Accum<St,Td>::cycles(){
	uint32_t prev = phaseI();
	nextPhase();
	return ~phaseI() & prev;
}

template<class St, class Td> inline uint32_t Accum<St,Td>::once(){
	uint32_t prev = phaseI();
	uint32_t c = cycle();
	if(c) mPhaseI = prev;
	return c;
}

template<class St, class Td> inline bool Accum<St,Td>::seq(uint32_t pat){
	uint32_t prev = phaseI();
	nextPhase();

	// Did any of the 5 MSBs change?
	if((phaseI() ^ prev) & 0xf8000000){
		return (pat >> (phaseI()>>27)) & 0x1;
	}
	return false;
}





//---- AccumPhase

template<class Tv, class Td>
AccumPhase<Tv, Td>::AccumPhase(Tv f, Tv p, Tv a)
:	mInc(f), mAmp(a), mFreqToInc(1)
{
	Td::refreshDomain();
	this->phase(p);
}

template<class Tv, class Td> inline Tv AccumPhase<Tv, Td>::mapFreq(Tv v) const { return v*mFreqToInc; }
template<class Tv, class Td> inline Tv AccumPhase<Tv, Td>::mapPhase(Tv v) const { return v*Tv(2)*mAmp; }

template<class Tv, class Td>
inline Tv AccumPhase<Tv, Td>::nextPhaseUsing(Tv inc){
	mPhase = scl::wrap(mPhase, mAmp,-mAmp); // guarantees that result is in [-A, A)
	Tv r = mPhase;
	mPhase += inc;
	return r;
}

template<class Tv, class Td> inline Tv AccumPhase<Tv, Td>::nextPhase(){
	return nextPhaseUsing(mInc);
}

template<class Tv, class Td> inline Tv AccumPhase<Tv, Td>::nextPhase(Tv frqMod){
	return nextPhaseUsing(mInc + mapFreq(frqMod));
}

template<class Tv, class Td> inline void AccumPhase<Tv, Td>::freq(Tv v){ mInc = mapFreq(v); }
template<class Tv, class Td> inline void AccumPhase<Tv, Td>::period(Tv v){ freq(Tv(1)/v); }
template<class Tv, class Td> inline void AccumPhase<Tv, Td>::phase(Tv v){ mPhase = mapPhase(v); }
template<class Tv, class Td> inline void AccumPhase<Tv, Td>::phaseAdd(Tv v){ mPhase += mapPhase(v); }
template<class Tv, class Td> inline void AccumPhase<Tv, Td>::amp(Tv v){
	Tv f = freq();
	mAmp=v;
	recache();
	freq(f);
}

template<class Tv, class Td> inline Tv AccumPhase<Tv, Td>::freq() const { return mInc/mFreqToInc; }
template<class Tv, class Td> inline Tv AccumPhase<Tv, Td>::period() const { return Tv(1)/freq(); }
template<class Tv, class Td> inline Tv AccumPhase<Tv, Td>::phase() const { return mPhase/(Tv(2)*mAmp); }
template<class Tv, class Td> inline Tv AccumPhase<Tv, Td>::amp() const { return mAmp; }

template<class Tv, class Td> void AccumPhase<Tv, Td>::onDomainChange(double r){
	Tv f=freq();
	recache();
	freq(f);
}
template<class Tv, class Td> void AccumPhase<Tv, Td>::recache(){
	mFreqToInc = Tv(Td::ups() * Tv(2)*mAmp);
}

template<class Tv, class Td> void AccumPhase<Tv, Td>::print(FILE * fp, const char * append){
	fprintf(fp, "%f %f %f%s", freq(), phase(), mInc, append);
}


//---- CSine

template<class Tv, class Td> CSine<Tv, Td>::CSine(Tv f, Tv a, Tv dcy60, Tv p)
	: val(a, 0), mAmp(a), mFreq(f), mDcy60(dcy60)
{
	Td::refreshDomain();
	this->phase(p);
}

template<class Tv, class Td> void CSine<Tv, Td>::amp(Tv v){
	if(scl::abs(mAmp) > Tv(0.000001)){
		Tv factor = v / mAmp;
		val *= factor;
	} else {
		val(v, Tv(0));
	}
	mAmp = v;
}

template<class Tv, class Td> void CSine<Tv, Td>::decay(Tv v){
	mDcy60 = v;
	mDcy = v > Tv(0) ? Tv(scl::t60(v * Td::spu())) : Tv(1);
	freq(mFreq);
}

template<class Tv, class Td> void CSine<Tv, Td>::freq(Tv v){
	mFreq = v;
	Tv phaseInc = v * Td::ups() * Tv(M_2PI);
	mInc.fromPolar(mDcy, phaseInc);
	//printf("%f %f %f %f\n", phaseInc, mDcy, c1, s1);
}

template<class Tv, class Td> void CSine<Tv, Td>::phase(Tv v){
	// set phase without changing current magnitude
	val.fromPolar(val.norm(), v*Tv(M_2PI));
}

template<class Tv, class Td> void CSine<Tv, Td>::reset(){ val(mAmp, Tv(0)); }

template<class Tv, class Td> void CSine<Tv, Td>::set(Tv frq, Tv phase, Tv amp, Tv dcy60){
	mFreq = frq;
	decay(dcy60);
	this->amp(amp);
	this->phase(phase);
}

template<class Tv, class Td> inline Complex<Tv> CSine<Tv, Td>::operator()(){
	complex c = val;
	val *= mInc;
	return c;
}

template<class Tv, class Td> void CSine<Tv, Td>::onDomainChange(double r){
	decay(mDcy60); // this sets frequency as well
}




//---- LFO
#define TLFO LFO<St,Td>
template<class St, class Td> TLFO::LFO(): Base(){ mod(0.5); }
template<class St, class Td> TLFO::LFO(float f, float p, float m): Base(f, p){ mod(m); }

template<class St, class Td> inline TLFO& TLFO::set(float f, float p, float m){ this->freq(f); this->phase(p); return mod(m); }
template<class St, class Td> inline TLFO& TLFO::mod(double v){ return modI(castIntRound(v*4294967296.)); }
template<class St, class Td> inline TLFO& TLFO::modI(uint32_t v){ mMod=v; return *this; }

template<class St, class Td> inline float TLFO::line2(){
	using namespace gam::scl;
	
//	// Starts at 1
//	float r1 = rampDown(phaseI());
//	float r2 = rampDown(phaseI() + mMod);

	// Starts at -1 (better for creating attack/decay like envelopes)
	uint32_t m = scl::clip<uint32_t>(mMod, 0xffefffff, 512); // avoid div by zero
	float r1 = rampDown(phaseI() - m);
	float r2 = rampDown(phaseI());
	float p  = rampUpU(m);

	float r = (r1*r1 - r2*r2)/(4.f*p*(1.f - p));
	nextPhase();
	return r;
}

template<class St, class Td> inline float TLFO::line2U(){ return line2()*0.5f+0.5f; }

#define DEF(name, exp) template<class St, class Td> inline float TLFO::name{ float r = exp; return r; }
//DEF(cos(),		tri(); r *= 0.5f * r*r - 1.5f)
//DEF(cos(),		up(); r=scl::abs(r*r) )//r = -1.f - scl::pow2(2.f*r)*(scl::abs(r)-1.5f) )
DEF(cos(),		up(); r = -1.f - r*r*(4.f*scl::abs(r)-6.f) )
DEF(down(),		scl::rampDown(nextPhase()))
DEF(even3(),	up(); static const float c=-1.50f*sqrtf(3.f); r *= (1.f-r*r)*c;)
DEF(even5(),	up(); static const float c=-1.25f*::powf(5.f,0.25f); r *= (1.f-scl::pow4(r))*c;)
DEF(imp(),		scl::pulseU(nextPhase(), this->freqI()) )
DEF(para(),		paraU()*1.5f - 0.5f)
DEF(pulse(),	scl::pulse(nextPhase(), mMod))
DEF(sinPara(),	scl::sinPara(nextPhase()))
DEF(stair(),	scl::stair(nextPhase(), mMod))
DEF(sqr(),		scl::square(nextPhase()))
DEF(tri(),		scl::triangle(nextPhase()))
DEF(up(),		scl::rampUp(nextPhase()))
DEF(up2(),		scl::rampUp2(nextPhase(), mMod))
DEF(cosU(),		tri(); r = scl::mapSinSU(r))
DEF(downU(),	scl::rampDownU(nextPhase()))
DEF(hann(),		tri(); r = r * (0.25f * r*r - 0.75f) + 0.5f)
DEF(paraU(),	up(); r*=r;)
DEF(pulseU(),	scl::pulseU(nextPhase(), mMod))
DEF(sqrU(),		scl::squareU(nextPhase()))
DEF(stairU(),	scl::stairU(nextPhase(), mMod))
DEF(triU(),		scl::triangleU(nextPhase()))
DEF(upU(),		scl::rampUpU(nextPhase()))
DEF(up2U(),		scl::rampUp2U(nextPhase(), mMod))
DEF(patU(),		scl::rampUpU(nextPhase() & mMod))

DEF(patU(uint32_t mul), scl::rampUpU((nextPhase() & mMod) * mul))

DEF(sineT9(),	up(); r = scl::sinT9(r * M_PI))
DEF(sineP9(),	up(); r = scl::sinP9(r))

#undef DEF

#undef TLFO


//---- Buzz

template<class Tv, class Td> Buzz<Tv,Td>::Buzz(Tv f, Tv p, Tv harmonics)
:	Base(f, p), mAmp(0), mPrev(Tv(0))
{
	onDomainChange(1);
	this->harmonics(harmonics);
}

template<class Tv, class Td> inline void Buzz<Tv,Td>::harmonics(Tv v){
	mN = mNDesired = scl::floor(v);
	setAmp();
	mNFrac = v - mN;
}

template<class Tv, class Td> inline void Buzz<Tv,Td>::harmonicsMax(){ harmonics(maxHarmonics()); }

template<class Tv, class Td> inline void Buzz<Tv,Td>::antialias(){
	float maxN = scl::floor(maxHarmonics());
	mN = mNDesired > maxN ? maxN : mNDesired;
	setAmp();
}

template<class Tv, class Td> inline Tv Buzz<Tv,Td>::maxHarmonics(){ return mSPU_2 / this->freq(); }

template<class Tv, class Td> inline void Buzz<Tv,Td>::setAmp(){
	// Normally, the amplitude is 1/(2N), but we will linearly interpolate
	// based on fractional harmonics to avoid sudden changes in amplitude to
	// the lower harmonics which is very noticeable.
	mAmp = (mN != Tv(0)) ? (Tv(0.5) / (mN+mNFrac)) : 0;
	//mAmp = (mN != Tv(0)) ? (Tv(0.5) / (mN)) : 0;
}

#define EPS 0.000001
template<class Tv, class Td> inline Tv Buzz<Tv, Td>::operator()(){
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

template<class Tv, class Td> inline Tv Buzz<Tv,Td>::odd(){
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

template<class Tv, class Td>
inline Tv Buzz<Tv,Td>::saw(Tv i){ return mPrev=(*this)()*0.125 + i*mPrev; }
template<class Tv, class Td>
inline Tv Buzz<Tv,Td>::square(Tv i){ return mPrev=odd()*0.125 + i*mPrev; }

template<class Tv, class Td> void Buzz<Tv,Td>::onDomainChange(double r){
	Base::onDomainChange(r);
	mSPU_2 = Tv(Td::spu() * 0.5);
}



//---- DSF

template<class Tv, class Td> DSF<Tv,Td>::DSF(Tv frq, Tv freqRatioA, Tv ampRatioA, Tv harmonicsA)
	: Base(frq)
{
	freq(frq);
	freqRatio(freqRatioA);
	harmonics(harmonicsA);
	ampRatio(ampRatioA);

	mBeta = 0.f;
}

template<class Tv, class Td> inline void DSF<Tv,Td>::freq(Tv v){
	Base::freq(v);
	updateBetaInc();	
}

template<class Tv, class Td> inline void DSF<Tv,Td>::freqRatio(Tv v){
	mFreqRatio = v;
	updateBetaInc();
}

template<class Tv, class Td> inline void DSF<Tv,Td>::ampRatio(Tv v){
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

template<class Tv, class Td> inline void DSF<Tv,Td>::harmonics(Tv v){
	if(v != mN){
		mN = mNDesired = v;
		updateAPow();
	}
}

template<class Tv, class Td> inline void DSF<Tv,Td>::harmonicsMax(){ harmonics(maxHarmonics()); }

template<class Tv, class Td> inline void DSF<Tv,Td>::antialias(){
	Tv maxN = maxHarmonics();
	if(mNDesired > maxN)	mN = maxN;
	else					mN = mNDesired;
	updateAPow();
}

template<class Tv, class Td> inline Tv DSF<Tv,Td>::ampRatio(){ return mA; }
template<class Tv, class Td> inline Tv DSF<Tv,Td>::freqRatio(){ return mFreqRatio; }
template<class Tv, class Td> inline Tv DSF<Tv,Td>::harmonics(){ return mN; }

template<class Tv, class Td> inline Tv DSF<Tv,Td>::maxHarmonics(){
	return scl::floor((Tv(this->spu()) * Tv(0.5)/this->freq() - Tv(1))/freqRatio() + Tv(1));
}

template<class Tv, class Td> inline void DSF<Tv,Td>::updateAPow(){ mAPow = ::pow(mA, mN); }
template<class Tv, class Td> inline void DSF<Tv,Td>::updateBetaInc(){ mBetaInc = this->mInc * mFreqRatio; }

// Generalized DSF formula:
// sum{k=0, N}( a^k sin(T + k B) )
//		=  sin(T) - a sin(T - B) - a^(N+1) [sin(T + (N+1) B) - a sin(T + N B))]
//			/ 1 + a^2 - 2a cos(B)

#define SIN scl::sinT7
#define COS scl::cosT8
//#define SIN sin
//#define COS cos
template<class Tv, class Td> inline Tv DSF<Tv,Td>::operator()(){
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

template<class Tv, class Td> void DSF<Tv,Td>::onDomainChange(double r){
	Base::onDomainChange(r);
	freq(Base::freq());
	harmonics(mNDesired);
}

} // gam::
#endif
