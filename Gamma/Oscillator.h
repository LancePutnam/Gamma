#ifndef GAMMA_OSCILLATOR_H_INC
#define GAMMA_OSCILLATOR_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include "Gamma/gen.h"
#include "Gamma/scl.h"
#include "Gamma/tbl.h"
#include "Gamma/Strategy.h"
#include "Gamma/Domain.h"
#include "Gamma/Types.h"

namespace gam{

/// Periodic waveforms to be used as sound or modulation sources

/// \defgroup Oscillator


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
/// \ingroup Oscillator
template <class Sp = phsInc::Loop, class Td = DomainObserver>
class Accum : public Td {
public:

	/// \param[in] frq		Frequency
	/// \param[in] phs		Phase in [0, 1)
	Accum(float frq=440, float phs=0);


	void freq(float v);				///< Set frequency
	void freqI(uint32_t v);			///< Set fixed-point frequency
	void freqAdd(float v);			///< Add value to frequency for 1 sample
	void freqMul(float v);			///< Multiply frequency by value for 1 sample
	void phase(float v);			///< Set phase from [0, 1) of one period
	void phaseMax();				///< Set phase to maximum value
	void phaseAdd(float v);			///< Add value to phase [0, 1)
	void period(float v);			///< Set period length

	void reset(){ mPhaseI=0; mSp.reset(); }	///< Reset phase accumulator
	void finish(){ phaseMax(); }	///< Set phase to end (maximum value)

	Sp& phsInc(){ return mSp; }		///< Get phase increment strategy

	bool done() const;				///< Returns true if done cycling
	bool cycled() const;			///< Returns whether phase cycled on last iteration

	float freq() const;				///< Get frequency
	uint32_t freqI() const;			///< Get fixed-point frequency
	float freqUnit() const;			///< Get frequency in [0, 1)
	float period() const;			///< Get period
	float phase() const;			///< Get phase in [0, 1)
	uint32_t phaseI() const;		///< Get fixed-point phase

	/// Iterates accumulator; \returns true on phase wrap, false otherwise
	bool operator()();

	uint32_t nextPhase();			///< Increment phase and return updated phase
	uint32_t nextPhase(float freqOffset);
	uint32_t cycles();				///< Get 1 to 0 transitions of all accumulator bits
	bool cycle();
	bool once();

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

	void onDomainChange(double r);

//protected:
private:
	float mFreq;		// Current frequency
	double mFreqToInt;
	uint32_t mPhaseI;	// Current fixed-point phase in [0, 2^32)
	uint32_t mFreqI;	// Current fixed-point frequency
	Sp mSp;

	uint32_t mapFreq(float v) const;
};

#define ACCUM_INHERIT\
	using Accum<Sp,Td>::phaseI;\
	using Accum<Sp,Td>::freqI;\
	using Accum<Sp,Td>::nextPhase;


/// Linear sweep in interval [0,1)

/// \tparam Sp	Phase increment strategy (e.g., phsInc::Loop, phsInc::Oneshot)
/// \tparam Td	Domain type
///\ingroup Oscillator
template <class Sp = phsInc::Loop, class Td = DomainObserver>
class Sweep : public Accum<Sp, Td> {
public:
	/// \param[in] frq		Frequency
	/// \param[in] phs		Phase in [0,1)
	Sweep(float frq=440, float phs=0)
	:	Accum<Sp, Td>(frq, phs){}

	float operator()(){
		float r = this->phase();
		this->nextPhase();
		return r;
	}
};


/// Floating-point phase accumulator with output in [-A, A)

/// \tparam Tv	Value (sample) type
/// \tparam Td	Domain type
///\ingroup Oscillator
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
	void freqAdd(Tv v);		///< Add value to frequency for 1 sample
	void freqMul(Tv v);		///< Multiply frequency by value for 1 sample
	void period(Tv v);		///< Set period length
	void phase(Tv v);		///< Set phase from [0, 1) of one period
	void phaseAdd(Tv v);	///< Add value to unit phase
	void amp(Tv v);			///< Set amplitude
	
	Tv freq() const;		///< Get frequency
	Tv freqUnit() const;	///< Get frequency in [0, 1)
	Tv period() const;		///< Get period
	Tv phase() const ;		///< Get normalized phase in [0, 1)
	Tv amp() const;			///< Get amplitude
	
	void onDomainChange(double r);
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
/// \ingroup Oscillator Envelope
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

	/// Default constructor does not allocate table memory
	Osc()
	:	ArrayPow2<Tv>(defaultArray<Tv>(), 1)
	{}

	/// Constructor that allocates an internal table

	/// \param[in]	frq			Frequency
	/// \param[in]	phs			Phase in [0, 1)
	/// \param[in]	size		Size of table (actual number is power of 2 ceiling)
	Osc(float frq, float phs=0, uint32_t size=512)
	:	Accum<Sp,Td>(frq, phs), ArrayPow2<Tv>(size, Tv())
	{}

	/// Constructor that references an external table

	/// \param[in]	frq			Frequency
	/// \param[in]	phs			Phase in [0, 1)
	/// \param[in]	src			A table to use as a reference
	Osc(float frq, float phs, ArrayPow2<Tv>& src)
	:	Accum<Sp,Td>(frq, phs), ArrayPow2<Tv>(src.elems(), src.size())
	{}


	/// Generate next sample
	Tv operator()(){
		Tv r = val();
		this->nextPhase();
		return r;
	}

	/// Get current value
	Tv val() const { return atPhaseI(this->phaseI()); }

	/// Get table value at fixed-point phase
	Tv atPhaseI(uint32_t v) const { return mIpol(table(), v); }
	
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

	/// Get reference to table
	ArrayPow2<Tv>& table(){ return *this; }
	const ArrayPow2<Tv>& table() const { return *this; }

protected:
	Si<Tv> mIpol;
};



/// Complex sinusoid oscillator

/// This oscillator outputs a (decaying) complex sinusoid whose real and
/// imaginary components correspond to a cosine and sine wave, respectively.
/// The complex sinusoid is generated by recursively multiplying two complex
/// numbers and thus requires only four multiplications and two additions per
/// iteration. One disadvantage of this oscillator is that it is expensive to
/// change its frequency. Another is that the amplitude may drift over time.
/// This, however, can be mostly resolved by using double precision.
///
/// Reference: Mathews, M., Smith, J. 2003. "Methods for synthesizing very high
/// Q parametrically well behaved two pole filters."
///
/// \tparam Tv	Value (sample) type
/// \tparam Td	Domain type
/// \ingroup Oscillator 
///  \sa Osc, TableSine, CSine, LFO, Sine, SineR
template<class Tv = gam::real, class Td = DomainObserver>
class CSine : public Td{
public:

	typedef Complex<Tv> complex;

	/// \param[in] frq		Frequency
	/// \param[in] amp		Amplitude
	/// \param[in] dcy		-60 dB decay length (negative for no decay)
	/// \param[in] phs		Phase in [0, 1)
	CSine(Tv frq=440, Tv amp=1, Tv dcy=-1, Tv phs=0);


	complex val;				///< Current value

	complex operator()();		///< Generate next sample

	void amp(Tv v);				///< Set amplitude
	void decay(Tv v);			///< Set -60 dB decay length (negative for no decay)
	void freq(Tv v);			///< Set frequency
	void freq(complex v){ mInc=v; }
	void phase(Tv v);			///< Set unit phase, in [0,1]
	void reset();				///< Reset amplitude and set phase to 0
	void set(Tv frq, Tv phs, Tv amp, Tv dcy);

	Tv amp() const {return mAmp;}		///< Get amplitude
	Tv decay() const {return mDcy60;}	///< Get decay length
	Tv freq() const {return mFreq;}		///< Get frequency

	void onDomainChange(double r);

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
/// \ingroup Oscillator 
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
/// \ingroup Oscillator 
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

	/// Set period
	void period(Tv v){ Base::freq(Td::ups()/v); }

	/// Set all control parameters
	void set(Tv frq, Tv amp, Tv phs=0){ Base::set(frq*Td::ups(), phs, amp); }


	void onDomainChange(double ratio){ Base::freq(Base::freq()/ratio); }

private:
	typedef gen::RSin<Tv> Base;
};



/// Multiple SineRs

/// This is a dynamically-sized bank of SineR objects.
///
/// \tparam Tv	Value (sample) type
/// \tparam Td	Domain type
/// \ingroup Oscillator 
/// \sa SineR (synthesizes a single sine)
template <class Tv = double, class Td = DomainObserver>
class SineRs : public Array<SineR<Tv, Domain1> >, Td{
public:

	SineRs(){}

	/// \param[in]	num		Number of resonators
	SineRs(unsigned num): Base(num){ onDomainChange(1); }

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


	void onDomainChange(double ratio){
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
/// \ingroup Oscillator 
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

	/// Set decay length
	void decay(Tv v){ Base::decay(decayFactor(v)); }

	/// Set frequency
	void freq(Tv v){ Base::freq(v*Td::ups()); }

	/// Set period
	void period(Tv v){ Base::freq(Td::ups()/v); }

	/// Set all control parameters

	/// \param[in]	frq		Frequency
	/// \param[in]	amp		Amplitude
	/// \param[in]	dcy		T60 decay length (negative == no decay)
	/// \param[in]	phs		Phase in [0, 1)
	void set(Tv frq, Tv amp, Tv dcy, Tv phs=0){
		Base::set(
			frq * Td::ups(),
			phs,
			decayFactor(dcy),
			amp
		);
	}

	void onDomainChange(double ratio){
		Base::freq(Base::freq()/ratio);
		//printf("%g\n", Base::decay());
		//printf("%g %g %g\n", Base::decay(), ratio, ::pow(Base::decay(), 1./ratio));
		// double radius60(double dcy, double ups){ return ::exp(M_LN001/dcy * ups); }
		//		same as (0.001)^(ups/dcy)
		Base::decay(::pow(Base::decay(), 1./ratio));
	}

private:
	typedef gen::RSin2<Tv> Base;
	Tv decayFactor(Tv length){
		return length > Tv(0) ? Tv(scl::radius60(length, Td::ups())) : Tv(1);
	}
};



/// Multiple SineDs

/// This is a dynamically-sized bank of SineD objects.
///
/// \tparam Tv	Value (sample) type
/// \tparam Td	Domain type
/// \ingroup Oscillator
/// \sa SineD
template <class Tv = double, class Td = DomainObserver>
class SineDs : public Array<SineD<Tv, Domain1> >, Td{
public:

	SineDs(){}

	/// \param[in]	num		Number of resonators
	SineDs(unsigned num): Base(num){
		onDomainChange(1);
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


	void onDomainChange(double ratio){
		for(unsigned i=0; i<this->size(); ++i){
			(*this)[i].onDomainChange(ratio);
		}
	}
private:
	typedef Array<SineD<Tv, Domain1> > Base;
};



/// Swept sinusoid with Gaussian envelope

/// This generates a sinusoid with a linear frequency sweep and Gaussian
/// envelope. Only two complex multiplies are required per sample.
/// Single-precision should be used only for very short envelopes due to
/// accumulation error.
///
/// \tparam Tv	Value (sample) type
/// \tparam Td	Domain type
/// \ingroup Oscillator
template<class Tv = double, class Td = DomainObserver>
class Chirplet : public Td{
public:

	typedef Complex<Tv> complex;

	/// \param[in] frq1		Start frequency
	/// \param[in] frq2		End frequency
	/// \param[in] amp		Amplitude
	/// \param[in] len		Length
	/// \param[in] phs		Phase in [0, 1)
	Chirplet(Tv frq1=440, Tv frq2=880, Tv amp=1, Tv len=1, Tv phs=0){
		freq(frq1, frq2);
		this->amp(amp);
		length(len);
	}


	/// Set start and end frequencies
	Chirplet& freq(double start, double end){
		mFreq1 = start;
		mFreq2 = end;
		mRGauss.mul1.arg(freqToRad(start));
		mRGauss.mul2.arg(freqToRad(end-start)/mRGauss.length());
		//printf("(%g, %g), (%g, %g)\n", mRGauss.mul2.mag(), mRGauss.mul2.arg()/2/M_PI, mRGauss.mul1.mag(), mRGauss.mul1.arg()/2/M_PI);
		return *this;
	}

	/// Set frequency
	Chirplet& freq(double v){
		mFreq1 = v;
		mFreq2 = v;
		mRGauss.mul1.arg(freqToRad(v));
		mRGauss.mul2.arg(0);
		return *this;
	}

	/// Set envelope length
	Chirplet& length(double v, const complex& offset=complex(0.01)){
		mLength = v;
		mRGauss.set(v*Td::spu(), offset);
		return *this;
	}

	/// Set amplitude
	Chirplet& amp(double v){
		mRGauss.mul1.mag(v);
		return *this;
	}

	/// Get next sample
	complex operator()(){
		return mRGauss();
	}

	/// Reset envelope
	Chirplet& reset(){ return length(mLength); }

	/// Returns true if envelope done
	bool done(float thresh=0.001) const { return mRGauss.done(); }

	/// Get length
	Tv length() const { return mLength; }

	void onDomainChange(double r){
		//Tv A = mRGauss.mul1.mag();
		length(mLength);
		freq(mFreq1, mFreq2);
		//amp(A);
	}

protected:
	gen::RGauss<complex> mRGauss;
	Tv mFreq1, mFreq2, mLength;
	double freqToRad(double v){ return v * 2*M_PI*Td::ups(); }
};



/// Low-frequency oscillator

/// This object generates various waveform types by mapping the output of a 
/// an accumulator through mathematical functions. The resulting waveforms are
/// non-band-limited.
///
/// \tparam Sp	Phase increment strategy (e.g., phsInc::Loop, phsInc::Oneshot)
/// \tparam Td	Domain type
/// \ingroup Oscillator 
/// \sa Osc, TableSine, CSine, Sine, SineR
template <class Sp = phsInc::Loop, class Td = DomainObserver>
class LFO : public Accum<Sp,Td>{
public:

	LFO();
	
	/// \param[in] frq		Frequency
	/// \param[in] phase	Phase in [0, 1)
	/// \param[in] mod		Modifier amount in [0, 1)
	LFO(double frq, double phase=0, double mod=0.5);


	/// Set frequency, phase and modifier amount
	LFO& set(float f, float p, float m);

	LFO& mod(double n);		///< Set modifier from unit value
	LFO& modI(uint32_t v);	///< Set modifier from integer

	/// Get modifier value
	uint32_t modI() const { return mMod; }
	double mod() const { return mMod / 4294967296.; }

	float cos();		///< Cosine-like wave based on 3rd order polynomial
	float down();		///< Downward ramp (1 to -1); S1 Clausen function
	float even3();		///< Even harmonic sine-like wave (3rd order); S3 Clausen function
	float even5();		///< Even harmonic sine-like wave (5th order)
	float imp();		///< Impulse train with aliasing reduction
	float line2();		///< 2-segment line; 'mod' changes wave from down to tri to up
	float para();		///< Parabolic wave; C2 Clausen function
	float pulse();		///< Pulse (up + down). 'mod' controls pulse width
	float pulseRange(); ///< Pulse (up + down). 'mod' controls pulse width. amplitude doesn't change with mod.
	float sinPara();	///< Sine-like wave constructed from parabolas; integral of triangle
	float stair();		///< Stair (square + square). 'mod' controls pulse width
	float sqr();		///< Square (-1 to 1)
	float tri();		///< Triangle (starts at 1 goes down to -1 then up to 1)
	float up();			///< Upward ramp
	float up2();		///< Dual upward ramp (up + up). 'mod' controls pulse width.

	float S1();			///< S1 Clausen function; sum_k sin(kt)/k^1
	float C2();			///< C2 Clausen function; sum_k cos(kt)/k^2
	float S3();			///< S3 Clausen function; sum_k sin(kt)/k^3
	float C4();			///< C4 Clausen function; sum_k cos(kt)/k^4
	float S5();			///< S5 Clausen function; sum_k sin(kt)/k^5

	float cosU();		///< Unipolar cosine-like wave based on 3rd order polynomial
	float downU();		///< Unipolar downward ramp
	float hann();		///< Hann-like window based on 3rd order polynomial
	float impU();		///< Unipolar impulse train
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



/// Differenced wave oscillator

/// This oscillator uses a difference between two waveforms to reduce aliasing.
/// Computation time is higher than that of an LFO unit generator, however,
/// the output exhibits far less aliasing at high frequencies.
template <class Sp = phsInc::Loop, class Td = DomainObserver>
class DWO : public Accum<Sp,Td>{
public:

	DWO();
	
	/// \param[in] frq		Frequency
	/// \param[in] phase	Phase in [0, 1)
	/// \param[in] mod		Modifier amount in [0, 1)
	DWO(float frq, float phase=0, float mod=0.5);


	DWO& mod(double n);		///< Set modifier from unit value
	DWO& modI(uint32_t v);	///< Set modifier from integer

	/// Get modifier value
	uint32_t modI() const { return mMod; }
	double mod() const { return mMod / 4294967296.; }

	/// Set frequency
	void freq(float v);
	void period(float v){ freq(1.f/v); }

	float up();				///< Upward saw
	float down();			///< Downward saw
	float sqr();			///< Square
	float para();			///< Parabolic
	float tri();			///< Triangle
	float pulse();			///< Pulse

	void onDomainChange(double r);

private:
	typedef Accum<Sp,Td> Base;
	uint32_t mMod;			// Modifier parameter
	float mGain;
	//float mPrev;
	//float diff(float v);
};



/// Sum of cosine waves

/// This produces a finite sum of harmonic, equi-amplitude cosine waves that 
/// approach the shape of a periodic impulse train. The algorithm uses a 
/// closed-form solution to the sum thus its computational complexity is always
/// O(1) (constant) regardless of the number of harmonics.
///
/// \tparam Tv	Value (sample) type
/// \tparam Td	Domain type
/// \ingroup Oscillator 
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
	void normalize(bool v);		///< Whether to normalize amplitude

	Tv operator()();			///< Returns next sample of all harmonic impulse
	Tv odd();					///< Returns next sample of odd harmonic impulse
	Tv saw(Tv intg=0.999);		///< Returns next sample of saw waveform
	Tv square(Tv intg=0.999);	///< Returns next sample of square waveform
	
	Tv maxHarmonics() const;	///< Get number of harmonics below Nyquist based on current settings

	void onDomainChange(double r);

protected:
	Tv mAmp;			// amplitude normalization factor
	Tv mN;				// actual number of harmonics
	Tv mNDesired;		// desired number of harmonics
	Tv mNFrac;		
	Tv mSPU_2;			// cached local
	Tv mPrev;			// previous output for integration
	bool mNormalize;
	void setAmp();
private: typedef AccumPhase<Tv,Td> Base;
};



/// Band-limited impulse train

/// This produces a Fourier representation of an impulse train where the number
/// of harmonics is adjusted automatically to prevent aliasing.
///
/// \tparam Tv	Value (sample) type
/// \tparam Td	Domain type
/// \ingroup Oscillator
/// \sa Buzz
template <class Tv = gam::real, class Td = DomainObserver>
struct Impulse : public Buzz<Tv,Td>{
public:

	/// \param[in] frq		Frequency
	/// \param[in] phs		Phase, in [0, 1)
	Impulse(Tv frq=440, Tv phs=0): Base(frq, phs){ onDomainChange(1); }

	/// Set frequency
	void freq(Tv v){ Base::freq(v); Base::harmonicsMax(); }

	void onDomainChange(double r){
		Base::onDomainChange(r);
		freq(AccumPhase<Tv,Td>::freq());
	}

	using Buzz<Tv,Td>::freq; // needed for getter

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
/// \ingroup Oscillator 
template <class Tv = gam::real, class Td = DomainObserver>
struct Saw : public Impulse<Tv,Td> {

	/// \param[in] frq		Frequency
	/// \param[in] phs		Phase, in [0, 1)
	Saw(Tv frq=440, Tv phs=0): Impulse<Tv, Td>(frq, phs){}

	/// Generate next sample
	
	/// \param[in] itg		Leaky integration factor
	///
	Tv operator()(Tv itg=0.999){ return Impulse<Tv,Td>::saw(itg); }
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
/// \ingroup Oscillator 
template <class Tv = gam::real, class Td = DomainObserver>
struct Square : public Impulse<Tv,Td> {

	/// \param[in] frq		Frequency
	/// \param[in] phs		Phase, in [0, 1)
	Square(Tv frq=440, Tv phs=0) : Impulse<Tv,Td>(frq, phs){}

	/// Generate next sample
	
	/// \param[in] itg		Leaky integration factor
	///
	Tv operator()(Tv itg=0.999){ return Impulse<Tv,Td>::square(itg); }
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
/// \ingroup Oscillator 
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

	Tv ampRatio() const;		///< Get amplitude ratio
	Tv freqRatio() const;		///< Get frequency ratio
	Tv harmonics() const;		///< Get current number of harmonics
	Tv maxHarmonics() const;	///< Get maximum number of harmonics for current settings
	
	void onDomainChange(double r);

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



/// Upsamples and interpolates a signal

/// This interpolates between the samples of a lower-rate signal. Use cases
/// include a sample-and-hold and a low-frequency noise generator.
///
/// \tparam Gen	Signal generator
/// \tparam Si	Sequence interpolation strategy
/// \tparam Sp	Phase increment strategy
/// \tparam Td	Domain type
/// \ingroup Oscillator
template <
	class Gen = gen::Default<>,
	template <typename> class Si = iplSeq::Linear,
	class Sp = phsInc::Loop,
	class Td = DomainObserver
>
class Upsample : public Accum<Sp,Td>{
public:

	typedef typename Gen::value_type value_type;

	/// \param[in] frq		Frequency of lower-rate signal
	Upsample(float frq=440)
	:	Accum<Sp,Td>(frq)
	{
		this->finish();
		mIpl.push(0);
	}

	// Upsample an external generator
	template <class SampleGen>
	value_type operator()(SampleGen& gen){
		if(this->cycle()){
			mIpl.push(gen());
		}
		return mIpl(this->phase());
	}

	// Sample and interpolate a signal
	value_type operator()(value_type in){
		return (*this)(gen::Val<value_type>(in));
	}

	// Upsample internal generator
	value_type operator()(){
		return (*this)(mGen);
	}

	/// Get internal generator
	Gen& gen(){ return mGen; }

private:
	Gen mGen;
	Si<value_type> mIpl;
};



// Implementation_______________________________________________________________

namespace{

// Convert unit floating-point to fixed-point integer.
// 32-bit float is good enough here since [0.f, 1.f) uses 29 bits.
inline uint32_t mapFI(float v){
	//return scl::unitToUInt(v);
	//return (uint32_t)(v * 4294967296.);
	return castIntRound(v * 4294967296.);
}

// Convert fixed-point integer to unit floating-point.
inline double mapIF(uint32_t v){
	return v/4294967296.;
	//return uintToUnit<float>(v); // not enough precision
}

};


//---- Accum
template<class Sp, class Td>
Accum<Sp,Td>::Accum(float f, float p)
:	mFreq(f), mFreqToInt(4294967296.), mFreqI(0)
{
	onDomainChange(1);
	phase(p);
}

template<class Sp, class Td>
inline uint32_t Accum<Sp,Td>::mapFreq(float v) const {
	//return mapFI(v * Td::ups());
	return castIntRound(v * mFreqToInt);
}

template<class Sp, class Td>
void Accum<Sp,Td>::onDomainChange(double /*r*/){ //printf("Accum: onDomainChange (%p)\n", this);
	mFreqToInt = 4294967296. / Td::spu();
	freq(mFreq);
}

template<class Sp, class Td>
inline void Accum<Sp,Td>::freq(float v){
	mFreq = v;
	mFreqI= mapFreq(v);
}

template<class Sp, class Td>
inline void Accum<Sp,Td>::freqI(uint32_t v){
	mFreqI= v;
	mFreq = mapIF(v) * Td::spu();
}

template<class Sp, class Td> inline void Accum<Sp,Td>::freqAdd(float v){ phaseAdd(v*Td::ups()); }
template<class Sp, class Td> inline void Accum<Sp,Td>::freqMul(float v){ freqAdd((v-1.f)*freq()); }

template<class Sp, class Td> inline void Accum<Sp,Td>::period(float v){ freq(1.f/v); }
template<class Sp, class Td> inline void Accum<Sp,Td>::phase(float v){ mPhaseI = mapFI(v); }
template<class Sp, class Td> inline void Accum<Sp,Td>::phaseAdd(float v){ mSp(mPhaseI, mapFI(v)); }
template<class Sp, class Td> void Accum<Sp,Td>::phaseMax(){ mPhaseI = 0xffffffff; }

template<class Sp, class Td>
inline bool Accum<Sp,Td>::done() const { return mSp.done(phaseI()); }

template<class Sp, class Td>
inline bool Accum<Sp,Td>::cycled() const {
	uint32_t prev = phaseI();
	Sp temp = mSp;
	temp(prev, ~freqI());
	return ((~phaseI() & prev) & 0x80000000) != 0;
}

template<class Sp, class Td> inline float Accum<Sp,Td>::freq() const { return mFreq; }
template<class Sp, class Td> inline uint32_t Accum<Sp,Td>::freqI() const { return mFreqI; }
template<class Sp, class Td> inline float Accum<Sp,Td>::freqUnit() const { return mapIF(mFreqI); }
template<class Sp, class Td> inline float Accum<Sp,Td>::period() const { return 1.f/freq(); }
template<class Sp, class Td> inline float Accum<Sp,Td>::phase() const { return mapIF(mPhaseI); }
template<class Sp, class Td> inline uint32_t Accum<Sp,Td>::phaseI() const { return mPhaseI; }

template<class Sp, class Td> inline uint32_t Accum<Sp,Td>::nextPhase(float frqOffset){
	uint32_t p = mPhaseI;
	mSp(mPhaseI, mFreqI + mapFreq(frqOffset)); // apply phase inc strategy
	return p;
}

template<class Sp, class Td> inline uint32_t Accum<Sp,Td>::nextPhase(){
	uint32_t p = mPhaseI;
	mSp(mPhaseI, mFreqI); // apply phase inc strategy
	return p;
}

template<class Sp, class Td> inline bool Accum<Sp,Td>::operator()(){ return cycle(); }

template<class Sp, class Td> inline bool Accum<Sp,Td>::cycle(){ return (cycles() & 0x80000000) != 0; }

//template<class Sp, class Td> inline uint32_t Accum<Sp,Td>::cycle(uint32_t mask){
//	return cycles() & mask;
//}

template<class Sp, class Td> inline uint32_t Accum<Sp,Td>::cycles(){
	uint32_t prev = phaseI();
	nextPhase();
	return ~phaseI() & prev;
}

template<class Sp, class Td> inline bool Accum<Sp,Td>::once(){
	uint32_t prev = phaseI();
	bool c = cycle();
	if(c) mPhaseI = prev;
	return c;
}

template<class Sp, class Td> inline bool Accum<Sp,Td>::seq(uint32_t pat){
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
	onDomainChange(1);
	this->phase(p);
}

template<class Tv, class Td> inline Tv AccumPhase<Tv, Td>::mapFreq(Tv v) const {
	return v*mFreqToInc;
}

template<class Tv, class Td> inline Tv AccumPhase<Tv, Td>::mapPhase(Tv v) const {
	return v*Tv(2)*mAmp;
}

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
template<class Tv, class Td> inline void AccumPhase<Tv, Td>::freqAdd(Tv v){
	nextPhaseUsing(mapFreq(v));
}
template<class Tv, class Td> inline void AccumPhase<Tv, Td>::freqMul(Tv v){
	nextPhaseUsing((v-Tv(1))*mInc);
}
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
template<class Tv, class Td> inline Tv AccumPhase<Tv, Td>::freqUnit() const { return freq()*this->ups(); }
template<class Tv, class Td> inline Tv AccumPhase<Tv, Td>::period() const { return Tv(1)/freq(); }
template<class Tv, class Td> inline Tv AccumPhase<Tv, Td>::phase() const { return mPhase/(Tv(2)*mAmp); }
template<class Tv, class Td> inline Tv AccumPhase<Tv, Td>::amp() const { return mAmp; }

template<class Tv, class Td> void AccumPhase<Tv, Td>::onDomainChange(double /*r*/){
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
	onDomainChange(1);
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

template<class Tv, class Td> void CSine<Tv, Td>::onDomainChange(double /*r*/){
	decay(mDcy60); // this sets frequency as well
}




//---- LFO
#define TLFO LFO<Sp,Td>
template<class Sp, class Td> TLFO::LFO(): Base(){ mod(0.5); }
template<class Sp, class Td> TLFO::LFO(double f, double p, double m): Base(f, p){ mod(m); }

template<class Sp, class Td> inline TLFO& TLFO::set(float f, float p, float m){
	this->freq(f);
	this->phase(p);
	return mod(m);
}
template<class Sp, class Td> inline TLFO& TLFO::mod(double v){
	return modI(castIntRound(v*4294967296.));
}
template<class Sp, class Td> inline TLFO& TLFO::modI(uint32_t v){
	mMod=v;
	return *this;
}

template<class Sp, class Td> inline float TLFO::line2(){
	uint32_t m = scl::clip<uint32_t>(mMod, 0xffefffff, 512); // avoid div by zero

	/* Starts at 1
	float r1 = scl::rampDown(phaseI());
	float r2 = scl::rampDown(phaseI() + m); //*/

	//* Starts at -1 (better for creating attack/decay like envelopes)
	float r1 = scl::rampDown(phaseI() - m);
	float r2 = scl::rampDown(phaseI()); //*/

	float p  = punUF(Expo2<float>() | (m >> 9)) - 2.f; // [0, 2);
	float r = (r1*r1 - r2*r2)/(p*(2.f - p));
	nextPhase();
	return r;
}

template<class Sp, class Td> inline float TLFO::line2U(){ return line2()*0.5f+0.5f; }

#define DEF(name, exp) template<class Sp, class Td> inline float TLFO::name{ float r = exp; return r; }
//DEF(cos(),		tri(); r *= 0.5f * r*r - 1.5f)
//DEF(cos(),		up(); r=scl::abs(r*r) )//r = -1.f - scl::pow2(2.f*r)*(scl::abs(r)-1.5f) )
DEF(cos(),		up(); r = -1.f - r*r*(4.f*scl::abs(r)-6.f) )
DEF(down(),		scl::rampDown(nextPhase()))
DEF(even3(),	up(); float c=-2.598076211353315;/*-1.50*pow(3,0.50)*/ r *= (1.f-r*r)*c;)
DEF(even5(),	up(); float c=-1.869185976526527;/*-1.25*pow(5,0.25)*/ r *= (1.f-scl::pow4(r))*c;)
DEF(imp(),		scl::pulse(nextPhase(), this->freqI()) )
DEF(para(),		paraU()*1.5f - 0.5f)
DEF(pulse(),	scl::pulse(nextPhase(), mMod))
DEF(sinPara(),	scl::sinPara(nextPhase()))
DEF(stair(),	scl::stair(nextPhase(), mMod))
DEF(sqr(),		scl::square(nextPhase()))
DEF(tri(),		scl::triangle(nextPhase()))
DEF(up(),		scl::rampUp(nextPhase()))
DEF(up2(),		scl::rampUp2(nextPhase(), mMod))
DEF(S1(),		up(); float c=          -M_PI /  2; r = c*r)
DEF(C2(),		up(); float c= scl::pow2(M_PI)/  4; r = c*(r*r - 1.f/3))
DEF(S3(),		up(); float c= scl::pow3(M_PI)/ 12; r = c*(r*r*r - r))
DEF(C4(),		up(); float c=-scl::pow4(M_PI)/ 48; float rr=r*r; r = c*(rr*(rr - 2.f) + 7.f/15))
DEF(S5(),		up(); float c=-scl::pow5(M_PI)/240; float rr=r*r; r = c*r*(rr*(rr - 10.f/3) + 7.f/3))
DEF(cosU(),		up(); r = r*r*(-2.f*scl::abs(r)+3.f))
DEF(downU(),	scl::rampDownU(nextPhase()))
//DEF(hann(),	tri(); r = r * (0.25f * r*r - 0.75f) + 0.5f)
DEF(hann(),		up(); r = 1.f + r*r*(2.f*scl::abs(r)-3.f))
DEF(impU(),		scl::pulseU(nextPhase(), this->freqI()) )
DEF(paraU(),	up(); r*=r;)
DEF(pulseU(),	scl::pulseU(nextPhase(), mMod))
DEF(pulseRange(), scl::pulseRange(nextPhase(), mMod))
DEF(sqrU(),		scl::squareU(nextPhase()))
DEF(stairU(),	scl::stairU(nextPhase(), mMod))
DEF(triU(),		scl::triangleU(nextPhase()))
DEF(upU(),		scl::rampUpU(nextPhase()))
DEF(up2U(),		scl::rampUp2U(nextPhase(), mMod))
DEF(patU(),		scl::rampUpU(nextPhase() & mMod))

DEF(patU(uint32_t mul), scl::rampUpU((nextPhase() & mMod) * mul))

/* The input domain for these is [-1,1] corresponding to [-pi,pi], but that will
give us an upside-down sine. We therefore use a downward ramp to flip the wave
on the time axis.*/
DEF(sineT9(),	down(); r = scl::sinT9(r * M_PI))
DEF(sineP9(),	down(); r = scl::sinP9(r))

#undef DEF
#undef TLFO



template <class Sp, class Td>
DWO<Sp,Td>::DWO()
//:	mPrev(0)
{
	mod(0.5);
}

template<class Sp, class Td>
DWO<Sp,Td>::DWO(float f, float p, float m)
//:	mPrev(0)
{
	freq(f);
	this->phase(p);
	mod(m);
}

template<class Sp, class Td>
inline DWO<Sp,Td>& DWO<Sp,Td>::mod(double v){
	return modI(castIntRound(v*4294967296.));
}
template<class Sp, class Td>
inline DWO<Sp,Td>& DWO<Sp,Td>::modI(uint32_t v){
	mMod=v;
	return *this;
}

template <class Sp, class Td>
inline void DWO<Sp,Td>::freq(float v){
	Base::freq(v);
	float freq1 = this->freqUnit();
	// Very low freq will produce quantization noise
	if(freq1 < 1e-5) freq1 = 1e-5;
	mGain = -0.25/freq1;
}

/* Ideally we would use a differencing filter, however, it has a serious shortcoming. Any sudden changes in phase (or frequency) will lead a large amplitude impulse in the output which becomes worse the lower the frequency. A related problem is what to initialize the filter's previous input sample to. Instead of a filter, we use an analytic approach which subtracts two phase-shifted waveforms.
*/
namespace{
	inline float para01(uint32_t p){
		float s = scl::rampUp(p);
		return s*s;
	}
	inline float triangle02(uint32_t p){
		p = Expo4<float>() | (p >> 9); // [4, 8)
		return scl::abs(punUF(p) - 6.f);
	}
}

template <class Sp, class Td>
inline float DWO<Sp,Td>::up(){
	/*
	float s = para01(this->nextPhase());
	return diff(s);//*/
	//*
	uint32_t p = this->nextPhase();
	float s = para01(p);
	float t = para01(p + this->freqI());
	return (s - t)*mGain;//*/
}

template <class Sp, class Td>
inline float DWO<Sp,Td>::down(){
	/*float s = para01(this->nextPhase());
	return diff(-s);*/
	uint32_t p = this->nextPhase();
	float s = para01(p);
	float t = para01(p + this->freqI());
	return (t - s)*mGain;
}

template <class Sp, class Td>
inline float DWO<Sp,Td>::sqr(){
	/*
	float s = triangle02(this->nextPhase());
	return diff(s);//*/
	//*
	uint32_t p = this->nextPhase();
	float s = triangle02(p);
	float t = triangle02(p + this->freqI());
	return (t - s)*mGain;//*/
}

template <class Sp, class Td>
inline float DWO<Sp,Td>::para(){
	static const float c = (-M_PI*M_PI*M_PI/12.)*0.5;
	uint32_t p = this->nextPhase();
	float s = scl::rampUp(p);
	s = s*s*s - s;
	float t = scl::rampUp(p + this->freqI());
	t = t*t*t - t;
	return (t - s)*c*mGain;
}

template <class Sp, class Td>
inline float DWO<Sp,Td>::tri(){
	/*
	float s = scl::sinPara(this->nextPhase())*0.5f;
	return diff(s);//*/
	//*
	uint32_t p = this->nextPhase();
	float s = scl::sinPara(p);
	float t = scl::sinPara(p + this->freqI());
	return (t - s)*-0.5f*mGain;//*/
}

template <class Sp, class Td>
inline float DWO<Sp,Td>::pulse(){
	/*
	uint32_t p = this->nextPhase();
	float s1 = para01(p);
	float s2 = para01(p + mMod);
	return diff((s1 - s2)*0.5f);//*/
	//*
	uint32_t p = this->nextPhase();
	float s1 = para01(p);
	float s2 = para01(p - mMod);
	float s  = s1 - s2;
	uint32_t q = p + this->freqI();
	float t1 = para01(q);
	float t2 = para01(q - mMod);
	float t  = t1 - t2;
	return (t - s)*0.5f*mGain;//*/
}

/*template <class Sp, class Td>
inline float DWO<Sp,Td>::diff(float v){
	float res = (v - mPrev)*mGain;
	mPrev = v;
	return res;
}*/

template<class Sp, class Td>
void DWO<Sp,Td>::onDomainChange(double r){
	Base::onDomainChange(r);
	freq(Base::freq());
}


//---- Buzz

template<class Tv, class Td> Buzz<Tv,Td>::Buzz(Tv f, Tv p, Tv harmonics)
:	Base(f, p), mAmp(0), mPrev(Tv(0)), mNormalize(true)
{
	onDomainChange(1);
	this->harmonics(harmonics);
}

template<class Tv, class Td> inline void Buzz<Tv,Td>::harmonics(Tv v){
	mNDesired = v;
	mN = scl::floor(v);
	mNFrac = v - mN;
	setAmp();
}

template<class Tv, class Td> inline void Buzz<Tv,Td>::harmonicsMax(){
	harmonics(maxHarmonics());
}

template<class Tv, class Td> inline void Buzz<Tv,Td>::antialias(){
	Tv newN = maxHarmonics();
	newN = mNDesired > newN ? newN : mNDesired;

	mN = scl::floor(newN);
	mNFrac = newN - mN;
	setAmp();
}

template<class Tv, class Td> void Buzz<Tv,Td>::normalize(bool v){
	mNormalize = v;
	setAmp();
}

template<class Tv, class Td> inline Tv Buzz<Tv,Td>::maxHarmonics() const {
	return mSPU_2 / this->freq();
}

template<class Tv, class Td> inline void Buzz<Tv,Td>::setAmp(){

	if(mNormalize){
		// Normally, the amplitude is 1/(2N), but we will linearly interpolate
		// based on fractional harmonics to avoid sudden changes in amplitude to
		// the lower harmonics which is very noticeable.
		mAmp = (mN != Tv(0)) ? (Tv(0.5) / (mN+mNFrac)) : 0;
	}
	else{
		mAmp = Tv(0.5);
	}
}

#define EPS 0.000001
template<class Tv, class Td> inline Tv Buzz<Tv, Td>::operator()(){
	/*        1   / sin((N+0.5)x)    \
	   f(x) = -- |  ------------- - 1 |
	          2N  \   sin(0.5x)      /
	*/
	Tv theta = this->nextPhase();
	Tv result;
	Tv denom = scl::sinT7(theta * Tv(0.5));

	// denominator goes to zero when theta is an integer multiple of 2 pi
	if(scl::abs(denom) < Tv(EPS)){
		result = Tv(2) * mN * mAmp;
		//printf("Buzz::operator(): oops\n");
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
		//printf("Buzz::odd(): oops\n");
	}
	else result = scl::sinT7(scl::wrapPhase(n2 * theta)) / (denom * (n2 + n2frac));
	
	return result;
}
#undef EPS

template<class Tv, class Td>
inline Tv Buzz<Tv,Td>::saw(Tv b){ return mPrev=(*this)() + b*mPrev; }

template<class Tv, class Td>
inline Tv Buzz<Tv,Td>::square(Tv b){ return mPrev=odd() + b*mPrev; }

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

template<class Tv, class Td> inline void DSF<Tv,Td>::harmonicsMax(){
	harmonics(maxHarmonics());
}

template<class Tv, class Td> inline void DSF<Tv,Td>::antialias(){
	Tv maxN = maxHarmonics();
	if(mNDesired > maxN)	mN = maxN;
	else					mN = mNDesired;
	updateAPow();
}

template<class Tv, class Td> inline Tv DSF<Tv,Td>::ampRatio() const {
	return mA;
}
template<class Tv, class Td> inline Tv DSF<Tv,Td>::freqRatio() const {
	return mFreqRatio;
}
template<class Tv, class Td> inline Tv DSF<Tv,Td>::harmonics() const {
	return mN;
}

template<class Tv, class Td> inline Tv DSF<Tv,Td>::maxHarmonics() const {
	return scl::floor((Tv(this->spu()) * Tv(0.5)/this->freq() - Tv(1))/freqRatio() + Tv(1));
}

template<class Tv, class Td> inline void DSF<Tv,Td>::updateAPow(){
	mAPow = ::pow(mA, mN);
}

template<class Tv, class Td> inline void DSF<Tv,Td>::updateBetaInc(){
	mBetaInc = this->mInc * mFreqRatio;
}

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
