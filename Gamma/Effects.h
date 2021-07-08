#ifndef GAMMA_EFFECTS_H_INC
#define GAMMA_EFFECTS_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include "Gamma/scl.h" // roundN
#include "Gamma/Delay.h"
#include "Gamma/Envelope.h"
#include "Gamma/Filter.h"
#include "Gamma/Noise.h"
#include "Gamma/Oscillator.h"

namespace gam{

///\defgroup Effects


/// Amplitude modulator

/// Amplitude modulation is based on the operation (M x d + 1) x C where
/// M is the modulator, C is the carrier, and d is the modulation depth. The
/// resulting signal contains the carrier signal plus the sums and differences
/// of all frequencies of the carrier and modulator.
/// Ring modulation is given by the operation M x C and thus, unlike amplitude
/// modulation, does not pass the carrier signal.
/// \ingroup Effects
template <class Tp=real>
class AM{
public:

	/// \param[in] modDepth		modulation depth
	AM(Tp modDepth = Tp(1))
	:	mDepth(modDepth){}

	/// Set depth of amplitude modulation
	void depth(Tp v){ mDepth=v; }

	/// Modulate amplitude of carrier by modulator
	template <class Tv>
	Tv operator()(Tv car, Tv mod){
		return (mDepth*mod)*car + car;
	}

private:
	Tp mDepth;
};


/// 3 biquad filters (of floats) in parallel

/// \ingroup Effects Filter
class Biquad3{
public:
	/// Constructor
	Biquad3(float f0, float f1, float f2, float q=8, FilterType type=BAND_PASS):
		bq0(f0,q,type), bq1(f1,q,type), bq2(f2,q,type){}

	/// Set center frequencies
	void freq(float f0, float f1, float f2){ bq0.freq(f0); bq1.freq(f1); bq2.freq(f2); }
	
	/// Return filtered sample
	float operator()(float i0){ return bq0(i0) + bq1(i0) + bq2(i0); }
	
	Biquad<> bq0, bq1, bq2;
};


/// Percussive noise burst consisting of resonant-filtered white noise with a rapidly decaying amplitude
class Burst{
public:
	Burst(float frq1=20000, float frq2=4000, float dec=0.1, float res=2) : 
		freq1(frq1), freq2(frq2), fil(frq1, res, BAND_PASS), env(dec)
	{}
	
	float operator()(){
		if(env.done()) return 0.f; // faster and avoids denormals
		fil.freq(freq2 + (freq1-freq2)*env());
		return fil(src()) * env.value();
	}
	
	void operator()(float frq1, float frq2, float dec, bool rst=true){
		freq1 = frq1; freq2 = frq2; env.decay(dec); if(rst) reset();
	}
	
	void reset(){
		env.reset(amp);
		fil.zero();
	}
	
	float freq1, freq2, amp=1.f;
	NoiseWhite<RNGMulCon> src;
	Biquad<> fil;
	Decay<float> env;
};



/// Sine wave with frequency and amplitude driven by an exponentially decaying envelope.
template <class T=gam::real>
class Chirp{
public:
	/// \param[in] freq1	start frequency
	/// \param[in] freq2	end frequency
	/// \param[in] decay60	length after which envelope decays by 60 dB
	Chirp(T freq1=220, T freq2=0, T decay60=0.2):
		osc(freq1), env(decay60), freq1(freq1), freq2(freq2)
	{}
	
	/// Generate next sample
	T operator()(){
		if(env.value() > 0.0001f){
			T e = env();
			osc.freq(freq2 + (freq1 - freq2) * e);
			return osc() * e;
		}
		return 0.f;
	}
	
	void decay(T v){ env.decay(v); }
	void freq(T start, T end){ freq1=start; freq2=end; }
	
	void operator()(T frq1, T frq2, T dcy, bool doReset=false){
		freq1 = frq1; freq2 = frq2; env.decay(dcy); if(doReset) reset();
	}
	
	/// Reset envelope
	void reset(){ osc.phase(0); env.reset(); }
	
	Sine<T> osc;		///< Sine oscillator
	Decay<T> env;		///< Envelope
	T freq1;			///< Start frequency
	T freq2;			///< End frequency
};



/// Nth order Chebyshev transfer function

/// This filter applies a Chebyshev polynomial to generate the 2nd through Nth 
/// cosine harmonics of the input signal which is presumed to be a unity gain 
/// sinusoid.
///\ingroup Filter
template <unsigned N, class T=gam::real> 
class ChebyN{
public:
	T c[N];			///< Harmonic coefficients
	
	ChebyN(const T& fundAmp = T(1)){ zero(); c[0]=fundAmp; }
	
	/// Returns filtered sample
	T operator()(T i0) const { return i0*c[0] + wet(i0); }
	
	/// Returns cosine overtones of sinusoidal input
	T wet(T i0) const {
//		T d1 = i0 * T(2);	// Chebyshev polynomial of 2nd kind
		T d1 = i0;			// Chebyshev polynomial of 1st kind
		T d2 = T(1);
		T b1 = T(2) * i0;
		
		T o0 = T(0);
		for(unsigned i=1; i<N; ++i){
			T d0 = b1 * d1 - d2;
			d2 = d1;
			d1 = d0;
			o0 += d0 * c[i];
		}
		return o0;
	}
	
	/// Get number of harmonics
	unsigned size() const { return N; }

	/// Get reference to harmonic coefficient
	T& coef(int i){ return c[i]; }

	/// Set harmonic amplitudes
	template <class V>
	ChebyN& set(const V* weights){
		for(unsigned i=0; i<N; ++i) c[i] = weights[i];
		return *this;
	}

	/// Zero all harmonic amplitudes
	ChebyN& zero(){
		for(unsigned i=0; i<N; ++i) c[i] = T(0);
		return *this;
	}
};



/// Dual delay-line chorus driven by quadrature sinusoid

///
/// \ingroup Effects
template <class T=gam::real>
class Chorus{
public:
	/// \param[in] delay	Delay interval
	/// \param[in] depth	Depth of delay-line modulation
	/// \param[in] freq		Frequency of modulation
	/// \param[in] ffd		Feedforward amount
	/// \param[in] fbk		Feedback amount
	Chorus(float delay=0.0021, float depth=0.002, float freq=1, float ffd=0.9, float fbk=0.1):
		comb1(delay + depth, delay, ffd, fbk),
		comb2(delay + depth, delay, ffd, fbk),
		mod(double(freq), double(depth)),
		mDelay(delay)
	{}

	Chorus& maxDelay(float v){ comb1.maxDelay(v); comb2.maxDelay(v); return *this; }
	Chorus& delay(float v){ mDelay=v; return *this; }
	Chorus& fbk(float v){ comb1.fbk(v); comb2.fbk(v); return *this; }
	Chorus& ffd(float v){ comb1.ffd(v); comb2.ffd(v); return *this; }
	Chorus& freq(float v){ mod.freq(v); return *this; }
	Chorus& depth(float v){ mod.amp(v); return *this; }

	/// Filter sample (mono-mono)
	T operator()(const T& v){
		modulate();
		return (comb1(v) + comb2(v)) * 0.5f;
	}
	
	/// Filter sample (mono-stereo)
	void operator()(const T& in, T& o1, T& o2){
		(*this)(in,in, o1,o2);
	}

	/// Filter samples (stereo-stereo)
	template <class V>
	Vec<2,V> operator()(const Vec<2,V>& v){
		modulate();
		return Vec<2,V>(comb1(v[0]), comb2(v[1]));
	}
	
	/// Filter samples (stereo-stereo)
	void operator()(const T& i1, const T& i2, T& o1, T& o2){
		modulate();
		o1=comb1(i1); o2=comb2(i2);
	}
	
	/// Perform delay modulation step (must manually step comb filters after use!)
	void modulate(){
		comb1.delay(mDelay + mod.val.r);
		comb2.delay(mDelay + mod.val.i);
		mod();
	}

	Comb<T, ipl::Cubic> comb1, comb2;		///< Comb filters
	CSine<double> mod;						///< Modulator

private:
	float mDelay; // Delay interval
};



/// Frequency shifter

/// This effect shifts all frequencies of an input signal by a constant amount.
/// It is also known as single-sideband modulation.
///\ingroup Effects
template <class T=gam::real>
class FreqShift{
public:
	/// \param[in] shift	frequency shift amount
	FreqShift(float shift=1): mod(shift){}

	/// Frequency shift input
	T operator()(T in){
		return (hil(in) * mod()).r;
	}
	
	/// Set frequency shift amount
	FreqShift& freq(T v){ mod.freq(v); return *this; }

	CSine<T> mod;
	Hilbert<T> hil;
};



/// Saw oscillator with sweepable filter.
class MonoSynth{
public:
	MonoSynth(float freq=440, float dur=0.8, float ctf1=1000, float ctf2=100, float res=3):
		osc(freq), filter(ctf1, res), env(dur), opEnv(100), ctf1(ctf1), ctf2(ctf2)
	{}

	float operator()(){
		if(env.done()) return 0.f;

		float e = opEnv(env());
		filter.freq(ctf2 + (ctf1 - ctf2) * e);
		float smp = osc() * e;
		return filter(smp);
	}
	
	void freq(float v){ osc.freq(v); }
	
	void reset(){ env.reset(); }

	Saw<float> osc;
	Biquad<> filter;
	Decay<float> env;
	OnePole<float> opEnv;
	float ctf1, ctf2;
};



/// Equal-power 2-channel panner

///
/// \ingroup Spatial
template <class T=gam::real>
class Pan{
public:

	/// \param[in] pos	Position, in [-1, 1]
	Pan(T pos=0){ this->pos(pos); }

	/// Filter sample (mono to stereo)
	Vec<2,T> operator()(T in){
		return Vec<2,T>(in*w1, in*w2);
	}	

	/// Filter sample (mono to stereo)
	template <class V>
	void operator()(T in, V& out1, V& out2){
		out1 = in*w1; out2 = in*w2;
	}

	/// Filter sample (stereo to stereo)
	template <class V>
	void operator()(T in1, T in2, V& out1, V& out2){
		out1 = in1*w1 + in2*w2;
		out2 = in1*w2 + in2*w1;
	}

	/// Set position (constant power law)

	/// This is a constant power pan where the sum of the squares of the two
	/// channel gains is always 1. A quadratic approximation is used to avoid
	/// expensive trig function calls. The approximation is good enough for most
	/// purposes and gives the exact result at positions of -1, 0, 1.
	///
	/// \param[in] v	Position, in [-1, 1]
	void pos(T v){
		static const T c0 = 1./sqrt(2);
		static const T c1 = 0.5 - c0;
		static const T c2 =-0.5/c1;
		v = scl::clip(v, T(1), T(-1));
		w1 = c1 * v * (v + c2) + c0;
		w2 = w1 + v;

		/* // exact formula, for reference
		v = scl::clip(v, (T)1, (T)-1);
		v = (v+T(1))*M_PI_4;	// put in [0, pi/2]
		w1 = cos(v);
		w2 = sin(v);*/
	}

	/// Set position and overall gain

	/// This is identical to pos(), but includes an overall gain.
	///
	/// \param[in] v		Position, in [-1, 1]
	/// \param[in] gain		Overall gain to apply to both channels
	void pos(T v, T gain){
		pos(v); w1*=gain; w2*=gain;
	}

	/// Set position (constant gain law)

	/// This uses a constant gain law where the sum of the channel gains is
	/// always 1.
	/// \param[in] v	Position, in [-1, 1]
	void posL(T v){
		v = scl::clip(v, T(1), T(-1));
		w1 = T(0.5) - v * T(0.5);
		w2 = w1 + v;
	}

protected:
	T w1, w2; // channel weights
};



/// Plucked string source/filter

/// \ingroup Oscillator
///
class Pluck{
public:
	Pluck(double freq=440, double decay=0.99)
	:	env(0.1), fil(3000, 1, LOW_PASS), comb(1./27.5, 1./freq, 1, decay)
	{}
	
	float operator()(){ return comb(fil(noise() * env())); }
	float operator()(float in){ return comb(fil(in * env())); }
	void reset(){ env.reset(); }
	void freq(float v){ comb.freq(v); }
	
	NoiseWhite<> noise;
	Decay<> env;
	Biquad<> fil;
	Comb<> comb;
};



/// Downsamples and quantizes amplitudes

/// This effect is also known as a bitcrusher.
///\ingroup Effects
template <class T=gam::real, class Td=GAM_DEFAULT_DOMAIN>
class Quantizer : public Td{
public:
	/// \param[in] freq		Frequency of sequence quantization
	/// \param[in] stepAmt	Step amount of amplitude quantization
	Quantizer(double freq=2000, float stepAmt=0)
	:	mPeriod(1./freq){
		step(stepAmt);
	}

	/// Set freqency of sequence quantization
	void freq(double v){ period(1./v); }

	/// Set period of sequence quantization
	void period(double v){
		mPeriod = v;
		mSamples = v * this->spu();
		if(mSamples < 1.) mSamples = 1.;
	}

	/// Set amplitude quantization amount
	void step(float v){
		mStep = v;
		mDoStep = mStep > 0.f;
		if(mDoStep) mStepRec = 1.f/mStep;
	}

	/// Return next filtered sample
	T operator()(T in){
		if(++mCount >= mSamples){
			mCount -= mSamples;
			mHeld = mDoStep ? scl::roundN(in, mStep, mStepRec) : in;
		}
		return mHeld;
	}
	
	virtual void onDomainChange(double r){
		period(mPeriod);
	}

private:
	T mHeld = T(0);
	// A float can represent integers up to 16,777,216.
	// Max period at 44.1kHz is thus 38 seconds.
	float mCount=0., mSamples, mPeriod;
	float mStep, mStepRec;
	bool mDoStep;
};

} // gam::
#endif
