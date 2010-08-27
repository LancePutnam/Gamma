#ifndef GAMMA_EFFECTS_H_INC
#define GAMMA_EFFECTS_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include "Gamma/Delay.h"
#include "Gamma/Envelope.h"
#include "Gamma/Noise.h"
#include "Gamma/Oscillator.h"

#define TEM template<class T>

namespace gam{


/// Amplitude envelope extractor
template <class Tv=gam::real, class Tp=gam::real>
struct AmpEnv{
	/// @param[in] freq		Cutoff frequency of smoothing filter
	AmpEnv(Tp freq=10):lpf(freq){}
	Tv operator()(Tv i0){ return lpf(scl::abs(i0)); }	///< Returns next sample
	Tv last(){ return lpf.last(); }
	OnePole<Tv, Tp> lpf;									///< Low-pass filter
};



/// 3 biquad filters in parallel
struct Biquad3{
	/// Constructor
	Biquad3(float f0, float f1, float f2, float q=8, int type=Filter::BP):
		bq0(f0,q,type), bq1(f1,q,type), bq2(f2,q,type){}

	/// Set center frequencies
	void freq(float f0, float f1, float f2){ bq0.freq(f0); bq1.freq(f1); bq2.freq(f2); }
	
	/// Return filtered sample
	float operator()(float i0){ return bq0(i0) + bq1(i0) + bq2(i0); }
	
	Biquad<> bq0, bq1, bq2;
};



struct Burst{
	Burst(float frq1=20000, float frq2=4000, float dec=0.1, float res=2) : 
		freq1(frq1), freq2(frq2), fil(frq1, res, Filter::BP), env(dec)
	{}
	
	float operator()(){
		if(env.done()) return 0.f;
		fil.freq(freq2 + (freq1-freq2)*env());
		return fil(src()) * env.value();
	}
	
	void operator()(float frq1, float frq2, float dec, bool rst=true){
		freq1 = frq1; freq2 = frq2; env.decay(dec); if(rst) reset();
	}
	
	void reset(){ env.reset(); }
	
	float freq1, freq2;
	NoiseWhite<RNGMulLinCon> src;
	Biquad<> fil;
	Decay<float> env;
};



/// Sine wave with frequency and amplitude driven by an exponentially decaying envelope.
template <class T=gam::real>
struct Chirp{
	/// @param[in] freq1	start frequency
	/// @param[in] freq2	end frequency
	/// @param[in] decay60	units to decay by 60 dB
	Chirp(T freq1=220, T freq2=0, T decay=0.2):
		osc(freq1), env(decay), freq1(freq1), freq2(freq2)
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



/// 8th order Chebyshev transfer function

/// This filter applies a Chebyshev polynomial to generate the 2nd through 8th 
/// harmonics of the input signal which is presumed to be a unity gain cosine.
template <class T=gam::real> 
struct Cheby8{
	T c[7];									///< 2nd-8th harmonics coefficients
	
	Cheby8(T h2=0, T h3=0, T h4=0, T h5=0, T h6=0, T h7=0, T h8=0){
		c[0]=h2; c[1]=h3; c[2]=h4; c[3]=h5; c[4]=h6; c[5]=h7; c[6]=h8;
	}
	
	/// Returns filtered sample
	T operator()(T i0){ return i0 + wet(i0); }
	
	/// Returns 2nd-8th harmonics of cosine input
	T wet(T i0){
		gen::RSin<T> rsin((T)1, i0, (T)2 * i0); T o0 = (T)0;
		#define DO(h) o0 += rsin() * c[h];
			DO(0) DO(1) DO(2) DO(3) DO(4) DO(5) DO(6)
		#undef DO
		return o0;
	}
};


/// Dual delay-line chorus driven by quadrature sinusoid.
struct Chorus{
	/// @param[in] delay	Delay interval
	/// @param[in] depth	Depth of delay-line modulation
	/// @param[in] freq		Frequency of modulation
	/// @param[in] ffd		Feedforward amount
	/// @param[in] fbk		Feedback amount
	Chorus(float delay, float depth=0.002, float freq=1, float ffd=0.9, float fbk=0.1):
		comb1(delay + depth, delay, ffd, fbk),
		comb2(delay + depth, delay, ffd, fbk),
		mod((double)freq, (double)depth),
		delay(delay)
	{}
	
	/// Filter sample (mono-mono)
	float operator()(float in){
		m(); return (comb1(in) + comb2(in)) * 0.5f;
	}
	
	/// Filter sample (mono-stereo)
	void operator()(float in, float& o1, float& o2){
		m(); o1=comb1(in); o2=comb2(in);
	}
	
	/// Filter samples (stereo-stereo)
	void operator()(float i1, float i2, float& o1, float& o2){
		m(); o1=comb1(i1); o2=comb2(i2);		
	}

	Comb<float, ipl::AllPass> comb1, comb2;		///< Comb filters
	//Comb<float, ipl::Linear> comb1, comb2;	///< Comb filters
	Quadra<double> mod;	///< Modulator
	float delay;		///< Delay interval
	
private: 
	void m(){ comb1.delay(delay + mod.val.r); comb2.delay(delay + mod.val.i); mod(); } // modulate
};



/// Group of 4 comb filters.

template <template<class> class Tipol=ipl::Linear>
struct Combs4{

	/// @param[in] d1		Delay length of filter 1
	/// @param[in] d2		Delay length of filter 2
	/// @param[in] d3		Delay length of filter 3
	/// @param[in] d4		Delay length of filter 4
	/// @param[in] ffd		Feedforward amount for all filters
	/// @param[in] fbk		Feedback amount for all filters
	Combs4(float d1, float d2, float d3, float d4, float ffd, float fbk)
	: c1(d1, ffd, fbk), c2(d2, ffd, fbk), c3(d3, ffd, fbk), c4(d4, ffd, fbk){}
	
	/// Set decay length for filters.
	void decay(float v, float end = 0.001f){
		c1.decay(v, end); c2.decay(v, end); c3.decay(v, end); c4.decay(v, end);
	}
	
	/// Set delay length for filters.
	void delay(float d1, float d2, float d3, float d4){
		c1.delay(d1); c2.delay(d2); c3.delay(d3); c4.delay(d4);
	}

	float nextP(float i0){ return c1(i0) + c2(i0) + c3(i0) + c4(i0); }
	
	float nextS(float i0){ return c4(c3(c2(c1(i0)))); }
	
	Comb<float, Tipol> c1, c2, c3, c4;
};




/// Diffuser using 4 parallel combs and 4 series comb allpass. 
struct Diffuser{

	/// param[in] decay		Decay length of parallel combs
	Diffuser(float decay = 1.f) :
		// amounts based on Schroeder model
		comb4P(0.0297, 0.0371, 0.0411, 0.0437 , 0, 0.85),
	
		// |feedback| < 0.85 (22 dB p2p)
		// gain should be fixed at 1/sqrt(2)
		comb4S(0.0137, 0.0127, 0.0101, 0.00773, -0.707, 0.707)
	{
		this->decay(decay);
	}
	
	/// Returns next processed sample.
	float operator()(float input){
		input = comb4P.nextP(input);	// needs precision to avoid beating
		return comb4S.nextS(input);
	}
	
	/// Set 60 dB decay length of parallel combs.
	void decay(float value){ comb4P.decay(value); }

	Combs4<ipl::AllPass> comb4P;	///< Parallel feedback comb filters
	Combs4<ipl::Trunc> comb4S;		///< Series allpass comb filters	
};



template <class T=gam::real>
struct Modulet{
	Modulet(T cfreq=1000, T q=10, T mfreq=1, T mphs=0, T depth=1)
	:	fil(cfreq, q, Filter::BPC), osc(mfreq, mphs), depth(depth){}
	
	T operator()(T in){ return fil(in) * scl::mapDepth(osc.cos(), depth); }

	Biquad<> fil;
	LFO<> osc;
	T depth;
};



/// Frequency shifter
template <class T=gam::real>
struct FreqShift{

	FreqShift(float shift=1): mod(shift){}

	T operator()(T in){
		return (hil(in) * mod()).r;
	}
	
	void freq(T v){ mod.freq(v); }

	Quadra<T> mod;
	Hilbert<T> hil;
};



// Saw oscillator with sweepable filter.
struct MonoSynth{
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



/// Equal-power 2-channel panner.
template <class T=gam::real>
class Pan{
public:

	/// @param[in] pos	Signed unit position in [-1, 1]
	Pan(T pos=0){ this->pos(pos); }
	
	/// Filter sample (mono-to-stereo)
	void operator()(T in, T& out1, T& out2){
		out1 = in * w1; out2 = in * w2;
	}

	/// Filter sample (stereo-to-stereo)
	void operator()(T in1, T in2, T& out1, T& out2){
		out1 = in1 * w1 + in2 * w2;
		out2 = in1 * w2 + in2 * w1;
	}

	/// Set position
//	void pos(T v){
//		v = scl::clip(v, (T)1, (T)-1);
//		v = (v+T(1))*M_PI_4;	// put in [0, ¹/2]
//		w1 = cos(v);
//		w2 = sin(v); 
//	}

	/// Set position using a quadratic approximation.
	void pos(T v){
		// gives correct result at -1, 0, and 1
		static const T c0 = 1./sqrt(2);
		static const T c1 = 0.5 - c0;
		static const T c2 =-0.5/c1;
		v = scl::clip(v, T(1), T(-1));
		//w1 = (T)-0.25 * v * (v + (T)2) + (T)0.75;
		w1 = c1 * v * (v + c2) + c0;
		w2 = w1 + v;
	}

	/// Set position using a linear approximation.
	void posL(T v){
		v = scl::clip(v, T(1), T(-1));
		w1 = -v * T(0.5) + T(0.5);
		w2 = w1 + v; 
	}

protected:
	T w1, w2; // channel weights
};



/// Plucked string source/filter
struct Pluck{
	Pluck() : env(0.1), fil(3000, 0.2, Filter::LP), comb(1./27.5, 1./440., 1, 0.99){}
	
	float operator()(){ return comb(fil(noise() * env())); }
	float operator()(float in){ return comb(fil(in * env())); }
	void reset(){ env.reset(); }
	void freq(float v){ comb.freq(v); }
	
	NoiseWhite<> noise;
	Decay<> env;
	Biquad<> fil;
	Comb<> comb;
};



/// Quantizes sequence sampling and element magnitudes.
template <class T=gam::real>
class Quantizer : public Synced{
public:
	/// @param[in] freq		Frequency of sequence quantization
	/// @param[in] step		Step size of amplitude quantization
	Quantizer(double freq, T step=0);

	void freq(double value);	///< Set freqency of sequence quantization
	void period(double value);	///< Set period of sequence quantization
	void step(T value);			///< Set amplitude quantization amount

	T operator()(T input);		///< Return next filtered sample
	
	virtual void onSyncChange();

private:
	T mHeld;
	double mCount, mSamples, mPeriod;
	T mStep, mStepRec;
	bool mDoStep;
};

TEM Quantizer<T>::Quantizer(double freq, T step){
	this->freq(freq);
	this->step(step);
}

TEM inline void Quantizer<T>::freq(double value){ period(1./value); }

TEM inline void Quantizer<T>::period(double value){
	mPeriod = value;
	mSamples = value * spu();
	if(mSamples < 1.) mSamples = 1.;
}

TEM inline void Quantizer<T>::step(T value){
	mStep = value;
	mDoStep = mStep > 0.;
	if(mDoStep) mStepRec = 1./mStep;
}

TEM inline T Quantizer<T>::operator()(T input){
	if(++mCount >= mSamples){
		mCount -= mSamples;
		mHeld = mDoStep ? scl::round(input, mStep, mStepRec) : input;
	}
	return mHeld;
}

TEM void Quantizer<T>::onSyncChange(){
	period(mPeriod);
}



/// Compares signal magnitude to a threshold

/// This filter compares the input magnitude to a threshold and returns 1 if 
/// it's greater than the threshold and 0 otherwise.  The output is sent through
/// a one-pole low-pass filter.
template <class T=gam::real>
struct Threshold{
	/// @param[in] thresh	Comparing threshold
	/// @param[in] freq		Cutoff frequency of output smoother
	Threshold(T thresh, T freq=10):lpf(freq), thresh(thresh){}
	
	inline T operator()(T i0){ return lpf(scl::abs(i0) > thresh ? (T)1 : (T)0); }	///< Returns next sample
	inline T        inv(T i0){ return lpf(scl::abs(i0) > thresh ? (T)0 : (T)1); }	///< Returns 1 if less than threshold, 0 otherwise.
	
	OnePole<T> lpf;	///< Output smoother
	T thresh;		///< Threshold value
};



/*
// Generates 7 octaves of a unitary amplitude cosine input.
template <class T=gam::real> 
struct Fract8{
	T c[7];									///< 2nd-8th octave coefficients
	
	Fract8(T o2=0, T o3=0, T o4=0, T o5=0, T o6=0, T o7=0, T o8=0){
		c[0]=o2; c[1]=o3; c[2]=o4; c[3]=o5; c[4]=o6; c[5]=o7; c[6]=o8;
	}
	
	/// Returns filtered sample
	T operator()(T i0){ return i0 + wet(i0); }
	
	/// Returns 2nd-8th octaves of cosine input
	T wet(T i0){	
		T t = i0, o0 = (T)0;
		#define DO(i) t = t*t*(T)2 - (T)1; o0 += t * c[i];
			DO(0) DO(1) DO(2) DO(3) DO(4) DO(5) DO(6)
		#undef DO
		return o0;
	}
};
*/

} // gam::
#undef TEM
#endif
