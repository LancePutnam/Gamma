#ifndef GAMMA_NOISE_H_INC
#define GAMMA_NOISE_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include "Gamma/rnd.h"
#include "Gamma/scl.h"

namespace gam{

/// Random signals

///\defgroup Noise


template <class RNG>
class NoiseBase{
public:
	typedef float value_type;

	NoiseBase(){}

	/// \param[in] seed		random number generator seed
	NoiseBase(uint32_t seed): rng(seed){}

	/// Set seed of RNG
	void seed(uint32_t v){ rng = v; }
	
	RNG rng;

protected:
	float U(){ return rnd::uni_float(rng); }
	float Us(){ return rnd::uniS_float(rng); }
	bool P(float p=0.5){ return U() < p; }
};


/// Brown noise

/// Brown noise has a power spectrum of 1/f^2. This corresponds to an amplitude
/// spectrum with a -6 dB/octave slope. Here it is produced by integrating
/// white noise. The output value is clipped within a specified interval.
/// \ingroup Noise
template <class RNG = RNGLinCon>
class NoiseBrown : public NoiseBase<RNG> {
public:

	/// \param[in] val		start value
	/// \param[in] step		accumulation step factor
	/// \param[in] min		minimum of clipping interval
	/// \param[in] max		maximum of clipping interval
	/// \param[in] seedi	random number generator seed; 0 generates a random seed
	NoiseBrown(float val=0, float step=0.04, float min=-1, float max=1, uint32_t seedi=0) 
	:	val(val), step(step), min(min), max(max)
	{	if(seedi) this->seed(seedi); }
	
	/// Generate next value
	float operator()(){
		val = scl::clip(val + this->Us() * step, max, min);
		return val;
	}
	
	float val, step, min, max;
};


/// Pink Noise

/// Pink noise has a power spectrum of 1/f. This corresponds to an amplitude
/// spectrum with a -3 dB/octave slope.
/// In this implementation, it is produced by summing together 12 octaves of
/// downsampled white noise.
/// \ingroup Noise
template <class RNG = RNGLinCon>
class NoisePink : public NoiseBase<RNG> {
public:

	NoisePink();

	/// \param[in] seed		random number generator seed
	NoisePink(uint32_t seed);

	/// Generate next value
	float operator()();
	
private:
	float mOctave[11];
	uint32_t mPhase;
	float mRunningSum;
	void init();
};


/// White noise

/// White noise has a uniform power spectrum.
/// \ingroup Noise
template <class RNG = RNGLinCon>
class NoiseWhite : public NoiseBase<RNG> {
public:

	/// Generate next value
	float operator()(){ return this->Us(); }
};


/// Violet noise

/// Violet noise has a power spectrum of f^2. This corresponds to an amplitude
/// spectrum with a +6 dB/octave slope. Here it is produced by differentiating
/// white noise.
/// \ingroup Noise
template <class RNG = RNGLinCon>
class NoiseViolet : public NoiseBase<RNG> {
public:
	
	/// Generate next value
	float operator()(){
		auto curr = punUF(Expo1<float>() | (this->rng()>>9)); // in [1,2)
		auto diff = curr - mPrev;
		mPrev = curr;
		return diff;
	}

private:
	float mPrev = 1.5f;
};


/// Binary noise

/// This noise flips randomly between -A and A. When A is very small, say 1e-20,
/// then this can be added to the input of filters to prevent denormals.
/// \ingroup Noise
template <class RNG = RNGLinCon>
class NoiseBinary : public NoiseBase<RNG> {
public:

	/// \param[in] ampi		peak-to-peak amplitude of noise
	/// \param[in] seedi	random number generator seedi; 0 generates a random seed
	NoiseBinary(float ampi = 1.f, uint32_t seedi = 0)
	:	amp(ampi)
	{	if(seedi) this->seed(seedi); }

	/// Generate next value
	float operator()(){
		return punUF((this->rng()&0x80000000) ^ punFU(amp));
	}

	float amp;
};


// Implementation_______________________________________________________________

template<class T>
NoisePink<T>::NoisePink(){
	init();
}
    
template<class T>
NoisePink<T>::NoisePink(uint32_t seed)
:	NoiseBase<T>(seed)
{
	init();
}

template<class T>
void NoisePink<T>::init(){
	mRunningSum = 0.f;
	for(uint32_t i=0; i<11; ++i){	/* init octaves with uniform randoms */
		auto r = this->Us();
		mOctave[i] = r;
		mRunningSum += r;
	}
	mPhase = 0;
}

template<class T>
inline float NoisePink<T>::operator()(){
	// phasor range:	[1, 2048]
	//					[0, 10]		trailing zeroes
	++mPhase;
	if(mPhase != 2048){
		uint32_t i = scl::trailingZeroes(mPhase);	// # trailing zeroes is octave
		auto r = this->Us();					// uniform random
		mRunningSum += r - mOctave[i];			// add new & subtract stored
		mOctave[i] = r;							// store new
	}
	else{
		mPhase = 0;
	}
	
	// add white noise every sample
	return (mRunningSum + this->Us()) * 0.083333333f;
}


/*
////////////////////////////////////////////////////////////////////////////////
Pink Noise
	Description:
	Infinite sum of octaves of white noise.  Since we cannot compute an infinite
	sum, we'll use 12 octaves. This gives a range of 10.77 - 44,100 Hz.

	12345678901234567890123456789012	octave  length
	********************************	0		0
	* * * * * * * * * * * * * * * *		1		2
	 *   *   *   *   *   *   *   *		2		4
	   *       *       *       *		3		8
	       *               *			4		16
		           *					5		32
	12131214121312151213121412131210

   0				0		00000	
	1		000		1		00001
	 2		001		2		00010
	1		000		3		00011
	  3		010		4		00100
	1		000		5		00101
	 2		001		6		00110
	1		000		7		00111
	    4   011		8		01000
	1		000		9		01001
	 2		001		10		01010
	1		000		11		01011
	   3	010		12		01100
	1		000		13		01101
	 2		001		14		01110
	1		000		15		01111
	     5  100		16		10000
	1		000		17		10001
	 2		001		18		10010
	1		000		19		10011
	  3		010		20		10100
	1		000		21		10101
	 2		001		22		10110
	1		000		23		10111
	    4   011		24		11000
	1		000		25		11001
	 2		001		26		11010
	1		000		27		11011
	   3	010		28		11100
	1		000		29		11101
	 2		001		30		11110
	1		000		31		11111
	
*/

} // gam::

#endif
