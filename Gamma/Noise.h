#ifndef GAMMA_NOISE_H_INC
#define GAMMA_NOISE_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include "Gamma/rnd.h"
#include "Gamma/scl.h"

namespace gam{


/// Brownian noise

/// Brownian noise has a power spectrum of 1/f^2.
/// It is produced by integrating white (uniform) noise.
template <class T=RNGLinCon>
class NoiseBrown{
public:

	NoiseBrown(float val=0, float step=0.04, float min=-1, float max=1, uint32_t seed=0) 
	:	val(val), step(step), min(min), max(max)
	{	if(seed) rng = seed; }
	
	/// Generate next value
	float operator()(){
		val = scl::clip(val + rnd::uniS_float(rng) * step, max, min);
		return val;
	}
	
	void seed(uint32_t v){ rng = v; }	///< Set seed value of RNG
	
	T rng;
	float val, step, min, max;
};


/// Pink Noise

/// Pink noise has a power spectrum of 1/f.
/// It is produced by summing together 12 octaves of white noise.
template <class T=RNGLinCon>
class NoisePink{
public:
	NoisePink();

	/// @param[in] seed Initial random number generator seed.
	NoisePink(uint32_t seed);
	
	float operator()();						///< Generate next value
	void operator()(float * arr, uint32_t len, float mul);
	void seed(uint32_t v){ rng = v; }		///< Set seed value of RNG
	
	T rng;
	
private:
	float mOctave[11];
	uint32_t mPhasor;
	float mRunningSum;
	void init();
};


/// White noise

/// White noise has a uniform power spectrum.
///
template <class T=RNGLinCon>
class NoiseWhite{
public:
	NoiseWhite(): rng(){}
	NoiseWhite(uint32_t seed) : rng(seed){}
	float operator()() const { return rnd::uniS_float(rng); }	///< Generate next value
	float operator[](uint32_t i) const { return (*this)(); }
	void seed(uint32_t v){ rng = v; }					///< Set seed value of RNG
	
	mutable T rng;
};




// Implementation_______________________________________________________________

#define TEM template<class T>

// NoisePink

TEM NoisePink<T>::NoisePink(): rng(){ init(); }
TEM NoisePink<T>::NoisePink(uint32_t seed): rng(seed){ init(); }

TEM void NoisePink<T>::init(){
	mRunningSum = 0.f;
	for(uint32_t i=0; i<11; ++i){	/* init octaves with uniform randoms */
		float r = rnd::uniS_float(rng);
		mOctave[i] = r;
		mRunningSum += r;
	}
	mPhasor = 0;
}

TEM inline float NoisePink<T>::operator()(){
	// phasor range:	[1, 2048]
	//					[0, 10]		trailing zeroes
	mPhasor++;
	if(mPhasor != 2048){
		uint32_t i = scl::trailingZeroes(mPhasor);	// # trailing zeroes is octave
		float r = rnd::uniS_float(rng);			// uniform random
		mRunningSum += r - mOctave[i];			// add new & subtract stored
		mOctave[i] = r;							// store new
	}
	else{
		mPhasor = 0;
	}
	
	// add white noise every sample
	return (mRunningSum + rnd::uniS_float(rng)) * 0.083333333f;
}

TEM inline void NoisePink<T>::operator()(float * arr, uint32_t len, float mul){	
	for(uint32_t i=0; i<len; ++i) arr[i] = (*this)() * mul;
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

#undef TEM

#endif
