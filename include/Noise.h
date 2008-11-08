#ifndef GAMMA_NOISE_H_INC
#define GAMMA_NOISE_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include "rnd.h"
#include "scl.h"
#include "MacroD.h"

namespace gam{


/// 1/f^2 noise
template <class T = RNGLinCon>
class NoiseBrown{
public:

	NoiseBrown(float val=0, float step=0.1, float min=-1, float max=1, uint32_t seed=0) 
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



/// 1/f noise

/// Summation of 12 octaves of white noise.
///
template <class T = RNGLinCon>
class NoisePink{
public:
	NoisePink();

	/// @param[in] seed Initial random number generator seed.
	NoisePink(uint32_t seed);
	
	float operator()();						///< Generate next value
	void operator()(float * arr, uint32_t len);
	void operator()(float * arr, uint32_t len, float mul);
	void seed(uint32_t v){ rng = v; }		///< Set seed value of RNG
	
	T rng;
	
private:
	float mOctave[11];
	uint32_t mPhasor;
	float mRunningSum;
	void init();
};



/// Uniform noise
template <class T = RNGLinCon>
class NoiseWhite{
public:
	NoiseWhite(): rng(){}
	NoiseWhite(uint32_t seed) : rng(seed){}
	float operator()() const { return rnd::uniS_float(rng); }	///< Generate next value
	float operator[](uint32_t i) const { return (*this)(); }
	void seed(uint32_t v){ rng = v; }					///< Set seed value of RNG
	
	mutable T rng;
};


// Qubit
template <class T = RNGLinCon>
class Qubit{
public:
	Qubit(float p): rng(), prob(p){}
	Qubit(float p, uint32_t seed) : rng(seed), prob(p){}
	bool operator()(){ return rnd::prob(rng, prob); }	///< Generate next value
	void seed(uint32_t v){ rng = v; }					///< Set seed value of RNG
	
	T rng;
	float prob;
};



//typedef uint32_t randomGen();
//typedef float randomFunc(randomGen &);

//template <class T = RNGLinCon>
//class RndUni{ public:
//	float operator()(T & rng){ return rnd::uni_float(rng); }
//	T rng;
//};
//
////
////template <randomFunc F, class T = RNGLinCon>
//template <class T = RNGLinCon>
//class Noise{
//public:
//	Noise(uint32_t seed = rnd::gen()) : rng(seed){}
//	//float operator()(){ return F(rng); }	///< Get next sample
//	
//	template <class Tf>
//	float operator()(Tf & fnc){ return fnc(rng); }
//	
//	void seed(uint32_t v){ rng = v; }			///< Set seed value of PRNG
//	
//	T rng;
//};



// Implementation_______________________________________________________________

// NoisePink

TEM NoisePink<T>::NoisePink(): rng(){ init(); }
TEM NoisePink<T>::NoisePink(uint32_t seed): rng(seed){ init(); }

TEM void NoisePink<T>::init(){
	mRunningSum = 0.f;
	LOOP(11,	/* init octaves with uniform randoms */
		float r = rnd::uniS_float(rng);
		mOctave[i] = r;
		mRunningSum += r;
	)
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

TEM inline void NoisePink<T>::operator()(float * arr, uint32_t len){
	LOOP_P(len, *arr++ = (*this)(); )
}

TEM inline void NoisePink<T>::operator()(float * arr, uint32_t len, float mul){
	LOOP_P(len,	*arr++ = (*this)() * mul; )
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

} // end namespace gam

#include "MacroU.h"

#endif

