#ifndef GAMMA_RND_H_INC
#define GAMMA_RND_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information

	File Description:
	Random variable classes
*/

#include <math.h>
#include <stdio.h>
#include <time.h>	/* for time() */

#include "Gamma/gen.h"
#include "Gamma/mem.h"
#include "Gamma/scl.h"
#include "Gamma/Conversion.h"
#include "Gamma/Types.h"

#define TEM template <class T>

namespace gam{


namespace rnd{
	namespace{
		static bool initSeed = false;
		static uint32_t mSeedPush[4];
		static gen::RMulAdd<uint32_t> seedGen(1664525, 1013904223);
	}
	
	/// Get a random seed
	static uint32_t getSeed(){
		if(!initSeed){ seedGen.val = time(NULL); initSeed = true; } 
		return seedGen();
	}

} // rnd::


/// Linear congruential uniform pseudo-random number generator.

///	This generator is very fast requiring only a single integer multiply and add per
/// iteration.  However, the least significant bits of the numbers are less 
/// random; the most extreme case being the LSB which at best flips between 0 and 1.
/// This generator also exhibits poor dimensional distribution, therefore it is
/// best to have a different generator for each dimension, rather than sharing one.
struct RNGLinCon : public gen::RMulAdd<uint32_t>{

	RNGLinCon(){ val=rnd::getSeed(); type(0); }

	/// @param[in] seed	Initial seed value
	RNGLinCon(uint32_t seed): gen::RMulAdd<uint32_t>(1,0,seed){ type(0); }

	/// Change the type of equation used.
	
	/// 0 - Knuth, Numerical Recipes in C\n
	/// 1 - BCPL
	void type(int v){
		switch(v){
		case 1:	mul = 2147001325; add =  715136305; break; // BCPL
		default:mul =    1664525; add = 1013904223;        // Knuth, Numerical Recipes in C
		}
	}
};



/// Multiplicative linear congruential uniform pseudo-random number generator.

///	This generator is a faster LCG requiring only a single integer multiply.
///
struct RNGMulLinCon : public gen::RMul<uint32_t>{

	RNGMulLinCon(){ val=rnd::getSeed(); type(0); }
	
	/// @param[in] seed	Initial seed value
	RNGMulLinCon(uint32_t seed): gen::RMul<uint32_t>(1,seed){ type(0); }
	
	/// Change the type of equation used.

	/// 0 - Marsaglia, Super-Duper\n
	///
	void type(int v){
		switch(v){
		default: mul = 69069;	// Super-duper
		}
	}
};



/// Combined Tausworthe uniform pseudo-random number generator.

/// This generator produces highly random numbers, but is more expensive than
/// than a linear congruential RNG.
/// It is based on the paper 
/// P. L'Ecuyer, "Maximally Equidistributed Combined Tausworthe Generators", 
/// Mathematics of Computation, 65, 213 (1996), 203--213.
/// http://www.iro.umontreal.ca/~lecuyer/papers.html
class RNGTaus{
public:
	RNGTaus(){ (*this) = rnd::getSeed(); }

	/// @param[in] seed		Initial seed value
	RNGTaus(uint32_t seed);
	
	uint32_t s1, s2, s3, s4;
	
	uint32_t operator()();				///< Generates uniform random unsigned integer in [0, 2^32).
	void operator = (uint32_t seed);	///< Set seed
	void seed(uint32_t s1, uint32_t s2, uint32_t s3, uint32_t s4); ///< Set seed

private:
	void iterate();
};



/// Random number functions

/// Unless specified, these operations use an internal Tausworthe RNG, favoring
/// quality over a slight deficiency in speed compared to a linear congruential
/// RNG.\n
///\n
/// The following template functions are available:\n
///\code
/// // Returns random number linearly mapped to [bound1, bound2)
/// T distr(T bound2, T bound1 = 0);
///
/// // Fills array with random numbers in [0, 1)
/// void distr(T * dst, uint32_t len);
///
/// // Fills array with random numbers linearly mapped to [bound1, bound2)
/// void distr(T * dst, uint32_t len, T bound2, T bound1 = 0);
///
/// // Copies random elements from 'src' to 'dst'.
/// void distr(T * dst, uint32_t len, T * src, uint32_t srcLen);
///\endcode
///
/// where "distr" is one of the following random number distributions:\n\n
/// add2  - sum of 2 uniform values (triangle) \n
/// add2I - inverse pdf of add2 \n
/// add3  - sum of 3 uniform values (2nd order normal curve) \n
/// binS  - uniform signed binary \n
///	lin   - linear decrease \n
/// mul2  - product of 2 uniform values \n
/// pow2  - uniform value squared \n
/// pow3  - uniform value cubed \n
/// uni   - uniform value \n
///\n
/// Since casts between floats and generic types are used in 
/// the ranged varieties, signed integer bounds will always be exclusive.\n\n

namespace rnd{

	/// Conditional probability with two values.
	
	/// 'pab' and 'pba' are the probabilities of a given b and b given a.
	///
	TEM T& cond(T& v, const T& va, const T& vb, float pab, float pba);

	/// Returns ith element in geometric series
	TEM T geom(int n, T mul, T start=1, float p=0.5);

	/// Returns value negated with probability, otherwise the value.
	TEM T neg(T val, float prob=0.5f);

	/// Returns val1 with probability, otherwise val2.
	TEM const T& pick(const T& val1, const T& val2, float prob=0.5f);

	/// Returns true with a probability of p.
	TEM bool prob(T& rng, float p=0.5f);
	
	/// Returns true with a probability of p.
	bool prob(float p=0.5f);
	bool prob(double p);
	
	/// Characters 0-8 return true with probability c/8.
	
	/// In addition to numeric characters, '/' returns true and '.' returns false.
	///
	bool prob(char c);

	#ifndef DOXYGEN_SHOULD_SKIP_THIS
	#define FUNCS(name)\
	TEM T name(T bound2, T bound1 = 0);\
	TEM void name(T * dst, uint32_t len);\
	TEM void name(T * dst, uint32_t len, T bound2, T bound1 = 0);\
	TEM void name(T * dst, uint32_t len, T * src, uint32_t srcLen);\
	TEM float name##_float(T & rng);
	//template <class Tv, class Tr> Tv name(Tr & rng);

	FUNCS(add2) FUNCS(add2I) FUNCS(add3) FUNCS(binS) FUNCS(lin) FUNCS(mul2)
	FUNCS(pow2) FUNCS(pow3) FUNCS(uni) FUNCS(uniS)

	#undef FUNCS
	#endif

	/// Push current RNG state onto stack (stack size = 1).
	
	/// After pushing, the current RNG is seeded with 'seed' unless 'seed' = 0.
	///
	void push(uint32_t seed=0);
	void pop();					///< Pop RNG state from stack.

	/// Randomly permutes (shuffles) elements in array.
	TEM void permute(T * arr, uint32_t len);

	/// Returns value in [0, 1) quantized by q divisions.
	float quan(uint32_t q);
	
	/// Returns value in [o, 2*o) quantized by q divisions.
	TEM T quanOct(uint32_t q, T o);

	/// Seed global PRNG. If seed is 0, then the system time (in seconds) is used.
	void seed(uint32_t value=0);

	/// Randomly set a certain amount of elements to a value.
	TEM void set(T * arr, uint32_t len, uint32_t num, T val=1);

	/// Zeroes elements according to a probability.
	TEM void thin(T * arr, uint32_t len, float prob=0.5f);
	
	/// Returns uniform random within interval [min, max) excluding 'exc' argument.
	TEM T uniExc(const T& exc, const T& max, const T& min=T(0));
	
	/// Returns random integer in [0, num) according to weights (a PDF).
	
	/// If the weights are not normalized, then the proper weightsSum must
	/// be passed in.
	TEM uint32_t weighted(T * weights, uint32_t num, T weightsSum=(T)1);
	
	static RNGTaus gen(rnd::getSeed());	///< Shared RNG
}


// Implementation_______________________________________________________________


//---- RNGTaus

inline RNGTaus::RNGTaus(uint32_t sd){ (*this) = sd; }

inline uint32_t RNGTaus::operator()(){
	iterate();
	return s1 ^ s2 ^ s3 ^ s4;
}

inline void RNGTaus::operator=(uint32_t s){
	gen::RMulAdd<uint32_t> g(1664525, 1013904223, s);
	g.val = s; g();
	seed(g(), g(), g(), g());
}

inline void RNGTaus::seed(uint32_t v1, uint32_t v2, uint32_t v3, uint32_t v4){
	//printf("%d %d %d %d\n", v1, v2, v3, v4);
	v1 & 0xffffffe ? s1 = v1 : s1 = ~v1;
	v2 & 0xffffff8 ? s2 = v2 : s2 = ~v2;
	v3 & 0xffffff0 ? s3 = v3 : s3 = ~v3;
	v4 & 0xfffff80 ? s4 = v4 : s4 = ~v4;
}
	
inline void RNGTaus::iterate(){
	s1 = ((s1 & 0xfffffffe) << 18) ^ (((s1 <<  6) ^ s1) >> 13);
	s2 = ((s2 & 0xfffffff8) <<  2) ^ (((s2 <<  2) ^ s2) >> 27);
	s3 = ((s3 & 0xfffffff0) <<  7) ^ (((s3 << 13) ^ s3) >> 21);
	s4 = ((s4 & 0xffffff80) << 13) ^ (((s4 <<  3) ^ s4) >> 12);
}


//---- rnd
namespace rnd{

#define R uni_float(rng)
TEM inline float add2_float(T& rng){ return (R + R) * 0.5f; }
TEM inline float add3_float(T& rng){ return (R + R + R) * 0.33333333333f; }
TEM inline float  lin_float(T& rng){ float r = R; float s = R; return r < s ? r : s; }
TEM inline float mul2_float(T& rng){ return R * R; }
TEM inline float pow2_float(T& rng){ float r = R; return r * r; }
TEM inline float pow3_float(T& rng){ float r = R; return r * r * r;	}

TEM inline float add2I_float(T & rng){
	float r = add2_float(rng) + 0.5f; // [0.5, 1.5)
	return (r < 1.f) ? r : r - 1.f;
}
#undef R

TEM inline float uni_float(T& rng){ return uintToUnit<float>(rng()); }
TEM inline float uniS_float(T& rng){ return uintToUnitS<float>(rng()); }

TEM inline float binS_float(T & rng){
	uint32_t r = rng() & MaskSign<float>() | Expo1<float>();
	return punUF(r);
}

TEM inline T & cond(T& v, const T& va, const T& vb, float pab, float pba){
	     if(v == va) v = pick(vb, va, pba);
	else if(v == vb) v = pick(va, vb, pab);
	return v;
}

template <class RNG>
float gaussian(RNG& rng = rnd::gen){
	float x1, x2, w, y1, y2;

	do{
		x1 = uniS_float(rng);
		x2 = uniS_float(rng);
		w = x1 * x1 + x2 * x2;
	} while( w >= 1.f );

	w = sqrt((-2.f * ::log(w)) / w);
	y1 = x1 * w;
	y2 = x2 * w;
	return y1;
}

TEM inline T geom(int n, T mul, T start, float p){
	for(; n>0; --n){
		if(prob(p)) break;
		start *= mul;
	}
	return start;
}

TEM inline T neg(T v, float p){ return prob(p) ? -v : v; }
TEM inline const T & pick(const T & v1, const T & v2, float p){ return prob(p) ? v1 : v2; }
TEM inline bool prob(T& rng, float p){ return uni_float(rng) < p; }
inline bool prob(float p){ return prob(gen, p); }

inline bool prob(double p){ return prob((float)p); }

inline bool prob(char c){
	if(c>'7'||c=='/') return true; 
	if(c<'1')         return false;
	return prob((c-48) * 0.125f);
}

#define LOOP(n) for(uint32_t i=0;i<n;++i)

TEM inline void set(T * arr, uint32_t len, uint32_t num, T val){
	LOOP(num){ arr[rnd::uni(len)] = val; }
}

TEM inline void permute(T * arr, uint32_t len){
	LOOP(len-1){ mem::swap(arr[i], arr[rnd::uni(len, i)]); }
}

inline float quan(uint32_t q){ return uni(q) / (float)q; }

TEM inline T quanOct(uint32_t q, T o){ return quan(q) * o + o; }

TEM inline void thin(T * arr, uint32_t len, float p){
	LOOP(len){ if(prob(p)) arr[i]=T(0); }
}

TEM inline T uniExc(const T& exc, const T& max, const T& min){
	T r=exc; while(exc==r){ r=uni(max,min); } return r;
}

#define DEF(rnd_t, fnc)\
	TEM inline T fnc(T b2, T b1){\
		return (T)( (rnd_t)b1 + fnc##_##rnd_t(gen) * (rnd_t)(b2-b1) );\
	}\
	TEM inline void fnc(T * dst, uint32_t len){\
		LOOP(len){ dst[i] = (T) fnc##_##rnd_t (gen); }\
	}\
	TEM inline void fnc(T * dst, uint32_t len, T bound2, T bound1){\
		rnd_t b1 = (rnd_t)bound1;\
		rnd_t df = (rnd_t)bound2 - b1;\
		LOOP(len){ dst[i] = (T)(b1 + fnc##_##rnd_t(gen) * df); }\
	}\
	TEM inline void fnc(T * dst, uint32_t len, T * src, uint32_t srcLen){\
		LOOP(len){ dst[i] = src[rnd::fnc(srcLen)]; }\
	}
	DEF(float, add2) DEF(float, add2I) DEF(float, add3) DEF(float, lin)
	DEF(float, mul2) DEF(float, pow2 ) DEF(float, pow3) DEF(float, uni)
	DEF(float, uniS)
#undef DEF

#undef LOOP

inline void push(uint32_t seedA){
	mSeedPush[0] = gen.s1;
	mSeedPush[1] = gen.s2;
	mSeedPush[2] = gen.s3;
	mSeedPush[3] = gen.s4;
	if(seedA) seed(seedA);
}
inline void pop(){
	gen.s1 = mSeedPush[0];
	gen.s2 = mSeedPush[1];
	gen.s3 = mSeedPush[2];
	gen.s4 = mSeedPush[3];
}

inline void seed(uint32_t value){
	gen = value ? value : getSeed();
}

TEM uint32_t weighted(T * weights, uint32_t num, T weightsSum){
	if(0 == num--) return 0;
	T probSum = weights[0];
	T rnd = (T)(uni_float(gen) * (float)weightsSum);
	
	if(rnd < probSum) return 0;
	
	for(uint32_t i=1; i<num; ++i){
		probSum += weights[i];
		if(rnd < probSum) return i;
	}
	return num;
}


/*

// some distributions from old code

inline float RandomLC::exp(){
	uint32_t rand = nextU() & 0x3fffffff;  // 00111111111111111111111111111111
	return AS_FLOAT(rand) * 0.5f;
}

inline float RandomLC::expS(){
	uint32_t rand = nextU() & 0xbfffffff;	// 10111111111111111111111111111111
	return AS_FLOAT(rand) * 0.5f;
}

inline float RandomLC::gaussian(){
	return sqrt(-2.f * log(next1())) * sin(next1() * 6.28318530718f);
}

*/

} // rnd::
} // gam::

#undef TEM

#endif
