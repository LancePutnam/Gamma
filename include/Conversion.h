#ifndef GAMMA_CONVERSION_H_INC
#define GAMMA_CONVERSION_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

/*	File description: 
	Functions/objects for converting between static data types.
*/

#include <stdio.h>
#include "Constants.h"
#include "Types.h"

namespace gam{


namespace{
	const double roundMagic = 6755399441055744.; // 2^52 * 1.5
	
	template<class T> const T roundEps();
	template<> inline const float  roundEps<float >(){ return 0.499999925f; }
	template<> inline const double roundEps<double>(){ return 0.499999985; }
}


/// Union for twiddling bits of floats
template<class T> struct Twiddle;

template<> struct Twiddle<float>{
	Twiddle(const float& v): f(v){}
	Twiddle(const uint32_t& v): u(v){}
	Twiddle(const int32_t& v): i(v){}
	union{ int32_t i; uint32_t u; float f; };
};

template<> struct Twiddle<double>{
	Twiddle(const double& v): f(v){}
	Twiddle(const uint64_t& v): u(v){}
	Twiddle(const int64_t& v): i(v){}
	union{ int64_t i; uint64_t u; double f; };
};


/// Convert 2-byte array to 16-bit unsigned integer.
uint16_t bytesToUInt16(const uint8_t * bytes2);

/// Convert 4-byte array to 32-bit unsigned integer.
uint32_t bytesToUInt32(const uint8_t * bytes4);

/// Cast value to signed integer using rounding.

/// This is much faster than using a standard C style cast.  Floor or
/// ceiling casts can be accomplished by subtracting or adding 0.5
/// from the input, respectively.
int32_t castIntRound(double value);

template <class T> long castIntTrunc(T value);

/// Returns biased decimal value of 32-bit float exponent field.

/// The true exponent is the return value minus 127.
/// For example, values in [0.5, 1) return 126 (01111110), so the true
///	exponent is 126 - 127 = -1.
uint32_t floatExponent(float value);

/// Returns mantissa field as float between [0, 1).
float floatMantissa(float value);

/// Cast float to long.

/// Reliable up to 2^24 (16777216)
///
long floatToInt(float value);

/// Cast float to unsigned long.

/// Reliable up to 2^24 (16777216)
///
uint32_t floatToUInt(float value);

/// Converts linear integer phase to fraction
float fraction(uint32_t bits, uint32_t phase);

/// Convert 16-bit signed integer to signed floating point normal [-1, 1).
float intToNormal(short value);

/// Convert floating point normal [0, 1) to unsigned long [0, 2^32)

/// This conversion is most accurate on an exponential scale.
///	Input values outside [-1, 1) return 0.
///	Values in [-1, 0] behave as positive values in [0, 1).
uint32_t normalToUInt(float normal);

/// Convert floating point normal [0, 1) to unsigned long [0, 2^32)

/// This conversion is most accurate on a linear scale.
/// Input values outside [0, 1) result in undefined behavior.
uint32_t normalToUInt2(float normal);

/// Maps a position in [-1, 1] to an index in [0, n). No boundary operations are performed.
inline int posToInd(float v, int n){ return n * (v*0.49999f + 0.5f); }

/// Type-pun 32-bit unsigned int to 32-bit float

/// This function uses a union to avoid problems with direct pointer casting
/// when fstrict-aliasing is on.
inline float punUF32(uint32_t v){ Twiddle<float> u(v); return u.f; }

/// Type-pun 32-bit float to 32-bit unsigned int

/// This function uses a union to avoid problems with direct pointer casting
/// when fstrict-aliasing is on.
inline uint32_t punFU32(float v){ Twiddle<float> u(v); return u.u; }

/// Get fractional and integer parts from a float.

/// Works reliably up to 2^24 == 16777216
/// Useful for linearly interpolated table lookups
float split(float value, long& intPart);

float splitInt512(uint32_t value, uint32_t& intPart);

/// Split integer accumulator into table index (size=1024) and interpolation fraction.
float splitInt1024(uint32_t value, uint32_t& intPart);

template<class T> T uintToNormal (uint32_t value);
template<class T> T uintToNormalS(uint32_t value);




// Implementation

inline int32_t castIntRound(double v){
	v += roundMagic;
	union{ double f; int32_t i[2]; } u; u.f = v;
	return u.i[endian]; // result in lsb
}

template <class T> inline long castIntTrunc(T v){
	return castIntRound( v + (v > (T)0 ? -roundEps<T>() : roundEps<T>()) );
}



/*
f32 range		u32 range		f32 exponent
[0.50,  1.00)	[2^31, 2^32)	01111110 (126)	
[0.25,  0.50)	[2^30, 2^31)	01111101 (125)	
[0.125, 0.25)	[2^29, 2^30)	01111100 (124)	

1. prepend 1 to fraction ( 'or' with 1<<24 (0x800000) )
2. shift left by 8
3. shift right according to exponent

0 01111110 00000000000000000000000

0 01111110 Fffffffffffffffffffffff	[1/2, 1/1)
1 Ffffffff fffffffffffffff00000000

0 01111101 Fffffffffffffffffffffff	[1/4, 1/2)
0 1Fffffff ffffffffffffffff0000000

0 01111100 Fffffffffffffffffffffff	[1/8, 1/4)
0 01Ffffff fffffffffffffffff000000

0 01111011 Fffffffffffffffffffffff	[1/16, 1/8)
0 001Fffff ffffffffffffffffff00000

effective  precision
9 x 2^24 + 2^23 + 2^22 + ... + 1
= 167,804,826

2^24
=  16,777,216

*/

inline uint32_t floatExponent(float v){
	return punFU32(v) << 1 >> 24;
}

inline float floatMantissa(float v){
	uint32_t frac = punFU32(v);
	frac = frac & MASK_F32_FRAC | 0x3f800000;
	return punUF32(frac) - 1.f;
}

inline float fraction(uint32_t bits, uint32_t phase){	
	phase = phase << bits >> 9 | 0x3f800000;
	return punUF32(phase) - 1.f;
}

inline float intToNormal(short v){
//	ULONG vu = ((ULONG)v) + 32768;
//	vu = vu << 7 | 0x40000000;
//	return *(float *)&vu - 3.f;

	uint32_t vu = (((uint32_t)v) + 0x808000) << 7; // set fraction in float [2, 4)
	return punUF32(vu) - 3.f;

	//return (float)v / 32768.f; // naive method
	//return (float)v * 0.000030517578125f; // less naive method
}

inline uint32_t normalToUInt(float v){
	uint32_t normalU = punFU32(v);
	uint32_t rbs = 126UL - (normalU >> 23UL);
//	printf("%x %lu\n", (normalU | 0x800000) << 8, rbs);
//	printf("%x\n", 0x80000000UL >> rbs);
	return ((normalU | 0x800000UL) << 8UL) >> rbs;
//Her00	
//float y = v + 1.f; 
//return ((unsigned long&)v) & 0x7FFFFF;      // last 23 bits 
}

inline uint32_t normalToUInt2(float v){
	v++;	// go into [1,2] range, FP fraction is now result
	return punFU32(v) << 9;
}

inline float splitInt512(uint32_t v, uint32_t& intPart){
	Twiddle<float> u(v & 0x007fffff | 0x3f800000);
	intPart = v >> 22;
	return u.f - 1.f;
}

inline float splitInt1024(uint32_t v, uint32_t& intPart){
	Twiddle<float> u((v<<1) & 0x007fffff | 0x3f800000);
	intPart = v >> 22;
	return u.f - 1.f;
}

template<> inline float uintToNormal<float>(uint32_t v){
	v = v >> 9 | 0x3f800000; 
	return punUF32(v) - 1.f;
}

template<> inline float uintToNormalS<float>(uint32_t v){
	v = v >> 9 | 0x40000000;
	return punUF32(v) - 3.f;
}

} // gam::

#endif
