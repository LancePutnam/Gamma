#ifndef GAMMA_CONVERSION_H_INC
#define GAMMA_CONVERSION_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information

	File Description:
	Functions/objects for converting amongst standard data types and strings.
*/

#include <stdio.h>
#include <iostream>
#include <sstream>		/* string conversion */
#include "Gamma/Constants.h"

namespace gam{

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



/// Convert a string of 1s and 0s to an integer.
uint32_t bits(const char * string);

/// Converts bit string to unsigned integer
uint32_t bitsToUInt(const char * bits);

/// Sets argument to zero if subnormal
void blockSubnormal(float& v);

/// Sets argument to zero if subnormal
void blockSubnormal(double& v);

/// Cast value to signed integer using rounding.

/// This is much faster than using a standard C style cast.  Floor or
/// ceiling casts can be accomplished by subtracting or adding 0.5
/// from the input, respectively.
int32_t castIntRound(double v);

template <class T> long castIntTrunc(T v);

/// Returns biased decimal value of 32-bit float exponent field.

/// The true exponent is the return value minus 127.
/// For example, values in [0.5, 1) return 126 (01111110), so the true
///	exponent is 126 - 127 = -1.
uint32_t floatExponent(float v);

/// Returns mantissa field as float between [0, 1).
float floatMantissa(float v);

/// Cast float to int

/// Reliable up to 2^24 (16777216)
///
int32_t floatToInt(float v);

/// Cast float to unsigned int

/// Reliable up to 2^24 (16777216)
///
uint32_t floatToUInt(float v);

/// Converts linear integer phase to fraction

///	2^bits is the effective size of the lookup table. \n
///	Note: the fraction only has 24-bits of precision.
float fraction(uint32_t bits, uint32_t phase);

/// Convert 16-bit signed integer to floating point in [-1, 1)
float intToUnit(int16_t v);

/// Type-pun 32-bit unsigned int to 32-bit float

/// This function uses a union to avoid problems with direct pointer casting
/// when the fstrict-aliasing compiler flag is on.
inline float punUF(uint32_t v){ Twiddle<float> u(v); return u.f; }

/// Type-pun 32-bit float to 32-bit unsigned int

/// This function uses a union to avoid problems with direct pointer casting
/// when the fstrict-aliasing compiler flag is on.
inline uint32_t punFU( float v){ Twiddle< float> u(v); return u.u; }
inline  int32_t punFI( float v){ Twiddle< float> u(v); return u.i; }

inline  int64_t punFI(  double v){ Twiddle<double> u(v); return u.i; }
inline uint64_t punFU(  double v){ Twiddle<double> u(v); return u.u; }
inline   double punUF(uint64_t v){ Twiddle<double> u(v); return u.f; }
inline   double punIF( int64_t v){ Twiddle<double> u(v); return u.f; }

/// Get fractional and integer parts from a float.

/// Works reliably up to 2^24 == 16777216
/// Useful for linearly interpolated table lookups
float split(float value, int32_t& intPart);

template<class T> T uintToUnit (uint32_t v);
template<class T> T uintToUnitS(uint32_t v);

/// Convert floating point in [0, 1) to unsigned long in [0, 2^32)

/// This conversion is most accurate on an exponential scale.
///	Input values outside [-1, 1) return 0.
///	Values in [-1, 0] behave as positive values in [0, 1).
uint32_t unitToUInt(float u);

/// Convert floating point in [0, 1) to unsigned long in [0, 2^32)

/// This conversion is most accurate on a linear scale.
/// Input values outside [0, 1) result in undefined behavior.
uint32_t unitToUInt2(float u);

/// Convert unit float in [0,1) to 8-bit unsigned int in [0, 256).
uint8_t unitToUInt8(float u);



// Implementation

inline void blockSubnormal(float& v){
	const uint32_t i = punFU(v);
	const uint32_t frac = i & MaskFrac<float>(); 
	const uint32_t expo = i & MaskExpo<float>(); 
	if(expo == 0 && frac != 0) v = 0.f;
}

inline void blockSubnormal(double& v){
	const uint64_t i = punFU(v);
	const uint64_t frac = i & MaskFrac<double>(); 
	const uint64_t expo = i & MaskExpo<double>(); 
	if(expo == 0 && frac != 0) v = 0.;
}


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
	return punFU(v) >> 23 & 0xff;
}

inline float floatMantissa(float v){
	uint32_t frac = punFU(v);
	frac = (frac & MaskFrac<float>()) | Expo1<float>();
	return punUF(frac) - 1.f;
}

inline float fraction(uint32_t bits, uint32_t phase){	
	phase = phase << bits >> 9 | Expo1<float>();
	return punUF(phase) - 1.f;
}

inline float intToUnit(int16_t v){
	uint32_t vu = (((uint32_t)v) + 0x808000) << 7; // set fraction in float [2, 4)
	return punUF(vu) - 3.f;
}

template<> inline float uintToUnit<float>(uint32_t v){
	v = v >> 9 | Expo1<float>(); 
	return punUF(v) - 1.f;
}

template<> inline float uintToUnitS<float>(uint32_t v){
	v = v >> 9 | 0x40000000;
	return punUF(v) - 3.f;
}

inline uint32_t unitToUInt(float v){
	// Convert to double in [1,2), then extract MSBs from mantissa
	// uint64_t -> uint32_t takes value modulo 2^32.
	//   (See https://en.cppreference.com/w/c/language/conversion)
	return punFU(double(v)+1.) >> 20;

	/* Old version: issues near zero
	uint32_t normalU = punFU(v);
	uint32_t rbs = 126UL - (normalU >> 23UL);
	return ((normalU | 0x800000) << 8UL) >> rbs;
	//*/
}

// Faster version, but with 23-bit precision
inline uint32_t unitToUInt2(float v){
	// Convert to float in [1,2), then extract MSBs from mantissa
	return punFU(v+1.f) << 9;
}

/*
inline uint32_t unitToUInt2(double v){
	// Convert to double in [1,2), then extract MSBs from mantissa
	return punFU(v+1.) << 12;
}//*/

inline uint8_t unitToUInt8(float u){
	++u;
	return uint8_t((punFU(u) >> 15) & MaskFrac<float>());
}

} // gam::

#endif
