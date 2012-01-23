/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include <ctype.h>
#include <string.h>
#include "Gamma/scl.h"

namespace gam{
namespace scl{

bool almostEqual(float a, float b, int maxUlps){
	// Make sure maxUlps is non-negative and small enough that the
	// default NAN won't compare as equal to anything.
	//assert(maxUlps > 0 && maxUlps < 4 * 1024 * 1024);

	int32_t ai = punFI(a);
	int32_t bi = punFI(b);

	// Make aInt and bInt lexicographically ordered as a twos-complement int
	if(ai < 0) ai = MaskSign<float>() - ai;
	if(bi < 0) bi = MaskSign<float>() - bi;

	return abs(ai - bi) <= maxUlps;
}


bool almostEqual(double a, double b, int maxUlps){
	// Make sure maxUlps is non-negative and small enough that the
	// default NAN won't compare as equal to anything.
	//assert(maxUlps > 0 && maxUlps < 4 * 1024 * 1024);

	int64_t ai = punFI(a);
	int64_t bi = punFI(b);

	// Make aInt and bInt lexicographically ordered as a twos-complement int
	if(ai < 0) ai = 0x8000000000000000ULL - ai;
	if(bi < 0) bi = 0x8000000000000000ULL - bi;

	return abs(ai - bi) <= maxUlps;
}


float clipMag(float value, float max, float min){
	Twiddle<float> v(value);
	uint32_t sign = v.u & MaskSign<float>();
	v.u |= MaskSign<float>();
	v.f = clip(v.f, max, min);
	v.u |= sign;
	return v.f;
}


double freq(const char * note){

	char c = *note++;
	if(within(c, 'a', 'g')){
		c -= 97;

		static char r[7] = {9,11,0,2,4,5,7};
		char result = r[(unsigned)c];
		
		c = *note++;
		     if(c == '+'){ result++; c = *note; }
		else if(c == '-'){ result--; c = *note; }
		else if(c == ' '){ c = *note; }
		
		return ::pow(2., (double)(result + (c-48)*12) / 12.) * 8.1757989157741;		
	}
	return 0.;
}

/*

exponent will not be modified
must quantize fraction

1. if	exponent is greater than 01111110, continue
   else return 0
2. zero last n bits in fraction

range		exp				zero lsb
[1/2, 1)	01111110 (126)	>= 32
[1, 2)		01111111 (127)	23
[2, 4)		10000000 (128)	22
[4, 8)		10000001 (129)	21

0 01111110 11111111111111111111111	just under 1

0 01111110 Fffffffffffffffffffffff	[1/2, 1)

0 01111111 Fffffffffffffffffffffff	[1, 2)

0 10000000 Fffffffffffffffffffffff	[2, 4)

0 01111111 10000000000000000000000	1.5
0 01111111 00000000000000000000000	1

0 10000000 01000000000000000000000	2.5
0 10000000 00000000000000000000000	2

0 10000000 11000000000000000000000	3.5
0 10000000 10000000000000000000000	3

0 00000000 00000000000000000000000	[0, 1)	huge
1 11111111 00000000000000000000000	[1, 2)	0
1 11111111 10000000000000000000000	[2, 4)	1
1 11111111 11000000000000000000000	[4, 8)	2

0 00000000 11111111111111111111111
1 11111111 00000000000000000000000

*/


// branchless
//float SclOp::fastFloor(float value){
//	ULONG valueU = *(ULONG *)&value;
////	ULONG shift = 23U - ((valueU >> 23U & 0xff) - 127U);
////	valueU = valueU >> shift << shift;
//	ULONG shift = (valueU >> 23U & 0xff) - 127U;
//	valueU &= ~(~(0xff800000 >> shift));
//	//printBinary((ULONG)~(0x7FFFFF >> shift)); printf("\n");
//	return *(float *)&valueU;
//}

//float SclOp::fastFloor(float value){
//	union { float f; ULONG u; } word;
//	word.f = value;
//
//	ULONG bits = word.u;
//	bits += 0x800000;
//
//	if(bits & 0x40000000){	// mag outside [0, 1)
//		int shift = 23U - (bits >> 23U & 0x7F);
//		bits = word.u;
//		bits = bits >> shift << shift;
//		word.u = bits;
//		return word.f;
//	}
//	else{
//		return 0.f;
//	}
//}

//inline float SclOp::truncFast(float value){
//
//	ULONG valueU = *(ULONG *)&value;
//	
//	// Is value outside [0, 1) ?
//	if((valueU & 0x7fffffff) > 0x3f7fffff){
//		ULONG shift = 150UL - (valueU >> 23UL & 0x7f);
//		valueU = valueU >> shift << shift;
//		return *(float *)&valueU;
//	}
//	
//	return 0.f;
//}

} // scl::
} // gam::

