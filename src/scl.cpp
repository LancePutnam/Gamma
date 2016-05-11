/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include <cctype> // tolower
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

double eqLoudAmp(double freq, double maxAmp){
	double ff = freq*freq;
	double c1 =    20.6; c1*=c1;
	double c2 =   107.7; c2*=c2;
	double c3 =	  737.9; c3*=c3;
	double c4 = 12200.0; c4*=c4;
	double n  = 1.258925411794167; // 10^(1/10); 2 dB offset to A-weight
	double A = ((ff + c1) * (::sqrt((ff+c2)*(ff+c3))) * (ff + c4)) / (n*ff*ff*c4);
	return A < maxAmp ? A : maxAmp;
}

double freq(const char * note){

	char c = std::tolower(*note++);
	if(within(c, 'a','g')){
		c -= 97;

		// get pitch class
		static char r[] = {9,11,0,2,4,5,7};
		char result = r[(unsigned)c];

		c = *note++;

		// apply accidental, if any
		     if(c == '+' || c == '#'){ ++result; c = *note; }
		else if(c == '-' || c == 'b'){ --result; c = *note; }
		else if(c == ' '){ c = *note; }

		// add octave
		result += (c-48)*12;

		return std::pow(2., double(result-9)/12.) * 27.5;
	}
	return 0.;
}

double nearest(double val, const char * intervals, long div){
	long vr = castIntRound(val);
	long numWraps = 0;
	long vm = wrap(vr, numWraps, div, 0L);
	long min = 0;

	struct F{
		static int base36To10(char v){
			v = std::tolower(v);
			if(v>='0' && v<='9') return v - '0';
			if(v>='a' && v<='z') return v - 'a' + 10;
			return 0;	// non-alphanumeric
		}
	};

	while(*intervals){
		long dia = F::base36To10(*intervals++);
		long max = min + dia;
		if(vm < max){	// are we within current interval?
			if(vm < (min + dia*0.5))	vm = min;
			else						vm = max;
			break;
		}
		min = max;
	}

	return double(vm + numWraps * div);
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

