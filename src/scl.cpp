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

} // scl::
} // gam::
