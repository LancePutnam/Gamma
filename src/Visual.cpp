/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include "Visual.h"

namespace gam{


#define LOOP_BITS(exp) for(int i=(msb-1); i>=0; --i){ exp }

void printBinary(uint32_t v, const char * zero, const char * one, int msb){
	LOOP_BITS(
		0 == ((v>>i) & 1) ? printf(zero) : printf(one);
	)
}

void printBinary(uint64_t v, const char * zero, const char * one, int msb){
	LOOP_BITS(
		0 == ((v>>i) & 1) ? printf(zero) : printf(one);
	)
}

void printBinary(float value, const char * zero, const char * one, int msb){
	FloatUInt<float> v(value);
	LOOP_BITS(
		0 == ((v.i>>i) & 1) ? printf(zero) : printf(one);
		if((i==31) || (i==23)) printf(" ");
	)
}

void printBinary(void * value32, const char * zero, const char * one, int msb){
	printBinary(*(uint32_t *)value32, zero, one, msb);
}

#undef LOOP_BITS

void printPlot(float value, uint32_t width, bool spaces, const char * point){
	int clipFlag;
	value = scl::clip(value, clipFlag, 1.f, -1.f);
	
	const char * pt = clipFlag != 0 ? "+" : point;
	
	uint32_t pos = castIntRound((value + 1.f) * 0.5f * (float)(width));
	uint32_t mid = width >> 1;
	uint32_t i=0;

	if(pos < mid){	// [-1, 0)
		for(; i<pos; ++i) printf(" ");
		printf(pt); ++i;
		for(; i<mid; ++i) printf("-");
		printf("|");
	}
	else{			// (0, 1]
		for(; i<mid; ++i) printf(" ");
		if(pos == mid){ printf(pt); goto end; }
		printf("|"); ++i;
		for(; i<pos; ++i) printf("-");
		printf(pt);
	}
	
	end: 
	if(spaces) for(; i<width; ++i) printf(" ");
}


} // gam::
