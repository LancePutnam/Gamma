/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include <stdlib.h>
#include "Gamma/Print.h"

namespace gam{

void colorHSV(float r, float g, float b, float &h, float &s, float &v){

	float min = r < g ? (r < b ? r : b) : (g < b ? g : b);
	float max = r > g ? (r > b ? r : b) : (g > b ? g : b);

	v = max;					// v
	float delta = max - min;	// delta RGB value

	if ( delta != 0.f && max != 0.f ){		// chromatic data...
		s = delta / max;		// s
		
		float hl;
		if     ( r == max )	hl =       ( g - b ) / delta;	// between yellow & magenta
		else if( g == max )	hl = 2.f + ( b - r ) / delta;	// between cyan & yellow
		else				hl = 4.f + ( r - g ) / delta;	// between magenta & cyan

		if( hl < 0.f ) hl += 6.f;

		h = hl * 0.166666667f;
	}
	else{				// this is a gray, no chroma...
	   h = 0.f;
	   s = 0.f;
	}
}


void colorRGB(float h, float s, float v, float &r, float &g, float &b){
	
	if( s == 0.f ) {
		r = g = b = v;	// achromatic (grey)
		return;
	}
	
	h *= 6.f;
										// sector 0 to 5
	unsigned int i = (unsigned int)(h);	// integer part of h
	float f = h - (float)i;				// fractional part of h
	float p = v * ( 1.f - s );
		
	float q;	// depends on hue section
	if(i & 1U)	q = v * ( 1.f - s * f );			// odd
	else		q = v * ( 1.f - s * ( 1.f - f ) );	// even

	switch( i ) {
		case 0: r = v; g = q; b = p; break;
		case 1:	r = q; g = v; b = p; break;
		case 2:	r = p; g = v; b = q; break;
		case 3:	r = p; g = q; b = v; break;
		case 4: r = q; g = p; b = v; break;
		default:r = v; g = p; b = q; break;
	}
}

#define LOOP_BITS(exp) for(int i=(msb-1); i>=0; --i){ exp }
void printBinary(uint32_t v, const char * zero, const char * one, int msb){
	LOOP_BITS(
		0 == ((v>>i) & 1) ? printf("%s", zero) : printf("%s", one);
	)
}

void printBinary(uint64_t v, const char * zero, const char * one, int msb){
	LOOP_BITS(
		0 == ((v>>i) & 1) ? printf("%s", zero) : printf("%s", one);
	)
}

void printBinary(float value, const char * zero, const char * one, int msb){
	Twiddle<float> v(value);
	LOOP_BITS(
		0 == ((v.u>>i) & 1) ? printf("%s", zero) : printf("%s", one);
		if((i==31) || (i==23)) printf(" ");
	)
}

void printBinary(void * value32, const char * zero, const char * one, int msb){
	printBinary(*(uint32_t *)value32, zero, one, msb);
}
#undef LOOP_BITS

void printHexArray(float * a, uint32_t len, uint32_t valuesPerLine){
	printf("{");
	for(uint32_t i=0; i<len; ++i){
		if((i % valuesPerLine) == 0) printf("\n\t");
		Twiddle<float> t(a[i]);
		printf("0x%08x%s", t.u, i == len-1 ? "\n};" : ",");
	}
}

void printPlot(float value, uint32_t width, bool spaces, const char * point){
	int clipFlag;
	value = scl::clip(value, clipFlag, 1.f, -1.f);
	
	const char * pt = clipFlag != 0 ? "+" : point;
	
	uint32_t pos = castIntRound((value + 1.f) * 0.5f * (float)(width));
	uint32_t mid = width >> 1;
	uint32_t i=0;

	if(pos < mid){	// [-1, 0)
		for(; i<pos; ++i) printf(" ");
		printf("%s", pt); ++i;
		for(; i<mid; ++i) printf("-");
		printf("|");
	}
	else{			// (0, 1]
		for(; i<mid; ++i) printf(" ");
		if(pos == mid){ printf("%s", pt); goto end; }
		printf("|"); ++i;
		for(; i<pos; ++i) printf("-");
		printf("%s", pt);
	}
	
	end: 
	if(spaces) for(; i<width; ++i) printf(" ");
}


void err(const char * msg, const char * src, bool exits){
	fprintf(stderr, "%s%serror: %s\n", src, src[0]?" ":"", msg);
	if(exits) exit(EXIT_FAILURE);
}

void warn(const char * msg, const char * src){
	fprintf(stderr, "%s%swarning: %s\n", src, src[0]?" ":"", msg);
}


} // gam::
