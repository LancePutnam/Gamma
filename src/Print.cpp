/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include <cstdlib> // exit
#include "Gamma/Print.h"

namespace gam{

std::string plotString(
	float value, unsigned width, bool spaces, bool sign, const char * point
){
	std::string res;

	float min, max;
	if(sign)	{ min=-1; max=1; }
	else		{ min= 0; max=1; }
	float dia = max-min;

	int clipFlag;
	value = scl::clip(value, clipFlag, max, min);
	
	const char * pt = clipFlag != 0 ? "+" : point;

	int imin = min/dia * width;	// normalize by diameter, then multiply by width
	int imax = imin + width;
	int pos = castIntRound(value / dia * width);

	if(!spaces) imax = pos<=0 ? 1 : pos;

	if(pos >= imax) pos = imax-1;
	
	for(int i=imin; i<imax; ++i){
		if(i == pos)	res += pt;
		else if(i==0)	res += "|";
		else if((i>pos && i<0) || (i<pos && i>0))	res += "-";
		else			res += " ";
	}
	return res;
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

void printBinary(const void * value32, const char * zero, const char * one, int msb){
	printBinary(*(uint32_t *)value32, zero, one, msb);
}
#undef LOOP_BITS

void printHexArray(const float * a, unsigned len,unsigned valuesPerLine){
	printf("{");
	for(unsigned i=0; i<len; ++i){
		if((i % valuesPerLine) == 0) printf("\n\t");
		Twiddle<float> t(a[i]);
		printf("0x%08x%s", t.u, i == len-1 ? "\n};" : ",");
	}
}

void printPlot(float value, unsigned width, bool spaces, bool sign, const char * point){
	fprintf(stdout, "%s", plotString(value,width,spaces,sign,point).c_str());
}


void err(const char * msg, const char * src, bool exits){
	fprintf(stderr, "%s%serror: %s\n", src, src[0]?" ":"", msg);
	if(exits) std::exit(EXIT_FAILURE);
}

void warn(const char * msg, const char * src){
	fprintf(stderr, "%s%swarning: %s\n", src, src[0]?" ":"", msg);
}


} // gam::
