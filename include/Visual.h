#ifndef GAMMA_VISUAL_H_INC
#define GAMMA_VISUAL_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

/*	File description: 
	Visual functions for printing and colors.
*/

#include <stdio.h>
#include "Constants.h"
#include "Types.h"
#include "scl.h"

namespace gam{


/// Returns an ASCII character most closely matching an intensity value in [0,1].
char intensityToASCII(float v);

template<class T> void print(const T& v, const char * post="", const char * pre="", FILE * fp=stdout);

// Binary printing methods
void printBinary(uint32_t value, const char * zero="0", const char * one="1", int msb=32);
void printBinary(uint64_t value, const char * zero="0", const char * one="1", int msb=64);
void printBinary(float value, const char * zero="0", const char * one="1", int msb=32);
void printBinary(void * value32, const char * zero="0", const char * one="1", int msb=32);

template<class T> void print2D(T* pix, int nx, int ny, FILE * fp=stdout);

/// Print signed normalized value on a horizontal plot.

/// @param[in]	value	Normalized value to plot
/// @param[in]	width	Character width of plot excluding center point
/// @param[in]	spaces	Print extra filling spaces to the right
/// @param[in]	point	The print character for points
void printPlot(float value, uint32_t width=50, bool spaces=true, const char * point="o");





// Implementation

inline char intensityToASCII(float v){
	static const char map[] =
	" .,;-~_+<>i!lI?/|)(1}{][rcvunxzjftLCJUYXZO0Qoahkbdpqwm*WMB8&%$#@";
	//"$@B%8&WM#*oahkbdpqwmZO0QLCJUYXzcvunxrjft/\|()1{}[]?-_+~<>i!lI;:,\"^`'. ";
//	 123456789.123456789.123456789.123456789.123456789.123456789.1234
	static const int N = sizeof(map)-1;
	return map[int(N*scl::clip(v,0.9999999f))];
}

#define DEF(type, spec)\
template<>\
inline void print<type>(const type& v, const char * post, const char * pre, FILE * fp){\
	fprintf(fp, "%s%"#spec"%s", pre, v, post);\
}
DEF(float, f) DEF(double, f) DEF(uint32_t, d) DEF(int, d)
#undef DEF

template<class T> void print2D(T* pix, int nx, int ny, FILE * fp=stdout){
	for(int j=0; j<nx; ++j){
	for(int i=0; i<ny; ++i){
		float v = pix[j*nx + i];
		fprintf(fp, "%c ", intensityToASCII(v));
	} printf("\n"); }
}

} // gam::

#endif
