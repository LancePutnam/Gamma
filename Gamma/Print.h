#ifndef GAMMA_VISUAL_H_INC
#define GAMMA_VISUAL_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information

	File Description: 
	Visual functions for printing and colors.
*/

#include <string>
#include <stdio.h>
#include "Gamma/Constants.h"
#include "Gamma/Types.h"
#include "Gamma/scl.h"

namespace gam{

void colorRGB(float h, float s, float v, float &r, float &g, float &b);
void colorHSV(float r, float g, float b, float &h, float &s, float &v);

/// Returns an ASCII character most closely matching an intensity value in [0,1].
char intensityToASCII(float v);

template<class T> void print(const T& v, const char * post="\n", const char * pre="", FILE * fp=stdout);

/// Prints 2D pixel array
template<class T> void print2D(T* pix, int nx, int ny, FILE * fp=stdout);

// Binary printing methods
void printBinary(uint32_t value, const char * zero="0", const char * one="1", int msb=32);
void printBinary(uint64_t value, const char * zero="0", const char * one="1", int msb=64);
void printBinary(float value, const char * zero="0", const char * one="1", int msb=32);
void printBinary(void * value32, const char * zero="0", const char * one="1", int msb=32);

/// Prints array as hexidecimal values.
void printHexArray(float * table, uint32_t len, uint32_t valuesPerLine);

/// Print signed unit value on a horizontal plot.

/// @param[in]	value	Normalized value to plot
/// @param[in]	width	Character width of plot excluding center point
/// @param[in]	spaces	Print extra filling spaces to the right
/// @param[in]	point	The print character for points
void printPlot(float value, uint32_t width=50, bool spaces=true, const char * point="o");

/// Prints error messge to stderr and optionally calls exit()
void err(const char * msg, const char * src="", bool exits=true);

/// Prints warning messge to stderr
void warn(const char * msg, const char * src="");




// Implementation

inline char intensityToASCII(float v){
	static const char map[] =
	" .,;-~_+<>i!lI?/|)(1}{][rcvunxzjftLCJUYXZO0Qoahkbdpqwm*WMB8&%$#@";
	//"$@B%8&WM#*oahkbdpqwmZO0QLCJUYXzcvunxrjft/\|()1{}[]?-_+~<>i!lI;:,\"^`'. ";
//	 123456789.123456789.123456789.123456789.123456789.123456789.1234
	static const int N = sizeof(map)-1;
	return map[int(N*scl::clip(v,0.9999999f))];
}

#define DEF(type, spec, val)\
template<>\
inline void print<type>(const type& v, const char * post, const char * pre, FILE * fp){\
	fprintf(fp, "%s%"#spec"%s", pre, val, post);\
}
DEF(float, g, v) DEF(double, g, v) DEF(uint32_t, d, v) DEF(int, d, v)
DEF(std::string, s, v.c_str())
#undef DEF

template<class T> void print2D(T* pix, int nx, int ny, FILE * fp){
	for(int j=0; j<nx; ++j){
	for(int i=0; i<ny; ++i){
		float v = pix[j*nx + i];
		fprintf(fp, "%c ", intensityToASCII(v));
	} printf("\n"); }
}

} // gam::

#endif
