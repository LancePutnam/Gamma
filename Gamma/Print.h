#ifndef GAMMA_PRINT_H_INC
#define GAMMA_PRINT_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information

	File Description: 
	Printing functions
*/

#include <stdio.h>
#include <string>
#include "Gamma/Constants.h"
#include "Gamma/scl.h"

namespace gam{

/// Returns an ASCII character most closely matching an intensity value in [0,1].
char intensityToASCII(float v);

std::string plotString(
	float value, unsigned width=50, bool spaces=true, bool sign=true, const char * point="o"
);


/// Prints 2D pixel array
template<class T>
void print2D(T* pix, unsigned nx, unsigned ny, FILE * fp=stdout);

// Binary printing methods
void printBinary(uint32_t value, const char * zero="0", const char * one="1", int msb=32);
void printBinary(uint64_t value, const char * zero="0", const char * one="1", int msb=64);
void printBinary(float value, const char * zero="0", const char * one="1", int msb=32);
void printBinary(const void * value32, const char * zero="0", const char * one="1", int msb=32);

/// Prints array as hexidecimal values
void printHexArray(const float * table, unsigned len, unsigned valuesPerLine);

/// Print signed unit value on a horizontal plot

/// \param[in]	value	Normalized value to plot
/// \param[in]	width	Character width of plot excluding center point
/// \param[in]	spaces	Print extra filling spaces to the right
/// \param[in]	sign	Whether plot is signed
/// \param[in]	point	The print character for points
void printPlot(float value, unsigned width=50, bool spaces=true, bool sign=true, const char * point="o");

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

template<class T>
void print2D(T* pix, unsigned nx, unsigned ny, FILE * fp){
	for(unsigned j=0; j<nx; ++j){
		for(unsigned i=0; i<ny; ++i){
			float v = pix[j*nx + i];
			fprintf(fp, "%c ", intensityToASCII(v));
		}
		fprintf(fp, "\n");
	}
}

} // gam::

#endif
