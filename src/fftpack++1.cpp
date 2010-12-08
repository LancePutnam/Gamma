/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include "fftpack++.h"
#include <math.h>

/* These definitions need to be changed depending on the floating-point precision */
#define POST(x)	x##f
#define FUNC(x)	x##1
typedef float real_t;

#include "fftpack++.inc"

#undef POST
#undef FUNC
