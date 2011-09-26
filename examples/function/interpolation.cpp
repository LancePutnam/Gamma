/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Example:		Function / Interpolation
	Description:	Demonstration of various interpolation functions.
*/

#include <stdio.h>
#include "Gamma/Gamma.h"
#include "Gamma/Print.h"

using namespace gam;

int main(){

	const int len = 16;
	gen::RAdd<> frac(1./len);

	#define LOOP for(int i=0; i<len; i++)
	#define DO(fnc)\
		frac = -frac.add; printf("\n"#fnc":\n");\
		LOOP{\
			float v = fnc;\
			printf("\t[%3d] % 6.3f  ", i, v); printPlot(v, 32); printf("\n");\
		}
	
	DO(ipl::nearest(frac(), -1., 1.))	
	DO(ipl::linear(frac(), -1., 1.))
	DO(ipl::quadratic(frac(), -1., 1., 1.))
	DO(ipl::cubic(frac(), 0., -1., 1., 0.))
	DO(ipl::cubic2(frac(), 0., -1., 1., 0.))
	
	return 0;
}

