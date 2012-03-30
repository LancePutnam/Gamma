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

	double seq[] = {0, -1, 1, -1, 0};
	double step = 1./16;
	
	#define DO(fnc)\
		printf("\n"#fnc":\n");\
		for(double phase=1; phase<3; phase+=step){\
			int i1 = int(phase);\
			double a = seq[i1-1];\
			double b = seq[i1];\
			double c = seq[i1+1];\
			double d = seq[i1+2];\
			a=a; d=d;\
			double fr = phase - i1;\
			double v = fnc;\
			printf("\t% 6.3f  ", v); printPlot(v, 32); printf("\n");\
		}\
	
	DO(ipl::nearest		(fr,   b,c));
	DO(ipl::linear		(fr,   b,c));
	DO(ipl::quadratic	(fr,   b,c,d));
	DO(ipl::cubic		(fr, a,b,c,d));
	return 0;
}

