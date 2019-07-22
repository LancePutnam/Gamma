/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example:		Generator / Oscillator Plots
Description:	This displays plots of the various oscillators.
*/

#include "Gamma/Noise.h"
#include "Gamma/Print.h"
using namespace gam;

int main(){

	const int N = 32;		// Number of samples per unit of position

	#define PLOT(Class, setup, func, description){\
		Class g;\
		setup;\
		printf("\n" #Class "::" #func ": %s\n", description);\
		for(int i=0; i<N; ++i){\
			float v = g.func();\
			printf("[%2u] % 6.3f  ", i, v);\
			printPlot(v, 32, false);\
			printf("\n");\
		}\
	}

	PLOT(NoiseWhite<>, , operator(), "")
	PLOT(NoisePink<>, , operator(), "")
	PLOT(NoiseBrown<>, , operator(), "")
	PLOT(NoiseViolet<>, , operator(), "")
	PLOT(NoiseBinary<>, , operator(), "")
}
