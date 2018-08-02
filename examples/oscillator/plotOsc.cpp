/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example:		Generator / Oscillator Plots
Description:	This displays plots of the various oscillators.
*/

#include "Gamma/Oscillator.h"
#include "Gamma/Print.h"
using namespace gam;

int main(){

	const int N = 32;		// Number of samples per unit of position
	sampleRate(1);

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

	PLOT(Accum<>, g.period(N/2); g.phase(0), operator(), "period=N/2, phase=0")
	PLOT(Accum<>, g.period(N/2); g.phase(0.9999), operator(), "period=N/2, phase=0.9999")

	PLOT(Sine<>, g.period(N);, operator(), "period=N, phase=0")
	PLOT(Sine<>, g.period(N); g.phase(0.25), operator(), "period=N, phase=0.25")

	PLOT(SineR<>, g.set(1./N, 1.0, 0);, operator(), "period=N")
	PLOT(SineR<>, g.set(1./N, 0.5, 0);, operator(), "period=N, amp=0.5")

	PLOT(SineD<>, g.set(4./N, 1.0, 4*N, 0);, operator(), "period=N/4, decay=4N")

	PLOT(Buzz<>, g.period(N); g.harmonics(1), operator(), "period=N, harmonics=1")
	PLOT(Buzz<>, g.period(N); g.harmonics(2), operator(), "period=N, harmonics=2")
	PLOT(Buzz<>, g.period(N); g.harmonics(4), operator(), "period=N, harmonics=4")
	PLOT(Buzz<>, g.period(N); g.harmonics(8), operator(), "period=N, harmonics=8")
	PLOT(Buzz<>, g.period(N); g.harmonics(8); g.normalize(false), operator(), "period=N, harmonics=8, normalize=false")

	PLOT(LFO<>, g.period(N/2);, cos, "period=N/2")
	PLOT(LFO<>, g.period(N/2);, tri, "period=N/2")
	PLOT(LFO<>, g.period(N/2);, sinPara, "period=N/2")
	PLOT(LFO<>, g.period(N/2);, para, "period=N/2")
	PLOT(LFO<>, g.period(N/2);, even3, "period=N/2")
	PLOT(LFO<>, g.period(N/2);, even5, "period=N/2")
	PLOT(LFO<>, g.period(N/2);, up, "period=N/2")
	PLOT(LFO<>, g.period(N/2);, sqr, "period=N/2")
	PLOT(LFO<>, g.period(N/2);, imp, "period=N/2")
	PLOT(LFO<>, g.period(N/2);, S1, "period=N/2")
	PLOT(LFO<>, g.period(N/2);, C2, "period=N/2")
	PLOT(LFO<>, g.period(N/2);, S3, "period=N/2")
	PLOT(LFO<>, g.period(N/2);, C4, "period=N/2")
	PLOT(LFO<>, g.period(N/2);, S5, "period=N/2")
	PLOT(LFO<>, g.period(N/2); g.mod(0.25), pulse, "period=N/2, pulse width=0.25")
	PLOT(LFO<>, g.period(N/2); g.mod(0.25), up2, "period=N/2, mod=0.25")
	PLOT(LFO<>, g.period(N/2); g.mod(0.25), line2, "period=N/2, pulse width=0.25")
	PLOT(LFO<>, g.period(N/2); g.mod(0.25), stair, "period=N/2, mod=0.25")

	PLOT(LFO<>, g.period(N/2);, hann, "period=N/2")
	PLOT(LFO<>, g.period(N/2);, triU, "period=N/2")
	PLOT(LFO<>, g.period(N/2);, paraU, "period=N/2")
	PLOT(LFO<>, g.period(N/2);, upU, "period=N/2")
	PLOT(LFO<>, g.period(N/2);, sqrU, "period=N/2")
	PLOT(LFO<>, g.period(N/2); g.mod(0.25), pulseU, "period=N/2, pulse width=0.25")
	PLOT(LFO<>, g.period(N/2); g.mod(0.25), line2U, "period=N/2, pulse width=0.25")

	PLOT(DWO<>, g.period(N/2);, up, "period=N/2")
	PLOT(DWO<>, g.period(N/2);, sqr, "period=N/2")
	PLOT(DWO<>, g.period(N/2);, para, "period=N/2")
	PLOT(DWO<>, g.period(N/2);, tri, "period=N/2")

}

