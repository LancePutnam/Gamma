/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example:		Generator / Envelope Plots
Description:	This displays plots of the various envelope generators.
*/

#include "Gamma/Envelope.h"
#include "Gamma/Print.h"
using namespace gam;

int main(){

	const int N = 32;		// Number of samples per unit of position
	sampleRate(1);

	#define PLOT(Class, setup, func, description){\
		Class g;\
		setup;\
		printf("\n" #Class ": %s\n", description);\
		for(int i=0; i<N; ++i){\
			float v = g.func();\
			printf("[%2u] % 6.3f  ", i, v);\
			printPlot(v, 32, false, false);\
			printf("\n");\
		}\
	}

	PLOT(Decay<>, g.decay(N);, operator(), "decay=N")
	PLOT(Decay<>, g.decay(N/2);, operator(), "decay=N/2")

	PLOT(AD<>, g.curve(-4); g.attack(N/4); g.decay(3*N/4);, operator(),
		"attack=N/4, decay=3N/4, curve=-4")

	PLOT(AD<>, g.curve(0); g.attack(N/4); g.decay(3*N/4);, operator(),
		"attack=N/4, decay=3N/4, curve=0")

	PLOT(AD<>, g.curve(4); g.attack(N/4); g.decay(3*N/4);, operator(),
		"attack=N/4, decay=3N/4, curve=4")

	PLOT(ADSR<>,
		g.curve(-4); g.attack(N/4); g.decay(N/4); g.sustain(0.5); g.release(N/2); g.sustainDisable(),
		operator(),
		"attack=N/4, decay=N/4, sustain=0.5, release=N/2, curve=-4")

	PLOT(ADSR<>,
		g.curve(0); g.attack(N/4); g.decay(N/2); g.sustain(1); g.release(N/4); g.sustainDisable(),
		operator(),
		"attack=N/4, decay=N/2, sustain=1, release=N/4, curve=0")
}

