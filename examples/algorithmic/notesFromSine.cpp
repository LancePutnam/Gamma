/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Example:		Algorithmic / Notes From A Sinusoid
	Description:	This is an example of using a high-frequency sine oscillator
					to determine the pitches of notes. This is technique is 
					inspired by the Strang moiré patterns described in:

					Richet, N. (1992). Strang’s strange figures. 
					The American Mathematical Monthly, 99(2):101-107.
*/

#include "../examples.h"

Sine<> src;
AD<> env;
double phase=0;
Accum<> tmr;

void audioCB(AudioIOData& io){

	tmr.period(0.15);
	env.attack(0.01).decay(0.14);

	while(io()){

		if(tmr()){
			// Add to the phase a low-integer division of 2π plus a small 
			// offset for variation.
			phase = scl::wrapPhase(phase + 2*M_PI / (rnd::prob(0.9)?3:2) + 0.012);

			float frq = sin(phase) * 8;
			
			// Map sine value to nearest pitch in C-major scale
			frq = scl::ratioET(scl::nearest(frq, "2212221")) * scl::freq("c4");

			src.freq(frq);
			env.reset();
		}

		float s = src() * env() * 0.2;

		io.out(0) = s;
		io.out(1) = s;
	}
}

RUN_AUDIO_MAIN
