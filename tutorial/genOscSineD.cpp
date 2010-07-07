/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Tutorial:		Generator / Oscillator / Damped Sine
	Description:	Decaying sinusoid
*/

#include "tutorial.h"

Accum<> tmr(2, 2);			// Timer for randomizing parameters
SineD<> osc(440, 0.2, 1);	// Decaying sinusoid

void audioCB(AudioIOData& io){

	for(uint32_t i=0; i<io.framesPerBuffer(); i++){
		
		// Randomize frequency and decay time every so often...
		if(tmr())
			osc.set(rnd::uni(10, 1)*50, 0.2, rnd::lin(2., 0.1));
		
		float s = osc();

		io.out(0)[i] = io.out(1)[i] = s;
	}
}


RUN(audioCB)
