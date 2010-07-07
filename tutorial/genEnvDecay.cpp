/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Tutorial:		Generator / Envelope
	Description:	Using an exponentially decaying envelope to control the
					amplitude of a noise source.
*/

#include "tutorial.h"

Accum<> tmr(1);			// Timer for resetting envelope
NoiseWhite<> src;		// Noise source
Decay<> env(1, 0.2);	// Envelope starting at 0.2 with duration of 1 second.
						// The duration is the amount of time it takes the
						// envelope to decay by -60 dB or 1/1000th of its 
						// starting value.

void audioCB(AudioIOData& io){

	for(uint32_t i=0; i<io.framesPerBuffer(); i++){
	
		if(tmr()){
			env.value(0.2);		// reset amplitude of envelope
			osc.phase(0);		// reset phase of oscillator to avoid click
		}

		float s = src() * env();

		io.out(0)[i] = io.out(1)[i] = s;
	}
}

RUN(audioCB);
