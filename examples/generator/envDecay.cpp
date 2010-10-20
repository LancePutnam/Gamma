/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Example:		Generator / Envelope
	Description:	Using an exponentially decaying envelope to control the
					amplitude of a noise source.
*/

#include "../examples.h"

Accum<> tmr(1);			// Timer for resetting envelope
NoiseWhite<> src;		// Noise source
Decay<> env(1, 0.2);	// Envelope starting at 0.2 with duration of 1 second.
						// The duration is the amount of time it takes the
						// envelope to decay by -60 dB or 1/1000th of its 
						// starting value.

void audioCB(AudioIOData& io){

	while(io()){
	
		if(tmr()){
			env.value(0.2);		// reset amplitude of envelope
		}

		float s = src() * env();

		io.out(0) = io.out(1) = s;
	}
}

RUN(audioCB);
