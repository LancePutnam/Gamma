/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Example:		Generator / Envelope
	Description:	Using an exponentially decaying envelope to control the
					amplitude of a noise source.
*/

#include "../examples.h"

Accum<> tmr(1,1);			// Timer for resetting envelope
NoiseWhite<> src;		// Noise source

SegExp<> env(	// Exponential envelope
	1,			// Length
	0.2,		// Curvature
	0.2,0		// Start/end values
);

float curvature = 10;

void audioCB(AudioIOData& io){

	while(io()){
	
		if(tmr()){
			// reset envelope to beginning
			env.reset();
			
			// set a new curvature value
			if(--curvature < -10) curvature+=20;
			env.curve(curvature);
			printf("curvature = % g\n", curvature);
		}

		float s = src() * env();

		io.out(0) = io.out(1) = s;
	}
}

RUN(audioCB);
