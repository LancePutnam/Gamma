/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Example:		Generator / Envelope
	Description:	Using an exponential AD (attack/decay) envelope to control
					the amplitude of an oscillator.
*/

#include "../examples.h"

Accum<> tmr(1/1.2);			// Timer for resetting envelope
NoiseWhite<> src;			// Noise source
AD<> env(0, 1);				// Attack/decay envelope
double tilt = 0;			// Tilt of envelope; 0=percussive, 1=reversive

void audioCB(AudioIOData& io){

	while(io()){
	
		if(tmr()){
			tilt += 0.1; if(tilt > 1) tilt=0;	// increment tilt amount
			env.set(tilt, 1-tilt);				// set new attack/decay times
			env.reset();						// reset amplitude of envelope
		}

		float s = src() * env() * 0.2;

		io.out(0) = io.out(1) = s;
	}
}

RUN(audioCB);
