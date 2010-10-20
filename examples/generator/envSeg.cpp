/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Example:		Generator / Envelope
	Description:	This uses a linear envelope segment to control the beating 
					of a simple additive synth.
*/

#include "../examples.h"

float fund = 440;					// fundamental frequency
Accum<> tmr(1/1.2, 2);				// timer for resetting envelope
Sine<> osc1(fund), osc2(fund*2);	// source oscillators
Seg<> env(1, fund*2, fund*2);		// beat envelope

void audioCB(AudioIOData& io){

	while(io()){
	
		if(tmr()){
			env = fund*2 + rnd::uni(10.0);	// set new target value of envelope
		}
		
		osc2.freq(env());					// modulate frequency of 2nd harmonic
		float s = (osc1() + osc2()) * 0.1;

		io.out(0) = io.out(1) = s;
	}
}


RUN(audioCB);
