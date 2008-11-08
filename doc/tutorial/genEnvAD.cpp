/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Tutorial:		Generator / Envelope
	Description:	Using an exponential AD (attack/decay) envelope to control
					the amplitude of an oscillator.
*/

#include "tutorial.h"

Accum<> tmr(1/1.2);			// timer for resetting envelope
Sine<> osc(440);			// source
AD<> env(0, 1);				// attack/decay envelope
double tilt = 0;			// tilt of envelope; 0=percussive, 1=reversive

void audioCB(AudioIOData& io){

	for(uint32_t i=0; i<io.numFrames(); i++){
	
		if(tmr()){
			tilt += 0.1; if(tilt > 1) tilt=0;	// increment tilt amount
			env.set(tilt, 1-tilt);				// set new attack/decay times
			env.reset();						// reset amplitude of envelope
		}
		
		float s = osc() * env() * 0.2;

		io.out(0)[i] = io.out(1)[i] = s;
	}
}


int main(int argc, char* argv[]){
	AudioIO io(256, 44100., audioCB, NULL, 2);
	Sync::master().spu(io.fps());
	
	io.start();
	printf("\nPress 'enter' to quit...\n"); getchar();
	return 0;
}
