/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Tutorial:		Generator / Envelope
	Description:	Using a linear envelope segment to control the beating of a 
					simple additive synth.
*/

#include "tutorial.h"

float fund = 220;					// fundamental frequency
Accum<> tmr(1/1.2, 2);				// timer for resetting envelope
Sine<> osc1(fund), osc2(fund*2);	// source oscillators
Seg<> env(1, fund*2, fund*2);		// beat envelope

void audioCB(AudioIOData& io){

	for(uint32_t i=0; i<io.numFrames(); i++){
	
		if(tmr()){
			env = fund*2 + rnd::uni(10.0);	// set new target value of envelope
		}
		
		osc2.freq(env());					// modulate frequency of 2nd harmonic
		float s = (osc1() + osc2()) * 0.1;

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
