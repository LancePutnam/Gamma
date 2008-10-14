/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Tutorial:		Generator/Oscillator/LFO
	Description:	Using an LFO as an amplitude envelope.	
*/

#include "tutorial.h"

Accum<> tmr(0.5);			// timer to switch between LFO types
NoiseWhite<> noise;			// source noise
LFO<> env(5,0,0.25);		// LFO as amplitude envelope
gen::Counter cnt(7);		// counter for LFO type


void audioCB(AudioIOData& io){

	for(uint32_t i=0; i<io.numFrames(); i++){
	
		if(tmr()) cnt();	// increment LFO type
	
		float s = 0.f;

		switch(cnt.val){			
			case  0: s = env.upU(); break;
			case  1: s = env.downU(); break;
			case  2: s = env.triU(); break;
			case  3: s = env.cosU(); break;
			case  4: s = env.sqrU(); break;
			case  5: s = env.pulseU(); break;
			case  6: s = env.imp(); break;
		}
		
		io.out(0)[i] = io.out(1)[i] = s * noise() * 0.4f;
	}
}


int main(int argc, char* argv[]){
	AudioIO io(256, 44100., audioCB, NULL, 2);
	Sync::master().spu(io.fps());
	
	io.start();
	printf("\nPress 'enter' to quit...\n"); getchar();
	return 0;
}
