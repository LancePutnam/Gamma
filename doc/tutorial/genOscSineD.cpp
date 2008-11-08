/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Tutorial:		Generator / Oscillator / Damped Sine
	Description:	Damped sinusoial oscillator
*/

#include "tutorial.h"

Accum<> tmr(2, 2);
SineD<> osc(440, 0.2, 1);

void audioCB(AudioIOData& io){

	for(uint32_t i=0; i<io.numFrames(); i++){
		
		// Randomize frequency and decay time every so often...
		if(tmr())
			osc.set(rnd::uni(10, 1)*50, 0.2, rnd::lin(2., 0.1));
		
		float s = osc();

		io.out(0)[i] = io.out(1)[i] = s;
	}
}


int main(int argc, char* argv[]){

	AudioIO io(256, 44100, audioCB, NULL, 2);
	Sync::master().spu(io.fps());

	io.start();
	printf("\nPress 'enter' to quit...\n"); getchar();
	return 0;
}
