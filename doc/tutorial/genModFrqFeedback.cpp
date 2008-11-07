/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Tutorial:		Generator / Modulation / Frequency
	Description:	Single-oscillator feedback frequency modulation.
*/

#include "tutorial.h"

Sine<> osc(220);			// source sine
LFO<> mod(0.1);		// modulator of modulation index (0.1 hz)
float prevSample = 0;

void audioCB(AudioIOData& io){

	for(uint32_t i=0; i<io.numFrames(); i++){
	
//		float s = 0.f;			// initialize current sample to zero
		osc.phaseAdd(prevSample*0.002*mod.hann());
		
		float s = osc();							// get current sample
		prevSample = s;
		io.out(0)[i] = io.out(1)[i] = s * 0.4f;
	}
}


int main(int argc, char* argv[]){
	AudioIO io(256, 44100., audioCB, NULL, 2);
	Sync::master().spu(io.fps());
	
	io.start();
	printf("\nPress 'enter' to quit...\n"); getchar();
	return 0;
}
