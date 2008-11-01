/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Tutorial:		Generator/Oscillator/LFO
	Description:	Using an LFO as an oscillator source.
*/

#include "tutorial.h"

Accum<> tmr(0.4);			// timer to switch between LFO types
LFO<> osc(220, 0, 0.25);	// source (220 hz, 0 phase, 0.25 mod)
LFO<> mod(0.5);				// modulator of source's modifier parameter
gen::Trigger cnt(12);		// counter for LFO type

void audioCB(AudioIOData& io){

	for(uint32_t i=0; i<io.numFrames(); i++){
	
		if(tmr()) cnt();		// increment LFO type

		osc.mod(mod.cosU());	// modulate modifier parameter with unipolar cosine wave
	
		float s = 0.f;			// initialize current sample to zero
		
		switch(cnt.val){ 
				
			// non-modifiable generators ordered from dull to bright:
			case  0: s = osc.cos();		break;		// Single harmonic
			case  1: s = osc.even3();	break;		// Even harmonic sine-like wave (3rd order)
			case  2: s = osc.even5();	break;		// Even harmonic sine-like wave (5th order)
			case  3: s = osc.tri();		break;		// 1/f^2 odd harmonics
			case  4: s = osc.sqr();		break;		// 1/f odd harmonics
			case  5: s = osc.down();	break;		// 1/f all harmonics
			case  6: s = osc.up();		break;		// 1/f all harmonics
			case  7: s = osc.imp();		break;		// flat spectrum all harmonics
			
			// modifiable generators ordered from dull to bright:
			case  8: s = osc.stair();	break;		// Mix between a square and impulse. 
			case  9: s = osc.pulse();	break;		// Mix between up and down
			case 10: s = osc.line2();	break;		// Mix between a saw and a triangle
			case 11: s = osc.up2();		break;		// Mix between two saws
		}
		
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
