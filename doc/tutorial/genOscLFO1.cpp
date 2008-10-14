/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Tutorial:		Generator/Oscillator/LFO
	Description:	Using an LFO as an oscillator source.
*/

#include "tutorial.h"

Accum<> tmr(0.4);			// timer to switch between LFO types
LFO<> osc(220, 0, 0.25);	// source
LFO<> mod(0.5);				// modulator of source's modifier parameter
gen::Counter cnt(11);		// counter for LFO type


void audioCB(AudioIOData& io){

	for(uint32_t i=0; i<io.numFrames(); i++){
	
		if(tmr()) cnt();		// increment LFO type

		osc.mod(mod.cosU());	// modulate modifier parameter
	
		float s = 0.f;

		switch(cnt.val){
			case  0: s = osc.cos(); break;		// non-modifiable ordered from dull to bright
			case  1: s = osc.even3(); break;
			case  2: s = osc.even5(); break;
			case  3: s = osc.tri(); break;
			case  4: s = osc.sqr(); break;
			case  5: s = osc.down(); break;
			case  6: s = osc.imp(); break;
			
			case  7: s = osc.stair(); break;	// modifiable ordered from dull to bright
			case  8: s = osc.pulse(); break;
			case  9: s = osc.line2(); break;
			case 10: s = osc.up2(); break;
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
