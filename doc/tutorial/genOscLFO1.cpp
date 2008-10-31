/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Tutorial:		Generator/Oscillator/LFO
	Description:	Using an LFO as an oscillator source.
*/

#include "tutorial.h"

Accum<> tmr(0.4);			// timer to switch between LFO types
LFO<> osc(220, 0, 0.25);	// source (220 hz, 0 phase, 0.25 mod)
LFO<> mod(0.5);				// modulator of source's modifier parameter
gen::Trigger cnt(11);		// counter for LFO type


void audioCB(AudioIOData& io){

	for(uint32_t i=0; i<io.numFrames(); i++){
	
		if(tmr()) cnt();		// increment LFO type

		osc.mod(mod.cosU());	// modulate modifier parameter
	
		float s = 0.f;

		// An LFO is a subclass of a unsigned 32 bit phase accumulator.
		// It's output is determined by a series of named generation
		// methods which map the internal integer phase to a floating
		// point value according to various strategies. This allows 
		// the user to quickly and efficiently switch different LFOs
		// generator methods at run-time using the same object instance.
		
		switch(cnt.val){ 
				
			// non-modifiable ordered from dull to bright:
			case  0: s = osc.cos(); break;		// Cosine based on 3rd order polynomial. When phase is zero, value is 1
			case  1: s = osc.even3(); break;    // Even harmonic sine-like wave (3rd order)
			case  2: s = osc.even5(); break;	// Even harmonic sine-like wave (5th order)
			case  3: s = osc.tri(); break;		// Triangle (starts at 1 going down to -1 then up to 1)
			case  4: s = osc.sqr(); break;		// Square (-1 to 1)
			case  5: s = osc.down(); break;		// Downward ramp (1 to -1)
			case  6: s = osc.imp(); break;		// Impulse (1 occurs at beginning of cycle, 0 otherwise)
			
			// modifiable ordered from dull to bright:
			case  7: s = osc.stair(); break;	
			case  8: s = osc.pulse(); break;	// Pulse (up + down). 'mod' controls pulse width
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
