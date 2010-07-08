/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Tutorial:		Generator / Oscillator / Impulse
	Description:	Band-limited impulse
*/

#include "tutorial.h"

Impulse<> impulse(55);			// Linearly sweeping phase
LFO<> mod(1./16);

void audioCB(AudioIOData& io){

	for(uint32_t i=0; i<io.framesPerBuffer(); i++){
	
		int nh = mod.triU() * 32;
		impulse.harmonics(nh);
	
		float s = impulse();

		io.out(0)[i] = io.out(1)[i] = s*0.2;
	}
}

RUN(audioCB);
