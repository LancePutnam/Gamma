/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Tutorial:		Generator / Oscillator / Additive
	Description:	Two-oscillator beating
*/

#include "tutorial.h"

float ff = 220;
float freqBeat = 1;

Sine<> osc1(ff);					// fundamental oscillator
Sine<> osc2(ff + freqBeat, 0.5);	// beat oscillator starting 180 out-of-phase
									// from fundamental

void audioCB(AudioIOData& io){

	for(uint32_t i=0; i<io.framesPerBuffer(); i++){

		float s = (osc1() + osc2()) * 0.1;

		io.out(0)[i] = io.out(1)[i] = s;
	}
}

RUN(audioCB);
