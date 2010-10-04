/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Example:		Generator / Oscillator / Buzz
	Description:	Variable harmonic impulse
*/

#include "../examples.h"

Buzz<> buzz(55);			// Linearly sweeping phase
LFO<> mod(1./16);

void audioCB(AudioIOData& io){

	for(int i=0; i<io.framesPerBuffer(); ++i){
	
		int nh = mod.triU() * 64;
		buzz.harmonics(nh);
	
		float s = buzz();

		io.out(0)[i] = io.out(1)[i] = s*0.2;
	}
}

RUN(audioCB);
