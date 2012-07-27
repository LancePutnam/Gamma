/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Example:		Generator / Oscillator / Buzz
	Description:	Buzz is a variable harmonic impulse train. Note that
					changing the number of harmonics during playback will 
					usually create a click.
*/

#include "../examples.h"

Buzz<> buzz(55);			// Linearly sweeping phase
LFO<> mod(1./16);

void audioCB(AudioIOData& io){

	while(io()){
	
		int nh = mod.triU() * 64;
		buzz.harmonics(nh);
	
		float s = buzz();

		io.out(0) = io.out(1) = s*0.2;
	}
}

RUN_AUDIO_MAIN
