/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Example:		Effect / Echo
	Description:	This demonstrates how to create an echo effect from a
					delay line. An echo is simply a feedback delay line.
*/

#include "../examples.h"

Delay<> delay(0.2);				// 200 ms delay
Accum<> tmr(1./4);
SineD<> src(1000, 0.1, 0.1);	// short sine plink

void audioCB(AudioIOData& io){

	while(io()){

		if(tmr()) src.reset();

		float s = src();

		// Here is the feedback delay:
		// We read the oldest sample from the delay line, scale it, and feed
		// it back into the delay.
		s += delay(s + delay()*0.7);
	
		io.out(0) = io.out(1) = s;
	}
}

RUN_AUDIO_MAIN
