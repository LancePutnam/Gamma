/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Example:		Effect / Envelope Follower
	Description:	This demonstrates the EnvFollow object which can be used to
					estimate the amplitude of a signal. The example uses an
					EnvFollow to reset a plucked string whenever its amplitude
					goes below a certain threshold.
*/

#include "../examples.h"

Pluck pluck(440);
EnvFollow<> envFollow;

void audioCB(AudioIOData& io){

	while(io()){

		float s = pluck() * 0.2;

		if(envFollow(s) < 0.001){
			pluck.reset();
			pluck.freq(rnd::uni(1, 20)*100);
		}

		io.out(0) = io.out(1) = s;
	}
}

RUN_AUDIO_MAIN
