/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Example:		Generator / Oscillator / Sine
	Description:	Plays a sinusoid at 440 Hz.
*/

#include "../examples.h"

using namespace gam;

Sine<> osc(440);

void audioCB(AudioIOData& io){
	while(io()){
		float s = osc() * 0.2;
		io.out(0) = s;
		io.out(1) = s;
	}
}

RUN_AUDIO_MAIN
