/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Example:		Generator / Oscillator / Additive
	Description:	This demonstrates the binaural beating effect produced
					from playing two sinusoids with slightly different
					frequencies in each ear.
*/

#include "../examples.h"

float ff = 220;
float freqBeat = 1;

Sine<> osc1(ff);				// "left" sinusoid
Sine<> osc2(ff + freqBeat);		// "right" sinusoid

void audioCB(AudioIOData& io){

	while(io()){

		io.out(0) = osc1() * 0.1;
		io.out(1) = osc2() * 0.1;
	}
}


RUN_AUDIO_MAIN
