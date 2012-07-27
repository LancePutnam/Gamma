/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Example:		Generator / Modulation / Frequency
	Description:	Transition from vibratory to timbral frequency modulation.
*/

#include "../examples.h"

float ff = 220;			// fundamental frequency
Sine<> oscC(ff);		// carrier
Sine<> oscM;			// modulator
LFO<> modFreq(1./20.);	// envelope for changing modulator frequency

void audioCB(AudioIOData& io){

	while(io()){

		oscM.freq(modFreq.hann() * 110 + 1);	// change modulator frequency
		oscC.freq(ff + oscM()*100);				// modulate frequency of carrier

		float s = oscC() * 0.2;
		
		io.out(0) = io.out(1) = s;
	}
}

RUN_AUDIO_MAIN
