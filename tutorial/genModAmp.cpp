/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Tutorial:		Generator / Modulation / Amplitude
	Description:	Transition from undulating to timbral amplitude modulation.
*/

#include "tutorial.h"

float ff = 440;			// fundamental frequency
Sine<> oscC(ff);		// carrier
Sine<> oscM;			// modulator
LFO<> modFreq(1./20.);	// envelope for changing modulator frequency

void audioCB(AudioIOData& io){

	for(uint32_t i=0; i<io.framesPerBuffer(); i++){

		oscM.freq(modFreq.hann() * 110 + 1);	// change modulator frequency

		float s = oscC() * (oscM() * 0.5 + 0.5) * 0.2;
		
		io.out(0)[i] = io.out(1)[i] = s;
	}
}

RUN(audioCB);
