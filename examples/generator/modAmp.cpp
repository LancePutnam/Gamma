/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Example:		Generator / Modulation / Amplitude
	Description:	Transition from undulating to timbral amplitude modulation.
*/

#include "../examples.h"

float ff = 440;			// fundamental frequency
Sine<> oscC(ff);		// carrier
Sine<> oscM;			// modulator
LFO<> modFreq(1./20.);	// envelope for changing modulator frequency

void audioCB(AudioIOData& io){

	while(io()){

		oscM.freq(modFreq.hann() * 110 + 1);	// change modulator frequency

		float s = oscC() * (oscM() * 0.5 + 0.5) * 0.2;
		
		io.out(0) = io.out(1) = s;
	}
}

RUN(audioCB);
