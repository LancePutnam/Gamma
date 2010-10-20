/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Example:		Generator / Modulation / Frequency
	Description:	Single carrier and modulator frequency modulation.
*/

#include "../examples.h"

float ff = 220;					// fundamental frequency
Accum<> tmr(1./2., 2);
Sine<> oscC(ff);				// carrier
Sine<> oscM;					// modulator
LFO<> modIndex(tmr.freq());		// envelope for changing modulator frequency

float modFreq = 1;				// modulator frequency as multiple of carrier frequency

void audioCB(AudioIOData& io){

	while(io()){

		if(tmr()){
			printf("c:m ratio = 1:% 5.2f\n", modFreq);
			oscM.freq(ff*modFreq + 1);		// change modulator frequency
			modFreq+=0.25;					// increment modulator frequency
			if(modFreq > 4) modFreq = 1;	// wrap modulator frequency
		}

		oscC.freq(ff + oscM()*modIndex.hann()*200);		// modulate frequency of carrier

		float s = oscC() * 0.2;
		
		io.out(0) = io.out(1) = s;
	}
}

RUN(audioCB);
