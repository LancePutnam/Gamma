/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Tutorial:		Filter / One-pole filtering
	Description:	Filtering with one-pole
*/

#include "tutorial.h"

LFO<> mod(1./8, 0.5);
NoiseWhite<> src;
OnePole<> filt(10);

void audioCB(AudioIOData& io){

	for(uint32_t i=0; i<io.framesPerBuffer(); i++){
		
		float cutoff = scl::pow3(mod.triU()) * 10000;
		
		filt.freq(cutoff);
		
		float s = filt(src());
			
		io.out(0)[i] = io.out(1)[i] = s * 0.2f;
	}
}

RUN(audioCB);
