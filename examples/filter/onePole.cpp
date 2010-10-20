/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Example:		Filter / One-pole filtering
	Description:	Filtering with one-pole
*/

#include "../examples.h"

LFO<> mod(1./8, 0.5);
NoiseWhite<> src;
OnePole<> filt(10);

void audioCB(AudioIOData& io){

	while(io()){
		
		float cutoff = scl::pow3(mod.triU()) * 10000;
		
		filt.freq(cutoff);
		
		float s = filt(src());
			
		io.out(0) = io.out(1) = s * 0.2f;
	}
}

RUN(audioCB);
