/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Example:		Filter / One-pole filtering
	Description:	This demonstrates the effect of a one-pole low-pass filter
					on a noise source. One-pole filters are very effective
					at controlling the brightness of sounds as they have an
					adjustable cutoff frequency and gentle slope.
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

RUN_AUDIO_MAIN
