/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Example:		Filter / One-pole Smoothing
	Description:	The demonstrates how a one-pole filter can be used to smooth
					out control signals--in this case random stepping of the 
					frequency of an oscillator.
*/

#include "../examples.h"

Accum<> tmr(4);
Sine<> src;
OnePole<> freq(10);

void audioCB(AudioIOData& io){

	while(io()){
		
		if(tmr()){
			freq = pow(2, rnd::uniS(1.))*440;
		}
		
		src.freq(freq());
		
		float s = src();
			
		io.out(0) = io.out(1) = s * 0.2f;
	}
}

RUN_AUDIO_MAIN
