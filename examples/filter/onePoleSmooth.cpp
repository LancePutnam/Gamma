/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Example:		Filter / One-pole Smoothing
	Description:	Smoothing frequency changes with a one-pole filter
*/

#include "../examples.h"

Accum<> tmr(4);
Sine<> src;
OnePole<> freq(10);

void audioCB(AudioIOData& io){

	for(int i=0; i<io.framesPerBuffer(); ++i){
		
		if(tmr()){
			freq = pow(2, rnd::uniS(1.))*440;
		}
		
		src.freq(freq());
		
		float s = src();
			
		io.out(0)[i] = io.out(1)[i] = s * 0.2f;
	}
}

RUN(audioCB);
