/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Tutorial:		Generator / Noise
	Description:	The different colors of noise.
*/

#include "tutorial.h"

Accum<> tmr(1./2);
NoiseWhite<> white;		// 1/f0 noise
NoisePink<> pink;		// 1/f1 noise
NoiseBrown<> brown;		// 1/f2 noise
int type = 0;			// Noise type


void audioCB(AudioIOData& io){

	for(uint32_t i=0; i<io.framesPerBuffer(); i++){
		
		if(tmr()) (++type)%=3;
		
		float s = 0;
		
		switch(type){
			case 0: s = white();	break;
			case 1: s = pink();		break;
			case 2: s = brown();	break;
		}

		io.out(0)[i] = io.out(1)[i] = s*0.2;
	}
}

RUN(audioCB);
