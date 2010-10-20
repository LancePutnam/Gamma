/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Example:		Generator / Noise
	Description:	The different colors of noise.
*/

#include "../examples.h"

Accum<> tmr(1./2);
NoiseWhite<> white;		// 1/f0 noise
NoisePink<> pink;		// 1/f1 noise
NoiseBrown<> brown;		// 1/f2 noise
int type = 0;			// Noise type


void audioCB(AudioIOData& io){

	while(io()){
		
		if(tmr()) (++type)%=3;
		
		float s = 0;
		
		switch(type){
			case 0: s = white()*0.4;break;
			case 1: s = pink();		break;
			case 2: s = brown();	break;
		}

		io.sum(s*0.2, 0,1);
	}
}

RUN(audioCB);
