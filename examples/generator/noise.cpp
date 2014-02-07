/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Example:		Generator / Noise
	Description:	The different colors of noise. Each of the noise types has
					a frequency-dependent power spectrum according to 
					power(f) = 1/f^n. For white noise, n=0 resulting in a flat
					spectrum. For pink and brown noise, n=1 and n=2, 
					respectively producing progressively darker sounding noise.
*/

#include "../examples.h"

NoiseWhite<> white;		// 1/f^0 noise
NoisePink<> pink;		// 1/f^1 noise
NoiseBrown<> brown;		// 1/f^2 noise
int type = 2;			// Noise type
Accum<> tmr(1./2);		// Timer for switching noise

void audioCB(AudioIOData& io){

	while(io()){
	
		if(tmr()) (++type)%=3;
		
		float s = 0;
		
		switch(type){
			case 0: s = white()*0.4;break;
			case 1: s = pink();		break;
			case 2: s = brown();	break;
		}

		s *= 0.2;

		io.out(0) = s;
		io.out(1) = s;
	}
}

RUN_AUDIO_MAIN
