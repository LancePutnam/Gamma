/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Tutorial:		Generator / Noise
	Description:	The different colors of noise.
*/

#include "tutorial.h"

Accum<> tmr(0.5);
NoiseWhite<> white;		// 1/f0 noise
NoisePink<> pink;		// 1/f1 noise
NoiseBrown<> brown;		// 1/f2 noise

gen::Trigger cnt(3);


void audioCB(AudioIOData& io){

	for(uint32_t i=0; i<io.framesPerBuffer(); i++){
		
		if(tmr()) cnt();
		
		float s = 0;
		
		switch(cnt.val){
			case 0: s = white();	break;
			case 1: s = pink();		break;
			case 2: s = brown();	break;
		}

		io.out(0)[i] = io.out(1)[i] = s*0.2;
	}
}


int main(int argc, char* argv[]){

	AudioIO io(256, 44100, audioCB, NULL, 2);
	Sync::master().spu(io.framesPerSecond());

	io.start();
	printf("\nPress 'enter' to quit...\n"); getchar();
	return 0;
}
