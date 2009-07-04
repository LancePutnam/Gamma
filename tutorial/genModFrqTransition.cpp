/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Tutorial:		Generator / Modulation / Frequency
	Description:	Transition from vibratory to timbral frequency modulation.
*/

#include "tutorial.h"

float ff = 220;			// fundamental frequency
Sine<> oscC(ff);		// carrier
Sine<> oscM;			// modulator
LFO<> modFreq(1./20.);	// envelope for changing modulator frequency

void audioCB(AudioIOData& io){

	for(uint32_t i=0; i<io.framesPerBuffer(); i++){

		oscM.freq(modFreq.hann() * 110 + 1);	// change modulator frequency
		oscC.freq(ff + oscM()*100);				// modulate frequency of carrier

		//oscC.phaseAdd(oscM()*0.002);

		float s = oscC() * 0.2;
		
		io.out(0)[i] = io.out(1)[i] = s;
	}
}


int main(int argc, char* argv[]){
	AudioIO io(256, 44100., audioCB, NULL, 2);
	Sync::master().spu(io.framesPerSecond());
	
	io.start();
	printf("\nPress 'enter' to quit...\n"); getchar();
	return 0;
}
