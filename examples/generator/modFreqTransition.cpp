/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Example:		Generator / Modulation / Frequency
	Description:	Transition from vibratory to timbral frequency modulation.
*/

#include "../examples.h"

float ff = 220;			// fundamental frequency
Sine<> oscC(ff);		// carrier
Sine<> oscM;			// modulator
LFO<> modFreq(1./20.);	// envelope for changing modulator frequency

void audioCB(AudioIOData& io){

	while(io()){

		oscM.freq(modFreq.hann() * 110 + 1);	// change modulator frequency
		oscC.freq(ff + oscM()*100);				// modulate frequency of carrier

		//oscC.phaseAdd(oscM()*0.002);

		float s = oscC() * 0.2;
		
		io.out(0) = io.out(1) = s;
	}
}


int main(int argc, char* argv[]){
	AudioIO io(256, 44100., audioCB, NULL, 2);
	Sync::master().spu(io.framesPerSecond());
	
	io.start();
	printf("\nPress 'enter' to quit...\n"); getchar();
	return 0;
}
