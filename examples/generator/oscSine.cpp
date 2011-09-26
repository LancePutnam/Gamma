/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Example:		Generator / Oscillator / Sine
	Description:	Plays a sinusoid
*/

#include <stdio.h>
#include "Gamma/AudioIO.h"
#include "Gamma/Oscillator.h"

using namespace gam;

Sine<> osc(440);

void audioCB(AudioIOData& io){
	while(io()){
		float s = osc() * 0.2;
		io.out(0) = s;
		io.out(1) = s;
	}
}

int main(){
	AudioIO io(256, 44100., audioCB);
	Sync::master().spu(io.framesPerSecond());
	io.start();
	printf("\nPress 'enter' to quit...\n"); getchar();
	return 0;
}
