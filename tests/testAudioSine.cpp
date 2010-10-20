#include <stdio.h>
#include "Gamma/AudioIO.h"
#include "Gamma/Oscillator.h"

using namespace gam;

Sine<> osc(440);
//SineR<> osc(440);
//SineD<> osc(440, 1, 10);

void audioCB(AudioIOData& io){
	while(io()){
		io.out(0) = io.out(1) = osc() * 0.2;
	}
}

int main(int argc, char* argv[]){
	AudioIO io(256, 44100., audioCB, NULL, 2);
	Sync::master().spu(io.framesPerSecond());
	io.start();
	printf("\nPress 'enter' to quit...\n"); getchar();
	return 0;
}
