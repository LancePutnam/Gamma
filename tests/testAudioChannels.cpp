#include <stdio.h>
#include "Gamma/Gamma.h"
#include "Gamma/AudioIO.h"
#include "Gamma/Delay.h"
#include "Gamma/Noise.h"
#include "Gamma/Oscillator.h"
#include "Gamma/Sync.h"

using namespace gam;

uint32_t chan=-1;
Accum<> tmr(2, 2);
NoisePink<> src;
SineD<> osc(2000, 0.1, 0.1);

void audioCB(AudioIOData & io){

	io.zeroOut();
	
	for(uint32_t i=0; i<io.framesPerBuffer(); ++i){

		if(tmr()){
			++chan;
			if(chan >= io.channelsOut()) chan=0;
			osc.ampPhase(0.1);
			printf("chan: %d\n", chan);
		}

		io.out(chan)[i] = osc();
	}
}

int main(int argc, char* argv[]){
	AudioIO io(256, 44100., audioCB, NULL, 8);
	Sync::master().spu(io.framesPerSecond());
	io.start(); io.print();
	
	printf("\nPress 'enter' to quit...\n"); getchar();
	return 0;
}
