#include <stdio.h>
#include "Gamma.h"
#include "AudioIO.h"
#include "Delay.h"
#include "Noise.h"
#include "Oscillator.h"
#include "Sync.h"

using namespace gam;

int chan=-1;
Accum<> tmr(2, 2);
NoisePink<> src;
SineD<> osc(2000, 0.1, 0.1);

void audioCB(AudioIOData & io){

	io.zeroOut();
	
	for(int i=0; i<io.numFrames(); ++i){

		if(tmr()){
			chan++;
			if(chan >= io.outChans()) chan=0;
			osc.ampPhase(0.1);
			printf("chan: %d\n", chan);
		}

		io.out(chan)[i] = osc();
	}
}

int main(int argc, char* argv[]){
	AudioIO io(256, 44100., audioCB, NULL, 8);
	Sync::master().spu(io.fps());
	io.start(); io.print();
	
	printf("\nPress 'enter' to quit...\n"); getchar();
	return 0;
}
