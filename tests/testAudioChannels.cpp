#include <stdio.h>
#include "Gamma/Gamma.h"
#include "Gamma/AudioIO.h"
#include "Gamma/Delay.h"
#include "Gamma/Noise.h"
#include "Gamma/Oscillator.h"
#include "Gamma/Sync.h"

using namespace gam;

int chan=-1;
Accum<> tmr(2, 2);
NoisePink<> src;
SineD<> osc(2000, 0.1, 0.1);

void audioCB(AudioIOData & io){
	
	while(io()){

		if(tmr()){
			++chan;
			if(chan >= io.channelsOut()) chan=0;
			osc.reset();
			printf("chan: %d\n", chan);
		}

		io.out(chan) = osc();
	}
}

int main(int argc, char* argv[]){
	
	int maxOChans = AudioDevice::defaultOutput().channelsOutMax();
	int maxIChans = AudioDevice::defaultOutput().channelsOutMax();
	//printf("%d %d\n", maxIChans, maxOChans);

	AudioIO io(256, 44100., audioCB, NULL, maxOChans, maxIChans);
	Sync::master().spu(io.framesPerSecond());
	io.start(); io.print();
	
	printf("\nPress 'enter' to quit...\n"); getchar();
	return 0;
}
