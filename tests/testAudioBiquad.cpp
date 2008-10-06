#include <stdio.h>
#include "Gamma.h"
#include "AudioIO.h"
#include "Delay.h"
#include "Noise.h"
#include "Oscillator.h"
#include "Sync.h"

using namespace gam;

NoiseWhite<> src;
Biquad<> bq(800, 4);
LFO<> lfoCtf(0.1);
Accum<> tmr(lfoCtf.freq());

void audioCB(AudioIOData & io){
	float * out0 = io.out(0);
	float * out1 = io.out(1);
	int numFrames = io.numFrames();

	for(int i=0; i<numFrames; ++i){

		// cycle through filter types
		if(tmr()) bq.type((bq.type() + 1) % Filter::AP);

		bq.freq(scl::pow2(lfoCtf.triU()) * 20000.f + 20.f);

		out0[i] = out1[i] = bq(src()) * 0.2f;
	}
}

int main(int argc, char* argv[]){
	AudioIO io(256, 44100., audioCB, NULL, 2);
	Sync::master().spu(io.fps());	
	io.start();
	printf("\nPress 'enter' to quit...\n"); getchar();
	return 0;
}
