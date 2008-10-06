#include <stdio.h>
#include "Gamma.h"
#include "AudioIO.h"
#include "Noise.h"
#include "Oscillator.h"

using namespace gam;

NoiseWhite<> white;
NoisePink<> pink;
NoiseBrown<> brown;
Accum<> tmr(0.5);
gen::Counter cnt(3);

void audioCB(AudioIOData & io){
	float * out0 = io.out(0);
	float * out1 = io.out(1);
	
	for(uint32_t i=0; i<io.numFrames(); i++){
	
		if(tmr()) cnt();
	
		float s = 0.f;
		
		switch(cnt.val){
			case 0: s = white(); break;
			case 1: s = pink();  break;
			case 2: s = brown(); break;
		}
	
		out0[i] = out1[i] = s * 0.2;
	}
}

int main(int argc, char* argv[]){
	AudioIO io(256, 44100., audioCB, NULL);
	Sync::master().spu(io.fps());
	io.start();
	printf("\nPress 'enter' to quit...\n"); getchar();
	return 0;
}
