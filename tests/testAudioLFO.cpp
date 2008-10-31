#include <stdio.h>
#include "Gamma.h"
#include "AudioIO.h"
#include "Noise.h"
#include "Oscillator.h"

using namespace gam;

Accum<> tmr(0.5);
NoisePink<> noise;
LFO<> lfo(2,0,0.25);
gen::Trigger trig(7);

void audioCB(AudioIOData & io){

	for(uint32_t i=0; i<io.numFrames(); i++){
	
		if(tmr()) trig();
	
		float s = 0.f;
		
		#define CS(v,f) case v: s = lfo.f(); break;
		switch(trig.val){
			CS(0, upU) CS(1, downU) CS(2, triU) CS(3, cosU) CS(4, sqrU)
			CS(5, pulseU) CS(6, imp)
		}
		
		io.out(0)[i] = io.out(1)[i] = s * noise() * 0.4f;
	}
}

int main(int argc, char* argv[]){
	AudioIO io(256, 44100., audioCB, NULL, 2);
	Sync::master().spu(io.fps());
	
	io.start();
	printf("\nPress 'enter' to quit...\n"); getchar();
	return 0;
}
