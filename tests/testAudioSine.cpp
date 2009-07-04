#include <stdio.h>
#include "AudioIO.h"
#include "Oscillator.h"

using namespace gam;

const unsigned long l2Size = 8;
float table[1<<l2Size];

// The many ways to skin a cat...
TableOsc<>   sineT1(table, l2Size, 440);	// table look-up (full-wave)
TableSine<>  sineT2(440);					// table look-up (quarter-wave)
Quadra<>     sineQd(440);					// complex mul
LFO<>        sineC1(440, -0.25f);			// computed 3rd order poly
Sine<>       sineC2(440);					// computed taylor
AccumPhase<> sineC3(440);					// direct math.h sin() (ground truth)

//SineRes< Multi<2> > sineRes(Multi<2>(220, 222));

void audioCB(AudioIOData & io){
	float * out0 = io.out(0);
	float * out1 = io.out(1);
	unsigned long numFrames = io.framesPerBuffer();

	for(unsigned long f=0; f<numFrames; f++){
		float s =	sineT1.nextL() +
					sineT2.nextL() +
					sineQd().i +
					sineC1.cos() +
					sineC2();

		s = s / 5.f - sin(sineC3.nextPhase());	// pass only the artifacts
		out0[f] = out1[f] = s * 0.2f;
		
//		out0[f] = sin(sineC3.nextPhase()) - sineT1.nextL();
	}

}

int main(int argc, char* argv[]){
	AudioIO io(256, 44100., audioCB, NULL, 2);
	Sync::master().spu(io.framesPerSecond());

	tbl::sine(table, 1<<l2Size);	// generate 1 sine cycle
	
	io.start();
	printf("\nPress 'enter' to quit...\n"); getchar();
	return 0;
}
