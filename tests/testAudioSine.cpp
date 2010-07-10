#include <stdio.h>
#include "Gamma/AudioIO.h"
#include "Gamma/Oscillator.h"

using namespace gam;

// Several ways to produce a sine oscillator
Osc<>        sineT1(440);			// Table look-up (full-wave)
TableSine<>  sineT2(440);			// Table look-up (quarter-wave)
Quadra<>     sineQd(440);			// Complex mul
LFO<>        sineC1(440, -0.25f);	// Computed 3rd order poly
Sine<>       sineC2(440);			// Computed taylor
AccumPhase<> sineC3(440);			// Direct math.h sin() (ground truth)

//SineRes< Multi<2> > sineRes(Multi<2>(220, 222));

void audioCB(AudioIOData & io){
	float * out0 = io.out(0);
	float * out1 = io.out(1);
	unsigned long numFrames = io.framesPerBuffer();

	for(unsigned long f=0; f<numFrames; f++){
		float s =	sineT1() +
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

	sineT1.addSine(1);
	
	io.start();
	printf("\nPress 'enter' to quit...\n"); getchar();
	return 0;
}
