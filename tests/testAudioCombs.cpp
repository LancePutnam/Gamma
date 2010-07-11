/*
	Resonating feedback combs flanged by a feedforward comb.
*/

#include <stdio.h>
#include "Gamma/Gamma.h"
#include "Gamma/AudioIO.h"
#include "Gamma/Delay.h"
#include "Gamma/Oscillator.h"

using namespace gam;

Comb<float, ipl::Cubic> cmbS1(1./220., 0, -0.97), cmbS2(1./330., 0, -0.97),
		cmbF(1./100., 1, 0.2), cmbE(0.28,0,0.9);
Biquad<> biq(14000, 0.2, Filter::LP);
LFO<> lfo(1./12.);

void audioCB(AudioIOData & io){
	float * out0 = io.out(0);
	float * out1 = io.out(1);
	unsigned long numFrames = io.framesPerBuffer();

	// Fill buffer with white noise
	rnd::uniS(out0, numFrames, rnd::prob(0.01) ? 0.15f : 0.f);

	for(unsigned long f=0; f<numFrames; f++){
		float smp = out0[f];
		
		// Harmonically resonate noise
		smp = (cmbS1(smp) + cmbS2(smp)) * 0.5f;
		
		// Flanger
		cmbF.delay(lfo.cosU() * 0.008f + 0.0001f);
		smp = cmbF(smp);
		
		// Low-pass echo
		smp += cmbE(smp, biq(cmbE())) * 0.5;
		
		out0[f] = out1[f] = smp;
	}
}

int main(int argc, char* argv[]){
	AudioIO io(256, 44100., audioCB, NULL, 2);
	Sync::master().spu(io.framesPerSecond());
	io.start();
	printf("\nPress 'enter' to quit...\n"); getchar();
	return 0;
}

