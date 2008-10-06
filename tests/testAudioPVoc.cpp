#include <stdio.h>
#include "AudioIO.h"
#include "Delay.h"
#include "DFT.h"
#include "Noise.h"
#include "Oscillator.h"
#include "Gamma.h"
#include "MacroD.h"

using namespace gam;

NoisePink<> src;
//Sine<float> src(440);
LFO<> lfoA1(1./1000., 0, 0.2), lfoA2(1./800., 0, 0.2), lfoF1(1./1.3), lfoF2(1./1.61);
Biquad<> bq0(1./400., 10, Filter::BP);
Biquad<> bq1(1./400., 10, Filter::BP);

ULONG hopSize = 512;
STFT stft(hopSize * 4, hopSize, 0, WinType::Hann, Bin::Polar);

void audioCB(AudioIOData & io){
	float * out0 = io.out(0);
	float * out1 = io.out(1);
	ULONG numFrames = io.numFrames();

	for(ULONG f=0; f<numFrames; ++f){
		float smp = src() * 0.5f;
		
		if( stft(smp) ){

			lfoA1.phase(lfoF1.upU());
			lfoA2.phase(lfoF2.downU());
			
			for(int i=0; i<stft.numBins(); ++i){
			
				stft.bins(i)[0] *= lfoA1.triU() - lfoA2.triU();
				//f0 *= lfoA1.impulse();
				//f0 = bq0(s0);
				//f1 = bq1(s1);
			}	
		}
		out0[f] = out1[f] = stft();
	}
}

int main(int argc, char* argv[]){
	AudioIO io(128, 44100., audioCB, NULL, 2);
	Sync::master().spu(io.fps());

	stft.syncHop()  << lfoF1 << lfoF2;
	stft.syncFreq() << lfoA1 << lfoA2 << bq0 << bq1;
		
	io.start();
	printf("\nPress 'enter' to quit...\n"); getchar();
	return 0;
}
