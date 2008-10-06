/*
Audio example of using a multiple table oscillator to produce bandlimited
waveforms.
*/

#include <stdio.h>
#include "AudioIO.h"
#include "Oscillator.h"

using namespace gam;

// Waveform LUTs
const uint32_t tableBits = 8;
const uint32_t numTables = 8;
const uint32_t numSamples = (1<<tableBits) * numTables;
float tables[numSamples];

MultiTableOsc<> osc(tables, tableBits, numTables, 110);
LFO<> lfoF(0.5);
Accum<> tmr(lfoF.freq(), 2);
gen::Counter cnt(3);

void audioCB(AudioIOData & io){
	float * out0 = io.out(0);
	float * out1 = io.out(1);

	for(unsigned long i=0; i<io.numFrames(); ++i){
		
		if(tmr()){
			mem::zero(tables, numSamples);
			switch(cnt.val){
				case 0:  tbl::multiWave(tables, 1<<tableBits, numTables, tbl::triangleSum); break;
				case 1:  tbl::multiWave(tables, 1<<tableBits, numTables, tbl::squareSum); break;
				case 2:  tbl::multiWave(tables, 1<<tableBits, numTables, tbl::sawSum); break;
				default:;
			}
			cnt();
		}
	
		osc.freqLL(scl::pow2( lfoF.upU() ) * 4000.f);
		out0[i] = out1[i] = osc.nextLL() * 0.1f;
	}
}

int main(int argc, char* argv[]){
	AudioIO io(64, 44100., audioCB, NULL, 2);
	Sync::master().spu(io.fps());	
	io.start();
	printf("\nPress 'enter' to quit...\n"); getchar();
	return 0;
}
