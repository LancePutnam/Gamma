/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Example:		Effect / Bitcrush
	Description:	This demonstrates a bitcrushing effect which adds a lo-fi,
					crunchy element to a sound. It is a combination of sample
					rate and bit reduction (quantization) of the input.
*/

#include "../examples.h"

SineDs<> src(3);			// Modal strike
Accum<> tmr(1,1);
LFO<> modSR(1./ 8, 0.0);	// Sample rate reduction modulator
LFO<> modQn(1./32, 0.5);	// Quantization modulator
Quantizer<> qnt;			// The bitcrush effect

void audioCB(AudioIOData& io){

	while(io()){

		if(tmr()){
			src.set(0,  220, 1, 2.0);
			src.set(1,  347, 1, 1.2);
			src.set(2, 1237, 1, 0.2);
			tmr.freq(rnd::uni(2., 1.));
		}

		qnt.freq(modSR.triU()*4000 + 1400);
		qnt.step(modQn.paraU()*0.25);

		float s = src() * 0.2;
		s = qnt(s);

		io.out(0) = io.out(1) = s;
	}
}

RUN(audioCB);
