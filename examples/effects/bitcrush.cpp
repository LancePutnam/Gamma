/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example:	Bitcrushing
Author:		Lance Putnam, 2012

Description:
This demonstrates a bitcrushing effect which adds a lo-fi, crunchy element to a
sound. It is a combination of sample rate and bit reduction (quantization) of
the input.
*/

#include "../AudioApp.h"
#include "Gamma/Effects.h"
#include "Gamma/Oscillator.h"
using namespace gam;

class MyApp : public AudioApp{
public:

	SineDs<> src;			// Modal strike
	Accum<> tmr;			// Timer for strike
	LFO<> modSR;			// Sample rate reduction modulator
	LFO<> modQn;			// Quantization modulator
	Quantizer<> qnt;		// The bitcrush effect

	MyApp(){
		src.resize(3);
		tmr.period(1);
		tmr.phaseMax();
		modSR.period(23);
		modQn.period(32);
	}

	void onAudio(AudioIOData& io){

		while(io()){

			if(tmr()){
				src.set(0,  220, 1, 2.0);
				src.set(1,  347, 1, 1.2);
				src.set(2, 1237, 1, 0.2);
				tmr.freq(rnd::uni(2., 1.));
			}

			// Produce modal strike
			float s = src();

			// Set sample rate quantization
			qnt.freq(modSR.triU()*4000 + 1400);

			// Set amplitude quantization
			qnt.step(modQn.paraU()*0.5);

			// Apply the bitcrush
			s = qnt(s);

			io.out(0) = io.out(1) = s * 0.2;
		}
	}
};

int main(){
	MyApp().start();
}

