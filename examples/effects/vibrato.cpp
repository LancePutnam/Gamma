/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Example:		Effect / Vibrato
	Description:	This demonstrates how to create a vibrato effect by
					slowly modulating the delay time of a delay line.
*/

#include "../examples.h"


struct Vibrato{

	Vibrato(float modAmount=1./400, float modFreq=5)
	:	modAmount(modAmount),
		delay(0.1, 0), mod(modFreq)
	{}
	
	float operator()(float i){
		delay.delay(mod.hann()*modAmount + 0.0001);
		return delay(i);
	}

	float modAmount;
	Delay<> delay;
	LFO<> mod;
};


LFO<> src(220);			// A rich source
Vibrato vibrato;		// Vibrato unit


void audioCB(AudioIOData& io){

	while(io()){

		float s = src.tri();
		s = vibrato(s);

		io.out(0) = io.out(1) = s;
	}
}

RUN(audioCB);
