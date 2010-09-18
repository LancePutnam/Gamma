/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Tutorial:		Generator / Envelope
	Description:	This demonstrates how to use a Seg to generate smooth
					noise. The result is used to control the pitch of an
					oscillator.
*/

#include "tutorial.h"

NoiseWhite<> noise;
Seg<float, iplSeq::Linear> seg(1/20.);
Sine<> osc;

void audioCB(AudioIOData& io){

	for(int i=0; i<io.framesPerBuffer(); ++i){
	
		float s = seg(noise);
		osc.freq(1000 + 400*s);
		s = osc()*0.1;

		io.out(0)[i] = io.out(1)[i] = s;
	}
}

RUN(audioCB);
