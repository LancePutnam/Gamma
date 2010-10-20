/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Example:		Generator / Envelope
	Description:	This demonstrates how to use a Seg to generate smooth
					noise. The result is used to control the pitch of an
					oscillator.
*/

#include "../examples.h"

NoiseWhite<> noise;
Seg<float, iplSeq::Linear> seg(1/20.);
Sine<> osc;

void audioCB(AudioIOData& io){

	while(io()){
	
		float s = seg(noise);
		osc.freq(1000 + 400*s);
		s = osc()*0.1;

		io.out(0) = io.out(1) = s;
	}
}

RUN(audioCB);
