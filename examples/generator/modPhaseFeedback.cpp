/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Example:		Generator / Modulation / Phase
	Description:	Single-oscillator feedback phase modulation.
*/

#include "../examples.h"

float prev=0;
Sine<> osc(220);			// source sine
LFO<> mod(1./10, 0.5);		// modulator of feedback amount

void audioCB(AudioIOData& io){

	while(io()){

		float fbk = prev*mod.hann()*0.2;

		osc.phaseAdd(fbk);

		prev = osc();
		io.out(0) = io.out(1) = prev * 0.2f;
		
		osc.phaseAdd(-fbk);
	}
}

RUN(audioCB);
