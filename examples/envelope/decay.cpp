/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example:	Exponential Decay
Author:		Lance Putnam, 2012

Description:
This shows how to use an exponentially decaying envelope to control the
amplitude of a sound source. Since this envelope never reaches zero, its length
is the amount of time it takes to decay by -60 dB or 1/1000th of its starting
value.
*/

#include "../AudioApp.h"
#include "Gamma/Envelope.h"
#include "Gamma/Oscillator.h"
using namespace gam;

class MyApp : public AudioApp{
public:

	Accum<> tmr;	// Timer for resetting envelope
	Sine<> src;		// Source to envelope
	Decay<> env;	// Exponentially decaying envelope

	MyApp(){
		tmr.period(1);
		tmr.phaseMax();
	}

	void onAudio(AudioIOData& io){

		while(io()){
			if(tmr()){
				env.decay(1);	// Set decay length to 1 second
				env.reset(0.2);	// Reset envelope and specify amplitude
			}

			// Envelope source
			float s = src() * env();

			io.out(0) = io.out(1) = s;
		}
	}
};

int main(){
	MyApp().start();
}
