/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example:	Attack-decay (AD) Envelope
Author:		Lance Putnam, 2012

Description:
This shows how to use an exponential AD (attack-decay) envelope to control the 
amplitude of a sound source.
*/

#include "../AudioApp.h"
#include "Gamma/Envelope.h"
#include "Gamma/Oscillator.h"
using namespace gam;

class MyApp : public AudioApp{
public:

	Accum<> tmr;		// Timer for resetting envelope
	Sine<> src;			// Source to envelope
	AD<> env;			// Attack-decay envelope
	float tilt;			// Tilt of envelope; 0=percussive, 1=reversive

	MyApp(){
		tmr.period(1.2);
		tmr.phaseMax();
		tilt = 0;
	}

	void onAudio(AudioIOData& io){

		while(io()){

			if(tmr()){
				env.attack(tilt);	// Set attack, in seconds
				env.decay(1-tilt);	// Set decay, in seconds
				env.amp(0.2);		// Set amplitude
				env.reset();		// Reset envelope

				// Increment tilt amount
				tilt += 0.1;
				if(tilt > 1) tilt=0;
			}

			float s = src() * env();

			io.out(0) = io.out(1) = s;
		}
	}
};

int main(){
	MyApp().start();
}
