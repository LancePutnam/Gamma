/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example:	Envelope Segment
Author:		Lance Putnam, 2012

Description:
This uses an envelope segment to make an oscillator randomly glide between
frequencies.
*/

#include "../AudioApp.h"
#include "Gamma/rnd.h"
#include "Gamma/Envelope.h"
#include "Gamma/Oscillator.h"
using namespace gam;

class MyApp : public AudioApp{
public:

	Accum<> tmr;	// Timer for resetting envelope
	Sine<> osc;		// Source oscillator
	Seg<> env;		// Envelope segment

	MyApp(){
		tmr.period(2);		// Set new target frequency every 2 seconds
		tmr.phaseMax();
		env.length(0.5);	// Set length, in seconds, of the segment
	}

	void onAudio(AudioIOData& io){

		while(io()){
			if(tmr()){
				// Set new target value of envelope
				env = 100. + rnd::uni(400.);
			}

			// Map envelope to oscillator frequency
			osc.freq(env());

			float s = osc() * 0.2;

			io.out(0) = io.out(1) = s;
		}
	}
};

int main(){
	MyApp().start();
}
