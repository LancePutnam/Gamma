/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example: Damped Sine Wave
Author: Lance Putnam, 2015

Description:
The uses SineD to generate an exponentially decaying (damped) sinusoid.
SineD is a "set-and-forget" unit generator in that it is very cheap to generate
samples, but expensive to set its parameters. For the most part, this means we
do not want to modulate its frequency.
*/

#include "../AudioApp.h"
#include "Gamma/rnd.h"
#include "Gamma/Oscillator.h"
using namespace gam;

class MyApp : public AudioApp{
public:

	Accum<> tmr;	// Timer for randomizing parameters
	SineD<> osc;	// Damped sinusoid

	MyApp(){
		tmr.finish();
	}

	void onAudio(AudioIOData& io){

		// Set period of timer, in seconds
		tmr.period(1);

		while(io()){

			// Retrigger damped sine periodically
			if(tmr()){
				osc.set(
					rnd::uni(10, 1)*50,	// frequency, in Hz
					0.2,				// amplitude
					rnd::lin(2., 0.1)	// decay length, in seconds
				);
			}
			
			float s = osc();

			io.out(0) = io.out(1) = s;
		}
	}

};

int main(){
	MyApp().start();
}
