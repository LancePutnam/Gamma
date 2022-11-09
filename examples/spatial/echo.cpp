/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example:	Echo
Author:		Lance Putnam, 2022

Description:
This demonstrates how to create an echo effect. The echo is implemented using
a feedback delay. A custom filter may be inserted in the feedback loop. The
default loop filter is a simple gain. A low-pass loop filter (Loop1P) can
emulate air absorption effects.
*/

#include "../AudioApp.h"
#include "Gamma/Oscillator.h"
#include "Gamma/Spatial.h"
using namespace gam;

class MyApp : public AudioApp{
public:

	SineD<> src;		// Decaying sine wave grain
	Accum<> tmr;		// Timer for firing grains
	unsigned seed = 1;	// RNG seed
	//Echo<> echo{0.5};	// Echo w/ simple gain loop filter: arg is (max) delay time
	Echo<float, ipl::Linear, Loop1P> echo{0.5};	// Echo w/ low-pass loop filter

	MyApp(){
		tmr.period(2);
		tmr.finish();

		// Set delay time
		echo.delay(0.11);

		// Set decay length, in seconds
		echo.decay(20);

		// Set high-frequency damping factor in [0, 1] (loop filter only)
		echo.damping(0.7);

		// Alternatively, negative values create a high-pass effect
		//echo.damping(-0.2);
	}

	void onAudio(AudioIOData& io){
		while(io()){

			// Trigger a new random grain on a timer
			if(tmr()){
				seed *= 69069;
				float f = (seed%2000) + 400;
				src.set(f, 1, 32./f);				
			}

			// Get next grain sample
			float s = src();

			// Pass grain through echo and mix with dry source
			s += echo(s) * 0.5;

			io.out(0) = s;
			io.out(1) = s;
		}
	}
};

int main(){
	MyApp().start();
}

