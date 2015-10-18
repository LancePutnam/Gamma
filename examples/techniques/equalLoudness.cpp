/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example:	Equal Loudness
Author:		Lance Putnam, 2015

Description:
This shows how to use the scl::eqLoudAmp function to synthesize a sine wave with
equal perceptual loudness at all frequencies. The effect is most noticeable at
lower frequencies. The function used is the reciprocal of an A-weighted equal
loudness curve based on the 40-phon level Fletcher-Munson curve.
*/

#include "../AudioApp.h"
#include "Gamma/scl.h"
#include "Gamma/Oscillator.h"
using namespace gam;

struct MyApp : public AudioApp{

	Sine<> osc;		// Sine oscillator
	LFO<> ramp;		// Used as frequency modulator
	bool doEqLoud;	// Whether to apply equal loudness

	MyApp(){
		ramp.period(12);
		doEqLoud = true;
	}
	
	void onAudio(AudioIOData& io){

		while(io()){

			// Set frequency of sine
			float frq = scl::pow3(ramp.triU()) * 10000 + 60;
			osc.freq(frq);

			if(ramp.cycled()){
				doEqLoud ^= true;
				printf("equal loudness %s\n", doEqLoud?"on":"off");
			}

			// Determine amplitude from frequency
			float amp = doEqLoud ? scl::eqLoudAmp(frq) : 1;

			// Generate sine and apply gain
			float s = osc() * amp * 0.2;

			io.out(0) = s;
			io.out(1) = s;
		}
	}
};

int main(){
	MyApp().start();
}
