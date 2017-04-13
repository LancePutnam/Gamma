/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example:	Oscillator drifting
Author:		Lance Putnam, 2012

Description:
This demonstrates how to add oscillator drift to produce a more human-like
sound or emulate the slight imprecisions of analog generators.
*/

#include "../AudioApp.h"
#include "Gamma/Noise.h"
#include "Gamma/Filter.h"
#include "Gamma/Oscillator.h"
using namespace gam;

class MyApp : public AudioApp{
public:

	Accum<> tmr{1./8};
	LFO<> osc;							// Oscillator
	Upsample<NoiseWhite<>> drift{1./2};	// Low-frequency noise
	bool applyDrift = false;

	void onAudio(AudioIOData& io){

		while(io()){

			if(tmr()) applyDrift^=true;
			if(tmr.cycled()) printf("Drift %s\n", applyDrift?"on":"off");

			float freq = 440;
			float amp = 1;
			if(applyDrift){
				float d = drift(); // in [-1,1]
				freq *= 1 + d*0.02;
				amp *= 1 + d*0.2;
			}

			osc.freq(freq);
			
			float s = osc.tri();
				
			io.out(0) = io.out(1) = s*0.2;
		}
	}
};

int main(){
	MyApp().start();
}
