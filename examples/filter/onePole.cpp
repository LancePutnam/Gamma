/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example:	One-pole Filter
Author:		Lance Putnam, 2012

Description:
This demonstrates the effect of a one-pole low-pass filter on a noise source.
One-pole filters are very effective at controlling the "tone" of sounds as they 
have an adjustable cutoff frequency and gentle slope.
*/

#include "../AudioApp.h"
#include "Gamma/Filter.h"
#include "Gamma/Noise.h"
#include "Gamma/Oscillator.h"
using namespace gam;

class MyApp : public AudioApp{
public:

	LFO<> mod;
	NoiseWhite<> src;
	OnePole<> onePole;

	MyApp(){
		mod.period(8);
		mod.phase(0.5);
	}

	void onAudio(AudioIOData& io){

		while(io()){

			float cutoff = scl::pow3(mod.triU()) * 10000;

			// Set one-pole cutoff frequency
			onePole.freq(cutoff);

			// Filter sample
			float s = onePole(src());

			io.out(0) = io.out(1) = s * 0.2f;
		}
	}
};

int main(){
	MyApp().start();
}
