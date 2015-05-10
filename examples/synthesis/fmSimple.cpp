/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example:	Simple Frequency Modulation
Author:		Lance Putnam, 2012

Description:
This demonstrates simple FM using two sine oscillators---a carrier and a 
modulator. The c:m ratio controls the multiplier between sinusoidal partials
---integer ratios produce harmonic sounds, non-integer ratios produce inharmonic
sounds. The index of modulation, I, determines roughly the number the sinusoidal
side-bands around the carrier frequency.
*/
#include "../AudioApp.h"
#include "Gamma/Oscillator.h"
using namespace gam;

class MyApp : public AudioApp{
public:

	Sine<> car;	// Carrier sine (gets its frequency modulated)
	Sine<> mod;	// Modulator sine (used to modulate frequency)

	void onAudio(AudioIOData& io){

		float fc = 220;		// Carrier frequency, e.g., the fundamental
		float ratio = 2./1;	// c:m ratio (described here as fm/fc)
		float I = 3;		// Modulation index

		while(io()){

			// Modulation frequency is carrier frequency over the c:m ratio
			float fm = fc*ratio;
			mod.freq(fm);

			// Compute frequency deviation signal
			float df = mod() * fm * I;
			
			// Set frequency of carrier
			car.freq(fc + df);

			// Most of FM is just mapping to the carrier frequency!

			// Here, we generate the carrier's next sample.
			float s = car() * 0.2;

			io.out(0) = s;
			io.out(1) = s;
		}
	}
};

int main(){
	MyApp().start();
}

