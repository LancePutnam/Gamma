/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example:	Frequency Modulation Transition
Author:		Lance Putnam, 2012

Description:
This demonstrates the transition from sub-audio- to audio-rate frequency
modulation. At sub-audio rate, we hear a vibrato effect. At audio rate, we hear
sidebands described by Bessel functions of the first kind.
*/
#include "../AudioApp.h"
#include "Gamma/Oscillator.h"
using namespace gam;

class MyApp : public AudioApp{
public:

	Sine<> car;	// Carrier sine (gets its frequency modulated)
	Sine<> mod;	// Modulator sine (used to modulate frequency)
	LFO<> env{1./30};

	void onAudio(AudioIOData& io){

		float fc = 220;		// Carrier frequency, e.g., the fundamental
		float I = 3;		// Modulation index

		while(io()){
			float ratio = env.hann();

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

