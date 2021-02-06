/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example:	Amplitude Modulation Transition
Author:		Lance Putnam, 2012

Description:
This demonstrates the transition from sub-audio- to audio-rate amplitude 
modulation. At sub-audio rate, we hear a tremolo effect. At audio rate, we hear
sidebands---the sum and difference tones.
*/
#include "../AudioApp.h"
#include "Gamma/Oscillator.h"
using namespace gam;

class MyApp : public AudioApp{
public:

	Sine<> oscC{440};		// Carrier
	Sine<> oscM;			// Modulator
	LFO<> modFreq{1./20};	// Envelope for changing modulator frequency

	void onAudio(AudioIOData& io){

		while(io()){

			// Set modulator frequency (time-varying)
			oscM.freq(modFreq.hann() * 110 + 1);

			// Do the amplitude modulation
			float s = oscC() * (oscM() + 1.) * 0.2;

			io.out(0) = io.out(1) = s;
		}
	}
};

int main(){
	MyApp().start();
}
