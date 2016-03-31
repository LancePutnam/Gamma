/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example:	Panning
Author:		Lance Putnam, 2012

Description:
A simple two-channel panning effect.
*/

#include "../AudioApp.h"
#include "Gamma/Effects.h"
#include "Gamma/Oscillator.h"
using namespace gam;

class MyApp : public AudioApp{
public:

	Sine<> src;		// Mono source
	Pan<> pan;		// Two-channel panner
	LFO<> panMod;	// Pan position modulator

	void onAudio(AudioIOData& io){

		src.freq(220);
		panMod.period(2);

		while(io()){

			// Generate our mono signal
			float s = src() * 0.2;

			// Modulate pan position between left (-1) and right (1) channels		
			pan.pos(panMod.tri());

			// The output is two floats (stereo)
			float2 s2 = pan(s);

			io.out(0) = s2[0];
			io.out(1) = s2[1];

		}
	}

};

int main(){
	MyApp().start();
}
