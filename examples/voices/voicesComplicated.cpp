/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example:	Voices: Doing More Complicated Stuff...
Author:		Lance Putnam, 2021

Description:
This demonstrates more advanced use of Voices, namely, how to tie it to a
keyboard controller for indefinitely sustained notes and modulation effects.
*/

#include "../AudioApp.h"
#include "Gamma/rnd.h"
#include "Gamma/Envelope.h"
#include "Gamma/Oscillator.h"
#include "Gamma/Voices.h"
#include <array>
using namespace gam;

class MyApp : public AudioApp{
public:

	// Our voice is in stereo, so we pass float2 as the Voice sample type
	struct MyVoice : public Voice<gam::float2> {
		Sine<> osc;
		ADSR<> env{0.02, 0.1, 0.7, 3.};
		float pan = 0.5;

		// Define what happens when the voice attacks
		void onAttack(float freq, float amp){
			osc.freq(freq);
			env.reset();
			env.amp(amp);
		}

		// Define what happens when the voice is released
		void onRelease(){
			env.release();
		}

		// Override this to give the parent Voices a better idea of when the
		// voice is done and can be released. If not overridden, a default
		// silence detection algorithm is used.
		bool done(){
			return env.done();
		}

		// value_type is the template parameter passed into Voice
		value_type operator()(){
			float s = osc() * env();
			return {s*(1.f-pan), s*pan};
		}
	};


	Voices<MyVoice, 16> voices; 	// A 16-voice polyphonic synth
	Sine<> vibOsc{5};				// Global vibrato oscillator
	AD<> vibEnv{1, 4, 0.01};		// Global vibrato envelope
	Accum<> tmr{1., 0.999};			// Timer for triggering new notes
	std::array<bool,4> keyStates; 	// State of keyboard keys

	MyApp(){

		// Iterate through all voices (for additional setup, etc.)
		for(auto& v : voices){
		}

		for(auto& v : keyStates) v=false;
	}


	void onAudio(AudioIOData& io){
		while(io()){

			if(tmr()){
				// Not a big fan of atonal music, but it's easier to generate...
				float f = pow(2, rnd::uni(36)/12.)*110;

				// Get a random keybord key
				int i = rnd::uni(keyStates.size());

				if(keyStates[i] ^= true){ // press key
					// Attack with the ID of a keyboard key / MIDI note
					// All attack functions return a reference to the new voice
					// which we can use for further configuration.
					auto& voice = voices.attackWithID(i, f, 0.2);

					// Delay onset of voice (in samples)
					voice.delay(io.fps() * rnd::uni(0.3));

					voice.pan = rnd::uni(1.);

				} else { // release key
					// Release the voice with the passed-in key ID
					voices.releaseWithID(i);
				}

				// Trigger a random vibrato modulation
				if(rnd::prob(0.1)) vibEnv.reset();

				// We can release all voices (in case of stuck notes)
				if(rnd::prob(0.01)) voices.releaseAll();

				voices.print(); // print out debugging info
			}

			// Apply a vibrato to all active voices
			float vib = vibOsc() * vibEnv();
			voices.traverseActive([&](unsigned i){
				voices.voice(i).osc.freqMul(1. + vib);
			});

			auto s = voices();

			io.out(0) = s[0];
			io.out(1) = s[1];
		}
	}
};

int main(){
	MyApp().start();
}

