/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example:	Voices: Basic Spawning
Author:		Lance Putnam, 2021

Description:
This demonstrates how to use the Voices class to implement a polyphonic synth.
First, create a subclass of Voice that implements a single voice of the synth.
Voice has attack and release methods that must be overridden. Second, plug the
custom Voice subclass into a Voices object to create the polyphonic synth.
*/

#include "../AudioApp.h"
#include "Gamma/rnd.h"
#include "Gamma/Envelope.h"
#include "Gamma/Oscillator.h"
#include "Gamma/Voices.h"
using namespace gam;

class MyApp : public AudioApp{
public:

	// Define a single voice in our polyphonic synth
	struct MyVoice : public Voice<float> {
		Sine<> osc;
		AD<> env{0.02, 6.}; // long decay to make polyphony evident

		// Define what happens when the voice attacks
		void onAttack(float freq, float amp){
			osc.freq(freq);
			env.reset();
			env.amp(amp);
		}

		// Define what happens when the voice is released
		// In our case, the voice is never explicitly released as it decays on
		// its own. However, it's good practice to define this anyway.
		void onRelease(){
			env.release();
		}

		// Generates the next sample of output
		// This MUST be overridden.
		float operator()(){
			return osc() * env();
		}
	};


	Voices<MyVoice, 16> voices; // A 16-voice polyphonic synth
	Accum<> tmr{1., 0.999};		// Timer for triggering new notes


	void onAudio(AudioIOData& io){
		while(io()){

			// Trigger a new voice on a timer
			if(tmr()){
				// Not a big fan of atonal music, but it's easier to generate...
				float f = pow(2, rnd::uni(36)/12.)*110;
				// This spawns a new voice calling onAttack
				voices.attack(f, 0.2);
			}

			// Get synth output (sum of all active voices)
			float s = voices();

			io.out(0) = s;
			io.out(1) = s;
		}
	}
};

int main(){
	MyApp().start();
}

