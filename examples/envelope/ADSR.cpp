/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example:	ADSR Envelope
Author:		Lance Putnam, 2015

Description:
This shows how to use an ADSR (attack-decay-sustain-release) envelope to create
a sound event. The envelope goes from 0 to 1 over the attack time, from 1 to the
sustain level over the decay time, and the sustain level to 0 over the release
time. The envelope will hold its value indefinitely at the sustain level until
told to release. A typical ADSR shape looks something like the following:

|   /\
|  /  \-------\ ..... S
| /            \
|/              \
+-------------------> time
  A  D         R
*/

#include "../AudioApp.h"
#include "Gamma/Envelope.h"
#include "Gamma/Oscillator.h"
using namespace gam;

class MyApp : public AudioApp{
public:

	Sine<> src;		// Source to envelope
	ADSR<> env;		// ADSR envelope
	float time;

	MyApp(){
		time = 0;
		env.attack(0.01);	// Set attack time in seconds
		env.decay(1);		// Set decay time in seconds
		env.sustain(0.5);	// Set sustain level
		env.release(5);		// Set release time in seconds
		env.amp(0.2);		// Set overall amplitude (default == 1)
	}

	void onAudio(AudioIOData& io){

		// Count seconds passed
		time += io.secondsPerBuffer();

		while(io()){

			// Enter release stage after 4 seconds
			if(time > 4){
				env.release();
			}

			// Apply envelope to source
			float s = src() * env();

			io.out(0) = io.out(1) = s;
		}
	}
};

int main(){
	MyApp().start();
}
