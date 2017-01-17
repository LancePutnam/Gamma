/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example:	Generator / Envelope
Author:		Lance Putnam, 2012

Description:
This demonstrates how to use an exponential segment to smooth out sudden changes
in frequency of an oscillator which, musically, results in a portamento.
*/

#include "../AudioApp.h"
#include "Gamma/Envelope.h"
#include "Gamma/Oscillator.h"
using namespace gam;

class MyApp : public AudioApp{
public:

	Accum<> tmr{1};		// Timer for selecting new frequency
	Sine<> osc{100};	// Test oscillator
	SegExp<> env{		// Exponential envelope used to smooth pitch
		0.2,		// Length in seconds
		-4,			// Curvature
		2000,100	// Start/end values
	};

	void onAudio(AudioIOData& io){
		while(io()){
		
			if(tmr()){
				// Change the frequency
				float freq = osc.freq();
				if(freq < 900){
					env = freq + 100;
				} else {
					env = 100;
				}
			}

			// Assign envelope value to oscillator frequency
			osc.freq(env());

			float s = osc() * 0.1;

			io.out(0) = io.out(1) = s;
		}
	}
};

int main(){
	MyApp().start();
}
