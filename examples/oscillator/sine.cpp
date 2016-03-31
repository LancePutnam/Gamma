/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example:	Sine Oscillator
Author:		Lance Putnam, 2015

Description:
This demonstrates how to use the Sine unit generator to create a sine wave at
440 Hz.
*/

#include "../AudioApp.h"
#include "Gamma/Oscillator.h"
using namespace gam;

class MyApp : public AudioApp{
public:

	// Declare a single sine wave oscillator
	Sine<> osc;

	void onAudio(AudioIOData& io){

		// Set the frequency of the sine oscillator, in Hz
		osc.freq(440);

		while(io()){

			// Generate next sine wave sample
			float s = osc();

			// Scale the sample down a bit for output
			s *= 0.2;

			io.out(0) = s;
			io.out(1) = s;
		}
	}

};

int main(){
	MyApp().start();
}
