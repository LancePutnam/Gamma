/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example: Chirplet grain
Author: Lance Putnam, 2015

Description:
Chirplet is a Gaussian enveloped sine wave produced from two complex multiplies.
It is a "set-and-forget" unit generator in that it is very cheap to generate
samples, but expensive to set its parameters. For the most part, this means we
do not want to modulate its frequency.
*/

#include "../AudioApp.h"
#include "Gamma/rnd.h"
#include "Gamma/Oscillator.h"
using namespace gam;

class MyApp : public AudioApp{
public:

	Chirplet<> osc;	// Gaussian enveloped sine

	void onAudio(AudioIOData& io){

		while(io()){

			// Retrigger chirplet when done
			if(osc.done()){
				// Set length of envelope
				float length = rnd::prob(0.8) ? osc.length() : rnd::lin(0.5, 0.01);
				osc.length(length);

				// Set start and end frequencies
				osc.freq(rnd::lin(2000, 40), rnd::lin(2000, 40));
				//osc.amp(0.2);
			}
			
			float s = osc().i*0.2;

			io.out(0) = io.out(1) = s;
		}
	}

};

int main(){
	MyApp().start();
}
