/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information

Example:	"Buzz" Oscillator
Author:		Lance Putnam, 2015

Description:
Buzz generates a finite sum of equal amplitude harmonic cosine waves.
Note it is not a good idea to change the number of harmonics dynamically as it
will generally produce a click. It is only done here for demonstration purposes.
*/

#include "../AudioApp.h"
#include "Gamma/Oscillator.h"
using namespace gam;

class MyApp : public AudioApp{
public:

	Buzz<> buzz;
	LFO<> mod;

	void onAudio(AudioIOData& io){

		// Set the frequency, in Hz
		buzz.freq(55);

		mod.period(16);

		while(io()){

			int nh = mod.triU() * 64;

			// Set number of harmonics
			buzz.harmonics(nh);
	
			float s = buzz() * 0.2;

			io.out(0) = s;
			io.out(1) = s;
		}
	}

};

int main(){
	MyApp().start();
}
