/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example:	Notes From Pink Noise
Author:		Lance Putnam, 2012

Description:
This is an example of using a pink noise generator to determine the pitches of 
notes. Compared to other noise types, like white and brown, pink noise 
supposedly produces the most musical results.
*/
#include "../AudioApp.h"
#include "Gamma/Envelope.h"
#include "Gamma/Noise.h"
#include "Gamma/Oscillator.h"
using namespace gam;

class MyApp : public AudioApp{
public:

	Sine<> src;
	AD<> env;
	NoisePink<> pink;
	Accum<> tmr;

	MyApp(){
		tmr.period(0.25);
		env.attack(0.01).decay(0.24);
	}

	void onAudio(AudioIOData& io){
		while(io()){
			if(tmr()){
				float frq = pink() * 12;
				
				// Map noise value to nearest pitch in C-major scale
				frq = scl::ratioET(scl::nearest(frq, "2212221")) * scl::freq("c4");

				src.freq(frq);
				env.reset();
			}

			float s = src() * env() * 0.2;

			io.out(0) = s;
			io.out(1) = s;
		}
	}
};

int main(){
	MyApp().start();
}
