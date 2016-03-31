/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example:	Notes From A Sinusoid
Author:		Lance Putnam, 2012

Description:
This is an example of using a high-frequency sine oscillator to determine the 
pitches of notes. This is technique is inspired by the Strang moiré patterns 
described in:

	Richet, N. (1992). Strang’s strange figures.
	The American Mathematical Monthly, 99(2):101-107.
*/

#include "../AudioApp.h"
#include "Gamma/rnd.h"
#include "Gamma/Envelope.h"
#include "Gamma/Oscillator.h"
using namespace gam;

class MyApp : public AudioApp{
public:

	Sine<> src;
	AD<> env;
	Accum<> tmr;
	double phase;

	MyApp(){
		tmr.period(0.15);
		env.attack(0.01).decay(0.14);
		phase = 0;
	}

	void onAudio(AudioIOData& io){
		while(io()){
			if(tmr()){
				// Add to the phase a low-integer division of 2π plus a small 
				// offset for variation.
				phase = scl::wrapPhase(phase + 2*M_PI / (rnd::prob(0.9)?3:2) + 0.012);

				float frq = sin(phase) * 8;
				
				// Map sine value to nearest pitch in C-major scale
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
