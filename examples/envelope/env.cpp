/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example:	Env
Author:		Lance Putnam, 2015

Description:
This demonstrates usage of the Env class, an n-stage exponential envelope.
In particular, we demonstrate how to construct an AHD (attack-hold-decay)
envelope using a 3-stage Env.

|  /----------\
| /            \
|/              \
+-------------------> time
  A     H      D

*/

#include "../AudioApp.h"
#include "Gamma/Envelope.h"
#include "Gamma/Oscillator.h"
using namespace gam;

class MyApp : public AudioApp{
public:

	Sine<> src;		// Source to envelope
	Env<3> env;		// 3-stage envelope

	MyApp(){

		// Set levels of 4 breakpoints
		env.levels(0, 0.2, 0.2, 0);

		// Set length, in seconds, of each stage
		env.lengths(2, 4, 2);

		// Set the curve factor of the segments
		env.curve(-4); // fast approach (default)
		//env.curve(0); // linear approach
		//env.curve(4); // slow approach

		// We can also specify to envelope to loop (default == false)
		//env.loop(true);
	}

	void onAudio(AudioIOData& io){

		while(io()){

			// Apply envelope to source
			float s = src() * env();

			io.out(0) = io.out(1) = s;
		}
	}
};

int main(){
	MyApp().start();
}
