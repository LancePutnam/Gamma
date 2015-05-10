/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example:	Delay-line
Author:		Lance Putnam, 2012

Description:
This demonstrates the multiple uses of a delay-line.
*/

#include "../AudioApp.h"
#include "Gamma/Delay.h"
#include "Gamma/Effects.h"
using namespace gam;

class MyApp : public AudioApp{
public:

	Delay<float, ipl::Trunc> delay;	// Delay-line with truncating interpolation
	//Delay<float, ipl::Linear> delay;	// Delay-line with linear interpolation
	//Delay<float, ipl::Cubic> delay;	// Delay-line with cubic interpolation
	//Delay<float, ipl::AllPass> delay;	// Delay-line with allpass interpolation
	Accum<> tmr;
	Burst burst;

	MyApp()
	:	burst(2e4,2e3, 0.1)
	{
		delay.maxDelay(0.4);	// Allocate delay-line memory for 400 ms
		delay.delay(0.2);		// Set delay time to 200 ms
		tmr.period(1);
		tmr.phaseMax();
	}

	void onAudio(AudioIOData& io){

		while(io()){
			if(tmr()) burst.reset();

			//float s = io.in(0);
			float s = burst()*0.5;

			// This short-hand method is convenient for simple delays.
			s += delay(s);

			// The delay-line can also be read from and written to using
			// separate calls.
			//delay.write(s);
			//s += delay();

			// We can create infinite echoes by feeding back a small amount of
			// the output back into the input on each iteration.
			//s += delay(s + delay()*0.2);

			// We can also create multitap delay-lines through multiple calls to
			// the read() method.
			//s += delay(s) + delay.read(0.15) + delay.read(0.39);

			// How about multitap feedback?
			//s += delay(s + delay.read(0.197)*0.3 + delay.read(0.141)*0.4 + delay.read(0.093)*0.2)*0.5;

			io.out(0) = io.out(1) = s;
		}
	}
};

int main(){
	MyApp().start();
}
