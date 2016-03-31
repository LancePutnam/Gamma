/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example:	Echo
Author:		Lance Putnam, 2012

Description:
This demonstrates how to create an echo effect from a delay line. The output of
the delay line is fed back into the input to create a series of exponentially
decaying echoes. A low-pass "loop" filter is inserted into the feedback path to
simulate high-frequency damping due to air absorption.
*/
#include "../AudioApp.h"
#include "Gamma/Delay.h"
#include "Gamma/Filter.h"
#include "Gamma/Oscillator.h"
using namespace gam;

class MyApp : public AudioApp{
public:

	Accum<> tmr;	// Timer for triggering sound
	SineD<> src;	// Sine grain
	Delay<> delay;	// Delay line
	OnePole<> lpf;

	MyApp(){
		// Allocate 200 ms in the delay line
		delay.maxDelay(0.2);

		tmr.period(4);
		tmr.phaseMax();

		// Configure a short cosine grain
		src.set(1000, 0.8, 0.04, 0.25);

		// Set up low-pass filter
		lpf.type(gam::LOW_PASS);
		lpf.freq(2000);
	}

	void onAudio(AudioIOData& io){
		while(io()){

			if(tmr()) src.reset();

			float s = src();

			// Read the end of the delay line to get the echo
			float echo = delay();

			// Low-pass filter and attenuate the echo
			echo = lpf(echo) * 0.8;

			// Write sum of current sample and echo to delay line
			delay(s + echo);

			// Finally output sum of dry and wet signals
			s += echo;
		
			io.out(0) = io.out(1) = s;
		}
	}
};

int main(){
	MyApp().start();
}
