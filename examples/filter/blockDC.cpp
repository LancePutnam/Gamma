/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example:	DC Blocker
Author:		Lance Putnam, 2012

Description:
This shows how a DC blocker can be used to eliminate unwanted zero/low 
frequencies.
*/

#include "../AudioApp.h"
#include "Gamma/Filter.h"
#include "Gamma/Oscillator.h"
using namespace gam;

class MyApp : public AudioApp{
public:

	Sine<> osc;			// The signal we want to hear
	Sine<> oscLow;		// Unwanted low-frequency component
	BlockDC<> blockDC;
	Accum<> tmr;
	bool blockingOn;

	MyApp(){
		osc.freq(440);
		oscLow.freq(0.1);
		tmr.period(4);
		tmr.phaseMax();
		blockingOn = true;
	}

	void onAudio(AudioIOData& io){

		while(io()){
			if(tmr()){
				blockingOn ^= true;
				printf("DC blocking %s\n", blockingOn?"on":"off");
			}
			
			// Generate a sine wave plus a very low frequency oscillation
			float s = osc() + oscLow();
			
			// Do DC blocking?
			if(blockingOn) s = blockDC(s);
			
			// Simulate clipping on the DAC
			s = scl::clipS(s) * 0.2;
				
			io.out(0) = io.out(1) = s;
		}
	}
};

int main(){
	MyApp().start();
}
