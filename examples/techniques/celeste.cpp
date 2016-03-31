/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example:	Celeste Oscillator
Author:		Lance Putnam, 2012

Description:
The demonstrates the old organ technique of voix celeste ("heavenly voice")
whereby two pipes are sounded together slightly off pitch. This produces a very
rich beating pattern with small detune amounts. Typically, the amount of 
detuning will track the base frequency so that higher frequencies beat faster
and lower frequencies beat slower.
*/

#include "../AudioApp.h"
#include "Gamma/rnd.h"
#include "Gamma/Filter.h"
#include "Gamma/Oscillator.h"
using namespace gam;

class MyApp : public AudioApp{
public:

	Accum<> tmr;		// Triggers new frequency value
	Osc<> osc1, osc2;	// Table oscillators (any oscillator will do)
	int count;

	MyApp(){
		tmr.period(4);
		tmr.phaseMax();
		osc1.resize(512);
		osc1.addSine(1, 0.25);
		osc1.addSine(2, 0.50);
		osc1.addSine(3, 0.25);
		osc2.source(osc1);
		count = 0;
	}

	void onAudio(AudioIOData& io){
		while(io()){

			if(tmr()){
				// Frequency needs to be divided by the middle harmonic number
				float f = 110 * pow(2, rnd::uni(0,12)/12.);

				// Select a new detune amount every 4 notes
				float detunings[] = {1.005, 1.01, 1.015, 1.002};
				float detune = detunings[(count/4)%4];
				++count;

				osc1.freq(f);			// Base frequency
				osc2.freq(f * detune);	// Detune sharp
			}
			
			float s = osc1() + osc2();
				
			io.out(0) = io.out(1) = s * 0.2;
		}
	}
};

int main(){
	MyApp().start();
}
