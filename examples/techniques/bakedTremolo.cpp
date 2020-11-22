/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example:	Baked Tremolo
Author:		Lance Putnam, 2012

Description:
The shows how to "bake" a tremolo into a wavetable to save on having an
additional modulator oscillator. The only caveat is that the tremolo rate is 
always an integer division of the oscillator frequency. However, this can
produce very rich modulation when oscillators are combined in parallel.
*/

#include "../AudioApp.h"
#include "Gamma/rnd.h"
#include "Gamma/Filter.h"
#include "Gamma/Oscillator.h"
using namespace gam;

class MyApp : public AudioApp{
public:

	Accum<> tmr;			// Triggers new frequency value
	Osc<> osc1, osc2, osc3;	// Table oscillators

	MyApp(){
		tmr.period(4);
		tmr.phaseMax();
		osc1.resize(1024);
		osc1.addSine(63, 0.2);	// Lower sideband
		osc1.addSine(64, 0.6);	// Carrier
		osc1.addSine(65, 0.2);	// Upper sideband
		osc2.source(osc1);
		osc3.source(osc1);
	}

	void onAudio(AudioIOData& io){
		while(io()){

			if(tmr()){
				// Frequency needs to be divided by the middle harmonic number
				float f = 440. / 64;

				// Create random chord
				osc1.freq(f*pow(2, rnd::uni(-2, 1)/12.));
				osc2.freq(f*pow(2, rnd::uni( 3, 5)/12.));
				osc3.freq(f*pow(2, rnd::uni( 5, 8)/12.));
			}
			
			float s = osc1() + osc2() + osc3();
				
			io.out(0) = io.out(1) = s * 0.2;
		}
	}
};

int main(){
	MyApp().start();
}
