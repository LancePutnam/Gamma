/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example:	Table Lookup Oscillator
Author:		Lance Putnam, 2015

Description:
This shows how to use Osc, a table lookup oscillator.
*/

#include "../AudioApp.h"
#include "Gamma/rnd.h"
#include "Gamma/Oscillator.h"
using namespace gam;

class MyApp : public AudioApp{
public:

	Accum<> tmr;
	Osc<> osc1, osc2, osc3;

	MyApp(){
		tmr.phase(0.999);
		osc1.resize(512);	// Set table size
		osc2.source(osc1);	// Have osc2 to use osc1's table
		osc3.source(osc1);	// Have osc3 to use osc1's table
	}

	void onAudio(AudioIOData& io){

		// Set period of timer, in seconds
		tmr.period(5);

		while(io()){
			if(tmr()){

				// Set oscillator frequencies
				float fund = 55 * pow(2, rnd::uni(12)/12.);
				osc1.freq(fund);
				osc2.freq(fund * 1.005);
				osc3.freq(fund * 0.991);
				//osc2.freq(fund * 1.498);
				//osc3.freq(fund * 1.260);

				osc1.zero();			// zero wavetable
				osc1.addSine(1,1);		// add the fundamental
			
				// add some 1/f overtones
				int skip = rnd::uni(5, 1);
				int maxh = rnd::uni(256, 16);
				for(int i=1+skip; i<=maxh; i+=skip){
					if(rnd::prob()) osc1.addSine(i, 1./i, 0);
				}
			}

			float s = (osc1() + osc2() + osc3()) * 0.2;
			
			io.out(0) = io.out(1) = s;
		}
	}

};

int main(){
	MyApp().start();
}
