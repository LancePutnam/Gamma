/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Example:		Generator / Oscillator / Osc
	Description:	Using a general-purpose table-based oscillator.
*/

#include "../examples.h"

float ff = 110;						// Fundamental frequency
Accum<> tmr(0.2, 2);				// Timer for modifying wavetable contents
Osc<> osc1(ff       , 0, 512);		// Oscillator owning a 512-element wavetable
Osc<> osc2(ff + 0.17, 0, osc1);		// Detuned oscillator sharing osc1's table
Osc<> osc3(ff + 0.03, 0, osc1);		// Detuned oscillator sharing osc1's table

void audioCB(AudioIOData& io){

	while(io()){
	
		if(tmr()){
			osc1.zero();			// zero wavetable
			osc1.addSine(1,1);		// add the fundamental
		
			// add some 1/f overtones
			int skip = rnd::uni(5, 1);
			for(int i=1+skip; i<=rnd::uni(65, 16); i+=skip){
				if(rnd::prob()) osc1.addSine(i, 1./i, 0);
			}
		}

		float s = (osc1() + osc2()/2 + osc3()/4) * 0.2;
		
		io.out(0) = io.out(1) = s;
	}
}

RUN(audioCB);
