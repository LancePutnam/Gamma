/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Example:		Generator / Oscillator / Accumulator
	Description:	Using an accumulator as a timer to trigger impulses
*/

#include "../examples.h"

Accum<> tmr(10);			// timer to switch between LFO types

void audioCB(AudioIOData& io){

	while(io()){
	
		float s = 0;
	
		if(tmr()){
			tmr.period(rnd::pick(0.4, 0.2));
			s = 0.2;
		}
		
		io.out(0) = io.out(1) = s;
	}
}

RUN(audioCB);
