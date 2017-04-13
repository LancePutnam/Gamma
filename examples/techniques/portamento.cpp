/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example:	Portamento
Author:		Lance Putnam, 2012

Description:
This demonstrates how to implement portamento using a one-pole filter.
One-pole filters can be used to smooth out control signals--in this case random 
stepping of the frequency of an oscillator.
*/

#include "../AudioApp.h"
#include "Gamma/rnd.h"
#include "Gamma/Filter.h"
#include "Gamma/Oscillator.h"
using namespace gam;

class MyApp : public AudioApp{
public:

	Accum<> tmr;		// Triggers new frequency value
	Sine<> src;			// Sine oscillator
	OnePole<> freq;

	void onAudio(AudioIOData& io){

		tmr.period(0.25);

		// Set time taken to reach new frequency value
		freq.lag(0.1);

		while(io()){

			if(tmr()){
				// Set new target frequency of one-pole
				freq = pow(2, rnd::uniS(1.))*440;
			}

			// Use smoothed output of one-pole for oscillator frequency
			src.freq(freq());
			
			float s = src();
				
			io.out(0) = io.out(1) = s * 0.2f;
		}
	}
};

int main(){
	MyApp().start();
}
