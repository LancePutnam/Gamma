/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example:	Periodic Timing
Author:		Lance Putnam, 2015

Description:
This demonstrates how to use a phase accumulator as a periodic timer to trigger
events.
*/

#include "../AudioApp.h"
#include "Gamma/Oscillator.h"
using namespace gam;

class MyApp : public AudioApp{
public:

	Accum<> tmr;	// Phase accumulator

	void onAudio(AudioIOData& io){

		// Set period of timer, in seconds
		tmr.period(1);

		while(io()){

			float s = 0;

			// Accum's function operator returns true after each period
			if(tmr()){
				s = 0.2;
				std::cout << "Click!\n";
			}
			
			io.out(0) = io.out(1) = s;
		}
	}

};

int main(){
	MyApp().start();
}
