/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example:	Hilbert Transform
Author:		Lance Putnam, 2012

Description:
This demonstrates how one can frequency shift a signal by using a Hilbert 
transform. The Hilbert transform converts a real signal into a complex signal.
If we multiply a complex signal by a complex sinusoid, we perform a frequency
shift or single-sideband amplitude modulation.
*/

#include "../AudioApp.h"
#include "Gamma/Filter.h"
#include "Gamma/Oscillator.h"
using namespace gam;

class MyApp : public AudioApp{
public:

	LFO<> osc;
	LFO<> shiftMod;
	Hilbert<> hil;
	CSine<> shifter;

	MyApp(){
		osc.freq(100);
		shiftMod.period(10);
		shifter.freq(100);
	}

	void onAudio(AudioIOData& io){

		while(io()){
			float s = osc.up()*0.1;

			// The Hilbert transform returns a complex number
			Complex<float> c = hil(s);
			
			// Perform a frequency shift
			shifter.freq(shiftMod.hann()*1000);
			c *= shifter();
			
			// Output the real and imaginary components
			float s0 = c.r;
			float s1 = c.i;
		
			io.out(0) = s0;
			io.out(1) = s1;
		}
	}
};

int main(){
	MyApp().start();
}

