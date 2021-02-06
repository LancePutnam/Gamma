/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example:	Beating
Author:		Lance Putnam, 2012

Description:
This demonstrates the beating effect produced by summing two sinusoids with 
equal amplitude and nearly equal frequency. The sum of two sinusoids with 
frequencies f1 and f2 is equivalent to the product of two sinusoids with 
frequencies (f1+f2)/2 and (f1-f2)/2.
*/
#include "../AudioApp.h"
#include "Gamma/Oscillator.h"
using namespace gam;

class MyApp : public AudioApp{
public:

	Sine<> osc1{220};
	Sine<> osc2{220 + 1};

	void onAudio(AudioIOData& io){

		while(io()){
			float s = (osc1() + osc2()) * 0.2;
			io.out(0) = io.out(1) = s;
		}
	}
};

int main(){
	MyApp().start();
}
