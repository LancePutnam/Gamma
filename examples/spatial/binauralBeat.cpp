/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example:	Binaural Beating
Author:		Lance Putnam, 2012

Description:
Demonstration of a binaural beating effect. This effect is produced by
outputting two sine waves with slightly different frequencies to the left and
right channels. With headphones on, the result sounds like a single sine wave
rotating around one's head at a rate that is the difference between the two sine
wave frequencies.
*/

#include "../AudioApp.h"
#include "Gamma/Oscillator.h"
using namespace gam;

class MyApp : public AudioApp{
public:

	Sine<> sineLeft, sineRight;

	void onAudio(AudioIOData& io){

		// Set frequencies to be 1 Hz apart
		sineLeft .freq(220);
		sineRight.freq(221);

		while(io()){

			// Output sine waves on left and right channels
			io.out(0) = sineLeft() * 0.2;
			io.out(1) = sineRight()* 0.2;
		}
	}

};

int main(){
	MyApp().start();
}
