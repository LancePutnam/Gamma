/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example:
Author:		Lance Putnam, 2012

Description:
This demonstrates how to use a Seg to generate smooth noise which is 
subsequently mapped to the frequency of an oscillator. The process is akin
to upsampling whereby we generate a noise sample at a lower rate and then use
the segment to produce an interpolated sample at the audio sample rate.
*/

#include "../AudioApp.h"
#include "Gamma/Envelope.h"
#include "Gamma/Noise.h"
#include "Gamma/Oscillator.h"
using namespace gam;

class MyApp : public AudioApp{
public:

	NoiseWhite<> noise;
	Seg<float, iplSeq::Linear> seg;
	Sine<> osc;

	void onAudio(AudioIOData& io){

		// Set frequency of interpolation
		seg.freq(20);

		while(io()){

			// seg will request a new sample from noise by calling its function
			// operator and then return the interpolated result.
			float s = seg(noise);

			// Map smooth noise to oscillator frequency
			osc.freq(1000 + 400*s);

			s = osc()*0.1;

			io.out(0) = io.out(1) = s;
		}
	}
};

int main(){
	MyApp().start();
}

