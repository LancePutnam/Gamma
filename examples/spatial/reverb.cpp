/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example:	Reverb
Author:		Lance Putnam, 2015

Description:
This demonstrates how to add reverberation, in the form of late reflections, to
a sound. The particular reverberation model used here is the Schroeder model
which consists of a parallel-comb, series-allpass configuration of delay lines.
*/

#include "../AudioApp.h"
#include "Gamma/Oscillator.h"
#include "Gamma/Spatial.h"
using namespace gam;

class MyApp : public AudioApp{
public:

	SineD<> src;		// Decaying sine wave grain
	Accum<> tmr;		// Timer for firing grains
	unsigned seed;		// RNG seed
	ReverbMS<> reverb;	// Schroeder reverberator
	//ReverbMS<float, Loop1P1Z> reverb; // Use a 1-pole, 1-zero loop filter

	MyApp(){
		seed = 1;
		tmr.period(2);
		tmr.finish();

		// Specify a reverb preset
		//reverb.resize(gam::JCREVERB); // Chowning (4-comb, 3-allpass)
		reverb.resize(gam::FREEVERB); // Jezar (8-comb, 4-allpass)

		// Or invent your own...
		//reverb.resize({1289, 2951, 2013}, {499, 951}); // pretty smooth
		//reverb.resize({3229, 3281, 4833}, {4487, 4097}); // crunchy
		//reverb.resize({3387, 3255, 3121}, {583, 491}); // colored
		//reverb.resize({3431, 4403, 5813}, {1645, 1299});

		//reverb.resize({1323, 1807, 1945, 2045}, {1559, 797, 297}); // pretty smooth
		//reverb.resize({1305, 1503, 1581, 1837}, {709, 535, 237}); // also smooth
		//reverb.resize({3013, 3799, 3479, 1799}, {335, 689, 907}); // less diffuse

		// Set decay length, in seconds
		reverb.decay(8);

		// Set high-frequency damping factor in [0, 1]
		reverb.damping(0.2);
	}

	void onAudio(AudioIOData& io){
		while(io()){

			// Trigger a new random grain on a timer
			if(tmr()){
				seed *= 69069;
				float f = (seed%4000) + 200;
				src.set(f, 1, 16./f);				
			}

			// Get next grain sample
			float s = src();

			// Pass grain through reverberator and mix with dry source
			s += reverb(s) * 0.2;

			io.out(0) = s;
			io.out(1) = s;
		}
	}
};

int main(){
	MyApp().start();
}

