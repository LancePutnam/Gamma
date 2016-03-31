/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example:	Chebyshev Waveshaping
Author:		Lance Putnam, 2012

Description:
This demonstrates how to waveshape a sine wave into a sum of harmonic cosine
waves using Chebyshev polynomials of the first kind.
*/

#include "../AudioApp.h"
#include "Gamma/Effects.h"
#include "Gamma/Oscillator.h"
using namespace gam;

class MyApp : public AudioApp{
public:

	Sine<> src;				// Sine wave
	ChebyN<12> cheby;		// Chebyshev waveshaper with 12 harmonics
	LFO<> mod;				// Harmonic coefficient modulator

	void onAudio(AudioIOData& io){

		src.freq(110);
		mod.freq((io.fps()*0.5 + 0.1) / cheby.size());

		while(io()){

			// Create harmonic traveling wave
			for(unsigned k=0; k<cheby.size(); ++k){

				// Set amplitude of kth harmonic
				cheby.coef(k) = mod.para();
			}

			// The waveshaper takes a sine wave in [-1, 1]
			float s = src();

			// Generate harmonics
			s = cheby(s);
			
			// Divide by number of harmonics to prevent clipping
			s /= cheby.size();

			io.out(0) = io.out(1) = s;
		}
	}

};

int main(){
	MyApp().start();
}
