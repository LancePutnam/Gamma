/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Example:		Effect / Chebyshev Waveshaping
	Description:	This demonstrates how to waveshape a sine wave into a sum
					of cosine waves using Chebyshev polynomials of the first 
					kind.
*/

#include "../examples.h"

Sine<> src(110);		// Sine wave
ChebyN<12> cheby;		// Chebyshev waveshaper
LFO<> mod;				// Harmonic coefficient modulator

void audioCB(AudioIOData& io){

	mod.freq((io.fps()*0.5 + 0.1) / cheby.size());

	while(io()){

		// Create harmonic traveling wave
		for(unsigned k=0; k<cheby.size(); ++k){
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

RUN(audioCB);
