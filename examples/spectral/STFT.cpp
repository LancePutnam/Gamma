/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example:	Short-time Fourier Transform (STFT)
Author:		Lance Putnam, 2012

Description:
This demonstrates how to use the STFT class to do frequency-domain analysis and 
processing and resynthesis back into the time domain.
*/

#include "../AudioApp.h"
#include "Gamma/DFT.h"
#include "Gamma/Oscillator.h"
using namespace gam;

class MyApp : public AudioApp{
public:

	// Short-time Fourier transform
	STFT stft{
		2048,		// Window size
		2048/4,		// Hop size; number of samples between transforms
		0,			// Pad size; number of zero-valued samples appended to window
		HANN,		// Window type: BARTLETT, BLACKMAN, BLACKMAN_HARRIS,
					//		HAMMING, HANN, WELCH, NYQUIST, or RECTANGLE
		COMPLEX		// Format of frequency samples:
					//		COMPLEX, MAG_PHASE, or MAG_FREQ
	};

	// Test tone
	Sine<> src;

	void onAudio(AudioIOData& io){
		while(io()){
			src.freq(220);
			float s = src();

			// Input next sample for analysis
			// When this returns true, then we have a new spectral frame
			if(stft(s)){
			
				// Loop through all the bins
				for(unsigned k=0; k<stft.numBins(); ++k){
					// Here we simply scale the complex sample
					stft.bin(k) *= 0.2;
				}
			}
		
			// Get next resynthesized sample
			s = stft();
		
			io.out(0) = s;
			io.out(1) = s;
		}
	}
};

int main(){
	MyApp().start();
}
