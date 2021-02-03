/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example:	Spectral Peak Detection
Author:		Lance Putnam, 2014

Description:
This demonstrates a simple peak detection algorithm in the frequency domain.
A peak is said to be present in a bin if its magnitude is larger than its two
nearest neighbors.
*/
#include "../AudioApp.h"
#include "Gamma/DFT.h"
#include "Gamma/Noise.h"
#include "Gamma/Oscillator.h"
using namespace gam;

class MyApp : public AudioApp{
public:

	// args: winSize, hopSize, padSize, winType, spectralFormat
	STFT stft{4096, 4096/4, 0, HANN, MAG_PHASE};
	NoiseWhite<> noise;
	Saw<> saw{220};

	void onAudio(AudioIOData& io){
		while(io()){
			// Our signal is a saw wave plus a little noise
			float s = saw()*0.3 + noise()*0.05;

			// Input next sample for analysis
			// When this returns true, then we have a new spectral frame
			if(stft(s)){

				// Loop through the bins
				for(unsigned k=1; k<stft.numBins()-1; ++k){
			
					// Get neighborhood of magnitudes
					float m0 = stft.bin(k-1)[0];
					float m1 = stft.bin(k  )[0];
					float m2 = stft.bin(k+1)[0];

					// Is the current bin a peak?
					if(m1 > m0 && m1 > m2){

						// Zero the peak (for demonstration purposes)
						stft.bin(k-1)[0] = 0;
						stft.bin(k  )[0] = 0;
						stft.bin(k+1)[0] = 0;

						k++; // Skip the next bin---it cannot be a peak
					}
				}
			}
		
			// Get next resynthesized sample
			s = stft();
		
			io.out(0) = io.out(1) = s;
		}
	}
};

int main(){
	MyApp().start();
}

