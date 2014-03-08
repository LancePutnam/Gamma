/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Example: Spectral / Magnitude Thresholding
	
	This demonstrates spectral magnitude thresholding which can be used to 
	denoise or "clean" a signal.
*/

#include "Gamma/AudioIO.h"
#include "Gamma/DFT.h"
#include "Gamma/Noise.h"
#include "Gamma/Oscillator.h"
using namespace gam;

// STFT(winSize, hopSize, padSize, winType, spectralFormat)
STFT stft(4096, 4096/4, 0, HANN, MAG_PHASE);
NoiseWhite<> noise;
Saw<> saw;


void audioCB(AudioIOData& io){

	saw.freq(220);

	while(io()){

		// Our signal is a saw wave plus a little noise
		float s = saw()*0.8 + noise()*0.01;

		if(stft(s)){

			for(unsigned k=0; k<stft.numBins(); ++k){

				// Get the bin magnitude (the first bin element)
				float mag = stft.bin(k)[0];

				// If the bin magnitude is less than our threshold, then we
				// zero its magnitude. The assumption here is that noisy bins
				// will have a relatively small magnitude.
				if(mag < 0.0001) stft.bin(k)[0] = 0;
				
				// If we flip the comparison, then we keep only the noise.
				//if(mag > 0.0001) stft.bin(k)[0] = 0;
			}
		}
		
		// Get next resynthesized sample
		// (Comment this out to hear the original signal with noise.)
		s = stft();
		
		io.out(0) = s;
		io.out(1) = s;
	}
}


int main(){
	AudioIO io(256, 44100, audioCB, NULL, 2);
	Domain::master().spu(io.framesPerSecond());
	io.start();
	printf("Press 'enter' to quit...\n"); getchar();
}
