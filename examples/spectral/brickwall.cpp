/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Example: Spectral / Brick-wall
	
	This shows how to create a brick-wall bandpass filter using an STFT.
*/

#include "Gamma/AudioIO.h"
#include "Gamma/DFT.h"
#include "Gamma/Noise.h"
using namespace gam;

// STFT(winSize, hopSize, padSize, winType, spectralFormat)
STFT stft(2048, 2048/4, 0, HANN, COMPLEX);
NoisePink<> src;

void audioCB(AudioIOData& io){

	while(io()){

		float s = src();

		if(stft(s)){
		
			// Define the band edges, in Hz
			float freqLo =  400;
			float freqHi = 1800;

			for(unsigned k=0; k<stft.numBins(); ++k){
				
				// Compute the frequency, in Hz, of this bin
				float freq = k*stft.binFreq();

				// If the bin frequency is outside of our band, then zero
				// the bin.
				if(freq < freqLo || freq > freqHi){
					stft.bin(k) = 0;
				}
			}
		}

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
