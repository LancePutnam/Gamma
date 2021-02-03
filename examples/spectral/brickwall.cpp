/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example:	Brick-wall Filter
Author:		Lance Putnam, 2014

Description:
This shows how to create a brick-wall bandpass filter using an STFT.
*/
#include "../AudioApp.h"
#include "Gamma/DFT.h"
#include "Gamma/Noise.h"
using namespace gam;

class MyApp : public AudioApp{
public:

	// args: winSize, hopSize, padSize, winType, spectralFormat)
	STFT stft{2048, 2048/4, 0, HANN, COMPLEX};
	NoisePink<> src;

	void onAudio(AudioIOData& io){
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
};

int main(){
	MyApp().start();
}
