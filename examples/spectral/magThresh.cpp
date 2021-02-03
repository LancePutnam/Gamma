/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example:	Magnitude Thresholding
Author:		Lance Putnam, 2014

Description:
This demonstrates spectral magnitude thresholding which can be used to denoise
or "clean" a signal. The procedure involves zeroing bin magnitudes which are
less than some threshold value. The assumption is that noise is present as low
amplitude energy across all frequencies.
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
			float s = saw()*0.3 + noise()*0.04;

			if(stft(s)){

				for(unsigned k=0; k<stft.numBins(); ++k){

					// Get the bin magnitude (the first bin element)
					float mag = stft.bin(k)[0];

					// If the bin magnitude is less than our threshold, then we
					// zero its magnitude. The assumption here is that noisy bins
					// will have a relatively small magnitude.
					if(mag < 0.0004) stft.bin(k)[0] = 0;
				
					// If we flip the comparison, then we keep only the noise.
					//if(mag > 0.0001) stft.bin(k)[0] = 0;
				}
			}
		
			// Get next resynthesized sample
			// (Comment this out to hear the original signal with noise.)
			s = stft();
		
			io.out(0) =  io.out(1) = s;
		}
	}
};

int main(){
	MyApp().start();
}
