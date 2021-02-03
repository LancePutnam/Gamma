/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example:	Convolution
Author:		Lance Putnam, 2014

Description:
This demonstrates how to perform convolution using two STFTs. The procedure is 
to take the Fourier transform of each signal, multiply their (complex-valued) 
spectra, and then perform an inverse transform on the product. The STFTs need to
be zero-padded by at least the window size to prevent circular convolution.
*/
#include "../AudioApp.h"
#include "Gamma/DFT.h"
#include "Gamma/Noise.h"
#include "Gamma/Oscillator.h"
using namespace gam;

class MyApp : public AudioApp{
public:

	static const int N = 4096;
	// We need two STFTs with the same settings
	STFT stft1{N, N/1, N, RECTANGLE, COMPLEX};
	STFT stft2{N, N/1, N, RECTANGLE, COMPLEX};
	NoisePink<> noise;
	Impulse<> imp1, imp2, imp3;

	MyApp(){
		imp1.normalize(false);
		imp2.normalize(false);
		imp3.normalize(false);
		imp1.freq(220 * ::pow(2, 0./12));
		imp2.freq(220 * ::pow(2, 4./12));
		imp3.freq(220 * ::pow(2, 7./12));
		stft1.inverseWindowing(false);
	}

	void onAudio(AudioIOData& io){
		while(io()){

			float s1 = noise()*2;
			float s2 = imp1() + imp2() + imp3();

			stft2(s2);

			if(stft1(s1)){

				// Loop through all the bins
				for(unsigned k=0; k<stft1.numBins(); ++k){
					// Multiply complex-valued spectra
					stft1.bin(k) *= stft2.bin(k);
				}
			}

			// Get next resynthesized sample
			float s = stft1();

			io.out(0) = s;
			io.out(1) = s;
		}
	}
};

int main(){
	MyApp().start();
}

