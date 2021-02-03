/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example:	Brick-wall Sweep
Author:		Lance Putnam, 2012

Description:
This shows how to create a sweeping brick-wall low-pass filter using a
modulator and an STFT. The modulator runs in the STFT hop rate domain since we
only need to compute its value when we get a new spectral frame.
*/

#include "../AudioApp.h"
#include "Gamma/DFT.h"
#include "Gamma/Noise.h"
#include "Gamma/Oscillator.h"
using namespace gam;

class MyApp : public AudioApp{
public:

	// args: winSize, hopSize, padSize, winType, spectralFormat
	STFT stft{2048, 2048/4, 0, HANN, COMPLEX};
	NoisePink<> src;
	LFO<> modCutoff{1./16, 0.5}; // args: freq, phase

	MyApp(){
		// Attach the cutoff frequency modulator to the hop-rate domain
		stft.domainHop() << modCutoff;
	}

	void onAudio(AudioIOData& io){
		while(io()){
			float s = src() * 0.4;

			if(stft(s)){

				// Cutoff frequency as fraction of sample rate
				float frac = scl::pow3( modCutoff.triU() );

				// Compute bin number at cutoff frequency
				int N = stft.numBins();
				int kCutoff = frac * N;
			
				// Zero bins higher than the cutoff bin
				for(int k=kCutoff; k<N; ++k){
					stft.bin(k) = 0;
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

