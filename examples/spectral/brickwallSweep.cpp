/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Example: Spectral / Brick-wall Sweep
	
	This shows how to create a sweeping brick-wall low-pass filter using an
	oscillator and an STFT.
*/

#include "Gamma/AudioIO.h"
#include "Gamma/DFT.h"
#include "Gamma/Noise.h"
#include "Gamma/Oscillator.h"
using namespace gam;

// STFT(winSize, hopSize, padSize, winType, spectralFormat)
STFT stft(2048, 2048/4, 0, HANN, COMPLEX);
NoisePink<> src;
LFO<> modCutoff(1./16, 0.5);


void audioCB(AudioIOData& io){

	while(io()){

		float s = src() * 0.2;

		if(stft(s)){
		
			float frac = scl::pow3( modCutoff.triU() );

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


int main(){
	// Attach the cutoff frequency modulator to the hop-rate domain
	stft.domainHop() << modCutoff;

	AudioIO io(256, 44100, audioCB, NULL, 2);
	Domain::master().spu(io.framesPerSecond());
	io.start();
	printf("Press 'enter' to quit...\n"); getchar();
}
