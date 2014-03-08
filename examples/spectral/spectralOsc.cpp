/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Example: Spectral / Spectral Oscillators
	
	This demonstrates how to use oscillators running across the frequency domain
	to create time-varying filtering patterns.
*/

#include "Gamma/AudioIO.h"
#include "Gamma/DFT.h"
#include "Gamma/Noise.h"
#include "Gamma/Oscillator.h"
using namespace gam;

NoisePink<> src;

// STFT(winSize, hopSize, padSize, winType, spectralFormat)
STFT stft(2048, 2048/4, 0, HANN, COMPLEX);

// Spectral magnitude oscillators
LFO<> oscMag1, oscMag2;

// LFOs used to vary the starting phase of the magnitude oscillators
LFO<> modPhase1, modPhase2;


void audioCB(AudioIOData & io){

	while(io()){
		float s = src() * 0.5f;
		
		if(stft(s)){

			// Set periods of phase modulators, in seconds
			modPhase1.period(1.33);
			modPhase2.period(1.61);
			
			// Set periods of magnitude oscillators, in Hz
			oscMag1.period(1000);
			oscMag2.period( 800);

			// Modulate phases according to raised-cosine wave
			oscMag1.phase(modPhase1.hann());
			oscMag2.phase(modPhase2.hann());
			
			for(unsigned i=0; i<stft.numBins(); ++i){
			
				// Compute the gain to apply to this bin
				float m  = oscMag1.triU() - oscMag2.triU();
			
				stft.bin(i) *= m;
			}	
		}

		io.out(0) = io.out(1) = stft();
	}
}

int main(){
	// Attach the phase modulation LFOs to the hop-rate domain
	stft.domainHop() << modPhase1 << modPhase2;

	// Attach objects to frequency domain
	stft.domainFreq() << oscMag1 << oscMag2;

	AudioIO io(128, 44100., audioCB, NULL, 2);
	Domain::master().spu(io.framesPerSecond());		
	io.start();
	printf("\nPress 'enter' to quit...\n"); getchar();
}
