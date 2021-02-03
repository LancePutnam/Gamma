/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example:	Spectral Plucking
Author:		Lance Putnam, 2012

Description:
This demonstrates a randomized bin plucking effect. The idea is to apply a
downward ramp function across all bins that "travels" towards DC. By itself,
this creates a sweeping band-pass filter effect, so the envelope samples are
scrambled bin-wise to eliminate the pattern.
*/
#include "../AudioApp.h"
#include "Gamma/DFT.h"
#include "Gamma/Noise.h"
#include "Gamma/Oscillator.h"
using namespace gam;

class MyApp : public AudioApp{
public:

	// args: winSize, hopSize, padSize, winType, spectralFormat, auxBufs
	STFT stft{2048, 2048/4, 0, RECTANGLE, COMPLEX, 1};
	NoisePink<> src;
	Sweep<> pluckEnvs;

	MyApp(){
		// Attach the pluck envelope accumulator to the hop-rate domain
		stft.domainHop() << pluckEnvs;
	}

	void onAudio(AudioIOData& io){
		while(io()){
			float s = src()*4;

			if(stft(s)){
				// Set period of pluck envelopes, in seconds.
				// This is negative to get a downward ramp going from 1 to 0.
				pluckEnvs.period(-60);

				// Create a moving ramp function in the aux buffer
				float amp = pluckEnvs();
				for(unsigned k=0; k<stft.numBins(); ++k){
					float A = scl::wrap(amp + float(k)/stft.numBins());
					// This gives the envelopes an exponential-like decay
					for(int j=0; j<5; ++j) A*=A;
					stft.aux(0)[k] = A;
				}

				// Do a seeded random shuffle of the envelopes so we don't hear
				// a pattern.
				rnd::push(14122);
				rnd::permute(stft.aux(0), stft.numBins());
				rnd::pop();

				for(unsigned k=0; k<stft.numBins(); ++k){
					stft.bin(k) *= stft.aux(0)[k];
				}
			}

			io.out(0) = io.out(1) = stft();
		}
	}
};

int main(){
	MyApp().start();
}
