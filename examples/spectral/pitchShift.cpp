/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example:	Pitch Shift
Author:		Lance Putnam, 2014

Description:
This shows how to shift the pitch of a sound without changing its duration.
For improved results, we use frequency estimates for resynthesis and reset
the phases whenever a transient is detected (defined here as a large positive
change in spectral flux).
*/
#include "../AudioApp.h"
#include "Gamma/DFT.h"
#include "Gamma/SamplePlayer.h"
using namespace gam;

class MyApp : public AudioApp{
public:

	// args: winSize, hopSize, padSize, winType, sampType, auxBufs
	STFT stft{4096, 4096/4, 0, gam::HAMMING, gam::MAG_FREQ, 3};
	SamplePlayer<> play;

	MyApp(){
		play.load("../../sounds/count.wav");
	}

	void onAudio(AudioIOData& io){
		float pshift = 1.7831; //pshift = 1./pshift;

		while(io()){
			float s = play(); play.loop();

			if(stft(s)){
		
				enum{
					PREV_MAG=0,
					TEMP_MAG,
					TEMP_FRQ
				};

				// Compute spectral flux (L^1 norm on positive changes)
				float flux = 0;
				for(unsigned k=0; k<stft.numBins(); ++k){
					float mcurr = stft.bin(k)[0];
					float mprev = stft.aux(PREV_MAG)[k];

					if(mcurr > mprev){
						flux += mcurr - mprev;
					}
				}

				//printf("%g\n", flux);
				//gam::printPlot(flux); printf("\n");

				// Store magnitudes for next frame
				stft.copyBinsToAux(0, PREV_MAG);

				// Given an onset, we would like the phases of the output frame
				// to match the input frame in order to preserve transients.
				if(flux > 0.2){
					stft.resetPhases();
				}

				// Initialize buffers to store pitch-shifted spectrum
				for(unsigned k=0; k<stft.numBins(); ++k){
					stft.aux(TEMP_MAG)[k] = 0.;
					stft.aux(TEMP_FRQ)[k] = k*stft.binFreq();
				}

				// Perform the pitch shift:
				// Here we contract or expand the bins. For overlapping bins,
				// we simply add the magnitudes and replace the frequency.
				// Reference:
				// http://oldsite.dspdimension.com/dspdimension.com/src/smbPitchShift.cpp
				if(pshift > 0){
					unsigned kmax = stft.numBins() / pshift;
					if(kmax >= stft.numBins()) kmax = stft.numBins()-1;
					for(unsigned k=1; k<kmax; ++k){
						unsigned j = k*pshift;
						stft.aux(TEMP_MAG)[j] += stft.bin(k)[0];
						stft.aux(TEMP_FRQ)[j] = stft.bin(k)[1]*pshift;
					}
				}

				// Copy pitch-shifted spectrum over to bins
				stft.copyAuxToBins(TEMP_MAG, 0);
				stft.copyAuxToBins(TEMP_FRQ, 1);
			}

			s = stft();

			io.out(0) = io.out(1) = s;
		}
	}
};

int main(){
	MyApp().start();
}

