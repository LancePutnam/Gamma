/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example:	Spectral Freezing
Author:		Lance Putnam, 2015

Description:
This shows how to create a spectral freezing effect with an STFT. We capture
one window of audio, perform a forward transform to get its spectrum and then
perform resynthesis multiple times on this spectrum. The STFT is configured to 
compute frequency estimates in each bin so upon resynthesis we obtain a more 
pleasing and natural (non-buzzy) continuation of the sound.
*/
#include "../AudioApp.h"
#include "Gamma/DFT.h"
#include "Gamma/Noise.h"
#include "Gamma/Oscillator.h"
using namespace gam;

class MyApp : public AudioApp{
public:

	// args: winSize, hopSize, padSize, winType, spectralFormat
	STFT stft{2048, 2048/4, 0, HANN, MAG_FREQ};
	NoisePink<> src;
	Accum<> tmr{1./2.};
	int captureCount = 0;

	void onAudio(AudioIOData& io){
		while(io()){

			if(tmr()){
				captureCount=0;
			}

			float s = src()*0.7;

			// Capture 4 hops to make up the size of one window
			if(captureCount<4 && stft(s)){
				captureCount++;
			}

			// Resynthesis based on internal spectral data
			s = stft();
		
			io.out(0) = io.out(1) = s;
		}
	}
};

int main(){
	MyApp().start();
}

