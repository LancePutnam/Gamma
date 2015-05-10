/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example:	Spectral Freezing
Author:		Lance Putnam, 2015

Description:
This shows how to create a spectral freezing effect with an STFT. We capture
one window of audio, perform a forward transform to get its spectrum and then
perform resynthesis multiple times on this spectrum. The STFT is configured to 
compute frequency estimates in each bin so upon resynthesis we obtain a more 
pleasing and natural continuation of the sound, i.e., it will not sound buzzy.
*/
#include "../AudioApp.h"
#include "Gamma/DFT.h"
#include "Gamma/Noise.h"
#include "Gamma/Oscillator.h"
using namespace gam;

class MyApp : public AudioApp{
public:

	STFT stft;
	NoisePink<> src;
	Accum<> tmr;
	int captureCount;

	MyApp()
	:	// STFT(winSize, hopSize, padSize, winType, spectralFormat)
		stft(2048, 2048/4, 0, HANN, MAG_FREQ)
	{
		captureCount = 0;
		tmr.period(2);
	}

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

