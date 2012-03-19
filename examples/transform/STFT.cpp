/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Example:		Transform / STFT
	Description:	Sweeping brick-wall filter
*/

#include "../examples.h"

STFT stft(
	2048,			// Window size
	2048/4,			// Hop size; interval between transforms
	0,				// Pad size; zero padding amount on window
	HANN,			// Window type: BARTLETT, BLACKMAN, BLACKMAN_HARRIS,
					// HAMMING, HANN, WELCH, NYQUIST, or RECTANGLE
	COMPLEX			// Format of frequency samples
);

NoisePink<> src;
LFO<> edge(1./16, 0.5);


void audioCB(AudioIOData& io){

	while(io()){

		float s = src();

		// Input next sample for analysis
		// When this returns true, then we have a new spectral frame
		if(stft(s)){
		
			float frac = scl::pow3( edge.triU() );
			int N = stft.numBins();
			
			for(int k=0; k<N; ++k){
				int indKnee = frac * N;
				stft.bin(k) *= k < indKnee ? 1:0;
			}
		}
		
		// Get next resynthesized sample
		s = stft() * 0.2;
		
		io.out(0) = s;
		io.out(1) = s;
	}
}


int main(){

	stft.syncHop() << edge;

	AudioIO io(256, 44100, audioCB, NULL, 2);
	Sync::master().spu(io.framesPerSecond());
	io.start();
	printf("Press 'enter' to quit...\n"); getchar();
}
