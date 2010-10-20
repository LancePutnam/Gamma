/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Example:		Transform / STFT
	Description:	Sweeping brick-wall filter
*/

#include "../examples.h"

STFT stft(
	2048,			// Window size
	2048/4,			// Hop size
	0,				// Pad size
	WinType::Hann,	// Window type: Bartlett, Blackman, BlackmanHarris, 
					// Hamming, Hann, Welch, Rectangle
	Bin::Rect		// Format of frequency samples
);

NoisePink<> src;
LFO<> edge(1./16, 0.5);


void audioCB(AudioIOData& io){

	while(io()){

		float s = src();

		if(stft(s)){
		
			float frac = scl::pow3( edge.triU() );
			int N = stft.numBins();
			
			for(int k=0; k<N; ++k){
				int indKnee = frac * N;
				stft.bins(k) *= k < indKnee ? 1:0;
			}
		}
		
		s = stft() * 0.2;
		
		io.out(0) = s;
		io.out(1) = s;
	}
}


int main(int argc, char* argv[]){

	stft.syncHop() << edge;

	AudioIO io(256, 44100, audioCB, NULL, 2);
	Sync::master().spu(io.framesPerSecond());
	io.start();
	printf("Press 'enter' to quit...\n"); getchar();
	return 0;
}
