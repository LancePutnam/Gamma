/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Tutorial:		Generator / Oscillator / Additive
	Description:	Wavepacket with dispersion.
*/

#include "tutorial.h"

float ff = 27.5;
SineRs<> oscs(512);

void audioCB(AudioIOData& io){

	// sample-based generation
	for(uint32_t i=0; i<io.framesPerBuffer(); i++){

		float s = oscs()*4;
		io.out(0)[i] = io.out(1)[i] = s;
	}
}


int main(int argc, char* argv[]){

	AudioIO io(256, 44100, audioCB, NULL, 2);
	Sync::master().spu(io.framesPerSecond());

	float k1 = 0.;					// velocity
	float k2 = 0.0001;				// dispersion
	float rs = 1./oscs.size();

	for(unsigned i=0; i<oscs.size(); ++i){
	
		float h = i+1 + i*k1 + i*i*k2;		
		float a = 16 * rs/h;
		float p = 0;
	
		oscs.set(i,  ff*h,a,p);
	}

	io.start();
	printf("\nPress 'enter' to quit...\n"); getchar();
	return 0;
}
