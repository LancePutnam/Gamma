/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example:	Noise Colors
Author:		Lance Putnam, 2015

Description:
This plays some different "colors" of noise. Each of the noise types has a 
frequency-dependent power spectrum according to power(f) = 1/f^n. For white 
noise, n=0 resulting in a flat spectrum. For pink and brown noise, n=1 and n=2, 
respectively producing progressively darker sounding noise. Violet noise has 
n=-2 giving it an bright, airy quality.
*/

#include "../AudioApp.h"
#include "Gamma/Noise.h"
#include "Gamma/Oscillator.h"
using namespace gam;

class MyApp : public AudioApp{
public:

	NoiseViolet<> violet;	//   f^2 noise generator
	NoiseWhite<> white;		// 1/f^0 noise generator
	NoisePink<> pink;		// 1/f^1 noise generator
	NoiseBrown<> brown;		// 1/f^2 noise generator
	Accum<> tmr{1./4};		// Timer for switching noise
	int type=0;				// Noise type

	void onAudio(AudioIOData& io){

		while(io()){
		
			if(tmr()) (++type)%=4;
			
			float s = 0;

			auto printType = [&](const char * msg){ if(tmr.cycled()) printf("%s\n",msg); };

			switch(type){
				case 0: s = violet();    printType("Violet noise"); break;
				case 1: s = white()*0.5; printType("White noise"); break;
				case 2: s = pink();      printType("Pink noise"); break;
				default:s = brown();     printType("Brown noise"); break;
			}

			io.out(0) = io.out(1) = s*0.2;
		}
	}
};

int main(){
	MyApp().start();
}
