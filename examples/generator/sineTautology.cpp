/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Example:		Generator / Sine Tautology
	Description:	Demonstrates multiple ways of generating sinusoids.
*/

#include "../examples.h"

Osc<>        sineT1(440);			// Table look-up (full-wave)
TableSine<>  sineT2(440);			// Table look-up (quarter-wave)
Quadra<>     sineQd(440);			// Complex mul
LFO<>        sineC1(440, -0.25f);	// Computed 3rd order poly
Sine<>       sineC2(440);			// Computed taylor
SineR<>      sineRs(440);			// Real-valued resonator
AccumPhase<> sineC3(440);			// Direct math.h sin() (ground truth)


void audioCB(AudioIOData& io){

	while(io()){
		float s	= sineT1()
				+ sineT2.nextL()
				+ sineQd().i
				+ sineC1.cos()
				+ sineC2()
				+ sineRs()
		;

		s = s/6 - sin(sineC3.nextPhase());	// pass only the artifacts
		
		io.out(0) = io.out(1) = s * 0.2f;
	}

}

int main(int argc, char* argv[]){
	AudioIO io(256, 44100., audioCB, NULL, 2);
	Sync::master().spu(io.framesPerSecond());

	sineT1.addSine(1);
	
	io.start();
	printf("\nPress 'enter' to quit...\n"); getchar();
	return 0;
}
