/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Example:		Generator / Sine Tautology
	Description:	Demonstrates multiple ways of generating sinusoids.
*/

#include "../examples.h"

Osc<>        sineT1(440);			// Table look-up
CSine<>      sineCS(440);			// Complex mul
LFO<>        sineC1(440, -0.25f);	// Computed 3rd order poly
Sine<>       sineC2(440);			// Computed taylor
SineR<>      sineRs(440);			// Real-valued resonator
AccumPhase<> sineC3(440);			// Direct math.h sin() (ground truth)


void audioCB(AudioIOData& io){

	while(io()){
		float s	= sineT1()
				+ sineCS().i
				+ sineC1.cos()
				+ sineC2()
				+ sineRs()
		;

		s = s/5. - sin(sineC3.nextPhase());	// pass only the artifacts
		
		io.out(0) = io.out(1) = s * 0.2f;
	}

}

int main(int argc, char* argv[]){
	sineT1.addSine(1);

	AudioIO io(256, 44100., audioCB, NULL, 2);
	Domain::master().spu(io.framesPerSecond());	
	io.start();
	printf("\nPress 'enter' to quit...\n"); getchar();
	return 0;
}
