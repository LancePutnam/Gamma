/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Example:		Generator / Oscillator / Additive
	Description:	This demonstrates the beating effect produced by summing two
					sinusoids with equal amplitude and nearly equal frequency.
					The sum of two sinusoids with frequencies f1 and f2 is 
					equivalent to the product of two sinusoids with frequencies
					(f1+f2)/2 and (f1-f2)/2.
*/

#include "../examples.h"

float ff = 220;
float freqBeat = 1;

Sine<> osc1(ff);				// sinusoid 1
Sine<> osc2(ff + freqBeat);		// sinusoid 2

void audioCB(AudioIOData& io){

	while(io()){

		float s = (osc1() + osc2()) * 0.1;

		io.out(0) = io.out(1) = s;
	}
}

RUN_AUDIO_MAIN
