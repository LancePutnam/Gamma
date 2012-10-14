/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Example:		Effect / Panning
	Description:	A simple two-channel panning effect.
*/

#include "../examples.h"

Sine<> src(220);		// Sine wave
Pan<> pan;
LFO<> panMod(0.5);

void audioCB(AudioIOData& io){

	while(io()){

		// Generate our mono signal
		float s = src() * 0.2;

		// Modulate pan position between left (-1) and right (1) channels		
		pan.pos(panMod.tri());

		// The output is two floats (stereo)
		float2 xy = pan(s);

		io.out(0) = xy.x;
		io.out(1) = xy.y;
	}
}

RUN_AUDIO_MAIN
