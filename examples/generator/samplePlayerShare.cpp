/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Example:		Generator / SamplePlayer
	Description:	This shows how two sample players can be made to share the
					same sample buffer loaded from a sound file.
*/

#include "../examples.h"

Accum<> tmr(1);
SamplePlayer<> player1("../../sounds/water3.wav");
SamplePlayer<> player2;

void audioCB(AudioIOData& io){

	// Assign player1's buffer to player2
	player2.buffer(player1);

	while(io()){
		if(tmr()){
			player1.rate(rnd::uni(0.5,1.));
			player2.rate(rnd::uni(0.5,1.));
			player1.reset();
			player2.reset();
		}

		float s = (player1() + player2()) * 0.2;
		io.out(0) = io.out(1) = s;
	}
}

RUN_AUDIO_MAIN
