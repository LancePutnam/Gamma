/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Example:		Generator / Player
	Description:	Demonstration of playing back a sound file at different
					rates.
*/

#include "../examples.h"

Accum<> tmr(1);
Player<float, ipl::Cubic, tap::Wrap> player("sounds/water3.wav");
//Player<float, ipl::Cubic, tap::Wrap> player("recording.aif");

void audioCB(AudioIOData& io){

	while(io()){
///*	
		if(tmr()){
			float r = pow(2, rnd::uniS(1.));
			player.range(rnd::uni(0.5), rnd::uni(1.));
			player.rate(rnd::neg(r, 0.3));
			player.reset();	
		}

		float s = player() * 0.2;
		io.out(0) = io.out(1) = s;
//*/
/*
		// stereo playback
		float s1 = player.read(0) * 0.2;
		float s2 = player.read(1) * 0.2;
		player.advance();
		io.out(0) = s1;
		io.out(1) = s2;
*/
	}
}

RUN(audioCB);
