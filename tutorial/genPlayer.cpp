/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Tutorial:		Generator / Player
	Description:	Demonstration of playing back a sound file at different
					rates.
*/

#include "tutorial.h"

Accum<> tmr(1);
Player<> player("sounds/water4.wav");
//Player<float, ipl::Trunc, tap::Wrap> player("sounds/water4.wav", 0.1);

void audioCB(AudioIOData& io){

	for(uint32_t i=0; i<io.framesPerBuffer(); i++){
	
		if(tmr()){
			float r = pow(2, rnd::uniS(1.));

			// play forward
			if(rnd::prob(0.7)){
				player.rate(r);
				player.phase(0);
			}
			
			// play backward
			else{
				player.rate(-r);
				player.phase(1);			
			}
		}

		float s = player() * 0.2;

		io.out(0)[i] = io.out(1)[i] = s;
	}
}

RUN(audioCB);
