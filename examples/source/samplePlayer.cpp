/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example:	Sample Playback
Author:		Lance Putnam, 2012

Description:
Demonstration of playing back a sound file at different rates.
*/

#include "../AudioApp.h"
#include "Gamma/rnd.h"
#include "Gamma/Oscillator.h"
#include "Gamma/SamplePlayer.h"
using namespace gam;

class MyApp : public AudioApp{
public:

	Accum<> tmr;
	SamplePlayer<> player;	// Uses linear interpolation
	//SamplePlayer<float, ipl::Trunc> player; // Uses no interpolation
	//SamplePlayer<float, ipl::Cubic> player; // Uses cubic interpolation

	MyApp(){
		tmr.period(1);

		// Load sound file into buffer
		player.load("../../sounds/water3.wav");
	}

	void onAudio(AudioIOData& io){

		while(io()){

			if(tmr()){
				// Set playback range
				player.range(rnd::uni(0.5), rnd::uni(1.));

				// Set playback rate; negative rates play in reverse
				float r = pow(2, rnd::uniS(1.));
				player.rate(rnd::neg(r, 0.3));

				// Start sample again from beginning
				player.reset();	
			}

			//* Mono (left-channel) playback
			float s = player() * 0.2;
			//player.loop(); // This makes the sample loop
			io.out(0) = io.out(1) = s;
			//*/

			/* Stereo playback
			float s1 = player.read(0) * 0.2;
			float s2 = player.read(1) * 0.2;
			player.advance();
			io.out(0) = s1;
			io.out(1) = s2;
			//*/
		}
	}
};

int main(){
	MyApp().start();
}
