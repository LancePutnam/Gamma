/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example:	Shared Sample Playback
Author:		Lance Putnam, 2014

Description:
This shows how two sample players can be made to share the same sample buffer 
loaded from a sound file.
*/

#include "../AudioApp.h"
#include "Gamma/SamplePlayer.h"
using namespace gam;

class MyApp : public AudioApp{
public:

	SamplePlayer<> player1, player2;

	MyApp(){
		// Load sound file (only once!) into buffer
		player1.load("../../sounds/water3.wav");

		// Assign player1's buffer to player2;
		// this can be called safely in the audio thread.
		player2.buffer(player1);

		// Make second playback rate slightly higher to create "phasing"
		player2.rate(1.005);
	}

	void onAudio(AudioIOData& io){
		while(io()){
			float s = (player1() + player2()) * 0.5;
			player1.loop();
			player2.loop();
			io.out(0) = io.out(1) = s;
		}
	}
};

int main(){
	MyApp().start();
}
