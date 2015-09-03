/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example:	Distance Coder
Author:		Lance Putnam, 2015

Description:
This shows how to use the Dist class to spatialize a sound source. Dist 
filters a sound based on its distance from a given receiver/listener. The
distance cues are interaural level difference (ILD), interaural time difference
(ITD), and high-frequency damping from air absorption.
*/

#include "../AudioApp.h"
#include "Gamma/Oscillator.h"
#include "Gamma/Spatial.h"
using namespace gam;

class MyApp : public AudioApp{
public:

	Accum<> tmr;		// Timer to trigger sound
	SineD<> src;		// Sound source (decaying sine wave)
	Dist<2> dist;		// Filters source at 2 destinations based on distance cues
	LFO<> pathx, pathy;	// Sound source path

	MyApp(){
		pathx.period(4);
		pathy.period(5);
		tmr.period(0.17);
		src.set(2000, 1, 0.01);

		// Set distances between which sound is heard; these are the defaults
		dist.near(0.1);
		dist.far(10);

		// Set speed of sound
		dist.speedOfSound(343.2); // air (default)
		//dist.speedOfSound(1500); // water
	}

	void onAudio(AudioIOData& io){
		while(io()){

			// Synthesize source sound
			if(tmr()) src.reset();
			float s = src();

			// Position of source w.r.t. origin
			float x = pathx.cos()*4;
			float y = pathy.cos()*4;

			// Set distances from source to destinations (ears)
			// The ears are set to be 20 cm apart---the average for a human head
			dist.dist(0, x - -0.1, y);
			dist.dist(1, x - +0.1, y);

			// Apply distance cues
			float2 ear = dist(s);

			io.out(0) = ear[0];
			io.out(1) = ear[1];
		}
	}
};

int main(){
	MyApp().start();
}

