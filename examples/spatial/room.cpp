/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example:	Room Spatialization
Author:		Lance Putnam, 2015

Description:
This demonstrates how to combine distance cues with late reverberation to
simulate a sound moving in a room. Two distance cues model the direct sound and
a third models the "walls" of the room.
The late reverb is delayed and low-pass filtered according to the size of the 
room. Assuming the sound and listener are both in the center of the room, then
the minimum distance the echoes must travel is twice the distance to the nearest
wall.
*/

#include "../AudioApp.h"
#include "Gamma/Oscillator.h"
#include "Gamma/SamplePlayer.h"
#include "Gamma/Spatial.h"
using namespace gam;

class MyApp : public AudioApp{
public:

	SamplePlayer<> play;	// Sound source
	LFO<> pathx, pathy;		// Path of sound source in space
	Dist<3> dist;			// Filters source at 3 destinations based on distance cues
	ReverbMS<> reverb[2]; 	// One reverb for each ear

	MyApp(){
		play.load("../../sounds/count.wav");
		pathx.period(30.13 / 3);
		pathy.period(30.13 / 5);

		for(int i=0; i<2; ++i){
			unsigned d = i*2;
			// The delay lengths will differ by 2 samples for each ear
			reverb[i].resize(gam::JCREVERB, d);

			// Set T60 decay time
			reverb[i].decay(4);

			// High-frequency damping due to air absorption
			reverb[i].damping(0.25);		
		}
	}

	void onAudio(AudioIOData& io){

		// The second parameter is the diameter of the room
		dist.dist(2, 32);

		while(io()){

			// Generate sound source sample
			float s = play(); play.loop();

			// Position of source w.r.t. origin
			float x = pathx.cos()*16;
			float y = pathy.cos()*16;

			// Set distances from source to destinations (ears)
			// The ears are set to be 20 cm apart---the average for a human head
			dist.dist(0, x - -0.1, y);
			dist.dist(1, x - +0.1, y);

			// Get source filtered at three destinations (2 ears + wall)
			float3 spat = dist(s);

			// The walls are not perfect reflectors, so pre-attenuate
			spat[2] *= 0.2;

			// Apply reverb to the wall signal
			float2 echoes(
				reverb[0](spat[2]),
				reverb[1](spat[2])
			);

			// At each ear, we hear the direct sound and room echoes
			float2 ear = spat.get(0,1) + echoes;

			io.out(0) = ear[0];
			io.out(1) = ear[1];
		}
	}
};

int main(){
	MyApp().start();
}

