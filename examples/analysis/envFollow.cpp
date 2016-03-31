/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example:	Envelope Follower
Author:		Lance Putnam, 2012

Description:
This shows how to use EnvFollow to estimate the amplitude of a signal.
The example uses an EnvFollow to reset a plucked string whenever its amplitude 
goes below a certain threshold.
*/

#include "../AudioApp.h"
#include "Gamma/Analysis.h"
#include "Gamma/Effects.h"
using namespace gam;

class MyApp : public AudioApp{
public:

	Pluck pluck;
	EnvFollow<> envFollow;

	void onAudio(AudioIOData& io){

		// You can also set the response time (lag) of the envelope follower.
		// Here, increasing the lag will increase the duration between notes.
		//envFollow.lag(2);

		while(io()){

			// Generate source signal
			float s = pluck() * 0.2;

			// Get amplitude estimate
			float ampEst = envFollow(s);

			if(ampEst < 0.001){
				pluck.reset();
				pluck.freq(rnd::uni(1, 20)*100);
			}

			io.out(0) = io.out(1) = s;
		}
	}
};

int main(){
	MyApp().start();
}

