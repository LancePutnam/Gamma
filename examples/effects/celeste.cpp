/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example:	Celeste Effect
Author:		Lance Putnam, 2025

Description:
This demonstrates a celeste (from Latin for "heavenly") effect where two
pitch-shifted versions of a signal are summed together to produce rich beating
effects. The effect presented here can operate on ANY signal, i.e., it is more
general than the similar technique of summing two oscillators with slightly 
different frequencies (detuning).
*/
#include "../AudioApp.h"
#include "Gamma/Effects.h"
#include "Gamma/SamplePlayer.h"
using namespace gam;


class MyApp : public AudioApp{
public:

	SamplePlayer<> player;	// Yes, we can apply this to an arbitrary signal!
	Celeste<> celeste;		// Effect unit
	Accum<> tmr{1./4};
	bool effectActive = false;

	MyApp(){
		player.load("../../sounds/count.wav");
		tmr.finish();
	}

	void onAudio(AudioIOData& io){
		while(io()){

			float s = player();
			player.loop();

			if(tmr()){
				effectActive^=true;
				printf("Celeste %s\n", effectActive?"on":"off");
			}

			if(effectActive){
				s = celeste(s)*0.71f; 		// minimal usage emits two detuned signals
				//s = (s+celeste(s))*0.5f;	// slightly fancier including original
			}

			io.out(0) = io.out(1) = s;
		}
	}
};

int main(){
	MyApp().start();
}
