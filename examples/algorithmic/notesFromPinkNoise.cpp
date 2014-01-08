/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Example:		Algorithmic / Notes From Pink Noise
	Description:	This is an example of using a pink noise generator to
					determine the pitches of notes. Compared to other noise
					types, like white and brown, pink noise supposedly produces 
					the most musical results.
*/

#include "../examples.h"

Sine<> src;
AD<> env;
NoisePink<> pink;
Accum<> tmr;

void audioCB(AudioIOData& io){

	tmr.period(0.25);
	env.attack(0.01).decay(0.24);

	while(io()){

		if(tmr()){
			float frq = pink() * 12;
			
			// Map noise value to nearest pitch in C-major scale
			frq = scl::ratioET(scl::nearest(frq, "2212221")) * scl::freq("c4");

			src.freq(frq);
			env.reset();
		}

		float s = src() * env() * 0.2;

		io.out(0) = s;
		io.out(1) = s;
	}
}

RUN_AUDIO_MAIN
