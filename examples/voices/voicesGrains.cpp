/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example:	Voices: Generating Grains
Author:		Lance Putnam, 2021

Description:
This demonstrates how to use the Voices class to implement a granular synth.
Voices is designed to work just fine at the sample level, so do not be afraid
of using it on micro time scales.
*/

#include "../AudioApp.h"
#include "Gamma/rnd.h"
#include "Gamma/Envelope.h"
#include "Gamma/Oscillator.h"
#include "Gamma/Voices.h"
using namespace gam;

class MyApp : public AudioApp{
public:

	struct SineGrain : public Voice<float> {
		Sine<> osc;
		AD<> env{0.01, 0.02};

		void onAttack(float freq, float amp){
			osc.freq(freq);
			env.reset();
			env.amp(amp);
		}

		void onRelease(){
			env.release();
		}

		float operator()(){
			return osc() * env();
		}
	};


	Voices<SineGrain, 32> grains;	// Up to 32 simultaneous grains
	Accum<> tmr{0.1, 0.999};		// Timer for triggering new grains
	float spread = 3;
	bool reverseEnv = false;

	void onAudio(AudioIOData& io){
		while(io()){

			// Trigger a new grain
			if(tmr()){
				float n = rnd::uni(spread); // log-freq "note"
				float f = pow(2, n)*110;
				if(rnd::prob(0.02)) reverseEnv^=true;
				float D = rnd::uni(0.1) + 0.01;
				float A = rnd::uni(0.01) + 0.002;

				float dur = 1./(1. + n*n);
				A *= dur; D *= dur;

				if(reverseEnv) std::swap(A,D);

				auto& grain = grains.attack(f, 0.2);
				grain.env.attack(A).decay(D);

				tmr.period(rnd::uni(0.05) + 0.01);
				if(rnd::prob(0.1)) spread = rnd::uni(6.);
			}

			float s = grains();

			io.out(0) = s;
			io.out(1) = s;
		}
	}
};

int main(){
	MyApp().start();
}

