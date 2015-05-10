/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example:	Bassline
Author:		Lance Putnam, 2012

Description:

*/
#include "../AudioApp.h"
#include "Gamma/Envelope.h"
#include "Gamma/Filter.h"
#include "Gamma/Oscillator.h"
using namespace gam;

class MyApp : public AudioApp{
public:

	Saw<> saw;			// Saw oscillator
	Biquad<> lpf;		// Resonant low-pass filter
	AD<> env;			// Envelope on amplitude and cutoff
	LFO<> modCutoff;	// Modulator on cutoff
	OnePole<> freq;		// Portamento filter
	Accum<> tmr;		// Sequencer timer
	int step;			// Sequencer step

	MyApp(){
		lpf.type(LOW_PASS);
		lpf.res(4);
		env.attack(0.01);
		env.decay(0.4);
		tmr.freq(120./60.*4.);
		tmr.phaseMax();
		modCutoff.period(40);
		modCutoff.phase(0.5);
		freq.lag(0.1);
		step=0;
	}

	void onAudio(AudioIOData& io){
		while(io()){

			if(tmr()){
				float pitches[] = {0,0,12,0,0,10,-5,0};
				float f = 55 * pow(2, pitches[step%8]/12.);
				++step;
				freq = f;
				env.resetSoft();
			}

			saw.freq(freq());

			float e = env();
			lpf.freq(e * (modCutoff.triU()*4000 + 1000) + 40);
			float s = saw() * 0.3;
			s = lpf(s) * e;
			io.out(0) = io.out(1) = s;
		}
	}
};

int main(){
	MyApp().start();
}
