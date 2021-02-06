/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example:	Bassline Pattern
Author:		Lance Putnam, 2012

Description:
This demonstrates how to create a bassline pattern with elementary unit
generators. The basic technique is pass a saw wave through a resonant low-pass
filter. A single attack-decay envelope is mapped to amplitude and cutoff 
frequency to produce a "pluck" with a natural-sounding decrease in
brightness over time. The filter cutoff is also modulated slowly over time for
interest and to demonstrate the range of different timbres possible.
*/
#include "../AudioApp.h"
#include "Gamma/Envelope.h"
#include "Gamma/Filter.h"
#include "Gamma/Oscillator.h"
using namespace gam;

class MyApp : public AudioApp{
public:

	DWO<> osc;			// Saw oscillator
	Biquad<> lpf;		// Resonant low-pass filter
	AD<> env;			// Envelope on amplitude and cutoff
	LFO<> modCutoff;	// Modulator on cutoff
	OnePole<> freq;		// Portamento filter
	Accum<> tmr;		// Sequencer timer
	int step = 0;		// Sequencer step

	MyApp(){
		lpf.type(LOW_PASS);		// Set filter to low-pass response
		lpf.res(4);				// Set resonance amount to emphasize filter
		env.attack(0.01);		// Set short (10 ms) attack
		env.decay(0.4);			// Set longer (400 ms) decay
		tmr.freq(120./60.*4.);	// Set timer frequency to 120 BPM
		tmr.phaseMax();			// Ensures timer triggers on first sample
		modCutoff.period(30);	// Set period of cutoff modulation
		modCutoff.phase(0.5);	// Start half-way through cycle
		freq.lag(0.1);			// Lag time of portamento effect
	}

	void onAudio(AudioIOData& io){
		while(io()){

			if(tmr()){
				// Our sequence of pitches
				float pitches[] = {0,0,12,0,0,10,-5,0};
				// Map pitch class to a frequency in Hz
				float f = 55 * pow(2, pitches[step]/12.);
				// Increment step counter
				step = (step + 1) % 8;
				// Set new target frequence of portamento
				freq = f;
				// Restart envelope using a soft reset (to avoid clicks)
				env.resetSoft();
			}

			// Set oscillator frequency from portamento filter
			osc.freq(freq());

			// Get next envelope value
			float e = env();
			// Map envelope value to cutoff frequency
			lpf.freq(e * (modCutoff.paraU()*6000 + 500) + 40);
			// Generate next saw sample
			float s = osc.saw() * 0.2;
			// Filter saw sample
			s = lpf(s) * e;
			// Send sample to DAC
			io.out(0) = io.out(1) = s;
		}
	}
};

int main(){
	MyApp().start();
}
