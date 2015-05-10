/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example:	LFO Modulation
Author:		Lance Putnam, 2012

Description:
This demonstrates use of an LFO as an amplitude and frequency modulator.
*/

#include "../AudioApp.h"
#include "Gamma/Noise.h"
#include "Gamma/Oscillator.h"
using namespace gam;

class MyApp : public AudioApp{
public:

	Accum<> tmr;			// Timer to switch between LFO types
	NoiseWhite<> noise;		// Carrier used for amplitude modulation
	LFO<> osc;				// Carrier used for frequency modulation
	LFO<> lfo;				// Modulator on amplitude or frequency
	LFO<> mod;				// Modulator on modifier parameter
	gen::Trigger cnt;		// Counter for changing LFO type
	bool modMode;			// true = amp mod, false = freq mod

	MyApp()
	:	cnt(10,10)
	{
		tmr.period(4);
		tmr.phaseMax();
		mod.period(2);
		lfo.freq(5);
		modMode = false;
	}

	void onAudio(AudioIOData& io){

		while(io()){

			if(tmr()){
				if(cnt()){
					modMode ^= true;			// increment LFO type, switch mod mode when wraps
					printf("\nMod mode: %s modulation\n", modMode? "Amp" : "Freq");
				}
			}
				
			lfo.mod(mod.cosU());	// modulate modifier parameter with unipolar cosine wave

			float s = 0.f;			// initialize current sample to zero

			switch(cnt.val){			
					
				// non-modifiable generators ordered from smooth to sharp:
				case  0: s = lfo.cosU();	break;	// unipolar cosine
				case  1: s = lfo.hann();	break;	// a computed hann window (inverted cosU)
				case  2: s = lfo.triU();	break;	// unipolar triangle
				case  3: s = lfo.upU();		break;  // unipolar upward ramp
				case  4: s = lfo.downU();	break;  // unipolar downward ramp
				case  5: s = lfo.sqrU();	break;	// unipolar square

				// modifiable generator ordered from smooth to sharp:
				case  6: s = lfo.pulseU();	break;	// Mix between upward ramp and downward ramp
				case  7: s = lfo.stairU();	break;	// Mix between a square and impulse.
				case  8: s = lfo.line2U();	break;	// Mix between a ramp and a triangle
				case  9: s = lfo.up2U();	break;  // Mix between two ramps
			}
			
			if(modMode){			// amplitude modulate noise with envelope
				s *= noise();
			}
			else{					   // frequency modulate oscillator with envelope
				osc.freq(s*400 + 200); // between 100 and 200 hz
				s = osc.cos();
			}
			
			io.out(0) = io.out(1) = s*0.2;
		}
	}
};

int main(){
	MyApp().start();
}
