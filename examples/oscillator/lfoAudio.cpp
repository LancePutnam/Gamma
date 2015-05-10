/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example:	LFO Audio Source
Author:		Lance Putnam, 2012

Description:
This plays through the different waveforms of the LFO used as an audio source.
*/

#include "../AudioApp.h"
#include "Gamma/Oscillator.h"
using namespace gam;

class MyApp : public AudioApp{
public:

	Accum<> tmr;	// Timer to switch between LFO types
	LFO<> osc;		// Low-frequency oscillator
	LFO<> mod;		// Modulator on modifier parameter
	int waveform;

	MyApp(){
		tmr.period(2);
		mod.period(2);
		osc.freq(220);
		waveform=0;
	}

	void onAudio(AudioIOData& io){

		while(io()){

			// Increment waveform type
			if(tmr()) (++waveform)%=13;

			// Modulate modifier parameter with raised cosine
			osc.mod(mod.cosU());
		
			float s = 0.f;
			
			switch(waveform){ 
					
				// non-modifiable generators ordered from dull to bright:
				case  0: s = osc.cos();		break;		// Cosine approximation: one harmonic
				case  1: s = osc.even3();	break;		// Even harmonic sine-like wave (3rd order)
				case  2: s = osc.even5();	break;		// Even harmonic sine-like wave (5th order)
				case  3: s = osc.tri();		break;		// Triangle wave: 1/f^2 odd harmonics
				case  4: s = osc.para();	break;		// Parabola train: 1/f^2 all harmonics
				case  5: s = osc.sqr()/4;	break;		// Square wave: 1/f odd harmonics
				case  6: s = osc.down()/4;	break;		// Downward saw wave: 1/f all harmonics
				case  7: s = osc.up()/4;	break;		// Upward saw wave: 1/f all harmonics
				case  8: s = osc.imp()/4;	break;		// Impulse train: flat spectrum all harmonics
				
				// modifiable generators ordered from dull to bright:
				case  9: s = osc.stair()/4;	break;		// Mix between a square and impulse
				case 10: s = osc.pulse()/4;	break;		// Mix between up and down
				case 11: s = osc.line2()/4;	break;		// Mix between a saw and a triangle
				case 12: s = osc.up2()/4;	break;		// Mix between two saws
			}
			
			io.out(0) = io.out(1) = s * 0.2f;
		}
	}
};

int main(){
	MyApp().start();
}
