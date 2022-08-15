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

	Accum<> tmr{1./2.};		// Timer to switch between LFO types
	LFO<> osc{220};			// Low-frequency oscillator
	LFO<> mod{1./2.};		// Modulator on modifier parameter
	int waveform = 0;

	void onAudio(AudioIOData& io){

		while(io()){

			// Increment waveform type
			if(tmr()) ++waveform;

			// Modulate modifier parameter with raised cosine
			osc.mod(mod.cosCubU());
		
			float s = 0.f;
			
			switch(waveform){ 
					
				// non-modifiable generators ordered from dull to bright:
				case  0: s = osc.sin();		break;		// Sine approximation: one harmonic
				case  1: s = osc.cosCub();	break;		// Cosine approximation from cubic poly
				case  2: s = osc.sinPara();	break;		// Sine approximation from parabolas
				case  3: s = osc.even3();	break;		// Even harmonic sine-like wave (3rd order)
				case  4: s = osc.even5();	break;		// Even harmonic sine-like wave (5th order)
				case  5: s = osc.tri();		break;		// Triangle wave: 1/f^2 odd harmonics
				case  6: s = osc.para();	break;		// Parabola train: 1/f^2 all harmonics
				case  7: s = osc.sqr()/4;	break;		// Square wave: 1/f odd harmonics
				case  8: s = osc.down()/4;	break;		// Downward saw wave: 1/f all harmonics
				case  9: s = osc.up()/4;	break;		// Upward saw wave: 1/f all harmonics
				case 10: s = osc.imp()/4;	break;		// Impulse train: flat spectrum all harmonics
				
				// modifiable generators ordered from dull to bright:
				case 11: s = osc.stair()/4;	break;		// Mix between a square and impulse
				case 12: s = osc.pulse()/4;	break;		// Mix between up and down
				case 13: s = osc.line2()/4;	break;		// Mix between a saw and a triangle
				case 14: s = osc.up2()/4;	break;		// Mix between two saws
				default: waveform = 0;
			}
			
			io.out(0) = io.out(1) = s * 0.2f;
		}
	}
};

int main(){
	MyApp().start();
}
