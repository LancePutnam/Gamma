/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example:	Feedback Phase Modulation
Author:		Lance Putnam, 2012

Description:
Single-oscillator feedback phase modulation.
*/

#include "../AudioApp.h"
#include "Gamma/Oscillator.h"
using namespace gam;

class MyApp : public AudioApp{
public:

	Sine<> osc;		// Source sine
	LFO<> mod;		// Modulator on feedback amount
	float prev;		// Previous sample

	MyApp(){
		osc.freq(220);
		mod.period(10);
		prev=0;
	}

	void onAudio(AudioIOData& io){
		while(io()){
			float fbk = prev*mod.hann()*0.2;

			// Add feedback to phase
			osc.phaseAdd(fbk);

			prev = osc();

			// Subtract feedback from phase to avoid changing pitch
			osc.phaseAdd(-fbk);

			io.out(0) = io.out(1) = prev * 0.2f;
		}
	}
};

int main(){
	MyApp().start();
}
