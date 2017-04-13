/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example:	Hard/soft sync
Author:		Lance Putnam, 2017

Description:
This demonstrates how to generate combed spectra by syncing the phase of a 
"formant" oscillator to a "sync" oscillator. The spectrum consists of peaks 
centered on the partials of the formant oscillator and has a fundamental 
frequency matching that of the sync oscillator. A "hard" sync simply resets the 
phase of the formant oscillator whenever the sync oscillator completes one 
cycle. This creates hard edges in the resulting waveform leading to harsh 
distortion and large changes in brightness when the formant oscillator frequency 
changes. A "soft" sync works just like a hard sync, except the resulting 
waveform is multiplied by a periodic window function that drops to zero at the 
sync points. This removes the hard edges and produces a more even sound as the 
formant oscillator frequency is changed.
*/
#include "../AudioApp.h"
#include "Gamma/Oscillator.h"
using namespace gam;

class MyApp : public AudioApp{
public:

	LFO<> sync{110};	// The sync osc determines the fundamental frequency
	LFO<> formant;		// The formant osc determines the formant locations
	LFO<> mod{1./8};	// Used to sweep the formant frequencies
	bool softSync = true;
	Accum<> syncMode{1./8, 0.999};

	void onAudio(AudioIOData& io){
		while(io()){
			if(syncMode()){
				softSync^=true;
				printf("%s sync\n", softSync?"soft":"hard");
			}

			// When the sync osc completes cycle, reset the phase of the formant osc
			float w = sync.up();
			if(sync.cycled()) formant.phase(0);

			// Set the pitch of the formants
			formant.freq(mod.triU()*2200+200);

			// The formant oscillator determines the spectral envelope
			float s = formant.tri();

			// Drop amp to zero at sync points?
			if(softSync){
				w = 1.f - scl::pow8(w);	// flattop window
				s *= w;					// apply window
			}

			io.out(0) = io.out(1) = s*0.2;
		}
	}
};

int main(){
	MyApp().start();
}
