/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example:	Biquad Filter
Author:		Lance Putnam, 2012

Description:
This shows how to use a Biquad to perform filtering with various frequency
responses.
*/

#include "../AudioApp.h"
#include "Gamma/Filter.h"
#include "Gamma/Noise.h"
#include "Gamma/Oscillator.h"
using namespace gam;

class MyApp : public AudioApp{
public:
	NoiseWhite<> src;	// Source to filter
	Biquad<> bq;		// Biquad filter
	LFO<> mod;			// Modulator on cutoff frequency
	Accum<> tmr;		// Timer to switch between filter types
	int cnt;			// Counter for filter type

	MyApp(){
		bq.res(4);		// Set resonance of filter
		bq.level(2);	// Set peak level (PEAKING type only)
		mod.period(5);
		tmr.period(5);
		tmr.phaseMax();
		cnt=0;
	}

	void onAudio(AudioIOData& io){

		while(io()){
			if(tmr()){
				switch(cnt){
					case 0: bq.type(LOW_PASS); printf("Low-pass\n"); break;
					case 1: bq.type(HIGH_PASS); printf("High-pass\n"); break;
					case 2: bq.type(BAND_PASS); printf("Band-pass\n"); break;
					case 3: bq.type(BAND_REJECT); printf("Band-reject\n"); break;
					case 4: bq.type(PEAKING); printf("Peaking\n"); break;
				}
				++cnt %= 5;
			}
			
			float cutoff = scl::pow3(mod.triU()) * 10000;

			// Set cutoff frequency of filter
			bq.freq(cutoff);

			// Filter source
			float s = bq(src());
		
			io.out(0) = io.out(1) = s * 0.1;
		}
	}
};

int main(){
	MyApp().start();
}
