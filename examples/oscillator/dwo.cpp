/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example:	Differenced wave oscillator (DWO)
Author:		Lance Putnam, 2017

Description:
This example compares alias-reduced waveforms using differencing with their 
naive aliasing counterparts. The aliasing oscillator (LFO) is played first and
then the alias-reduced oscillator (DWO) for comparison.
*/

#include "../AudioApp.h"
#include "Gamma/Oscillator.h"
using namespace gam;

class MyApp : public AudioApp{
public:

	Accum<> tmr{1./4};	// Timer to switch between LFO types
	DWO<> dwo;			// Alias-reduced oscillator
	LFO<> lfo;			// "Naive" aliasing oscillator
	int type=0;

	void onAudio(AudioIOData& io){

		lfo.mod(0.2);
		dwo.mod(0.2);

		while(io()){

			// Increment waveform type
			if(tmr()) (++type)%=10;

			float freq = tmr.phase() * 2000 + 440;

			lfo.freq(freq);
			dwo.freq(freq);
		
			float s = 0;

			auto printType = [&](const char * msg){ if(tmr.cycled()) printf("%s\n",msg); };

			switch(type){
				case 0: s = lfo.down();     printType("LFO saw wave"); break;
				case 1: s = dwo.down();     printType("DWO saw wave"); break;
				case 2: s = lfo.sqr();      printType("LFO square wave"); break;
				case 3: s = dwo.sqr();      printType("DWO square wave"); break;
				case 4: s = lfo.pulse();    printType("LFO pulse wave"); break;
				case 5: s = dwo.pulse();    printType("DWO pulse wave"); break;
				case 6: s = lfo.tri()*1.5;  printType("LFO triangle wave"); break;
				case 7: s = dwo.tri()*1.5;  printType("DWO triangle wave"); break;
				case 8: s = lfo.para()*1.5; printType("LFO parabolic wave"); break;
				case 9: s = dwo.para()*1.5; printType("DWO parabolic wave"); break;
			}
			
			io.out(0) = io.out(1) = s*0.2;
		}
	}
};

int main(){
	MyApp().start();
}
