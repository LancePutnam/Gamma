/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example:	Comb Filter
Author:		Lance Putnam, 2012

Description:
This demonstrates various frequency responses of a comb filter.
*/

#include "../AudioApp.h"
#include "Gamma/Delay.h"
#include "Gamma/Oscillator.h"
using namespace gam;

class MyApp : public AudioApp{
public:

	int feedType;
	Accum<> tmr;	// Switch between flanging types
	Saw<> src;
	//NoiseWhite<> src;
	LFO<> mod;
	Comb<float, ipl::Switchable> comb;

	MyApp(){
		feedType=0;
		tmr.phaseMax();
		tmr.period(4);
		mod.period(4);
		src.freq(100);
		comb.maxDelay(1./100);
	}

	void onAudio(AudioIOData& io){

		while(io()){

			if(tmr()){
				switch(feedType){
				case 0:	printf("Low-pass feedforward\n");
					comb.feeds( 1,0); break;
				case 1:	printf("High-pass feedforward\n");
					comb.feeds(-1,0); break;
				case 2:	printf("Low-pass feedback\n");
					comb.feeds(0,0.7); break;
				case 3:	printf("High-pass feedback\n");
					comb.feeds(0,-0.7); break;
				case 4:	printf("Low-pass dual-feed\n");
					comb.feeds(0.7,0.7); break;
				case 5:	printf("High-pass dual-feed\n");
					comb.feeds(-0.7,-0.7); break;
				case 6:	printf("All-pass 1\n");
					comb.feeds(0.9,-0.9); break;
				case 7:	printf("All-pass 2\n");
					comb.feeds(-0.9,0.9); break;
				}
				(++feedType) %= 8;
			}
			
			float s = src()*0.1;
			
			comb.ipolType(ipl::ROUND);
			comb.delay(mod.triU() * 1./400 + 1./10000);
			s = comb(s);

			io.out(0) = io.out(1) = s;
		}
	}
};

int main(){
	MyApp().start();
}
