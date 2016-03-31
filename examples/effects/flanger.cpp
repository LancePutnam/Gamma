/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example:	Flanger Effect
Author:		Lance Putnam, 2012

Description:
This demonstrates how to create a flanging effect by slowly modulating the delay
time of a delay line.
*/
#include "../AudioApp.h"
#include "Gamma/Delay.h"
#include "Gamma/Oscillator.h"
using namespace gam;

class Flanger{
public:
	Flanger(float delay=1./500, float modAmount=1./1000, float modFreq=0.5)
	:	delay(delay), modAmount(modAmount),
		comb(1./20, delay, 1, 0), mod(modFreq)
	{}
	
	float operator()(float i){
		comb.delay(delay + mod.cos()*modAmount);
		return comb(i);
	}
	
	void feeds(float ffd, float fbk){
		comb.ffd(ffd); comb.fbk(fbk);
	}

	float delay, modAmount;
	Comb<float, ipl::AllPass> comb;
	LFO<> mod;
};

class MyApp : public AudioApp{
public:

	Accum<> tmr;		// Timer for switching flanging types
	LFO<> src;			// Harmonically rich source
	Flanger flanger;	// Flanger unit
	int flangeType;		// Flanging type

	MyApp(){
		flangeType = 0;
		tmr.period(5);
		tmr.phaseMax();
		src.freq(110);
	}

	void onAudio(AudioIOData& io){
		while(io()){

			if(tmr()){
				switch(flangeType){
				case 0:	printf("Low-pass feedforward\n");
					flanger.feeds(0.7, 0); break;
				case 1:	printf("High-pass feedforward\n");
					flanger.feeds(-0.7,0); break;
				case 2:	printf("Low-pass feedback\n");
					flanger.feeds(0,0.7); break;
				case 3:	printf("High-pass feedback\n");
					flanger.feeds(0,-0.7); break;
				case 4:	printf("Low-pass dual-feed\n");
					flanger.feeds(0.7,0.7); break;
				case 5:	printf("High-pass dual-feed\n");
					flanger.feeds(-0.7,-0.7); break;
				case 6:	printf("All-pass 1\n");
					flanger.feeds(0.9,-0.9); break;
				case 7:	printf("All-pass 2\n");
					flanger.feeds(-0.9,0.9); break;
				}
				(++flangeType) %= 8;
			}
			
			float s = src.up();
			s = flanger(s)*0.1;

			io.out(0) = io.out(1) = s;
		}
	}
};

int main(){
	MyApp().start();
}
