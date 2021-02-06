/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example:	Plucked String
Author:		Lance Putnam, 2012

Description:
Simulation of a plucked string using noise and a feedback delay-line.
*/

#include "../AudioApp.h"
#include "Gamma/Delay.h"
#include "Gamma/Envelope.h"
#include "Gamma/Filter.h"
#include "Gamma/Noise.h"
#include "Gamma/Oscillator.h"
using namespace gam;


class PluckedString{
public:
	PluckedString(float frq=440)
	:	env(0.1), fil(2), delay(1./27.5, 1./frq){}

	float operator()(){
		return (*this)(noise()*env());
	}

	float operator()(float in){
		return delay(
			fil( delay() + in )
		);
	}

	void reset(){ env.reset(); }
	void freq(float v){ delay.freq(v); }
	
	NoiseWhite<> noise;
	Decay<> env;
	MovingAvg<> fil;
	Delay<float, ipl::Trunc> delay;
};


class MyApp : public AudioApp{
public:

	Accum<> tmr{1./0.1};
	PluckedString
		pluck1{scl::freq("d6")},
		pluck2{scl::freq("g5")},
		pluck3{scl::freq("a4")},
		pluck4{scl::freq("d3")};

	void onAudio(AudioIOData& io){

		while(io()){

			if(tmr()){
				if(rnd::prob(0.1/1)) pluck1.reset();
				if(rnd::prob(0.1/2)) pluck2.reset();
				if(rnd::prob(0.1/3)) pluck3.reset();
				if(rnd::prob(0.1/4)) pluck4.reset();
			}
			
			float s = pluck1() + pluck2() + pluck3() + pluck4();
				
			io.out(0) = io.out(1) = s * 0.2f;
		}
	}
};

int main(){
	MyApp().start();
}

