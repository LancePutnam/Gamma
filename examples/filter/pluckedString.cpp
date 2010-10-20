/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Example:		Filter / Plucked String
	Description:	Simulation of a plucked string with noise and a feedback 
					delay-line.
*/

#include "../examples.h"

struct PluckedString{
	PluckedString(float frq=440) : env(0.1), fil(2), delay(1./27.5, 1./frq){}

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


Accum<> tmr(10);
PluckedString pluck1(scl::freq("d6"));
PluckedString pluck2(scl::freq("g5"));
PluckedString pluck3(scl::freq("a4"));
PluckedString pluck4(scl::freq("d3"));

void audioCB(AudioIOData& io){

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

RUN(audioCB);
