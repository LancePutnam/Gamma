/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example:	Saw Cluster
Author:		Lance Putnam, 2012

Description:
A set of harmonically tuned saw waves are fed into slowly modulated feedback
delay-lines and filters.
*/

#include "../AudioApp.h"
#include "Gamma/rnd.h"
#include "Gamma/Delay.h"
#include "Gamma/Envelope.h"
#include "Gamma/Filter.h"
#include "Gamma/Oscillator.h"
using namespace gam;

class MyApp : public AudioApp{
public:

	Accum<> tmr;
	Seg<> env0;
	LFO<> osc0, osc1, osc2, osc3;
	//Sine<> mod0(0.63,0.5), mod1(0.71,0.75), mod2(0.87,1), mod3(1,1);
	Sine<> mod0, mod1, mod2, mod3;
	Comb<> del0, del1, ech0, ech1;
	OnePole<> lag, frq0, frq1;
	AllPass1<> ap0, ap1;
	Reson<> res0, res1;

	MyApp()
	:	del0(1, -0., 0.9), del1(1, -0., 0.9), ech0(0.3,0,0.7), ech1(0.32,0,0.7)
	{
		tmr.period(2); tmr.phaseMax();
		env0.period(2);
		osc0.freq(440);
		osc1.freq(330.04);
		osc2.freq(180.04);
		osc3.freq(90.03);
		mod0.freq(0.63);
		mod1.freq(0.71);
		mod2.freq(0.87);
		mod3.freq(1);
		lag.lag(1./400);
		frq0.lag(10);
		frq1.lag(10);
		ap0.freq(2000);
		ap1.freq(2000);
		res0.freq(400); res0.width(2000);
		res1.freq(800); res1.width(2000);
	}

	void onAudio(AudioIOData& io){
		while(io()){

			if(tmr()){
				env0 = rnd::uni(0.4, 0.39);
				if(rnd::prob(0.8)){
					float r = rnd::uni(1.);
					tmr.period(r * 4);
					env0.period(r * 4);
				}
				
				int a = rnd::pick(8,6, 0.7);
				if(rnd::prob(0.2)) osc0.freq(rnd::quanOct(a, 440.));
				if(rnd::prob(0.1)) osc1.freq(rnd::quanOct(a, 220.));
				if(rnd::prob(0.1)) osc2.freq(rnd::quanOct(a, 110.));
				if(rnd::prob(0.1)) osc3.freq(rnd::quanOct(a,  55.));
				if(rnd::prob(0.2)) frq0 = rnd::lin(8000, 400);
				if(rnd::prob(0.2)) frq1 = rnd::lin(8000, 400);//printf("d");
			}
			
			float e = lag(env0());
			del0.delay(e);
			del1.delay(e * 0.9);
			
			float s	= osc0.up() * mod0()
					+ osc1.up() * mod1()
					+ osc2.up() * mod2()
					+ osc3.up() * mod3();
			s *= 0.05;
			res0.freq(frq0()); res1.freq(frq1());
			s = res0(s) + res1(s);

			float sl = ech0(del0(s), ap0(ech0()));
			float sr = ech1(del1(s), ap1(ech1()));

			io.out(0) = sl;
			io.out(1) = sr;
		}
	}
};

int main(){
	MyApp().start();
}
