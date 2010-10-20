/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Example:		Algorithmic / Saw Cluster
	Description:	A set of harmonically tuned saw waves are fed into slowly
					modulated feedback delay-lines and filters.
*/

#include "../examples.h"

// Sounds a bit like...
Accum<> tmr(0.5, 2);
Seg<> env0(2);
LFO<> osc0(440), osc1(330.04), osc2(180.04), osc3(90.03);
Sine<> mod0(0.63,0.5), mod1(0.71,0.75), mod2(0.87,1), mod3(1,1);
Comb<> del0(1, -0., 0.9), del1(1, -0., 0.9), ech0(0.3,0,0.7), ech1(0.32,0,0.7);
OnePole<> lag(400), frq0(0.1), frq1(0.1);
AllPass1<> ap0(2000), ap1(2000);
Reson<> res0(400, 2000), res1(800, 2000);

void audioCB(AudioIOData& io){

	while(io()){
		using namespace gam::rnd;

		if(tmr()){
			env0 = uni(0.4, 0.39);
			if(prob(0.8)){
				float r = uni(1.);
				tmr.period(r * 4);
				env0.period(r * 4);
			}
			
			int a = pick(8,6, 0.7);
			if(prob(0.2)) osc0.freq(quanOct(a, 440.));
			if(prob(0.1)) osc1.freq(quanOct(a, 220.));
			if(prob(0.1)) osc2.freq(quanOct(a, 110.));
			if(prob(0.1)) osc3.freq(quanOct(a,  55.));
			if(prob(0.2)) frq0 = lin(8000, 400);
			if(prob(0.2)) frq1 = lin(8000, 400);//printf("d");
		}
		
		float e = lag(env0());
		del0.delay(e);
		del1.delay(e * 0.9);
		
		float s = (osc0.up() * mod0() + osc1.up() * mod1() + osc2.up() * mod2() + osc3.up() * mod3()) * 0.05;
		res0.freq(frq0()); res1.freq(frq1());
		s = res0(s) + res1(s);

		float sl = ech0(del0(s), ap0(ech0()));
		float sr = ech1(del1(s), ap1(ech1()));

		io.out(0) = sl;
		io.out(1) = sr;
	}
	
}

RUN(audioCB)
