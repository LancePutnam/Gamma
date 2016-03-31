/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example:	Chirp Tank
Author:		Lance Putnam, 2012

Description:
A single decaying chirp is fed into multiple feedback delay-lines. The delay 
times and other parameters are randomly modulated over time to avoid
steady-state patterns and to allow interesting developments of the sound.
*/
#include "../AudioApp.h"
#include "Gamma/Effects.h"
using namespace gam;

// Four comb filters connected in series
template <class T = gam::real>
struct Combs4{
	Combs4(float d1, float d2, float d3, float d4, float ffd, float fbk)
	: c1(d1, ffd, fbk), c2(d2, ffd, fbk), c3(d3, ffd, fbk), c4(d4, ffd, fbk){}

	T operator()(T v){ return c4(c3(c2(c1(v)))); }
	
	Comb<T, ipl::Linear> c1, c2, c3, c4;
};

class MyApp : public AudioApp{
public:

	Accum<> tmr;				// Periodic timer for parameter randomization
	Chirp<> src;				// Sine chirp
	LFO<> lfoD;					// Chirp decay time modulator
	Comb<> cmb1, cmb2, cmb3;	// Parallel combs
	Combs4<> cmbs4;				// Series combs
	Chorus<> chr;				// Chorus effect
	OnePole<> opC1, opC2, opC3; // Comb delay interpolation curves
	OnePole<> opSpread;			// Stereo spread interpolation curve
	float dcyMul;				

	MyApp()
	:	src(100., 1, 0.1), chr(0.5, 0.3, 0.01, 1, 0.7),
		cmb1(1, -0.9, 0.9), cmb2(1, -0.9, 0.9), cmb3(1, -0.9, 0.9),
		cmbs4(0.11, 0.17, 0.23, 0.31, -0.9,0.9)
	{
		tmr.freq(16);
		tmr.phaseMax();
		lfoD.freq(2.51);
		opC1.freq(7);
		opC2.freq(7);
		opC3.freq(7);
		opSpread.lag(15);
		dcyMul = 0.1f;
	}

	void onAudio(AudioIOData& io){
		while(io()){
			using namespace gam::rnd;

			float dcy60 = lfoD.upU();
			
			if(tmr()){
				// Randomize chirp parameters		
				if(prob(0.1  )){
					float f = scl::round(pow3(12000.f, 150.f), 300.f);
					src.freq(f, 0);
				}
				if(prob(0.001)) dcyMul = lin(0.4f, 0.01f);
				if(prob(0.3  )){ src.reset(); src.decay(dcy60 * dcyMul); }
			
				// Randomize comb delay times
				if(prob(0.008)){
					float m = 1./add2I(40., 1.);
					opC1 = m / uni(2., 1.);
					opC2 = m / uni(4., 2.);
					opC3 = m / uni(8., 4.);
				}

				// Randomize chorusing rate
				chr.mod.freq(pick(0.001, lin(0.01, 0.1), 0.98));

				// Randomize stereo spread
				if(prob(0.01)) opSpread = add2I(1.f);

				// Randomize timer frequency
				if(prob(0.005)) tmr.freq(add2I(32., 4.));
			}

			// Synthesize chirp
			float t0 = src() * 0.1;

			// Set parallel comb delays
			cmb1.delay(opC1());
			cmb2.delay(opC2());
			cmb3.delay(opC3());

			// Pass chirp through comb filter network
			t0 = cmb1(t0) + cmb2(t0) + cmb3(t0);
			t0 += cmbs4(t0)*0.1;

			// Apply chorus effect
			float t1;
			chr(t0, t0, t1);

			// Apply stereo spread
			scl::mix2(t0, t1, opSpread());
			
			io.out(0) = t0;
			io.out(1) = t1;
		}
	}
};

int main(){
	MyApp().start();
}
