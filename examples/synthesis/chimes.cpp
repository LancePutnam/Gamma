/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example:	Chimes
Author:		Lance Putnam, 2014

Description:
Wind chimes created from an additive model of a struck bar.
*/

#include "../AudioApp.h"
#include "Gamma/Effects.h"
#include "Gamma/Oscillator.h"
using namespace gam;

// Returns frequency ratio of a mode of a bar clamped at one end
float barClamp(float mode){
	float res = mode - 0.5;
	return 2.81f*res*res;
}

// Returns frequency ratio of a mode of a freely vibrating bar
float barFree(float mode){
	float res = mode + 0.5;
	return 0.441f*res*res;
}

class MyApp : public AudioApp{
public:

	static const int Nc = 9; // # of chimes
	static const int Nm = 5; // # of modes
	SineDs<> src{Nc * Nm};
	Accum<> tmr;
	Chorus<> chr1{0.10}, chr2{0.11}; // chorusing for more natural beating

	MyApp(){
		tmr.finish();
	}

	void onAudio(AudioIOData& io){
		while(io()){
			if(tmr()){
				static double freqs[Nc] = {
					scl::freq("c4"),
					scl::freq("f4"), scl::freq("g4"), scl::freq("a4"), scl::freq("d5"),
					scl::freq("f5"), scl::freq("g5"), scl::freq("a5"), scl::freq("d6"),
				};

				int i = rnd::uni(Nc);
				float f0 = freqs[i];
				float A = rnd::uni(0.1,1.);

				//      osc #   frequency      amplitude               length
				src.set(i*Nm+0, f0*1.000,      A*1.0*rnd::uni(0.8,1.), 1600.0/f0);
				src.set(i*Nm+1, f0*barFree(2), A*0.5*rnd::uni(0.5,1.), 1200.0/f0);
				src.set(i*Nm+2, f0*barFree(3), A*0.4*rnd::uni(0.5,1.),  800.0/f0);
				src.set(i*Nm+3, f0*barFree(4), A*0.3*rnd::uni(0.5,1.),  400.0/f0);
				src.set(i*Nm+4, f0*barFree(5), A*0.2*rnd::uni(0.5,1.),  200.0/f0);

				tmr.period(rnd::uni(0.01,1.5));
			}

			float2 s = src() * 0.1;

			s.x += chr1(s.x)*0.2;
			s.y += chr2(s.y)*0.2;

			io.out(0) = s.x;
			io.out(1) = s.y;
		}
	}
};

int main(){
	MyApp().start();
}
