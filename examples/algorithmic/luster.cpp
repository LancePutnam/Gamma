/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example:	Luster
Author:		Lance Putnam, 2012

Description:
Thinning the spectrum of a drone-like sound with a STFT
*/
#include "../AudioApp.h"
#include "Gamma/DFT.h"
#include "Gamma/Effects.h"
#include "Gamma/Oscillator.h"
using namespace gam;

class MyApp : public AudioApp{
public:

	int seed;
	STFT stft;
	Accum<> tmr;
	LFO<> osc1, osc2;
	LFO<> oscA, oscB;
	LFO<> mix;
	FreqShift<> fshift1, fshift2;
	LFO<> modfs1, modfs2;
	Chorus<> chrA1, chrA2, chrA3, chrB1, chrB2, chrB3;

	MyApp()
	:	stft(4096, 4096/4, 0, HANN, COMPLEX),
		chrA1(0.31, 0.002, 0.20111, -0.7, 0.9),
		chrA2(0.22, 0.002, 0.10151, -0.7, 0.9),
		chrA3(0.13, 0.002, 0.05131, -0.7, 0.9),
		chrB1(0.31, 0.002, 0.20141, -0.7, 0.9),
		chrB2(0.22, 0.002, 0.10171, -0.7, 0.9),
		chrB3(0.13, 0.002, 0.05111, -0.7, 0.9)
	{
		tmr.period(10);
		osc1.freq(40);
		osc2.freq(40.003);
		oscA.freq(62);
		oscB.freq(62.003);
		mix.period(60);
		modfs1.period(101);
		modfs2.period(102);
	}

	void onAudio(AudioIOData& io){
		while(io()){

			// Generate new seed?
			if(tmr()){ if(rnd::prob()) seed++; }

			// Modulated mix between pulse waves
			float mx = mix.hann();
			float s = (oscA.up() - oscB.up())*mx + (osc1.up() - osc2.up())*(1-mx);

			// Add frequency shifted versions of signal to get barberpole combs
			fshift1.freq(modfs1.tri()*2);
			fshift2.freq(modfs2.tri()*2);
			s += (fshift1(s) + fshift2(s))*0.5;

			if(stft(s)){
				// Apply spectral thin
				rnd::push(seed);
				float prb = rnd::uni(0.3, 0.1);
				for(unsigned k=0; k<stft.numBins(); ++k){
					//float frac = double(k)/stft.numBins();
					float m = rnd::pick(rnd::pick(2,1, 0.1),0, prb);
					stft.bin(k) *= m;
				}
				rnd::pop();
				stft.zeroEnds();
			}

			s = stft()*0.5;

			// "Spatialize" with modulated echoes
			float s0 = chrA3(chrA2(chrA1(s)));
			float s1 = chrB3(chrB2(chrB1(s)));

			io.out(0) = s0;
			io.out(1) = s1;
		}
	}
};

int main(){
	MyApp().start();
}
