/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Example:		Algorithmic / Luster
	Description:	Thinning the spectrum of a drone-like sound with a STFT
*/

#include "../examples.h"

int rng = 1;
const int dftSize = 2048*2;
STFT stft(dftSize, dftSize/4, 0, WinType::Hann, Bin::Rect, 1);
Accum<> tmr(1./10);
LFO<> osc1(40), osc2(40.003);
LFO<> oscA(62), oscB(62.003);
LFO<> mix(1./60);
FreqShift<> fshift1, fshift2;
LFO<> modfs1(1./101), modfs2(1./102);
Chorus chrA1(0.31, 0.002, 0.20111, -0.7, 0.9);
Chorus chrA2(0.22, 0.002, 0.10151, -0.7, 0.9);
Chorus chrA3(0.13, 0.002, 0.05131, -0.7, 0.9);
Chorus chrB1(0.31, 0.002, 0.20141, -0.7, 0.9);
Chorus chrB2(0.22, 0.002, 0.10171, -0.7, 0.9);
Chorus chrB3(0.13, 0.002, 0.05111, -0.7, 0.9);

void audioCB(AudioIOData& io){
	using namespace gam::rnd;
	while(io()){

		if(tmr()){ if(prob()) rng++; }

		float mx = mix.hann();
		float s = (oscA.up() - oscB.up())*mx + (osc1.up() - osc2.up())*(1-mx);

		fshift1.freq(modfs1.tri()*2);
		fshift2.freq(modfs2.tri()*2);
		s = ((fshift1(s) + fshift2(s))/2 + s)/1;

		if(stft(s)){
		
			rnd::push(rng);
			float prb = rnd::uni(0.3, 0.1);
			for(unsigned k=0; k<stft.numBins(); ++k){
				//float frac = double(k)/stft.numBins();
				float m = pick(pick(2,1, 0.1),0, prb);
				stft.bins(k) *= m;
			}
			rnd::pop();
			stft.zeroEnds();
		}
			
		s = stft()*0.5;
		//float s0=s, s1=s;
		float s0 = chrA3(chrA2(chrA1(s)));
		float s1 = chrB3(chrB2(chrB1(s)));

		io.out(0) = s0;
		io.out(1) = s1;
	}
}

RUN(audioCB);
