/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Example:		Filter / Biquad
	Description:	Filtering with a multimode biquad filter
*/

#include "../examples.h"

LFO<> mod(0.2, 0.5);
NoiseWhite<> src;
Biquad<> filt(10, 4);

Accum<> tmr(0.2,1);			// Timer to switch between filter types
int cnt=0;					// Counter for filter type

void audioCB(AudioIOData& io){

	while(io()){
		
		if(tmr()){
			switch(cnt){
				case 0: filt.type(Filter::LP); printf("Low-pass\n"); break;
				case 1: filt.type(Filter::HP); printf("High-pass\n"); break;
				case 2: filt.type(Filter::BP); printf("Band-pass\n"); break;
				case 3: filt.type(Filter::BR); printf("Band-reject\n"); break;
			}
			++cnt %= 4;
		}
		
		float cutoff = scl::pow3(mod.triU()) * 10000;

		filt.res(4);
		filt.freq(cutoff);
		
		float s = filt(src());
	
		io.out(0) = io.out(1) = s * 0.1f;
	}
}

RUN(audioCB);
