/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Example:		Algorithmic / Singer
	Description:	Demonstration of singing voice synthesis

*/

#include "../examples.h"

struct VowelFilter{

	VowelFilter(): params(1){}

	float operator()(float v){
		params();
		float r=0;
		for(int i=0; i<N; ++i){
			filters[i].freq(params.last()[2*i+0]);
			r += filters[i](v)*params.last()[2*i+1];
		}
		return r;
	}
	
	void set(int i, float freq, float amp){
		filters[i].type(Filter::BP);
		filters[i].res(1./0.02);
		params.stored()[i*2+0] = freq;
		params.stored()[i*2+1] = amp;
	}
	
	int size() const { return N; }
	
	enum{N=3};
	Biquad<> filters[N];
	OnePole<Vec<N*2, float> > params; // array of (cfreq, amp) pairs
	float amps[N];
};

Accum<> tmr(1./2, 2);
ValWrap<int> ivowel(Vowel::NUM_PHONEMES);
VowelFilter filt;
LFO<> osc;
Sine<> vib(5);

void audioCB(AudioIOData& io){

	while(io()){

		if(tmr()){
			for(int i=0; i<filt.size(); ++i){
				filt.set(i,
					Vowel::freq(Vowel::WOMAN, (Vowel::Phoneme)ivowel(), i),
					Vowel::amp(Vowel::WOMAN, (Vowel::Phoneme)ivowel(), i)
				);
			}
			++ivowel;
		}

		osc.freq(330.5 + vib()*3);
		osc.mod(0.2);

		float s = osc.pulse();
		
		s = filt(s);

		io.out(0) = io.out(1) = s * 0.1f;
	}
}

RUN(audioCB);
