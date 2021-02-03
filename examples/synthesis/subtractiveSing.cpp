/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example:	Singing Voice
Author:		Lance Putnam, 2012

Description:
Demonstration of singing voice synthesis using subtractive synthesis. We start
with a rich harmonic source and then pass it through a bank or parallel
bandpass filters to simulate formants. The formant center frequencies are
changed to produce different vowel sounds.
*/

#include "../AudioApp.h"
#include "Gamma/rnd.h"
#include "Gamma/Delay.h"
#include "Gamma/Filter.h"
#include "Gamma/FormantData.h"
#include "Gamma/Oscillator.h"
using namespace gam;

struct VowelFilter{

	VowelFilter(){ params.lag(1); }

	float operator()(float v){
		Vec<N*2, float> frqAmp = params();
		float r=0;
		for(int i=0; i<N; ++i){
			filters[i].freq(frqAmp[2*i+0]);
			r += filters[i](v)*frqAmp[2*i+1];
		}
		return r;
	}
	
	void set(int i, float freq, float amp){
		filters[i].type(RESONANT);
		filters[i].res(1./0.08);
		params.stored()[i*2+0] = freq;
		params.stored()[i*2+1] = amp;
	}
	
	int size() const { return N; }
	
	enum{N=3};
	Biquad<> filters[N];
	OnePole<Vec<N*2, float> > params; // array of (cfreq, amp) pairs
};

class MyApp : public AudioApp{
public:

	Accum<> tmr;
	unsigned ivowel = 0;
	VowelFilter filt;
	OnePole<> pitch;
	DWO<> osc;
	Sine<> vib{5};
	Delay<> echo{0.27};

	MyApp(){
		tmr.period(2);
		tmr.phaseMax();
		pitch.lag(1);
	}

	void onAudio(AudioIOData& io){
		while(io()){
			if(tmr()){
				// Set new vowel
				for(int i=0; i<filt.size(); ++i){
					filt.set(i,
						Vowel::freq(Vowel::WOMAN, Vowel::Phoneme(ivowel), i),
						Vowel::amp(Vowel::WOMAN, Vowel::Phoneme(ivowel), i)
					);
				}

				// Set new pitch
				if(0 == (ivowel%4)){
					pitch = pow(2, rnd::uni(7)/12.) * 220;
				}

				(++ivowel) %= Vowel::NUM_PHONEMES;
			}

			osc.freq(pitch() * (1 + vib()*0.01));

			// Generate pulse wave
			osc.mod(0.1);
			float s = osc.pulse()*0.1;

			// Apply formant filter
			s = filt(s);

			// Add some echo
			s += echo(s + echo()*0.7)*0.3;

			io.out(0) = io.out(1) = s;
		}
	}
};

int main(){
	MyApp().start();
}
