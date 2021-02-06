/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example:	Fourier Synthesis
Author:		Lance Putnam, 2012

Description:
Example of building waveforms using Fourier series (sums of harmonics).
*/
#include "../AudioApp.h"
#include "Gamma/Oscillator.h"
using namespace gam;

class MyApp : public AudioApp{
public:

	Accum<> timer{1};
	Osc<> osc{220, 0, 512}; // freq, phase, table size
	static const int K = 16; // max harmonics
	int hNum; // harmonic number
	int wave = 0; // waveform

	MyApp(){
		reset();
	}

	void reset(){
		hNum = 0;
		osc.zero();
	}

	void onAudio(AudioIOData& io){
		while(io()){

			if(timer()){
				++hNum;
				if(hNum == K){
					reset();
					++wave;
				}
				else{
					// addSine(h, amp, phase) adds a sine wave to the table
					// where h is the harmonic number, amp is the amplitude, and
					// phase is the phase of the sine wave (0.25 for a cosine).
					switch(wave){
					case 0:
						if(hNum==1) printf("Saw Wave\n");
						osc.addSine(hNum, 1./hNum);
						break;
					case 1:
						if(hNum==1) printf("Square Wave\n");
						osc.addSine(2*hNum-1, 1./(2*hNum-1));
						break;
					case 2:
						if(hNum==1) printf("Parabolic Wave\n");
						osc.addSine(hNum, 1./(hNum*hNum), 0.25);
						break;
					case 3:
						if(hNum==1) printf("Triangle Wave\n");
						osc.addSine(2*hNum-1, 1./((2*hNum-1)*(2*hNum-1)), 0.25);
						break;
					default:;
					}
				}
			}

			float s = osc() * 0.25;

			io.out(0) = s;
			io.out(1) = s;
		}
	}
};

int main(){
	MyApp().start();
}
