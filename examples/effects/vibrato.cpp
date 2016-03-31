/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example:	Vibrato Effect
Author:		Lance Putnam, 2012

Description:
This demonstrates how to create a vibrato effect by slowly modulating the delay
time of a delay line. Since the vibrato uses a delay line, it can be applied to
any sound source.
*/
#include "../AudioApp.h"
#include "Gamma/Delay.h"
#include "Gamma/Oscillator.h"
using namespace gam;

class Vibrato{
public:
	Vibrato(float modAmount=1./400, float modFreq=5)
	:	modAmount(modAmount),
		delay(0.1, 0), mod(modFreq)
	{}
	
	float operator()(float i){
		delay.delay(mod.hann()*modAmount + 0.0001);
		return delay(i);
	}

	float modAmount;
	Delay<> delay;
	LFO<> mod;
};

class MyApp : public AudioApp{
public:

	LFO<> src1, src2;	// Sound sources
	Vibrato vibrato;	// Vibrato unit

	void onAudio(AudioIOData& io){

		src1.freq(220);
		src2.freq(220 * 1.26);

		while(io()){
			float s = (src1.tri() + src2.tri()) * 0.2;

			// Apply the vibrato
			s = vibrato(s);

			io.out(0) = io.out(1) = s;
		}
	}
};

int main(){
	MyApp().start();
}
