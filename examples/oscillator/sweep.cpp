/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information

Example:	Sine Oscillator
Author:		Lance Putnam, 2015

Description:
Here we use Sweep, a periodic ramp, to scan through a table to containing a 
sine wave.
*/

#include "../AudioApp.h"
#include "Gamma/Oscillator.h"
using namespace gam;

class MyApp : public AudioApp{
public:

	static const int N=2048;
	float table[N]; // Table to store sine wave
	Sweep<> sweep;

	MyApp(){
		// Fill sine table
		for(int i=0; i<N; ++i)
			table[i] = sin(float(i)/N*M_2PI);
	}

	void onAudio(AudioIOData& io){

		sweep.freq(440);

		while(io()){
			float s = sweep();		// upward ramp between 0 and 1
			s = table[int(s*N)];

			io.out(0) = io.out(1) = s*0.2;
		}
	}

};

int main(){
	MyApp().start();
}
