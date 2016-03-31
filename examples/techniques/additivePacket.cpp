/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example:	Wavepacket with Dispersion
Author:		Lance Putnam, 2012

Description:

*/
#include "../AudioApp.h"
#include "Gamma/Oscillator.h"
using namespace gam;

class MyApp : public AudioApp{
public:

	SineRs<> oscs;

	MyApp(){
		oscs.resize(512);

		float ff = 27.5;			// Fundamental frequency
		float k1 = 0.;				// Velocity
		float k2 = 0.00001;			// Dispersion
		float rs = 1./oscs.size();

		for(unsigned i=0; i<oscs.size(); ++i){
			float h = i+1 + i*k1 + i*i*k2;		
			float a = 16 * rs/h;
			float p = 0;
			oscs.set(i,  ff*h,a,p);
		}
	}

	void onAudio(AudioIOData& io){
		while(io()){
			io.out(0) = io.out(1) = oscs()*4;
		}
	}
};

int main(){
	MyApp().start();
}
