/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Example:		IO / Audio Channel Test
	Description:	Plays a noise burst through each output channel.
*/

#include <stdio.h>
#include "Gamma/Gamma.h"
#include "Gamma/AudioIO.h"

using namespace gam;

struct TestSound{
	double env;

	TestSound(): env(1){}

	float operator()(){
		float r = rnd::uniS(1.);
		r *= env;
		env *= 0.9997;
		return r;
	}
	
	void reset(){ env = 1; }
	bool done(){ return env < 0.001; }
};

TestSound src;
int chan=-1;


void audioCB(AudioIOData & io){
	
	while(io()){

		if(src.done()){
			++chan;
			if(chan >= io.channelsOut()) chan=0;
			src.reset();
			printf("chan: %d\n", chan);
		}

		io.out(chan) = src();
	}
}


int main(){
	
	int maxOChans = AudioDevice::defaultOutput().channelsOutMax();
	int maxIChans = AudioDevice::defaultOutput().channelsOutMax();
	printf("Max input channels:  %d\n", maxIChans);
	printf("Max output channels: %d\n", maxOChans);

	// To open the maximum number of channels, we can hand in the queried values...
	//AudioIO io(256, 44100., audioCB, NULL, maxOChans, maxIChans);
	
	// ... or just use -1
	AudioIO io(256, 44100., audioCB, NULL, -1, -1);

	//Sync::master().spu(io.framesPerSecond());
	io.start();
	//io.print();
	
	printf("\nPress 'enter' to quit...\n"); getchar();
	return 0;
}
