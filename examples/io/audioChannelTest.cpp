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

	AudioDevice adevi = AudioDevice::defaultInput();
	AudioDevice adevo = AudioDevice::defaultOutput();
	//AudioDevice adevi = AudioDevice("Microphone (AUREON");
	//AudioDevice adevo = AudioDevice("Speakers (AUREON");
	//AudioDevice adevi = AudioDevice(11);
	//AudioDevice adevo = AudioDevice(16);
		
	int maxIChans = adevi.channelsInMax();
	int maxOChans = adevo.channelsOutMax();
	printf("Max input channels:  %d\n", maxIChans);
	printf("Max output channels: %d\n", maxOChans);

	AudioIO io(256, 44100., audioCB);
	
	// Set i/o devices
	io.deviceIn(adevi);
	io.deviceOut(adevo);
	
	// Use the maximum number of channels
	io.channelsIn(-1);
	io.channelsOut(-1);

	//gam::sampleRate(io.framesPerSecond());
	if(io.start()){
		printf("start successful\n");
		io.print();
		printf("\nPress 'enter' to quit...\n"); getchar();
	}
	else{
		printf("start failed\n");
	}
}
