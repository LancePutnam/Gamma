#include <stdio.h>
#include "AudioIO.h"

using namespace gam;

// Audio callback definition
void audioCB(AudioIOData & io){
	const float * in0 = io.in(0);
	float * out0 = io.out(0);
	float * out1 = io.out(1);
	
	for(unsigned long i=0; i<io.numFrames(); ++i){
		out0[i] = out1[i] = -in0[i];
	}
}

int main(int argc, char* argv[]){
	
	AudioIO::initDevices();
	AudioIO::printDevices();
	
	//AudioIO io(128, 44100., audioCB, NULL, 16, 16);
	AudioIO io(128, 44100., audioCB, NULL, 2, 1);
	io.start();
	io.print();

	printf("\nPress 'enter' to quit...\n"); getchar();
	return 0;
}

