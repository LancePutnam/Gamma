/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Description:	Include files for Gamma examples	
*/

#include <stdio.h>				// for printing to stdout
#define GAMMA_H_INC_ALL			// define this to include all header files
#include "Gamma/Gamma.h"
using namespace gam;

#define RUN_AUDIO_MAIN \
int main(int argc, char* argv[]){\
	AudioIO io(256, 44100, audioCB, NULL, 2);\
	Domain::master().spu(io.framesPerSecond());\
	io.start();\
	printf("Press 'enter' to quit...\n"); getchar();\
}

