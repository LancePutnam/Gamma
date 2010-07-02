/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Description:	Include files for Gamma tutorials.	
*/

#include <stdio.h>				// for printing to stdout
#include "Gamma/Gamma.h"		// core functions
#include "Gamma/Access.h"
#include "Gamma/AudioIO.h"
#include "Gamma/Delay.h"
#include "Gamma/Effects.h"
#include "Gamma/Envelope.h"
#include "Gamma/Noise.h"
#include "Gamma/Oscillator.h"
#include "Gamma/Print.h"
#include "Gamma/SoundFile.h"
#include "Gamma/Types.h"

using namespace gam;

#define RUN(audioCB)\
int main(int argc, char* argv[]){\
	AudioIO io(256, 44100, audioCB, NULL, 2);\
	Sync::master().spu(io.framesPerSecond());\
	io.start();\
	printf("\nPress 'enter' to quit...\n"); getchar();\
	return 0;\
}
