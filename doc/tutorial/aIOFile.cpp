/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Tutorial:		Audio / IO
	Description:	Reading and writing a sound file.
*/

#include "tutorial.h"

int main(int argc, char* argv[]){

	SoundFile sf("./test.aiff");
	
	sf.format(SF_FORMAT_AIFF | SF_FORMAT_PCM_16);
	sf.channels(2);
	sf.frameRate(44100);
	sf.openWrite();
	
	const float sampeRate = 44100;
	const float lenSec = 2;
	const int numFrames = sampeRate * lenSec;
	
	float freq = 440;
	float buf[numFrames*2];
	
	for(int i=0; i<numFrames; i++){
		
		float s = i/float(numFrames) * lenSec;	// compute time in seconds
		float p1 = s * freq * M_2PI;
		float p2 = s * (freq+4) * M_2PI;
		
		buf[i*2  ] = sin(p1)*0.5;
		buf[i*2+1] = sin(p2)*0.5;
	}
	
	sf.write(buf, numFrames);
	
	return 0;
}
