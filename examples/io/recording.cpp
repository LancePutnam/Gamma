/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Example:		IO / Recording
	Description:	This demonstrates how to record live audio to a sound file.
*/

#include "../examples.h"

SoundFile sf("recording.aif");
Recorder rec(2);
Sine<> src1(440), src2(441);

void audioCB(AudioIOData& io){

	while(io()){		
		float s1 = src1();
		float s2 = src2();
		rec.write(s1, s2);
		io.sum(s1*0.2, 0);
		io.sum(s2*0.2, 1);
	}
}


int main(){

	sf.channels(2);
	sf.frameRate(44100);
	sf.openWrite();

	AudioIO io(256, 44100, audioCB, 0, 2);
	Sync::master().spu(io.framesPerSecond());
	io.start();

	// Write samples from ring buffer into sound file.
	// This will typically be done in a separate lower-priority thread, not main.
	int i=200;
	while(--i){
		float * buf;
		int n = rec.read(buf);
		sf.write(buf, n);
		//printf("%d\n", n);
		sleepSec(0.01);
	}

	//printf("\nPress 'enter' to quit...\n"); getchar();
	return 0;
}
