/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example:	Real-time Recording
Author:		Lance Putnam, 2012

Description:
This demonstrates how to record audio to a sound file in real-time.
*/
#include "../AudioApp.h"
#include "Gamma/Oscillator.h"
#include "Gamma/SoundFile.h"
#include "Gamma/Recorder.h"
#include "Gamma/Timer.h"
using namespace gam;

class MyApp : public AudioApp{
public:

	gam::SoundFile sf;		// The sound file of the recording
	gam::Recorder rec;		// Create a recorder
	gam::Sine<> src1, src2;	// Sine waves to record

	MyApp(){
		sf.channels(2);
		sf.frameRate(44100);
		sf.path("recording.wav");
		sf.openWrite();

		rec.resize(
			2,		// number of channels
			8192	// size of buffer, in frames
		);

		src1.freq(440);
		src2.freq(441);

		start(false);

		// Write samples from ring buffer into sound file.
		// This will typically be done in a separate lower-priority thread, 
		// definitely not here...
		int i=200;
		while(--i){
			float * buf;
			int n = rec.read(buf);
			sf.write(buf, n);
			//printf("consumed %d frames\n", n);
			gam::sleepSec(0.01);
		}
	}

	void onAudio(gam::AudioIOData& io){
		while(io()){
			float s1 = src1()*0.25;
			float s2 = src2()*0.25;
			rec.write(s1, s2);

			io.out(0) = s1;
			io.out(1) = s2;
		}
	}

};

int main(){
	MyApp();
}

