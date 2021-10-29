#ifndef GAMMA_AUDIO_APP_H_INC
#define GAMMA_AUDIO_APP_H_INC
/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Description: A very basic audio application abstraction
*/

#include <stdio.h>
#include "Gamma/AudioIO.h"
#include "Gamma/Domain.h"

namespace gam{

class AudioApp : public AudioCallback{
public:

	AudioApp(){
		mAudioIO.append(*this);
		initAudio(44100);
	}

	void initAudio(
		double framesPerSec, unsigned framesPerBuffer=128,
		unsigned outChans=2, unsigned inChans=0
	){
		mAudioIO.framesPerSecond(framesPerSec);
		mAudioIO.framesPerBuffer(framesPerBuffer);
		mAudioIO.channelsOut(outChans);
		mAudioIO.channelsIn(inChans);
		sampleRate(framesPerSec);
	}

	AudioIO& audioIO(){ return mAudioIO; }

	void start(bool block=true){
		mAudioIO.start();
		if(block){		
			printf("Press 'enter' to quit...\n");
			fflush(stdout);
			getchar();
		}
	}

private:
	AudioIO mAudioIO;
};

}

#endif
