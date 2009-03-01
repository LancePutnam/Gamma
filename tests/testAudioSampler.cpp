#include "rnd.h"
#include "AudioIO.h"
#include "Sampler.h"

using namespace gam;

typedef Player<float, ipl::Linear, tap::Clip> Player_t;
const int numSmps = 4;

#define SF_ROOT "../../../../sounds/"
Player_t smps[numSmps] = {
	Player_t(SF_ROOT"water1.wav"), Player_t(SF_ROOT"water2.wav"),
	Player_t(SF_ROOT"water3.wav"), Player_t(SF_ROOT"water4.wav")
};

Player_t * smpA = smps;
bool loops = false;

void audioCB(AudioIOData & io){
	float * out0 = io.out(0);
	float * out1 = io.out(1);
	unsigned long numFrames = io.numFrames();
	
	using namespace gam::rnd;
	if(prob(0.008)){						// Randomize sampler
		smpA = &smps[rnd::uni(numSmps)];		
		smpA->pos(0.);
		smpA->rate(pow(2., rnd::uni(-2., 1.)));
		
		if(prob(0.1)){
			smpA->rate(-smpA->rate());
			smpA->phase(1.);
		}		
		loops = prob(0.2);
	}

	for(unsigned long f=0; f<numFrames; f++){
		out0[f] = loops ? (*smpA)() : (*smpA)();
		//out0[f] = scl::intToNormal(loops ? smpA->nextN() : smpA->nextClipN());
		out1[f] = out0[f];
	}
}

int main(int argc, char* argv[]){
	AudioIO io(256, 44100., audioCB, NULL, 2, 1);
	Sync::master().spu(io.fps());

	//if(!smpA.hasSource()){ printf("Can't open sound file...\n"); return -1; }
	
	io.start();

	printf("\nPress 'enter' to quit...\n"); getchar();
	return 0;
}
