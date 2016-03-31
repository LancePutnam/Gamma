/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Example:		
	Description:	
*/

#include "examples.h"

class SineEnv : public Process<AudioIOData> {
public:

	SineEnv(double startTime=0)
	{
		dt(startTime);
		set (6.5, 260, 0.3, 1, 2);
		mAmpEnv.curve(0); // make segments lines
		mAmpEnv.levels(0,1,1,0);
	}

	SineEnv& freq(float v){ mOsc.freq(v); return *this; }
	SineEnv& amp(float v){ mAmp=v; return *this; }
	SineEnv& attack(float v){
		mAmpEnv.lengths()[0] = v;
		return *this;
	}
	SineEnv& decay(float v){
		mAmpEnv.lengths()[2] = v;
		return *this;
	}
	SineEnv& dur(float v){ mDur=v; return *this; }
	SineEnv& pan(float v){ mPan.pos(v); return *this; }
	SineEnv& set(
		float a, float b, float c, float d, float e, float f=0
	){
		return dur(a).freq(b).amp(c).attack(d).decay(e).pan(f);
	}

	//
	void onProcess(AudioIOData& io){

		mAmpEnv.totalLength(mDur, 1);

		while(io()){
			float s1 = mOsc() * mAmpEnv() * mAmp;
			float s2;
			mPan(s1, s1,s2);
			io.out(0) += s1;
			io.out(1) += s2;
		}
		if(mAmpEnv.done()) free();
	}

protected:
	float mAmp;
	float mDur;
	Pan<> mPan;
	Sine<> mOsc;
	Env<3> mAmpEnv;
};



int main(){

	Scheduler s;
	s.add<SineEnv>( 0  ).set(0.5, 260, 0.3, .011, .2);
	s.add<SineEnv>( 0  ).set(0.5, 510, 0.3, .011, .2);
	s.add<SineEnv>( 3.5).set(0.5, 233, 0.3, .011, .2);
	s.add<SineEnv>( 3.5).set(0.5, 340, 0.3, .011, .2);
	s.add<SineEnv>( 3.5).set(7.5, 710, 0.3, 1, 2);
	s.add<SineEnv>( 6.5).set(1.5,  60, 0.3, 1, 2);
	s.add<SineEnv>( 8.0).set(3.5, 230, 0.3, 1, 2);
	s.add<SineEnv>(11.5).set(7.5, 450, 0.3, 1, 2);
	s.add<SineEnv>(11.5).set(3.0, 240, 0.3, 1, 2);
	s.add<SineEnv>(11.5).set(1.0,  90, 0.3, 1, 2);
	s.add<SineEnv>(12.5).set(2.0, 120, 0.3, 1, 2);
	s.add<SineEnv>(14.5).set(4.5, 380, 0.3, 1, 2);
	s.add<SineEnv>(14.5).set(3.0, 240, 0.3, 1, 2);
	s.add<SineEnv>(17.5).set(1.5, 110, 0.3, 1, 2);
	s.add<SineEnv>(19.0).set(8.0, 940, 0.3, 1, 2);
	s.add<SineEnv>(19.0).set(8.0, 780, 0.3, 1, 2);
	s.add<SineEnv>(19.0).set(8.0, 640, 0.3, 1, 2);
	s.add<SineEnv>(27.0).freq(233).amp(0.1);
	s.add<SineEnv>(27.0).freq(234).amp(0.1);
	s.add<SineEnv>(27.0).freq(235).amp(0.1);
	s.add<SineEnv>(27.0).freq(232).amp(0.1);
	s.add<SineEnv>(27.0).freq(231).amp(0.1);
	s.add<SineEnv>(27.0).freq(230).amp(0.1);
	s.add<SineEnv>(31.0).freq(329.6).amp(0.1);
	s.add<SineEnv>(31.0).freq(330.6).amp(0.1);
	s.add<SineEnv>(31.0).freq(331.6).amp(0.1);
	s.add<SineEnv>(31.0).freq(334.6).amp(0.1);
	s.add<SineEnv>(31.0).freq(335.6).amp(0.1);
	s.add<SineEnv>(31.0).freq(337.6).amp(0.1);
	s.add<SineEnv>(35.0).freq(233).amp(0.1);
	s.add<SineEnv>(35.0).freq(234).amp(0.1);
	s.add<SineEnv>(35.0).freq(235).amp(0.1);
	s.add<SineEnv>(35.0).freq(232).amp(0.1);
	s.add<SineEnv>(35.0).freq(231).amp(0.1);
	s.add<SineEnv>(35.0).freq(230).amp(0.1);
	s.add<SineEnv>(39.0).freq(440).amp(0.1);
	s.add<SineEnv>(39.0).freq(441).amp(0.1);
	s.add<SineEnv>(39.0).freq(442).amp(0.1);
	s.add<SineEnv>(39.0).freq(443).amp(0.1);
	s.add<SineEnv>(39.0).freq(445).amp(0.1);
	s.add<SineEnv>(39.0).freq(446).amp(0.1);
	s.add<SineEnv>(45.0).set(8.0, 940, 0.1, 1, 2);
	s.add<SineEnv>(45.0).set(8.0, 941, 0.1, 1, 2);
	s.add<SineEnv>(45.0).set(8.0, 942, 0.1, 1, 2);
	s.add<SineEnv>(47.0).set(8.0, 780, 0.1, 1, 2);
	s.add<SineEnv>(47.0).set(8.0, 781, 0.1, 1, 2);
	s.add<SineEnv>(47.0).set(8.0, 782, 0.1, 1, 2);
	s.add<SineEnv>(49.0).set(8.0, 640, 0.1, 1, 2);
	s.add<SineEnv>(49.0).set(8.0, 641, 0.1, 1, 2);
	s.add<SineEnv>(49.0).set(8.0, 642, 0.1, 1, 2);
	s.add<SineEnv>(51.5).set(7.5, 450, 0.3, 1, 2);
	s.add<SineEnv>(51.5).set(3.0, 240, 0.3, 1, 2);
	s.add<SineEnv>(51.5).set(1.0,  90, 0.3, 1, 2);
	s.add<SineEnv>(52.5).set(2.0, 120, 0.3, 1, 2);
	s.add<SineEnv>(54.5).set(4.5, 380, 0.3, 1, 2);
	s.add<SineEnv>(54.5).set(3.0, 240, 0.3, 1, 2);
	s.add<SineEnv>(57.5).set(1.5, 110, 0.3, 1, 2);
	s.add<SineEnv>( 68  ).set(6.5, 260, 0.3, 1, 2);
	s.add<SineEnv>( 68  ).set(3.5, 510, 0.3, 1, 2);
	s.add<SineEnv>( 65.5).set(3.0, 233, 0.3, 1, 2);
	s.add<SineEnv>( 65.5).set(4.5, 340, 0.3, 1, 2);
	s.add<SineEnv>( 65.5).set(7.5, 710, 0.3, 1, 2);
	s.add<SineEnv>( 62.5).set(1.5,  60, 0.3, 1, 2);
	s.add<SineEnv>( 60.0).set(3.5, 230, 0.3, 1, 2);
	
	AudioIO io(256, 44100., Scheduler::audioCB, &s);
	gam::sampleRate(io.fps());
	io.start();
	printf("\nPress 'enter' to quit...\n"); getchar();
}
