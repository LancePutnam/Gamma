/*	Description:	
*/

#include "examples.h"

class Sub : public Process<AudioIOData> {
public:

	Sub(double startTime=0)
	:	mAmp(1), mDur(2)
	{
		dt(startTime);
		mAmpEnv.curve(0); // linear segments
		mAmpEnv.levels(0,1,1,0);
		hnum(12);
	}

	Sub& freq(float v){ mOsc.freq(v); return *this; }
	Sub& amp(float v){ mAmp=v; return *this; }
	Sub& attack(float v){ mAmpEnv.lengths()[0]=v; return *this; }
	Sub& decay (float v){ mAmpEnv.lengths()[2]=v; return *this; }
	Sub& dur(float v){ mDur=v; return *this; }
	Sub& pan(float v){ mPan.pos(v); return *this; }

	Sub& noise(float v){ mNoiseMix=v; return *this; }
	Sub& hnum(int i){ mOsc.harmonics(i); return *this; }
	Sub& hamp(float v){ mOsc.ampRatio(v); return *this; }
	Sub& cf1(float v){ mCFEnv.levels(v, mCFEnv.levels()[1], v); return *this; }
	Sub& cf2(float v){ mCFEnv.levels()[1]=v; return *this; }
	Sub& cfRise(float v){ mCFEnv.lengths(v, 1-v); return *this; }
	Sub& bw1(float v){ mBWEnv.levels(v, mBWEnv.levels()[1], v); return *this; }
	Sub& bw2(float v){ mBWEnv.levels()[1]=v; return *this; }
	Sub& bwRise(float v){ mBWEnv.lengths(v, 1-v); return *this; }
	//Sub& hlow(int i){ return *this; }

	Sub& set(
		float a, float b, float c, float d, float e, float f,
		float g, float h, float i,
		float j, float k, float l,
		float m, float n,
		float o=0
	){	return
		dur(a).freq(b).amp(c).attack(d).decay(e).noise(f)
		.cf1(g).cf2(h).cfRise(i)
		.bw1(j).bw2(k).bwRise(l)
		.hnum(m).hamp(n)
		.pan(o);
	}

	//
	void onProcess(AudioIOData& io){

		mAmpEnv.totalLength(mDur, 1);
		mCFEnv.totalLength(mDur);
		mBWEnv.totalLength(mDur);

		while(io()){
			// mix oscillator with noise
			float s1 = mOsc()*(1-mNoiseMix) + mNoise()*mNoiseMix;

			// apply resonant filter
			mRes.set(mCFEnv(), mBWEnv());
			s1 = mRes(s1);

			// appy amplitude envelope
			s1 *= mAmpEnv() * mAmp;

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

	DSF<> mOsc;
	NoiseWhite<> mNoise;
	float mNoiseMix;
	Reson<> mRes;
	Env<3> mAmpEnv;
	Env<2> mCFEnv;
	Env<2> mBWEnv;
};



int main(){

	Scheduler s;
	s.add<Sub>( 0).set(5, 220, 0.3, 0.1, 0.1, 0,  50,5000,1.0,  20,  20,0.5, 40,1.0);
	s.add<Sub>( 6).set(5, 220, 0.3, 0.1, 0.1, 0,  50,  50,1.0,  20,2000,0.5, 10,1.0);
	s.add<Sub>(12).set(5,   0, 0.8, 0.1, 0.1, 1, 262, 262,1.0, 100,5000,0.5, 10,1.0);
	s.add<Sub>(18).set(5,   0, 0.8, 0.1, 0.1, 1, 100,1000,1.0,  10,  10,1.0, 10,1.0);
	s.add<Sub>(24).set(5, 120, 0.8, 1.0, 1.0, 0,  50,  50,1.0,  20,2000,0.5, 20,0.5);
	s.add<Sub>(30).set(5,   0, 0.8, 1.0, 1.0, 1, 100, 100,1.0, 100,5000,0.5, 20,0.5);

	AudioIO io(256, 44100., Scheduler::audioCB, &s);
	gam::sampleRate(io.fps());
	io.start();
	printf("\nPress 'enter' to quit...\n"); getchar();
}
