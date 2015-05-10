/*	Description:
 
 */

#include "examples.h"

class OscTrm : public Process<AudioIOData> {
public:
	
	OscTrm(double startTime=0)
	:	mAmp(1), mDur(2)
	{
		dt(startTime);
		set(10, 262, 0.5, 0.1,2,0.8, 0.4,4,8,0.5, mOsc, 0.8);
		mAmpEnv.levels(0,1,1,0);
	}
	
	OscTrm& freq(float v){ mOsc.freq(v); return *this; }
	OscTrm& amp(float v){ mAmp=v; return *this; }
	OscTrm& dur(float v){ mDur=v; return *this; }
	OscTrm& attack(float v){ mAmpEnv.lengths()[0]=v; return *this; }
	OscTrm& decay(float v){  mAmpEnv.lengths()[2]=v; return *this; }
	OscTrm& sus(float v){ mAmpEnv.levels()[2]=v; return *this; }
	
	OscTrm& trm1(float v){ mTrmEnv.levels(v, mTrmEnv.levels()[1], v); return *this; }
	OscTrm& trm2(float v){ mTrmEnv.levels()[1]=v; return *this; }
	OscTrm& trmDepth(float v){ mTrmDepth=v; return *this; }
	OscTrm& trmRise(float v){ mTrmEnv.lengths(v,1-v); return *this; }
	
	OscTrm& table(ArrayPow2<float>& v){ mOsc.source(v); return *this; }
	
	OscTrm& pan(float v){ mPan.pos(v); return *this; }
	
	
	OscTrm& set(
		float a, float b, float c, float d, float e, float f,
		float g, float h, float i, float j,
		ArrayPow2<float>& k, float l
	){
		return dur(a).freq(b).amp(c).attack(d).decay(e).sus(f)
			.trmDepth(g).trm1(h).trm2(i).trmRise(j)
			.table(k).pan(l);
	}
	
	void onProcess(AudioIOData& io){
		
		mAmpEnv.totalLength(mDur, 1);
		mTrmEnv.totalLength(mDur);
		
		while(io()){
			
			mTrm.freq(mTrmEnv());
			//float trmAmp = mAmp - mTrm()*mTrmDepth; // Replaced with line below
			float trmAmp = (mTrm()*0.5+0.5)*mTrmDepth + (1-mTrmDepth); // Corrected 
			float s1 = mOsc() * mAmpEnv() * trmAmp * mAmp;
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
	
	Sine<> mTrm;
	float mTrmDepth;
	Env<2> mTrmEnv;
	Osc<> mOsc;
	Env<3> mAmpEnv;
};


int main(){
	
	ArrayPow2<float> tbSin(2048), tbSqr(2048), tbSaw(2048);
	addWave(tbSin, SINE);
	addWave(tbSqr, SQUARE);
	addWave(tbSaw, SAW);
	
	Scheduler s;
	s.add<OscTrm>( 0).set(10, 262, 0.5, 0.1,2,0.8, 0.4,4,8,0.5, tbSin, 0.8); 
	
	s.add<OscTrm>(11).set(10, 262, 0.5, 0.1,2,0.8, 0.4,4,8,0.5, tbSqr, 0.8); 
	
	s.add<OscTrm>(20).freq(120).table(tbSaw);
	
	s.add<OscTrm>(30).set(20, 262, 0.5, 0.1,2,0.8, 0.4,4,80,0.2, tbSin, 0.8); 
	
	s.add<OscTrm>(50).set(10, 262, 0.5, 0.1,2,0.8, 0.8,4,8,0.5, tbSin, 0.8); 
	
	AudioIO io(256, 44100., Scheduler::audioCB, &s);
	gam::sampleRate(io.fps());
	io.start();
	printf("\nPress 'enter' to quit...\n"); getchar();
}