/*	Description:

*/

#include "examples.h"

class OscAM : public Process<AudioIOData> {
public:

	OscAM(double startTime=0)
	:	mAmp(1), mDur(2)
	{
		dt(startTime);
		mAmpEnv.levels(0,1,1,0);
		mAMEnv.curve(0);
	}

	OscAM& freq(float v){ mOsc.freq(v); return *this; }
	OscAM& amp(float v){ mAmp=v; return *this; }
	OscAM& dur(float v){ mDur=v; return *this; }
	OscAM& attack(float v){ mAmpEnv.lengths()[0]=v; return *this; }
	OscAM& decay(float v){  mAmpEnv.lengths()[2]=v; return *this; }
	OscAM& sus(float v){ mAmpEnv.levels()[2]=v; return *this; }
	
	OscAM& am1(float v){ mAMEnv.levels(v, mAMEnv.levels()[1], v); return *this; }
	OscAM& am2(float v){ mAMEnv.levels()[1]=v; return *this; }
	OscAM& amRatio(float v){ mAMRatio=v; return *this; }
	OscAM& amRise(float v){ mAMEnv.lengths(v,1-v); return *this; }
	OscAM& amFunc(ArrayPow2<float>& v){ mAM.source(v); return *this; }
	
	OscAM& pan(float v){ mPan.pos(v); return *this; }


	OscAM& set(
		float a, float b, float c, float d, float e, float f,
		float g, float h, float i, float j, ArrayPow2<float>& k,
		float l=0
	){
		return dur(a).freq(b).amp(c).attack(d).decay(e).sus(f)
			.am1(g).am2(h).amRise(i).amRatio(j).amFunc(k)
			.pan(l);
	}

	void onProcess(AudioIOData& io){

		mAmpEnv.totalLength(mDur, 1);
		mAMEnv.totalLength(mDur);

		while(io()){
		
			mAM.freq(mOsc.freq()*mAMRatio);			// set AM freq according to ratio
			float amAmt = mAMEnv();					// AM amount envelope
		
			float s1 = mOsc();						// non-modulated signal
			s1 = s1*(1-amAmt) + (s1*mAM())*amAmt;	// mix modulated and non-modulated
			
			s1 *= mAmpEnv() *mAmp;
			
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

	Osc<> mAM;
	float mAMRatio;
	Env<2> mAMEnv;
	Sine<> mOsc;
	Env<3> mAmpEnv;
};


int main(){

	ArrayPow2<float> tbSin(2048), tbSqr(2048), tbPls(2048), tbDin(2048);
	addWave(tbSin, SINE);
	addWave(tbSqr, SQUARE);
	addWave(tbPls, IMPULSE, 4);

	// inharmonic partials
	{	float A[] = {1, 0.7, 0.45, 0.3, 0.15, 0.08};
		float C[] = {10, 27, 54, 81, 108, 135};
		addSines(tbDin, A,C,6);
	}

	Scheduler s;
	s.add<OscAM>( 0).set(5, 262, 0.5, 0.1,0.08,0.8, 0.2,0.8,0.5, 1.4870, tbSin);
	s.add<OscAM>( 5).set(5, 262, 0.5, 0.1,0.08,0.8, 0.8,0.2,0.5, 2.0001, tbSin);
	s.add<OscAM>(10).set(5, 262, 0.5, 0.1,0.08,0.8, 0.2,0.8,1.0, 2.0001, tbPls);
	s.add<OscAM>(15).set(5, 262, 0.5, 0.1,0.08,0.8, 0.2,0.8,1.0, 2.0001/10., tbDin);

	AudioIO io(256, 44100., Scheduler::audioCB, &s);
	gam::sampleRate(io.fps());
	io.start();
	printf("\nPress 'enter' to quit...\n"); getchar();
}
