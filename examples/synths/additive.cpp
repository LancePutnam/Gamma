/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Example:		
	Description:	
*/

#include "examples.h"

class AddSyn : public Process<AudioIOData> {
public:

	AddSyn(double startTime=0)
	{
		dt(startTime);
		set(6.2,155.6,0.01,0.5,0.1,0.1,0.8,0.5,0.001,0.1,0.8,0.6,0.01,0.075,0.9,1,2.001,3,4.00009,5.0002,6,7,8,9);
		mEnvStri.curve(-4); // make segments lines
		mEnvStri.levels(0,1,1,0);
		mEnvLow.curve(-4); // make segments lines
		mEnvLow.levels(0,1,1,0);
		mEnvUp.curve(-4); // make segments lines
		mEnvUp.levels(0,1,1,0);
	}

	
	AddSyn& freq(float v){ 
		mOscFrq=v; 

		mOsc1.freq(mfreqStri1*mOscFrq); 
		mOsc2.freq(mfreqStri2*mOscFrq); 
		mOsc3.freq(mfreqStri3*mOscFrq);
		mOsc4.freq(mfreqLow1*mOscFrq);
		mOsc5.freq(mfreqLow2*mOscFrq);
		mOsc6.freq(mfreqUp1*mOscFrq);
		mOsc7.freq(mfreqUp2*mOscFrq); 
		mOsc8.freq(mfreqUp3*mOscFrq);
		mOsc9.freq(mfreqUp4*mOscFrq);		
		return *this; 
	}
	
	AddSyn& freqStri1(float v){ mfreqStri1=v; mOsc1.freq(v*mOscFrq); return *this; }
	AddSyn& freqStri2(float v){ mfreqStri2=v; mOsc2.freq(v*mOscFrq); return *this; }
	AddSyn& freqStri3(float v){ mfreqStri3=v; mOsc3.freq(v*mOscFrq); return *this; }
	AddSyn& freqLow1(float v){ mfreqLow1=v; mOsc4.freq(v*mOscFrq); return *this; }
	AddSyn& freqLow2(float v){ mfreqLow2=v; mOsc5.freq(v*mOscFrq); return *this; }
	AddSyn& freqUp1(float v){ mfreqUp1=v; mOsc6.freq(v*mOscFrq); return *this; }
	AddSyn& freqUp2(float v){ mfreqUp2=v; mOsc7.freq(v*mOscFrq); return *this; }
	AddSyn& freqUp3(float v){ mfreqUp3=v; mOsc8.freq(v*mOscFrq); return *this; }
	AddSyn& freqUp4(float v){ mfreqUp4=v; mOsc9.freq(v*mOscFrq); return *this; }
	AddSyn& freqUp(float v1, float v2, float v3, float v4) { 
		mOsc6.freq(v1*mOscFrq);
		mOsc7.freq(v2*mOscFrq);
		mOsc8.freq(v3*mOscFrq);
		mOsc9.freq(v4*mOscFrq);
		return *this;
	}
	
	AddSyn& freq(float v1, float v2) {
		dur(v1);
		freq(v2);
		return *this;
	}
	
	AddSyn& amp(float v){ mAmp=v; return *this; }
	AddSyn& ampStri(float v){ mAmpStri=v; return *this; }
	AddSyn& attackStri(float v){ mEnvStri.lengths()[0] = v; return *this; }
	AddSyn& decayStri(float v){ mEnvStri.lengths()[2] = v; return *this; }
	AddSyn& susStri(float v){ mEnvStri.levels()[2]=v; return *this; }
	AddSyn& ampLow(float v){ mAmpLow=v; return *this; }
	AddSyn& attackLow(float v){ mEnvLow.lengths()[0] = v; return *this; }
	AddSyn& decayLow(float v){ mEnvLow.lengths()[2] = v; return *this; }
	AddSyn& susLow(float v){ mEnvLow.levels()[2]=v; return *this; }
	AddSyn& ampUp(float v){ mAmpUp=v; return *this; }
	AddSyn& attackUp(float v){ mEnvUp.lengths()[0] = v; return *this; }
	AddSyn& decayUp(float v){ mEnvUp.lengths()[2] = v; return *this; }
	AddSyn& susUp(float v){ mEnvUp.levels()[2]=v; return *this; }
	
	
	AddSyn& dur(float v){ mDur=v; return *this; }

	AddSyn& pan(float v){ mPan.pos(v); return *this; }

	AddSyn& set(
		float a, float b, float c, float d, float e,float f, float g, float h, float i, float j, 
				float k, float l, float m, float n, float o, float p, float q, float r, float s, float t,
				float u, float v, float w, float x, float y=0
	){
		return dur(a).freq(b).amp(c).ampStri(d).attackStri(e).decayStri(f).susStri(g).ampLow(h)
		.attackLow(i).decayLow(j).susLow(k).ampUp(l).attackUp(m).decayUp(n).susUp(o)
		.freqStri1(p).freqStri2(q).freqStri3(r).freqLow1(s).freqLow2(t).freqUp1(u)
		.freqUp2(v).freqUp3(w).freqUp4(x).pan(y);
	}

	//
	void onProcess(AudioIOData& io){

		mEnvStri.totalLength(mDur, 1);
		mEnvLow.totalLength(mDur, 1);
		mEnvUp.totalLength(mDur, 1);
		
		while(io()){
			float s1 = (mOsc1() + mOsc2() + mOsc3()) * mEnvStri() * mAmpStri;
			s1 += (mOsc4() + mOsc5()) * mEnvLow() * mAmpLow;
			s1 += (mOsc6() + mOsc7() + mOsc8() + mOsc9()) * mEnvUp() * mAmpUp;
			s1 *= mAmp;
			float s2;
			mPan(s1, s1,s2);
			io.out(0) += s1;
			io.out(1) += s2;
		}
		if(mEnvStri.done()) free();
	}


protected:
	float mAmp;
	float mAmpStri;
	float mAmpLow;
	float mAmpUp;
	float mDur;
	float mOscFrq;
	
	float mfreqStri1, mfreqStri2, mfreqStri3, mfreqLow1, mfreqLow2, mfreqUp1, mfreqUp2, mfreqUp3, mfreqUp4;
	
	
	Pan<> mPan;
	Sine<> mOsc;
	Sine<> mOsc1;
	Sine<> mOsc2;
	Sine<> mOsc3;
	Sine<> mOsc4;
	Sine<> mOsc5;
	Sine<> mOsc6;
	Sine<> mOsc7;
	Sine<> mOsc8;
	Sine<> mOsc9;
	Env<3> mEnvStri;
	Env<3> mEnvLow;
	Env<3> mEnvUp;
};

class Chimes : public AddSyn {
public:
	Chimes(double dt=0): AddSyn(dt) {
		set (6.2,440,0.1,0.5,0.0001,3.8,0.3,0.4,0.0001,6.0,0.99,0.3,0.0001,6.0,0.9,2,3,4.07,0.56,0.92,1.19,1.7,2.75,3.36);
	}
};


void addOneNote(Scheduler &s) {
	s.add<Chimes>(0.0).freq(6.2,622.2);
}


float harmonicSeriesScale[20];

void initScaleToHarmonicSeries() {
	for (int i=0;i<20;++i) {
		harmonicSeriesScale[i] = 100*i;
	}
}
float randomFromHarmonicSeries() {
	int index = rnd::uni(0,20);
	// std::cout << "index " << index << " is " << myScale[index] << std::endl;
	return harmonicSeriesScale[index];
}


float halfStepScale[20];
float halfStepInterval = 1.05946309; // 2^(1/12)
void initScaleTo12TET(float lowest) {
	float f = lowest;
	for (int i=0;i<20;++i) {
		halfStepScale[i] = f;
		f *= halfStepInterval;
	}
}

float randomFrom12TET() {
	int index = rnd::uni(0,20);
	// std::cout << "index " << index << " is " << myScale[index] << std::endl;
	return halfStepScale[index];
}


float myScale[6];
void initScaleToMyPitches() {
	myScale[0] = 100;
	myScale[1] = 121.245;
	myScale[2] = 145.876234;
	myScale[3] = 167.4786234;
	myScale[4] = 367.4786234;
	myScale[5] = 1367.4786234;
}

float randomFromMyScale() {
	int index = rnd::uni(0,6);
	std::cout << "index " << index << " is " << myScale[index] << std::endl;
	return myScale[index];
}


// Chooses frequencies from a random uniform distribution
void fillTime(Scheduler &s, float from, float to, float minNoteDur, float maxNoteDur, float minFreq, float maxFreq) {
	while (from <= to) {
		float nextDur = rnd::uni(minNoteDur,maxNoteDur);
		s.add<Chimes>(from).dur(nextDur).freq(rnd::uni(minFreq,maxFreq)).amp(.3);
		std::cout << "old from " << from << " plus nextDur " << nextDur << std::endl;
		from += nextDur;
	}
}

void fillTimeWith12TET(Scheduler &s, float from, float to, float minNoteDur, float maxNoteDur) {
	while (from <= to) {
		float nextDur = rnd::uni(minNoteDur,maxNoteDur);
		float f = randomFrom12TET();
		s.add<Chimes>(from).dur(nextDur).freq(f).amp(.3);
		std::cout << "12 old from " << from << " plus nextDur " << nextDur << std::endl;
		from += nextDur;
	}
}





int main(){
	// initScaleTo12TET(110);
	initScaleToMyPitches();
	initScaleTo12TET(110);
	
	Scheduler s;
	
	s.add<AddSyn>(0).set(6.2,155.6,0.1,0.5,0.0001,3.8,0.3,0.4,0.0001,6.0,0.99,0.3,0.0001,6.0,0.9,2,3,4.07,0.56,0.92,1.19,1.7,2.75,3.36);
	s.add<AddSyn>(7.0).set(6.2,622.2,0.1,0.5,0.0001,6.1,0.99,0.4,0.0005,6.1,0.99,0.3,0.0005,6.1,0.9,2,3,4.07,0.56,0.92,1.19,1.7,2.75,3.36);
	s.add<AddSyn>(14.0).set(6.2,155.6,0.01,0.5,0.1,0.1,0.8,0.5,0.001,0.1,0.8,0.6,0.01,0.075,0.9,1,2.001,3,4.00009,5.0002,6,7,8,9); 
	s.add<AddSyn>(21.0).set(6.2,77.78,0.01,0.5,0.1,0.4,0.8,0.5,0.001,0.4,0.8,0.6,0.01,0.4,0.5,1,2.0001,3,4.00009,5.0002,6,7,8,9); 
	s.add<AddSyn>(28.0).set(6.2,311.1,0.01,0.5,0.1,0.4,0.8,0.5,0.001,0.4,0.8,0.6,0.01,0.4,0.5,1,1.0001,3,3.0009,5.0002,5,7,7.0009,9);
	s.add<AddSyn>(35.0).set(6.2,1245,0.01,0.5,0.0001,6.1,0.99,0.4,0.0005,6.2,0.99,0.3,0.0005,6.2,0.9,1,3,4.07,.56,.92,1.19,1.7,2.74,3.36);
	// s.add<Chimes>(0).freq(rnd::uni(.2,.5), rnd::uni(50.0f,500.0f));
	fillTime(s,40,44, .1, .5, 1000, 2000);
	// fillTime(s,4,7, .5, 1.5, 200, 500);
	fillTimeWith12TET(s,45,58,.100,.400);
	s.add<AddSyn>(59.0).freq(6.2, 155.6); 
	s.add<AddSyn>(65.0).freq(6.2, randomFromMyScale()); 
	s.add<AddSyn>(71.0).freq(6.2,311.1).freqUp(2., 3., 4., 5.).decayLow(456);
	s.add<Chimes>(77.0).freq(6.2,1245);
	

	// A bunch of random-freq notes over 10 seconds: 
	//for (float t = 0; t<10; t += 0.5) {
	//	s.add<Chimes>(t).dur(0.4).freq(randomFromHarmonicSeries());
	//}

	
	AudioIO io(256, 44100., Scheduler::audioCB, &s);
	gam::sampleRate(io.fps());
	io.start();
	printf("\nPress 'enter' to quit...\n"); getchar();
}
