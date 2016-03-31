/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Example:		Filter / Plucked String
	Description:	Simulation of a plucked string with noise and a feedback 
					delay-line.
*/

#include "examples.h"

class PluckedString : public Process<AudioIOData> {
public:

	PluckedString(double startTime=0, float frq=440)
	:	mAmp(1), mDur(2), delay1(0.4,	0.2),
		env(0.1), fil(2), delay(1./27.5, 1./frq)
	{
		dt(startTime);
		decay(1.0);
		mAmpEnv.curve(4); // make segments lines
		mAmpEnv.levels(1,1,0);
	}
		
	PluckedString& freq(float v){delay.freq(v); return *this; }
	PluckedString& amp(float v){ mAmp=v; return *this; }
	PluckedString& dur(float v){ 
		mAmpEnv.lengths()[0] = v;
		return *this; }
	PluckedString& decay(float v){
		mAmpEnv.lengths()[1] = v;
		return *this;
	}
	
	PluckedString& pan(float v){ mPan.pos(v); return *this; }
	void reset(){ env.reset(); }
	
	PluckedString& set(
		float a, float b, float c, float d, float e=0
	
	){
		return dur(a).freq(b).amp(c).decay(d).pan(e);
	}

	float operator() (){
		return (*this)(noise()*env());
	}
	
	float operator() (float in){
		return delay(
					 fil( delay() + in )
					 );
	}
	
	void onProcess(AudioIOData& io){
	
		while(io()){
			float s =  (*this)() * mAmpEnv() * mAmp;
			// This short-hand method is convenient for simple delays.
			//float s1 = s += delay1(s);
		
			// We can create infinite echoes by feeding back a small amount of the 
			// output back into the input on each iteration.
			float s1 = s += delay1(s + delay1()*0.2);
		
			// We can also create mult-tap delay-lines through multiple calls to 
			// the read() method.
			//float s1 = s += delay1(s) + delay1.read(0.15) + delay1.read(0.39);
	
			// How about multi-tap feedback?
			//float s1 = s += delay1(s + delay.read(0.197)*0.3 + delay1.read(0.141)*0.4 + delay1.read(0.093)*0.2)*0.5;
			
			float s2;
			mEnvFollow(s1);
			mPan(s1, s1,s2);
			io.out(0) += s1;
			io.out(1) += s2;
		}
		if(mAmpEnv.done() && (mEnvFollow.value() < 0.001)) free();
	}

protected:
	float mAmp;
	float mDur;
	Pan<> mPan;
	NoiseWhite<> noise;
	Decay<> env;
	MovingAvg<> fil;
	Delay<float, ipl::Trunc> delay, delay1;
	Env<2> mAmpEnv;
	EnvFollow<> mEnvFollow;
};

int main(){

	Scheduler s;
	s.add<PluckedString>( 0  ).set(6.5, 110,  0.3, .005, -1);
	s.add<PluckedString>( 3.5).set(6.5, 233,  0.3, .1, 0);
	PluckedString &thirdPluck = s.add<PluckedString>( 6.5).set(6.5, 329,  0.7, .0001, 1);
	s.add(Func(thirdPluck, &PluckedString::freq, 440)).dt(8);
	
	AudioIO io(256, 44100., Scheduler::audioCB, &s);
	gam::sampleRate(io.fps());
	io.start();
	printf("\nPress 'enter' to quit...\n"); getchar();
}
