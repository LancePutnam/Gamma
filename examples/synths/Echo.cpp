/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
    By Professor Phil Conrad 
*/

#include "examples.h"

int globalCounter = 0;


struct EchoEffect : public Process<AudioIOData>{
  EchoEffect(): echo(0.4, 0.323, 0, 0.8){}
  //Echo(const Process& p): Process(p), echo(0.4, 0.323, 0, 0.8){}
	
	void onProcess(AudioIOData& io){
    while(io()){
      float2 s = float2(io.out(0), io.out(1));
      s = echo(s)*0.5;
      io.out(0) += s[0];
      io.out(1) += s[1];
    }
  }
  Comb<float2> echo;
};

/// Dual delay-line chorus driven by quadrature sinusoid

struct ChorusEffect : public Process<AudioIOData> {
  ChorusEffect(): chorus(0.0021, 0.001, 1, 0.9, 0.1){}
 // ChorusEffect(const Process& p): Process(p), chorus(0.0021, 0.001, 1, 0.9, 0.1){}
	
	void onProcess(AudioIOData& io){
    while(io()){
      float2 s = float2(io.out(0), io.out(1));
      s = chorus(s)*0.5;
      io.out(0) += s[0];
      io.out(1) += s[1];
    }
  }
  Chorus<float2> chorus;
};
	

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


// #include "wellTempered.h"

const double hs = pow(2.0,1/12.0);
const double doh = 220.0;
const double re = doh * hs * hs;
const double me = re * hs;
const double mi = re * hs * hs;
const double fa = mi * hs;
const double sol = fa * hs * hs;
const double le = sol * hs;	  
const double ti = doh / hs;

int main(){

	Scheduler s;

	ProcessNode& effects = s.add<ProcessNode>(); // add effects group
	//s.add<Echo>(effects); // add Echo to effects
	//s.add<Chorus1>(effects); // add chorus to effects

	s.add<ChorusEffect>(effects); // add chorus to effects
	s.add<EchoEffect>(effects); // add Echo to effects

	double attack = 0.1;
	double decay = 0.1;
	double pan = 0.0; 

	double clip = 0.05;

	double q=0.25 - clip;
	// double e=0.125 - clip;
	double h=0.5 - clip;
	//  double w=1.0 - clip;

	double start = 0.0;
	double aug = 1.0;
	double dyn = 0.2;

	s.add<SineEnv>( start + 0     ).set( aug * q, doh, dyn, attack, decay, -1.0*pan);
	s.add<SineEnv>( start + (0.5*aug)  ).set( aug * q, me, dyn, attack, decay, -0.8*pan);
	s.add<SineEnv>( start + (1.0*aug)   ).set( aug * q, sol, dyn, attack, decay, -0.6*pan);
	s.add<SineEnv>( start + (1.5*aug)  ).set( aug * q, le, dyn, attack, decay, -0.4*pan);
	s.add<SineEnv>( start + (2.00*aug)  ).set( aug * h+q, ti, dyn, attack, decay, -0.2*pan);


	AudioIO io(256, 44100., Scheduler::audioCB, &s);
	gam::sampleRate(io.fps());
	io.start();
	printf("\nPress 'enter' to quit...\n"); getchar();
}
