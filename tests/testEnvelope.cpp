#include <stdio.h>
#include "Gamma/scl.h"
#include "Gamma/Envelope.h"
#include "Gamma/Noise.h"
#include "Gamma/Print.h"

using namespace gam;

int main(int argc, char* argv[]){

	const int n = 32;
	Sync::master().spu(n);

	AD<> ad(0.2, 0.8, -3, 3);
	Curve<> curve(n, -3);
	Seg<> segLin(1, 1, -1);
	Seg<float, iplSeq::Cosine> segCos(0.25, 0, 0, 1);
	Seg<float, iplSeq::Cubic> segCub(0.25, 0, 0, 1);
	SegExp<> segExp(1, -3, 1, -1);
	
	NoiseWhite<> noiseWhite;
	gen::Nyquist<> nyq;
	
	#define DO(fnc)\
		printf("\n%s:\n", #fnc);\
		for(int i=0; i<n+4; ++i){\
			double v = fnc;\
			printf("[%2d] % 6.4f ", i, v);\
			printPlot(v, 32); printf("\n");\
		}\
	
	DO(ad())
	curve.set(n,-3); DO(curve()) curve.reset();
	curve.set(n, 3); DO(curve()) curve.reset();
	DO(segExp()) segExp = 1;
	DO(segExp())
	DO(segLin()) segLin = 1;
	DO(segLin())
	
	segLin.period(0.25);
	DO(segLin(nyq))
	DO(segLin(noiseWhite))
	DO(segCos(nyq))
	DO(segCos(noiseWhite))
	DO(segCub(nyq))
	DO(segCub(noiseWhite))
	


//	int iterations = 44100 * 8;
//	Decay<double> dd(iterations);
//	Decay<float>  df(iterations);
//
//	printf("num      end(64b) end(32b)  err(64b)    err(32b) avg err/sample (32b)\n");
//	for(int j=1; j<16; j++){
//		int num = 1280000*j;
//		dd.decay60(num); dd.value1();
//		df.decay60(num); df.value1();
//		double avgErr = 0;
//		for(int i=0; i<num; i++){
//			dd();
//			df();
//			avgErr += dd.value() - (double)df.value();
//		}
//		printf("%8i %8.6f %8.6f % 10.3g % 10.3g % f\n", 
//			num, dd.value(), df.value(), dd.value() - 0.001, df.value() - 0.001f, avgErr / (double)num);
//	}

//	float dur = 64.f;
//	Curve<float> curve(dur, -4.f);
//	
//	PRINT(next);
//	
//	curve.reset();
//	curve.scaleOffset(-1.f, 1.f);
//	
//	//PRINT(nextSO);
//
//	printf("Checking curve = 0\n");
//	curve.set(dur, 0.f);
//	
//	PRINT(next);

//	int num = 64;
//	Decay<double> d1(num), d2(num);
//	//printf("%f\n", d1.decay());
//	double mul = 1. / d1.decay();
//	double att = 0.001 * mul;
//
//	for(int i=0; i<num; i++){
//		//double v = 1. - d1();
//		double v = att;
//		printf("%6.4f ", v);
//		SclOp::printPlot(v); printf("\n");
//		att *= mul;
//	}


}

