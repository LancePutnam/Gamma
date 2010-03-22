/*
 *  testRandom.cpp
 *  dsp
 *
 *  Created by Lance Putnam on 10/15/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */

#include <stdio.h>
#include "Gamma/arr.h"
#include "Gamma/gen.h"
#include "Gamma/rnd.h"
#include "Gamma/scl.h"
#include "Gamma/Access.h"
#include "Gamma/Print.h"

using namespace gam;

int main(int argc, char* argv[]){
	
	const int iterations = 24;
	#define LOOP for(int i=0; i<iterations; i++)

	float floats[iterations];
	int ints[iterations];

	RNGLinCon lc;
	RNGMulLinCon mlc;
	RNGTaus tw;

	printf("\n\nRNG LinCon float:\n");
		LOOP printf("%4.2f ", rnd::uni_float(lc));

	printf("\n\nRNG MulLinCon float:\n");
		LOOP printf("%4.2f ", rnd::uni_float(mlc));

	printf("\n\nRNG Tausworthe float:\n");
		LOOP printf("%4.2f ", rnd::uni_float(tw));
		
	printf("\n\nRNG shared float:\n");
		LOOP printf("%4.2f ", rnd::uni(1.f));

	printf("\n\npermute:\n");
		for(int i=0; i<16; i++){
			slice(ints,iterations) = gen::RAdd1<int>();
			rnd::permute(ints, iterations);
			LOOP printf("%2i ", ints[i]);
			printf("\n");
		}

	printf("\n\ncond(0.2, 0.2):\n");
		char v = '|';
		for(int i=0; i<256; ++i) printf("%c", rnd::cond(v, '.', '|', 0.2, 0.2));

	printf("\n\ncond(0.8, 0.8):\n");
		for(int i=0; i<256; ++i) printf("%c", rnd::cond(v, '.', '|', 0.8, 0.8));

	printf("\n\nprob(0.2):\n");
		for(int i=0; i<256; ++i) printf("%c", rnd::prob(0.2) ? '|' : '.');
	
	printf("\n\nprob(0.8):\n");
		for(int i=0; i<256; ++i) printf("%c", rnd::prob(0.8) ? '|' : '.');

	printf("\n\nthin(0.2):\n");
		slice(floats, iterations) = gen::Val<>(1);
		rnd::thin(floats, iterations, 0.2);
		LOOP printf("%c", floats[i] != 0.f ? '|' : '.');

	printf("\n\nthin(0.8):\n");
		slice(floats, iterations) = gen::Val<>(1);
		rnd::thin(floats, iterations, 0.8);
		LOOP printf("%c", floats[i] != 0.f ? '|' : '.');

	printf("\n\nuni [0.f, 1.f)\n");
		rnd::uni(floats, iterations);
		LOOP printf("%4.2f ", floats[i]);
	
	printf("\n\nuni [-1.f, 1.f)\n");
		rnd::uni(floats, iterations, -1.f, 1.f);
		LOOP printf("%4.2f ", floats[i]);

	printf("\n\nuni (10, 1)\n");
		rnd::uni(ints, iterations, 10, 1);
		LOOP printf("%i ", ints[i]);

	printf("\n\nuni (-2, 2)\n");
		rnd::uni(ints, iterations, -2, 2);
		LOOP printf("%i ", ints[i]);

	printf("\n\nuniS (1, 4)\n");
		rnd::uniS(ints, iterations, 4);
		LOOP printf("%i ", ints[i]);

	printf("\n\npow3 (0, 10)\n");
		rnd::pow3(ints, iterations, 10);
		LOOP printf("%i ", ints[i]);

	printf("\n\nweighted [0.7, 0.1, 0.1, 0.1]\n");
		float weights[4] = {0.7, 0.1, 0.1, 0.1};
		LOOP printf("%d ", rnd::weighted(weights, 4));

	printf("\n\nPDFs:\n");
	const uint32_t len = 8192 * 8;
	const uint32_t numBins = 16;
	float samples[len];
	float interval[2] = {-1.f, 1.f};
	float normFactor = 1.f;
	uint32_t bins[numBins];
	
	#define HIST(fnc)\
		printf("\n%s\n", #fnc);\
		rnd::fnc(samples, len, interval[1], interval[0]);\
		mem::zero(bins, numBins);\
		arr::histogram(samples, len, bins, numBins, (float)numBins/2, (float)numBins/2);\
		normFactor = (float)len/(float)bins[arr::max(bins, numBins)];\
		for(unsigned long i=0; i<numBins; i++){\
			float amt = (float)bins[i] / (float)len;\
			printf("[% 4.2f] %5.3f ", scl::mapLin((float)i, 0.f, (float)numBins, interval[0], interval[1]), amt);\
			printPlot(amt * normFactor, 30); printf("\n");\
		}\
	
	HIST(add2) HIST(add2I) HIST(add3) HIST(lin) HIST(mul2) HIST(pow2) HIST(pow3) HIST(uni) //HIST(uniS)

}

