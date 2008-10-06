#include <stdio.h>
#include <math.h>
#include <float.h>
#include "arr.h"

#define LINE	arr::lineSlope1(arrE, lenE);\
				arr::lineSlope1(arrO, lenO);
#define ONES	mem::set(arrE, gen::Val<>(1), Loop(lenE));\
				mem::set(arrO, gen::Val<>(1), Loop(lenO));
#define ZEROS	mem::zero(arrE, lenE);\
				mem::zero(arrO, lenO);

#define PRINT_BOTH \
	for(i=0; i<lenE; i++){\
		printf("\t[%2i] % 5.2f  % 5.2f\n", i, arrO[i], arrE[i]);\
	} printf("\t[%2i] % 5.2f\n", lenO-1, arrO[lenO-1]);

#define PRINT_EVEN \
	for(i=0; i<lenE; i++){\
		printf("\t[%2i] % 5.2f\n", i, arrE[i]);\
	}
	
#define ASSERT_GUARDS \
	if(arrE[lenE] != 100.f){ printf("**arrE guard modified to % 5.2f\n", arrE[lenE]); return 1;}\
	if(arrO[lenO] != 100.f){ printf("**arrO guard modified to % 5.2f\n", arrO[lenO]); return 1;}

using namespace gam;

int main(int argc, char* argv[]){
	int i=0;

	const uint32_t lenE = 8;
	const uint32_t lenO = lenE + 1;
	const uint32_t lenE2 = lenE>>1;
	const uint32_t lenO2 = lenO>>1;
	
	float arrE[lenE + 1];
	float arrO[lenO + 1];
	float tempE[lenE];
	float tempO[lenO];

	// Set guard elements
	arrE[lenE] = 100.f;
	arrO[lenO] = 100.f;

	printf("\nAdd flipped 1:\n");
	ZEROS
	//mem::set(arrE + 1, lenE, 2.f, 2);
	arr::addFlip(arrE, lenE, 1.f);
	PRINT_EVEN ASSERT_GUARDS

	printf("\nClip to [1, -1]:\n");
	printf("Input:\n");
	arrE[0] = 1.f / FLT_MAX;
	arrE[1] = -1.f / FLT_MAX;
	arrE[2] = sqrt(-1.);
	arrE[3] = -sqrt(-1.);	
	arrE[4] = 0.f;
	arrE[5] = -0.f;
	arrE[6] = 1.f;
	arrE[7] = -1.f;
	PRINT_BOTH
	
	printf("Output:\n");
	arr::clip1(arrE, lenE);
	arr::clip1(arrO, lenO);
	PRINT_BOTH

	printf("\nDifferentiate:\n");
	LINE
	float prev = -1.f; arr::differentiate(arrE, lenE, prev);
	      prev = -1.f; arr::differentiate(arrO, lenO, prev);
	PRINT_BOTH ASSERT_GUARDS

	printf("\nExponentiate base 2:\n");
	LINE
	arr::expBase(arrE, lenE, 2.f);
	arr::expBase(arrO, lenO, 2.f);
	PRINT_BOTH ASSERT_GUARDS

	printf("\nLinear to dB (thresh = -12):\n");
	printf("Output:\n");
	arr::linToDB(arrE, lenE, -12.f);
	PRINT_EVEN

	printf("\nMirror (dp):\n");
	LINE
	arr::mirror_dp(arrE, lenE);
	arr::mirror_dp(arrO, lenO);
	PRINT_BOTH ASSERT_GUARDS

	printf("\nMirror (dq):\n");
	LINE
	arr::mirror_dq(arrE, lenE);
	arr::mirror_dq(arrO, lenO);
	PRINT_BOTH ASSERT_GUARDS

	printf("\nMultiply 1s by Bartlett:\n");
	ONES
	arr::mulBartlett(arrE, lenE);
	arr::mulBartlett(arrO, lenO);
	PRINT_BOTH

	printf("\nNormalize:\n");
	printf("Input:\n");
	arrE[0] = 0.f;
	arrE[1] = 0.1f;
	arrE[2] = -0.1f;
	arrE[3] = 0.2f;	
	arrE[4] = -0.2f;
	arrE[5] = 0.3f;
	arrE[6] = -0.3f;
	arrE[7] = -0.4f;
	PRINT_EVEN
	float norm = arr::normalize(arrE, lenE);
	printf("Output (norm = %f):\n", norm);
	PRINT_EVEN
	
	return 0;
}








