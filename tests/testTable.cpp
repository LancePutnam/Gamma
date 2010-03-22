#include <stdio.h>
#include <string.h>

#include "Gamma/arr.h"
#include "Gamma/gen.h"
#include "Gamma/tbl.h"
#include "Gamma/Access.h"
#include "Gamma/Constants.h"
#include "Gamma/Conversion.h"
#include "Gamma/Print.h"

#define TEST_WINDOWS
#define TEST_TABLEOSC
//#define TEST_CURVES
//#define TEST_ACCESS

#define PRINT_TABLE \
	for(unsigned long i=0; i<len; i++){\
		float v = arr[i];\
		printf("[%4lu] % 5.3f %8lx  ", i, v, (long unsigned int)punFU(v));\
		printPlot(v); printf("\n");\
	}

using namespace gam;

int main(int argc, char* argv[]){
	
	const uint32_t len = 32;
	float arr[len];
	unsigned long numHarmonics = len >> 1;
	
	#ifdef TEST_TABLEOSC

	printf("\n\nCosine...\n");
		mem::zero(arr, len);
		tbl::cosine(arr, len);
		PRINT_TABLE

	printf("\n\nDecay (order=-4)...\n");
		mem::zero(arr, len);
		tbl::decay(arr, len, -4.);
		PRINT_TABLE
	
	printf("\n\nSine...\n");
		mem::zero(arr, len);
		tbl::sine(arr, len);
		PRINT_TABLE

	printf("\n\nSinusoid (phase=pi/2, periods=3)...\n");
		mem::zero(arr, len);
		tbl::sinusoid(arr, len, M_PI_2, 3.);
		PRINT_TABLE
	
	printf("\n\nImpulse...\n");
		mem::zero(arr, len);
		tbl::impulseSum(arr, len, 1, (numHarmonics>>0)-0);
		//ArrOp::normalize(arr, len);
		PRINT_TABLE

	printf("\n\nImpulse (full harmonics)...\n");
		mem::zero(arr, len);
		tbl::impulseSum(arr, len);
		PRINT_TABLE

	printf("\n\nSaw...\n");
		mem::zero(arr, len);
		tbl::sawSum(arr, len, 1, numHarmonics);
		PRINT_TABLE

	printf("\n\nSquare...\n");
		mem::zero(arr, len);
		tbl::squareSum(arr, len, 1, numHarmonics);
		PRINT_TABLE

	printf("\n\nTriangle...\n");
		mem::zero(arr, len);
		tbl::triangleSum(arr, len, 1, numHarmonics);
		PRINT_TABLE
	
	#endif
	
	
	#ifdef TEST_WINDOWS
	
	#define WIN(name) printf("\n\n"#name":\n"); tbl::name(arr, len); PRINT_TABLE
	
	WIN(blackman) WIN(blackmanHarris) WIN(hamming) WIN(hann) WIN(nyquist)
	WIN(bartlett) WIN(welch)
	
	#endif
	
	
	#ifdef TEST_CURVES
	printf("\n\nDecay...\n");
		tbl::decay(arr, len, -4.f);
		PRINT_TABLE
		
	#endif
	
	#ifdef TEST_ACCESS

	#define Q_NUM_BITS 3
	float qSine[1<<Q_NUM_BITS + 1];
	unsigned long qLength = 1<<Q_NUM_BITS;

	tbl::sinusoid(qSine, qLength, 0.f, 0.25f);
	qSine[qLength] = 1.f;	// set the end sample
	printf("\n\nQuarter sine table values...\n");
	ArrOp::printHex(qSine, qLength+1);
	
	unsigned long bits = 30 - Q_NUM_BITS;
	unsigned long phaseInc = 1 << bits;
	
	printf("\nQuarter sine table full period access...\n");
	for(unsigned int i=0; i<(qLength*4); i++){
		float value = tbl::atQ(qSine, bits, i * phaseInc);
		printf("[%3i] % 5.3f %x\n", i, value, value);
	}

	#define H_NUM_BITS 4
	float hSine[1<<H_NUM_BITS];
	unsigned long hLength = 1<<H_NUM_BITS;

	tbl::sinusoid(hSine, hLength, 0.f, 0.5f);
	printf("\n\nHalf sine table values...\n");
	ArrOp::printHex(hSine, hLength);
	
	bits = 31 - H_NUM_BITS;
	phaseInc = 1 << bits;
	
	printf("\nHalf sine table full period access...\n");
	for(unsigned int i=0; i<(hLength*2); i++){
		float value = tbl::atH(hSine, bits, i * phaseInc);
		printf("[%3i] % 5.3f %x\n", i, value, value);
	}	

//	unsigned int dirBit = 0x40000000;
//	//unsigned int dirBit = 0x00000000;
//	dirBit >>= 30;
//	for(unsigned int i=0; i<32; i++){
//	//	printf("%u\t%u\n", i, i | (dirBit<<1));
//	//	printf("%u\t%u\n", i, (-i) & 31U);
//	//	printf("%u\t%u\n", i, (~i + 1) & 31U);
//	//	printf("%u\t%u\n", i, ((i^0xffffffff) + 1) & 31U);
//	//	printf("%u\t%u\n", i, ((i^0xffffffff) + (dirBit>>30)) & 31U);
//		printf("%u\t%u\n", i, ((i^(-dirBit)) + dirBit) & 31U);
//	}

	// 1 -> 0xffffffff	0001 -> 1111
	// 0 -> 0x00000000	0000 -> 0000

	printf("\n\nqSine as C hex array:\n");
	tbl::printHexArray(qSine, qLength+1, 4);

	printf("\n\nhSine as C hex array:\n");
	tbl::printHexArray(hSine, hLength, 4);

	#endif

//	for(unsigned long i=0; i<32; i++){
//		printf("[%3i] %f\n", i, tbl::fraction(i, 1<<(32-i-1)));
//	}
	{ using namespace gam::gen;
		
		#define GEN(fnc) slice(arr,len) = fnc; printf("\n\n"#fnc":\n"); PRINT_TABLE
		
		GEN(Val<>(1))
		GEN(RMul<>(0.9))
		GEN(RMul<>(1.1, 0.04))
		GEN(RAdd<>(1./(double)len))
		GEN(RAdd<>(-1./(double)len, 1.))
		GEN(Recip<>(1))
		GEN(Sin<>(M_2PI/(double)len, 0))
		GEN(Sin<>(M_2PI/(double)len, M_PI_2))
	}
	
	return 0;
}

