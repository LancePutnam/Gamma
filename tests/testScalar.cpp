#include <stdio.h>
#include "Gamma/rnd.h"
#include "Gamma/scl.h"
#include "Gamma/Print.h"

using namespace gam;

int main(int argc, char* argv[]){

	const int len = 12;
	float array[len];

	array[0] = INFINITY;
	array[1] = -INFINITY;
	array[2] = NAN;
	array[3] = -NAN;
	array[4] = 2.;
	array[5] = 2.25;
	array[6] = 2.75;
	array[7] = 3.25;
	array[8] = 1.;
	array[9] = 0.75;
	array[10] = 0.25;
	array[11] = -0.25;


	float f1, f2, f3, f4;
	double d1, d2;

	#define LOOP for(int i=0; i<len; i++)

	printf("\nbase36To10:\n");
	for(char c='0'; c<='Z'; c++) printf("%2c ", c); printf("\n");
	for(char c='0'; c<='Z'; c++) printf("%2d ", base36To10(c)); printf("\n");

	printf("\nbase10To36:\n");
	for(int i=0; i<36; i++) printf("%2d ", i); printf("\n");
	for(int i=0; i<36; i++) printf("%2c ", base10To36(i)); printf("\n");

	printf("\nCeiling (power of 2):\n");
	for(int i=1; i<33; i++){
		unsigned long val = (1<<i) - i;
		printf("of %lu = %d\n", val, scl::ceilPow2(val));
	}

	printf("\nFloor\n");
		for(double d = -8; d < 8; d+=1.25){
			printf("% 5.2f -> % 5.2f\n", d, scl::floor(d));
		}

	printf("\nFloor (power of 2):\n");
	for(int i=1; i<33; i++){
		uint32_t val = (1<<i) - i;
		printf("of %d = %d\n", val, scl::floorPow2(val));
	}

	printf("\nFold:\n");
		LOOP{
			printf("% 5.2f in [1, 2) -> % 5.2f\n", array[i], scl::fold(array[i], 2.f, 1.f));
		}

	printf("\nLinear Map:\n");
		float scale, offset;
		f1 = -1; f2 = 1; f3 = 0; f4 = 4;
		scl::mapLin(f1, f2, f3, f4, scale, offset);
		printf("[%.2f, %.2f] -> [%.2f, %.2f] = scale: %.2f, offset: %.2f\n", f1, f2, f3, f4, scale, offset);		

	printf("\nLog (base 2):\n");
	for(int i=1; i<33; i++){
		uint32_t val = (1<<i) - i;
		printf("of %d = %d\n", val, scl::log2(val));
	}

	printf("\nMax:\n");
		f1 = 0.f; f2 = 1.f;
		printf("%f, %f = %f\n", f1, f2, scl::max(f1, f2));
		f1 = -1.f; f2 = 0.f;
		printf("%f, %f = %f\n", f1, f2, scl::max(f1, f2));
	
		d1 = 0.; d2 = 1.;
		printf("%f, %f = %f\n", d1, d2, scl::max(d1, d2));

	printf("\nMin:\n");
		f1 = 0.f; f2 = INFINITY;
		printf("%f, %f = %f\n", f1, f2, scl::min(f1, f2));

	printf("\nMul Complex:\n");
		f1 = 1.f; f2 = 0.f;
		for(int i=0; i<17; i++){
			printf("% 9.6f % 9.6f\n", f1, f2);
			float ang = M_2PI/(float)16;
			scl::mulComplex(f1, f2, cosf(ang), sinf(ang));
		}
		
	printf("\nNearest:\n");
		for(float i=0.f; i<12.f; i+=0.77f){
			printf("% 6.2f -> % 6.2f\n", i, scl::nearest(i));
		}

	printf("\nNext multiple:\n");
		printf("of 2 of %d = %d\n", 3, scl::nextMultiple(3U, 2U));		
	
	printf("\nrecSqrtFast()\n");
		for(int i=0; i<16; ++i){
			float v = i+1;
			printf("% 7.3f % 9.7f % 9.7f\n", v, 1./sqrt(v), scl::invSqrt<1>(v));
		}
	
	printf("\nRound (0.5)\n");
		LOOP{
			printf("% 5.2f -> % 5.2f\n", array[i], scl::round(array[i], 0.5f));
		}

	printf("\nRound1\n");
		for(double d = -8; d < 8; d+=1.31){
			printf("% 5.2f -> % 5.2f\n", d, scl::round(d));
		}

	printf("\nSlope:\n");
		f1 = 1; f2 = 1; f3 = 2; f4 = 3;
		printf("(%.2f, %.2f), (%.2f, %.2f) = %f\n", f1, f2, f3, f4, scl::slope(f1, f2, f3, f4));

	printf("\nsqrt<1>()\n");
		for(int i=0; i<16; ++i){
			float v = i+1;
			float real = sqrt(v);
			float apx = scl::sqrt<1>(v);
			printf("% 7.3f % 9.7f % 9.7f, error= % 7.4g\n", v, real, apx, scl::error(real,apx)*100);
		}

	printf("\nsqrt<3>()\n");
		for(int i=0; i<16; ++i){
			double v = i+1;
			double real = sqrt(v);
			double apx = scl::sqrt<3>(v);
			printf("% 7.3f % 9.7f % 9.7f, error= % 7.4g\n", v, real, apx, scl::error(real,apx)*100);
		}

	printf("\nTrunc:\n");
		for(double d = -8; d < 8; d+=1.31){
			printf("% 5.2f -> % 5.2f\n", d, scl::trunc(d));
		}
	
	printf("\nWrap:\n");
		LOOP{
			printf("% 5.2f in [1, 2) -> % 5.2f\n", array[i], scl::wrap(array[i], 2.f, 1.f));
		}

	printf("\nWrap:\n");
		LOOP{
			long numWraps;
			float val = scl::wrap(array[i], numWraps, 2.f, 1.f);
			printf("% 5.2f in [1, 2) -> % 5.2f (%li wraps)\n", array[i], val, numWraps);
		}
		
	printf("\nWrap1:\n");
		LOOP{
			printf("% 5.2f in [0, 1) -> % 5.2f\n", array[i], scl::wrap(array[i]));
		}
		
	printf("%f\n", scl::freq("a5"));
	
	for(int i=0; i<25; ++i){
		printf("%2d -> %f\n", i, scl::nearest((float)i, "22122212"));
	}
	
	//float tx[2], ty[2];
	//scl::cross(tx, ty);
	
//	for(int i=0; i<32; ++i){
//		float p = (i/32.)*2-1;
//		printf("% 6.4f: % 6.4f (% 6.4f)\n", p, sinFast(p*M_PI), sin(p*M_PI));
//	}

}

