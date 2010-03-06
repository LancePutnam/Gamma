/*
Test of the conversion functions.
*/

#include <stdio.h>
#include "Gamma/arr.h"
#include "Gamma/ipl.h"
#include "Gamma/scl.h"

using namespace gam;

inline float map_lin_2_exp (float val, float k) { 
	const float a = (k - 1) / (k + 1); 
	const float b = (4 - 2*k) / (k + 1);	// 1 - 3*a
	const float c = 2*a;
	val = (a * val + b) * val + c;
	return (val);
}

inline float fast_log2 (float val)
{
	// assert (val > 0);

	int * const  exp_ptr = reinterpret_cast <int *> (&val);
	int          x = *exp_ptr;
	const int    log_2 = ((x >> 23) & 255) - 128;
	x &= ~(255 << 23);
	x += 127 << 23;
	*exp_ptr = x;

	val = map_lin_2_exp(val, 0.5f);

	return (val + log_2);
}

int main(int argc, char* argv[]){
	
//	arr::conversionInit();
//	
//	printf("2 10987654 32109876543210987654321\n");
//	scl::printBinary(1.f);
//	printf("\n\n");
//
//	for(float f=-64.f; f<64.f; f+=1.1f){
//		printf("% 6.2f ", f);
//
//		scl::printBinary(f);
//
//		long intPart;
//		float frac;
//		frac = gam::split(f, intPart);
//		printf(" %li % 8.6f", intPart, frac);
//
//		printf(" %d", floatToUInt(f));
//		printf(" %li", floatToInt(f));
//		printf(" %g",  scl::trunc(f));
//
//		printf("\n");
//	}
//
//	printf(" %d\n", floatToUInt(16777218.f));
//	//printf(" %lu\n", normalToIndex(1.f, 23));
//
//
//	uint32_t ul = 0x80000000;
//	int32_t sl = 0x80000000;
//	printf("\n"); scl::printBinary(ul); printf(" [>>31] "); scl::printBinary(ul >> 31);
//	//printf("\n"); scl::printBinary(sl); printf(" [>>31] "); scl::printBinary(sl >> 31);
//	printf("\n");
//
//	printf("\nCast double to int (round)\n");
//	for(double i=-1.5; i<1.5; i+=0.125){
//		printf("%6.3f -> %d\n", i, castIntRound(i));
//	}
//
//	printf("\nCast float to int (round)\n");
//	for(float i=-1.5f; i<1.5f; i+=0.125f){
//		printf("%6.3f -> %d\n", i, castIntRound(i));
//	}
//
//	printf("\nFloat exponent\n");
//	for(float i=-4.f; i<4.f; i+=0.25f){
//		printf("%6.3f -> %d\n", i, floatExponent(i));
//	}
//	
//	printf("\nFloat mantissa\n");
//	for(float i=-4.f; i<4.f; i+=0.25f){
//		printf("%6.3f -> %.3f\n", i, floatMantissa(i));
//	}
//
//	printf("\nInt16 to unit\n");
//	for(long i=-32768; i<32767; i+=4096){
//		printf("% 6li -> % f\n", i, intToUnit(i));
//	}
//	printf("% 6d -> % f\n", 0x7fff, intToUnit(0x7fff));
//		
//	printf("\nNumber\tTrailing Zeroes\n");
//	for(unsigned long i=0; i<32; i++){
//		unsigned long num = 1<<i;
//		printf("%10lu%4d\n", num, scl::trailingZeroes(num));
//	}
//	
//	printf("\nNormal to UInt\n");
//	for(double i=-1.5; i<1.5; i+=0.125){
//		printf("%6.3f -> %d\n", i, unitToUInt(i));
//		//printf("%6.3f -> %lu\n", i, (unsigned long)(((double)i) * 4294967296.));
//	}
//	
////	printf("\nNormal to UInt (2)\n");
////	for(double i=-1.5; i<1.5; i+=0.125){
////		printf("%6.3f -> %lu\n", i, scl::unitToUInt2(i));
////	}
//
//	printf("\nNormal to UInt (using castIntRound(double))\n");
//	for(double i=-1.5; i<1.5; i+=0.125)
//		printf("%6.3f -> %lu\n", i, (unsigned long)castIntRound(i * 4294967296.));
//	
//	printf("\nNormal to UInt (using castIntRound(float))\n");
//	for(float i=-1.5f; i<1.5f; i+=0.125f)
//		printf("%6.3f -> %lu\n", i, (unsigned long)castIntRound(i * 4294967296.f));
//	
//	long x = (1<<23 - 1);
//	printf("%ld %d\n", x, castIntRound((float)x));
//	
//	printf("\n");
//	ul = 0x1234567;
//	printf("value              = "); scl::printBinary(ul); printf("\n");
//	printf("~value             = "); scl::printBinary(~ul); printf("\n");
//	printf("0xffffffff - value = "); scl::printBinary(0xffffffff - ul); printf("\n");
	//printf("value/0            = "); scl::printBinary(ul/0); printf("\n");
	
	/*
	const unsigned long len = 16;
	float arr1[len], arr2[len];
	
	printf("\nPolar to Rect:\n");
		mem::set(arr1, gen::Val<>(1), Loop(len));
		arr::line(arr2, len, 0.f, (float)M_2PI);
		
		arr::polarToRect(arr1, arr2, len);
		arr::print(arr1, arr2, len);

	printf("\nPolar to Rect (Fast):\n");
		mem::set(arr1, gen::Val<>(1), Loop(len));
		arr::line(arr2, len, 0.f, (float)M_2PI);
		
		arr::polarToRectFast(arr1, arr2, len);
		arr::print(arr1, arr2, len);
	
	scl::printBinary(0.f); printf("\n");
	printf("%d\n", scl::unitToUInt(0.f));
	
	for(float i=-32.f; i<=32.f; i++){
		float v = i/16.f;
		printf("log2(%6.3f) = % 8.3f % 8.3f (actual=% 8.3f)\n", v, scl::log2Fast(v), fast_log2(v), log2(v));
	}
	*/
	printf("\nendian: %d\n", endian);
	
//	int N=5000;
//	for(int i=0; i<N; ++i){
//		float v = float(i)/N;
//		uint32_t m = unitToUInt(v);
//		uint32_t a = v*4294967295UL;
//		if((m-a)!=0) printf("m=%d, a=%d\n", m,a);
//	}
	
	return 0;
}

