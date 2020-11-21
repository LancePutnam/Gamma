/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include <cmath>
#include "Gamma/arr.h"
#include "Gamma/mem.h" // deepCopy
#include "Gamma/Constants.h"

#define LOOP(n,s) for(unsigned i=0; i<n; i+=s)

namespace gam{
namespace arr{

void linToDB(float * arr, unsigned len, float minDB){
	float normFactor = 20.f / minDB;

	LOOP(len,1){
		uint32_t * arrI = (uint32_t *)arr;
		uint32_t sign = (*arrI) & MaskSign<float>();
		
		float val = std::fabs(*arr);
		
		if(val == 0.f){
			*arr++ = 0.f;
		}
		else{
			val = 1.f - normFactor * log10(val);
			if(val < 0.f) val = 0.f;
			*arr++ = val;
			*arrI |= sign;
		}
	}
}

//inline T scl::linToDB(T v){ return (T)log10(v) * (T)20; }
//inline T scl::dBToLin(T v){ return pow(10., v * 0.05); }

void clip1(float * arr, unsigned len, unsigned str){
	//uint32_t * arrUI = (uint32_t *)arr;
	LOOP(len, str){
		uint32_t val = punFU(arr[i]);
		uint32_t sign = val & MaskSign<float>();
		val &= 0x7fffffff;	// mask off the sign to catch infs and NaNs
		uint32_t above = (val + 0x00800000) >> 30;
		if(above) val = Expo1<float>();
		arr[i] = punUF(val | sign);
	}
}

void compact(float * dst, const float * src, unsigned len, unsigned chunkSize){

	if(chunkSize < 2){
		mem::deepCopy(dst, src, len);
		return;
	}
	else if(chunkSize > len){
		chunkSize = len;
	}

//	for(unsigned i=0; i<len; i+=chunkSize){
//		unsigned min, max;
//		extrema(src, chunkSize, &min, &max);
//		
//		if(min < max){
//			*dst++ = src[min];
//			*dst++ = src[max];
//		}
//		else{
//			*dst++ = src[max];
//			*dst++ = src[min];
//		}
//
//		src += chunkSize;
//	}
	for(unsigned i=0; i<len; i+=chunkSize){
		unsigned max = indexOfMaxNorm(src, chunkSize);
		*dst++ = src[max];
		src += chunkSize;
	}
}


//unsigned zeroCross(const float * src, unsigned len, float prevVal){
//	
//	unsigned count = 0;
//	uint32_t * srcI = (uint32_t *)src;
//	uint32_t prev = *(uint32_t *)(&prevVal);
//
//	LOOP(len,
//		uint32_t now = *srcI++;
//
//		if((now > 0 && prev <= 0) || (now < 0 && prev >= 0)) count++;
//		//count += (now ^ prev)>>31;
//
//		prev = now;
//	)
//
//	return count;
//}

unsigned zeroCross(const float * src, unsigned len, float prevVal){
	unsigned count = 0;
	float prev = prevVal;
	LOOP(len,1){
		float curr = *src++;
		count += unsigned((curr > 0.f && prev <= 0.f) || (curr < 0.f && prev >= 0.f));
		prev = curr;
	}
	return count;
}

unsigned zeroCrossFirst(const float * src, unsigned len){
	uint32_t * srcI = (uint32_t *)src;
	uint32_t prev = *srcI++;

	for(unsigned i=0; i<len-1; ++i){
		uint32_t now = *srcI++;
		if((now ^ prev)>>31){
			return i-1;
		}
		prev = now;
	}
	return 0;
}

unsigned zeroCrossN(const float * src, unsigned len, float prevVal){
	
	unsigned count = 0;
	uint32_t * srcI = (uint32_t *)src;
	uint32_t prev = Twiddle<float>(prevVal).u;

	LOOP(len,1){
		uint32_t now = *srcI++;

		//if((now > 0) && (prev <= 0))		count++;
		//else if((now < 0) && (prev >= 0))	count++;

		//if((now > 0) && (prev <= 0))		count++;

		count += ((now & ~prev)>>31);

		prev = now;
	}

	return count;
}

} // arr::
} // gam::

#undef LOOP
