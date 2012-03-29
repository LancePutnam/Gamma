/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */
 
#include <string.h>
#include <math.h>
#include <stdio.h>

#include "Gamma/Constants.h"
#include "Gamma/arr.h"
#include "Gamma/ipl.h"
#include "Gamma/mem.h"
#include "Gamma/tbl.h"

#define LOOP(n,s) for(uint32_t i=0; i<n; i+=s)

namespace gam{
namespace arr{

void linToDB(float * arr, uint32_t len, float minDB){
	float normFactor = 20.f / minDB;
	
	LOOP(len,1){
		uint32_t * arrI = (uint32_t *)arr;
		uint32_t sign = (*arrI) & MaskSign<float>();
		
		float val = fabs(*arr);
		
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

void clip1(float * arr, uint32_t len, uint32_t str){
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

void compact(float * dst, const float * src, uint32_t len, uint32_t chunkSize){

	if(chunkSize < 2){
		mem::deepCopy(dst, src, len);
		return;
	}
	else if(chunkSize > len){
		chunkSize = len;
	}

//	for(uint32_t i=0; i<len; i+=chunkSize){
//		uint32_t min, max;
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
	for(uint32_t i=0; i<len; i+=chunkSize){
		uint32_t max;
		max = indexOfMaxNorm(src, chunkSize);

		*dst++ = src[max];
		src += chunkSize;
	}
}


//uint32_t zeroCross(const float * src, uint32_t len, float prevVal){
//	
//	uint32_t count = 0;
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

uint32_t zeroCross(const float * src, uint32_t len, float prevVal){
	uint32_t count = 0;
	float prev = prevVal;
	LOOP(len,1){
		float curr = *src++;
		if((curr > 0.f && prev <= 0.f) || (curr < 0.f && prev >= 0.f)) count++;
		prev = curr;
	}
	return count;
}

uint32_t zeroCrossFirst(const float * src, uint32_t len){
	uint32_t * srcI = (uint32_t *)src;
	uint32_t prev = *srcI++;

	for(uint32_t i=0; i<len-1; ++i){
		uint32_t now = *srcI++;
		if((now ^ prev)>>31){
			return i-1;
		}
		prev = now;
	}
	return 0;
}

uint32_t zeroCrossN(const float * src, uint32_t len, float prevVal){
	
	uint32_t count = 0;
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
