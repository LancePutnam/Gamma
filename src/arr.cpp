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


static const uint32_t LUTSize = 2048 - 1; // power of two - 1
static const uint32_t LUTMask = LUTSize;
static const uint32_t LUTSize2 = LUTSize >> 1;
static const uint32_t LUTSize4 = LUTSize >> 2;
static const float LUTSize2F = (float)LUTSize2;
static const float phaseToIndex = ((float)LUTSize) * M_1_2PI;	

static float * atanLUT = 0;
static float * magLUT = 0;
static ArrayPow2<float> * sinLUT = 0;


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
		max = maxNorm(src, chunkSize);

		*dst++ = src[max];
		src += chunkSize;
	}
}

uint32_t fundHPS(float * tmp, const float * mag, uint32_t len, uint32_t downSample){
	hps(tmp, mag, len, downSample);
	return max(tmp, len);
}

void hps(float * dst, const float * src, uint32_t len, uint32_t downSample){
	
	memcpy(dst, src, len * sizeof(float));
	
	for(uint32_t d=2; d <= downSample; d++){
		for(uint32_t i=0; i<len; i++){
			dst[i/d] *= src[i];
		}
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


void magFrqToPolar(float * frq, float * phsAccum, uint32_t len, float factorUnwrap){

//	float binFreq = fundFreq;
//	float expPhaseDiff = fundRadians;
//	
//	LOOP(len,
//		float phaseDiff = (*frq - binFreq) * factorUnwrap;
//		float phaseNew = *phsAccum;
//
//		// original DM style
//		//phaseNew += phaseDiff;			// UNWRAP phase differences
//
//		// SB style
//		phaseNew += phaseDiff + expPhaseDiff;			// UNWRAP phase differences
//		
//		// reduced SB
//		//phaseNew += *frq * factorUnwrap;
//		
//		phaseNew = scl::wrapPhase(phaseNew);		// wrap to [-pi, pi) to retain precision
//		//phaseNew = scl::wrap(phaseNew, 128.f * (float)M_PI, -128.f * (float)M_PI);
//		*frq++ = *phsAccum++ = phaseNew;
//		
//		binFreq += fundFreq;
//		expPhaseDiff += fundRadians;
//	)
	
	LOOP(len,1){
		float phaseNew = *phsAccum + *frq * factorUnwrap;
		phaseNew = scl::wrapPhase(phaseNew);		// wrap to [-pi, pi) to retain precision
		*frq++ = *phsAccum++ = phaseNew;
	}
	
}


/*
factorWrap  := spu / (sizeHop * M_2PI)
fundFreq    := spu / sizeDFT
fundRadians := M_2PI * sizeHop / sizeDFT

unitsHop := sizeHop * ups;
*/
void polarToMagFrq(float * p0, float * p1, uint32_t len, float factorWrap, float fundFreq, float fundRadians){
	
	float binFreq = fundFreq;			// center frequency of bin
	float expPhaseDiff = fundRadians;	// expected phase difference based on bin number
	
	LOOP(len,1){
				
		// Wrap phase difference into range [-pi, pi)		
		float dp = scl::wrapPhase(*p0 - *p1 - expPhaseDiff);	// SB
		//float dp = scl::wrapPhase(*p0 - *p1);					// DM

		*p1++ = *p0;						// Store this frame's phase for next frame
		*p0++ = dp * factorWrap + binFreq;	// Compute instantaneous frequency
		
		binFreq += fundFreq;
		expPhaseDiff += fundRadians;
	}
}

//TEM void phaseToFreq(T * p0, T * p1, uint32_t len, T ups){
//	T factor = (T)1 / (M_2PI * ups);
//	LOOP_P(len,
//		T dp = scl::wrapPhase(*p0 - *p1);		// wrap phase into [-pi, pi)
//		*p1++ = *p0;							// prev phase = curr phase
//		*p0++ = dp * factor;
//	)
//}

void polarToRect(float * mag, float * phs, uint32_t len){
	LOOP(len,1){
		float m = *mag;
		float p = *phs;
		*mag++ = m * cos(p);
		*phs++ = m * sin(p);
	}
}


void polarToRectFast(float * magA, float * phsA, uint32_t len){
	
	uint32_t bits = sinLUT->log2Size();
	uint32_t fracBits = sinLUT->fracBits();
	uint32_t quarter = 0x40000000;
	uint32_t oneIndex = sinLUT->oneIndex();
	float * sinTable = sinLUT->elems();

	LOOP(len,1){
		float mag = *magA;
		float phs = *phsA;
		
		phs = scl::wrap(phs, (float)M_2PI);
		phs *= (float)M_1_2PI;
		
		//uint32_t index = scl::unitToUInt(phs);
		uint32_t index = unitToUInt2(phs);	// better cuz phase is on a linear scale
		
//		float sn = MemOpF::at(sinTable, fracBits, index);
//		float cs = MemOpF::at(sinTable, fracBits, index + quarter);

		float frac = fraction(bits, index);
		float sn = ipl::linear(frac,
			mem::at(sinTable, fracBits, index),
			mem::at(sinTable, fracBits, index + oneIndex)
		);
		
		index += quarter;
		float cs = ipl::linear(frac,
			mem::at(sinTable, fracBits, index),
			mem::at(sinTable, fracBits, index + oneIndex)
		);
		
		*magA++ = mag * cs;
		*phsA++ = mag * sn;
	}
}


void rectToPolar(float * re, float * im, uint32_t len, uint32_t str){
	LOOP(len,str){ scl::rectToPolar(re[i], im[i]); }
}


void rectToPolarFast(float * realA, float * imagA, uint32_t len){
	
	LOOP(len, 1){
		float real = *realA;
		float imag = *imagA;

		float realAbs = fabs(real);
		float imagAbs = fabs(imag);
		
		if(realAbs > imagAbs){
			float ratio = imag / real;	// (-1, 1)
			//float ratio = 0.5f;
			
			uint32_t index = (uint32_t)(LUTSize2F + LUTSize2F * ratio);
			//uint32_t index = 10;
			//uint32_t index = normalToIndex(ratio * 0.5f + 0.5f, 11);
			//int index = roundFloatToInt(LUTSize2F + LUTSize2F * ratio);
			
			if(real > 0.f)  imag = atanLUT[index];				// -pi/4, pi/4
			else			imag = M_PI + atanLUT[index];		// 3pi/4, 5pi/4
			real = realAbs * magLUT[index];
		}
		
		else if((imagAbs == 0.f) && (realAbs == 0.f)){
			real = 0.f;
			imag = 0.f;
		}
		
		else{
			float ratio = real / imag;
			//float ratio = 0.5f;
			
			uint32_t index = (uint32_t)(LUTSize2F + LUTSize2F * ratio);
			//uint32_t index = 10;
			//uint32_t index = normalToIndex(ratio * 0.5f + 0.5f, 11);
			//int index = roundFloatToInt(LUTSize2F + LUTSize2F * ratio);

			if(imag > 0.f)  imag = M_PI_2 - atanLUT[index];		// pi/4, 3pi/4
			else			imag = M_3PI_2 - atanLUT[index];	// 5pi/4, 7pi/4
			real = imagAbs * magLUT[index];
		}
		
		*realA++ = real;
		*imagA++ = imag;
	}
}



void conversionInit(){


	
//	ArrayPow2<float> sinLUT(11);
//	float atanLUT[LUTSize + 1], magLUT[LUTSize + 1];

//	ArrayPow2<float> * sinLUT;
//	float * atanLUT, * magLUT;
//	const uint32_t LUTSize;
//	
//	const uint32_t LUTMask;
//	const uint32_t LUTSize2;
//	const uint32_t LUTSize4;
//	const float LUTSize2F;
//	const float phaseToIndex;

	if(0 == atanLUT){
		atanLUT = new float[LUTSize + 1];
		magLUT  = new float[LUTSize + 1];
		sinLUT  = new ArrayPow2<float>(2048);
		
		//printf("arr:conversionInit(): %d\n", sinLUT->size());

		tbl::sine(sinLUT->elems(), sinLUT->size());
			
		double LUTSizeRec = 1. / (double)LUTSize;
		double ramp2 = -1.;
		double ramp2Inc = 2. * LUTSizeRec;
		
		for(uint32_t i=0; i<LUTSize; ++i){
			
			double angle = atan(ramp2);
			atanLUT[i] = angle;
			magLUT[i] = 1. / cos(angle);
			
			ramp2 += ramp2Inc;
		}
	}
}

} // arr::
} // gam::

#undef LOOP
