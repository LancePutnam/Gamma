/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */
 
#include <string.h>
#include <math.h>
#include <stdio.h>
#include "Constants.h"

#include "arr.h"
#include "ipl.h"
#include "mem.h"
#include "tbl.h"

#include "MacroD.h"

namespace gam{
namespace arr{

void linToDB(float * arr, ULONG len, float minDB){
	float normFactor = 20.f / minDB;
	
	LOOP(len,
		ULONG * arrI = (ULONG *)arr;
		ULONG sign = (*arrI) & MASK_F32_SIGN;
		
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
	)
}

//inline T scl::linToDB(T v){ return (T)log10(v) * (T)20; }
//inline T scl::dBToLin(T v){ return pow(10., v * 0.05); }

void compact(float * dst, const float * src, ULONG len, ULONG chunkSize){

	if(chunkSize < 2){
		mem::copy(dst, src, len);
		return;
	}
	else if(chunkSize > len){
		chunkSize = len;
	}

//	for(ULONG i=0; i<len; i+=chunkSize){
//		ULONG min, max;
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
	for(ULONG i=0; i<len; i+=chunkSize){
		ULONG max;
		max = maxAbs(src, chunkSize);

		*dst++ = src[max];
		src += chunkSize;
	}
}

ULONG fundHPS(float * tmp, const float * mag, ULONG len, ULONG downSample){
	hps(tmp, mag, len, downSample);
	return max(tmp, len);
}

void hps(float * dst, const float * src, ULONG len, ULONG downSample){
	
	memcpy(dst, src, len * sizeof(float));
	
	for(ULONG d=2; d <= downSample; d++){
		for(ULONG i=0; i<len; i++){
			dst[i/d] *= src[i];
		}
	}
}






//ULONG zeroCross(const float * src, ULONG len, float prevVal){
//	
//	ULONG count = 0;
//	ULONG * srcI = (ULONG *)src;
//	ULONG prev = *(ULONG *)(&prevVal);
//
//	LOOP(len,
//		ULONG now = *srcI++;
//
//		if((now > 0 && prev <= 0) || (now < 0 && prev >= 0)) count++;
//		//count += (now ^ prev)>>31;
//
//		prev = now;
//	)
//
//	return count;
//}

ULONG zeroCross(const float * src, ULONG len, float prevVal){
	ULONG count = 0;
	float prev = prevVal;
	LOOP(len,
		float curr = *src++;
		if((curr > 0.f && prev <= 0.f) || (curr < 0.f && prev >= 0.f)) count++;
		prev = curr;
	)
	return count;
}

ULONG zeroCrossFirst(const float * src, ULONG len){
	ULONG * srcI = (ULONG *)src;
	ULONG prev = *srcI++;

	for(ULONG i=0; i<len-1; ++i){
		ULONG now = *srcI++;
		if((now ^ prev)>>31){
			return i-1;
		}
		prev = now;
	}
	return 0;
}

ULONG zeroCrossN(const float * src, ULONG len, float prevVal){
	
	ULONG count = 0;
	ULONG * srcI = (ULONG *)src;
	ULONG prev = *(ULONG *)(&prevVal);

	LOOP(len,
		ULONG now = *srcI++;

		//if((now > 0) && (prev <= 0))		count++;
		//else if((now < 0) && (prev >= 0))	count++;

		//if((now > 0) && (prev <= 0))		count++;

		count += ((now & ~prev)>>31);

		prev = now;
	)

	return count;
}


void magFrqToPolar(float * frq, float * phsAccum, ULONG len, float factorUnwrap){

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
	
	LOOP(len,
		float phaseNew = *phsAccum + *frq * factorUnwrap;
		phaseNew = scl::wrapPhase(phaseNew);		// wrap to [-pi, pi) to retain precision
		*frq++ = *phsAccum++ = phaseNew;
	)
	
}


/*
factorWrap  := spu / (sizeHop * M_2PI)
fundFreq    := spu / sizeDFT
fundRadians := M_2PI * sizeHop / sizeDFT

unitsHop := sizeHop * ups;
*/
void polarToMagFrq(float * p0, float * p1, ULONG len, float factorWrap, float fundFreq, float fundRadians){
	
	float binFreq = fundFreq;			// center frequency of bin
	float expPhaseDiff = fundRadians;	// expected phase difference based on bin number
	
	LOOP(len,
				
		// Wrap phase difference into range [-pi, pi)		
		float dp = scl::wrapPhase(*p0 - *p1 - expPhaseDiff);	// SB
		//float dp = scl::wrapPhase(*p0 - *p1);					// DM

		*p1++ = *p0;						// Store this frame's phase for next frame
		*p0++ = dp * factorWrap + binFreq;	// Compute instantaneous frequency
		
		binFreq += fundFreq;
		expPhaseDiff += fundRadians;
	)
}

//TEM void phaseToFreq(T * p0, T * p1, ULONG len, T ups){
//	T factor = (T)1 / (M_2PI * ups);
//	LOOP_P(len,
//		T dp = scl::wrapPhase(*p0 - *p1);		// wrap phase into [-pi, pi)
//		*p1++ = *p0;							// prev phase = curr phase
//		*p0++ = dp * factor;
//	)
//}

void polarToRect(float * mag, float * phs, ULONG len){
	LOOP_P(len,
		float m = *mag;
		float p = *phs;
		*mag++ = m * cos(p);
		*phs++ = m * sin(p);
	)
}


void polarToRectFast(float * magA, float * phsA, ULONG len){
	
	ULONG bits = sinLUT->log2Size();
	ULONG fracBits = sinLUT->fracBits();
	ULONG quarter = 0x40000000;
	ULONG oneIndex = sinLUT->oneIndex();
	float * sinTable = sinLUT->elems();

	LOOP(len,
		float mag = *magA;
		float phs = *phsA;
		
		phs = scl::wrap(phs, (float)M_2PI);
		phs *= (float)M_1_2PI;
		
		//ULONG index = scl::normalToUInt(phs);
		ULONG index = scl::normalToUInt2(phs);	// better cuz phase is on a linear scale
		
//		float sn = MemOpF::at(sinTable, fracBits, index);
//		float cs = MemOpF::at(sinTable, fracBits, index + quarter);

		float frac = tbl::fraction(bits, index);
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
	)
}


void rectToPolar(float * r, float * i, ULONG len){
	LOOP_P(len, scl::rectToPolar(*r, *i); r++; i++; )
}


void rectToPolarFast(float * realA, float * imagA, ULONG len){
	
	LOOP(len, 
		float real = *realA;
		float imag = *imagA;

		float realAbs = fabs(real);
		float imagAbs = fabs(imag);
		
		if(realAbs > imagAbs){
			float ratio = imag / real;	// (-1, 1)
			//float ratio = 0.5f;
			
			ULONG index = (ULONG)(LUTSize2F + LUTSize2F * ratio);
			//ULONG index = 10;
			//ULONG index = normalToIndex(ratio * 0.5f + 0.5f, 11);
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
			
			ULONG index = (ULONG)(LUTSize2F + LUTSize2F * ratio);
			//ULONG index = 10;
			//ULONG index = normalToIndex(ratio * 0.5f + 0.5f, 11);
			//int index = roundFloatToInt(LUTSize2F + LUTSize2F * ratio);

			if(imag > 0.f)  imag = M_PI_2 - atanLUT[index];		// pi/4, 3pi/4
			else			imag = M_3PI_2 - atanLUT[index];	// 5pi/4, 7pi/4
			real = imagAbs * magLUT[index];
		}
		
		*realA++ = real;
		*imagA++ = imag;
		
	)
}


void conversionInit(){
	if(0 == atanLUT){
		atanLUT = new float[LUTSize + 1];
		magLUT  = new float[LUTSize + 1];
		sinLUT  = new ArrayPow2<float>(2048);
		
		//printf("arr:conversionInit(): %d\n", sinLUT->size());

		tbl::sine(sinLUT->elems(), sinLUT->size());
			
		double LUTSizeRec = 1. / (double)LUTSize;
		double ramp2 = -1.;
		double ramp2Inc = 2. * LUTSizeRec;
		
		for(ULONG i=0; i<LUTSize; ++i){
			
			double angle = atan(ramp2);
			atanLUT[i] = angle;
			magLUT[i] = 1. / cos(angle);
			
			ramp2 += ramp2Inc;
		}
	}
}


void print(const float * src, ULONG len){
	for(ULONG i=0; i<len; i++) printf("[%4d]\t% f\n", (int)i, *src++);
}

void print(const float * src1, const float * src2, ULONG len){
	for(ULONG i=0; i<len; i++) printf("[%4d]\t% f % f\n", (int)i, *src1++, *src2++);
}

void printHex(const float * src, ULONG len){
	for(ULONG i=0; i<len; i++){
		float v = *src++;
		printf("[%4d] % 5.3f %8lx\n", (int)i, v, (unsigned long)*(ULONG *)&v);
	}
}

} // arr::
} // gam::

#include "MacroU.h"

