/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include "Conversion.h"

namespace gam{

uint32_t bitsToUInt(const char * bits){
	uint32_t i=0, r=0;
	for(; bits[i] && i<32; ++i) r |= ((bits[i]=='1'?1:0) << (31-i));
	return r>>(32-i);
}


uint32_t bytesToUInt32(const uint8_t * bytes4){
	uint32_t word = 0;

	if(0 == endian){
		word  = bytes4[3] << 24;
		word |= (bytes4[2] & 0xff) << 16;
		word |= (bytes4[1] & 0xff) << 8;
		word |= bytes4[0] & 0xff;
	}
	else{
		word  = bytes4[0] << 24;
		word |= (bytes4[1] & 0xff) << 16;
		word |= (bytes4[2] & 0xff) << 8;
		word |= bytes4[3] & 0xff;
	}
	
	return word;
}


uint16_t bytesToUInt16(const uint8_t * bytes2){
	uint16_t word = 0;

	if(0 == endian){
		word  = bytes2[0] & 0xff;
		word |= (bytes2[1] & 0xff) << 8;
	}
	else{
		word  = bytes2[1] & 0xff;
		word |= (bytes2[0] & 0xff) << 8;
	}
	
	return word;
}


uint32_t floatToUInt(float value){
	Twiddle<float> u(value);
	u.u += 0x800000;

	if(u.u & 0x40000000){	// mag outside [0, 1)		
		uint32_t shift = (u.u >> 23) & 0x7F;	
		return (1<<shift) | ((u.u & MASK_F32_FRAC) >> (23 - shift));
	}
	else{
		return 0;
	}
}


int32_t floatToInt(float value){
	Twiddle<float> u(value);
	u.u = (u.u + 0x800000);

	if(u.u & 0x40000000){	// mag outside [0, 1)
		int32_t shift = ((u.u)>>23) & 0x7F;
		int32_t sign = u.u & 0x80000000;
		int32_t result = (1<<shift) | ((u.u & MASK_F32_FRAC)>>(23-shift));
		
		if(sign){	// negative number
			result = ~result + 1;	// 2's complement
		}
		return result;
	}
	else{
		return 0;
	}
}


float split(float value, int32_t& intPart){
	Twiddle<float> u(value);

	u.u = (u.u + 0x800000);

	if(u.u & 0x40000000){
		int32_t shift = ((u.u)>>23) & 0x7F;
		intPart = (1<<shift) | ((u.u & MASK_F32_FRAC)>>(23-shift));
		u.u = 0x3F800000 | ((u.u << shift) & MASK_F32_FRAC);
		return u.f - 1.f;
	}
	else{
		intPart = 0;
		return value;
	}
}

//ULONG SclOp::floatToUInt(float value){
//	union { float f; ULONG u; } word;
//	word.f = value;
//
//	word.u = (word.u + 0x800000);
//
//	if(word.u & 0x40000000){	// mag outside [0, 1)
//		//int shift = ((word.u)>>23) & 0x7F;
//		ULONG shift = ((word.u)>>23) & 0x7F;
//		return (1<<shift) | ((word.u & MASK_F32_FRAC)>>(23-shift));
//	}
//	else{
//		return 0;
//	}
//}

} // gam::
