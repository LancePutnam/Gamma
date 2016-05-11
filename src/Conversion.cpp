/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include "Gamma/Conversion.h"
#include <cstring> // strlen

namespace gam{

uint32_t bits(const char * string){
	uint32_t v=0; int n = std::strlen(string);
	for(int i=0; i<n; ++i) if(string[i] == '1') v |= 1<<(n-1-i);
	return v;
}

uint32_t bitsToUInt(const char * bits){
	uint32_t i=0, r=0;
	for(; bits[i] && i<32; ++i) r |= ((bits[i]=='1'?1:0) << (31-i));
	return r>>(32-i);
}

uint32_t floatToUInt(float value){
	Twiddle<float> u(value);
	u.u += 0x800000;

	if(u.u & 0x40000000){	// mag outside [0, 1)		
		uint32_t shift = (u.u >> 23) & 0x7F;	
		return (1<<shift) | ((u.u & MaskFrac<float>()) >> (23 - shift));
	}
	else{
		return 0;
	}
}

int32_t floatToInt(float value){
	Twiddle<float> u(value);
	u.u += 0x800000;

	if(u.u & 0x40000000){	// mag outside [0, 1)
		int32_t shift = (u.u>>23) & 0x7F;
		int32_t sign = u.u & MaskSign<float>();
		int32_t result = (1<<shift) | ((u.u & MaskFrac<float>())>>(23-shift));
		
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
	u.u += 0x800000;

	if(u.u & 0x40000000){
		int32_t shift = (u.u>>23) & 0x7F;
		intPart = (1<<shift) | ((u.u & MaskFrac<float>())>>(23-shift));
		u.u = Expo1<float>() | ((u.u << shift) & MaskFrac<float>());
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
