#ifndef GAMMA_MEM_H_INC
#define GAMMA_MEM_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "Gamma/Access.h"

#define LOOP(n,s) for(uint32_t i=0; i<n; i+=s)
#define TEM template <class T>
#define TEM2 template <class T1, class T2>

namespace gam{

/// Memory functions

///	These operations only involve moving data positions in memory and thus are 
/// completely object type independent, i.e. generic.  Some operations require 
/// the objects to understand the = and == operators.
/// Operations post-fixed with 'L' or 'R' are referring to the direction of
/// the action either 'to the left' or 'to the right', respectively.
namespace mem{

/// Copies every 'chunkSize'th element to 'dst'; akin to downsampling.
TEM void compact(T * dst, const T * src, uint32_t len, uint32_t chunkSize);

//////
////// Ring buffer copying functions.
//////
////// Elements are copied starting at ring + ringSize and elements past the end of 
////// ring are wrapped back to the beginning of ring.  No bounds checking is performed
////// on the ring taps.  length must be <= ringSize

/// Copies elements into a ring buffer. 

/// Returns the next tap index.  This will not guaranteed to be in the range
/// [0, ringSize) and therefore should be wrapped into the proper range.
TEM uint32_t copyToRing(T * ring, uint32_t ringSize, uint32_t ringTap, const T * src, uint32_t len);

/// Copies a number of element from a ring buffer.
TEM void copyFromRing(const T * ring, uint32_t ringSize, uint32_t ringTap, T * dst, uint32_t len);

/// Copies all elements from a ring buffer.

/// @param[in] ring			pointer to source ring buffer
/// @param[in] ringSize		number of elements in source ring
/// @param[in] ringTap		ring location to start copy from, usually location of oldest element
/// @param[in] dst			destination array (size must at least match ring size)
TEM void copyAllFromRing(const T * ring, uint32_t ringSize, uint32_t ringTap, T * dst);

/// Copy source 'src' array bytes to destination 'dst'.
TEM void deepCopy(T * dst, const T * src, uint32_t len){ memcpy(dst, src, len*sizeof(T)); }

/// Returns whether array bytes are equal (==).
TEM bool deepEqual(const T * src1, const T * src2, uint32_t len){
	return 0 == memcmp(src1, src2, len*sizeof(T));
}

/// Move source 'src' array bytes to destination 'dst'. Arrays can overlap.
TEM void deepMove(T * dst, const T * src, uint32_t len){ memmove(dst, src, len*sizeof(T)); }

/// Set elements' bytes to zero.
TEM void deepZero(T * arr, uint32_t len){ memset(arr, 0, len*sizeof(T)); }

/// Deinterleave an array of elements with any number of channels.

/// numFrames = length / numChannels
///
TEM void deinterleave(T * dst, const T * src, uint32_t numFrames, uint32_t numChannels);

/// Deinterleave an array of elements with 2 channels.

/// numFrames = length / numChannels \n
/// Example: abababab -> aaaabbbb
TEM void deinterleave2(T * dst, const T * src, uint32_t numFrames);

/// Expands elements from 'src' to 'dst'.

/// Elements are copied contiguously from 'src' to strided locations in 'dst'.
/// 'dst' must have room for 'lenSrc' times 'amount' elements. \n
/// Ex.: 1234 -> 1.2.3.4.     (amount = 2)
/// Ex.: 1234 -> 1..2..3..4.. (amount = 3)
TEM void expand(T * dst, const T * src, uint32_t lenSrc, uint32_t amount);

/// Like standard free, but checks if pointer is valid (!=0) and sets it to zero aftyer freeing it.

/// This uses C-style memory management. No destructors will be called on
/// class-type objects.
TEM void free(T *& ptr);

/// Finds index 'ind' of first element in 'src' equal to 'elem'.

/// Returns whether a match was found.
///
TEM bool indexOf(const T& elem, uint32_t& ind, const T * src, uint32_t len, uint32_t str=1);

/// Interleave an array of elements with any number of channels.

/// numFrames = length / numChannels
///
TEM void interleave(T * dst, const T * src, uint32_t numFrames, uint32_t numChannels);

/// Interleave an array of elements with 2 channels.

/// numFrames = length / numChannels \n	
/// Example: aaaabbbb -> abababab
TEM void interleave2(T * dst, const T * src, uint32_t numFrames);

/// Keeps every Nth element; the rest are zeroed.

/// @param[in]	arr		Array to operate on.
/// @param[in]	len		Number of elements in array.
/// @param[in]	stride	Spacing between kept elements.
/// @param[in]	offset	Offset of spacing from start of array.
TEM void keep(T * arr, uint32_t len, uint32_t stride, uint32_t offset=0);

/// Perform kth permutation on array.

///	The length of the array must be <= 20.
///
TEM void permute(T * arr, uint32_t len, uint64_t k);

/// Pivot elements around 'index' element.
TEM void pivot(uint32_t index, T * arr, uint32_t len);

/// Reflect right half of array to left half.

/// The first len/2 elements will be overwritten. \n
/// Example: ....1234 -> 43211234
TEM void reflectLeft(T * arr, uint32_t len, uint32_t str=1);

/// Reflect left half of array to right half.

/// The last len/2 elements will be overwritten. \n
/// Example: 1234.... -> 12344321
TEM void reflectRight(T * arr, uint32_t len, uint32_t str=1);

/// Repeats first 'chunkSize' elements until end of array.

/// Example: 1234.... -> 12341234
///
TEM void repeat(uint32_t chunkSize, T * arr, uint32_t len);

/// Replace instances of one element with another
TEM void replace(const T& val, const T& with, T * arr, uint32_t len, uint32_t str=1){
	LOOP(len,str){ if(arr[i] == val) arr[i] = with; }
}

/// Replace multiple instances of one element with another
TEM void replace(const T * val, const T * with, uint32_t num, T * arr, uint32_t len, uint32_t str=1){
	LOOP(num,1){ replace(val[i], with[i], arr, len,str); }
}

/// Resizes array.  Returns true if resized, false otherwise.

/// This uses C-style memory management. No constructors or destructors will
/// be called on class-type objects.
TEM bool resize(T *& arr, uint32_t sizeNow, uint32_t sizeNew);

/// Reverse elements' order in array.

///	Example: 1234 -> 4321
///
TEM void reverse(T * arr, uint32_t len, uint32_t str=1);

/// Reverse every two elements.

/// Example: 1234 -> 2143
///
TEM void reverse2(T * arr, uint32_t len, uint32_t str=1);

/// Rotate elements half the length of the array.

/// Works only for even length arrays.
///
TEM void rotateHalf(T * arr, uint32_t len, uint32_t str=1);

/// Rotate elements left by 1.

/// Example: 1234 -> 2341
///
TEM void rotateLeft1(T * arr, uint32_t len, uint32_t str=1);

/// Rotate elements left by 'order' elements.
TEM void rotateLeft(uint32_t order, T * arr, uint32_t len);

/// Rotate elements right by 1.

/// Example: 1234 -> 4123
///
TEM void rotateRight1(T * arr, uint32_t len, uint32_t str=1);

/// Copies elements from 'src' to fractionally strided locations in 'dst'. 

/// The destination indices are rounded up.
/// 'dst' must have room for 'lenSrc' times 'amount' elements. \n
/// Ex.: 12345678 -> 12.34.56.78    (amount = 1.5) \n
/// Ex.: 1234 -> 1..2..3..4.. (amount = 3)
TEM void scale(T * dst, const T * src, uint32_t lenSrc, float stride);

/// Copies elements from 'src' to fractionally strided locations in 'dst'. 

/// The destination indices are rounded up with overrun elements being 
/// cropped, i.e. not copied.\n
/// Ex.: 12345678 -> 12.34.56 (amount = 1.5) \n
/// Ex.: 1234 -> 1..2 (amount = 3)
TEM void scaleCrop(T * dst, const T * src, uint32_t len, float stride);

/// Returns (n+1)th argument.
TEM inline const T& select(int n, const T& v1, const T& v2, const T& v3){
	return (n<=0)?v1 : (n==1)?v2 : v3;
}

/// Returns (n+1)th argument.
TEM inline const T& select(int n, const T& v1, const T& v2, const T& v3, const T& v4){
	return (n<=0)?v1 : (n==1)?v2 : (n==2)?v3 : v4;
}

/// Stretches array by duplicating every element 'amount' times.

/// 'len' must be an integer multiple of 'amount'. \n
/// Example: 1234.... -> 11223344 
TEM void stretch(T * arr, uint32_t len, uint32_t amount);

/// Stretches array by duplicating every element 'amount' times.

/// 'dst' must have room for 'lenSrc' times 'amount' elements. \n
/// Example: 1234 -> 11223344	
TEM void stretch(T * dst, const T * src, uint32_t lenSrc, uint32_t amount);

/// Swaps two elements in memory.
TEM void swap(T & elem1, T & elem2);

/// Transpose a 2 x (len/2) matrix.

/// Example: 12121212 -> 11112222
///
TEM void transpose2(T * arr, uint32_t len);


/// Returns true if all elements are different
TEM bool unique(const T * src, uint32_t len, uint32_t str=1);
TEM bool unique(const T * src, uint32_t len, uint32_t str){
	for(uint32_t j=0; j<len-str; j+=str){
		const T& v = src[j];
		for(uint32_t i=j+str; i<len; i+=str){
			if(v == src[i]) return false;
		}
	}
	return true;
}


/// Sets elements to zero
TEM inline void zero(T * arr, uint32_t len, uint32_t str=1){
	LOOP(len,str){ arr[i]=T(0); }
}

//TEM inline void zero(T * arr, unsigned int len){ zero(arr, (uint)len); }
//
//template <class To0> void zero(To0& o0){ zero(&o0[0], (uint)o0.size()); }
//
///// o0[i] = 0
//template <class To0, class Tind>
//void zero(To0& o0, const Tind& ind){ LOOP_IND(o0[i] = 0;) }

/// Get value from a power-of-two array.

///	'src':		waveform \n
///	'fbits':	number of bits in fractional part of phase \n
///	'phase':	fixed-point phase of lookup (full period is [0, 2^32))
TEM T at(const T * src, uint32_t fbits, uint32_t phase);

/// Set value in a power-of-two array.

///	'dst':		waveform \n
///	'fbits':	number of bits in fractional part of phase \n
///	'phase':	fixed-point phase of lookup (full period is [0, 2^32))	
TEM void put(T * dst, uint32_t fbits, uint32_t phase, T value);



// Implementation_______________________________________________________________

TEM inline void compact(T * dst, const T * src, uint32_t len, uint32_t chunkSize){
	if(chunkSize < 2){			deepCopy(dst, src, len); return;	}
	else if(chunkSize > len){	*dst = *src; return; }
	
	for(uint32_t i=0; i<len; i+=chunkSize){
		*dst++ = src[i];
	}
}

TEM void copy(T * dst, const T * src, const uint32_t * indices, uint32_t numIndices){
	LOOP(numIndices,1){
		uint32_t index = *indices++;
		dst[index] = src[index];
	}
}

TEM uint32_t copyToRing(T * ring, uint32_t ringSize, uint32_t ringTap, const T * src, uint32_t len){

	uint32_t endTap = ringTap + len;

	if(endTap <= ringSize){		// haven't gone past end
		deepCopy(ring + ringTap, src, len);
	}
	else{						// have gone past end, do wrapped copy
		uint32_t under	= ringSize - ringTap;
		uint32_t over	= endTap - ringSize;
		deepCopy(ring + ringTap, src, under);
		deepCopy(ring, src + under, over);
	}

	return endTap;
}

TEM void copyFromRing(const T * ring, uint32_t ringSize, uint32_t ringTap, T * dst, uint32_t len){

	uint32_t endTap = ringTap + len;

	if(endTap <= ringSize){
		deepCopy(dst, ring + ringTap, len);
	}
	else{
		uint32_t under	= ringSize - ringTap;
		uint32_t over	= endTap - ringSize;
		deepCopy(dst, ring + ringTap, under);
		deepCopy(dst + under, ring, over);
	}

}

TEM inline void copyAllFromRing(const T * ring, uint32_t ringSize, uint32_t ringTap, T * dst){
	uint32_t under = ringSize - ringTap;
	deepCopy(dst, ring + ringTap, under);
	deepCopy(dst + under, ring, ringTap);
}

TEM inline void expand(T * dst, const T * src, uint32_t lenSrc, uint32_t amount){	
	LOOP(lenSrc,1){ *dst = *src++; dst += amount; }
}

TEM inline void free(T *& ptr){
	if(ptr){ ::free(ptr); ptr=0; }
}

TEM inline bool indexOf(const T& v, uint32_t& ind, const T * src, uint32_t len, uint32_t str){	
	LOOP(len,str){ if(src[i] == v){ ind=i; return true; } }
	return false;
}

// s=2, o=0		1.1.1.1.
// s=2, o=1		.1.1.1.1
// s=3, o=0		1..1..1.
// s=3, o=1		.1..1..1
// s=3, o=2		..1..1..

TEM inline void keep(T * arr, uint32_t len, uint32_t stride, uint32_t offset){	
	uint32_t c = offset % stride;
	LOOP(len,1){
		if(0 == c){	arr++; c = stride; }
		else		*arr++ = (T)0;
		--c;
	}
}

TEM void permute(T * arr, uint32_t len, uint64_t k) {
     uint64_t factorial = 1;
     for(uint32_t i=1; i<len; i++){
        factorial *= i;
		uint32_t j = i - ((k / factorial) % (i + 1));
		swap(arr[i], arr[j]);
     }
}

TEM inline void pivot(uint32_t index, T * arr, uint32_t len){
	
	uint32_t lt = index - 1;
	uint32_t rt = index + 1;
	
	LOOP((len-1)>>1,1){
		if(lt >= len) lt = len - 1;
		if(rt >= len) rt = 0;
		swap(arr[lt], arr[rt]);
		lt--;
		rt++;
	}
}

TEM inline void reflectLeft(T * arr, uint32_t len, uint32_t str){	
	uint32_t end = indexLast(len,str);
	LOOP(len>>1, str){
		arr[i] = arr[end-i];
	}
}

TEM inline void reflectRight(T * arr, uint32_t len, uint32_t str){
	uint32_t end = indexLast(len,str);
	LOOP(len>>1, str){
		arr[end-i] = arr[i];
	}	
}	

TEM inline void repeat(uint32_t chunkSize, T * arr, uint32_t len){
	uint32_t rd = 0;
	uint32_t wr = chunkSize;
	LOOP(len - chunkSize, 1){
		if(rd == chunkSize) rd = 0;
		arr[wr] = arr[rd];
		wr++;
		rd++;
	}
}

/*	void * realloc(void *ptr, size_t size);
	 
     The realloc() function tries to change the size of the allocation pointed
     to by ptr to size, and return ptr.  If there is not enough room to
     enlarge the memory allocation pointed to by ptr, realloc() creates a new
     allocation, copies as much of the old data pointed to by ptr as will fit
     to the new allocation, frees the old allocation, and returns a pointer to
     the allocated memory.  realloc() returns a NULL pointer if there is an
     error, and the allocation pointed to by ptr is still valid.
*/
TEM bool resize(T *& arr, uint32_t sizeNow, uint32_t sizeNew){
	if((sizeNow != sizeNew) && (0 != sizeNew)){
		T * ptr = (T *)realloc(arr, sizeNew * sizeof(T));
		if(0 != ptr){ arr = ptr; return true; }	// successful resize
	}
	return false;
}

TEM inline void reverse(T * arr, uint32_t len, uint32_t str){
	uint32_t end = indexLast(len,str);
	LOOP(len>>1, str){ swap(arr[i], arr[end-i]); }
}

TEM inline void reverse2(T * arr, uint32_t len, uint32_t str){	
	LOOP(len, 2){ swap(arr[i], arr[i+1]); }
}

TEM inline void rotateHalf(T * arr, uint32_t len, uint32_t str){	
	T * half = arr + (((len>>1)/str)*str);
	LOOP(len>>1, str){ swap(arr[i], half[i]); }
}

TEM void rotateLeft(uint32_t order, T * arr, uint32_t len){
	
	order %= len;

	uint32_t numSwaps = len>>1;
	
	if( (len & 1) == 0 && (order & 1) == 1 ) numSwaps--;

	uint32_t rt = (order + 1)>>1;
	uint32_t lt = rt - 1 - (order & 1);
	
	LOOP(numSwaps, 1){
		if(lt >= len) lt = len - 1;
		if(rt >= len) rt = 0;
		swap(arr[lt--], arr[rt++]);
	}
	
	reverse(arr, len);
}

TEM inline void rotateLeft1(T * arr, uint32_t len, uint32_t str){
	LOOP(len-str, str){	swap(arr[i], arr[i+str]); }
}

TEM inline void rotateRight1(T * arr, uint32_t len, uint32_t str){
	uint32_t end = indexLast(len,str);
	for(uint32_t i=end; i>=str; i-=str){
		swap(arr[i], arr[i-str]);
	}
}

TEM inline void scale(T * dst, const T * src, uint32_t lenSrc, float amount){
	*dst = *src++;
	float dstIndex = 0.5f;							// add offset to round
	lenSrc--;

	LOOP(lenSrc, 1){
		dstIndex += amount;
		dst[(uint32_t)dstIndex] = *src++;
		//printf("%d\n", (uint32_t)dstIndex);
	}
}

TEM inline void scaleCrop(T * dst, const T * src, uint32_t len, float amount){
	uint32_t temp = (uint32_t)((float)len / amount);	// theoretical # to copy from src
	if(temp < len) len = temp;					// actual # to copy from src

	scale(dst, src, len, amount);
}

TEM inline void stretch(T * arr, uint32_t len, uint32_t amount){

	T * end = arr + len;
	
	len /= amount;
	arr += len;

	LOOP(len, 1){
		end -= amount;
		T v = *--arr;
		for(uint32_t k=0; k<amount; ++k) end[k]=v;
	}
}

TEM inline void stretch(T * dst, const T * src, uint32_t lenSrc, uint32_t amount){
	LOOP(lenSrc, 1){
		T v = *src++;
		for(uint32_t k=0; k<amount; ++k) dst[k]=v;
		dst += amount;
	}
}

TEM inline void swap(T& a, T& b){
	T t=a; a=b; b=t;
}

TEM inline void transpose2(T * arr, uint32_t len){
	for(uint32_t i=2; i <= len-2; i+=2){
		reverse2(arr + (i>>1), len - i);
	}
}

TEM inline void deinterleave(T * dst, const T * src, uint32_t numFrames, uint32_t numChannels){
	uint32_t numSamples = numFrames * numChannels;
	for(uint32_t c=0; c < numChannels; c++){
		for(uint32_t i=c; i < numSamples; i+=numChannels){
			*dst++ = src[i];
		}
	}
}

TEM inline void deinterleave2(T * dst, const T * src, uint32_t numFrames){	
	T * dst2 = dst + numFrames;
	LOOP(numFrames, 1){
		*dst2++ = *src++;
		*dst++  = *src++;
	}
}

TEM inline void interleave(T * dst, const T * src, uint32_t numFrames, uint32_t numChannels){
	uint32_t numSamples = numFrames * numChannels;
	for(uint32_t c=0; c < numChannels; c++){
		for(uint32_t i=c; i < numSamples; i+=numChannels){
			dst[i] = *src++;
		}
	}
}

TEM inline void interleave2(T * dst, const T * src, uint32_t numFrames){	
	const T * src2 = src + numFrames;
	LOOP(numFrames, 1){
		*dst++ = *src2++;
		*dst++ = *src++;
	}
}

TEM inline T at(const T * src, uint32_t fbits, uint32_t phase){
	return src[phase >> fbits];
}

TEM inline void put(T * dst, uint32_t fbits, uint32_t phase, T value){
	dst[phase >> fbits] = value;
}

/// Duplicates strided elements in array
/// a_a_a_a_ -> aaaaaaaa	numFrames = 4, stride = 2
/// b___b___ -> bbbbbbbb	numFrames = 2, stride = 4
inline void arrayStridedDup(float * dst, const float * src, uint32_t numFrames, uint32_t stride){
	for(uint32_t i=1; i<stride; i++){
		for(uint32_t j=0; j<numFrames * stride; j+=stride){
			dst[i + j] = src[j];
		}
	}
}

} // mem::
} // gam::

#undef LOOP
#undef TEM
#undef SWAP

#endif

