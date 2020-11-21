#ifndef GAMMA_MEM_H_INC
#define GAMMA_MEM_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include <cstring> // memcpy, memmove, etc.
#include <cstdlib> // free, realloc
#include "Gamma/Access.h" // indexLast
#include "Gamma/std_swap.h"

#define LOOP(n,s) for(unsigned i=0; i<n; i+=s)

namespace gam{

/// Memory functions

///	These operations only involve moving data positions in memory and thus are 
/// completely object type independent, i.e. generic.  Some operations require 
/// the objects to understand the = and == operators.
/// Operations post-fixed with 'L' or 'R' are referring to the direction of
/// the action either 'to the left' or 'to the right', respectively.
namespace mem{

/// Copies every 'chunkSize'th element to 'dst'; akin to downsampling.
template <class T>
void compact(T * dst, const T * src, unsigned len, unsigned chunkSize);

//////
////// Ring buffer copying functions.
//////
////// Elements are copied starting at ring + ringSize and elements past the end of 
////// ring are wrapped back to the beginning of ring.  No bounds checking is performed
////// on the ring taps.  length must be <= ringSize

/// Copies elements into a ring buffer. 

/// Returns the next tap index.  This will not guaranteed to be in the range
/// [0, ringSize) and therefore should be wrapped into the proper range.
template <class T>
unsigned copyToRing(T * ring, unsigned ringSize, unsigned ringTap, const T * src, unsigned len);

/// Copies a number of element from a ring buffer.
template <class T>
void copyFromRing(const T * ring, unsigned ringSize, unsigned ringTap, T * dst, unsigned len);

/// Copies all elements from a ring buffer.

/// \param[in] ring			pointer to source ring buffer
/// \param[in] ringSize		number of elements in source ring
/// \param[in] ringTap		ring location to start copy from, usually location of oldest element
/// \param[in] dst			destination array (size must at least match ring size)
template <class T>
void copyAllFromRing(const T * ring, unsigned ringSize, unsigned ringTap, T * dst);

/// Copy source 'src' array bytes to destination 'dst'.
template <class T>
void deepCopy(T * dst, const T * src, unsigned len){ std::memcpy(dst, src, len*sizeof(T)); }

/// Returns whether array bytes are equal (==).
template <class T>
bool deepEqual(const T * src1, const T * src2, unsigned len){
	return 0 == std::memcmp(src1, src2, len*sizeof(T));
}

/// Move source 'src' array bytes to destination 'dst'. Arrays can overlap.
template <class T>
void deepMove(T * dst, const T * src, unsigned len){ std::memmove(dst, src, len*sizeof(T)); }

/// Set elements' bytes to zero.
template <class T>
void deepZero(T * arr, unsigned len){ std::memset(arr, 0, len*sizeof(T)); }

/// Deinterleave an array of elements with any number of channels.

/// numFrames = length / numChannels
///
template <class T>
void deinterleave(T * dst, const T * src, unsigned numFrames, unsigned numChannels);

/// Deinterleave an array of elements with 2 channels.

/// numFrames = length / numChannels \n
/// Example: abababab -> aaaabbbb
template <class T>
void deinterleave2(T * dst, const T * src, unsigned numFrames);

/// Expands elements from 'src' to 'dst'.

/// Elements are copied contiguously from 'src' to strided locations in 'dst'.
/// 'dst' must have room for 'lenSrc' times 'amount' elements. \n
/// Ex.: 1234 -> 1.2.3.4.     (amount = 2)
/// Ex.: 1234 -> 1..2..3..4.. (amount = 3)
template <class T>
void expand(T * dst, const T * src, unsigned lenSrc, unsigned amount);

/// Like standard free, but checks if pointer is valid (!=0) and sets it to zero after freeing it.

/// This uses C-style memory management. No destructors will be called on
/// class-type objects.
template <class T>
void free(T *& ptr);

/// Finds index 'ind' of first element in 'src' equal to 'elem'.

/// Returns whether a match was found.
///
template <class T>
bool indexOf(const T& elem, unsigned& ind, const T * src, unsigned len, unsigned str=1);

/// Interleave an array of elements with any number of channels.

/// numFrames = length / numChannels
///
template <class T>
void interleave(T * dst, const T * src, unsigned numFrames, unsigned numChannels);

/// Interleave an array of elements with 2 channels.

/// numFrames = length / numChannels \n	
/// Example: aaaabbbb -> abababab
template <class T>
void interleave2(T * dst, const T * src, unsigned numFrames);

/// Keeps every Nth element; the rest are zeroed.

/// \param[in]	arr		Array to operate on.
/// \param[in]	len		Number of elements in array.
/// \param[in]	stride	Spacing between kept elements.
/// \param[in]	offset	Offset of spacing from start of array.
template <class T>
void keep(T * arr, unsigned len, unsigned stride, unsigned offset=0);

/// Perform kth permutation on array.

///	The length of the array must be <= 20.
///
template <class T>
void permute(T * arr, unsigned len, unsigned long long k);

/// Pivot elements around 'index' element.
template <class T>
void pivot(unsigned index, T * arr, unsigned len);

/// Reflect right half of array to left half.

/// The first len/2 elements will be overwritten. \n
/// Example: ....1234 -> 43211234
template <class T>
void reflectLeft(T * arr, unsigned len, unsigned str=1);

/// Reflect left half of array to right half.

/// The last len/2 elements will be overwritten. \n
/// Example: 1234.... -> 12344321
template <class T>
void reflectRight(T * arr, unsigned len, unsigned str=1);

/// Repeats first 'chunkSize' elements until end of array.

/// Example: 1234.... -> 12341234
///
template <class T>
void repeat(unsigned chunkSize, T * arr, unsigned len);

/// Replace instances of one element with another
template <class T>
void replace(const T& val, const T& with, T * arr, unsigned len, unsigned str=1){
	LOOP(len,str){ if(arr[i] == val) arr[i] = with; }
}

/// Replace multiple instances of one element with another
template <class T>
void replace(const T * val, const T * with, unsigned num, T * arr, unsigned len, unsigned str=1){
	LOOP(num,1){ replace(val[i], with[i], arr, len,str); }
}

/// Resizes array.  Returns true if resized, false otherwise.

/// This uses C-style memory management. No constructors or destructors will
/// be called on class-type objects.
template <class T>
bool resize(T *& arr, unsigned sizeNow, unsigned sizeNew);

/// Reverse elements' order in array.

///	Example: 1234 -> 4321
///
template <class T>
void reverse(T * arr, unsigned len, unsigned str=1);

/// Reverse every two elements.

/// Example: 1234 -> 2143
///
template <class T>
void reverse2(T * arr, unsigned len, unsigned str=1);

/// Rotate elements half the length of the array.

/// Works only for even length arrays.
///
template <class T>
void rotateHalf(T * arr, unsigned len, unsigned str=1);

/// Rotate elements left by 1.

/// Example: 1234 -> 2341
///
template <class T>
void rotateLeft1(T * arr, unsigned len, unsigned str=1);

/// Rotate elements left by 'order' elements.
template <class T>
void rotateLeft(unsigned order, T * arr, unsigned len);

/// Rotate elements right by 1.

/// Example: 1234 -> 4123
///
template <class T>
void rotateRight1(T * arr, unsigned len, unsigned str=1);

/// Rotate elements right by 'order' elements.
template <class T>
void rotateRight(unsigned order, T * arr, unsigned len);

/// Copies elements from 'src' to fractionally strided locations in 'dst'. 

/// The destination indices are rounded up.
/// 'dst' must have room for 'lenSrc' times 'amount' elements. \n
/// Ex.: 12345678 -> 12.34.56.78 (amount = 1.5) \n
/// Ex.: 1234 -> 1..2..3..4.. (amount = 3)
template <class T>
void scale(T * dst, const T * src, unsigned lenSrc, float stride);

/// Copies elements from 'src' to fractionally strided locations in 'dst'. 

/// The destination indices are rounded up with overrun elements being 
/// cropped, i.e. not copied.\n
/// Ex.: 12345678 -> 12.34.56 (amount = 1.5) \n
/// Ex.: 1234 -> 1..2 (amount = 3)
template <class T>
void scaleCrop(T * dst, const T * src, unsigned len, float stride);

/// Returns (n+1)th argument.
template <class T>
inline const T& select(int n, const T& v1, const T& v2, const T& v3){
	return (n<=0)?v1 : (n==1)?v2 : v3;
}

/// Returns (n+1)th argument.
template <class T>
inline const T& select(int n, const T& v1, const T& v2, const T& v3, const T& v4){
	return (n<=0)?v1 : (n==1)?v2 : (n==2)?v3 : v4;
}

/// Stretches array by duplicating every element 'amount' times.

/// 'len' must be an integer multiple of 'amount'. \n
/// Example: 1234.... -> 11223344 
template <class T>
void stretch(T * arr, unsigned len, unsigned amount);

/// Stretches array by duplicating every element 'amount' times.

/// 'dst' must have room for 'lenSrc' times 'amount' elements. \n
/// Example: 1234 -> 11223344	
template <class T>
void stretch(T * dst, const T * src, unsigned lenSrc, unsigned amount);

/// Transpose a 2 x (len/2) matrix.

/// Example: 12121212 -> 11112222
///
template <class T>
void transpose2(T * arr, unsigned len);


/// Returns true if all elements are different
template <class T>
bool unique(const T * src, unsigned len, unsigned str=1);

template <class T>
bool unique(const T * src, unsigned len, unsigned str){
	for(unsigned j=0; j<len-str; j+=str){
		const T& v = src[j];
		for(unsigned i=j+str; i<len; i+=str){
			if(v == src[i]) return false;
		}
	}
	return true;
}


/// Sets elements to zero
template <class T>
inline void zero(T * arr, unsigned len, unsigned str=1){
	LOOP(len,str){ arr[i]=T(0); }
}



// Implementation_______________________________________________________________

template <class T>
inline void compact(T * dst, const T * src, unsigned len, unsigned chunkSize){
	if(chunkSize < 2){			deepCopy(dst, src, len); return;	}
	else if(chunkSize > len){	*dst = *src; return; }
	
	for(unsigned i=0; i<len; i+=chunkSize){
		*dst++ = src[i];
	}
}

template <class T>
void copy(T * dst, const T * src, const unsigned * indices, unsigned numIndices){
	LOOP(numIndices,1){
		unsigned index = *indices++;
		dst[index] = src[index];
	}
}

template <class T>
unsigned copyToRing(T * ring, unsigned ringSize, unsigned ringTap, const T * src, unsigned len){

	unsigned endTap = ringTap + len;

	if(endTap <= ringSize){		// haven't gone past end
		deepCopy(ring + ringTap, src, len);
	}
	else{						// have gone past end, do wrapped copy
		unsigned under	= ringSize - ringTap;
		unsigned over	= endTap - ringSize;
		deepCopy(ring + ringTap, src, under);
		deepCopy(ring, src + under, over);
	}

	return endTap;
}

template <class T>
void copyFromRing(const T * ring, unsigned ringSize, unsigned ringTap, T * dst, unsigned len){

	unsigned endTap = ringTap + len;

	if(endTap <= ringSize){
		deepCopy(dst, ring + ringTap, len);
	}
	else{
		unsigned under	= ringSize - ringTap;
		unsigned over	= endTap - ringSize;
		deepCopy(dst, ring + ringTap, under);
		deepCopy(dst + under, ring, over);
	}

}

template <class T>
inline void copyAllFromRing(const T * ring, unsigned ringSize, unsigned ringTap, T * dst){
	unsigned under = ringSize - ringTap;
	deepCopy(dst, ring + ringTap, under);
	deepCopy(dst + under, ring, ringTap);
}

template <class T>
inline void expand(T * dst, const T * src, unsigned lenSrc, unsigned amount){
	LOOP(lenSrc,1){ *dst = *src++; dst += amount; }
}

template <class T>
inline void free(T *& ptr){
	if(ptr){ std::free(ptr); ptr=0; }
}

template <class T>
inline bool indexOf(const T& v, unsigned& ind, const T * src, unsigned len, unsigned str){
	LOOP(len,str){ if(src[i] == v){ ind=i; return true; } }
	return false;
}

// s=2, o=0		1.1.1.1.
// s=2, o=1		.1.1.1.1
// s=3, o=0		1..1..1.
// s=3, o=1		.1..1..1
// s=3, o=2		..1..1..

template <class T>
inline void keep(T * arr, unsigned len, unsigned stride, unsigned offset){
	unsigned c = offset % stride;
	LOOP(len,1){
		if(0 == c){	arr++; c = stride; }
		else		*arr++ = (T)0;
		--c;
	}
}

template <class T>
void permute(T * arr, unsigned len, unsigned long long k) {
	unsigned long long factorial = 1;
	for(unsigned i=1; i<len; i++){
		factorial *= i;
		unsigned j = i - ((k / factorial) % (i + 1));
		std::swap(arr[i], arr[j]);
	}
}

template <class T>
inline void pivot(unsigned index, T * arr, unsigned len){
	
	unsigned lt = index - 1;
	unsigned rt = index + 1;
	
	LOOP((len-1)>>1,1){
		if(lt >= len) lt = len - 1;
		if(rt >= len) rt = 0;
		std::swap(arr[lt], arr[rt]);
		lt--;
		rt++;
	}
}

template <class T>
inline void reflectLeft(T * arr, unsigned len, unsigned str){
	unsigned end = indexLast(len,str);
	LOOP(len>>1, str){
		arr[i] = arr[end-i];
	}
}

template <class T>
inline void reflectRight(T * arr, unsigned len, unsigned str){
	unsigned end = indexLast(len,str);
	LOOP(len>>1, str){
		arr[end-i] = arr[i];
	}	
}	

template <class T>
inline void repeat(unsigned chunkSize, T * arr, unsigned len){
	unsigned rd = 0;
	unsigned wr = chunkSize;
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
template <class T>
bool resize(T *& arr, unsigned sizeNow, unsigned sizeNew){
	if((sizeNow != sizeNew) && (0 != sizeNew)){
		T * ptr = (T *)std::realloc(arr, sizeNew * sizeof(T));
		if(0 != ptr){ arr = ptr; return true; }	// successful resize
	}
	return false;
}

template <class T>
inline void reverse(T * arr, unsigned len, unsigned str){
	unsigned end = indexLast(len,str);
	LOOP(len>>1, str){ std::swap(arr[i], arr[end-i]); }
}

template <class T>
inline void reverse2(T * arr, unsigned len, unsigned str){
	LOOP(len, 2){ std::swap(arr[i], arr[i+1]); }
}

template <class T>
inline void rotateHalf(T * arr, unsigned len, unsigned str){
	T * half = arr + (((len>>1)/str)*str);
	LOOP(len>>1, str){ std::swap(arr[i], half[i]); }
}

template <class T>
void rotateLeft(unsigned order, T * arr, unsigned len){
	
	order %= len;

	unsigned numSwaps = len>>1;
	
	if( (len & 1) == 0 && (order & 1) == 1 ) numSwaps--;

	unsigned rt = (order + 1)>>1;
	unsigned lt = rt - 1 - (order & 1);
	
	LOOP(numSwaps, 1){
		if(lt >= len) lt = len - 1;
		if(rt >= len) rt = 0;
		std::swap(arr[lt--], arr[rt++]);
	}
	
	reverse(arr, len);
}

template <class T>
inline void rotateLeft1(T * arr, unsigned len, unsigned str){
	LOOP(len-str, str){	std::swap(arr[i], arr[i+str]); }
}

template <class T>
void rotateRight(unsigned order, T * arr, unsigned len){
	rotateLeft(len - order, arr, len);
}

template <class T>
inline void rotateRight1(T * arr, unsigned len, unsigned str){
	unsigned end = indexLast(len,str);
	for(unsigned i=end; i>=str; i-=str){
		std::swap(arr[i], arr[i-str]);
	}
}

template <class T>
inline void scale(T * dst, const T * src, unsigned lenSrc, float amount){
	*dst = *src++;
	float dstIndex = 0.5f;							// add offset to round
	lenSrc--;

	LOOP(lenSrc, 1){
		dstIndex += amount;
		dst[(unsigned)dstIndex] = *src++;
		//printf("%d\n", (unsigned)dstIndex);
	}
}

template <class T>
inline void scaleCrop(T * dst, const T * src, unsigned len, float amount){
	unsigned temp = (unsigned)((float)len / amount);	// theoretical # to copy from src
	if(temp < len) len = temp;					// actual # to copy from src

	scale(dst, src, len, amount);
}

template <class T>
inline void stretch(T * arr, unsigned len, unsigned amount){

	T * end = arr + len;
	
	len /= amount;
	arr += len;

	LOOP(len, 1){
		end -= amount;
		T v = *--arr;
		for(unsigned k=0; k<amount; ++k) end[k]=v;
	}
}

template <class T>
inline void stretch(T * dst, const T * src, unsigned lenSrc, unsigned amount){
	LOOP(lenSrc, 1){
		T v = *src++;
		for(unsigned k=0; k<amount; ++k) dst[k]=v;
		dst += amount;
	}
}

template <class T>
inline void transpose2(T * arr, unsigned len){
	for(unsigned i=2; i <= len-2; i+=2){
		reverse2(arr + (i>>1), len - i);
	}
}

template <class T>
inline void deinterleave(T * dst, const T * src, unsigned numFrames, unsigned numChannels){
	unsigned numSamples = numFrames * numChannels;
	for(unsigned c=0; c < numChannels; c++){
		for(unsigned i=c; i < numSamples; i+=numChannels){
			*dst++ = src[i];
		}
	}
}

template <class T>
inline void deinterleave2(T * dst, const T * src, unsigned numFrames){
	T * dst2 = dst + numFrames;
	LOOP(numFrames, 1){
		*dst2++ = *src++;
		*dst++  = *src++;
	}
}

template <class T>
inline void interleave(T * dst, const T * src, unsigned numFrames, unsigned numChannels){
	unsigned numSamples = numFrames * numChannels;
	for(unsigned c=0; c < numChannels; c++){
		for(unsigned i=c; i < numSamples; i+=numChannels){
			dst[i] = *src++;
		}
	}
}

template <class T>
inline void interleave2(T * dst, const T * src, unsigned numFrames){
	const T * src2 = src + numFrames;
	LOOP(numFrames, 1){
		*dst++ = *src2++;
		*dst++ = *src++;
	}
}

/// Duplicates strided elements in array
/// a_a_a_a_ -> aaaaaaaa	numFrames = 4, stride = 2
/// b___b___ -> bbbbbbbb	numFrames = 2, stride = 4
inline void arrayStridedDup(float * dst, const float * src, unsigned numFrames, unsigned stride){
	for(unsigned i=1; i<stride; i++){
		for(unsigned j=0; j<numFrames * stride; j+=stride){
			dst[i + j] = src[j];
		}
	}
}

} // mem::
} // gam::

#undef LOOP
#undef SWAP

#endif

