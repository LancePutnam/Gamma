#ifndef GAMMA_CONTAINERS_H_INC
#define GAMMA_CONTAINERS_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include <stdlib.h>
#include <vector>

#include "Types.h"
#include "gen.h"
#include "mem.h"
#include "scl.h"

#include "MacroD.h"

namespace gam{


/// Size type for ArrayPow2
struct SizeArrayPow2{
	SizeArrayPow2(uint32_t size){ (*this)(size); }
	uint32_t operator()() const { return (1<<mBitsI) & 0xfffffffe/*avoids 1*/; }
	void operator()(uint32_t v){ mBitsI = scl::log2(convert(v)); mBitsF = 32U - mBitsI; /*printf("%d %d\n", mBitsI, mBitsF);*/ }
	uint32_t convert(uint32_t v){ v=scl::ceilPow2(v); return v!=1 ? v : 2; }	// should return 0,2,4,8,16,...
	uint32_t mBitsI;	// integer portion # bits
	uint32_t mBitsF;	// fraction portion # bits
};


/// Size type for Array
struct SizeArray{
	SizeArray(uint32_t size): mSize(size){}
	uint32_t operator()() const { return mSize; }
	void operator()(uint32_t v){ mSize = v; }
	uint32_t convert(uint32_t v){ return v; }
	uint32_t mSize;
};


/// Abstract base class for array types

/// When the array is resized, if the elements are class-types, then their
/// default constructors are called and if the elements are non-class-types,
/// then they are left uninitialized.
template <class T, class S>
class ArrayBase{
public:
	ArrayBase();
	ArrayBase(uint32_t size);
	ArrayBase(T * src, uint32_t size);
	ArrayBase(const ArrayBase<T,S>& src);

	virtual ~ArrayBase();

	T& operator[](uint32_t i);
	const T& operator[](uint32_t i) const;
	
	ArrayBase& operator=(const T& v){ for(uint32_t i=0; i<size(); ++i) (*this)[i] = v; return *this; }
	
	T * elems() const;				///< Returns pointer to first array element.
	bool owner() const;				///< Returns whether element memory is owned.
	uint32_t size() const;			///< Returns number of elements in array.

	void freeElements();			///< Frees memory

	/// Ensures ownership of elements.
	
	/// If the array is not the owner, new memory is allocated and the previously
	/// referenced array elements are copied.
	void own();
	
	void resize(uint32_t newSize);			///< Resizes number of elements in array
	void source(const ArrayBase<T,S>& src);	///< Sets source of array elements to another array
	void source(T * src, uint32_t size);	///< Sets source of array elements to another array
	//void zero();							// Zeroes array elements

	virtual void onResize(){}

protected:
	T * mElems;
	S mSize;
	bool mOwner;
};



/// Container for storing arrayed elements.
template <class T>
class Array : public ArrayBase<T, SizeArray>{
public:
	typedef ArrayBase<T, SizeArray> super;

	Array(): super(){}
	Array(uint32_t size): super(size){}
	Array(T * src, uint32_t size): super(src, size){}
	Array(const Array<T>& src): super(src){}
	
	Array& operator=(const T& v){ super::operator=(v); return *this; }

	virtual ~Array(){}
};



/// Container for storing a power-of-two number of arrayed elements.
template <class T>
class ArrayPow2 : public ArrayBase<T, SizeArrayPow2>{
public:
	typedef ArrayBase<T, SizeArrayPow2> super;

	ArrayPow2(): super(){}
	ArrayPow2(uint32_t size): super(size){}
	ArrayPow2(T * src, uint32_t size): super(src, size){}
	ArrayPow2(const ArrayPow2<T>& src): super(src){}

	virtual ~ArrayPow2(){}

	uint32_t fracBits() const;			///< Returns number of bits for fraction (32 - bits()).
	float fraction(uint32_t phase) const;
	uint32_t index(uint32_t phase) const;
	uint32_t log2Size() const;			///< Returns log base-2 of the number of array elements.
	uint32_t oneIndex() const;			///< Returns 32-bit phase acc increment for one element index.

	T atPhase(uint32_t phase) const;
	void putPhase(uint32_t phase, T v);
	
private:
	using super::mSize;
};



/// Container for sharing arrays between threads.
template <class T>
class Buffer : public Array<T> {
public:

	/// @param[in]	size		Number of elements in buffer.
	Buffer(uint32_t size);
	virtual ~Buffer(){}
	
	/// Writes 'src' to buffer and locks buffer.
	
	/// The buffer will only be written to if it is unlocked.  'src' must have
	/// at least this.size() elements.
	void writeLock(const T * src);

	/// Writes 'src' to a portion of the buffer and locks buffer.
	
	/// @param[in]	src			Source array to write into buffer.
	/// @param[in]	numRead		Number of elements to read from 'src'.
	/// @param[in]	writeOffset	Starting write position in buffer.
	/// The buffer will only be written to if it is unlocked.  'src' must have
	/// at least 'numRead' elements.  'writePos' + 'numRead' should not exceed
	/// this.size().
	void writeLock(const T * src, uint32_t numRead, uint32_t writeOffset=0);
	
	void lock();		///< Locks buffer to prevent writing.
	void unlock();		///< Unlocks buffer to permit writing.
	
	bool isLocked();	///< Returns whether the buffer is locked.
	bool isUnlocked();	///< Returns whether the buffer is unlocked.
	
protected:
	bool mLocked;
};




/// Container for reading and writing arrays using front and back buffers.
template <class T>
class DoubleBuffer{
public:
	
	/// @param[in]	singleBufSize	Number of elements in single buffer (allocated size will be twice this).
	DoubleBuffer(uint32_t singleBufSize);
	~DoubleBuffer();

	/// Writes single buffer to back and performs swap().
	void write(const T * singleBuffer);
										
	void swap();	///< Swaps the front and back buffers.

	T * front();	///< Returns a pointer to the front (newer) buffer.
	T * back();		///< Returns a pointer to the back (older) buffer.
	void copyAll(T * dst);	///< Copy contents to dst, back buffer first then front.

protected:
	T * mBuffer;
	uint32_t mSizeSingle;
	uint32_t mTap;
};




/// Ring buffer
template <class T>
class Ring : public Array<T> {
public:

	typedef Array<T> super; using super::elems; using super::size;

	/// @param[in]	size		Number of elements in ring.
	Ring(uint32_t size);
	
	Ring& operator=(const T& v){ super::operator=(v); return *this; }

	const T& atPrev(uint32_t ago) const { return (*this)[indexPrev(ago)]; }
	
	/// Copies len elements starting from element pos() - delay into dst.
	void copy(T * dst, uint32_t len, uint32_t delay) const;
	
	/// Copy elements starting from last in into dst unwrapping from ring
	void copyUnwrap(T * dst, uint32_t len) const;
	
	uint32_t pos() const;			///< Return absolute index of last written element
	uint32_t indexPrev(uint32_t ago) const;	///< Returns absolute index of a previously written element

	void operator()(const T& v);	///< Write new element
	void pos(uint32_t index);		///< Set absolute buffer index of writer.
	void reset();					///< Reset write position to beginning
	void writeClip(const T& v);		///< Writes element unless at end of buffer.
	
protected:
	uint32_t mPos;
	void incPos();
};




/// Double buffered ring-buffer
template <class T>
struct DoubleRing : public Ring<T>{
	DoubleRing(uint32_t size): Ring<T>(size), read(size){}
	
	/// Copy elements into read buffer unwrapping from ring
	
	/// Returns a pointer to the read buffer
	///
	T * copyUnwrap(){ Ring<T>::copyUnwrap(read.elems(), read.size()); return read.elems(); }
	
	/// Copy elements into read buffer "as is" from ring
	
	/// Returns a pointer to the read buffer
	///
	T * copy(){ mem::copy(read.elems(), Ring<T>::elems(), read.size()); return read.elems(); }
	
	Array<T> read;
};




/// Fixed N-sample delay
template <class T>
struct DelayN: public Ring<T>{
	using Ring<T>::incPos; using Ring<T>::pos;

	/// @param[in]	size		Number of elements to delay.
	DelayN(uint32_t size): Ring<T>(size){}
	
	DelayN& operator=(const T& v){ Ring<T>::operator=(v); return *this; }

	/// Write new element, return oldest.
	T operator()(T input){
		incPos();				// inc write pos
		T dly = (*this)[pos()];	// read oldest element
		(*this)[pos()] = input;	// write new element
		return dly;				//	... write pos left at last written element
	}
};




/// Fixed-sized array with a sequence generator
template <uint32_t N, class T=gam::real, class G=gen::RAdd1<uint32_t> >
class Seq: public Multi<N,T>{
public:

	Seq(const T& value){ mem::set(this->elems, gen::val(value), N); }
	Seq(const T * values){ mem::copy(this->elems, values, N); }

	/// Generate next array element
	T operator()(){ return (*this)[((uint32_t)mTap())%N]; }

	/// Get reference to index generator
	G& tap(){ return mTap; }
	
private:
	G mTap;
};






// Implementation_______________________________________________________________

// ArrayBase

#undef TEM2
#define TEM2 template <class T, class S>
#define ARRAYBASE_INIT mElems(0), mSize(0)

TEM2 ArrayBase<T,S>::ArrayBase()
:	ARRAYBASE_INIT, mOwner(true){}

TEM2 ArrayBase<T,S>::ArrayBase(uint32_t size)
:	ARRAYBASE_INIT, mOwner(true)
{	resize(size); }

TEM2 ArrayBase<T,S>::ArrayBase(T * src, uint32_t size)
:	ARRAYBASE_INIT, mOwner(false)
{	source(src, size); }

TEM2 ArrayBase<T,S>::ArrayBase(const ArrayBase<T,S>& src)
:	mSize(src.size()), mElems(src.elems()), mOwner(false){}

#undef ARRAYBASE_INIT

TEM2 ArrayBase<T,S>::~ArrayBase(){ freeElements(); }


TEM2 inline T& ArrayBase<T,S>::operator[](uint32_t i){ return elems()[i]; }
TEM2 inline const T&  ArrayBase<T,S>::operator[](uint32_t i) const { return elems()[i]; }

TEM2 inline T * ArrayBase<T,S>::elems() const { return mElems; }

TEM2 void ArrayBase<T,S>::freeElements(){ //printf("ArrayBase::freeElements(): mElems=%p\n", mElems);
	//if(owner()){ mem::free(mElems); mSize(0); }
	//if(owner()){ delete[] mElems; mElems=0; mSize(0); }	// TODO: delete[] is causing crash
	
	if(owner()){
		//printf("(%p) ArrayBase::freeElements(): mElems=%p\n", this, mElems);
		delete [] mElems; 
		mElems=0; mSize(0);
	}
}

TEM2 void ArrayBase<T,S>::own(){
	if(!owner()){
		mOwner = true;
		uint32_t refSize = size();
		T * refElems = elems();
		mElems = 0;
		mSize(0);
		resize(refSize);						// allocate new memory
		mem::copy(elems(), refElems, size());	// copy referenced elements
	}
}

TEM2 inline bool ArrayBase<T,S>::owner() const { return mOwner; }

TEM2 void ArrayBase<T,S>::resize(uint32_t newSize){
	newSize = mSize.convert(newSize);
	//if(owner() && mem::resize(mElems, size(), newSize)){
	if(owner() && (newSize != size())){
		freeElements();
		mElems = new T[newSize];
		mSize(newSize);
		//zero();
		onResize();
	}
	//printf("ArrayBase::resize(): mElems=%p, size=%d\n", mElems, size());
}

TEM2 inline uint32_t ArrayBase<T,S>::size() const { return mSize(); }

TEM2 void ArrayBase<T,S>::source(const ArrayBase<T,S>& src){
	source(src.elems(), src.size());
}

TEM2 void ArrayBase<T,S>::source(T * src, uint32_t size){
	freeElements();
	mOwner = false;
	mElems = src;
	mSize(size);
}

//TEM2 void ArrayBase<T,S>::zero(){ if(owner()) mem::zero(elems(), size()); }

#undef TEM2



// ArrayPow2

TEM inline uint32_t ArrayPow2<T>::oneIndex() const { return 1<<fracBits(); }
TEM inline uint32_t ArrayPow2<T>::log2Size() const { return mSize.mBitsI; }
TEM inline uint32_t ArrayPow2<T>::fracBits() const { return mSize.mBitsF; }
TEM inline uint32_t ArrayPow2<T>::index(uint32_t phase) const { return phase >> fracBits(); }

TEM inline T ArrayPow2<T>::atPhase(uint32_t phase) const { return (*this)[index(phase)]; }
TEM inline void ArrayPow2<T>::putPhase(uint32_t phase, T v){ (*this)[index(phase)] = v; }

TEM inline float ArrayPow2<T>::fraction(uint32_t phase) const{	
	phase = phase << log2Size() >> 9 | 0x3f800000;
	return scl::punUF32(phase) - 1.f;
}



// Buffer

TEM Buffer<T>::Buffer(uint32_t size)
	: Array<T>(size)
{
	unlock();
}

TEM inline void Buffer<T>::writeLock(const T * src){
	if(isLocked()) return;
	lock();
	memcpy(this->elems(), src, this->size() * sizeof(T));
}

TEM inline void Buffer<T>::writeLock(const T * src, uint32_t numRead, uint32_t writeOffset){
	if(isLocked()) return;
	lock();
	memcpy(this->elems() + writeOffset, src, numRead * sizeof(T));
}

TEM inline void Buffer<T>::unlock(){ mLocked = false; }
TEM inline void Buffer<T>::lock(){ mLocked = true; }
TEM inline bool Buffer<T>::isLocked(){ return mLocked; }
TEM inline bool Buffer<T>::isUnlocked(){ return !mLocked; }



//---- Ring

TEM Ring<T>::Ring(uint32_t size) : Array<T>(size), mPos(size-1){}

TEM inline void Ring<T>::operator()(const T& v){
	incPos();				// inc write pos; do first to avoid out-of-bounds access
	(*this)[pos()] = v;		// write new element
}

TEM void Ring<T>::copy(T * dst, uint32_t len, uint32_t delay) const{
	// pos() points to most recently written slot
	//uint32_t tap = (pos() - delay) % size();
	uint32_t tap = (uint32_t)scl::wrap((int32_t)pos() - (int32_t)delay, (int32_t)size());

	// this ensures that we don't copy across the ring tap boundary
	// we add one to maxLen because of a fence post anomaly
	uint32_t maxLen = (tap < pos() ? (pos() - tap) : (pos() + (size() - tap))) + 1;
	len = scl::min(len, maxLen);
	
	mem::copyFromRing(elems(), size(), tap, dst, len);
}

TEM void Ring<T>::copyUnwrap(T * dst, uint32_t len) const { copy(dst, len, size() - 1); }

TEM inline uint32_t Ring<T>::indexPrev(uint32_t v) const {
	return scl::wrapOnce<int>(pos() - v, size());
}

TEM inline uint32_t Ring<T>::pos() const { return mPos; }

TEM inline void Ring<T>::incPos(){ if(++mPos >= size()) mPos = 0; }
TEM void Ring<T>::pos(uint32_t index){ mPos = index; }

TEM void Ring<T>::reset(){ pos(size()-1); }

TEM	inline void Ring<T>::writeClip(const T& v){
	if(mPos < size()){
		(*this)[mPos] = v;
		mPos++;
	}
}


// DoubleBuffer

TEM DoubleBuffer<T> :: DoubleBuffer(uint32_t singleBufSize){
	mTap = 0;
	mSizeSingle = singleBufSize;
	mBuffer = (T *)malloc((singleBufSize<<1) * sizeof(T));
}

TEM DoubleBuffer<T> :: ~DoubleBuffer(){
	free(mBuffer);
}

TEM inline void DoubleBuffer<T>::write(const T * singleBuffer){
	memcpy(mBuffer + mTap, singleBuffer, mSizeSingle * sizeof(T));
	swap();
}

TEM inline void DoubleBuffer<T>::swap(){
	mTap = mSizeSingle - mTap;
}

TEM inline T * DoubleBuffer<T>::front(){
	return mBuffer + (mSizeSingle - mTap);
}

TEM inline T * DoubleBuffer<T>::back(){
	return mBuffer + mTap;
}

TEM inline void DoubleBuffer<T>::copyAll(T * dst){
	if(0 == mTap){
		memcpy(dst, mBuffer, (mSizeSingle<<1) * sizeof(T));
	}
	else{
		memcpy(dst, mBuffer + mSizeSingle, mSizeSingle * sizeof(T));
		memcpy(dst + mSizeSingle, mBuffer, mSizeSingle * sizeof(T));
	}
}

} // end namespace gam


#include "MacroU.h"

#endif

