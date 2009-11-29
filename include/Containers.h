#ifndef GAMMA_CONTAINERS_H_INC
#define GAMMA_CONTAINERS_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information

	File Description:
	Dynamically sized generic containers.
*/

#include <stdlib.h>
#include <vector>
#include <map>

#include "Conversion.h"
#include "mem.h"
#include "scl.h"

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
	ArrayBase(uint32_t size, const T& initial);
	ArrayBase(T * src, uint32_t size);
	ArrayBase(const ArrayBase<T,S>& src);

	virtual ~ArrayBase();

	/// Writable element access
	T& operator[](uint32_t i);
	
	/// Read-only element access
	const T& operator[](uint32_t i) const;
	
	/// Sets all elements to argument
	ArrayBase& operator=(const T& v){ for(uint32_t i=0; i<size(); ++i) (*this)[i] = v; return *this; }
	
	T * elems() const;				///< Returns pointer to first array element.
	uint32_t size() const;			///< Returns number of elements in array.

	void freeElements();			///< Frees memory

	/// Ensures ownership of elements.
	
	/// If the array is not the owner, new memory is allocated and the previously
	/// referenced array elements are copied.
	void own();
	
	void resize(uint32_t newSize);			///< Resizes number of elements in array
	void source(const ArrayBase<T,S>& src);	///< Sets source of array elements to another array
	void source(T * src, uint32_t size);	///< Sets source of array elements to another array

	virtual void onResize(){}

	static int references(T * const m){ return managing(m) ? refCount()[m] : 0; }

protected:
	T * mElems;
	S mSize;

	typedef std::map<T *, int> RefCount;

	static RefCount& refCount(){
		static RefCount * o = new RefCount;
		return *o;
	}
	
	// is memory being managed automatically?
	static bool managing(T* const m){ return refCount().count(m) != 0; }
};



/// Container for storing arrayed elements.
template <class T>
class Array : public ArrayBase<T, SizeArray>{
public:
	typedef ArrayBase<T, SizeArray> super;

	Array(): super(){}
	Array(uint32_t size): super(size){}
	Array(uint32_t size, const T& initial): super(size, initial){}
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
	ArrayPow2(uint32_t size, const T& initial): super(size, initial){}
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
	/// @param[in]	value		Initial value of all elements.
	Ring(uint32_t size, const T& value=T());
	
	/// Sets all elements to value
	Ring& operator=(const T& v){ super::operator=(v); return *this; }

	/// Returns reference to backmost element
	T& atBack(){ return (*this)[indexBack()]; }
	const T& atBack() const { return (*this)[indexBack()]; }
	
	/// Returns reference to frontmost element
	T& atFront(){ return (*this)[indexFront()]; }
	const T& atFront() const { return (*this)[indexFront()]; }
	
	/// Returns reference to element 'ago' indices from front
	T& atPrev(uint32_t ago){ return (*this)[indexPrev(ago)]; }
	const T& atPrev(uint32_t ago) const { return (*this)[indexPrev(ago)]; }

	/// Copies len elements starting from element pos() - delay into dst.
	void copy(T * dst, uint32_t len, uint32_t delay) const;
	
	/// Copy elements starting from last in into dst unwrapping from ring
	void copyUnwrap(T * dst, uint32_t len) const;
	
	uint32_t pos() const;			///< Return absolute index of frontmost (newest) element
	uint32_t indexBack() const;		///< Returns absolute index of backmost (oldest) element
	uint32_t indexFront() const;	///< Returns absolute index of frontmost (newest) element
	uint32_t indexPrev(uint32_t ago) const;	///< Returns absolute index of a previously written element

	void operator()(const T& v);	///< Write new element
	void pos(uint32_t index);		///< Set absolute buffer index of writer.
	void reset();					///< Reset write position to beginning
	void writeClip(const T& v);		///< Writes element unless at end of buffer.
	
protected:
	uint32_t mPos;
	void incPos();
};


/// Ring buffer that keeps track of its fill.
template <class T>
class RingFill : public Ring<T> {
public:
	typedef Ring<T> super;
	
	/// @param[in]	size		Number of elements in ring.
	/// @param[in]	value		Initial value of all elements.
	RingFill(uint32_t size, const T& value=T()): super(size, value), mFill(0){}

	void operator()(const T& v){
		this->super::operator()(v);
		if(mFill < super::size()) ++mFill;
	}
	
	/// Reset fill to zero.
	void reset(){ mFill=0; super::reset(); }
	
	/// Returns current fill of buffer.
	uint32_t fill() const { return mFill; }
	
protected:
	uint32_t mFill;
	virtual void onResize(){ reset(); }
};



/// Double buffered ring-buffer
template <class T>
struct DoubleRing : public Ring<T>{

	/// @param[in]	size		Number of elements in ring.
	/// @param[in]	value		Initial value of all elements.
	DoubleRing(uint32_t size, const T& value=T()): Ring<T>(size, value), read(size){}
	
	/// Copy elements into read buffer unwrapping from ring
	
	/// Returns a pointer to the read buffer
	///
	T * copyUnwrap(){ Ring<T>::copyUnwrap(read.elems(), read.size()); return read.elems(); }
	
	/// Copy elements into read buffer "as is" from ring
	
	/// Returns a pointer to the read buffer
	///
	T * copy(){ mem::deepCopy(read.elems(), Ring<T>::elems(), read.size()); return read.elems(); }
	
	Array<T> read;
};




/// Fixed N-sample delay
template <class T>
struct DelayN: public Ring<T>{
	using Ring<T>::incPos; using Ring<T>::pos;

	/// @param[in]	size		Number of elements to delay.
	/// @param[in]	value		Initial value of all elements.
	DelayN(uint32_t size, const T& value=T()): Ring<T>(size, value){}
	
	DelayN& operator=(const T& v){ Ring<T>::operator=(v); return *this; }

	/// Write new element, return oldest.
	T operator()(T input){
		incPos();				// inc write pos
		T dly = (*this)[pos()];	// read oldest element
		(*this)[pos()] = input;	// write new element
		return dly;				//	... write pos left at last written element
	}
};






// Implementation_______________________________________________________________

#define TEM template<class T>

// ArrayBase

#undef TEM2
#define TEM2 template <class T, class S>
#define ARRAYBASE_INIT mElems(0), mSize(0)

TEM2 ArrayBase<T,S>::ArrayBase()
:	ARRAYBASE_INIT{}

TEM2 ArrayBase<T,S>::ArrayBase(uint32_t size)
:	ARRAYBASE_INIT
{	resize(size); }

TEM2 ArrayBase<T,S>::ArrayBase(uint32_t size, const T& initial)
:	ARRAYBASE_INIT
{	resize(size); for(uint32_t i=0;i<this->size();++i) (*this)[i] = initial; }

TEM2 ArrayBase<T,S>::ArrayBase(T * src, uint32_t size)
:	ARRAYBASE_INIT
{	source(src, size); }

TEM2 ArrayBase<T,S>::ArrayBase(const ArrayBase<T,S>& src)
:	ARRAYBASE_INIT
{	source(src); }

#undef ARRAYBASE_INIT

TEM2 ArrayBase<T,S>::~ArrayBase(){ freeElements(); }


TEM2 inline T& ArrayBase<T,S>::operator[](uint32_t i){ return elems()[i]; }
TEM2 inline const T&  ArrayBase<T,S>::operator[](uint32_t i) const { return elems()[i]; }

TEM2 inline T * ArrayBase<T,S>::elems() const { return mElems; }

TEM2 void ArrayBase<T,S>::freeElements(){ //printf("ArrayBase::freeElements(): mElems=%p\n", mElems);
	
//	if(owner()){
//		//printf("(%p) ArrayBase::freeElements(): mElems=%p\n", this, mElems);
//		
//		delete [] mElems; 
//		mElems=0; mSize(0);
//	}
	
	// we are managing this memory
	if(mElems && managing(mElems)){
		int& c = refCount()[mElems];
		--c;
		if(0 == c){
			refCount().erase(mElems);
			delete[] mElems;
		}
		mElems=0; mSize(0);
	}
}

TEM2 void ArrayBase<T,S>::own(){
//	if(!owner()){
//		mOwner = true;
//		uint32_t refSize = size();
//		T * refElems = elems();
//		mElems = 0;
//		mSize(0);
//		resize(refSize);						// allocate new memory
//		mem::copy(elems(), refElems, size());	// copy referenced elements
//	}
	
	T * oldElems = elems();
	if(references(oldElems) != 1){
		uint32_t oldSize = size();
		freeElements();
		resize(oldSize);
		for(uint32_t i=0; i<size(); ++i) (*this)[i] = oldElems[i];
	}
}

TEM2 void ArrayBase<T,S>::resize(uint32_t newSize){
	newSize = mSize.convert(newSize);

	if(newSize != size()){
		freeElements();
		mElems = new T[newSize];
		refCount()[mElems] = 1;
		mSize(newSize);
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
	if(managing(src)){
		++refCount()[src];
	}
	mElems = src;
	mSize(size);
}

#undef TEM2



// ArrayPow2

TEM inline uint32_t ArrayPow2<T>::oneIndex() const { return 1<<fracBits(); }
TEM inline uint32_t ArrayPow2<T>::log2Size() const { return mSize.mBitsI; }
TEM inline uint32_t ArrayPow2<T>::fracBits() const { return mSize.mBitsF; }
TEM inline uint32_t ArrayPow2<T>::index(uint32_t phase) const { return phase >> fracBits(); }

TEM inline T ArrayPow2<T>::atPhase(uint32_t phase) const { return (*this)[index(phase)]; }
TEM inline void ArrayPow2<T>::putPhase(uint32_t phase, T v){ (*this)[index(phase)] = v; }

TEM inline float ArrayPow2<T>::fraction(uint32_t phase) const{		
	return gam::fraction(log2Size(), phase);
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

TEM Ring<T>::Ring(uint32_t size, const T& v) : Array<T>(size,v), mPos(size-1){}

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

TEM inline uint32_t Ring<T>::indexBack() const {
	uint32_t i = pos() + 1;
	return (i != size()) ? i : 0;
}

TEM inline uint32_t Ring<T>::indexFront() const { return pos(); }

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

#undef TEM

#endif
