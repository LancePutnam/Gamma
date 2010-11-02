#ifndef GAMMA_CONTAINERS_H_INC
#define GAMMA_CONTAINERS_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information

	File Description:
	Dynamically sizable generic containers.
*/

#include <stdlib.h>
#include <vector>
#include <map>

#include "Gamma/Allocator.h"
#include "Gamma/Conversion.h"
#include "Gamma/mem.h"
#include "Gamma/scl.h"

namespace gam{


/// Size type for ArrayPow2
struct SizeArrayPow2{
	SizeArrayPow2(uint32_t size){ (*this)(size); }
	uint32_t operator()() const { return (1<<mBitsI) & 0xfffffffe/*avoids 1*/; }
	void operator()(uint32_t v){ mBitsI = scl::log2(convert(v)); mBitsF = 32U - mBitsI; /*printf("%d %d\n", mBitsI, mBitsF);*/ }
	static uint32_t convert(uint32_t v){ v=scl::ceilPow2(v); return v!=1 ? v : 2; }	// should return 0,2,4,8,16,...
	uint32_t mBitsI;	// integer portion # bits
	uint32_t mBitsF;	// fraction portion # bits
};


/// Size type for Array
struct SizeArray{
	SizeArray(uint32_t size): mSize(size){}
	uint32_t operator()() const { return mSize; }
	void operator()(uint32_t v){ mSize = v; }
	static uint32_t convert(uint32_t v){ return v; }
	uint32_t mSize;
};


/// Abstract base class for array types

/// When the array is resized, if the elements are class-types, then their
/// default constructors are called and if the elements are non-class-types,
/// then they are left uninitialized.
template <class T, class S, class A=gam::Allocator<T> >
class ArrayBase : private A{
public:
	ArrayBase();
	explicit ArrayBase(uint32_t size);
	ArrayBase(uint32_t size, const T& init);
	ArrayBase(T * src, uint32_t size);
	explicit ArrayBase(const ArrayBase<T,S,A>& src);

	virtual ~ArrayBase();

	/// Set element
	T& operator[](uint32_t i);
	
	/// Get element
	const T& operator[](uint32_t i) const;
	
	/// Sets all elements to value
	ArrayBase& assign(const T& v);

	/// Sets linear slice of elements to value
	
	/// @param[in] v		value to be copied as new content
	/// @param[in] end		end index (exclusive)
	/// @param[in] stride	index stride amount
	/// @param[in] start	start index (inclusive)
	ArrayBase& assign(const T& v, uint32_t end, uint32_t stride=1, uint32_t start=0);
	
	T * elems() const;				///< Returns pointer to first array element.
	uint32_t size() const;			///< Returns number of elements in array.

	/// Destroys all elements and frees memory
	void clear();

	/// Ensures ownership of elements.
	
	/// If the array is not the owner, new memory is allocated and the previously
	/// referenced array elements are copied.
	void own();
	
	/// Resizes number of elements in array
	
	/// If the new size is less than the old size, then elements are truncated.
	/// If the new size is greater than the old size, then the argument value
	/// is copied into the additional elements.
	void resize(uint32_t newSize, const T& c=T());
	void source(const ArrayBase<T,S,A>& src);	///< Sets source of array elements to another array
	void source(T * src, uint32_t size);	///< Sets source of array elements to another array

	virtual void onResize(){}

	/// Returns number of pointers to memory address being managed
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

private: ArrayBase& operator=(const ArrayBase& v);
};



/// Resizable array
template <class T, class A=gam::Allocator<T> >
class Array : public ArrayBase<T, SizeArray, A>{
public:
	typedef ArrayBase<T, SizeArray, A> Base;

	Array(): Base(){}
	explicit Array(uint32_t size): Base(size){}
	Array(uint32_t size, const T& init): Base(size, init){}
	Array(T * src, uint32_t size): Base(src, size){}
	Array(const Array& src): Base(src){}

	virtual ~Array(){}

private: Array& operator=(const Array& v);
};



/// Resizable array with a power-of-2 number of elements
template <class T, class A=gam::Allocator<T> >
class ArrayPow2 : public ArrayBase<T, SizeArrayPow2, A>{
public:
	typedef ArrayBase<T, SizeArrayPow2, A> Base;

	ArrayPow2(): Base(){}
	explicit ArrayPow2(uint32_t size): Base(size){}
	ArrayPow2(uint32_t size, const T& initial): Base(size, initial){}
	ArrayPow2(T * src, uint32_t size): Base(src, size){}
	explicit ArrayPow2(const ArrayPow2& src): Base(src){}

	virtual ~ArrayPow2(){}

	uint32_t fracBits() const;				///< Returns number of bits in fraction (32 - bits())
	float fraction(uint32_t phase) const;	///< Get floating-point fractional part of fixed-point phase
	uint32_t index(uint32_t phase) const;	///< Get integer part of fixed-point phase
	uint32_t log2Size() const;				///< Returns log base-2 of the number of array elements
	uint32_t oneIndex() const;				///< Returns 32-bit phase increment for one element index

	const T& atPhase(uint32_t phase) const;	///< Get element at truncated fixed-point phase
	void putPhase(uint32_t phase, T v);		///< Set element at truncated fixed point phase
	
private:
	using Base::mSize;
private: ArrayPow2& operator=(const ArrayPow2& v);
};



/// Container for sharing arrays between threads.
template <class T, class A=gam::Allocator<T> >
class Buffer : public Array<T,A> {
public:

	/// @param[in]	size		Number of elements in buffer.
	explicit Buffer(uint32_t size);
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
template <class T, class A=gam::Allocator<T> >
class DoubleBuffer : public Array<T,A>{
public:
	
	/// @param[in]	singleBufSize	Number of elements in single buffer (allocated size will be twice this).
	explicit DoubleBuffer(uint32_t singleBufSize);

	/// Set front buffer element
	T& operator[](uint32_t i){ return front()[i]; }
	
	/// Get front buffer element
	const T& operator[](uint32_t i) const { return front()[i]; }

	T * back();				///< Returns a pointer to the back (older) buffer
	T * front();			///< Returns a pointer to the front (newer) buffer

	void swap();			///< Swaps the front and back buffers

	void copyAll(T * dst);	///< Copy contents to dst, back buffer first then front

protected:
	uint32_t mSizeSingle;
	T * mFront, * mBack;
	
	typedef Array<T,A> Base;
	
	virtual void onResize();
	bool backIsFirst();

	/// Writes single buffer to back and performs swap().
	//void write(const T * singleBuffer);
};



/// Ring buffer
template <class T, class A=gam::Allocator<T> >
class Ring : public Array<T,A> {
public:

	typedef Array<T,A> Base; using Base::elems; using Base::size;

	/// @param[in]	size		Number of elements in ring.
	/// @param[in]	value		Initial value of all elements.
	explicit Ring(uint32_t size, const T& value=T());

	/// Returns reference to backmost (oldest) element
	T& atBack(){ return (*this)[indexBack()]; }
	const T& atBack() const { return (*this)[indexBack()]; }
	
	/// Returns reference to frontmost (newest) element
	T& atFront(){ return (*this)[indexFront()]; }
	const T& atFront() const { return (*this)[indexFront()]; }
	
	/// Returns reference to element 'ago' indices behind front
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
	void pos(uint32_t index);		///< Set absolute buffer index of writer
	void reset();					///< Reset write position to beginning
	void writeClip(const T& v);		///< Writes element unless at end of buffer
	
protected:
	uint32_t mPos;
	void incPos();
};


/// Ring buffer that keeps track of its fill amount.
template <class T, class A=gam::Allocator<T> >
class RingFill : public Ring<T,A> {
public:
	typedef Ring<T,A> Base;
	
	/// @param[in]	size		Number of elements in ring.
	/// @param[in]	value		Initial value of all elements.
	explicit RingFill(uint32_t size, const T& value=T()): Base(size, value), mFill(0){}

	/// Write new element
	void operator()(const T& v){
		this->Base::operator()(v);
		if(mFill < Base::size()) ++mFill;
	}
	
	/// Reset internal tap and fill amount to zero
	void reset(){ mFill=0; Base::reset(); }
	
	/// Returns current fill amount of buffer
	
	/// The fill amount is the number of elements written to the buffer since
	/// the last call to reset(). The fill amount has a maximum value equal
	/// to the size of the buffer.
	uint32_t fill() const { return mFill; }
	
protected:
	uint32_t mFill;
	virtual void onResize(){ reset(); }
};



/// Double buffered ring-buffer
template <class T, class A=gam::Allocator<T> >
class DoubleRing : public Ring<T,A>{
public:
	/// @param[in]	size		Number of elements in ring.
	/// @param[in]	value		Initial value of all elements.
	explicit DoubleRing(uint32_t size, const T& value=T()): Ring<T>(size, value), mRead(size){}
	
	/// Copy elements into read buffer unwrapping from ring
	
	/// Returns a pointer to the read buffer
	///
	T * copyUnwrap(){ Ring<T,A>::copyUnwrap(mRead.elems(), mRead.size()); return mRead.elems(); }
	
	/// Copy elements into read buffer "as is" from ring
	
	/// Returns a pointer to the read buffer
	///
	T * copy(){
		mem::deepCopy(mRead.elems(), Ring<T,A>::elems(), mRead.size());
		//for(uint32_t i=0; i<read.size(); ++i) construct(read.elems()+i, (*this)[i]);
		return mRead.elems();
	}
	
	// Returns reference to the read buffer
	const Array<T,A>& read() const { return mRead; }

protected:	
	Array<T,A> mRead;
};




/// N-sample delay
template <class T, class A=gam::Allocator<T> >
struct DelayN: public Ring<T,A>{
	using Ring<T,A>::incPos; using Ring<T,A>::pos;

	/// @param[in]	size		Number of elements to delay.
	/// @param[in]	value		Initial value of all elements.
	explicit DelayN(uint32_t size, const T& value=T()): Ring<T,A>(size, value){}

	/// Write new element and return oldest
	T operator()(const T& input){
		incPos();				// inc write pos
		T dly = (*this)[pos()];	// read oldest element
		(*this)[pos()] = input;	// write new element
		return dly;				//	... write pos left at last written element
	}
};






// Implementation_______________________________________________________________


// ArrayBase

#define TEM3 template <class T, class S, class A>
#define ARRAYBASE_INIT mElems(0), mSize(0)

TEM3 ArrayBase<T,S,A>::ArrayBase()
:	ARRAYBASE_INIT{}

TEM3 ArrayBase<T,S,A>::ArrayBase(uint32_t sz)
:	ARRAYBASE_INIT
{	resize(sz); }

TEM3 ArrayBase<T,S,A>::ArrayBase(uint32_t sz, const T& initial)
:	ARRAYBASE_INIT
{	resize(sz); for(uint32_t i=0;i<this->size();++i) (*this)[i] = initial; }

TEM3 ArrayBase<T,S,A>::ArrayBase(T * src, uint32_t sz)
:	ARRAYBASE_INIT
{	source(src, sz); }

TEM3 ArrayBase<T,S,A>::ArrayBase(const ArrayBase<T,S,A>& src)
:	ARRAYBASE_INIT
{	source(src); }

#undef ARRAYBASE_INIT

TEM3 ArrayBase<T,S,A>::~ArrayBase(){ clear(); }

TEM3 inline ArrayBase<T,S,A>& ArrayBase<T,S,A>::assign(const T& v){
	return assign(v, size());
}

TEM3 inline ArrayBase<T,S,A>& ArrayBase<T,S,A>::assign(
	const T& v, uint32_t end, uint32_t stride, uint32_t start
){
	for(uint32_t i=start; i<end; i+=stride) A::construct(mElems+i, v);
	return *this;
}

TEM3 inline T& ArrayBase<T,S,A>::operator[](uint32_t i){ return elems()[i]; }
TEM3 inline const T& ArrayBase<T,S,A>::operator[](uint32_t i) const { return elems()[i]; }

TEM3 inline T * ArrayBase<T,S,A>::elems() const { return mElems; }

TEM3 void ArrayBase<T,S,A>::clear(){ //printf("ArrayBase::clear(): mElems=%p\n", mElems);
	
	// Is this memory being ref counted?
	// NOTE: this check is only necessary because mElems can be assigned from
	// a raw external pointer. Maybe having this functionality is not a good
	// idea in the first place???
	// The culprit method is source(T * src, uint32_t size).
	if(mElems && managing(mElems)){
		int& c = refCount()[mElems];
		--c;
		if(0 == c){
			refCount().erase(mElems);
			for(uint32_t i=0; i<size(); ++i) A::destroy(mElems+i);
			A::deallocate(mElems, size());
		}
		mElems=0; mSize(0);
	}
}

TEM3 void ArrayBase<T,S,A>::own(){	
	T * oldElems = elems();
	
	// Check to see if we are already the owner
	if(references(oldElems) != 1){
		uint32_t oldSize = size();
		clear();
		resize(oldSize);
		for(uint32_t i=0; i<size(); ++i) A::construct(mElems+i, oldElems[i]);
	}
}

TEM3 void ArrayBase<T,S,A>::resize(uint32_t newSize, const T& c){
//printf("ArrayBase::resize() %p\n", this);
	newSize = mSize.convert(newSize);

	if(0 == newSize){
		clear();
	}

	if(newSize != size()){
		
		T * newElems = A::allocate(newSize);

		// If successful allocation...
		if(newElems){

			uint32_t nOldToCopy = newSize<size() ? newSize : size();
		
			// Copy over old elements
			for(uint32_t i=0; i<nOldToCopy; ++i){
				A::construct(newElems+i, (*this)[i]);
			}
			
			// Copy argument into any additional elements
			for(uint32_t i=nOldToCopy; i<newSize; ++i){
				A::construct(newElems+i, c);
			}
		
			clear();
			mElems = newElems;
		
			refCount()[mElems] = 1;
			mSize(newSize);
			onResize();
		}
	}
	//printf("ArrayBase::resize(): mElems=%p, size=%d\n", mElems, size());
}

TEM3 inline uint32_t ArrayBase<T,S,A>::size() const { return mSize(); }

TEM3 void ArrayBase<T,S,A>::source(const ArrayBase<T,S,A>& src){
	source(src.elems(), src.size());
}

TEM3 void ArrayBase<T,S,A>::source(T * src, uint32_t size){
	clear();
	if(managing(src)){
		++refCount()[src];
	}
	mElems = src;
	mSize(size);
}

#undef TEM3


#define TEM template<class T, class A>

// ArrayPow2

TEM inline uint32_t ArrayPow2<T,A>::oneIndex() const { return 1<<fracBits(); }
TEM inline uint32_t ArrayPow2<T,A>::log2Size() const { return mSize.mBitsI; }
TEM inline uint32_t ArrayPow2<T,A>::fracBits() const { return mSize.mBitsF; }
TEM inline uint32_t ArrayPow2<T,A>::index(uint32_t phase) const { return phase >> fracBits(); }

TEM inline const T& ArrayPow2<T,A>::atPhase(uint32_t phase) const { return (*this)[index(phase)]; }
TEM inline void ArrayPow2<T,A>::putPhase(uint32_t phase, T v){ (*this)[index(phase)] = v; }

TEM inline float ArrayPow2<T,A>::fraction(uint32_t phase) const{		
	return gam::fraction(log2Size(), phase);
}


// Buffer

TEM Buffer<T,A>::Buffer(uint32_t size)
	: Array<T>(size)
{
	unlock();
}

TEM inline void Buffer<T,A>::writeLock(const T * src){
	if(isLocked()) return;
	lock();
	memcpy(this->elems(), src, this->size() * sizeof(T));
}

TEM inline void Buffer<T,A>::writeLock(const T * src, uint32_t numRead, uint32_t writeOffset){
	if(isLocked()) return;
	lock();
	memcpy(this->elems() + writeOffset, src, numRead * sizeof(T));
}

TEM inline void Buffer<T,A>::unlock(){ mLocked = false; }
TEM inline void Buffer<T,A>::lock(){ mLocked = true; }
TEM inline bool Buffer<T,A>::isLocked(){ return mLocked; }
TEM inline bool Buffer<T,A>::isUnlocked(){ return !mLocked; }



//---- Ring

TEM Ring<T,A>::Ring(uint32_t size, const T& v) : Array<T,A>(size,v), mPos(size-1){}

TEM inline void Ring<T,A>::operator()(const T& v){
	incPos();				// inc write pos; do first to avoid out-of-bounds access
	(*this)[pos()] = v;		// write new element
}

TEM void Ring<T,A>::copy(T * dst, uint32_t len, uint32_t delay) const{
	// pos() points to most recently written slot
	//uint32_t tap = (pos() - delay) % size();
	uint32_t tap = (uint32_t)scl::wrap((int32_t)pos() - (int32_t)delay, (int32_t)size());

	// this ensures that we don't copy across the ring tap boundary
	// we add one to maxLen because of a fence post anomaly
	uint32_t maxLen = (tap < pos() ? (pos() - tap) : (pos() + (size() - tap))) + 1;
	len = scl::min(len, maxLen);
	
	mem::copyFromRing(elems(), size(), tap, dst, len);
}

TEM void Ring<T,A>::copyUnwrap(T * dst, uint32_t len) const { copy(dst, len, size() - 1); }

TEM inline uint32_t Ring<T,A>::indexBack() const {
	uint32_t i = pos() + 1;
	return (i != size()) ? i : 0;
}

TEM inline uint32_t Ring<T,A>::indexFront() const { return pos(); }

TEM inline uint32_t Ring<T,A>::indexPrev(uint32_t v) const {
	return scl::wrapOnce<int>(pos() - v, size());
}

TEM inline uint32_t Ring<T,A>::pos() const { return mPos; }

TEM inline void Ring<T,A>::incPos(){ if(++mPos >= size()) mPos = 0; }
TEM void Ring<T,A>::pos(uint32_t index){ mPos = index; }

TEM void Ring<T,A>::reset(){ pos(size()-1); }

TEM	inline void Ring<T,A>::writeClip(const T& v){
	if(mPos < size()){
		(*this)[mPos] = v;
		mPos++;
	}
}


// DoubleBuffer

TEM DoubleBuffer<T,A>::DoubleBuffer(uint32_t n)
: Array<T>(n*2), mSizeSingle(n)
{
	mBack  = Base::mElems;
	mFront = Base::mElems + n;
}

TEM inline bool DoubleBuffer<T,A>::backIsFirst(){ return mBack==Base::mElems; }

TEM inline T * DoubleBuffer<T,A>::back(){ return mBack; }
TEM inline T * DoubleBuffer<T,A>::front(){ return mFront; }
TEM inline void DoubleBuffer<T,A>::swap(){
	if(backIsFirst()){
		mFront= Base::mElems;
		mBack = Base::mElems + mSizeSingle;
	}
	else{
		mBack = Base::mElems;
		mFront= Base::mElems + mSizeSingle;	
	}
}

TEM void DoubleBuffer<T,A>::copyAll(T * dst){
	if(backIsFirst()){
		memcpy(dst, mBack, (mSizeSingle<<1) * sizeof(T));
	}
	else{
		memcpy(dst, mBack, mSizeSingle * sizeof(T));
		memcpy(dst + mSizeSingle, mFront, mSizeSingle * sizeof(T));
	}
}

TEM void DoubleBuffer<T,A>::onResize(){
	mSizeSingle = Base::size()/2;
	mBack  = Base::mElems;
	mFront = Base::mElems + mSizeSingle;
}

//TEM inline void DoubleBuffer<T>::write(const T * singleBuffer){
//	memcpy(mBuffer + mTap, singleBuffer, mSizeSingle * sizeof(T));
//	swap();
//}

} // gam::

#undef TEM

#endif
