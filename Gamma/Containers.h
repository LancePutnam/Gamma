#ifndef GAMMA_CONTAINERS_H_INC
#define GAMMA_CONTAINERS_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information

	File Description:
	Dynamically sizable generic containers.
*/

/// \defgroup Containers

#include <stdlib.h>
#include <vector>
#include <map>

#include "Gamma/Allocator.h"
#include "Gamma/Conversion.h"
#include "Gamma/mem.h"
#include "Gamma/scl.h"

namespace gam{


///Size functor for ArrayPow2
    
///\ingroup Containers
struct SizeArrayPow2{
	SizeArrayPow2(uint32_t size){ (*this)(size); }
	uint32_t operator()() const { return (1<<mBitsI) & 0xfffffffe/*avoids 1*/; }
	void operator()(uint32_t v){ mBitsI = scl::log2(convert(v)); mBitsF = 32U - mBitsI; /*printf("%d %d\n", mBitsI, mBitsF);*/ }
	static uint32_t convert(uint32_t v){ v=scl::ceilPow2(v); return v!=1 ? v : 2; }	// should return 0,2,4,8,16,...
	uint32_t mBitsI;	// integer portion # bits
	uint32_t mBitsF;	// fraction portion # bits
};


/// Size functor for Array
    
/// \ingroup Containers
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
    
/// \ingroup Containers
    
template <class T, class S, class A=gam::Allocator<T> >
class ArrayBase : private A{
public:

	/// Default constructor that does not allocate memory
	ArrayBase();

	/// Constructor that allocates memory, but does not initialize elements

	/// \param[in] size		number of elements to allocate
	///
	explicit ArrayBase(uint32_t size);

	/// \param[in] size		number of elements to allocate
	/// \param[in] init		value to initialize all elements to
	ArrayBase(uint32_t size, const T& init);

	/// \param[in] src		external array to reference
	/// \param[in] size		size of external array
	ArrayBase(T * src, uint32_t size);

	/// \param[in] src		external array to reference
	explicit ArrayBase(ArrayBase<T,S,A>& src);

	virtual ~ArrayBase();


	/// Get write reference to element
	T& operator[](uint32_t i);
	
	/// Get read-only reference to element
	const T& operator[](uint32_t i) const;
	
	/// Sets all elements to value
	ArrayBase& assign(const T& v);

	/// Assign elements from another array
	template <class Arr>
	ArrayBase& assign(const Arr& src);

	/// Sets linear slice of elements to value
	
	/// \param[in] v		value to be copied as new content
	/// \param[in] end		end index (exclusive)
	/// \param[in] stride	index stride amount
	/// \param[in] start	start index (inclusive)
	ArrayBase& assign(const T& v, uint32_t end, uint32_t stride=1, uint32_t start=0);


	T * elems();					///< Get writable pointer to elements	
	const T * elems() const;		///< Get read-only pointer to elements
	uint32_t size() const;			///< Returns number of elements in array

	/// Destroys all elements and frees memory
	void clear();

	/// Ensures ownership of elements
	
	/// If the array is not already the sole owner, new memory is allocated and 
	/// the previously referenced array elements are copied.
	void own();
	
	/// Returns true if we are the sole owner of data internally allocated
	bool isSoleOwner() const;
	
	bool usingExternalSource() const;
	
	/// Resizes number of elements in array
	
	/// If the new size is less than the old size, then elements are truncated.
	/// If the new size is greater than the old size, then the argument value
	/// is copied into the additional elements.
	void resize(uint32_t newSize, const T& c=T());
	void source(ArrayBase<T,S,A>& src);		///< Sets source of array elements to another array
	void source(T * src, uint32_t size);	///< Sets source of array elements to another array

	/// Called whenever the size changes
	virtual void onResize(){}

	/// Returns number of pointers to memory address being managed
	static int references(T* m){
		typename RefCount::const_iterator it = refCount().find(m);
		return it != refCount().end() ? it->second : 0;
	}

protected:
	T * mElems;
	S mSize;

	typedef std::map<T *, int> RefCount;

	static RefCount& refCount(){
		static RefCount * o = new RefCount;
		return *o;
	}
	
	// is memory being managed automatically?
	static bool managing(T* m){ return references(m) != 0; }

private: ArrayBase& operator=(const ArrayBase& v);
};



/// Resizable array
    
/// \ingroup Containers
template <class T, class A=gam::Allocator<T> >
class Array : public ArrayBase<T, SizeArray, A>{
public:
	typedef ArrayBase<T, SizeArray, A> Base;

	Array(): Base(){}
	explicit Array(uint32_t size): Base(size){}
	Array(uint32_t size, const T& init): Base(size, init){}
	Array(T * src, uint32_t size): Base(src, size){}
	Array(Array& src): Base(src){}

	virtual ~Array(){}

private: Array& operator=(const Array& v);
};



///Resizable array with a power-of-2 number of elements
    
///\ingroup Containers
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



/// Ring buffer
    
/// \ingroup Containers
template <class T, class A=gam::Allocator<T> >
class Ring : public Array<T,A> {
public:

	typedef Array<T,A> Base; using Base::elems; using Base::size;

	/// \param[in]	size		Number of elements in ring.
	/// \param[in]	value		Initial value of all elements.
	explicit Ring(uint32_t size=0, const T& value=T());

	/// Returns reference to backmost (oldest) element
	T& readBack(){ return (*this)[indexBack()]; }
	const T& readBack() const { return (*this)[indexBack()]; }
	
	/// Returns reference to frontmost (newest) element
	T& readFront(){ return (*this)[indexFront()]; }
	const T& readFront() const { return (*this)[indexFront()]; }
	
	/// Returns reference to element 'ago' indices behind front
	T& read(uint32_t ago){ return (*this)[indexPrev(ago)]; }
	const T& read(uint32_t ago) const { return (*this)[indexPrev(ago)]; }

	/// Copies len elements starting from element pos() - delay into dst.
	void copy(T * dst, uint32_t len, uint32_t delay) const;
	
	/// Copy elements starting from last in into dst unwrapping from ring
	void copyUnwrap(T * dst, uint32_t len) const;
	
	uint32_t pos() const;			///< Return absolute index of frontmost (newest) element
	bool reachedEnd() const;		///< Returns whether the last element written was at the end of the array
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



/// Double buffered ring-buffer

/// This is a two-part buffer consisting of a ring buffer for writing and
/// a standard (absolute indexed) array for reading.
template <class T, class A=gam::Allocator<T> >
class DoubleRing : public Ring<T,A>{
public:
	/// @param[in]	size		Number of elements in ring.
	/// @param[in]	value		Initial value of all elements.
	explicit DoubleRing(uint32_t size=0, const T& value=T())
	:	Ring<T>(size, value), mRead(size)
	{}

	/// Returns reference to the reading buffer
	const Array<T,A>& readBuf() const { return mRead; }

	/// Copy elements into read buffer unwrapping from ring
	
	/// \returns a pointer to the read buffer
	///
	T * copyUnwrap(){ Ring<T,A>::copyUnwrap(mRead.elems(), mRead.size()); return mRead.elems(); }
	
	/// Copy elements into read buffer "as is" from ring
	
	/// \returns a pointer to the read buffer
	///
	T * copy(){
		mem::deepCopy(mRead.elems(), Ring<T,A>::elems(), mRead.size());
		//for(uint32_t i=0; i<read.size(); ++i) construct(read.elems()+i, (*this)[i]);
		return mRead.elems();
	}

	/// Resize buffers
	void resize(int n){ Ring<T,A>::resize(n); mRead.resize(n); }

protected:	
	Array<T,A> mRead;
};



/// N-sample delay
    
/// \ingroup Containers
template <class T, class A=gam::Allocator<T> >
struct DelayN: public Ring<T,A>{
	using Ring<T,A>::incPos; using Ring<T,A>::pos;

	/// \param[in]	size		Delay size, greater than 0
	/// \param[in]	value		Initial value of all elements
	explicit DelayN(uint32_t size, const T& value=T())
	:	Ring<T,A>(size, value)
	{}

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

#define ARRAYBASE_INIT mElems(0), mSize(0)

template <class T, class S, class A>
ArrayBase<T,S,A>::ArrayBase()
:	ARRAYBASE_INIT{}

template <class T, class S, class A>
ArrayBase<T,S,A>::ArrayBase(uint32_t sz)
:	ARRAYBASE_INIT
{	resize(sz); }

template <class T, class S, class A>
ArrayBase<T,S,A>::ArrayBase(uint32_t sz, const T& initial)
:	ARRAYBASE_INIT
{	resize(sz); assign(initial); }

template <class T, class S, class A>
ArrayBase<T,S,A>::ArrayBase(T * src, uint32_t sz)
:	ARRAYBASE_INIT
{	source(src, sz); }

template <class T, class S, class A>
ArrayBase<T,S,A>::ArrayBase(ArrayBase<T,S,A>& src)
:	ARRAYBASE_INIT
{	source(src); }

#undef ARRAYBASE_INIT

template <class T, class S, class A>
ArrayBase<T,S,A>::~ArrayBase(){ clear(); }

template <class T, class S, class A>
ArrayBase<T,S,A>& ArrayBase<T,S,A>::assign(const T& v){
	return assign(v, size());
}

template <class T, class S, class A>
template <class Arr>
ArrayBase<T,S,A>& ArrayBase<T,S,A>::assign(const Arr& src){
	unsigned N = src.size();
	if(N > size()) N = size();
	for(unsigned i=0; i<N; ++i) A::construct(mElems+i, T(src[i]));
	return *this;
}

template <class T, class S, class A>
ArrayBase<T,S,A>& ArrayBase<T,S,A>::assign(
	const T& v, uint32_t end, uint32_t stride, uint32_t start
){
	for(uint32_t i=start; i<end; i+=stride) A::construct(mElems+i, v);
	return *this;
}

template <class T, class S, class A>
inline T& ArrayBase<T,S,A>::operator[](uint32_t i){ return elems()[i]; }
template <class T, class S, class A>
inline const T& ArrayBase<T,S,A>::operator[](uint32_t i) const { return elems()[i]; }

template <class T, class S, class A>
inline T * ArrayBase<T,S,A>::elems(){ return mElems; }
template <class T, class S, class A>
inline const T * ArrayBase<T,S,A>::elems() const { return mElems; }

template <class T, class S, class A>
void ArrayBase<T,S,A>::clear(){ //printf("ArrayBase::clear(): mElems=%p\n", mElems);
	
	// We will only attempt to deallocate the data if it exists and is being 
	// managed (reference counted) by ArrayBase.
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

template <class T, class S, class A>
void ArrayBase<T,S,A>::own(){
	T * oldElems = elems();
	
	// If we are not the sole owner, do nothing...
	if(!isSoleOwner()){
		uint32_t oldSize = size();
		clear();
		resize(oldSize);
		for(uint32_t i=0; i<size(); ++i) A::construct(mElems+i, oldElems[i]);
	}
}

template <class T, class S, class A>
bool ArrayBase<T,S,A>::isSoleOwner() const {
	return references((T*)elems()) == 1;
}

template <class T, class S, class A>
bool ArrayBase<T,S,A>::usingExternalSource() const {
	return elems() && !managing((T*)elems());
}

template <class T, class S, class A>
void ArrayBase<T,S,A>::resize(uint32_t newSize, const T& c){
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

template <class T, class S, class A>
inline uint32_t ArrayBase<T,S,A>::size() const { return mSize(); }

template <class T, class S, class A>
void ArrayBase<T,S,A>::source(ArrayBase<T,S,A>& src){
	source(src.elems(), src.size());
}

template <class T, class S, class A>
void ArrayBase<T,S,A>::source(T * src, uint32_t size){
	if(src == mElems) return; // check for self assignment
	clear();
	if(managing(src)){
		++refCount()[src];
	}
	mElems = src;
	mSize(size);
	onResize();
}

// ArrayPow2

template<class T, class A>
inline uint32_t ArrayPow2<T,A>::oneIndex() const { return 1<<fracBits(); }

template<class T, class A>
inline uint32_t ArrayPow2<T,A>::log2Size() const { return mSize.mBitsI; }

template<class T, class A>
inline uint32_t ArrayPow2<T,A>::fracBits() const { return mSize.mBitsF; }
    
template<class T, class A>
inline uint32_t ArrayPow2<T,A>::index(uint32_t phase) const { return phase >> fracBits(); }

template<class T, class A>
inline const T& ArrayPow2<T,A>::atPhase(uint32_t phase) const { return (*this)[index(phase)]; }

template<class T, class A>
inline void ArrayPow2<T,A>::putPhase(uint32_t phase, T v){ (*this)[index(phase)] = v; }

template<class T, class A>
inline float ArrayPow2<T,A>::fraction(uint32_t phase) const{
	return gam::fraction(log2Size(), phase);
}



//---- Ring

template<class T, class A>
Ring<T,A>::Ring(uint32_t size, const T& v) : Array<T,A>(size,v), mPos(size-1){}

template<class T, class A>
inline void Ring<T,A>::operator()(const T& v){
	incPos();				// inc write pos; do first to avoid out-of-bounds access
	(*this)[pos()] = v;		// write new element
}

template<class T, class A>
void Ring<T,A>::copy(T * dst, uint32_t len, uint32_t delay) const{
	// pos() points to most recently written slot
	//uint32_t tap = (pos() - delay) % size();
	uint32_t tap = (uint32_t)scl::wrap((int32_t)pos() - (int32_t)delay, (int32_t)size());

	// this ensures that we don't copy across the ring tap boundary
	// we add one to maxLen because of a fence post anomaly
	uint32_t maxLen = (tap < pos() ? (pos() - tap) : (pos() + (size() - tap))) + 1;
	len = scl::min(len, maxLen);
	
	mem::copyFromRing(elems(), size(), tap, dst, len);
}

template<class T, class A>
void Ring<T,A>::copyUnwrap(T * dst, uint32_t len) const { copy(dst, len, size() - 1); }

template<class T, class A>
inline uint32_t Ring<T,A>::indexBack() const {
	uint32_t i = pos() + 1;
	return (i != size()) ? i : 0;
}

template<class T, class A>
inline uint32_t Ring<T,A>::indexFront() const { return pos(); }

template<class T, class A>
inline uint32_t Ring<T,A>::indexPrev(uint32_t v) const {
	return scl::wrapOnce<int>(pos() - v, size());
}

template<class T, class A>
inline uint32_t Ring<T,A>::pos() const { return mPos; }
    
template<class T, class A>
inline bool Ring<T,A>::reachedEnd() const { return pos() == (size()-1); }

template<class T, class A>
inline void Ring<T,A>::incPos(){ if(++mPos >= size()) mPos = 0; }
    
template<class T, class A>
void Ring<T,A>::pos(uint32_t index){ mPos = index; }

template<class T, class A>
void Ring<T,A>::reset(){ pos(size()-1); }

template<class T, class A>
inline void Ring<T,A>::writeClip(const T& v){
	if(mPos < size()){
		(*this)[mPos] = v;
		mPos++;
	}
}

} // gam::
#endif
