#ifndef GAMMA_ALLOCATOR_H_INC
#define GAMMA_ALLOCATOR_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information

	File Description:
	Interface for and default implementation of memory allocator
*/

#include <stdlib.h>		/* size_t */

namespace gam{


template <class T> class Allocator;

// specialize for void:
template<> class Allocator<void> {
public:
  typedef void*       pointer;
  typedef const void* const_pointer;
  // reference to void members are impossible.
  typedef void value_type;
  template <class U> struct rebind { typedef Allocator<U> other; };
};

template <class T> class Allocator{
public:
	typedef size_t    size_type;
	typedef ptrdiff_t difference_type;
	typedef T*        pointer;
	typedef const T*  const_pointer;
	typedef T&        reference;
	typedef const T&  const_reference;
	typedef T         value_type;
	template <class U> struct rebind { typedef Allocator<U> other; };

public:
	explicit Allocator(){}
	explicit Allocator(const Allocator&){}
	template <class U> explicit Allocator(const Allocator<U>&){}
	~Allocator(){}

	pointer address(reference x) const { return &x; }
	const_pointer address(const_reference x) const { return &x; }

	pointer allocate(size_type n, Allocator<void>::const_pointer hint = 0){
		return reinterpret_cast<pointer>(::operator new(n * sizeof(T)));
	}

	void deallocate(pointer p, size_type n){ ::operator delete(p); }

	size_type max_size() const {
		return static_cast<size_type>(-1) / sizeof(T);
	}

	void construct(pointer p, const T& val){ new(p) T(val); }
	void destroy(pointer p){ p->~T(); }
};

template <class T1, class T2>
bool operator==(const Allocator<T1>&, const Allocator<T2>&){ return true; }

template <class T1, class T2>
bool operator!=(const Allocator<T1>&, const Allocator<T2>&){ return false; }



//template <class T, class Alloc=Allocator<T> >
//class Buffer : private Alloc{
//
//	explicit Buffer(int n, const T& v=T(), const Alloc& a=Alloc())
//	: Alloc(a), mMemBegin(0), mMemEnd(0)
//	{
//		resize(n);
//	}
//	
//	~Buffer(){
//		clear();
//	}
//	
//	
//	
//	
//	void clear(){
//	
//		
//	
////		if(mMemBegin != mMemEnd){
////			for(T * i = mMemBegin; i<mMemEnd; ++i)
////				Alloc::destroy(i);
////			Alloc::deallocate(mMemBegin, mMemEnd - mMemBegin);
////		}	
//	}
//	
//	void resize(int n){
//		if((mMemEnd - mMemBegin) < n){
//			mMemBegin = Alloc::allocate(n);
//			for(int i=0; i<n; ++i){
//				Alloc::construct(mMemBegin + i, v);
//			}		
//		}
//	}
//	
//	void reserve(int n){
//		
//	}
//
//protected:
//	T * mMemBegin;		// First element in memory block
//	T * mMemEnd;		// 1 past last allocated element
//	T * mEnd;			// 1 past last element in buffer
//	
//};

// vector


/*

Full ISO C++ standard (20.4.1)

namespace std {
	template <class T> class allocator;

	// specialize for void:
	template <> class allocator<void> {
	public:
		typedef void*       pointer;
		typedef const void* const_pointer;
		// reference to void members are impossible.
		typedef void value_type;
		template <class U> struct rebind { typedef allocator<U> other; };
	};

	template <class T> class allocator {
	public:
		typedef size_t    size_type;
		typedef ptrdiff_t difference_type;
		typedef T*        pointer;
		typedef const T*  const_pointer;
		typedef T&        reference;
		typedef const T&  const_reference;
		typedef T         value_type;
		template <class U> struct rebind { typedef allocator<U> other; };

		allocator() throw();
		allocator(const allocator&) throw();
		template <class U> allocator(const allocator<U>&) throw();
		~allocator() throw();

		pointer address(reference x) const;
		const_pointer address(const_reference x) const;

		pointer allocate(size_type, allocator<void>::const_pointer hint = 0);
		void deallocate(pointer p, size_type n);
		size_type max_size() const throw();

		void construct(pointer p, const T& val);
		void destroy(pointer p);
	};
	
	template <class T1, class T2>
	bool operator==(const allocator<T1>&, const allocator<T2>&) throw();

	template <class T1, class T2>
	bool operator!=(const allocator<T1>&, const allocator<T2>&) throw();
	
}
  
  
  
  
*/

} //gam::

#endif
