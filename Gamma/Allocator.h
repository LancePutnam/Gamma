#ifndef GAMMA_ALLOCATOR_H_INC
#define GAMMA_ALLOCATOR_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information

	File Description:
	Interface for and default implementation of memory allocator
*/

#include <cstddef> // ptrdiff_t
#include <cstdlib> // size_t

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
	typedef std::size_t		size_type;
	typedef std::ptrdiff_t	difference_type;
	typedef T*				pointer;
	typedef const T*		const_pointer;
	typedef T&				reference;
	typedef const T&		const_reference;
	typedef T				value_type;
	template <class U> struct rebind { typedef Allocator<U> other; };

public:
	explicit Allocator(){}
	explicit Allocator(const Allocator&){}
	template <class U> explicit Allocator(const Allocator<U>&){}
	~Allocator(){}

	pointer address(reference x) const { return &x; }
	const_pointer address(const_reference x) const { return &x; }

	pointer allocate(size_type n){
		return reinterpret_cast<pointer>(::operator new(n * sizeof(T)));
	}

	void deallocate(pointer p, size_type /*n*/){ ::operator delete(p); }

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

} //gam::

#endif
