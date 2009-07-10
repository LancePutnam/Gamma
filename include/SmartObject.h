#ifndef GAMMA_SMARTOBJECT_H_INC
#define GAMMA_SMARTOBJECT_H_INC

#include <list>

namespace gam{

/*	
We would like to prevent these types of situations from happening:
	Vec3 * v1 = new Vec3;
	Vec3 * v2 = v1;
	delete v1;
	// v2 is now dangling!

and:
	Vec3 * v1 = new Vec3;
	Vec3 * v2 = new Vec3;
	
	v1 = v2;
	// v1's memory leaked
*/

/// Mixin class for doing advanced memory management

/// The template parameter should be the base class in the derivation chain.
///
template <class BaseClass>
struct SmartObject{

	typedef std::list<void *> NewObjectList;

	SmartObject(): mDynamicAlloc(false){
		//printf("%x: ctor\n", this);
		// Check new list to see if this was dynamically allocated.
		// If so, then remove from list and set internal flag.
		NewObjectList& l = newObjects();
		for(NewObjectList::iterator it=l.begin(); it!=l.end(); ++it){
			
			// check if we lie within the memory footprint of the base class
			if(withinFootprint(*it)){
				mDynamicAlloc=true;
				l.erase(it);
				break;
			}
		}
	}

	void * operator new(size_t sz){
		void * m = malloc(sz);	// this will point to the base class
		newObjects().push_back(m);
		//printf("%x: new\n", m);
		return m;
	}

	void operator delete(void * m){
		if(m) free(m);
	}
	
	/// Returns true if the object was created dynamically with the new operator, false otherwise.
	bool dynamicAlloc() const { return mDynamicAlloc; }

private:
	bool mDynamicAlloc;
	
	// Returns true
	bool withinFootprint(void * m){
		BaseClass * b = (BaseClass *)m;
		return (this >= b) && (this < (b + sizeof(BaseClass)));
	}
	
	// Temporarily list of 'new' objects
	static NewObjectList& newObjects(){
		static NewObjectList * l = new NewObjectList;
		return *l;
	}
};


/*
// class for counted reference semantics
// - deletes the object to which it refers when the last CountedPtr
//   that refers to it is destroyed

template <class T>
class CountedPtr{
public:
	// initialize pointer with existing pointer
	// - requires that the pointer p is a return value of new
	explicit CountedPtr(T* p=0)
	: mPtr(p), mCount(new long(1)){}
	
	// copy pointer (one more owner)
	CountedPtr(const CountedPtr<T>& p)
	: mPtr(p.mPtr), mCount(p.mCount)
	{
		++*mCount;
	}
	
	// destructor (delete value if this was the last owner)
	~CountedPtr(){ dispose(); }
	
	// assignment (unshare old and share new value)
	CountedPtr<T>& operator= (const CountedPtr<T>& p){
		if(this != &p){
			dispose();
			mPtr = p.mPtr;
			mCount = p.mCount;
			++*mCount;
		}
		return *this;
	}
	
	// access the value to which the pointer refers
	T& operator *() const { return *mPtr; }
	T* operator->() const { return  mPtr; }
	
private:
	T* mPtr;        // pointer to the value
	long* mCount;   // shared number of owners

	void dispose() {
		if(0 == --*mCount){
			delete mCount;
			delete mPtr;
		}
	}
};
*/


} // gam::

#endif
