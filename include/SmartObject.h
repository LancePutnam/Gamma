#ifndef GAMMA_SMARTOBJECT_H_INC
#define GAMMA_SMARTOBJECT_H_INC

#include <list>

namespace gam{

/// Mixin class for doing advanced memory management

/// The template parameter should be the base class in the derivation chain.
///
template <class BaseClass>
struct SmartObject {

	typedef std::list<void *> NewObjectList;

	SmartObject(): mDynamicAlloc(false){
		//printf("%x: ctor\n", this);
		// Check new list to see if this was dynamically allocated.
		// If so, then remove from list and set internal flag.
		NewObjectList& l = newObjects();
		for(NewObjectList::iterator it=l.begin(); it!=l.end(); ++it){
			
			// check if we lie within the memory footprint of the base class
			if(withinFootPrint(*it)){
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
		free(m);
	}
	
	/// Returns true if the object was created dynamically with the new operator, false otherwise.
	bool dynamicAlloc() const { return mDynamicAlloc; }

private:
	bool mDynamicAlloc;
	
	// Returns true
	bool withinFootPrint(void * m){
		BaseClass * b = (BaseClass *)m;
		return (this >= b) && (this < (b + sizeof(BaseClass)));
	}
	
	// Temporarily list of 'new' objects
	static NewObjectList& newObjects(){
		static NewObjectList * l = new NewObjectList;
		return *l;
	}
};


} // gam::

#endif
