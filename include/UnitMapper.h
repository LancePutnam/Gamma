#ifndef GAMMA_UNITMAPPER_H_INC
#define GAMMA_UNITMAPPER_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include <math.h>
#include <map>
#include <vector>
#include "scl.h"
#include "MacroD.h"

namespace gam{

namespace MapType{
	/// Mapping types.
	enum {	Lin,		/**< b0 + u * (b1-b0)		*/
			Pow,		/**< b0 + u^p * (b1-b0)		*/
			Exp2		/**< p 2^[b0 + u * (b1-b0)]	*/
	};
};

// Maps a normalized value to a warped, ranged value.
template <class T>
class UnitMapper{
public:
	UnitMapper();
	UnitMapper(T bound1, T bound0=0., T p1=1., int type = MapType::Pow, bool clip=true);

	T bound0, bound1, p1;
	int type;
	bool clip;

	void set(T bound1, T bound0=0., T p1=1., int type = MapType::Pow, bool clip=true);
	
	T map(T normal);			///< Map a unit value
	T unmap(T value);			///< Unmap a value to a unit value
	
private:
	T mapLin (T u);
	T mapPow (T u);	// Map normal directly using power function
	T mapExp2(T u);	// Map normal directly using exponentiation function
	
	// clip input normal (or not)
	void doClip(T& nrm){ if(clip) nrm = scl::clip(nrm); }
};




// Implementation ______________________________________________________________

TEM UnitMapper<T>::UnitMapper(){
	set((T)1);
}

TEM UnitMapper<T>::UnitMapper(T bound1, T bound0, T p1, int type, bool clip){
	set(bound1, bound0, p1, type, clip);
}

TEM void UnitMapper<T>::set(T bound1, T bound0, T p1, int type, bool clip){
	this->bound1 = bound1;
	this->bound0 = bound0;
	this->p1 = p1;
	this->type = type;
	this->clip = clip;
}

TEM inline T UnitMapper<T>::map(T u){
	switch(type){
	case MapType::Lin:	return mapLin(u);
	case MapType::Pow:	return mapPow(u);
	case MapType::Exp2:	return mapExp2(u);
	default:;
	}
}

TEM T UnitMapper<T>::unmap(T value){
	switch(type){
	case MapType::Lin:
		return scl::mapLin(value, bound0, bound1, (T)0, (T)1);

	case MapType::Pow:
		value = scl::mapLin(value, bound0, bound1, (T)0, (T)1);
		return pow(value, 1. / p1);
		//return mapPow(pow(value, 1. / p1));

	case MapType::Exp2:
		value = log2(value / p1);
		return scl::mapLin(value, bound0, bound1, (T)0, (T)1);
		//value = scl::mapLin(value, bound0, bound1, (T)0, (T)1);
		//return mapExp2(value);
	default: return 0;
	}
}


TEM T UnitMapper<T>::mapLin(T u){
	doClip(u); return bound0 + u * (bound1 - bound0);
}

TEM T UnitMapper<T>::mapPow(T u){
	doClip(u); return (T)scl::mapNormal(u, bound1, bound0, p1);
}

TEM T UnitMapper<T>::mapExp2(T u){
	doClip(u);
	return (T)(pow(2., scl::mapNormal(u, bound1, bound0, 1.)) * p1);
}

} // gam::

#include "MacroU.h"

#endif
