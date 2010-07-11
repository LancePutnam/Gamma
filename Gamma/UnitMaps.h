#ifndef GAMMA_UNITMAPS_H_INC
#define GAMMA_UNITMAPS_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include <math.h>
#include <map>
#include <vector>
#include "Gamma/scl.h"
#include "Gamma/Access.h"
#include "Gamma/Strategy.h"
#include "Gamma/Types.h"

namespace gam{

namespace MapType{
	/// Mapping types.
	enum {	Lin,		/**< b0 + u * (b1-b0)		*/
			Pow,		/**< b0 + u^p * (b1-b0)		*/
			Exp2		/**< p 2^[b0 + u * (b1-b0)]	*/
	};
};



/// Tabulated function with real number lookup.
template <class T, template<class> class Sipl=ipl::Linear, class Sacc=acc::Wrap, class A=gam::Allocator<T> >
class FunctionTable : public Array<T,A>{
	typedef Array<T,A> Base;
public:
	using Base::elems; using Base::size;

	//explicit FunctionTable(): Base(){}

	/// Constructor that allocates an internal table

	/// @param[in] size		Number of elements (actual number is power of 2 ceiling)
	/// @param[in] init		Initializes all elements to this value
	explicit FunctionTable(uint32_t size, const T& init=T(0))
	:	Base(size, init)
	{
		endpoints(0, size);
	}
	
	virtual ~FunctionTable(){}
	

	/// Returns f(x) where x lies in the function domain, [0,1).
	T operator()(double x) const {

		x = scl::wrap(x);
		double f;
//		index_t i1 = mIndMap(x, f);
//		return mIpl(mAcc, *this, i1,f, size()-1);
		index_t i1 = mIndMap(x, f) + mInterval.min();
		return mIpl(mAcc, *this, i1,f, mInterval.max(), mInterval.min());

		//int i2 = i1+1; if(i2==size()) i2=0;
		//return (*this)[i1]*(1.f-f) + (*this)[i2]*f;
	}

	/// Sums generator stream with table elements
	template <class Gen>
	FunctionTable& operator+=(Gen& g){
		for(uint32_t i=0; i<size(); ++i) (*this)[i] += g();
		return *this;
	}

	template <class Gen>
	FunctionTable& operator+=(const Gen& g){
		for(uint32_t i=0; i<size(); ++i) (*this)[i] += g();
		return *this;
	}
	
	/// Set indexing interval for look-up [min, max)
	FunctionTable& endpoints(index_t min, index_t max){
		mInterval.endpoints(min, max-1); // interpolator max index is inclusive
		mIndMap.max(max-min, 1.);
		return *this;
	}

protected:

	virtual void onResize(){ mIndMap.max(size(), 1.); }
	
	IndexMap<double> mIndMap;
	Interval<index_t> mInterval;

	Sipl<T> mIpl;
	Sacc mAcc;
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

template <class T> UnitMapper<T>::UnitMapper(){
	set((T)1);
}

template <class T> UnitMapper<T>::UnitMapper(T bound1, T bound0, T p1, int type, bool clip){
	set(bound1, bound0, p1, type, clip);
}

template <class T> void UnitMapper<T>::set(T bound1, T bound0, T p1, int type, bool clip){
	this->bound1 = bound1;
	this->bound0 = bound0;
	this->p1 = p1;
	this->type = type;
	this->clip = clip;
}

template <class T> inline T UnitMapper<T>::map(T u){
	switch(type){
	case MapType::Lin:	return mapLin(u);
	case MapType::Pow:	return mapPow(u);
	case MapType::Exp2:	return mapExp2(u);
	default:;
	}
}

template <class T> T UnitMapper<T>::unmap(T value){
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


template <class T> T UnitMapper<T>::mapLin(T u){
	doClip(u); return bound0 + u * (bound1 - bound0);
}

template <class T> T UnitMapper<T>::mapPow(T u){
	doClip(u); return (T)scl::mapPower(u, bound1, bound0, p1);
}

template <class T> T UnitMapper<T>::mapExp2(T u){
	doClip(u);
	return (T)(pow(2., scl::mapPower(u, bound1, bound0, 1.)) * p1);
}

} // gam::

#endif
