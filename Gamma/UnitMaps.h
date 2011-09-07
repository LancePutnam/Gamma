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



// Fixed-size power of 2 table supporting fixed point lookup

// This table minimizes memory usage and table look-up speed at the expense
// of having a fixed size that is a power of two.
template <uint32_t B, class T>
class TablePow2{
public:

	enum{ N = 1<<B };

	const T& operator[](unsigned i) const { return mElems[i]; }
	T& operator[](unsigned i){ return mElems[i]; }

	/// Read value using truncating interpolation
	
	/// @param[in] phase	phase value in [0,1)
	///
	const T& read(double phase) const {
		return read(phaseR2I(phase));
	}

	/// Read value using truncating interpolation
	const T& read(uint32_t phase) const {
		return (*this)[phase >> shift()];
	}

	/// Read value using linear interpolation
	T readL(uint32_t phase) const {
		const T& lo = read(phase);
		const T& hi = read(phase + oneIndex());
		float fr = gam::fraction(bits(), phase);
		return lo + (hi - lo)*fr;
	}

	static uint32_t mask(){ return N-1; }
	static uint32_t size(){ return N; }
	static uint32_t bits(){ return B; }
	static uint32_t shift(){ return 32U-B; }
	static uint32_t oneIndex(){ return 1<<shift(); }

protected:
	T mElems[N];
	static uint32_t phaseR2I(double v){ return static_cast<uint32_t>(v * 4294967296.); }
};



// Complex sinusoid table
// B is the log2 size of each table
// D is the number of tables
// The effective table size is (2^B)^D
template <unsigned B=10, unsigned D=2, class T=double>
class CSinTable{
public:
	typedef TablePow2<B, Complex<T> > Arc;

	CSinTable(){ init(); }
	
	/// Get sinusoidal value at unit phase. No bounds checking is performed.
	Complex<T> operator()(double phase){
		return (*this)(uint32_t(phase * 4294967296.));
	}
	
	/// Get value from fixed-point phase in interval [0, 2^(B*D))
	Complex<T> operator()(uint32_t p){

		p >>= shift();

		// start with finest sample
		Complex<T> r = arc(D-1)[p & Arc::mask()];

		// iterate remaining arcs (fine to coursest)
		for(unsigned i=2; i<D+1; ++i){
			p >>= B;
			r *= arc(D-i)[p & Arc::mask()];
		}
		return r;
	}
	
	static uint32_t shift(){ return 32U - (B*D); }

private:
	enum{ M = Arc::N };

	// Get pointer to complex arc tables
	static Arc * arcs(){
		static Arc tables[D]; // higher index = finer resolution
		return tables;
	}
	
	static Arc& arc(unsigned res){ return arcs()[res]; }
	
	static void init(){
		static bool tabulate=true;
		if(tabulate){
			tabulate=false;
			unsigned long long N = M;
			
			for(unsigned j=0; j<D; ++j){		// iterate resolution (course to fine)
				for(unsigned i=0; i<M; ++i){	// iterate arc of unit circle
					double p = (i*M_PI)/(N>>1);
					arc(j)[i] = Polar<T>(1, p);
				}
				N *= M;
			}
			
		}
	}
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
