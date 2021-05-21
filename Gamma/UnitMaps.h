#ifndef GAMMA_UNITMAPS_H_INC
#define GAMMA_UNITMAPS_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include <functional>
#include <initializer_list>
#include <map>
#include <vector>
#include "Gamma/scl.h"
#include "Gamma/Access.h"
#include "Gamma/Constants.h"
#include "Gamma/Strategy.h"
#include "Gamma/Types.h"

namespace gam{


/// Maps value in unit interval to a tabulated function
template<
	class T,
	template<class> class Sipl=ipl::Linear,
	class Sacc=acc::Wrap,
	class A=gam::Allocator<T>
>
class LookupTable : public Array<T,A>{
	typedef Array<T,A> Base;
public:
	using Base::elems; using Base::size;


	/// Constructor that allocates an internal table

	/// \param[in] size		Number of elements (actual number is power of 2 ceiling)
	/// \param[in] init		Initial value of elements
	explicit LookupTable(unsigned size=2048, const T& init=T(0))
	:	Base(size, init)
	{
		endpoints(0, size);
	}
	
	virtual ~LookupTable(){}
	
	/// Get array
	Array<T,A>& array(){ return *this; }
	const Array<T,A>& array() const { return *this; }

	/// Returns f(x) where x lies in the function domain, [0,1)
	T operator()(double x) const {

		x = scl::wrap(x);
		double f;
//		index_t i1 = mIndMap(x, f);
//		return mIpl(mAcc, *this, i1,f, size()-1);
		index_t i1 = mIndMap(x, f) + mIMin;
		return mIpl(mAcc, elems(), i1,f, mIMax, mIMin);

		//int i2 = i1+1; if(i2==size()) i2=0;
		//return (*this)[i1]*(1.f-f) + (*this)[i2]*f;
	}


	/// Set indexing interval for look-up [min, max)
	LookupTable& endpoints(index_t min, index_t max){
		if(min < max){
			mIMin = min;
			mIMax = max-1; // make index inclusive
		}else{
			mIMin = max;
			mIMax = min-1; // make index inclusive
		}
		mIndMap.max((mIMax+1)-mIMin, 1.);
		return *this;
	}


	/// Sums generator stream with table elements
	template <class Gen>
	LookupTable& operator+=(Gen& g){
		for(unsigned i=0; i<size(); ++i) (*this)[i] += g();
		return *this;
	}

	template <class Gen>
	LookupTable& operator+=(const Gen& g){
		for(unsigned i=0; i<size(); ++i) (*this)[i] += g();
		return *this;
	}

protected:
	IndexMap<double> mIndMap;
	index_t mIMin, mIMax; // min and max indices, inclusive
	Sipl<T> mIpl;
	Sacc mAcc;

	virtual void onResize(){ mIndMap.max(size(), 1.); }
};



/// Fixed-size power of 2 table supporting fixed point lookup

/// This table minimizes memory usage and table look-up speed at the expense
/// of having a fixed size that is a power of two.
///
/// \tparam B	log2 size of table (size == 2^B)
/// \tparam T	element type
template <uint32_t B, class T>
class TablePow2{
public:

	enum{ N = 1<<B };

	/// Construct and leave elements uninitialized
	TablePow2(){}
	
	/// Construct and assign value to all elements
	TablePow2(const T& fillValue){
		for(auto& v : *this) v = fillValue;
	}
	
	/// Construct and fill table according to function
	
	/// @param[in] x01	Unit position in table, in [0,1)
	///
	TablePow2(const std::function<T(float x01)>& fillFunc){
		assign(fillFunc);
	}

	/// Construct from initializer list
	TablePow2(std::initializer_list<T> vals){
		int len = vals.size() < N ? vals.size() : N;
		for(int i=0; i<len; ++i)
			mElems[i] = vals.begin()[i];
	}


	/// Assign elements from function
	
	/// @param[in] x01	Unit position in table, in [0,1)
	///
	TablePow2& assign(const std::function<T(float x01)>& f){
		for(int i=0; i<N; ++i) mElems[i] = f(float(i)/N);
		return *this;
	}

	/// Assign elements from another array
	template <class U>
	TablePow2& assign(const U * src){
		for(int i=0; i<N; ++i) mElems[i] = src[i];
		return *this;
	}

	const T& operator[](unsigned i) const { return mElems[i]; }
	T& operator[](unsigned i){ return mElems[i]; }

	const T * begin() const { return mElems; }
	T * begin(){ return mElems; }
	const T * end() const { return mElems+N; }
	T * end(){ return mElems+N; }


	/// Read value using truncating interpolation
	
	/// \param[in] phase	phase value in [0,1)
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

	static constexpr uint32_t mask(){ return N-1; }
	static constexpr uint32_t size(){ return N; }
	static constexpr uint32_t bits(){ return B; }
	static constexpr uint32_t shift(){ return 32U-B; }
	static constexpr uint32_t oneIndex(){ return 1<<shift(); }

protected:
	T mElems[N];
	static uint32_t phaseR2I(double v){ return static_cast<uint32_t>(v * 4294967296.); }
};



/// Complex sinusoid lookup table

/// This retrieves a complex sinusoid e^it given a phase t. D tables each with a
/// size of 2^B are used to compute the value. D-1 tables are used for recursive
/// interpolation. This makes the lookup equivalent to a non-interpolating
/// lookup in a table of size (2^B)^D.
///
/// \tparam B	log2 size of each table (size == 2^B)
/// \tparam D	the number of tables
template <unsigned B=10, unsigned D=2, class TComplex=Complex<double> >
class CSinTable{
public:
	typedef TablePow2<B, TComplex> Arc;

	CSinTable(){ init(); }
	
	/// Get sinusoidal value at unit phase. No bounds checking is performed.
	TComplex operator()(double phase){
		return (*this)(uint32_t(phase * 4294967296.));
	}
	
	/// Get value from fixed-point phase in interval [0, 2^(B*D))
	TComplex operator()(uint32_t p){

		p >>= shift();

		// start with finest sample
		TComplex r(arc(D-1)[p & Arc::mask()]);

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
					arc(j)[i].real() = cos(p);
					arc(j)[i].imag() = sin(p);
				}
				N *= M;
			}
			
		}
	}
};



/// Maps a unit value through an invertible function
template <class T>
class UnitMapper{
public:

	/// Mapping types
	enum MapType{
		MAP_LIN,		/**< b0 + u * (b1-b0)		*/
		MAP_POW,		/**< b0 + u^p * (b1-b0)		*/
		MAP_EXP2		/**< p 2^[b0 + u * (b1-b0)]	*/
	};


	T min, max, p1;
	MapType type;
	bool clip;


	UnitMapper();
	
	/// \param[in] max		upper endpoint of interval
	/// \param[in] min		lower endpoint of interval
	/// \param[in] p1		mapping function function parameter
	/// \param[in] type		mapping function
	/// \param[in] clip		whether to clip values to interval
	UnitMapper(T max, T min=0., T p1=1., MapType type = MAP_POW, bool clip=true);


	/// Set all attributes
	UnitMapper& set(T max, T min=0., T p1=1., MapType type = MAP_POW, bool clip=true);
	
	T map(T unit);			///< Map a unit value
	T unmap(T value);		///< Unmap a value to a unit value
	
private:
	T mapLin (T u);
	T mapPow (T u);	// Map normal directly using power function
	T mapExp2(T u);	// Map normal directly using exponentiation function
	
	// clip input normal (or not)
	void doClip(T& nrm){ if(clip) nrm = scl::clip(nrm); }
};




// Implementation ______________________________________________________________

template <class T> UnitMapper<T>::UnitMapper(){
	set(T(1));
}

template <class T> UnitMapper<T>::UnitMapper(T max, T min, T p1, MapType type, bool clip){
	set(max, min, p1, type, clip);
}

template <class T> UnitMapper<T>& UnitMapper<T>::set(T max, T min, T p1, MapType type, bool clip){
	this->max = max;
	this->min = min;
	this->p1 = p1;
	this->type = type;
	this->clip = clip;
	return *this;
}

template <class T> inline T UnitMapper<T>::map(T u){
	switch(type){
	case MAP_LIN:	return mapLin(u);
	case MAP_POW:	return mapPow(u);
	case MAP_EXP2:	return mapExp2(u);
	default:;
	}
}

template <class T> T UnitMapper<T>::unmap(T v){
	switch(type){
	case MAP_LIN:
		return scl::mapLin(v, min, max, T(0), T(1));

	case MAP_POW:
		v = scl::mapLin(v, min, max, T(0), T(1));
		return pow(v, 1. / p1);
		//return mapPow(pow(v, 1. / p1));

	case MAP_EXP2:
		v = log(v / p1) * M_LOG2E; // log2(x)
		return scl::mapLin(v, min, max, T(0), T(1));
		//v = scl::mapLin(v, min, max, T(0), T(1));
		//return mapExp2(value);
	default:
		return 0;
	}
}


template <class T> T UnitMapper<T>::mapLin(T u){
	doClip(u);
	return min + u * (max - min);
}

template <class T> T UnitMapper<T>::mapPow(T u){
	doClip(u);
	return (T)scl::mapPower(u, max, min, p1);
}

template <class T> T UnitMapper<T>::mapExp2(T u){
	doClip(u);
	return (T)(pow(2., scl::mapPower(u, max, min, 1.)) * p1);
}

} // gam::

#endif
