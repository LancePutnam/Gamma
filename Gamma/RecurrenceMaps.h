#ifndef GAMMA_RECURRENCEMAPS_H_INC
#define GAMMA_RECURRENCEMAPS_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include "Gamma/mem.h"
#include "Gamma/scl.h"

namespace gam{


/// Base recursive map
template <uint32_t Nv, uint32_t Nc, class Tv>
struct MapBase{

	MapBase(){
		mem::deepZero(v, numVals());
		mem::deepZero(a, numCoefs());
	}
	
	Tv v[Nv];								///< Values
	Tv a[Nc];								///< Coefficients

	static uint32_t numCoefs(){ return Nc; }	///< Returns number of coefficients
	static uint32_t numVals (){ return Nv; }	///< Returns number of values
};


#define INHERIT(der, nv, nc)\
using MapBase<nv,nc,T>::v; using MapBase<nv,nc,T>::a;\
der() : MapBase<nv,nc,T>(){}


/// 1-D linear map.

/// Iterated map of function:
/// \code
/// v0 = a0 + a1 v0
/// \endcode
template <class T=gam::real>
struct MapLin : public MapBase<1,2,T>{

	INHERIT(MapLin, 1, 2)

	MapLin(T v0, T a0=0, T a1=0){ v[0]=v0; coef(a0,a1); }

	void coef(T a0, T a1){ a[0] = a0; a[1] = a1; }

	T operator()(T max, T min=0){ return v[0] = scl::wrap((*this)(), max, min); }

	T operator()(){
		return v[0] = a[0] + v[0] * a[1];
	}
};


/// 2-D linear map.

/// Iterated map of function:
/// \code
/// v0 = a0 + a1 v0 + a2 v1
/// v1 = a3 + a4 v0 + a5 v1
/// \endcode
template <class T=gam::real>
struct MapLin2 : public MapBase<2,6,T>{

	INHERIT(MapLin2, 2, 6)

	MapLin2(T v0, T v1){ v[0]=v0; v[1]=v1; }

	T operator()(){
		T v0 = v[0];
		v[0] = a[0] + v0 * a[1] + v[1] * a[2];
		v[1] = a[3] + v0 * a[4] + v[1] * a[5];
		return v[0];
	}
};


/// 3-D linear map.

/// Iterated map of function:
/// \code
/// v0 = a0 + a1 v0 +  a2 v1 +  a3 v2
/// v1 = a4 + a5 v0 +  a6 v1 +  a7 v2
/// v2 = a8 + a9 v0 + a10 v1 + a11 v2
/// \endcode
template <class T=gam::real>
struct MapLin3 : public MapBase<3,12,T>{

	INHERIT(MapLin3, 3, 12)

	MapLin3(T v0, T v1, T v2){ v[0]=v0; v[1]=v1; v[2]=v2; }

	void operator()(){
		T v0 = v[0];
		T v1 = v[1];
		v[0] = a[ 0] + v0 * a[ 1] + v1 * a[ 2] + v[2] * a[ 3];
		v[1] = a[ 4] + v0 * a[ 5] + v1 * a[ 6] + v[2] * a[ 7];
		v[2] = a[ 8] + v0 * a[ 9] + v1 * a[10] + v[2] * a[11];
	}
};


/// 4-D linear map.

/// Iterated map of function:
/// \code
/// v0 =  a0 +  a1 v0 +  a2 v1 +  a3 v2 +  a4 v3
/// v1 =  a5 +  a6 v0 +  a7 v1 +  a8 v2 +  a9 v3
/// v2 = a10 + a11 v0 + a12 v1 + a13 v2 + a14 v3
/// v3 = a15 + a16 v0 + a17 v1 + a18 v2 + a19 v3
/// \endcode
template <class T=gam::real>
struct MapLin4 : public MapBase<4,20,T>{

	INHERIT(MapLin4, 4, 20)

	MapLin4(T v0, T v1, T v2, T v3){ v[0]=v0; v[1]=v1; v[2]=v2; v[3]=v3; }

	void operator()(){
		T v0=v[0]; T v1=v[1]; T v2=v[2];
		v[0] = a[ 0] + v0 * a[ 1] + v1 * a[ 2] + v2 * a[ 3] + v[3] * a[ 4];
		v[1] = a[ 5] + v0 * a[ 6] + v1 * a[ 7] + v2 * a[ 8] + v[3] * a[ 9];
		v[2] = a[10] + v0 * a[11] + v1 * a[12] + v2 * a[13] + v[3] * a[14];
		v[3] = a[15] + v0 * a[16] + v1 * a[17] + v2 * a[18] + v[3] * a[19];
	}
};


/// 5-D linear map.

/// Iterated map of function:
/// \code
/// v0 =  a0 +  a1 v0 +  a2 v1 +  a3 v2 +  a4 v3 +  a5 v4
/// v1 =  a6 +  a7 v0 +  a8 v1 +  a9 v2 + a10 v3 + a11 v4
/// v2 = a12 + a13 v0 + a14 v1 + a15 v2 + a16 v3 + a17 v4
/// v3 = a18 + a19 v0 + a20 v1 + a21 v2 + a22 v3 + a23 v4
/// v4 = a24 + a25 v0 + a26 v1 + a27 v2 + a28 v3 + a29 v4
/// \endcode
template <class T=gam::real>
struct MapLin5 : public MapBase<5,30,T>{

	INHERIT(MapLin5, 5, 30)

	MapLin5(T v0, T v1, T v2, T v3){ v[0]=v0; v[1]=v1; v[2]=v2; v[3]=v3; }

	void operator()(){
		T v0=v[0]; T v1=v[1]; T v2=v[2]; T v3=v[3];
		v[0] = a[ 0] + v0 * a[ 1] + v1 * a[ 2] + v2 * a[ 3] + v3 * a[ 4] + v[4] * a[ 5];
		v[1] = a[ 6] + v0 * a[ 7] + v1 * a[ 8] + v2 * a[ 9] + v3 * a[10] + v[4] * a[11];
		v[2] = a[12] + v0 * a[13] + v1 * a[14] + v2 * a[15] + v3 * a[16] + v[4] * a[17];
		v[3] = a[18] + v0 * a[19] + v1 * a[20] + v2 * a[21] + v3 * a[22] + v[4] * a[23];
		v[4] = a[24] + v0 * a[25] + v1 * a[26] + v2 * a[27] + v3 * a[28] + v[4] * a[29];
	}
};


/// 1-D quadratic map.

/// Iterated map of function:
/// \code
/// v0 = a0 + a1 v0 + a2 v0^2
/// \endcode
template <class T=gam::real>
struct MapQuad : public MapBase<1,3,T>{

	INHERIT(MapQuad, 1, 3)

	MapQuad(T v0){ v[0]=v0; }

	void operator()(){
		v[0] = a[0] + v[0] * (a[1] + v[0] * a[2]);
	}
};


/// 2-D quadratic map.

/// Iterated map of function:
/// \code
/// v0 = a0 + a1 v0 + a2 v0^2 + a3 v0 v1 +  a4 v1  +  a5 v1^2
/// v1 = a6 + a7 v0 + a8 v0^2 + a9 v0 v1 + a10 v1  + a11 v1^2
/// \endcode
template <class T=gam::real>
struct MapQuad2 : public MapBase<2,12,T>{

	INHERIT(MapQuad2, 2, 12)

	MapQuad2(T v0, T v1){ v[0]=v0; v[1]=v1; }

	void operator()(){
		T v0 = v[0];
		v[0] = a[0] + v0 * (a[1] + v0 * a[2] + v[1] * a[3]) + v[1] * (a[ 4] + v[1] * a[ 5]);
		v[1] = a[6] + v0 * (a[7] + v0 * a[8] + v[1] * a[9]) + v[1] * (a[10] + v[1] * a[11]);
	}
};


/// 3-D quadratic map.

/// Iterated map of function:
/// \code
/// v0 =  a0 +  a1 v0 +  a2 v0^2 +  a3 v0 v1 +  a4 v0 v2 +  a5 v1 +  a6 v1^2 +  a7 v1 v2 +  a8 v2 +  a9 v2^2
/// v1 = a10 + a11 v0 + a12 v0^2 + a13 v0 v1 + a14 v0 v2 + a15 v1 + a16 v1^2 + a17 v1 v2 + a18 v2 + a19 v2^2
/// v2 = a20 + a21 v0 + a22 v0^2 + a23 v0 v1 + a24 v0 v2 + a25 v1 + a26 v1^2 + a27 v1 v2 + a28 v2 + a29 v2^2
/// \endcode
template <class T=gam::real>
struct MapQuad3 : public MapBase<3,30,T>{

	INHERIT(MapQuad3, 3, 30)

	MapQuad3(T v0, T v1, T v2){ v[0]=v0; v[1]=v1; v[2]=v2; }

	void operator()(){
		T v0 = v[0];
		T v1 = v[1];
		v[0] = a[ 0] + v0 * (a[ 1] + v0 * a[ 2] + v1 * a[ 3] + v[2] * a[ 4]) + v1 * (a[ 5] + v1 * a[ 6] + v[2] * a[ 7]) + v[2] * (a[ 8] + v[2] * a[ 9]);
		v[1] = a[10] + v0 * (a[11] + v0 * a[12] + v1 * a[13] + v[2] * a[14]) + v1 * (a[15] + v1 * a[16] + v[2] * a[17]) + v[2] * (a[18] + v[2] * a[19]);
		v[2] = a[20] + v0 * (a[21] + v0 * a[22] + v1 * a[23] + v[2] * a[24]) + v1 * (a[25] + v1 * a[26] + v[2] * a[27]) + v[2] * (a[28] + v[2] * a[29]);
	}
};

#undef INHERIT

} // gam::

#endif
