#ifndef GAMMA_ACCESS_H_INC
#define GAMMA_ACCESS_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information

	File Description: 
	Functions/objects for accessing and indexing arrays.
*/

/// \defgroup access Array access

#include "Gamma/pstdint.h"

namespace gam{

typedef int32_t index_t;

/// Returns last index of an arithmetic iteration starting from 0
inline uint32_t indexLast(uint32_t len, uint32_t str){ return ((len-1)/str)*str; }

/// Maps a position in [-1, 1] to an index in [0, n). No boundary operations are performed.
inline index_t posToInd(float v, index_t n){ return n * (v*0.49999f + 0.5f); }


// Neighbor accessing strategies
// Valid index interval is [mn, mx]. The max value is inclusive to speed up
// checks of closest neighbors, which is done most often.
namespace acc{

	struct None{
		static index_t mapM1(index_t v, index_t mx, index_t mn){ return v; }
		static index_t mapP1(index_t v, index_t mx, index_t mn){ return v; }
		static index_t map  (index_t v, index_t mx, index_t mn){ return v; }
	};

	struct Wrap{
		static index_t mapM1(index_t v, index_t mx, index_t mn){ return v<mn ? mx : v; }
		static index_t mapP1(index_t v, index_t mx, index_t mn){ return v>mx ? mn : v; }
		static index_t map  (index_t v, index_t mx, index_t mn){ return v<mn ? v+mx+1-mn : (v>mx ? v-(mx+1-mn) : v); }
	};

	struct Clip{
		static index_t mapM1(index_t v, index_t mx, index_t mn){ return v>mn ? v : mn; }
		static index_t mapP1(index_t v, index_t mx, index_t mn){ return v<mx ? v : mx; }
		static index_t map  (index_t v, index_t mx, index_t mn){ return v>mn ? (v<mx ? v : mx) : mn; }
	};

//	struct Fold{
//		static index_t checkM1(index_t i, index_t mx, index_t mn){ return i>=mn ? i : mn+1; }
//		static index_t mapP1(index_t i, index_t mx, index_t mn){ return i<=mx ? i : mx-1; }
//	};
};



/// Maps a real number in [0, pmax) to an integer in [0, imax).
    
///\ingroup access    
template <class T>
class IndexMap{

public:
	IndexMap(index_t idxMax=1, const T& posMax=T(1)){ max(idxMax, posMax); }
	
	index_t operator()(const T& x) const { return cast(x*mMul); }
	
	index_t operator()(const T& x, T& f) const {
		f = x*mMul;
		index_t i = cast(f);
		f -= cast(i); 
		return i;
	}
	
	T operator()(index_t i) const { return cast(i) * mRec; }

	void max(index_t idxMax, const T& posMax){ mMul=idxMax/posMax; mRec=1/mMul; }

private:
	T mMul, mRec;
	//index_t cast(const T& v) const { return castIntTrunc(v); }
	index_t cast(const T& v) const { return index_t(v); }	// use native cast, usually optimized
	T cast(index_t v) const { return T(v); }
};




#define L1 for(int32_t i=0;i<count();++i)
#define L2 int32_t n=minCount(v); for(int32_t i=0;i<n;++i)


/// Uniformly strided section of an array

/// For operations between different slices, the minimum count between the
/// two slices will be used for iteration.
/// All operations requiring elements to be copied perform a shallow copy
/// (i.e. use '=' operator) and therefore are safe to use with objects.
///  \ingroup access
template <class T>
class Slice{
public:

	/// \param[in] src		pointer to array elements
	/// \param[in] count_	how many elements to iterate over
	/// \param[in] stride_	stride increment through array
	/// \param[in] offset_	absolute offset into array, -1 is last, -2 penultimate, etc.
	Slice(T * src, int32_t count_, int32_t stride_=1, int32_t offset_=0)
	:	A(src), C(count_), S(stride_)
	{	offset(offset_);	}

	/// Returns new sub-slice
	Slice operator()(int32_t cnt, int32_t str=1, int32_t off=0) const { return Slice(A, cnt,str,off); }
	
	/// Returns ith count element
	T& operator[](int32_t i) const { return B[i*S]; }

	template <class Gen>
	const Slice& operator  = (const Gen& v) const { L1{ (*this)[i] =v(); } return *this; }
	const Slice& operator  = (const   T& v) const { L1{ (*this)[i] =v  ; } return *this; }

	template <class Gen>
	bool operator == (const Gen& v) const { L1{ if(v() != (*this)[i]) return false; } return true; }
	bool operator == (const   T& v) const { L1{ if(v   != (*this)[i]) return false; } return true; }

	template <class U>
	const Slice& operator += (const Slice<U>& v) const { L2{ (*this)[i]+=T(v[i]); } return *this; }

	template <class Gen>
	const Slice& operator += (const Gen& v) const { L1{ (*this)[i]+=v(); } return *this; }
	const Slice& operator += (const   T& v) const { L1{ (*this)[i]+=v  ; } return *this; }

	template <class U>
	const Slice& operator -= (const Slice<U>& v) const { L2{ (*this)[i]-=T(v[i]); } return *this; }

	template <class Gen>
	const Slice& operator -= (const Gen& v) const { L1{ (*this)[i]-=v(); } return *this; }
	const Slice& operator -= (const   T& v) const { L1{ (*this)[i]-=v  ; } return *this; }

	template <class U>
	const Slice& operator *= (const Slice<U>& v) const { L2{ (*this)[i]*=T(v[i]); } return *this; }

	template <class Gen>
	const Slice& operator *= (const Gen& v) const { L1{ (*this)[i]*=v(); } return *this; }
	const Slice& operator *= (const   T& v) const { L1{ (*this)[i]*=v  ; } return *this; }

	template <class U>
	const Slice& operator /= (const Slice<U>& v) const { L2{ (*this)[i]/=T(v[i]); } return *this; }

	template <class Gen>
	const Slice& operator /= (const Gen& v) const { L1{ (*this)[i]/=v(); } return *this; }
	const Slice& operator /= (const   T& v) const { L1{ (*this)[i]/=v  ; } return *this; }

	/// Copy elements from another slice.
	
	/// Source elements are statically cast to the type of the destination slice.
	///
	template <class U>
	const Slice& copy(const Slice<U>& v) const { L2{ (*this)[i]=v[i]; } return *this; }

	/// Apply filter in-place
	template <class Fil>
	const Slice& filter(const Fil& v) const { L1{ (*this)[i]=v((*this)[i]); } return *this; }
	
	/// Apply C-style unary function in-place, x = func(x, a1)
	template <class R, class X, class A1>
	const Slice& filter(R (* const func)(X, A1), const A1& a1){
		L1{ (*this)[i] = func((*this)[i], a1); }
		return *this;
	}

	/// Apply C-style binary function in-place, x = func(x, a1, a2)
	template <class R, class X, class A1, class A2>
	const Slice& filter(R (* const func)(X, A1,A2), const A1& a1, const A2& a2){
		L1{ (*this)[i] = func((*this)[i], a1,a2); }
		return *this;
	}

	/// Reverse slice
	const Slice& reverse(){ B=B+(C-1)*S; S=-S; return *this; }
	
	/// Returns reversed slice
	Slice reversed() const { Slice r=*this; r.B=B+(C-1)*S; r.S=-r.S; return r; }

	/// Set all elements to argument
	const Slice& set(const T& v=T()) const { return (*this = v); }

	/// Swaps elements
	template <class U>
	const Slice& swap(const Slice<U>& v) const {
		L2{ T t=(*this)[i]; (*this)[i]=v[i]; v[i]=t; }
		return *this;
	}
	
	/// Returns mean of elements
	T mean() const { return sum()/C; }

	/// Returns sum of elements in slice
	T sum() const { T r=T(0); L1{ r+=(*this)[i]; } return r; }

	int32_t count() const { return C; }
	int32_t offset() const { return B-A; }
	int32_t stride() const { return S; }
	int32_t N() const { return (S>0?S:-S)*C; }

	Slice& count(int32_t v){ C=v; return *this; }
	Slice& offset(int32_t v){ B=A+(v<0 ? N()+v : v); return *this; }
	Slice& stride(int32_t v){ S=v; return *this; }

protected:
	T * A, * B;		// absolute, relative pointers
	int32_t C,S;	// count, stride
	int32_t minCount(const Slice& o) const { return count()<o.count() ? count() : o.count(); }
};

#undef L1
#undef L2

/// Slice object function
template <class T>
Slice<T> slice(T * src, int32_t cnt, int32_t str=1, int32_t off=0){ return Slice<T>(src,cnt,str,off); }


/*
1 2 3 4 5 6 7 8
1 5 2 6 3 7 4 8		d = 4 + 1/2
1 3 5 7 2 4 6 8		d = 2 + 1/4

1 2 3 4 5 6 7 8 9
1 4 7 2 5 8 3 6 9	d = 3 + 1/3

*/

} // gam::

#endif
