#ifndef GAMMA_ACCESS_H_INC
#define GAMMA_ACCESS_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information

	File Description: 
	Functions/objects for accessing and indexing arrays.
*/

#include "pstdint.h"

namespace gam{


/// Compute 2-D array indices from 1-D array index
template <class T>
inline void index1to2(T index1, T sizeX, T& x, T& y){
	y = index1 / sizeX; x = index1 % sizeX;
}

/// Compute 1-D array index from 3-D array indices

/// The x indices move fastest followed by y, then z.
///
template <class T>
inline T index3to1(T x, T y, T z, T sizeX, T sizeY){
	return x + sizeX * (y + sizeY * z);
}

/// Returns last index of an arithmetic iteration starting from 0
inline uint32_t indexLast(uint32_t len, uint32_t str){ return ((len-1)/str)*str; }

inline void neighborsWrap(int x1, int N, int& x0, int& x2){
	x0 = x1-1; if(x1< 0) x1+=N;
	x2 = x1+1; if(x2>=N) x2-=N;
}

/// Maps a position in [-1, 1] to an index in [0, n). No boundary operations are performed.
inline int posToInd(float v, int n){ return n * (v*0.49999f + 0.5f); }


// Neighbor accessing strategies
// Valid index interval is [mn, mx]. The max value is inclusive to speed up
// checks of closest neighbors, which is done most often.
namespace acc{

	struct None{
		static int mapM1(int v, int mx, int mn){ return v; }
		static int mapP1(int v, int mx, int mn){ return v; }
		static int map  (int v, int mx, int mn){ return v; }
	};

	struct Wrap{
		static int mapM1(int v, int mx, int mn){ return v<mn ? mx : v; }
		static int mapP1(int v, int mx, int mn){ return v>mx ? mn : v; }
		static int map  (int v, int mx, int mn){ return v<mn ? v+mx+1-mn : (v>mx ? v-(mx+1-mn) : v); }
	};

	struct Clip{
		static int mapM1(int v, int mx, int mn){ return v>mn ? v : mn; }
		static int mapP1(int v, int mx, int mn){ return v<mx ? v : mx; }
		static int map  (int v, int mx, int mn){ return v>mn ? (v<mx ? v : mx) : mn; }
	};

//	struct Fold{
//		static int checkM1(int i, int mx, int mn){ return i>=mn ? i : mn+1; }
//		static int mapP1(int i, int mx, int mn){ return i<=mx ? i : mx-1; }
//	};
};


// First iteration will give (begin-1, begin, begin+1) with specified 
// neighbor accessing strategy.
template <class Tacc=acc::Wrap>
class AccessStream1{
public:

	/// @param[in] size		size of dimension
	/// @param[in] min		minimum index
	/// @param[in] begin	relative beginning index in [0, size)
	AccessStream1(int size, int min=0, int stride=1, int begin=0)
	:	i0(0), i1(0), i2(0), mMax((size-1)*stride + min), mMin(min), mStride(stride)
	{
		begin += min;
		i1=mapM1(begin-mStride);
		i2=begin;
	}

	void operator()(){ i0=i1; i1=i2; i2=mapP1(i2+mStride); }
	void operator()(int& a0, int& a1, int& a2){	(*this)(); a0=i0; a1=i1; a2=i2; }
	bool valid(int i){ return (i<mMin) || (i>mMax); }

	int i0, i1, i2;

private:
	int mMax, mMin, mStride;
	
	int mapM1(int i){ return Tacc::mapM1(i, mMax, mMin); }
	int mapP1(int i){ return Tacc::mapP1(i, mMax, mMin); }
};



struct BoundaryClip{
	int operator()(int ix, int sizeX) const{
		return ix < 0 ? 0 : ix >= sizeX ? sizeX-1 : ix;
	}

	int operator()(int ix, int iy, int sizeX, int sizeY) const{		
		return (*this)(ix, sizeX) + (*this)(iy, sizeY)*sizeX;
	}
};

struct BoundaryWrap{
	int operator()(int ix, int sizeX) const{
		return ix < 0 ? sizeX+ix : ix >= sizeX ? ix-sizeX : ix;
	}

	int operator()(int ix, int iy, int sizeX, int sizeY) const{
		return (*this)(ix, sizeX) + (*this)(iy, sizeY)*sizeX;
	}
};



/// Maps a real number in [0, pmax) to an integer in [0, imax).
template <class T>
class IndexMap{
public:
	IndexMap(int indMax, const T& posMax=T(1)){ max(indMax, posMax); }
	
	int operator()(const T& x) const { return cast(x*mMul); }
	
	int operator()(const T& x, T& f) const {
		f = x*mMul;
		int i = cast(f);
		f -= cast(i); 
		return i;
	}
	
	T operator()(int i) const { return cast(i) * mRec; }

	void max(int indMax, const T& posMax){ mMul=indMax/posMax; mRec=1/mMul; }

private:
	T mMul, mRec;
	//int cast(const T& v) const { return castIntTrunc(v); }
	int cast(const T& v) const { return int(v); }	// use native cast, usually optimized
	T cast(int v) const { return T(v); }
};



template <int N=1, class Sbounds=BoundaryWrap>
struct Neighbors1D {

	Neighbors1D(int sizeX): sizeX(sizeX){}

	void operator()(int ix){
		for(int i=0; i<N; ++i){
			l[i] = bounds(ix-(i+1), sizeX);
			r[i] = bounds(ix+(i+1), sizeX);
		}
	}
	
	template <class T>
	T access(const T * src, int ix){ return src[bounds(ix, sizeX)]; }
	
	int l[N], r[N];

private:
	Sbounds bounds;
	int sizeX;
};



template <class Sneigh, class Sbounds=BoundaryWrap>
struct Neighbors2D : public Sneigh{

	Neighbors2D(int sizeX, int sizeY): sizeX(sizeX), sizeY(sizeY){}

	void operator()(int ix, int iy){
		(*this)(ix, iy, sizeX, sizeY, mBounds);
	}

private:
	Sbounds mBounds;
	int sizeX, sizeY;
};



template <class SBounds>
struct NeighborsCross2D{

	void operator()(int ix, int iy, int sizeX, int sizeY, const SBounds& sBounds){
		l = sBounds(ix-1,iy, sizeX, sizeY);
		r = sBounds(ix+1,iy, sizeX, sizeY);
		t = sBounds(ix,iy-1, sizeX, sizeY);
		b = sBounds(ix,iy+1, sizeX, sizeY);
	}

	union{
		struct{ int l, r, t, b; }; 
		int indices[4];
	};
};



template <class SBounds>
struct NeighborsDiag2D{

	void operator()(int ix, int iy, int sizeX, int sizeY, const SBounds& sBounds){
		tl = sBounds(ix-1,iy-1, sizeX, sizeY);
		br = sBounds(ix+1,iy+1, sizeX, sizeY);
		tr = sBounds(ix+1,iy-1, sizeX, sizeY);
		bl = sBounds(ix-1,iy+1, sizeX, sizeY);
	}

	union{
		struct{ int tl, br, tr, bl; }; 
		int indices[4];
	};
};





struct Scan1{
	Scan1(int size): n(size),i(-1){}

	bool operator()(){
		return ((++i)<n);
	}
	
	int flat() const { return i; }
	double frac() const { return double(i)/n; }
	double fracS() const { return double(i<<1)/n - 1.; }

	operator int() const { return flat(); }

	int n, i;
};

struct Scan2{
	Scan2(int size): s1(size),s2(size),i(-1),j(0){}
	Scan2(int size1, int size2): s1(size1),s2(size2),i(-1),j(0){}

	bool operator()(){
		if(++i==s1){
			i=0;
			if(++j==s2) return false;
		}
		return true;
	}
	
	int flat() const { return flat(i,j); }
	int flat(int i1, int i2) const { return i1 + i2*s1; }
	double frac1() const { return double(i)/s1; }
	double frac2() const { return double(j)/s2; }
	double frac1S() const { return double(i<<1)/s1 - 1.; }
	double frac2S() const { return double(j<<1)/s2 - 1.; }

	operator int() const { return flat(); }

	int s1, s2;
	int i, j;
};


struct Scan3{
	Scan3(int size): s1(size),s2(size),s3(size),i(-1),j(0),k(0){}
	Scan3(int size1, int size2, int size3): s1(size1),s2(size2),s3(size3),i(-1),j(0),k(0){}

	bool operator()(){
		if(++i==s1){
			i=0;
			if(++j==s2){
				j=0;
				if(++k==s3)	return false;
			}
		}
		return true;
	}
	
	int flat() const { return flat(i,j,k); }
	int flat(int i1, int i2, int i3) const { return index3to1(i1,i2,i3, s1,s2); }
	double frac1() const { return double(i)/s1; }
	double frac2() const { return double(j)/s2; }
	double frac3() const { return double(k)/s3; }
	double frac1S() const { return double(i<<1)/s1 - 1.; }
	double frac2S() const { return double(j<<1)/s2 - 1.; }
	double frac3S() const { return double(k<<1)/s3 - 1.; }

	operator int() const { return flat(); }

	int s1, s2, s3;
	union{
		struct{ int i, j, k; };
		struct{ int x, y, z; };
	};
};






/// Unbounded no-op index map
class Indexer{
public:
	Indexer(int32_t count_, double stride_=1.)
	:	mCount(count_)
	{}

	void count(int32_t v){ mCount=v; }
	
	int32_t operator()(int32_t i, int32_t max=-1) const { return i; }

	int32_t count() const { return mCount; }

protected:
	int32_t mCount;
};


/// Unbounded integer strided index map
class IndexerInt : public Indexer {
public:
	IndexerInt(int32_t count_, double stride_=1.)
	:	Indexer(count_), mStride(stride_)
	{}

	void stride(int32_t v){ mStride=v; }

	int32_t operator()(int32_t i, int32_t max=-1) const { return i*stride(); }

	int32_t stride() const { return mStride; }

protected:
	int32_t mStride;
};


/// Bounded real strided index map
class IndexerReal : public Indexer {
public:
	IndexerReal(int32_t count_, double stride_=1.)
	:	Indexer(count_), mEps(0)
	{ stride(stride_); }
	
	// This should be half of the smallest fractional component of stride
	void eps(double v){ mEps=v; }

	void stride(double v){ mStride=v; eps(0.); }
	
	void strides(int32_t s1, int32_t s2){ mStride=s1 + 1./s2; eps(0.5/s2); }
	void strides(int32_t s1, int32_t s2, int32_t s3){ mStride=s1 + 1./s2 + 1./(s2*s1); eps(0.5/(s2*s1)); }

	int32_t operator()(int32_t i, int32_t max) const {	
		// conversion can fail on integer boundaries because fraction is not exact
		// add offset equal to 1/2 the fractional part of the stride...
		int32_t r = (i*stride() + mEps);
		
		// i*n1 + i/n2

		int32_t d = max;
		if(r>=max){ r-=max; if(r>=max){ r-=d*(r/d); } }
		else if(r<0){ r+=max; if(r<0){ r+=d*(1 - r/d); } }
		return r;
	}

	double stride() const { return mStride; }

protected:
	double mStride, mEps;
};




#define L1 for(int32_t i=0;i<count();++i)
#define L2 int32_t n=minCount(v); for(int32_t i=0;i<n;++i)


/// Uniformly strided section of an array

/// For operations between different slices, the minimum count between the
/// two slices will be used for iteration.
/// All operations requiring elements to be copied perform a shallow copy
/// (i.e. use '=' operator) and therefore are safe to use with objects.
template <class T>
class Slice{
public:

	/// @param[in] src		pointer to array elements
	/// @param[in] count	how many elements to iterate over
	/// @param[in] stride	stride increment through array
	/// @param[in] offset	absolute offset into array, -1 is last, -2 penultimate, etc.
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

/*

int sx = 64;
int sy = 32;
float data[sx*sy];

Neighbors2D<NeighborsCross2D, BoundaryWrap> n2d(sx, sy);

for(int j=0; j<sy; ++j){
	for(int i=0; i<sx; ++i){
		n2d(i, j);
		float d1xc = (data[n2d.r] - data[n2d.l])*0.5;
		float d1yc = (data[n2d.b] - data[n2d.t])*0.5;
	}
}

*/



//class Indices : public Indexer{
//public:
//
//	Indices(uint32_t maxSize)
//	:	Indexer(0), mMaxSize(maxSize)
//	{	mElems = new uint32_t[maxSize]; }
//
//	~Indices(){ delete[] mElems; }
//
//	// add value to indices
//	void operator+= (uint32_t v){
//		for(uint32_t i=0; i<end(); ++i)
//			mElems[i] = (mElems[i] + v) % mMaxSize;
//	}
//
//	// multiply indices by value
//	void operator*= (double v){
//		for(uint32_t i=0; i<end(); ++i)
//			mElems[i] = ((uint32_t)((double)mElems[i] * v + 0.5)) % mMaxSize;
//	}
//
//	// add new index
//	Indices& operator<< (uint32_t index){
//		if(end() < mMaxSize && index < mMaxSize) mElems[mEnd++] = index;
//		return *this;
//	}
//
//	void clear(){ mEnd=0; }
//	uint32_t * elems() const { return mElems; }
//
//	uint32_t index(uint32_t i) const { return mElems[i]; }
//
//private:
//	uint32_t mMaxSize;
//	uint32_t * mElems;
//};



// experimental index set
// similar to vector, but without auto-memory management
//class Indices{
//public:
//
//	Indices(uint maxSize)
//	:	mSize(0), mMaxSize(maxSize)
//	{	mElems = new uint[maxSize]; }
//
//	~Indices(){ delete[] mElems; }
//
//	// add value to indices
//	void operator+= (uint v){
//		for(uint i=0; i<size(); ++i)
//			mElems[i] = (mElems[i] + v) % mMaxSize;
//	}
//
//	// multiply indices by value
//	void operator*= (double v){
//		for(uint i=0; i<size(); ++i)
//			mElems[i] = ((uint)((double)mElems[i] * v + 0.5)) % mMaxSize;
//	}
//
//	// add new index
//	Indices& operator<< (uint index){
//		if(mSize < mMaxSize && index < mMaxSize) mElems[mSize++] = index;
//		return *this;
//	}
//
//	void clear(){ mSize=0; }
//	uint * elems() const { return mElems; }
//	uint size() const { return mSize; }
//
//	uint begin() const { return 0; }
//	uint end() const { return mSize; }
//	uint index(uint i) const { return mElems[i]; }
//	uint stride() const { return 1; }
//
//private:
//	uint mSize, mMaxSize;
//	uint * mElems;
//};



} // gam::

#endif
