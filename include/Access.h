#ifndef GAMMA_ACCESS_H_INC
#define GAMMA_ACCESS_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information

	File Description: 
	Functions/objects for accessing and indexing arrays.
*/

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

inline void neighborsWrap(int x1, int N, int& x0, int& x2){
	x0 = x1-1; if(x1< 0) x1+=N;
	x2 = x1+1; if(x2>=N) x2-=N;
}

/// Maps a position in [-1, 1] to an index in [0, n). No boundary operations are performed.
inline int posToInd(float v, int n){ return n * (v*0.49999f + 0.5f); }


// Neighbor accessing strategies
namespace acc{

	struct None{
		static int checkM1(int i, int mx, int mn){ return i; }
		static int checkP1(int i, int mx, int mn){ return i; }
	};

	struct Wrap{
		static int checkM1(int i, int mx, int mn){ return i<mn ? mx : i; }
		static int checkP1(int i, int mx, int mn){ return i>mx ? mn : i; }
	};

	struct Clip{
		static int checkM1(int i, int mx, int mn){ return i>mn ? i : mn; }
		static int checkP1(int i, int mx, int mn){ return i<mx ? i : mx; }
	};

//	struct Fold{
//		static int checkM1(int i, int mx, int mn){ return i>=mn ? i : mn+1; }
//		static int checkP1(int i, int mx, int mn){ return i<=mx ? i : mx-1; }
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
		i1=checkM1(begin-mStride);
		i2=begin;
	}

	void operator()(){ i0=i1; i1=i2; i2=checkP1(i2+mStride); }
	void operator()(int& a0, int& a1, int& a2){	(*this)(); a0=i0; a1=i1; a2=i2; }
	bool valid(int i){ return (i<mMin) || (i>mMax); }

	int i0, i1, i2;

private:
	int mMax, mMin, mStride;
	
	int checkM1(int i){ return Tacc::checkM1(i, mMax, mMin); }
	int checkP1(int i){ return Tacc::checkP1(i, mMax, mMin); }
};



template <class T>
class IndexMap{
public:
	IndexMap(int indMax, const T& posMax){ max(indMax, posMax); }
	
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
	int cast(const T& v) const { return castIntTrunc(v); }
	T cast(int v) const { return T(v); }
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



// 1-D arithmetic series index generator
// for(Series1 i(N); i();)
//struct Series1{
//	Series1(int end, int stride=1, int start=0)
//	:	i(start-stride), n(end-start), mEnd(end), mStride(stride){}
//
//	bool operator()(){ i += mStride; return i < mEnd; }
//	
//	int flat() const { return i; }
//	double frac() const { return double(i)/n; }
//	double fracS() const { return double(i<<1)/n - 1.; }
//	
//	operator int(){ return flat(); }
//	
//	int i, n;
//private:
//	int mEnd, mStride;
//};

//struct Series3{
//	Series1(int end1, int end2, int end3, int stride=1)
//	:	i(start-stride), n(end-start), mEnd(end), mStride(stride){}
//
//	bool operator()(){ i += mStride; return i < mEnd; }
//	
//	int flat() const { return i; }
//	double frac() const { return double(i)/n; }
//	double fracS() const { return double(i<<1)/n - 1.; }
//	
//	operator int(){ return flat(); }
//	
//private:
//	int mEnd, mStride;
//};


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


} // gam::

#endif
