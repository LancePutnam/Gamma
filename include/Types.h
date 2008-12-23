#ifndef GAMMA_TYPES_H_INC
#define GAMMA_TYPES_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include "pstdint.h"	/* for uint32_t, uint16_t, etc... */
#include <math.h>
#include <stdlib.h>

namespace gam{

typedef unsigned long	uint;	// default natural number type
typedef float			real;	// default real number type



/// Complex number
template <class T=gam::real>
struct Complex{

	typedef Complex<T> C;

	struct Polar{
		Polar(const T& mg, const T& ph): m(mg), p(ph){}
		T m, p;
	};

	union{
		struct{ T r, i; };
		//struct{ T x, y; };
		T elems[2];
	};
	
	Complex(const Complex& v): r(v.r), i(v.i){}
	Complex(const Polar& v){ *this = v; }

	Complex(const T& re=(T)1, const T& im=(T)0): r(re), i(im){}
	Complex(const T& m, const T& p, int fromPol){ (*this) = Polar(m,p); }

	//static Complex polar(const T& m, const T& p){ return C(Polar(m,p)); }

	C& fromPhase(const T& p){ r=cos(p); i=sin(p); return *this; }
	C& fromPolar(const T& m, const T& p){ return (*this)(Polar(m,p)); }

	C& operator()(const T& vr, const T& vi){ r=vr; i=vi; return *this; }
	C& operator()(const Polar& p){ return *this = p; }
	T& operator[](uint32_t ind){ return elems[ind];}
	const T& operator[](uint32_t i) const { return elems[i]; }
	
	C& operator = (const Polar& v){ r=v.m*cos(v.p); i=v.m*sin(v.p); return *this; }
	C& operator = (const C& v){ r=v.r; i=v.i; return *this; }
	C& operator = (const T& v){ r=v;   i=v;   return *this; }
	C  operator - () const { return C(-r, -i); }
	C  operator - (const C& v) const { return C(r-v.r, i-v.i); }
	C  operator - (const T& v) const { return C(r-v,   i-v); }
	C& operator -=(const C& v){ r-=v.r; i-=v.i; return *this; }
	C& operator -=(const T& v){ r-=v;   i-=v;   return *this; }
	C  operator + (const C& v) const { return C(r+v.r, i+v.i); }
	C  operator + (const T& v) const { return C(r+v,   i+v); }
	C& operator +=(const C& v){ r+=v.r; i+=v.i; return *this; }
	C& operator +=(const T& v){ r+=v;   i+=v;   return *this; }
	C  operator * (const C& v) const { C c(*this); return mul(c, v); }
	C  operator * (const T& v) const { return C(r*v,   i*v); }
	C& operator *=(const C& v){ return mul(*this, v); }
	C& operator *=(const T& v){ r*=v; i*=v; return *this; }
	C  operator / (const C& v) const { C c(*this); return div(c, v); }
	C  operator / (const T& v) const { return C(r/v, i/v); }
	C& operator /=(const C& v){ return div(*this, v); }
	C& operator /=(const T& v){ r/=v; i/=v; return *this; }

	C  conj() const { return C(r,-i); }

	T dot() const { return dot(*this); }
	T dot(const C& v) const { return r*v.r + i*v.i; }
	T mag() const { return sqrt(dot()); }
	C mul2(const C& v) const { return C(r*v.r, i*v.i); }
	C& normalize(){ return *this /= mag(); }
	T phase() const { return atan2(i, r); }
	C recip() const { T m=1./dot(); return C(r*m, -i*m); }
	
	bool operator < (const C& c) const { return dot() < c.dot(); }

	static C& mul(C& a, const C& b){
		return a(a.r*b.r - a.i*b.i, a.i*b.r + a.r*b.i);
	}

	static C& div(C& a, const C& b){
		T den = (T)1 / b.dot();
		return a((a.r*b.r + a.i*b.i) * den, (a.i*b.r - a.r*b.i) * den);
	}
};

typedef Complex<float > Complexf;
typedef Complex<double> Complexd;





/// Quaternion
template <class T=gam::real>
struct Quat{
	typedef Quat<T> Q;

	struct Unit3{
		Unit3(const T& vx, const T& vy, const T& vz, const T& scale=1)
		:	x(vx), y(vy), z(vz)
		{	T m = scale / sqrt(x*x + y*y + z*z); x*=m; y*=m; z*=m; }
		
		T x,y,z;
	};
	
	union{
		struct{ T r, i, j, k; };
		//struct{ T w, x, y, z; };
		T elems[4];
	};
	
	//Quat(const Quat& q): r(q.r), i(q.i), j(q.j), k(q.k){}
	Quat(const T& r=(T)1, const T& i=(T)0, const T& j=(T)0, const T& k=(T)0): r(r), i(i), j(j), k(k){}
	Quat(const T& a, const Unit3& u){ fromAxis(a,u); }

	Q& operator ()(const T& vr, const T& vi, const T& vj, const T& vk){ r=vr; i=vi; j=vj; k=vk; return *this; }
	T& operator[](uint32_t i){ return elems[i];}
	const T& operator[](uint32_t i) const { return elems[i]; }
	
//	Q& operator = (const Q& v){ r=v.r; i=v.i; j=v.j; k=v.k; return *this; }
//	Q& operator = (const T& v){ r=v;   i=v;   j=v;   k=v;   return *this; }
//	Q  operator - (const Q& v) const { return Q(r-v.r, i-v.i, j-v.j, k-v.k); }
	Q  operator - () const { return Q(-r, -i, -j, -k); }
	Q  operator - (const Q& v) const { return Q(r-v.r, i-v.i, j-v.j, k-v.k); }
	Q& operator -=(const Q& v){ r-=v.r; i-=v.i; j-=v.j; k-=v.k; return *this; }
	Q  operator + (const Q& v) const { return Q(r+v.r, i+v.i, j+v.j, k+v.k); }
	Q& operator +=(const Q& v){ r+=v.r; i+=v.i; j+=v.j; k+=v.k; return *this; }
	Q  operator * (const Q& v) const { Q q(*this); return mul(q, v); }
	Q  operator * (const T& v) const { return Q(r*v, i*v, j*v, k*v); }
	Q& operator *=(const Q& v){ return mul(*this, v); }
	Q& operator *=(const T& v){ r*=v; i*=v; j*=v; k*=v; return *this; }
//	Q  operator / (const Q& v) const { Q c(*this); return div(c, v); }
	Q  operator / (const T& v) const { return Q(r/v, i/v, j/v, k/v); }
//	Q& operator /=(const Q& v){ return div(*this, v); }
	Q& operator /=(const T& v){ r/=v; i/=v; j/=v; k/=v; return *this; }

	Q conj() const { return Q(r,-i,-j,-k); }
	
	T dot() const { return r*r + i*i + j*j + k*k; }
	
	// Set from angle (radians) and unit vector (x,y,z)
	Q& fromAxis(T a, T x, T y, T z){
		a *= (T)0.5;
		T s2 = sin(a);
		return (*this)(cos(a), x*s2, y*s2, z*s2);
	}
	
	Q& fromAxis(T a, const Unit3& u){ return fromAxis(a, u.x, u.y, u.z); }
	
	Q& fromEuler(T a, T b, T c){
		a*=(T)0.5; b*=(T)0.5; c*=(T)0.5;
		T c1=cos(a), s1=sin(a), c2=cos(b), s2=sin(b), c3=cos(c), s3=sin(c);
		T c1c2 = c1*c2;
		T s1s2 = s1*s2;
		return (*this)(c1c2*c3 - s1s2*s3, c1c2*s3 + s1s2*c3, s1*c2*c3 + c1*s2*s3, c1*s2*c3 - s1*c2*s3);
	}
	
	Q& identity(){ return (*this)(1,0,0,0); }
	
	T mag() const { return sqrt(dot()); }
	Q& normalize(){ return *this /= mag(); }
	
	/// Rotate a vector by current quaternion
	void rotate(T& x, T& y, T& z) const {
		Q p(-i*x - j*y - k*z, r*x + j*z - k*y, r*y - i*z + k*x,	r*z + i*y - j*x);
		p *= conj(); x=p.i; y=p.j; z=p.k;
	}
	
	void toAxis(T& a, T& x, T& y, T& z) const {
		a = (T)2 * acos(r);
		T s = 1./sqrt(i*i + j*j + k*k);
		x = i*s; y = j*s; z = k*s;
	}
	
	void toVectorX(T& x, T& y, T& z) const {
		x = (j*j + k*k) * (T)-2 + (T)1;
		y = (i*j + k*r) * (T) 2;
		z = (i*k - j*r) * (T) 2;	
	}

	void toVectorY(T& x, T& y, T& z) const {
		x = (i*j - k*r) * (T) 2;
		y = (i*i + k*k) * (T)-2 + (T)1;
		z = (j*k + i*r) * (T) 2;
	}

	void toVectorZ(T& x, T& y, T& z) const {
		x = (i*k + j*r) * (T) 2;
		y = (j*k - i*r) * (T) 2;
		z = (i*i + j*j) * (T)-2 + (T)1;
	}

	static Q& mul(Q& a, const Q& b){
		return a(
			a.r*b.r - a.i*b.i - a.j*b.j - a.k*b.k,
			a.r*b.i + a.i*b.r + a.j*b.k - a.k*b.j,
			a.r*b.j - a.i*b.k + a.j*b.r + a.k*b.i,
			a.r*b.k + a.i*b.j - a.j*b.i + a.k*b.r
		);
	}
	
//	static Q& mul3(Q& a, const Q& b){
//		return a(
//			a.r*b.r - a.i*b.i - a.j*b.j,
//			a.r*b.i + a.i*b.r,
//			a.r*b.j           + a.j*b.r,
//			          a.i*b.j - a.j*b.i
//		);
//	}
};





/// Multi-element container

/// This is a fixed size array to enable better loop unrolling optimizations
/// by the compiler and to avoid an extra 'size' data member for small-sized
/// arrays. It also lacks a constructor to allow C-style struct initializations.
template <uint32_t N, class T>
struct Multi{

//	Multi(){}
//	Multi(const T& e ){ mem::set(elems, N, e); }
//	Multi(const T* es){ mem::copy(elems, es, N); }

	T elems[N];
	
	/// Set element at index with no bounds checking
	T& operator[](uint32_t i){ return elems[i];}
	
	/// Get element at index with no bounds checking
	const T& operator[](uint32_t i) const { return elems[i]; }

	/// Returns size of array
	static uint32_t size(){ return N; }

	/// Zeros all elements.
	void zero(){ memset(elems, 0, N * sizeof(T)); }
};


template <class T>
struct Multi3: public Multi<3,T>{
	Multi3(T v1=0, T v2=0, T v3=0){ (*this)[0]=v1; (*this)[1]=v2; (*this)[2]=v3; }
};




// 3-tuple
template <class T1, class T2, class T3>
struct Tup3{
	Tup3(T1 v1=0, T2 v2=0, T3 v3=0): v1(v1), v2(v2), v3(v3){}
	
	void set(T1 a1, T2 a2, T3 a3){ v1=a1; v2=a2; v3=a3; }
	
	T1 v1; T2 v2; T3 v3;
};

typedef Tup3<float, float, float> Tup3f;

template <class T1, class T2, class T3> 
inline Tup3<T1,T2,T3> tup(T1 v1, T2 v2, T3 v3){ return Tup3<T1,T2,T3>(v1,v2,v3); }




///< Fixed size vector.
template <uint32_t N, class T>
struct Vec : public Multi<N,T> {

	typedef Vec<N,T> V;

	Vec(const T& v=T(0)){ (*this) = v; }

	#define DO for(uint32_t i=0; i<N; ++i)

	bool operator !=(const V& v){ DO{ if((*this)[i] == v[i]) return false; } return true; }
	bool operator !=(const T& v){ DO{ if((*this)[i] == v   ) return false; } return true; }
	V& operator = (const V& v){ DO (*this)[i] = v[i]; return *this; }
	V& operator = (const T& v){ DO (*this)[i] = v;    return *this; }
	bool operator ==(const V& v){ DO{ if((*this)[i] != v[i]) return false; } return true; }
	bool operator ==(const T& v){ DO{ if((*this)[i] != v   ) return false; } return true; }
	V  operator * (const V& v) const { V r; DO r[i] = (*this)[i] * v[i]; return r; }
	V  operator * (const T& v) const { V r; DO r[i] = (*this)[i] * v;    return r; }
	V& operator *=(const V& v){ DO (*this)[i] *= v[i]; return *this; }
	V& operator *=(const T& v){ DO (*this)[i] *= v;    return *this; }
	V  operator / (const V& v) const { V r; DO r[i] = (*this)[i] / v[i]; return r; }
	V  operator / (const T& v) const { V r; DO r[i] = (*this)[i] / v;    return r; }
	V& operator /=(const V& v){ DO (*this)[i] /= v[i]; return *this; }
	V& operator /=(const T& v){ DO (*this)[i] /= v;    return *this; }
	V  operator - (          ) const { V r; DO r[i] = -(*this)[i]; return r; }
	V  operator - (const V& v) const { V r; DO r[i] = (*this)[i] - v[i]; return r; }
	V  operator - (const T& v) const { V r; DO r[i] = (*this)[i] - v;    return r; }
	V& operator -=(const V& v){ DO (*this)[i] -= v[i]; return *this; }
	V& operator -=(const T& v){ DO (*this)[i] -= v;    return *this; }
	V  operator + (const V& v) const { V r; DO r[i] = (*this)[i] + v[i]; return r; }
	V  operator + (const T& v) const { V r; DO r[i] = (*this)[i] + v;    return r; }
	V& operator +=(const V& v){ DO (*this)[i] += v[i]; return *this; }
	V& operator +=(const T& v){ DO (*this)[i] += v;    return *this; }

	T dot() const { return dot(*this); }
	T dot(const V& v) const { T r=(T)0; DO r+=(*this)[i]*v[i]; return r; }

	#undef DO
};



///< Two element vector
template <class T>
struct Vec2 : public Vec<2, T> {
	using Vec<2,T>::operator=;
	Vec2(const Vec<2, T>& v){ (*this) = v; }
	Vec2(T v=(T)0){ (*this)(v, v); }
	Vec2(T v1, T v2){ (*this)(v1, v2); }
	void operator()(T v1, T v2){ (*this)[0]=v1; (*this)[1]=v2; }
};


///< Three element vector
template <class T>
struct Vec3 : public Vec<3, T> {
	using Vec<3,T>::operator=;
	Vec3(const Vec<3, T>& v){ (*this) = v; }
	Vec3(T v=(T)0){ (*this)(v, v, v); }
	Vec3(T v1, T v2, T v3=(T)0){ (*this)(v1, v2, v3); }
	void operator()(T v1, T v2, T v3){ (*this)[0]=v1; (*this)[1]=v2; (*this)[2]=v3; }
	
	T dot() const { return dot(*this); }
	T dot(const Vec3& v) const { return v[0]*(*this)[0] + v[1]*(*this)[1] + v[2]*(*this)[2]; }
};






// We want this thing to be able to provide a standard interface to indexing
// classes without relying heavily (or at all) on virtual functions.

class Indexer{
public:
	Indexer(uint end=1, uint stride=1, uint begin=0)
	:	mEnd(end), mBegin(begin), mStride(stride)
	{}
	
	virtual ~Indexer(){}
	
	uint begin() const { return mBegin; }			// Begin index (inclusive)
	bool cond(uint i) const { return i < end(); }	// Continuation condition
	uint end() const { return mEnd; }				// End index (exclusive)
	virtual uint index(uint i) const { return i; }	// Iteration to index map 
	uint stride() const { return mStride; }			// Index stride amount

protected:
	uint32_t mEnd, mBegin, mStride;
};

typedef Indexer Loop;


class Indices : public Indexer{
public:

	Indices(uint maxSize)
	:	Indexer(0), mMaxSize(maxSize), mElems(new uint[maxSize])
	{}

	~Indices(){ delete[] mElems; }

	// add value to indices
	void operator+= (uint v){
		for(uint i=0; i<end(); ++i)
			mElems[i] = (mElems[i] + v) % mMaxSize;
	}

	// multiply indices by value
	void operator*= (double v){
		for(uint i=0; i<end(); ++i)
			mElems[i] = ((uint)((double)mElems[i] * v + 0.5)) % mMaxSize;
	}

	// add new index
	Indices& operator<< (uint index){
		if(end() < mMaxSize && index < mMaxSize) mElems[mEnd++] = index;
		return *this;
	}

	void clear(){ mEnd=0; }
	uint * elems() const { return mElems; }

	uint index(uint i) const { return mElems[i]; }

private:
	uint mMaxSize;
	uint * mElems;
};




/*
// Upward counting indexer
class Loop{
public:
	Loop(uint end=1, uint stride=1, uint begin=0)
	:	mEnd(end), mBegin(begin), mStride(stride)
	{}
	
	Loop& operator()(uint end, uint stride=1, uint begin=0){
		mEnd = end; mBegin = begin; mStride = stride; return *this;
	}
	
	uint begin() const { return mBegin; }			// Begin index (inclusive)
	bool cond(uint i) const { return i < end(); }	// Continuation condition
	uint end() const { return mEnd; }				// End index (exclusive)
	uint index(uint i) const { return i; }			// Iteration to index map 
	uint stride() const { return mStride; }			// Index stride amount

private:
	uint mEnd, mBegin, mStride;
};
*/


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


}

#endif
