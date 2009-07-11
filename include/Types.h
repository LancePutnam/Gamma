#ifndef GAMMA_TYPES_H_INC
#define GAMMA_TYPES_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information 

	File Description: 
	Static (fixed size) types
*/


#include "pstdint.h"			// for cross-platform uint32_t, uint16_t, etc...
#include <math.h>
#include <stdlib.h>

namespace gam{

typedef unsigned long	uint;	// default natural number type
typedef float			real;	// default real number type


/// Integer-based bit field
template <class T=uint32_t>
class Bits{
public:

	Bits(): mVal(T(0)){}
	
	/// Returns current value
	T operator()() const { return mVal; }
	
	/// Returns mask of enabled bits in input mask
	T enabled(T bitMask) const { return mVal & bitMask; }

	T triggered(T bitMask){
		T r = enabled(bitMask);
		disable(bitMask);
		return r;
	}
	
	/// Disable bits using mask
	Bits& disable(T bitMask){ mVal&=~bitMask; return *this; }
	
	/// Enable bits using mask
	Bits& enable(T bitMask){ mVal|=bitMask; return *this; }
	
	/// Set bits using mask
	Bits& set(T bitMask, bool v){ v ? enable(bitMask) : disable(bitMask); return *this; }
	
	/// Toggle bits using mask
	Bits& toggle(T bitMask){ mVal^=bitMask; return *this; }
	
	/// Zero all bits
	Bits& zero(){ mVal=T(0); return *this; }

	/// Get mask from specific bit number where 0 is the LSB.
	static T mask(T bitNum){ return 1<<bitNum; }
	
	/// Get mask from two bit numbers where 0 is the LSB.
	static T mask(T bitNum1, T bitNum2){ return (1<<bitNum1) | (1<<bitNum2); }

private:
	T mVal;
};



/// Complex number
template <class T=gam::real>
struct Complex{

	typedef Complex<T> C;

	struct Polar{
		Polar(const T& m, const T& p): m(m), p(p){}
		T m, p;
	};

	union{
		struct{ T r, i; };
		//struct{ T x, y; };
		T elems[2];
	};
	
	Complex(const Complex& v): r(v.r), i(v.i){}
	Complex(const Polar& v){ *this = v; }

	Complex(const T& r=(T)1, const T& i=(T)0): r(r), i(i){}
	Complex(const T& m, const T& p, int fromPolar){ (*this) = Polar(m,p); }

	C& fromPhase(const T& p){ r=::cos(p); i=::sin(p); return *this; }
	C& fromPolar(const T& m, const T& p){ return (*this)(Polar(m,p)); }

	C& operator()(const T& vr, const T& vi){ r=vr; i=vi; return *this; }
	C& operator()(const Polar& p){ return *this = p; }
	T& operator[](uint32_t i){ return elems[i];}
	const T& operator[](uint32_t i) const { return elems[i]; }

	bool operator ==(const C& v) const { return (r==v.r) && (i==v.i); }		///< Returns true if all components are equal
	bool operator !=(const C& v) const { return (r!=v.r) || (i!=v.i); }		///< Returns true if any components are not equal
	bool operator > (const C& v) const { return norm2() > v.norm2(); }		///< Returns true if norm is greater than argument's norm
	bool operator < (const C& c) const { return norm2() < c.norm2(); }		///< Returns true if norm is less than argument's norm

	C& operator = (const Polar& v){ r=v.m*::cos(v.p); i=v.m*::sin(v.p); return *this; }
	C& operator = (const C& v){ r=v.r; i=v.i; return *this; }
	C& operator = (const T& v){ r=v;   i=T(0); return *this; }
	C& operator -=(const C& v){ r-=v.r; i-=v.i; return *this; }
	C& operator -=(const T& v){ r-=v; return *this; }
	C& operator +=(const C& v){ r+=v.r; i+=v.i; return *this; }
	C& operator +=(const T& v){ r+=v; return *this; }
	C& operator *=(const C& v){ return (*this)(r*v.r - i*v.i, i*v.r + r*v.i); }
	C& operator *=(const T& v){ r*=v; i*=v; return *this; }
	C& operator /=(const C& v){ return (*this) *= v.recip(); }
	C& operator /=(const T& v){ r/=v; i/=v; return *this; }

	const C operator - () const { return C(-r, -i); }
	const C operator - (const C& v) const { return C(*this) -= v; }
	const C operator - (const T& v) const { return C(*this) -= v; }
	const C operator + (const C& v) const { return C(*this) += v; }
	const C operator + (const T& v) const { return C(*this) += v; }
	const C operator * (const C& v) const { return C(*this) *= v; }
	const C operator * (const T& v) const { return C(*this) *= v; }
	const C operator / (const C& v) const { return C(*this) /= v; }
	const C operator / (const T& v) const { return C(*this) /= v; }
	
	T arg() const { return atan2(i, r); }					///< Returns argument (angle)
	C conj() const { return C(r,-i); }						///< Returns conjugate, z*
	T dot(const C& v) const { return r*v.r + i*v.i; }		///< Returns vector dot product
	const C exp() const { return Polar(::exp(r), i); }		///< Returns e^z
	const C log() const { return Complex<T>(T(0.5)*::log(norm2()), arg()); } ///< Returns log(z)
	//C mul2(const C& v) const { return C(r*v.r, i*v.i); }
	T norm() const { return sqrt(norm2()); }				///< Returns norm (radius), |z|
	T norm2() const { return dot(*this); }					///< Returns square of norm, |z|^2
	C& normalize(){ return *this /= norm(); }				///< Sets norm (radius) to 1, |z|=1
	const C pow(const C& v) const { return ((*this).log()*v).exp(); }	///< Returns z^v
	const C pow(const T& v) const { return ((*this).log()*v).exp(); }	///< Returns z^v
	const C recip() const { return conj()/norm2(); }		///< Return multiplicative inverse, 1/z
	const C sgn() const { return C(*this).normalize(); }	///< Returns signum, z/|z|, the closest point on unit circle
	const C sqr() const { return C(r*r-i*i, T(2)*r*i); }	///< Returns square

	const C cos() const { return C(::cos(r)*::cosh(i),-::sin(r)*::sinh(i)); } ///< Returns cos(z)
	const C sin() const { return C(::sin(r)*::cosh(i), ::cos(r)*::sinh(i)); } ///< Returns sin(z)

	T abs() const { return norm(); }						///< Returns norm (radius), |z|
	T mag() const { return norm(); }						///< Returns norm (radius), |z|
	T phase() const { return arg(); }						///< Returns argument (angle)
};

typedef Complex<float > Complexf;
typedef Complex<double> Complexd;

#define TEM template <class T>
TEM const Complex<T> exp(const Complex<T>& c){ return c.exp(); }
TEM const Complex<T> log(const Complex<T>& c){ return c.log(); }
TEM const Complex<T> pow(const Complex<T>& b, const Complex<T>& e){ return b.pow(e); }
TEM const Complex<T> operator + (T r, const Complex<T>& c){ return  c+r; }
TEM const Complex<T> operator - (T r, const Complex<T>& c){ return -c+r; }
TEM const Complex<T> operator * (T r, const Complex<T>& c){ return  c*r; }
TEM const Complex<T> operator / (T r, const Complex<T>& c){ return  c.conj()*(r/c.norm()); }
#undef TEM



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
	Quat(const T& r=T(1), const T& i=T(0), const T& j=T(0), const T& k=T(0)): r(r), i(i), j(j), k(k){}
	Quat(const T& a, const Unit3& u){ fromAxis(a,u); }

	Q& operator()(const T& vr, const T& vi, const T& vj, const T& vk){ r=vr; i=vi; j=vj; k=vk; return *this; }
	T& operator[](uint32_t i){ return elems[i];}
	const T& operator[](uint32_t i) const { return elems[i]; }

	bool operator ==(const Q& v) const { return (r==v.r) && (i==v.i) && (j==v.j) && (k==v.k); } ///< Returns true if all components are equal
	bool operator !=(const Q& v) const { return (r!=v.r) || (i!=v.i) || (j!=v.j) || (k!=v.k); } ///< Returns true if any components are not equal
	
	Q& operator = (const Q& v){ r=v.r; i=v.i; j=v.j; k=v.k; return *this; }
	Q& operator = (const T& v){ r=v;   i=T(0);j=T(0);k=T(0);return *this; }
	Q& operator -=(const Q& v){ r-=v.r; i-=v.i; j-=v.j; k-=v.k; return *this; }
	Q& operator -=(const T& v){ r-=v; return *this; }
	Q& operator +=(const Q& v){ r+=v.r; i+=v.i; j+=v.j; k+=v.k; return *this; }
	Q& operator +=(const T& v){ r+=v; return *this; }
	Q& operator *=(const Q& v){ return (*this)(
		r*v.r - i*v.i - j*v.j - k*v.k,
		r*v.i + i*v.r + j*v.k - k*v.j,
		r*v.j - i*v.k + j*v.r + k*v.i,
		r*v.k + i*v.j - j*v.i + k*v.r);
	}
	Q& operator *=(const T& v){ r*=v; i*=v; j*=v; k*=v; return *this; }
	Q& operator /=(const Q& v){ return (*this) *= v.recip(); }
	Q& operator /=(const T& v){ r/=v; i/=v; j/=v; k/=v; return *this; }

	Q operator - () const { return Q(-r, -i, -j, -k); }
	Q operator - (const Q& v) const {return Q(*this)-=v;}
	Q operator - (const T& v) const {return Q(*this)-=v;}
	Q operator + (const Q& v) const {return Q(*this)+=v;}
	Q operator + (const T& v) const {return Q(*this)+=v;}
	Q operator * (const Q& v) const {return Q(*this)*=v;}
	Q operator * (const T& v) const {return Q(*this)*=v;}
	Q operator / (const Q& v) const {return Q(*this)/=v;}
	Q operator / (const T& v) const {return Q(*this)/=v;}

	Q conj() const { return Q(r,-i,-j,-k); }				///< Returns conjugate
	T dot(const Q& v) const { return r*v.r + i*v.i + j*v.j + k*v.k; }		///< Returns vector dot product
	Q& identity(){ return (*this)(1,0,0,0); }
	T norm() const { return sqrt(norm2()); }				///< Returns norm (length)
	T norm2() const { return dot(*this); }					///< Returns square of norm
	Q& normalize(){ return *this /= norm(); }				///< Sets norm to 1
	const Q recip() const { return conj()/norm2(); }		///< Return multiplicative inverse
	const Q sgn() const { return Q(*this).normalize(); }	///< Returns signum, q/|q|, the closest point on unit 3-sphere
	
	// Set from angle (radians) and unit vector (x,y,z)
	Q& fromAxis(T a, T x, T y, T z){
		a *= T(0.5);
		return fromAxis(cos(a), sin(a), x,y,z);
	}

	// Set from cosine and sine of half-angle (radians) and unit vector (x,y,z)
	Q& fromAxis(T cosA2, T sinA2, T x, T y, T z){
		return (*this)(cosA2, x*sinA2, y*sinA2, z*sinA2);
	}
	
	Q& fromAxis(T a, const Unit3& u){ return fromAxis(a, u.x, u.y, u.z); }
	
	Q& fromEuler(T a, T b, T c){
		a*=(T)0.5; b*=(T)0.5; c*=(T)0.5;
		T c1=cos(a), s1=sin(a), c2=cos(b), s2=sin(b), c3=cos(c), s3=sin(c);
		T c1c2 = c1*c2;
		T s1s2 = s1*s2;
		return (*this)(c1c2*c3 - s1s2*s3, c1c2*s3 + s1s2*c3, s1*c2*c3 + c1*s2*s3, c1*s2*c3 - s1*c2*s3);
	}
	

	void rotate(T& x, T& y, T& z) const {
		Q p(-i*x - j*y - k*z, r*x + j*z - k*y, r*y - i*z + k*x,	r*z + i*y - j*x);
		p *= conj(); x=p.i; y=p.j; z=p.k;
	}
	
	/// Rotate a vector by current quaternion
	template <class V3>
	void rotate(V3& v) const {
		Q p(-i*v[0] - j*v[1] - k*v[2], r*v[0] + j*v[2] - k*v[1], r*v[1] - i*v[2] + k*v[0], r*v[2] + i*v[1] - j*v[0]);
		p *= conj(); v[0]=p.i; v[1]=p.j; v[2]=p.k;
	}

	// Rotates a vector assuming input z component is zero
	template <class V3>
	void rotateXY(V3& v) const {
		Q p(-i*v[0] - j*v[1], r*v[0] - k*v[1], r*v[1] + k*v[0], i*v[1] - j*v[0]);
		p *= conj(); v[0]=p.i; v[1]=p.j; v[2]=p.k;
	}

	// Rotates a vector assuming input y and z components are zero
	template <class V3>
	void rotateX(V3& v) const {
		Q p(-i*v[0], r*v[0], k*v[0], -j*v[0]);
		p *= conj(); v[0]=p.i; v[1]=p.j; v[2]=p.k;
	}
	
	// Rotates a vector assuming input x and y components are zero
	template <class V3>
	void rotateZ(V3& v) const {
		Q p(-k*v[2], j*v[2], -i*v[2], r*v[2]);
		p *= conj(); v[0]=p.i; v[1]=p.j; v[2]=p.k;
	}
	
	/// Perform spherical linear interpolation from this to target, q.
	Q slerp(const Q& q, T f){
		T dot = q.dot(*this);
		if(dot > 0.9995) return (*this + f*(q - *this)).normalize(); // linear interpolation
		T theta = acos(dot < T(-1) ? T(-1) : (dot > T(1) ? T(1) : dot)*f);
		return (*this)*::cos(theta) + (q - (*this)*dot).normalize()*::sin(theta);
	}
	
	void toAxis(T& a, T& x, T& y, T& z) const {
		a = T(2) * acos(r);
		T s = T(1)/sqrt(i*i + j*j + k*k);
		x = i*s; y = j*s; z = k*s;
	}

	void toVectorX(T& x, T& y, T& z) const {
		x = (j*j + k*k) * T(-2) + T(1);
		y = (i*j + k*r) * T( 2);
		z = (i*k - j*r) * T( 2);	
	}

	void toVectorY(T& x, T& y, T& z) const {
		x = (i*j - k*r) * T( 2);
		y = (i*i + k*k) * T(-2) + T(1);
		z = (j*k + i*r) * T( 2);
	}

	void toVectorZ(T& x, T& y, T& z) const {
		x = (i*k + j*r) * T( 2);
		y = (j*k - i*r) * T( 2);
		z = (i*i + j*j) * T(-2) + T(1);
	}

	/// Set quaternion from rotation matrix
	
	/// Code taken from van Waveren (2005). 'From Quaternion to Matrix and Back'
	///
	template <class V3>
	void fromMat3(const V3& vx, const V3& vy, const V3& vz){
		
		const T &x0=vx[0], &x1=vx[1], &x2=vx[2];
		const T &y0=vy[0], &y1=vy[1], &y2=vy[2];
		const T &z0=vz[0], &z1=vz[1], &z2=vz[2];
		
		if(x0 + y1 + z2 > T(0)){
			T t = x0 + y1 + z2 + T(1);
			T s = T(0.5)/sqrt(t);
			r = t * s;
			i = (y2 - z1) * s;
			j = (z0 - x2) * s;
			k = (x1 - y0) * s;
		}
		else if(x0 > y1 && x0 > z2){
			T t = x0 - y1 - z2 + T(1);
			T s = T(0.5)/sqrt(t);
			r = (y2 - z1) * s;
			i = t * s;
			j = (x1 + y0) * s;
			k = (x2 + z0) * s;
		} 
		else if(y1 > z2){
			T t =-x0 + y1 - z2 + T(1);
			T s = T(0.5)/sqrt(t);
			r = (z0 - x2) * s;
			i = (x1 + y0) * s;
			j = t * s;
			k = (y2 + z1) * s;
		} 
		else{
			T t =-x0 - y1 + z2 + T(1);
			T s = T(0.5)/sqrt(t);
			r = (x1 - y0) * s;
			i = (x2 + z0) * s;
			j = (y2 + z1) * s;
			k = t * s;
		}
	}
	
	/*
		Quaternion mult can be broken down into 4 complex mult:
		
		A = Ari + Ajk
		B = Bri + Bjk
		
		Ari Bri		-Ajk   Bjk*
		Ari Bjk		 Ajk* -Bri
		
		'-' rotates 180 and '*' (conj) reflects around first axis
	*/
	
//	static Q& mul3(Q& a, const Q& b){
//		return a(
//			a.r*b.r - a.i*b.i - a.j*b.j,
//			a.r*b.i + a.i*b.r,
//			a.r*b.j           + a.j*b.r,
//			          a.i*b.j - a.j*b.i
//		);
//	}
};

#define TEM template <class T>
TEM const Quat<T> operator + (T r, const Quat<T>& q){ return  q+r; }
TEM const Quat<T> operator - (T r, const Quat<T>& q){ return -q+r; }
TEM const Quat<T> operator * (T r, const Quat<T>& q){ return  q*r; }
TEM const Quat<T> operator / (T r, const Quat<T>& q){ return  q.conj()*(r/q.norm2()); }
#undef TEM

typedef Quat<float> Quatf;
typedef Quat<double> Quatd;



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
	T norm() const { return sqrt(dot()); }

	V sgn() const { V(*this) /= norm(); }

	#undef DO
};



///< Two element vector
template <class T>
struct Vec2 : public Vec<2, T> {
	using Vec<2,T>::operator=;
	Vec2(const Vec<2, T>& v){ (*this)=v; }
	Vec2(const T& v=T(0)){ (*this)(v,v); }
	Vec2(const T& v1, const T& v2){ (*this)(v1,v2); }
	Vec2& operator()(const T& v1, const T& v2){ (*this)[0]=v1; (*this)[1]=v2; return *this; }
};


///< Three element vector
template <class T>
struct Vec3 : public Vec<3, T> {
	using Vec<3,T>::operator=;

	Vec3(const Vec<3, T>& v){ (*this)=v; }
	Vec3(const T& v=T(0)){ (*this)(v,v,v); }
	Vec3(const T& v1, const T& v2, const T& v3=T(0)){ (*this)(v1,v2,v3); }

	Vec3& operator()(const T& v1, const T& v2, const T& v3){ (*this)[0]=v1; (*this)[1]=v2; (*this)[2]=v3; return *this; }
	
	T dot() const { return dot(*this); }
	T dot(const Vec3& v) const { return v[0]*(*this)[0] + v[1]*(*this)[1] + v[2]*(*this)[2]; }
	Vec3 cross(const Vec3& v) const {
		Vec3 r; const Vec3& t = *this;
		r[0] = t[1]*v[2] - t[2]*v[1];
		r[1] = t[2]*v[0] - t[0]*v[2];
		r[2] = t[0]*v[1] - t[1]*v[0];
		return r;
	}
};


///< Four element vector
template <class T>
struct Vec4 : public Vec<4, T> {
	using Vec<4,T>::operator=;

	Vec4(const Vec<4, T>& v){ (*this)=v; }
	Vec4(const T& v=T(0)){ (*this)(v,v,v,v); }
	Vec4(const T& v1, const T& v2, const T& v3, const T& v4){ (*this)(v1,v2,v3,v4); }

	Vec4& operator()(const T& v1, const T& v2, const T& v3, const T& v4){
		(*this)[0]=v1; (*this)[1]=v2; (*this)[2]=v3; (*this)[3]=v4; return *this; }
	
	T dot() const { return dot(*this); }
	T dot(const Vec4& v) const { return v[0]*(*this)[0] + v[1]*(*this)[1] + v[2]*(*this)[2] + v[3]*(*this)[3]; }
	
	Vec3<T> xyz() const { return Vec3<T>((*this)[0], (*this)[1], (*this)[2]); }
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
	:	Indexer(0), mMaxSize(maxSize)
	{	mElems = new uint[maxSize]; }

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




// Trying to abstract too many things. Hard to write subclasses of templates
//template <class T>
//struct ElemData2{
//	union{
//		struct{ T r,i; };
//		struct{ T x,y; };
//		T elems[2];
//	};
//};
//
//
//template <uint32_t N, class T, template <class> class E>
//struct ElemN : public E<T>{
//	using E<T>::elems;
//
//	/// Set element at index (no bounds checking)
//	T& operator[](uint32_t i){ return elems[i]; }
//	
//	/// Get element at index (no bounds checking)
//	const T& operator[](uint32_t i) const { return elems[i]; }
//
//	/// Returns size of array
//	static uint32_t size(){ return N; }
//};
//
//template <class T>
//struct Elem2 : public ElemN<2, T, ElemData2>{
//	using ElemN<2, T, ElemData2>::elems;
//	using ElemData2<T>::x;
//	//Elem2(const T& v1, const T& v2): x(v1), y(v2){}
//	Elem2(const T& v1, const T& v2){ elems[0]=v1; elems[1]=v2; }
//};



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
