#ifndef GAMMA_TYPES_H_INC
#define GAMMA_TYPES_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information 

	File Description: 
	Static (fixed size) POD types including complex numbers, quaternions, 
	2,3,4-vectors, 3x3 matrix, and fixed-size array.
	
*/


#include "Gamma/pstdint.h"		// for cross-platform uint32_t, uint16_t, etc...
#include <math.h>
#include <stdlib.h>

namespace gam{

typedef float			real;	// default real number type


template<class T> class Complex;
template<class T> class Quat;
template<class T> class Mat3;
template<class T> class Vec2;
template<class T> class Vec3;
template<class T> class Vec4;


//typedef Polar<float > Polarf;
//typedef Polar<double> Polard;
//typedef Complex<float > Complexf;
//typedef Complex<double> Complexd;
//typedef Mat3<float> Mat3f;
//typedef Mat3<double> Mat3d;
//typedef Vec2<float> Vec2f;
//typedef Vec2<double> Vec2d;
//typedef Vec3<float> Vec3f;
//typedef Vec3<double> Vec3d;
//typedef Vec4<float> Vec4f;
//typedef Vec4<double> Vec4d;
//typedef Quat<float> Quatf;
//typedef Quat<double> Quatd;


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



/// Polar number
template <class T>
struct Polar{
	Polar(const T& p): m(1.), p(p){}
	Polar(const T& m, const T& p): m(m), p(p){}
	Polar(const Complex<T>& v){ *this = v; }

	Polar& operator = (const Complex<T>& v){ m=v.norm(); p=v.arg(); return *this; }

	T m, p;
};


/// Complex number
template <class T=gam::real>
struct Complex{

	typedef Complex<T> C;

	union{
		struct{ T r, i; };
		T elems[2];
	};
	
	Complex(const Complex& v): r(v.r), i(v.i){}
	Complex(const Polar<T>& v){ *this = v; }
	Complex(const T& r=(T)1, const T& i=(T)0): r(r), i(i){}
	Complex(const T& m, const T& p, int fromPolar){ (*this) = Polar<T>(m,p); }

	
	C& arg(const T& v){ return fromPolar(norm(), v); }					///< Set phase leaving magnitude the same
	C& fromPhase(const T& v){ r=::cos(v); i=::sin(v); return *this; }	///< Set phase and normalize
	C& fromPolar(const T& m, const T& p){ return (*this)(Polar<T>(m,p)); }	///< Set magnitude and phase
	C& norm(const T& v){ return fromPolar(v, arg()); }					///< Set magnitude leaving phase the same

	C& operator()(const T& vr, const T& vi){ r=vr; i=vi; return *this; }
	C& operator()(const Polar<T>& p){ return *this = p; }
	T& operator[](uint32_t i){ return elems[i];}
	const T& operator[](uint32_t i) const { return elems[i]; }

	bool operator ==(const C& v) const { return (r==v.r) && (i==v.i); }		///< Returns true if all components are equal
	bool operator ==(const T& v) const { return (r==v  ) && (i==T(0));}		///< Returns true if real and equals value
	bool operator !=(const C& v) const { return (r!=v.r) || (i!=v.i); }		///< Returns true if any components are not equal
	bool operator > (const C& v) const { return norm2() > v.norm2(); }		///< Returns true if norm is greater than argument's norm
	bool operator < (const C& c) const { return norm2() < c.norm2(); }		///< Returns true if norm is less than argument's norm

	C& operator = (const Polar<T>& v){ r=v.m*::cos(v.p); i=v.m*::sin(v.p); return *this; }
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

	C operator - () const { return C(-r, -i); }
	C operator - (const C& v) const { return C(*this) -= v; }
	C operator - (const T& v) const { return C(*this) -= v; }
	C operator + (const C& v) const { return C(*this) += v; }
	C operator + (const T& v) const { return C(*this) += v; }
	C operator * (const C& v) const { return C(*this) *= v; }
	C operator * (const T& v) const { return C(*this) *= v; }
	C operator / (const C& v) const { return C(*this) /= v; }
	C operator / (const T& v) const { return C(*this) /= v; }
	
	T arg() const { return atan2(i, r); }					///< Returns argument (angle)
	C conj() const { return C(r,-i); }						///< Returns conjugate, z*
	T dot(const C& v) const { return r*v.r + i*v.i; }		///< Returns vector dot product
	C exp() const { return Polar<T>(::exp(r), i); }			///< Returns e^z
	C log() const { return Complex<T>(T(0.5)*::log(norm2()), arg()); } ///< Returns log(z)
	T norm() const { return ::sqrt(norm2()); }				///< Returns norm (radius), |z|
	T norm2() const { return dot(*this); }					///< Returns square of norm, |z|^2
	C& normalize(){ return *this /= norm(); }				///< Sets norm (radius) to 1, |z|=1
	C pow(const C& v) const { return ((*this).log()*v).exp(); }	///< Returns z^v
	C pow(const T& v) const { return ((*this).log()*v).exp(); }	///< Returns z^v
	C recip() const { return conj()/norm2(); }				///< Return multiplicative inverse, 1/z
	C sgn() const { return C(*this).normalize(); }			///< Returns signum, z/|z|, the closest point on unit circle
	C sqr() const { return C(r*r-i*i, T(2)*r*i); }			///< Returns square

	/// Returns square root
	C sqrt() const {
		static const T c = T(1)/::sqrt(T(2));
		T n = norm();
		T a = ::sqrt(n+r) * c;
		T b = ::sqrt(n-r) * (i<T(0) ? -c : c);		
		return C(a,b);
	}

	C cos() const { return C(::cos(r)*::cosh(i),-::sin(r)*::sinh(i)); } ///< Returns cos(z)
	C sin() const { return C(::sin(r)*::cosh(i), ::cos(r)*::sinh(i)); } ///< Returns sin(z)

	T abs() const { return norm(); }						///< Returns norm (radius), |z|
	T mag() const { return norm(); }						///< Returns norm (radius), |z|
	T phase() const { return arg(); }						///< Returns argument (angle)
};

#define TEM template <class T>
TEM Complex<T> exp(const Complex<T>& c){ return c.exp(); }
TEM Complex<T> log(const Complex<T>& c){ return c.log(); }
TEM Complex<T> pow(const Complex<T>& b, const Complex<T>& e){ return b.pow(e); }
TEM Complex<T> operator + (T r, const Complex<T>& c){ return  c+r; }
TEM Complex<T> operator - (T r, const Complex<T>& c){ return -c+r; }
TEM Complex<T> operator * (T r, const Complex<T>& c){ return  c*r; }
TEM Complex<T> operator / (T r, const Complex<T>& c){ return  c.conj()*(r/c.norm()); }
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
		a*=T(0.5); b*=T(0.5); c*=T(0.5);
		T c1=cos(a), s1=sin(a), c2=cos(b), s2=sin(b), c3=cos(c), s3=sin(c);
		T c1c2 = c1*c2;
		T s1s2 = s1*s2;
		T c1s2 = c1*s2;
		T s1c2 = s1*c2;
		return (*this)(c1c2*c3 - s1s2*s3, c1c2*s3 + s1s2*c3, s1c2*c3 + c1s2*s3, c1s2*c3 - s1c2*s3);
		
//		Q qa, qb, qc;
//		qa.fromAxis(c, 1,0,0);
//		qb.fromAxis(a, 0,1,0);
//		qc.fromAxis(b, 0,0,1);
//		return *this = qa*qb*qc;
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

	void rotate(T& x, T& y, T& z) const {
		Q p(-i*x - j*y - k*z, r*x + j*z - k*y, r*y - i*z + k*x,	r*z + i*y - j*x);
		p *= conj(); x=p.i; y=p.j; z=p.k;
	}
	
	/// Rotate a vector by current quaternion
	template <class V3>
	void rotate(V3& v) const {
		Q p(
			-i*v[0] - j*v[1] - k*v[2],
			 r*v[0] + j*v[2] - k*v[1],
			 r*v[1] - i*v[2] + k*v[0],
			 r*v[2] + i*v[1] - j*v[0]
		);
		p *= conj();
		v[0]=p.i; v[1]=p.j; v[2]=p.k;
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
		T s = sqrt(T(1) - r*r);
		if(s > T(0.00000001)){ x = i*s; y = j*s; z = k*s; }
		else{ x=T(1); y=z=T(0); }
	}
	
	/// Convert to 3x3 rotation matrix
	void toMatrix(Mat3<T>& m) const {
		static const T _2 = T(2);
		static const T _1 = T(1);		
		T _2r=_2*r, _2i=_2*i, _2j=_2*j;
		T ri=_2r*i, rj=_2r*j, rk=_2r*k, ii=_2i*i, ij=_2i*j, ik=_2i*k, jj=_2j*j, jk=_2j*k, kk=_2*k*k;
		m(	(-kk-jj)+_1,( ij-rk),	( rj+ik),
			( rk+ij),	(-kk-ii)+_1,( jk-ri),
			( ik-rj),	( ri+jk),	(-jj-ii)+_1);
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
TEM Quat<T> operator + (T r, const Quat<T>& q){ return  q+r; }
TEM Quat<T> operator - (T r, const Quat<T>& q){ return -q+r; }
TEM Quat<T> operator * (T r, const Quat<T>& q){ return  q*r; }
TEM Quat<T> operator / (T r, const Quat<T>& q){ return  q.conj()*(r/q.norm2()); }
#undef TEM




template <class T>
struct Mat3{
	
	Mat3(
		const T& v00=T(0), const T& v01=T(0), const T& v02=T(0),
		const T& v10=T(0), const T& v11=T(0), const T& v12=T(0),
		const T& v20=T(0), const T& v21=T(0), const T& v22=T(0)
	)
	:	a00(v00), a01(v01), a02(v02),
		a10(v10), a11(v11), a12(v12),
		a20(v20), a21(v21), a22(v22)
	{}
	
	Mat3& operator()(
		const T& v00, const T& v01, const T& v02,
		const T& v10, const T& v11, const T& v12,
		const T& v20, const T& v21, const T& v22
	){
		a00=v00; a01=v01; a02=v02; a10=v10; a11=v11; a12=v12; a20=v20; a21=v21; a22=v22;
		return *this;		
	}
	
	/// Set element at index with no bounds checking
	T& operator[](uint32_t i){ return elems[i];}

	/// Get element at index with no bounds checking
	const T& operator[](uint32_t i) const { return elems[i]; }
	
	/// Set element at row i, column j
	T& operator()(uint32_t i, uint32_t j){ return elems2[i][j]; }
	
	/// Get element at row i, column j
	const T& operator()(uint32_t i, uint32_t j) const { return elems2[i][j]; }

	Mat3& operator *= (const Mat3& v){
		return (*this)(
			x[0]*v.x[0] + x[1]*v.y[0] + x[2]*v.z[0],
			x[0]*v.x[1] + x[1]*v.y[1] + x[2]*v.z[1],
			x[0]*v.x[2] + x[1]*v.y[2] + x[2]*v.z[2],
			
			y[0]*v.x[0] + y[1]*v.y[0] + y[2]*v.z[0],
			y[0]*v.x[1] + y[1]*v.y[1] + y[2]*v.z[1],
			y[0]*v.x[2] + y[1]*v.y[2] + y[2]*v.z[2],
			
			z[0]*v.x[0] + z[1]*v.y[0] + z[2]*v.z[0],
			z[0]*v.x[1] + z[1]*v.y[1] + z[2]*v.z[1],
			z[0]*v.x[2] + z[1]*v.y[2] + z[2]*v.z[2]
		);
	}
	
	Mat3 operator * (const Mat3& v) const { return Mat3(*this) *= v; }

	Vec3<T> col0() const { return Vec3<T>(a00, a10, a20); }
	Vec3<T> col1() const { return Vec3<T>(a01, a11, a21); }
	Vec3<T> col2() const { return Vec3<T>(a02, a12, a22); }

	template <class S>
	void cols(Vec3<S>& c0, Vec3<S>& c1, Vec3<S>& c2) const {
		c0(a00,a10,a20); c1(a01,a11,a21); c2(a02,a12,a22);
	}

	Vec3<T> row0() const { return Vec3<T>(a00, a01, a02); }
	Vec3<T> row1() const { return Vec3<T>(a10, a11, a12); }
	Vec3<T> row2() const { return Vec3<T>(a20, a21, a22); }

	void rows(Vec3<T>& r0, Vec3<T>& r1, Vec3<T>& r2) const {
		r0(a00,a01,a02); r1(a10,a11,a12); r2(a20,a21,a22);
	}

	union{
		T elems[9];
		T elems2[3][3];
		struct { T x[3], y[3], z[3]; };
		struct{ 
			T a00, a01, a02;
			T a10, a11, a12;
			T a20, a21, a22;
		};
	};
};


/// Returns product of matrix multiplied by column vector
template <class T1, class T2>
Vec3<T2> operator * (const Mat3<T1>& m, const Vec3<T2>& v){
	return Vec3<T2>(
		m(0,0)*v[0] + m(0,1)*v[1] + m(0,2)*v[2],
		m(1,0)*v[0] + m(1,1)*v[1] + m(1,2)*v[2],
		m(2,0)*v[0] + m(2,1)*v[1] + m(2,2)*v[2]
	);
}

/// Returns product of row vector multiplied by matrix
template <class T1, class T2>
Vec3<T1> operator * (const Vec3<T1>& v, const Mat3<T2>& m){
	return Vec3<T1>(
		v[0]*m(0,0) + v[1]*m(1,0) + v[2]*m(2,0),
		v[0]*m(0,1) + v[1]*m(1,1) + v[2]*m(2,1),
		v[0]*m(0,2) + v[1]*m(1,2) + v[2]*m(2,2)
	);
}


/// Multi-element container

/// This is a fixed size array to enable better loop unrolling optimizations
/// by the compiler and to avoid an extra 'size' data member for small-sized
/// arrays. It also lacks a constructor to allow C-style struct initializations.
template <uint32_t N, class T>
struct Multi{
	typedef Multi M;
//	Multi(){}
//	Multi(const T& e ){ mem::set(elems, N, e); }
//	Multi(const T* es){ mem::copy(elems, es, N); }

	T elems[N];
	
	/// Set element at index with no bounds checking
	T& operator[](uint32_t i){ return elems[i];}
	
	/// Get element at index with no bounds checking
	const T& operator[](uint32_t i) const { return elems[i]; }

	#define DO for(uint32_t i=0; i<N; ++i)

	bool operator !=(const M& v){ DO{ if((*this)[i] == v[i]) return false; } return true; }
	bool operator !=(const T& v){ DO{ if((*this)[i] == v   ) return false; } return true; }
	M& operator   = (const M& v){ DO{ (*this)[i] = v[i]; } return *this; }
	M& operator   = (const T& v){ DO{ (*this)[i] = v;    } return *this; }
	bool operator ==(const M& v){ DO{ if((*this)[i] != v[i]) return false; } return true; }
	bool operator ==(const T& v){ DO{ if((*this)[i] != v   ) return false; } return true; }

	#undef DO

	/// Returns size of array
	static uint32_t size(){ return N; }

	/// Zeros all elements.
	void zero(){ memset(elems, 0, N * sizeof(T)); }
};


template <class T>
struct Multi3: public Multi<3,T>{
	Multi3(T v1=0, T v2=0, T v3=0){ (*this)[0]=v1; (*this)[1]=v2; (*this)[2]=v3; }
};


/// A value in the form: base^expo
template <class T=double>
struct ValPower{
	ValPower(const T& base=2, const T& expo=0){ (*this)(base, expo); }
	
	T operator()() const { return mVal; }
	T base() const { return mBase; }
	T expo() const { return mExpo; }
	
	void operator()(const T& base, const T& expo){ mBase=base; mExpo=expo; computeVal(); }
	void base(const T& v){ mBase=v; computeVal(); }
	void expo(const T& v){ mExpo=v; computeVal(); }
	void expoAdd(const T& v){ expo(mExpo+v); }

private:
	T mBase, mExpo, mVal;
	void computeVal(){ mVal=::pow(mBase, mExpo); }
};



/// A closed interval [min, max]

/// An interval is a connected region of the real line. Geometrically, it
/// describes a 0-sphere. Order is strongly enforced so that the endpoints will
/// always satisfy min <= max.
template <class T>
class Interval{
public:

	Interval()
	:	mMin(0), mMax(1){}

	Interval(const T& min, const T& max)
	{ endpoints(min,max); }

	T center() const { return (max()+min())/T(2); }	///< Returns center point

	/// Test is point is contained exclusively within interval
	bool contains(const T& v) const { return v>=min() && v<=max(); }

	bool degenerate() const { return min()==max(); }///< Returns true if diameter is zero
	T diameter() const { return max()-min(); }		///< Returns absolute difference of endpoints
	const T& max() const { return mMax; }			///< Get maximum endpoint
	const T& min() const { return mMin; }			///< Get minimum endpoint
	bool proper() const { return min()!=max(); }	///< Returns true if diameter is non-zero
	T radius() const { return diameter()/T(2); }	///< Returns one-half the diameter

	/// Linearly map point in interval to point in the unit interval
	T toUnit(const T& v) const { return (v-min())/diameter(); }
	
	template <class U>
	bool operator == (const Interval<U>& v){ return min()==v.min() && max()==v.max(); }

	template <class U>
	bool operator != (const Interval<U>& v){ return !(*this == v); }
	
	template <class U>
	Interval& operator +=(const Interval<U>& v){ endpoints(min()+v.min(), max()+v.max()); return *this; }

	template <class U>
	Interval& operator -=(const Interval<U>& v){ endpoints(min()-v.max(), max()-v.min()); return *this; }
	
	template <class U>
	Interval& operator *=(const Interval<U>& v){
		T a=min()*v.min(), b=min()*v.max(), c=max()*v.min(), d=max()*v.max();
		mMin = min(min(a,b),min(c,d));
		mMax = max(max(a,b),max(c,d));
		return *this;
	}

	template <class U>
	Interval& operator /=(const Interval<U>& v){
		T a=min()/v.min(), b=min()/v.max(), c=max()/v.min(), d=max()/v.max();
		mMin = min(min(a,b),min(c,d));
		mMax = max(max(a,b),max(c,d));
		return *this;
	}

	/// Set center point preserving diameter
	Interval& center(const T& v){ return centerDiameter(v, diameter()); }

	/// Set diameter (width) preserving center
	Interval& diameter(const T& v){ return centerDiameter(center(), v); }

	/// Set center and diameter
	Interval& centerDiameter(const T& c, const T& d){
		mMin = c - d*T(0.5);
		mMax = mMin + d;
		return *this;
	}

	/// Set the endpoints
	Interval& endpoints(const T& min, const T& max){
		mMax=max; mMin=min;
		if(mMin > mMax){ T t=mMin; mMin=mMax; mMax=t; }
		return *this;
	}

	/// Translate interval by fixed amount
	Interval& translate(const T& v){ mMin+=v; mMax+=v; return *this; }

	/// Set maximum endpoint
	Interval& max(const T& v){ return endpoints(min(), v); }
	
	/// Set minimum endpoint
	Interval& min(const T& v){ return endpoints(v, max()); }

private:
	T mMin, mMax;

	const T& min(const T& a, const T& b){ return a<b?a:b; }
	const T& max(const T& a, const T& b){ return a>b?a:b; }
};


/// A value wrapped to an interval [min, max)

/// Mathematical correctness is strongly enforced. The value will always lie in 
/// the specified interval.
template <class T>
class ValWrap : public Interval<T>{
public:
	typedef Interval<T> I;
	typedef ValWrap<T> V;

	ValWrap(): I(), mVal(0){}

	ValWrap(const T& max, const T& min=T(0), const T& v=T(0))
	: I(min, max), mVal(v)
	{}
	
	ValWrap& operator= (const T& v){ return val(v); }				///< Set value
	ValWrap& operator+=(const T& v){ return val(mVal+v); }			///< Add value
	//ValWrap& operator+=(const ValWrap& v){ return (*this)+=v(); }	///< Add wrapped value
	ValWrap& operator-=(const T& v){ return val(mVal-v); }			///< Subtract value
	//ValWrap& operator-=(const ValWrap& v){ return (*this)-=v(); }	///< Subtract wrapped value

	ValWrap& operator*=(const T& v){ return val(mVal*v); }			///< Multiply value
	//ValWrap& operator*=(const ValWrap& v){ return (*this)*=v(); }	///< Multiply wrapped value
	ValWrap& operator/=(const T& v){ return val(mVal/v); }			///< Divide value
	//ValWrap& operator/=(const ValWrap& v){ return (*this)/=v(); }	///< Divide wrapped value
	
	ValWrap  operator+ (const T& v) const { return V(*this)+=v; }
	ValWrap  operator- (const T& v) const { return V(*this)-=v; }
	ValWrap  operator* (const T& v) const { return V(*this)*=v; }
	ValWrap  operator/ (const T& v) const { return V(*this)/=v; }
	
	ValWrap  operator++(int){ V r=*this; ++(*this); return r; }		///< Postfix increment value
	ValWrap  operator--(int){ V r=*this; --(*this); return r; }		///< Postfix decrement value
	ValWrap& operator++(){ return val(++mVal); }					///< Prefix increment value
	ValWrap& operator--(){ return val(--mVal); }					///< Prefix decrement value

	/// Set wrapping interval
	ValWrap& endpoints(const T& min, const T& max){
		I::endpoints(min,max);
		return val(val());
	}

	ValWrap& max(const T& v){ I::max(v); return val(val()); }
	ValWrap& min(const T& v){ I::min(v); return val(val()); }

	ValWrap& val(T v){
	
		if(I::proper()){
			if(!(v < I::max())){
				T d = I::diameter();
				v -= d;
				if(!(v < I::max())) v -= d * (T)(uint32_t)((v - I::min())/d);
			}
			else if(v < I::min()){
				T d = I::diameter();
				v += d;
				if(v < I::min()) v += d * (T)(uint32_t)(((I::min() - v)/d) + 1);
			}
			mVal = v;
		}
		else{
			mVal = I::min();
		}
		return *this;
	}


	const T& operator()() const { return mVal; }

	/// Returns positive unit fractional position in interval
	double fraction() const {
		return I::proper() ? double(val())/I::diameter() : 0.;
	}

	const T& val() const { return mVal; }
	
	operator T() const { return mVal; }
	
protected:
	T mVal;
};


/// Fixed size shift buffer
template <int N, class T>
struct ShiftBuffer : public Multi<N,T>{

	typedef Multi<N,T> base;
	using base::elems;
	using base::operator=;

	ShiftBuffer(const T& v=T()){ *this = v; }

	/// Push new element onto buffer. Newest element is at index 0.
	void operator()(const T& v){
		for(int i=N-1; i>0; --i) elems[i] = elems[i-1];
		elems[0]=v;
	}
};




///< Fixed size vector.
template <uint32_t N, class T>
struct Vec : public Multi<N,T> {

	typedef Vec<N,T> V;

	Vec(const T& v=T(0)){ (*this) = v; }
	Vec(const V& v){ (*this) = v; }

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

	T dot(const V& v) const { T r=(T)0; DO r+=(*this)[i]*v[i]; return r; }
	T norm() const { return sqrt(norm2()); }
	T norm2() const { return dot(*this); }
	V& normalize(){ return *this /= norm(); }

	V sgn() const { return V(*this).normalize(); }

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
	using Vec<3,T>::operator*;
	using Vec<3,T>::operator*=;

	Vec3(const Vec<3, T>& v){ (*this)=v; }
	Vec3(const T& v=T(0)){ (*this)(v,v,v); }
	Vec3(const T& v1, const T& v2, const T& v3=T(0)){ (*this)(v1,v2,v3); }

	Vec3& operator()(const T& v1, const T& v2, const T& v3){ (*this)[0]=v1; (*this)[1]=v2; (*this)[2]=v3; return *this; }
	
	Vec3& operator*= (const Mat3<T>& m){
		Vec3& t = *this;
		return t(
			t[0]*m[0] + t[1]*m[1] + t[2]*m[2],
			t[0]*m[3] + t[1]*m[4] + t[2]*m[5],
			t[0]*m[6] + t[1]*m[7] + t[2]*m[8]			
		);
	}
	
	Vec3 operator * (const Mat3<T>& m) const { return Vec3(*this) *= m; }

	Vec3 cross(const Vec3& v) const {
		Vec3 r; const Vec3& t = *this;
		r[0] = t[1]*v[2] - t[2]*v[1];
		r[1] = t[2]*v[0] - t[0]*v[2];
		r[2] = t[0]*v[1] - t[1]*v[0];
		return r;
	}

	T dot() const { return dot(*this); }
	T dot(const Vec3& v) const { return v[0]*(*this)[0] + v[1]*(*this)[1] + v[2]*(*this)[2]; }
	
	Vec3& rotateX(const Complex<T>& v){ Complex<T> t(y(),z()); t*=v; return (*this)(x(),t.r,t.i); }
	Vec3& rotateY(const Complex<T>& v){ Complex<T> t(z(),x()); t*=v; return (*this)(t.i,y(),t.r); }
	Vec3& rotateZ(const Complex<T>& v){ Complex<T> t(x(),y()); t*=v; return (*this)(t.r,t.i,z()); }
	
	const T& x() const { return (*this)[0]; }
	const T& y() const { return (*this)[1]; }
	const T& z() const { return (*this)[2]; }
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


/// Rotate vector towards perpendicular vector by angle using right-hand rule.
template <class T>
Vec3<T> rotate(const Vec3<T>& v, const Vec3<T>& p, const Complex<T>& a){
	return v*a.r + p*a.i;
}

template <class VecN, class T>
VecN rotate(const VecN& v, const VecN& p, const Complex<T>& a){
	return v*a.r + p*a.i;
}

/// Rotates two vectors by angle in plane formed from bivector v1 ^ v2
template <class VecN, class T>
void rotatePlane(VecN& v1, VecN& v2, const Complex<T>& a){
	VecN t = gam::rotate(v1, v2, a);
	v2 = gam::rotate(v2, VecN(-v1), a);
	v1 = t;
}


template <class T>
Vec3<T> rotateX(const Vec3<T>& v, const Complex<T>& a){
	Complex<T> t(v[1], v[2]); t*=a; return Vec3<T>(v[0], t[0], t[1]);
}

template <class T>
Vec3<T> rotateY(const Vec3<T>& v, const Complex<T>& a){
	Complex<T> t(v[2], v[0]); t*=a; return Vec3<T>(t[1], v[1], t[0]);
}

template <class T>
Vec3<T> rotateZ(const Vec3<T>& v, const Complex<T>& a){
	Complex<T> t(v[0], v[1]); t*=a; return Vec3<T>(t[0], t[1], v[2]);
}



template<class T> inline T norm(const Complex<T>& v){ return v.norm(); }
template<class T> inline T normCompare(const Complex<T>& v){ return v.norm2(); }
template<class T> inline T norm(const Quat<T>& v){ return v.norm(); }
template<class T> inline T normCompare(const Quat<T>& v){ return v.norm2(); }
template<int N,class T> inline T norm(const Vec<N,T>& v){ return v.norm(); }
template<int N,class T> inline T normCompare(const Vec<N,T>& v){ return v.norm2(); }


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


} // gam::

#endif
