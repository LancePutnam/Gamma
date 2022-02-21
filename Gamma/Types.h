#ifndef GAMMA_TYPES_H_INC
#define GAMMA_TYPES_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information 

	File Description: 
	Complex numbers and n-vectors.
*/

#include <cmath>
#include <initializer_list>

namespace gam{

using std::cos;
using std::sin;

template<class T> class Complex;
template<unsigned N, class T> class Vec;

typedef float real;				///< Default real number type
typedef Vec<2,float > float2;	///< Vector of 2 floats
typedef Vec<2,double> double2;	///< Vector of 2 doubles
typedef Vec<3,float > float3;	///< Vector of 3 floats
typedef Vec<3,double> double3;	///< Vector of 3 doubles
typedef Vec<4,float > float4;	///< Vector of 4 floats
typedef Vec<4,double> double4;	///< Vector of 4 doubles

#define IT(n) for(unsigned i=0; i<(n); ++i)


/// Polar number with argument in radians
template <class T>
struct Polar{

	typedef T value_type;

	union{
		struct{ T m, p; };	///< Magnitude and phase values
		T elems[2];			///< Component 2-vector
	};

	Polar(const T& p=T(0)): m(T(1)), p(p){}
	Polar(const T& m, const T& p): m(m), p(p){}
	Polar(const Complex<T>& v){ *this = v; }

	T& operator[](unsigned i){ return elems[i];}
	const T& operator[](unsigned i) const { return elems[i]; }

	Polar& operator = (const Complex<T>& v){ m=v.norm(); p=v.arg(); return *this; }
};


/// Complex number
template <class T=gam::real>
class Complex{
public:

	typedef Complex<T> C;
	typedef T value_type;

	union{
		struct{ T r, i; };	///< Real and imaginary values
		T elems[2];			///< Component 2-vector
	};

	Complex(const Complex& v): r(v.r), i(v.i){}
	Complex(const Polar<T>& v){ *this = v; }
	Complex(const T& r=T(0), const T& i=T(0)): r(r), i(i){}
	Complex(const T& m, const T& p, int fromPolar){ (*this) = Polar<T>(m,p); }


	// Accessors compatible with std::complex
	T& real(){return r;}
	const T& real() const {return r;}
	T& imag(){return i;}
	const T& imag() const {return i;}


	C& arg(const T& v){ return fromPolar(norm(), v); }					///< Set argument leaving norm the same
	C& norm(const T& v){ return fromPolar(v, arg()); }					///< Set norm leaving argument the same
	C& mag(const T& v){ return norm(v); }

	C& fromPhase(const T& v){ r=cos(v); i=sin(v); return *this; }		///< Set phase and normalize
	C& fromPolar(const T& m, const T& p){ return (*this)(Polar<T>(m,p)); }	///< Set magnitude and phase

	template <class U>
	C& operator = (const Polar<U>& v){ r=v.m*cos(v.p); i=v.m*sin(v.p); return *this; }

	template <class U>
	C& operator = (const Complex<U>& v){ r=v.r; i=v.i; return *this; }

	C& operator = (const T& v){ r=v; i=T(0); return *this; }

	template <class U>
	C& set(const Complex<U>& v){ return *this = v; }

	C& set(const T& vr, const T& vi){ r=vr; i=vi; return *this; }

	template <class U>
	C& set(const Polar<U>& p){ return *this = p; }


	T& operator[](unsigned i){ return elems[i];}
	const T& operator[](unsigned i) const { return elems[i]; }

	bool operator ==(const C& v) const { return (r==v.r) && (i==v.i); }		///< Returns true if all components are equal
	bool operator ==(const T& v) const { return (r==v  ) && (i==T(0));}		///< Returns true if real and equals value
	bool operator !=(const C& v) const { return (r!=v.r) || (i!=v.i); }		///< Returns true if any components are not equal
	bool operator > (const C& v) const { return normSqr() > v.normSqr(); }	///< Returns true if norm is greater than argument's norm
	bool operator < (const C& c) const { return normSqr() < c.normSqr(); }	///< Returns true if norm is less than argument's norm

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
	T argUnit() const { T r=arg()/(2*M_PI); return r>0 ? r : r+1; }	///< Return argument in unit interval [0, 1)
	C conj() const { return C(r,-i); }						///< Returns conjugate, z*
	T dot(const C& v) const { return r*v.r + i*v.i; }		///< Returns vector dot product
	T norm() const { return sqrt(normSqr()); }				///< Returns norm (radius), |z|
	T normSqr() const { return dot(*this); }				///< Returns square of norm, |z|^2
	C& normalize(T m=T(1)){ return *this *= (m/norm()); }	///< Sets norm (radius) to 1, |z|=1
	C recip() const { return conj()/normSqr(); }			///< Return multiplicative inverse, 1/z
	C sgn(T m=T(1)) const { return C(*this).normalize(m); }	///< Returns signum, z/|z|, the closest point on unit circle
	C sqr() const { return C(r*r-i*i, T(2)*r*i); }			///< Returns square

	T mag() const { return norm(); }						///< Returns norm (radius), |z|
	T magSqr() const { return normSqr(); }					///< Returns magnitude squared, |z|^2
	T phase() const { return arg(); }						///< Returns argument (angle)


	// deprecated
	C& operator()(const T& vr, const T& vi){ return set(vr,vi); }
	C& operator()(const Polar<T>& p){ return set(p); }
};

#define TEM template <class T>
TEM Complex<T> abs(const Complex<T>& z){ return z.norm(); }
TEM Complex<T> cos(const Complex<T>& z){ return Complex<T>(cos(z.r)*cosh(z.i),-sin(z.r)*sinh(z.i)); }
TEM Complex<T> sin(const Complex<T>& z){ return Complex<T>(sin(z.r)*cosh(z.i), cos(z.r)*sinh(z.i)); }
TEM Complex<T> exp(const Complex<T>& z){ return Polar<T>(exp(z.r), z.i); }
TEM Complex<T> log(const Complex<T>& z){ return Complex<T>(T(0.5)*log(z.normSqr()), z.arg()); }
TEM Complex<T> pow(const Complex<T>& b, const Complex<T>& e){ return exp(log(b)*e); }
TEM Complex<T> pow(const Complex<T>& b, const T& e){ return exp(log(b)*e); }
TEM Complex<T> sqrt(const Complex<T>& z){
	static const T c = T(1)/sqrt(T(2));
	T n = z.norm();
	T a = sqrt(n + z.r) * c;
	T b = sqrt(n - z.r) * (z.i<T(0) ? -c : c);
	return Complex<T>(a,b);
}
TEM Complex<T> operator + (T r, const Complex<T>& c){ return  c+r; }
TEM Complex<T> operator - (T r, const Complex<T>& c){ return -c+r; }
TEM Complex<T> operator * (T r, const Complex<T>& c){ return  c*r; }
TEM Complex<T> operator / (T r, const Complex<T>& c){ return  c.conj()*(r/c.normSqr()); }
#undef TEM


template <unsigned N, class T> struct NamedElems{ union{ T x; T mElems[N]; }; };
template<class T> struct NamedElems<0,T>{ static T x; };
template<class T> struct NamedElems<1,T>{ T x; };
template<class T> struct NamedElems<2,T>{ T x,y; };
template<class T> struct NamedElems<3,T>{ T x,y,z; };
template<class T> struct NamedElems<4,T>{ T x,y,z,w; };



/// N-vector or fixed-size array

/// This is fixed in size to enable better loop unrolling optimizations and to 
/// avoid an extra 'size' data member for small sizes.
template <unsigned N, class T>
class Vec : public NamedElems<N,T> {
public:

	typedef T value_type;

	Vec(const T& v=T()){ set(v); }
	Vec(const T& v1, const T& v2){ set(v1,v2); }
	Vec(const T& v1, const T& v2, const T& v3){ set(v1,v2,v3); }
	Vec(const T& v1, const T& v2, const T& v3, const T& v4){ set(v1,v2,v3,v4); }

	template <class U>
	Vec(const U * src){ set(src); }

	template <unsigned N2, class T2>
	Vec(const Vec<N2, T2>& v){ set(v); }

	template <class Tv, class Ts>
	Vec(const Vec<N-1, Tv>& v, Ts s){ set(v,s);}


	/// Returns size of vector
	static constexpr unsigned size(){ return N; }

	T * elems(){ return &(this->x); }
	const T * elems() const { return &(this->x); }

	T * begin(){ return elems(); }
	const T * begin() const { return elems(); }
	T * end(){ return begin()+N; }
	const T * end() const { return begin()+N; }

	/// Set element at index (no bounds checking)
	T& operator[](unsigned i){ return elems()[i];}

	/// Get element at index (no bounds checking)
	const T& operator[](unsigned i) const { return elems()[i]; }

	/// Get a vector comprised of indexed elements
	Vec<2,T> get(int i0, int i1) const {
		return Vec<2,T>((*this)[i0], (*this)[i1]); }

	/// Get a vector comprised of indexed elements
	Vec<3,T> get(int i0, int i1, int i2) const {
		return Vec<3,T>((*this)[i0], (*this)[i1], (*this)[i2]); }

	/// Get a vector comprised of indexed elements
	Vec<4,T> get(int i0, int i1, int i2, int i3) const {
		return Vec<4,T>((*this)[i0], (*this)[i1], (*this)[i2], (*this)[i3]); }

	/// Get a subvector
	template <int M, int Begin=0>
	const Vec<M,T>& sub() const {
		return const_cast<Vec*>(this)->sub<M,Begin>();
	}
	template <int M, int Begin=0>
	Vec<M,T>& sub(){
		static_assert((Begin+M)<=N, "Invalid subvector range");
		return *(Vec<M,T>*)(elems()+Begin);
	}

	bool operator !=(const Vec& v){ IT(N){ if((*this)[i] == v[i]) return false; } return true; }
	bool operator !=(const   T& v){ IT(N){ if((*this)[i] == v   ) return false; } return true; }
	Vec& operator = (const Vec& v){ IT(N) (*this)[i] = v[i]; return *this; }
	Vec& operator = (const   T& v){ IT(N) (*this)[i] = v; return *this; }
	bool operator ==(const Vec& v){ IT(N){ if((*this)[i] != v[i]) return false; } return true; }
	bool operator ==(const   T& v){ IT(N){ if((*this)[i] != v   ) return false; } return true; }

	Vec  operator * (const Vec& v) const { Vec r; IT(N) r[i] = (*this)[i] * v[i]; return r; }
	Vec  operator * (const   T& v) const { Vec r; IT(N) r[i] = (*this)[i] * v;    return r; }
	Vec& operator *=(const Vec& v){ IT(N) (*this)[i] *= v[i]; return *this; }
	Vec& operator *=(const   T& v){ IT(N) (*this)[i] *= v;    return *this; }
	Vec  operator / (const Vec& v) const { Vec r; IT(N) r[i] = (*this)[i] / v[i]; return r; }
	Vec  operator / (const   T& v) const { Vec r; IT(N) r[i] = (*this)[i] / v;    return r; }
	Vec& operator /=(const Vec& v){ IT(N) (*this)[i] /= v[i]; return *this; }
	Vec& operator /=(const   T& v){ IT(N) (*this)[i] /= v;    return *this; }
	Vec  operator - (          ) const { Vec r; IT(N) r[i] = -(*this)[i]; return r; }
	Vec  operator - (const Vec& v) const { Vec r; IT(N) r[i] = (*this)[i] - v[i]; return r; }
	Vec  operator - (const   T& v) const { Vec r; IT(N) r[i] = (*this)[i] - v;    return r; }
	Vec& operator -=(const Vec& v){ IT(N) (*this)[i] -= v[i]; return *this; }
	Vec& operator -=(const   T& v){ IT(N) (*this)[i] -= v;    return *this; }
	Vec  operator + (const Vec& v) const { Vec r; IT(N) r[i] = (*this)[i] + v[i]; return r; }
	Vec  operator + (const   T& v) const { Vec r; IT(N) r[i] = (*this)[i] + v;    return r; }
	Vec& operator +=(const Vec& v){ IT(N) (*this)[i] += v[i]; return *this; }
	Vec& operator +=(const   T& v){ IT(N) (*this)[i] += v;    return *this; }

	template <class V, class Func, class... Args>
	Vec<N,V> map(Func func, Args... args) const {
		Vec<N,V> r;
		for(int i=0; i<size(); ++i)
			r[i] = func((*this)[i], args...);
		return r;
	}

	/// Map elements through function into new vector

	/// @param[in] func		Function taking old value and returning new value
	/// @param[in] args		Extra function arguments
	template <class Func, class... Args>
	Vec map(Func func, Args... args) const {
		return map<T>(func, args...);
	}

	/// Reduce elements into scalar

	/// @param[in] prev		Initial previous value
	/// @param[in] func		Function taking previous and current values as first 
	///						two arguments and returning new value
	/// @param[in] args		Extra function arguments
	template <class Func, class... Args>
	T reduce(const T& prev, Func func, Args... args) const {
		T r = prev;
		for(auto& v : *this) r = func(r, v, args...);
		return r;
	}

	/// Zeros all elements
	void zero(){ memset(elems(), 0, N * sizeof(T)); }

	T dot(const Vec& v) const { T r=T(0); IT(N) r+=(*this)[i]*v[i]; return r; }
	T sum() const { T r=T(0); IT(N) r+=(*this)[i]; return r; }
	T mag() const { return std::sqrt(magSqr()); }
	T magSqr() const { return dot(*this); }
	Vec normalized() const { return Vec(*this).normalize(); }

	Vec& normalize(){
		T msqr = magSqr();
		if(msqr > 0)	return (*this) /= sqrt(msqr);
		else			return setIdentity();
	}

	template <unsigned N2, class T2>
	Vec& set(const Vec<N2, T2>& v){ IT(N<N2?N:N2){ (*this)[i] = T(v[i]); } return *this; }

	template <class Tv, class Ts>
	Vec& set(const Vec<N-1, Tv>& v, Ts s){ (*this)[N-1]=s; return set(v); }

	/// Set all elements to the same value
	Vec& set(const T& v){ return (*this = v); }

	/// Set first 2 elements
	Vec& set(const T& v1, const T& v2){
		return set(v1,v2,v1,v1,v1,v1); }

	/// Set first 3 elements
	Vec& set(const T& v1, const T& v2, const T& v3){
		return set(v1,v2,v3,v1,v1,v1); }

	/// Set first 4 elements
	Vec& set(const T& v1, const T& v2, const T& v3, const T& v4){
		return set(v1,v2,v3,v4,v1,v1); }

	/// Set first 5 elements
	Vec& set(const T& v1, const T& v2, const T& v3, const T& v4, const T& v5){
		return set(v1,v2,v3,v4,v5,v1); }

	/// Set first 6 elements
	Vec& set(const T& v1, const T& v2, const T& v3, const T& v4, const T& v5, const T& v6){		
		switch(N){
		default:(*this)[5] = v6;
		case 5: (*this)[4] = v5;
		case 4: (*this)[3] = v4;
		case 3: (*this)[2] = v3;
		case 2: (*this)[1] = v2;
		case 1: (*this)[0] = v1;
		}
		return *this;
	}

	/// Set elements to values from C array
	template <class U>
	Vec& set(const U * src){ IT(N){ (*this)[i]=src[i]; } return *this; }

	/// Set to identity, i.e., {1, 0, ..., 0}
	Vec& setIdentity(){
		(*this)[0] = T(1);
		for(unsigned i=1; i<N; ++i) (*this)[i] = T(0);
		return *this;
	}
};


template <unsigned N, class T, class S>
inline Vec<N,T> operator + (const S& s, const Vec<N,T>& v){ return  v+s; }

template <unsigned N, class T, class S>
inline Vec<N,T> operator - (const S& s, const Vec<N,T>& v){ return -v+s; }

template <unsigned N, class T, class S>
inline Vec<N,T> operator * (const S& s, const Vec<N,T>& v){ return  v*s; }

template <unsigned N, class T, class S>
inline Vec<N,T> operator / (const S& s, const Vec<N,T>& v){
	Vec<N,T> r; IT(N){ r[i] = s/v[i]; } return r;
}

#define DEF_CVOP(op)\
template <unsigned N, class T, class U>\
inline Vec<N,T> operator op (const Vec<N,T>& v, const Complex<U>& c){\
	Vec<N,T> r; for(auto i:{0,1}) r[i]=v[i] op c[i]; return r;\
}\
template <unsigned N, class T, class U>\
inline Complex<U> operator op (const Complex<U>& c, const Vec<N,T>& v){\
	Complex<U> r; for(auto i:{0,1}) r[i]=c[i] op v[i]; return r;\
}\
template <unsigned N, class T, class U>\
inline Vec<N,T>& operator op##= (Vec<N,T>& v, const Complex<U>& c){\
	for(auto i:{0,1}) v[i] op##= c[i]; return v;\
}\
template <unsigned N, class T, class U>\
inline Complex<U>& operator op##= (Complex<U>& c, const Vec<N,T>& v){\
	for(auto i:{0,1}) c[i] op##= v[i]; return c;\
}
DEF_CVOP(+) DEF_CVOP(-)
#undef DEF_CVOP

namespace scl{

template<unsigned N, class T>
inline Vec<N,T> abs(Vec<N,T> a){
	Vec<N,T> r;
	for(unsigned i=0; i<N; ++i) r[i] = abs(a[i]);
	return r;
}

template<unsigned N, class T, class U>
inline Vec<N,T> max(Vec<N,T> a, Vec<N,U> b){
	Vec<N,T> r;
	for(unsigned i=0; i<N; ++i) r[i] = max(a[i], b[i]);
	return r;
}

template<unsigned N, class T, class U>
inline Vec<N,T> min(Vec<N,T> a, Vec<N,U> b){
	Vec<N,T> r;
	for(unsigned i=0; i<N; ++i) r[i] = min(a[i], b[i]);
	return r;
}

}


/// Returns cross product, a x b
template <class T>
inline Vec<3,T> cross(const Vec<3,T>& a, const Vec<3,T>& b){
	return Vec<3,T>(
		a[1]*b[2] - a[2]*b[1],
		a[2]*b[0] - a[0]*b[2],
		a[0]*b[1] - a[1]*b[0]
	);
}

/// Returns spherical Euler product (ZXZ convention)
template <class T>
Vec<3,T> productZXZ(const Complex<T>& a, const Complex<T>& b, const Complex<T>& c){
	return Vec<3,T>(
		a.r*b.i - a.i*b.r*c.i,
		a.i*b.i + a.r*b.r*c.i,
		b.r*c.r
	);
}

/// Rotate vector towards perpendicular vector by angle using right-hand rule.
template <int N, class T>
Vec<N,T> rotate(const Vec<N,T>& v, const Vec<N,T>& p, const Complex<T>& a){
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
Vec<3,T> rotateX(const Vec<3,T>& v, const Complex<T>& a){
	Complex<T> t(v[1], v[2]); t*=a; return Vec<3,T>(v[0], t[0], t[1]);
}

template <class T>
Vec<3,T> rotateY(const Vec<3,T>& v, const Complex<T>& a){
	Complex<T> t(v[2], v[0]); t*=a; return Vec<3,T>(t[1], v[1], t[0]);
}

template <class T>
Vec<3,T> rotateZ(const Vec<3,T>& v, const Complex<T>& a){
	Complex<T> t(v[0], v[1]); t*=a; return Vec<3,T>(t[0], t[1], v[2]);
}

template<class T> inline T magSqr(const Complex<T>& v){ return v.magSqr(); }
template<class T> inline double norm(const Complex<T>& v){ return v.norm(); }
template<class T> inline double normCompare(const Complex<T>& v){ return v.normSqr(); }

template<int N,class T> inline T magSqr(const Vec<N,T>& v){ return v.magSqr(); }
template<int N,class T> inline double norm(const Vec<N,T>& v){ return v.mag(); }
template<int N,class T> inline double normCompare(const Vec<N,T>& v){ return v.magSqr(); }

#undef IT

} // gam::

#endif
