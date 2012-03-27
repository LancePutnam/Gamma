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
template<uint32_t N, class T> class Vec;

typedef Vec<2,float > float2;
typedef Vec<2,double> double2;

//typedef Polar<float > Polarf;
//typedef Polar<double> Polard;
//typedef Complex<float > Complexf;
//typedef Complex<double> Complexd;
//typedef Vec2<float> Vec2f;
//typedef Vec2<double> Vec2d;
//typedef Vec3<float> Vec3f;
//typedef Vec3<double> Vec3d;
//typedef Vec4<float> Vec4f;
//typedef Vec4<double> Vec4d;


/// Polar number with argument in radians
template <class T>
struct Polar{

	union{
		struct{ T m, p; };	///< Magnitude and phase values
		T elems[2];			///< Component 2-vector
	};

	Polar(const T& p=0): m(1.), p(p){}
	Polar(const T& m, const T& p): m(m), p(p){}
	Polar(const Complex<T>& v){ *this = v; }

	T& operator[](uint32_t i){ return elems[i];}
	const T& operator[](uint32_t i) const { return elems[i]; }

	Polar& operator = (const Complex<T>& v){ m=v.norm(); p=v.arg(); return *this; }
};


/// Complex number
template <class T=gam::real>
struct Complex{

	typedef Complex<T> C;

	union{
		struct{ T r, i; };	///< Real and imaginary values
		T elems[2];			///< Component 2-vector
	};
	
	Complex(const Complex& v): r(v.r), i(v.i){}
	Complex(const Polar<T>& v){ *this = v; }
	Complex(const T& r=(T)1, const T& i=(T)0): r(r), i(i){}
	Complex(const T& m, const T& p, int fromPolar){ (*this) = Polar<T>(m,p); }


	// Accessors compatible with std::complex
	T& real(){return r;}
	const T& real() const {return r;}
	T& imag(){return i;}
	const T& imag() const {return i;}


	C& arg(const T& v){ return fromPolar(norm(), v); }					///< Set phase leaving magnitude the same
	C& fromPhase(const T& v){ r=::cos(v); i=::sin(v); return *this; }	///< Set phase and normalize
	C& fromPolar(const T& m, const T& p){ return (*this)(Polar<T>(m,p)); }	///< Set magnitude and phase
	C& mag(const T& v){ return norm(v); }
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
	T magSqr() const { return norm2(); }					///< Returns magnitude squared, |z|^2
	T phase() const { return arg(); }						///< Returns argument (angle)
};

#define TEM template <class T>
TEM Complex<T> abs(const Complex<T>& c){ return c.abs(); }
TEM Complex<T> exp(const Complex<T>& c){ return c.exp(); }
TEM Complex<T> log(const Complex<T>& c){ return c.log(); }
TEM Complex<T> pow(const Complex<T>& b, const Complex<T>& e){ return b.pow(e); }
TEM Complex<T> operator + (T r, const Complex<T>& c){ return  c+r; }
TEM Complex<T> operator - (T r, const Complex<T>& c){ return -c+r; }
TEM Complex<T> operator * (T r, const Complex<T>& c){ return  c*r; }
TEM Complex<T> operator / (T r, const Complex<T>& c){ return  c.conj()*(r/c.norm()); }
#undef TEM



/// Multi-element container

/// This is a fixed size array to enable better loop unrolling optimizations
/// by the compiler and to avoid an extra 'size' data member for small-sized
/// arrays. It intentionally lacks a constructor to allow C-style struct 
/// initializations.
template <uint32_t N, class T>
struct Multi{
	typedef Multi M;

	T elems[N];
	
	/// Set element at index with no bounds checking
	T& operator[](uint32_t i){ return elems[i];}
	
	/// Get element at index with no bounds checking
	const T& operator[](uint32_t i) const { return elems[i]; }

	#define IT for(uint32_t i=0; i<N; ++i)

	bool operator !=(const M& v){ IT{ if((*this)[i] == v[i]) return false; } return true; }
	bool operator !=(const T& v){ IT{ if((*this)[i] == v   ) return false; } return true; }
	M& operator   = (const M& v){ IT{ (*this)[i] = v[i]; } return *this; }
	M& operator   = (const T& v){ IT{ (*this)[i] = v;    } return *this; }
	bool operator ==(const M& v){ IT{ if((*this)[i] != v[i]) return false; } return true; }
	bool operator ==(const T& v){ IT{ if((*this)[i] != v   ) return false; } return true; }

	#undef IT

	/// Returns size of array
	static uint32_t size(){ return N; }

	/// Zeros all elements.
	void zero(){ memset(elems, 0, N * sizeof(T)); }
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




/// Fixed-size vector
template <uint32_t N, class T>
struct Vec : public Multi<N,T> {

	typedef Vec<N,T> V;

	Vec(const T& v=T()){ set(v); }
	Vec(const T& v1, const T& v2){ set(v1,v2); }
	Vec(const T& v1, const T& v2, const T& v3){ set(v1,v2,v3); }
	Vec(const T& v1, const T& v2, const T& v3, const T& v4){ set(v1,v2,v3,v4); }

	template <uint32_t N2, class T2>
	Vec(const Vec<N2, T2>& v){ set(v); }

	/// Get a vector comprised of indexed elements
	Vec<2,T> get(int i0, int i1) const {
		return Vec<2,T>((*this)[i0], (*this)[i1]); }

	/// Get a vector comprised of indexed elements
	Vec<3,T> get(int i0, int i1, int i2) const {
		return Vec<3,T>((*this)[i0], (*this)[i1], (*this)[i2]); }

	/// Get a vector comprised of indexed elements
	Vec<4,T> get(int i0, int i1, int i2, int i3) const {
		return Vec<4,T>((*this)[i0], (*this)[i1], (*this)[i2], (*this)[i3]); }

	#define IT(n) for(uint32_t i=0; i<n; ++i)

	bool operator !=(const V& v){ IT(N){ if((*this)[i] == v[i]) return false; } return true; }
	bool operator !=(const T& v){ IT(N){ if((*this)[i] == v   ) return false; } return true; }
	V& operator = (const V& v){ IT(N) (*this)[i] = v[i]; return *this; }
	V& operator = (const T& v){ IT(N) (*this)[i] = v;    return *this; }
	bool operator ==(const V& v){ IT(N){ if((*this)[i] != v[i]) return false; } return true; }
	bool operator ==(const T& v){ IT(N){ if((*this)[i] != v   ) return false; } return true; }
	V  operator * (const V& v) const { V r; IT(N) r[i] = (*this)[i] * v[i]; return r; }
	V  operator * (const T& v) const { V r; IT(N) r[i] = (*this)[i] * v;    return r; }
	V& operator *=(const V& v){ IT(N) (*this)[i] *= v[i]; return *this; }
	V& operator *=(const T& v){ IT(N) (*this)[i] *= v;    return *this; }
	V  operator / (const V& v) const { V r; IT(N) r[i] = (*this)[i] / v[i]; return r; }
	V  operator / (const T& v) const { V r; IT(N) r[i] = (*this)[i] / v;    return r; }
	V& operator /=(const V& v){ IT(N) (*this)[i] /= v[i]; return *this; }
	V& operator /=(const T& v){ IT(N) (*this)[i] /= v;    return *this; }
	V  operator - (          ) const { V r; IT(N) r[i] = -(*this)[i]; return r; }
	V  operator - (const V& v) const { V r; IT(N) r[i] = (*this)[i] - v[i]; return r; }
	V  operator - (const T& v) const { V r; IT(N) r[i] = (*this)[i] - v;    return r; }
	V& operator -=(const V& v){ IT(N) (*this)[i] -= v[i]; return *this; }
	V& operator -=(const T& v){ IT(N) (*this)[i] -= v;    return *this; }
	V  operator + (const V& v) const { V r; IT(N) r[i] = (*this)[i] + v[i]; return r; }
	V  operator + (const T& v) const { V r; IT(N) r[i] = (*this)[i] + v;    return r; }
	V& operator +=(const V& v){ IT(N) (*this)[i] += v[i]; return *this; }
	V& operator +=(const T& v){ IT(N) (*this)[i] += v;    return *this; }


	T dot(const V& v) const { T r=T(0); IT(N) r+=(*this)[i]*v[i]; return r; }
	T sum() const { T r=T(0); IT(N) r+=(*this)[i]; return r; }
	T mag() const { return sqrt(magSqr()); }
	T magSqr() const { return dot(*this); }
	
	//T norm() const { return sqrt(norm2()); }
	//T norm2() const { return dot(*this); }
	Vec sgn() const { return V(*this).normalize(); }

	Vec& normalize(){
		T msqr = magSqr();
		if(msqr > 0) return (*this) /= sqrt(msqr);
		return Vec().setIdentity();
	}

	template <int N2, class T2>
	Vec& set(const Vec<N2, T2> &v){ IT(N<N2?N:N2){ (*this)[i] = T(v[i]); } return *this; }

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
	
	Vec& setIdentity(){
		(*this)[0] = T(1);
		for(uint32_t i=1; i<N; ++i) (*this)[i] = T(0);
		return *this;
	}

	#undef IT
};

namespace scl{

template<uint32_t N, class T>
inline Vec<N,T> abs(Vec<N,T> a){
	Vec<N,T> r;
	for(uint32_t i=0; i<N; ++i) r[i] = abs(a[i]);
	return r;
}

template<uint32_t N, class T, class U>
inline Vec<N,T> max(Vec<N,T> a, Vec<N,U> b){
	Vec<N,T> r;
	for(uint32_t i=0; i<N; ++i) r[i] = max(a[i], b[i]);
	return r;
}

template<uint32_t N, class T, class U>
inline Vec<N,T> min(Vec<N,T> a, Vec<N,U> b){
	Vec<N,T> r;
	for(uint32_t i=0; i<N; ++i) r[i] = min(a[i], b[i]);
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



template<class T> inline double norm(const Complex<T>& v){ return v.norm(); }
template<class T> inline double normCompare(const Complex<T>& v){ return v.norm2(); }
template<int N,class T> inline double norm(const Vec<N,T>& v){ return v.norm(); }
template<int N,class T> inline double normCompare(const Vec<N,T>& v){ return v.norm2(); }

} // gam::

#endif
