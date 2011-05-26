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
template<class T> class Vec2;
template<class T> class Vec3;
template<class T> class Vec4;


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




/// Fixed size vector.
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



// Two element vector
template <class T>
struct Vec2 : public Vec<2, T> {
	using Vec<2,T>::operator=;
	Vec2(const Vec<2, T>& v){ (*this)=v; }
	Vec2(const T& v=T(0)){ (*this)(v,v); }
	Vec2(const T& v1, const T& v2){ (*this)(v1,v2); }
	Vec2& operator()(const T& v1, const T& v2){ (*this)[0]=v1; (*this)[1]=v2; return *this; }
};


// Three element vector
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


// Four element vector
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

} // gam::

#endif
