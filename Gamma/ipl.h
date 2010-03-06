#ifndef GAMMA_IPL_H_INC
#define GAMMA_IPL_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information

	File Description:
	Interpolation functions.
*/

#include "Gamma/Constants.h"

namespace gam{

/// Interpolation functions.

/// The naming convention for values is that their alphabetical order is
/// equivalent to their sequence order (oldest to newest).
/// The interpolated value lies 'frac' distance between 'x' and 'y'
///	where 'frac' lies in [0, 1).
namespace ipl{

/// First order allpass interpolation.

/// Best for delay lines, not good for random access.
///
template <class Tf, class Tv>
Tv allpass(Tf frac, const Tv& x, const Tv& y, Tv& o1);

/// First order allpass interpolation with warped fraction.
template <class Tf, class Tv>
Tv allpassFixed(Tf frac, const Tv& x, const Tv& y, Tv& o1);

/// Bezier curve, 3-point quadratic.

/// 'frac' [0, 1) is the value on the curve btw x2 and x0
///
template <class T> T bezier(T frac, T x2, T x1, T x0);

/// Bezier curve, 4-point cubic.

/// 'frac' [0, 1) is the value on the curve btw x3 and x0
///
template <class T> T bezier(T frac, T x3, T x2, T x1, T x0);

template <class Tp, class Tv>
Tv hermite(Tp f, const Tv& w, const Tv& x, const Tv& y, const Tv& z, Tp tension, Tp bias);


/// Computes FIR coefficients for Waring-Lagrange interpolation.

///		'h' are the FIR coefficients and should be of size ('order' + 1). \n
///		'delay' is a fractional delay in samples. \n
///		As order increases, this converges to sinc interpolation.
template <class T> void lagrange(T * h, T delay, int order);

/// Optimized lagrange() for first order.
template <class T> void lagrange1(T * h, T delay);

/// Optimized lagrange() for second order.
template <class T> void lagrange2(T * h, T delay);

/// Optimized lagrange() for third order.
template <class T> void lagrange3(T * h, T delay);

/// Simplified parabolic interpolation of 3 points.

/// This assumes the points are spaced evenly on the x axis from [-1, 1].
/// The output is an offset from 0.
template <class T> T parabolic(T xm1, T x, T xp1);

// Various functions to perform Waring-Lagrange interpolation.
//		These are much faster than using a general purpose FIR filter since
//		the coefs are computed directly and nested multiplication is used
//		rather than directly evaluating the polynomial (FIR).

/// Cubic interpolation

///	This is also known as a Catmull-Rom spline or Cardinal spline with a=-0.5.
///
template <class Tf, class Tv>
Tv cubic(Tf frac, const Tv& w, const Tv& x, const Tv& y, const Tv& z);

template <class Tf, class Tv>
Tv cubic2(Tf d, const Tv& w, const Tv& x, const Tv& y, const Tv& z);

/// Linear interpolation.  Identical to first order Lagrange.
template <class Tf, class Tv>
Tv linear(Tf frac, const Tv& x, const Tv& y);

/// Linear interpolation between three elements.
template <class Tf, class Tv>
Tv linear(Tf frac, const Tv& x, const Tv& y, const Tv& z);

template <class T> void linear(T * dst, const T * xs, const T * xp1s, uint32_t len, T frac);

/// Nearest neighbor interpolation.
template <class Tf, class Tv>
Tv nearest(Tf frac, const Tv& x, const Tv& y);

/// Quadratic interpolation
template <class Tf, class Tv>
Tv quadratic(Tf frac, const Tv& x, const Tv& y, const Tv& z); 


/// Trilinear interpolation between eight corners of a cube.
template <class Tf3, class Tv>
inline Tv trilinear(
	const Tf3& f, 
	const Tv& v000, const Tv& v100,
	const Tv& v010, const Tv& v110,
	const Tv& v001, const Tv& v101,
	const Tv& v011, const Tv& v111
);




// Implementation_______________________________________________________________

template <class Tf, class Tv>
inline Tv allpass(Tf f, const Tv& x, const Tv& y, Tv& o1){
	//f = f * 0.87 - 0.05;	// avoid pole near z = -1
	return o1 = (y - o1) * f + x;
}


template <class Tf, class Tv>
inline Tv allpassFixed(Tf f, const Tv& x, const Tv& y, Tv& o1){
//	f = 1.f-f;
//	f = allpassCoef(f);						// compute allpass coefficient
//	return allpass(f, x, xp1, ym1);	// apply filter
	
	//f = f / (2.f - f);	// warp down
	f = (Tf(2) * f) / (Tf(1) + f);	// warp up
	//f = (1.5f * f) / (0.5f + f);
	f -= Tf(0.1);				// avoid pole near z = -1
	return allpass(f, x, y, o1);	// apply filter
}
//
//inline float Ipol::allpassCoef(float f, float offset){
//	return (1.f - f) / (1.f + f) + offset;
//}

template <class T> inline T bezier(T d, T x2, T x1, T x0){
	T d2 = d * d;
	T dm1 = T(1) - d;
	T dm12 = dm1 * dm1;
	return x2 * dm12 + T(2) * x1 * dm1 * d + x0 * d2;

//	x2 (1-d)(1-d) + 2 x1 (1-d) d + x0 d d
//	x2 - d 2 x2 + d d x2 + d 2 x1 - d d 2 x1 + d d x0
//	x2 - d (2 x2 + d x2 + 2 x1 - d 2 x1 + d x0)
//	x2 - d (2 (x2 + x1) + d (x2 - 2 x1 + x0))

//	float c2 = x2 - 2.f * x1 + x0;
//	float c1 = 2.f * (x2 + x1);
//	return x2 - (d * c2 + c1) * d;
}


template <class T> inline T bezier(T d, T x3, T x2, T x1, T x0){
	T c1 = T(3) * (x2 - x3);
	T c2 = T(3) * (x1 - x2) - c1;
	T c3 = x0 - x3 - c1 - c2;
	return ((c3 * d + c2) * d + c1) * d + x3;
}


template <class Tf, class Tv>
inline Tv cubic(Tf f, const Tv& w, const Tv& x, const Tv& y, const Tv& z){	
//	Tv c3 = (x - y)*(Tf)1.5 + (z - w)*(Tf)0.5;
//	Tv c2 = w - x*(Tf)2.5 + y*(Tf)2. - z*(Tf)0.5;
//	Tv c1 = (y - w)*(Tf)0.5;
//	return ((c3 * f + c2) * f + c1) * f + x;

	// -w + 3x - 3y + z	
	// 2w - 5x + 4y - z
	// c2 = w - 2x + y - c3

//	Tv c3 = (x - y)*(Tf)3 + z - w;
//	Tv c2 = w - x*(Tf)2 + y - c3;
//	Tv c1 = y - w;
//	return (((c3 * f + c2) * f + c1)) * f * (Tf)0.5 + x;
	
//	Tv c3 = (x - y)*(Tf)1.5 + (z - w)*(Tf)0.5;
//	Tv c2 = (y + w)*(Tf)0.5 - x - c3;
//	Tv c1 = (y - w)*(Tf)0.5;	
//	return ((c3 * f + c2) * f + c1) * f + x;

	Tv c1 = (y - w)*Tf(0.5);
	Tv c3 = (x - y)*Tf(1.5) + (z - w)*Tf(0.5);
	Tv c2 = c1 + w - x - c3;
	return ((c3 * f + c2) * f + c1) * f + x;
}

template <class T>
void cubic(T * dst, const T * xm1s, const T * xs, const T * xp1s, const T * xp2s, uint32_t len, T f){
	for(uint32_t i=0; i<len; ++i) dst[i] = cubic(f, xm1s[i], xs[i], xp1s[i], xp2s[i]);
}


// From http://astronomy.swin.edu.au/~pbourke/other/interpolation/ (Paul Bourke)
template <class Tf, class Tv>
inline Tv cubic2(Tf f, const Tv& w, const Tv& x, const Tv& y, const Tv& z){
	Tv c3 = z - y - w + x;
	Tv c2 = w - x - c3;
	Tv c1 = y - w;
	return ((c3 * f + c2) * f + c1) * f + x;
}


// From http://astronomy.swin.edu.au/~pbourke/other/interpolation/ (Paul Bourke)
/*
   Tension: 1 is high, 0 normal, -1 is low
   Bias: 0 is even,
         positive is towards first segment,
         negative towards the other
*/
template <class Tp, class Tv>
inline Tv hermite(Tp f,
	const Tv& w, const Tv& x, const Tv& y, const Tv& z,
	Tp tension, Tp bias)
{
	tension = (Tp(1) - tension)*Tp(0.5);

	// compute endpoint tangents
	//Tv m0 = ((x-w)*(1+bias) + (y-x)*(1-bias))*tension;
	//Tv m1 = ((y-x)*(1+bias) + (z-y)*(1-bias))*tension;
	Tv m0 = ((x*Tv(2) - w - y)*bias + y - w)*tension;
	Tv m1 = ((y*Tv(2) - x - z)*bias + z - x)*tension;
	
//	x - w + x b - w b + y - x - y b + x b
//	-w + 2x b - w b + y - y b
//	b(2x - w - y) + y - w			
//	
//	y - x + y b - x b + z - y - z b + y b
//	-x + 2y b - x b + z - z b
//	b(2y - x - z) + z - x

	Tp f2 = f  * f;
	Tp f3 = f2 * f;

	// compute hermite basis functions
	Tp a3 = Tp(-2)*f3 + Tp(3)*f2;
	Tp a0 = Tp(1) - a3;
	Tp a2 = f3 - f2;
	Tp a1 = f3 - Tp(2)*f2 + f;

	return x*a0 + m0*a1 + m1*a2 + y*a3;
}


template <class T> void lagrange(T * a, T delay, int order){
	for(uint32_t i=0; i<=order; i++){
		T coef = T(1);
		T i_f = T(i); 
		for(uint32_t j=0; j<=order; j++){
			if(j != i){
				T j_f = (T)j;
				coef *= (delay - j_f) / (i_f - j_f);
			}
		}
		*a++ = coef;
	}
}


template <class T> inline void lagrange1(T * h, T d){
	h[0] = T(1) - d;
	h[1] = d;
}


template <class T> inline void lagrange2(T * h, T d){
	h[0] =      (d - T(1)) * (d - T(2)) * T(0.5);
	h[1] = -d              * (d - T(2))         ;
	h[2] =  d * (d - T(1))              * T(0.5);
}


template <class T> inline void lagrange3(T * h, T d){
	T d1 = d - T(1);
	T d2 = d - T(2);
	T d3 = d - T(3);
	h[0] =     -d1 * d2 * d3 * T(1./6.);
	h[1] =  d      * d2 * d3 * T(0.5);
	h[2] = -d * d1      * d3 * T(0.5);
	h[3] =  d * d1 * d2      * T(1./6.);
}

/*
x1 (1 - d) + x0 d
x1 - x1 d + x0 d
x1 + (x0 - x1) d

x2 (d - 1) (d - 2) /2 - x1 d (d - 2) + x0 d (d - 1) /2
d d /2 x2 - d 3/2 x2 + x2 - d d x1 + d 2 x1 + d d /2 x0 - d /2 x0
d d /2 x2 - d d x1 + d d /2 x0 - d 3/2 x2 + d 2 x1 - d /2 x0 + x2
d (d (/2 x2 - x1 + /2 x0) - 3/2 x2 + 2 x1 - /2 x0) + x2
*/

template <class Tf, class Tv>
inline Tv linear(Tf f, const Tv& x, const Tv& y){
	return (y - x) * f + x;
}

template <class Tf, class Tv>
inline Tv linear(Tf frac, const Tv& x, const Tv& y, const Tv& z){
	frac *= Tf(2);
	if(frac<Tf(1)) return ipl::linear(frac, x,y);
	return ipl::linear(frac-Tf(1), y,z);
}

template <class T>
void linear(T * dst, const T * xs, const T * xp1s, uint32_t len, T f){
	for(uint32_t i=0; i<len; ++i) dst[i] = linear(f, xs[i], xp1s[i]);
}


template <class Tf, class Tv>
inline Tv nearest(Tf f, const Tv& x, const Tv& y){
	return (f < Tf(0.5)) ? x : y;
}


template <class T> inline T parabolic(T xm1, T x, T xp1){
	T numer = xm1 - xp1;
	T denom = x - xp1 + x - xm1;
	return T(-0.5) * numer / denom;
}


template <class Tf, class Tv>
inline Tv quadratic(Tf f, const Tv& x, const Tv& y, const Tv& z){
	Tv c2 = (x + z)*Tf(0.5) - y;
	//Tv c1 = x*(Tf)-1.5 + y*(Tf)2 - z*(Tf)0.5;
	Tv c1 = -x + y - c2;
	return (c2 * f + c1) * f + x;
}


template <class F3, class Tv>
inline Tv trilinear(
	const F3& f, 
	const Tv& v000, const Tv& v100,
	const Tv& v010, const Tv& v110,
	const Tv& v001, const Tv& v101,
	const Tv& v011, const Tv& v111
){
	float f1=f[0], f2=f[1];

	return
	linear(f[2],
		linear(f2,
			linear(f1, v000, v100),
			linear(f1, v010, v110)
		),
		linear(f2,
			linear(f1, v001, v101),
			linear(f1, v011, v111)			
		)
	);
}

} // ipl::
} // gam::
#endif
