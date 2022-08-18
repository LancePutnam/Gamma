#ifndef GAMMA_IPL_H_INC
#define GAMMA_IPL_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information

	File Description:
	Functions for interpolating sample data. Given a set of sample values,
	interpolation estimates values between the samples using some continuous
	function. It is a form of upsampling.
*/

namespace gam{

/// Interpolation functions.

/// The naming convention for values is that their alphabetical order is
/// equivalent to their sequence order (oldest to newest).
/// The interpolated value lies 'frac' distance between 'x' and 'y'
///	where 'frac' lies in [0, 1).
namespace ipl{


/// Interpolation types
enum Type{
	TRUNC=0,	/**< Truncating interpolation */
	ROUND,		/**< Nearest neighbor interpolation */
	MEAN2,		/**< Mean of two nearest neighbors */
	LINEAR,		/**< Linear interpolation */
	CUBIC,		/**< Cubic interpolation */
	ALLPASS		/**< Allpass interpolation */
};


/// Nearest neighbor interpolation
template <class Tf, class Tv>
Tv nearest(Tf frac, const Tv& x, const Tv& y){
	return (frac < Tf(0.5)) ? x : y;
}

/// Linear interpolation (first order Lagrange)
template <class Tf, class Tv>
Tv linear(Tf frac, const Tv& x, const Tv& y){
	return (y - x) * frac + x;
}

/// Linear interpolation between three elements
template <class Tf, class Tv>
Tv linear(Tf frac, const Tv& x, const Tv& y, const Tv& z){
	frac *= Tf(2);
	if(frac<Tf(1)) return ipl::linear(frac, x,y);
	return ipl::linear(frac-Tf(1), y,z);
}

/// Linear interpolation (lookup) within array
template <class Tf, class Tv>
Tv linear(Tf pos01, const Tv * src, unsigned len){
	if(pos01 >= Tf(1)) return src[len-1];
	Tf idxf = pos01 * (len-1);
	unsigned i0 = unsigned(idxf);
	unsigned i1 = i0+1; if(i1>=len) i1=len-1;
	float f = idxf - i0;
	return linear(f, src[i0], src[i1]);
}

template <class Tf, class Array>
typename Array::value_type linear(Tf pos01, const Array& src){
	return linear(pos01, &src[0], src.size());
}

/// Linear interpolation between arrays
template <class T>
void linear(T * dst, const T * xs, const T * xp1s, unsigned len, T frac){
	for(unsigned i=0; i<len; ++i) dst[i] = linear(frac, xs[i], xp1s[i]);
}

/// Cyclic linear interpolation between three elements
template <class Tf, class Tv>
Tv linearCyclic(Tf frac, const Tv& x, const Tv& y, const Tv& z){
	frac *= Tf(3);
	if(frac <= Tf(1))		return ipl::linear(frac, x,y);
	else if(frac >= Tf(2))	return ipl::linear(frac-Tf(2), z,x);
							return ipl::linear(frac-Tf(1), y,z);
}

/// Cyclic linear interpolation (lookup) within array
template <class Tf, class Tv>
Tv linearCyclic(Tf pos01, const Tv * src, unsigned len){
	Tf idxf = pos01 * len;
	unsigned i0 = unsigned(idxf);
	unsigned i1 = i0+1; if(i1>=len) i1=0;
	float f = idxf - i0;
	return linear(f, src[i0], src[i1]);
}

template <class Tf, class Array>
typename Array::value_type linearCyclic(Tf pos01, const Array& src){
	return linearCyclic(pos01, &src[0], src.size());
}

/// Trapezoidal interpolation
template <class Tf, class Tv>
Tv trapz(Tf frac, const Tv& x, const Tv& y){
	return trapz2(frac*Tf(2), x,y);
}

/// Trapezoidal interpolation; fraction in [0,2] to save a multiply
template <class Tf, class Tv>
Tv trapz2(Tf frac02, const Tv& x, const Tv& y){
	if(frac02 < Tf(1))	return x + y*frac02;
	else				return x*(Tf(2)-frac02) + y;
}

/// Quadratic interpolation
template <class Tf, class Tv>
Tv quadratic(Tf frac, const Tv& x, const Tv& y, const Tv& z){
	Tv c2 = (x + z)*Tf(0.5) - y;
	//Tv c1 = x*(Tf)-1.5 + y*(Tf)2 - z*(Tf)0.5;
	Tv c1 = -x + y - c2;
	return (c2 * frac + c1) * frac + x;
}

/// Cubic interpolation

///	This is a Cardinal spline with a tension of 0 (AKA a Catmull-Rom spline).
/// The resulting value will never exceed the range of the interpolation points.
template <class Tf, class Tv>
Tv cubic(Tf frac, const Tv& w, const Tv& x, const Tv& y, const Tv& z){
	/*
	Tv c3 = (x - y)*Tf(1.5) + (z - w)*Tf(0.5);
	Tv c2 = w - x*Tf(2.5) + y*Tf(2.) - z*Tf(0.5);
	Tv c1 = (y - w)*Tf(0.5);
	return ((c3 * frac + c2) * frac + c1) * frac + x;//*/

	// -w + 3x - 3y + z	
	// 2w - 5x + 4y - z
	// c2 = w - 2x + y - c3

	/*
	Tv c3 = (x - y)*Tf(3) + z - w;
	Tv c2 = w - x*Tf(2) + y - c3;
	Tv c1 = y - w;
	return (((c3 * frac + c2) * frac + c1)) * frac * Tf(0.5) + x;//*/

	/*
	Tv c3 = (x - y)*Tf(1.5) + (z - w)*Tf(0.5);
	Tv c2 = (y + w)*Tf(0.5) - x - c3;
	Tv c1 = (y - w)*Tf(0.5);	
	return ((c3 * frac + c2) * frac + c1) * frac + x; //*/

	//*
	Tv c1 = (y - w)*Tf(0.5);
	Tv c3 = (x - y)*Tf(1.5) + (z - w)*Tf(0.5);
	Tv c2 = c1 + w - x - c3;
	return ((c3 * frac + c2) * frac + c1) * frac + x;//*/
}

/// Cubic interpolation between arrays
template <class T>
void cubic(T * dst, const T * xm1s, const T * xs, const T * xp1s, const T * xp2s, unsigned len, T frac){
	for(unsigned i=0; i<len; ++i) dst[i] = cubic(frac, xm1s[i], xs[i], xp1s[i], xp2s[i]);
}

/// Cyclic cubic interpolation (lookup) within array
template <class Tf, class Tv>
Tv cubicCyclic(Tf pos01, const Tv * src, unsigned len){
	Tf idxf = pos01 * len;
	unsigned i0 = unsigned(idxf);
	unsigned im1= i0>0 ? i0-1 : len-1;
	unsigned i1 = i0+1; if(i1>=len) i1-=len;
	unsigned i2 = i0+2; if(i2>=len) i2-=len;
	float f = idxf - i0;
	return cubic(f, src[im1], src[i0], src[i1], src[i2]);
}

template <class Tf, class Array>
typename Array::value_type cubicCyclic(Tf pos01, const Array& src){
	return cubicCyclic(pos01, &src[0], src.size());
}

/// Simplified parabolic interpolation of 3 points.

/// This assumes the points are spaced evenly on the x axis from [-1, 1].
/// The output is an offset from 0.
template <class T>
T parabolic(const T& xm1, const T& x, const T& xp1){
	T numer = xm1 - xp1;
	T denom = x - xp1 + x - xm1;
	return T(-0.5) * numer / denom;
}

/// First order allpass interpolation

/// Allpass interpolation preserves the spectral magnitude of the interpolated
/// sequence at the expense of phase delay. It should only be used for
/// interpolating at a regular speed across a sequence. It will fail miserably
/// at random access.
template <class Tf, class Tv>
Tv allpass(Tf frac, const Tv& w, const Tv& x, const Tv& y, Tv& o1){
	// y[n]	= a x[n] + x[n-1] - a y[n-1]
	//		= a (x[n] - y[n-1]) + x[n-1]

	// The delay is shifted to minimize transients when changing the delay time.
	// This avoids having a pole on the unit circle. See:
	// Valimaki, V., Laakso, T. I., and Mackenzie, J. (1995). Elimination of transients in time-varying allpass fractional delay filters with application to digital waveguide modeling. In Proceedings of the International Computer Music Conference, pages 327–334.
	const Tf shift = 0.5;

	frac = Tf(1)+shift-frac;	// convert to sample delay

	auto i1 = x;
	auto i0 = y;

	if(frac >= Tf(1)){ // delay must be in [0,1], so jump back 1 sample of input
		--frac;
		i1=w;
		i0=x;
	}

	//Tf a = (Tf(1)-frac)/(Tf(1)+frac);

	// Approx above to avoid division
	//Tf a = Tf(1)-frac; // apx
	//Tf a = Tf(0.5)*(Tf(1)-frac)*(Tf(2)-frac); // apx
	Tf a = (Tf(1)-frac); a = Tf(0.4545454545)*a*(a*a + Tf(1.2));

	return o1 = (i0 - o1) * a + i1;
}

/// Computes FIR coefficients for Waring-Lagrange interpolation.

/// This function is useful if the interpolation fraction is constant.
///
/// \param[in] h		FIR coefficients; must be of size 'order' + 1.
///	\param[in] delay	Fractional delay in samples
///	\param[in] order	As order increases, this converges to sinc interpolation
template <class T> void lagrange(T * h, T delay, unsigned order){
	for(unsigned i=0; i<=order; ++i){
		T coef = T(1);
		T i_f = T(i); 
		for(unsigned j=0; j<=order; ++j){
			if(j != i){
				T j_f = (T)j;
				coef *= (delay - j_f) / (i_f - j_f);
			}
		}
		*h++ = coef;
	}
}

/// Optimized lagrange() for first order.
template <class T> void lagrange1(T * h, T delay){
	auto& d = delay;
	h[0] = T(1) - d;
	h[1] = d;
}

/// Optimized lagrange() for second order.
template <class T> void lagrange2(T * h, T delay){
	auto& d = delay;
	h[0] =      (d - T(1)) * (d - T(2)) * T(0.5);
	h[1] = -d              * (d - T(2))         ;
	h[2] =  d * (d - T(1))              * T(0.5);
}

/// Optimized lagrange() for third order.
template <class T> void lagrange3(T * h, T delay){
	auto& d = delay;
	T d1 = d - T(1);
	T d2 = d - T(2);
	T d3 = d - T(3);
	h[0] =     -d1 * d2 * d3 * T(1./6.);
	h[1] =  d      * d2 * d3 * T(0.5);
	h[2] = -d * d1      * d3 * T(0.5);
	h[3] =  d * d1 * d2      * T(1./6.);
}

} // ipl::
} // gam::
#endif
