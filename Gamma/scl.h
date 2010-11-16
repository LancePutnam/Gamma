#ifndef GAMMA_SCL_H_INC
#define GAMMA_SCL_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information

	File Description:
	Mathematical scalar operations.
*/

#include <math.h>
#include <stdlib.h>				/* labs(long) */
#include "Gamma/Conversion.h"
#include "Gamma/Types.h"

#define TEM template<class T>

// undefine macros in windows.h
#ifdef max
#undef max
#endif
#ifdef min
#undef min
#endif
#ifdef sinc
#undef sinc
#endif

namespace gam{

/// Returns a positive length associated with argument
TEM double norm(const T& v);

/// Returns a positive number valid for comparing norms
TEM double normCompare(const T& v);


/// Scalar rank functions for numerical types.
namespace scl{


template<int N, class T, template<class> class F> struct NewtonIterator{
	NewtonIterator(T& v, T v0){
		F<T>(v,v0);
		NewtonIterator<N-1,T,F>(v,v0); // this just iterates
	}
};
template<class T, template<class> class F> struct NewtonIterator<0,T,F>{ NewtonIterator(T& v, T v0){} };

template<class T> struct NewtonSqrtMap{ NewtonSqrtMap(T& v, T v0){ v=T(0.5)*(v+v0/v); } };
template<int N, class T> struct SqrtNewton
:	public NewtonIterator<N,T, NewtonSqrtMap>{
	SqrtNewton(T& v, T v0): NewtonIterator<N,T, NewtonSqrtMap>(v,v0){}
};


template<class T> struct NewtonInvSqrtMap{ NewtonInvSqrtMap(T& v, T v0_2){ v *= T(1.5)-v0_2*v*v; } };
template<int N, class T> struct InvSqrtNewton
:	public NewtonIterator<N,T, NewtonInvSqrtMap>{
	InvSqrtNewton(T& v, T v0): NewtonIterator<N,T, NewtonInvSqrtMap>(v,v0){}
};

template<class T> const Twiddle<T> invSqrtMagic();
template<> inline const Twiddle<float > invSqrtMagic(){ return Twiddle<float >(0x5f3759df); }
template<> inline const Twiddle<double> invSqrtMagic(){ return Twiddle<double>(UINT64_C(0x5fe6ec85e7de30da)); }


/// Approximate square root using a quick log base-2 method.
inline float sqrtLog2(float v){
	Twiddle<float> u(v);
	u.u=(1<<29) + (u.u>>1) - (1<<22);
	return u.f;
}

/// Approximate square root using a quick log base-2 method.
inline double sqrtLog2(double v){
	Twiddle<double> u(v);
	u.u=((uint64_t(1))<<61) + (u.u>>1) - ((uint64_t(1))<<51);
	return u.f;
}

/// Approximate square root using Newton's method.
template<int N, class T> void sqrtNewton(T& v, T v0){ SqrtNewton<N,T>(v,v0); }

/// Approximate square root using log base-2 and Newton methods.

/// 'N' determines the accuracy of the approximation. For N=0, a quick and dirty
/// log base-2 approximation is performed. For N>0, N-1 Newton iterations
/// are applied to improve the result.
template<uint32_t N, class T>
inline T sqrt(T v){
	T r=sqrtLog2(v);
	sqrtNewton<N>(r,v);
	return r;
}

/// Approximate inverse square root using Newton's method.
template<int N, class T> void invSqrtNewton(T& v, T v0){ InvSqrtNewton<N,T>(v,v0); }

/// Approximate inverse square root using a quick log base-2 method.
template <class T>
inline T invSqrtLog2(T v){
	Twiddle<T> u(v);
	u.u = invSqrtMagic<T>().u - (u.u>>1);
	return u.f;
}

/// Approximate inverse square root using log base-2 and Newton methods.

/// 'N' determines the accuracy of the approximation. For N=0, a quick and dirty
/// log base-2 approximation is performed. For N>0, N-1 Newton iterations
/// are applied to improve the result.
template<uint32_t N, class T>
inline T invSqrt(T v){
	T r=invSqrtLog2(v);
	invSqrtNewton<N>(r, v*T(0.5));
	return r;
}



/// Returns absolute value.
TEM T abs(T value);

/// Tests whether two values are close in value

/// 'maxULP' (maximum units in the last place) is the maximum difference allowed
/// between the values punned as integers.
/// Code taken from Bruce Dawson, "Comparing floating point numbers."
bool almostEqual(float a, float b, int maxULP=10);
bool almostEqual(double a, double b, int maxUlps=10);

/// Fast approximation to atan2().

// Author: Jim Shima, http://www.dspguru.com/comp.dsp/tricks/alg/fxdatan2.htm.
// |error| < 0.01 rad
TEM T atan2Fast(T y, T x);

/// Returns value clipped ouside of range [-eps, eps]
TEM T atLeast(T v, T eps);

/// Returns bump function, psi(x)
TEM T bump(T x);

/// Returns floating point value rounded to next highest integer.
TEM T ceil(T val);
TEM T ceil(T val, T step);
TEM T ceil(T val, T step, T recStep);

/// Returns power of two ceiling of value.

/// This uses an algorithm devised by Sean Anderson, Sep. 2001.
/// From "Bit Twiddling Hacks", http://graphics.stanford.edu/~seander/bithacks.html.
uint32_t ceilPow2(uint32_t value);

/// Returns even number ceiling.
uint32_t ceilEven(uint32_t value);

/// Returns value clipped to [lo, hi].
TEM T clip(T value, T hi=T(1), T lo=T(0));

/// Returns value clipped to [-hi, hi].
TEM T clipS(T value, T hi=T(1));

/// Returns value clipped to [lo, hi] and signifies clipping behavior.

/// clipFlag signifies if and where clipping occured.  0 means no clipping
/// occured, -1 means clipping occured at the lower bound, and 1 means
/// clipping at the upper bound.
TEM T clip(T value, int & clipFlag, T hi, T lo);

/// Returns value whose magnitude is clipped to [min, max].
float clipMag(float value, float min, float max);

/// Third order polynomial approximation to first half of cosine.

/// 'u' must be in the range [0, 0.5] which corresponds to the first half
/// of the cosine.
TEM T cosP3(T u);

/// 8th order Taylor series approximation to a cosine.

/// 'radians' must be in [-pi, pi].
///
TEM T cosT8(T radians);

/// Compute curvature around point b of three successive points a, b, and c.
template <class T, template <class> class V>
T curvature(const V<T>& a, const V<T>& b, const V<T>& c);

/// Convert decibels to amplitude
template <class T>
inline T dBToAmp(const T& db){ return ::pow(10, db/20.); }

/// Convert amplitude to decibels
template <class T>
inline T ampTodB(const T& amp){ return 20*::log(amp); }

/// Returns two element dot product x1 * y1 + x2 * y2.
TEM T dot2(T x1, T x2, T y1, T y2);

/// Returns relative error between a true and measured value.
TEM inline double error(const T& truth, const T& measured){ return (measured-truth)/double(truth);}

/// Returns weights for linear fade.
TEM void fadeLin(T& weight1, T& weight2, T fade);

/// Returns weights for triangular window fade.

/// The weights returned are from two overlapping normalized unipolar 
/// triangular windows in the fade interval [0, 2/3] and [1/3, 1].\n
/// fade	weight1		weight2 \n
/// 0.25	1			0       \n
/// 0.5		0.5			0.5		\n
/// 0.75	0			1		
TEM void fadeTri(T& weight1, T& weight2, T fade);

uint32_t factorial12(uint32_t value);				///< Returns factorial of value in [0, 12].

/// Return Fejer weighting factor for kth harmonic in a Fourier series of length n.

/// The function is a line from (0,1) to (n,0). This is used for reducing the
/// Gibbs effect from a Fourier series summation.
TEM T fejer(T k, T n){ return (n-k)/n; }

/// Returns floor of floating point value.
TEM T floor(T val);
TEM T floor(T val, T step);
TEM T floor(T val, T step, T recStep);

/// Returns power of two floor of value.

/// This uses an algorithm devised by Sean Anderson, Sep. 2001.
/// From "Bit Twiddling Hacks", http://graphics.stanford.edu/~seander/bithacks.html.
uint32_t floorPow2(uint32_t value);

/// Returns value folded into [lo, hi].

/// For out-of-range values, the boundaries act like mirrors reflecting
/// the value into the range. For an even number of periods out of the range
/// this is identical to a wrap().
TEM T fold(T value, T hi=T(1), T lo=T(0));
TEM T fold(T v, long& numFolds, T hi=T(1), T lo=T(0));

/// Returns value folded into [lo, hi] one time.
TEM T foldOnce(T value, T hi=T(1), T lo=T(0));

/// Returns frequency in Hz from a 12-TET note string.

/// Notes are specified by a letter in [a, g], followed optionally by one of 
/// '+', '-', ' ', to specify a sharp, flat or natural, and finally an integer
/// in [0,9] representing the octave number.  For example, middle C is specified
/// as "c5" or "c 5" and the A sharp below that as "a+4".
double freq(const char * note);

/// Convert domain frequency to radians clipped to interval [0, pi).
TEM inline T freqToRad(T freq, double ups){ return scl::clip(freq * ups, 0.499) * M_PI; }

/// Compute Frenet frame (tangent, normal) from 1st difference
template <class V2>
void frenet(const V2& d1, V2& t, V2& n);


/// Compute Frenet frame (tangent, normal, binormal) from 1st and 2nd differences
template <class V3>
void frenet(const V3& d1, const V3& d2, V3& t, V3& n, V3& b);

/// Compute Frenet frame (tangent, normal, binormal) from 3 consecutive points
template <class V3>
void frenet(const V3& p2, const V3& p1, const V3& p0, V3& t, V3& n, V3& b);


/// Returns e^(-v*v)
TEM T gaussian(T v){ return ::exp(-v*v); }

/// Return greatest common divisor of two arguments
TEM T gcd(const T& x, const T& y);

TEM T hypot(T x, T y);

/// Generalized Laguerre polynomial L{n,k}
///
/// http://en.wikipedia.org/wiki/Laguerre_polynomials
double laguerre(int n, int k, double x);

/// Returns least common multiple
TEM inline T lcm(const T& x, const T& y){ return (x*y)/gcd(x,y); }

/// Associated Legendre polynomial
///
/// P_l^m(cos(t)) = (-1)^{l+m} / (2^l l!) sin^m(t) (d/d cos(t))^{l+m} sin^{2l}(t)
///
/// @param[in]	l	degree s.t. l>=0
/// @param[in]	m	order s.t. -l<=m <=l
/// @param[in]	t	variable
///
/// http://comp.cs.ehime-u.ac.jp/~ogata/nac/index.html
double legendre(int l, int m, double t);
double legendre(int l, int m, double ct, double st);

/// Convert linear value to log2 in range [0, 1]
TEM T linLog2(T v, T recMin);

/// Returns base 2 logarithm of value.

/// If the value is not an exact power of two, the logarithm of the next
/// highest power of two will taken.
/// This uses an algorithm devised by Eric Cole, Jan. 2006.
/// From "Bit Twiddling Hacks", http://graphics.stanford.edu/~seander/bithacks.html.
uint32_t log2(uint32_t v);

/// Fast base 2 logarithm.  For value <= 0, behavior is undefined.
float log2Fast(float v);

/// Maps value from [-1,1] to [depth, 1].
TEM T mapDepth(T v, T depth){ return (v - T(1)) * T(0.5) * depth + T(1);  }

/// Inverse 2nd power mapping for interval [0, 1].
TEM T mapInvPow2(T v);

/// Computes scale and offset necessary to map value from [i0, i1] to [o0, o1].
TEM void mapLin(T i0, T i1, T o0, T o1, T & scale, T & offset);

/// Linearly maps value from [i0, i1] to [o0, o1].
TEM T mapLin(T v, T i0, T i1, T o0, T o1);

/// Returns v^power linearly mapped to [bound0, bound1].
double mapPower(double v, double bound1, double bound0, double power=1.);

/// Mixes two values together (1 = thru, 0.5 = mono, 0 = swap).
TEM void mix2(T& io1, T& io2, T mix);

/// Perform complex multiplication, c1 = c1 c2.
TEM void mulComplex(T& r1, T& i1, const T& r2, const T& i2);

TEM T nearest(T v, const char * interval="2212221", long div=12);

/// Returns nearest integer division of one value to another
TEM T nearestDiv(T of, T to);

/// Smooth negative map.

/// The return value is close to 1 if v < 0 and close to 0 if v > 0.
/// The smoothness is controlled with the bw argument.
TEM T negative(T v, T bw);

/// Same as negative(T,T), but 'amt' controls positive level (0,1) -> (-1,1)
TEM T negative(T v, T bw, T amt);

/// Returns the next representable floating-point or integer value following x in the direction of y
TEM T nextAfter(T x, T y);

/// Returns the number of digits in the integer portion
TEM T numInt(const T& v){ return scl::floor(::log10(v)) + 1; }

/// Returns pole radius given a bandwidth and sampling interval
TEM	inline T poleRadius(T bw, double ups){ return ::exp(-M_PI * bw * ups); }
//return (T)1 - (M_2PI * bw * ups); // linear apx for fn < ~0.02

/// Evaluates polynomial a0 + a1 x + a2 x^2
TEM T poly(T x, T a0, T a1, T a2);

/// Evaluates polynomial a0 + a1 x + a2 x^2 + a3 x^3
TEM T poly(T x, T a0, T a1, T a2, T a3);

/// Smooth positive map

/// The return value is close to 1 if v > 0 and close to 0 if v < 0.
/// The smoothness is controlled with the bw argument.
TEM T positive(T v, T bw);

TEM T pow2(T v);			///< Returns value to the 2nd power.
TEM T pow2S(T v);			///< Returns value to the 2nd power preserving sign.
TEM T pow3(T v);			///< Returns value to the 3rd power.
TEM T pow3Abs(T v);			///< Returns absolute value to the 3rd power.
TEM T pow4(T v);			///< Returns value to the 4th power.
TEM T pow5(T v);			///< Returns value to the 5th power.
TEM T pow6(T v);			///< Returns value to the 6th power.
TEM T pow8(T v);			///< Returns value to the 8th power.
TEM T pow16(T v);			///< Returns value to the 16th power.
TEM T pow64(T v);			///< Returns value to the 64th power.

unsigned char prime(uint32_t n);	///< Returns (n+1)th prime number up to n=53.
TEM T prime(uint32_t n, T mul);	///< Returns scaled (n+1)th prime number up to n=53.

/// Returns spherical Euler product (ZXZ convention)
TEM Vec3<T> productZXZ(const Complex<T>& a, const Complex<T>& b, const Complex<T>& c);

/// Returns pole radius given a T60 decay length and units/sample
inline double radius60(double dcy, double ups){ return ::exp(M_LN001/dcy * ups); } // u/s * 1/u

/// Returns equal temperament ratio- octave^(pc/divisions)
TEM T ratioET(T pc, T divisions=12, T octave=2);

/// Returns the value r such that r = x - n*y.
TEM T remainder(const T& x, const T& y);

/// Returns floating point value rounded to nearest integer.
TEM T round(T v);

/// Returns floating point value rounded to nearest integer multiple of 'step'.
TEM T round(T v, T step);

/// Returns floating point value rounded to nearest integer multiple of 'step'. Faster version to avoid 1/step divide.
TEM T round(T v, T step, T recStep);

/// Returns value rounded to nearest integer away from zero.
TEM T roundAway(T v);

/// Returns value rounded to nearest to nearest integer multiple of 'step' away from zero.
TEM T roundAway(T v, T step);

/// Returns the section 'v' lies in in [0,num] divided into 'div' sections.
inline int section(int v, int num, int divs){ return (v*divs)/double(num); }

/// Signum function for real numbers.
TEM inline T sgn(const T& v, const T& norm=T(1)){ return v<T(0) ? -norm : norm; }

/// Returns value of spherical harmonic Y{l,m}(theta, phi).
TEM Complex<T> sharm(int l, int m, T theta, T phi);

#define DEF(name) TEM inline Complex<T> name(T ct, T st, T cp, T sp)

// Spherical harmonics of l and m. For m!=0, the +m harmonic is returned.
// The input angles are in terms of cosine and sine of theta and phi.
// The phi arguments should be of an angle m times the base phi angle.
// For odd, negative m, the result should by multiplied by -1.

DEF(sharm00){ static const T c= 0.50*::sqrt(    M_1_PI ); return Complex<T>(c, 0); }
DEF(sharm10){ static const T c= 0.50*::sqrt( 3.*M_1_PI ); return Complex<T>(c*ct, 0); }
DEF(sharm11){ static const T c=-0.50*::sqrt( 3.*M_1_2PI); return Complex<T>(cp,sp)*c*st; }
DEF(sharm20){ static const T c= 0.25*::sqrt( 5.*M_1_PI ); return Complex<T>(c*(3.*ct*ct - 1.), 0); }
DEF(sharm21){ static const T c=-0.50*::sqrt(15.*M_1_2PI); return Complex<T>(cp,sp)*c*ct*st; }
DEF(sharm22){ static const T c= 0.25*::sqrt(15.*M_1_2PI); return Complex<T>(cp,sp)*c*st*st; }

#undef DEF

/// Continuous sign map.

/// The return value is close to 1 if v > 0 and close to -1 if v < 0.
/// 'bw' controls the width of the transition region.
TEM T sign(T v, T bw);

//
TEM T sinFast(T radians);

/// 7th order minimax polynomial approximation to a sine.
TEM T sinP7(T normal);

/// 9th order minimax polynomial approximation to a sine.
TEM T sinP9(T normal);

/// 7th order Taylor series approximation to a sine.

/// 'radians' must be in [-pi, pi].
///
TEM T sinT7(T radians);

/// 9th order Taylor series approximation to a sine.

/// 'radians' must be in [-pi, pi].
///
TEM T sinT9(T radians);

/// Unnormalized sinc function
TEM T sinc(T radians, T eps=T(0.0001));

/// Sort values so that value1 <= value2.
TEM void sort(T& value1, T& value2);

/// Returns spherical product of two complex numbers
TEM Vec3<T> spherical(const Complex<T>& a, const Complex<T>& b){ return Vec3<T>(a.r*b.r, a.i*b.r, b.i); }

/// Sum of integers squared from 1 to n.
TEM T sumOfSquares(T n){
	static const T c1_6 = 1/T(6);
	static const T c2_6 = c1_6*T(2);
	return n*(n+1)*(c2_6*n+c1_6);
}


/// Truncates floating point value at decimal.
TEM T trunc(T value);

/// Truncates floating point value to step.
TEM T trunc(T value, T step);

/// Truncates floating point value to step. Faster version to avoid 1/step divide.
TEM T trunc(T value, T step, T recStep);

/// Returns multiplicaton factor for reaching -60 dB after 'samples' iterations.
double t60(double samples);

/// Warp a value in [-1,1] to a sine-like value with domain [-pi, pi] and range [-1,1]
TEM T warpSinSS(T v);

/// Warp a value in [-1,1] to a sine-like value with domain [-pi, pi] and range [0,1]
TEM T warpSinSU(T v);

/// Warp a value in [0,1] to a sine-like value with domain [-pi, pi] and range [-1,1]
TEM T warpSinUS(T v);

/// Warp a value in [0,1] to a sine-like value with domain [-pi, pi] and range [0,1]
TEM T warpSinUU(T v);

/// Returns value wrapped in [lo, hi).
TEM T wrap(T value, T hi=(T)1, T lo=(T)0);

/// Returns value wrapped in [lo, hi).

/// 'numWraps' reports how many wrappings occured where the sign, + or -,
/// signifies above 'hi' or below 'lo', respectively.
TEM T wrap(T value, long & numWraps, T hi=(T)1, T lo=(T)0);

/// Returns value incremented by 1 and wrapped into interval [0, max).
TEM T wrapAdd1(T v, T max){ ++v; return v == max ? 0 : v; }

/// Like wrap(), but only adds or subtracts 'hi' once from value.
TEM T wrapOnce(T value, T hi=(T)1);

TEM T wrapOnce(T value, T hi, T lo);

TEM T wrapPhase(T radians);			///< Returns value wrapped in [-pi, pi).
TEM T wrapPhaseOnce(T radians);		///< Like wrapPhase(), but only wraps once.

/// Returns a continuous value measuring how close to zero the value is.

/// The graph of this function resembles a resonant peak. The function uses the
/// square of the value for evaluation.
TEM T zero(T v, T bw);

/// Same as zero(), but takes a unary function to be applied to the value.
template <class T, class F>
T zero(T v, T bw, F f);


TEM T notch(T v, T bw){ return v/(v+bw); }
TEM T notch1(T v, T bw){ return notch(scl::abs(v), bw); }

// Same as notch1, but 'amt' controls notch depth.
TEM T notch1(T v, T bw, T amt){ return (T)1 - amt*peak1(v, bw); }
TEM T notch2(T v, T bw){ return notch(v*v, bw); }

/// Peak function.

/// When v is the output of y=x^2,
/// this formula approximates the true formula for a resonant peak
/// 1/(1 - 2rcos(theta) + r^2). The argument bw is equivalent to (1-r)^2.
/// In general, the approximation has a slightly
/// smaller bandwidth than the true response. Also, the true response is
/// periodic, while this one is not.
TEM T peak(T v, T bw){ return bw/(v+bw); }
TEM T peak1(T v, T bw){ return bw/(scl::abs(v)+bw); }



//
// Analysis
//	

/// Returns number of bits set to 1.

/// From "Bit Twiddling Hacks", 
/// http://graphics.stanford.edu/~seander/bithacks.html#CountBitsSetParallel
uint32_t bitsSet(uint32_t v);

/// Returns whether or not an integer value is even.
TEM bool even(T v);

/// Returns whether the absolute value is less than an epsilon.
TEM bool lessAbs(T v, T eps=(T)0.000001);

/// Returns maximum of two values.
TEM T max(T v1, T v2);

/// Returns maximum of three values.
TEM T max(T v1, T v2, T v3);

/// Returns mean of two values.
TEM T mean(T v1, T v2);

/// Returns minimum of two values.
TEM T min(T v1, T v2);

/// Returns minimum of three values.
TEM T min(T v1, T v2, T v3);

/// Returns next largest value of 'val' that is a multiple of 'multiple'.
TEM T nextMultiple(T val, T multiple);

/// Returns whether or not an integer value is odd.
TEM bool odd(T v);

/// Returns whether the value is a power of two.
bool powerOf2(int v);

/// Returns slope of line passing through two points.
TEM T slope(T x1, T y1, T x2, T y2);

/// Returns number of trailing zeros in 32-bit int

/// This implements an algorithm from the paper 
/// "Using de Bruijn Sequences to Index 1 in a Computer Word"
/// by Charles E. Leiserson, Harald Prokof, and Keith H. Randall.
uint32_t trailingZeroes(uint32_t v);

/// Returns whether value is within [lo, hi].
TEM bool within(T v, T lo, T hi);
TEM bool within3(T v1, T v2, T v3, T lo, T hi);

/// Returns whether value is within [lo, hi).
TEM bool withinIE(T v, T lo, T hi);

/// Returns whether a positive zero crossing occured.
bool zeroCrossP(float prev, float now);



/// Convert polar to rectangular coordinates
TEM void polarToRect(T mag, T phs, T& real, T& imag);

/// Convert rectangular coordinates to polar.
TEM void rectToPolar(T& r, T& i);

/// In-place spherical to cartesian conversion for floating point types.
TEM void sphericalToCart(T & rho, T & phi, T & theta);




//---- Waveform generators

//---- Bipolar waveforms [-1, 1)
float rampDown	(uint32_t phase);	///< Returns value of bipolar downward ramp function.
float rampUp	(uint32_t phase);	///< Returns value of bipolar upward ramp function.
float square	(uint32_t phase);	///< Returns value of bipolar square function.
float triangle	(uint32_t phase);	///< Returns value of bipolar triangle function.

/// Returns value of bipolar pulse function (rampDown() + rampUp()).
float pulse		(uint32_t phase, uint32_t width);

/// Returns value of bipolar stair function (square() + square()).
float stair		(uint32_t phase, uint32_t width);

/// Returns value of bipolar dual upward ramp function (rampUp() + rampUp()).
float rampUp2	(uint32_t phase, uint32_t width);	// rampUp + rampUp

//---- Unipolar waveforms [0, 1)
float pulseU	(uint32_t phase, uint32_t width);	///< Returns value of unipolar pulse function.
float rampUpU	(uint32_t phase);	///< Returns value of unipolar downward ramp function.
float rampUp2U	(uint32_t phase);	///< Returns value of unipolar upward ramp2 function.
float rampDownU	(uint32_t phase);	///< Returns value of unipolar upward ramp function.
float squareU	(uint32_t phase);	///< Returns value of unipolar square function.
float stairU(uint32_t phase, uint32_t width); ///< Returns value of unipolar stair function.
float triangleU	(uint32_t phase);	///< Returns value of unipolar triangle function.

/// Dirichlet kernel, an impulse with n harmonics.

/// This is a closed form solution to the summation:\n
/// Dn(x)	= 1 + 2 * sum(1,n){ cos( kx ) } \n
///	Dn(x)	= sin( nx + x/2 ) / sin( x/2 )
float dirichlet(float phase, float n);

float dirichlet2(float phase, float n);

TEM T bartlett(T nphase);				///< Bartlett window. nphase => [-1, 1)
TEM T blackman(T phase);				///< Blackman window function.
TEM T blackmanHarris(T phase);			///< Blackman-Harris window function.
TEM T hamming(T phase);					///< Hamming window function.
TEM T hann(T phase);					///< von Hann window function.
TEM T raisedCosine(T phase, T a, T b);	///< Raised cosine f(x) = a - b cos(x).
TEM T welch(T nphase);					///< Welch window function. nphase => [-1, 1)


// internal
namespace{

//	const float mFactorial12f[13] = {
//		1.f, 1.f, 2.f, 6.f, 24.f, 120.f, 720.f, 5040.f, 40320.f, 
//		362880.f, 3628800.f, 39916800.f, 479001600.f
//	};

	const uint32_t mFactorial12u[13] = {
		1, 1, 2, 6, 24, 120, 720, 5040, 40320, 
		362880, 3628800, 39916800, 479001600
	};

	const uint32_t deBruijnBitPosition[32] = {
		 0,  1, 28,  2, 29, 14, 24,  3, 30, 22, 20, 15, 25, 17,  4,  8, 
		31, 27, 13, 23, 21, 19, 16,  7, 26, 12, 18,  6, 11,  5, 10,  9
	};
	
	inline uint32_t deBruijn(uint32_t v){
		return deBruijnBitPosition[(uint32_t(v * 0x077CB531UL)) >> 27];
	}

	const unsigned char mPrimes54[54] = {
	/*	  0    1    2    3    4    5    6    7    8    9   */
		  2,   3,   5,   7,  11,  13,  17,  19,  23,  29, // 0
		 31,  37,  41,  43,  47,  53,  59,  61,	 67,  71, // 1
		 73,  79,  83,  89,  97, 101, 103, 107, 109, 113, // 2
		127, 131, 137, 139, 149, 151, 157, 163, 167, 173, // 3
		179, 181, 191, 193, 197, 199, 211, 223, 227, 229, // 4
		233, 239, 241, 251								  // 5
	};
	
	TEM T taylorFactor3(T vv, T c1, T c2, T c3);
	TEM T taylorFactor4(T vv, T c1, T c2, T c3, T c4);
	TEM T taylorFactor5(T vv, T c1, T c2, T c3, T c4, T c5);
}


// Implementation_______________________________________________________________

#define GEN(t, f) template<> inline t abs<t>(t v){ return f(v); }
GEN(int, ::abs) GEN(long, labs) GEN(long long, llabs) GEN(float, fabsf) GEN(double, fabs)
TEM inline T abs(T v){ return v < T(0) ? -v : v; }
#undef GEN


TEM T atan2Fast(T y, T x){
	
	T r, angle;
	T ay = scl::abs(y) + (T)1e-10;      // kludge to prevent 0/0 condition
	
	if(x < (T)0){
		r = (x + ay) / (ay - x);
		angle = (T)M_3PI_4;
	}
	else{
		r = (x - ay) / (x + ay);
		angle = (T)M_PI_4;
	}
	
	angle += ((T)0.1963*r*r - (T)0.9817)*r;
	return y < (T)0 ? -angle : angle;
}

TEM inline T atLeast(T v, T e){	return v > (T)0 ? max(v, e) : min(v, -e); }

TEM inline T bump(T x){ return abs(x)<T(1) ? ::exp(T(-1)/(T(1) - x*x)) : T(0); }

TEM inline T ceil(T v){ return round(v + roundEps<T>()); }
TEM inline T ceil(T v, T s){ return ceil(v/s)*s; }
TEM inline T ceil(T v, T s, T r){ return ceil(v*r)*s; }

inline uint32_t ceilEven(uint32_t v){ return v += v & 1UL; }

inline uint32_t ceilPow2(uint32_t v){
	v--;
	v |= v >> 1;
	v |= v >> 2;
	v |= v >> 4;
	v |= v >> 8;
	v |= v >> 16;
	return ++v;
}

TEM inline T clip(T v, T hi, T lo){
	     if(v < lo) return lo;
	else if(v > hi)	return hi;
	return v;
}

TEM inline T clip(T v, int & clipFlag, T hi, T lo){
	clipFlag = 0;
	     if(v < lo){ clipFlag = -1; return lo; }
	else if(v > hi){ clipFlag =  1; return hi; }
	return v;
}

TEM inline T clipS(T v, T hi){ return clip(v, hi, -hi); }

template <class T, template <class> class V>
T curvature(const V<T>& a, const V<T>& b, const V<T>& c){

	V<T> d1b = b-a;				// first backward difference
	V<T> d1f = c-b;				// first forward difference
	V<T> d1  = (d1f+d1b) * 0.5;	// first mid difference
	
	V<T> d2  = d1f - d1b;		// second difference
	
	T d1n = d1.norm();
	
	return (d1.cross(d2)).norm() / (d1n*d1n*d1n);
}

TEM inline T dot2(T x1, T x2, T y1, T y2){ return x1 * y1 + x2 * y2; }

TEM inline T equals(T v1, T v2, T bw){ return equals(v1, v2, bw, pow2<T>); }

template <class T, class F>
inline T equals(T v1, T v2, T bw, F f){ return bw/(f(v1-v2) + bw); }

TEM inline void fadeLin(T & w1, T & w2, T f){ w1 = (T)1 - f; w2 = f; }

TEM inline void fadeTri(T & w1, T & w2, T f){
	if(f < 0.25){
		w1 = f * 4;
		w2 = 0;
	}
	else if( (f >= 0.25) && (f < 0.75) ){
		fadeLin(w1, w2, (T)(f * 2 - 0.5));
	}
	else{
		w1 = 0;
		w2 = (1 - f) * 4;
	}
}

//inline float factorial12(float v){ return mFactorial12f[(uint32_t)v]; }
inline uint32_t factorial12(uint32_t v){ return mFactorial12u[v]; }
TEM inline T floor(T v){ return round(v - roundEps<T>()); }
TEM inline T floor(T v, T s){ return floor(v/s)*s; }
TEM inline T floor(T v, T s, T r){ return floor(v*r)*s; }

inline uint32_t floorPow2(uint32_t v){
	v |= v >> 1;
	v |= v >> 2;
	v |= v >> 4;
	v |= v >> 8;
	v |= v >> 16;
	return (v >> 1) + 1;
}

TEM inline T fold(T v, T hi, T lo){
	long t;
	return fold(v,t,hi,lo);
}

TEM inline T fold(T v, long& numFolds, T hi, T lo){
	long numWraps;
	v = wrap(v, numWraps, hi, lo);
	if(numWraps & 1) v = hi + lo - v;
	numFolds = numWraps;
	return v;
}

TEM inline T foldOnce(T v, T hi, T lo){
	if(v > hi) return hi + (hi - v);
	if(v < lo) return lo + (lo - v);
	return v;
}

template <class V2>
void frenet(const V2& d1, V2& t, V2& n){
	t = d1;
	t *= invSqrt<2>(t.dot());
	n(-t[1], t[0]);	// normal according to right-hand rule
}

template <class V3>
void frenet(const V3& d1, const V3& d2, V3& t, V3& n, V3& b){	
	b = cross(d2, d1);
	n = cross(d1, b);
	t = d1 * invSqrt<2>(d1.magSqr());
	b *= invSqrt<2>(b.magSqr());
	n *= invSqrt<2>(n.magSqr());
}

template <class V3>
void frenet(const V3& p2, const V3& p1, const V3& p0, V3& t, V3& n, V3& b){
	//const V3 d1 = p0 - p1, d2 = d1 - (p1 - p2);
	const V3 d1 = (p0 - p2)*0.5;
	const V3 d2 = (d1 - p1)*2.0; // p0 - 2*p1 + p2 = p0 - p2 - 2*p1 = 2*d1 - 2*p1
	frenet(d1,d2, t,n,b);
}

TEM T gcd(const T& x, const T& y){
	if(y==T(0)) return x;
	return gcd(y, remainder(x,y));
}

TEM inline T hypot(T x, T y){ return ::sqrt(x*x + y*y); }

TEM inline T linLog2(T v, T recMin){
	v = log2Fast(scl::abs(v) + (T)0.000001);	// offset to avoid -inf
	return scl::max(v * recMin, (T)-1) + (T)1;
}

inline uint32_t log2(uint32_t v){ return deBruijn(ceilPow2(v)); }

inline float log2Fast(float v){
	Twiddle<float> u(v);
	return (float)((u.i - int32_t(Expo1<float>()))) * 0.0000001192092896f;// / 8388608.f;
}

TEM inline T mapInvPow2(T v){ return v*(T(2)-v); }

TEM inline void mapLin(T i0, T i1, T o0, T o1, T& scale, T& offset){
	scale = slope(i0, o0, i1, o1);
	offset = o0 - scale * i0;
}

TEM inline T mapLin(T v, T i0, T i1, T o0, T o1){
	float scale = slope(i0, o0, i1, o1);
	return (v - i0) * scale + o0;
}

inline double mapPower(double v, double b1, double b0, double p){
	if(p != 1.) v = ::pow(v, p);
	return b0 + (b1 - b0) * v;
}

TEM inline void mix2(T& io1, T& io2, T mix){
	T t1 = (io1 - io2) * mix;
	T t2 = io1 - t1;
	io1 = t1 + io2;
	io2 = t2;
	//io1 = io1 * mix + io2 * ((T)1 - mix); 
	//io2 = io2 * mix + io1 * ((T)1 - mix);
}

TEM inline void mulComplex(T & r1, T & i1, const T & r2, const T & i2){
	T rt = r1;
	r1 = r1 * r2 - i1 * i2;
	i1 = i1 * r2 + rt * i2;
}

TEM T nearest(T val, const char * interval, long div){
	long vr = castIntRound(val);
	long numWraps = 0;
	long vm = wrap(vr, numWraps, div, 0L);
	long sum = 0;

	while(*interval){
		if(vm <= sum){ vm = sum; break; }
		sum += (long)(*interval++ - 48);
	}
	
	return (T)(vm + numWraps * div);
}

TEM inline T nearestDiv(T of, T to){ return to / round(to/of); }
TEM inline T negative(T v, T bw){ return T(0.5) - sign(v, bw)*T(0.5); }
TEM inline T negative(T v, T bw, T a){ return a - (T(1)-a)*sign(v, bw); }
TEM inline T nextAfter(T x, T y){ return x<y ? x+1 : x-1; }
template<> inline float nextAfter(float x, float y){ return nextafterf(x,y); }
template<> inline double nextAfter(double x, double y){ return nextafter(x,y); }
template<> inline long double nextAfter(long double x, long double y){ return nextafterl(x,y); }
TEM inline T poly(T v, T a0, T a1, T a2){ return a0 + v*(a1 + v*a2); }
TEM inline T poly(T v, T a0, T a1, T a2, T a3){ return a0 + v*(a1 + v*(a2 + v*a3)); }
TEM inline T positive(T v, T bw){ return T(0.5) + sign(v, bw)*T(0.5); }
TEM inline T pow2 (T v){ return v*v; }
TEM inline T pow2S(T v){ return v*scl::abs(v); }
TEM inline T pow3 (T v){ return v*v*v; }
TEM inline T pow3Abs(T v){ return scl::abs(pow3(v)); }
TEM inline T pow4 (T v){ return pow2(pow2(v)); }
TEM inline T pow5 (T v){ return v * pow4(v); }
TEM inline T pow6 (T v){ return pow3(pow2(v)); }
TEM inline T pow8 (T v){ return pow4(pow2(v)); }
TEM inline T pow16(T v){ return pow4(pow4(v)); }
TEM inline T pow64(T v){ return pow8(pow8(v)); }

inline unsigned char prime(uint32_t n){ return mPrimes54[n]; }
TEM inline T prime(uint32_t n, T mul){ return (T)prime(n) * mul; }

TEM inline T ratioET(T pc, T divisions, T interval){
	return (T)::pow((double)interval, (double)pc / (double)divisions);
}


template<> inline float remainder<float>(const float& x, const float& y){ return ::remainderf(x,y); }
template<> inline double remainder<double>(const double& x, const double& y){ return ::remainder(x,y); }
template<> inline long double remainder<long double>(const long double& x, const long double& y){ return ::remainderl(x,y); }
TEM inline T remainder(const T& x, const T& y){ return x-(x/y)*y; }

//TEM inline T round(T v){ return (v + roundMagic<T>()) - roundMagic<T>(); }
TEM inline T round(T v){ double r=v; return (r + roundMagic) - roundMagic; }
TEM inline T round(T v, T s){ return round<double>(v/s) * s; }
TEM inline T round(T v, T s, T r){ return round<T>(v * r) * s; }
TEM inline T roundAway(T v){ return v<T(0) ? floor(v) : ceil(v); }
TEM inline T roundAway(T v, T s){ return v<T(0) ? floor(v,s) : ceil(v,s); }

//inline float round(float val, float step){
//	union { float f; unsigned long i; } u;
//	u.f = val;
//	u.i = u.i & 0x80000000 | 0x3F000000;		
//	//val = fastFloor(val * step + u.f) * stepRec;
//	return float(long(val/step + u.f)) * step;
//}
//
//inline float round(float val, float step, float stepRec){
//	union { float f; unsigned long i; } u;
//	u.f = val;
//	u.i = u.i & 0x80000000 | 0x3F000000;		
//	//val = fastFloor(val * step + u.f) * stepRec;
//	return float(long(val * stepRec + u.f)) * step;
//}

TEM Complex<T> sharm(int l, int m, T theta, T phi){
	T c = T(factorial12(l-m)) / T(factorial12(l+m));
	c = ::sqrt((2*l + 1) / M_4PI * c);
	c *= legendre(l, m, theta);
	
	phi *= m;
	return Complex<T>(c*cos(phi), c*sin(phi));
}

TEM inline T sign(T v, T bw){ return v/(scl::abs(v) + bw); }

TEM inline T taylorFactor3(T vv, T c1, T c2, T c3){
	return c1 * vv * (c2 - vv * (c3 - vv));
}
TEM inline T taylorFactor4(T vv, T c1, T c2, T c3, T c4){
	return c1 * vv * (c2 - vv * (c3 - vv * (c4 - vv)));
}
TEM inline T taylorFactor5(T vv, T c1, T c2, T c3, T c4, T c5){
	return c1 * vv * (c2 - vv * (c3 - vv * (c4 - vv * (c5 - vv))));
}


TEM inline T cosP3(T n){
	return (T)1 - (T)32 * n * n * ((T)0.75 - n);
}


#define t84 56.
#define t83 1680.
#define t82 20160.
#define t81 2.4801587302e-05
#define t73 42.
#define t72 840.
#define t71 1.9841269841e-04

TEM inline T cosT8(T r){
	
	if(r < (T)M_PI_4 && r > (T)-M_PI_4){
		float rr = r*r;
		return (T)1 - rr * (T)t81 * ((T)t82 - rr * ((T)t83 - rr * ((T)t84 - rr)));
	}
	else if(r > (T)0){
		r -= (T)M_PI_2;
		T rr = r*r;
		return -r * ((T)1 - (T)t71 * rr * ((T)t72 - rr * ((T)t73 - rr)));
	}
	else{
		r += (T)M_PI_2;
		float rr = r*r;
		return r * ((T)1 - (T)t71 * rr * ((T)t72 - rr * ((T)t73 - rr)));
	}
}


TEM inline T sinFast(const T& r){
    const T B = 4 / M_PI, C = -4 / (M_PI*M_PI);
    T y = B * r + C * r * scl::abs(r);
	const T P = 0.225; // const float Q = 0.775;
	return P * (y * scl::abs(y) - y) + y;   // Q * y + P * y * abs(y)
}

TEM inline T sinP7(T n){
	T nn = n*n;
	return n * ((T)3.138982 + nn * ((T)-5.133625 + nn * ((T)2.428288 - nn * (T)0.433645)));
}

TEM inline T sinP9(T n){	
	T nn = n*n;
	return n * ((T)3.1415191 + nn * ((T)-5.1662729 + nn * ((T)2.5422065 + nn * ((T)-0.5811243 + nn * (T)0.0636716))));
}

TEM inline T sinT7(T r){

//	if(r < (T)M_PI_4 && r > (T)-M_PI_4){
//		return r * ((T)1 - taylorFactor3<T>(r*r, t71, t72, t73));
//	}
//	else if(r > (T)0){
//		r -= (T)M_PI_2;
//		return (T)1 - taylorFactor4<T>(r*r, t81, t82, t83, t84);
//	}
//	else{
//		r += (T)M_PI_2;
//		return (T)-1 + taylorFactor4<T>(r*r, t81, t82, t83, t84);
//	}
	
	if(r < (T)M_PI_4 && r > (T)-M_PI_4){
		T rr = r*r;
		return r * ((T)1 - (T)t71 * rr * ((T)t72 - rr * ((T)t73 - rr)));
	}
	else if(r > (T)0){
		r -= (T)M_PI_2;
		float rr = r*r;
		return (T)1 - rr * (T)t81 * ((T)t82 - rr * ((T)t83 - rr * ((T)t84 - rr)));
	}
	else{
		r += (T)M_PI_2;
		float rr = r*r;
		return (T)-1 + rr * (T)t81 * ((T)t82 - rr * ((T)t83 - rr * ((T)t84 - rr)));
	}
}
#undef t84
#undef t83
#undef t82
#undef t81
#undef t73
#undef t72
#undef t71

#define ta5 90.
#define ta4 5040.
#define ta3 151200.
#define ta2 1814400.
#define ta1 2.7557319224e-07
#define t94 72.
#define t93 3024.
#define t92 60480.
#define t91 2.7557319224e-06

TEM inline T sinT9(T r){
	if(r < (T)M_PI_4 && r > (T)-M_PI_4){
		T rr = r*r;
		return r * ((T)1 - (T)t91 * rr * ((T)t92 - rr * ((T)t93 - rr * ((T)t94 - rr))));
	}
	else if(r > 0.){
		r -= (T)M_PI_2;
		T rr = r*r;
		return (T)1 - rr * (T)ta1 * ((T)ta2 - rr * ((T)ta3 - rr * ((T)ta4 - rr * ((T)ta5 - rr))));
	}
	else{
		r += (T)M_PI_2;
		T rr = r*r;
		return (T)-1 + rr * (T)ta1 * ((T)ta2 - rr * ((T)ta3 - rr * ((T)ta4 - rr * ((T)ta5 - rr))));
	}
}
#undef ta5
#undef ta4
#undef ta3
#undef ta2
#undef ta1
#undef t94
#undef t93
#undef t92
#undef t91

TEM T sinc(T r, T eps){ return (scl::abs(r) > eps) ? sin(r)/r : cos(r); }

TEM inline void sort(T& v1, T& v2){ if(v1>v2){ T t=v1; v1=v2; v2=t; } } 

inline double t60(double samples){ return ::pow(0.001, 1./samples); }

TEM inline T trunc(T v){ return round( (v > (T)0) ? v-roundEps<T>() : v+roundEps<T>() ); }
TEM inline T trunc(T v, T s){ return trunc(v/s)*s; }
TEM inline T trunc(T v, T s, T r){ return trunc(v*r)*s; }

/*

1. if	exponent is greater than 011111110
   then	continue
2. convert last n fraction bits to int

0 01111110 Fffffffffffffffffffffff	[1/2, 1)

0 01111111 Fffffffffffffffffffffff	[1, 2)

0 10000000 Fffffffffffffffffffffff	[2, 4)
1 00000000 111
*/

//inline float wrap1(float value){
//	uint32_t valueU = *(uint32_t *)&value;
//	
//	if((valueU & 0x7fffffff) > 0x3f7fffff){
//		uint32_t shift = floatExponent(value) - 127;
//		valueU = (valueU << shift) & 0x007fffff | 0x3f800000;
//		return *(float *)&valueU - 1.f;
//	}
//	else{
//		return value;
//	}
//	
//	
//}

TEM inline T warpSinSS(T v){ return v*((T)1.50 - v*v*(T)0.50); }
TEM inline T warpSinSU(T v){ return v*((T)0.75 - v*v*(T)0.25) + (T)0.5; }
TEM inline T warpSinUS(T v){ return v*v*((T)6 - v*(T)4) - (T)1; }
TEM inline T warpSinUU(T v){ return v*v*(v*(T)-2 + (T)3); }

TEM inline T wrap(T v, T hi, T lo){
	if(lo == hi) return lo;
	
	//if(v >= hi){
	if(!(v < hi)){
		T diff = hi - lo;
		v -= diff;
		if(!(v < hi)) v -= diff * (T)(uint32_t)((v - lo)/diff);
	}
	else if(v < lo){
		T diff = hi - lo;
		v += diff;	// this might give diff if range is too large, so check at end of block...
		if(v < lo) v += diff * (T)(uint32_t)(((lo - v)/diff) + 1);
		if(v==diff) return nextAfter(v, lo);
	}
	return v;
}

TEM inline T wrap(T v, long& numWraps, T hi, T lo){
	if(lo == hi){ numWraps = 0xFFFFFFFF; return lo; }
	
	T diff = hi - lo;
	numWraps = 0;
	
	if(v >= hi){
		v -= diff;
		if(v >= hi){
			numWraps = (long)((v - lo)/diff);
			v -= diff * T(numWraps);
		}
		++numWraps;
	}
	else if(v < lo){
		v += diff;
		if(v < lo){
			numWraps = (long)((v - lo)/diff) - 1;
			v -= diff * T(numWraps);
		}
		--numWraps;
	}
	return v;
}

TEM inline T wrapOnce(T v, T hi){
	     if(v >= hi ) return v - hi;
	else if(v < (T)0) return v + hi;
	return v;
}

TEM inline T wrapOnce(T v, T hi, T lo){
	     if(v >= hi) return v - hi + lo;
	else if(v <  lo) return v + hi - lo;
	return v;
}

TEM inline T wrapPhase(T r){
	if(r >= (T)M_PI){
		r -= (T)M_2PI;
		if(r < (T)M_PI) return r;
	}
	else if (r < (T)-M_PI){
		r += (T)M_2PI;
		if(r >= (T)-M_PI) return r;
	}
	else return r;
	
	return r - (T)M_2PI * scl::floor<T>((r + (T)M_PI) * (T)M_1_2PI);
}

TEM inline T wrapPhaseOnce(T r){
	if(r >= (T)M_PI)		return r - (T)M_2PI;
	else if(r < (T)-M_PI)	return r + (T)M_2PI;
	return r;
}

TEM inline T zero(T v, T bw){ return zero(v, bw, scl::pow2<T>); }

template <class T, class F>
inline T zero(T v, T bw, F f){ return bw/(f(v) + bw); }



inline uint32_t bitsSet(uint32_t v){
	v = v - ((v >> 1) & 0x55555555);                    // reuse input as temporary
	v = (v & 0x33333333) + ((v >> 2) & 0x33333333);     // temp
	return ((v + ((v >> 4) & 0xF0F0F0F)) * 0x1010101) >> 24; // count
}

TEM inline bool even(T v){ return 0 == odd(v); }

TEM inline bool lessAbs(T v, T eps){ return scl::abs(v) < eps; }
TEM inline T max(T v1, T v2){ return v1<v2?v2:v1; }
TEM inline T max(T v1, T v2, T v3){ return max(max(v1,v2),v3); }
TEM inline T mean(T v1, T v2){ return (v1 + v2) * (T)0.5; }
TEM inline T min(T v1, T v2){ return v1<v2?v1:v2; }
TEM inline T min(T v1, T v2, T v3){ return min(min(v1,v2),v3); }

TEM inline T nextMultiple(T v, T m){
	uint32_t div = (uint32_t)(v / m);	
	return (T)(div + 1) * m;
}

TEM inline bool odd(T v){ return v & T(1); }

TEM inline T slope(T x1, T y1, T x2, T y2){ return (y2 - y1) / (x2 - x1); }

inline uint32_t trailingZeroes(uint32_t v){ return deBruijn(v & -v); }

TEM inline bool within  (T v, T lo, T hi){ return !((v < lo) || (v > hi)); }
TEM inline bool withinIE(T v, T lo, T hi){ return (!(v < lo)) && (v < hi); }

TEM inline bool within3(T v1, T v2, T v3, T lo, T hi){
	return within(v1,lo,hi) && within(v2,lo,hi) && within(v3,lo,hi);
}

inline bool zeroCrossP(float prev, float now){
	union{ float f; int32_t i; } u1, u0;
	u1.f = prev; u0.f = now;
	return (u1.i <= 0) && (u0.i > 0);
	//(prevI <= 0) && (nowI > 0)
//	prevI += 0x00800000;
//	nowI += 0x00800000;
//	(now & ~prev)>>31
}


TEM inline void polarToRect(T m, T p, T& r, T& i){
	r = m * cos(p);
	i = m * sin(p);
	//printf("%f %f %f %f\n", m, p, r, i);
}

TEM Vec3<T> productZXZ(const Complex<T>& a, const Complex<T>& b, const Complex<T>& c){
	return Vec3<T>(
		a.r*b.i - a.i*b.r*c.i,
		a.i*b.i + a.r*b.r*c.i,
		b.r*c.r
	);
}

TEM inline void rectToPolar(T& r, T& i){
//	i = atan2(i, r);
//	r /= cos(i);
	T m = scl::hypot(i, r);
	i = atan2(i, r);
	r = m;
}

TEM inline void sphericalToCart(T& r, T& p, T& t){
	T rsinp = r * sin(p);
	T tt = t;
	t = r * cos(p);
	r = rsinp * cos(tt);
	p = rsinp * sin(tt);
}



// [1, 0.5, 0, -0.5]
// Freq precision:	32 bits
// Amp precision:	24 bits
inline float rampDown(uint32_t p){
	p = (p >> 9) | 0x40000000;
	return 3.f - punUF(p);
}

// [-1, -0.5, 0, 0.5]
// Freq precision:	32 bits
// Amp precision:	24 bits
inline float rampUp(uint32_t p){
	p = (p >> 9) | 0x40000000;
	return punUF(p) - 3.f;
}

// [1, 1, -1, -1] 
// Freq precision:	31 bits
// Amp precision:	NA
inline float square(uint32_t p){
//	phase = 0x3f800000 | (phase & 0x80000000);
//	return *(float *)&phase;
	//return p & 0x80000000 ? -1.f : 1.f;
	return p & MaskSign<float>() ? -1.f : 1.f;
}

// [-1, 0, 1, 0]
// Freq precision:	31 bits
// Amp precision:	25 bits
inline float triangle(uint32_t p){
	uint32_t dir = p >> 31;
	p = ((p^(-dir)) + dir);
	p = (p >> 8) | 0x40000000;
	return punUF(p) - 3.f;
}

// Just another triangle wave algorithm
//inline float triangle(uint32_t phase){	
//	uint32_t dir = phase & 0x80000000;
//	dir |= 0x40000000;
//	phase = (phase << 1 >> 9) | dir;
//	dir |= 0x400000;	// make it +/-3
//	return *(float *)&phase - *(float *)&dir;
//}

// and another...
//inline float triangle(uint32_t phase){
//	return rampUp(phase<<1) * square(phase);
//}

// Freq precision:	32 bits
// Amp precision:	24 bits
// Width precision:	32 bits
inline float pulse(uint32_t p, uint32_t w){
	// output floating point exponent should be [1, 2)
	uint32_t saw1 = ((p-w) >> 9) | Expo1<float>();
	uint32_t saw2 = ( p    >> 9) | Expo1<float>();
	return punUF(saw1) - punUF(saw2);
}

inline float stair(uint32_t p, uint32_t w){
	uint32_t sqr1 = 0x3f000000 | ( p    & MaskSign<float>());
	uint32_t sqr2 = 0x3f000000 | ((p+w) & MaskSign<float>());
	return punUF(sqr1) + punUF(sqr2);
}

inline float stairU(uint32_t p, uint32_t w){
	return ((p & MaskSign<float>()) ? 0.5f : 0.f) + (((p+w) & MaskSign<float>()) ? 0.5f : 0.f);
}
	
inline float pulseU(uint32_t p, uint32_t w){
	return p > w ? 0.f : 1.f;
}

inline float rampUp2(uint32_t p, uint32_t w){
	uint32_t saw1 = ( p    >> 9) | Expo1<float>();
	uint32_t saw2 = ((p+w) >> 9) | Expo1<float>();
	return punUF(saw1) + punUF(saw2) - 3.f;
}

// [0, 0.25, 0.5, 0.75]
// Freq precision:	32 bits
// Amp precision:	24 bits
inline float rampUpU(uint32_t p){
	p = (p >> 9) | Expo1<float>();
	return punUF(p) - 1.f;
}
	
inline float rampUp2U(uint32_t p, uint32_t w){
	uint32_t saw1 = ( p    >> 9) | 0x3F000000;
	uint32_t saw2 = ((p+w) >> 9) | 0x3F000000;
	return punUF(saw1) + punUF(saw2) - 1.f;
}

// [1, 0.75, 0.5, 0.25]
// Freq precision:	32 bits
// Amp precision:	24 bits
inline float rampDownU(uint32_t p){
	p = (p >> 9) | 0xbf800000;
	return punUF(p) + 2.f;
}

inline float squareU(uint32_t p){
//	phase = (phase & 0x80000000) >> 1;
//	return *(float *)&phase * 0.5f;
	return p & MaskSign<float>() ? 0.f : 1.f;
}

// [1, 0.5, 0, 0.5]
// Freq precision:	32 bits
// Amp precision:	24 bits
inline float triangleU(uint32_t p){
	union{ float f; uint32_t i; } u;
	u.i = (p >> 9) | 0x40000000;
	u.f -= 3.f;
	u.i &= 0x7fffffff;
	return u.f;
}

#define EPS 0.0000001
inline float dirichlet(float p, float n){
	float den = sin(p * 0.5f);
	if( scl::abs(den) < EPS )	return n;
	return sin(p * (n + 0.5f)) / den;
}

inline float dirichlet2(float p, float n){
	float den = sin(p);
	if( scl::abs(den) < EPS ){
		p = scl::wrapPhase(p);		
		return (p > -M_PI_2 && p < M_PI_2) ? n : -n;
	}
	return sin(p * n) / den;
}
#undef EPS


TEM inline T bartlett(T n){	return (T)1 - scl::abs(n); }

TEM inline T blackman(T r){
	return (T)0.42 + (T)0.08 * cos((T)2. * r) - (T)0.5 * cos(r);	// prevents -0s
}
TEM inline T blackmanHarris(T r){
	return (T)0.35875 - (T)0.48829 * cos(r) + (T)0.14128 * cos((T)2. * r) - (T)0.01168 * cos((T)3. * r);
}
TEM inline T hamming(T r){ return raisedCosine(r, (T)0.53836, (T)0.46164); }
TEM inline T hann(T r){ return raisedCosine(r, (T)0.5, (T)0.5); }
TEM inline T raisedCosine(T r, T a, T b){ return a - b * cos(r); }
TEM inline T welch(T n){ return (T)1 - n*n; }

} // scl::

TEM inline double norm(const T& v){ return scl::abs(v); }
TEM inline double normCompare(const T& v){ return norm(v); }

} // gam::

#undef TEM
#endif
