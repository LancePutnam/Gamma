#ifndef GAMMA_SCL_H_INC
#define GAMMA_SCL_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include <math.h>
#include <stdio.h>				/* printf */
#include <stdlib.h>				/* labs(long) */
#include "Constants.h"
#include "mem.h"

#include "MacroD.h"

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

/// Scalar rank functions for numerical types.
namespace scl{

const float justUnder1f = 0.99999997f; //0x3f7fffff;

/// Returns absolute value.
TEM T abs(T value);

/// Fast approximation to atan2().

// Author: Jim Shima, http://www.dspguru.com/comp.dsp/tricks/alg/fxdatan2.htm.
// |error| < 0.01 rad
TEM T atan2Fast(T y, T x);

/// Returns value clipped ouside of range [-eps, eps]
TEM T atLeast(T v, T eps);

/// Convert a string of 1s and 0s to an integer.
uint32_t bits(const char * string);

/// Returns floating point value rounded to next highest integer.
TEM T ceil(T val);
TEM T ceil(T val, T step);
TEM T ceil(T val, T step, T recStep);

/// Returns power of two ceiling of value.

/// This uses an algorithm devised by Sean Anderson, Sep. 2001.
/// From "Bit Twiddling Hacks", http://graphics.stanford.edu/~seander/bithacks.html.
ULONG ceilPow2(ULONG value);

ULONG ceilEven(ULONG value);		///< Returns even number ceiling.

/// Returns value clipped to [lo, hi].
TEM T clip(T value, T hi=(T)1, T lo=(T)0);

/// Returns value clipped to [-hi, hi].
TEM T clipS(T value, T hi=(T)1);

/// Returns value clipped to [lo, hi] and signifies clipping behavior.

/// clipFlag signifies if and where clipping occured.  0 means no clipping
/// occured, -1 means clipping occured at the lower bound, and 1 means
/// clipping at the upper bound.
TEM T clip(T value, int & clipFlag, T hi, T lo);

/// Returns value whose magnitude is clipped to [min, max].
float clipMag(float value, float min, float max);

void colorRGB(float h, float s, float v, float &r, float &g, float &b);
void colorHSV(float r, float g, float b, float &h, float &s, float &v);

/// Third order polynomial approximation to first half of cosine.

/// 'normal' must be in the range [0, 0.5] which corresponds to the first half
/// of the cosine.
TEM T cosP3(T normal);

/// 8th order Taylor series approximation to a cosine.

/// 'radians' must be in [-pi, pi].
///
TEM T cosT8(T radians);

/// Cross product of 3-element arrays.
template <class T1, class T2, class T3>
void cross(const T1& a, const T2& b, T3& r);

/// Returns two element dot product x1 * y1 + x2 * y2.
TEM T dot2(T x1, T x2, T y1, T y2);

/// Returns weights for linear fade.
TEM void fadeLin(T & weight1, T & weight2, T fade);

/// Returns weights for triangular window fade.

/// The weights returned are from two overlapping normalized unipolar 
/// triangular windows in the fade interval [0, 2/3] and [1/3, 1].\n
/// fade	weight1		weight2 \n
/// 0.25	1			0       \n
/// 0.5		0.5			0.5		\n
/// 0.75	0			1		
TEM void fadeTri(T & weight1, T & weight2, T fade);

float factorial12(float value);				///< Returns factorial of value in [0, 12].
ULONG factorial12(ULONG value);				///< Returns factorial of value in [0, 12].

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
ULONG floorPow2(ULONG value);

/// Returns value folded into [lo, hi].

/// For out-of-range values, the boundaries act like mirrors reflecting
/// the value into the range. For an even number of periods out of the range
/// this is identical to a wrap().
TEM T fold(T value, T hi=(T)1, T lo=(T)0);

/// Returns value folded into [lo, hi] one time.
TEM T foldOnce(T value, T hi=(T)1, T lo=(T)0);

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
TEM T gaussian(T v){ return exp(-v*v); }

TEM T hypot(T x, T y);

/// Convert linear value to log2 in range [0, 1]
TEM T linLog2(T v, T recMin);

/// Returns base 2 logarithm of value.

/// If the value is not an exact power of two, the logarithm of the next
/// highest power of two will taken.
/// This uses an algorithm devised by Eric Cole, Jan. 2006.
/// From "Bit Twiddling Hacks", http://graphics.stanford.edu/~seander/bithacks.html.
ULONG log2(ULONG value);

/// Fast base 2 logarithm.  For value <= 0, behavior is undefined.
float log2Fast(float value);

/// Maps value from [-1,1] to [depth, 1].
TEM T mapDepth(T v, T depth){ return (v - (T)1) * (T)0.5 * depth + (T)1;  }

/// Computes scale and offset necessary to map value from [i0, i1] to [o0, o1].
TEM void mapLin(T i0, T i1, T o0, T o1, T & scale, T & offset);

/// Linearly maps value from [i0, i1] to [o0, o1].
TEM T mapLin(T value, T i0, T i1, T o0, T o1);
//TEM T mapLin(T val, T i0, T i1=1, T o0=0, T o1=1);

/// Returns normal^power linearly mapped to [bound0, bound1].
double mapNormal(double normal, double bound1, double bound0, double power=1.);

/// Mixes two values together (1 = thru, 0.5 = mono, 0 = swap).
TEM void mix2(T& io1, T& io2, T mix);

/// Perform complex multiplication, c1 = c1 c2.
TEM void mulComplex(T& r1, T& i1, const T& r2, const T& i2);

/// Perform quaternion multiplication, q1 = q1 q2.
TEM void mulQuat(T& r1, T& i1, T& j1, T& k1, T r2, T i2, T j2, T k2);

TEM T nearest(T val, const char * interval = "2212221", long div=12);

/// Returns nearest integer division of one value to another
TEM T nearestDiv(T of, T to);

/// Returns pole radius given a bandwidth and sampling interval
TEM	inline T poleRadius(T bw, double ups){ return exp(-M_PI * bw * ups); }
//return (T)1 - (M_2PI * bw * ups); // linear apx for fn < ~0.02

/// Evaluates polynomial a0 + a1 x + a2 x^2.
TEM T poly(T x, T a0, T a1, T a2);

/// Evaluates polynomial a0 + a1 x + a2 x^2 + a3 x^3.
TEM T poly(T x, T a0, T a1, T a2, T a3);

/// Continuous positive map.

/// The return value is close to 1 if v > 0 and close to 0 if v < 0.
///
TEM T positive(T v, T bw);

TEM T pow2(T value);			///< Returns value to the 2nd power.
TEM T pow3(T value);			///< Returns value to the 3rd power.
TEM T pow4(T value);			///< Returns value to the 4th power.
TEM T pow5(T value);			///< Returns value to the 5th power.
TEM T pow6(T value);			///< Returns value to the 6th power.
TEM T pow8(T value);			///< Returns value to the 8th power.
TEM T pow16(T value);			///< Returns value to the 16th power.
TEM T pow2S(T value);			///< Returns value to the 2nd power preserving sign.
unsigned char prime(ULONG n);	///< Returns (n+1)th prime number up to n=53.
TEM T prime(ULONG n, T mul);	///< Returns scaled (n+1)th prime number up to n=53.

/// Returns pole radius given a T60 decay length and units/sample
inline double radius60(double dcy, double ups){ return exp(M_LN001/dcy * ups); } // u/s * 1/u

/// Returns equal temperament ratio- octave^(pc/divisions)
TEM T ratioET(T pc, T divisions=12, T octave=2);

/// Fast reciprocal of square root of value.

/// Lomont, C. 2003. "Fast inverse square root."
///
float recSqrtFast(float value);

/// Returns floating point value rounded to nearest integer.
TEM T round(T value);

/// Returns floating point value rounded to nearest step multiple.
TEM T round(T value, T step);

/// Returns floating point value rounded to nearest step multiple. Faster version to avoid 1/step divide.
TEM T round(T value, T step, T recStep);

/// Continuous sign map.

/// The return value is close to 1 if v > 0 and close to -1 if v < 0.
///
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
TEM T sinc(T radians, T eps=(T)0.0001);

/// Sort values so that value1 <= value2.
TEM void sort2(T & value1, T & value2);

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

/// Peak function.

/// When v is the output of y=x^2,
/// this formula approximates the true formula for a resonant peak
/// 1/(1 - 2rcos(theta) + r^2). The argument bw is equivalent to (1-r)^2.
/// In general, the approximation has a slightly
/// smaller bandwidth than the true response. Also, the true response is
/// periodic, while this one is not.
TEM T peak(T v, T bw){ return bw/(v+bw); }




//
// Analysis
//	

/// Returns number of bits set to 1.

/// From "Bit Twiddling Hacks", 
/// http://graphics.stanford.edu/~seander/bithacks.html#CountBitsSetParallel
ULONG bitsSet(ULONG v);

/// Returns whether or not the value is even.
bool even(ULONG v);

/// Returns whether the absolute value is less than an epsilon.
TEM bool lessAbs(T v, T eps=(T)0.000001);

/// Returns maximum of two values.
TEM T max(T v1, T v2);

/// Returns mean of two values.
TEM T mean(T v1, T v2);

/// Returns minimum of two values.
TEM T min(T v1, T v2);

/// Returns minimum of three values.
TEM T min(T v1, T v2, T v3);

/// Returns next largest value of 'val' that is a multiple of 'multiple'.
TEM T nextMultiple(T val, T multiple);

/// Returns whether the value is a power of two.
bool powerOf2(int v);

/// Returns slope of line passing through two points.
TEM T slope(T x1, T y1, T x2, T y2);

/// Returns number of trailing zeros in 32-bit int

/// This implements an algorithm from the paper 
/// "Using de Bruijn Sequences to Index 1 in a Computer Word"
/// by Charles E. Leiserson, Harald Prokof, and Keith H. Randall.
ULONG trailingZeroes(ULONG v);

/// Returns whether value is within [lo, hi].
TEM bool within(T v, T lo, T hi);
TEM bool within3(T v1, T v2, T v3, T lo, T hi);

/// Returns whether value is within [lo, hi).
TEM bool withinIE(T v, T lo, T hi);

/// Returns whether a positive zero crossing occured.
bool zeroCrossP(float prev, float now);


//
// Conversion
//

/// Convert decimal integer to ascii base-36 character
char base10To36(int dec10);

/// Convert ascii base-36 character to decimal integer 
int base36To10(char ascii36);

/// Convert 2-byte array to 16-bit unsigned integer.
uint16_t bytesToUInt16(const uint8_t * bytes2);

/// Convert 4-byte array to 32-bit unsigned integer.
uint32_t bytesToUInt32(const uint8_t * bytes4);

/// Cast value to signed integer using rounding.

/// This is much faster than using a standard C style cast.  Floor or
/// ceiling casts can be accomplished by subtracting or adding 0.5
/// from the input, respectively.
int32_t castIntRound(double value);

TEM long castIntTrunc(T value);

/// Returns biased decimal value of exponent field.

/// The true exponent is the return value minus 127. \n
/// For example, values in [0.5, 1) return 126 (01111110), so the true
///	exponent is 126 - 127 = -1.
ULONG floatExponent(float value);

/// Returns mantissa field as float between [0, 1).
float floatMantissa(float value);

/// Cast float to unsigned long.

/// Reliable up to 2^24 (16777216)
///
ULONG floatToUInt(float value);

/// Cast float to long.

/// Reliable up to 2^24 (16777216)
///
long floatToInt(float value);

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

/// Convert 16-bit signed integer to signed floating point normal [-1, 1).
float intToNormal(short value);

/// Convert floating point normal [0, 1) to unsigned long [0, 2^32)

/// This conversion is most accurate on an exponential scale.
///	Input values outside [-1, 1) return 0.
///	Values in [-1, 0] behave as positive values in [0, 1).
ULONG normalToUInt(float normal);

/// Convert floating point normal [0, 1) to unsigned long [0, 2^32)

/// This conversion is most accurate on a linear scale.
/// Input values outside [0, 1) result in undefined behavior.
ULONG normalToUInt2(float normal);

TEM void polarToRect(T mag, T phs, T& real, T& imag);

/// Maps a position in [-1, 1] to an index in [0, n). No boundary operations are performed.
inline int posToInd(float v, int n){ return n * (v*0.49999f + 0.5f); }

/// Type-pun 32-bit unsigned int to 32-bit float

/// This function uses a union to avoid problems with direct pointer casting
/// when fstrict-aliasing is on.
inline float punUF32(uint32_t v){ union{float f; uint32_t i;} u; u.i=v; return u.f; }

/// Type-pun 32-bit float to 32-bit unsigned int

/// This function uses a union to avoid problems with direct pointer casting
/// when fstrict-aliasing is on.
inline uint32_t punFU32(float v){ union{float f; uint32_t i;} u; u.f=v; return u.i; }

/// Convert rectangular coordinates to polar.
TEM void rectToPolar(T& r, T& i);

/// Get fractional and integer parts from a float.

/// Works reliably up to 2^24 == 16777216
/// Useful for linearly interpolated table lookups
float split(float value, long & intPart);

float splitInt512(ULONG value, ULONG & intPart);

/// Split integer accumulator into table index (size=1024) and interpolation fraction.
float splitInt1024(ULONG value, ULONG & intPart);

/// In-place spherical to cartesian conversion for floating point types.
TEM void sphericalToCart(T & rho, T & phi, T & theta);

TEM T uintToNormal (uint32_t value);
TEM T uintToNormalS(uint32_t value);

//---- Waveform generators

//---- Bipolar waveforms [-1, 1)
float rampDown	(ULONG phase);	///< Returns value of bipolar downward ramp function.
float rampUp	(ULONG phase);	///< Returns value of bipolar upward ramp function.
float square	(ULONG phase);	///< Returns value of bipolar square function.
float triangle	(ULONG phase);	///< Returns value of bipolar triangle function.

/// Returns value of bipolar pulse function (rampDown() + rampUp()).
float pulse		(ULONG phase, ULONG width);

/// Returns value of bipolar stair function (square() + square()).
float stair		(ULONG phase, ULONG width);

/// Returns value of bipolar dual upward ramp function (rampUp() + rampUp()).
float rampUp2	(ULONG phase, ULONG width);	// rampUp + rampUp

//---- Unipolar waveforms [0, 1)
float pulseU	(ULONG phase, ULONG width);	///< Returns value of unipolar pulse function.
float rampUpU	(ULONG phase);	///< Returns value of unipolar downward ramp function.
float rampUp2U	(ULONG phase);	///< Returns value of unipolar upward ramp2 function.
float rampDownU	(ULONG phase);	///< Returns value of unipolar upward ramp function.
float squareU	(ULONG phase);	///< Returns value of unipolar square function.
float stairU(ULONG phase, ULONG width); ///< Returns value of unipolar stair function.
float triangleU	(ULONG phase);	///< Returns value of unipolar triangle function.

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

//
// I/O utilities
//

/// Returns an ASCII character most closely matching an intensity value in [0,1].
inline char intensityToASCII(float v){
	static const char map[] =
	" .,;-~_+<>i!lI?/|)(1}{][rcvunxzjftLCJUYXZO0Qoahkbdpqwm*WMB8&%$#@";
	//"$@B%8&WM#*oahkbdpqwmZO0QLCJUYXzcvunxrjft/\|()1{}[]?-_+~<>i!lI;:,\"^`'. ";
//	 123456789.123456789.123456789.123456789.123456789.123456789.1234
	const int N  = sizeof(map)-1;
	return map[int((N*scl::clip(v, 0.9999999f)))];
}

TEM void print(T & v, const char * post="", const char * pre="", FILE * fp=stdout);

// Binary printing methods
void printBinary(ULONG value, const char * zero="0", const char * one="1", int msb=32);
void printBinary(unsigned long long value, const char * zero="0", const char * one="1", int msb=64);
void printBinary(float value, const char * zero="0", const char * one="1", int msb=32);
void printBinary(void * value32, const char * zero="0", const char * one="1", int msb=32);


TEM void print2D(T* pix, int nx, int ny, FILE * fp=stdout){
	for(int j=0; j<nx; ++j){
	for(int i=0; i<ny; ++i){
		float v = pix[j*nx + i];
		fprintf(fp, "%c ", scl::intensityToASCII(v));
	} printf("\n"); }
}

/// Print signed normalized value on a horizontal plot.

/// @param[in]	value	Normalized value to plot
/// @param[in]	width	Character width of plot excluding center point
/// @param[in]	spaces	Print extra filling spaces to the right
/// @param[in]	point	The print character for points
void printPlot(float value, ULONG width=50, bool spaces=true, const char * point="o");





// internal
namespace{

	const float mFactorial12f[13] = {
		0.f, 1.f, 2.f, 6.f, 24.f, 120.f, 720.f, 5040.f, 40320.f, 
		362880.f, 3628800.f, 39916800.f, 479001600.f
	};

	const ULONG mFactorial12u[13] = {
		0, 1, 2, 6, 24, 120, 720, 5040, 40320, 
		362880, 3628800, 39916800, 479001600
	};

	const ULONG multiplyDeBruijnBitPosition[32] = {
	  0, 1, 28, 2, 29, 14, 24, 3, 30, 22, 20, 15, 25, 17, 4, 8, 
	  31, 27, 13, 23, 21, 19, 16, 7, 26, 12, 18, 6, 11, 5, 10, 9
	};

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

TEM const T roundEps();

inline bool isLittleEndian(){ return endian == 0; }




// Implementation_______________________________________________________________


namespace{
	const double roundMagic = 6755399441055744.; // 2^52 * 1.5
}

	template<> inline const float  roundEps<float >(){ return 0.499999925f; }
	template<> inline const double roundEps<double>(){ return 0.499999985; }




#define GEN(t, f) template<> inline t abs<t>(t v){ return f(v); }
GEN(int, abs) GEN(long, labs) GEN(float, fabsf) GEN(double, fabs)
TEM inline T abs(T v){ return v < (T)0 ? -v : v; }
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

inline uint32_t bits(const char * string){
	uint32_t v=0; int n = strlen(string);
	for(int i=0; i<n; ++i) if(string[i] == '1') v |= 1<<(n-1-i);
	return v;
}

TEM inline T ceil(T v){ return round(v + roundEps<T>()); }
TEM inline T ceil(T v, T s){ return ceil(v/s)*s; }
TEM inline T ceil(T v, T s, T r){ return ceil(v*r)*s; }

inline ULONG ceilEven(ULONG v){ return v += v & 1UL; }

inline ULONG ceilPow2(ULONG v){
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

template <class T1, class T2, class T3>
inline void cross(const T1& a, const T2& b, T3& r){
	r[0] = a[1] * b[2] - a[2] * b[1];
	r[1] = a[2] * b[0] - a[0] * b[2];
	r[2] = a[0] * b[1] - a[1] * b[0];
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

inline float factorial12(float v){ return mFactorial12f[(ULONG)v]; }
inline ULONG factorial12(ULONG v){ return mFactorial12u[v]; }
TEM inline T floor(T v){ return round(v - roundEps<T>()); }
TEM inline T floor(T v, T s){ return floor(v/s)*s; }
TEM inline T floor(T v, T s, T r){ return floor(v*r)*s; }

inline ULONG floorPow2(ULONG v){
	v |= v >> 1;
	v |= v >> 2;
	v |= v >> 4;
	v |= v >> 8;
	v |= v >> 16;
	return (v >> 1) + 1;
}

TEM inline T fold(T v, T hi, T lo){
	long numWraps;
	v = wrap(v, numWraps, hi, lo);
	if(numWraps & 1) v = hi + lo - v;
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
	t *= recSqrtFast(t.dot());
	n(-t[1], t[0]);	// normal according to right-hand rule
}

template <class V3>
void frenet(const V3& d1, const V3& d2, V3& t, V3& n, V3& b){	
	t = d1;
	cross(d2, d1, b);
	cross(d1,  b, n);
	
	t *= recSqrtFast(t.dot());
	b *= recSqrtFast(b.dot());
	n *= recSqrtFast(n.dot());	
}

template <class V3>
void frenet(const V3& p2, const V3& p1, const V3& p0, V3& t, V3& n, V3& b){
	//const V3 d1 = p0 - p1, d2 = d1 - (p1 - p2);
	const V3 d1 = (p0 - p2)*0.5;
	const V3 d2 = (d1 - p1)*2.0; // p0 - 2*p1 + p2 = p0 - p2 - 2*p1 = 2*d1 - 2*p1
	frenet(d1,d2, t,n,b);
}

TEM inline T hypot(T x, T y){ return sqrt(x*x + y*y); }

TEM inline T linLog2(T v, T recMin){
	v = log2Fast(scl::abs(v) + (T)0.000001);	// offset to avoid -inf
	return scl::max(v * recMin, (T)-1) + (T)1;
}

inline ULONG log2(ULONG v){
	v = ceilPow2(v);
	return multiplyDeBruijnBitPosition[(v * 0x077CB531UL) >> 27];
}

inline float log2Fast(float v){
	union{ float f; int32_t i; } u; u.f=v;
	return (float)((u.i - 0x3f800000)) * 0.0000001192092896f;// / 8388608.f;
}

TEM inline void mapLin(T i0, T i1, T o0, T o1, T& scale, T& offset){
	scale = slope(i0, o0, i1, o1);
	offset = o0 - scale * i0;
}

TEM inline T mapLin(T v, T i0, T i1, T o0, T o1){
	float scale = slope(i0, o0, i1, o1);
	return (v - i0) * scale + o0;
}

inline double mapNormal(double v, double b1, double b0, double power){
	if(power != 1.) v = pow(v, power);
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

TEM inline void mulQuat(T & r1, T & i1, T & j1, T & k1, T r2, T i2, T j2, T k2){
	T rt = r1, it = i1, jt = j1;
	r1 = r1 * r2 - i1 * i2 - j1 * j2 - k1 * k2;
	i1 = rt * i2 + it * r2 + jt * k2 - k1 * j2;
	j1 = rt * j2 - it * k2 + jt * r2 + k1 * i2;
	k1 = rt * k2 + it * j2 - jt * i2 + k1 * r2;
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

TEM inline T poly(T v, T a0, T a1, T a2){ return a0 + v*(a1 + v*a2); }
TEM inline T poly(T v, T a0, T a1, T a2, T a3){ return a0 + v*(a1 + v*(a2 + v*a3)); }
TEM inline T positive(T v, T bw){ return sign(v, bw)*0.5 + 0.5; }
TEM inline T pow2 (T v){ return v*v; }
TEM inline T pow3 (T v){ return v*v*v; }
TEM inline T pow4 (T v){ return pow2(pow2(v)); }
TEM inline T pow5 (T v){ return v * pow4(v); }
TEM inline T pow6 (T v){ return pow3(pow2(v)); }
TEM inline T pow8 (T v){ return pow4(pow2(v)); }
TEM inline T pow16(T v){ return pow4(pow4(v)); }
TEM inline T pow64(T v){ v*=v; v*=v; v*=v; v*=v; v*v; return v*v; }

TEM inline T pow2S(T v){ return v*scl::abs(v); }

TEM inline T pow3Abs(T v){ return scl::abs(pow3(v)); }

inline unsigned char prime(ULONG n){ return mPrimes54[n]; }
TEM inline T prime(ULONG n, T mul){ return (T)prime(n) * mul; }

TEM inline T ratioET(T pc, T divisions, T interval){
	return (T)pow((double)interval, (double)pc / (double)divisions);
}

//TEM inline T round(T v){ return (v + roundMagic<T>()) - roundMagic<T>(); }
TEM inline T round(T v){ double r=v; return (r + roundMagic) - roundMagic; }
TEM inline T round(T v, T s){ return round<double>(v/s) * s; }
TEM inline T round(T v, T s, T r){ return round<T>(v * r) * s; }

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

TEM T sinc(T r, T eps=(T)0.0001){ return (scl::abs(r) > eps) ? sin(r)/r : cos(r); }

TEM inline void sort2(T & v1, T & v2){ if(v1 > v2) mem::swap(v1, v2); }

inline float recSqrtFast(float v){
	
	union{float f; uint32_t i; } u;
	
	float v2 = v * 0.5f;
	u.f = v;						// get bits for floating value 
	u.i = 0x5f375a86 - (u.i>>1);	// gives initial guess y0 
	v = u.f;						// convert bits back to float
	
	// Newton step, repeating increases accuracy 
	v *= 1.5f - v2 * v * v;
	
	return v; 
} 

inline double t60(double samples){ return pow(0.001, 1./samples); }

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
//	ULONG valueU = *(ULONG *)&value;
//	
//	if((valueU & 0x7fffffff) > 0x3f7fffff){
//		ULONG shift = floatExponent(value) - 127;
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
		if(!(v < hi)) v -= diff * (T)(ULONG)((v - lo)/diff);
	}
	else if(v < lo){
		T diff = hi - lo;
		v += diff;
		if(v < lo) v += diff * (T)(ULONG)(((lo - v)/diff) + 1);
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
			v -= diff * (T)numWraps;
		}
		numWraps++;
	}
	else if(v < lo){
		v += diff;
		if(v < lo){
			numWraps = (long)((v - lo)/diff) - 1;
			v -= diff * (T)numWraps;
		}
		numWraps--;
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



inline ULONG bitsSet(ULONG v){
	v = v - ((v >> 1) & 0x55555555);                    // reuse input as temporary
	v = (v & 0x33333333) + ((v >> 2) & 0x33333333);     // temp
	return ((v + (v >> 4) & 0xF0F0F0F) * 0x1010101) >> 24; // count
}

inline bool even(ULONG v){ return 0 == (v & 1); }

TEM inline bool lessAbs(T v, T eps){ return scl::abs(v) < eps; }
TEM inline T max(T v1, T v2){ return v1<v2?v2:v1; }
TEM inline T mean(T v1, T v2){ return (v1 + v2) * (T)0.5; }
TEM inline T min(T v1, T v2){ return v1<v2?v1:v2; }
TEM inline T min(T v1, T v2, T v3){ return (v1<v2 && v1<v3) ? v1 : scl::min(v2,v3); }

TEM inline T nextMultiple(T v, T m){
	ULONG div = (ULONG)(v / m);	
	return (T)(div + 1) * m;
}

TEM inline T slope(T x1, T y1, T x2, T y2){ return (y2 - y1) / (x2 - x1); }

inline ULONG trailingZeroes(ULONG v){
	return multiplyDeBruijnBitPosition[((v & -v) * 0x077CB531UL) >> 27];
}

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


inline int32_t castIntRound(double v){
	v += roundMagic;
	union{ double f; int32_t i[2]; } u; u.f = v;
	return u.i[endian]; // result in lsb
}

TEM inline long castIntTrunc(T v){
	return castIntRound( v + (v > (T)0 ? -roundEps<T>() : roundEps<T>()) );
}

inline ULONG floatExponent(float v){

	return scl::punFU32(v) << 1 >> 24;

//	PUN_F2I(v,u) ULONG exp = u.i; 
//	return exp << 1 >> 24;
}

inline float floatMantissa(float v){
	uint32_t frac = scl::punFU32(v);
	frac = frac & MASK_F32_FRAC | 0x3f800000;
	return scl::punUF32(frac) - 1.f;
}

inline float intToNormal(short v){
//	ULONG vu = ((ULONG)v) + 32768;
//	vu = vu << 7 | 0x40000000;
//	return *(float *)&vu - 3.f;

	uint32_t vu = (((uint32_t)v) + 0x808000) << 7; // set fraction in float [2, 4)
	return scl::punUF32(vu) - 3.f;

	//return (float)v / 32768.f; // naive method
	//return (float)v * 0.000030517578125f; // less naive method
}

/*
f32 range		u32 range		f32 exponent
[0.50,  1.00)	[2^31, 2^32)	01111110 (126)	
[0.25,  0.50)	[2^30, 2^31)	01111101 (125)	
[0.125, 0.25)	[2^29, 2^30)	01111100 (124)	

1. prepend 1 to fraction ( 'or' with 1<<24 (0x800000) )
2. shift left by 8
3. shift right according to exponent

0 01111110 00000000000000000000000

0 01111110 Fffffffffffffffffffffff	[1/2, 1/1)
1 Ffffffff fffffffffffffff00000000

0 01111101 Fffffffffffffffffffffff	[1/4, 1/2)
0 1Fffffff ffffffffffffffff0000000

0 01111100 Fffffffffffffffffffffff	[1/8, 1/4)
0 01Ffffff fffffffffffffffff000000

0 01111011 Fffffffffffffffffffffff	[1/16, 1/8)
0 001Fffff ffffffffffffffffff00000

effective  precision
9 x 2^24 + 2^23 + 2^22 + ... + 1
= 167,804,826

2^24
=  16,777,216

*/

inline ULONG normalToUInt(float v){
	ULONG normalU = punFU32(v);
	ULONG rbs = 126UL - (normalU >> 23UL);
//	printf("%x %lu\n", (normalU | 0x800000) << 8, rbs);
//	printf("%x\n", 0x80000000UL >> rbs);
	return (normalU | 0x800000UL) << 8UL >> rbs;

//Her00	
//float y = v + 1.f; 
//return ((unsigned long&)v) & 0x7FFFFF;      // last 23 bits 
}

inline ULONG normalToUInt2(float v){
	v++;	// go into [1,2] range, FP fraction is now result
	return punFU32(v) << 9;
}

TEM inline void polarToRect(T m, T p, T& r, T& i){
	r = m * cos(p);
	i = m * sin(p);
	//printf("%f %f %f %f\n", m, p, r, i);
}

TEM inline void rectToPolar(T& r, T& i){
//	i = atan2(i, r);
//	r /= cos(i);
	T m = scl::hypot(i, r);
	i = atan2(i, r);
	r = m;
}

template<> inline float uintToNormal<float>(ULONG v){
	v = v >> 9 | 0x3f800000; 
	return punUF32(v) - 1.f;
}

template<> inline float uintToNormalS<float>(ULONG v){
	v = v >> 9 | 0x40000000;
	return punUF32(v) - 3.f;
}

inline float splitInt512(ULONG v, ULONG & intPart){
	union{ float f; ULONG i; }u;
	u.i = v & 0x007fffff | 0x3f800000;
	intPart = v >> 22;
	return u.f - 1.f;
}

inline float splitInt1024(ULONG v, ULONG & intPart){
	union{ float f; ULONG i; }u;
	//u.i = v << 10 >> 10 | 0x3F800000;
	u.i = v << 1 & 0x007fffff | 0x3f800000;
	intPart = v >> 22;
	return u.f - 1.f;
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
inline float rampDown(ULONG p){
	p = (p >> 9) | 0x40000000;
	return 3.f - punUF32(p);
}

// [-1, -0.5, 0, 0.5]
// Freq precision:	32 bits
// Amp precision:	24 bits
inline float rampUp(ULONG p){
	p = (p >> 9) | 0x40000000;
	return punUF32(p) - 3.f;
}

// [1, 1, -1, -1] 
// Freq precision:	31 bits
// Amp precision:	NA
inline float square(ULONG p){
//	phase = 0x3f800000 | (phase & 0x80000000);
//	return *(float *)&phase;
	return p & 0x80000000 ? 1.f : -1.f;
}

// [-1, 0, 1, 0]
// Freq precision:	31 bits
// Amp precision:	25 bits
inline float triangle(ULONG p){
	ULONG dir = p >> 31;
	p = ((p^(-dir)) + dir);
	p = (p >> 8) | 0x40000000;
	return scl::punUF32(p) - 3.f;
}

// Just another triangle wave algorithm
//inline float triangle(ULONG phase){	
//	ULONG dir = phase & 0x80000000;
//	dir |= 0x40000000;
//	phase = (phase << 1 >> 9) | dir;
//	dir |= 0x400000;	// make it +/-3
//	return *(float *)&phase - *(float *)&dir;
//}

// and another...
//inline float triangle(ULONG phase){
//	return rampUp(phase<<1) * square(phase);
//}

// Freq precision:	32 bits
// Amp precision:	24 bits
// Width precision:	32 bits
inline float pulse(ULONG p, ULONG w){
	// output floating point exponent should be [1, 2)
	ULONG saw1 = (p >> 9) | 0x3F800000;
	ULONG saw2 = ((p+w) >> 9) | 0x3F800000;
	return scl::punUF32(saw1) - scl::punUF32(saw2);
}

inline float stair(ULONG p, ULONG w){
	ULONG sqr1 = 0x3f000000 | (p & 0x80000000);
	ULONG sqr2 = 0x3f000000 | ((p+w) & 0x80000000);
	return punUF32(sqr1) + punUF32(sqr2);
}

inline float stairU(ULONG p, ULONG w){
	return ((p & 0x80000000) ? 0.5f : 0.f) + (((p+w) & 0x80000000) ? 0.5f : 0.f);
}
	
inline float pulseU(ULONG p, ULONG w){
	return p > w ? 0.f : 1.f;
}

inline float rampUp2(ULONG p, ULONG w){
	ULONG saw1 = (p >> 9) | 0x3f800000;
	ULONG saw2 = ((p+w) >> 9) | 0x3f800000;
	return punUF32(saw1) + punUF32(saw2) - 3.f;
}

// [0, 0.25, 0.5, 0.75]
// Freq precision:	32 bits
// Amp precision:	24 bits
inline float rampUpU(ULONG p){
	p = (p >> 9) | 0x3f800000;
	return punUF32(p) - 1.f;
}
	
inline float rampUp2U(ULONG p, ULONG w){
	ULONG saw1 = (p >> 9) | 0x3F000000;
	ULONG saw2 = ((p+w) >> 9) | 0x3F000000;
	return punUF32(saw1) + punUF32(saw2) - 1.f;
}

// [1, 0.75, 0.5, 0.25]
// Freq precision:	32 bits
// Amp precision:	24 bits
inline float rampDownU(ULONG p){
	p = (p >> 9) | 0xbf800000;
	return punUF32(p) + 2.f;
}

inline float squareU(ULONG p){
//	phase = (phase & 0x80000000) >> 1;
//	return *(float *)&phase * 0.5f;
	return p & 0x80000000 ? 1.f : 0.f;
}

// [1, 0.5, 0, 0.5]
// Freq precision:	32 bits
// Amp precision:	24 bits
inline float triangleU(ULONG p){
	union{ float f; ULONG i; } u;
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

#define DEF(type, spec)\
template<>\
inline void print<type>(type & v, const char * post, const char * pre, FILE * fp){\
	fprintf(fp, "%s%"#spec"%s", pre, v, post);\
}
DEF(float, f) DEF(double, f) DEF(uint32_t, d) DEF(int, d)
#undef DEF

} // scl::
} // gam::

#undef PUN_F2I
#undef PUN_I2F
#undef ROUND_EPS

#include "MacroU.h"

#endif

