#ifndef GAMMA_SCL_H_INC
#define GAMMA_SCL_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information

	File Description:
	This file defines some commonly needed scalar functions.
*/

#include "Gamma/Config.h"

#if GAM_WINDOWS
	// Undefine macros in windows.h
	#ifdef max
	#undef max
	#endif
	#ifdef min
	#undef min
	#endif
#endif

#include <cmath>
#include <type_traits> // is_fundamental
#include "Gamma/Conversion.h"


namespace gam{

using std::exp;
using std::log;
using std::pow;
using std::sqrt;

/// Returns a positive length associated with argument
template<class T> double norm(const T& v);

/// Returns a positive number valid for comparing norms
template<class T> double normCompare(const T& v);


/// Scalar rank functions for numerical types
namespace scl{

#define DEF(T)\
inline T max(T v1, T v2){ return v1<v2?v2:v1; }\
inline T min(T v1, T v2){ return v1<v2?v1:v2; }
DEF(float) DEF(double)
DEF(long long) DEF(unsigned long long)
DEF(int) DEF(unsigned) 
DEF(short) DEF(unsigned short)
DEF(char) DEF(unsigned char)
#undef DEF

/// Absolute value
template<class T> T abs(T v){ return std::abs(v); }

/// Fast approximation to atan2

/// Author: Jim Shima, http://www.dspguru.com/comp.dsp/tricks/alg/fxdatan2.htm.
/// |error| < 0.01 rad
template<class T> T atan2Fast(T y, T x);

/// Returns floating point value rounded to next highest integer.
template<class T> T ceil(T val);
template<class T> T ceil(T val, T step);
template<class T> T ceil(T val, T step, T recStep);

/// Returns power of two ceiling of value.

/// This uses an algorithm devised by Sean Anderson, Sep. 2001.
/// From "Bit Twiddling Hacks", http://graphics.stanford.edu/~seander/bithacks.html.
uint32_t ceilPow2(uint32_t value);

/// Returns value clipped to [lo, hi].
template<class T> T clip(T value, T hi=T(1), T lo=T(0));

/// Returns value clipped to [-hi, hi].
template<class T> T clipS(T value, T hi=T(1));

/// Returns value clipped to [lo, hi] and signifies clipping behavior.

/// clipFlag signifies if and where clipping occured.  0 means no clipping
/// occured, -1 means clipping occured at the lower bound, and 1 means
/// clipping at the upper bound.
template<class T> T clip(T value, int & clipFlag, T hi, T lo);

/// Returns value whose magnitude is clipped to [min, max].
float clipMag(float value, float max=1.f, float min=0.f);

/// Third order polynomial approximation to first half of cosine.

/// 'u' must be in the range [0, 0.5] which corresponds to the first half
/// of the cosine.
template<class T> T cosP3(T u);

/// 8th order Taylor series approximation to a cosine.

/// 'radians' must be in [-pi, pi].
///
template<class T> T cosT8(T radians);

/// Convert decibels to amplitude
template <class T>
inline T dBToAmp(const T& db){ return pow(10, db/20.); }

/// Convert amplitude to decibels
template <class T>
inline T ampTodB(const T& amp){ return 20*log(amp); }

/// Returns an equal-loudness amplitude for a given frequency, in Hz.

/// The curve used is an A-weighting curve which approximates a 40-phon equal
/// loudness Fletcher-Munson contour.
/// \param[in] freq		frequency, in Hz
/// \param[in] maxAmp	maximum amplitude to clamp return value to
double eqLoudAmp(double freq, double maxAmp=4.);

/// Returns weights for linear fade.
template<class T> void fadeLin(T& weight1, T& weight2, T fade);

/// Returns weights for triangular window fade.

/// The weights returned are from two overlapping normalized unipolar
/// triangular windows in the fade interval [0, 2/3] and [1/3, 1].\n
/// fade	weight1		weight2 \n
/// 0.25	1			0       \n
/// 0.5		0.5			0.5		\n
/// 0.75	0			1
template<class T> void fadeTri(T& weight1, T& weight2, T fade);

/// Return Fejer weighting factor for kth harmonic in a Fourier series of length n.

/// The function is a line from (0,1) to (n,0). This is used for reducing the
/// Gibbs effect from a Fourier series summation.
template<class T> T fejer(T k, T n){ return (n-k)/n; }

/// Returns floor of floating point value.
template<class T> T floor(T val);
template<class T> T floor(T val, T step);
template<class T> T floor(T val, T step, T recStep);

/// Returns value folded into [lo, hi].

/// For out-of-range values, the boundaries act like mirrors reflecting
/// the value into the range. For an even number of periods out of the range
/// this is identical to a wrap().
template<class T> T fold(T value, T hi=T(1), T lo=T(0));
template<class T> T fold(T v, long& numFolds, T hi=T(1), T lo=T(0));

/// Returns value folded into [lo, hi] one time.
template<class T> T foldOnce(T value, T hi=T(1), T lo=T(0));

/// Returns base-2 log linearly mapped to [0, 1]

/// This function applies a log2 mapping and then normalizes the result to the
/// interval [0, 1]. Since the log(0) = -inf, the result of the log is clipped
/// below -1/recMin. The equation used is\n
///
/// v -> (log2(|v|) recMin) + 1
///
/// \param[in] v		input value
/// \param[in] recMin	negative of reciprocal of minimum log value
/// \returns log base-2 linearly mapped to [0, 1]
template<class T> T linLog2(T v, T recMin = T(1./16));

/// Returns base-2 logarithm of value.

/// If the value is not an exact power of two, the logarithm of the next
/// highest power of two will taken.
/// This uses an algorithm devised by Eric Cole, Jan. 2006.
/// From "Bit Twiddling Hacks", http://graphics.stanford.edu/~seander/bithacks.html.
uint32_t log2(uint32_t v);

/// Compile-time log2 (for powers of 2)
template <uint32_t PowerOf2> constexpr uint32_t log2();

/// Fast base-2 logarithm. For value <= 0, behavior is undefined.
float log2Fast(float v);

/// Maps value from [-1,1] to [1-depth, 1].
template<class T>
T mapDepth(T v, T depth){ return (v - T(1)) * T(0.5) * depth + T(1);  }

/// Inverse 2nd power mapping for interval [0, 1].
template<class T> T mapInvPow2(T v);

/// Computes scale and offset necessary to map value from [i0, i1] to [o0, o1].
template<class T>
void mapLin(T i0, T i1, T o0, T o1, T & scale, T & offset);

/// Linearly maps value from [i0, i1] to [o0, o1].
template<class T>
T mapLin(T v, T i0, T i1, T o0, T o1);

/// Returns v^power linearly mapped to [bound0, bound1].
double mapPower(double v, double bound1, double bound0, double power=1.);

/// Map a value in [-1,1] to a cubic approximating sin(pi/2 x)
template<class T> T mapSinSS(T v);

/// Map a value in [-1,1] to a cubic approximating 1/2 sin(pi/2 x) + 1/2
template<class T> T mapSinSU(T v);

/// Map a value in [ 0,1] to a cubic approximating -cos(pi x)
template<class T> T mapSinUS(T v);

/// Map a value in [ 0,1] to a cubic approximating -1/2 cos(pi x) + 1/2
template<class T> T mapSinUU(T v);

/// Mixes two values together (1 = thru, 0.5 = mono, 0 = swap).
template<class T> void mix2(T& io1, T& io2, T mix);

/// Returns the next representable floating-point or integer value following x in the direction of y
template<class T> T nextAfter(T x, T y);

template<class T> T pow2(T v);			///< Returns value to the 2nd power
template<class T> T pow3(T v);			///< Returns value to the 3rd power
template<class T> T pow4(T v);			///< Returns value to the 4th power
template<class T> T pow5(T v);			///< Returns value to the 5th power
template<class T> T pow6(T v);			///< Returns value to the 6th power
template<class T> T pow7(T v);			///< Returns value to the 7th power
template<class T> T pow8(T v);			///< Returns value to the 8th power
template<class T> T pow9(T v);			///< Returns value to the 9th power

/// Returns pole radius given a T60 decay length and units/sample
inline double radius60(double dcy, double ups){ return exp(M_LN001/dcy * ups); } // u/s * 1/u

/// Returns floating point value rounded to nearest integer.
template<class T> T round(T v);

/// Returns floating point value rounded to nearest integer multiple of 'step'.
template<class T> T round(T v, T step);

/// Returns floating point value rounded to nearest integer multiple of 'step'. Faster version to avoid 1/step divide.
template<class T> T round(T v, T step, T recStep);

/// Returns value rounded to nearest integer away from zero.
template<class T> T roundAway(T v);

/// Returns value rounded to nearest to nearest integer multiple of 'step' away from zero.
template<class T> T roundAway(T v, T step);

/// Returns the section 'v' lies in in [0,num] divided into 'div' sections.
inline int section(int v, int num, int divs){ return int((v*divs)/double(num)); }

/// Returns sign of a number
template<class T> int sgn(T v);

/// Fast sine function approximation

/// \param[in] x	angle in [-2, 2] corresponding to [-pi, pi]
///
template<class T> T sinFast(T x);

/// 7th order minimax polynomial approximation to sin(pi x).

/// Error is spread evenly across domain.
/// \param[in] normal	phase, in [-1, 1] (corresponding to [-pi, pi] radians)
template<class T> T sinP7(T normal);

/// 9th order minimax polynomial approximation to sin(pi x).

/// Error is spread evenly across domain.
/// \param[in] normal	phase, in [-1, 1] (corresponding to [-pi, pi] radians)
template<class T> T sinP9(T normal);

/// 7th order Taylor series approximation to a sine.

/// \param[in] radians	phase, in [-pi, pi]
///
template<class T> T sinT7(T radians);

/// 9th order Taylor series approximation to a sine.

/// \param[in] radians	phase, in [-pi, pi]
///
template<class T> T sinT9(T radians);

/// Smooth negative map

/// The return value is close to 1 if v < 0 and close to 0 if v > 0.
/// The smoothness is controlled with the bw argument.
template<class T> T smoothNeg(T v, T bw);

/// Same as smoothNeg(T,T), but 'amt' controls positive level (0,1) -> (-1,1)
template<class T> T smoothNeg(T v, T bw, T amt);

/// Smooth positive map

/// The return value is close to 1 if v > 0 and close to 0 if v < 0.
/// The smoothness is controlled with the bw argument.
template<class T> T smoothPos(T v, T bw);

template<class T> T smoothNotchFunc(T v, T bw);
template<class T> T smoothNotch(T v, T bw);

// Same as smoothNotch1, but 'amt' controls notch depth.
template<class T> T smoothNotch(T v, T bw, T amt);
template<class T> T smoothNotch2(T v, T bw);

/// Peak function.

/// When v is the output of y=x^2,
/// this formula approximates the true formula for a resonant peak
/// 1/(1 - 2rcos(theta) + r^2). The argument bw is equivalent to (1-r)^2.
/// In general, the approximation has a slightly
/// smaller bandwidth than the true response. Also, the true response is
/// periodic, while this one is not.
template<class T> T smoothPeak(T v, T bw);
template<class T> T smoothPeak1(T v, T bw);

/// Continuous sign map

/// The return value is close to 1 if v > 0 and close to -1 if v < 0.
/// 'bw' controls the width of the transition region. If 'bw' is 0, then this
/// turns into a signum function.
template<class T> T smoothSign(T v, T bw);

/// Same as smoothZero(), but takes a unary function to be applied to the value
template <class T, class F> T smoothZero(T v, T bw, F f);

/// Returns a continuous value measuring how close to zero the value is

/// The graph of this function resembles a resonant peak. The function uses the
/// square of the value for evaluation.
template<class T> T smoothZero(T v, T bw);

/// Approximate square root using log base-2 and Newton methods.

/// 'N' determines the accuracy of the approximation. For N=0, a quick and dirty
/// log base-2 approximation is performed. For N>0, N Newton iterations
/// are applied to improve the result.
template<unsigned N, class T> T sqrt(T v);

/// Approximate inverse square root using log base-2 and Newton methods.

/// 'N' determines the accuracy of the approximation. For N=0, a quick and dirty
/// log base-2 approximation is performed. For N>0, N Newton iterations
/// are applied to improve the result.
template<unsigned N, class T> T invSqrt(T v);

/// Truncates floating point value at decimal.
template<class T> T trunc(T value);

/// Truncates floating point value to step.
template<class T> T trunc(T value, T step);

/// Truncates floating point value to step. Faster version to avoid 1/step divide.
template<class T> T trunc(T value, T step, T recStep);

/// Returns multiplicaton factor for reaching -60 dB after 'samples' iterations.
double t60(double samples);

/// Returns value wrapped in [lo, hi).
template<class T> T wrap(T value, T hi=T(1), T lo=T(0));

/// Returns value wrapped in [lo, hi).

/// 'numWraps' reports how many wrappings occured where the sign, + or -,
/// signifies above 'hi' or below 'lo', respectively.
template<class T> T wrap(T value, long & numWraps, T hi=T(1), T lo=T(0));

/// Returns value incremented by 1 and wrapped into interval [0, max).
template<class T> T wrapAdd1(T v, T max){ ++v; return v == max ? 0 : v; }

/// Like wrap(), but only adds or subtracts 'hi' once from value.
template<class T> T wrapOnce(T value, T hi=T(1));

template<class T> T wrapOnce(T value, T hi, T lo);

template<class T> T wrapPhase(T radians);			///< Returns value wrapped in [-pi, pi)
template<class T> T wrapPhaseOnce(T radians);		///< Like wrapPhase(), but only wraps once




//---- ANALYSIS

/// Tests whether two values are close in value

/// 'maxULP' (maximum units in the last place) is the maximum difference allowed
/// between the values punned as integers.
/// Code taken from Bruce Dawson, "Comparing floating point numbers."
bool almostEqual(float a, float b, int maxULP=10);
bool almostEqual(double a, double b, int maxUlps=10);

/// Returns whether or not an integer value is even.
template<class T> bool even(T v);

/// Returns maximum of three values.
template<class T> T max(T v1, T v2, T v3);

/// Returns minimum of three values.
template<class T> T min(T v1, T v2, T v3);

/// Returns whether or not an integer value is odd.
template<class T> bool odd(T v);

/// Returns slope of line passing through two points.
template<class T> T slope(T x1, T y1, T x2, T y2);

/// Returns number of trailing zeros in 32-bit int

/// This implements an algorithm from the paper
/// "Using de Bruijn Sequences to Index 1 in a Computer Word"
/// by Charles E. Leiserson, Harald Prokof, and Keith H. Randall.
uint32_t trailingZeroes(uint32_t v);

/// Returns whether value is within [lo, hi].
template<class T> bool within(T v, T lo, T hi);

/// Returns whether a positive zero crossing occured.
bool zeroCrossP(float prev, float now);



//---- WAVEFORM GENERATORS

/// Returns value quantized to multiples of 2^q
uint32_t quantizePow2(uint32_t value, uint32_t q);

//---- Bipolar waveforms [-1, 1)
float rampDown	(uint32_t phase);	///< Returns value of bipolar downward ramp function.
float rampUp	(uint32_t phase);	///< Returns value of bipolar upward ramp function.
float square	(uint32_t phase);	///< Returns value of bipolar square function.
float triangle	(uint32_t phase);	///< Returns value of bipolar triangle function.

/// Returns value of fast sine function approximation
float sinFast	(uint32_t phase);

/// Returns value of sine-like function constructed from parabolic sections
float sinPara	(uint32_t phase);

/// Returns value of bipolar pulse function (rampDown() + rampUp()).
float pulse		(uint32_t phase, uint32_t width);

/// Returns value of bipolar stair function (square() + square()).
float stair		(uint32_t phase, uint32_t width);

/// Returns value of bipolar dual upward ramp function (rampUp() + rampUp()).
float rampUp2	(uint32_t phase, uint32_t width);	// rampUp + rampUp

//---- Unipolar waveforms [0, 1)
float pulseU	(uint32_t phase, uint32_t width);	///< Returns value of unipolar pulse function.
float pulseRange(uint32_t phase, uint32_t width, float highVal = 1.0, float lowVal = -1.0); ///< Returns value of ranged pulse function.
float rampUpU	(uint32_t phase);	///< Returns value of unipolar downward ramp function.
float rampUp2U	(uint32_t phase);	///< Returns value of unipolar upward ramp2 function.
float rampDownU	(uint32_t phase);	///< Returns value of unipolar upward ramp function.
float squareU	(uint32_t phase);	///< Returns value of unipolar square function.
float stairU(uint32_t phase, uint32_t width); ///< Returns value of unipolar stair function.
float triangleU	(uint32_t phase);	///< Returns value of unipolar triangle function.

template<class T> T bartlett(T nphase);								///< Bartlett window. nphase => [-1, 1)
template<class T> T blackman(T phase);								///< Blackman window function.
template<class T> T blackmanHarris(T phase);					///< Blackman-Harris window function.
template<class T> T blackmanNuttall(T phase);					///< Blackman-Nuttall window function.
template<class T> T flatTop(T phase);									///< Flat-Top window function.
template<class T> T hamming(T phase);									///< Hamming window function.
template<class T> T hann(T phase);										///< von Hann window function.
template<class T> T raisedCosine(T phase, T a, T b);	///< Raised cosine f(x) = a - b cos(x).
template<class T> T welch(T nphase);									///< Welch window function. nphase => [-1, 1)



//---- NOTE-BASED FUNCTIONS

/// Returns frequency in Hz from a 12-TET note string.

/// Notes are specified by a letter in [A, G] or [a, g], followed by an 
/// accidental '+' or '#' for a sharp, '-' or 'b' for a flat, or ' ' or nothing 
/// for a natural, and finally an integer in [0, 9] representing the octave 
/// number. For example, middle C is specified as "c4" or "c 4" and the A sharp 
/// below that as "a+3" or a#3.
double freq(const char * note);

/// Returns nearest note number within a pitch class set

/// \param[in] v			the note number to match
/// \param[in] intervals	sequence of base-36 intervals 
/// \param[in] mod			modulo amount
///
/// The sum of the values in the interval array should be equal to 'mod'.
double nearest(double v, const char * intervals="2212221", long mod=12);

/// Returns equal temperament ratio- octave^(pitch/divs)

/// \param[in] pitch	pitch class
/// \param[in] divs		number of equally tempered divisions in octave
/// \param[in] octave	base multiplier of (pseudo) octave
double ratioET(double pitch, double divs=12, double octave=2);




// Implementation ______________________________________________________________

// internal
namespace{
	inline uint32_t deBruijn(uint32_t v){
		static const unsigned char deBruijnBitPosition[32] = {
			 0,  1, 28,  2, 29, 14, 24,  3, 30, 22, 20, 15, 25, 17,  4,  8,
			31, 27, 13, 23, 21, 19, 16,  7, 26, 12, 18,  6, 11,  5, 10,  9
		};
		// Note: this is basically a hash function
		return deBruijnBitPosition[(uint32_t(v * 0x077CB531UL)) >> 27];
	}
}

//#define GEN(t, f) template<> inline t abs<t>(t v){ return f(v); }
//GEN(int, ::abs) GEN(long, labs) GEN(long long, llabs) GEN(float, fabsf) GEN(double, fabs)
////template<class T> inline T abs(T v){ return v < T(0) ? -v : v; } // only allow specializations
//#undef GEN

template<class T> T atan2Fast(T y, T x){

	T r, angle;
	T ay = scl::abs(y) + T(1e-10); // kludge to prevent 0/0 condition

	if(x < T(0)){
		r = (x + ay) / (ay - x);
		angle = T(M_3PI_4);
	}
	else{
		r = (x - ay) / (x + ay);
		angle = T(M_PI_4);
	}

	angle += (T(0.1963)*r*r - T(0.9817))*r;
	return y < T(0) ? -angle : angle;
}

template<class T> inline T ceil(T v){ return round(v + roundEps<T>()); }
template<class T> inline T ceil(T v, T s){ return ceil(v/s)*s; }
template<class T> inline T ceil(T v, T s, T r){ return ceil(v*r)*s; }

inline uint32_t ceilPow2(uint32_t v){
	v--;
	v |= v >> 1;
	v |= v >> 2;
	v |= v >> 4;
	v |= v >> 8;
	v |= v >> 16;
	return ++v;
}

template<class T> inline T clip(T v, T hi, T lo){
	     if(v < lo) return lo;
	else if(v > hi)	return hi;
	return v;
}

template<class T> inline T clip(T v, int & clipFlag, T hi, T lo){
	clipFlag = 0;
	     if(v < lo){ clipFlag = -1; return lo; }
	else if(v > hi){ clipFlag =  1; return hi; }
	return v;
}

template<class T>
inline T clipS(T v, T hi){ return clip(v, hi, -hi); }

template<class T>
inline T equals(T v1, T v2, T bw){ return equals(v1, v2, bw, pow2<T>); }

template <class T, class F>
inline T equals(T v1, T v2, T bw, F f){ return bw/(f(v1-v2) + bw); }

template<class T>
inline void fadeLin(T& w1, T& w2, T f){ w1 = T(1) - f; w2 = f; }

template<class T>
inline void fadeTri(T& w1, T& w2, T f){
	if(f < 0.25){
		w1 = f * 4;
		w2 = 0;
	}
	else if( (f >= 0.25) && (f < 0.75) ){
		fadeLin(w1, w2, T(f * 2 - 0.5));
	}
	else{
		w1 = 0;
		w2 = (1 - f) * 4;
	}
}

template<class T> inline T floor(T v){ return round(v - roundEps<T>()); }
template<class T> inline T floor(T v, T s){ return floor(v/s)*s; }
template<class T> inline T floor(T v, T s, T r){ return floor(v*r)*s; }

template<class T> inline T fold(T v, T hi, T lo){
	long t;
	return fold(v,t,hi,lo);
}

template<class T> inline T fold(T v, long& numFolds, T hi, T lo){
	long numWraps;
	v = wrap(v, numWraps, hi, lo);
	if(numWraps & 1) v = hi + lo - v;
	numFolds = numWraps;
	return v;
}

template<class T> inline T foldOnce(T v, T hi, T lo){
	if(v > hi) return hi + (hi - v);
	if(v < lo) return lo + (lo - v);
	return v;
}

template<class T> inline T linLog2(T v, T recMin){
	v = log2Fast(scl::abs(v) + T(0.000001));	// offset to avoid -inf
	return scl::max(v * recMin, T(-1)) + T(1);
}

inline uint32_t log2(uint32_t v){ return deBruijn(ceilPow2(v)); }

template<> constexpr uint32_t log2<   1>(){ return  0; }
template<> constexpr uint32_t log2<   2>(){ return  1; }
template<> constexpr uint32_t log2<   4>(){ return  2; }
template<> constexpr uint32_t log2<   8>(){ return  3; }
template<> constexpr uint32_t log2<  16>(){ return  4; }
template<> constexpr uint32_t log2<  32>(){ return  5; }
template<> constexpr uint32_t log2<  64>(){ return  6; }
template<> constexpr uint32_t log2< 128>(){ return  7; }
template<> constexpr uint32_t log2< 256>(){ return  8; }
template<> constexpr uint32_t log2< 512>(){ return  9; }
template<> constexpr uint32_t log2<1024>(){ return 10; }
template<> constexpr uint32_t log2<2048>(){ return 11; }
template<> constexpr uint32_t log2<4096>(){ return 12; }
template<> constexpr uint32_t log2<8192>(){ return 13; }

inline float log2Fast(float v){
	Twiddle<float> u(v);
	return (float)((u.i - int32_t(Expo1<float>()))) * 0.0000001192092896f;// / 8388608.f;
}

template<class T> inline T mapInvPow2(T v){ return v*(T(2)-v); }

template<class T>
inline void mapLin(T i0, T i1, T o0, T o1, T& scale, T& offset){
	scale = slope(i0, o0, i1, o1);
	offset = o0 - scale * i0;
}

template<class T>
inline T mapLin(T v, T i0, T i1, T o0, T o1){
	float scale = slope(i0, o0, i1, o1);
	return (v - i0) * scale + o0;
}

inline double mapPower(double v, double b1, double b0, double p){
	if(p != 1.) v = ::pow(v, p);
	return b0 + (b1 - b0) * v;
}

template<class T> inline T mapSinSS(T v){ return v*(T(1.5 ) - v*v*T(0.50)); }
template<class T> inline T mapSinSU(T v){ return v*(T(0.75) - v*v*T(0.25)) + T(0.5); }
template<class T> inline T mapSinUS(T v){ return v*v*(T(6) - v*T(4)) - T(1); }
template<class T> inline T mapSinUU(T v){ return v*v*(T(3) - v*T(2)); }

template<class T> inline void mix2(T& io1, T& io2, T mix){
	T t1 = (io1 - io2) * mix;
	T t2 = io1 - t1;
	io1 = t1 + io2;
	io2 = t2;
	//io1 = io1 * mix + io2 * (T(1) - mix);
	//io2 = io2 * mix + io1 * (T(1) - mix);
}

#ifdef __MSYS__
// MSYS2 does not support C++11 std::nextafter like it should
template<class T>
inline T nextAfter(T x, T y){ return x<y ? x+1 : x-1; }
template<>
inline float nextAfter(float x, float y){ return nextafterf(x,y); }
template<>
inline double nextAfter(double x, double y){ return nextafter(x,y); }
#else
template<class T>
inline T nextAfter(T x, T y){ return std::nextafter(x,y); }
#endif

template<class T> inline T pow2(T v){ return v*v; }
template<class T> inline T pow3(T v){ return v*v*v; }
template<class T> inline T pow4(T v){ return pow2(v*v); }
template<class T> inline T pow5(T v){ return pow4(v)*v; }
template<class T> inline T pow6(T v){ return pow3(pow2(v)); }
template<class T> inline T pow7(T v){ return pow6(v)*v; }
template<class T> inline T pow8(T v){ return pow4(v*v); }
template<class T> inline T pow9(T v){ return pow3(pow3(v)); }
template<class T> inline T powN(T v, int power){
	switch(power){
	case  0: return T(1);
	case  1: return v;
	case  2: return pow2(v);
	case  3: return pow3(v);
	case  4: return pow4(v);
	case  5: return pow5(v);
	case  6: return pow6(v);
	case  7: return pow7(v);
	case  8: return pow8(v);
	case  9: return pow9(v);
	case 10: return pow2(pow5(v));
	case 11: {auto vv=v*v; return pow3(vv*v)*vv;}
	case 12: return pow3(pow4(v));
	case 13: return pow3(pow4(v))*v;
	default: return std::pow(v,power);
	}
}

inline double ratioET(double pc, double divs, double ival){
	return pow(ival, pc/divs);
}

//template<class T> inline T round(T v){ return (v + roundMagic<T>()) - roundMagic<T>(); }
template<class T>
inline T round(T v){ double r=v; return T((r + roundMagic) - roundMagic); }
template<class T>
inline T round(T v, T s){ return round<double>(v/s) * s; }
template<class T>
inline T round(T v, T s, T r){ return round<T>(v * r) * s; }
template<class T>
inline T roundAway(T v){ return v<T(0) ? floor(v) : ceil(v); }
template<class T>
inline T roundAway(T v, T s){ return v<T(0) ? floor(v,s) : ceil(v,s); }

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

namespace detail{
	template <class T, bool isFund> struct RoundN;

	template <class T>
	struct RoundN<T,true>{
		T operator()(const T& v, float step, float stepRec){
			return round<T>(v,step,stepRec);
		}
	};

	template <class T>
	struct RoundN<T,false>{
		T operator()(const T& v, float step, float stepRec){
			auto res = v;
			for(auto& e : res) e = round<typename T::value_type>(e,step,stepRec);
			return res;
		}
	};
}

template <class T> T roundN(const T& v, float step, float stepRec){
	return detail::RoundN<T, std::is_fundamental<T>::value>()(v,step,stepRec);
}

template <class T> T roundN(const T& v, float step){
	return detail::RoundN<T, std::is_fundamental<T>::value>()(v,step,1.f/step);
}

template <class T> T roundN(const T& v){
	return detail::RoundN<T, std::is_fundamental<T>::value>()(v,1.f,1.f);
}

template<class T>
inline int sgn(T v){ return (T(0) < v) - (v < T(0)); }

template<class T>
inline T smoothNeg		(T v, T bw){ return T(0.5) - smoothSign(v, bw)*T(0.5); }
template<class T>
inline T smoothNeg		(T v, T bw, T a){ return a - (T(1)-a)*smoothSign(v, bw); }
template<class T>
inline T smoothPos		(T v, T bw){ return T(0.5) + smoothSign(v, bw)*T(0.5); }
template<class T>
inline T smoothNotchFunc(T v, T bw){ return v/(v+bw); }
template<class T>
inline T smoothNotch	(T v, T bw){ return smoothNotchFunc(scl::abs(v), bw); }
template<class T>
inline T smoothNotch	(T v, T bw, T amt){ return T(1) - amt*smoothPeak1(v, bw); }
template<class T>
inline T smoothNotch2	(T v, T bw){ return smoothNotchFunc(v*v, bw); }
template<class T>
inline T smoothPeak		(T v, T bw){ return bw/(v+bw); }
template<class T>
inline T smoothPeak1	(T v, T bw){ return bw/(scl::abs(v)+bw); }
template<class T>
inline T smoothSign		(T v, T bw){ return v/(scl::abs(v) + bw); }
template <class T, class F>
inline T smoothZero		(T v, T bw, F f){ return bw/(f(v) + bw); }
template<class T>
inline T smoothZero		(T v, T bw){ return smoothZero(v, bw, scl::pow2<T>); }


template<unsigned N, class T, template<class> class F>
struct NewtonIterator{
	NewtonIterator(T& v, T v0){
		F<T>(v,v0);
		NewtonIterator<N-1,T,F>(v,v0); // this just iterates
	}
};

template<class T, template<class> class F>
struct NewtonIterator<0,T,F>{
	NewtonIterator(T& /*v*/, T /*v0*/){}
};

template<class T> struct NewtonSqrtMap{
	NewtonSqrtMap(T& v, T v0){ v=T(0.5)*(v+v0/v); }
};

template<unsigned N, class T>
struct SqrtNewton : public NewtonIterator<N,T, NewtonSqrtMap> {
	SqrtNewton(T& v, T v0): NewtonIterator<N,T, NewtonSqrtMap>(v,v0)
	{}
};


template<class T> struct NewtonInvSqrtMap{
	NewtonInvSqrtMap(T& v, T v0_2){ v *= T(1.5)-v0_2*v*v; }
};

template<unsigned N, class T>
struct InvSqrtNewton : public NewtonIterator<N,T, NewtonInvSqrtMap>{
	InvSqrtNewton(T& v, T v0): NewtonIterator<N,T, NewtonInvSqrtMap>(v,v0)
	{}
};

template<class T> const Twiddle<T> invSqrtMagic();
template<> inline const Twiddle<float > invSqrtMagic(){ return Twiddle<float >(uint32_t(0x5f3759df)); }
template<> inline const Twiddle<double> invSqrtMagic(){ return Twiddle<double>(uint64_t(0x5fe6ec85e7de30daULL)); }

inline float sqrtLog2(float v){
	Twiddle<float> u(v);
	u.u=(1<<29) + (u.u>>1) - (1<<22);
	return u.f;
}

inline double sqrtLog2(double v){
	Twiddle<double> u(v);
	u.u=((uint64_t(1))<<61) + (u.u>>1) - ((uint64_t(1))<<51);
	return u.f;
}

template<unsigned N, class T> void sqrtNewton(T& v, T v0){
	SqrtNewton<N,T>(v,v0);
}

template<unsigned N, class T>
inline T sqrt(T v){
	T r=sqrtLog2(v);
	sqrtNewton<N>(r,v);
	return r;
}

template<unsigned N, class T> void invSqrtNewton(T& v, T v0){
	InvSqrtNewton<N,T>(v,v0);
}

template <class T>
inline T invSqrtLog2(T v){
	Twiddle<T> u(v);
	u.u = invSqrtMagic<T>().u - (u.u>>1);
	return u.f;
}

template<unsigned N, class T>
inline T invSqrt(T v){
	T r=invSqrtLog2(v);
	invSqrtNewton<N>(r, v*T(0.5));
	return r;
}


// sin/cos Taylor polynomial coefs
namespace{ static const double
	t73 = 42.,
	t72 = 840.,
	t71 = 1.9841269841e-04,
	t84 = 56.,
	t83 = 1680.,
	t82 = 20160.,
	t81 = 2.4801587302e-05,
	t94 = 72.,
	t93 = 3024.,
	t92 = 60480.,
	t91 = 2.7557319224e-06,
	ta5 = 90.,
	ta4 = 5040.,
	ta3 = 151200.,
	ta2 = 1814400.,
	ta1 = 2.7557319224e-07;
}

template<class T> inline T cosT8(T r){
	if(r < T(M_PI_4) && r > T(-M_PI_4)){
		T rr = r*r;
		return T(1) - rr * T(t81) * (T(t82) - rr * (T(t83) - rr * (T(t84) - rr)));
	}
	else if(r > T(0)){
		r -= T(M_PI_2);
		T rr = r*r;
		return -r * (T(1) - T(t71) * rr * (T(t72) - rr * (T(t73) - rr)));
	}
	else{
		r += T(M_PI_2);
		T rr = r*r;
		return r * (T(1) - T(t71) * rr * (T(t72) - rr * (T(t73) - rr)));
	}
}

template<class T> inline T sinT7(T r){
	if(r < T(M_PI_4) && r > T(-M_PI_4)){
		T rr = r*r;
		return r * (T(1) - T(t71) * rr * (T(t72) - rr * (T(t73) - rr)));
	}
	else if(r > T(0)){
		r -= T(M_PI_2);
		T rr = r*r;
		return T(1) - rr * T(t81) * (T(t82) - rr * (T(t83) - rr * (T(t84) - rr)));
	}
	else{
		r += T(M_PI_2);
		T rr = r*r;
		return T(-1) + rr * T(t81) * (T(t82) - rr * (T(t83) - rr * (T(t84) - rr)));
	}
}

template<class T> inline T sinT9(T r){
	if(r < T(M_PI_4) && r > T(-M_PI_4)){
		T rr = r*r;
		return r * (T(1) - T(t91) * rr * (T(t92) - rr * (T(t93) - rr * (T(t94) - rr))));
	}
	else if(r > T(0)){
		r -= T(M_PI_2);
		T rr = r*r;
		return T(1) - rr * T(ta1) * (T(ta2) - rr * (T(ta3) - rr * (T(ta4) - rr * (T(ta5) - rr))));
	}
	else{
		r += T(M_PI_2);
		T rr = r*r;
		return T(-1) + rr * T(ta1) * (T(ta2) - rr * (T(ta3) - rr * (T(ta4) - rr * (T(ta5) - rr))));
	}
}

// Input is in [-2, 2] corresponding to [-pi, pi]
template<class T> inline T sinFast(T x){
	T y = x * (T(2) - gam::scl::abs(x));
	return y * (T(0.775) + T(0.225) * gam::scl::abs(y));
}

// Input is in [0, 2] corresponding to [0, pi]
template <class T> inline T sinFast02(T x){
	T y = x * (T(2) - x);
	return y * (T(0.775) + T(0.225) * y);
}

template<class T> inline T cosP3(T n){
	return T(1) - T(32) * n * n * (T(0.75) - n);
}

template<class T> inline T sinP7(T n){
	T nn = n*n;
	return n * (T(3.138982) + nn * (T(-5.133625) + nn * (T(2.428288) - nn * T(0.433645))));
}

template<class T> inline T sinP9(T n){
	T nn = n*n;
	return n * (T(3.1415191) + nn * (T(-5.1662729) + nn * (T(2.5422065) + nn * (T(-0.5811243) + nn * T(0.0636716)))));
}


inline double t60(double samples){ return ::pow(0.001, 1./samples); }

template<class T>
inline T trunc(T v){ return round( (v > T(0)) ? v-roundEps<T>() : v+roundEps<T>() ); }

template<class T>
inline T trunc(T v, T s){ return trunc(v/s)*s; }

template<class T>
inline T trunc(T v, T s, T r){ return trunc(v*r)*s; }

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

template<class T> inline T wrap(T v, T hi, T lo){
	if(lo == hi) return lo;

	//if(v >= hi){
	if(!(v < hi)){
		T diff = hi - lo;
		v -= diff;
		if(!(v < hi)) v -= diff * (T)(unsigned)((v - lo)/diff);
	}
	else if(v < lo){
		T diff = hi - lo;
		v += diff;	// this might give diff if range is too large, so check at end of block...
		if(v < lo) v += diff * (T)(unsigned)(((lo - v)/diff) + 1);
		if(v==diff) return nextAfter(v, lo);
	}
	return v;
}

template<class T> inline T wrap(T v, long& numWraps, T hi, T lo){
	if(lo == hi){ numWraps = 0xFFFFFFFF; return lo; }

	T diff = hi - lo;
	numWraps = 0;

	if(v >= hi){
		v -= diff;
		if(v >= hi){
			numWraps = long((v - lo)/diff);
			v -= diff * T(numWraps);
		}
		++numWraps;
	}
	else if(v < lo){
		v += diff;
		if(v < lo){
			numWraps = long((v - lo)/diff) - 1;
			v -= diff * T(numWraps);
		}
		--numWraps;
	}
	return v;
}

template<class T> inline T wrapOnce(T v, T hi){
	     if(v >= hi ) return v - hi;
	else if(v < T(0)) return v + hi;
	return v;
}

template<class T> inline T wrapOnce(T v, T hi, T lo){
	     if(v >= hi) return v - hi + lo;
	else if(v <  lo) return v + hi - lo;
	return v;
}

template<class T> inline T wrapPhase(T r){
	// The result is		[r+pi - 2pi floor([r+pi] / 2pi)] - pi
	// which simplified is	r - 2pi floor([r+pi] / 2pi) .
	if(r >= T(M_PI)){
		r -= T(M_2PI);
		if(r < T(M_PI)) return r;
		return r - T(long((r+M_PI)*M_1_2PI)  )*T(M_2PI);
	}
	else if (r < T(-M_PI)){
		r += T(M_2PI);
		if(r >= T(-M_PI)) return r;
		return r - T(long((r+M_PI)*M_1_2PI)-1)*T(M_2PI);
	}
	else return r;
}

template<class T> inline T wrapPhaseOnce(T r){
	if(r >= T(M_PI))		return r - T(M_2PI);
	else if(r < T(-M_PI))	return r + T(M_2PI);
	return r;
}




template<class T> inline bool even(T v){ return 0 == odd(v); }

//template<class T> inline T max(T v1, T v2){ return v1<v2?v2:v1; }
template<class T>
inline T max(T v1, T v2, T v3){ return max(max(v1,v2),v3); }
//template<class T> inline T min(T v1, T v2){ return v1<v2?v1:v2; }
template<class T> inline T min(T v1, T v2, T v3){ return min(min(v1,v2),v3); }

template<class T> inline bool odd(T v){ return v & T(1); }

template<class T> inline T slope(T x1, T y1, T x2, T y2){ return (y2 - y1) / (x2 - x1); }

inline uint32_t trailingZeroes(uint32_t v){ return deBruijn(v & -int32_t(v)); }

template<class T> inline bool within(T v, T lo, T hi){ return !((v < lo) || (v > hi)); }

inline bool zeroCrossP(float prev, float now){
	union{ float f; int32_t i; } u1, u0;
	u1.f = prev; u0.f = now;
	return (u1.i <= 0) && (u0.i > 0);
	//(prevI <= 0) && (nowI > 0)
//	prevI += 0x00800000;
//	nowI += 0x00800000;
//	(now & ~prev)>>31
}


inline uint32_t quantizePow2(uint32_t v, uint32_t q){
	return v & (0xFFFFFFFF << q);
}


inline float pulse(uint32_t p, uint32_t w){
	// output floating point exponent should be [1, 2)
	uint32_t saw1 = Expo1<float>() | ((p-w) >> 9);
	uint32_t saw2 = Expo1<float>() | ( p    >> 9);
	return punUF(saw1) - punUF(saw2);
}

inline float pulseU(uint32_t p, uint32_t w){	
	return p < w ? 1.f : 0.f;
}

inline float pulseRange(uint32_t p, uint32_t w, float highVal, float lowVal){
	return p < w ? highVal : lowVal;
}

// [1, 0.5, 0, -0.5]
inline float rampDown(uint32_t p){
	p = Expo2<float>() | (p >> 9);
	return 3.f - punUF(p);
}

// [1, 0.75, 0.5, 0.25]
inline float rampDownU(uint32_t p){
	p = ExpoNeg1<float>() | (p >> 9);
	return punUF(p) + 2.f;
}

// [-1, -0.5, 0, 0.5]
inline float rampUp(uint32_t p){
	p = Expo2<float>() | (p >> 9);
	return punUF(p) - 3.f;
}

// [0, 0.25, 0.5, 0.75]
inline float rampUpU(uint32_t p){
	p = Expo1<float>() | (p >> 9);
	return punUF(p) - 1.f;
}

inline float rampUp2(uint32_t p, uint32_t w){
	uint32_t saw1 = (p    )>>10;
	uint32_t saw2 = (p + w)>>10;
	return punUF(Expo2<float>() | (saw1 + saw2)) - 3.f;
}

inline float rampUp2U(uint32_t p, uint32_t w){
	uint32_t saw1 = (p    )>>10;
	uint32_t saw2 = (p + w)>>10;
	return punUF(Expo1<float>() | (saw1 + saw2)) - 1.f;
}

inline float sinFast(uint32_t p){
	return sinFast(6.f - punUF(Expo4<float>() | p>>9));
}

inline float sinPara(uint32_t p){
	/*uint32_t saw = ((p)                   >> 9) | Expo4<float>(); // [4, 8]
	uint32_t tri = ((p+MaskSign<float>()) >> 9) | Expo4<float>();
	return (6.f - punUF(saw)) * abs(6.f - punUF(tri));//*/
	auto x = punUF(Expo4<float>() | (p>>9)) - 6.f;
	return x*(scl::abs(x) - 2.f);
}

// [1, 1,-1,-1]
inline float square(uint32_t p){
	// use MSB to set sign of 1.f
	p = (p & MaskSign<float>()) | Expo1<float>();
	return punUF(p);
	// branching version
	//return p & MaskSign<float>() ? -1.f : 1.f;
}

inline float squareU(uint32_t p){
	// sign shift
//	p = (~p & MaskSign<float>()) >> 1; // 2, 0, 2, 0 ...
//	return punUF(p) * 0.5f;
	// half-amp square plus 1/2
	p = (p & MaskSign<float>()) | Expo1_2<float>();
	return punUF(p) + 0.5f;
	// branching version
	//return p & MaskSign<float>() ? 0.f : 1.f;
}

inline float stair(uint32_t p, uint32_t w){
	uint32_t sqr1 = ( p    & MaskSign<float>()) | Expo1_2<float>();
	uint32_t sqr2 = ((p-w) & MaskSign<float>()) | Expo1_2<float>();
	return punUF(sqr1) - punUF(sqr2);
}

inline float stairU(uint32_t p, uint32_t w){
	return ((p & MaskSign<float>()) ? 0.5f : 0.f) + (((p+w) & MaskSign<float>()) ? 0.5f : 0.f);
}

// [1, 0, -1, 0]; abs on ramp down
inline float triangle(uint32_t p){
	p = Expo4<float>() | (p >> 9); // [4, 8]
	return scl::abs(gam::punUF(p) - 6.f) - 1.f;
}

// [1, 0.5, 0, 0.5]
inline float triangleU(uint32_t p){
	return scl::abs(rampDown(p));
}

template<class T> inline T bartlett(T n){
	return T(1) - scl::abs(n);
}

template<class T> inline T blackman(T r){
	return	T(0.42)
		+ T(0.08) * cos(T(2)*r)
		- T(0.5 ) * cos(     r); // prevents -0s
}
template<class T> inline T blackmanHarris(T r){
	return T(0.35875)
		- T(0.48829) * cos(     r)
		+ T(0.14128) * cos(T(2)*r)
		- T(0.01168) * cos(T(3)*r);
}
template<class T> inline T blackmanNuttall(T r){
	return T(0.3635819)
		- T(0.4891775) * cos(     r)
		+ T(0.1365995) * cos(T(2)*r)
		- T(0.0106411) * cos(T(3)*r);
}
template<class T> inline T nuttall(T r){
	return T(0.355768)
		- T(0.487396) * cos(     r)
		+ T(0.144232) * cos(T(2)*r)
		- T(0.012604) * cos(T(3)*r);
}
template<class T> inline T flatTop(T r){
	return T(1.0)
		- T(1.93 ) * cos(     r)
		+ T(1.29 ) * cos(T(2)*r)
		- T(0.388) * cos(T(3)*r)
		+ T(0.028) * cos(T(4)*r);
}
template<class T> inline T hamming(T r){ return raisedCosine(r, T(0.53836), T(0.46164)); }
template<class T> inline T hann(T r){ return raisedCosine(r, T(0.5), T(0.5)); }
template<class T> inline T raisedCosine(T r, T a, T b){ return a - b * cos(r); }
template<class T> inline T welch(T n){ return T(1) - n*n; }

} // scl::

template<class T> inline T magSqr(const T& v){ return v*v; }

template<class T> inline double norm(const T& v){ return scl::abs(v); }
template<class T> inline double normCompare(const T& v){ return norm(v); }

} // gam::

#endif
