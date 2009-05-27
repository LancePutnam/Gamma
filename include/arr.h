#ifndef GAMMA_ARR_H_INC
#define GAMMA_ARR_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include <math.h>
#include <string.h>
#include "Containers.h"
#include "ipl.h"
#include "mem.h"
#include "scl.h"
#include "MacroD.h"

namespace gam{

/// Array rank functions for numerical types.

/// The following are names commonly used for input parameters. \n\n
///		'src' - array that will only be read from \n
///		'dst' - array that will be written to \n
///		'arr' - array that will be read from and written to, i.e. in-place \n
///		'len' - length of array \n
///		'stride' - how many elements to step through the array by
/// \n
/// Many of these functions are numerically generic meaning that the object
/// is required to understand one or more of arithmetic operators +, -, *, \,
/// and/or relational operators <, <=, >, >=.
namespace arr{

// These functions take input, output, and/or filtered indexable objects along
// with an indexer object.
// f*  -- filtered (i/o) indexable
// i*  -- input indexable
// o*  -- output indexable
// ind -- indexer object (begin, end, stride, index)

/// f0[i] += i0[i]
template <class Tf0, class Ti0>
Tf0& add(Tf0& f0, const Ti0& i0, const Indexer& ind){
	LOOP_IND(f0[i] += i0[i];) return f0;
}

/// o0[i] = i0[i] + i1[i]
template <class To0, class Ti0, class Ti1>
To0& add(To0& o0, const Ti0& i0, const Ti1& i1, const Indexer& ind){
	LOOP_IND(o0[i] = i0[i] + i1[i];) return o0;
}

/// f0[i] = -i0[i]
template <class To0, class Ti0>
To0& invert(To0& o0, const Ti0& i0, const Indexer& ind){
	LOOP_IND(o0[i] = -i0[i];) return o0;
}

/// f0[i] *= i0[i]
template <class Tf0, class Ti0>
Tf0& mul(Tf0& f0, const Ti0& i0, const Indexer& ind){
	LOOP_IND(f0[i] *= i0[i];) return f0;
}

/// o0[i] = i0[i] * i1[i]
template <class To0, class Ti0, class Ti1>
To0 & mul(To0& o0, const Ti0& i0, const Ti1& i1, const Indexer& ind){
	LOOP_IND(o0[i] = i0[i] * i1[i];) return o0;
}

/// f0[i] -= i0[i]
template <class Tf0, class Ti0>
Tf0 & sub(Tf0& f0, const Ti0& i0, const Indexer& ind){
	LOOP_IND(f0[i] -= i0[i];) return f0;
}

/// o0[i] = i0[i] - i1[i]
template <class To0, class Ti0, class Ti1>
To0 & sub(To0& o0, const Ti0& i0, const Ti1& i1, const Indexer& ind){
	LOOP_IND(o0[i] = i0[i] - i1[i];) return o0;
}


//---- Below are old the school pointer type function prototypes
//---- TODO: Convert these to interfaced object types as above.

/// Add flipping value to even-sized array (i.e. a Nyquist amount).

/// For instance if value is 1, it will add 1, -1, 1, -1, ...
///
TEM void addFlip(T * arr, uint32_t len, T value);

/// Add then multiply array values by fixed amounts.
TEM void addMul(T * arr, uint32_t len, T add, T mul);

/// Sum elements from src into ring-buffer ring.

/// Returns the next tap index.  This will not guaranteed to be in the range [0, ringSize).
///
TEM uint32_t addToRing(T * ring, uint32_t ringSize, uint32_t ringTap, const T * src, uint32_t len);

/// Differentiate array.
TEM void differentiate(T * arr, uint32_t len, T & prev);

/// Exponentiates base using array values as exponents.
TEM void expBase(T * arr, uint32_t len, double base = 2.);

/// One-pole filter.

/// Performs equation: arr[i] = b0 * arr[i] + a1 * src[i]
///
TEM void filter1P(T * arr, const T * src, uint32_t len, T b0, T a1);

/// Integrate array.
TEM void integrate(T * arr, uint32_t len, T & prev);

/// Linearly maps array values from [i0, i1) to [o0, o1).
TEM void mapLin(T * arr, uint32_t len, T i0, T i1, T o0, T o1);

//	/// Applies mirror isometry sequence [dbqp] from first quarter of array.
//	
//	/// The sequence of mirror isometries are identity (d), reflection (b),
//	/// glide reflection (q), and rotation (p). The array should hold the first
//	/// len/4 + 1 elements of the signal.\n
//	/// Ex.: [ 1, 2, 3, x, x, x, x, x] -> [ 1, 2, 3, 2, -1,-2,-3,-2]
//	/// Ex.: [ 1, 2, x, x, x, x, x, x] -> [ 1, 2, 2, 1, -1,-2,-2,-1]
//	TEM void mirror_dbqp(T * arr, uint32_t len);

/// Applies mirror isometry sequence [dp] from first half of array.

/// The sequence of mirror isometries are identity (d) and rotation (p).
/// The first len/2 elements of the array are mirrored.\n
/// Ex.: [ 1, 2, 3, 4, x, x, x, x] -> [ 1, 2, 3, 4,-4,-3,-2,-1]
TEM void mirror_dp(T * arr, uint32_t len);

/// Applies mirror isometry sequence [dq] from first half of array.

/// The sequence of mirror isometries are identity (d) and glide relfection (q).
/// The first len/2 elements of the array are mirrored.\n
/// Ex.: [ 1, 2, 3, 4, x, x, x, x] -> [ 1, 2, 3, 4,-1,-2,-3,-4]
TEM void mirror_dq(T * arr, uint32_t len);

/// Multiply then add to array values by fixed amounts.
TEM void mulAdd(T * arr, uint32_t len, T mul, T add);

/// Multiply array by a Bartlett (triangle) window.

/// Works only for even sized arrays.
///
TEM void mulBartlett(T * arr, uint32_t len);

TEM void mulComplex(T * arrR, T * arrI, const T * srcR, const T * srcI, uint32_t len);

/// Multiply 'arr' by 'src' where 'src' is the first 'len'/2 + 1 elements
/// of a symmetric window.
TEM void mulHalfWindow(T * arr, const T * src, uint32_t len);

/// Multiply array by a line from [start, end]
TEM void mulLine(T * arr, uint32_t len, T start, T end, bool includeEnd = false);

/// Uniformly scale array values to fit in [-1, 1].

/// Returns normalization factor.
///
TEM T normalize(T * arr, uint32_t len);

TEM void overlapAdd(T * arr, const T * src, uint32_t len, uint32_t hop);

TEM void pow2(T * arr, uint32_t len);	///< Raise array values to 2nd power
TEM void pow3(T * arr, uint32_t len);	///< Raise array values to 3rd power
TEM void pow4(T * arr, uint32_t len);	///< Raise array values to 4th power

/// Weighted averaging of values between two arrays.

/// Implements formula: arr[i] = arr[i] * weight + src[i] * (1 - weight)
TEM void smooth(T * arr, const T * src, uint32_t len, T weight);

TEM void smooth(T * arr, uint32_t len, T weight, T prev=0);

//
// Filling functions
//

/// Fills array with line from [start, end) or [start, end].
TEM void line(T * dst, uint32_t len, T start, T end, bool includeEnd = false);

/// Fills array with line from [start, start + 1).
TEM void line1(T * dst, uint32_t len, T start=0);

/// Fills array with line from [start, start + (slope * len)).
TEM void lineSlope(T * dst, uint32_t len, T start, T slope);

/// Fills array with line from [start, start + len).
TEM void lineSlope1(T * dst, uint32_t len, T start = 0);


//
// Non-linear transformations
//

/// Make array values positive.
template <class Tf0>
inline void abs(Tf0& f0, const Indexer& ind){
	LOOP_IND(f0[i] = scl::abs(f0[i]); )
}

/// Clip array values to range.
template <class Tf0, class Tp>
inline void clip(Tf0& f0, const Indexer& ind, const Tp& hi, const Tp& lo){
	LOOP_IND(f0[i] = scl::clip<Tp>(f0[i], hi, lo); )
}

/// Clip array values between [-1, 1].
void clip1(float * arr, uint32_t len);

/// Mapping from linear range [-1, 1] to normalized dB range [-1, 1].
void linToDB(float * arr, uint32_t len, float minDB);

/// Perform max operation element-wise on array values with scalar value.
TEM void max(T * arr, uint32_t len, T val);

/// Perform min operation element-wise on array values with scalar value.
TEM void min(T * arr, uint32_t len, T val);

/// Round array values to integer multiples.
TEM void round(T * arr, uint32_t len);

/// Round array values to integer multiple of 'step'.
TEM void round(T * arr, uint32_t len, T step);

/// Wrap array values into range [lo, hi).
TEM void wrap(T * arr, uint32_t len, T hi=(T)1, T lo=(T)0);

/// Zero array values greater than 'threshold'.
TEM void zeroAbove(T * arr, uint32_t len, T threshold);

/// Zero array values less than 'threshold'.
TEM void zeroBelow(T * arr, uint32_t len, T threshold);



//
// Analysis
//

/// Finds elements that are within a threshold of their nearest neighbors.

/// @param[in]  src			Source array of elements
/// @param[in]  indices		Index array used for iterating the source array.
/// @param[out] indices		Cluster indices.
/// @param[in]  numIndices	Number of source indices.
/// @param[out] numIndices	Number of cluster indices.
/// @param[in]  threshold	Magnitude threshold of cluster.
TEM void cluster(const T * src, uint32_t * indices, uint32_t & numIndices, T threshold);

void compact(float * dst, const float * src, uint32_t len, uint32_t chunkSize);

/// Measures distances from 2-d points to a single source point.
void distance2(float * dst, const float * xFrom, const float * yFrom, uint32_t len, float xTo, float yTo);

/// Returns dot-product of two arrays.
TEM T dot(const T * src1, const T * src2, uint32_t len);

/// Returns dot-product of two arrays of length 4.
TEM T dot4(const T * src1, const T * src2);

/// Returns energy of array values.

/// Energy is the sum of values squared.
///
TEM T energy(const T * src, uint32_t len);

/// Get indices of min and max values.
TEM void extrema(const T * src, uint32_t len, uint32_t & indexMin, uint32_t & indexMax);

/// Perform linear least squares fitting of array.

/// The independent axis is the array indices, i.  The best fit line
/// equation is y = inter + slope * i.
template <class T1, class T2, class T3>
void fitLine(const T1 * src, uint32_t len, T2& slope, T3& inter);

/// Estimate the fundamental frequency of a spectrum using HPS.

/// Returns index of detected fundamental.
///
uint32_t fundHPS(float * tmp, const float * mag, uint32_t len, uint32_t downSample=4);

//	uint32_t fundLowPeak


/// Compute histogram of 'src'.

/// Values in 'src' are tallied and placed in 'bins', where the index of the
/// bin is the integer part of the source values.  Source values greater than 
/// the number of bins are ignored.  The scale and offset parameters can be
/// used to put the src values into the proper range.
template <class Ts, class Tb>
void histogram(const Ts * src, uint32_t len, Tb * bins, uint32_t numBins, Ts scale=1);

template <class Ts, class Tb>
void histogram(const Ts * src, uint32_t len, Tb * bins, uint32_t numBins, Ts scale, Ts offset);

/// Computes harmonic product spectrum.

/// 'downSample' specifies how many times to downsample and multiply
///
void hps(float * dst, const float * src, uint32_t len, uint32_t downSample);

/// Returns index of maximum value.
TEM uint32_t max(const T * src, uint32_t len);

/// Returns index of maximum absolute value (i.e. magnitude).
TEM uint32_t maxAbs(const T * src, uint32_t len);

/// Locates local maxima and writes their indices into 'dst'.

///	Returns number of maxima found.
///
TEM uint32_t maxima(const T * src, uint32_t len, uint32_t * dst);

/// Returns the mean (average) value of array values.
TEM T mean(const T * src, uint32_t len);

/// Returns the mean absolute value of array values.
TEM T meanAbs(const T * src, uint32_t len);

/// Returns mean absolute difference of array values.
TEM T meanAbsDiff(const T * src, uint32_t len);

/// Returns weighted mean of array values.

/// Weights must be positive.
///
TEM T meanWeighted(const T * src, T * weights, uint32_t len);

/// Returns weighted mean in [0, len) of indices of weights.

/// Weights must be positive.
///		Can be used to compute centroid of spectrum.
TEM T meanWeightedIndex(const T * weights, uint32_t len);

/// Returns index of minimum value in array.
TEM uint32_t min(const T * src, uint32_t len);

TEM void minimaRemove(const T * src, uint32_t * indices, uint32_t & numIndices);

/// Returns norm of array values.

/// The norm is the square root of energy.
///
TEM T norm(const T * src, uint32_t len);

/// Returns taxicab norm of array values (sum of absolute values).
TEM T normTaxi(const T * src, uint32_t len);

/// Returns unnormalized Nyquist value for use with DFT.
TEM T nyquist(const T * src, uint32_t len);

/// Returns root mean square- the normalized norm.
TEM T rms(const T * src, uint32_t len);

/// Returns index of absolute maximum slope in array.
TEM uint32_t slopeAbsMax(const T * src, uint32_t len);

/// Returns index of maximum slope in array.
TEM uint32_t slopeMax(const T * src, uint32_t len);

/// Insertion sort of elements.

/// Elements are sorted from lowest to highest.
/// This sort is fastest for small length arrays and mostly sorted sets.
TEM void sortInsertion(T * arr, uint32_t len);

/// Insertion sort of indexed elements.

/// Elements are sorted from lowest to highest.
/// This sort is fastest for small length arrays and mostly sorted sets.
TEM void sortInsertion(const T * src, uint32_t * indices, uint32_t numIndices);

/// Quick sort of elements.

/// Elements are sorted from lowest to highest.
///
TEM void sortQuick(const T * src, uint32_t * indices, long beg, long end);

/// Returns sum of values in array (e.g. DC amount).
TEM T sum(const T * src, uint32_t len);

/// Variance (deviation from mean).
TEM T variance(const T * src, uint32_t len);

/// Returns number of values within [-threshold, theshold].
TEM uint32_t within(const T * src, uint32_t len, T threshold);

/// Returns number of values within [lo, hi].
TEM uint32_t within(const T * src, uint32_t len, T lo, T hi);

/// Returns number of values that equal zero.
TEM uint32_t zeroCount(const T * src, uint32_t len);

/// Returns number of zero-crossings in array.

/// 'prev' is the last value from the previous buffer.
///
uint32_t zeroCross(const float * src, uint32_t len, float prev);

TEM void zeroCross(const T * src, uint32_t len, uint32_t & nzc, uint32_t & pzc);

/// Returns index of first zero-crossing or 0 if none detected.
uint32_t zeroCrossFirst(const float * src, uint32_t len);

TEM uint32_t zeroCrossMax(const T * src, uint32_t len);

/// Returns # of negative slope zero-crossings.

/// 'prev' is the last value from the previous buffer.
///
uint32_t zeroCrossN(const float * src, uint32_t len, float prev);


//
// Conversion
//

/// Generates tables for fast conversion methods.
void conversionInit();

/// Sets indices [numIndices, maxNumIndices) to complement indices.

/// Indices must be sorted from low to high.
///
void indicesComplement(uint32_t * indices, uint32_t numIndices, uint32_t maxNumIndices);

/// In-place magnitude-frequency to polar conversion.
void magFrqToPolar(float * frq, float * phsAccum, uint32_t len, float factorUnwrap);

/// Compute frequencies based on phase differences (in-place).

/// Upon completion, p1 holds the current phases and p0 holds the computed 
/// frequencies.\n\n
/// The basic algorithm is:\n
/// p0[i] = p1[i]\n
/// p1[i] = freq\n
TEM void phaseToFreq(T * p0, T * p1, uint32_t len, T ups);

/// Compute frequencies based on phase differences.

/// @param[out]	frq		frequency values
/// @param[in]	phs0	current phases
/// @param[in]	phs1	previous phases
/// @param[in]	len		length of arrays
/// @param[in]	ups		number of units between phase samplings
TEM void phaseToFreq(T * frq, const T * phs0, const T * phs1, uint32_t len, T ups);

/// In-place polar to magnitude-frequency conversion.
void polarToMagFrq(float * phs, float * phsPrev, uint32_t len, float factorWrap, float fundFreq, float fundRadians);

/// In-place polar to rectangular conversion.
void polarToRect(float * mag, float * phs, uint32_t len);

/// Fast in-place polar to rectangular conversion.

/// Call conversionInit() before using.
///
void polarToRectFast(float * real, float * imag, uint32_t len);

/// In-place rectangular to polar conversion.
void rectToPolar(float * real, float * imag, uint32_t len);

/// Fast in-place rectangular to polar conversion.

/// Call conversionInit() before using.
///
void rectToPolarFast(float * real, float * imag, uint32_t len);


//
// Printing
//

void print(const float * src, uint32_t len);
void print(const float * src1, const float * src2, uint32_t len);
void printHex(const float * src, uint32_t len);

namespace{
	
	const uint32_t LUTSize = 2048 - 1; // power of two - 1
	const uint32_t LUTMask = LUTSize;
	const uint32_t LUTSize2 = LUTSize >> 1;
	const uint32_t LUTSize4 = LUTSize >> 2;
	const float LUTSize2F = (float)LUTSize2;
	const float phaseToIndex = ((float)LUTSize) * M_1_2PI;	

	float * atanLUT = 0;
	float * magLUT = 0;
	ArrayPow2<float> * sinLUT = 0;
	
//	ArrayPow2<float> sinLUT(11);
//	float atanLUT[LUTSize + 1], magLUT[LUTSize + 1];

//	ArrayPow2<float> * sinLUT;
//	float * atanLUT, * magLUT;
//	const uint32_t LUTSize;
//	
//	const uint32_t LUTMask;
//	const uint32_t LUTSize2;
//	const uint32_t LUTSize4;
//	const float LUTSize2F;
//	const float phaseToIndex;
}



// Implementation_______________________________________________________________

TEM inline void addFlip(T * arr, uint32_t len, T value){
	len >>= 1;
	LOOP_P(len,
		*arr++ += value;
		*arr++ -= value;
	)
}

TEM inline void addMul(T * arr, uint32_t len, T add, T mul){
	LOOP_P(len, *arr = (*arr + add) * mul; arr++; )
}

TEM uint32_t addToRing(T * ring, uint32_t ringSize, uint32_t ringTap, const T * src, uint32_t len){
	uint32_t endTap = ringTap + len;

	if(endTap <= ringSize){		// haven't gone past end
		add(ring + ringTap, src, len);
	}
	else{						// have gone past end, do wrapped addition
		uint32_t samplesUnder	= ringSize - ringTap;
		uint32_t samplesOver	= endTap - ringSize;
		add(ring + ringTap, src, samplesUnder);
		add(ring, src + samplesUnder, samplesOver);
	}
	return endTap;
}

TEM inline void smooth(T * arr, const T * src, uint32_t len, T weight){
	LOOP(len, *arr = ipl::linear(weight, *src++, *arr); arr++; )
}

TEM inline void smooth(T * arr, uint32_t len, T weight, T prevValue){
	T prev = prevValue;
	LOOP(len,
		T curr = ipl::linear(weight, *arr, prev);
		*arr++ = curr;
		prev = curr;
	)
}

TEM inline void differentiate(T * arr, uint32_t len, T & prev){
	LOOP(len,
		T curr = *arr;
		*arr++ = curr - prev;
		prev = curr;
	)
}

TEM inline void expBase(T * arr, uint32_t len, double base){
	LOOP(len, arr[i] = ::pow(base, arr[i]);)
}

TEM inline void filter1P(T * arr, const T * src, uint32_t len, T b0, T a1){
	LOOP_P(len, *arr = scl::dot2(*arr, *src, b0, a1); arr++; src++; )
}

TEM inline void integrate(T * arr, uint32_t len, T & prev){
	LOOP_P(len,
		T curr = *arr;
		*arr++ = curr + prev;
		prev = curr;
	)
}

TEM inline void mapLin(T * arr, uint32_t len, T i0, T i1, T o0, T o1){
	LOOP(len, arr[i] = scl::mapLin(arr[i], i0, i1, o0, o1);)
}

//TEM inline void mirror_dbqp(T * arr, uint32_t len){
//	
//	T * arr3 = arr + (len>>1);	// 3rd quad start
//	T * arr2 = arr3 - 2;		// 2nd quad end
//	T * arr4 = arr + len - 2;	// 4th quad end
//	
//	//for(uint32_t j=1; j<=(len>>2); ++j){
//	for(uint32_t j=0; j<len; j+=4){
//		float val = *arr++;
//		*arr2-- =  val;
//		*arr3++ = -val;
//		*arr4-- = -val;
//	}
//}

TEM inline void mirror_dp(T * arr, uint32_t len){
	T * arr2 = arr + len - 1;	// 2nd half end
	LOOP_S(len, 2, *arr2-- = -*arr++; )
}

TEM inline void mirror_dq(T * arr, uint32_t len){
	T * arr2 = arr + (len>>1);	// 2nd half begin
	LOOP_S(len, 2, *arr2++ = -*arr++; )
}

TEM inline void mulAdd(T * arr, uint32_t len, T mul, T add){
	LOOP_P(len, *arr = *arr * mul + add; arr++; )
}

TEM inline void mulBartlett(T * arr, uint32_t len){
	T * end = arr + len - 1;
	uint32_t len_2 = len >> 1;
	const T slope = (T)1. / (T)len_2;
	T line = slope;
	
	*arr++ = (T)0;
	LOOP(len_2 - 1,
		*arr++ *= line;
		*end-- *= line;
		line += slope;
	)
}

TEM inline void mulComplex(T * arrR, T * arrI, const T * srcR, const T * srcI, uint32_t len){
	LOOP_P(len, scl::mulComplex(*arrR++, *arrI++, *srcR++, *srcI++); )
}

TEM inline void mulHalfWindow(T * arr, const T * src, uint32_t len){
	T * end = arr + len - 1;
	len >>= 1;
	
	*arr++ *= *src++;
	LOOP(len - 1,
		const T val = *src++;
		*arr++ *= val;
		*end-- *= val;
	)
	*arr *= *src;
}

TEM inline void mulLine(T * arr, uint32_t len, T start, T end, bool includeEnd){
	T slope = (end - start)/(T)(len - (includeEnd ? 1 : 0));
	LOOP_P(len, *arr++ *= start; start += slope;)
}

TEM T normalize(T * arr, uint32_t len){
	T max = scl::abs(arr[maxAbs(arr, len)]);
	T normFactor = T(1)/max;
	if(max != T(0)){ mul(arr, gen::Val<T>(normFactor), len); }
	return normFactor;
}

TEM inline void overlapAdd(T * arr, const T * src, uint32_t len, uint32_t hop){
	uint32_t lap = len - hop;
	add(arr, arr + hop, src, lap);
	mem::copy(arr + lap, src + lap, hop);
}

TEM inline void pow2(T * arr, uint32_t len){
	LOOP_P(len, *arr *= *arr; arr++; )
}

TEM inline void pow3(T * arr, uint32_t len){
	LOOP_P(len, *arr = scl::pow3(*arr); arr++; )
}

TEM inline void pow4(T * arr, uint32_t len){
	LOOP_P(len, *arr = scl::pow4(*arr); arr++; )
}

TEM inline void line(T * dst, uint32_t len, T start, T end, bool includeEnd){
	T slope = (end - start)/(T)(len - (includeEnd ? 1 : 0));
	LOOP_P(len, *dst++ = start; start += slope;)
}

TEM inline void line1(T * dst, uint32_t len, T start){
	T slope = 1/(T)len;
	LOOP_P(len, *dst++ = start; start += slope;)
}

TEM inline void lineSlope(T * dst, uint32_t len, T start, T slope){
	LOOP_P(len, *dst++ = start; start += slope;)
}

TEM inline void lineSlope1(T * dst, uint32_t len, T start){
	T val = start;
	LOOP_P(len, *dst++ = val; val++;)
}

inline void clip1(float * arr, uint32_t len){
	uint32_t * arrUI = (uint32_t *)arr;
	LOOP(len,
		uint32_t val = *arrUI;
		uint32_t sign = val & 0x80000000;
		val &= 0x7fffffff;	// mask off the sign to catch infs and NaNs
		uint32_t above = (val + 0x00800000) >> 30;
		
		if(above) val = 0x3f800000;
		
		//val = -above | 
		
		*arrUI++ = val | sign; 
	)
}

TEM inline void max(T * arr, uint32_t len, T val){
	LOOP_P(len, *arr = scl::max(*arr, val); arr++; )
}

TEM inline void min(T * arr, uint32_t len, T val){
	LOOP_P(len, *arr = scl::min(*arr, val); arr++; )
}

TEM inline void round(T * arr, uint32_t len){
	LOOP_P(len, *arr = scl::round(*arr); arr++; )
}

TEM inline void round(T * arr, uint32_t len, T step){
	T stepRec = (T)1 / step;
	LOOP_P(len, *arr = scl::round(*arr, step, stepRec); arr++; )
}

TEM inline void wrap(T * arr, uint32_t len, T hi, T lo){
	LOOP_P(len, *arr = scl::wrap(*arr, hi, lo); arr++; )
}

TEM inline void zeroAbove(T * arr, uint32_t len, T threshold){
	LOOP_P(len, if(*arr > threshold){ *arr = (T)0; } arr++; )
}

TEM inline void zeroBelow(T * arr, uint32_t len, T threshold){
	LOOP_P(len, if(*arr < threshold){ *arr = (T)0; } arr++; )
}


TEM void cluster(const T * src, uint32_t * indices, uint32_t & numIndices, T threshold){
	
	if(numIndices == 0) return;
	
	uint32_t * newIndices = indices;
	uint32_t newNumIndices = 0;
	
	T prev = src[*indices++];
	bool inCluster = false;
	
	LOOP(numIndices - 1,
		uint32_t index = *indices++;
		T curr = src[index];
	
		if( fabs(curr - prev) <= threshold ){
			
			if(!inCluster){
				// Add previous index
				*newIndices++ = indices[-2];
				newNumIndices++;
			}
			
			// Add current index
			*newIndices++ = index;
			newNumIndices++;
			
			inCluster = true;
		}
		else{
			inCluster = false;
		}
		
		prev = curr;
	)
	
	numIndices = newNumIndices;
}

inline void distance2(float * dst, const float * xFrom, const float * yFrom, uint32_t len, float xTo, float yTo){
	LOOP_P(len, *dst++ = hypot(*xFrom - xTo, *yFrom - yTo); xFrom++; yFrom++; )
}

TEM inline T dot(const T * src1, const T * src2, uint32_t len){
	T sum = (T)0;
	LOOP_P(len, sum += *src1 * *src2; src1++; src2++; )
	return sum;
}

TEM inline T dot4(const T * src1, const T * src2){
	return	src1[0] * src2[0]
		+	src1[1] * src2[1]
		+	src1[2] * src2[2]
		+	src1[3] * src2[3];
}

TEM inline T energy(const T * src, uint32_t len){
	T nrg = (T)0;
	LOOP_P(len, nrg += *src * *src; src++;)
	return nrg;
}

TEM void extrema(const T * src, uint32_t len, uint32_t & indexMin, uint32_t & indexMax){
	indexMin = 0;
	indexMax = 0;
	
	T min = *src++;
	T max = min;
	
	for(uint32_t i = 1; i < len; i++){
		T val = *src++;
		if(val > max){
			max = val;
			indexMax = i;
		}
		else if(val < min){
			min = val;
			indexMin = i;
		}
	}
}

/*
	slope = cov(x, y) / var(x)
	inter = x_mean - slope * y_mean
*/

template <class T, class T2, class T3>
inline void fitLine(const T * src, uint32_t len, T2& slope, T3& inter){
	T lenT  = T(len);
	T meanX = (lenT - T(1)) * T(0.5);	// mean of independent variables (the indices)
	T meanY = sum(src, len) / lenT;		// mean of dependent variables
	T varX  = T(2)*scl::sumOfSquares(meanX); // variance of x

	T cov  = T(0);	// covariance
	T dx   =-meanX;	// current distance from point to center along x

 	LOOP_P(len,
		cov += dx++ * (*src++ - meanY);
 	)

	slope = cov / varX;
	inter = meanY - slope * meanX;
}


/*

Histogram using multiplicity characteristic of a multiset

indices.set(buf, len);
arr::add(bins, indices, 1);

void Indices::set(const T * src, uint len){
	clear();
	for(uint i=0; i<len; ++i){
		long ind = (long)src[i];
		if(ind >= 0) (*this) << (uint)ind;
	}
}

void add(T * arr, Indices & ind, T val){
	for(uint i=ind.begin(); i<ind.end(); i+=ind.stride()) arr[i] += val;
}

*/

template <class Ts, class Tb>
inline void histogram(const Ts * src, uint32_t len, Tb * bins, uint32_t numBins, Ts scale){
	LOOP_P(len,
		int32_t index = (int32_t)(*src++ * scale);
		if(index >= 0 && index < (int32_t)numBins) bins[index]++;
	)
}

template <class Ts, class Tb>
inline void histogram(const Ts * src, uint32_t len, Tb * bins, uint32_t numBins, Ts scale, Ts offset){
	LOOP_P(len,
		int32_t index = (int32_t)(*src++ * scale + offset);
		if(index >= 0 && index < (int32_t)numBins) bins[index]++;
	)
}

TEM uint32_t max(const T * src, uint32_t len){
	uint32_t index = 0;
	T max = src[0];

	for(uint32_t i=1; i<len; i++){
		T val = src[i];
		if(val > max){
			max = val;
			index = i;
		}
	}
	return index;
}

TEM uint32_t maxAbs(const T * src, uint32_t len){
	uint32_t index = 0;
	T max = scl::abs(*src++);

	for(uint32_t i=1; i<len; i++){
		T val = scl::abs(*src++);
		if(val > max){
			max = val;
			index = i;
		}
	}
	return index;
}

TEM uint32_t maxima(const T * src, uint32_t len, uint32_t * dst){
	
	T prev = *src++;
	T curr = *src++;
	
	uint32_t numPeaks = 0;
	
	for(uint32_t i=1; i<len-1; ++i){
		T next = *src++;

		if(curr > prev && curr > next){
			*dst++ = i;
			numPeaks++;
			i++;
			prev = next;
			curr = *src++;
		}
		else{
			prev = curr;
			curr = next;
		}
	}
	return numPeaks;
}

TEM inline T mean(const T * src, uint32_t len){
	return sum(src, len) / (T)len;
}

TEM inline T meanAbs(const T * src, uint32_t len){
	T mean = (T)0;
	T rec = (T)1 / (T)len;
	LOOP_P(len, mean += scl::abs(*src++); )
	return mean * rec;
}

TEM T meanAbsDiff(const T * src, uint32_t len){
	T sum = (T)0;
	T mean = (T)0;
	
	T prev = *src++;
    LOOP(len - 1,
		T now = *src++;
		T diff = now - prev;
		T abs = scl::abs(now);
		sum += scl::abs(diff);
		mean += abs;
		prev = now;
    )
	//if(mean < 0.0001f) mean = 0.0001f;
	return mean < T(0.25) ? T(-1) : sum * T(0.5) / mean;

}

TEM T meanWeighted(const T * src, T * weights, uint32_t len){
	T mean = (T)0;

	// One loop: faster, but less accurate
//	LOOP(len,
//		mean += *src++ * *weight;
//		normFactor += *weight++;
//	)
//	return mean /= normFactor;

	// Two loops: slower, but avoids overflow
	T normFactor = (T)1 / sum(weights, len);
	LOOP(len, mean += *src++ * *weights++ * normFactor;)
	return mean;
}

TEM T meanWeightedIndex(const T * weights, uint32_t len){
	T mean = (T)0;
	T normFactor = sum(weights, len);
	if(normFactor == (T)0) return (T)0;
	normFactor = (T)1 / normFactor;
	T weightFactor = normFactor;
	LOOP(len,
		mean += *weights++ * weightFactor;
		weightFactor += normFactor;
	)
	return mean;
}


TEM uint32_t min(const T * src, uint32_t len){
	uint32_t index = 0;
	T min = src[0];

	for(uint32_t i=1; i<len; i++){
		T val = src[i];
		if(val < min){
			min = val;
			index = i;
		}
	}
	return index;
}

TEM void minimaRemove(const T * src, uint32_t * indices, uint32_t & numIndices){

	if(numIndices < 3) return;

	uint32_t * newIndices = indices + 1;	// first index is never a minima
	uint32_t newNumIndices = 2;			// always keep first and last indices

	T prev = src[*indices++];
	T curr = src[*indices++];

	for(uint32_t i=1; i<numIndices-1; i++){
		uint32_t index = *indices++;
		T next = src[index];

		if(curr >= prev || curr >= next){		// it's not a minima
			*newIndices++ = indices[-2];
			newNumIndices++;
		}
		
		prev = curr;
		curr = next;
	}

	*newIndices = *--indices;		// add last index
	numIndices = newNumIndices;
}

//TEM void minimaRemove(const T * src, uint32_t * indices, uint32_t & numIndices){
//
//	if(numIndices < 3) return;
//
//	uint32_t * newIndices = indices + 1;	// first index is never a minima
//	uint32_t newNumIndices = 2;			// always keep first and last indices
//
//	T prev = src[*indices++];
//	T curr = src[*indices++];
//
//	numIndices -= 2;
//
//	for(uint32_t i=0; i<numIndices; i++){
//		uint32_t index = *indices++;
//		T next = src[index];
//
//		if(curr < prev && curr < next){		// it's a minima
//			if(++i == numIndices) break;
//			index = *indices++;
//			*newIndices++ = index;
//			newNumIndices++;			
//			prev = next;
//			curr = src[index];
//		}
//		else{
//			*newIndices++ = index;
//			newNumIndices++;
//			prev = curr;
//			curr = next;
//		}
//	}
//
//	*newIndices = *--indices;		// add last index
//	numIndices = newNumIndices;
//}

TEM inline T norm(const T * src, uint32_t len){
	return (T)sqrt((double)energy(src, len));
}

TEM T normTaxi(const T * src, uint32_t len){
	T v = (T)0;
	LOOP(len, v += scl::abs(src[i]);)
	return v;
}

TEM inline T nyquist(const T * src, uint32_t len){
	T sum = (T)0;
	LOOP_S(len, 2, sum += src[i] - src[i+1];)
	return sum;
}

TEM inline T rms(const T * src, uint32_t len){
	return norm(src, len) / (T)len;	
}

TEM uint32_t slopeAbsMax(const T * src, uint32_t len){
	uint32_t index = 0;
	T prev = *src++;
	T maxSlope = (T)0;

	LOOP(len-1,
		T now = *src++;
		T slope = scl::abs(now - prev);
		if(slope > maxSlope){
			maxSlope = slope;
			index = i;
		}
		prev = now;
	)
	return index;
}

TEM uint32_t slopeMax(const T * src, uint32_t len){
	uint32_t index = 0;
	T prev = *src++;
	T maxSlope = (T)0;

	for(uint32_t i=0; i<len-1; ++i){
		T now = *src++;
		T slope = now - prev;
		if(slope > maxSlope){
			maxSlope = slope;
			index = i;
		}
		prev = now;
	}
	return index;
}

//TEM void sort3(T * arr){
//
//	T v1 = *arr;
//	T v2 = arr[1];
//	T v3 = arr[2];
//	
//	
//
//}

TEM void sortInsertion(T * arr, uint32_t len){                                              
	for(uint32_t i = 1; i < len; i++){
		T val = arr[i];
		uint32_t j = i - 1;
		for(; (j < len) && (arr[j] > val); j--){
			arr[j + 1] = arr[j];                                      
		}
		arr[j + 1] = val;
	} 
}

TEM void sortInsertion(const T * src, uint32_t * indices, uint32_t numIndices){                                              
	for(uint32_t i = 1; i < numIndices; i++)
	{
		uint32_t index = indices[i];
		T val = src[index];
		uint32_t j = i - 1;
		while((j < numIndices) && (val < src[indices[j]]))
		{
			indices[j + 1] = indices[j]; 
			j--;                                         
		}
		indices[j + 1] = index;
	} 
}

TEM void sortQuick(const T * src, uint32_t * indices, long beg, long end){
	// must be at least 1 element to sort
	if(end > beg + 1) {
		long piv = indices[beg], l = beg + 1, r = end;
		T pivVal = src[piv];
		while (l < r){
			if ( src[indices[l]] <= pivVal ) 
				l++;
			else 
				mem::swap(indices[l], indices[--r]);
		}
		mem::swap(indices[--l], indices[beg]);
		sortQuick(src, indices, beg, l);
		sortQuick(src, indices, r, end);
	}
}

TEM inline T sum(const T * src, uint32_t len){
	T result = (T)0;
	LOOP_P(len, result += *src++;)
	return result;
}

TEM inline T variance(const T * src, uint32_t len){
	T rec = (T)1/(T)len;
	T avg = sum(src, len) * rec;
	T var = (T)0;
	
	LOOP_P(len,
		T diff = *src++ - avg;
		var += diff * diff;
	)
	return var * rec;
}

TEM uint32_t within(const T * src, uint32_t len, T threshold){
	uint32_t count = 0;
	LOOP(len,
		T val = src[i];
		if((val <= threshold) && (val >= -threshold)) count++;
	)
	return count;
}

TEM uint32_t within(const T * src, uint32_t len, T lo, T hi){
	uint32_t count = 0;
	LOOP(len,
		T val = src[i];
		if((val <= hi) && (val >= lo)) count++;
	)
	return count;
}

TEM inline uint32_t zeroCount(const T * src, uint32_t len){
	uint32_t zeros = 0;
	LOOP_P(len, if(*src == (T)0){ zeros++; } src++; )
	return zeros;
}

TEM void zeroCross(const T * src, uint32_t len, uint32_t& nzc, uint32_t& pzc){
	pzc = 0;
	nzc = 0;
	T prev = *src++;
	LOOP(len - 1,
		T curr = src[i];
		if		(curr > (T)0 && prev <= (T)0) pzc++;
		else if	(curr <= (T)0 && prev > (T)0) nzc++;
		prev = curr;
	)
}

TEM uint32_t zeroCrossMax(const T * src, uint32_t len){
	T prev = *src++;
	T max = 0;
	uint32_t ind = 0;
	
	LOOP(len - 1,
		T curr = src[i];
		if(curr > (T)0 && prev <= (T)0){
			T slope = curr - prev;
			if(slope > max){
				max = slope;
				ind = i;
			}
		}
		prev = curr;
	)
	return ind;
}

/*

7243....
23470156

*/
inline void indicesComplement(uint32_t * indices, uint32_t numIndices, uint32_t maxNumIndices){
	uint32_t * comp = indices + numIndices;
	for(uint32_t i=0; i<maxNumIndices; i++){
		if(*indices == i)	indices++;
		else				*comp++ = i;
	}
}

TEM void phaseToFreq(T * p0, T * p1, uint32_t len, T ups){
	T factor = (T)1 / (M_2PI * ups);
	LOOP_P(len,
		T dp = scl::wrapPhase(*p0 - *p1);		// wrap phase into [-pi, pi)
		*p1++ = *p0;							// prev phase = curr phase
		*p0++ = dp * factor;
	)
}

TEM void phaseToFreq(T * f, const T * p0, const T * p1, uint32_t len, T ups){
	T factor = (T)1 / (M_2PI * ups);
	LOOP_P(len,
		T dp = scl::wrapPhase(*p0 - *p1);		// wrap phase into [-pi, pi)
		*f++ = dp * factor;
		p0++; p1++;
	)
//	subtract(f, p0, p1, len);
//	mul(f, len, factor);
//	wrap(f, len, -0.5 / ups, 0.5 / ups);
}


/*
00		00		0	0
01		10		1	2
10		01		2	1
11		11		3	3

000		000		0	0	
001		100		1	4	x
010		010		2	2	x
011		110		3	6	x
100		001		4	1	
101		101		5	5	x
110		011		6	3	
111		111		7	7

0000	0000	0	0	
0001	1000	1	8	x
0010	0100	2	4	x
0011	1100	3	12	x
0100	0010	4	2	
0101	1010	5	10	x
0110	0110	6	6	x
0111	1110	7	14	x
1000	0001	8	1	
1001	1001	9	9	x
1010	0101	10	5
1011	1101	11	13	x
1100	0011	12	3
1101	1011	13	11	
1110	0111	14	7
1111	1111	15	15
*/
// In-place bit reversed reordering of array.
//unsigned int bitReverse(float * arr, uint32_t len, uint32_t numBits){
//	for(uint32_t i=1; i<len; i++){
//	
//	}
//}

} // arr::

} // gam::

#include "MacroU.h"

#endif

