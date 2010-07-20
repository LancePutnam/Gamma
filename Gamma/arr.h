#ifndef GAMMA_ARR_H_INC
#define GAMMA_ARR_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information

	File Description:
	Functions for processing arrays of data.
*/

#include <math.h>
#include <string.h>
#include "Gamma/ipl.h"
#include "Gamma/mem.h"
#include "Gamma/scl.h"
#include "Gamma/Containers.h"

#define LOOP(n,s) for(uint32_t i=0; i<n; i+=s)
#define TEM template <class T>

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

/// Add source array to destination array
TEM void add(T * dst, const T * src, uint32_t len, uint32_t str=1){
	LOOP(len,str){ dst[i] += src[i]; }
}

/// Add two source arrays into destination array
TEM void add(T * dst, const T * src1, const T * src2, uint32_t len, uint32_t str=1){
	LOOP(len,str){ dst[i] = src1[i] + src2[i]; }
}

/// Sum elements from src into ring-buffer ring.

/// Returns the next tap index.  This will not guaranteed to be in the range [0, ringSize).
///
TEM uint32_t addToRing(T * ring, uint32_t ringSize, uint32_t ringTap, const T * src, uint32_t len);

/// Clip array values between [-1, 1].
void clip1(float * arr, uint32_t len, uint32_t str=1);

/// Finds elements that are within a threshold of their nearest neighbors.

/// @param[in]  src			Source array of elements
/// @param[in]  indices		Index array used for iterating the source array.
/// @param[out] indices		Cluster indices.
/// @param[in]  numIndices	Number of source indices.
/// @param[out] numIndices	Number of cluster indices.
/// @param[in]  threshold	Magnitude threshold of cluster.
TEM void cluster(const T * src, uint32_t * indices, uint32_t& numIndices, T threshold);

void compact(float * dst, const float * src, uint32_t len, uint32_t chunkSize);

/// Returns dot-product of two arrays
TEM inline T dot(const T * src1, const T * src2, uint32_t len, uint32_t str=1){
	T r=T(0); LOOP(len, str){ r += src1[i]*src2[i]; } return r;
}

/// Returns dot-product of two arrays of length 4.
TEM T dot4(const T * src1, const T * src2);

/// Get indices of min and max values.
TEM void extrema(const T * src, uint32_t len, uint32_t& indexMin, uint32_t& indexMax);

/// Perform linear least squares fitting of array.

/// The independent axis is the array indices, i.  The best fit line
/// equation is y = inter + slope * i.
template <class T1, class T2, class T3>
void fitLine(const T1 * src, uint32_t len, T2& slope, T3& inter);

/// Compute histogram of 'src'.

/// Values in 'src' are tallied and placed in 'bins', where the index of the
/// bin is the integer part of the source values.  Source values greater than 
/// the number of bins are ignored.  The scale and offset parameters can be
/// used to put the src values into the proper range.
template <class Ts, class Tb>
void histogram(const Ts * src, uint32_t len, Tb * bins, uint32_t numBins, Ts scale=1);

template <class Ts, class Tb>
void histogram(const Ts * src, uint32_t len, Tb * bins, uint32_t numBins, Ts scale, Ts offset);

/// Mapping from linear range [-1, 1] to normalized dB range [-1, 1].
void linToDB(float * arr, uint32_t len, float minDB);

/// Returns index of maximum value.
TEM uint32_t max(const T * src, uint32_t len, uint32_t str=1);

/// Returns index of maximum normed value (i.e., magnitude).
TEM uint32_t maxNorm(const T * src, uint32_t len, uint32_t str=1);

/// Locates local maxima and writes their indices into 'dst'.

///	Returns number of maxima found.
///
TEM uint32_t maxima(uint32_t * dst, const T * src, uint32_t len, uint32_t str=1);

/// Returns the mean (average) value of array values
TEM inline T mean(const T * src, uint32_t len, uint32_t str=1){
	return sum(src, len, str) / T(len/str);
}

/// Returns the mean norm of array values.
TEM inline T meanNorm(const T * src, uint32_t len, uint32_t str=1){
	T r=T(0); LOOP(len,str){ r += gam::norm(src[i]); } return r / T(len/str);
}

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

TEM void minimaRemove(const T * src, uint32_t * indices, uint32_t& numIndices);

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

/// Multiply destination array by source array
TEM void mul(T * dst, const T * src, uint32_t len, uint32_t str=1){
	LOOP(len,str){ dst[i] *= src[i]; }
}

/// Multiply array by a Bartlett (triangle) window.

/// Works only for even sized arrays.
///
TEM void mulBartlett(T * arr, uint32_t len);

/// Multiply 'arr' by 'src' where 'src' is the first 'len'/2 + 1 elements
/// of a symmetric window.
TEM void mulHalfWindow(T * arr, const T * src, uint32_t len);

/// Uniformly scale array values to fit in [-1, 1].

/// Returns the applied normalization multiplication factor.
///
TEM double normalize(T * arr, uint32_t len, double scale=1);

/// Returns norm of array values.
TEM double norm(const T * src, uint32_t len, uint32_t str=1){
	return sqrt((double)sumSquares(src, len,str));
}

/// Returns taxicab norm of array values (sum of absolute values).
TEM double normTaxi(const T * src, uint32_t len, uint32_t str=1){
	double r=0; LOOP(len,str){ r+=gam::norm(src[i]); } return r;
}

/// Returns unnormalized Nyquist value for use with DFT.
TEM inline T nyquist(const T * src, uint32_t len, uint32_t str=1){
	T r=T(0); LOOP(len,(str<<1)){ r += src[i] - src[i+str]; } return r;
}

TEM void overlapAdd(T * arr, const T * src, uint32_t len, uint32_t hop);

/// Returns root mean square- the normalized norm.
TEM inline T rms(const T * src, uint32_t len, uint32_t str=1){
	return norm(src, len,str) / T(len/str);	
}

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

/// Returns sum of values
TEM inline T sum(const T * src, uint32_t len, uint32_t str=1){
	T r=T(0); LOOP(len, str){ r += src[i]; } return r;
}

/// Returns sum of values squared
TEM inline T sumSquares(const T * src, uint32_t len, uint32_t str=1){
	return dot(src,src,len,str);
}

/// Variance (deviation from mean).
TEM T variance(const T * src, uint32_t len, uint32_t str=1);

/// Returns number of values within [-threshold, theshold].
TEM uint32_t within(const T * src, uint32_t len, T threshold);

/// Returns number of values within [lo, hi].
TEM uint32_t within(const T * src, uint32_t len, T lo, T hi);

/// Returns number of values that equal zero.
TEM uint32_t zeroCount(const T * src, uint32_t len, uint32_t str=1);

/// Returns number of zero-crossings in array.

/// 'prev' is the last value from the previous buffer.
///
uint32_t zeroCross(const float * src, uint32_t len, float prev);

TEM void zeroCross(const T * src, uint32_t len, uint32_t& nzc, uint32_t& pzc);

/// Returns index of first zero-crossing or 0 if none detected.
uint32_t zeroCrossFirst(const float * src, uint32_t len);

TEM uint32_t zeroCrossMax(const T * src, uint32_t len);

/// Returns # of negative slope zero-crossings.

/// 'prev' is the last value from the previous buffer.
///
uint32_t zeroCrossN(const float * src, uint32_t len, float prev);




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
void rectToPolar(float * real, float * imag, uint32_t len, uint32_t str=1);

/// Fast in-place rectangular to polar conversion.

/// Call conversionInit() before using.
///
void rectToPolarFast(float * real, float * imag, uint32_t len);


//--- TODO: consider for removal or refactoring

/// Estimate the fundamental frequency of a spectrum using HPS.

/// Returns index of detected fundamental.
///
uint32_t fundHPS(float * tmp, const float * mag, uint32_t len, uint32_t downSample=4);

/// Computes harmonic product spectrum.

/// 'downSample' specifies how many times to downsample and multiply
///
void hps(float * dst, const float * src, uint32_t len, uint32_t downSample);




// Implementation_______________________________________________________________

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
	LOOP(len, 2){ *arr2-- = -*arr++; }
}

TEM inline void mirror_dq(T * arr, uint32_t len){
	T * arr2 = arr + (len>>1);	// 2nd half begin
	LOOP(len, 2){ *arr2++ = -*arr++; }
}

TEM inline void mulBartlett(T * arr, uint32_t len){
	T * end = arr + len - 1;
	uint32_t len_2 = len >> 1;
	const T slope = (T)1. / (T)len_2;
	T line = slope;
	
	*arr++ = (T)0;
	LOOP(len_2 - 1, 1){
		*arr++ *= line;
		*end-- *= line;
		line += slope;
	}
}

TEM inline void mulHalfWindow(T * arr, const T * src, uint32_t len){
	T * end = arr + len - 1;
	len >>= 1;
	
	*arr++ *= *src++;
	LOOP(len - 1, 1){
		const T val = *src++;
		*arr++ *= val;
		*end-- *= val;
	}
	*arr *= *src;
}

TEM double normalize(T * arr, uint32_t len, double scale){
	double max = gam::norm(arr[maxNorm(arr, len)]);
	double normFactor = scale/max;
	if(max != 0.){ for(uint32_t i=0; i<len; ++i){ arr[i]*=normFactor; } }
	return normFactor;
}

TEM inline void overlapAdd(T * arr, const T * src, uint32_t len, uint32_t hop){
	uint32_t lap = len - hop;
	add(arr, arr + hop, src, lap);
	mem::copy(arr + lap, src + lap, hop);
}

TEM void cluster(const T * src, uint32_t * indices, uint32_t & numIndices, T threshold){
	
	if(numIndices == 0) return;
	
	uint32_t * newIndices = indices;
	uint32_t newNumIndices = 0;
	
	T prev = src[*indices++];
	bool inCluster = false;
	
	LOOP(numIndices - 1, 1){
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
	}
	
	numIndices = newNumIndices;
}

TEM inline T dot4(const T * a, const T * b){
	return a[0]*b[0] + a[1]*b[1] + a[2]*b[2] + a[3]*b[3];
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

 	LOOP(len, 1){
		cov += dx++ * (*src++ - meanY);
 	}

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
	LOOP(len, 1){
		int32_t j = (int32_t)(src[i] * scale);
		if(j >= 0 && j < (int32_t)numBins) bins[j]++;
	}
}

template <class Ts, class Tb>
inline void histogram(const Ts * src, uint32_t len, Tb * bins, uint32_t numBins, Ts scale, Ts offset){
	LOOP(len, 1){
		int32_t j = (int32_t)(src[i] * scale + offset);
		if(j >= 0 && j < (int32_t)numBins) bins[j]++;
	}
}

TEM uint32_t max(const T * src, uint32_t len, uint32_t str){
	uint32_t r=0;
	T max = src[0];
	for(uint32_t i=str; i<len; i+=str){
		const T& v = src[i];
		if(v > max){ max=v; r=i; }
	}
	return r;
}

TEM uint32_t maxNorm(const T * src, uint32_t len, uint32_t str){
	uint32_t r = 0;
	double max = normCompare(src[0]);
	for(uint32_t i=str; i<len; i+=str){
		double v = normCompare(src[i]);
		if(v > max){ max=v; r=i; }
	}
	return r;
}

TEM uint32_t maxima(uint32_t * dst, const T * src, uint32_t len, uint32_t str){
	T prev = src[0];
	T curr = src[str];
	
	uint32_t numPeaks = 0;
	
	for(uint32_t i=(str<<1); i<len; i+=str){
		T next = src[i];

		if(curr > prev && curr > next){
			*dst++ = i-str;
			numPeaks++;
			i+=str;
			if(i<len){
				prev = next;
				curr = src[i];
			}
		}
		else{
			prev = curr;
			curr = next;
		}
	}
	return numPeaks;
}



TEM T meanAbsDiff(const T * src, uint32_t len){
	T sum = (T)0;
	T mean = (T)0;
	
	T prev = *src++;
    LOOP(len - 1, 1){
		T now = *src++;
		T diff = now - prev;
		T abs = scl::abs(now);
		sum += scl::abs(diff);
		mean += abs;
		prev = now;
    }
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
	LOOP(len,1){ mean += *src++ * *weights++ * normFactor; }
	return mean;
}

TEM T meanWeightedIndex(const T * weights, uint32_t len){
	T mean = (T)0;
	T normFactor = sum(weights, len);
	if(normFactor == (T)0) return (T)0;
	normFactor = (T)1 / normFactor;
	T weightFactor = normFactor;
	LOOP(len,1){
		mean += *weights++ * weightFactor;
		weightFactor += normFactor;
	}
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


TEM uint32_t slopeAbsMax(const T * src, uint32_t len){
	uint32_t index = 0;
	T prev = *src++;
	T maxSlope = (T)0;

	LOOP(len-1,1){
		T now = *src++;
		T slope = scl::abs(now - prev);
		if(slope > maxSlope){
			maxSlope = slope;
			index = i;
		}
		prev = now;
	}
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

TEM inline T variance(const T * src, uint32_t len, uint32_t str){
	T rec = T(1)/T(len/str);
	T avg = sum(src, len,str) * rec;
	T r = T(0);
	
	LOOP(len, str){
		T d = src[i] - avg;
		r += d*d;
	}
	return r * rec;
}

TEM uint32_t within(const T * src, uint32_t len, T threshold){
	uint32_t count = 0;
	LOOP(len, 1){
		T val = src[i];
		if((val <= threshold) && (val >= -threshold)) count++;
	}
	return count;
}

TEM uint32_t within(const T * src, uint32_t len, T lo, T hi){
	uint32_t count = 0;
	LOOP(len, 1){
		T val = src[i];
		if((val <= hi) && (val >= lo)) count++;
	}
	return count;
}

TEM inline uint32_t zeroCount(const T * src, uint32_t len, uint32_t str){
	uint32_t r=0; LOOP(len,str){ if(src[i] == T(0)) r++; } return r;
}

TEM void zeroCross(const T * src, uint32_t len, uint32_t& nzc, uint32_t& pzc){
	pzc = 0;
	nzc = 0;
	T prev = *src++;
	LOOP(len - 1, 1){
		T curr = src[i];
		if		(curr > (T)0 && prev <= (T)0) pzc++;
		else if	(curr <= (T)0 && prev > (T)0) nzc++;
		prev = curr;
	}
}

TEM uint32_t zeroCrossMax(const T * src, uint32_t len){
	T prev = *src++;
	T max = 0;
	uint32_t ind = 0;
	
	LOOP(len-1, 1){
		T curr = src[i];
		if(curr > (T)0 && prev <= (T)0){
			T slope = curr - prev;
			if(slope > max){
				max = slope;
				ind = i;
			}
		}
		prev = curr;
	}
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
	LOOP(len, 1){
		T dp = scl::wrapPhase(*p0 - *p1);		// wrap phase into [-pi, pi)
		*p1++ = *p0;							// prev phase = curr phase
		*p0++ = dp * factor;
	}
}

TEM void phaseToFreq(T * f, const T * p0, const T * p1, uint32_t len, T ups){
	T factor = (T)1 / (M_2PI * ups);
	LOOP(len, 1){
		T dp = scl::wrapPhase(*p0 - *p1);		// wrap phase into [-pi, pi)
		*f++ = dp * factor;
		p0++; p1++;
	}
//	subtract(f, p0, p1, len);
//	mul(f, len, factor);
//	wrap(f, len, -0.5 / ups, 0.5 / ups);
}

} // arr::
} // gam::

#undef LOOP
#undef TEM

#endif
