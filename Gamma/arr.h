#ifndef GAMMA_ARR_H_INC
#define GAMMA_ARR_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information

	File Description:
	Functions for processing arrays of data.
*/

#include <cmath>
#include "Gamma/scl.h" // norm
#include "Gamma/std_swap.h"

#define LOOP(n,s) for(unsigned i=0; i<n; i+=s)

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


// These functions are defined at the top since they are used by other functions
// defined in this file...

/// Returns dot-product of two arrays
template <class T> 
inline T dot(const T * src1, const T * src2, unsigned len, unsigned str=1){
	T r=T(0); LOOP(len, str){ r += src1[i]*src2[i]; } return r;
}

/// Returns sum of values squared
template <class T>
inline T sumSquares(const T * src, unsigned len, unsigned str=1){
	return dot(src,src,len,str);
}




/// Add source array to destination array
template <class T>
void add(T * dst, const T * src, unsigned len, unsigned str=1){
	LOOP(len,str){ dst[i] += src[i]; }
}

/// Add two source arrays into destination array
template <class T>
void add(T * dst, const T * src1, const T * src2, unsigned len, unsigned str=1){
	LOOP(len,str){ dst[i] = src1[i] + src2[i]; }
}

/// Sum elements from src into ring-buffer ring.

/// Returns the next tap index.  This will not guaranteed to be in the range [0, ringSize).
///
template <class T>
unsigned addToRing(T * ring, unsigned ringSize, unsigned ringTap, const T * src, unsigned len);

/// Clip array values between [-1, 1].
void clip1(float * arr, unsigned len, unsigned str=1);

/// Finds elements that are within a threshold of their nearest neighbors.

/// \param[in]  src			Source array of elements
/// \param[in]  indices		Index array used for iterating the source array.
/// \param[out] indices		Cluster indices.
/// \param[in]  numIndices	Number of source indices.
/// \param[out] numIndices	Number of cluster indices.
/// \param[in]  threshold	Magnitude threshold of cluster.
template <class T, class Index>
void cluster(const T * src, Index * indices, unsigned& numIndices, T threshold);

void compact(float * dst, const float * src, unsigned len, unsigned chunkSize);

/// Returns dot-product of two arrays of length 4.
template <class T>
T dot4(const T * src1, const T * src2);

/// Get indices of min and max values.
template <class T, class Index>
void extrema(const T * src, unsigned len, Index& indexMin, Index& indexMax);

/// Compute histogram of 'src'.

/// Values in 'src' are tallied and placed in 'bins', where the index of the
/// bin is the integer part of the source values.  Source values greater than 
/// the number of bins are ignored.  The scale and offset parameters can be
/// used to put the src values into the proper range.
template <class Ts, class Tb>
void histogram(const Ts * src, unsigned len, Tb * bins, unsigned numBins, Ts scale=1);

template <class Ts, class Tb>
void histogram(const Ts * src, unsigned len, Tb * bins, unsigned numBins, Ts scale, Ts offset);

/// Returns index of maximum value
template <class T>
unsigned indexOfMax(const T * src, unsigned len, unsigned str=1);

/// Returns index of maximum normed value (i.e., magnitude)
template <class T>
unsigned indexOfMaxNorm(const T * src, unsigned len, unsigned str=1);

/// Returns index of minimum value
template <class T>
unsigned indexOfMin(const T * src, unsigned len, unsigned str=1);

/// Sets indices [numIndices, maxNumIndices) to complement indices.

/// Indices must be sorted from low to high.
///
template <class Index>
void indicesComplement(Index * indices, unsigned numIndices, unsigned maxNumIndices);

/// Perform linear least squares fitting of array.

/// This computes the linear equation 'y = slope * x + inter' that best fits the
/// array values in a least squares sense. Here, the independent variable, x, is
/// the array index and the dependent variable, y, is the array element value.
template <class T1, class T2, class T3>
void lineFit(const T1 * src, unsigned len, T2& slope, T3& inter);

/// Mapping from linear range [-1, 1] to normalized dB range [-1, 1].
void linToDB(float * arr, unsigned len, float minDB);

/// Locates local maxima and writes their indices into 'dst'.

///	Returns number of maxima found.
///
template <class Index, class T>
unsigned maxima(Index * dst, const T * src, unsigned len, unsigned str=1);

/// Returns the mean (average) value of array values
template <class T>
T mean(const T * src, unsigned len, unsigned str=1);

/// Returns the mean norm of array values.
template <class T>
inline T meanNorm(const T * src, unsigned len, unsigned str=1){
	T r=T(0); LOOP(len,str){ r += gam::norm(src[i]); } return r / T(len/str);
}

/// Returns mean absolute difference of array values.
template <class T>
T meanAbsDiff(const T * src, unsigned len);

/// Returns weighted mean of array values.

/// Weights must be positive.
///
template <class T>
T meanWeighted(const T * src, const T * weights, unsigned len);

/// Returns weighted mean in [0, len) of indices of weights.

/// Weights must be positive.
///		Can be used to compute centroid of spectrum.
template <class T>
T meanWeightedIndex(const T * weights, unsigned len);

/// Computes min, max and mean values of array
template <class T>
void minMaxMean(T& min, T& max, T& mean, const T * src, unsigned len, unsigned str=1);

template <class T, class Index>
void minimaRemove(const T * src, Index * indices, unsigned& numIndices);

//	/// Applies mirror isometry sequence [dbqp] from first quarter of array.
//	
//	/// The sequence of mirror isometries are identity (d), reflection (b),
//	/// glide reflection (q), and rotation (p). The array should hold the first
//	/// len/4 + 1 elements of the signal.\n
//	/// Ex.: [ 1, 2, 3, x, x, x, x, x] -> [ 1, 2, 3, 2, -1,-2,-3,-2]
//	/// Ex.: [ 1, 2, x, x, x, x, x, x] -> [ 1, 2, 2, 1, -1,-2,-2,-1]
//	template <class T> void mirror_dbqp(T * arr, unsigned len);

/// Applies mirror isometry sequence [dp] from first half of array.

/// The sequence of mirror isometries are identity (d) and rotation (p).
/// The first len/2 elements of the array are mirrored.\n
/// Ex.: [ 1, 2, 3, 4, x, x, x, x] -> [ 1, 2, 3, 4,-4,-3,-2,-1]
template <class T>
void mirror_dp(T * arr, unsigned len);

/// Applies mirror isometry sequence [dq] from first half of array.

/// The sequence of mirror isometries are identity (d) and glide relfection (q).
/// The first len/2 elements of the array are mirrored.\n
/// Ex.: [ 1, 2, 3, 4, x, x, x, x] -> [ 1, 2, 3, 4,-1,-2,-3,-4]
template <class T>
void mirror_dq(T * arr, unsigned len);

/// Multiply destination array by source array
template <class T>
void mul(T * dst, const T * src, unsigned len, unsigned str=1){
	LOOP(len,str){ dst[i] *= src[i]; }
}

/// Multiply array by a Bartlett (triangle) window.

/// Works only for even sized arrays.
///
template <class T>
void mulBartlett(T * arr, unsigned len);

/// Multiply 'arr' by 'src' where 'src' is the first 'len'/2 + 1 elements
/// of a symmetric window.
template <class T>
void mulHalfWindow(T * arr, const T * src, unsigned len);

/// Uniformly scale array values to fit in [-1, 1].

/// Returns the applied normalization multiplication factor.
///
template <class T>
double normalize(T * arr, unsigned len, double scale=1);

template<class T, template<class> class ArrayType>
double inline normalize(ArrayType<T>& arr, double scale=1){
	return normalize(&arr[0], arr.size(), scale);
}

/// Returns norm of array values.
template <class T>
double norm(const T * src, unsigned len, unsigned str=1){
	return sqrt((double)sumSquares(src, len,str));
}

/// Returns taxicab norm of array values (sum of absolute values).
template <class T>
double normTaxi(const T * src, unsigned len, unsigned str=1){
	double r=0; LOOP(len,str){ r+=gam::norm(src[i]); } return r;
}

/// Returns unnormalized Nyquist value for use with DFT.
template <class T>
inline T nyquist(const T * src, unsigned len, unsigned str=1){
	T r=T(0); LOOP(len,(str<<1)){ r += src[i] - src[i+str]; } return r;
}

/// Returns root mean square- the normalized norm.
template <class T>
inline T rms(const T * src, unsigned len, unsigned str=1){
	return norm(src, len,str) / T(len/str);	
}

/// Returns index of absolute maximum slope in array.
template <class T>
unsigned slopeAbsMax(const T * src, unsigned len);

/// Returns index of maximum slope in array.
template <class T>
unsigned slopeMax(const T * src, unsigned len);

/// Insertion sort of elements.

/// Elements are sorted from lowest to highest.
/// This sort is fastest for small length arrays and mostly sorted sets.
template <class T>
void sortInsertion(T * arr, unsigned len);

/// Insertion sort of indexed elements.

/// Elements are sorted from lowest to highest.
/// This sort is fastest for small length arrays and mostly sorted sets.
template <class T, class Index>
void sortInsertion(const T * src, Index * indices, unsigned numIndices);

/// Quick sort of elements.

/// Elements are sorted from lowest to highest.
///
template <class T, class Index>
void sortQuick(const T * src, Index * indices, long beg, long end);

/// Returns sum of values
template <class T>
T sum(const T * src, unsigned len, unsigned str=1);

/// Variance (deviation from mean).
template <class T>
T variance(const T * src, unsigned len, unsigned str=1);

/// Returns number of values within [-threshold, theshold].
template <class T>
unsigned within(const T * src, unsigned len, T threshold);

/// Returns number of values within [lo, hi].
template <class T>
unsigned within(const T * src, unsigned len, T lo, T hi);

/// Returns number of values that equal zero.
template <class T>
unsigned zeroCount(const T * src, unsigned len, unsigned str=1);

/// Returns number of zero-crossings in array.

/// 'prev' is the last value from the previous buffer.
///
unsigned zeroCross(const float * src, unsigned len, float prev);

template <class T, class Index>
void zeroCross(const T * src, unsigned len, Index& nzc, Index& pzc);

/// Returns index of first zero-crossing or 0 if none detected.
unsigned zeroCrossFirst(const float * src, unsigned len);

template <class T>
unsigned zeroCrossMax(const T * src, unsigned len);

/// Returns # of negative slope zero-crossings.

/// 'prev' is the last value from the previous buffer.
///
unsigned zeroCrossN(const float * src, unsigned len, float prev);




// Implementation_______________________________________________________________

template <class T>
unsigned addToRing(T * ring, unsigned ringSize, unsigned ringTap, const T * src, unsigned len){
	unsigned endTap = ringTap + len;

	if(endTap <= ringSize){		// haven't gone past end
		add(ring + ringTap, src, len);
	}
	else{						// have gone past end, do wrapped addition
		unsigned samplesUnder	= ringSize - ringTap;
		unsigned samplesOver	= endTap - ringSize;
		add(ring + ringTap, src, samplesUnder);
		add(ring, src + samplesUnder, samplesOver);
	}
	return endTap;
}

//template <class T> inline void mirror_dbqp(T * arr, unsigned len){
//	
//	T * arr3 = arr + (len>>1);	// 3rd quad start
//	T * arr2 = arr3 - 2;		// 2nd quad end
//	T * arr4 = arr + len - 2;	// 4th quad end
//	
//	//for(unsigned j=1; j<=(len>>2); ++j){
//	for(unsigned j=0; j<len; j+=4){
//		float val = *arr++;
//		*arr2-- =  val;
//		*arr3++ = -val;
//		*arr4-- = -val;
//	}
//}

template <class T>
inline void mirror_dp(T * arr, unsigned len){
	T * arr2 = arr + len - 1;	// 2nd half end
	LOOP(len, 2){ *arr2-- = -*arr++; }
}

template <class T>
inline void mirror_dq(T * arr, unsigned len){
	T * arr2 = arr + (len>>1);	// 2nd half begin
	LOOP(len, 2){ *arr2++ = -*arr++; }
}

template <class T>
inline void mulBartlett(T * arr, unsigned len){
	T * end = arr + len - 1;
	unsigned len_2 = len >> 1;
	const T slope = (T)1. / (T)len_2;
	T line = slope;
	
	*arr++ = (T)0;
	LOOP(len_2 - 1, 1){
		*arr++ *= line;
		*end-- *= line;
		line += slope;
	}
}

template <class T>
inline void mulHalfWindow(T * arr, const T * src, unsigned len){
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

template <class T>
double normalize(T * arr, unsigned len, double scale){
	double max = gam::norm(arr[indexOfMaxNorm(arr, len)]);
	double normFactor = scale/max;
	if(max != 0.){
		for(unsigned i=0; i<len; ++i){
			arr[i]*=normFactor;
		}
	}
	return normFactor;
}

template <class T, class Index>
void cluster(const T * src, Index * indices, unsigned& numIndices, T threshold){

	if(numIndices == 0) return;

	auto * newIndices = indices;
	unsigned newNumIndices = 0;

	T prev = src[*indices++];
	bool inCluster = false;

	LOOP(numIndices - 1, 1){
		auto index = *indices++;
		T curr = src[index];

		if( std::abs(curr - prev) <= threshold ){
			if(!inCluster){
				// Add previous index
				*newIndices++ = indices[-2];
				++newNumIndices;
			}

			// Add current index
			*newIndices++ = index;
			++newNumIndices;

			inCluster = true;
		}
		else{
			inCluster = false;
		}
		prev = curr;
	}
	numIndices = newNumIndices;
}

template <class T>
inline T dot4(const T * a, const T * b){
	return a[0]*b[0] + a[1]*b[1] + a[2]*b[2] + a[3]*b[3];
}

template <class T, class Index>
void extrema(const T * src, unsigned len, Index& idxMin, Index& idxMax){
	idxMin = 0;
	idxMax = 0;

	T min = *src++;
	T max = min;

	for(unsigned i = 1; i < len; i++){
		T val = *src++;
		if(val > max){
			max = val;
			idxMax = i;
		}
		else if(val < min){
			min = val;
			idxMin = i;
		}
	}
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
inline void histogram(const Ts * src, unsigned len, Tb * bins, unsigned numBins, Ts scale){
	LOOP(len, 1){
		auto j = int(src[i] * scale);
		if(j >= 0 && j < int(numBins)) bins[j]++;
	}
}

template <class Ts, class Tb>
inline void histogram(const Ts * src, unsigned len, Tb * bins, unsigned numBins, Ts scale, Ts offset){
	LOOP(len, 1){
		auto j = int(src[i] * scale + offset);
		if(j >= 0 && j < int(numBins)) bins[j]++;
	}
}

template <class T>
unsigned indexOfMax(const T * src, unsigned len, unsigned str){
	unsigned r = 0;
	T mx = src[0];
	for(unsigned i=str; i<len; i+=str){
		const T& v = src[i];
		if(v > mx){ mx=v; r=i; }
	}
	return r;
}

template <class T>
unsigned indexOfMaxNorm(const T * src, unsigned len, unsigned str){
	unsigned r = 0;
	double max = normCompare(src[0]);
	for(unsigned i=str; i<len; i+=str){
		double v = normCompare(src[i]);
		if(v > max){ max=v; r=i; }
	}
	return r;
}

template <class T>
unsigned indexOfMin(const T * src, unsigned len, unsigned str){
	unsigned r = 0;
	T mn = src[0];
	for(unsigned i=str; i<len; i+=str){
		const T& v = src[i];
		if(v < mn){ mn=v; r=i; }
	}
	return r;
}

/*
	slope = cov(x, y) / var(x)
	inter = x_mean - slope * y_mean
*/
template <class T, class T2, class T3>
void lineFit(const T * src, unsigned len, T2& slope, T3& inter){
	T lenT  = T(len);
	T meanX = (lenT - T(1)) * T(0.5);	// mean of independent variables (the indices)
	T meanY = sum(src, len) / lenT;		// mean of dependent variables
	T sumSqrs = meanX*(meanX+T(1))*(meanX*T(2./6) + T(1./6));
	T varX  = T(2)*sumSqrs; // variance of x

	T cov  = T(0);	// covariance
	T dx   =-meanX;	// current distance from point to center along x

	LOOP(len, 1){
		cov += dx++ * (*src++ - meanY);
	}

	slope = cov / varX;
	inter = meanY - slope * meanX;
}

template <class Index, class T>
unsigned maxima(Index * dst, const T * src, unsigned len, unsigned str){
	T prev = src[0];
	T curr = src[str];
	
	unsigned numPeaks = 0;
	
	for(unsigned i=(str<<1); i<len; i+=str){
		T next = src[i];

		if(curr > prev && curr > next){
			*dst++ = i-str;
			++numPeaks;
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


template <class T>
inline T mean(const T * src, unsigned len, unsigned str){
	return arr::sum(src, len, str) / T(len/str);
}


template <class T>
T meanAbsDiff(const T * src, unsigned len){
	T sum = (T)0;
	T mean = (T)0;
	
	T prev = *src++;
	LOOP(len - 1, 1){
		T now = *src++;
		T diff = now - prev;
		sum += std::abs(diff);
		mean += std::abs(now);
		prev = now;
	}
	//if(mean < 0.0001f) mean = 0.0001f;
	return mean < T(0.25) ? T(-1) : sum * T(0.5) / mean;

}

template <class T>
T meanWeighted(const T * src, const T * weights, unsigned len){
	T mean = T(0);

	// One loop: faster, but less accurate
//	LOOP(len,
//		mean += *src++ * *weight;
//		normFactor += *weight++;
//	)
//	return mean /= normFactor;

	// Two loops: slower, but avoids overflow
	T normFactor = T(1) / sum(weights, len);
	LOOP(len,1){ mean += *src++ * *weights++ * normFactor; }
	return mean;
}

template <class T>
T meanWeightedIndex(const T * weights, unsigned len){
	T mean = T(0);
	T normFactor = sum(weights, len);
	if(normFactor == T(0)) return T(0);
	normFactor = T(1) / normFactor;
	T weightFactor = T(0);
	LOOP(len,1){
		mean += *weights++ * weightFactor;
		weightFactor += normFactor;
	}
	return mean;
}

template <class T>
void minMaxMean(T& min, T& max, T& mean, const T * src, unsigned len, unsigned str){
	min = max = src[0];
	mean = T(0);
	LOOP(len,str){
		const auto& v = src[i];
		mean += v;
		if(v < min) min = v;
		else if (v > max) max = v;
	}
	mean /= (len/str);
}

template <class T, class Index>
void minimaRemove(const T * src, Index * indices, unsigned& numIndices){

	if(numIndices < 3) return;

	auto * newIndices = indices + 1;	// first index is never a minima
	unsigned newNumIndices = 2;			// always keep first and last indices

	T prev = src[*indices++];
	T curr = src[*indices++];

	for(unsigned i=1; i<numIndices-1; i++){
		auto index = *indices++;
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

template <class T>
unsigned slopeAbsMax(const T * src, unsigned len){
	unsigned index = 0;
	T prev = *src++;
	T maxSlope = (T)0;

	LOOP(len-1,1){
		T now = *src++;
		T slope = std::abs(now - prev);
		if(slope > maxSlope){
			maxSlope = slope;
			index = i;
		}
		prev = now;
	}
	return index;
}

template <class T>
unsigned slopeMax(const T * src, unsigned len){
	unsigned index = 0;
	T prev = *src++;
	T maxSlope = (T)0;

	for(unsigned i=0; i<len-1; ++i){
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

template <class T>
void sortInsertion(T * arr, unsigned len){
	for(unsigned i = 1; i < len; i++){
		T val = arr[i];
		unsigned j = i - 1;
		for(; (j < len) && (arr[j] > val); j--){
			arr[j + 1] = arr[j];
		}
		arr[j + 1] = val;
	} 
}

template <class T, class Index>
void sortInsertion(const T * src, Index * indices, unsigned numIndices){
	for(unsigned i = 1; i < numIndices; i++)
	{
		auto index = indices[i];
		T val = src[index];
		unsigned j = i - 1;
		while((j < numIndices) && (val < src[indices[j]]))
		{
			indices[j + 1] = indices[j]; 
			--j;
		}
		indices[j + 1] = index;
	} 
}

template <class T, class Index>
void sortQuick(const T * src, Index * indices, long beg, long end){
	// must be at least 1 element to sort
	if(end > beg + 1) {
		long piv = indices[beg], l = beg + 1, r = end;
		T pivVal = src[piv];
		while (l < r){
			if ( src[indices[l]] <= pivVal ) 
				l++;
			else 
				std::swap(indices[l], indices[--r]);
		}
		std::swap(indices[--l], indices[beg]);
		sortQuick(src, indices, beg, l);
		sortQuick(src, indices, r, end);
	}
}

template <class T>
inline T sum(const T * src, unsigned len, unsigned str){
	T r=T(0);
	LOOP(len, str){ r += src[i]; }
	return r;
}

template <class T>
inline T variance(const T * src, unsigned len, unsigned str){
	T rec = T(1)/T(len/str);
	T avg = sum(src, len,str) * rec;
	T r = T(0);
	
	LOOP(len, str){
		T d = src[i] - avg;
		r += d*d;
	}
	return r * rec;
}

template <class T>
inline unsigned within(const T * src, unsigned len, T threshold){
	return within(src,len, -threshold, threshold);
}

template <class T>
unsigned within(const T * src, unsigned len, T lo, T hi){
	unsigned count = 0;
	LOOP(len, 1){
		T val = src[i];
		count += unsigned((val <= hi) && (val >= lo));
	}
	return count;
}

template <class T>
inline unsigned zeroCount(const T * src, unsigned len, unsigned str){
	unsigned r=0;
	LOOP(len,str){
		r += unsigned(src[i] == T(0));
	}
	return r;
}

template <class T, class Index>
void zeroCross(const T * src, unsigned len, Index& nzc, Index& pzc){
	pzc = 0;
	nzc = 0;
	T prev = *src++;
	LOOP(len - 1, 1){
		T curr = src[i];
		if		(curr > (T)0 && prev <= (T)0) ++pzc;
		else if	(curr <= (T)0 && prev > (T)0) ++nzc;
		prev = curr;
	}
}

template <class T>
unsigned zeroCrossMax(const T * src, unsigned len){
	T prev = *src++;
	T max = 0;
	unsigned ind = 0;
	
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

template <class Index>
void indicesComplement(Index * indices, unsigned numIndices, unsigned maxNumIndices){
	auto * comp = indices + numIndices;
	for(unsigned i=0; i<maxNumIndices; i++){
		if(*indices == i)	++indices;
		else				*comp++ = i;
	}
}

} // arr::
} // gam::

#undef LOOP

#endif
