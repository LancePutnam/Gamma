#ifndef GAMMA_TBL_H_INC
#define GAMMA_TBL_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */
 
#include <stdlib.h>
#include "Gamma/arr.h"
#include "Gamma/mem.h"
#include "Gamma/scl.h"
#include "Gamma/Constants.h"

#define TEM template<class T>
#define LOOP(n,s) for(uint32_t i=0; i<n; i+=s)

namespace gam{


/// Various window types.
class WinType{
public:

	/// Window types.
	enum type{ 	
		Bartlett,			/**< Bartlett (Triangle)*/
		Blackman,			/**< Blackman */
		BlackmanHarris,		/**< Blackman-Harris */
		Hamming,			/**< Hamming */
		Hann,				/**< von Hann */
		Welch,				/**< Welch */
		Nyquist,			/**< Nyquist */
		Rectangle			/**< Rectangle (1s) */
	};
	
	/// Returns human readable string of window type.
	inline static const char * string(WinType::type type){
		#define CS(name) case WinType::name: return #name;
		switch(type){
			CS(Bartlett) CS(Blackman) CS(BlackmanHarris) CS(Hamming) CS(Hann)
			CS(Welch) CS(Nyquist) CS(Rectangle)
			default: return "Unknown";
		}
		#undef CS
	}
	
};

/// Table functions
namespace tbl{

/// Fills array with one period of a cosine wave.
TEM void cosine(T * dst, uint32_t len);

/// Fills table with section of an exponential decay.

///	Values are normalized to descend from 1 to 0. \n
///	Negative 'order' curves downward and positive 'order' curves upward. \n
///	If 'order' is 0, you get a line.
TEM void decay(T * dst, uint32_t len, double order);

/// Fills array with one period of a sine wave.
TEM void sine(T * dst, uint32_t len);

/// Fills array with arbitrary phase and length sinusoid.
TEM void sinusoid(T * dst, uint32_t len, double phase, double periods);

// Max harmonics
// [len, -1, 0, -1, ..., 0, -1]
// Max harmonics - 1
// [len - 1, 0, -1, 0, ..., -1, 0]
TEM void impulseSum(T * dst, uint32_t len);
	
/// Sums band-limited impulse wave into array.

/// The waveform includes harmonics in the range [hrmLo, hrmHi].
/// The amplitude of the waveform will not be normalized.
/// The ideal waveform shape is [4, -1, 0, -1, 0, -1, 0, -1 ]
TEM void impulseSum(T * dst, uint32_t len, uint32_t hrmLo, uint32_t hrmHi);

/// Sums band-limited saw wave into array.

/// The waveform includes harmonics in the range [hrmLo, hrmHi].
/// The ideal waveform shape is [1, 0.75, 0.5, 0.25, 0, -0.25, -0.5, -0.75]
TEM void sawSum(T * dst, uint32_t len, uint32_t hrmLo, uint32_t hrmHi);

/// Sums band-limited square wave into array.

/// The waveform includes harmonics in the range [hrmLo, hrmHi].
///	The ideal waveform shape is [ 1, 1, 1, 1, -1, -1, -1, -1].
TEM void squareSum(T * dst, uint32_t len, uint32_t hrmLo, uint32_t hrmHi);

/// Sums band-limited triangle wave into array.
	
/// The waveform includes harmonics in the range [hrmLo, hrmHi].
///	The ideal waveform shape is [ 0, 0.5, 1, 0.5, 0, -0.5, -1, -0.5].
TEM void triangleSum(T * dst, uint32_t len, uint32_t hrmLo, uint32_t hrmHi);

/// 
TEM void multiWave(T * dst, uint32_t len, uint32_t order, void (* func)(T *, uint32_t, uint32_t, uint32_t));

/// Returns maximum number of harmonics that will fit in array.
inline uint32_t maxHarmonics(uint32_t len){ return len>>1; }

/// Fills array with specified window type.
TEM void window			(T * dst, uint32_t len, WinType::type type);
TEM void bartlett		(T * dst, uint32_t len); ///< Fills array with Bartlett window.
TEM void blackman		(T * dst, uint32_t len); ///< Fills array with Blackman window.
TEM void blackmanHarris	(T * dst, uint32_t len); ///< Fills array with Blackman-Harris window.
TEM void hamming		(T * dst, uint32_t len); ///< Fills array with Hamming window.
TEM void hann			(T * dst, uint32_t len); ///< Fills array with von Hann window.
TEM void welch			(T * dst, uint32_t len); ///< Fills array with Welch window.
TEM void rectangle		(T * dst, uint32_t len); ///< Fills array with Rectangle window.
TEM void nyquist		(T * dst, uint32_t len, uint32_t str=1); ///< Fills array with Nyquist window.


//
// Accessing
//

// Return value from a table with the first half of a dq-symmetric 
// waveform.  The table size must be a power of two.
//
//	'table':	first half of waveform
//	'fbits':	= 31 - (# bits in table)
//	'phase':	phase of lookup (full waveform period is [0, 2^32))
//
//	'phase' bit format (b = 'fbits'):
//	32:			sign bit (0 = positive, 1 = negative)
//	[31, b]:	phase integer part
//	[ b, 0]:	phase fractional part
float atH(const float * src, uint32_t fbits, uint32_t phase);

// Return value from a table with the first quarter of a dbqp-symmetric 
// waveform.  The table size must be a power of two plus one.
//
//	'table':	first quarter of waveform
//	'fbits':	= 30 - (# bits in table)
//	'phase':	phase of lookup (full waveform period is [0, 2^32))
//
//	'phase' bit format (b = 'fbits'):
//	32:			sign bit (0 = positive, 1 = negative)
//	31:			direction bit (0 = forward, 1 = backward)
//	[30, b]:	phase integer part
//	[ b, 0]:	phase fractional part
float atQ(const float * src, uint32_t fbits, uint32_t phase);

// Get interpolation data from integer phasor.
//float getIpol2(uint32_t bits, uint32_t phase, uint32_t &i, uint32_t &ip1);

/// Returns phase increment factor.

///	Multiply by desired frequency in Hz to get integer phase increment.
///
float phaseIncFactor(double framesPerSec);



// Implementation_______________________________________________________________

TEM void cosine(T * dst, uint32_t len){
	double inc = M_2PI / (double)len;
	double phs = inc;
	len >>= 1;
	
	T * dst2 = dst + len;
	
	*dst++  = (T) 1;
	*dst2++ = (T)-1;

	len -= 1;
	LOOP(len, 1){
		T val = (T)cos(phs);
		*dst++  =  val;
		*dst2++ = -val;
		phs += inc;
	}
}

TEM void decay(T * dst, uint32_t len, double order){
	double final = ::pow(2., order);
	double lambda = ::log(final) / (float)len;
	double time = 1.;
	double scale = 1. / (1. - final);
	double offset = -final;
	*dst++ = (T)1;
	len -= 1;
	LOOP(len, 1){
		*dst++ = (T)(::exp(lambda * time) + offset) * scale;
		time++;
	}
}

TEM void sine(T * dst, uint32_t len){
	double inc = M_2PI / (double)len;
	double phs = inc;
	len >>= 1;
	
	T * dst2 = dst + len;
	
	*dst++  = (T)0;
	*dst2++ = (T)0;
	
	len -= 1;
	LOOP(len, 1){
		T val = (T)sin(phs);
		*dst++  =  val;
		*dst2++ = -val;
		phs += inc;
	}
}

// VERY accurate, but not so fast
TEM void sinusoid(T * dst, uint32_t len, double phase, double periods){
	
	double inc = M_2PI * periods / (double)len;

	for(double i = 0.; i < (double)len; i++){
		*dst++ = (T)sin(inc * i + phase);
	}
}

TEM void impulseSum(T * dst, uint32_t len){
	uint32_t harmonics = (len>>1) - 1;
	*dst++ = (T)harmonics;
	*dst++ = (T)0;
	LOOP(harmonics,1){
		*dst++ = (T)-1;
		*dst++ = (T) 0;
	}
}

TEM void impulseSum(T * dst, uint32_t len, uint32_t hrmLo, uint32_t hrmHi){
	double inc = M_2PI / (double)len;
	uint32_t hLen = len >> 1;
	
	for(uint32_t k = hrmLo; k <= hrmHi; ++k){
		double phaseInc = (double)k * inc;
		double phs = 0.;
		
		T * dst1 = dst;
		
		LOOP(hLen+1, 1){
			*dst1++ += (T)(cos(phs));
			phs += phaseInc;
		}
	}
	
	// Extrapolate remaining from [db] symmetry
	mem::reflectRight(dst + 1, len - 1);
}

TEM void sawSum(T * dst, uint32_t len, uint32_t hrmLo, uint32_t hrmHi){

	static const double sawFactor = 2.0 / M_PI;
	double inc = M_2PI / (double)len;
	uint32_t hLen = len >> 1;
	
	dst++;
	
	for(uint32_t i = hrmLo; i <= hrmHi; ++i){
		double h = (double)i;
		double phaseInc = h * inc;
		double phs = phaseInc;
		double amp = sawFactor / h;
		
		T * dst1 = dst;
		
		for(uint32_t j=1; j<hLen; ++j){
			*dst1++ += (T)(amp * sin(phs));
			phs += phaseInc;
		}
	}
	
	// Extrapolate remaining from [dp] symmetry
	arr::mirror_dp(dst, len-1);	
}

TEM void squareSum(T * dst, uint32_t len, uint32_t hrmLo, uint32_t hrmHi){

	static const double sqrFactor = 4.0 / M_PI;
	double inc = M_2PI / (double)len;
	uint32_t qLen = len >> 2;
	
	dst++;
	
	hrmLo |= 1;	// next highest odd if even
	
	// Calculate first quadrant
	for(uint32_t i = hrmLo; i <= hrmHi; i+=2){

		double h = (double)i;
		double phaseInc = h * inc;
		double phs = phaseInc;
		double amp = sqrFactor / h;
		
		T * dst1 = dst;
		
		for(uint32_t j=1; j<=qLen; ++j){
			*dst1++ += (T)(amp * sin(phs));
			phs += phaseInc;
		}
	}

	// Extrapolate remaining from [dbqp] symmetry	
	mem::reflectRight(dst, (len >> 1) - 1);
	arr::mirror_dq(--dst, len);
}


TEM void triangleSum(T * dst, uint32_t len, uint32_t hrmLo, uint32_t hrmHi){

	static const double triFactor = 8.0 / (M_PI * M_PI);
	double inc = M_2PI / (double)len;
	uint32_t qLen = len >> 2;
	
	dst++;
	
	hrmLo |= 1;	// next highest odd if even
	
	double factor = hrmLo & 0x2 ? -triFactor : triFactor;
	
	// Calculate first quadrant
	for(uint32_t i = hrmLo; i <= hrmHi; i+=2){

		double h = (double)i;
		double phaseInc = h * inc;
		double phs = phaseInc;
		double amp = factor / (h * h);
		factor = -factor;
		
		T * dst1 = dst;
		
		for(uint32_t j=1; j<=qLen; ++j){
			*dst1++ += (T)(amp * sin(phs));
			phs += phaseInc;
		}
	}

	// Extrapolate remaining from [dbqp] symmetry	
	mem::reflectRight(dst, (len >> 1) - 1);
	arr::mirror_dq(--dst, len);
}

/*	1	1
	2	2
	3	4
	5	8
	9	16	
*/
TEM void multiWave(T * dst, uint32_t len, uint32_t order, void (* func)(T *, uint32_t, uint32_t, uint32_t)){

	dst += len * (order - 1);

	func(dst, len, 1, 1);
	
	uint32_t hrmLo = 2;
	uint32_t hrmHi = 2;

	for(uint32_t o=0; o<order-1; o++){
	
		T * dstPrev = dst;
		dst -= len;
		mem::deepCopy(dst, dstPrev, len);
		func(dst, len, hrmLo, hrmHi);
	
		hrmLo = hrmHi + 1;
		hrmHi <<= 1;
	}
}

TEM void window(T * dst, uint32_t len, WinType::type type){
	switch(type){
		case WinType::Bartlett:			bartlett(dst, len);			break;
		case WinType::Blackman:			blackman(dst, len);			break;
		case WinType::BlackmanHarris:	blackmanHarris(dst, len);	break;
		case WinType::Hamming:			hamming(dst, len);			break;
		case WinType::Hann:				hann(dst, len);				break;
		case WinType::Welch:			welch(dst, len);			break;
		case WinType::Nyquist:			nyquist(dst, len);			break;
		default:						rectangle(dst, len);
	};
}

#define SYM_WIN(period, phs0, eqn) \
	double inc = period / (double)(len);\
	double phs = phs0;\
	T * dst2 = dst + len - 1;\
	*dst++ = (T)eqn;\
	LOOP(len>>1, 1){\
		phs += inc;\
		T val = (T)eqn;\
		*dst++  = val;\
		*dst2-- = val;\
	}
	
TEM void bartlett      (T * dst, uint32_t len){ SYM_WIN(2.   , 0., phs) }
TEM void blackman      (T * dst, uint32_t len){ SYM_WIN(M_2PI, 0., scl::blackman(phs)) }
TEM void blackmanHarris(T * dst, uint32_t len){ SYM_WIN(M_2PI, 0., scl::blackmanHarris(phs)) }
TEM void hamming       (T * dst, uint32_t len){ SYM_WIN(M_2PI, 0., scl::hamming(phs)) }
TEM void hann          (T * dst, uint32_t len){ SYM_WIN(M_2PI, 0., scl::hann(phs)) }
TEM void welch         (T * dst, uint32_t len){ SYM_WIN(2.   ,-1., scl::welch(phs)) }

#undef SYM_WIN

TEM void rectangle(T * dst, uint32_t len){
	for(uint32_t i=0; i<len; ++i) dst[i]=T(1);
}

TEM void nyquist(T * dst, uint32_t len, uint32_t str){
	LOOP(len>>1, str){
		dst[i    ] = (T) 1;
		dst[i+str] = (T)-1;
	}
}



//void TblOp::decay(float * dst, uint32_t len, float order){
//	float final = pow(2.f, order);
//	float lambda = log(final) / (float)len;
//	float time = 1.f;
//	float scale = 1.f / (1.f - final);
//	float offset = -final;
//	*dst++ = 1.f;
//	LOOP(len - 1,
//		*dst++ = (exp(lambda * time) + offset) * scale;
//		time++;
//	)
//}


// Return value from a table containing the first quarter of a sine wave.
// The table size must be a power of two.
//
//	'qsin':		first quarter of a sine wave
//	'bits':		= 30 - (# bits in table)
//	'mask':		= (table size) - 1
//	'phase':	phase of lookup (2^32 is one period of the sine wave)
//
//	'phase' bit format (b = 'bits'):
//	32:			sign bit (0 = positive, 1 = negative)
//	31:			direction bit (0 = forward, 1 = backward)
//	[30, b]:	phase integer part
//	[ b, 0]:	phase fractional part

//inline float TblOp::atQ(float * qsin, uint32_t bits, uint32_t mask, uint32_t phase){
//	union{float f; uint32_t i;} val;
//	uint32_t sign = phase & 0x80000000;
//	uint32_t dir  = phase & 0x40000000;
//
//	phase >>= bits;			// integer part of index
//	
//	// Check direction and reverse index
//	if(dir){
//		phase = -phase; /* index = (len - index);  */
//		if((phase & mask) == 0){
//			val.i = 0x3f800000 | sign;
//			return val.f;
//		}
//	}
//	
//	val.f = qsin[phase & mask];
//	val.i |= sign;	// sign bit
//	return val.f;
//}

// 1 branch
//inline float TblOp::atQ(float * qsin, uint32_t bits, uint32_t mask, uint32_t phase){
//	union{float f; uint32_t i;} val;
//	uint32_t sign = phase & 0x80000000;
//	uint32_t dir  = (phase & 0x40000000) >> 30;	
//
//	phase >>= bits;	// integer part of index
//	
//	/*
//	uint32_t phase2[2];
//	phase2[0] = phase;
//	phase2[1] = -phase;
//	phase = phase2[dir] & mask;
//	*/
//
//	phase = ((phase^-dir) + dir) & mask;	// 2s complement gives flipped phase
//
//	if(phase == 0){
//		val.i = 0x3f800000 & (-dir) | sign;
//		return val.f;
//	}
//	
//	val.f = qsin[phase];
//	val.i |= sign;	// sign bit
//	return val.f;
//}

// 0 branches
// requires n + 1 table with qsin[n] = 1
//inline float TblOp::atQ(float * qsin, uint32_t bits, uint32_t mask, uint32_t phase){
//inline float TblOp::atQ(float * table, uint32_t bits, uint32_t phase){
//	union{float f; uint32_t i;} val;
//	uint32_t sign = phase & 0x80000000;
//	uint32_t dir  = (phase & 0x40000000) >> 30;	// 0 = fwd or 1 = bwd
//
//	//phase >>= bits;	// integer part of index
//	
//	//uint32_t tmp = phase >> 30;
//	//printBinaryUserInt32((phase^-tmp) + tmp, ".", "1");
//	
//	// ulong complementor(ulong v, ulong c){ return v^(-c) + c; }
//	//phase = ((phase^-dir) + dir) & ((mask << 1) | 1);
//	//phase = ((phase^(-dir)) + dir) & mask); // if mask = n<<1 - 1
//
//	phase = (((phase^-dir) + (dir<<bits)) & 0x7fffffff) >> bits;
//
//	val.f = table[phase];
//	val.i |= sign;	// sign bit
//	return val.f;
//}
///*
inline float atQ(const float * src, uint32_t fbits, uint32_t phase){
	uint32_t sign = phase & MaskSign<float>();
	uint32_t dir  = (phase & 0x40000000) >> 30;	// 0 = fwd or 1 = bwd
	Twiddle<float> v(src[(((phase^-dir) + (dir<<fbits)) & 0x7fffffff) >> fbits]);
	v.i |= sign;	// sign bit
	return v.f;
}
//*/
/*
inline float TblOp::atQ(float * table, uint32_t bits, uint32_t phase){
	switch(phase>>30){
	case 0: return  table[( phase & 0x3fffffff) >> bits]; break;
	case 1: return  table[(-phase & 0x7fffffff) >> bits]; break;
	case 2: return -table[( phase & 0x3fffffff) >> bits]; break;
	default:return -table[(-phase & 0x7fffffff) >> bits]; break;
	}
}*/

inline float atH(const float * src, uint32_t bits, uint32_t phase){
	Twiddle<float> v(src[(phase & 0x7fffffff) >> bits]);
	v.i |= phase & MaskSign<float>();				// sign bit
	return v.f;
}

// i: 0 1 2 3 4 5 6 7
// o: 0 1 2 3 0 1 2 3

//TEM inline T at_dq(const T * src, uint32_t len_2, uint32_t i){
//	return i < len_2 ? src[i] : -src[i - len_2];
//}

inline float phaseIncFactor(double framesPerSec){
	//return float(4294967296. / framesPerSec);
	return float(65536. / framesPerSec) * 65536.;
}

} // tbl::
} // gam::

#undef TEM
#undef LOOP

#endif
