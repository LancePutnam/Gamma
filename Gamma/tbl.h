#ifndef GAMMA_TBL_H_INC
#define GAMMA_TBL_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include "Gamma/arr.h"
#include "Gamma/mem.h"
#include "Gamma/scl.h"
#include "Gamma/Constants.h"
#include "Gamma/Types.h"

namespace gam{

/// Window types
enum WindowType{
	BARTLETT,			/**< Bartlett (Triangle) */
	BLACKMAN,			/**< Blackman */
	BLACKMAN_HARRIS,	/**< Blackman-Harris */
	BLACKMAN_NUTTALL,	/**< Blackman-Nuttall */
	FLATTOP,			/**< Flat-Top */
	HAMMING,			/**< Hamming */
	HANN,				/**< von Hann */
	NUTTALL,			/**< Nuttall */
	WELCH,				/**< Welch */
	NYQUIST,			/**< Nyquist */
	RECTANGLE			/**< Rectangle (no window) */
};


/// Waveform types
enum WaveformType{
	SINE,				/**< Sine wave */
	COSINE,				/**< Cosine wave */
	TRIANGLE,			/**< Triangle wave */
	PARABOLIC,			/**< Parabolic wave */
	SQUARE,				/**< Square wave */
	SAW,				/**< Saw wave */
	IMPULSE				/**< Impulse wave */
};


/// Returns human readable string of window type
inline static const char * toString(gam::WindowType v){
	#define CS(name) case name: return #name;
	switch(v){
		CS(BARTLETT) CS(BLACKMAN) CS(BLACKMAN_HARRIS) CS(HAMMING) CS(HANN)
		CS(WELCH) CS(NYQUIST) CS(RECTANGLE)
		default: return "Unknown";
	}
	#undef CS
}



/// Add sine wave to array

/// \param[out] dst		destination array
/// \param[in] len		length of array
/// \param[in] cycles	number of cycles of sine wave, must be integer for periodic waves
/// \param[in] amp		amplitude of sine wave
/// \param[in] phs		phase of sine wave, in [0,1]
template <class T>
void addSine(T * dst, unsigned len, double cycles=1, double amp=1, double phs=0);

/// Add sine wave to array

/// \param[out] dst		destination array
/// \param[in] cycles	number of cycles of sine wave, must be integer for periodic waves
/// \param[in] amp		amplitude of sine wave
/// \param[in] phs		phase of sine wave, in [0,1]
template <class T, class Alloc, template<class,class> class ArrayType>
void inline addSine(
	ArrayType<T,Alloc>& dst,
	double cycles=1, double amp=1, double phs=0
){
	addSine(&dst[0], dst.size(), cycles, amp, phs);
}


/// Add harmonic series to array with specified amplitudes

/// \param[out] dst		destination array
/// \param[in] len		length of destination array
/// \param[in] amps		amplitudes of harmonic series, size must be numh
/// \param[in] numh		total number of harmonics
/// \param[in] hmul		harmonic number multiplication factor
/// \param[in] hshf		harmonic number shift amount
/// \param[in] hphs		phase of sine wave, in [0,1]
template <class T, class A>
void inline addSines(
	T * dst, unsigned len, const A * amps, int numh,
	double hmul=1, double hshf=1, double hphs=0
){
	for(int i=0; i<numh; ++i){
		if(A(0) != amps[i]) addSine(dst,len, i*hmul+hshf, amps[i], hphs);
	}
}

/// Add harmonic series to array with specified amplitudes

/// \param[out] dst		destination array
/// \param[in] amps		amplitudes of harmonic series, size must be numh
/// \param[in] numh		total number of harmonics
/// \param[in] hmul		harmonic number multiplication factor
/// \param[in] hshf		harmonic number shift amount
/// \param[in] hphs		phase of sine wave, in [0,1]
template <class T, class Alloc, template<class,class> class ArrayType, class A>
void inline addSines(
	ArrayType<T,Alloc>& dst, const A * amps, int numh,
	double hmul=1, double hshf=1, double hphs=0
){
	addSines(&dst[0],dst.size(),amps,numh,hmul,hshf,hphs);
}


/// Add harmonics to array with specified amplitudes and harmonic numbers

/// \param[out] dst		destination array
/// \param[in] len		length of destination array
/// \param[in] amps		harmonic amplitudes of series, size must be numh
/// \param[in] cycs		harmonic numbers of series, size must be numh
/// \param[in] numh		total number of harmonics
/// \param[in] hphs		phase of sine wave, in [0,1]
template <class T, class A, class C>
void addSines(
	T * dst, unsigned len, const A * amps, const C * cycs, int numh, double hphs=0
){
	for(int i=0; i<numh; ++i) addSine(dst,len,cycs[i],amps[i],hphs);
}

/// Add harmonics to array with specified amplitudes and harmonic numbers

/// \param[out] dst		destination array
/// \param[in] amps		harmonic amplitudes of series, size must be numh
/// \param[in] cycs		harmonic numbers of series, size must be numh
/// \param[in] numh		total number of harmonics
/// \param[in] hphs		phase of sine wave, in [0,1]
template <class T, class Alloc, template<class,class> class ArrayType, class A, class C>
void inline addSines(
	ArrayType<T,Alloc>& dst, const A * amps, const C * cycs, int numh, double hphs=0)
{
	addSines(&dst[0],dst.size(), amps,cycs,numh,hphs);
}


/// Add sine waves to array using inverse power law for amplitudes

/// \tparam InvPower	amplitudes will be set to 1 / h^InvPower
/// \param[out] dst		destination array
/// \param[in] len		length of destination array
/// \param[in] numh		total number of harmonics
/// \param[in] hmul		harmonic number multiplication factor
/// \param[in] hshf		harmonic number shift amount
/// \param[in] amp		overall amplitude scaling factor
/// \param[in] hphs		phase of (sine) harmonics, in [0,1]
/// \param[in] wphs		phase of composite waveform, in [0,1]
/// \param[in] cycles	number of cycles of wave, must be integer for periodic waves
template <int InvPower, class T>
void addSinesPow(
	T * dst, unsigned len, int numh,
	double hmul=1, double hshf=1, double amp=1, double hphs=0, double wphs=0, double cycles=1
);

/// Add sine waves to array using inverse power law for amplitudes

/// \tparam InvPower	amplitudes will be set to 1 / h^InvPower
/// \param[out] dst		destination array
/// \param[in] numh		total number of harmonics
/// \param[in] hmul		harmonic number multiplication factor
/// \param[in] hshf		harmonic number shift amount
/// \param[in] amp		overall amplitude scaling factor
/// \param[in] hphs		phase of (sine) harmonics, in [0,1]
/// \param[in] wphs		phase of composite waveform, in [0,1]
/// \param[in] cycles	number of cycles of wave, must be integer for periodic waves
template <int InvPower, class T, class Alloc, template<class,class> class ArrayType>
inline void addSinesPow(
	ArrayType<T,Alloc>& dst, int numh,
	double hmul=1, double hshf=1, double amp=1, double hphs=0, double wphs=0, double cycles=1
){
	addSinesPow<InvPower>(&dst[0], dst.size(), numh,hmul,hshf,amp,hphs,wphs,cycles);
}


/// Add predefined waveform to array

/// The produced waveforms are not normalized; the fundamental always has a 
/// unit amplitude.
/// \param[out] dst		destination array
/// \param[in] len		length of destination array
/// \param[in] type		waveform type
/// \param[in] numh		total number of harmonics
/// \param[in] amp		amplitude of waveform
/// \param[in] phs		phase of waveform, in [0,1]
/// \param[in] hshf		harmonic number shift amount
/// \param[in] cycles	number of cycles of wave, must be integer for periodic waves
template <class T>
void addWave(
	T * dst, unsigned len, gam::WaveformType type,
	int numh=32, double amp=1, double phs=0, double hshf=1, double cycles=1
);

/// Add predefined waveform to array

/// The produced waveforms are not normalized; the fundamental always has a 
/// unit amplitude.
/// \param[out] dst		destination array
/// \param[in] type		waveform type
/// \param[in] numh		total number of harmonics
/// \param[in] amp		amplitude of waveform
/// \param[in] phs		phase of waveform, in [0,1]
/// \param[in] hshf		harmonic number shift amount
/// \param[in] cycles	number of cycles of wave, must be integer for periodic waves
template <class T, class Alloc, template<class,class> class ArrayType>
void inline addWave(
	ArrayType<T,Alloc>& dst, gam::WaveformType type,
	int numh=32, double amp=1, double phs=0, double hshf=1, double cycles=1
){
	addWave(&dst[0],dst.size(), type,numh,amp,phs,hshf,cycles);
}


/// Get Fourier series normalization constant for a waveform
template<gam::WaveformType W> double normConstant();



/// Table functions
namespace tbl{


/// Fills array with one period of a cosine wave.
template<class T>
void cosine(T * dst, unsigned len);

/// Fills array with one period of a sine wave.
template<class T>
void sine(T * dst, unsigned len);

/// Fills array with arbitrary phase and length sinusoid.
template<class T>
void sinusoid(T * dst, unsigned len, double phase, double periods);

/// Sums band-limited impulse wave into multi-wavetable array

/// The waveform includes harmonics in the range [hrmLo, hrmHi].
/// The amplitude of the waveform will not be normalized.
/// The ideal waveform shape is [4, -1, 0, -1, 0, -1, 0, -1 ]
template<class T>
void multiImpulse(T * dst, unsigned len, unsigned hrmLo, unsigned hrmHi);

/// Sums band-limited saw wave into multi-wavetable array

/// The waveform includes harmonics in the range [hrmLo, hrmHi].
/// The ideal waveform shape is [1, 0.75, 0.5, 0.25, 0, -0.25, -0.5, -0.75]
template<class T>
void multiSaw(T * dst, unsigned len, unsigned hrmLo, unsigned hrmHi);

/// Sums band-limited square wave into multi-wavetable array

/// The waveform includes harmonics in the range [hrmLo, hrmHi].
///	The ideal waveform shape is [ 1, 1, 1, 1, -1, -1, -1, -1].
template<class T>
void multiSquare(T * dst, unsigned len, unsigned hrmLo, unsigned hrmHi);

/// Sums band-limited triangle wave into multi-wavetable array
	
/// The waveform includes harmonics in the range [hrmLo, hrmHi].
///	The ideal waveform shape is [ 0, 0.5, 1, 0.5, 0, -0.5, -1, -0.5].
template<class T>
void multiTriangle(T * dst, unsigned len, unsigned hrmLo, unsigned hrmHi);

/// Create multi-wavetable
template<class T>
void multiWave(T * dst, unsigned len, unsigned order, void (* func)(T *, unsigned, unsigned, unsigned));

/// Returns maximum number of harmonics that will fit in array.
inline unsigned maxHarmonics(unsigned len){ return len>>1; }

/// Fills array with specified window type
template<class T>
void window			(T * dst, unsigned len, WindowType type);

template<class T>
void bartlett		(T * dst, unsigned len); ///< Fills array with Bartlett window

template<class T>
void blackman		(T * dst, unsigned len); ///< Fills array with Blackman window

template<class T>
void blackmanHarris	(T * dst, unsigned len); ///< Fills array with Blackman-Harris window

template<class T>
void blackmanNuttall(T * dst, unsigned len); ///< Fills array with Blackman-Nuttall window

template<class T>
void flatTop		(T * dst, unsigned len); ///< Fills array with Flat-Top window

template<class T>
void hamming		(T * dst, unsigned len); ///< Fills array with Hamming window

template<class T>
void hann			(T * dst, unsigned len); ///< Fills array with von Hann window

template<class T>
void nuttall		(T * dst, unsigned len); ///< Fills array with Nuttall window

template<class T>
void welch			(T * dst, unsigned len); ///< Fills array with Welch window

template<class T>
void rectangle		(T * dst, unsigned len); ///< Fills array with Rectangle window

template<class T>
void nyquist		(T * dst, unsigned len, unsigned str=1); ///< Fills array with Nyquist window



/// Get value from a power-of-two array.

///	'src':		waveform \n
///	'fbits':	number of bits in fractional part of phase \n
///	'phase':	fixed-point phase of lookup (full period is [0, 2^32))
template <class T>
T at(const T * src, uint32_t fbits, uint32_t phase);

/// Set value in a power-of-two array.

///	'dst':		waveform \n
///	'fbits':	number of bits in fractional part of phase \n
///	'phase':	fixed-point phase of lookup (full period is [0, 2^32))	
template <class T>
void put(T * dst, uint32_t fbits, uint32_t phase, T value);


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

#define LOOP(n,s) for(unsigned i=0; i<n; i+=s)

template<class T>
void cosine(T * dst, unsigned len){
	double inc = M_2PI / (double)len;
	double phs = inc;
	len >>= 1;
	
	T * dst2 = dst + len;
	
	*dst++  = T( 1);
	*dst2++ = T(-1);

	len -= 1;
	LOOP(len, 1){
		T val = T(cos(phs));
		*dst++  =  val;
		*dst2++ = -val;
		phs += inc;
	}
}

template<class T>
void sine(T * dst, unsigned len){
	double inc = M_2PI / (double)len;
	double phs = inc;
	len >>= 1;
	
	T * dst2 = dst + len;
	
	*dst++  = T(0);
	*dst2++ = T(0);
	
	--len;
	LOOP(len, 1){
		T val = sin(phs);
		*dst++  =  val;
		*dst2++ = -val;
		phs += inc;
	}
}

// VERY accurate, but not so fast
template<class T>
void sinusoid(T * dst, unsigned len, double phase, double periods){
	double inc = M_2PI * periods / len;
	for(unsigned i=0; i<len; ++i){
		*dst++ = sin(inc * i + phase);
	}
}


template<class T>
void multiImpulse(T * dst, unsigned len, unsigned hrmLo, unsigned hrmHi){
	double inc = M_2PI / (double)len;
	unsigned hLen = len >> 1;
	
	for(unsigned k = hrmLo; k <= hrmHi; ++k){
		double phaseInc = (double)k * inc;
		double phs = 0.;
		
		T * dst1 = dst;
		
		LOOP(hLen+1, 1){
			*dst1++ += T(cos(phs));
			phs += phaseInc;
		}
	}
	
	// Extrapolate remaining from [db] symmetry
	mem::reflectRight(dst + 1, len - 1);
}

template<class T>
void multiSaw(T * dst, unsigned len, unsigned hrmLo, unsigned hrmHi){

	static const double sawFactor = 2.0 / M_PI;
	double inc = M_2PI / (double)len;
	unsigned hLen = len >> 1;
	
	dst++;
	
	for(unsigned i = hrmLo; i <= hrmHi; ++i){
		double h = (double)i;
		double phaseInc = h * inc;
		double phs = phaseInc;
		double amp = sawFactor / h;
		
		T * dst1 = dst;
		
		for(unsigned j=1; j<hLen; ++j){
			*dst1++ += T(amp * sin(phs));
			phs += phaseInc;
		}
	}
	
	// Extrapolate remaining from [dp] symmetry
	arr::mirror_dp(dst, len-1);	
}

template<class T>
void multiSquare(T * dst, unsigned len, unsigned hrmLo, unsigned hrmHi){

	static const double sqrFactor = 4.0 / M_PI;
	double inc = M_2PI / (double)len;
	unsigned qLen = len >> 2;
	
	dst++;
	
	hrmLo |= 1;	// next highest odd if even
	
	// Calculate first quadrant
	for(unsigned i = hrmLo; i <= hrmHi; i+=2){

		double h = (double)i;
		double phaseInc = h * inc;
		double phs = phaseInc;
		double amp = sqrFactor / h;
		
		T * dst1 = dst;
		
		for(unsigned j=1; j<=qLen; ++j){
			*dst1++ += T(amp * sin(phs));
			phs += phaseInc;
		}
	}

	// Extrapolate remaining from [dbqp] symmetry	
	mem::reflectRight(dst, (len >> 1) - 1);
	arr::mirror_dq(--dst, len);
}


template<class T>
void multiTriangle(T * dst, unsigned len, unsigned hrmLo, unsigned hrmHi){

	static const double triFactor = 8.0 / (M_PI * M_PI);
	double inc = M_2PI / (double)len;
	unsigned qLen = len >> 2;
	
	dst++;
	
	hrmLo |= 1;	// next highest odd if even
	
	double factor = hrmLo & 0x2 ? -triFactor : triFactor;
	
	// Calculate first quadrant
	for(unsigned i = hrmLo; i <= hrmHi; i+=2){

		double h = (double)i;
		double phaseInc = h * inc;
		double phs = phaseInc;
		double amp = factor / (h * h);
		factor = -factor;
		
		T * dst1 = dst;
		
		for(unsigned j=1; j<=qLen; ++j){
			*dst1++ += T(amp * ::sin(phs));
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
template<class T>
void multiWave(T * dst, unsigned len, unsigned order, void (* func)(T *, unsigned, unsigned, unsigned)){

	dst += len * (order - 1);

	func(dst, len, 1, 1);

	unsigned hrmLo = 2;
	unsigned hrmHi = 2;

	for(unsigned o=0; o<order-1; o++){
		T * dstPrev = dst;
		dst -= len;
		mem::deepCopy(dst, dstPrev, len);
		func(dst, len, hrmLo, hrmHi);
	
		hrmLo = hrmHi + 1;
		hrmHi <<= 1;
	}
}

template<class T>
void window(T * dst, unsigned len, WindowType type){
	switch(type){
		case BARTLETT:			bartlett(dst, len);			break;
		case BLACKMAN:			blackman(dst, len);			break;
		case BLACKMAN_HARRIS:	blackmanHarris(dst, len);	break;
		case BLACKMAN_NUTTALL:	blackmanNuttall(dst, len);	break;
		case FLATTOP:			flatTop(dst, len);			break;
		case HAMMING:			hamming(dst, len);			break;
		case HANN:				hann(dst, len);				break;
		case NUTTALL:			nuttall(dst, len);			break;
		case WELCH:				welch(dst, len);			break;
		case NYQUIST:			nyquist(dst, len);			break;
		default:				rectangle(dst, len);
	};
}

#define SYM_WIN(period, phs0, eqn) \
	double inc = period / (double)(len);\
	double phs = phs0;\
	T * dst2 = dst + len - 1;\
	*dst++ = T(eqn);\
	LOOP(len>>1, 1){\
		phs += inc;\
		T val = T(eqn);\
		*dst++  = val;\
		*dst2-- = val;\
	}
	
template<class T>
void bartlett		(T * dst, unsigned len){ SYM_WIN(2.   , 0., phs) }

template<class T>
void blackman		(T * dst, unsigned len){ SYM_WIN(M_2PI, 0., scl::blackman(phs)) }

template<class T>
void blackmanHarris	(T * dst, unsigned len){ SYM_WIN(M_2PI, 0., scl::blackmanHarris(phs)) }

template<class T>
void blackmanNuttall(T * dst, unsigned len){ SYM_WIN(M_2PI, 0., scl::blackmanNuttall(phs)) }

template<class T>
void flatTop		(T * dst, unsigned len){ SYM_WIN(M_2PI, 0., scl::flatTop(phs)) }

template<class T>
void hamming		(T * dst, unsigned len){ SYM_WIN(M_2PI, 0., scl::hamming(phs)) }

template<class T>
void hann			(T * dst, unsigned len){ SYM_WIN(M_2PI, 0., scl::hann(phs)) }

template<class T>
void nuttall		(T * dst, unsigned len){ SYM_WIN(M_2PI, 0., scl::nuttall(phs)) }

template<class T>
void welch			(T * dst, unsigned len){ SYM_WIN(2.   ,-1., scl::welch(phs)) }

#undef SYM_WIN

template<class T> void rectangle(T * dst, unsigned len){
	for(unsigned i=0; i<len; ++i) dst[i]=T(1);
}

template<class T>
void nyquist(T * dst, unsigned len, unsigned str){
	LOOP(len, str*2){
		dst[(i+0)*str] = T( 1);
		dst[(i+1)*str] = T(-1);
	}
}

#undef LOOP

template <class T>
inline T at(const T * src, uint32_t fbits, uint32_t phase){
	return src[phase >> fbits];
}

template <class T>
inline void put(T * dst, uint32_t fbits, uint32_t phase, T value){
	dst[phase >> fbits] = value;
}

inline float atH(const float * src, uint32_t bits, uint32_t phs){
	Twiddle<float> v(src[(phs & 0x7fffffff) >> bits]);
	v.i |= phs & MaskSign<float>();				// sign bit
	return v.f;
}

inline float atQ(const float * src, uint32_t fbits, uint32_t phs){
	uint32_t sign = phs & MaskSign<float>();
	uint32_t dir = (phs >> 30) & 1; // 0 = fwd or 1 = bwd
	Twiddle<float> v(src[(((phs^-int32_t(dir)) + (dir<<fbits)) & 0x7fffffff) >> fbits]);
	v.i |= sign;	// sign bit
	return v.f;
}

/*inline float atQ(const float * src, uint32_t bits, uint32_t phs){
	switch(phs>>30){
	case 0: return  src[( phs & 0x3fffffff) >> bits]; break;
	case 1: return  src[(-phs & 0x7fffffff) >> bits]; break;
	case 2: return -src[( phs & 0x3fffffff) >> bits]; break;
	default:return -src[(-phs & 0x7fffffff) >> bits]; break;
	}
}*/

// i: 0 1 2 3 4 5 6 7
// o: 0 1 2 3 0 1 2 3

//template<class T> inline T at_dq(const T * src, uint32_t len_2, uint32_t i){
//	return i < len_2 ? src[i] : -src[i - len_2];
//}

inline float phaseIncFactor(double framesPerSec){
	//return float(4294967296. / framesPerSec);
	return float(65536. / framesPerSec) * 65536.0f;
}

} // tbl::


template<WaveformType W> inline double normConstant(){ return 1.; }
template<> inline double normConstant<TRIANGLE	>(){ return 8/(M_PI*M_PI); }
template<> inline double normConstant<PARABOLIC	>(){ return 2/ M_PI; }
template<> inline double normConstant<SQUARE	>(){ return 4/ M_PI; }
template<> inline double normConstant<SAW		>(){ return 2/ M_PI; }

namespace{
	// This enables support for complex-valued tables
	template<class T>
	T getSin(double p);

	template<class T>
	T getSin(double p){
		return T(sin(p));
	}

	template<>
	Complex<float> getSin<Complex<float> >(double p){
		return Complex<float>(float(cos(p)), float(sin(p)));
	}

	template<>
	Complex<double> getSin<Complex<double> >(double p){
		return Complex<double>(cos(p), sin(p));
	}
};

template <class T>
void addSine(T * dst, unsigned len, double cycles, double amp, double phs){
	double f = cycles/len;
	for(unsigned i=0; i<len; ++i){
		double p = (f*i + phs)*M_2PI;
		dst[i] += getSin<T>(p) * amp;
	}
}

template <int InvPower, class T>
void addSinesPow(
	T * dst, unsigned len, int numh,
	double hmul, double hshf, double amp, double hphs, double wphs, double cycles
){
	using std::pow;
	const double inc1 = (M_2PI / len) * cycles;

	for(int j=0; j<numh; ++j){
		const double h = j*hmul + hshf;
		if(InvPower && 0==h) continue; // don't bother with DC component

		double A = amp;		
		switch(InvPower){
		case 0: break;
		case 1: A /= h; break;
		case 2: A /= h*h; break;
		case 3: A /= h*h*h; break;
		default:A *= pow(h, -InvPower);
		}

		const double P = (hphs + h*wphs) * M_2PI;
		const double inch = inc1 * h;

		for(unsigned i=0; i<len; ++i){
			dst[i] += A*getSin<T>(inch*i + P);
		}
	}
}


template <class T>
void addWave(
	T * dst, unsigned len, WaveformType type,
	int numh, double amp, double phs, double hshf, double cycles
){
//	static const double ctri = normConstant<TRIANGLE>();
//	static const double csaw = normConstant<SAW>();
//	static const double csqr = normConstant<SQUARE>();
	static const double ctri = 1;
	static const double csaw = 1;
	static const double csqr = 1;

	switch(type){
	case SINE:		addSine(dst,len, hshf*cycles,amp,phs); break;
	case COSINE:	addSine(dst,len, hshf*cycles,amp,phs+0.25); break;
	case TRIANGLE:	addSinesPow<2>(dst,len, numh,2,hshf,amp*ctri,0.25,phs-0.25,cycles); break;
	case PARABOLIC:	addSinesPow<2>(dst,len, numh,1,hshf,amp*csaw,0.25,phs,cycles); break;
	case SQUARE:	addSinesPow<1>(dst,len, numh,2,hshf,amp*csqr,0.00,phs,cycles); break;
	case SAW:		addSinesPow<1>(dst,len, numh,1,hshf,amp*csaw,0.00,phs,cycles); break;
	case IMPULSE:	addSinesPow<0>(dst,len, numh,1,hshf,amp     ,0.25,phs,cycles); break;
	default:;
	}
}


} // gam::

#endif
