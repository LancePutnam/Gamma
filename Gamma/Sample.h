#ifndef GAMMA_SAMPLE_H_INC
#define GAMMA_SAMPLE_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include <cmath> // abs, max
#include <cstdint> // int32_t
#include <vector>

namespace gam{

/// Convert sample from one type to another

/// The input/output ranges are (expected to be)
///
///   double: [-1, 1]\n
///   float:  [-1, 1]\n
///   short:  [-32767, 32767]\n
///   char:   [-127, 127]\n
///
/// Signed n-bit int types should be encoded in [-2^(n-1)+1, 2^(n-1)-1].
/// The minimum value should be dropped to have an equal # samples around 0.
template <class To, class From> To sampleTo(From v);

template<> inline double sampleTo<double,double>(double v){ return v; }
template<> inline double sampleTo<double,float>(float v){ return double(v); }
template<> inline double sampleTo<double,int32_t>(int32_t v){ return double(v)*(1./2147483647.); }
template<> inline double sampleTo<double,short>(short v){ return double(v)*(1./32767.); }
template<> inline double sampleTo<double,char>(char v){ return double(v)*(1./127.); }
template<> inline float sampleTo<float,double>(double v){ return float(v); }
template<> inline float sampleTo<float,float>(float v){ return v; }
template<> inline float sampleTo<float,int32_t>(int32_t v){ return sampleTo<double>(v); }
template<> inline float sampleTo<float,short>(short v){ return float(v)*(1.f/32767.f); }
template<> inline float sampleTo<float,char>(char v){ return float(v)*(1.f/127.f); }
namespace{
	template <class I, class R>
	inline I roundTo(R v){ return I(v + (v<R(0)?R(-0.5):R(0.5))); }
}
template<> inline int32_t sampleTo<int32_t,double>(double v){ return roundTo<int32_t>(v*2147483647.); }
template<> inline int32_t sampleTo<int32_t,float>(float v){ return sampleTo<int32_t>(double(v)); }
template<> inline int32_t sampleTo<int32_t,int32_t>(int32_t v){ return v; }
template<> inline int32_t sampleTo<int32_t,short>(short v){ return int32_t(v)*65538; }
template<> inline int32_t sampleTo<int32_t,char>(char v){ return int32_t(v)*16909320; }
template<> inline short sampleTo<short,double>(double v){ return roundTo<short>(v*32767.); }
template<> inline short sampleTo<short,float>(float v){ return roundTo<short>(v*32767.f); }
template<> inline short sampleTo<short,int32_t>(int32_t v){ return v/65536; }
template<> inline short sampleTo<short,short>(short v){ return v; }
template<> inline short sampleTo<short,char>(char v){ return short(v)*258; }
template<> inline char sampleTo<char,double>(double v){ return roundTo<char>(v*127.); }
template<> inline char sampleTo<char,float>(float v){ return roundTo<char>(v*127.f); }
template<> inline char sampleTo<char,int32_t>(int32_t v){ return v/16777216; }
template<> inline char sampleTo<char,short>(short v){ return v/256; }
template<> inline char sampleTo<char,char>(char v){ return v; }

/// Sample data
template <class T>
struct Sample{
	typedef T value_type;

	std::vector<value_type> data;	///< Sample data array
	double frameRate = 1.;			///< Frame rate (AKA sample rate)
	int channels = 1;				///< Channels per frame
	float gain = 1.;				///< Nominal playback gain (e.g. for equal loudness)
	bool interleaved=true;			///< Whether samples are interleaved in array

	/// Get number of (multi-channel) frames
	int frames() const { return data.size()/channels; }

	/// Get total length, in seconds
	float length() const { return frames()/frameRate; }

	/// Returns whether sample array is empty
	bool empty() const { return data.empty(); }

	/// Get peak value (extremum)
	T peak() const {
		T mx = T(0);
		for(auto& v : data) mx = std::max<T>(std::abs(v),mx);
		return mx;
	}

	/// Fit all samples within (unit) range

	/// \returns normalization factor
	///
	float normalize(float peakGain = 1.){
		float nrm = sampleTo<float>(peak());
		if(nrm > 0.f){
			nrm = peakGain/nrm;
			for(auto& v : data) v = value_type(float(v)*nrm);
		} else {
			nrm = 1.f;
		}
		return nrm;
	}
};

} // gam::

#endif
