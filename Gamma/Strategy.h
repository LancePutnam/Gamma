#ifndef GAMMA_STRATEGY_H_INC
#define GAMMA_STRATEGY_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include <limits> // numeric_limits
#include "Gamma/Access.h"
#include "Gamma/Containers.h"
#include "Gamma/ipl.h"
#include "Gamma/scl.h"

namespace gam{

/// Function objects representing algorithms

/// \defgroup Strategy

// \defgroup ipl Interpolation

/// Gamma supports several interpolation strategies.  These can be used, for example,
/// to make a delay line whose delay amount is a non-integer number of samples.
/// Julius Smith's <A HREF="https://ccrma.stanford.edu/~jos/pasp/Delay_Line_Interpolation.html">
/// Delay-Line Interpolation page</A>
/// \defgroup ipl Random-access Interpolation Strategies

namespace ipl{

/// Truncating random-access interpolation strategy

/// \ingroup Strategy, ipl
template <class T>
struct Trunc{

	ipl::Type type() const { return TRUNC; }
	void type(ipl::Type v){}

	/// Return interpolated element from power-of-2 array
	T operator()(const ArrayPow2<T>& a, uint32_t phase) const{
		return a.atPhase(phase);
	}

	/// Return interpolated element from array

	/// \tparam AccessStrategy	access strategy type (\sa access)
	///
	/// \param[in] acc			access strategy
	/// \param[in] src			source array
	/// \param[in] iInt			integer part of index
	/// \param[in] iFrac		fractional part of index, in [0, 1)
	/// \param[in] max			maximum index for accessing
	/// \param[in] min			minimum index for accessing
	/// \param[in] str			array stride (post-multiplier on indices)
	template <class AccessStrategy>
	T operator()(const AccessStrategy& acc, const T * src, index_t iInt, double iFrac, index_t max, index_t min=0, index_t str=1) const{
		return src[iInt*str];
	}

	T operator()(const T * src, index_t iInt, double iFrac, index_t max, index_t min=0, index_t str=1) const{
		return (*this)(acc::Wrap(), src, iInt, iFrac, max, min, str);
	}
};


/// Nearest neighbor random-access interpolation strategy

/// \ingroup Strategy, ipl
template <class T>
struct Round{

	ipl::Type type() const { return ROUND; }
	void type(ipl::Type v){}

	/// Return interpolated element from power-of-2 array
	T operator()(const ArrayPow2<T>& src, uint32_t phase) const{
		// accessing normally truncates, so add half fraction to round
		return src.atPhase(phase + (src.oneIndex()>>1));
	}

	/// Return interpolated element from array

	/// \tparam AccessStrategy	access strategy type (\sa access)
	///
	/// \param[in] acc			access strategy
	/// \param[in] src			source array
	/// \param[in] iInt			integer part of index
	/// \param[in] iFrac		fractional part of index, in [0, 1)
	/// \param[in] max			maximum index for accessing
	/// \param[in] min			minimum index for accessing
	/// \param[in] str			array stride (post-multiplier on indices)
	template <class AccessStrategy>
	T operator()(const AccessStrategy& acc, const T * src, index_t iInt, double iFrac, index_t max, index_t min=0, index_t str=1) const{
		return ipl::nearest(
			iFrac,
			src[iInt*str],
			src[acc.mapP1(iInt+1, max, min)*str]
		);
	}

	T operator()(const T * src, index_t iInt, double iFrac, index_t max, index_t min=0, index_t str=1) const{
		return (*this)(acc::Wrap(), src, iInt, iFrac, max, min, str);
	}
};


/// Sum of two nearest neighbors random-access interpolation strategy

/// This is equivalent to truncating interpolation with a 1-zero low-pass
/// filter.
/// \ingroup Strategy, ipl
template <class T>
struct Mean2{

	ipl::Type type() const { return MEAN2; }
	void type(ipl::Type v){}

	/// Return interpolated element from power-of-2 array
	T operator()(const ArrayPow2<T>& a, uint32_t phase) const{
		return (a.atPhase(phase) + a.atPhase(phase + a.oneIndex()))*0.5;
	}

	/// Return interpolated element from array

	/// \tparam AccessStrategy	access strategy type (\sa access)
	///
	/// \param[in] acc			access strategy
	/// \param[in] src			source array
	/// \param[in] iInt			integer part of index
	/// \param[in] iFrac		fractional part of index, in [0, 1)
	/// \param[in] max			maximum index for accessing
	/// \param[in] min			minimum index for accessing
	/// \param[in] str			array stride (post-multiplier on indices)
	template <class AccessStrategy>
	T operator()(const AccessStrategy& acc, const T * src, index_t iInt, double iFrac, index_t max, index_t min=0, index_t str=1) const{
		return (src[iInt*str] + src[acc.mapP1(iInt+1, max, min)*str])*0.5;
	}

	T operator()(const T * src, index_t iInt, double iFrac, index_t max, index_t min=0, index_t str=1) const{
		return (*this)(acc::Wrap(), src, iInt, iFrac, max, min, str);
	}
};


/// Linear random-access interpolation strategy

/// \ingroup Strategy, ipl
template <class T>
struct Linear{

	ipl::Type type() const { return LINEAR; }
	void type(ipl::Type v){}

	/// Return interpolated element from power-of-2 array
	T operator()(const ArrayPow2<T>& a, uint32_t phase) const{
		return ipl::linear(
			a.fraction(phase),
			a.atPhase(phase),
			a.atPhase(phase + a.oneIndex())
		);
	}

	/// Return interpolated element from array

	/// \tparam AccessStrategy	access strategy type (\sa access)
	///
	/// \param[in] acc			access strategy
	/// \param[in] src			source array
	/// \param[in] iInt			integer part of index
	/// \param[in] iFrac		fractional part of index, in [0, 1)
	/// \param[in] max			maximum index for accessing
	/// \param[in] min			minimum index for accessing
	/// \param[in] str			array stride (post-multiplier on indices)
	template <class AccessStrategy>
	T operator()(const AccessStrategy& acc, const T * src, index_t iInt, double iFrac, index_t max, index_t min=0, index_t str=1) const{
		return ipl::linear(
			iFrac,
			src[iInt*str],
			src[acc.mapP1(iInt+1, max, min)*str]
		);
	}

	T operator()(const T * src, index_t iInt, double iFrac, index_t max, index_t min=0, index_t str=1) const{
		return (*this)(acc::Wrap(), src, iInt, iFrac, max, min, str);
		// TODO: wrapping access maybe not correct for one-shot playback
	}
};


/// Cubic random-access interpolation strategy

/// \ingroup Strategy, ipl
template <class T>
struct Cubic{

	ipl::Type type() const { return CUBIC; }
	void type(ipl::Type v){}

	/// Return interpolated element from power-of-2 array
	T operator()(const ArrayPow2<T>& a, uint32_t phase) const{
		//*
		const auto oneIndex = a.oneIndex();
		return ipl::cubic(
			a.fraction(phase),
			a.atPhase(phase - oneIndex),
			a.atPhase(phase),
			a.atPhase(phase + oneIndex),
			a.atPhase(phase + (oneIndex<<1))
		);//*/
		/* faster
		unsigned i = a.index(phase);
		unsigned m = a.size() - 1;
		return ipl::cubic(
			a.fraction(phase),
			a[(i-1)&m],
			a[ i ],
			a[(i+1)&m],
			a[(i+2)&m]
		);//*/
	}

	/// Return interpolated element from array

	/// \tparam AccessStrategy	access strategy type (\sa access)
	///
	/// \param[in] acc			access strategy
	/// \param[in] src			source array
	/// \param[in] iInt			integer part of index
	/// \param[in] iFrac		fractional part of index, in [0, 1)
	/// \param[in] max			maximum index for accessing
	/// \param[in] min			minimum index for accessing
	/// \param[in] str			array stride (post-multiplier on indices)
	template <class AccessStrategy>
	T operator()(const AccessStrategy& acc, const T * src, index_t iInt, double iFrac, index_t max, index_t min=0, index_t str=1) const{
		return ipl::cubic(
			iFrac,
			src[acc.mapM1(iInt-1, max, min)*str],
			src[iInt*str],
			src[acc.mapP1(iInt+1, max, min)*str],
			src[acc.map  (iInt+2, max, min)*str]
		);
	}

	T operator()(const T * src, index_t iInt, double iFrac, index_t max, index_t min=0, index_t str=1) const{
		return (*this)(acc::Wrap(), src, iInt, iFrac, max, min, str);
	}
};


/// Allpass random-access interpolation strategy

/// \ingroup Strategy, ipl
template <class T>
struct AllPass{

	AllPass(T prev=0): prev(prev){}

	ipl::Type type() const { return ALLPASS; }
	void type(ipl::Type v){}

	/// Return interpolated element from power-of-2 array
	T operator()(const ArrayPow2<T>& a, uint32_t phase) const{
		const auto oneIndex = a.oneIndex();
		return ipl::allpass(
			a.fraction(phase),
			a.atPhase(phase - oneIndex),
			a.atPhase(phase), 
			a.atPhase(phase + oneIndex),
			prev
		);
	}

	/// Return interpolated element from array

	/// \tparam AccessStrategy	access strategy type (\sa access)
	///
	/// \param[in] acc			access strategy
	/// \param[in] src			source array
	/// \param[in] iInt			integer part of index
	/// \param[in] iFrac		fractional part of index, in [0, 1)
	/// \param[in] max			maximum index for accessing
	/// \param[in] min			minimum index for accessing
	/// \param[in] str			array stride (post-multiplier on indices)
	template <class AccessStrategy>
	T operator()(const AccessStrategy& acc, const T * src, index_t iInt, double iFrac, index_t max, index_t min=0, index_t str=1) const{
		return ipl::allpass(
			iFrac,
			src[acc.mapM1(iInt-1, max, min)*str],
			src[iInt*str],
			src[acc.mapP1(iInt+1, max, min)*str],
			prev
		);
	}

	T operator()(const T * src, index_t iInt, double iFrac, index_t max, index_t min=0, index_t str=1) const{
		return (*this)(acc::Wrap(), src, iInt, iFrac, max, min, str);
	}
	
	mutable T prev;
};


/// Dynamically switchable random-access interpolation strategy

/// \ingroup Strategy, ipl
template <class T>
struct Switchable{

	Switchable(): mType(TRUNC){}
	
	ipl::Type type() const { return mType; }
	void type(ipl::Type v){ mType=v; }

	/// Return interpolated element from power-of-2 array
	T operator()(const ArrayPow2<T>& a, uint32_t phase) const{		
		switch(mType){
			case ROUND:		return round	(a, phase);
			case LINEAR:	return linear	(a, phase);
			case CUBIC:		return cubic	(a, phase);
			case ALLPASS:	return allpass	(a, phase);
			default:		return trunc	(a, phase);
		}
	}

	template <class AccessStrategy>
	T operator()(const AccessStrategy& s, const Array<T>& a, index_t iInt, double iFrac, index_t max, index_t min=0, index_t str=1) const{
		switch(mType){
			case ROUND:		return round	(s,a,iInt,iFrac,max,min,str);
			case LINEAR:	return linear	(s,a,iInt,iFrac,max,min,str);
			case CUBIC:		return cubic	(s,a,iInt,iFrac,max,min,str);
			case ALLPASS:	return allpass	(s,a,iInt,iFrac,max,min,str);
			default:		return trunc	(s,a,iInt,iFrac,max,min,str);
		}
	}

	T operator()(const Array<T>& a, index_t iInt, double iFrac, index_t max, index_t min=0, index_t str=1) const{
		switch(mType){
			case ROUND:		return round	(a,iInt,iFrac,max,min,str);
			case LINEAR:	return linear	(a,iInt,iFrac,max,min,str);
			case CUBIC:		return cubic	(a,iInt,iFrac,max,min,str);
			case ALLPASS:	return allpass	(a,iInt,iFrac,max,min,str);
			default:		return trunc	(a,iInt,iFrac,max,min,str);
		}
	}

protected:
	ipl::Type mType;
	Trunc<T> trunc;
	Round<T> round;
	Linear<T> linear;
	Cubic<T> cubic;
	AllPass<T> allpass;
};

} // ipl::



/// Sequence (stream-based) interpolation strategies

// interface:
// 
// method		desc
// ()			return value at fraction
// push			push new value into interpolation window

/// \defgroup iplSeq Sequence Interpolation Strategies
namespace iplSeq{

	/// Base class for sequence interpolation strategies

	/// \ingroup Strategy, iplSeq
	template <unsigned N, class T>
	struct Base{
		Base(const T& v=0){ set(v); }
	
		/// Push a new value onto sequence
		void push(T va){ for(unsigned i=N-1; i>0; --i) v[i]=v[i-1]; v[0]=va; }
		
		/// Set sequence history to value
		void set(T va){ for(unsigned i=0; i<N; ++i) v[i]=va; }	
		
		/// Get current sequence value
		T val() const { return v[0]; }
		
		/// Set current sequence value
		void val(const T& va){ v[0]=va; }
			
		T v[N];	///< Value buffer, 0 is newest, N-1 is oldest
	};

	/// Non-interpolating sequence interpolation strategy

	/// \ingroup iplSeq
	template <class T>
	struct None : public Base<1,T>{
		using Base<1,T>::v;
		None(const T& v=0): Base<1,T>(v){}
		T operator()(float f) const { return v[0]; }
	};

	/// Truncating sequence interpolation strategy

	/// \ingroup iplSeq
	template <class T>
	struct Trunc : public Base<2,T>{
		using Base<2,T>::v;
		Trunc(const T& v=0): Base<2,T>(v){}
		T operator()(float f) const { return v[1]; }
	};

	/// Round half up sequence interpolation strategy

	/// \ingroup iplSeq
	template <class T>
	struct Round : public Base<2,T>{
		using Base<2,T>::v;
		Round(const T& v=0): Base<2,T>(v){}
		T operator()(float f) const { return ipl::nearest(f, v[1], v[0]); }
	};

	/// Linear sequence interpolation strategy

	/// \ingroup Strategy, iplSeq
	template <class T>
	struct Linear : public Base<2,T>{
		using Base<2,T>::v;
		Linear(const T& v=0): Base<2,T>(v){}
		T operator()(float f) const { return ipl::linear(f, v[1], v[0]); }
	};

	/// Cubic sequence interpolation strategy

	/// \ingroup Strategy, iplSeq
	template <class T>
	struct Cubic : public Base<4,T>{
		using Base<4,T>::v;
		Cubic(const T& v=0): Base<4,T>(v){}
		T operator()(float f) const { return ipl::cubic(f, v[3], v[2], v[1], v[0]); }
		T val() const { return v[1]; }
		void val(const T& va){ v[1]=va; }
	};

	/// Cosine sequence interpolation strategy

	/// \ingroup Strategy, iplSeq
	template <class T>
	struct Cosine : public Base<2,T>{
		using Base<2,T>::v;
		Cosine(const T& v=0): Base<2,T>(v){}
		T operator()(float f) const { return ipl::linear(scl::mapSinUU(f), v[1], v[0]); }
	};

} // iplSeq::



/// Phase increment strategies

// The expected interface is:
//	uint32_t operator()(uint32_t& pos, uint32_t inc);	// fixed-point tap increment
//	bool done(uint32_t pos);							// fixed-point tap done reading
//	T operator()(T v, T inc, T max, T min);				// float tap increment
//	void reset();										// reset internal state, if any

/// \defgroup phsInc Phase Increment Strategies
namespace phsInc{

	// Increment and clip to closed-open interval [min, max)
	template <class T>
	T incClip(T v, T inc, T max, T min){
		T res = v + inc;
		if(res >= max) return v < max ? v : scl::nextAfter(max, min);
		if(res <  min) return min;
		return res;
	}

	/// Loop waveform indefinitely.

	/// \ingroup Strategy, phsInc
	struct Loop{
		void reset(){}
	
		uint32_t operator()(uint32_t& pos, uint32_t inc){ return pos+=inc; }
		bool done(uint32_t /*pos*/) const { return false; }
		
		template <class T>
		T operator()(T v, T inc, T max, T min){ return scl::wrap(v+inc, max, min); }
	};


	/// Play then hold waveform
	struct Pulse{

		Pulse(){
			pulse(2,4);
			reset();
		}

		/// Set width and period of pulse

		/// \param[in] width	pulse width in cycles
		/// \param[in] period	total pulse period in cycles
		Pulse& pulse(uint32_t width, uint32_t period){
			mWidth = width;
			mPeriod = period;
			return *this;
		}

		void reset(){ mPhase=0; mCycle=0; setOn(); }

		uint32_t operator()(uint32_t& pos, uint32_t inc){
			uint32_t prev = mPhase;
			mPhase += inc;

			// Check for phase wrap
			if((prev > mPhase) ^ (inc >> 31)){
				//if(++mCycle >= mPeriod) mCycle=0;
				mCycle = (mCycle+1) % mPeriod;
				setOn();
			}

			if(mOn){
				pos = mPhase;
			}

			return pos;
		}

		bool done(uint32_t /*pos*/) const { return false; }

	private:
		uint32_t mPhase;
		uint32_t mCycle;
		uint32_t mWidth;
		uint32_t mPeriod;
		uint8_t mOn;
		void setOn(){ mOn = mCycle<mWidth; }
	};


	/// Play waveform one cycle, then hold at the end. A one-shot.

	/// \ingroup Strategy, phsInc
	struct OneShot{
		void reset(){}
	
		uint32_t operator()(uint32_t& pos, uint32_t inc){
			uint32_t prev = pos;
			pos += inc;
			if(~pos & prev & 0x80000000) pos = 0xffffffff;
			// correct version
			/*if((prev>pos) ^ (inc>>31)){
				pos = (inc>>31) ? 0 : -1;
			}*/
			// pos3:	1101	inc = 0001
			// pos2:	1110
			// pos1:	1111
			// pos0:	0000	msb goes from 1 to 0
			return pos;
		}
		bool done(uint32_t pos) const { return pos == 0xffffffff; }
		
		template <class T>
		T operator()(T v, T inc, T max, T min){ return incClip(v,inc,max,min); }
	};


	/// Repeat waveform a fixed number of times, then hold at the end.

	/// \ingroup Strategy, phsInc
	struct NShot{
		NShot(){ repeats(1); reset(); }
		
		void reset(){ mCount=0; }

		uint32_t operator()(uint32_t& pos, uint32_t inc){
			uint32_t prev = pos;
			pos += inc;
			
			// Check MSB goes from 1 to 0
			// TODO: works only for positive increments and non-zero mNumber
			if((~pos & prev) & 0x80000000){
				if(++mCount >= mRepeats) pos = 0xffffffff;
			}
			return pos;
		}
		
		bool done(uint32_t pos) const { return (mCount >= mRepeats) && (pos == 0xffffffff); }
		
		template <class T>
		T operator()(T v, T inc, T max, T min){
			v += inc;
			if(v >= max || v < min) ++mCount;
			return mCount < mRepeats ? scl::wrap(v, max, min) : incClip(v, inc, max, min);
		}
		
		/// Set number of repetitions
		NShot& repeats(uint32_t v){ mRepeats=v; return *this; }
		NShot& number(uint32_t v){ mRepeats=v; return *this; }

	private:
		uint32_t mRepeats;
		uint32_t mCount;
	};


	/// Play waveform forward, backwards, forward, etc.  Like Wrap, loops indefinitely.

	/// \ingroup Strategy, phsInc
	struct PingPong{
		PingPong(): dir(0){}
	
		void reset(){ dir=0; }
	
		uint32_t operator()(uint32_t& pos, uint32_t inc){
			uint32_t prev = pos;
			pos += dir ? -int32_t(inc) : inc;
			if(~pos & prev & 0x80000000) dir^=1;
			return pos;
		}

		bool done(uint32_t /*pos*/) const { return false; }
		
		template <class T>
		T operator()(T v, T inc, T max, T min){		
			v += dir ? -inc : inc;
			long n;
			v = scl::fold(v, n, max, min);
			dir ^= n!=0;
			return v;
		}
		
		uint32_t dir;
	};


	/// Plays and holds waveform according to binary repeating pattern.

	/// \ingroup Strategy, phsInc
	struct Rhythm{
	
		Rhythm(){
			pattern(2,2);
			reset();
		}
	
		void reset(){ mPhase=0; mIndex=0; setOn(); }

		uint32_t operator()(uint32_t& pos, uint32_t inc){
			uint32_t prev = mPhase;
			mPhase += inc;

			// Check for phase wrap
			if((prev>mPhase) ^ (inc>>31)){
				mIndex = (mIndex+1) % mSize;
				setOn();
			}

			if(mOn){
				pos = mPhase;
			}

			return pos;
		}

		bool done(uint32_t /*pos*/) const { return false; }
		
		Rhythm& pattern(uint64_t bits, uint8_t size){
			mPattern=bits;
			mSize=size;
			return *this;
		}
		
		/// Specify on/off pattern as string, e.g., "/./....."
		Rhythm& pattern(const char* bits, char offChar='.'){
			mSize = uint8_t(strlen(bits));
			mPattern = 0;
			for(int i=0; i<mSize; ++i){
				mPattern |= (bits[i]!=offChar) << (mSize-1-i);
			}
			return *this;
		}

		/// Set pattern of play then hold cycles

		/// \param[in] pulseWidth	pulse width in cycles
		/// \param[in] length		total pulse length in cycles
		Rhythm& pulse(uint8_t pulseWidth, uint8_t length){
			mSize = length;
			auto ones = std::numeric_limits<uint64_t>::max();
			mPattern = (ones>>(64-mSize)) & (ones<<(length-pulseWidth));			
			return *this;
		}

	private:
		uint64_t mPattern;
		uint32_t mPhase;
		uint8_t mOn;
		uint8_t mIndex;
		uint8_t mSize;
		void setOn(){ mOn = (mPattern >> (mSize-1 - mIndex)) & uint64_t(1); }
	};

	typedef Loop Wrap;
	typedef OneShot Clip;
	typedef NShot Rep;
	typedef PingPong Fold;
	typedef Rhythm Pat;

} // phsInc::

namespace tap = phsInc;

} // gam::
#endif
