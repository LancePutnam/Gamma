#ifndef GAMMA_STRATEGY_H_INC
#define GAMMA_STRATEGY_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

// Function objects representing algorithms
/// \defgroup Strategies Strategies

// \defgroup ipl Interpolation

#include "Gamma/Access.h"
#include "Gamma/Containers.h"
#include "Gamma/ipl.h"
#include "Gamma/scl.h"

namespace gam{

/// Gamma supports several interpolation strategies.  These can be used, for example,
/// to make a delay line whose delay amount is a non-integer number of samples.
/// Julius Smith's <A HREF="https://ccrma.stanford.edu/~jos/pasp/Delay_Line_Interpolation.html">
/// Delay-Line Interpolation page</A>
/// \defgroup ipl Random-access Interpolation Strategies
namespace ipl{

/// Truncating random-access interpolation strategy

/// \ingroup Strategies, ipl
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
	template <class AccessStrategy>
	T operator()(const AccessStrategy& acc, const T * src, index_t iInt, double iFrac, index_t max, index_t min=0) const{
		return src[iInt];
	}

	T operator()(const T * src, index_t iInt, double iFrac, index_t max, index_t min=0) const{
		return (*this)(acc::Wrap(), src, iInt, iFrac, max, min);
	}
};


/// Nearest neighbor random-access interpolation strategy
    
/// \ingroup Strategies, ipl
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
	template <class AccessStrategy>
	T operator()(const AccessStrategy& acc, const T * src, index_t iInt, double iFrac, index_t max, index_t min=0) const{
		return ipl::nearest(
			iFrac,
			src[iInt],
			src[acc.mapP1(iInt+1, max, min)]
		);
	}

	T operator()(const T * src, index_t iInt, double iFrac, index_t max, index_t min=0) const{
		return (*this)(acc::Wrap(), src, iInt, iFrac, max, min);
	}
};


/// Linear random-access interpolation strategy
    
/// \ingroup Strategies, ipl
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
	template <class AccessStrategy>
	T operator()(const AccessStrategy& acc, const T * src, index_t iInt, double iFrac, index_t max, index_t min=0) const{
		return ipl::linear(
			iFrac,
			src[iInt],
			src[acc.mapP1(iInt+1, max, min)]
		);
	}

	T operator()(const T * src, index_t iInt, double iFrac, index_t max, index_t min=0) const{
		return (*this)(acc::Wrap(), src, iInt, iFrac, max, min);
	}
};


/// Cubic random-access interpolation strategy
    
/// \ingroup Strategies, ipl
template <class T>
struct Cubic{

	ipl::Type type() const { return CUBIC; }
	void type(ipl::Type v){}

	/// Return interpolated element from power-of-2 array
	T operator()(const ArrayPow2<T>& a, uint32_t phase) const{
		uint32_t one = a.oneIndex();
		return ipl::cubic(
			a.fraction(phase),
			a.atPhase(phase - one),
			a.atPhase(phase),
			a.atPhase(phase + one),
			a.atPhase(phase + (one<<1))
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
	template <class AccessStrategy>
	T operator()(const AccessStrategy& acc, const T * src, index_t iInt, double iFrac, index_t max, index_t min=0) const{
		return ipl::cubic(
			iFrac,
			src[acc.mapM1(iInt-1, max, min)],
			src[iInt],
			src[acc.mapP1(iInt+1, max, min)],
			src[acc.map  (iInt+2, max, min)]
		);
	}

	T operator()(const T * src, index_t iInt, double iFrac, index_t max, index_t min=0) const{
		return (*this)(acc::Wrap(), src, iInt, iFrac, max, min);
	}
/*
	// TODO: is it worth trying to support strided arrays?
	template <class AccessStrategy>
	T operator()(const AccessStrategy& s, const Array<T>& a, index_t iInt, double iFrac, index_t max, index_t min, index_t str) const{		
		return ipl::cubic(
			iFrac,
			a[s.map(iInt-str, max, min)],
			a[iInt],
			a[s.map(iInt+str, max, min)],
			a[s.map(iInt+(str<<1), max, min)]
		);
	}

	T operator()(const Array<T>& a, index_t iInt, double iFrac, index_t max, index_t min, index_t str) const{
		return (*this)(acc::Wrap(), a, iInt,iFrac, max,min,str);
	}
*/
};


/// Allpass random-access interpolation strategy
    
/// \ingroup Strategies, ipl
template <class T>
struct AllPass{

	AllPass(T prev=0): prev(prev){}

	ipl::Type type() const { return ALLPASS; }
	void type(ipl::Type v){}

	/// Return interpolated element from power-of-2 array
	T operator()(const ArrayPow2<T>& a, uint32_t phase) const{
		return ipl::allpass(
			a.fraction(phase), 
			a.atPhase(phase), 
			a.atPhase(phase + a.oneIndex()),
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
	template <class AccessStrategy>
	T operator()(const AccessStrategy& acc, const T * src, index_t iInt, double iFrac, index_t max, index_t min=0) const{
		return ipl::allpass(
			iFrac,
			src[iInt],
			src[acc.mapP1(iInt+1, max, min)],
			prev
		);
	}

	T operator()(const T * src, index_t iInt, double iFrac, index_t max, index_t min=0) const{
		return (*this)(acc::Wrap(), src, iInt, iFrac, max, min);
	}
	
	mutable T prev;
};


/// Dynamically switchable random-access interpolation strategy
    
/// \ingroup Strategies, ipl
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
	T operator()(const AccessStrategy& s, const Array<T>& a, index_t iInt, double iFrac, index_t max, index_t min=0) const{
		switch(mType){
			case ROUND:		return round	(s,a,iInt,iFrac,max,min);
			case LINEAR:	return linear	(s,a,iInt,iFrac,max,min);
			case CUBIC:		return cubic	(s,a,iInt,iFrac,max,min);
			case ALLPASS:	return allpass	(s,a,iInt,iFrac,max,min);
			default:		return trunc	(s,a,iInt,iFrac,max,min);
		}
	}

	T operator()(const Array<T>& a, index_t iInt, double iFrac, index_t max, index_t min=0) const{
		switch(mType){
			case ROUND:		return round	(a,iInt,iFrac,max,min);
			case LINEAR:	return linear	(a,iInt,iFrac,max,min);
			case CUBIC:		return cubic	(a,iInt,iFrac,max,min);
			case ALLPASS:	return allpass	(a,iInt,iFrac,max,min);
			default:		return trunc	(a,iInt,iFrac,max,min);
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
    
    /// \ingroup Strategies, iplSeq
	template <uint32_t N, class T>
	struct Base{
		Base(const T& v=0){ set(v); }
	
		/// Push a new value onto sequence
		void push(T va){ for(uint32_t i=N-1; i>0; --i) v[i]=v[i-1]; v[0]=va; }
		
		/// Set sequence history to value
		void set(T va){ for(uint32_t i=0; i<N; ++i) v[i]=va; }	
		
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
    
    /// \ingroup Strategies, iplSeq
	template <class T>
	struct Linear : public Base<2,T>{
		using Base<2,T>::v;
		Linear(const T& v=0): Base<2,T>(v){}
		T operator()(float f) const { return ipl::linear(f, v[1], v[0]); }
	};

	/// Cubic sequence interpolation strategy
    
    /// \ingroup Strategies, iplSeq
	template <class T>
	struct Cubic : public Base<4,T>{
		using Base<4,T>::v;
		Cubic(const T& v=0): Base<4,T>(v){}
		T operator()(float f) const { return ipl::cubic(f, v[3], v[2], v[1], v[0]); }
		T val() const { return v[1]; }
		void val(const T& va){ v[1]=va; }
	};

	/// Cosine sequence interpolation strategy
    
    /// \ingroup Strategies, iplSeq
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
//	T operator()(T v, T max, T min);					// float tap post increment check
//	void reset();										// reset internal state, if any
    
/// \defgroup phsInc Phase Increment Strategies
namespace phsInc{

	// Increment and clip to closed-open interval [min, max)
	template <class T>
	T incClip(T v, T inc, T max, T min){
		T res = v + inc;
		if(res >= max) return v;
		if(res <  min) return min;
		return res;
	}

	/// Loop waveform indefinitely.
            
    /// \ingroup Strategies, phsInc
	struct Loop{
		void reset(){}
	
		uint32_t operator()(uint32_t& pos, uint32_t inc){ return pos+=inc; }
		bool done(uint32_t pos) const { return false; }
		
		template <class T>
		T operator()(T v, T inc, T max, T min){ return scl::wrap(v+inc, max, min); }
	};


	/// Play waveform one cycle, then hold at the end. A one-shot.
    
    /// \ingroup Strategies, phsInc
	struct OneShot{
		void reset(){}
	
		uint32_t operator()(uint32_t& pos, uint32_t inc){
			uint32_t prev = pos;
			pos += inc;
			if(~pos & prev & 0x80000000) pos = 0xffffffff;
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
    
    /// \ingroup Strategies, phsInc
	struct NShot{
		NShot(){ number(1); reset(); }
		
		void reset(){ mCount=0; }

		uint32_t operator()(uint32_t& pos, uint32_t inc){
			uint32_t prev = pos;
			pos += inc;
			
			// Check MSB goes from 1 to 0
			// TODO: works only for positive increments and non-zero mNumber
			if((~pos & prev) & 0x80000000){
				if(++mCount >= mNumber) pos = 0xffffffff;
			}
			return pos;
		}
		
		bool done(uint32_t pos) const { return (mCount >= mNumber) && (pos == 0xffffffff); }
		
		template <class T>
		T operator()(T v, T inc, T max, T min){
			v += inc;
			if(v >= max || v < min) ++mCount;
			return mCount < mNumber ? scl::wrap(v, max, min) : incClip(v, inc, max, min);
		}
		
		// Set number of repetitions
		NShot& number(uint32_t v){ mNumber=v; return *this; }

	private:
		uint32_t mNumber;
		uint32_t mCount;
	};


	/// Play waveform forward, backwards, forward, etc.  Like Wrap, loops indefinitely.
    
    /// \ingroup Strategies, phsInc
	struct PingPong{
		PingPong(): dir(0){}
	
		void reset(){ dir=0; }
	
		uint32_t operator()(uint32_t& pos, uint32_t inc){
			uint32_t prev = pos;
			pos += dir ? -inc : inc;
			if(~pos & prev & 0x80000000) dir^=1;
			return pos;
		}

		bool done(uint32_t pos) const { return false; }
		
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
    
    /// \ingroup Strategies, phsInc
	struct Rhythm{
	
		Rhythm(){
			pattern(0,0);
			reset();
		}
	
		void reset(){ mPhase=0; mIndex=0; }

		uint32_t operator()(uint32_t& pos, uint32_t inc){
			uint32_t prev = mPhase;
			mPhase += inc;

			// Check MSB goes from 1 to 0
			// TODO: works only for positive increments
			if((~mPhase & prev) & 0x80000000){
				if(++mIndex >= mSize) mIndex=0;
			}

			uint32_t bit = (mPattern >> (mSize-1 - mIndex)) & 1UL;
			if(bit){
				pos = mPhase;
			}
			return pos;
		}
		
		Rhythm& pattern(uint32_t bits, uint16_t size){
			mPattern=bits;
			mSize=size;
			return *this;
		}
		
		Rhythm& pattern(const char* bits){
			mSize = strlen(bits);
			mPattern = 0;
			for(int i=0; i<mSize; ++i){
				mPattern |= (bits[i]!='.') << (mSize-1-i);
			}
			return *this;
		}
		
	private:
		uint32_t mPhase;
		uint32_t mPattern;
		uint16_t mIndex;
		uint16_t mSize;
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
