#ifndef GAMMA_STRATEGY_H_INC
#define GAMMA_STRATEGY_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include "Access.h"
#include "Containers.h"
#include "ipl.h"
#include "scl.h"

namespace gam{


namespace ipl{

/// Truncating random-access interpolation strategy
struct Trunc{

	/// Return element from power-of-2 array
	template <class T>
	T operator()(const ArrayPow2<T>& a, uint32_t phase) const{
		return a.atPhase(phase);
	}

	template <class T, class AccessStrategy>
	T operator()(const AccessStrategy& s, const Array<T>& a, uint32_t iInt, double iFrac, uint32_t max, uint32_t min=0) const{
		return a[iInt];
	}

	template <class T>
	T operator()(const Array<T>& a, uint32_t iInt, double iFrac, uint32_t max, uint32_t min=0) const{
		return (*this)(acc::Wrap(), a,iInt,iFrac, max,min);
	}
};


/// Nearest neighbor random-access interpolation strategy
struct Round{

	/// Return element from power-of-2 array
	template <class T>
	T operator()(const ArrayPow2<T>& a, uint32_t phase) const{

		// accessing normally truncates, so add half fraction to round
		return a.atPhase(phase + (a.oneIndex()>>1));
	}
	
	template <class T, class AccessStrategy>
	T operator()(const AccessStrategy& s, const Array<T>& a, uint32_t iInt, double iFrac, uint32_t max, uint32_t min=0) const{
		return ipl::nearest(
			iFrac,
			a[iInt],
			a[s.mapP1(iInt+1, max, min)]
		);
	}

	template <class T>
	T operator()(const Array<T>& a, uint32_t iInt, double iFrac, uint32_t max, uint32_t min=0) const{
		return (*this)(acc::Wrap(), a, iInt, iFrac, max, min);
	}
};


/// Linear random-access interpolation strategy
struct Linear{

	/// Return element from power-of-2 array
	template <class T>
	T operator()(const ArrayPow2<T>& a, uint32_t phase) const{
		return ipl::linear(
			a.fraction(phase),
			a.atPhase(phase),
			a.atPhase(phase + a.oneIndex())
		);
	}

	template <class T, class AccessStrategy>
	T operator()(const AccessStrategy& s, const Array<T>& a, uint32_t iInt, double iFrac, uint32_t max, uint32_t min=0) const{		
		return ipl::linear(
			iFrac,
			a[iInt],
			a[s.mapP1(iInt+1, max, min)]
		);
	}

	template <class T>
	T operator()(const Array<T>& a, uint32_t iInt, double iFrac, uint32_t max, uint32_t min=0) const{
//		return ipl::linear(
//			iFrac,
//			a[iInt],
//			a[scl::wrapOnce(iInt + 1, max, min)]
//		);
		return (*this)(acc::Wrap(), a, iInt, iFrac, max, min);
	}	

};


/// Cubic random-access interpolation strategy
struct Cubic{

	/// Return element from power-of-2 array
	template <class T>
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
	
	template <class T, class AccessStrategy>
	T operator()(const AccessStrategy& s, const Array<T>& a, uint32_t iInt, double iFrac, uint32_t max, uint32_t min=0) const{		
		return ipl::cubic(
			iFrac,
			a[s.mapM1(iInt-1, max, min)],
			a[iInt],
			a[s.mapP1(iInt+1, max, min)],
			a[s.map  (iInt+2, max, min)]
		);
	}
	
	template <class T>
	T operator()(const Array<T>& a, uint32_t iInt, double iFrac, uint32_t max, uint32_t min=0) const{
		return (*this)(acc::Wrap(), a, iInt,iFrac, max,min);
	}
};


/// Allpass random-access interpolation strategy
template <class T>
struct AllPass{

	AllPass(T prev=0): prev(prev){}

	/// Return element from power-of-2 array
	T operator()(const ArrayPow2<T>& a, uint32_t phase) const{
		return ipl::allpass(					// Standard fraction
		//return Ipol::allpassFixed(				// Fixed fractional delay
			a.fraction(phase), 
			a.atPhase(phase), 
			a.atPhase(phase + a.oneIndex()),
			prev
		);
	}
	
	template <class AccessStrategy>
	T operator()(const AccessStrategy& s, const Array<T>& a, uint32_t iInt, double iFrac, uint32_t max, uint32_t min=0) const{		
		return ipl::allpass(
			iFrac,
			a[iInt],
			a[s.mapP1(iInt+1, max, min)],
			prev
		);
	}

	T operator()(const Array<T>& a, uint32_t iInt, double iFrac, uint32_t max, uint32_t min=0) const{
		return (*this)(acc::Wrap(), a, iInt,iFrac, max,min);
	}
	
	mutable T prev;
};


} // ipl::



// Sequence interpolation strategies

// interface:
// 
// method		desc
// ()			return value at fraction
// push			push new value into interpolation window 
namespace iplSeq{

	template <uint32_t N, class T>
	struct Base{
		Base(const T& v=0){ set(v); }
	
		void push(T va){ for(uint32_t i=N-1; i>0; --i) v[i]=v[i-1]; v[0]=va; }
		void set(T va){ for(uint32_t i=0; i<N; ++i) v[i]=va; }		
		T val() const { return v[0]; }
		void val(const T& va){ v[0]=va; }
			
		T v[N];	// value buffer, 0 == newest, N-1 == oldest
	};

	template <class T>
	struct Cubic : public Base<4, T>{
		using Base<4, T>::v;
		Cubic(const T& v=0): Base<4, T>(v){}
		T operator()(float f) const { return ipl::cubic(f, v[3], v[2], v[1], v[0]); }
		T val() const { return v[1]; }
		void val(const T& va){ v[1]=va; }
	};
	
	template <class T>
	struct Linear : public Base<2, T>{
		using Base<2, T>::v;
		Linear(const T& v=0): Base<2,T>(v){}
		T operator()(float f) const { return ipl::linear(f, v[1], v[0]); }
	};

	template <class T>
	struct Cosine : public Base<2, T>{
		using Base<2, T>::v;
		Cosine(const T& v=0): Base<2,T>(v){}
		T operator()(float f) const { return ipl::linear(scl::warpSinUU(f), v[1], v[0]); }
	};

} // iplSeq::



// Read tap strategies.

// The expected strategy interface is:
//	void operator()(uint32_t& pos, uint32_t inc);	// integer tap increment
//	bool done(uint32_t pos);						// integer tap done reading
//	T operator()(T v, T max, T min);				// float tap post increment check
namespace tap{

	struct Wrap{
		void operator()(uint32_t& pos, uint32_t inc){ pos += inc; }
		bool done(uint32_t pos){ return false; }
		
		template <class T>
		T operator()(T v, T max, T min){ return scl::wrap(v, max, min); }
	};

	struct Clip{
		void operator()(uint32_t& pos, uint32_t inc){
			uint32_t prev = pos;
			pos += inc;
			if(~pos & prev & 0x80000000) pos = 0xffffffff;

			// pos3:	1101	inc = 0001
			// pos2:	1110
			// pos1:	1111
			// pos0:	0000	msb goes from 1 to 0
		}
		bool done(uint32_t pos){ return pos == 0xffffffff; }
		
		template <class T>
		T operator()(T v, T max, T min){ return scl::clip(v, max, min); }
	};

	template <uint32_t N>
	struct Rep{

		Rep(): count(0){}

		void operator()(uint32_t& pos, uint32_t inc){
			uint32_t prev = pos;
			pos += inc;
			if(~pos & prev & 0x80000000) ++count < N ?: pos = 0xffffffff;
		}
		
		bool done(uint32_t pos){ return (count >= N) && (pos == 0xffffffff); }
		
		template <class T>
		T operator()(T v, T max, T min){
			if(v >= max || v < min) count++;
			return count < N ? scl::wrap(v, max, min) : scl::clip(v, max, min);
		}
		
		uint32_t count;
	};

} // tap::



// Functors for warping values in domain/range [0,1]
namespace warp{

	struct None{
		template <class T>
		T operator()(T v){ return v; }
	};

	struct S{
		template <class T>
		T operator()(T v){ return scl::warpSinUU(v); }
	};

} // warp::

} // gam::

#endif
