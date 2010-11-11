#ifndef GAMMA_STRATEGY_H_INC
#define GAMMA_STRATEGY_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include "Gamma/Access.h"
#include "Gamma/Containers.h"
#include "Gamma/ipl.h"
#include "Gamma/scl.h"

namespace gam{


namespace ipl{

enum{
	TRUNC=0,
	ROUND,
	LINEAR,
	CUBIC,
	ALLPASS
};


/// Truncating random-access interpolation strategy
template <class T>
struct Trunc{

	int type() const { return TRUNC; }
	void type(int v){}

	/// Return element from power-of-2 array
	T operator()(const ArrayPow2<T>& a, uint32_t phase) const{
		return a.atPhase(phase);
	}

	template <class AccessStrategy>
	T operator()(const AccessStrategy& s, const Array<T>& a, index_t iInt, double iFrac, index_t max, index_t min=0) const{
		return a[iInt];
	}

	T operator()(const Array<T>& a, index_t iInt, double iFrac, index_t max, index_t min=0) const{
		return (*this)(acc::Wrap(), a,iInt,iFrac, max,min);
	}
};


/// Nearest neighbor random-access interpolation strategy
template <class T>
struct Round{

	int type() const { return ROUND; }
	void type(int v){}

	/// Return element from power-of-2 array
	T operator()(const ArrayPow2<T>& a, uint32_t phase) const{

		// accessing normally truncates, so add half fraction to round
		return a.atPhase(phase + (a.oneIndex()>>1));
	}
	
	template <class AccessStrategy>
	T operator()(const AccessStrategy& s, const Array<T>& a, index_t iInt, double iFrac, index_t max, index_t min=0) const{
		return ipl::nearest(
			iFrac,
			a[iInt],
			a[s.mapP1(iInt+1, max, min)]
		);
	}

	T operator()(const Array<T>& a, index_t iInt, double iFrac, index_t max, index_t min=0) const{
		return (*this)(acc::Wrap(), a, iInt, iFrac, max, min);
	}
};


/// Linear random-access interpolation strategy
template <class T>
struct Linear{

	int type() const { return LINEAR; }
	void type(int v){}

	/// Return element from power-of-2 array
	T operator()(const ArrayPow2<T>& a, uint32_t phase) const{
		return ipl::linear(
			a.fraction(phase),
			a.atPhase(phase),
			a.atPhase(phase + a.oneIndex())
		);
	}

	template <class AccessStrategy>
	T operator()(const AccessStrategy& s, const Array<T>& a, index_t iInt, double iFrac, index_t max, index_t min=0) const{		
		return ipl::linear(
			iFrac,
			a[iInt],
			a[s.mapP1(iInt+1, max, min)]
		);
	}

	T operator()(const Array<T>& a, index_t iInt, double iFrac, index_t max, index_t min=0) const{
//		return ipl::linear(
//			iFrac,
//			a[iInt],
//			a[scl::wrapOnce(iInt + 1, max, min)]
//		);
		return (*this)(acc::Wrap(), a, iInt, iFrac, max, min);
	}	

};


/// Cubic random-access interpolation strategy
template <class T>
struct Cubic{

	int type() const { return CUBIC; }
	void type(int v){}

	/// Return element from power-of-2 array
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
	
	template <class AccessStrategy>
	T operator()(const AccessStrategy& s, const Array<T>& a, index_t iInt, double iFrac, index_t max, index_t min=0) const{		
		return ipl::cubic(
			iFrac,
			a[s.mapM1(iInt-1, max, min)],
			a[iInt],
			a[s.mapP1(iInt+1, max, min)],
			a[s.map  (iInt+2, max, min)]
		);
	}

	T operator()(const Array<T>& a, index_t iInt, double iFrac, index_t max, index_t min=0) const{
		return (*this)(acc::Wrap(), a, iInt,iFrac, max,min);
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
template <class T>
struct AllPass{

	AllPass(T prev=0): prev(prev){}

	int type() const { return ALLPASS; }
	void type(int v){}

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
	T operator()(const AccessStrategy& s, const Array<T>& a, index_t iInt, double iFrac, index_t max, index_t min=0) const{		
		return ipl::allpass(
			iFrac,
			a[iInt],
			a[s.mapP1(iInt+1, max, min)],
			prev
		);
	}

	T operator()(const Array<T>& a, index_t iInt, double iFrac, index_t max, index_t min=0) const{
		return (*this)(acc::Wrap(), a, iInt,iFrac, max,min);
	}
	
	mutable T prev;
};


//
template <class T>
struct Any{

	Any(): mType(TRUNC){}
	
	int type() const { return mType; }
	void type(int v){ mType=v; }

	/// Return element from power-of-2 array
	T operator()(const ArrayPow2<T>& a, uint32_t phase) const{		
		switch(mType){
			case ROUND:		return round(a, phase);
			case LINEAR:	return linear(a, phase);
			case CUBIC:		return cubic(a, phase);
			case ALLPASS:	return allpass(a, phase);
			default:		return trunc(a, phase);
		}
	}

	template <class AccessStrategy>
	T operator()(const AccessStrategy& s, const Array<T>& a, index_t iInt, double iFrac, index_t max, index_t min=0) const{
		switch(mType){
			case ROUND:		return round(s,a,iInt,iFrac,max,min);
			case LINEAR:	return linear(s,a,iInt,iFrac,max,min);
			case CUBIC:		return cubic(s,a,iInt,iFrac,max,min);
			case ALLPASS:	return allpass(s,a,iInt,iFrac,max,min);
			default:		return trunc(s,a,iInt,iFrac,max,min);
		}
	}

	T operator()(const Array<T>& a, index_t iInt, double iFrac, index_t max, index_t min=0) const{
		switch(mType){
			case ROUND:		return round(a,iInt,iFrac,max,min);
			case LINEAR:	return linear(a,iInt,iFrac,max,min);
			case CUBIC:		return cubic(a,iInt,iFrac,max,min);
			case ALLPASS:	return allpass(a,iInt,iFrac,max,min);
			default:		return trunc(a,iInt,iFrac,max,min);
		}
	}

protected:
	int mType;
	Trunc<T> trunc;
	Round<T> round;
	Linear<T> linear;
	Cubic<T> cubic;
	AllPass<T> allpass;
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
			
		T v[N];	// value buffer, 0 is newest, N-1 is oldest
	};

	template <class T>
	struct Trunc : public Base<1,T>{
		using Base<1,T>::v;
		Trunc(const T& v=0): Base<1,T>(v){}
		T operator()(float f) const { return v[0]; }
	};
	
	template <class T>
	struct Linear : public Base<2,T>{
		using Base<2,T>::v;
		Linear(const T& v=0): Base<2,T>(v){}
		T operator()(float f) const { return ipl::linear(f, v[1], v[0]); }
	};

	template <class T>
	struct Cubic : public Base<4,T>{
		using Base<4,T>::v;
		Cubic(const T& v=0): Base<4,T>(v){}
		T operator()(float f) const { return ipl::cubic(f, v[3], v[2], v[1], v[0]); }
		T val() const { return v[1]; }
		void val(const T& va){ v[1]=va; }
	};

	template <class T>
	struct Cosine : public Base<2,T>{
		using Base<2,T>::v;
		Cosine(const T& v=0): Base<2,T>(v){}
		T operator()(float f) const { return ipl::linear(scl::warpSinUU(f), v[1], v[0]); }
	};

} // iplSeq::



// Read tap strategies.

// The expected strategy interface is:
//	void operator()(uint32_t& pos, uint32_t inc);	// fixed-point tap increment
//	bool done(uint32_t pos);						// fixed-point tap done reading
//	T operator()(T v, T max, T min);				// float tap post increment check
//	void reset();									// reset internal state, if any
namespace tap{

	struct Clip{
		void reset(){}
	
		uint32_t& operator()(uint32_t& pos, uint32_t inc){
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
		T operator()(T v, T inc, T max, T min){ return scl::clip(v+inc, max, min); }
	};

	struct Fold{
		Fold(): dir(0){}
	
		void reset(){ dir=0; }
	
		uint32_t& operator()(uint32_t& pos, uint32_t inc){
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

	template <uint32_t N>
	struct Rep{
		Rep(): count(0){}
		
		void reset(){ count=0; }

		uint32_t& operator()(uint32_t& pos, uint32_t inc){
			uint32_t prev = pos;
			pos += inc;
			if(~pos & prev & 0x80000000) ++count < N ?: pos = 0xffffffff;
			return pos;
		}
		
		bool done(uint32_t pos) const { return (count >= N) && (pos == 0xffffffff); }
		
		template <class T>
		T operator()(T v, T inc, T max, T min){
			v += inc;
			if(v >= max || v < min) ++count;
			return count < N ? scl::wrap(v, max, min) : scl::clip(v, max, min);
		}
		
		uint32_t count;
	};

	struct Wrap{
		void reset(){}
	
		uint32_t& operator()(uint32_t& pos, uint32_t inc){ return pos+=inc; }
		bool done(uint32_t pos) const { return false; }
		
		template <class T>
		T operator()(T v, T inc, T max, T min){ return scl::wrap(v+inc, max, min); }
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
