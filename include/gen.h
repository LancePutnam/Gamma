#ifndef GAMMA_GEN_H_INC
#define GAMMA_GEN_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include "scl.h"
#include "Constants.h"
#include "Types.h"
//#include "MacroD.h"

namespace gam{


///	Generator function objects

/// A generator is a lightweight object that generates a sequence of elements.
/// Generators have a standard interface specified by the Val class. The array
/// access operator, [], is overloaded so generators can be treated like
/// arrays. A const qualified generator only means that its generating function 
/// parameters are held constant; its current value can change.
/// A convention used is that generators' current value is the one most recently
/// generated. This means one must be careful when initializing a generator so 
/// that it generates the first element properly. Usually this means setting its
/// value to what would have been the previously generated value. For instance, 
/// to get the sequence 0,1,2,... from RAdd1, its value must be initialized to 
/// -1. Constructors 'rewind' the generator, so that the initial value argument
/// is returned on the next generate call.

namespace gen{

/// Single value generator
template <class T=gam::real>
struct Val{
	Val(): val((T)0){}								///< Constructor
	Val(T v): val(v){}								///< Constructor
	Val& operator = (T v){ val=v; return *this; }	///< Set value
	T operator()() const { return val; }			///< Generate next value
	T& operator[](uint i)      { return val; }		///< Array set; sets current value 
	T  operator[](uint i) const{ return (*this)(); }///< Array get; generates next element
	mutable T val;									///< Value
	// Since this is a generator, we will allow its value to be modified if 
	// it's a const.
};


// This is needed since templates are not always smart about inheriting super members.
#define INHERIT\
	using Val<T>::val; using Val<T>::operator=;\
	T   operator[](uint i) const { return (*this)(); }


template<class T = gam::real>
struct Impulse : public Val<T>{ INHERIT;
	Impulse(T val=1): Val<T>(val){}						///< Constructor
	T operator()() const {T t=val; val=0; return t;}	///< Generate next value
};

/// Nyquist sequence generator
template<class T = gam::real>
struct Nyquist : public Val<T>{ INHERIT;
	Nyquist(T val=1): Val<T>(-val){}					///< Constructor
	T operator()() const { return val = -val; }			///< Generate next value
};

/// Reciprocal sequence generator
template <class T=gam::real>
struct Recip : public Val<T>{ INHERIT;
	Recip(T val=1): Val<T>(val){}						///< Constructor
	T operator()() const { return (T)1/val++; }			///< Generate next value
};


/// Cosine generator based on recursive formula x0 = c x1 - x2
template <class T=gam::real>
struct RCos : public Val<T>{ INHERIT;

	/// Constructor
	RCos(T frq=0, T amp=1): val2(0), c1(0){ set(frq,amp); }

	/// Generate next value.
	T operator()() const {
		T v0 = val*c1 - val2;
		val2 = val; val = v0;
		return val2;
	}
	
	T freq() const { return acos(c1*0.5) * M_1_2PI; }
	
	/// Set parameters from unit freq, phase, and amplitude.
	RCos& set(T frq, T amp=(T)1){
		c1  = T(cos(frq*M_2PI));
		val2= c1*amp;
		val = amp;
		c1 *= T(2);
		return *this;
	}

	mutable T val2;		///< 2-previous value
	T c1;				///< 1-previous coefficient
};


/// Sinusoidal generator based on recursive formula x0 = c x1 - x2
template <class T=gam::real>
struct RSin : public Val<T>{ INHERIT;

	/// Constructor
	RSin(T frq=0, T phs=0, T amp=1): val2(0), mul(0){ set(frq,phs,amp); }

	/// Generate next value.
	T operator()() const {
		T v0 = mul * val - val2;
		val2 = val;
		return val = v0;
	}
	

	void next3(T& o2, T& o1, T& o0) const {
		T v0 = o0 = mul * val  - val2;
		o2 = val2 = mul * v0   - val;
		o1 = val  = mul * val2 - v0;
	}

	
	T freq() const { return acos(mul*0.5) * M_1_2PI; }
	
	/// Set parameters from unit freq, phase, and amplitude.
	RSin& set(T frq, T phs, T amp=(T)1){
		frq *= M_2PI; phs *= M_2PI;
		mul  = (T)2 * (T)cos(frq);
		val2 = (T)sin(phs - frq * (T)2) * amp;
		val  = (T)sin(phs - frq       ) * amp;
		return *this;
	}

	mutable T val2;
	T mul;			///< Multiplication factor. [-2, 2] range gives stable sinusoids.
};


template <class T=gam::real>
struct RSin2 : public Val<T>{ INHERIT;

	/// Constructor
	RSin2(T frq=0, T phs=0, T dcy=1, T amp=1){ set(frq,phs,dcy,amp); }

	/// Generate next value
	T operator()() const {
		T v0 = mul1 * val + mul2 * val2;
		val2 = val;
		return val = v0;
	}
	
	/// Filter next value
	T operator()(const T& v) const {
		T v0 = mul1 * val + mul2 * val2 + v;
		val2 = val;
		return val = v0;
	}
	
	T freq() const { return acos(mul1*0.5/decay()) * M_1_2PI; }
	T decay() const{ return sqrt(-mul2); }
	
	/// Set parameters from freq (rad/unit), phase (rad), decay, and amplitude.
	RSin2& set(T frq, T phs, T dcy, T amp=(T)1){
		frq *= M_2PI; phs *= M_2PI;
		dcy = scl::atLeast(dcy, (T)0.00000001);
		
		mul1 = (T)2 * dcy * cos(frq);
		mul2 = -dcy*dcy;
		T rdcy = 1./dcy;
		val2 = (T)sin(phs - frq * (T)2) * amp * rdcy * rdcy;
		val  = (T)sin(phs - frq       ) * amp * rdcy;
		return *this;
	}

	mutable T val2;
	T mul2, mul1;			///< Multiplication factor. [-2, 2] range gives stable sinusoids.
};


/// Recursive add generator
template <class T=gam::real>
struct RAdd: public Val<T>{ INHERIT;
	RAdd(T add=1, T val=0): Val<T>(val-add), add(add){}	///< Constructor
	T operator()() const { return val += add; }			///< Generate next value
	T add;												///< Addition amount
};

/// Recursive add 1 generator (iota function)
template <class T=gam::real>
struct RAdd1: public Val<T>{ INHERIT;
	RAdd1(T val=0): Val<T>(val-1){}						///< Constructor
	T operator()() const { return ++val; }				///< Generate next value
};

/// Recursive add integer generator
template <int N=1, class T=gam::real>
struct RAddN: public Val<T>{ INHERIT;
	RAddN(T val=0): Val<T>(val-N){}						///< Constructor
	T operator()() const { return val+=N; }				///< Generate next value
};

/// Recursive multiply generator
template <class T=gam::real>
struct RMul: public Val<T>{ INHERIT;
	RMul(T mul=1, T val=1): Val<T>(val/mul), mul(mul){}	///< Constructor
	T operator()() const { return val *= mul; }			///< Generate next value
	T mul;												///< Multiplication amount
};

/// Recursive multiply-add generator
template <class T=gam::real>
struct RMulAdd: public Val<T>{ INHERIT;
	/// Constructor
	RMulAdd(T mul=1, T add=0, T val=0): Val<T>((val-add)/mul), mul(mul), add(add){} 
	T operator()() const { return val=val*mul+add; }	///< Generate next value
	T mul;												///< Multiplication amount
	T add;												///< Addition amount
};


/// Sawtooth wave in interval [0, max)
template <class T=gam::real>
struct Saw: public RAdd<T>{ INHERIT;
	Saw(T add, T val=0, T max=1): RAdd<T>(add, val), max(max){}
	T operator()() const {
		RAdd<T>::operator()();
		if(val >= max) val -= max;
		return val;
	} 
	T max;
};


/// Sinusoid sequence generator
template <class T=gam::real>
struct Sin: public RAdd<T>{ INHERIT;
	Sin(T add, T val=0, T amp=1): RAdd<T>(add, val), amp(amp){}		///< Constructor
	T operator()() const { return (T)sin(RAdd<T>::operator()()) * amp; }	///< Generate next value
	T amp;
};


// Temp object functions
// Note: we cannot use default arguments with two different template types
//
#define OF1(Ob, fu, a)\
template<class T> inline Ob<T> fu(T v1=a){ return Ob<T>(v1); }

//#define OF2(Ob, fu, a, b)
//template<class T, class U> inline Ob<T> fu(T v1=a,U v2=b){ return Ob<T>(v1,v2); }

#define OF2(Ob, fu, a, b)\
template<class T> inline Ob<T> fu(T v1=a,T v2=b){ return Ob<T>(v1,v2); }

#define OF3(Ob, fu, a, b, c)\
template<class T> inline Ob<T> fu(T v1=a,T v2=b,T v3=c){ return Ob<T>(v1,v2,v3); }

OF1(Val,		val,		0)
OF1(Nyquist,	nyquist,	0)
OF1(RAdd1,		rAdd1,		0)
OF1(Recip,		recip,		1)

OF2(RAdd,		rAdd,		1,0)
OF2(RCos,		rCos,		0,1)
OF2(RMul,		rMul,		1,1)

OF3(RMulAdd,	rMulAdd,	1,0,0)
OF3(RSin,		rSin,		0,0,1)

#undef OF1
#undef OF2

/// Triggers after a specified number of iterations and then resets
struct Trigger{
	Trigger(uint32_t num, uint32_t val=0) : val(val), num(num){}
	
	/// Returns (triggers) true upon reset
	bool operator()(){
		if(++val >= num){ val = 0; return true; }
		return false;
	}
	uint32_t val;		///< Value
	uint32_t num;		///< Maximum value
};


struct OnOff{
	OnOff(uint32_t max, uint32_t ons) : max(max), ons(ons), cnt(0){}
	
	bool operator()(){
		cnt++;
		if(cnt <= ons) return true;
		if(cnt >= max) cnt = 0;
		return ons >= max;
	}
	
	void set(uint32_t max, uint32_t ons, uint32_t cnt){
		this->max = max; this->ons = ons; this->cnt = cnt;
	}
	
	uint32_t max, ons, cnt;
};

#undef INHERIT

//template <class T=uint32_t>
//struct Wrapper{
//	Wrapper(T max, T inc=1, T val=0) : val(val), max(max){}
//
//	bool operator()(){
//		if((val += inc) >= max){ scl::wrap(val, max); return true; }
//		return false;
//	}
//	T val, inc, max;
//};

//template <class Ta=gam::real, class Tv=gam::real>
//struct TAdd{
//	TAdd(Ta add, Tv val=0): val(val), add(add){}		
//	Tv operator()(){ Tv r = val; val+=add(); return r; }
//	Tv val;
//	Ta add;
//};
//
//template <class Ta=gam::real, class Tv=gam::real>
//struct TSin : public TAdd<Ta, Tv>{
//	TSin(Ta add, Tv val=0): TAdd<Ta, Tv>(add, val){}		
//	Tv operator()(){ return (Tv)sin( TAdd<Ta, Tv>::operator()() ); }
//};
//
//template <class Tg1, class Tg2, typename Tr=gam::real>
//struct TMulOp{
//	
//	TMulOp(Tg1 gen1, Tg2 gen2): gen1(gen1), gen2(gen2){}
//
//	Tr operator()(){ return gen1() * gen2(); }
//	
//	Tg1 gen1;
//	Tg2 gen2;
//};
//
//template <class Tg1, class Tg2>
//float mul(Tg1 & gen1, Tg2 & gen2){ return (float)(gen1() * gen2()); }


//template <class T=gam::real>
//struct Test : public Val<T>{ INHERIT;
//
//	/// Constructor
//	Test(T val3=0, T val2=0, T val1=0, T mul3=1, T mul2=-1, T mul1=0)
//	:	Val<T>(val1), val3(val3), val2(val2), mul3(mul3), mul2(mul2), mul1(mul1){}
//
//	///< Generate next value.
//	T operator()() const {
//		T v0 = mul1*val + mul2*val2 + mul3*val3;
//		val3 = val2;
//		val2 = val;
//		return val = v0;
//	}
//
//	mutable T val3, val2;
//	T mul3, mul2, mul1;			///< Multiplication factor. [-2, 2] range gives stable sinusoids.
//};


} // gen::
} // gam::

//#include "MacroU.h"

#endif

