#ifndef GAMMA_GEN_H_INC
#define GAMMA_GEN_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information

	File Description:

	Generator function objects

	A generator is a lightweight object that generates a sequence of elements.
	Generators have a standard interface specified by the Val class. The array
	access operator, [], is overloaded so generators can be treated like
	arrays. A const qualified generator only means that its generating function 
	parameters are held constant; its current value can change.
	A convention used is that generators' current value is the one most recently
	generated. This means one must be careful when initializing a generator so 
	that it generates the first element properly. Usually this means setting its
	value to what would have been the previously generated value. For instance, 
	to get the sequence 0,1,2,... from RAdd1, its value must be initialized to 
	-1. Constructors 'rewind' the generator, so that the initial value argument
	is returned on the next generate call.
*/

#include "Gamma/scl.h"
#include "Gamma/Constants.h"
#include "Gamma/Types.h"

namespace gam{

/// Generator function objects
namespace gen{

/// Single value generator
template <class T=gam::real>
struct Val{
	Val(): val(T(0)){}										///< Constructor
	Val(const T& v): val(v){}								///< Constructor
	Val& operator = (const T& v){ val=v; return *this; }	///< Set value
	T operator()() const { return val; }					///< Generate next value
	T& operator[](uint32_t i)      { return val; }			///< Array set; sets current value 
	T  operator[](uint32_t i) const{ return (*this)(); }	///< Array get; generates next element
	mutable T val;											///< Value
	// Since this is a generator, we will allow its value to be modified if 
	// it's a const.
};


// This is needed since templates are not always smart about inheriting super members.
#define INHERIT\
	using Val<T>::val; using Val<T>::operator=;\
	T   operator[](uint32_t i) const { return (*this)(); }

template<class T=gam::real>
struct Impulse : public Val<T>{ INHERIT;
	Impulse(const T& val=T(1)): Val<T>(val){}			///< Constructor
	T operator()() const {T t=val; val=0; return t;}	///< Generate next value
};

/// Nyquist sequence generator
template<class T=gam::real>
struct Nyquist : public Val<T>{ INHERIT;
	Nyquist(const T& val=T(1)): Val<T>(-val){}			///< Constructor
	T operator()() const { return val = -val; }			///< Generate next value
};

/// Reciprocal sequence generator
template <class T=gam::real>
struct Recip : public Val<T>{ INHERIT;
	Recip(const T& val=T(1)): Val<T>(val){}				///< Constructor
	T operator()() const { return T(1)/val++; }			///< Generate next value
};


/// Cosine generator based on recursive formula x0 = c x1 - x2
template <class T=double>
struct RCos : public Val<T>{ INHERIT;

	/// Constructor
	RCos(const T& frq=T(0), const T& amp=T(1))
	:	val2(0), c1(0){ set(frq,amp); }

	/// Generate next value.
	T operator()() const {
		T v0 = val*c1 - val2;
		val2 = val; val = v0;
		return val2;
	}
	
	T freq() const { return acos(c1*0.5) * M_1_2PI; }
	
	/// Set parameters from unit freq, phase, and amplitude.
	RCos& set(const T& frq, const T& amp=T(1)){
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
template <class T=double>
struct RSin : public Val<T>{ INHERIT;

	/// Constructor
	RSin(const T& frq=T(0), const T& phs=T(0), const T& amp=T(1))
	:	val2(0), mul(0){ set(frq,phs,amp); }

	/// Generate next value.
	T operator()() const {
		T v0 = mul * val - val2;
		val2 = val;
		return val = v0;
	}

	/// Get amplitude
	T amp() const {	return mAmp; }

	/// Get unit frequency
	T freq() const { return mFreq; }

	/// Get unit phase
	T phase() const { return mPhase; }

	/// Set amplitude
	RSin& amp(const T& v){ return set(freq(), phase(), v); }
	
	/// Set unit frequency
	RSin& freq(const T& v){	return set(v, phase(), amp()); }

	/// Set unit phase
	RSin& phase(const T& v){ return set(freq(), v, amp()); }

	/// Reset state from stored parameters
	RSin& reset(){ set(freq(), phase(), amp()); return *this; }

	/// Set parameters from unit freq, phase, and amplitude.
	RSin& set(const T& frq, const T& phs, const T& amp=T(1)){
//		printf("%g %g %g\n", frq, phs, amp);
		mFreq = frq;
		mAmp = amp;
		mPhase = phs;

		T f=frq*M_2PI, p=phs*M_2PI;
		mul  = 2 * cos(f);
		val2 = sin(p - f * T(2))*amp;
		val  = sin(p - f       )*amp;
//		printf("%g %g %g\n", freq(), phase(), this->amp());
		return *this;
	}

// Note: these methods compute parameters directly from coefficients, but is buggy...
//	/// Get amplitude
//	T amp() const {	T a,p; ampPhase(a,p); return a; }
//	
//	/// Get amplitude and unit phase
//	void ampPhase(T& a, T& p) const {
//		p = phase();
//
//		const T eps = 1e-8;
//		if(p > eps && p < (1-eps) && scl::abs(p-0.5) > eps)
//			a = val /sin(p * M_2PI);
//		else
//			a = val2/sin((p - freq())*M_2PI);
//		return;
//	}
//
//	/// Get unit frequency
//	T freq() const { return acos(mul*0.5) * M_1_2PI; }
//
//	/// Get unit phase
//	T phase() const {
//		if(val == val2) return 0;
//		T f = freq()*M_2PI;
//		T y = val * sin(f);
//		T x = val * cos(f) - val2;
//		T r = atan2(y, x) * M_1_2PI;
//		if(r < 0) r += 1;
//		return r;
//	}
//
//	/// Set amplitude
//	RSin& amp(const T& v){ return set(freq(), phase(), v); }
//	
//	/// Set unit frequency
//	RSin& freq(const T& v){	T a,p; ampPhase(a,p); return set(v,a,p); }
//
//	/// Set unit phase
//	RSin& phase(const T& v){ return set(freq(), v, amp()); }

	mutable T val2;
	T mul;			///< Multiplication factor

protected:
	T mFreq, mAmp, mPhase;
};


template <class T=double>
struct RSin2 : public Val<T>{ INHERIT;

	/// Constructor
	RSin2(const T& frq=T(0), const T& phs=T(0), const T& dcy=T(1), const T& amp=T(1))
	{ set(frq,phs,dcy,amp); }

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

	/// Get amplitude
	T amp() const {	return mAmp; }

	/// Get decay factor
	//T decay() const{ return mDecay; }
	T decay() const{ return sqrt(-mul2); }

	/// Get unit frequency
	T freq() const { return mFreq; }
	//T freq() const { return acos(mul1*0.5/decay()) * M_1_2PI; }

	/// Get unit phase
	T phase() const { return mPhase; }


	/// Set amplitude
	RSin2& amp(const T& v){ return set(freq(), phase(), decay(), v); }

	/// Set decay
	RSin2& decay(const T& v){ return set(freq(), phase(), v, amp()); }
	
	/// Set unit frequency
	RSin2& freq(const T& v){ return set(v, phase(), decay(), amp()); }

	/// Set unit phase
	RSin2& phase(const T& v){ return set(freq(), v, decay(), amp()); }

	/// Reset state from stored parameters
	RSin2& reset(){ set(freq(), phase(), decay(), amp()); return *this; }

	/// Set parameters from freq (rad/unit), phase (rad), decay, and amplitude.
	RSin2& set(T frq, T phs, T dcy, T amp=T(1)){
//		printf("%g %g %g %g\n", frq, phs, dcy, amp);
		mFreq = frq;
		mAmp = amp;
		mPhase = phs;
		mDecay = dcy;

		frq *= M_2PI; phs *= M_2PI;
		//dcy = scl::atLeast(dcy, (T)0.00000001);
		
		mul1 = (T)2 * dcy * cos(frq);
		mul2 = -dcy*dcy;
		T rdcy = 1./dcy;
		val2 = (T)sin(phs - frq * (T)2) * amp * rdcy * rdcy;
		val  = (T)sin(phs - frq       ) * amp * rdcy;
//		printf("%g %g %g %g\n", freq(), phase(), decay(), this->amp());
		return *this;
	}

	mutable T val2;
	T mul2, mul1;			///< Multiplication factors

protected:
	T mFreq, mAmp, mPhase, mDecay;
};


/// Recursive add generator
template <class T=gam::real>
struct RAdd: public Val<T>{ INHERIT;
	RAdd(const T& add=T(1), const T& val=T(0))
	: Val<T>(val-add), add(add){}						///< Constructor
	T operator()() const { return val += add; }			///< Generate next value
	T add;												///< Addition amount
};

/// Recursive add 1 generator (iota function)
template <class T=gam::real>
struct RAdd1: public Val<T>{ INHERIT;
	RAdd1(const T& val=T(0)): Val<T>(val-1){}			///< Constructor
	T operator()() const { return ++val; }				///< Generate next value
};

/// Recursive add integer generator
template <int N=1, class T=gam::real>
struct RAddN: public Val<T>{ INHERIT;
	RAddN(const T& val=T(0)): Val<T>(val-N){}			///< Constructor
	T operator()() const { return val+=N; }				///< Generate next value
};

/// Recursive addition wrapped in interval [min, max)
template <class T=gam::real>
struct RAddWrap: public RAdd<T>{ INHERIT;
	RAddWrap(const T& add, const T& val=T(0), const T& max=T(1), const T& min=T(0))
	:	RAdd<T>(add, val), max(max), min(min){}			///< Constructor
	T operator()() const {								///< Generate next value
		RAdd<T>::operator()();
		return val=scl::wrap(val,max,min);
	} 
	T max, min;
};

/// Recursive multiply generator
template <class T=gam::real>
struct RMul: public Val<T>{ INHERIT;
	RMul(const T& mul=T(1), const T& val=T(1))
	:	Val<T>(val/mul), mul(mul){}						///< Constructor
	T operator()() const { return val *= mul; }			///< Generate next value
	T mul;												///< Multiplication amount
};

/// Recursive multiply-add generator
template <class T=gam::real>
struct RMulAdd: public Val<T>{ INHERIT;
	/// Constructor
	RMulAdd(const T& mul=T(1), const T& add=T(0), const T& val=T(0))
	:	Val<T>((val-add)/mul), mul(mul), add(add){}		///< Constructor
	T operator()() const { return val=val*mul+add; }	///< Generate next value
	T mul;												///< Multiplication amount
	T add;												///< Addition amount
};

/// Sinusoid sequence generator
template <class T=gam::real>
struct Sin: public RAdd<T>{ INHERIT;
	Sin(const T& frq, const T& phs=T(0), const T& amp=T(1))
	:	RAdd<T>(frq, phs), amp(amp){}					///< Constructor

	T operator()() const
	{ return (T)sin(RAdd<T>::operator()()) * amp; }		/// Generate next value
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
#undef OF3
#undef INHERIT


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


struct OneOff{
	OneOff(bool v=true): mVal(v) {}
	bool operator()(){ bool r=mVal; mVal=false; return r; }
	void set(){ mVal=true; }
	
private:
	bool mVal;
};


/// Fixed-sized array with a sequence generator
template <uint32_t N, class T=gam::real, class G=gen::RAdd1<uint32_t> >
class Seq: public Multi<N,T>{
public:

	Seq(const T& val){ for(int i=0; i<N; ++i){ this->elems[i]=val; } }
	Seq(const T * vals){ for(int i=0; i<N; ++i){ this->elems[i]=vals[i]; } }

	/// Generate next array element
	T operator()(){ return (*this)[((uint32_t)mTap())%N]; }

	/// Get reference to index generator
	G& tap(){ return mTap; }
	
private:
	G mTap;
};


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


} // gen::
} // gam::

#endif
