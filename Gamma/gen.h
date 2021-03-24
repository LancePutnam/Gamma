#ifndef GAMMA_GEN_H_INC
#define GAMMA_GEN_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
*/

#include "Gamma/scl.h"			// abs, wrap
#include "Gamma/Constants.h"	// M_2PI, M_1_2PI
#include "Gamma/Types.h"		// Complex

namespace gam{

/// Generator function objects

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

/// Generates the default value of its associated type
template <class T=gam::real>
struct Default{
	typedef T value_type;
	T operator()() const { return T(); }					///< Generate next value
};


/// Single value generator
template <class T=gam::real>
struct Val : public Default<T>{
	Val(): val(T(0)){}										///< Constructor
	Val(const T& v): val(v){}								///< Constructor
	Val& operator = (const T& v){ val=v; return *this; }	///< Set value
	T operator()() const { return val; }					///< Generate next value
	T& operator[](unsigned i)      { return val; }			///< Array set; sets current value 
	T  operator[](unsigned i) const{ return (*this)(); }	///< Array get; generates next element

	template<class U> bool operator> (const U& v) const { return val> v; }
	template<class U> bool operator>=(const U& v) const { return val>=v; }
	template<class U> bool operator< (const U& v) const { return val< v; }
	template<class U> bool operator<=(const U& v) const { return val<=v; }

	mutable T val;											///< Value
	// Since this is a generator, we will allow its value to be modified if 
	// it's a const.
};


// This is needed since templates are not always smart about inheriting super members.
#define INHERIT\
	using Val<T>::val; using Val<T>::operator=;\
	T   operator[](unsigned i) const { return (*this)(); }

template<class T=gam::real>
struct Impulse : public Val<T>{ INHERIT;
	Impulse(const T& val=T(1)): Val<T>(val){}			///< Constructor
	T operator()() const {T t=val; val=0; return t;}	///< Generate next value
};


/// Generates a Nyquist signal, i.e., -1, 1, -1, 1, …
template<class T=gam::real>
struct Nyquist : public Val<T>{ INHERIT;
	Nyquist(const T& val=T(1)): Val<T>(-val){}			///< Constructor
	T operator()() const { return val = -val; }			///< Generate next value
};


///Reciprocal sequence generator

///Given a type that can be initialized by passing
///the integer 1 to the constructor (let the value be “x”), 
///it generates the sequence x/1, x/2, x/3, x/4, etc.
///http://www.britannica.com/EBchecked/topic/1500010/harmonic-sequence
template <class T=gam::real>
struct Recip : public Val<T>{ INHERIT;
	Recip(const T& val=T(1)): Val<T>(val){}				///< Constructor
	T operator()() const { return T(1)/val++; }			///< Generate next value
};


/// Cosine generator based on recursive formula x0 = c x1 - x2

/// This has an optimal initial state configuration (set) requiring only one
/// trig call. The general sinusoid requires three trig calls.
template <class T=double>
struct RCos : public Val<T>{ INHERIT;

	/// Default constructor
	RCos(){ zero(); }

	/// Constructor
	RCos(const T& frq, const T& amp=T(1)){
		set(frq,amp);
	}

	/// Generate next value
	T operator()() const {
		T v0 = val*c1 - val2;
		val2 = val; val = v0;
		return val2;
	}

	/// Get frequency
	T freq() const { return acos(c1*0.5) * M_1_2PI; }
	
	/// Set parameters from unit freq and amplitude
	RCos& set(const T& frq, const T& amp=T(1)){
		c1  = T(cos(frq*M_2PI));
		val2= c1*amp;
		val = amp;
		c1 *= T(2);
		return *this;
	}

	RCos& constant(const T& amp = T(1)){
		c1 = T(2);
		val2 = amp;
		val = amp;
		return *this;
	}

	RCos& zero(){ return constant(T(0)); }

	mutable T val2;		///< 2-previous value
	T c1;				///< 1-previous coefficient
};


// TODO: zero-overhead version of this with only delays and coef

/// Sinusoidal generator based on recursive formula x0 = c x1 - x2
template <class T=double>
struct RSin : public Val<T>{ INHERIT;

	/// Constructor
	RSin(const T& frq=T(0), const T& amp=T(1), const T& phs=T(0))
	:	val2(0), mul(0){ set(frq,amp,phs); }

	/// Generate next value
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
	RSin& amp(const T& v){ return set(freq(), v, phase()); }
	
	/// Set unit frequency
	RSin& freq(const T& v){	return set(v, amp(), phase()); }

	/// Set unit phase
	RSin& phase(const T& v){ return set(freq(), amp(), v); }

	/// Reset state from stored parameters
	RSin& reset(){ set(freq(), amp(), phase()); return *this; }

	/// Set parameters from unit freq, amplitude and phase
	RSin& set(const T& frq, const T& amp=T(1), const T& phs=T(0)){
		//printf("%g %g %g\n", frq, phs, amp);
		mFreq = frq;
		mAmp = amp;
		mPhase = phs;

		T f=frq*M_2PI, p=phs*M_2PI;
		mul  = 2 * cos(f);
		val2 = sin(p - f * T(2))*amp;
		val  = sin(p - f       )*amp;
		//printf("%g %g %g\n", freq(), phase(), this->amp());
		return *this;
	}

/*	Note: these methods compute parameters directly from coefs, but are buggy...
	
	/// Get amplitude and unit phase
	void ampPhase(T& a, T& p) const {
		p = phase();

		const T eps = 1e-8;
		if(p > eps && p < (1-eps) && scl::abs(p-0.5) > eps)
			a = val /sin(p * M_2PI);
		else
			a = val2/sin((p - freq())*M_2PI);
		return;
	}

	/// Get unit frequency
	T freq() const { return acos(mul*0.5) * M_1_2PI; }

	/// Get amplitude
	T amp() const {	T a,p; ampPhase(a,p); return a; }

	/// Get unit phase
	T phase() const {
		if(val == val2) return 0;
		T f = freq()*M_2PI;
		T y = val * sin(f);
		T x = val * cos(f) - val2;
		T r = atan2(y, x) * M_1_2PI;
		if(r < 0) r += 1;
		return r;
	}
	
	/// Set unit frequency
	RSin& freq(const T& v){	T a,p; ampPhase(a,p); return set(v,a,p); }
*/
	mutable T val2;
	T mul;			///< Multiplication factor

protected:
	T mFreq, mAmp, mPhase;
};


template <class T=double>
struct RSin2 : public Val<T>{ INHERIT;

	/// Constructor
	RSin2(const T& frq=T(0), const T& amp=T(1), const T& dcy=T(1), const T& phs=T(0))
	{ set(frq,amp,dcy,phs); }

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
	RSin2& amp(const T& v){ return set(freq(), v, decay(), phase()); }

	/// Set decay
	RSin2& decay(const T& v){ return set(freq(), amp(), v, phase()); }
	
	/// Set unit frequency
	RSin2& freq(const T& v){ return set(v, amp(), decay(), phase()); }

	/// Set unit phase
	RSin2& phase(const T& v){ return set(freq(), amp(), decay(), v); }

	/// Reset state from stored parameters
	RSin2& reset(){ set(freq(), amp(), decay(), phase()); return *this; }

	/// Set parameters from freq (rad/unit), amplitude, decay and phase
	RSin2& set(T frq, T amp=T(1), T dcy=T(1), T phs=T(0)){
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


/// Recursive Gaussian window generator
template <class T=gam::real>
struct RGauss{
	RGauss(): val(1), mul1(1), mul2(1)
	{}

	/// (Re)set window function

	/// \param[in] length	length of window in samples
	/// \param[in] offset	start/end offset of window over its length
	///						(must be greater than 0)
	void set(double length, T offset=0.001){
		using namespace std;
		length *= 0.5;
		mul2= pow(offset, T(2.)/(length*length)); // b^2
		mul1= pow(T(1.)/mul2, T(length+0.5)); // rewind from center of window
		val = offset;
	}

	/// Generate next value
	T operator()(){
		T res = val;
		mul1 *= mul2;
		val  *= mul1;
		return res;
	}

	/// Returns whether envelope is done
	bool done(float thresh=0.001) const {
		// slope negative and value < threshold
		return (gam::magSqr(mul1) < 1.) && (gam::magSqr(val) < thresh*thresh);
	}

	double length() const {
		/*m1 = (1/m2)^(L+0.5)
		log(m1) = (L+0.5)log(1/m2)
		L = log(m1)/log(1/m2) - 0.5
		L = log(m1)/(-log(m2)) - 0.5*/
		return (log(gam::norm(mul1))/-log(gam::norm(mul2)) - 0.5)*2;
	}

	T offset() const {
		/*m2 = O^(2/LL)
		m2^(LL/2) = O*/
		double L = length();
		return pow(mul2, T(L*L*0.125));
	}

	T val, mul1, mul2;
};


/// Recursive add generator that generates lines
template <class T=gam::real>
struct RAdd: public Val<T>{ INHERIT;

	/// \param[in] add	addition amount
	/// \param[in] val	current value
	RAdd(const T& add=T(1), const T& val=T(0))
	: Val<T>(val-add), add(add){}
	
	/// Generate next value
	const T& operator()() const { return val += add; }

	/// Go back one step
	const T& recede() const { return val -= add; }

	/// Set to generate line between points (0, begin) and (length, end)
	void line(T begin, T end, T length){
		val = begin;
		add = (end - begin)/length;
		recede();
	}

	/// Set to generate line between points (0, val) and (length, end)
	void line(T end, T length){ line(val, end, length); }

	/// Set to generate constant value
	void constant(T v){ val=v; add=T(0); }

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
	:	Val<T>((val-add)/mul), mul(mul), add(add){}

	T operator()() const { return val=val*mul+add; }	///< Generate next value

	/// Go back one step
	const T& recede() const { return val = (val-add)/mul; }

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



/// Complex resonator

/// A complex resonator consists of two complex numbers- one is an absolute
/// phase and amplitude and the other is a relative (differential) frequency 
/// and decay/grow factor.
template <class T=double>
class CReson : public Complex<T>{
public:
	typedef Complex<T> C;
	using C::operator();
	using C::operator=;

	/// \param[in] frq	unit frequency
	/// \param[in] amp	amplitude
	/// \param[in] phs	unit phase
	/// \param[in] dec	unit decay/grow factor
	CReson(const T& frq=T(0), const T& amp=T(1), const T& phs=T(0), const T& dec=T(1)){
		set(frq, amp, phs);
	}

	/// Advance one iteration and return value
	const C& operator()(){ return (*this)*=mFactor; }

	/// Filter input
	const C& operator()(const C& v){ return (*this) = (*this)*mFactor + v; }
	const C& operator()(const T& v){ return (*this) = (*this)*mFactor + v; }
	
	/// Recede one iteration and return value
	const C& recede(){ return (*this)/=mFactor; }

	/// Set amplitude
	void amp(const T& v){ (*this).fromPolar(v, this->arg()); }
	
	/// Set complex amplitude
	template <class U>
	void amp(const Complex<U>& v){ (*this)(v.r, v.i); }
	
	/// Set amplitude decay/grow factor after N iterations
	void decay(const T& target, const T& N=T(1)){
		// NOTE: this handles negative decays, thought better to leave this to frequency component
		//factor(freq(), (1==N) ? target : (pow(scl::abs(target), 1./N)*scl::sgn(target)));
		factor(freq(), (1==N) ? target : (pow(scl::abs(target), 1./N)));
	}

	/// Set recursive multiplication factor (frequency and decay/growth factor)
	void factor(const T& frq, const T& dec=T(1)){
		mFactor.fromPolar(dec, frq*M_2PI);
	}
	
	void factor(const Complex<T>& v){ mFactor=v; }

	/// Set unit frequency
	void freq(const T& v){ factor(v, decay()); }

	/// Set unit frequency, amplitude, unit phase, and decay/grow factor

	/// The phase state will be rewound 1 iteration so the first function call
	/// will return a complex number at the desired phase.
	void set(const T& frq, const T& amp, const T& phs, const T& dec=T(1)){
		this->fromPolar(amp, phs*M_2PI);
		factor(frq, dec);
		recede();
	}

	void set(const T& frq, const Complex<T>& phs){
		(*this)(phs.r, phs.i); freq(frq);
	}

	/// Get value one iteration ahead of current state
	C ahead() const { return (*this)*mFactor; }
	
	/// Get value one iteration behind current state
	C behind() const { return (*this)/mFactor; }

	/// Get unit decay
	T decay() const { return mFactor.mag(); }
	
	/// Get unit frequency
	T freq() const { return mFactor.phase()*M_1_2PI; }

	const C& factor() const { return mFactor; }
	C& factor(){ return mFactor; }

protected:
	C mFactor;

	// Set 60 dB decay interval
	//void decay(const Tv& v){ width(T(2.198806796637603 /* -ln(0.001)/pi */)/v); }

	// Set unit bandwidth
	//void width(const Tv& v){ mDecay=::exp(-M_PI*v); freq(freq()); }
};

typedef CReson<float>	CResonf;
typedef CReson<double>	CResond;



struct OnOff{
	OnOff(unsigned _max, unsigned _ons) : max(_max), ons(_ons), cnt(0){}

	bool operator()(){
		cnt++;
		if(cnt <= ons) return true;
		if(cnt >= max) cnt = 0;
		return ons >= max;
	}

	void set(unsigned _max, unsigned _ons, unsigned _cnt){
		this->max = _max; this->ons = _ons; this->cnt = _cnt;
	}

	unsigned max, ons, cnt;
};


/// Returns true once, then false until reset
class OneOff{
public:
	OneOff(bool v=true): mVal(v) {}

	/// Get next value
	bool operator()(){ bool r=mVal; mVal=false; return r; }

	/// Reset trigger
	void reset(){ mVal=true; }

private:
	bool mVal;
};


/// Fixed-sized array with a sequence generator
template <unsigned N, class T=gam::real, class G=gen::RAdd1<unsigned> >
class Seq: public Vec<N,T>{
public:

	Seq(const T& val){ this->set(val); }
	Seq(const T * vals){ this->set(vals); }

	/// Generate next array element
	T operator()(){ return (*this)[((unsigned)mTap())%N]; }

	/// Get reference to index generator
	G& tap(){ return mTap; }

private:
	G mTap;
};


/// Triggers after a specified number of iterations and then resets.

/// Outputs true on every nth sample and false on the rest.
/// Argument "num" determines the length of the sequence.
/// Argument "val", with a default value of zero can be set by the user
/// to adjust the location of the triggering sample within the sequence.
struct Trigger{
	Trigger(unsigned num, unsigned val=0) : val(val), num(num){}

	/// Returns (triggers) true upon reset
	bool operator()(){
		if(++val >= num){ val = 0; return true; }
		return false;
	}

	unsigned val;		///< Value
	unsigned num;		///< Maximum value
};


} // gen::
} // gam::

#endif
