#ifndef GAMMA_DELAY_H_INC
#define GAMMA_DELAY_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

//#include "Gamma/arr.h"
#include "Gamma/ipl.h"
#include "Gamma/mem.h"
#include "Gamma/scl.h"
//#include "Gamma/tbl.h"

#include "Gamma/Containers.h"
#include "Gamma/Strategy.h"
#include "Gamma/Sync.h"
#include "Gamma/Types.h"

namespace gam{

/// \defgroup Delays

/// Variable length delay line

/// A delay line is an all-pass filter that shifts samples forward along the 
/// sampling domain. E.g., if the sampling domain is time, then it shifts 
/// samples into the future. This delay line operates using one write head and 
/// one read head by default. More general read and write methods are provided 
/// so that the delay line can be used as a multi-tap delay line.
/// Ideally, a delay line is an all-pass filter which means it does not modify 
/// the magnitude spectrum of the signal, only its phase. However, when reading 
/// samples from a delay line using Lagrange interpolation (linear, cubic, etc.)
/// the output signal will be colored as the interpolation acts like a low-pass
/// filter. To avoid coloration, one must use either no interpolation or 
/// all-pass interpolation. However, a major caveat of both no and all-pass 
/// interpolation is that they introduce undesirable transients when the delay 
/// length is dynamically varied.
///
/// \tparam Tv	value (sample) type
/// \tparam Tp	parameter type
/// \tparam Si	interpolation strategy
/// \tparam Ts	sync type
/// \ingroup Delays    
template <class Tv=gam::real, template<class> class Si=ipl::Linear, class Ts=Synced>
class Delay : public ArrayPow2<Tv>, Ts{
public:

	/// Default constructor. Does not allocate memory.
	Delay();

	/// \param[in]	delay		Delay length
	/// The size of the delay buffer will be the smallest possible power of two.
	Delay(float delay);

	/// \param[in]	maxDelay	Maximum delay length
	/// \param[in]	delay		Delay length
	/// The size of the delay buffer will be the smallest possible power of two.
	Delay(float maxDelay, float delay);


	void delay(float v);						///< Set delay length
	void delaySamples(uint32_t);				///< Set delay length in samples
	void delayUnit(float u);					///< Set delay as (0, 1) of buffer size
	void freq(float v);							///< Set natural frequency (1/delay())
	void ipolType(ipl::Type v){ mIpol.type(v); }///< Set interpolation type
	void maxDelay(float v);						///< Set maximum delay length
	void zero();								///< Sets all delay elements to zero

	Tv operator()(const Tv& v);					///< Returns next filtered value
	Tv operator()() const;						///< Reads delayed element from buffer
	Tv read(float ago);							///< Returns element 'ago' units ago
	void write(const Tv& v);					///< Writes element into buffer. Tap is post-incremented.
	void writePre(const Tv& v);					///< Writes element into buffer. Tap is pre-incremented.
	
	float delay() const;						///< Get current delay length
	uint32_t delaySamples() const;				///< Get current delay length in samples
	float delayUnit() const;					///< Get unit delay (relative to max delay)
	uint32_t delayIndex(uint32_t delay) const;	///< Get index of delayed element	
	float freq() const { return 1.f/delay(); }	///< Get frequency of delay line
	uint32_t indexBack() const;					///< Get index of backmost element
	float maxDelay() const;						///< Get maximum delay length units

	virtual void onResize();
	virtual void onResync(double r);

	void print();

protected:
	Si<Tv> mIpol;

	float mMaxDelay;				// maximum delay length
	float mDelayFactor;				// multiplication factor when setting delay
	float mDelayLength;				// current delay length
	uint32_t mPhase;				// write tap
	uint32_t mPhaseInc;				// phase increment
	uint32_t mDelay;				// read tap as delay from write tap

	void incPhase();				// increment phase
	void refreshDelayFactor();
	uint32_t delayFToI(float v);	// convert f.p. delay to fixed-point
};



/// Variable delay-line with multiple read taps
    
/// \ingroup Delays
///
template <class Tv=gam::real, template <class> class Si=ipl::Linear, class Ts=Synced>
class Multitap : public Delay<Tv,Si,Ts> {
public:

	/// \param[in]	delay		Delay length. The size of the delay line will 
	///							be the smallest possible power of two.
	/// \param[in]	numTaps		Number of reader taps
	Multitap(float delay, uint32_t numTaps)
	:	Delay<Tv,Si,Ts>(delay)
	{	taps(numTaps); }

	/// Get number of read taps
	uint32_t taps() const { return mDelays.size(); }	

	/// Read sample from tap
	Tv read(uint32_t tap) const {
		return mIpol(*this, this->mPhase - mDelays[tap]);
	}

	/// Set delay length
	void delay(float length, uint32_t tap){
		mDelays[tap] = this->delayFToI(length);
	}

	/// Set number of read taps
	void taps(uint32_t numTaps){ mDelays.resize(numTaps); }

protected:
	std::vector<uint32_t> mDelays;
};



/// Fixed-size delay that uses memory-shifting.
    
/// Where N is the number of elements in the delay, insertion is O(N) which is 
/// slower than that of the average ring buffer at O(1). Access, however, will 
/// be optimal (a direct array access) versus that of the ring buffer which 
/// requires an additional conversion of a relative index into an absolute index.
///
/// \tparam N	size of delay
/// \tparam T	value (sample) type
/// \ingroup Delays    
template <uint32_t N, class T>
class DelayShift{
public:
	#define IT for(uint32_t i=0; i<N; ++i)

	/// \param[in] v	Initial value of elements
	DelayShift(const T& v=T()){ IT mElems[i]=v; }

	/// Set nth delayed element
	T& operator[](uint32_t i){ return mElems[i]; }
	
	/// Get nth delayed element
	const T& operator[](uint32_t i) const { return mElems[i]; }

	/// Get elements
	T * elems(){ return mElems; }
	const T * elems() const { return mElems; }


	/// Input element and return Nth delayed element
	T operator()(const T& v) const {
		const T r = mElems[N-1];
		for(uint32_t i=N-1; i>0; --i) mElems[i] = mElems[i-1];
		mElems[0]=v;
		return r;
	}

	/// Get size of delay
	static uint32_t size(){ return N; }

	#undef IT

protected:
	mutable T mElems[N];
};


/// One sample delay. Returns last input sample.

/// \ingroup Delays
///  
template<class T=gam::real> 
class Delay1 : public DelayShift<1,T>{
public:

	/// \param[in] v	Initial value of elements
	Delay1(const T& v=T()): DelayShift<1,T>(v){}
};


/// Two sample delay. Returns second to last input sample.
    
/// \ingroup Delays
///
template<class T=gam::real> 
class Delay2 : public DelayShift<2,T>{
public:

	/// \param[in] v	Initial value of elements
	Delay2(const T& v=T()): DelayShift<2,T>(v){}
	
	/// \param[in] v2	Initial value of 2nd delayed element
	/// \param[in] v1	Initial value of 1st delayed element
	Delay2(const T& v2, const T& v1){ (*this)[1]=v2; (*this)[0]=v1; }
};



/// Variable length delay-line with feedback and/or feedforward.

/// The general comb filter transfer function provides N evenly spaced poles
/// and/or zeroes around the unit circle. Feedback and feedforward produce
/// evenly spaced resonances and notches, respectively, in the frequency
/// response. Positive feeds result in even harmonics and negative feeds give
/// odd harmonics. If the feedback and feedforward amounts are inverses of each
/// other, an Nth order all-pass filter results. Comb filters are stable as
/// long as |feedback| < 1.
///
/// \tparam Tv	value type
/// \tparam Si	interpolation strategy
/// \tparam Tp	parameter type
/// \tparam Ts	sync type
/// \ingroup Delays   
/// \ingroup Filters
// H(z) = (ffd + z^-m) / (1 - fbk z^-m)
// y[n] = ffd x[n] + x[n-m] + fbk y[n-m]
template<
	class Tv=gam::real,
	template <class> class Si=ipl::Linear,
	class Tp=gam::real,
	class Ts=Synced
>
class Comb : public Delay<Tv,Si,Ts> {

private:
	typedef Delay<Tv,Si,Ts> Base;

public:
	using Base::operator();

	/// Default constructor. Does not allocate memory.
	Comb();

	/// \param[in]	delay		Delay length. The size of the delay line will 
	///							be the smallest possible power of two.
	/// \param[in]	ffd			Feedforward amount, in [-1, 1]
	/// \param[in]	fbk			Feedback amount, in (-1, 1)
	Comb(float delay, const Tp& ffd = Tp(0), const Tp& fbk = Tp(0));
	
	/// \param[in]	maxDelay	Maximum delay length. The size of the delay line 
	///							will be the smallest possible power of two.
	/// \param[in]	delay		Delay length
	/// \param[in]	ffd			Feedforward amount, in [-1, 1]
	/// \param[in]	fbk			Feedback amount, in (-1, 1)
	Comb(float maxDelay, float delay, const Tp& ffd, const Tp& fbk);

	
	/// Set number of units for response to decay to end value
	
	/// The sign of the decay length determines the sign of the feedback coefficient.
	/// The default end value of 0.001 (-60 dB) is the reverberation time of 
	/// the filter.  Setting the decay amount effects only the feedback value.
	/// The decay must be updated whenever the delay length of the filter changes.
	void decay(float units, float end = 0.001f);

	/// Sets feedback to argument and feedforward to argument negated
	void allPass(const Tp& v);

	void fbk(const Tp& v);					///< Set feedback amount, in (-1, 1)
	void ffd(const Tp& v);					///< Set feedforward amount [-1, 1]
	void feeds(const Tp& fwd, const Tp& bwd){ ffd(fwd); fbk(bwd); }

	void set(float delay, const Tp& ffd, const Tp& fbk); ///< Set several parameters

	Tv operator()(const Tv& i0);				///< Returns next filtered value
	Tv operator()(const Tv& i0, const Tv& oN);	///< Circulate filter with ffd & fbk
	Tv circulateFbk(const Tv& i0, const Tv& oN);///< Circulate filter with fbk only	

	/// Filters sample (feedback only).
	Tv nextFbk(const Tv& i0);
	
	float norm() const;				///< Get unity gain scale factor
	float normFbk() const;			///< Get unity gain scale factor due to feedback
	float normFfd() const;			///< Get unity gain scale factor due to feedforward
	Tp ffd() const;					///< Get feedforward amount
	Tp fbk() const;					///< Get feedback amount

protected:
	Tp mFFD, mFBK;
};




// Implementation_______________________________________________________________

#define TM1 template <class Tv, template <class> class Ti, class Ts>
#define TM2 Tv,Ti,Ts

#define DELAY_INIT mMaxDelay(0), mDelayFactor(0), mDelayLength(0), mPhase(0), mPhaseInc(0), mDelay(0)

TM1 Delay<TM2>::Delay()
:	ArrayPow2<Tv>(), DELAY_INIT
{	Ts::initSynced(); }

TM1 Delay<TM2>::Delay(float maxDly, float dly)
:	ArrayPow2<Tv>(), DELAY_INIT
{
	Ts::initSynced();
	maxDelay(maxDly);
	zero();
	delay(dly);
}

TM1 Delay<TM2>::Delay(float dly)
:	ArrayPow2<Tv>(), DELAY_INIT
{	//printf("Delay::Delay(float)\n");
	Ts::initSynced();
	maxDelay(dly);
	zero();
	delay(dly);
}

#undef DELAY_INIT

TM1 void Delay<TM2>::maxDelay(float length){ //printf("Delay::maxDelay(%f)\n", length);
	mMaxDelay = length;
	if(Ts::sync() && Ts::sync()->hasBeenSet()){
		//printf("Delay::maxDelay(): resize to %d\n", (uint32_t)(mMaxDelay * spu()));
		
		// This will trigger onResize() -> onResync(double r) calls ONLY if
		// the size changes to prevent infinite recursion.
		this->resize(uint32_t(mMaxDelay * Ts::spu()));
	}
}

TM1 void Delay<TM2>::zero(){ this->assign(Tv(0)); }

TM1 inline Tv Delay<TM2>::operator()() const{ return mIpol(*this, mPhase - mDelay); }

TM1 inline Tv Delay<TM2>::operator()(const Tv& i0){
//	writePre(i0);
//	return (*this)();
	Tv o0 = (*this)();	// read delayed element
	write(i0);			// write input element
	return o0;
}

TM1 inline uint32_t Delay<TM2>::delayFToI(float v){
	return castIntRound((v * mDelayFactor) * 4294967296.);
	//return scl::unitToUInt(v * mDelayFactor);
}

TM1 inline void Delay<TM2>::incPhase(){ mPhase += mPhaseInc; }

TM1 void Delay<TM2>::onResize(){ //printf("Delay::onResize %d elements\n", this->size());
	mPhaseInc = this->oneIndex();
	//for(uint32_t i=0; i<this->size(); ++i) (*this)[i] = Tv(0);
	if(this->isSoleOwner()) zero();
	onResync(1);
}

TM1 void Delay<TM2>::onResync(double r){ //printf("Delay::onSyncChange\n");
	if(this->usingExternalSource()){
		mMaxDelay = this->size() * Ts::ups();
	}
	else{
		maxDelay(mMaxDelay);
	}
	refreshDelayFactor();
	delay(mDelayLength);
}

TM1 inline Tv Delay<TM2>::read(float ago){ return mIpol(*this, mPhase - delayFToI(ago)); }

TM1 void Delay<TM2>::refreshDelayFactor(){ mDelayFactor = 1. / maxDelay(); }

TM1 inline void Delay<TM2>::write(const Tv& v){
	mem::put(this->elems(), this->fracBits(), mPhase, v);
	incPhase();
}

TM1 inline void Delay<TM2>::writePre(const Tv& v){
	incPhase();
	mem::put(this->elems(), this->fracBits(), mPhase, v);
}

TM1 inline void Delay<TM2>::delay(float v){
	mDelayLength = v;
	mDelay = delayFToI(v);
}

TM1 inline float Delay<TM2>::delay() const { return mDelayLength; }

TM1 inline void Delay<TM2>::delaySamples(uint32_t v){ mDelay = v << this->fracBits(); }

TM1 inline uint32_t Delay<TM2>::delaySamples() const { return mDelay >> this->fracBits(); }

TM1 inline void Delay<TM2>::delayUnit(float v){ mDelay = unitToUInt(v); }

TM1 inline float Delay<TM2>::delayUnit() const {
	return uintToUnit<float>(mDelay);
}

TM1 inline uint32_t Delay<TM2>::delayIndex(uint32_t delay) const {
	return this->index(mPhase - (delay << this->fracBits()));
}

TM1 inline void Delay<TM2>::freq(float v){ delay(1.f/v); }

TM1 inline uint32_t Delay<TM2>::indexBack() const {
	return this->index(mPhase + this->oneIndex());
}

TM1 inline float Delay<TM2>::maxDelay() const { return this->size() * Ts::ups(); }

TM1 void Delay<TM2>::print(){
	printf("SPU:       %f\n", Ts::spu());
	printf("Buffer:    %p\n", this->elems());
	printf("BufBits:   %d\n", this->log2Size());
	printf("FracBits:  %d\n", this->fracBits());
	printf("Phase:     %d\n", mPhase);
	printf("PhaseInc:  %d\n", mPhaseInc);
	printf("Delay:     %d, %f\n", mDelay, mDelayLength);
	printf("Max Delay: %f\n",  mMaxDelay);
	printf("DlyFactor: %f\n",  mDelayFactor);
}

#undef TM1
#undef TM2




#define TM1 template<class Tv, template<class> class Si, class Tp, class Ts>
#define TM2 Tv,Si,Tp,Ts
TM1 Comb<TM2>::Comb(): Delay<Tv,Si,Ts>(), mFFD(0), mFBK(0){}

TM1 Comb<TM2>::Comb(float delay, const Tp& ffd, const Tp& fbk)
:	Delay<Tv,Si,Ts>(delay), mFFD(ffd), mFBK(fbk)
{}

TM1 Comb<TM2>::Comb(float delayMax, float delay, const Tp& ffd, const Tp& fbk)
:	Delay<Tv,Si,Ts>(delayMax, delay), mFFD(ffd), mFBK(fbk)
{}

TM1 inline Tv Comb<TM2>::operator()(const Tv& i0){
	return (*this)(i0, (*this)()); }

TM1 inline Tv Comb<TM2>::operator()(const Tv& i0, const Tv& oN){
	Tv t = i0 + oN * mFBK;
	this->write(t);
	return oN + t * mFFD;
}

TM1 inline Tv Comb<TM2>::circulateFbk(const Tv& i0, const Tv& oN){
	this->write(i0 + oN * mFBK);
	return oN;
}

TM1 inline Tv Comb<TM2>::nextFbk(const Tv& i0){
	return circulateFbk(i0, (*this)()); }

TM1 inline void Comb<TM2>::decay(float units, float end){
	mFBK = ::pow(end, this->delay() / scl::abs(units));
	if(units < 0.f) mFBK = -mFBK;
}

TM1 inline void Comb<TM2>::allPass(const Tp& v){ fbk(v); ffd(-v); }
TM1 inline void Comb<TM2>::fbk(const Tp& v){ mFBK=v; }
TM1 inline void Comb<TM2>::ffd(const Tp& v){ mFFD=v; }
TM1 inline void Comb<TM2>::set(float d, const Tp& ff, const Tp& fb){ this->delay(d); ffd(ff); fbk(fb); }

TM1 inline Tp Comb<TM2>::fbk() const { return mFBK; }
TM1 inline Tp Comb<TM2>::ffd() const { return mFFD; }
TM1 inline float Comb<TM2>::norm() const { return (1.f - scl::abs(fbk()))/(1.f + scl::abs(ffd())); }
TM1 inline float Comb<TM2>::normFbk() const { return 1.f - scl::abs(fbk()); }
TM1 inline float Comb<TM2>::normFfd() const { return 1.f/(1.f + scl::abs(ffd())); }

#undef TM1
#undef TM2

} // gam::
#endif
