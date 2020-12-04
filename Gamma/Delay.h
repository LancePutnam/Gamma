#ifndef GAMMA_DELAY_H_INC
#define GAMMA_DELAY_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include "Gamma/ipl.h"
#include "Gamma/scl.h"
#include "Gamma/tbl.h"
#include "Gamma/Containers.h"
#include "Gamma/Domain.h"
#include "Gamma/Strategy.h"
#include "Gamma/Types.h"

namespace gam{

/// Filters that keep a history of previous input samples

/// \defgroup Delay


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
/// \tparam Tv	Value (sample) type
/// \tparam Si	Interpolation strategy
/// \tparam Td	Domain type
/// \ingroup Delay
template<
	class Tv = gam::real,
	template<class> class Si = ipl::Linear,
	class Td = GAM_DEFAULT_DOMAIN
>
class Delay : public ArrayPow2<Tv>, public Td{
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
	void delaySamples(uint32_t v);				///< Set delay length in samples
	void delaySamplesR(float v);				///< Set delay length in samples (real-valued)
	void delayUnit(float u);					///< Set delay as (0, 1) of buffer size
	void freq(float v);							///< Set natural frequency (1/delay())
	void ipolType(ipl::Type v){mIpol.type(v);}	///< Set interpolation type
	void maxDelay(float v, bool setDelay=true);	///< Set maximum delay length

	Tv operator()(const Tv& v);					///< Returns next filtered value
	Tv operator()() const;						///< Reads delayed element from buffer
	Tv read(float ago) const;					///< Returns element 'ago' units ago
	template <template <class> class InterpolationStrategy> Tv read(float ago) const;
	void write(const Tv& v);					///< Writes new element into buffer

	/// Copy delay elements to another array

	/// \param[out] dst		array to copy element to
	/// \param[ in] len		number of elements to copy
	/// \param[ in] end		copy begins at (len+end) elements ago and
	///							ends end elements ago
	template <class V>
	void read(V * dst, unsigned len, unsigned end=0) const;

	float delay() const;						///< Get current delay length
	uint32_t delaySamples() const;				///< Get current delay length in samples
	float delaySamplesR() const;				///< Get current delay length in samples (real-valued)
	float delayUnit() const;					///< Get unit delay (relative to max delay)
	uint32_t delayIndex(uint32_t delay) const;	///< Get index of delayed element	
	float freq() const { return 1.f/delay(); }	///< Get frequency of delay line
	uint32_t indexBack() const;					///< Get index of backmost element
	float maxDelay() const;						///< Get maximum delay length units

	virtual void onResize();
	void onDomainChange(double r);

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
	uint32_t delayFToI(float v) const; // convert f.p. delay to fixed-point
};



/// Variable delay-line with multiple read taps

/// \ingroup Delay
///
template <
	class Tv = gam::real,
	template <class> class Si = ipl::Linear,
	class Td = GAM_DEFAULT_DOMAIN
>
class Multitap : public Delay<Tv,Si,Td> {
public:

	/// \param[in]	delay		Delay length. The size of the delay line will 
	///							be the smallest possible power of two.
	/// \param[in]	numTaps		Number of reader taps
	Multitap(float delay, unsigned numTaps)
	:	Delay<Tv,Si,Td>(delay)
	{	taps(numTaps); }

	/// Get number of read taps
	unsigned taps() const { return mDelays.size(); }	

	/// Read sample from tap
	Tv read(unsigned tap) const {
		return this->mIpol(*this, this->mPhase - mDelays[tap]);
	}

	/// Set delay length
	void delay(float length, unsigned tap){
		mDelays[tap] = this->delayFToI(length);
	}
	
	/// Set a tap's delay length as a frequency
	void freq(float v, unsigned tap){
		delay(1.f/v, tap);
	}

	/// Set number of read taps
	void taps(unsigned numTaps){ mDelays.resize(numTaps); }

protected:
	std::vector<unsigned> mDelays;
};



/// Fixed-size delay that uses memory-shifting.

/// This delay employs a shift buffer rather than the ring buffer used in a
/// typical delay. The primary advantage of a shift buffer is that elements
/// remain sorted in the buffer chronologically, leading to optimal access.
/// Access is also simple making this delay a good choice for implementing an
/// small-sized FIR filter. The disadvantage of a shift buffer is that for every
/// new element added, all elements must be shifted (moved) over by one array 
/// slot (i.e. insertion is O(N)).
///
/// \tparam N	size of delay
/// \tparam T	value (sample) type
/// \ingroup Delay
template <unsigned N, class T>
class DelayShift{
public:

	/// \param[in] v	Initial value of elements
	DelayShift(const T& v=T()){ for(auto& e:mElems) e=v; }

	/// Set nth delayed element
	T& operator[](unsigned i){ return mElems[i]; }
	
	/// Get nth delayed element
	const T& operator[](unsigned i) const { return mElems[i]; }

	/// Get elements
	T * elems(){ return mElems; }
	const T * elems() const { return mElems; }

	/// Write new element
	void write(const T& v){
		for(unsigned i=N-1; i>0; --i) mElems[i] = mElems[i-1];
		mElems[0]=v;
	}

	/// Write new element and return Nth delayed element
	T operator()(const T& v){
		T r = mElems[N-1];
		write(v);
		return r;
	}

	/// Get size of delay
	static unsigned size(){ return N; }

protected:
	T mElems[N];
};


/// One sample delay. Returns last input sample.

/// \ingroup Delay
///  
template<class T = gam::real>
class Delay1 : public DelayShift<1,T>{
public:

	/// \param[in] v	Initial value of elements
	Delay1(const T& v=T()): DelayShift<1,T>(v){}
};


/// Two sample delay. Returns second to last input sample.

/// \ingroup Delay
///
template<class T = gam::real>
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
/// \tparam Tv	Value (sample) type
/// \tparam Si	Interpolation strategy
/// \tparam Tp	Parameter type
/// \tparam Td	Domain type
/// \ingroup Delay Filter
// H(z) = (ffd + z^-m) / (1 - fbk z^-m)
// y[n] = ffd x[n] + x[n-m] + fbk y[n-m]
template<
	class Tv = gam::real,
	template <class> class Si = ipl::Linear,
	class Tp = gam::real,
	class Td = GAM_DEFAULT_DOMAIN
>
class Comb : public Delay<Tv,Si,Td> {

private:
	typedef Delay<Tv,Si,Td> Base;

public:

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

	Tv operator()();
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

#define TM1 template <class Tv, template <class> class Ti, class Td>
#define TM2 Tv,Ti,Td

#define DELAY_INIT mMaxDelay(0), mDelayFactor(0), mDelayLength(0), mPhase(0), mPhaseInc(0), mDelay(0)

TM1 Delay<TM2>::Delay()
:	DELAY_INIT
{
	onDomainChange(1);
}

TM1 Delay<TM2>::Delay(float maxDly, float dly)
:	DELAY_INIT
{
	onDomainChange(1);
	maxDelay(maxDly, false);
	this->zero();
	delay(dly);
}

TM1 Delay<TM2>::Delay(float dly)
:	DELAY_INIT
{	//printf("Delay::Delay(float)\n");
	onDomainChange(1);
	maxDelay(dly, false);
	this->zero();
	delay(dly);
}

#undef DELAY_INIT

TM1 void Delay<TM2>::maxDelay(float length, bool setDelay){
	//printf("Delay::maxDelay(%f)\n", length);
	mMaxDelay = length;
	if(Td::domain() && Td::domain()->hasBeenSet()){
		//printf("Delay::maxDelay(): resize to %d\n", (unsigned)(mMaxDelay * spu()));
		
		unsigned maxDelayInSamples = unsigned(mMaxDelay * Td::spu());

		// If writing before reading,
		// we must add 1 to support delays at exactly the max delay length.
		//maxDelayInSamples += 1;

		// This will trigger onResize() -> onDomainChange(double r) calls ONLY if
		// the size changes to prevent infinite recursion.
		this->resize(maxDelayInSamples);
	}
	if(setDelay) delay(length);
}

TM1 inline Tv Delay<TM2>::operator()() const {
	return mIpol(*this, mPhase - mDelay);
}

TM1 inline Tv Delay<TM2>::operator()(const Tv& i0){
	// Read, then write.
	Tv res = (*this)();
	write(i0);
	return res;

	// We must write before reading to support delay lengths of zero.
	/*write(i0);
	return (*this)();*/
}

TM1 inline uint32_t Delay<TM2>::delayFToI(float v) const {
	return castIntRound((v * mDelayFactor) * 4294967296.);
	//return scl::unitToUInt(v * mDelayFactor);
}

TM1 inline void Delay<TM2>::incPhase(){ mPhase += mPhaseInc; }

TM1 void Delay<TM2>::onResize(){ //printf("Delay::onResize %d elements\n", this->size());
	mPhaseInc = this->oneIndex();
	//for(uint32_t i=0; i<this->size(); ++i) (*this)[i] = Tv(0);
	if(this->isSoleOwner()) this->zero();
	onDomainChange(1);
}

TM1 void Delay<TM2>::onDomainChange(double /*r*/){ //printf("Delay::onDomainChange\n");
	if(this->usingExternalSource()){
		mMaxDelay = float(this->size() * Td::ups());
	}
	else{
		maxDelay(mMaxDelay, false);
	}

	float currDelay = delay();
	refreshDelayFactor();
	delay(currDelay);
}

TM1 inline Tv Delay<TM2>::read(float ago) const {
	return mIpol(*this, mPhase - delayFToI(ago));
}

TM1 template <template <class> class InterpolationStrategy> inline Tv Delay<TM2>::read(float ago) const {
	return InterpolationStrategy<Tv>()(*this, mPhase - delayFToI(ago));
}

TM1
template <class V>
void Delay<TM2>::read(V * dst, unsigned len, unsigned end) const {
	unsigned mask = this->size()-1;
	unsigned begin = this->index(mPhase) - (end + len);
	for(unsigned i=0; i<len; ++i){
		dst[i] = (*this)[(begin+i)&mask];
	}
}

TM1 void Delay<TM2>::refreshDelayFactor(){ mDelayFactor = 1.0f/maxDelay(); }

TM1 inline void Delay<TM2>::write(const Tv& v){
	//incPhase();
	tbl::put(this->elems(), this->fracBits(), mPhase, v);
	incPhase();
}

TM1 inline void Delay<TM2>::delay(float v){
	mDelayLength = v;
	mDelay = delayFToI(v);
}

TM1 float Delay<TM2>::delay() const { return mDelayLength; }

TM1 inline void Delay<TM2>::delaySamples(uint32_t v){
	mDelay = v << this->fracBits();
	mDelayLength = v * Td::ups();
}

TM1 inline void Delay<TM2>::delaySamplesR(float v){
	delay(v * Td::ups());
}

TM1 uint32_t Delay<TM2>::delaySamples() const {
	return mDelay >> this->fracBits();
}

TM1 float Delay<TM2>::delaySamplesR() const {
	return delay() * Td::spu();
}

TM1 inline void Delay<TM2>::delayUnit(float v){
	mDelay = unitToUInt(v);
}

TM1 float Delay<TM2>::delayUnit() const {
	return uintToUnit<float>(mDelay);
}

TM1 inline uint32_t Delay<TM2>::delayIndex(uint32_t delay) const {
	return this->index(mPhase - (delay << this->fracBits()));
}

TM1 inline void Delay<TM2>::freq(float v){ delay(1.f/v); }

TM1 inline uint32_t Delay<TM2>::indexBack() const {
	return this->index(mPhase + this->oneIndex());
}

TM1 float Delay<TM2>::maxDelay() const { return float(this->size() * Td::ups()); }

TM1 void Delay<TM2>::print(){
	printf("SPU:       %f\n", Td::spu());
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




#define TM1 template<class Tv, template<class> class Si, class Tp, class Td>
#define TM2 Tv,Si,Tp,Td
TM1 Comb<TM2>::Comb()
:	mFFD(0), mFBK(0)
{}

TM1 Comb<TM2>::Comb(float delay, const Tp& ffd, const Tp& fbk)
:	Delay<Tv,Si,Td>(delay), mFFD(ffd), mFBK(fbk)
{}

TM1 Comb<TM2>::Comb(float delayMax, float delay, const Tp& ffd, const Tp& fbk)
:	Delay<Tv,Si,Td>(delayMax, delay), mFFD(ffd), mFBK(fbk)
{}

TM1 inline Tv Comb<TM2>::operator()(){
	return Base::operator()();
}

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
	mFBK = pow(end, this->delay() / scl::abs(units));
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
