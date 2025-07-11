#ifndef GAMMA_SPATIAL_H_INC
#define GAMMA_SPATIAL_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information

	File Description: 
	Sound spatialization
*/

#include <cmath> // min, max, sqrt
#include <initializer_list>
#include <vector>
#include "Gamma/scl.h"
#include "Gamma/Types.h"
#include "Gamma/Delay.h"
#include "Gamma/Filter.h"
#include "Gamma/Ramped.h"


namespace gam{

/// Simulate aspects of sound propagation through space

/// \defgroup Spatial


/// Gain loop filter (all-pass)
template <typename T>
class LoopGain{
public:

	/// Set filter gain
	LoopGain& gain(float v){
		mA0 = v;
		return *this;
	}

	LoopGain& damping(float v){ return *this; }

	T operator()(T in){
		return in*mA0;
	}

	/// Get filter gain
	float gain() const { return mA0; }

private:
	float mA0=0.f;
};


/// One-pole loop filter (low-/high-pass)
template <typename T>
class Loop1P{
public:

	/// Set filter gain
	Loop1P& gain(float v){
		mA0 = (1.f - scl::abs(mB1))*v;
		return *this;
	}

	/// Set damping amount

	/// \param[in] v	Damping amount in (-1,1) where
	///					zero is a no-op,
	///					positive values produce a low-pass and
	///					negative values produce a high-pass.
	Loop1P& damping(float v){
		float g = gain();
		const float vMax = 0.99999f;
		v = std::min(std::max(v, -vMax), vMax); // limit to (-1,1)
		mB1 = v;
		return gain(g);
	}

	T operator()(T in){
		return mO1 = in*mA0 + mO1*mB1;
	}

	/// Get filter gain
	float gain() const {
		return mA0 / (1.f - scl::abs(mB1));
	}

private:
	float mA0=0.f, mB1=0.f;
	T mO1=T(0);
};


/// One-pole/one-zero loop filter (low-pass)
template <typename T>
class Loop1P1Z{
public:

	/// Set filter gain
	Loop1P1Z& gain(float v){
		mA0 = (mB1*0.5f + 0.5f)*v;
		return *this;
	}

	/// Set damping amount in [0, 1)
	Loop1P1Z& damping(float v){
		float g = gain();
		const float vMax = 0.99999f;
		v = std::min(std::max(v, 0.f), vMax); // limit to [0, 1)
		mB1 = 1.f-2.f*v; // [0, 1) -> [1, -1)
		return gain(g);
	}

	T operator()(T in){
		mO1 = (in + mI1)*mA0 - mO1*mB1; // low-pass
		//mO1 = (in - mI1)*mA0 - mO1*mB1; // high-pass
		mI1 = in;
		return mO1;
	}

	/// Get filter gain
	float gain() const {
		return 2.f * mA0 / (mB1 + 1.f);
	}

private:
	float mA0=0.f, mB1=0.f;
	T mI1=T(0), mO1=T(0);
};



/// Recursive echo with loop filter

/// \tparam Tv			Value (sample) type
/// \tparam Si			Interpolation strategy
/// \tparam LoopFilter	Filter to insert in feedback loop
/// \tparam Td			Domain type
/// \ingroup Spatial
template<
	typename Tv = gam::real,
	template<typename> class Si = ipl::Linear,
	template<typename> class LoopFilter = LoopGain,
	class Td = GAM_DEFAULT_DOMAIN
>
class Echo : public Delay<Tv,Si,Td> {
public:

	typedef Delay<Tv,Si,Td> Base;


	Echo();

	/// \param[in] delay	delay length
	/// \param[in] decay	decay length
	/// \param[in] damp		damping factor of loop filter (if applicable)
	Echo(float delay, float decay=1.f, float damp=0.f);


	/// Set decay length
	Echo& decay(float v);

	/// Set feedback coefficient, in [-1,1]
	Echo& fbk(float v);

	/// Set damping factor of loop filter
	Echo& damping(float v);

	/// Process next sample
	Tv operator()(Tv in);

	/// Process next sample

	/// The output tap is just after the loop filter after the delay, thus the
	/// first echo is filtered.
	/// \returns wet (delayed) sample.
	Tv nextPost(Tv in);

	/// Process next sample; output tap is before loop filter

	/// The output tap is just after the delay before the loop filter, thus the
	/// first echo is not filtered.
	/// \returns wet (delayed) sample.
	Tv nextPre(Tv in);

	/// Get decay length
	float decay() const { return mDecay; }

	/// Get loop filter
	LoopFilter<Tv>& loopFilter(){ return mFilter; }

protected:
	float mDecay=1.f;
	LoopFilter<Tv> mFilter;

	virtual void onDomainChange(double r){
		Delay<Tv,Si,Td>::onDomainChange(r);
		decay(mDecay);
	}
};


/// Echo with damped complex sinusoidal response
template<
	typename Tv = gam::real,
	template<typename> class Si = ipl::Linear,
	template<typename> class LoopFilter = LoopGain,
	class Td = GAM_DEFAULT_DOMAIN
>
class EchoCSine : public Delay<Complex<Tv>, Si, Td> {
public:

	typedef Delay<Complex<Tv>, Si, Td> Base;


	EchoCSine();

	/// \param[in] delay	delay length (and max delay)
	EchoCSine(double delay);


	/// Set feedback factor
	EchoCSine& fbk(float amt, float ang=0);

	/// Set gain factor
	EchoCSine& gain(float amt, float ang=0);

	/// Set decay length
	EchoCSine& decay(float units);

	/// Set feedback oscillation frequency
	EchoCSine& fbkFreq(float frq, float addCycle=0);


	/// Filter next sample
	Complex<Tv> operator()(Tv real, Tv imag);
	Complex<Tv> operator()(Tv v);


	/// Get feedback factor
	Complex<float> fbk() const { return mB; }

	/// Get gain factor
	Complex<float> gain() const { return mG; }

private:
	Complex<Tv> mB{0.5,0}, mG{1,0};
	float mFbkFreq=0.f;

	virtual void onDomainChange(double r){
		Delay<Complex<Tv>,Si,Td>::onDomainChange(r);
		if(mFbkFreq>0.f) fbkFreq(mFbkFreq);
	}
};



enum ReverbFlavor{
	JCREVERB,	// John Chowning (4-comb, 3-allpass)
	FREEVERB,	// Jezar's Freeverb (8-comb, 4-allpass)
};


/// Schroeder reverberator

/// This simulates the late reflections in a reverberant space using an
/// algorithm devised by Manfred Schroeder. The network consists of a group of
/// parallel comb filters followed by a series of allpass comb filters.
/// The comb filters produce initial echoes and the allpass combs diffuse the
/// echoes by smearing transients and further increasing echo density.
/// Note that using long decay lengths and/or highly "tonal" inputs can lead to 
/// very metallic sounding responses due to the fixed resonances of the comb 
/// filters.
///
/// \tparam Tv			Value (sample) type
/// \tparam LoopFilter	Filter to insert in comb feedback loop
/// \tparam Si			Interpolation strategy
/// \tparam Td			Domain type
/// \ingroup Spatial
template<
	typename Tv = gam::real,
	template<typename> class LoopFilter = Loop1P,
	template<typename> class Si = ipl::Trunc,
	class Td = GAM_DEFAULT_DOMAIN
>
class ReverbMS : public Td {
public:

	typedef std::vector<Echo<Tv, Si, LoopFilter, Domain1>> Combs;
	typedef std::vector<Comb<Tv, Si, float, Domain1>> Allpasses;

	ReverbMS();

	ReverbMS(ReverbFlavor flavor, unsigned offset=0);

	ReverbMS(
		std::initializer_list<unsigned> combDelays,
		std::initializer_list<unsigned> allpassDelays,
		unsigned offset = 0
	);


	/// Resize delay lines based on a particular flavor of reverb

	/// \param[in] flavor	reverb flavor/algorithm
	/// \param[in] offset	offset to add to nomimal delay lengths, in samples
	ReverbMS& resize(ReverbFlavor flavor, unsigned offset=0);

	/// Resize delay lines

	/// \param[in] combDelays		comb delay sizes in samples
	/// \param[in] allpassDelays	allpass delay sizes in samples
	/// \param[in] offset			offset to add to delay lengths, in samples
	ReverbMS& resize(
		std::initializer_list<unsigned> combDelays,
		std::initializer_list<unsigned> allpassDelays,
		unsigned offset = 0
	);

	/// Set comb delay sizes in samples
	ReverbMS& resizeComb(std::initializer_list<unsigned> delays, unsigned offset=0);

	/// Set allpass delay sizes in samples
	ReverbMS& resizeAllpass(std::initializer_list<unsigned> delays, unsigned offset=0);

	/// Set decay length
	ReverbMS& decay(float v);

	/// Set damping factor
	ReverbMS& damping(float v);

	/// Filter next sample
	Tv operator()(Tv in);

	/// Get sum of comb delay taps
	Tv read(std::initializer_list<unsigned> delays) const;


	/// Get decay length
	float decay() const { return mDecay; }

	Combs& combs(){ return mCombs; }
	Allpasses& allpasses(){ return mAllpasses; }

	/// Check if valid for processing
	operator bool() const { return mCombs.size(); }

	void print() const;

private:
	float mDecay;
	Combs mCombs;
	Allpasses mAllpasses;

	virtual void onDomainChange(double r){
		decay(decay());
	}
};



/// Spatializes a source at one or more destinations

/// This effectively samples the wave field produced by a single source at
/// multiple points in space via distance cues---amplitude attenuation, time
/// delay, and high-frequency damping (air absorption).
///
/// \tparam Ndest	number of destination/sample points
/// \tparam T		value type to process
template <int Ndest=2, class T = gam::real>
class Dist{
public:

	/// \param[in] maxDelay		maximum delay interval
	/// \param[in] near			near clipping distance
	/// \param[in] far			far clipping distance
	Dist(float maxDelay=0.2, float near=0.1, float far=10);

	/// Set number of samples over which to smooth coefficients
	Dist& blockSize(int v){ mBlockSize=v; return *this; }

	/// Resets internal coefficient smoothers

	/// Call to re-initialize the object to avoid interpolating from the old 
	/// state to the new state, like when acquiring the object from a pool.
	Dist& reset(){ mReset=true; return *this; }

	/// Set near clipping distance
	Dist& near(float v){ mNear=v; return updateCoefs(); }
	float near() const { return mNear; }

	/// Set far clipping distance
	Dist& far(float v){ mFar=v; return updateCoefs(); }
	float far() const { return mFar; }

	/// Set gain at far clipping distance (must be less than 1)
	Dist& farGain(float v){ mFarGain=v; return updateCoefs(); }
	float farGain() const { return mFarGain; }

	/// Set whether attenuation goes down to zero at far clip

	/// If true, far gain acts as roll-off where [0,1) -> [impulse, line).
	///
	Dist& toZero(bool v){ mToZero=v; return updateCoefs(); }

	/// Set all attenuation function parameters (more efficient)
	Dist& set(float n, float f, float fgain){
		mNear=n; mFar=f; mFarGain=fgain;
		return updateCoefs();
	}
	Dist& set(float n, float f, float fgain, bool toZero){
		mToZero=toZero; return set(n,f, fgain);
	}

	/// Set speed of sound
	Dist& speedOfSound(float v){ mInvSpeedOfSound=1./v; return *this; }
	float speedOfSound() const { return 1./mInvSpeedOfSound; }

	/// Set distance from source to a destination, in meters
	Dist& dist(int dest, float d);

	/// Set distance vector from source to a destination, in meters
	Dist& dist(int dest, float x, float y, float z=0);

	/// Set distances from source to destinations, in meters
	template <typename V>
	Dist& dist(const Vec<Ndest,V>& d);

	/// Get distances from source to destinations, in meters
	const Vec<Ndest,float>& dist() const { return mDist; }

	/// Filter source sound
	Vec<Ndest, T> operator()(T in);

	/// Get delay line
	const Delay<T>& delayLine() const { return mDelay; }

	/// Returns attenuation for distance, based on current settings
	float inverse(float dist) const;

private:
	Delay<T> mDelay;
	float mNear, mFar;
	float mFarGain = 0.25;
	float mA, mB, mC; // coefs for attenuation function
	float mInvSpeedOfSound = 1./343.2;
	float mBlockSize = 64.;

	Vec<Ndest,float> mDist;
	Ramped<T> mDly[Ndest];
	Ramped<T> mAmp[Ndest];
	OnePole<T> mLPF[Ndest];
	bool mToZero = false;
	bool mReset = true; // used to reset ramps

	Dist& updateCoefs();
};



// Implementation_______________________________________________________________

namespace detail{

template <typename T>
inline T decayToFbk(T decay, T delay){
	return pow(0.001, delay/decay);
}

template <typename T>
inline T fbkToDecay(T fbk, T delay){
	return T(-3)*delay / log(scl::abs(fbk));
}

}

#define TDEC\
	typename Tv,\
	template<typename> class Si,\
	template<typename> class LoopFilter,\
	class Td

#define TARG Tv,Si,LoopFilter,Td

template<TDEC>
Echo<TARG>::Echo(){}

template<TDEC>
Echo<TARG>::Echo(float delay, float decay, float damp)
:	Base(delay)
{
	this->decay(decay);
	this->damping(damp);
}

template<TDEC>
Echo<TARG>& Echo<TARG>::decay(float v){
	mDecay = v;
	float fbk = detail::decayToFbk(mDecay, this->delay());
	mFilter.gain(fbk);
	return *this;
}

template<TDEC>
Echo<TARG>& Echo<TARG>::fbk(float v){
	mDecay = detail::fbkToDecay(v, this->delay());
	mFilter.gain(v);
	return *this;
}

template<TDEC>
Echo<TARG>& Echo<TARG>::damping(float v){
	mFilter.damping(v);
	return *this;
}

template<TDEC>
inline Tv Echo<TARG>::nextPost(Tv in){
	auto out = Base::read();
	out = mFilter(out);
	Base::write(in + out);
	return out;
}

template<TDEC>
inline Tv Echo<TARG>::nextPre(Tv in){
	auto out = Base::read();
	Base::write(in + mFilter(out));
	return out;
}

template<TDEC>
inline Tv Echo<TARG>::operator()(Tv in){
	return nextPre(in);
}



template<TDEC>
EchoCSine<TARG>::EchoCSine(){}

template<TDEC>
EchoCSine<TARG>::EchoCSine(double delay)
:	Base(delay)
{}

template<TDEC>
EchoCSine<TARG>& EchoCSine<TARG>::gain(float amt, float ang){
	mG.fromPolar(amt, ang*2*M_PI);
	return *this;
}

template<TDEC>
EchoCSine<TARG>& EchoCSine<TARG>::fbk(float amt, float ang){
	mFbkFreq = 0.f;
	mB.fromPolar(amt, ang*2*M_PI);
	return *this;
}

template<TDEC>
EchoCSine<TARG>& EchoCSine<TARG>::decay(float units){
	Tv fbk = detail::decayToFbk(units, this->delay());
	mB.mag(fbk);
	return *this;
}

template<TDEC>
EchoCSine<TARG>& EchoCSine<TARG>::fbkFreq(float frq, float addCycle){
	mFbkFreq = frq;
	/*
	d		seconds/echo
	frq		cycles/second
	*/
	mB.arg((frq*this->delay() + addCycle)*2*M_PI);
	return *this;
}

template<TDEC>
inline Complex<Tv> EchoCSine<TARG>::operator()(Tv real, Tv imag){
	Complex<Tv> echo = Base::operator()()*mB;
	return Base::operator()(Complex<Tv>(real,imag) + echo) * mG;
}

template<TDEC>
inline Complex<Tv> EchoCSine<TARG>::operator()(Tv v){
	return (*this)(v, Tv(0));
}

#undef TDEC
#undef TARG

#define TDEC\
	typename Tv,\
	template<typename> class LoopFilter,\
	template<typename> class Si,\
	class Td

#define TARG Tv,LoopFilter,Si,Td

template<TDEC>
ReverbMS<TARG>::ReverbMS(){
	decay(1);
}

template<TDEC>
ReverbMS<TARG>::ReverbMS(ReverbFlavor flavor, unsigned offset)
:	ReverbMS()
{
	resize(flavor, offset);
}

template<TDEC>
ReverbMS<TARG>::ReverbMS(
	std::initializer_list<unsigned> combDelays,
	std::initializer_list<unsigned> allpassDelays,
	unsigned offset
)
:	ReverbMS()
{
	resize(combDelays, allpassDelays, offset);
}

template<TDEC>
ReverbMS<TARG>& ReverbMS<TARG>::resizeComb(std::initializer_list<unsigned> delays, unsigned offset){
	if(delays.size() != mCombs.size()){
		mCombs.resize(delays.size());
	}
	for(unsigned i=0; i<mCombs.size(); ++i){
		unsigned d = delays.begin()[i] + offset;
		mCombs[i].maxDelay(d);
		mCombs[i].delay(d);
	}
	decay(decay());
	return *this;
}

template<TDEC>
ReverbMS<TARG>& ReverbMS<TARG>::resizeAllpass(std::initializer_list<unsigned> delays, unsigned offset){
	if(delays.size() != mAllpasses.size()){
		mAllpasses.resize(delays.size());
	}
	for(unsigned i=0; i<mAllpasses.size(); ++i){
		unsigned d = delays.begin()[i] + offset;
		mAllpasses[i].maxDelay(d);
		mAllpasses[i].delay(d);
		mAllpasses[i].allPass(0.71); // all use the same feedback amount
	}
	return *this;
}

template<TDEC>
ReverbMS<TARG>& ReverbMS<TARG>::resize(
	std::initializer_list<unsigned> combDelays,
	std::initializer_list<unsigned> allpassDelays,
	unsigned offset
){
	return resizeComb(combDelays,offset).resizeAllpass(allpassDelays,offset);
}

template<TDEC>
ReverbMS<TARG>& ReverbMS<TARG>::resize(ReverbFlavor flavor, unsigned d){
	switch(flavor){
	case JCREVERB:
		return resize(
			{1307+d, 1637+d, 1811+d, 1931+d},
			{1051+d, 337+d, 113+d}
		);
	case FREEVERB:
		return resize(
			{1116+d, 1188+d, 1277+d, 1356+d, 1422+d, 1491+d, 1557+d, 1617+d},
			{556+d, 441+d, 341+d, 225+d}
		);
	default:;
	}
	return *this;
}

template<TDEC>
ReverbMS<TARG>& ReverbMS<TARG>::decay(float v){
	mDecay = v;
	float decaySamples = v * Td::spu();
	for(unsigned i=0; i<mCombs.size(); ++i) mCombs[i].decay(decaySamples);
	return *this;
}

template<TDEC>
ReverbMS<TARG>& ReverbMS<TARG>::damping(float v){
	for(unsigned i=0; i<mCombs.size(); ++i) mCombs[i].damping(v);
	return *this;
}

template<TDEC>
inline Tv ReverbMS<TARG>::operator()(Tv in){

	// Series allpasses
	for(unsigned i=0; i<mAllpasses.size(); ++i){
		in = mAllpasses[i](in);
	}

	// Parallel combs
	Tv res = Tv(0);
	for(unsigned i=0; i<mCombs.size(); ++i){
		res += mCombs[i](in);
	}

	return res;
}

template<TDEC>
inline Tv ReverbMS<TARG>::read(std::initializer_list<unsigned> delays) const {
	Tv res = Tv(0);
	for(unsigned i=0; i<mCombs.size(); ++i){
		res += mCombs[i].read(delays.begin()[i]);
	}
	return res;
}

template<TDEC>
void ReverbMS<TARG>::print() const {
	unsigned Nc = mCombs.size();
	unsigned Na = mAllpasses.size();
	printf("comb delays = {");
	for(unsigned i=0; i<Nc; ++i)
		printf("%u%s", mCombs[i].delaySamples(), i!=(Nc-1)?", ":"");
	printf("} samples\n");
	printf("allpass delays = {");
	for(unsigned i=0; i<Na; ++i)
		printf("%u%s", mAllpasses[i].delaySamples(), i!=(Na-1)?", ":"");
	printf("} samples\n");
}

#undef TDEC
#undef TARG

#define TDEC int Ndest, class T
#define TARG Ndest, T

template<TDEC>
Dist<TARG>::Dist(float maxDelay, float near, float far)
:	mDelay(maxDelay), mNear(near), mFar(far)
{
	for(int i=0; i<Ndest; ++i){
		mDist[i] = 1e8;
		mAmp[i] = 0;
		mDly[i] = 0.001;
	}
	updateCoefs();
}

template<TDEC>
Dist<TARG>& Dist<TARG>::dist(int i, float d){
	mDist[i] = d;
	auto amp = inverse(d);
	auto dly = std::min(d * mInvSpeedOfSound, mDelay.maxDelay());
	if(mReset){
		mReset = false;
		mAmp[i] = amp;
		mDly[i] = dly;
	}
	mAmp[i].target(amp, mBlockSize);
	mDly[i].target(dly, mBlockSize);
	mLPF[i].freq(22000 * amp + 20.); // low-pass gate
	return *this;
}

template<TDEC>
Dist<TARG>& Dist<TARG>::dist(int i, float x, float y, float z){
	auto d = std::sqrt(x*x+y*y+z*z);
	return dist(i, d);
}

template<TDEC>
template <typename V>
Dist<TARG>& Dist<TARG>::dist(const Vec<Ndest,V>& d){
	for(unsigned i=0; i<Ndest; ++i) dist(i, d[i]);
	return *this;
}

template<TDEC>
inline Vec<Ndest, T> Dist<TARG>::operator()(T in){
	mDelay.write(in);
	Vec<Ndest, T> res;
	for(int i=0; i<Ndest; ++i)
		res[i] = mLPF[i](mDelay.read(mDly[i]())) * mAmp[i]();
	return res;
}

template<TDEC>
Dist<TARG>& Dist<TARG>::updateCoefs(){
	/*
	A(d) = 
		n / [n + r (d - n)]
		n / [n + r d - r n]
		(n/r) / (n/r - n + d)
	*/
	auto rollOff = (mNear/mFarGain - mNear) / (mFar - mNear);
	mA = mNear / rollOff;
	mB = mA - mNear;

	if(mToZero){
		auto toZeroGain = 1./(1.-mFarGain);
		mA *= toZeroGain;
		mC = mFarGain * toZeroGain;
	} else {
		mC = 0.;
	}
	return *this;
}

template<TDEC>
float Dist<TARG>::inverse(float dist) const {
	if(dist <= mNear) return 1.f;
	return std::max(mA / (mB + dist) - mC, 0.f);
}

#undef TDEC
#undef TARG

} // gam::
#endif

