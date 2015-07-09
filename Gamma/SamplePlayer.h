#ifndef GAMMA_SAMPLE_PLAYER_H_INC
#define GAMMA_SAMPLE_PLAYER_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include <stdio.h>
#include "Gamma/Containers.h"	// Array
#include "Gamma/ipl.h"
#include "Gamma/scl.h"
#include "Gamma/SoundFile.h"
#include "Gamma/Strategy.h"
#include "Gamma/Domain.h"

namespace gam{

/// Sample buffer player

/// This streams a sequence of frames from a n-channel buffer according to a 
/// specified playback rate. Minimum and maximum endpoint frames are supported
/// for playing back a subinterval of the buffer.
///	The number of frames in the sample should not exceed 2^32. This equates
///	to 27 hours at 44.1 kHz.
///
/// \tparam T	Value (sample) type
/// \tparam Si	Interpolation strategy
/// \tparam Sp	Phase increment strategy
template<
	class T = float,
	template<class> class Si = ipl::Trunc,
	class Sp = phsInc::OneShot
>
class SamplePlayer: public DomainObserver, public Array<T>{
public:
	using Array<T>::size;
	using Array<T>::elems;


	SamplePlayer();

	/// \param[in] src		Another SamplePlayer to read data from
	/// \param[in] rate		Playback rate
	explicit SamplePlayer(SamplePlayer<T>& src, double rate=1);

	/// \param[in] src		Sample array to reference
	/// \param[in] smpRate	Sample rate of samples
	/// \param[in] rate		Playback rate
	SamplePlayer(Array<T>& src, double smpRate, double rate=1);

	/// \param[in] pathToSoundFile		Path to sound file
	/// \param[in] rate					Playback rate
	explicit SamplePlayer(const char * pathToSoundFile, double rate=1);


	/// Load a sound file into internal sample buffer
	
	/// \returns whether the sound file loaded properly
	///
	bool load(const char * pathToSoundFile);


	/// Increment read tap
	void advance();

	/// Returns sample at current position on specified channel and increments phase
	T operator()(int channel=0);

	/// Returns sample at current position on specified channel (without incrementing phase)
	T read(int channel) const;

	/// Set sample buffer
	
	/// \param[in] src		Sample buffer (if multichannel, must be deinterleaved)
	/// \param[in] frmRate	Frame rate of sample buffer.
	///						If the sample is a wavetable, then this should be
	///						its period, in frames.
	/// \param[in] channels	Number of channels in sample buffer
	void buffer(Array<T>& src, double frmRate, int channels);

	/// Set sample buffer
	
	/// \param[in] src		A source SamplePlayer from which to use the same 
	///						samples, sample rate, and channel count
	void buffer(SamplePlayer& src);

	void free();							///< Free sample buffer (if owner)

	void freq(double v){ rate(v); }			///< Set frequency if sample buffer is a wavetable
	void max(double v);						///< Set playback interval maximum frame (open)
	void min(double v);						///< Set playback interval minimum frame (closed)
	void pos(double v);						///< Set playback position, in frames
	void phase(double v);					///< Set playback position, in [0, 1)
	void rate(double v);					///< Set playback rate scalar
	void range(double phs, double period);	///< Set playback interval start phase and period

	void reset();							///< Reset playback head
	void finish();							///< Set playback head to end

	/// Loop playback head if it's past an endpoint
	
	/// This is only applicable for one-shot playback (phsInc::OneShot) and is 
	/// here to provide run-time switchable looping behavior.
	/// \returns whether the playback head was looped
	bool loop();

	/// Apply linear fade-in/-out envelope(s) to frame buffer

	/// \param[in] fadeOutFrames	length of fade out, in frames; <2 for no fade
	/// \param[in] fadeInFrames		length of fade in, in frames; <2 for no fade
	void fade(int fadeOutFrames=4, int fadeInFrames=2);


	/// Returns whether sample playback has completed (non-looping only)
	bool done() const;

	int channels() const { return mChans; }	///< Get number of channels
	int frames() const { return size()/channels(); } ///< Get number of frames (samples divided by channels)
	double frameRate() const { return mFrameRate; } ///< Get frame rate of sample buffer
	double freq() const { return rate(); }	///< Get frequency if sample buffer is a wavetable
	double max() const { return mMax; }		///< Get playback interval maximum frame (open)
	double min() const { return mMin; }		///< Get playback interval minimum frame (closed)
	double period() const;					///< Get total period, in seconds, of sample data
	double pos() const { return mPos; }		///< Get playback position, in frames
	double posInInterval(double frac) const;///< Get position from fraction within interval
	double rate() const { return mRate; }	///< Get playback rate


	void onDomainChange(double r){ frameRate(mFrameRate); }

protected:
	Si<T> mIpol;
	Sp mPhsInc;

	double mPos, mInc;			// real index position and increment
	double mFrameRate;			// frame rate of array data
	int mChans;					// number of channels
	double mRate;				// playback rate factor
	double mMin, mMax;			// [min, max) playback interval, in frames
	
	void frameRate(double v){
		mFrameRate = v;
		rate(mRate);
	}

	T& sample(int idx, int chan){
		return (*this)[frames()*chan + idx];
	}
};



#define PRE template <class T, template<class> class Si, class Sp>
#define CLS SamplePlayer<T,Si,Sp>

PRE CLS::SamplePlayer()
:	Array<T>(defaultArray<T>(), 1),
	mPos(0), mInc(0),
	mFrameRate(1), mChans(1),
	mRate(1), mMin(0), mMax(1)
{}


PRE CLS::SamplePlayer(SamplePlayer<T>& src, double rate)
:	mPos(0), mInc(1), mRate(rate)
{
	buffer(src);
}

PRE CLS::SamplePlayer(Array<T>& src, double smpRate, double rate)
:	mPos(0), mInc(1), mRate(rate)
{
	buffer(src, smpRate, 1);
}


PRE CLS::SamplePlayer(const char * path, double rate)
:	Array<T>(), mPos(0), mInc(1), mChans(1), mRate(rate), mMin(0), mMax(1)
{	
	if(!load(path)){
		this->source(defaultArray<T>(), 1);
	}
}

PRE bool CLS::load(const char * pathToSoundFile){
	SoundFile sf(pathToSoundFile);
	
	if(sf.openRead()){
		Array<T>::resize(sf.samples());
		sf.readAllD(elems());
		frameRate(sf.frameRate());
		mChans = sf.channels();
		mMin = 0;
		mMax = frames();
		mPos = 0;
		sf.close();
		return true;
	}

	fprintf(stderr, 
		"gam::SamplePlayer: couldn't load sound file \"%s\"\n",
		pathToSoundFile);

	return false;
}


PRE inline void CLS::advance(){
	mPos = mPhsInc(pos(), mInc, max(), min()); // update read position, in frames
}

PRE inline T CLS::operator()(int channel){
	T r = read(channel);
	advance();
	return r;
}

PRE inline T CLS::read(int channel) const {
	int posi = int(pos());
	int Nframes= frames();
	int offset = channel*Nframes;
	return mIpol(elems(), posi+offset, pos()-posi, offset+Nframes-1, offset);
}


PRE void CLS::buffer(Array<T>& src, double smpRate, int channels){
	this->source(src);
	frameRate(smpRate);	// sets mFrameRate, mRate, and mInc
	mChans = channels;
	mMin = 0;
	mMax = frames();
	mPos = 0;
}

PRE void CLS::buffer(SamplePlayer& src){
	buffer(src, src.frameRate(), src.channels());
}

PRE inline void CLS::pos(double v){	mPos = v; }

PRE inline void CLS::phase(double v){ pos(v * frames()); }

PRE void CLS::min(double v){ mMin = scl::clip<double>(v, mMax, 0.); }	

PRE void CLS::max(double v){ mMax = scl::clip<double>(v, frames(), mMin); }

PRE void CLS::free(){ this->freeElements(); }

PRE inline void CLS::rate(double v){
	mRate = v;
	mInc = v * frameRate() * ups();
}

PRE void CLS::range(double posn, double period){
	phase(posn);
	min(pos());
	max(pos() + period * spu());
}

PRE void CLS::reset(){
	pos(rate()<0 ? max() : min());
	mPhsInc.reset();
}

PRE void CLS::finish(){
	pos(rate()<0 ? min() : max());
}

PRE inline bool CLS::done() const{
	// The trigger points are based on the logic of phsInc::OneShot
	if(rate() >= 0.){
		return pos() >= (max() - mInc);
	}
	else{
		return pos() <= min();
	}
}

PRE bool CLS::loop(){
	if(done()){
		mPos = phsInc::Loop()(mPos, mInc, max(), min());
		return true;
	}
	return false;
}

PRE void CLS::fade(int fadeOutFrames, int fadeInFrames){
	if(fadeInFrames > 0){
		double amp;
		double slope = 1./(fadeInFrames-1);
		for(int c=0; c<channels(); ++c){
			amp = 0.;
			for(int i=0; i<fadeInFrames; ++i){
				sample(i,c) *= amp;
				amp += slope;	
			}
		}
	}

	if(fadeOutFrames > 0){
		double amp;
		double slope =-1./(fadeOutFrames-1);
		for(int c=0; c<channels(); ++c){
			amp = 1;
			for(int i=frames()-1-fadeOutFrames; i<frames(); ++i){
				sample(i,c) *= amp;
				amp += slope;	
			}
		}
	}	
}

PRE inline double CLS::period() const { return frames() * ups(); }

PRE inline double CLS::posInInterval(double frac) const {
	return min() + (max() - min()) * frac;
}

#undef PRE
#undef CLS

} // gam::

#endif

