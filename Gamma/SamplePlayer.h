#ifndef GAMMA_SAMPLE_PLAYER_H_INC
#define GAMMA_SAMPLE_PLAYER_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include "Gamma/Containers.h"
#include "Gamma/ipl.h"
#include "Gamma/scl.h"
#include "Gamma/SoundFile.h"
#include "Gamma/Strategy.h"
#include "Gamma/Sync.h"
#include "Gamma/Types.h"

namespace gam{


/// Sample buffer player

///	The number of frames in the sample should not exceed 2^32. This equates
///	to 27 hours at 44.1 kHz.
///
/// \tparam T	Value (sample) type
/// \tparam Si	Interpolation strategy
/// \tparam Sp	Phase increment strategy
template<
	class T = real,
	template<class> class Si = ipl::Trunc,
	class Sp = phsInc::Clip
>
class SamplePlayer: public Synced, public Array<T>{
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
	template<class Char>
	explicit SamplePlayer(const Char * pathToSoundFile, double rate=1);


	/// Load a sound file into internal sample buffer
	
	/// \returns whether the sound file loaded properly
	///
	template<class Char>
	bool load(const Char * pathToSoundFile);


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
	void max(double v);						///< Set playback interval max frame (open)
	void min(double v);						///< Set playback interval min frame (closed)
	void pos(double v);						///< Set current read position, in frames
	void phase(double v);					///< Set current read position [0, 1)
	void rate(double v);					///< Set playback rate scalar
	void range(double phs, double period);	///< Set interval start phase and period
	void reset();							///< Reset playback head


	/// Whether sample playback has completed (non-looping only)
	bool done() const;

	int channels() const { return mChans; }	///< Get number of channels
	double frameRate() const { return mFrameRate; } ///< Get frame rate of sample buffer
	double freq() const { return rate(); }	///< Get frequency if sample buffer is a wavetable
	double max() const { return mMax; }		///< Get playback interval max frame (open)
	double min() const { return mMin; }		///< Get playback interval min frame (closed)
	double period() const;					///< Get total period of sample data
	double pos() const { return mPos; }		///< Get current read position, in frames
	double posInInterval(double frac) const;///< Get position from fraction within interval
	double rate() const { return mRate; }	///< Get playback rate

	/// Get whether the sample buffer is valid for playback
	bool valid() const;

	virtual void onResync(double r){ frameRate(mFrameRate); }

protected:	
	static T * defaultBuffer(){
		static T v(0);
		return &v;
	}

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

	int frames() const { return size()/channels(); }
};



#define PRE template <class T, template<class> class Si, class Sp>
#define CLS SamplePlayer<T,Si,Sp>

PRE CLS::SamplePlayer()
:	Array<T>(defaultBuffer(), 1),
	mPos(0), mInc(0),
	mFrameRate(1), mChans(1),
	mRate(1), mMin(0), mMax(1)
{}


PRE CLS::SamplePlayer(SamplePlayer<T>& src, double rate)
:	Array<T>(src), 
	mPos(0), mInc(1), 
	mFrameRate(src.frameRate()), mChans(src.channels()), 
	mRate(rate), mMin(0), mMax(src.size())
{	initSynced(); }

PRE CLS::SamplePlayer(Array<T>& src, double smpRate, double rate)
:	Array<T>(src),
	mPos(0), mInc(1),
	mFrameRate(smpRate), mChans(1),
	mRate(rate), mMin(0), mMax(src.size())
{
	initSynced();
	frameRate(smpRate);
}

PRE
template<class Char>
CLS::SamplePlayer(const Char * path, double rate)
:	Array<T>(), mPos(0), mInc(1), mChans(1), mRate(rate), mMin(0), mMax(1)
{	
	if(!load(path)){
		this->source(defaultBuffer(), 1);
	}
}

PRE
template<class Char>
bool CLS::load(const Char * pathToSoundFile){
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
	uint32_t posi = uint32_t(pos());
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

PRE inline bool CLS::done() const{
	if(rate() >= 0.){
		return pos() >= (max() - 1);
	}
	else{
		return pos() <= min();
	}
}

PRE inline double CLS::period() const { return frames() * ups(); }

PRE inline double CLS::posInInterval(double frac) const {
	return min() + (max() - min()) * frac;
}

PRE bool CLS::valid() const {
	return this->mElems && (this->mElems != defaultBuffer());
}

#undef PRE
#undef CLS

} // gam::

#endif

