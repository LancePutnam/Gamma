#ifndef GAMMA_DFT_H_INC
#define GAMMA_DFT_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include <math.h>
#include "mem.h"
#include "scl.h"
#include "tbl.h"
#include "Sync.h"
#include "Constants.h"
#include "Containers.h"
#include "FFT.h"
#include "Types.h"


namespace gam{

namespace Bin{

	/// Bin format types.
	enum t{
		Rect=0,		/**< Rectangular (Cartesian)*/
		Polar,		/**< Polar */
		MagFreq		/**< Magnitude/Frequency */
	};

	/// Returns bin type as human readable string.
	inline const char * string(int type){
		switch(type){
		case Bin::Rect:		return "Rectangular";
		case Bin::Polar:	return "Polar";
		case Bin::MagFreq:	return "Magnitude/Frequency";
		default:			return "Unknown\n";
		}	
	}

};



/// Sliding window for analysis
template <class T=gam::real>
class SlidingWindow{
public:
	SlidingWindow(uint winSize, uint hopSize);
	~SlidingWindow();

	void resize(uint winSize, uint hopSize);
	void sizeHop(uint size);
	void sizeWin(uint size);
	
	uint sizeHop() const;
	uint sizeWin() const;
	
	/// Returns pointer to internal sample window
	
	/// The samples returned cannot be modified directly since they point to an
	/// internal delay line.
	const T * window();
	const T * operator()();
	
	/// Returns true when sample window is ready to be processed
	bool operator()(T input);
	
	/// Returns true when sample window is ready to be processed.
	
	/// Upon returning true, the sample buffer is copied into 'dst'.
	/// 'dst' must have a size of at least sizeWin().
	bool operator()(T * dst, T input);

protected:
	void slide();	// Slides samples in window left by hop num samples
	// After processing the window samples, a call to slide() must be made
	// to prepare for the next window.

protected:
	T * mBuf;
	uint mSizeWin, mSizeHop;
	uint mTapW;	// current index to write to
	uint mHopCnt;	// counts samples for hop
	
private:
	uint hopStart() const;
};


/*
if(oa()){
	oa(buf);
}


TEM inline bool OverlapAdd<T>::operator()(T * input){	
	if(++mHopCnt == sizeHop()){
		overlapAdd(input);
		mHopCnt = 0;
	}
}

TEM void OverlapAdd<T>::overlapAdd(T * input){
	ArrOp::overlapAdd(mBuf, input, sizeWin(), sizeHop());
}

*/

template <class T=gam::real>
class DFTBase : public Synced{
public:
	DFTBase();
	virtual ~DFTBase();

	T * aux(uint num);		///< Returns pointer to an auxiliary buffer
//	T *& bins0();			///< Returns pointer to first elements of spectral bins
//	T& bins0(uint i);		///< Returns ith element of first spectral type
//	T * bins1();			///< Returns pointer to second elements of spectral bins
	
	Complex<T> * bins() const { return mBins; }
	Complex<T>& bins(uint32_t i) const { return mBins[i]; }
	
	double binFreq() const;	///< Returns width of frequency bins
	uint numBins() const;	///< Returns number of frequency bins
	uint sizeDFT() const;	///< Returns size of forward transform
	Sync& syncFreq();		///< Returns frequency domain synchronizer
	
	void numAux(uint num);	///< Sets number of auxilliary buffers, each of size numBins()

	virtual void onResync(double r);

protected:
	uint mSizeDFT, mNumAux;
	union{
		T * mBuf;		// FFT buffer
		Complex<T> * mBins;
	};
	
	T * mAux;		// aux buffers
	Sync mSyncFreq;
	T normForward() const;	// Returns norm factor for forward transform values
};


/// Discrete Fourier transform.

/// Requires FFTW library (http://www.fftw.org/) in single-precision (32-bit) mode.
///
class DFT : public DFTBase<float>{
public:
	/// Constructor.

	/// @param[in]	winSize		Number of samples in window.
	/// @param[in]	padSize		Number of zeros to append to window.
	/// @param[in]	complexType	Form of complex spectral data.
	/// @param[in]	numAux		Number of auxilliary buffers of size numBins() to create
	DFT(uint winSize, uint padSize=0, Bin::t t = Bin::Rect, uint numAux=0);
	virtual ~DFT();

	void binType(Bin::t t);						///< Sets format of complex spectrum data.
	DFT& precise(bool whether);					///< Whether to use precise (but slower) math. Default is off.
	void resize(uint windowSize, uint padSize);	///< Sets size parameters of transform.

	float freqRes() const;		///< Returns frequency resolution of analysis.
	float overlap() const;		///< Returns degree of transform overlap
	bool overlapping() const;	///< Whether the xform is overlapping
	uint sizeHop() const;		///< Returns size of hop
	uint sizePad() const;		///< Returns size of zero-padding
	uint sizeWin() const;		///< Returns size of window

	Sync& syncHop();	///< Hop rate synchronizer

	/// Reads next sample in for a forward transform.
	
	/// Returns true when sizeDFT() samples are read in and, subsequently,
	/// the forward DFT is performed.  Returns false otherwise.
	bool operator()(float input);


	/// Returns next sample from inverse transform.
	
	/// The inverse transform is performed every sizeWin() samples.
	///
	float operator()();


	/// Performs forward transform on a window of samples.
	
	/// 'src' must have at least sizeWin() number of elements.
	///
	void forward(const float * src);	

	
	/// Performs inverse transform on internal spectrum.
	
	///	The resynthesized samples are copied into 'dst'.  The destination
	/// array must have room for at least sizeDFT() number of elements. If 'dst'
	/// equals 0, then the resynthesized samples are not copied, but instead
	/// held in the internal transform buffer.
	virtual void inverse(float * dst);
	
	/// Returns true if next call to inverse() will perform the inverse transform.
	
	/// This method is used for doing inverse-only transforms.
	/// Basically, it tells you when you should set the frequency samples.
	bool inverseOnNext();

	void spctToRect();		// convert spectrum to rectangle format
	void spctToPolar();		// convert spectrum to polar format
	void zero();			///< Zeroes internal frequency bins
	void zeroEnds();		///< Zeroes DC and Nyquist bins

	virtual void onResync(double r);
	virtual void print(FILE * fp=stdout, const char * append="\n");
	
protected:
	uint mSizeWin;					// samples in analysis window
	uint mSizeHop;					// samples between forward transforms (= winSize() for DFT)
	Bin::t mSpctFormat;				// complex format of spectrum
	//FFTInfo mInfoFFT, mInfoIFFT;	// info for FFT
	RFFT mFFT;
	Sync mSyncHop;
	
	// Buffers
	float * mPadOA;			// Overlap-add buffer (alloc'ed only if zero-padded)
	float * mBufInv;		// Pointer to inverse sample buffer
	uint mTapW, mTapR;		// DFT i/o read/write taps
	bool mPrecise;

	// Magnitude normalization for inverse transform
	float normInverse() const;
};




/// Short-time Fourier transform.

/// Requires FFTW library (http://www.fftw.org/) in single-precision (32-bit) mode.
///
class STFT : public DFT {
public:

	/// Constructor.
	
	/// The default complex data type is rectangular and the default window
	/// is none.
	/// @param[in]	winSize		Number of samples to window.
	/// @param[in]	hopSize		Number of samples between successive windows.
	/// @param[in]	padSize		Number of zeros to append to window.
	/// @param[in]	winType		Type of forward transform window.
	/// @param[in]	complexType	Complex form of spectral data.
	/// @param[in]	createAux	Whether to create an auxillary spectral buffer (see createAux()).
	STFT(uint winSize, uint hopSize, uint padSize=0,
		WinType::type winType = WinType::Rectangle,
		Bin::t t = Bin::Rect,
		uint numAux=0);
	
	virtual ~STFT();
	
	using DFT::operator();

	bool  operator()(float input);
	//float operator()();
	
	void forward(float * input);
	virtual void inverse(float * dst);
	
	void resize(uint winSize, uint padSize);					///< Sets size parameters of transform.
	void sizeHop(uint size);
	void winType(WinType::type type);
	void rotateForward(bool v){ mRotateForward = v; }
	//void resetPhaseAccum(){  }

	float unitsHop();
	
	float * phases();
	
	virtual void print(FILE * fp=stdout, const char * append="\n");

protected:
	using DFT::sizeHop;
	
	void computeInvWinMul();	// compute inverse normalization factor (due to overlap-add)
	
	SlidingWindow<float> mSlide;
	float * mFwdWin;			// forward transform window
	float * mPhases;			// copy of current phases (mag-freq mode)
	WinType::type mWinType;		// type of analysis window used
	float mFwdWinMul, mInvWinMul;
	bool mWindowInverse;		// whether to window inverse samples
	bool mRotateForward;
};




/// Sliding discrete Fourier transform
template <class T>
class SDFT : public DFTBase<T> {
public:
	SDFT(uint sizeDFT, uint binLo, uint binHi);
	
	void forward(T input);
	void range(uint binLo, uint binHi);
	void resize(uint sizeDFT, uint binLo, uint binHi);
		
protected:
	uint mBinLo, mBinHi;
	DelayN<T> mDelay;
	T * mBufI;				// alias to imaginaries
	double mC0, mS0;		// phasor rotators
	double mCL, mSL;		// low bin phasor
	T mNorm;				// fwd transform normalization
};



// Implementation_______________________________________________________________

//---- SlidingWindow
#define TEM template<class T>
TEM SlidingWindow<T>::SlidingWindow(uint winSize, uint hopSize)
	: mBuf(0), mSizeWin(0), mSizeHop(0), mHopCnt(0)
{
	resize(winSize, hopSize);
	mem::zero(mBuf, sizeWin());
}

TEM SlidingWindow<T>::~SlidingWindow(){
	if(mBuf){ free(mBuf); mBuf = 0; }
}

TEM void SlidingWindow<T>::resize(uint winSize, uint hopSize){
	sizeWin(winSize);
	sizeHop(hopSize);
	//mTapW = hopStart();	// for single-buffer slide mode
	mTapW = 0;				// for single-buffer rotate mode
}

TEM void SlidingWindow<T>::sizeWin(uint size){
	if(mem::resize(mBuf, sizeWin(), size)) mSizeWin = size;
}

TEM void SlidingWindow<T>::sizeHop(uint size){
	mSizeHop = scl::clip(size, sizeWin(), (uint)1);
}

TEM inline uint SlidingWindow<T>::sizeHop() const { return mSizeHop; }
TEM inline uint SlidingWindow<T>::sizeWin() const { return mSizeWin; }
TEM inline const T * SlidingWindow<T>::window(){ return mBuf; }
TEM inline const T * SlidingWindow<T>::operator()(){ return mBuf; }

TEM inline bool SlidingWindow<T>::operator()(T input){

// faster, but requires explicit call to slide() by user which is error prone
//	mBuf[mTapW] = input;
//	
//	if(++mTapW == sizeWin()){
//		mTapW = hopStart();
//		return true;
//	}
//	return false;
	
	mBuf[mTapW] = input;
	
	if(++mTapW >= sizeHop()){
		mTapW = 0;
		mem::rotateL(mBuf, sizeWin(), sizeHop());
		return true;
	}
	return false;
}

TEM inline bool SlidingWindow<T>::operator()(T * output, T input){
	mBuf[mTapW] = input;
	if(++mTapW == sizeWin()) mTapW = 0;

	if(++mHopCnt == sizeHop()){
		mem::copyAllFromRing(mBuf, sizeWin(), mTapW, output);
		mHopCnt = 0;
		return true;
	}
	return false;
}

TEM inline void SlidingWindow<T>::slide(){
	mem::move(mBuf, mBuf + sizeHop(), hopStart());
}

TEM inline uint SlidingWindow<T>::hopStart() const { return sizeWin() - sizeHop(); }



//---- DFTBase

TEM DFTBase<T>::DFTBase() : mSizeDFT(0), mNumAux(0), mBuf(0), mAux(0){ initSynced(); }

TEM DFTBase<T>::~DFTBase(){ //printf("~DFTBase\n");
	mem::free(mBuf);
	mem::free(mAux);
}

TEM inline T *		DFTBase<T>::aux(uint num){ return mAux + numBins() * num; }
TEM inline double	DFTBase<T>::binFreq() const { return spu() / (double)sizeDFT(); }
TEM inline uint	    DFTBase<T>::numBins() const { return (sizeDFT() + 2)>>1; }
TEM inline uint	    DFTBase<T>::sizeDFT() const { return mSizeDFT; }
TEM inline Sync&	DFTBase<T>::syncFreq(){ return mSyncFreq; }
TEM inline T		DFTBase<T>::normForward() const { return (T)2 / (T)sizeDFT(); }

TEM void DFTBase<T>::numAux(uint num){
	if(mem::resize(mAux, mNumAux * numBins(), num * numBins())) mNumAux = num;
}

TEM void DFTBase<T>::onResync(double r){
	mSyncFreq.ups(binFreq());
}



//---- DFT

inline float DFT::freqRes() const { return spu() / (float)sizeWin(); }
inline float DFT::overlap() const { return (float)sizeWin() / (float)sizeHop(); }
inline bool DFT::overlapping() const { return sizeHop() < sizeWin(); }
inline uint DFT::sizeHop() const { return mSizeHop; }
inline uint DFT::sizePad() const { return mSizeDFT - mSizeWin; }
inline uint DFT::sizeWin() const { return mSizeWin; }
inline Sync& DFT::syncHop(){ return mSyncHop; }
//inline float DFT::normInverse() const { return 0.5f * (float)sizeHop() / (float)sizeWin(); } /* o-a factor depends on window */
inline float DFT::normInverse() const { return 0.5f; }

inline bool DFT::operator()(float input){
	mBuf[mTapW] = input;

	if(++mTapW >= sizeHop()){
		forward(mBuf);
		mTapW = 0;
		return true;
	}
	return false;
}

inline float DFT::operator()(){
	if(++mTapR >= sizeHop()){
		inverse(0);	// this is a virtual method
		mTapR = 0;
	}
	return mBufInv[mTapR];
}

inline bool DFT::inverseOnNext(){ return mTapR == (sizeHop() - 1); }
inline void DFT::zero(){ mem::zero(mBuf, sizeDFT() + 2); }
inline void DFT::zeroEnds(){ 
	//bins0()[0]=0; bins0()[numBins()-1]=0; 
	//bins1()[0]=0; bins1()[numBins()-1]=0;
	bins()[0](0,0);
	bins()[numBins()-1](0,0);
}


//---- STFT

inline bool STFT::operator()(float input){
	if(mSlide(mBuf, input)){
		forward(mBuf);
		return true;
	}
	return false;
}

inline float STFT::unitsHop(){ return (float)DFT::sizeHop() * ups(); }
inline float * STFT::phases(){ return mPhases; }



//---- SDFT

TEM SDFT<T>::SDFT(uint sizeDFT, uint binLo, uint binHi)
	: DFTBase<T>(), mBinLo(0), mBinHi(0), mDelay(0)
{
	resize(sizeDFT, binLo, binHi);
}

TEM void SDFT<T>::resize(uint sizeDFT, uint binLo, uint binHi){
	// may be able to keep these smaller?
	mem::resize(this->mBuf, this->mSizeDFT + 2, sizeDFT + 2);
	mem::zero(this->mBuf, sizeDFT + 2);
	
	mDelay.resize(sizeDFT);
	mDelay.zero();
	
	this->mSizeDFT = sizeDFT;
	mBufI = this->mBuf + this->numBins();
	
	range(binLo, binHi);
	this->onSyncChange();
}

#define COS cos
#define SIN sin
TEM void SDFT<T>::range(uint binLo, uint binHi){
	mBinLo = binLo;
	mBinHi = binHi;
	
	double theta = M_2PI / (double)this->sizeDFT();
	mC0 = COS(theta); 
	mS0 = SIN(theta);
	theta *= (double)mBinLo;
	mCL = COS(theta);
	mSL = SIN(theta);
	mNorm = (T)2 / (T)this->sizeDFT();
}
#undef COS
#undef SIN

TEM inline void SDFT<T>::forward(T input){
	T smp = (input - mDelay(input)) * mNorm;	// ffd comb zeroes
	double cs = mCL; double sn = mSL;			// phasor at low bin
	
	for(uint k=mBinLo; k<mBinHi; ++k){
		scl::mulComplex(DFTBase<T>::mBuf[k], mBufI[k], (T)cs, (T)sn);
		scl::mulComplex(cs, sn, mC0, mS0);
		DFTBase<T>::mBuf[k] += smp;
	}
}

//TEM inline T inverse(){}



/*
TEM class Counter{
public:

	bool operator()(){
		if(++mVal == mMax){
			onWrap();
			mVal = (T)0;
			return true;
		}
		return false;
	}

	virtual void onWrap();

protected:
	T mVal, mMax;
};
*/

/*
class STFT : public Synced {
public:
	/// Constructor.
	
	/// The default complex data type is rectangular and the default window
	/// is none.
	/// @param[in]	windowSize	Number of samples to window.
	/// @param[in]	hopSize		Number of samples between successive windows.
	/// @param[in]	padSize		Number of zeros to append to window.
	/// @param[in]	winType		Type of forward transform window.
	/// @param[in]	format		Form of complex spectral data.
	/// @param[in]	createAux	Whether to create an auxillary spectral buffer (see createAux()).
	STFT(ULONG windowSize, ULONG hopSize, ULONG padSize=0,
		WinType::type winType = WinType::Rectangle,
		ComplexType::type format = ComplexType::Rect,
		bool createAux = false);
	~STFT();

	// Setters
	
	/// Whether or not to apply tapering window during resynthesis.
	
	/// A tapering window is generally applied to the resynthesized samples
	/// to avoid phase discontinuities in the output signal.  The window used
	/// is a Bartlett (triangle) window.  Tapering is turned on by default.
	void doesTapering(bool val);
	void format(ComplexType::type format);						///< Sets format of complex spectrum data.
	void hopSize(ULONG newHopSize);								///< Sets hop size.
	void size(ULONG windowSize, ULONG hopSize, ULONG padSize);	///< Sets size parameters of transform.
	void windowType(WinType::type type);						///< Sets window type.

	// Getters
	float binFreq();	///< Returns fundamental bin frequency.
	ULONG dftSize();	///< Returns size of DFT.  This equals the window plus zero-padding size.
	float freqRes();	///< Returns frequency resolution of analysis.
	float hopDur();		///< Returns hop size in seconds.
	ULONG hopSize();	///< Returns size of analysis hop.
	ULONG numBins();	///< Returns number of frequency samples (bins) per frame.
	float overlap();	///< Returns overlap factor (windowSize / hopSize).
	ULONG padSize();	///< Returns size of zero-padding.
	ULONG windowSize();	///< Returns size of window.
	
	/// Returns pointer to a row of spectral data.
	
	/// Row 0 is either the reals or mags and row 1 is either imags, phases, or
	/// instantaneous frequencies.
	float * spct(ULONG row);
	
	/// Returns pointer to phase accumulation buffer used in mag-freq form.
	float * acc ();
	
	float * aux ();			///< Returns pointer to auxiliary buffer (see 'createAux()').
	float * fftBuffer();	///< Returns pointer to FFT buffer.
	
	/// Performs forward transform of sliding window DFT.
	void forward(const float * hopBuffer);
	void forward();
	
	///	Performs inverse transform using overlap-add.
	
	///	The resynthesized samples are copied into 'hopBuffer'.  The destination
	/// array must have room for at least hopSize number of elements.
	void inverse(float * hopBuffer);

	///	Performs inverse transform using overlap-add.
	
	/// Returns pointer to output hop samples.
	///
	float * inverse();
	
	/// Allocates space for an auxiliary buffer of size equal to 2 * numBins.
	
	/// Useful for analyses/transformations requiring a temp buffer.
	/// Access with aux().
	void createAux();
	
	void resetPhaseAccum();

	void spctToRect();		// convert spectrum to rectangle format
	void spctToPolar();		// convert spectrum to polar format

	virtual void print();

protected:

	float * mWindow;		// analysis window table
	Spectrum * mSpct;		// deinterleaved complex spectrum

	ULONG mWindowSize;		// number of samples in analysis window
	ULONG mPadSize;			// number of samples to zero pad for spectral interpolation
	ULONG mHopSize;			// number of samples between analyses
	WinType::type mWindowType;		// type of analysis window used
	ComplexType::type mSpctFormatT2F;	// format of incoming spectrum (time to freq)
	bool mDoesTapering;		// whether to do tapering during resynthesis

		// Magnitude normalization between time (t) and frequency (f) domains
	float t2fNormFactor();
	float f2tNormFactor();

		// Instantaneous frequency variables
	float factorWrap();
	float factorUnwrap();
	float fundRadians();

private:
	float * mFFTBuffer;	// buffer used by FFT and IFFT
	float * timeIn;		// ring buffer for time-domain input
	float * timeOut;	// ring buffer for time-domain output (overlap added)
	FFTInfo fft;		// info for FFT
	FFTInfo ifft;		// info for IFFT
	float * phase;		// phases of current spectral frame
	float * phaseAccum;	// phase accumulation buffer for mag/freq to polar conversion
	float * mAux;		// auxiliary buffer (created on demand)
	
	ULONG mTapIn, mTapOut;		// input/output ring buffer taps.

	// These methods can possibly be factored out
	
		// Copy samples from external buffer into local ring buffer
		void copyHopSamples(const float * hopBuffer);
	
	void wrapTapIn();
	void wrapTapOut();
};

inline ULONG STFT::numBins(){ return (dftSize() + 2)>>1; }
inline ULONG STFT::dftSize(){ return mWindowSize + mPadSize; }
inline ULONG STFT::hopSize(){ return mHopSize; }
inline ULONG STFT::padSize(){ return mPadSize; }
inline ULONG STFT::windowSize(){ return mWindowSize; }
inline float STFT::overlap(){ return (float)windowSize() / (float)mHopSize; }
inline float STFT::binFreq(){ return spu() / (float)dftSize(); }
inline float STFT::freqRes(){ return spu() / (float)windowSize(); }
inline float STFT::hopDur(){ return (float)mHopSize * ups(); }
inline float * STFT::spct(ULONG row){ return mSpct->row(row); }
inline float * STFT::acc (){ return phaseAccum; }
inline float * STFT::aux (){ return mAux; }
inline float * STFT::fftBuffer(){ return mFFTBuffer; }
inline float STFT::t2fNormFactor(){ return 2.f / (float)dftSize(); }
inline float STFT::f2tNormFactor(){ return 0.5f * (float)mHopSize / (float)mWindowSize; }
inline float STFT::factorWrap(){ return spu() / ((float)mHopSize * M_2PI); }
inline float STFT::factorUnwrap(){ return M_2PI * (float)mHopSize * ups(); }
inline float STFT::fundRadians(){ return M_2PI * (float)mHopSize / (float)dftSize(); }
inline void STFT::inverse(float * hopBuffer){ mem::copy(hopBuffer, inverse(), mHopSize); }

inline void STFT::wrapTapIn(){ if(mTapIn >= mWindowSize) mTapIn %= mWindowSize; }
inline void STFT::wrapTapOut(){ if(mTapOut >= dftSize()) mTapOut %= dftSize(); }
//*/


#undef TEM

} // gam::

#endif

