/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include "Gamma/DFT.h"

#include "Gamma/arr.h"
#include "Gamma/gen.h"
#include "Gamma/mem.h"
#include "Gamma/scl.h"
#include "Gamma/Access.h"
#include "Gamma/Sync.h"


#define RECT_POLAR\
	if(mPrecise){\
		for(uint32_t i=1; i<numBins()-1; ++i){\
			Complex<float> &c = bins()[i];\
			c(c.mag(), c.phase());\
		}\
	}\
	else{\
		for(uint32_t i=1; i<numBins()-1; ++i){\
			Complex<float> &c = bins()[i];\
			float m = c.norm2();\
			float p = scl::atan2Fast(c.i,c.r);\
			c(scl::sqrt<1>(m), p);\
		}\
	}

#define POLAR_RECT\
	if(mPrecise){\
		for(uint32_t i=1; i<numBins()-1; ++i){\
			Complex<float> &c = bins()[i];\
			c(c[0]*cos(c[1]), c[0]*sin(c[1]));\
		}\
	}\
	else{\
		for(uint32_t i=1; i<numBins()-1; ++i){\
			Complex<float> &c = bins()[i];\
			float p = scl::wrapPhase(c[1]);\
			c(c[0]*scl::cosT8(p), c[0]*scl::sinT9(p));\
		}\
	}\


namespace gam{

DFT::DFT(uint32_t winSize, uint32_t padSize, Bin::t binT, uint32_t numAuxA) :
	mSizeWin(0), mSizeHop(0),
	mFFT(0),
	mPadOA(0), mTapW(0), mTapR(0), mPrecise(false)
{
	//printf("DFT::DFT\n");
	resize(winSize, padSize);
	binType(binT);	
	numAux(numAuxA);
	
	arr::conversionInit();
}

DFT::~DFT(){
	//printf("~DFT\n");
	if(mBufInv != mBuf) mem::free(mBufInv);
	mem::free(mPadOA);
}

void DFT::resize(uint32_t newWindowSize, uint32_t newPadSize){ //printf("DFT::resize()\n");

	if(0 == newWindowSize && 0 == newPadSize) return;

	uint32_t oldDFTSize = sizeDFT();
	uint32_t newDFTSize = newWindowSize + newPadSize;
	uint32_t oldFrqSize = oldDFTSize+2;		// 2 extra for DC/Nyquist imaginary
	uint32_t newFrqSize = newDFTSize+2;		// "

	if(mem::resize(mBuf, oldFrqSize, newFrqSize)){
		mBufInv = mBuf;
		if(mAux) mem::resize(mAux, oldFrqSize, newFrqSize);
		mFFT.resize(newDFTSize);
		mem::deepZero(mBuf, newFrqSize);
	}
	
	mem::resize(mPadOA, sizePad(), newPadSize);
	mem::deepZero(mPadOA, newPadSize);
	
	mSizeDFT = newDFTSize;
	mSizeWin = newWindowSize;
	mSizeHop = mSizeWin;

	onResync(1);
}

void DFT::onResync(double r){
	DFTBase<float>::onResync(r);
	mSyncHop.ups((double)sizeHop() * ups());
	//printf("[%p] hop: %d, ups: %f\n", this, sizeHop(), spu());
}

void DFT::binType(Bin::t type){ mSpctFormat = type; }
DFT& DFT::precise(bool w){ mPrecise=w; return *this; }

void DFT::forward(const float * window){ //printf("DFT::forward(const float *)\n");

	// this attempts to avoid a memory move, but doesn't work when window == mBuf
//	if(window != mBuf) mem::deepCopy(mBuf+1, window, sizeWin());
//	mem::deepZero(mBuf+1 + sizeWin(), sizePad());	// zero pad
//	// TODO: zero-phase window (rotate buffer 180)
//
//	mFFT.forward(mBuf+1, true);
//	
//	// re-arrange DC and Nyquist bins
//	mBuf[0] = mBuf[1];
//	mBuf[1] = 0.f;
//	mBuf[numBins()*2-1] = 0.f;

	// old code which did a mem move, rather offsetting input buffer pointer...
	if(window != mBuf) mem::deepCopy(mBuf, window, sizeWin());
	mem::deepZero(mBuf + sizeWin(), sizePad());	// zero pad
	//... do zero-phase window (rotate buffer 180)

	mFFT.forward(mBuf, true);
	
	// re-arrange DC and Nyquist bins
	mem::deepMove(mBuf+2, mBuf+1, sizeDFT()-1);
	mBuf[1] = 0.f;
	mBuf[numBins()*2-1] = 0.f;
	
	switch(mSpctFormat){
	case Bin::Polar:
	case Bin::MagFreq:
		RECT_POLAR
		//arr::mul(bins0(), gen::val(normForward()), nbins);
		break;
	
	default:;	// rectangular
		//arr::mul(bins0(), gen::val(normForward()), nbins<<1);
	}
}

void DFT::inverse(float * output){
	//printf("DFT::inverse(float *)\n");
	
	switch(mSpctFormat){
	case Bin::Polar:
	case Bin::MagFreq:
		//arr::mul(bins0(), gen::val(normInverse()), nbins);
		POLAR_RECT
		break;
	
	default:;
		//arr::mul(bins0(), gen::val(normInverse()), nbins<<1);
	}	

	// arrange/scale bins for inverse xfm
	// TODO: can we avoid this move by pointer offsetting?
	mem::deepMove(mBuf+1, mBuf+2, sizeDFT()-1);
	slice(mBuf+1, sizeDFT()-2) *= 0.5f;

	mFFT.inverse(mBuf);
	
	// o.a. inverse window with prev spill
	if(sizePad() > 0){
		arr::add(mBuf, mPadOA, scl::min(sizePad(), sizeWin()));	// add prev spill
	
		if(sizePad() <= sizeWin()){	// no spill overlap
			mem::deepCopy(mPadOA, mBuf + sizeWin(), sizePad());	// save current spill
		}
		else{						// spill overlaps
			// add and save current spill to previous
			arr::add(mPadOA, mBuf + sizeWin(), mPadOA + sizeWin(), sizePad() - sizeWin());
			mem::deepCopy(mPadOA + sizePad() - sizeWin(), mBuf + sizePad(), sizeWin());
		}
	}

	if(output) mem::deepCopy(output, mBuf, sizeWin());
}

void DFT::spctToRect(){
	switch(mSpctFormat){
	case Bin::Polar: POLAR_RECT break;
	default:;
	}
	mSpctFormat = Bin::Rect;
}

void DFT::spctToPolar(){	
	switch(mSpctFormat){
	case Bin::Rect:	RECT_POLAR break;
	default:;
	}
	mSpctFormat = Bin::Polar;
}

void DFT::print(FILE * f, const char * a){
	fprintf(f, "DFT, Win, Hop: %d, %d, %d samples\n", (int)sizeDFT(), (int)sizeWin(), (int)sizeHop());
	fprintf(f, "# bins:        %d\n", (int)numBins());
	fprintf(f, "Freq res:      %f units/sample\n", freqRes());
	fprintf(f, "Bin freq:      %f units\n", binFreq());
	fprintf(f, "Data format:   %s\n", Bin::string(mSpctFormat));
	fprintf(f, "Precise:       %s\n", mPrecise ? "true" : "false");	
	fprintf(f, "Aux buffers:   %d\n", (int)mNumAux);
	//printf("buf, pad, inv, aux: %p %p %p %p", mBuf, mPadOA, mBufInv, mAux);
	fprintf(f, "%s", a);
}


//---- STFT

STFT::STFT(uint32_t winSize, uint32_t hopSize, uint32_t padSize, WinType::type windowType, Bin::t binT, uint32_t numAuxA) :
	DFT(0, 0, binT, numAuxA),
	mSlide(winSize, hopSize), mFwdWin(0), mPhases(0),
	mWinType(windowType),
	mWindowInverse(true), mRotateForward(false)
{
	mBufInv = 0;
	resize(winSize, padSize);
	sizeHop(hopSize);
}

STFT::~STFT(){ //printf("~STFT\n");
	mem::free(mFwdWin);
}


void STFT::computeInvWinMul(){

	// compute sum of overlapping elements of forward window
	if(overlapping()){
		float sum = 0.f;
		for(uint32_t i=0; i<sizeWin(); i+=sizeHop()){
			sum += mFwdWin[i] * (mWindowInverse ? scl::bartlett(2*i/(float)sizeWin() - 1.f): 1.f);
		}
		mInvWinMul = 1/sum;
	}
	
	// if no overlap, then do not scale output
	else{
		mInvWinMul = 1;
	}

	//printf("mInvWinMul: %f\n", mInvWinMul);
}


void STFT::winType(WinType::type type){
	tbl::window(mFwdWin, sizeWin(), type);
	mFwdWinMul = 1.f / arr::mean(mFwdWin, sizeWin());	// compute mul factor for normalization
	
	// scale forward window?
	//slice(mFwdWin, sizeWin()) *= mFwdWinMul;
	
	computeInvWinMul();
	mWinType = type;
}


void STFT::resize(uint32_t winSize, uint32_t padSize){
	float * tmp = mBufInv;
	uint32_t oldWinSize = sizeWin();
	uint32_t oldNumBins = numBins();
	
	// resize DFT buffers
	DFT::resize(winSize, padSize);
	
	mBufInv = tmp;	// DFT::resize changes this, so change it back
	
	// resize STFT-specific buffers
	mSlide.sizeWin(winSize);
	mem::resize(mFwdWin, oldWinSize, winSize);
	mem::resize(mBufInv, oldWinSize, winSize);
	mem::resize(mPhases, oldNumBins, numBins());
	
	mem::deepZero(mPhases, numBins());
	mem::deepZero(mBufInv, winSize);
	
	// re-compute fwd window
	winType(mWinType);
}


void STFT::sizeHop(uint32_t size){
	mSlide.sizeHop(size);
	mSizeHop = mSlide.sizeHop();
	computeInvWinMul();
	onResync(1);
}


// input is sizeWin
void STFT::forward(float * input){ //printf("STFT::forward(float *)\n");

	arr::mul(input, mFwdWin, sizeWin());	// apply forward window
	if(mRotateForward) mem::rotateHalf(input, sizeWin()); // do zero-phase windowing rotation?
	DFT::forward(input);					// do forward transform
	
	// compute frequency estimates?
	if(Bin::MagFreq == mSpctFormat){
		
		// This will effectively subtract the expected phase difference from the computed.
		// This extra step seems to give more precise frequency estimates.
		slice(mPhases, numBins()) += float(M_2PI * sizeHop()) / sizeDFT();
		
		// compute relative frequencies
		//arr::phaseToFreq(phs, mPhases, numBins(), unitsHop());
		float factor = 1.f / (M_2PI * unitsHop());
		for(uint32_t i=1; i<numBins()-1; ++i){
			float dp = scl::wrapPhase(bins()[i][1] - mPhases[i]);	// wrap phase into [-pi, pi)
			mPhases[i] = bins()[i][1];							// prev phase = curr phase
			bins()[i][1] = dp*factor;
		}
		
		// compute absolute frequencies by adding respective bin center frequency
		slice(mBuf, numBins(), 2) += gen::RAdd<float>(binFreq());
	}
}


void STFT::inverse(float * dst){
	//printf("STFT::inverse(float *)\n");
	if(Bin::MagFreq == mSpctFormat){
		//mem::copy(bins1(), mPhases, numBins()); // not correct, need to unwrap frequencies
		for(uint32_t i=1; i<numBins()-1; ++i) bins()[i] = mPhases[i];
	}	
	
	DFT::inverse(0);	// result goes into mBuf
	
	// undo zero-phase windowing rotation?
	if(mRotateForward) mem::rotateHalf(mBuf, sizeWin());
	
	// apply secondary window to smooth ends?
	if(mWindowInverse){
		arr::mulBartlett(mBuf, sizeWin());
	}

	if(overlapping()){	//inverse windows overlap?
	
		// scale inverse so overlap-add is normalized
		//arr::mul(mBuf, gen::val(mInvWinMul), sizeWin());
		slice(mBuf, sizeWin()) *= mInvWinMul;
	
		// shift old output left while adding new output
		arr::add(mBufInv, mBuf, mBufInv + sizeHop(), sizeWin() - sizeHop());
	}

	// copy remaining non-overlapped portion of new output
	uint32_t sizeOverlap = sizeWin() - sizeHop();
	mem::deepCopy(mBufInv + sizeOverlap, mBuf + sizeOverlap, sizeHop());

	// copy output if external buffer provided
	if(dst) mem::deepCopy(dst, mBufInv, sizeWin());	
}


void STFT::print(FILE * f, const char * a){
	DFT::print(f, "");
	fprintf(f, "Window type:   %s\n", WinType::string(mWinType));
	fprintf(f, "Inv. window:   %s\n", mWindowInverse ? "true" : "false");
	fprintf(f, "%s", a);
}

} // gam::

#undef RECT_POLAR
#undef POLAR_RECT

