/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include "Gamma/DFT.h"
#include "Gamma/arr.h"
#include "Gamma/mem.h"
#include "Gamma/scl.h"
#include "Gamma/Sync.h"

#define CART_TO_POL\
	if(mPrecise){\
		for(uint32_t i=1; i<numBins()-1; ++i){\
			Complex<float> &c = bins()[i];\
			c(c.mag(), c.phase());\
		}\
	}\
	else{\
		for(uint32_t i=1; i<numBins()-1; ++i){\
			Complex<float> &c = bins()[i];\
			float m = c.normSqr();\
			float p = scl::atan2Fast(c.i,c.r);\
			c(scl::sqrt<1>(m), p);\
		}\
	}

#define POL_TO_CART\
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

DFT::DFT(uint32_t winSize, uint32_t padSize, SpectralType specT, uint32_t numAuxA)
:	mSizeWin(0), mSizeHop(0),
	mFFT(0),
	mPadOA(0), mTapW(0), mTapR(0), mPrecise(false)
{
	//printf("DFT::DFT\n");
	resize(winSize, padSize);
	spectrumType(specT);	
	numAux(numAuxA);
}

DFT::~DFT(){
	//printf("~DFT\n");
	if(mBufInv != bufPos()) mem::free(mBufInv);
	mem::free(mPadOA);
}

void DFT::resize(uint32_t newWinSize, uint32_t newPadSize){ //printf("DFT::resize()\n");

	if(0 == newWinSize && 0 == newPadSize) return;

	uint32_t oldDFTSize = sizeDFT();
	uint32_t newDFTSize = newWinSize + newPadSize;
	uint32_t oldFrqSize = oldDFTSize+2;		// 2 extra for DC/Nyquist imaginary
	uint32_t newFrqSize = newDFTSize+2;		// "

	if(mem::resize(mBuf, oldFrqSize, newFrqSize)){
		mBufInv = bufPos();
		if(mAux) mem::resize(mAux, oldFrqSize, newFrqSize);
		mFFT.resize(newDFTSize);
		mem::deepZero(mBuf, newFrqSize);
	}
	
	mem::resize(mPadOA, sizePad(), newPadSize);
	mem::deepZero(mPadOA, newPadSize);
	
	mSizeDFT = newDFTSize;
	mSizeWin = newWinSize;
	mSizeHop = mSizeWin;

	onResync(1);
}

void DFT::onResync(double r){
	DFTBase<float>::onResync(r);
	mSyncHop.ups((double)sizeHop() * ups());
	//printf("[%p] hop: %d, ups: %f\n", this, sizeHop(), spu());
}

void DFT::forward(const float * src){ //printf("DFT::forward(const float *)\n");

	if(src != bufPos()) mem::deepCopy(bufPos(), src, sizeWin());
	mem::deepZero(bufPos() + sizeWin(), sizePad());	// zero pad

	mFFT.forward(mBuf, true, true); // complex buffer and normalize

//	// do a mem move rather than offsetting input buffer pointer...
//	if(src != mBuf) mem::deepCopy(mBuf, src, sizeWin());
//	mem::deepZero(mBuf + sizeWin(), sizePad());	// zero pad
//	//... do zero-phase window (rotate buffer 180)
//
//	mFFT.forward(mBuf);
//	
//	// re-arrange DC and Nyquist bins
//	mem::deepMove(mBuf+2, mBuf+1, sizeDFT()-1);
//	mBuf[1] = 0.f;
//	mBuf[numBins()*2-1] = 0.f;

	switch(mSpctFormat){
	case COMPLEX: break;
	case MAG_PHASE:
	case MAG_FREQ:
		CART_TO_POL
		break;
	default:;
	}
}

void DFT::inverse(float * dst){
	//printf("DFT::inverse(float *)\n");
	
	switch(mSpctFormat){
	case COMPLEX: break;
	case MAG_PHASE:
	case MAG_FREQ:
		POL_TO_CART
		break;
	}	

	mFFT.inverse(mBuf, true);

	// overlap-add inverse window with prev spill
	if(sizePad() > 0){
		arr::add(bufPos(), mPadOA, scl::min(sizePad(), sizeWin()));	// add prev spill
	
		if(sizePad() <= sizeWin()){	// no spill overlap
			mem::deepCopy(mPadOA, bufPos() + sizeWin(), sizePad());	// save current spill
		}
		else{						// spill overlaps
			// add and save current spill to previous
			arr::add(mPadOA, bufPos() + sizeWin(), mPadOA + sizeWin(), sizePad() - sizeWin());
			mem::deepCopy(mPadOA + sizePad() - sizeWin(), bufPos() + sizePad(), sizeWin());
		}
	}

	if(dst) mem::deepCopy(dst, bufPos(), sizeWin());

//	// arrange/scale bins for inverse xfm
//	// TODO: can we avoid this move by pointer offsetting?
//	mem::deepMove(mBuf+1, mBuf+2, sizeDFT()-1);
//	//slice(mBuf+1, sizeDFT()-2) *= 0.5f;
//
//	mFFT.inverse(mBuf);
//	
//	// o.a. inverse window with prev spill
//	if(sizePad() > 0){
//		arr::add(mBuf, mPadOA, scl::min(sizePad(), sizeWin()));	// add prev spill
//	
//		if(sizePad() <= sizeWin()){	// no spill overlap
//			mem::deepCopy(mPadOA, mBuf + sizeWin(), sizePad());	// save current spill
//		}
//		else{						// spill overlaps
//			// add and save current spill to previous
//			arr::add(mPadOA, mBuf + sizeWin(), mPadOA + sizeWin(), sizePad() - sizeWin());
//			mem::deepCopy(mPadOA + sizePad() - sizeWin(), mBuf + sizePad(), sizeWin());
//		}
//	}
//
//	if(dst) mem::deepCopy(dst, mBuf, sizeWin());
}

void DFT::spctToRect(){
	switch(mSpctFormat){
	case MAG_PHASE: POL_TO_CART break;
	default:;
	}
	mSpctFormat = COMPLEX;
}

void DFT::spctToPolar(){
	switch(mSpctFormat){
	case COMPLEX:	CART_TO_POL break;
	default:;
	}
	mSpctFormat = MAG_PHASE;
}


static const char * toString(SpectralType v){
	switch(v){
	case COMPLEX:	return "Complex";
	case MAG_PHASE:	return "Magnitude/Phase";
	case MAG_FREQ:	return "Magnitude/Frequency";
	default:		return "Unknown";
	}	
}

void DFT::print(FILE * f, const char * a){
	fprintf(f, "DFT, Win, Hop: %d, %d, %d samples\n", (int)sizeDFT(), (int)sizeWin(), (int)sizeHop());
	fprintf(f, "# bins:        %d\n", (int)numBins());
	fprintf(f, "Freq res:      %f units/sample\n", freqRes());
	fprintf(f, "Bin freq:      %f units\n", binFreq());
	fprintf(f, "Data format:   %s\n", toString(mSpctFormat));
	fprintf(f, "Precise:       %s\n", mPrecise ? "true" : "false");	
	fprintf(f, "Aux buffers:   %d\n", (int)mNumAux);
	//printf("buf, pad, inv, aux: %p %p %p %p", mBuf, mPadOA, mBufInv, mAux);
	fprintf(f, "%s", a);
}


//---- STFT

STFT::STFT(uint32_t winSize, uint32_t hopSize, uint32_t padSize, WindowType winType, SpectralType specType, uint32_t numAuxA)
:	DFT(0, 0, specType, numAuxA),
	mSlide(winSize, hopSize), mFwdWin(0), mPhases(0),
	mWinType(winType),
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
			float invWin =  mWindowInverse ? scl::bartlett(2*i/(float)sizeWin() - 1.f) : 1.f;
			sum += mFwdWin[i] * invWin;
		}
		mInvWinMul = 1/sum;
	}
	
	// if no overlap, then do not scale output
	else{
		mInvWinMul = 1;
	}
	//printf("mFwdWinMul: %f\n", mFwdWinMul);
	//printf("mInvWinMul: %f\n", mInvWinMul);
}


STFT& STFT::windowType(WindowType v){
	mWinType = v;
	tbl::window(mFwdWin, sizeWin(), mWinType);
	
	// compute forward normalization factor
	mFwdWinMul = 1.f / arr::mean(mFwdWin, sizeWin());
	
	// scale forward window?
	//slice(mFwdWin, sizeWin()) *= mFwdWinMul;
	
	computeInvWinMul();
	return *this;
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
	windowType(mWinType);
}


STFT& STFT::sizeHop(uint32_t size){
	mSlide.sizeHop(size);
	mSizeHop = mSlide.sizeHop();
	computeInvWinMul();
	onResync(1);
	return *this;
}


// input is sizeWin
void STFT::forward(const float * src){ //printf("STFT::forward(float *)\n");

	if(src != bufPos()) mem::deepCopy(bufPos(), src, sizeWin());

	// apply forward window
	arr::mul(bufPos(), mFwdWin, sizeWin());
	
	// do zero-phase windowing rotation?
	if(mRotateForward) mem::rotateHalf(bufPos(), sizeWin());

	DFT::forward(bufPos());
	
	// compute frequency estimates?
	if(MAG_FREQ == mSpctFormat){
		
//		// This will effectively subtract the expected phase difference from the computed.
//		// This extra step seems to give more precise frequency estimates.
//		{
//			float expdp = float(M_2PI * sizeHop()) / sizeDFT();
//			for(uint32_t i=0; i<numBins(); ++i){
//				mPhases[i] += expdp;
//			}
//		}
		
		// compute frequency estimates
		float factor = 1.f / (M_2PI * unitsHop()); // converts phase difference from radians to Hz
		
		// expected phase difference of fundamental
		float expdp1 = float(sizeHop())/sizeWin() * M_2PI;

		for(uint32_t i=1; i<numBins()-1; ++i){
			float ph = bin(i)[1];						// current phase
			float dp = ph - mPhases[i] - i*expdp1;		// compute phase difference
			dp = scl::wrapPhase(dp);					// wrap back into [-pi, pi)
			mPhases[i] = ph;							// save current phase
			bin(i)[1] = dp*factor;						// convert phase diff to freq
		}
		
		// compute absolute frequencies by adding respective bin center frequency
		for(uint32_t i=0; i<numBins(); ++i){
			bins()[i][1] += binFreq()*i;
		}
	}
}


void STFT::inverse(float * dst){
	//printf("STFT::inverse(float *)\n");
	if(MAG_FREQ == mSpctFormat){
		// TODO:
		for(uint32_t i=1; i<numBins()-1; ++i) bin(i)[1] = mPhases[i];
	}	
	
	DFT::inverse(0);	// result goes into bufPos()
	
	// undo zero-phase windowing rotation?
	if(mRotateForward) mem::rotateHalf(bufPos(), sizeWin());

	// apply secondary window to smooth ends?
	if(mWindowInverse){
		arr::mulBartlett(bufPos(), sizeWin());
	}

	if(overlapping()){	//inverse windows overlap?
	
		// scale inverse so overlap-add is normalized
		//slice(mBuf, sizeWin()) *= mInvWinMul;
		for(unsigned i=0; i<sizeWin(); ++i){
			bufPos()[i] *= mInvWinMul;
		}
	
		// shift old output left while adding new output
		arr::add(mBufInv, bufPos(), mBufInv + sizeHop(), sizeWin() - sizeHop());
	}

	// copy remaining non-overlapped portion of new output
	uint32_t sizeOverlap = sizeWin() - sizeHop();
	mem::deepCopy(mBufInv + sizeOverlap, bufPos() + sizeOverlap, sizeHop());

	// copy output if external buffer provided
	if(dst) mem::deepCopy(dst, mBufInv, sizeWin());	
}


void STFT::print(FILE * f, const char * a){
	DFT::print(f, "");
	fprintf(f, "Window type:   %s\n", toString(mWinType));
	fprintf(f, "Inv. window:   %s\n", mWindowInverse ? "true" : "false");
	fprintf(f, "%s", a);
}

} // gam::

#undef CART_TO_POL
#undef POL_TO_CART

