#include <stdio.h>
#include "Gamma/Recorder.h"

namespace gam{

Recorder::Recorder(int channels, int frames){
	resize(frames, channels);
}

//int write(const float * buf, int numFrames){
//
//	if(numFrames > frames()) numFrames = frames();
//
//	int Nr = mRing.size();
//	int Nw = numFrames * channels();
//	
//	if((mIW+Nw) > Nr){
//		memcpy(&mRing[mIW], buf, Nr-mIW);
//		memcpy(&mRing[  0], buf + Nr-mIW, mIW+Nw - Nr);
//		mIW += Nw;
//		mIW -= mRing.size();
//	}
//	else{
//		memcpy(&mRing[mIW], buf, Nw);
//		mIW += Nw;
//	}
//	return numFrames;
//}

int Recorder::read(float *& buf){
	if(mIR == mIW) return 0;

	int N = size();
	int iw = mIW;
	int behind = 0;

	/*
	01234567
	r---w		case 1
	 r---w		case 1
	   r---w	case 1
	w   r---	case 2
	-w   r--	case 2
	*/
	
	if(mIR < iw){	// case 1: contiguous copy
		behind = iw - mIR;
		for(int i=0; i<behind; ++i){
			mRead[i] = mRing[mIR+i];
		}
//		printf("w: %4d, r: %4d, d: %4d\n", (int)mIW, (int)mIR, behind);			
	}
	else{			// case 2: two fragment copy
		int d1 = N - mIR;
		for(int i=0; i<d1; ++i){
			mRead[i] = mRing[mIR+i];
		}
		
		int d2 = mIW;
		for(int i=d1; i<d1+d2; ++i){
			mRead[i] = mRing[i-d1];
		}
//		printf("w: %4d, r: %4d, d: %4d + %4d\n", mIW, mIR, d1, d2);
		behind = d1+d2;
	}
	mIR += behind;
	if(mIR >= N) mIR -= N;
	
	buf = &mRead[0];
	return behind/channels();
}

void Recorder::resize(int frames, int chans){
	mChans = chans;
	mRing.resize(frames*mChans);
	mRead.resize(frames*mChans);
}

} // gam::
