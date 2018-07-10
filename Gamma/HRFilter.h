#ifndef GAM_INC_HR_FILTER_H
#define GAM_INC_HR_FILTER_H

#include "Gamma/Filter.h"
#include "Gamma/Spatial.h"
#include "Gamma/Types.h"

namespace gam{

template <class Vec3>
float headShadow(const Vec3& earDir, const Vec3& srcDir){
	return earDir.dot(srcDir)*0.25+0.75;
}

/// Head-related filter

/// This is a parametric HRTF based largely on:
/// Iida, K. et al. (2007). Median plane localization using a parametric model of
/// the head-related transfer function based on spectral cues. 
/// Applied Acoustics, 68, 835-850.
/// Iida, K. and Ishii, Y. (2018). Effects of adding a spectral peak generated 
/// by the second pinna resonance to a parametric model of head-related transfer
/// functions on upper median plane sound localization.
/// Applied Acoustics, 129, 239-247.
///
/// Default distances are assumed to be in meters.
class HRFilter{
public:

	HRFilter(){
		mBackShelf.set(4000, 0.707, gam::HIGH_SHELF);
		mPinnaPeak1.set(5000, 1.000, gam::PEAKING); mPinnaPeak1.level(2); // concha resonance
		mPinnaPeak2.set(8000, 4.000, gam::PEAKING);
	}

	Dist<3>& dist(){ return mDist; }

	/// Set ear distance (measured from center of head)
	HRFilter& earDist(float v){ mEarDist = v; return *this; }

	/// Set position of source

	/// @param[in] sourcePos	Position of sound source
	/// @param[in] headPose		Head pose represented as a 4x4 matrix.
	///							The matrix must be right-handed and indexable 
	///							with column-major layout.
	template <class Vec3, class Mat4>
	HRFilter& pos(const Vec3& sourcePos, const Mat4& headPose){
		auto headToSrc = sourcePos - Vec3(headPose[12], headPose[13], headPose[14]);
		Vec3 headR(headPose[ 0], headPose[ 1], headPose[ 2]);
		Vec3 headU(headPose[ 4], headPose[ 5], headPose[ 6]);
		Vec3 headB(headPose[ 8], headPose[ 9], headPose[10]); // pose stores back vector since right-handed
		auto earR = headR*mEarDist;
		const auto& S = headToSrc;

		mDist.near(mEarDist);
		mDist.dist(0, (S + earR).mag()); // left channel
		mDist.dist(1, (S - earR).mag()); // right channel
		mDist.dist(2, mRoomSize); // diameter of room

		// TODO: apply for each ear?
		auto Su = S; Su.normalize();
		//float hsLvl = headB.dot(Su)*-0.15 + 0.85; // [0.7, 1] -> [back,front]
		float hsLvl = headB.dot(Su)*-0.25 + 0.75; // [0.5, 1] -> [back,front]
		mBackShelf.level(hsLvl);

		// Pinna filtering
		auto elVec = Su - (headR*Su.dot(headR)); // projection of source dir onto up-back plane
		elVec.normalize();
		float srcUp = headU.dot(elVec); // how much source is up, [-1,1], sine of elevation angle
		float pinnaFrq1 = srcUp*2000. + 8000.; // [6k, 10k] -> [below, above]
		// Notch weak above, and strong in front (-10 dB / 0.3)
		// Following a half-period raised cosine
		float pinnaLvl1 = srcUp*srcUp*0.7+0.3; // [0.3, 1]
		mPinnaNotch1.set(pinnaFrq1, 24., gam::PEAKING);
		mPinnaNotch1.level(pinnaLvl1);

		float srcB = headB.dot(elVec);
		float pinnaFrq2 = srcB*1000. + 10000.; // [9k, 11k] -> [front, above, behind]
		// Notch weak above, and strong in front (-10 dB / 0.3)
		// Following a half-period raised cosine
		float pinnaLvl2 = srcUp*srcUp*0.8+0.2; // [0.2, 1]
		mPinnaNotch2.set(pinnaFrq2, 48., gam::PEAKING);
		mPinnaNotch2.level(pinnaLvl2);

		float above = srcUp>0. ? srcUp : 0.;
		//mPinnaPeak2.level(srcUp*srcUp*1.+1.);
		mPinnaPeak2.level(above + 1.);

		mShadows[0] = headShadow(-headR, (S + earR).normalize());
		mShadows[1] = headShadow( headR, (S - earR).normalize());

		/* Torso filtering (only useful if body orientation available!)
		mTorsoAmt = al::abs(srcUp)*0.1;
		mTorsoDelay = 1./1000. * mTorsoAmt/0.3;
		//*/

		//printf("front/back high-shelf level: %f\n", mBackShelf.level());
		//printf("pinna notch 1: %.2f Hz @ %.2f\n", mPinnaNotch1.freq(), mPinnaNotch1.level());
		//printf("pinna notch 2: %.2f Hz @ %.2f\n", mPinnaNotch2.freq(), mPinnaNotch2.level());
		//printf("head shadow on left ear: %f\n", headShadow(-headR, (S + earR).normalize()));
		return *this;

		//return posOld(&headPose[0], &headToSrc[0]);
	}

	/// Return spatialized sample as (left, right, room)
	float3 operator()(float src){
		// Note: the correct way would be to filter at each ear
		src = mBackShelf(src);
		src = mPinnaNotch1(src);
		src = mPinnaNotch2(src);
		src = mPinnaPeak1(src);
		src = mPinnaPeak2(src);
		//src += mDist.delayLine().read(mTorsoDelay) * mTorsoAmt;
		auto res = mDist(src);
		for(int i=0;i<2;++i) res[i] *= mShadows[i];
		return res;
	}

/*
	// @param[in] headPose		pointer to elements of column-major 4x4 matrix
	//								using a right-handed coordinate system
	// @param[in] headToSrcVec	vector pointing from head to source position
	HRFilter& posOld(const float * headPose, const float * headToSrcVec){
		return pos();
	}
*/

public:
	Dist<2+1> mDist;
	Biquad<> mBackShelf;
	Biquad<> mPinnaNotch1, mPinnaNotch2, mPinnaPeak1, mPinnaPeak2;
	float mTorsoAmt=0., mTorsoDelay=0.;
	float mShadows[2] = {1,1}; // TODO: these should be LPFs
	float mEarDist = 0.07; // about half the average bitragion breadth
	float mRoomSize = 3;
};


/// Scene with multiple sources based on HRFilter
template <int Nsrc>
class HRScene{
public:
	typedef HRFilter Source;

	HRScene(){
		for(int i=0; i<2; ++i){
			mReverbs[i].resize(gam::JCREVERB, i*2);
			mReverbs[i].decay(4);
			mReverbs[i].damping(0.25);		
		}
		far(0.5);
	}

	/// Get number of sources
	int numSources() const { return Nsrc; }

	/// Get sound source
	Source& source(int i){ return mSources[i]; }

	/// Set/get source sample
	float& sample(int i){ return mSamples[i]; }

	/// Return next spatialized sample as (left, right, room)
	float2 operator()(){

		float3 spat(0,0,0);
		for(int i=0; i<Nsrc; ++i){
			spat += mSources[i](mSamples[i]);
		}

		spat[2] *= mWallAtten;

		float2 echoes(
			mReverbs[0](spat[2]),
			mReverbs[1](spat[2])
		);
		return spat.get(0,1) + echoes;
	}

	HRScene& far(float v){ for(auto& s:mSources) s.dist().far(v); return *this; }
	HRScene& reverbDecay(float v){ for(auto& r:mReverbs) r.decay(v); return *this; }
	HRScene& reverbDamping(float v){ for(auto& r:mReverbs) r.damping(v); return *this; }
	HRScene& wallAtten(float v){ mWallAtten=v; return *this; }

private:
	Source mSources[Nsrc];
	float mSamples[Nsrc] = {0.f};	// source input samples
	ReverbMS<> mReverbs[2]; 		// one reverb for each ear
	float mWallAtten = 0.1;
};

} // gam::

#endif
