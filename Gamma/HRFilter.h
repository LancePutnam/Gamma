#ifndef GAM_INC_HR_FILTER_H
#define GAM_INC_HR_FILTER_H

#include <cmath>
#include <algorithm>
#include "Gamma/Filter.h"
#include "Gamma/Noise.h"
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

	HRFilter(){}

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
		Vec3 Ur(headPose[ 0], headPose[ 1], headPose[ 2]);
		Vec3 Uu(headPose[ 4], headPose[ 5], headPose[ 6]);
		Vec3 Ub(headPose[ 8], headPose[ 9], headPose[10]); // pose stores back vector since right-handed
		
		//auto earR = Ur*mEarDist;
		//auto headToSrc = sourcePos - Vec3(headPose[12], headPose[13], headPose[14]);
		//const auto& S = headToSrc;

		mDist.near(mEarDist);
		//mDist.dist(0, (S + earR).mag()); // left channel
		//mDist.dist(1, (S - earR).mag()); // right channel
		mDist.dist(2, mRoomSize); // diameter of room

		for(int i=0; i<2; ++i){
			auto& e = mEarFilters[i];
			auto earVec = Ur*(i==0?-mEarDist:mEarDist);
			auto earToSource = sourcePos - (Vec3(headPose[12], headPose[13], headPose[14]) + earVec);
			auto dist = earToSource.mag() + 1e-8;
			mDist.dist(i, dist);
			auto Su = earToSource / dist; // ear to source direction vector

			float rightness = Su.dot(Ur); // how much source is to right of ear, in [-1,1]
			float pinnaStrength = std::abs(rightness)*0.5+0.5;

			//float hsLvl = Su.dot(Ub)*-0.15 + 0.85; // [0.7, 1] -> [back,front]
			float hsLvl = Su.dot(Ub)*-0.25 + 0.75; // [0.5, 1] -> [back,front]
			e.backShelf.level(hsLvl);
	
			// Pinna filtering
			auto elVec = Su - (Ur*rightness); // projection of source dir onto up-back plane
			elVec.normalize();
			float srcUp = Uu.dot(elVec); // how much source is up, [-1,1], sine of elevation angle
			float pinnaFrq1 = srcUp*2000. + 8000.; // [6k, 10k] -> [front, above]
			//float pinnaFrq1 = std::max(srcUp, 0.f)*4000. + 6000.;
			// Notch weak above, and strong in front (-10 dB / 0.3)
			// Following a half-period raised cosine
			float pinnaLvl1 = srcUp*srcUp*0.7+0.3; // [0.3, 1]
			pinnaLvl1 = (pinnaLvl1-1.)*pinnaStrength + 1.;
			e.pinnaNotch1.set(pinnaFrq1, 24./1., gam::PEAKING);
			e.pinnaNotch1.level(pinnaLvl1);

			float srcB = Ub.dot(elVec);
			float pinnaFrq2 = srcB*1000. + 10000.; // [9k, 11k] -> [front, above, behind]
			// Notch weak above, and strong in front (-10 dB / 0.3)
			// Following a half-period raised cosine
			float pinnaLvl2 = srcUp*srcUp*0.8+0.2; // [0.2, 1]
			pinnaLvl2 = (pinnaLvl2-1.)*pinnaStrength + 1.;
			e.pinnaNotch2.set(pinnaFrq2, 48./1., gam::PEAKING);
			e.pinnaNotch2.level(pinnaLvl2);

			float above = srcUp>0. ? srcUp : 0.;
			//e.pinnaPeak2.level(srcUp*srcUp*1.+1.);
			e.pinnaPeak2.level(above + 1.);

			e.shadow = (i==0?-rightness:rightness)*0.25+0.75;
		}
		
		//mShadows[0] = headShadow(-Ur, (S + earR).normalize()); // left
		//mShadows[1] = headShadow( Ur, (S - earR).normalize()); // right

		/* Torso filtering (only useful if body orientation available!)
		mTorsoAmt = al::abs(srcUp)*0.1;
		mTorsoDelay = 1./1000. * mTorsoAmt/0.3;
		//*/

		//printf("front/back high-shelf level: %f\n", mBackShelf.level());
		//printf("pinna notch 1: %.2f Hz @ %.2f\n", mPinnaNotch1.freq(), mPinnaNotch1.level());
		//printf("pinna notch 2: %.2f Hz @ %.2f\n", mPinnaNotch2.freq(), mPinnaNotch2.level());
		//printf("head shadow on left ear: %f\n", headShadow(-Ur, (S + earR).normalize()));
		return *this;
	}

	/// Return spatialized sample as (left, right, room)
	float3 operator()(float src){
		//src += mDist.delayLine().read(mTorsoDelay) * mTorsoAmt;
		auto res = mDist(src);
		for(int i=0;i<2;++i) res[i] = mEarFilters[i](res[i]);
		return res;
	}

public:
	Dist<2+1> mDist;

	struct EarFilter{
		Biquad<>
			backShelf {4000, 0.707, gam::HIGH_SHELF},
			pinnaPeak1{4000, 1.000, gam::PEAKING}, // concha resonance
			pinnaPeak2{8000, 4.000, gam::PEAKING}, // for above localization
			pinnaNotch1, pinnaNotch2;
		float shadow = 1.; // TODO: should be LPF

		EarFilter(){
			pinnaPeak1.level(2); 
		}

		float operator()(float s){
			s = backShelf(s);
			s = pinnaNotch1(s);
			s = pinnaNotch2(s);
			//s = pinnaPeak1(s); // not convinced this helps
			s = pinnaPeak2(s);
			return s * shadow;
		}
	};

	EarFilter mEarFilters[2];

	float mTorsoAmt=0., mTorsoDelay=0.;
	float mEarDist = 0.07; // about half the average bitragion breadth
	float mRoomSize = 3;
};


/// Scene with multiple sources based on HRFilter
template <int Nsrc>
class HRScene{
public:
	typedef HRFilter Source;

	HRScene(){
		for(int i=0; i<Nsrc; ++i){
			mSamples[i] = 0.f;
			mActive[i] = true;
		}
		for(int i=0; i<2; ++i){
			mReverbs[i].resize(gam::JCREVERB, i*2);
			mReverbs[i].decay(4);
			mReverbs[i].damping(0.25);		
		}
		far(0.5);
	}

	/// Get number of sources
	int numSources() const { return Nsrc; }
	int size() const { return Nsrc; }

	/// Get sound source
	Source& source(int i){ return mSources[i]; }

	/// Set source sample
	float& sample(int i){ return mSamples[i]; }

	/// Get source sample
	const float& sample(int i) const { return mSamples[i]; }

	/// Set all source samples to zero
	HRScene& zeroSamples(){ for(auto& v : mSamples) v=0.f; return *this; }

	/// Set all samples to zero and make all sources inactive
	HRScene& clear(){
		for(int i=0; i<Nsrc; ++i){
			mSamples[i] = 0.f;
			mActive[i] = false;
		}
		return *this;
	}

	/// Set whether a source is active
	HRScene& active(int i, bool v){ mActive[i]=v; return *this; }
	bool active(int i) const { return mActive[i]; }

	/// Return next spatialized sample as (left, right, room)
	float2 operator()(){

		// Any zero-valued input will devastate the CPU due to denormals 
		// emerging and propagating through all the IIR filters. We add some 
		// inaudible noise to all input samples to prevent this.
		float noise = mNoise();

		float3 spat(0,0,0);
		for(int i=0; i<Nsrc; ++i){
			if(mActive[i]) spat += mSources[i](mSamples[i] + noise);
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
	float mSamples[Nsrc];			// source input samples
	bool mActive[Nsrc];
	ReverbMS<> mReverbs[2]; 		// one reverb for each ear
	float mWallAtten = 0.1;
	NoiseBinary<RNGMulCon> mNoise{1e-20, 17};
};

} // gam::

#endif
