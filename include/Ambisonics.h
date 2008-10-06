#ifndef AMBISONICS_H_INC
#define AMBISONICS_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

/*
Modified from the CSL Ambisonic code written by Florian Hollerweger, 
Graham Wakefield, and Jorge Castellanos, 2005.
*/

#include <stdio.h>
#include "pstdint.h"			/* for uint32_t, uint16_t, etc... */
#include "MacroD.h"

#define USE_GAMMA

#ifdef USE_GAMMA
	#include "scl.h"
	#define COS			gam::scl::cosT8
	#define SIN			gam::scl::sinT7
	#define WRAP(phase)	phase = gam::scl::wrapPhase(phase)
#else
	#include <math.h>
	#define COS			cos
	#define SIN			sin
	#define WRAP(phase)
#endif

#ifndef TEM
#define TEM template <class T>
#endif
#ifndef SAFE_FREE
#define SAFE_FREE(ptr) if(ptr){ free(ptr); ptr=0; }
#endif
//#define MAX_ORDER 3

namespace gam{

/// Ambisonic base class
class AmbiBase{
public:
	AmbiBase(uint32_t dim, uint32_t order);
	virtual ~AmbiBase(){}

	uint32_t numChannels() const;	///< Returns total number of Ambisonic domain channels
	void order(uint32_t order);		///< Set the order.
	
	virtual void onChannelsChange(){}	///< Called whenever the number of Ambisonic channels changes.

	void print(FILE * fp = stdout, const char * append = "\n");
	
	static uint32_t channelsToUniformOrder(uint32_t channels);
	
	TEM static void encodeWeightsFuMa(T * weights, uint32_t dim, uint32_t order, T azimuth, T elevation);
	
	/// Brute force 3rd order.  Weights must be of size 16.
	TEM static void encodeWeightsFuMa16(T * weights, T azimuth, T elevation);
	TEM static void encodeWeightsFuMa16(T * ws, T xn, T yn, T zn);
	
	static uint32_t orderToChannels(uint32_t dim, uint32_t order);
	static uint32_t orderToChannelsH(uint32_t orderH);
	static uint32_t orderToChannelsV(uint32_t orderV);
	
protected:
	uint32_t mDim;			// dimensions - 2d or 3d
	uint32_t mOrder;		// order - 0th, 1st, 2nd, or 3rd
	uint32_t mNumChannels;	// cached for efficiency
	
	static const double c1_sqrt2, c8_11, c40_11;
};


/// Higher Order Ambisonic Decoding class
TEM class AmbiDecode : public AmbiBase{
public:
	AmbiDecode(uint32_t dim, uint32_t order, uint32_t numSpeakers, uint32_t flavor=1);
	virtual ~AmbiDecode();

	T decode(uint32_t speakerNum);	///< Decode speaker's sample from stored ambisonic frame.

	void decode(T ** dec, const T ** enc, uint32_t numDecFrames){
		
		// iterate speakers
		for(uint32_t s=0; s<numSpeakers(); ++s){
			T * out = dec[s];
			
			// iterate ambi channels
			for(uint32_t c=0; c<numChannels(); c++){
				T * in = enc[c];
				float w = decodeWeight(s, c);
				
				for(uint32_t i=0; i<numDecFrames; ++i) out[i] += in[i] * w;		
			}		
		}
	}

	void flavor(uint32_t type);
	void numSpeakers(uint32_t num);	///< Set number of speakers.  Positions are zeroed upon resize.
	void setSpeaker(uint32_t index, T azimuth, T elevation=0);
	void setSpeakerDegrees(uint32_t index, T azimuth, T elevation=0);
	void zero();				///< Zeroes out internal ambisonic frame.

	T * azimuths();				///< Returns pointer to speaker azimuths.
	T * elevations();			///< Returns pointer to speaker elevations.
	T * frame() const;			///< Returns pointer to ambisonic channel frame used by decode(uint32_t)
	uint32_t numSpeakers();		///< Returns number of speakers.
	uint32_t flavor();			///< Returns decode flavor.

	virtual void onChannelsChange();
	
	void print(FILE * fp = stdout, const char * append = "\n");
	
	T decodeWeight(uint32_t speaker, uint32_t channel){ 
		return mWChan[channel] * mDecodeMatrix[speaker * mNumChannels + channel];
	}

protected:
	uint32_t mNumSpeakers;
	uint32_t mFlavor;		// decode flavor

	T * mDecodeMatrix;		// deccoding matrix for each ambi channel & speaker
							// cols are channels and rows are speakers								
	T mWOrder[5];			// weights for each order
	T * mWChan;				// weights for each ambi channel
	T * mPositions;			// speakers' azimuths + elevations
	T * mFrame;				// an ambisonic channel frame used for decode(uint32_t)

	void updateChanWeights();
	void resizeArrays(uint32_t numChannels, uint32_t numSpeakers);

	T decode(T * encFrame, uint32_t encNumChannels, uint32_t speakerNum);	// is this useful?
	
	static T flavorWeights[4][5][5];
};


/// Higher Order Ambisonic encoding class
TEM class AmbiEncode : public AmbiBase{
public:
	AmbiEncode(uint32_t dim, uint32_t order);
	virtual ~AmbiEncode();
	
	void encode   (const AmbiDecode<T> &dec, T input);	///< Encode input sample and set decoder frame.
	void encodeAdd(const AmbiDecode<T> &dec, T input);	///< Encode input sample and add to decoder frame.
	
	//void encodeAdd(AmbiDecode<T> &dec, T input, );
	
	void position(T azimuth, T elevation);			///< Set position of source to be encoded.
	//void position(T x, T y, T z);
	
	T * weights();									///< Returns pointer to encoding weights.
	
	virtual void onChannelsChange();
	
protected:
	T * mWeights;			// encoding weights
	uint32_t mSizeWeights;
	
	void weightsResize(uint32_t numChannels);
};



// Implementation ______________________________________________________________

// AmbiBase

inline uint32_t AmbiBase::numChannels() const { return mNumChannels; }

inline void AmbiBase::order(uint32_t o){
	if(o != mOrder){
		mOrder = o;
		mNumChannels = orderToChannels(mDim, mOrder);
		onChannelsChange();
	}
}

inline uint32_t AmbiBase::channelsToUniformOrder(uint32_t channels){
	// M = floor(sqrt(N) - 1)
	return (uint32_t)floor(sqrt((double)channels) - 1);
}

TEM void AmbiBase::encodeWeightsFuMa
	(T * ws, uint32_t dim, uint32_t order, T az, T el)
{	
	*ws++ = (T)c1_sqrt2;							// W = 1/sqrt(2)
	
	if(order > 0){
		WRAP(az);
		WRAP(el);
	
		T cosel = COS(el);
		T x = COS(az) * cosel;
		T y = SIN(az) * cosel;
		T x2 = x * x;
		T y2 = y * y;
		
		*ws++ = x;							// X = cos(A)cos(E)	
		*ws++ = y;							// Y = sin(A)cos(E)
		
		if(order > 1){
			x2 = x*x;
			y2 = y*y;
			
			*ws++ = x2 - y2;					// U = cos(2A)cos2(E) = xx-yy
			*ws++ = (T)2 * x * y;				// V = sin(2A)cos2(E) = 2xy
			
			if(order > 2){
				*ws++ = x * (x2 - (T)3 * y2);		// P = cos(3A)cos3(E) = X(X2-3Y2)
				*ws++ = y * (y2 - (T)3 * x2);		// Q = sin(3A)cos3(E) = Y(3X2-Y2)
			}
		}
		
		if(dim == 3){
			T z = SIN(el);
			*ws++ = z;							// Z = sin(E)
			
			if(order > 1){
				T z2 = z*z;
		
				*ws++ = (T)2 * z * x;				// S = cos(A)sin(2E) = 2zx
				*ws++ = (T)2 * z * y;				// T = sin(A)sin(2E) = 2yz
				*ws++ = ((T)1.5 * z2) - (T)0.5;		// R = 1.5sin2(E)-0.5 = 1.5zz-0.5
				
				if(order > 2){
					T pre = (T)c40_11 * z2 - (T)c8_11;
					
					*ws++ = z * (x2-y2) * (T)0.5;		// N = cos(2A)sin(E)cos2(E) = Z(X2-Y2)/2
					*ws++ = x * y * z;					// O = sin(2A)sin(E)cos2(E) = XYZ
					*ws++ = pre * x;					// L = 8cos(A)cos(E)(5sin2(E) - 1)/11 = 8X(5Z2-1)/11
					*ws++ = pre * y;					// M = 8sin(A)cos(E)(5sin2(E) - 1)/11 = 8Y(5Z2-1)/11
					*ws   = z * ((T)2.5 * z2 - (T)1.5);	// K = sin(E)(5sin2(E) - 3)/2 = Z(5Z2-3)/2
				}
			}
		}
	}
}

// [x, y, z] is the normalized direction vector
TEM void AmbiBase::encodeWeightsFuMa16(T * ws, T x, T y, T z){
	T x2 = x * x;
	T y2 = y * y;
	T z2 = z * z;
	T pre = (T)c40_11 * z2 - (T)c8_11;

	ws[ 0] = (T)c1_sqrt2;					// W channel, shouldn't it be already defined?
	ws[ 1] = x;								// X = cos(A)cos(E)	
	ws[ 2] = y;								// Y = sin(A)cos(E)
	ws[ 3] = z;								// Z = sin(E)
	ws[ 4] = x2 - y2;						// U = cos(2A)cos2(E) = xx-yy
	ws[ 5] = (T)2. * x * y;					// V = sin(2A)cos2(E) = 2xy
	ws[ 6] = (T)2. * z * x;					// S = cos(A)sin(2E) = 2zx
	ws[ 7] = (T)2. * z * y;					// T = sin(A)sin(2E) = 2yz
	ws[ 8] = (T)1.5 * z2 - (T)0.5;			// R = 1.5sin2(E)-0.5 = 1.5zz-0.5
	ws[ 9] = x * (x2 - (T)3. * y2);			// P = cos(3A)cos3(E) = X(X2-3Y2)
	ws[10] = y * (y2 - (T)3. * x2);			// Q = sin(3A)cos3(E) = Y(3X2-Y2)
	ws[11] = z * (x2 - y2) * (T)0.5;		// N = cos(2A)sin(E)cos2(E) = Z(X2-Y2)/2
	ws[12] = x * y * z;						// O = sin(2A)sin(E)cos2(E) = XYZ
	ws[13] = pre * x;						// L = 8cos(A)cos(E)(5sin2(E) - 1)/11 = 8X(5Z2-1)/11
	ws[14] = pre * y;						// M = 8sin(A)cos(E)(5sin2(E) - 1)/11 = 8Y(5Z2-1)/11
	ws[15] = z * ((T)2.5 * z2 - (T)1.5);	// K = sin(E)(5sin2(E) - 3)/2 = Z(5Z2-3)/2	
}


TEM void AmbiBase::encodeWeightsFuMa16(T * ws, T az, T el){

	WRAP(az);
	WRAP(el);
	
	T cosel = COS(el);
	T x = COS(az) * cosel;
	T y = SIN(az) * cosel;
	T z = SIN(el);
	
	encodeWeightsFuMa16(ws, x, y, z);
}

inline uint32_t AmbiBase::orderToChannels(uint32_t dim, uint32_t order){
	uint32_t chans = orderToChannelsH(order);
	return dim == 2 ? chans : chans + orderToChannelsV(order);
}

inline uint32_t AmbiBase::orderToChannelsH(uint32_t orderH){ return (orderH << 1) + 1; }
inline uint32_t AmbiBase::orderToChannelsV(uint32_t orderV){ return orderV * orderV; }


// AmbiEncode

TEM AmbiEncode<T>::AmbiEncode(uint32_t dim, uint32_t order)
	: AmbiBase(dim, order), mWeights(0), mSizeWeights(0)
{
	onChannelsChange();
}

TEM AmbiEncode<T>::~AmbiEncode(){
	SAFE_FREE(mWeights);
}

TEM inline void AmbiEncode<T>::position(T az, T el){
	AmbiBase::encodeWeightsFuMa(mWeights, mDim, mOrder, az, el);
}

TEM inline void AmbiEncode<T>::encode(const AmbiDecode<T> &dec, T input){	
	for(uint32_t c=0; c<dec.numChannels(); ++c) dec.frame()[c] = weights()[c] * input;
}

TEM inline void AmbiEncode<T>::encodeAdd(const AmbiDecode<T> &dec, T input){
	for(uint32_t c=0; c<dec.numChannels(); ++c) dec.frame()[c] += weights()[c] * input;
}

TEM inline T * AmbiEncode<T>::weights(){ return mWeights; }

TEM void AmbiEncode<T>::weightsResize(uint32_t numChannels){
	if(numChannels != mSizeWeights){
		mSizeWeights = numChannels;
		mWeights = (T *)realloc(mWeights, mSizeWeights * sizeof(T));
		memset(mWeights, 0, mSizeWeights * sizeof(T));
	}
}

TEM void AmbiEncode<T>::onChannelsChange(){
	weightsResize(numChannels());
}


// AmbiDecode

TEM T AmbiDecode<T>::flavorWeights[4][5][5] = {
	{	// none:
		{1,		1,		1,		1,		1    }, // n = 0, M = 0, 1, 2, 3, 4
		{0,		1,		1,		1,		1    }, // n = 1, M = 0, 1, 2, 3, 4
		{0,		0,		1,		1,		1    }, // n = 2, M = 0, 1, 2, 3, 4
		{0,		0,		0,		1,		1    }, // n = 3, M = 0, 1, 2, 3, 4
		{0,		0,		0,		0,		1    }  // n = 4, M = 0, 1, 2, 3, 4
	},{	// default:
		{1,		0.707,	0.707,	0.707,	0.707},	// n = 0, M = 0, 1, 2, 3, 4
		{0,		1    ,	0.75 ,	0.75 ,	0.75 },	// n = 1, M = 0, 1, 2, 3, 4
		{0,		0    ,	0.5  ,	0.5  ,	0.5  },	// n = 2, M = 0, 1, 2, 3, 4
		{0,		0    ,	0    ,	0.3  ,	0.3  },	// n = 3, M = 0, 1, 2, 3, 4
		{0,		0    ,	0    ,	0    ,	0.1  } 	// n = 4, M = 0, 1, 2, 3, 4
	},{	// in phase
		{1,		1,		1,		1,		1    },	// n = 0, M = 0, 1, 2, 3, 4
		{0,		0.333,	0.5,	0.6,	0.667},	// n = 1, M = 0, 1, 2, 3, 4
		{0,		0,		0.1,	0.2,	0.286},	// n = 2, M = 0, 1, 2, 3, 4
		{0,		0,		0,		0.029,	0.071},	// n = 3, M = 0, 1, 2, 3, 4
		{0,		0,		0,		0,		0.008}	// n = 4, M = 0, 1, 2, 3, 4	
	},{	// max-rE
		{1,		1,		1,		1,		1    },	// n = 0, M = 0, 1, 2, 3, 4
		{0,		0.577,	0.775,	0.861,	0.906},	// n = 1, M = 0, 1, 2, 3, 4
		{0,		0,		0.4,	0.612,	0.732},	// n = 2, M = 0, 1, 2, 3, 4
		{0,		0,		0,		0.305,	0.501},	// n = 3, M = 0, 1, 2, 3, 4
		{0,		0,		0,		0,		0.246}	// n = 4, M = 0, 1, 2, 3, 4
	}
};

TEM AmbiDecode<T>::AmbiDecode(uint32_t dim, uint32_t order, uint32_t numSpeakers, uint32_t flavor)
	: AmbiBase(dim, order),
	mNumSpeakers(0), mDecodeMatrix(0), mWChan(0), mPositions(0), mFrame(0)
{
	resizeArrays(numChannels(), numSpeakers);
	this->flavor(flavor);
}

TEM AmbiDecode<T>::~AmbiDecode(){
	SAFE_FREE(mDecodeMatrix);
	SAFE_FREE(mWChan);
	SAFE_FREE(mPositions);
}

TEM inline T AmbiDecode<T>::decode(T * encFrame, uint32_t encNumChannels, uint32_t speakerNum){
	T smp = (T)0;
	T * dec = mDecodeMatrix + speakerNum * mNumChannels;
	T * wc = mWChan;
	for(; encNumChannels>0; encNumChannels--) smp += *dec++ * *wc++ * *encFrame++;
	return smp;
}

TEM inline T AmbiDecode<T>::decode(uint32_t speakerNum){
	return decode(mFrame, numChannels(), speakerNum);
}

TEM void AmbiDecode<T>::flavor(uint32_t type){
	if(type < 4){
		mFlavor = type;
		for(int i=0; i<4 ; i++) mWOrder[i] = flavorWeights[type][i][mOrder];
		updateChanWeights();
	} 
}

TEM void AmbiDecode<T>::numSpeakers(uint32_t num){
	resizeArrays(numChannels(), num);
}

TEM void AmbiDecode<T>::zero(){ memset(mFrame, 0, numChannels() * sizeof(T)); }

TEM inline T * AmbiDecode<T>::azimuths(){ return mPositions; }
TEM inline T * AmbiDecode<T>::elevations(){ return mPositions + mNumSpeakers; }
TEM inline T * AmbiDecode<T>::frame() const { return mFrame; }
TEM inline uint32_t AmbiDecode<T>::flavor(){ return mFlavor; }
TEM inline uint32_t AmbiDecode<T>::numSpeakers(){ return mNumSpeakers; }

TEM void AmbiDecode<T>::setSpeaker(uint32_t index, T az, T el){
	if (index < numSpeakers()) {	// verify speaker index
		azimuths()  [index] = az;	// update speaker location	
		elevations()[index] = el;
		
		// update encoding weights
		//mDecodeMatrix[index][0] *= AmbiBase::c1_sqrt2;
		encodeWeightsFuMa(mDecodeMatrix + index * numChannels(), mDim, mOrder, az, el);
	}	
}

TEM void AmbiDecode<T>::setSpeakerDegrees(uint32_t index, T az, T el){
	setSpeaker(index, az * (T)0.01745329252, el * (T)0.01745329252);
}

TEM void AmbiDecode<T>::updateChanWeights(){
	T * wc = mWChan;
	*wc++ = mWOrder[0];	
	
	if (mOrder > 0) {
		*wc++ = mWOrder[1];			// X
		*wc++ = mWOrder[1];			// Y
		if (mOrder > 1) {
			*wc++ = mWOrder[2];		// U
			*wc++ = mWOrder[2];		// V
			if (mOrder > 2) {
				*wc++ = mWOrder[3];	// P
				*wc++ = mWOrder[3];	// Q
			}
		}
		
		if(3 == mDim){
			*wc++ = mWOrder[1];			// Z
			if (mOrder > 1) {
				*wc++ = mWOrder[2];		// S
				*wc++ = mWOrder[2];		// T
				*wc++ = mWOrder[2];		// R
				if (mOrder > 2) {
					*wc++ = mWOrder[3];	// N
					*wc++ = mWOrder[3];	// O
					*wc++ = mWOrder[3];	// L
					*wc++ = mWOrder[3];	// M
					*wc   = mWOrder[3];	// K
				}
			}
		}
	}
}

TEM void AmbiDecode<T>::resizeArrays(uint32_t numChannels, uint32_t numSpeakers){

	uint32_t oldSize = mNumChannels * mNumSpeakers;
	uint32_t newSize = numChannels * numSpeakers;
	
	if(oldSize != newSize){
	
		// resize decoding matrix
		uint32_t bytesSize = newSize * sizeof(T);
		mDecodeMatrix = (T *)realloc(mDecodeMatrix, bytesSize);
		memset(mDecodeMatrix, 0, bytesSize);
		
		// resize ambi channel weights
		uint32_t bytesChan = numChannels * sizeof(T);
		mWChan = (T *)realloc(mWChan, bytesChan);
		memset(mWChan, 0, bytesChan);

		mFrame = (T *)realloc(mFrame, bytesChan);
		memset(mFrame, 0, bytesChan);
		
		// resize number of speakers (?)
		if(numSpeakers != mNumSpeakers){
			mNumSpeakers = numSpeakers;
			uint32_t bytesPos = numSpeakers * 2 * sizeof(T);
			mPositions = (T *)realloc(mPositions, bytesPos);
			memset(mPositions, 0, bytesPos);
		}
		
		// recompute decode matrix weights
		for(uint32_t i=0; i<numSpeakers; i++){
			setSpeaker(i, azimuths()[i], elevations()[i]);
		}
		
		updateChanWeights();
	}
	
	
	mNumChannels = numChannels;
}

TEM void AmbiDecode<T>::onChannelsChange(){
	resizeArrays(numChannels(), mNumSpeakers);
}

TEM void AmbiDecode<T>::print(FILE * fp, const char * append){
	AmbiBase::print(stdout, ", ");
	fprintf(fp, "s:%3d%s", mNumSpeakers, append);
}

} // gam::

#include "MacroU.h"

#undef TEM
#undef SAFE_FREE

#undef WRAP
#undef COS
#undef SIN

#endif

