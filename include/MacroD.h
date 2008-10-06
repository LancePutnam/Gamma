#ifndef GAMMA_MACROS
#define GAMMA_MACROS

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include "Types.h"

#define ULONG	uint32_t

#define TEM					template <class T>						/* single type */
#define TEM2				template <class T1, class T2>			/* double type */
#define TEM3				template <class T1, class T2, class T3>	/* triple type */
#define LOOP_P(n, exp)		{for(; n>0; --n){ exp }}				/* pointer loop down */
#define LOOP_1(n, exp)		{for(ULONG i=0; i<(n); ++i) { exp }}	/* array index loop up */
#define LOOP_2(n, exp)		{for(ULONG i=0; i<(n); i+=2){ exp exp }}/* unrolled by 2 */
#define LOOP_S(n, s, exp)	{for(ULONG i=0; i<(n); i+=s){ exp }}	/* strided loop */
#define LOOP				LOOP_1

// Indexer loop
#define LOOP_IND(exp)\
	for(gam::uint j=ind.begin(); ind.cond(j); j+=ind.stride()){ gam::uint i=ind.index(j); exp }

#define processFreq(dft, exp) {\
	float * bins0 = (dft).bins0();\
	float * bins1 = (dft).bins1();\
	for(ULONG b=0; b<(dft).numBins(); ++b){\
		float f0 = *bins0;\
		float f1 = *bins1;\
		exp\
		*bins0++ = f0;\
		*bins1++ = f1;\
	}}
	
#define AUDIO_LOOP for(unsigned long f=0; f<numFrames; f++)

// 2-channel (stereo) frame loop
// The output samples are floats called 's0' and 's1' 
// which get summed into the output buffers.
#define AUDIO_LOOP_2(exp)\
	AUDIO_LOOP{\
		float t0 = 0.f, t1 = 0.f;\
		exp;\
		o0[f] = t0; o1[f] = t1;\
	}

#define AUDIO_PROCESS(name, exp)\
extern "C"{\
	void name(AudioIOData & io){\
		float * o0 = io.out(0);\
		float * o1 = io.out(1);\
		unsigned long numFrames = io.numFrames();\
		exp\
	}\
}

#define processTime(name, exp)\
	AUDIO_PROCESS(name, AUDIO_LOOP_2(exp))


#define beginAudio(name)\
extern "C"{\
	void name(AudioIOData & io){\
		float * out0 = io.out(0);\
		float * out1 = io.out(1);\
		uint32_t numFrames = io.numFrames();

#define endAudio() }}

#define loopTime { for(uint32_t i=0; i<numFrames; ++i)

#define loopFreq(dft){\
	float * bin0 = dft.bins0();\
	float * bin1 = dft.bins1();\
	for(uint32_t i=0; i<dft.numBins(); ++i)
	

//#define SAFE_FREE(p) if(p){ free(p); p=0; }

/*
#define DUFF_DEVICE_8(aCount, aAction) \
{ \
int count_ = (aCount); \
int times_ = (count_ + 7) >> 3; \
switch (count_ & 7){ \
case 0: do { aAction; \
case 7: aAction; \
case 6: aAction; \
case 5: aAction; \
case 4: aAction; \
case 3: aAction; \
case 2: aAction; \
case 1: aAction; \
} while (--times_ > 0); \
} \
}
*/

#endif
