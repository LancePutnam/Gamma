#ifndef GAMMA_MACROS_H_INC
#define GAMMA_MACROS_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include "Types.h"

#define TEM					template <class T>						/* single type */
#define TEM2				template <class T1, class T2>			/* double type */
#define TEM3				template <class T1, class T2, class T3>	/* triple type */
#define LOOP_P(n, exp)		{for(; n>0; --n){ exp }}				/* pointer loop down */
#define LOOP_1(n, exp)		{for(uint32_t i=0; i<(n); ++i) { exp }}	/* array index loop up */
#define LOOP_2(n, exp)		{for(uint32_t i=0; i<(n); i+=2){ exp exp }}/* unrolled by 2 */
#define LOOP_S(n, s, exp)	{for(uint32_t i=0; i<(n); i+=s){ exp }}	/* strided loop */
#define LOOP				LOOP_1

// Indexer loop
#define LOOP_IND(exp)\
	for(gam::uint j=ind.begin(); ind.cond(j); j+=ind.stride()){ gam::uint i=ind.index(j); exp }

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
