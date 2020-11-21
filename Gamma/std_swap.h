#ifndef GAMMA_STD_SWAP_H_INC
#define GAMMA_STD_SWAP_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information

	File Description:
	C++ version aware std::swap include
*/

#include <cmath>
#if __cplusplus <= 199711L
	#include <algorithm> // swap
#else
	#include <utility> // swap
#endif

#endif
