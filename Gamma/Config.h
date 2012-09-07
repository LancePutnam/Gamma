#ifndef GAMMA_CONFIG_H_INC
#define GAMMA_CONFIG_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information

	File Description:
	Defines platform-dependent macros
*/

#if defined(__APPLE__) && defined(__MACH__)
	#define GAM_OSX 1

#elif defined(_WIN32) || defined(WIN32) || defined(_WIN64) || defined(WIN64)
	#define GAM_WINDOWS 1

#else
	#define GAM_LINUX 1

#endif

#endif
