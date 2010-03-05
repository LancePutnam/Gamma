#ifndef GAMMA_GAMMA_H_INC
#define GAMMA_GAMMA_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information

	File Description:
	Main Gamma includes
*/

#define GAMMA_VERSION 0.9.4
//#define GAMMA_H_INC_ALL

// Core Functions
// Everything else depends on these so always include them.
#include "arr.h"
#include "gen.h"
#include "fil.h"
#include "ipl.h"
#include "mem.h"
#include "scl.h"
#include "tbl.h"
#include "rnd.h"

#include "Containers.h"
#include "Strategy.h"

// Optional includes
#ifdef GAMMA_H_INC_ALL

	// State Functions
	#include "AudioIO.h"
	#include "Chaos.h"
	#include "Delay.h"
	#include "DFT.h"
	#include "Envelope.h"
	#include "Noise.h"
	#include "Oscillator.h"
	#include "Sampler.h"
	#include "SoundFile.h"
	#include "Sync.h"

	// Composite Objects
	#include "Effects.h"

#endif

#endif

