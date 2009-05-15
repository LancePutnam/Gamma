#ifndef GAMMA_GAMMA_H_INC
#define GAMMA_GAMMA_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#define GAMMA_VERSION 0.9.2
//#define GAMMA_H_INC_ALL

// Core Functions
// Everything else depends on these so always include them.
#include "arr.h"
#include "gen.h"
#include "ipl.h"
#include "mem.h"
#include "scl.h"
#include "tbl.h"
#include "rnd.h"

// Optional includes
#ifdef GAMMA_H_INC_ALL

	#include "Containers.h"

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

