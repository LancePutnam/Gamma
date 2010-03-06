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
#include "Gamma/arr.h"
#include "Gamma/gen.h"
#include "Gamma/fil.h"
#include "Gamma/ipl.h"
#include "Gamma/mem.h"
#include "Gamma/scl.h"
#include "Gamma/tbl.h"
#include "Gamma/rnd.h"

#include "Gamma/Containers.h"
#include "Gamma/Strategy.h"

// Optional includes
#ifdef GAMMA_H_INC_ALL

	// State Functions
	#include "Gamma/AudioIO.h"
	#include "Gamma/Chaos.h"
	#include "Gamma/Delay.h"
	#include "Gamma/DFT.h"
	#include "Gamma/Envelope.h"
	#include "Gamma/Noise.h"
	#include "Gamma/Oscillator.h"
	#include "Gamma/Sampler.h"
	#include "Gamma/SoundFile.h"
	#include "Gamma/Sync.h"

	// Composite Objects
	#include "Gamma/Effects.h"

#endif

#endif

