#ifndef GAMMA_GAMMA_H_INC
#define GAMMA_GAMMA_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information

	File Description:
	Main Gamma includes
*/

/*! \mainpage Gamma - Generic synthesis library

	\section intro_sec About

	Gamma is a cross-platform, C++ library for doing generic synthesis and 
	filtering of numerical data. It contains numerous mathematical functions, 
	common algebraic types, such as vectors, complex numbers, and quaternions, 
	an assortment of sequence generators and many objects for signal processing. 
	It is oriented towards real-time sound and graphics rendering, but is 
	equally useful for non-real-time tasks.

*/

#define GAMMA_VERSION 0.9.4x
//#define GAMMA_H_INC_ALL

// Core Functions
// Everything else depends on these so always include them.
#include "Gamma/Containers.h"
#include "Gamma/Strategy.h"
#include "Gamma/Types.h"

#include "Gamma/arr.h"
#include "Gamma/gen.h"
#include "Gamma/ipl.h"
#include "Gamma/mem.h"
#include "Gamma/scl.h"
#include "Gamma/tbl.h"
#include "Gamma/rnd.h"

// Optional includes
#ifdef GAMMA_H_INC_ALL

	// System/Utility
	#include "Gamma/AudioIO.h"
	#include "Gamma/Conversion.h"
	#include "Gamma/File.h"
	#include "Gamma/Print.h"
	#include "Gamma/Timer.h"

	// Generators/Filters
	#include "Gamma/Access.h"
	#include "Gamma/Delay.h"
	#include "Gamma/DFT.h"
	#include "Gamma/Envelope.h"
	#include "Gamma/FFT.h"
	#include "Gamma/Filter.h"
	#include "Gamma/FormantData.h"
	#include "Gamma/Noise.h"
	#include "Gamma/Oscillator.h"
	#include "Gamma/Player.h"
	#include "Gamma/Recorder.h"
	#include "Gamma/SoundFile.h"
	#include "Gamma/Sync.h"
	#include "Gamma/UnitMaps.h"

	// Composite Objects
	#include "Gamma/Effects.h"

#endif

#endif

