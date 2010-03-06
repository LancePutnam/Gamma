#include <math.h>
//#include <stdint.h>			/* for uint32_t, uint16_t, etc... */
#include <stdio.h>
#include "Gamma/Ambisonics.h"

namespace gam{

// AmbiBase

const double AmbiBase::c1_sqrt2 = 1. / sqrt(2.);
const double AmbiBase::c8_11    = 8. / 11.;
const double AmbiBase::c40_11   = 40. / 11.;

AmbiBase::AmbiBase(uint32_t dim, uint32_t order) : mDim(dim){
	this->order(order);
}

void AmbiBase::print(FILE * fp, const char * append){
	fprintf(fp, "o:%2d, d:%2d, c:%3d%s", mOrder, mDim, mNumChannels, append);
}

} // gam::
