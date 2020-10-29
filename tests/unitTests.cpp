#undef NDEBUG
#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <complex>
#define GAMMA_H_INC_ALL
#include "../Gamma/Gamma.h"

using namespace gam;

// #defined in Windef.h!
#ifdef near
#undef near
#endif
inline bool near(double a, double b, double eps=1e-6){
	return scl::abs(a-b) < eps;
}

namespace gam{
namespace scl{
template <class T>
bool almostEqual(const Complex<T>& a, const Complex<T>& b, int maxULP=10){
	return scl::almostEqual(a.r,b.r,maxULP) && scl::almostEqual(a.i,b.i,maxULP);
}
template <class T>
bool almostEqual(const Complex<T>& a, const Polar<T>& b, int maxULP=10){
	Complex<T> t = b;
	return scl::almostEqual(a.r,t.r,maxULP) && scl::almostEqual(a.i,t.i,maxULP);
}
}}


int main(int argc, char* argv[]){

	// Unit tests are ordered from the least to most dependent functions/objects
	// in order to catch errors in base functionality.
	
	// LookupTable
	{
		const int N=4;
		LookupTable<double, ipl::Linear, acc::Wrap> ft(N);
		for(int i=0; i<N; ++i) ft[i]=i;
		
		assert(ft(0./N) == 0);
		assert(ft(1./N) == 1);
		
		assert(ft(0.5/N) == 0.5);
		assert(ft(1.5/N) == 1.5);
	}

	#include "ut/utTypes.cpp"
	#include "ut/utConversion.cpp"
	#include "ut/utContainers.cpp"
	#include "ut/utAccess.cpp"
	#include "ut/utStrategy.cpp"

	#include "ut/utFFT.cpp"

	#include "ut/ut_gen.cpp"
	#include "ut/ut_mem.cpp"
	#include "ut/ut_scl.cpp"
	#include "ut/ut_ipl.cpp"
	#include "ut/ut_arr.cpp"

	#include "ut/utDomain.cpp"

	#include "ut/utDelay.cpp"
	#include "ut/utEnvelope.cpp"
	#include "ut/utFilter.cpp"
	#include "ut/utGenerators.cpp"
	#include "ut/utAnalysis.cpp"

//	printf("Unit testing succeeded.\n");

	return 0;
}
