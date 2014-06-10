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

//template<class A, class B>
//inline bool near(A a, B b, double eps=1e-6){
//	return scl::abs(a-b) < eps;
//}

int main(int argc, char* argv[]){

	// Unit tests are ordered from the least to most dependent functions/objects
	// in order to catch errors in base functionality.

	// File I/O
	{
		const char * path = "test.txt";
		File f(path, "w");
		assert(f.open());
		
		char buf[] = {'H','e','l','l','o',' ','W','o','r','l','d','!'};
		assert(f.write(buf, 1, sizeof(buf)) == sizeof(buf));
		f.close();
		assert(!f.opened());
		
		assert(File::exists(path));
		
		f.mode("r");
		assert(f.open());
		
		char * bufR = f.readAll();
		for(int i=0; i<f.size(); ++i) assert(buf[i] == bufR[i]);
		f.close();
	}
	
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

	#include "ut/utFFT.cpp"

	#include "ut/ut_gen.cpp"
	#include "ut/ut_mem.cpp"
	#include "ut/ut_scl.cpp"
	#include "ut/ut_ipl.cpp"
	#include "ut/ut_arr.cpp"

	#include "ut/utDelay.cpp"
	#include "ut/utEnvelope.cpp"
	#include "ut/utFilter.cpp"
	#include "ut/utGenerators.cpp"

//	printf("Unit testing succeeded.\n");

	return 0;
}
