#include "Gamma.h"
#include "Access.h"
#include "scl.h"
#include "fil.h"

using namespace gam;

int main(int argc, char* argv[]){
	
	const int N = 8;
	double ds = 1./N;
	double eps = 1e-12;
	double A[N];
	
	#define ASSERT(A, a,b,c,d,e,f,g,h) assert(A[0]==a && A[1]==b && A[2]==c && A[3]==d && A[4]==e && A[5]==f && A[6]==g && A[7]==h);

	fil::Resonator<double> res(ds, 1, 0, 1);
	
	assert(res.freq() == ds);
	assert(res.decay() == 1);
	assert(res.ahead() == 1);
	assert(scl::almostEqual(res.behind().i, sin(-2*ds*M_2PI)));
	for(int i=0; i<N; ++i) assert(scl::abs(sin(i*ds*M_2PI) - res().i) < eps);

	res.set(0, 1, 0, 1);
	res.decay(0.5, N);
	for(int i=0; i<N; ++i) res();
	assert(scl::almostEqual(res.r, 0.5) );

	res.set(ds, 0, 0, 1);
	for(int i=0; i<N; ++i) assert(scl::abs(sin(i*ds*M_2PI) - res(i?0:1).i) < eps);

	return 0;
}
