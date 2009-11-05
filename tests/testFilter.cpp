#include "Gamma.h"
#include "Access.h"
#include "scl.h"
#include "fil.h"
#include "Visual.h"

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

	fil::TransferFunc fr;
	//fr.addX(0.1, 0).addY( 0.9, 1);	// one-pole lo-pass
	//fr.addX(0.1, 0).addY(-0.9, 1);	// one-pole hi-pass
	//fr.addX(0.5, 0).addX( 0.5, 1);	// one-zero lo-pass
	//fr.addX(0.5, 0).addX(-0.5, 1);	// one-zero hi-pass
	//fr.addX(0.25,0).addX(-0.5, 1).addY(0.5, 1).gain(2); // 1st order all-pass
	fr.addX(0.25,0).addX(-0.5, 4).addY(0.5, 4).gain(2); // 4th order all-pass comb
	//fr.addX(0.5, 0).addX( 0.5, 8);	// 4-notch lo-pass comb
	//fr.addX(0.5, 0).addX(-0.5, 8);	// 4-notch hi-pass comb
	//fr.addX(0.2, 0).addY( 0.8, 8);	// 4-peak lo-pass comb
	//fr.addX(0.2, 0).addY(-0.8, 8);	// 4-peak hi-pass comb

	for(int i=0; i<32; ++i){
		float f = float(i)/32 * 0.5;
		Complexd::Polar p = fr(f);
		printf("% 6.2f % 6.2f ", p.m, p.p*M_1_PI);
		printPlot(p.m);
		printPlot(p.p*M_1_PI);
		printf("\n");
	}

	return 0;
}
