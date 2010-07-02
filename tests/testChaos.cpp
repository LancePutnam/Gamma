#include "Gamma/RecurrenceMaps.h"
#include "Gamma/scl.h"
#include "Gamma/Print.h"

using namespace gam;

int main(int argc, char ** argv){

//	MapLin<double> ln1(1.);
//	
//	ln1.coef(-1., -0.5);
//	
//	for(int i=0; i<32; i++){
//		//ln1.x = SclOp::wrap(ln1.x);
//		//printf("[% 3i] %.12f\n", i, ln1.x);
//		SclOp::printPlot(ln1.x); printf("\n");
//		ln1.next();
//	}

	MapLin2<double> l2(1., 1.);
	
	//l2.coef(0.2, 0.8, -0.2, 0.2, -0.2, 0.8);
	
	for(int i=0; i<32; i++){
		//ln1.x = SclOp::wrap(ln1.x);
		//printf("[% 3i] %.12f\n", i, ln1.x);
		printPlot(l2.v[0]); printf("\n");
		l2();
	}

	//TMapLin4<double> l4;
	MapLin4<double> l4(1,2,3,4);
	MapQuad<double> q1(1);
	MapQuad2<double> q2(1,2);
	MapQuad3<double> q3(1,2,3);

	return 0;
}

