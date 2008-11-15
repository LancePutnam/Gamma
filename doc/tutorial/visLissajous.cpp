/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Tutorial:		Visual
	Description:	
*/

#include "tutorial.h"

int main(int argc, char* argv[]){

	const int NP=512, N1=51, N2=N1*N1;
	float pixels[N2];
	mem::zero(pixels, N2);
	
	// draw axes
	for(int i=0; i<N1; i++) pixels[i + N1*(N1/2)] = 0.02;
	for(int i=N1/2; i<N2; i+=N1) pixels[i] = 0.02;

//	Quadra<> qua1(1/7., 0.8,-NP*2);
//	Quadra<> qua2(1/3., 0.1,-NP*2);

	float shift = 0.0004;
	Quadra<> qua1( 1./NP + shift, 0.50/1.,-NP*2);
	Quadra<> qua2( 8./NP + shift, 0.10/1.,-NP*2);

	for(int i=0; i<NP; ++i){
		float x = qua1.val[0] + qua2.val[0];
		float y = qua1.val[1] + qua2.val[1];
		int ix = scl::posToInd( x, N1);
		int iy = scl::posToInd(-y, N1);
		pixels[iy*N1 + ix] += 0.125;
		qua1();
		qua2();
	}

	scl::print2D(pixels, N1, N1);	

	return 0;
}
