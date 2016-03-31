/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Example:		Curve / Circle
	Description:	Draws a circle using complex multiplication.
*/

#include "Gamma/Gamma.h"
#include "Gamma/Print.h"

int main(){
	using namespace gam;

	const int NP=256;				// Number of points in curve
	const int N1=31, N2=N1*N1;		// Size of graph
	float pixels[N2];				// Accumulation buffer

	Complex<double> phase(1,0), freq;
	

	printf("\nUnit circle\n");
	mem::zero(pixels, N2);
	freq.fromPhase(M_2PI/NP);
	for(int i=0; i<NP; ++i){
		int ix = posToInd(phase[0], N1);
		int iy = posToInd(phase[1], N1);
		phase *= freq;
		pixels[iy*N1 + ix] += 0.125;
	}
	print2D(pixels, N1, N1);


	printf("\nHalf circle\n");
	mem::zero(pixels, N2);
	phase(0.5, 0);
	freq.fromPhase(M_2PI/NP);
	for(int i=0; i<NP; ++i){
		int ix = posToInd(phase[0], N1);
		int iy = posToInd(phase[1], N1);
		phase *= freq;
		pixels[iy*N1 + ix] += 0.125;
	}
	print2D(pixels, N1, N1);	
}
