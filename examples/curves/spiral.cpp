/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Example:		Curve / Logarithmic Spiral
	Description:	Drawing 2D logarithmic spirals using a damped harmonic 
					oscillator.
*/

#include "Gamma/Gamma.h"
#include "Gamma/Print.h"
#include "Gamma/Oscillator.h"

int main(){
	using namespace gam;

	const int NP=256;				// Number of points in curve
	const int N1=31, N2=N1*N1;		// Size of graph
	float pixels[N2];				// Accumulation buffer
	
	CSine<> csin;

	// Plot several winding and damping amounts
	for(int k=1; k<4; ++k){
	for(int j=1; j<4; ++j){

		printf("\nwinding=%d, decay=%d\n\n", k, j);
		mem::zero(pixels, N2);
		
		csin.decay(NP*j);			// Number of iterations for amp to decay by 0.0001
		csin.freq(k/(float)NP);		// Winding number
		csin.reset();				// Reset oscillator to (1,0)

		for(int i=0; i<NP; ++i){
			Complex<> c = csin();
			int ix = posToInd( c.r, N1);	// map real component to x
			int iy = posToInd(-c.i, N1);	// map imag component to y
			pixels[iy*N1 + ix] += 0.1;
		}

		print2D(pixels, N1, N1);
		printf("\n");
	}}
}
