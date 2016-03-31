/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Example:		Curve / Botanic
	Description:	A botanic curve made with an amplitude modulated complex 
					oscillator.
*/

#include "Gamma/Gamma.h"
#include "Gamma/Print.h"
#include "Gamma/Oscillator.h"

int main(){
	using namespace gam;

	const int NP=256;				// Number of points in curve
	const int N1=17, N2=N1*N1;		// Size of graph
	float pixels[N2];				// Accumulation buffer

	CSine<> osc(1./NP);
	SineR<> oscAM;

	for(int k=2; k<5; ++k){
	for(int j=0; j<3; ++j){
	
		mem::zero(pixels, N2);

		float w = j*0.3 + 0.2;		// "petalness" factor
		oscAM.set(k/(float)NP, w);
		
		for(int i=0; i<NP; ++i){
			Complex<> c = osc() * (1-w + oscAM());
			int ix = posToInd( c.r, N1);
			int iy = posToInd(-c.i, N1);
			pixels[iy*N1 + ix] += 0.1;
		}

		print2D(pixels, N1, N1);
		printf("\n");
	}}
}
