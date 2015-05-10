/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example:	Complex to complex Fourier Transform
Author:		Lance Putnam, 2012

Description:
This demonstrates how to do a complex fast Fourier transform using the CFFT
class.
*/

#include "Gamma/FFT.h"
#include "Gamma/Types.h"
#include "Gamma/Print.h"
using namespace gam;

int main(){

	const int N = 16;		// Transform size
	Complex<float> sig[N];	// Time-/position-domain signal
	Complex<float> buf[N];	// Transform buffer
	CFFT<float> fft(N);		// Complex FFT

	// Create signal
	for(int i=0; i<N; ++i){
		float p = float(i)/N*M_2PI;
		sig[i] = 1;
		sig[i]+= Polar<float>(p);
		sig[i]+= Polar<float>(2*p);
	}
	
	// CFFT operates in-place, so we copy our input signal to another buffer
	for(int i=0; i<N; ++i) buf[i] = sig[i];

	// Perform forward transform (time/position to frequency)
	fft.forward(buf);

	// Print out frequency-domain samples
	for(int i=0; i<N; ++i){
		Complex<float> c = buf[i];
		
		printf("[%2d] ", i);
		printf("% 5.2f % 5.2f ", c[0], c[1]);
		printPlot(c[0], 32);
		printPlot(c[1], 32);
		printf("\n");
	}

	// Perform inverse transform (frequency to time/position)
	fft.inverse(buf);

	// Print out original signal versus forward/inverse transformed signal
	printf("\n");
	for(int i=0; i<N; ++i){
		float sr = sig[i].r;
		float si = sig[i].i;
		float br = buf[i].r;
		float bi = buf[i].i;
		printf("[%2d] ", i);
		printf("(% 5.2f, % 5.2f) (% 5.2f, % 5.2f) ", sr, si, br, bi);
		printPlot(sr/3, 16); printPlot(si/3, 16);
		printPlot(br/3, 16); printPlot(bi/3, 16);
		printf("\n");
	}
}
