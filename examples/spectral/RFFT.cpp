/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
Example:	Fast Fourier Transform
Author:		Lance Putnam, 2012

Description:
This demonstrates how to do a real-to-complex fast Fourier transform using the
RFFT class. Note that if you are interested in doing Fourier transforms on
real-time audio streams, see the DFT and STFT classes.
*/

#include "Gamma/FFT.h"
#include "Gamma/Types.h"
#include "Gamma/Print.h"
using namespace gam;

int main(){

	const int N = 32;	// Transform size
	float sig[N];		// Time-/position-domain signal
	float buf[N];		// Transform buffer
	RFFT<float> fft(N);	// Real-to-complex FFT
	
	// Create signal
	for(int i=0; i<N; ++i){
		float p = float(i)/N*M_2PI;
		sig[i] = 1 + cos(p) + cos(2*p) + sin(3*p);
		sig[i]+= i&1?-1:1;
	}
	
	// CFFT operates in-place, so we copy our input signal to another buffer
	for(int i=0; i<N; ++i) buf[i] = sig[i];

	// Perform real-to-complex forward transform (time/position to frequency)
	fft.forward(buf);

	// Print out frequency-domain samples
	int numBins = N/2 + 1;

	for(int i=0; i<numBins; ++i){
		Complex<float> c;
		if		(  0 == i)	c(buf[  0], 0);
		else if	(N/2 == i)	c(buf[N-1], 0);
		else				c(buf[i*2-1], buf[i*2]);
		
		printf("[%2d] ", i);
		printf("% 5.2f % 5.2f ", c[0], c[1]);
		printPlot(c[0], 32);
		printPlot(c[1], 32);
		printf("\n");
	}

	// Perform complex-to-real inverse transform (frequency to time/position)
	fft.inverse(buf);

	// Print out original signal versus forward/inverse transformed signal
	printf("\n");
	for(int i=0; i<N; ++i){
		printf("[%2d] ", i);
		printf("% 5.2f % 5.2f ", sig[i], buf[i]);
		printPlot(sig[i]/5, 32);
		printPlot(buf[i]/5, 32);
		printf("\n");
	}
}
