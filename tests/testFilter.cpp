#include <stdio.h>
#include "Gamma/Gamma.h"
#include "Gamma/FFT.h"
#include "Gamma/Filter.h"
#include "Gamma/Print.h"

using namespace gam;

int main(int argc, char* argv[]){

	const int N = 64;		// Number of samples per unit of position
	float buf[N+2];
	RFFT<float> fft(N);

	#define FREQ_RESP(f, description)\
		f.zero();\
		printf("\n%s:\n", description);\
		for(int i=0; i<N; ++i){ buf[i+1]=f(i?0:1); }\
		fft.forward(buf, false, true);\
		for(int i=0; i<N+2; i+=2){\
			Complex<float> c(buf[i], buf[i+1]);\
			float m = c.mag();\
			float p = c.arg() * M_1_PI;\
			printf("[%2d] % 6.3f %6.3f ", i/2, m, p);\
			printPlot(m*0.8, 48);\
			printPlot(p, 48);\
			printf("\n");\
		}

	{
		IIR1<> f;
	}

	{
		IIR2<> f;
	}

	{
		IIRButter<> f;

		f.order(1);
//		f.freq(1./4 * 0.5); FREQ_RESP(f, "Order 1 Butterworth low-pass at 1/4 band")
		f.freq(2./4 * 0.5); FREQ_RESP(f, "Order 1 Butterworth low-pass at 1/2 band")

		f.order(2);
		f.freq(2./4 * 0.5); FREQ_RESP(f, "Order 2 Butterworth low-pass at 1/2 band")

		f.order(3);
		f.freq(2./4 * 0.5); FREQ_RESP(f, "Order 3 Butterworth low-pass at 1/2 band")

		f.order(4);
		f.freq(2./4 * 0.5); FREQ_RESP(f, "Order 4 Butterworth low-pass at 1/2 band")

		f.order(5);
		f.freq(2./4 * 0.5); FREQ_RESP(f, "Order 5 Butterworth low-pass at 1/2 band")

		printf("\n");
		f.order(6);
		f.freq(2./4 * 0.5); FREQ_RESP(f, "Order 6 Butterworth low-pass at 1/2 band")
	}

	{
		IIRCheby<> f;

		f.set(1./4, 1);

		f.order(1);		
//		f.freq(1./4 * 0.5); FREQ_RESP(f, "Order 1 Chebyshev low-pass at 1/4 band")
		f.freq(2./4 * 0.5); FREQ_RESP(f, "Order 1 Chebyshev low-pass at 1/2 band")

		f.order(2);
		f.freq(2./4 * 0.5); FREQ_RESP(f, "Order 2 Chebyshev low-pass at 1/2 band")

		f.order(3);
		f.freq(2./4 * 0.5); FREQ_RESP(f, "Order 3 Chebyshev low-pass at 1/2 band")

		f.order(4);
		f.freq(2./4 * 0.5); FREQ_RESP(f, "Order 4 Chebyshev low-pass at 1/2 band")

		f.order(5);
		f.freq(2./4 * 0.5); FREQ_RESP(f, "Order 5 Chebyshev low-pass at 1/2 band")

		f.order(6);
		f.freq(2./4 * 0.5); FREQ_RESP(f, "Order 6 Chebyshev low-pass at 1/2 band")
	}

	return 0;
}

