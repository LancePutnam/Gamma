#include <stdio.h>
#include "Gamma/Gamma.h"
#include "Gamma/DFT.h"
#include "Gamma/Filter.h"
#include "Gamma/Print.h"

using namespace gam;

int main(int argc, char* argv[]){

	const int N = 32;		// Number of samples per unit of position
	DFT dft(N, 0, Bin::Polar); dft.precise(true);

	#define FREQ_RESP(f, description)\
		f.zero();\
		printf("\n%s:\n", description);\
		for(uint32_t i=0; i<dft.sizeWin(); ++i){ float v=f(i?0:1); dft(v); }\
		dft.bins(0)[0] *= 2; dft.bins(dft.numBins()-1)[0] *= 2;\
		for(uint32_t i=0; i<dft.numBins(); ++i){\
			float m = dft.bins(i)[0] * N * 0.5;\
			float p = dft.bins(i)[1] * M_1_PI;\
			printf("[%2d] % 6.3f %6.3f ", i, m, p);\
			printPlot(m*0.7, 32);\
			printPlot(p, 32);\
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
		f.freq(1./4 * 0.5); FREQ_RESP(f, "Order 1 Butterworth low-pass at 1/4 band")
		f.freq(2./4 * 0.5); FREQ_RESP(f, "Order 1 Butterworth low-pass at 1/2 band")

		f.order(2);
		f.freq(2./4 * 0.5); FREQ_RESP(f, "Order 2 Butterworth low-pass at 1/2 band")

		f.order(3);
		f.freq(2./4 * 0.5); FREQ_RESP(f, "Order 3 Butterworth low-pass at 1/2 band")

		f.order(4);
		f.freq(2./4 * 0.5); FREQ_RESP(f, "Order 4 Butterworth low-pass at 1/2 band")

		f.order(5);
		f.freq(2./4 * 0.5); FREQ_RESP(f, "Order 5 Butterworth low-pass at 1/2 band")
	}

	{
		IIRCheby<> f;

		f.order(1);
		f.set(1./4, 6);
		f.freq(1./4 * 0.5); FREQ_RESP(f, "Order 1 Chebyshev low-pass at 1/4 band")
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

