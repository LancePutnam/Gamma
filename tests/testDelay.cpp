#include <stdio.h>
#include "Gamma/Delay.h"
#include "Gamma/DFT.h"
#include "Gamma/Oscillator.h"
#include "Gamma/gen.h"
#include "Gamma/scl.h"
#include "Gamma/Print.h"

using namespace gam;

int main(int argc, char* argv[]){

	const int n = 16;
	Sync::master().spu(n);
	
	#define SIZE(o) printf("%s\t%d\n", #o, sizeof(o));
//	SIZE(AllPass1<>)
//	typedef AllPass1<float, float, Synced1> AllPass11; SIZE(AllPass11)
//	SIZE(Hilbert<>)
//	SIZE(Hilbert2<>)

//	Hilbert<> hilbert;
//	printf("%f\n", hilbert.cf0.freq());
//	printf("%f\n", hilbert.cf1.freq());
//
//	Hilbert2<> hilbert2;
//	printf("%f\n", hilbert2.cf0.freq());
//	printf("%f\n", hilbert2.cf1.freq());
	
/*	
	printf("\nDelay:\n");
	Delay<float, ipl::Linear> delay(1, 5./n);
	gen::Sin<> genSin(M_2PI/n, M_PI_2);

	for(int i=0; i<2*n; ++i){
		float i0 = genSin();
		float o0 = delay(i0);
		printf("[%2d]", i);
		printf("  % 6.3f  % 6.3f  ", i0, o0);
		scl::printPlot(i0, 16, true); printf("  "); scl::printPlot(o0, 16);
		printf("\n");
	}
	
	
	printf("\nAllPass1:\n");
	LFO lfo(2);
	AllPass1<> allPass1(3);

	for(int i=0; i<2*n; ++i){
		float i0 = lfo.up();
		float o0 = allPass1.high(i0 * 0.5);
		printf("[%2d]", i);
		printf("  % 6.3f  % 6.3f  ", i0, o0);
		scl::printPlot(i0, 32, true); printf("  "); scl::printPlot(o0, 32);
		printf("\n");
	}
*/	

	DFT dft(32, 0, Bin::Polar); dft.precise(true);
	gen::Impulse<> sig;

	#define FREQ_RESP(f)\
		printf("\n%s:\n", #f);\
		sig = 1;\
		/*f.zero();*/\
		for(uint32_t i=0; i<dft.sizeWin(); ++i){ float v=sig(); v=f; dft(v); }\
		dft.bins(0)[0] *= 2; dft.bins(dft.numBins()-1)[0] *= 2;\
		for(uint32_t i=0; i<dft.numBins(); ++i){\
			float m = dft.bins(i)[0];\
			float p = dft.bins(i)[1] * M_1_PI;\
			printf("% 6.3f %6.3f ", m, p);\
			printPlot(m*n*0.5, 32);\
			printPlot(p, 32);\
			printf("\n");\
		}
	
	AllPass1<> allPass1;
	AllPass2<> allPass2(4, 4);
	Biquad<> biquad(4);
	BlockDC<> blockDC(0.5);
	BlockNyq<> blockNyq(0.5);
	Delay<> delay(1);
	Comb<> comb(1, 1, 0);
	Hilbert<> hilbert;
	OnePole<> onePole;
	MovingAvg<> movingAvg(4);
	Notch<> notch(2, 0.5);
	Reson<> reson(2, 0.5);
	
	allPass1.zero(); allPass1.freq(2); FREQ_RESP(allPass1(v))
	allPass1.zero(); allPass1.freq(4); FREQ_RESP(allPass1(v))
	allPass1.zero(); allPass1.freq(6); FREQ_RESP(allPass1(v))
	allPass1.zero(); allPass1.freq(2); FREQ_RESP(allPass1.low(v))
	allPass1.zero(); allPass1.freq(2); FREQ_RESP(allPass1.high(v))

	allPass2.zero(); allPass2.freq(2); FREQ_RESP(allPass2(v))
	allPass2.zero(); allPass2.freq(4); FREQ_RESP(allPass2(v))
	allPass2.zero(); allPass2.freq(6); FREQ_RESP(allPass2(v))

	biquad.zero(); biquad.freq(2); FREQ_RESP(biquad(v))
	biquad.zero(); biquad.freq(4); FREQ_RESP(biquad(v))
	biquad.zero(); biquad.freq(6); FREQ_RESP(biquad(v))

	
	FREQ_RESP(hilbert(v).i)

	onePole.zero(); onePole.freq(2); FREQ_RESP(onePole(v))
	onePole.zero(); onePole.freq(4); FREQ_RESP(onePole(v))
	onePole.zero(); onePole.freq(6); FREQ_RESP(onePole(v))

	notch.zero(); notch.freq(2); FREQ_RESP(notch(v))
	notch.zero(); notch.freq(4); FREQ_RESP(notch(v))
	notch.zero(); notch.freq(6); FREQ_RESP(notch(v))

	reson.zero(); reson.freq(2); FREQ_RESP(reson(v))
	reson.zero(); reson.freq(4); FREQ_RESP(reson(v))
	reson.zero(); reson.freq(6); FREQ_RESP(reson(v))
	
	delay.zero(); delay.delay(1./n); FREQ_RESP(delay(v))
	delay.zero(); delay.delay(2./n); FREQ_RESP(delay(v))
	
	comb.zero(); comb.set(1./n, 1, 0  ); FREQ_RESP(comb(v)*comb.normFfd())
	comb.zero(); comb.set(4./n, 1, 0  ); FREQ_RESP(comb(v)*comb.normFfd())
	comb.zero(); comb.set(1./n,-1, 0  ); FREQ_RESP(comb(v)*comb.normFfd())
	comb.zero(); comb.set(4./n,-1, 0  ); FREQ_RESP(comb(v)*comb.normFfd())
	comb.zero(); comb.set(4./n, 0, 0.5); FREQ_RESP(comb(v)*comb.normFbk())
	comb.zero(); comb.set(4./n, 0,-0.5); FREQ_RESP(comb(v)*comb.normFbk())
	comb.zero(); comb.set(4./n, 0.5,-0.5); FREQ_RESP(comb(v))
	comb.zero(); comb.set(4./n,-0.5, 0.5); FREQ_RESP(comb(v))
	
	FREQ_RESP(blockDC(v))
	FREQ_RESP(blockNyq(v))
	FREQ_RESP(movingAvg(v))
	return 0;
}

