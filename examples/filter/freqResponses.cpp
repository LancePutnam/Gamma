/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Example:		Transform / Filter Frequency Responses
	Description:	This prints out the frequency responses of several different
					filters with various settings. The left plot shows the
					complex magnitudes and the right plot shows the complex
					phases.
*/

#include "../examples.h"

int main(){

	const int N = 32;		// Number of samples per unit of position
	Sync::master().spu(1);
	DFT dft(N, 0, MAG_PHASE); dft.precise(true);

	#define FREQ_RESP(f, description)\
		printf("\n%s:\n", description);\
			for(uint32_t i=0; i<dft.sizeWin(); ++i){ float v=i?0:1; v=f; dft(v); }\
		for(uint32_t i=0; i<dft.numBins(); ++i){\
			float m = dft.bin(i)[0] * N;\
			float p = dft.bin(i)[1] * M_1_PI;\
			printf("% 6.3f %6.3f ", m, p);\
			printPlot(m*0.7, 32);\
			printPlot(p, 32);\
			printf("\n");\
		}
	
	AllPass1<> allPass1;
	AllPass2<> allPass2(4, 4);
	Biquad<> biquad(4);
	BlockDC<> blockDC(0.5/N);
	BlockNyq<> blockNyq(0.5/N);
	Delay<> delay(N);
	Comb<> comb(N, 1, 0);
	Hilbert<> hilbert;
	OnePole<> onePole;
	MovingAvg<> movingAvg(4);
	Notch<> notch(2, 0.5/N);
	Reson<> reson(2, 0.5/N);
	
	allPass1.zero(); allPass1.freq(1./4*0.5); FREQ_RESP(allPass1(v), "1st-order all-pass at 1/4 band")
	allPass1.zero(); allPass1.freq(2./4*0.5); FREQ_RESP(allPass1(v), "1st-order all-pass at 1/2 band")
	allPass1.zero(); allPass1.freq(3./4*0.5); FREQ_RESP(allPass1(v), "1st-order all-pass at 3/4 band")
	allPass1.zero(); allPass1.freq(1./4*0.5); FREQ_RESP(allPass1.low(v), "1st-order low-pass at 1/4 band")
	allPass1.zero(); allPass1.freq(1./4*0.5); FREQ_RESP(allPass1.high(v), "1st-order high-pass at 1/4 band")

	allPass2.width(0.1);
	allPass2.zero(); allPass2.freq(1./4*0.5); FREQ_RESP(allPass2(v), "2nd-order all-pass at 1/4 band")
	allPass2.zero(); allPass2.freq(2./4*0.5); FREQ_RESP(allPass2(v), "2nd-order all-pass at 1/2 band")
	allPass2.zero(); allPass2.freq(3./4*0.5); FREQ_RESP(allPass2(v), "2nd-order all-pass at 3/4 band")

	// Note: res=1 will give us the same phase response as a 1st-order all-pass
	biquad.res(1);
	biquad.zero(); biquad.freq(1./4*0.5); FREQ_RESP(biquad(v), "Biquad low-pass at 1/4 band")
	biquad.zero(); biquad.freq(2./4*0.5); FREQ_RESP(biquad(v), "Biquad low-pass at 1/2 band")
	biquad.zero(); biquad.freq(3./4*0.5); FREQ_RESP(biquad(v), "Biquad low-pass at 3/4 band")

	
	FREQ_RESP(hilbert(v).i, "Hilbert filter (90 degree phase shift)")

	onePole.zero(); onePole.freq(1./8*0.5); FREQ_RESP(onePole(v), "One-pole at 1/8 band")
	onePole.zero(); onePole.freq(1./4*0.5); FREQ_RESP(onePole(v), "One-pole at 1/4 band")
	onePole.zero(); onePole.freq(3./8*0.5); FREQ_RESP(onePole(v), "One-pole at 3/8 band")

	notch.zero(); notch.freq(1./4*0.5); FREQ_RESP(notch(v), "Two-zero notch at 1/4 band")
	notch.zero(); notch.freq(2./4*0.5); FREQ_RESP(notch(v), "Two-zero notch at 1/2 band")
	notch.zero(); notch.freq(3./4*0.5); FREQ_RESP(notch(v), "Two-zero notch at 3/4 band")

	reson.zero(); reson.freq(1./4*0.5); FREQ_RESP(reson(v), "Two-pole reson at 1/4 band")
	reson.zero(); reson.freq(2./4*0.5); FREQ_RESP(reson(v), "Two-pole reson at 1/2 band")
	reson.zero(); reson.freq(3./4*0.5); FREQ_RESP(reson(v), "Two-pole reson at 3/4 band")
	
	printf("%d\n", delay.size());
	delay.zero(); delay.delay(1); FREQ_RESP(delay(v), "Delay of 1 sample")
	delay.zero(); delay.delay(2); FREQ_RESP(delay(v), "Delay of 2 samples")
	
	comb.zero(); comb.set(1, 1  , 0  ); FREQ_RESP(comb(v)*comb.normFfd(), "Low-pass FFD comb at full band")
	comb.zero(); comb.set(2, 1  , 0  ); FREQ_RESP(comb(v)*comb.normFfd(), "Low-pass FFD comb at 1/2 band")
	comb.zero(); comb.set(1,-1  , 0  ); FREQ_RESP(comb(v)*comb.normFfd(), "High-pass FFD comb at full band")
	comb.zero(); comb.set(2,-1  , 0  ); FREQ_RESP(comb(v)*comb.normFfd(), "High-pass FFD comb at 1/2 band")
	comb.zero(); comb.set(2, 0  , 0.5); FREQ_RESP(comb(v)*comb.normFbk(), "Low-pass FBK comb at 1/2 band")
	comb.zero(); comb.set(2, 0  ,-0.5); FREQ_RESP(comb(v)*comb.normFbk(), "High-pass FBK comb at 1/2 band")
	comb.zero(); comb.set(2, 0.5,-0.5); FREQ_RESP(comb(v), "All-pass comb at 1/2 band")
	comb.zero(); comb.set(2,-0.5, 0.5); FREQ_RESP(comb(v), "All-pass comb at 1/2 band")
	
	FREQ_RESP(blockDC(v), "DC blocker")
	FREQ_RESP(blockNyq(v), "Nyquist blocker")
	movingAvg.resize(2); movingAvg.zero(); FREQ_RESP(movingAvg(v), "Moving average (N=2)")
	movingAvg.resize(3); movingAvg.zero(); FREQ_RESP(movingAvg(v), "Moving average (N=3)")
	movingAvg.resize(4); movingAvg.zero(); FREQ_RESP(movingAvg(v), "Moving average (N=4)")
	return 0;
}

