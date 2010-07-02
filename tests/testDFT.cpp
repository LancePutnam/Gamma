#include <stdlib.h>
#include <stdio.h>
#include "Gamma/DFT.h"
#include "Gamma/FFT.h"
#include "Gamma/Sync.h"
#include "Gamma/rnd.h"
#include "Gamma/scl.h"
#include "Gamma/tbl.h"
#include "Gamma/Print.h"

//#define TEST_RFFT
//#define TEST_CFFT
#define TEST_DFT
//#define TEST_STFT
	#define SIMPLE_TEST
	//#define OVERLAP_TEST


#define PRINT_BUF(src, len)\
	printf("\nTD samples...\n");\
	for(uint32_t i=0; i<(uint32_t)len; i++){\
		printf("[%2d] % 8.6f ", (int)i, src[i]);\
		printPlot(src[i], 20); printf("\n");\
	}
	
#define PRINT_FREQ(obj)\
	printf("\nFD samples...\n");\
	for(uint32_t i=0; i<obj.numBins(); i++)\
		printf("% f\t % f\n", obj.bins(i)[0], obj.bins(i)[1]);	
	
using namespace gam;

int main(int argc, char* argv[]){

	const int winSize = 16;
	//const int hopSize = winSize/4;
	const int padSize = winSize * 0;
	const int tableSize = winSize * 4;

	float table[tableSize], output[tableSize];
	tbl::sinusoid(table, tableSize, M_PI_2, 32.f);

	Sync::master().spu(44100);	// init manually since there is no timer




#ifdef TEST_RFFT
	{	printf("\nRFFT\n");
	
		const uint32_t size = 16;
		float src[size], buf[size];
		
		//rnd::uniS(src, size);
		tbl::sinusoid(src, size, M_PI_2, 0);
		
		mem::copy(buf, src, size);

		RFFT rfft(size);
		
		rfft.forward(buf);
		
		for(int i=0; i<rfft.size(); i+=2){
			float re=buf[i], im=buf[i+1];\
			scl::printPlot(re, 16); scl::printPlot(im, 16); printf("\n");\
		}
		
		rfft.inverse(buf);
		
		float err=0;
		for(int i=0; i<size; ++i){
			err += scl::abs(buf[i] - src[i]);
		}
		printf("\terror = %f %%\n", err*100);
	}
#endif


#ifdef TEST_CFFT
	{
		printf("\nCFFT\n");

		const uint32_t size = 16;
		Complexf src[size], buf[size];
		
		mem::set(src, 
			gen::rMul(Complexf(1, M_2PI/size, 1), Complexf(1,0)), 
			size);
			
		mem::copy(buf, src, size);
		
		CFFT<float> cfft(size/2);

		#define PRINT\
		for(int i=0; i<size; ++i){\
			float re=buf[i].r, im=buf[i].i;\
			printf("\t[%2d] % 6.3f % 6.3f ", i, re, im);\
			scl::printPlot(re, 16); scl::printPlot(im, 16); printf("\n");\
		} printf("\n");
		PRINT cfft.forward(buf+4);
		PRINT cfft.inverse(buf+4);
		PRINT
		#undef PRINT
		
		Complexf err=0;
		for(int i=0; i<size; ++i){
			err += Complex<float>(scl::abs(buf[i].r - src[i].r), scl::abs(buf[i].i - src[i].i));
		}
		err *= 100;
		printf("\terror = (%f, %f) %%\n", err[0], err[1]);
	}
#endif




#ifdef TEST_DFT
	printf("\nDFT:\n");
	DFT dft(winSize, padSize, Bin::Polar);
	dft.precise(true).print();

	// Scalar-based
	for(int i=0; i<winSize*2; i++){
		float si = table[i];	printf("[%2d] % 8.6f ", i, si);
		if(dft(si)){}
		float so = dft();		printf("% 8.6f\n", so);
	}
	
	// Vector-based
							PRINT_BUF(table, winSize)
	dft.forward(table);		PRINT_FREQ(dft);
	dft.inverse(output);	PRINT_BUF(output, winSize)
#endif



#ifdef TEST_STFT
	printf("\nSTFT:\n");
	//STFT stft(winSize, hopSize, padSize, WinType::Hann, Bin::MagFreq);
	STFT stft(winSize, hopSize, padSize, WinType::Hann, Bin::Rect);
	stft.precise(true).print();

	// Scalar-based
	for(int i=0; i<winSize*2; i++){
		float si = table[i];	printf("[ri:%2d] % 8.6f ", i, si);
		if(stft(si)){	printf("\n");
			for(uint32_t i=0; i<stft.numBins(); ++i)
				printf("[c: %2d] % 6.4f % 8.6f\n", i, stft.bins(i)[0], stft.bins(i)[1]);
		}
		float so = stft();		printf("[ro:%2d] % 8.6f\n", i, so);
	}
	
	// Vector-based
//							PRINT_BUF(table, winSize)
//	stft.forward(table);	PRINT_FREQ(stft);
//	for(int i=0; i<stft.numBins(); ++i){ printf("%f\n", stft.mPhases[i]); }
//	stft.inverse(output);	PRINT_BUF(output, winSize)

/*
	{	printf("\nOverlapped window test:\n");

		const uint32_t nt = 16, nw = 8, nh = 4;
		float t[nt]; mem::zero(t, nt);
		float w[nw]; tbl::hann(w, nw);

		for(uint32_t j=0; j<nt; j+=nh){

			for(uint32_t i=0; i<nw; ++i){
				uint32_t ind = i+j;
				if(ind < nt) t[ind] += w[i];
			}
			
			printf("\n");
			for(uint32_t k=0; k<nt; ++k){
				printf("\t[%2d] % 8.6f", k, t[k]); scl::printPlot(t[k]/4.); printf("\n");
			}
		}
		
		gen::Val<float> max;
		arr::add(max, w, Loop(nw, nh));
		printf("\n\tmax sum = %f\n", max());
	}
*/

	{	printf("\nSTFT overlap resynthesis test:\n");
	
		STFT stft(8, 4, 0, WinType::Rectangle);
	
		for(int i=0; i<16; ++i){
			stft(1.f);
			float v = stft();
			printf("\t[%2d] % 8.6f", i, v); scl::printPlot(v); printf("\n");
		}
	}


/*
	STFT * stft = new STFT(dftSize, hopSize, padSize);
	
	#ifdef  SIMPLE_TEST
	printf("\nInput TD samples...\n");
		for(unsigned long i=0; i<dftSize; i++){
			printf("[%2lu] % 8.6f ", i, table[i]);
			SclOp::printPlot(table[i], 20); printf("\n");
		}

	printf("\nFD samples...\n");
		stft->format(ComplexType::Polar);
		stft->forward(table);

		float * mag = stft->spct(0);
		float * phs = stft->spct(1);

		for(unsigned long i=0; i<stft->numBins(); i++){
			printf("% f\t % f\n", mag[i], phs[i]);
		}

	printf("\nOutput TD samples...\n");
		float * out = stft->inverse();
		
		for(unsigned long i=0; i<dftSize; i++){
			printf("[%2lu] % 8.6f ", i, out[i]);
			SclOp::printPlot(out[i], 20); printf("\n");
		}
	#endif
	
	#ifdef  OVERLAP_TEST
	float * hop = table;
	
	//stft->format(SPCT_POLAR);
	stft->format(ComplexType::MagFreq);
	
	for(unsigned long h=0; h<tableSize/hopSize; h++){

		printf("\nInput TD samples (%lu)...\n", h);
			for(unsigned long i=0; i<hopSize; i++){
				printf("% .16f\n", hop[i]);
			}
		
		printf("\nFD samples...\n");
			stft->forward(hop);
			
			float * mag = stft->spct(0);
			float * frq = stft->spct(1);
			
			for(unsigned long i=0; i<stft->numBins(); i++){
				printf("% f\t % f\n", mag[i], frq[i]);
			}	printf("\n");
		
		printf("\nOutput TD samples...\n");
			float * out = stft->inverse();
			for(unsigned long i=0; i<hopSize; i++){
				printf("% .16f\n", out[i]);
			}		
		
		hop += hopSize;
	}
	

	#endif
*/
#endif




/*
	printf("\nSliding Window:\n");
	const ULONG sizeSW = 8;
	SlidingWindow<ULONG> sw(sizeSW, 2);

	printf("\nSingle-buffer:\n");
	LOOP(16,
		if(sw(i + 1)){
			{ LOOP(sw.sizeWin(), printf("%2d ", sw.window()[i]); ) }
			printf("\n");
			//sw.slide();
		}
	)

	printf("\nDouble-buffer:\n");
	ULONG buf[sizeSW];
	
	LOOP(16,
		if(sw(buf, i + 1)){
			{ LOOP(sw.sizeWin(), printf("%2d ", buf[i]); ) }
			printf("\n");
		}
	)
*/

///*
	SDFT<float> sdft(winSize, 0, winSize/2 + 1);
	//MemOp::set(table, tableSize, 1.f);
	//MemOp::zero(table, tableSize); table[0] = 1.f;

	for(int i=0; i<tableSize; ++i){
		sdft.forward(table[i]);
		
		{ for(unsigned i=0; i<sdft.numBins(); ++i){
			Complex<float>& c = sdft.bins(i);
			printf("[%d] ", i);
			printf("% 5.3f ", c.r); printf("% 5.3f ", c.i);
			printPlot(c.r, 20); printf(" ");
			printPlot(c.i, 20); printf("\n");
		}} printf("\n");
		
	}
//*/

	return 0;
}

