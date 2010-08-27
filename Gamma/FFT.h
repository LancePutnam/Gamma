#ifndef GAMMA_FFT_H_INC
#define GAMMA_FFT_H_INC

#include "Gamma/Containers.h"


namespace gam{


/// One-dimensional real-to-complex fast Fourier transform

/// The complex sequence format for forward and inverse transforms is
/// [2r0, r1, i1, ... , r(n/2-1), i(n/2-1), 2r(n/2)]
template <class T>
class RFFT{
public:

	/// @param[in] size		size of real input sequence; most efficient when a product of small primes
	RFFT(int size=0);
	
	~RFFT();

	/// Get size of transform
	int size() const;
	
	/// Perform real-to-complex forward transform in-place
	void forward(T * buf, bool normalize=true);
	
	/// Perform complex-to-real inverse transform in-place
	void inverse(T * buf, bool normalized=true);

	/// Set size of transform
	void resize(int n);
	

private:
	class Impl; Impl * mImpl;
};




/// One-dimensional complex fast Fourier transform

/// The complex sequence format for forward and inverse transforms is
/// [r0, i0, r1, i1, ... , r(n-1), i(n-1)]
template <class T>
class CFFT{
public:

	/// @param[in] size		size of complex input sequence; most efficient when a product of small primes
	CFFT(int size=0);
	
	~CFFT();

	/// Get size of transform
	int size() const;

	/// Perform forward transform in-place
	void forward(T * buf, bool normalize=true);

	/// Perform inverse transform in-place
	void inverse(T * buf);

	/// Set size of transform
	void resize(int n);

	template <template <class> class ComplexType>
	void forward(ComplexType<T> * buf, bool normalize=true){
		forward((T*)buf, normalize);
	}

	template <template <class> class ComplexType>
	void inverse(ComplexType<T> * buf){ inverse((T*)buf); }

private:
	class Impl; Impl * mImpl;
};





/// One-dimensional, stridable complex fast Fourier transform

/// Sizes must be powers or two. Transforms can operate on strided arrays.
/// Modified from:
/// http://www.jjj.de/fft/cplxfft.h
template <class T>
class CFFT2{
public:

	typedef Complex<T> complex;

    CFFT2(uint32_t size=0,						// size is power of 2
		bool doBitRev=true,
		T scalef1 = 0.5, T scalef2 = 1.0,	// fwd transform scalings
		T scalei1 = 1.0, T scalei2 = 1.0	// rev xform
	);

    ~CFFT2();
	
	/// Perform in-place forward transform
    void forward(complex * buf, uint32_t stride=1){ fft(buf, 0, stride); }
	
	/// Perform in-place inverse transform
    void inverse(complex * buf, uint32_t stride=1){ fft(buf, 1, stride); }
	
	/// Return size of transform
    uint32_t size() const { return mN; }
	
	/// Set size of transform. Return true upon success, false otherwise.
	bool size(uint32_t n);

    // used to fill in last half of complex spectrum of real signal
    // when the first half is already there.
    //
    void hermitian(complex * buf) const;

private:
    uint32_t mN, mLog2N;	// these define size of FFT buffer
    complex * mW;			// array [N/2] of cos/sin values
    uint32_t * mBitRev;		// bit-reversing table, in 0..N
    T mFScales[2];			// f-transform scalings
    T mIScales[2];			// i-transform scales
	bool mDoBitRev;			// whether to do bit-reversal (convolution does not need it)
	
    void fft(complex * buf, int iflag, uint32_t stride) const;
	void freeMem();
};



//////////////////////////////  CFFT2 methods //////////////////////////////

/*
 * constructor takes an int, power-of-2.
 * scalef1,scalef2, are the post-pass and post-transform
 * scalings for the forward transform; scalei1 and scalei2 are
 * the same for the inverse transform.
 */
template <class T>
CFFT2<T>::CFFT2(uint32_t n, bool doBitRev, T scalef1, T scalef2, T scalei1, T scalei2)
:	mN(0), mLog2N(0), mW(0), mBitRev(0), mDoBitRev(doBitRev)
{
    mFScales[0] = scalef1;
    mFScales[1] = scalef2;
    mIScales[0] = scalei1;
    mIScales[1] = scalei2;
	size(n);
}


template <class T>
CFFT2<T>::~CFFT2(){ freeMem(); }


template <class T>
void CFFT2<T>::freeMem(){
	if(mBitRev){ delete[] mBitRev; mBitRev=0; }
    if(mW     ){ delete[] mW;      mW     =0; }
}


/*
 * hermitian() assumes the array has been filled in with values
 * up to the center and including the center point. It reflects these
 * conjugate-wise into the last half.
 */
template <class T>
void CFFT2<T>::hermitian(complex * buf) const {

    if(size() <= 2) return;		// nothing to do
    int i = (size()>>1)-1;		// input
    int j = i+2;				// output

    while(i > 0){
		buf[j] = buf[i].conj();
		--i; ++j;
    }
}



template <class T>
bool CFFT2<T>::size(uint32_t n){

	if(size() == n) return true;

	freeMem();

	// Convert size to power of two
	uint32_t k=0;
    for(; ; ++k){
		if((1UL<<k) == n) break;
		if(k==14 || (1UL<<k) > n){
			//throw "CFFT2: size not power of 2";
			return false;
		}
	}
	mN = 1<<k;
	mLog2N = k; //printf("N:%d log2(N):%d\n", mN, mLog2N);

	mBitRev = new uint32_t[size()];
	mW = k>0 ? new complex[size()>>1] : 0;

	if(mBitRev == 0 || ((k>0) && mW == 0)){
		//throw "CFFT2: out of memory";
		return false;
	}


	// compute bit-rev table
	mBitRev[0] = 0;

	for(uint32_t j=1; j<size(); j<<=1){
		for(uint32_t i=0; i<j; ++i){
			mBitRev[i] <<= 1;
			mBitRev[i+j] = mBitRev[i]+1;
		}
	}


    // prepare the cos/sin table. This is bit-reversed, and goes
    // like this: 0, 90, 45, 135, 22.5 ...  for N/2 entries.
    if(k > 0){
		k = (1<<(k-1));
		for(uint32_t i=0; i<k; ++i){
			T t = T(mBitRev[i<<1]) * M_PI / T(k);
			complex ww(cos(t), sin(t));
			mW[i] = ww.conj();		// force limiting of imag part if applic.
		}
    }
	
	return true;
}



/*
 * CFFT2::fft_func(buf,0) performs a forward fft on the data in the buffer specified.
 * CFFT2::fft_func(buf,1) performs an inverse fft on the data in the buffer specified.
 */
 /*
template <class T>
void CFFT2<T>::fft(complex * buf, int iflag) const {

	const T * sp = iflag ? mIScales : mFScales;
	T s = sp[0];		// per-pass scale

    if(mLog2N == 0){		// only 1 element !
		buf[0] *= sp[1];	// final scale only
		return;
    }

    // first pass:
    //  1st element  = sum of 1st & middle, middle element = diff.
    // repeat N/2 times.

    uint32_t k = size()>>1;

    if(mLog2N == 1) s *= sp[1];	// final scale
	
    complex * buf2 = buf + k;			// half-way through input buffer
    for(uint32_t i=0; i<k; ++i){		// first pass is faster
		complex z1(buf[i] + buf2[i]);
		complex z2(buf[i] - buf2[i]);
		buf [i] = z1 * s;
		buf2[i] = z2 * s;
    }
    if(mLog2N == 1) return;	// only 2!

    k>>=1;							// k is N/4 now
    complex * bufe = buf+size();	// past end
    for( ; k; k>>=1){
		if(k == 1){			// last pass - include final scale 
			s *= sp[1];		// final scale
		}
		
		complex * buf0 = buf;
		for(uint32_t j=0; buf0<bufe; ++j){
			
			complex zw(iflag ? mW[j].conj() : mW[j]);

			buf2 = buf0+k;
			for(uint32_t i=0; i<k; ++i){	// a butterfly
				complex z1(zw * buf2[i]);
				complex z2(buf0[i] + z1);
				buf2[i] = (buf0[i] - z1)*s;
				buf0[i] = z2 * s;
			}
			buf0 += (k<<1);
		}
    }
	
    // bitrev the sucker
	if(mDoBitRev){
		for(uint32_t i=0; i<size(); ++i){
			uint32_t j = mBitRev[i];
			if(i <= j) continue;		// don't do these
			complex z1(buf[i]);			// swap values
			buf[i] = buf[j];
			buf[j] = z1;
		}
	}
}
*/

template <class T>
void CFFT2<T>::fft(complex * buf, int iflag, uint32_t stride) const {

	const T * sp = iflag ? mIScales : mFScales;
	T s = sp[0];		// per-pass scale

    if(mLog2N == 0){		// only 1 element !
		buf[0] *= sp[1];	// final scale only
		return;
    }

    // first pass:
    //  1st element  = sum of 1st & middle, middle element = diff.
    // repeat N/2 times.
	
	uint32_t N = size()*stride;
    uint32_t k = N>>1;

    if(mLog2N == 1) s *= sp[1];				// final scale
	
    complex * buf2 = buf + k;				// half-way through input buffer
    for(uint32_t i=0; i<k; i+=stride){		// first pass is faster
		complex z1(buf[i] + buf2[i]);
		complex z2(buf[i] - buf2[i]);
		buf [i] = z1 * s;
		buf2[i] = z2 * s;
		//printf("%d ", i);
    }
    if(mLog2N == 1) return;	// only 2!

    k>>=1;									// k is N/4 now
    complex * bufe = buf + N;				// past end
    for( ; k>=stride; k>>=1){
		if(k == stride) s *= sp[1];			// last pass - include final scale 
		
		complex * buf0 = buf;
		complex * w = mW;
		for(uint32_t j=0; buf0<bufe; j+=stride){
			
			complex zw(iflag ? (*w).conj() : *w); w++;

			buf2 = buf0 + k;
			for(uint32_t i=0; i<k; i+=stride){	// a butterfly
				//printf("%d ", i);
				complex z1(zw * buf2[i]);
				complex z2(buf0[i] + z1);
				buf2[i] = (buf0[i] - z1)*s;
				buf0[i] = z2 * s;
			}
			buf0 += (k<<1);
		}
    }
	
    // bitrev the sucker
	if(mDoBitRev){
	
		uint32_t * bitRev = mBitRev;
	
		for(uint32_t i=0; i<N; i+=stride){
			uint32_t j = *bitRev++; j *= stride;
			if(i <= j) continue;		// don't do these
			complex z1(buf[i]);			// swap values
			buf[i] = buf[j];
			buf[j] = z1;
			//printf("%d ", i);
		}
	}
}


/*
// real-to-complex
template <class T>
void CFFT2<T>::fft(T * buf, int iflag) const {

	const T * sp = iflag ? mIScales : mFScales;
	T s = sp[0];		// per-pass scale

    if(mLog2N == 0){		// only 1 element !
		buf[0] *= sp[1];	// final scale only
		return;
    }

    // first pass:
    //  1st element  = sum of 1st & middle, middle element = diff.
    // repeat N/2 times.

    uint32_t k = size()>>1;

    if(mLog2N == 1){
		s *= sp[1];	// final scale
	}
	
    T * buf2 = buf + k;					// half-way through input buffer
    for(uint32_t i=0; i<k; ++i){		// first pass is faster
		T z1 = buf[i] + buf2[i];
		T z2 = buf[i] - buf2[i];
		buf [i] = z1 * s;
		buf2[i] = z2 * s;
    }
    if(mLog2N == 1) return;	// only 2!

    k>>=1;							// k is N/4 now
    T * bufe = buf+size();	// past end
    for( ; k; k>>=1){
		if(k == 1){			// last pass - include final scale 
			s *= sp[1];		// final scale
		}
		
		T * buf0 = buf;
		for(uint32_t j=0; buf0<bufe; ++j){
			
			complex zw(iflag ? mW[j].conj() : mW[j]);

			buf2 = buf0+k;
			for(uint32_t i=0; i<k; ++i){	// a butterfly
				complex z1( zw * buf2[i]);
				complex z2( z1 + buf0[i]);
				buf2[i] = (-z1 + buf0[i])*s;
				buf0[i] = z2 * s;
			}
			buf0 += (k<<1);
		}
    }
	
    // bitrev the sucker 
    for(uint32_t i=0; i<size(); ++i){
		uint32_t j = mBitRev[i];
		if(i <= j) continue;		// don't do these
		T z1 = buf[i];			// swap values
		buf[i] = buf[j];
		buf[j] = z1;
    }
}
*/

/*

Recursive even-odd splitting of input results in a bit-reversed scrambling.

000		000
001		100
010		010
011		110
100		001
101		101
110		011
111		111

01234567
0246 1357
04 26 15 37

0123456789abcdef
02468ace 13579bdf

x = ((x & 0x00ff00ff)<<8) | ((x & 0xff00ff00)>>8);
x = x<<16 | x>>16;

void fft(complex * buf, int n){

	int k = log2(n);
	
	

}


*/


} // gam::

#endif
