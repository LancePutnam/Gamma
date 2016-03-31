#include "Gamma/FFT.h"
#include "fftpack++.h"

namespace gam{


template <class T>
class CFFT<T>::Impl{
public:
	Impl(int sz): n(0), wsave(0){
		resize(sz);
	}
	
	~Impl(){ freeMem(); }
	
	void resize(int size){
		if(size != n){
			n = size;
			freeMem();
			wsave = new T[4*n+15];
			fftpack::cffti(&n, wsave, ifac);
		}
	}

	void freeMem(){ if(wsave){ delete[] wsave; wsave=0;} }

	int n;
	int ifac[sizeof(int) /*bytes/int*/ * 8 /*bits/byte*/ - 1];
	T * wsave;				// work array
};


template <class T>
CFFT<T>::CFFT(int size)
:	mImpl(new Impl(size))
{}


template <class T>
CFFT<T>::~CFFT(){
	if(mImpl){ delete mImpl; mImpl=0; }
}

template <class T>
void CFFT<T>::forward(T * buf, bool normalize, T nrmGain){
	fftpack::cfftf(&mImpl->n, buf, mImpl->wsave, mImpl->ifac);
	
	if(normalize){
		T m = nrmGain/size();
		for(int i=0; i<size()*2; ++i) buf[i] *= m;
	}
}
	
template <class T>
void CFFT<T>::inverse(T * buf){
	fftpack::cfftb(&mImpl->n, buf, mImpl->wsave, mImpl->ifac);
}

template <class T>
void CFFT<T>::resize(int n){ mImpl->resize(n); }

template <class T>
int CFFT<T>::size() const{ return mImpl->n; }




template <class T>
class RFFT<T>::Impl{
public:
	Impl(int sz): n(0), wsave(0){
		resize(sz);
	}
	
	~Impl(){ freeMem(); }
	
	void resize(int size){
		if(size != n){
			n = size;
			freeMem();
			wsave = new T[2*n+15];
			fftpack::rffti(&n, wsave, ifac);
		}
	}

	void freeMem(){ if(wsave){ delete[] wsave; wsave=0;} }

	int n;
	int ifac[sizeof(int) /*bytes/int*/ * 8 /*bits/byte*/ - 1];
	T * wsave;				// work array
};


template <class T>
RFFT<T>::RFFT(int size)
:	mImpl(new Impl(size))
{}


template <class T>
RFFT<T>::~RFFT(){
	if(mImpl){ delete mImpl; mImpl=0; }
}

template <class T>
void RFFT<T>::forward(T * iobuf, bool complexBuf, bool normalize, T nrmGain){

	T * buf = complexBuf ? iobuf+1 : iobuf;

	fftpack::rfftf(&mImpl->n, buf, mImpl->wsave, mImpl->ifac);

	if(normalize){
		const T m = nrmGain/size();
		for(int i=0; i<size(); ++i) buf[i] *= m;
	}

	if(complexBuf){
		iobuf[       0] = buf[0];
		iobuf[       1] = T(0);
		iobuf[size()+1] = T(0);
	}
}
	
template <class T>
void RFFT<T>::inverse(T * iobuf, bool complexBuf){

	T * buf = iobuf;

	if(complexBuf){
		buf++;
		buf[0] = iobuf[0];
	}

	fftpack::rfftb(&mImpl->n, buf, mImpl->wsave, mImpl->ifac);
}

template <class T>
void RFFT<T>::resize(int n){ mImpl->resize(n); }

template <class T>
int RFFT<T>::size() const{ return mImpl->n; }







// Explicit template instantiations
template class RFFT<float>;
template class RFFT<double>;

template class CFFT<float>;
template class CFFT<double>;

} // gam::


