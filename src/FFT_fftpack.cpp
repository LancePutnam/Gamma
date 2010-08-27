#include "Gamma/FFT.h"
#include "fftpack++.h"

namespace gam{


template <class T>
struct CFFT<T>::Impl{

	Impl(int sz): n(0), wsave(0)
	{
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
void CFFT<T>::forward(T * buf, bool normalize){
	fftpack::cfftf(&mImpl->n, buf, mImpl->wsave, mImpl->ifac);
	
	if(normalize){
		T m = T(1)/size();
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
struct RFFT<T>::Impl{

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
void RFFT<T>::forward(T * buf, bool normalize){
	fftpack::rfftf(&mImpl->n, buf, mImpl->wsave, mImpl->ifac);
	
	if(normalize){
		T m = T(2)/size();
		for(int i=1; i<size()-1; ++i) buf[i] *= m;
		m *= T(0.5);
		buf[0] *= m;
		buf[size()-1] *= m;
	}
}
	
template <class T>
void RFFT<T>::inverse(T * buf, bool normalized){
	if(normalized){
		T m = 1./T(2);
		for(int i=1; i<size()-1; ++i) buf[i] *= m;	
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


