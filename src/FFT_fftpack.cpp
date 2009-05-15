#include "FFT.h"
#include "fftpack.h"

namespace gam{

struct RFFTImpl{

	RFFTImpl(uint32_t sz)
	:	n(0)/*, ifac(0)*/, wsave(0)
	{
		//ifac = new integer_t[sizeof(integer_t)*8];
		resize(sz);
	}
	
	~RFFTImpl(){
		//if(ifac){ delete ifac; ifac=0; }
		freeMem();
	}
	
	void resize(integer_t size){
		if(size != n){
			n = size;
			freeMem();
			wsave = new real_t[2*n+15];
			rffti(&n, wsave, ifac);
		}
	}

	void freeMem(){ if(wsave){ delete wsave; wsave=0;} }

	integer_t n;
	//integer_t * ifac;
	integer_t ifac[sizeof(integer_t)*8];
	real_t * wsave;
};





RFFT::RFFT(uint32_t size)
:	mImpl(0)
{
	mImpl = new RFFTImpl(size);
}


RFFT::~RFFT(){
	if(mImpl){ delete mImpl; mImpl=0; }
}


void RFFT::forward(float * buf, bool normalize){
	rfftf(&mImpl->n, buf, mImpl->wsave, mImpl->ifac);
	
	if(normalize){
		float m = 2./size();
		for(uint32_t i=1; i<size()-1; ++i) buf[i] *= m;
		m *= 0.5f;
		buf[0] *= m;
		buf[size()-1] *= m;
	}
}
	

void RFFT::inverse(float * buf){
	rfftb(&mImpl->n, buf, mImpl->wsave, mImpl->ifac);
}


void RFFT::resize(uint32_t n){ mImpl->resize(n); }

uint32_t RFFT::size() const{ return mImpl->n; }



} // gam::
