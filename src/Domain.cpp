/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include "Gamma/Domain.h"

namespace gam{

// Note:	We can't zero the head node's links because it might be constructed
//			after nodes have been appended to it!
//
//			However, if no observers are added and notifyObservers() is called we get a crash!
Domain::Domain()
:	mSPU(1.), mUPS(1.), mHeadObserver(true, *this), mHasBeenSet(false)
{}

Domain::Domain(double spuA)
:	mHeadObserver(true, *this)
{ spu(spuA); }

Domain::~Domain(){
//	DomainObserver * s = mHeadObserver.nodeR;
//
//	while(s){
//		s->sync(0);
//		s = s->nodeR;
//	}
}

Domain& Domain::operator<< (DomainObserver& obs){ obs.domain(*this); return(*this); }

void Domain::attach(DomainObserver& obs){
	//printf("%p: ", &synced); mHeadSynced.print();
	if(&obs != &mHeadObserver){
		obs.nodeInsertR(mHeadObserver);
	}
	//printf("%p: ", &synced); mHeadSynced.print();
}

void Domain::notifyObservers(double r){
	DomainObserver * s = mHeadObserver.nodeR;

	while(s){	//printf("Domain %p: Notifying %p\n", this, s);
		s->scaleSPU(r);	// this will call onDomainChange()
		s = s->nodeR;
	}
}

void Domain::spu(double v){
	mHasBeenSet = true;
	if(v != mSPU){
		double r = v/mSPU;
		mSPU = v;
		mUPS = 1. / v;
		notifyObservers(r);	// calls onDomainChange() of each observer
	}
}

void Domain::ups(double val){ spu(1./val); }




DomainObserver& DomainObserver::operator= (const DomainObserver& rhs){
	if(this != &rhs){
		if(rhs.mSubject){
			domain(*rhs.mSubject);
		}
	}
	return *this;
}

void DomainObserver::scaleSPU(double v){
	mSPU *= v;
	mUPS = 1. / mSPU;
	onDomainChange(v);
}

void DomainObserver::scaleUPS(double v){ scaleSPU(1./v); }

void DomainObserver::domain(Domain& newSubject){
	if(&newSubject != mSubject){
		if(mSubject) nodeRemove();
		newSubject.attach(*this);
		double r = newSubject.spu() / (mSubject ? mSubject->spu() : 1.);
		mSubject = &newSubject;
		scaleSPU(r);	// calls onDomainChange()
	}
}


void sampleRate(double samplesPerSecond){
	Domain::master().spu(samplesPerSecond);
}

double sampleRate(){
	return Domain::master().spu();
}


} // gam::
