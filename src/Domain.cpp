/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include "Gamma/Domain.h"

namespace gam{

#define DOM_OBS_INIT mSubject(0), mSPU(1), mUPS(1)
DomainObserver::DomainObserver()
:	Node2<DomainObserver>(), DOM_OBS_INIT
{
	//printf("Synced::Synced() - %p\n", this);
	domain(Domain::master());
}

DomainObserver::DomainObserver(const DomainObserver& rhs)
:	Node2<DomainObserver>(), DOM_OBS_INIT
{
	//printf("DomainObserver::DomainObserver(const DomainObserver&) - %p\n", this);
	Domain& s = rhs.mSubject ? *rhs.mSubject : Domain::master();
	domain(s);
}
#undef DOM_OBS_INIT

DomainObserver::~DomainObserver(){
	if(mSubject && mSubject->mHeadObserver == this){
		if(nodeR){
			mSubject->mHeadObserver = nodeR;
		}
		else{
			mSubject->mHeadObserver = NULL;
		}
	}
}

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

void DomainObserver::scaleUPS(double v){
	scaleSPU(1./v);
}

void DomainObserver::domain(Domain& newSubject){
	if(&newSubject != mSubject){
		if(mSubject) nodeRemove();
		newSubject.attach(*this);
		double r = newSubject.spu() / (mSubject ? mSubject->spu() : 1.);
		mSubject = &newSubject;
		scaleSPU(r);	// calls onDomainChange()
	}
}



Domain::Domain()
:	mSPU(1.), mUPS(1.), mHeadObserver(NULL), mHasBeenSet(false)
{}

Domain::Domain(double spuA)
:	mHeadObserver(NULL)
{
	spu(spuA);
}

Domain::~Domain(){
//	DomainObserver * s = mHeadObserver.nodeR;
//
//	while(s){
//		s->sync(0);
//		s = s->nodeR;
//	}
}

Domain& Domain::operator<< (DomainObserver& obs){
	obs.domain(*this);
	return *this;
}

void Domain::attach(DomainObserver& obs){
	if(mHeadObserver){
		// insert new observer onto front of list
		obs.nodeInsertL(*mHeadObserver);
	}
	
	// The head observer is always the most recently attached
	mHeadObserver = &obs;

	//printf("%p: ", &obs); mHeadObserver->print();
}

void Domain::notifyObservers(double r){
	DomainObserver * s = mHeadObserver;

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

/*static*/ Domain& Domain::master(){
	static Domain * s = new Domain;
	return *s;
}

/*static*/ void sampleRate(double samplesPerSecond){
	Domain::master().spu(samplesPerSecond);
}

/*static*/ double sampleRate(){
	return Domain::master().spu();
}


} // gam::
