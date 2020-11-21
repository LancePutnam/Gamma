/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include <stdio.h>
#include "Gamma/Domain.h"

namespace gam{


DomainObserver::DomainObserver()
:	mSubject(0)
{
	//printf("DomainObserver::DomainObserver() - %p\n", this);
	domain(Domain::master());
}

DomainObserver::DomainObserver(const DomainObserver& rhs)
:	mSubject(0)
{
	//printf("DomainObserver::DomainObserver(const DomainObserver&) - %p\n", this);
	Domain& s = rhs.mSubject ? *rhs.mSubject : Domain::master();
	domain(s);
}

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

void DomainObserver::domain(Domain& newSubject){
	if(&newSubject != mSubject){
		if(mSubject){
			// If head of list, then set head to right node
			if(mSubject->mHeadObserver == this){
				mSubject->mHeadObserver = this->nodeR;
			}

			nodeRemove();
		}
		newSubject.attach(*this);
		double r = newSubject.spu() / (mSubject ? mSubject->spu() : 1.);
		mSubject = &newSubject;
		onDomainChange(r);
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
		s->onDomainChange(r);
		s = s->nodeR;
	}
}

void Domain::spu(double v){ //printf("[%p] Domain::spu(%g)\n", this, v);
	mHasBeenSet = true;
	if(v != mSPU){
		double r = v/mSPU;
		mSPU = v;
		mUPS = 1. / v;
		notifyObservers(r);	// calls onDomainChange() of each observer
	}
}

void Domain::ups(double val){ spu(1./val); }

void Domain::print() const {
	printf("Domain %p:\n\tspu = %f, ups = %f\n", this, spu(), ups());

	DomainObserver * o = mHeadObserver;
	unsigned numObs = 0;
	while(o){
		++numObs;
		o = o->nodeR;
	}
	
	printf("\t %u observers%s", numObs, numObs ? ": " : "\n");

	if(numObs){
		o = mHeadObserver;
		while(o){
			printf("%p ", o);
			o = o->nodeR;
		}
		printf("\n");
	}
}

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
