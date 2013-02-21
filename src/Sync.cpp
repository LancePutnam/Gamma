/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include "Gamma/Sync.h"

namespace gam{

// Sync

// Note:	We can't zero the head node's links because it might be constructed
//			after nodes have been appended to it!
//
//			However, if no Synceds are added and notifyObservers() is called we get a crash!
Sync::Sync()
:	mSPU(1.), mUPS(1.), mHeadObserver(true, *this), mHasBeenSet(false)
{}

Sync::Sync(double spuA)
:	mHeadObserver(true, *this)
{ spu(spuA); }

Sync::~Sync(){
//	Synced * s = mHeadSynced.nodeR;
//
//	while(s){
//		s->sync(0);
//		s = s->nodeR;
//	}
}

Sync& Sync::operator<< (Synced& obs){ obs.sync(*this); return(*this); }

void Sync::addObserver(Synced& obs){
	//printf("%p: ", &synced); mHeadSynced.print();
	if(&obs != &mHeadObserver){
		obs.nodeInsertR(mHeadObserver);
	}
	//printf("%p: ", &synced); mHeadSynced.print();
}

void Sync::notifyObservers(double r){
	Synced * s = mHeadObserver.nodeR;

	while(s){	//printf("Sync %p: Notifying %p\n", this, s);
		s->scaleSPU(r);	// this will call onResync()
		s = s->nodeR;
	}
}

void Sync::spu(double v){
	mHasBeenSet = true;
	if(v != mSPU){
		double r = v/mSPU;
		mSPU = v;
		mUPS = 1. / v;
		notifyObservers(r);	// call onResync() of my Synceds
	}
}

void Sync::ups(double val){ spu(1./val); }




// Synced

Synced& Synced::operator= (const Synced& rhs){
	if(this != &rhs){
		if(rhs.mSubject){
			sync(*rhs.mSubject);
		}
	}
	return *this;
}

void Synced::scaleSPU(double v){
	mSPU *= v;
	mUPS = 1. / mSPU;
	onResync(v);
}

void Synced::scaleUPS(double v){ scaleSPU(1./v); }

void Synced::sync(Sync& newSubject){
	if(&newSubject != mSubject){
		if(mSubject) nodeRemove();
		newSubject.addObserver(*this);
		double r = newSubject.spu() / (mSubject ? mSubject->spu() : 1.);
		mSubject = &newSubject;
		scaleSPU(r);	// calls onResync()
	}
}

} // gam::
