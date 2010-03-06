/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include "Gamma/Sync.h"

namespace gam{

// Sync

Sync * Sync::mMaster=0;


// Note:	We can't zero the head node's links because it might be constructed
//			after nodes have been appended to it!
//
//			However, if no Synceds are added and notifySynceds() is called we get a crash!
Sync::Sync()
:	mSPU(1.), mUPS(1.), mHeadSynced(true, *this), mHasBeenSet(false)
{}

Sync::Sync(double spuA)
:	mHeadSynced(true, *this)
{ spu(spuA); }

Sync::~Sync(){
//	Synced * s = mHeadSynced.nodeR;
//
//	while(s){
//		s->sync(0);
//		s = s->nodeR;
//	}
}

Sync& Sync::operator<< (Synced& synced){ synced.sync(*this); return(*this); }

void Sync::addSynced(Synced& synced){
	//printf("%p: ", &synced); mHeadSynced.print();
	if(&synced != &mHeadSynced){
		synced.nodeInsertR(mHeadSynced);
	}
	//printf("%p: ", &synced); mHeadSynced.print();
}

void Sync::notifySynceds(double r){
	Synced * s = mHeadSynced.nodeR;

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
		notifySynceds(r);	// call onResync() of my Synceds
	}
}

void Sync::ups(double val){ spu(1./val); }




// Synced

void Synced::scaleSPU(double v){
	mSPU *= v;
	mUPS = 1. / mSPU;
	onResync(v);
}

void Synced::scaleUPS(double v){ scaleSPU(1./v); }

void Synced::sync(Sync& src){
	if(&src != mSync){
		if(mSync) nodeRemove();
		src.addSynced(*this);
		double r = src.spu() / (mSync ? mSync->spu() : 1.);
		mSync = &src;
		scaleSPU(r);	// calls onResync()
	}
}



} // end namespace gam
