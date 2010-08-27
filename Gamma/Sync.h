#ifndef GAMMA_SYNC_H_INC
#define GAMMA_SYNC_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include "Gamma/Node.h"

namespace gam{

class Sync;


// Unit
class Synced1{
public:

	virtual ~Synced1(){}

	double spu() const { return 1.; }		///< Returns local samples/unit.
	double ups() const { return 1.; }		///< Returns local units/sample.

	///	Notification that subject's samples/unit has changed.

	/// Any instance state that depends on the samples/unit ratio should be 
	/// updated here. The ratio of the new to the old samples/unit is passed in.
	virtual void onResync(double ratioSPU){}
	
	void spu(double val){}			///< Set local samples/unit.	
	void ups(double val){}			///< Set local units/sample.

protected:
	void initSynced(){ onResync(1); }
};



/// Synchronized unit sampler.

/// A Synced will attempt to copy its local variables from the master Sync
/// upon construction. If the master Sync has not been constructed, it will
/// use its default values of 1. 
/// This class has a reference to a Sync and its own local scaling factors.
/// By default, the reference Sync is Sync::master.
class Synced : public Node2<Synced> {
public:
	/// Constructor
	Synced();
	
	/// Constructor
	Synced(bool zeroLinks, Sync& src);
	virtual ~Synced(){}

	double scaleSPU() const;	///< Returns ratio of my SPU to my Sync's SPU
	double spu() const;			///< Returns local samples/unit.
	double ups() const;			///< Returns local units/sample.
	const Sync * sync() const;	///< Returns reference to my Sync.

	/// Called by my Sync reference after it changes its value.
	
	///	Any instance state that depends on the sampling length should be updated here.
	/// TODO: Maybe this should be pure virtual?
	virtual void onResync(double ratioSPU){}

	void scaleSPU(double v);	///< Scales samples/unit by factor.
	void scaleUPS(double v);	///< Scales units/sample by factor.
	void spu(double v);			///< Set local samples/unit.
	void sync(Sync& src);		///< Set absolute Sync source.
	void ups(double v);			///< Set local units/sample.

protected:
	void initSynced();

private:
	Sync * mSync;	// Reference to my Sync.
	double mSPU;	// Local samples/unit.
	double mUPS;	// Local units/sample.
};





/// Synchronizing unit sampler.
class Sync{
public:

	Sync();

	/// @param[in]	spu		Initial samples/unit.
	Sync(double spu);

	~Sync();

	Sync& operator<< (Synced& synced);	///< Move a Synced into list.
	void notifySynceds(double r);		///< Calls onSyncChange() for all my Synceds.
	void spu(double v);					///< Set samples/unit.  If changed, will notify Synceds.
	void ups(double v);					///< Set units/sample.  If changed, will notify Synceds.

	bool hasBeenSet() const;			///< Returns true if spu has been set at least once.
	double spu() const;					///< Returns samples/unit, i.e. sample rate.
	double ups() const;					///< Returns units/sample, i.e. sample interval.	

	/// Master sync. By default, Synceds will be synced to this.
	static Sync& master(){
	
		// Note: This uses a Construct On First Use Idiom to avoid unpredictable
		// static initialization order.  The memory allocated will get freed from 
		// the heap when the program exits.		
		// Since we can't predict when a class variable is initialized, we
		// will assume that uninitialized variables are set to 0.
		// TODO: there may be a compiler flag to ensure this is the case
		if(!mMaster){
			mMaster = new Sync; //printf("new master %p\n", mMaster);
		}
		return *mMaster;
	}

protected:
	static Sync * mMaster;
	double mSPU, mUPS;
	Synced mHeadSynced;		// Head of Synced doubly-linked list.
	bool mHasBeenSet;

friend class Synced;
	void addSynced(Synced& synced);
};




// Implementation_______________________________________________________________

// Sync

inline bool Sync::hasBeenSet() const { return mHasBeenSet; }
inline double Sync::spu() const { return mSPU; }
inline double Sync::ups() const { return mUPS; }


// Synced

#define SYNCED_INIT mSync(0), mSPU(1), mUPS(1)

inline Synced::Synced()
:	Node2<Synced>(), SYNCED_INIT{
	sync(Sync::master());
}

// This ctor is used to create the head Synced in a Sync
inline Synced::Synced(bool zeroLinks, Sync& src)
:	Node2<Synced>(zeroLinks), SYNCED_INIT{
	sync(src);
}

inline void Synced::initSynced(){ if(sync()) spu(sync()->spu()); }
inline double Synced::scaleSPU() const { return sync() ? spu() / sync()->spu() : 1; }
inline void Synced::spu(double v){ scaleSPU(v * ups()); }
inline void Synced::ups(double v){ scaleUPS(v * spu()); }
inline double Synced::spu() const { return mSPU; }
inline double Synced::ups() const { return mUPS; }

inline const Sync * Synced::sync() const { return mSync; }

} // gam::

#endif

