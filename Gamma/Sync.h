#ifndef GAMMA_SYNC_H_INC
#define GAMMA_SYNC_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include "Gamma/Node.h"

namespace gam{

class Sync;


/// Normalized unit synchronization observer
class Synced1{
public:

	virtual ~Synced1(){}

	double spu() const { return 1.; }		///< Returns local samples/unit
	double ups() const { return 1.; }		///< Returns local units/sample
	const Sync * sync() const { return subjectSingleton(); } 	///< Returns reference to my subject

	///	Notification that subject's samples/unit has changed.

	/// Any instance state that depends on the samples/unit ratio should be 
	/// updated here. The ratio of the new to the old samples/unit is passed in.
	virtual void onResync(double ratioSPU){}
	
	void spu(double val){}			///< Set local samples/unit	
	void ups(double val){}			///< Set local units/sample

protected:
	void initSynced(){ onResync(1); }
	
	static Sync * subjectSingleton();
};



/// Unit synchronization observer

/// A Synced will attempt to copy its local variables from the master Sync
/// upon construction. If the master Sync has not been constructed, it will
/// use its default values of 1. 
/// This class has a reference to a Sync and its own local scaling factors.
/// By default, the reference Sync is Sync::master.
class Synced : public Node2<Synced> {
public:

	Synced();
	
	/// Copy constructor
	
	/// If the argument has a subject, then attach this as an observer to the
	/// same subject. Otherwise, attach as observer of Sync::master().
	Synced(const Synced& rhs);

	virtual ~Synced(){}

	double scaleSPU() const;	///< Returns ratio of my SPU to my Sync's SPU
	double spu() const;			///< Returns local samples/unit
	double ups() const;			///< Returns local units/sample
	const Sync * sync() const;	///< Returns reference to my Sync

	/// Called by my Sync reference after it changes its value
	
	///	Any instance state that depends on the sampling length should be 
	/// updated here.
	virtual void onResync(double ratioSPU){}

	void scaleSPU(double v);	///< Scales samples/unit by factor
	void scaleUPS(double v);	///< Scales units/sample by factor
	void spu(double v);			///< Set local samples/unit
	void sync(Sync& src);		///< Set absolute Sync source
	void ups(double v);			///< Set local units/sample

	Synced& operator= (const Synced& rhs);

protected:
	void initSynced(); ///< To be called from the constructor(s) of derived classes

private:
	friend class Sync;
	Synced(bool zeroLinks, Sync& src);

	Sync * mSubject;	// My subject
	double mSPU;		// Local samples/unit
	double mUPS;		// Local units/sample
};





/// Unit synchronization subject
class Sync{
public:

	Sync();

	/// \param[in]	spu		samples/unit
	Sync(double spu);

	~Sync();

	Sync& operator<< (Synced& obs);		///< Attach observer
	void notifyObservers(double r);		///< Notify observers (\see Synced::onResync)
	void spu(double v);					///< Set samples/unit and notify observers
	void ups(double v);					///< Set units/sample and notify observers

	bool hasBeenSet() const;			///< Returns true if spu has been set at least once
	double spu() const;					///< Returns samples/unit, i.e. sample rate
	double ups() const;					///< Returns units/sample, i.e. sample interval

	/// Master sync. By default, Synceds will be synced to this.
	static Sync& master(){
		static Sync * s = new Sync;
		return *s;
	}

protected:
	double mSPU, mUPS;
	Synced mHeadObserver;	// Head of observer doubly-linked list
	bool mHasBeenSet;

friend class Synced;
	void addObserver(Synced& obs);
};




// Implementation_______________________________________________________________

inline Sync * Synced1::subjectSingleton(){
	static Sync * s = new Sync(1.);
	return s;
}


#define SYNCED_INIT mSubject(0), mSPU(1), mUPS(1)
inline Synced::Synced()
:	Node2<Synced>(), SYNCED_INIT
{
	//printf("Synced::Synced() - %p\n", this);
	sync(Sync::master());
}

// This ctor is used to create the head Synced in a Sync
inline Synced::Synced(bool zeroLinks, Sync& src)
:	Node2<Synced>(zeroLinks), SYNCED_INIT
{
	sync(src);
}

inline Synced::Synced(const Synced& rhs)
:	Node2<Synced>(), SYNCED_INIT
{
	//printf("Synced::Synced(const Synced&) - %p\n", this);
	Sync& s = rhs.mSubject ? *rhs.mSubject : Sync::master();
	sync(s);
}
#undef SYNCED_INIT

inline void Synced::initSynced(){ if(sync()) spu(sync()->spu()); }
inline double Synced::scaleSPU() const { return sync() ? spu() / sync()->spu() : 1; }
inline void Synced::spu(double v){ scaleSPU(v * ups()); }
inline void Synced::ups(double v){ scaleUPS(v * spu()); }
inline double Synced::spu() const { return mSPU; }
inline double Synced::ups() const { return mUPS; }

inline const Sync * Synced::sync() const { return mSubject; }



inline bool Sync::hasBeenSet() const { return mHasBeenSet; }
inline double Sync::spu() const { return mSPU; }
inline double Sync::ups() const { return mUPS; }
} // gam::

#endif

