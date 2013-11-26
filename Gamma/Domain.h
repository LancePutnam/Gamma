#ifndef GAMMA_DOMAIN_H_INC
#define GAMMA_DOMAIN_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include "Gamma/Node.h"

namespace gam{

class Domain;


/// Domain with normalized sampling frequency and interval
class Domain1{
public:

	virtual ~Domain1(){}

	double spu() const {return 1.;}	///< Get local samples/unit
	double ups() const {return 1.;}	///< Get local units/sample
	const Domain1 * domain() const	///< Get pointer to my subject domain (myself)
		{return this;}

	bool hasBeenSet() const { return true; }

	///	Called when subject domain's samples/unit changes

	/// Any instance state that depends on the samples/unit ratio should be 
	/// updated here. The ratio of the new to the old samples/unit is passed in.
	virtual void onDomainChange(double ratioSPU){}
	
	void spu(double val){}			///< Set local samples/unit	
	void ups(double val){}			///< Set local units/sample

protected:
	void refreshDomain(){ onDomainChange(1); }
};



/// Domain observer

/// This observer will attempt to copy its local variables from the default
/// subject upon construction. If the default subject has not been constructed, 
/// it will use its default values of 1. 
/// This class has a reference to a subject and its own local scaling factors.
/// By default, the reference subject is Domain::master.
class DomainObserver : public Node2<DomainObserver> {
public:

	DomainObserver();
	
	/// Copy constructor
	
	/// If the argument has a subject, then attach this as an observer to the
	/// same subject. Otherwise, attach as observer of Domain::master().
	DomainObserver(const DomainObserver& rhs);

	virtual ~DomainObserver(){}

	double scaleSPU() const;		///< Get ratio of my SPU to my domain's SPU
	double spu() const;				///< Get local samples/unit
	double ups() const;				///< Get local units/sample
	const Domain * domain() const;	///< Get pointer to my subject domain

	///	Called when subject domain's samples/unit changes

	/// Any instance state that depends on the samples/unit ratio should be 
	/// updated here. The ratio of the new to the old samples/unit is passed in.
	virtual void onDomainChange(double ratioSPU){}

	void scaleSPU(double v);		///< Scales samples/unit by factor
	void scaleUPS(double v);		///< Scales units/sample by factor
	void spu(double v);				///< Set local samples/unit
	void domain(Domain& src);		///< Set domain subject
	void ups(double v);				///< Set local units/sample

	DomainObserver& operator= (const DomainObserver& rhs);

protected:
	void refreshDomain(); ///< To be called from the constructor(s) of derived classes

private:
	friend class Domain;
	DomainObserver(bool zeroLinks, Domain& src);

	Domain * mSubject;	// Pointer to my subject
	double mSPU;		// Local samples/unit
	double mUPS;		// Local units/sample
};





/// Domain subject
class Domain{
public:

	Domain();

	/// \param[in]	spu		samples/unit
	Domain(double spu);

	~Domain();

	Domain& operator<< (DomainObserver& obs);	///< Attach observer
	void notifyObservers(double r);		///< Notify observers (\sa DomainObserver::onDomainChange)
	void spu(double v);					///< Set samples/unit and notify observers
	void ups(double v);					///< Set units/sample and notify observers

	bool hasBeenSet() const;			///< Returns true if spu has been set at least once
	double spu() const;					///< Returns samples/unit, i.e. sample rate
	double ups() const;					///< Returns units/sample, i.e. sample interval

	/// Master domain. By default, all observers will be attached to this.
	static Domain& master(){
		static Domain * s = new Domain;
		return *s;
	}

protected:
	double mSPU, mUPS;
	DomainObserver mHeadObserver;		// Head of observer doubly-linked list
	bool mHasBeenSet;

friend class DomainObserver;
	void attach(DomainObserver& obs);
};


/// Set master sample rate
void sampleRate(double samplesPerSecond);

/// Get master sample rate
double sampleRate();



// Implementation_______________________________________________________________

#define DOM_OBS_INIT mSubject(0), mSPU(1), mUPS(1)
inline DomainObserver::DomainObserver()
:	Node2<DomainObserver>(), DOM_OBS_INIT
{
	//printf("Synced::Synced() - %p\n", this);
	domain(Domain::master());
}

// This ctor is used to create the head observer in the subject
inline DomainObserver::DomainObserver(bool zeroLinks, Domain& src)
:	Node2<DomainObserver>(zeroLinks), DOM_OBS_INIT
{
	domain(src);
}

inline DomainObserver::DomainObserver(const DomainObserver& rhs)
:	Node2<DomainObserver>(), DOM_OBS_INIT
{
	//printf("DomainObserver::DomainObserver(const DomainObserver&) - %p\n", this);
	Domain& s = rhs.mSubject ? *rhs.mSubject : Domain::master();
	domain(s);
}
#undef DOM_OBS_INIT

inline void DomainObserver::refreshDomain(){ if(domain()) spu(domain()->spu()); }
inline double DomainObserver::scaleSPU() const { return domain() ? spu() / domain()->spu() : 1; }
inline void DomainObserver::spu(double v){ scaleSPU(v * ups()); }
inline void DomainObserver::ups(double v){ scaleUPS(v * spu()); }
inline double DomainObserver::spu() const { return mSPU; }
inline double DomainObserver::ups() const { return mUPS; }

inline const Domain * DomainObserver::domain() const { return mSubject; }



inline bool Domain::hasBeenSet() const { return mHasBeenSet; }
inline double Domain::spu() const { return mSPU; }
inline double Domain::ups() const { return mUPS; }
} // gam::

#include "Gamma/Sync.h"

#endif
