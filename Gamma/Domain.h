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

	double spu() const {return 1.;}	///< Get samples/unit
	double ups() const {return 1.;}	///< Get units/sample
	const Domain1 * domain() const	///< Get pointer to my subject domain (myself)
		{return this;}

	bool hasBeenSet() const { return true; }

	///	Called when subject domain's samples/unit changes

	/// Any instance state that depends on the samples/unit ratio should be 
	/// updated here. The ratio of the new to the old samples/unit is passed in.
	virtual void onDomainChange(double ratioSPU){}
	
	void spu(double val){}			///< Set samples/unit
	void ups(double val){}			///< Set units/sample

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

	virtual ~DomainObserver();

	double spu() const;				///< Get samples/unit
	double ups() const;				///< Get units/sample
	const Domain * domain() const;	///< Get pointer to my subject domain

	///	Called when subject domain's samples/unit changes

	/// Any instance state that depends on the samples/unit ratio should be 
	/// updated here. The ratio of the new to the old samples/unit is passed in.
	virtual void onDomainChange(double ratioSPU){}

	/// Set domain subject
	
	/// If the object was already attached to another Domain, it will be
	/// detached.
	void domain(Domain& src);

	DomainObserver& operator= (const DomainObserver& rhs);

protected:
	/// Forces call to onDomainChange
	
	/// This should be called from the constructor(s) of derived classes
	/// since base classes cannot correctly call a virtual function in their
	/// constructor.
	void refreshDomain();

private:
	friend class Domain;

	Domain * mSubject;	// Pointer to my subject
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
	static Domain& master();

protected:
	double mSPU, mUPS;
	DomainObserver * mHeadObserver;	// Head of observer doubly-linked list
	bool mHasBeenSet;

friend class DomainObserver;
	void attach(DomainObserver& obs);
};


/// Set master sample rate
void sampleRate(double samplesPerSecond);

/// Get master sample rate
double sampleRate();



// Implementation_______________________________________________________________

inline void DomainObserver::refreshDomain(){ onDomainChange(1); }
inline double DomainObserver::spu() const { return domain()->spu(); }
inline double DomainObserver::ups() const { return domain()->ups(); }
inline const Domain * DomainObserver::domain() const { return mSubject; }


inline bool Domain::hasBeenSet() const { return mHasBeenSet; }
inline double Domain::spu() const { return mSPU; }
inline double Domain::ups() const { return mUPS; }
} // gam::

#include "Gamma/Sync.h"

#endif
