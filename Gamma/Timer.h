#ifndef GAMMA_TIMER_H_INC
#define GAMMA_TIMER_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

namespace gam{

// We use an 8-byte signed integer to represent time in nanoseconds.  This
// allows us represent times accurately up to +/- 292.5 years.

// Windows
#ifdef WIN32
	#include <windows.h>
	typedef __int64 nsec_t;		///< nanoseconds type
	
// Posix (Mac, Linux)
#else
	#include <sys/time.h>
	#include <time.h>
	typedef long long nsec_t;	///< nanoseconds type
	
#endif /* platform specific */


/// Suspend thread execution for dt nsec
static void sleep(nsec_t dt);

/// Suspend thread execution for a specified amount of seconds
static void sleepSec(double sec);

/// Suspend thread execution until absolute time, t. Returns ns slept.
static nsec_t sleepUntil(nsec_t t);

/// Get current time from OS
static nsec_t timeNow();

/// Convert nsec to sec
static double toSec(nsec_t nsec);

/// Convert sec to nsec
static nsec_t toNsec(double sec);


/// Timer
class Timer {
public:
	Timer();
	~Timer();

	nsec_t elapsed();		///< Returns nsec between start() and stop() calls
	double elapsedSec();	///< Returns  sec between start() and stop() calls
	double elapsedMSec();	///< Returns msec between start() and stop() calls
	void start();			///< Set start time as current time
	void stop();			///< Set stop time as current time

private:
	nsec_t mStart, mStop;	// start and stop times
};



// Implementation_______________________________________________________________

inline void sleepSec(double t){ sleep(toNsec(t)); }

inline nsec_t sleepUntil(nsec_t t){
	nsec_t now = timeNow();
	if(t > now) sleep(t - now);
	return t - now;
}

inline double toSec(nsec_t nsec){ return ((double)nsec) * 1e-9; }
inline nsec_t toNsec(double sec){ return (nsec_t)(sec * 1e9); }

inline void Timer::start(){ mStart = timeNow(); }
inline void Timer::stop (){ mStop  = timeNow(); }
inline nsec_t Timer::elapsed(){ return mStop - mStart; }
inline double Timer::elapsedSec() { return toSec(elapsed()); }
inline double Timer::elapsedMSec(){ return ((double)elapsed() * 1e-6); }


// platform specific

// Windows
#ifdef WIN32

	inline nsec_t timeNow(){
		return (nsec_t)timeGetTime() * (nsec_t)1e6;
	}

	inline void sleep(nsec_t dt){
		Sleep((DWORD)(dt / (nsec_t)1e6));
	}

	inline Timer::Timer(){
		timeBeginPeriod(1);		// adjust timer precision (1 millisecond)
	}

	inline Timer::~Timer(){
		timeEndPeriod(1);		// must call on program exit!
	}

// Posix (Mac, Linux)
#else

	#define NS_S ((nsec_t)1e9)

	inline nsec_t timeNow(){
		timeval t;
		gettimeofday(&t, NULL);	
		return ((nsec_t)t.tv_sec) * NS_S + (nsec_t)(t.tv_usec * 1e3);
	}

	inline void sleep(nsec_t t){
		time_t sec = (time_t)(t / NS_S);
		timespec tspec = { sec, (long)(t - ((nsec_t)sec * NS_S)) }; // { sec, nsec }
		nanosleep(&tspec, NULL);
	}

	#undef NS_S

	inline Timer::Timer(){}
	inline Timer::~Timer(){}
	
#endif	/* platform specific */

} // gam::

#endif

