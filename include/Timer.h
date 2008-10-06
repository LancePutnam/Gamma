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

/// Timer
class Timer {
public:
	Timer();
	~Timer();

	nsec_t elapsed();					///< Returns nsec between start() and stop() calls
	double elapsedSec();				///< Returns  sec between start() and stop() calls
	double elapsedMSec();				///< Returns msec between start() and stop() calls
	void start();						///< Set start time as current time
	void stop();						///< Set stop time as current time

	static nsec_t currentTime();		///< Get current time from OS
	static void sleep(nsec_t dt);		///< Suspend thread execution for dt nsec
	static void sleepSec(double dt);	///< Suspend thread execution for dt sec
	static nsec_t sleepUntil(nsec_t t);	///< Suspend thread execution until absolute time, t. Returns ns slept.

	static double sec(nsec_t nsec);		///< Convert nsec to sec
	static nsec_t nsec(double sec);		///< Convert sec to nsec

private:
	nsec_t mStart, mStop;	// start and stop times
};



// Implementation_______________________________________________________________

//inline double Timer::sec(nsec_t nsec){ return (double)nsec * 0.000000001; }
inline double Timer::sec(nsec_t nsec){ return (double)nsec * 1e-9; }
inline nsec_t Timer::nsec(double sec){ return (nsec_t)(sec * 1e9); }

inline void Timer::start(){ mStart = currentTime(); }
inline void Timer::stop (){ mStop  = currentTime(); }
inline nsec_t Timer::elapsed(){ return mStop - mStart; }
inline double Timer::elapsedSec() { return sec(elapsed()); }
inline double Timer::elapsedMSec(){ return ((double)elapsed() * 0.000001); }

inline void Timer::sleepSec(double t){ Timer::sleep(nsec(t)); }

inline nsec_t Timer::sleepUntil(nsec_t t){
	nsec_t now = Timer::currentTime();
	if(t > now){
		Timer::sleep(t - now);
	}
	return t - now;
}


// platform specific

// Windows
#ifdef WIN32

	inline Timer::Timer(){
		timeBeginPeriod(1);		// adjust timer precision (1 millisecond)
	}

	inline Timer::~Timer(){
		timeEndPeriod(1);		// must call on program exit!
	}

	inline nsec_t Timer::currentTime(){
		return (nsec_t)timeGetTime() * (nsec_t)1000000;
	}

	inline void Timer::sleep(nsec_t dt){
		Sleep((DWORD)(dt / (nsec_t)1000000));
	}

// Posix (Mac, Linux)
#else

	inline Timer::Timer(){}
	inline Timer::~Timer(){}

	#define NS_S (nsec_t)1e9

	inline nsec_t Timer::currentTime(){
		timeval t;
		gettimeofday(&t, NULL);	
		return ((nsec_t)t.tv_sec) * NS_S + (nsec_t)(t.tv_usec * 1000);
	}

	inline void Timer::sleep(nsec_t t){
		time_t sec = (time_t)(t / NS_S);
		timespec tspec = { sec, (long)(t - ((nsec_t)sec * NS_S)) }; // { sec, nsec }
		nanosleep(&tspec, NULL);
	}

	#undef NS_S
	
#endif	/* platform specific */

} // end namespace gam

#endif

