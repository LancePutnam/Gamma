#ifndef GAMMA_TIMER_H_INC
#define GAMMA_TIMER_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

namespace gam{

typedef long long nsec_t; ///< Nanoseconds type; good up to +/- 292.5 years


/// Suspend thread execution for dt nsec
void sleep(nsec_t dt);

/// Suspend thread execution for a specified amount of seconds
void sleepSec(double sec);

/// Suspend thread execution until absolute time, t. Returns ns slept.
nsec_t sleepUntil(nsec_t t);

/// Get current time from OS
nsec_t timeNow();

/// Convert nsec to sec
double toSec(nsec_t nsec);

/// Convert sec to nsec
nsec_t toNSec(double sec);


/// Timer
class Timer {
public:

	nsec_t elapsed();		///< Returns nsec between start() and stop() calls
	double elapsedSec();	///< Returns  sec between start() and stop() calls
	double elapsedMSec();	///< Returns msec between start() and stop() calls
	void start();			///< Set start time as current time
	void stop();			///< Set stop time as current time

private:
	nsec_t mStart, mStop;	// start and stop times
};



// Implementation_______________________________________________________________

inline double toSec(nsec_t nsec){ return double(nsec) * 1e-9; }
inline nsec_t toNSec(double sec){ return nsec_t(sec * 1e9); }

} // gam::

#endif

