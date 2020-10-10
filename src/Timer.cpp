/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include "Gamma/Timer.h"
#include "Gamma/Config.h"

// Windows
#if GAM_WINDOWS
	#define WIN32_LEAN_AND_MEAN
	#include <windows.h>
	
// Posix (Mac, Linux)
#else
	#include <sys/time.h>
	#include <time.h>
	
#endif /* platform specific */

namespace gam{

// Windows
#ifdef GAM_WINDOWS

	nsec_t steady_time_us(){
		// Windows 10 and above
		//PULONGLONG time;
		//QueryInterruptTimePrecise(&time);

		LARGE_INTEGER freq; // ticks/second
		LARGE_INTEGER time; // tick count
		QueryPerformanceFrequency(&freq);
		QueryPerformanceCounter(&time);
		// convert ticks to microseconds
		// As long as time.QuadPart < 9.0e15, this will work:
		return nsec_t(time.QuadPart / double(freq.QuadPart) * 1.0e6);
	}

	nsec_t timeNow(){
		return steady_time_us() * nsec_t(1000);
	}

	void sleep(nsec_t dt){
		Sleep((DWORD)(dt / (nsec_t)1e6));
	}

// Posix (Mac, Linux)
#else

	#define NS_S ((nsec_t)1e9)

	nsec_t timeNow(){
		timeval t;
		gettimeofday(&t, NULL);	
		return ((nsec_t)t.tv_sec) * NS_S + (nsec_t)(t.tv_usec * 1e3);
	}

	void sleep(nsec_t t){
		time_t sec = (time_t)(t / NS_S);
		timespec tspec = { sec, (long)(t - ((nsec_t)sec * NS_S)) }; // { sec, nsec }
		nanosleep(&tspec, NULL);
	}

	#undef NS_S
	
#endif	/* platform specific */

void sleepSec(double t){
	sleep(toNSec(t));
}

nsec_t sleepUntil(nsec_t t){
	nsec_t now = timeNow();
	if(t > now) sleep(t - now);
	return t - now;
}

void Timer::start(){ mStart = timeNow(); }
void Timer::stop (){ mStop  = timeNow(); }
nsec_t Timer::elapsed(){ return mStop - mStart; }
double Timer::elapsedSec() { return toSec(elapsed()); }
double Timer::elapsedMSec(){ return ((double)elapsed() * 1e-6); }

} // gam::
