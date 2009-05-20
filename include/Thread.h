#ifndef GAMMA_TIMER_H_INC
#define GAMMA_TIMER_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

/***************************************************/
/*! \class Thread
    \brief STK thread class.

    This class provides a uniform interface for cross-platform
    threads.  On unix systems, the pthread library is used.  Under
    Windows, the C runtime threadex functions are used.

    Each instance of the Thread class can be used to control a single
    thread process.  Routines are provided to signal cancelation
    and/or joining with a thread, though it is not possible for this
    class to know the running status of a thread once it is started.

    For cross-platform compatability, thread functions should be
    declared as follows:

    THREAD_FUNCTION(thread_function_name)

    by Perry R. Cook and Gary P. Scavone, 1995 - 2005.
*/
/***************************************************/

namespace gam{

#define USE_PTHREAD		(defined (__APPLE__) || defined (OSX) || defined (__LINUX__) || defined (__UNIX__))
#define USE_THREADEX	(defined(WIN32))

#if USE_PTHREAD
	#include <pthread.h>
	typedef pthread_t ThreadHandle;
	typedef void * (*ThreadFunction)(void *);
	#define THREAD_FUNCTION(name) void * name(void * user)

#elif USE_THREADEX
	#include <windows.h>
	#include <process.h>
	typedef unsigned long ThreadHandle;
	typedef unsigned (__stdcall *ThreadFunction)(void *);
	#define THREAD_FUNCTION(name) unsigned _stdcall * name(void * user)

#endif



class Thread{
public:
	//! Default constructor.
	Thread();
	
	Thread(ThreadFunction routine, void * ptr = NULL);

	//! The class destructor does not attempt to cancel or join a thread.
	~Thread();

	//! Begin execution of the thread \e routine.  Upon success, true is returned.
	/*!
	A data pointer can be supplied to the thread routine via the
	optional \e ptr argument.  If the thread cannot be created, the
	return value is false.
	*/
	bool start(ThreadFunction routine, void * ptr = NULL);

	//! Signal cancellation of a thread routine, returning \e true on success.
	/*!
	This function only signals thread cancellation.  It does not
	wait to verify actual routine termination.  A \e true return value
	only signifies that the cancellation signal was properly executed,
	not thread cancellation.  A thread routine may need to make use of
	the testCancel() function to specify a cancellation point.
	*/
	bool cancel();

	//! Block the calling routine indefinitely until the thread terminates.
	/*!
	This function suspends execution of the calling routine until the thread has 
	terminated.  It will return immediately if the thread was already 
	terminated.  A \e true return value signifies successful termination. 
	A \e false return value indicates a problem with the wait call.
	*/
	bool wait();

	//! Create a cancellation point within a thread routine.
	/*!
	This function call checks for thread cancellation, allowing the
	thread to be terminated if a cancellation request was previously
	signaled.
	*/
	void testCancel();

protected:
	ThreadHandle mHandle;
};




// Implementation

inline Thread::Thread() : mHandle(0){}

inline Thread::Thread(ThreadFunction routine, void * ptr) :
	mHandle(0)
{
	start(routine, ptr);
}

inline Thread::~Thread(){}

#if USE_PTHREAD

inline bool Thread::start( ThreadFunction routine, void * ptr ){
	if ( mHandle )	return false;
	return 0 == pthread_create(&mHandle, NULL, *routine, ptr);
}

inline bool Thread::cancel(){
	return 0 == pthread_cancel(mHandle);
}

inline bool Thread::wait(){
	if ( pthread_join(mHandle, NULL) == 0 ) {
		mHandle = 0;
		return true;
	}
	return false;
}

inline void Thread::testCancel(void){
	pthread_testcancel();
}


#elif USE_THREADEX

inline bool Thread::start( ThreadFunction routine, void * ptr ){
	if ( mHandle ) return false;
	unsigned thread_id;
	mHandle = _beginthreadex(NULL, 0, routine, ptr, 0, &thread_id);
	if ( mHandle ) return true;
	return false;
}

inline bool Thread::cancel(){
	TerminateThread((HANDLE)mHandle, 0);
	return true;
}

inline bool Thread::wait(){
	long retval = WaitForSingleObject( (HANDLE)mHandle, INFINITE );
	if ( retval == WAIT_OBJECT_0 ) {
		CloseHandle( (HANDLE)mHandle );
		mHandle = 0;
		return true;
	}
	return false;
}

void Thread::testCancel(void){
}

#endif

} // gam::

#endif

