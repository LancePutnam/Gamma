#ifndef GAMMA_THREAD_H_INC
#define GAMMA_THREAD_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

/***************************************************/
/* \class Thread
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

//#define GAM_USE_PTHREAD		(defined (__APPLE__) || defined (OSX) || defined (__LINUX__) || defined (__UNIX__))
//#define GAM_USE_THREADEX	(defined(WIN32) || defined(_WIN32) || defined(WIN64))

#include "Gamma/Config.h"

#define GAM_USE_PTHREAD		(GAM_OSX || GAM_LINUX)
#define GAM_USE_THREADEX	(GAM_WINDOWS)

#if GAM_USE_PTHREAD
	#include <pthread.h>
#elif GAM_USE_THREADEX
	#include <windows.h>
	#include <process.h>
#endif

namespace gam{

class Thread{
public:

	typedef void * (*Function)(void * user);

	#if GAM_USE_PTHREAD
		typedef pthread_t		Handle;
	#elif GAM_USE_THREADEX
		typedef unsigned long	Handle;
	#endif

	Thread()
	:	mHandle(0), mJoinOnDestroy(false)
	{}
	
	Thread(Function func, void * user = NULL)
	:	mHandle(0), mJoinOnDestroy(false)
	{	start(func, user); }

	~Thread(){
		if(mJoinOnDestroy) join();
	}

	/// Begin execution of the thread routine.  Upon success, true is returned.

	/// A data pointer can be supplied to the thread routine via the
	/// optional \e ptr argument.  If the thread cannot be created, the
	/// return value is false.
	bool start(Function func, void * user = NULL);

	/// Signal cancellation of a thread routine, returning true on success.
	
	/// This function only signals thread cancellation.  It does not
	/// wait to verify actual routine termination.  A true return value
	/// only signifies that the cancellation signal was properly executed,
	/// not thread cancellation.  A thread routine may need to make use of
	/// the testCancel() function to specify a cancellation point.
	bool cancel();

	/// Block the calling routine indefinitely until the thread terminates.
	
	/// This function suspends execution of the calling routine until the thread 
	/// has terminated.  It will return immediately if the thread was already 
	/// terminated.  A \e true return value signifies successful termination. 
	/// A false return value indicates a problem with the call.
	bool join();


	/// Set whether thread will automatically join upon destruction
	Thread& joinOnDestroy(bool v){ mJoinOnDestroy=v; return *this; }

protected:
	Handle mHandle;
	bool mJoinOnDestroy;
};




// Implementation

#if GAM_USE_PTHREAD

inline bool Thread::start(Thread::Function func, void * user){
	if(mHandle) return false;
	return 0 == pthread_create(&mHandle, NULL, *func, user);
}

inline bool Thread::cancel(){
	return 0 == pthread_cancel(mHandle);
}

inline bool Thread::join(){
	if(pthread_join(mHandle, NULL) == 0){
		mHandle = 0;
		return true;
	}
	return false;
}


#elif GAM_USE_THREADEX

inline bool Thread::start(Thread::Function func, void * user){
	if(mHandle) return false;
	
//	struct F{
//		Thread::Function func;
//		void * user;
//
//		static unsigned _stdcall * call(void * user){
//			F& f = *(F*)user;
//			f.func(f.userData);
//			return 0;
//		}
//	} f = { func, user };
//	
//	unsigned thread_id;
//	mHandle = _beginthreadex(NULL, 0, F::call, &f, 0, &thread_id);
//	if(mHandle) return true;
//	return false;


	struct F{
		Thread::Function func;
		void * userData;

		static unsigned _stdcall call(void * user){
			F *pF = reinterpret_cast<F*>(user);
			(*(pF->func))(pF->userData);
			delete pF;
			return 0;
		}
	};

	struct F* f = new F;
	f->func = func;
	f->userData = user;

	unsigned thread_id;
	mHandle = _beginthreadex(NULL, 0, F::call, f, 0, &thread_id);
	if(mHandle) return true;
	return false;
}

inline bool Thread::cancel(){
	TerminateThread((HANDLE)mHandle, 0);
	return true;
}

inline bool Thread::join(){
	long retval = WaitForSingleObject((HANDLE)mHandle, INFINITE);
	if(WAIT_OBJECT_0 == retval){
		CloseHandle((HANDLE)mHandle);
		mHandle = 0;
		return true;
	}
	return false;
}

#endif

} // gam::

#endif

