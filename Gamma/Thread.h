#ifndef GAMMA_THREAD_H_INC
#define GAMMA_THREAD_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include "Gamma/Config.h"

//#define GAM_USE_STD_THREAD	1
#define GAM_USE_PTHREAD		(GAM_OSX || GAM_LINUX)
#define GAM_USE_WINTHREAD	(GAM_WINDOWS)

#if GAM_USE_STD_THREAD
	#include <thread>
#elif GAM_USE_PTHREAD
	#include <pthread.h>
#elif GAM_USE_WINTHREAD
	#define WIN32_LEAN_AND_MEAN
	#include <windows.h>
	#ifdef far
	#undef far
	#endif
	#ifdef near
	#undef near
	#endif
#endif

namespace gam{

class Thread{
public:

	typedef void * (*Function)(void * user);

	#if GAM_USE_STD_THREAD
		typedef std::thread	Handle;
	#elif GAM_USE_PTHREAD
		typedef pthread_t	Handle;
	#elif GAM_USE_WINTHREAD
		typedef HANDLE		Handle;
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

#if GAM_USE_STD_THREAD

inline bool Thread::start(Thread::Function func, void * user){
	mHandle = std::thread(func, user);
}

inline bool Thread::join(){
	mHandle.join();
	return true;
}

#elif GAM_USE_PTHREAD

inline bool Thread::start(Thread::Function func, void * user){
	if(mHandle) return false;
	return 0 == pthread_create(&mHandle, NULL, *func, user);
}

inline bool Thread::join(){
	if(pthread_join(mHandle, NULL) == 0){
		mHandle = 0;
		return true;
	}
	return false;
}


#elif GAM_USE_WINTHREAD

namespace{
struct ThreadFunctor{
	Thread::Function func;
	void * userData;

	static DWORD WINAPI call(void * user){
		ThreadFunctor * pF = reinterpret_cast<ThreadFunctor*>(user);
		(*(pF->func))(pF->userData);
		delete pF;
		return 0;
	}
};
}

inline bool Thread::start(Thread::Function func, void * user){
	if(mHandle) return false;

	struct ThreadFunctor* f = new ThreadFunctor;
	f->func = func;
	f->userData = user;

	DWORD thread_id;
	// _beginthreadex should be used if the C run-time library is used in
	// the thread function. However, it is not available in MSYS2...
	//mHandle = _beginthreadex(NULL, 0, ThreadFunctor::call, f, 0, &thread_id);
	mHandle = CreateThread(NULL, 0, ThreadFunctor::call, f, 0, &thread_id);
	if(mHandle) return true;
	return false;
}

inline bool Thread::join(){
	long retval = WaitForSingleObject(mHandle, INFINITE);
	if(WAIT_OBJECT_0 == retval){
		CloseHandle(mHandle);
		mHandle = 0;
		return true;
	}
	return false;
}

#endif

} // gam::

#endif

