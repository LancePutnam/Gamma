#ifndef GAMMA_THREAD_H_INC
#define GAMMA_THREAD_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include <thread>

namespace gam{

class Thread{
public:

	typedef void * (*Function)(void * user);

	Thread(){}

	Thread(Function func, void * user = NULL){
		start(func, user);
	}

	~Thread(){
		if(mJoinOnDestroy) join();
	}

	/// Begin execution of the thread routine.  Upon success, true is returned.

	/// A data pointer can be supplied to the thread routine via the
	/// optional \e ptr argument.  If the thread cannot be created, the
	/// return value is false.
	bool start(Function func, void * user = NULL){
		mHandle = std::thread(func, user);
		return true;
	}

	/// Block the calling routine indefinitely until the thread terminates.
	
	/// This function suspends execution of the calling routine until the thread 
	/// has terminated.  It will return immediately if the thread was already 
	/// terminated.  A \e true return value signifies successful termination. 
	/// A false return value indicates a problem with the call.
	bool join(){
		mHandle.join();
		return true;
	}


	/// Set whether thread will automatically join upon destruction
	Thread& joinOnDestroy(bool v){ mJoinOnDestroy=v; return *this; }

protected:
	std::thread mHandle;
	bool mJoinOnDestroy = false;
};

} // gam::

#endif
