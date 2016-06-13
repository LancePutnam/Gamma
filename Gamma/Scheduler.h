#ifndef INC_GAM_SCHEDULER_H
#define INC_GAM_SCHEDULER_H

#include <cstdlib> // exit
#include <cstring> // memcpy, size_t
#include <queue>
#include <list>

#include "Gamma/Node.h"
#include "Gamma/Print.h"
#include "Gamma/Thread.h"

namespace gam{

class Scheduler;

#ifndef GAM_FUNC_MAX_DATA_SIZE
	#define GAM_FUNC_MAX_DATA_SIZE 64
#endif


/// Deferrable function
class Func{
public:
	typedef void (* func_t)(void * data);

	template <class R>
	Func(R (*fnc)()){
		struct Data{
			R (*fnc)();
			static void call(void * data){
				Data& d = *(Data *)data;
				d.fnc();
			}
		} data = {fnc};
		set(Data::call, (void *)&data, sizeof(Data));
	}

	template <class R, class A, class L>
	Func(R (*fnc)(A), L l){
		struct Data{
			R (*fnc)(A);
			L l;
			static void call(void * data){
				Data& d = *(Data *)data;
				d.fnc(d.l);
			}
		} data = {fnc,l};
		set(Data::call, (void *)&data, sizeof(Data));
	}

	template <class R, class A, class B, class L, class M>
	Func(R (*fnc)(A,B), L l, M m){
		struct Data{
			R (*fnc)(A,B);
			L l; M m;
			static void call(void * data){
				Data& d = *(Data *)data;
				d.fnc(d.l,d.m);
			}
		} data = {fnc,l,m};
		set(Data::call, (void *)&data, sizeof(Data));
	}

	template <class R, class A, class B, class C, class L, class M, class N>
	Func(R (*fnc)(A,B,C), L l, M m, N n){
		struct Data{
			R (*fnc)(A,B,C);
			L l; M m; N n;
			static void call(void * data){
				Data& d = *(Data *)data;
				d.fnc(d.l,d.m,d.n);
			}
		} data = {fnc,l,m,n};
		set(Data::call, (void *)&data, sizeof(Data));
	}

	template <class R, class A, class B, class C, class D, class L, class M, class N, class O>
	Func(R (*fnc)(A,B,C,D), L l, M m, N n, O o){
		struct Data{
			R (*fnc)(A,B,C,D);
			L l; M m; N n; O o;
			static void call(void * data){
				Data& d = *(Data *)data;
				d.fnc(d.l,d.m,d.n,d.o);
			}
		} data = {fnc,l,m,n,o};
		set(Data::call, (void *)&data, sizeof(Data));
	}


	template <class Obj1, class Obj2, class R>
	Func(Obj1& obj, R (Obj2::*mth)()){
		struct Data{
			Obj1& obj;
			R (Obj2::*mth)();
			static void call(void * data){
				Data& d = *(Data *)data;
				(d.obj.*d.mth)();
			}
		} data = {obj,mth};
		set(Data::call, (void *)&data, sizeof(Data));
	}

	template <class Obj1, class Obj2, class R, class A, class L>
	Func(Obj1& obj, R (Obj2::*mth)(A), L l){
		struct Data{
			Obj1& obj;
			R (Obj2::*mth)(A);
			L l;
			static void call(void * data){
				Data& d = *(Data *)data;
				(d.obj.*d.mth)(d.l);
			}
		} data = {obj,mth,l};
		set(Data::call, (void *)&data, sizeof(Data));
	}

	template <class Obj1, class Obj2, class R, class A, class B, class L, class M>
	Func(Obj1& obj, R (Obj2::*mth)(A,B), L l, M m){
		struct Data{
			Obj1& obj;
			R (Obj2::*mth)(A,B);
			L l; M m;
			static void call(void * data){
				Data& d = *(Data *)data;
				(d.obj.*d.mth)(d.l,d.m);
			}
		} data = {obj,mth,l,m};
		set(Data::call, (void *)&data, sizeof(Data));
	}

	template <class Obj1, class Obj2, class R, class A, class B, class C, class L, class M, class N>
	Func(Obj1& obj, R (Obj2::*mth)(A,B,C), L l, M m, N n){
		struct Data{
			Obj1& obj;
			R (Obj2::*mth)(A,B,C);
			L l; M m; N n;
			static void call(void * data){
				Data& d = *(Data *)data;
				(d.obj.*d.mth)(d.l,d.m,d.n);
			}
		} data = {obj,mth,l,m,n};
		set(Data::call, (void *)&data, sizeof(Data));
	}

	template <class Obj1, class Obj2, class R, class A, class B, class C, class D, class L, class M, class N, class O>
	Func(Obj1& obj, R (Obj2::*mth)(A,B,C,D), L l, M m, N n, O o){
		struct Data{
			Obj1& obj;
			R (Obj2::*mth)(A,B,C,D);
			L l; M m; N n; O o;
			static void call(void * data){
				Data& d = *(Data *)data;
				(d.obj.*d.mth)(d.l,d.m,d.n,d.o);
			}
		} data = {obj,mth,l,m,n,o};
		set(Data::call, (void *)&data, sizeof(Data));
	}

	/// Execute stored function
	void operator()(){ mFunc(mData); }

	const char * data() const { return mData; }
	const void * obj() const { return mObj; }

private:
	union{
		char mData[GAM_FUNC_MAX_DATA_SIZE];
		void * mObj;
	};
	func_t mFunc;

	void set(func_t f, void * data, int size){
		mFunc = f;
		#ifndef NDEBUG
		if(size > GAM_FUNC_MAX_DATA_SIZE){
			fprintf(stderr, "Func maximum data size exceeded. "
				"Attempt to use %d bytes with maximum size set to %d.\n",
				size, GAM_FUNC_MAX_DATA_SIZE);
			std::exit(-1);
		}
		#endif
		std::memcpy(mData, data, size);
		//printf("Func::set() with %d bytes\n", size);		
	}
};



/// Audio I/O data structure used by real-time scheduling system
struct SchedulerAudioIOData{

	SchedulerAudioIOData()
	:	buffersIn(NULL), buffersOut(NULL),
		framesPerSecond(1), framesPerBuffer(0), channelsIn(0), channelsOut(0),
		startFrame(0),
		mUserData(NULL), mUserDataTypeID(0)
	{}


	const float * buffersIn;	///< Non-interleaved input buffers
	float * buffersOut;			///< Non-interleaved output buffers
	double framesPerSecond;		///< Frames per second (frame rate)
	unsigned framesPerBuffer;	///< Frames per buffer (block size)
	unsigned channelsIn;		///< Number of input channels
	unsigned channelsOut;		///< Number of output channels
	unsigned startFrame;		///< Start frame to begin processing


	/// Set user data
	template <class T>
	void userData(T * user){
		mUserData = user;
		mUserDataTypeID = typeID<T>();
	}

	/// Get user data
	void * userData() const { return mUserData; }

	/// Get number of seconds in one buffer
	double secondsPerBuffer() const{
		return framesPerBuffer / framesPerSecond; }


	/// Map external to internal audio I/O data

	/// Use this method to map an external data structure holding audio I/O
	/// data (buffers, frame rate, block size, number of channels, etc.). This method
	/// exists primarily to make the Scheduler easy to use in conjunction with
	/// gam::AudioIOData.
	///
	/// \tparam TAudioIOData	A class with the same interface as gam::AudioIOData
	/// \param[in] externalIO	External audio I/O data
	template <class TAudioIOData>
	void mapAudioIOData(TAudioIOData& externalIO);

	/// Unmap audio I/O data (must be matched with a call to mapAudioIOData)

	/// \tparam TAudioIOData	A class with the same interface as gam::AudioIOData
	///
	template <class TAudioIOData>
	TAudioIOData& unmapAudioIOData();

private:
	template <class T>
	static std::size_t typeID(){
		static int x;
		return std::size_t(&x);
	}

	void * mUserData;				// User data (usually other audio I/O data)
	std::size_t mUserDataTypeID;	// Use for safe casting
};



// A block-rate processing node in the audio graph
class ProcessNode : public Node3<ProcessNode>{
public:

	ProcessNode(double delay=0.);

	virtual ~ProcessNode();


	/// Called whenever this node must process audio
	virtual void onProcessNode(SchedulerAudioIOData& io){}

	/// Called whenever this node is "reset"
	virtual void onReset(){}


	/// Set starting time offset, in seconds
	ProcessNode& dt(double v){ mDelay=v; return *this; }

	/// Flag self (and consequently all descendents) for deletion
	ProcessNode& free();
	
	/// Set whether processor is active
	
	/// If true, then the processor is executed in the synthesis network.
	/// If false, the processor and its descendents are skipped.
	ProcessNode& active(bool v);

	ProcessNode& reset();

	bool deletable() const { return mDeletable; }
	bool done() const { return DONE==mStatus; }
	bool active() const { return ACTIVE==mStatus; }
	bool inactive() const { return INACTIVE==mStatus; }

	void print();

protected:
	friend class Scheduler;
	
	enum{
		INACTIVE=0,		// node and descendents are not executed
		ACTIVE,			// node and descendents are executed
		DONE			// processing done, node and descendents can be removed
	};
	
	int mStatus;
	double mDelay;
	bool mDeletable;

	// Call my own processing algorithm, onProcess()
	ProcessNode * update(const ProcessNode * top, SchedulerAudioIOData& io);

	// Process all descendents recursively
	ProcessNode * process(const ProcessNode * top, SchedulerAudioIOData& io, int frameStart=0);
};



/// ProcessNode with callback using a gam::AudioIOData-like interface
template <class TAudioIOData>
class Process : public ProcessNode{
public:

	virtual void onProcess(TAudioIOData& io) = 0;

private:
	void onProcessNode(SchedulerAudioIOData& io){
		onProcess(io.unmapAudioIOData<TAudioIOData>());
	}
};


/// A function that can be delayed and/or repeated periodically
class ControlFunc{
public:
	ControlFunc(const Func& f, double dt=0)
	:	mFunc(f), mDelay(dt), mPeriod(0), mObjDel(0)
	{}
	
	ControlFunc& dt(double v){ mDelay=v; return *this; }
	ControlFunc& period(double v){ mPeriod=v; return *this; }

	void operator()(){ mFunc(); }

protected:
	friend class Scheduler;
	Func mFunc;
	double mDelay;
	double mPeriod;
	int mObjDel;
};



/// Schedules real-time audio processes

/// Before starting the scheduler, you must map your application's audio buffers 
/// and other information to the scheduler's SchedulerAudioIOData (accessed with 
/// the io() method).
class Scheduler : public ProcessNode{
public:

	typedef std::queue<ProcessNode *> FreeList;
	typedef std::list<ControlFunc> Funcs;

	Scheduler();
	~Scheduler();


	/// Test whether the synthesis graph is empty
	bool empty() const;
	
	/// Check free list for finished events and reclaim their memory
	
	/// \returns number of events deleted
	///
	int reclaim();

	/// Get internal audio I/O data structure
	const SchedulerAudioIOData& io() const { return mIO; }
	SchedulerAudioIOData& io(){ return mIO; }
	

	/// Add dynamically allocated process as first child of root node
	template <class AProcess>
	AProcess& add(){
		AProcess * v = new AProcess;
		cmdAdd(v); return *v;
	}

	template <class AProcess, class A>
	AProcess& add(const A& a){
		AProcess * v = new AProcess(a);
		cmdAdd(v); return *v;
	}

	template <class AProcess, class A, class B>
	AProcess& add(const A& a, const B& b){
		AProcess * v = new AProcess(a,b);
		cmdAdd(v); return *v;
	}

	template <class AProcess, class A, class B, class C>
	AProcess& add(const A& a, const B& b, const C& c){
		AProcess * v = new AProcess(a,b,c);
		cmdAdd(v); return *v;
	}

	template <class AProcess, class A, class B, class C, class D>
	AProcess& add(const A& a, const B& b, const C& c, const D& d){
		AProcess * v = new AProcess(a,b,c,d);
		cmdAdd(v); return *v;
	}

	template <class AProcess, class A, class B, class C, class D, class E>
	AProcess& add(const A& a, const B& b, const C& c, const D& d, const E& e){
		AProcess * v = new AProcess(a,b,c,d,e);
		cmdAdd(v); return *v;
	}

	template <class AProcess, class A, class B, class C, class D, class E, class F>
	AProcess& add(const A& a, const B& b, const C& c, const D& d, const E& e, const F& f){
		AProcess * v = new AProcess(a,b,c,d,e,f);
		cmdAdd(v); return *v;
	}


	/// Add dynamically allocated process as first child of specified node
	template <class AProcess>
	AProcess& add(ProcessNode& parent){
		AProcess * v = new AProcess;
		pushCommand(Command::ADD_FIRST_CHILD, &parent, v);
		return *v;
	}


	/// Add deferred function call
	ControlFunc& add(const Func& f){
		mFuncs.push_back(f);
		return mFuncs.back();
	}

	/// Execute all audio processes in execution tree

	/// This should be called at the audio block rate. 
	/// The latency of events will be determined by the block size.
	void update();

	/// Map external audio I/O and then update

	/// \tparam TAudioIOData	A class with the same interface as gam::AudioIOData
	///
	template <class TAudioIOData>
	void update(TAudioIOData& externalIO){
		io().mapAudioIOData(externalIO);
		update();
	}

	/// Set time period between low-priority actions
	Scheduler& period(float v);

	/// Start scheduler
	void start();

	/// Stop scheduler
	void stop();

	/// Record output to sound file in non-real-time

	/// \param[in] soundFilePath	path to sound file
	/// \param[in] durSec			duration, in seconds, of recording
	//void recordNRT(GAM_SCHEDULER_IO_DATA& io, const char * soundFilePath, double durSec);
	void recordNRT(const char * soundFilePath, double durSec);

	/// Record output to sound file in non-real-time

	/// \param[in] aio				audio i/o data to map to internal audio i/o data
	/// \param[in] soundFilePath	path to sound file
	/// \param[in] durSec			duration, in seconds, of recording
	template <class TAudioIOData>
	void recordNRT(TAudioIOData& aio, const char * soundFilePath, double durSec){
		io().mapAudioIOData(aio);
		recordNRT(soundFilePath, durSec);
	}

//	void print(){
//		printf("%d events\n", (int)events().size());
//		Events::iterator it = events().begin();
//		for(; it!=events().end(); ++it){
//			printf("\t"); (*it)->print();
//		}
//	}


	/// A static audio callback function that updates the Scheduler

	/// \tparam TAudioIOData	A class sharing the same interface as gam::AudioIOData
	///
	template <class TAudioIOData>
	static void audioCB(TAudioIOData& aio);

protected:

	struct Command{
		enum Type{
			ADD_FIRST_CHILD,
			ADD_LAST_CHILD,
			REMOVE_CHILD
		};
		
		//double time;
		Type type;
		ProcessNode * object;
		ProcessNode * other;

//		struct Compare{
//			bool operator()(const Command& a, const Command& b){
//				return a.time < b.time;
//			}
//		};
	};

//	std::priority_queue<Command, std::vector<Command>, Command::Compare> 
//		mCommandQueue;

	// LPT:  low-priority thread
	// HPT: high-priority thread
	std::queue<Command> mAddCommands;	// items newly allocated in LPT to be added to tree in HPT
	FreeList mFreeList;		// items removed from tree in HPT to be deleted in LPT
	Funcs mFuncs;
	Thread mLPThread;		// low-priority thread for garbage collection, etc.
	float mPeriod;
	double mTime;			// scheduler's time, in seconds
	SchedulerAudioIOData mIO;
	bool mRunning;
	
	static void * cLPThreadFunc(void * user);

	void pushCommand(Command::Type c, ProcessNode * object, ProcessNode * other);
	void cmdAdd(ProcessNode * v);

	// Execute pending graph manipulation commands from HPT.
	// This is to be called from the same thread that is executing the
	// processing within nodes (i.e., from the audio thread).
	void hpUpdateTree();

	void hpUpdateControlFuncs(double dt);
	
	// Moves branches marked as being done to a free list for cleanup by a 
	// lower priority thread.
	void hpUpdateFreeList();

	// TODO: are these needed???
	// Reclaims memory and returns number of events playing
	bool check();
};




// IMPLEMENTATION ______________________________________________________________

template <class TAudioIOData>
void SchedulerAudioIOData::mapAudioIOData(TAudioIOData& externalIO){
	userData(&externalIO);
	buffersIn = externalIO.inBuffer();
	buffersOut = externalIO.outBuffer();
	framesPerSecond = externalIO.framesPerSecond();
	framesPerBuffer = externalIO.framesPerBuffer();
	channelsIn = externalIO.channelsIn();
	channelsOut = externalIO.channelsOut();
}

template <class TAudioIOData>
TAudioIOData& SchedulerAudioIOData::unmapAudioIOData(){
	if(!mUserData /*|| !audioIODataType<TAudioIOData>(false)*/){
		gam::err(
			"member 'userData' is NULL. Did you make a matching call to mapAudioIOData?",
			"gam::SchedulerAudioIOData::unmapAudioIOData()"
		);
	}
	else if(typeID<TAudioIOData>() != mUserDataTypeID){
		gam::err(
			"Type mismatch between member 'userData' and template parameter.",
			"gam::SchedulerAudioIOData::unmapAudioIOData()"
		);
	}
	TAudioIOData& res = *(TAudioIOData *)mUserData;
	res.frame(startFrame);
	return res;
}


template <class TAudioIOData>
void Scheduler::audioCB(TAudioIOData& aio){
	if(!aio.user()){
		gam::err(
			"AudioIOData user data is NULL. Should be the address of the Scheduler.",
			"gam::Scheduler::audioCB()");
	}
	Scheduler& s = *(Scheduler*)aio.user();
	s.update(aio);
}

} //gam::

#endif
