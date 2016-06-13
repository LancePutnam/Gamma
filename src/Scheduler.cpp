#include <cstring> // memset
#include "Gamma/Scheduler.h"
#include "Gamma/SoundFile.h"
#include "Gamma/Timer.h"

/*
The audio thread reads through the event list and calls back into objects to
produce audio. If an object is finished making sound, then it is flagged as
'done'.

The main thread iterates through the event list periodically to check which
events are marked as done. When it finds a completed event, it removes it from
the event list and deletes its memory.

There is a potential problem if two separate threads are accessing the same
container. The audio thread will (or should) not modify the container, but
may access its items. The main thread will change the container when it removes
events. We therefore have the audio thread accessing a container that may be
concurrently having its items removed by another thread.

Solution 1:
Audio thread:
) lock event list
) add new events queued from main thread
) iterate events
) move finished events from list to somewhere else for main thread reclamation
) unlock event list

Responsibilities of threads:
Low-priority thread:
	Allocate/deallocate ProcessNode objects

High-priority thread:
	Add nodes:
		Move node from new list to processing tree
	Remove nodes:
		Move node from processing tree to free list


What happens when we add an event in the LPT?

In LPT:
1. Get current time, in seconds, from Scheduler start
2. Add event to queue with time plus event delta

In HPT:
1. Iterate through time-sorted event list
	If event time in block, then remove from list and add to process tree.
	Otherwise, end iteration.
	
2. Iterate through LP->HP event queue
	If event time in block, then add to process tree
	Otherwise, add to time-sorted event list


*/

/*
A ProcessNode is a node in a processing tree. Its placement in the tree determines 
when it is executed in relation to other nodes in the tree.

Consider the processing tree below.

	top
	|
	synths -> effects
	|         |
	|         chorus -> reverb
	|
	grain1 -> grain2

Execution of the tree occurs in a depth-first fashion resulting in the
processing sequence: top, synths, grain1, grain2, effects, chorus, reverb.
*/


/*
Scheduler s;
MySynth& synth = s.add<MySynth>(440, 0.2);
s.set(synth, &MySynth::freq, 330);
*/

namespace gam{

ProcessNode::ProcessNode(double delay)
:	mStatus(ACTIVE), mDelay(delay), mDeletable(false)
{}

ProcessNode::~ProcessNode(){
	removeFromParent();
	
	// delete children
	if(child){
	
		// remove all child's siblings first...
		while(child->sibling){
			if(child->sibling->deletable()) delete child->sibling;
			else child->sibling->removeFromParent();
		}
		
		// then remove the child
		if(child->deletable()) delete child;
		else child->removeFromParent();
	}
}

ProcessNode& ProcessNode::free(){
	mStatus = DONE;
	return *this;
}

ProcessNode& ProcessNode::active(bool v){ mStatus = v ? ACTIVE : INACTIVE; return *this; }

ProcessNode& ProcessNode::reset(){ onReset(); return *this; }

ProcessNode * ProcessNode::update(const ProcessNode * top, SchedulerAudioIOData& io){
	double dt = io.framesPerBuffer / io.framesPerSecond;
	unsigned frame = 0;
	if(mDelay >= dt){
		mDelay -= dt;
		return nextBreadth(top);
	}
	else if(mDelay > 0){	// 0 < delay <= delta
		frame = mDelay * io.framesPerSecond;
		mDelay=0;
		if(frame >= io.framesPerBuffer) return nextBreadth(top);
	}
	return process(top, io, frame);
}

ProcessNode * ProcessNode::process(const ProcessNode * top, SchedulerAudioIOData& io, int frameStart){
	if(active()){
		io.startFrame = frameStart;
		onProcessNode(io);
		if(active()) return next(top);
	}
	return nextBreadth(top);
}

void ProcessNode::print(){ printf("%p: %g sec, stat=%d\n", this, mDelay, mStatus); }



Scheduler::Scheduler()
:	mPeriod(1./10), mTime(0), mRunning(false)
{
	mDeletable = false;
}

Scheduler::~Scheduler(){
	stop();
}

bool Scheduler::empty() const {
	return (0==child) && mFreeList.empty() && mAddCommands.empty();
}

bool Scheduler::check(){
	reclaim(); return !empty();
}

int Scheduler::reclaim(){
	int r=0;
	while(!mFreeList.empty()){
		FreeList::value_type& v = mFreeList.front();
		
		Funcs::iterator it = mFuncs.begin();
		while(it != mFuncs.end()){
			const void * funcObj = it->mFunc.obj();
			if(funcObj == v){
				fprintf(stderr, "gam::Scheduler: reclaiming ProcessNode with active functions assigned\n");
				//printf("does %p == %p ?\n", (void *)v, funcObj);
				//it->mObjDel = 1;
				//while(it->mObjDel != 2){}
			}
			++it;
		}
		
		//printf("Scheduler: reclaiming %p\n", v);
		delete v;
		mFreeList.pop();
		++r;
	}
	return r;
}


void Scheduler::update(){

	double blockPeriod = io().framesPerBuffer / io().framesPerSecond;

	hpUpdateTree();
	hpUpdateControlFuncs(blockPeriod);

	// traverse tree
	ProcessNode * v = this;
	do{
		v = v->update(this, io());
	} while(v);
	
	// put nodes marked as 'done' into free list
	hpUpdateFreeList();
	
	mTime += blockPeriod;
}

Scheduler& Scheduler::period(float v){
	mPeriod=v;
	return *this;
}


void * Scheduler::cLPThreadFunc(void * user){
	Scheduler& s = *(Scheduler*)user;
	while(s.mRunning){
		//double t = gam::toSec(gam::timeNow());
		s.reclaim();
		//double dt = toSec(timeNow()) - t;
		//printf("%g\n", dt);
		::gam::sleepSec(s.mPeriod);
	}
	return NULL;
}

void Scheduler::start(){
	mLPThread.joinOnDestroy(true);
	mRunning = true;
	mLPThread.start(cLPThreadFunc, this);
}

void Scheduler::stop(){
	mRunning=false;
}


void Scheduler::cmdAdd(ProcessNode * v){
	pushCommand(Command::ADD_FIRST_CHILD, this, v);
}

void Scheduler::pushCommand(Command::Type type, ProcessNode * object, ProcessNode * other){
	other->mDeletable=true;
	Command c = { type, object, other };
	mAddCommands.push(c);
}

void Scheduler::hpUpdateControlFuncs(double dt){
	Funcs::iterator it = mFuncs.begin();
	while(it != mFuncs.end()){
		Funcs::value_type& f = *it;
		if(f.mDelay >= 0){
			f.mDelay -= dt;
			if(f.mDelay < 0){
				f();
				f.mDelay += f.mPeriod;
			}
			//printf("%g\n", f.mDelay);
			++it;
		}
		else{
			//printf("remove %p\n", &f);
			mFuncs.erase(it++);
		}
	}
}

void Scheduler::hpUpdateTree(){

	/* TODO:
	The tree should only contain processes that are active during the 
	current block. ProcessNodes still in the future should be kept in a 
	separate list and added to the tree when their time comes.
	In order to avoid scanning the whole list of future events, the events
	should be sorted according to their activation time. It is probably
	better to use absolute times rather than delta times so we don't
	have to perform any arithmetic to update timing status of the events.
	
	We should also avoid allocating memory for future ProcessNodes until when they
	are actually used. This will also allow us to use memory pools.
	*/

	while(!mAddCommands.empty()){
		Command& c = mAddCommands.front();
		
		switch(c.type){
		case Command::ADD_FIRST_CHILD:
			c.object->addFirstChild(c.other);
			break;
		case Command::ADD_LAST_CHILD:
			c.object->addLastChild(c.other);
			break;
		default:;
		}
		mAddCommands.pop();
	}	
}

void Scheduler::hpUpdateFreeList(){
	
	ProcessNode * p = this;
	ProcessNode * v = p->next(this);

	while(v){
		if(v->done()){
			v->removeFromParent();
			if(v->deletable()){	mFreeList.push(v); }
			v=p; // backtrack one node in traversal
		}
		else{
			p=v;
			v=v->next(this);
		}
	}		
}


void Scheduler::recordNRT(const char * soundFilePath, double durationSec){
	int numFrames = io().framesPerBuffer;
	int numChans  = io().channelsOut;

	SoundFile sf(soundFilePath);
	sf	.encoding(SoundFile::FLOAT)
		.channels(numChans)
		.frameRate(io().framesPerSecond)
	;
	if(sf.openWrite()){
		double  t = 0;
		double dt = io().secondsPerBuffer();
		
		// create buffer for interleaved samples
		float * buf = new float[numFrames*numChans];
		
		while(t < durationSec){

			//io.zeroOut();
			std::memset(io().buffersOut, 0, numChans*numFrames*sizeof(float));
			update();
			
			// interleave channel data
			for(int j=0; j<numChans; ++j){
				float * dst = buf + j;
				const float * src = io().buffersOut + j*numFrames;
				for(int i=0; i<numFrames; ++i){
					//printf("%d\n", i*numChans + j);
					dst[i*numChans] = src[i];
				}
			}

			// write to file
			sf.write(buf, numFrames);
			
			t += dt;
		}
		
		delete[] buf;
	}
}


} //gam::
