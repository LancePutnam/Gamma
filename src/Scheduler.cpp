#include "Gamma/Scheduler.h"

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


Low-priority thread:
Allocate/deallocate Process objects

High-priority thread:
	Add nodes:
		Move node from new list to processing tree
	Remove nodes:
		Move node from processing tree to free list

*/
/*
A Process is a node in a processing tree. Its placement in the tree determines 
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

Process::Process(double delay)
:	mStatus(ACTIVE), mDelay(delay), mDeletable(false)
{}

Process::~Process(){
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

Process& Process::free(){
	mStatus = DONE;
	return *this;
}

Process& Process::active(bool v){ mStatus = v ? ACTIVE : INACTIVE; return *this; }

Process& Process::reset(){ onReset(); return *this; }

Process * Process::update(const Process * top, AudioIOData& io){
	double dt = io.secondsPerBuffer();
	int frame = 0;
	if(mDelay >= dt){
		mDelay-=dt;
		return nextBreadth(top);
	}
	else if(mDelay > 0){	// 0 < delay <= delta
		frame = mDelay * io.framesPerSecond();
		mDelay=0;
		if(frame >= io.framesPerBuffer()) return nextBreadth(top);
	}
	return process(top, io, frame);
}

void Process::print(){ printf("%p: %g sec, stat=%d\n", this, mDelay, mStatus); }

Process * Process::process(const Process * top, AudioIOData& io, int frameStart){
	if(active()){
		io.frame(frameStart);
		onProcess(io);
		if(active()) return next(top);
	}
	return nextBreadth(top);
}




Scheduler::Scheduler()
:	mPeriod(1./10), mRunning(false)
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
				printf("reclaiming Process with active functions\n");
				//printf("does %p == %p ?\n", (void *)v, funcObj);
				//it->mObjDel = 1;
				//while(it->mObjDel != 2){}
			}
			++it;
		}
		
		//printf("reclaim %p\n", v);
		delete v;
		mFreeList.pop();
		++r;
	}
	return r;
}

void Scheduler::update(AudioIOData& io){
	updateTree();
	updateControlFuncs(io.secondsPerBuffer());
	
	// traverse tree
	Process * v = this;
	do{
		v = v->update(this, io);
	} while(v);
	
	// put nodes marked as 'done' into free list
	updateFreeList();
}

Scheduler& Scheduler::period(float v){
	mPeriod=v;
	return *this;
}


void * Scheduler::cThreadFunc(void * user){
	Scheduler& s = *(Scheduler*)user;
	while(s.mRunning){
		s.mTime = toSec(timeNow());
		s.reclaim();
		//double dt = toSec(timeNow()) - s.mTime;
		//printf("%g\n", dt);
		gam::sleepSec(s.mPeriod);
	}
	return NULL;
}

void Scheduler::start(){
	mThread.joinOnDestroy(true);
	mRunning = true;
	mThread.start(cThreadFunc, this);
}

void Scheduler::stop(){
	mRunning=false;
}

	
void Scheduler::updateControlFuncs(double dt){
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


//void Scheduler::cmdAdd(Process * v){
//	v->mDeletable=true;
//	Command c = { Command::ADD_FIRST_CHILD, this, v };
//	mAddCommands.push(c);
//}

void Scheduler::cmdAdd(Process * v){
	pushCommand(Command::ADD_FIRST_CHILD, this, v);
}

void Scheduler::pushCommand(Command::Type type, Process * object, Process * other){
	other->mDeletable=true;
	Command c = { type, object, other };
	mAddCommands.push(c);
}

void Scheduler::updateTree(){
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

void Scheduler::updateFreeList(){
	
	Process * p = this;
	Process * v = p->next(this);

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

} //gam::
