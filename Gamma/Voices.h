#ifndef GAMMA_VOICES_H_INC
#define GAMMA_VOICES_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include <map>
#include <cstdint> // uint64_t
#include <cmath>
#include <functional>
#include "Gamma/Types.h" // magSqr

namespace gam{

/// Index pool

/// This is a non-intrusive mechanism that can be used in conjunction with an
/// array to create an object pool. Obtaining and recycling an index is O(1).
/// For efficiency, the maximum supported size is 64.
class IndexPool{
public:

	typedef uint64_t bits_t;
	typedef unsigned index_t;
	
	static constexpr index_t npos = ~index_t(0);
	static constexpr index_t maxSize = sizeof(bits_t)*8;


	/// @param[in]	Maximum number of indices (up to 64)
	IndexPool(unsigned size = 8)
	:	mSize(size), mAutoPlay(~bits_t(0))
	{
		recycleAll();
	}

	/// \returns maximum number of indices
	unsigned size() const { return mSize; }

	/// \returns number of indices left in pool
	unsigned left() const { return mLeft; }

	/// \returns number of indices being used in pool
	unsigned used() const { return mSize-mLeft; }

	bool active(index_t i) const { return !!(mActive & bit(i)); }
	
	bool playing(index_t i) const { return !!(mPlaying & bit(i)); }

	bits_t activeBits() const { return mActive; }
	
	bits_t playingBits() const { return mPlaying; }


	/// Set whether a newly obtained index is automatically played
	void autoPlay(bool v){ mAutoPlay = -bits_t(v); }

	/// Set whether an index is playing
	void playing(index_t i, bool v){
		if(v)	mPlaying |=  bit(i);
		else	mPlaying &= ~bit(i);
	}

	/// Obtain the next free index
	
	/// \returns an object index or IndexPool::npos if there are no available
	/// indices.
	index_t obtain(){
		if(left()){
			auto idx = trailingZeroes(~mActive);
			obtain(idx);
			//printf("IndexPool: obtain %d\n", idx);
			return idx;
		}
		return npos;
	}

	/// Obtain the next free index and assign to it an integer ID
	index_t obtainWithID(unsigned ID){
		auto idx = obtain();
		if(idx != npos) mIDToIndex[ID] = idx;
		return idx;
	}

	// Set whether an index is active
	void active(index_t i, bool v){
		auto currActive = active(i);
		if(v && !currActive){
			obtain(i);
		} else if(!v && currActive){
			recycle(i);
		}
	}

	/// Free an index
	void recycle(index_t i){
		auto mask = bit(i);
		mLeft += (mActive & mask) >> i;
		mActive &= ~mask;
		mPlaying &= ~mask;
	}

	/// Free all indices
	void recycleAll(){
		mActive = 0;
		mPlaying = 0;
		mLeft = size();
	}

	/// Get index with given integer ID
	index_t indexWithID(unsigned ID){
		auto it = mIDToIndex.find(ID);
		if(it != mIDToIndex.end()){
			return it->second;
		}
		else{
			return npos;
		}
	}

	std::map<unsigned, index_t>& idToIndex(){
		return mIDToIndex;
	}


	/// Forward iterator through active indices
	struct iterator{

		iterator(bits_t bits)
		:	mBits(bits), mTemp(mBits & -mBits){}
		
		iterator(const iterator& mit)
		:	mBits(mit.mBits), mTemp(mit.mTemp){}

		iterator& operator++(){ // ++it
			mBits ^= mTemp;
			mTemp = mBits & -int64_t(mBits);
			return *this;
		}

		iterator operator++(int){ // it++
			iterator tmp(*this);
			operator++();
			return tmp;
		}

		bool operator==(const iterator& rhs){return mBits == rhs.mBits;}
		bool operator!=(const iterator& rhs){return !((*this)==rhs);}
		
		index_t operator*(){
			return deBruijnBitPosition(mTemp);
		}

	private:
		bits_t mBits;
		bits_t mTemp;
	};
	
	iterator begin() const { return iterator(mActive & mPlaying); }
	iterator end() const { return iterator(0); }


	void print(FILE * fp = stderr) const {
		fprintf(fp, "active (%3d): ", size()-mLeft);
		for(unsigned i=0; i<size(); ++i){
			fprintf(fp, "%c", active(i)?(playing(i)?'|':':'):'.');
		}
		fprintf(fp, "\n");
	}

private:
	bits_t mActive;		// bit array of active objects
	bits_t mPlaying;	// bit array of playing objects
	unsigned mLeft;		// number of free slots left
	unsigned mSize;		// maximum number of objects
	std::map<unsigned, index_t> mIDToIndex;
	bits_t mAutoPlay;

	void obtain(index_t i){
		auto mask = bit(i);
		if(0 == (mActive & mask)){
			mActive |= mask;
			mPlaying |= mask & mAutoPlay;
			--mLeft;
		}
	}

	static unsigned deBruijnBitPosition(uint64_t v){
		static const unsigned char deBruijn64[64] = {
			 0,  1,  2,  7,  3, 13,  8, 19, //  0- 7
			 4, 25, 14, 28,  9, 34, 20, 40, //  8-15
			 5, 17, 26, 38, 15, 46, 29, 48, // 16-23
			10, 31, 35, 54, 21, 50, 41, 57, // 24-31
			63,  6, 12, 18, 24, 27, 33, 39, // 32-39
			16, 37, 45, 47, 30, 53, 49, 56, // 40-47
			62, 11, 23, 32, 36, 44, 52, 55, // 48-55
			61, 22, 43, 51, 60, 42, 59, 58, // 56-63
		};

		return deBruijn64[(v*0x0218a392cd3d5dbfull) >> 58];
	}
	
	static unsigned trailingZeroes(uint64_t v){
		return deBruijnBitPosition(v & -int64_t(v));
	}

	static bits_t bit(unsigned i){ return bits_t(1)<<i; }
};



// These functions are used for measuring amplitude of different sample types.
namespace{
	template<unsigned N, class T>
	T ampSqr(const gam::Vec<N,T>& v){ return v.magSqr(); }

	template<class T>
	T ampSqr(const gam::Complex<T>& v){ return v.magSqr(); }

	template <class T>
	static T ampSqr(T v){ return float(v*v); }
};


class VoicesBase;
//template<class,unsigned> class Voices;


/// Representation of a single synthesizer voice

/// This class is to be subclassed by the user for a specific synth.
/// The function call operator must be overridden and return the next sample.
/// The onAttack and onRelease functions are overridden to trigger the attack
/// and release phases of a note. The 'done' function should be overriden to
/// indicate when the note is completely silent and eligible for placement back
/// in the voice pool. If 'done' is not overridden, then a basic silence
/// detection algorithm is used.
///
/// \tparam Tv			Value (sample) type of the voice generator
template <class Tv = float>
class Voice{
public:
	typedef Tv value_type;

	Voice(){
		silenceThreshold(0.001);
	}

	/// Generates next sample (must be overridden)
	Tv operator()() = delete;

	/// Called when the voice is started (key down)
	void onAttack(){}

	/// Called when the voice is released (key up)
	void onRelease(){}

	/// Called to check if the voice is done playing
	
	/// Unless overridden, this uses a simple silence detection algorithm
	/// based on power thresholding.
	bool done(){
		return detectSilence();
	}


	/// Get age of voice, in samples
	long age() const { return mAge; }

	/// Set delay before voice activates on attack
	void delay(long samples){
		mAge = -samples;
	}

	/// Set magnitude threshold for silence detection
	void silenceThreshold(float v){
		mSilenceThresh = v*v; // since we compare to magnitude squared
	}
	
	unsigned silenceCount() const { return mSilenceCount; }

	/// Get global voice parameter
	float getParam(unsigned idx) const;
	float getParamSafe(unsigned idx) const;

private:
	template<class,unsigned> friend class Voices;
	//friend class VoicesBase;

	Tv mOutput = Tv(0);
	float mSilenceThresh;
	unsigned short mSilenceCount = 0;
	unsigned short mSilenceCountMax = 12000;
	long mAge = 0; // max 27 hours @ 44.1 kHz
	unsigned mIndex = 0;
	VoicesBase * mParent = nullptr;

	// Returns true if the input has been below the threshold for the
	// maximum number of samples.
	bool detectSilence(){
		if(ampSqr(mOutput) < mSilenceThresh){
			++mSilenceCount;
			return mSilenceCount >= mSilenceCountMax;
		}
		mSilenceCount = 0;
		return false;
	}
	
	void voiceReset(){
		mSilenceCount = 0;
		mAge = 0;
	}
	
	void voiceUpdate(const Tv& newOutput){
		mOutput = newOutput;
		++mAge;
	}
};



/// A synthesis parameter with exponential (one-pole low-pass) smoothing
class Param{
public:

	Param(float initVal=0, float smooth=0.7)
	:	mValue(initVal), mTarget(initVal), mSmooth(smooth)
	{}


	float value() const { return mValue; }

	Param& operator= (float target){
		mTarget = target;
		return *this;
	}

	Param& init(float v){
		mValue = mTarget = v;
		return *this;
	}

	Param& smooth(float v){
		mSmooth = v;
		return *this;
	}

	void update(){
		//mValue = mValue*mSmooth + mTarget*(1-mSmooth);
		mValue = mTarget + mSmooth*(mValue - mTarget);
	}

private:
	float mValue, mTarget;
	float mSmooth;
};




class VoicesBase {
public:
	/// Set number of global voice parameters
	void numParams(unsigned size){
		mParams.resize(size);
	}

	/// Get number of global voice parameters
	unsigned numParams() const {
		return unsigned(mParams.size());
	}

	/// Get global voice parameter (object)
	Param& param(unsigned idx){ return mParams[idx]; }

	/// Get global voice parameter
	float getParam(unsigned idx) const { return mParams[idx].value(); }

	/// Initialize global voice parameter
	void initParam(unsigned idx, float val){
		if(idx >= mParams.size()) numParams(idx+1);
		mParams[idx].init(val);
	}

	/// Initialize global voice parameter
	void initParam(unsigned idx, float val, float smooth){
		initParam(idx, val);
		mParams[idx].smooth(smooth);
	}

	/// Set global voice parameter
	void setParam(unsigned idx, float val){ mParams[idx]=val; }

	void updateParams(){
		for(unsigned i=0; i<mParams.size(); ++i){
			mParams[i].update();
		}
	}

protected:
	std::vector<Param> mParams;
};


/// Fixed-size voice pool for polyphonic synthesis

/// This is a fixed-size voice pool that can used to create a polyphonic
/// synth or implement granular synthesis. Since the size is fixed, there are
/// no memory allocations at run-time. This means, however, that if the maximum
/// number of voices is exceeeded, then a voice stealing algorithm takes effect.
///
/// \tparam VoiceGen	A voice generator function object (a subclass of Voice)
/// \tparam Nvoices		Maximum number of active voices (up to IndexPool::maxSize)
template <class VoiceGen, unsigned Nvoices = 16>
class Voices : public VoicesBase {
public:

	static_assert(Nvoices <= IndexPool::maxSize, "Max size of index pool exceeded");
	
	typedef typename VoiceGen::value_type Tv;

	/// Voice stealing policy
	enum StealPolicy{
		OLDEST=0,		/** Steal oldest voice */
		NEWEST,			/** Steal newest voice */
		//QUIETEST		/** Steal quietest voice */
		//PRIORITY,		/* User-defined priority */
		//SAME_ID,		/* Steal voice with same ID */
		NONE,			/* No stealing */
	};


	Voices()
	:	mIndexPool(Nvoices), mStealIdx(0), mStealPolicy(OLDEST)
	{
		mIndexPool.autoPlay(false);
		for(unsigned i=0; i<Nvoices; ++i){
			mVoiceGens[i].mParent = this;
		}
	}


	/// Get maximum number of voices
	unsigned size() const { return mIndexPool.size(); }

	/// Get number of active voices
	unsigned numActive() const { return mIndexPool.size() - mIndexPool.left(); }

	/// Set whether voices should start playback immediately upon acquisition
	void autoPlay(bool v){ mIndexPool.autoPlay(v); }

	/// Flag a voice to start playback
	void play(VoiceGen& v){
		mIndexPool.playing(v.mIndex, true);
	}


	/// Get an array of the voices
	VoiceGen * voices(){ return mVoiceGens; }

	/// Get a reference to the voice at the specified index
	VoiceGen& voice(unsigned i){ return mVoiceGens[i]; }

	VoiceGen * begin(){ return mVoiceGens; }
	const VoiceGen * begin() const { return mVoiceGens; }
	VoiceGen * end(){ return mVoiceGens + size(); }
	const VoiceGen * end() const { return mVoiceGens + size(); }

	/// Get a reference to the voice with the given ID
	VoiceGen& voiceWithID(unsigned ID){
		auto idx = mIndexPool.indexWithID(ID);
		if(idx != mIndexPool.npos){
			return voice(idx);
		}
		else{ // ERROR!
			// Steal a voice...
			//updateStealIndex();
			//return mVoiceGens[mStealIdx];
			// Return a dummy voice
			return dummy();
		}
	}


	/// Obtain a new voice index
	unsigned obtainIndex(){
		auto idx = mIndexPool.obtain();
		
		// steal voice
		if(idx == mIndexPool.npos){
			updateStealIndex();
			idx = mStealIdx;
		}

		voice(idx).mIndex = idx;
		voice(idx).voiceReset(); // reset age, etc.
		return idx;
	}

	/// Obtain new voice without playing it
	VoiceGen& obtain(){
		return voice(obtainIndex());
	}

	/// Start new voice

	/// This calls Voice::onAttack.
	///
	VoiceGen& attack(){
		auto& v = obtain();
		v.onAttack(); // user-defined attack
		play(v);
		return v;
	}

	/// Start new voice

	/// The arguments are user-defined variables that are forwarded to
	/// Voice::onAttack. Typical parameters are frequency and amplitude.
	VoiceGen& attack(double a, double b){
		auto& v = obtain();
		v.onAttack(a,b); // user-defined attack
		play(v);
		return v;
	}

	/// Start new voice with given ID
	VoiceGen& attackWithID(unsigned ID){
		auto& v = attack();
		mIndexPool.idToIndex()[ID] = v.mIndex;
		return v;
	}

	/// Start new voice with given ID
	VoiceGen& attackWithID(unsigned ID, double a, double b){
		auto& v = attack(a,b);
		mIndexPool.idToIndex()[ID] = v.mIndex;
		return v;
	}

	/// Release voice with given ID
	VoiceGen& releaseWithID(unsigned ID){
		auto& v = voiceWithID(ID);
		v.onRelease();
		return v;
	}

	/// Releases all voices (to kill stuck notes)
	Voices& releaseAll(){
		for(auto& v : mVoiceGens) v.onRelease();
		return *this;
	}

	void updateStealIndex(){
		unsigned stealIdx = mStealIdx;

		switch(mStealPolicy){
		case OLDEST:{
			long stealCount = -2147483648;
			// Find voice with largest age value
			for(unsigned i=0; i<Nvoices; ++i){
				if(mIndexPool.active(i) && voice(i).age() > stealCount){
					stealCount = voice(i).age();
					stealIdx = i;
				}
			}
			}
			break;
		case NEWEST:{
			long stealCount = 2147483647;
			// Find voice with smallest age value
			for(unsigned i=0; i<Nvoices; ++i){
				if(mIndexPool.active(i) && voice(i).age() < stealCount){
					stealCount = voice(i).age();
					stealIdx = i;
				}
			}
			}
			break;
		/*case QUIETEST:
			for(unsigned i=0; i<Nvoices; ++i){
				if(mIndexPool.active(i) && mVoiceGens[i].silenceCount() > stealCount){
					stealCount = mStealCounts[i].count();
					stealIdx = i;
				}
			}
			break;*/
		default: // No stealing; set index to dummy voice
			stealIdx = Nvoices;
		}

		mStealIdx = stealIdx;
	}

	
	/// Returns next sum of all active voices
	Tv operator()(){

		updateParams();

		Tv res = Tv(0);

		for(auto idx : mIndexPool){

			auto& v = voice(idx);

			if(v.age() >= 0){
				Tv s = v(); // compute next sample
				res += s;
				v.voiceUpdate(s);

				// Check if voice is done
				if(v.done()){
					mIndexPool.recycle(idx);
				}
			}

			else{
				++v.mAge;
			}
		}
		return res;
	}

	/* Generate next buffer (surprisingly, no faster than single-sample!)
	void operator()(Tv * buf, unsigned size){
		for(int i=0; i<size; ++i) buf[i] = Tv(0);
		for(auto idx : mIndexPool){
			for(int i=0; i<size; ++i){
				Tv s = mVoiceGens[idx]();
				mVoiceGens[idx].voiceUpdate(s);
				buf[i] += s;
			}
			// Check if voice is done
			if(mVoiceGens[idx].done()){
				mIndexPool.recycle(idx);
				continue;
			}
		}
	}
	//*/

	/// Traverse all active voices
	Voices& traverseActive(const std::function<void(unsigned i)>& onVisit){
		for(auto idx : mIndexPool) onVisit(idx);
		return *this;
	}

	// Call member function on all active voices
	template <class Obj, class R>
	Voices& call(R (Obj::*func)()){
		return traverseActive([&,this](unsigned i){ (voice(i).*func)(); });
	}

	// Call member function on all active voices
	template <class Obj, class R, class A, class L>
	Voices& call(R (Obj::*func)(A), L l){
		return traverseActive([&,this](unsigned i){ (voice(i).*func)(l); });
	}

	// Call member function on all active voices
	template <class Obj, class R, class A, class B, class L, class M>
	Voices& call(R (Obj::*func)(A,B), L l, M m){
		return traverseActive([&,this](unsigned i){ (voice(i).*func)(l,m); });
	}


	/// Set the voice stealing policy
	Voices& stealPolicy(StealPolicy v){
		mStealPolicy = v; return *this;
	}


	void print() const {
		mIndexPool.print();
	}

private:
	VoiceGen mVoiceGens[Nvoices+1];
	IndexPool mIndexPool;
	StealPolicy mStealPolicy;
	unsigned mStealIdx;
	
	// Returned when voice stealing is off or encounters an error
	VoiceGen& dummy(){ return mVoiceGens[Nvoices]; }
};


template <class Tv>
inline float Voice<Tv>::getParam(unsigned idx) const {
	return mParent->getParam(idx);
}

template <class Tv>
float Voice<Tv>::getParamSafe(unsigned idx) const {
	if(mParent && (idx < mParent->numParams())){
		return mParent->getParam(idx);
	}
	else{
		return 0.f;
	}
}

} // gam::
#endif
