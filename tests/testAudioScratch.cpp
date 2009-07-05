//#include "AudioIO.h"
//#include "Sync.h"
//
//using namespace gam;
//
//void audioCB(AudioIOData & io){
//
//}
//
//
//int main(int argc, char* argv[]){
//	AudioIO io(256, 44100., audioCB, NULL, 2, 1);
//	Sync::master().spu(io.fps());
//	
//	io.start();
//	io.print();
//
//	printf("\nPress 'enter' to quit...\n"); getchar();
//	return 0;
//}


#include "Envelope.h"
#include "Oscillator.h"

namespace gam{

template <class T=gam::real>
struct Grain{

	Grain(T freq, T dur, T mul, T phase)
	:	osc(freq, phase), env(1./dur, 0, table), amp(mul)
	{}

	/// Generate next sample
	T operator()(){
		return osc() * env() * amp;  // doesn't work...
		//return osc() * amp;        // does work
	}

	// ensure static window table gets assigned to envelope member.
	// this must be called once after initTable()
	void assignTable(){
		env.source(table); 
	}

	void reset(){osc.phase(0); env.reset();}


	/// initialize the hann window
	static void fillTable(){
		tbl::hann(&table[0], table.size());	///< fill table with hann window
	}
	
	// force static table to allocate itself.
	// this will only allocate memory once (equal to the power of 2 given)
	// as long as the size does not change...
	static void initTable(){
		table.resize(10);
	}


	/// data members
	Sine<T> osc;                                    ///< Sine oscillator
	Osc<T, ipl::Linear, tap::Clip> env;             ///< Envelope
	T amp;	 ///< Grain Amplitude
	static ArrayPow2<T> table;	 ///< Envelope table
};

// call default constructor, memory not allocated yet...
template <class T> ArrayPow2<T> Grain<T>::table;


} // end namespace gam



#include <stdio.h>
#include "AudioIO.h"

using namespace gam;


Grain<>	grain(220, 2, 1.0, 0);


void audioCB(AudioIOData & io){

	for(uint32_t i=0; i<io.framesPerBuffer(); ++i){
		float s = grain();
		io.out(0)[i] = io.out(1)[i] = s * 0.5f;
	}
}



int main(int argc, char* argv[]){

	AudioIO io(256, 44100., audioCB, NULL, 2);
	Sync::master().spu(io.framesPerSecond());

	Grain<>::initTable();
	Grain<>::fillTable();

	grain.assignTable();

	io.start();
	printf("\nPress 'enter' to quit...\n"); getchar();
	return 0;
}