/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Example:		Filter / Delay-line
	Description:	This demonstrates the multiple uses of a delay-line.
*/

#include "../examples.h"

// Delay args: max delay, initial delay
Delay<float, ipl::Trunc> delay(0.4,	0.2);
Accum<> tmr(1);
Burst burst(2e4,2e3, 0.1);

void audioCB(AudioIOData& io){

	while(io()){
		
		if(tmr()) burst.reset();
		
		//float s = io.in(0)[i];
		float s = burst()*2;
		
		// This short-hand method is convenient for simple delays.
		s += delay(s);
		
		// The delay-line can also be read from and written to using separate calls.
		//delay.write(s);
		//s += delay();
		
		// We can create infinite echoes by feeding back a small amount of the 
		// output back into the input on each iteration.
		//s = delay(s + delay()*0.2);
		
		// We can also create mult-tap delay-lines through multiple calls to 
		// the read() method.
		//s += delay(s) + delay.read(0.15) + delay.read(0.39);

		// How about multi-tap feedback?
		//s += delay(s + delay.read(0.197)*0.3 + delay.read(0.141)*0.4 + delay.read(0.093)*0.2)*0.5;
	
		io.out(0) = io.out(1) = s;
	}
}

int main(){
	AudioIO io(128, 44100, audioCB, 0, 2, 1);
	Sync::master().spu(io.framesPerSecond());
	io.start();
	printf("\nPress 'enter' to quit...\n"); getchar();
	return 0;
}
