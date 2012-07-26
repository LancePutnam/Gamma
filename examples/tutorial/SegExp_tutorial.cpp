/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Example:		Tutorial / SegExp_tutorial by Josh Dickinson, 2012
	Description:	SegExp used to interpolate between two numbers
 */

#include "../examples.h"

// Timer to reset the curve
Accum<> tmr(0.25);

//SegExp is an exponential curve which asymptotically approaches a target value.  It is initialized with the arguments (length, curve, end, start)
SegExp<> mySegExp(2,-3.0, 0.0, 1.0);

//sound source
NoisePink<> pinkNoise;

int frameCount = 0;

void audioCB(AudioIOData& io){
    
    while(io()){
        //every frame we increment the internal clock of the attack timer.  It only returns true every 2.0 seconds. 
		if(tmr()){
            //mySegExp.set(10000,0,1.0,100.0);
            mySegExp.reset();
        }
        
        float tempValue = mySegExp();
        
        //print the current value of the curve every 1000 frames
        if((frameCount %1000) == 0)
        	std::cout<< "current SegExp value: " << tempValue <<std::endl;
        
        //multiply constant white noise times our curve every frame
		float s = pinkNoise()*0.2*tempValue;

        //left channel = right channel = s
		io.out(0) = io.out(1) = s;
		
		frameCount++;
	}
}

int main(){
    AudioIO io(256, 44100, audioCB, NULL, 2);
    Sync::master().spu(io.framesPerSecond());
    io.start();
    printf("Press 'enter' to quit...\n"); getchar();
    return 0;
}

