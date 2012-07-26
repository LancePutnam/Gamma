/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Example:		Tutorial / SegExp_tutorial by Josh Dickinson, 2012
	Description:	SegExp used to interpolate between two numbers
 */

#include "../examples.h"

//length, curve, end, start
//Curve<> mySegExp(10000,1, 50.0, 100.0);//,1.0,100.0);
SegExp<> mySegExp(10000,-0.01, 0.0, -100.0);//,1.0,100.0);

void audioCB(AudioIOData& io){
    
    while(io()){
        //every frame we increment the internal clock of the attack timer.  It only returns true every 2.0 seconds. 
		if(mySegExp.done()){
            //mySegExp.set(10000,0,1.0,100.0);
            mySegExp.reset();
        }
        
        float tempValue = mySegExp();
        
        //print the current value of the curve every frame
        std::cout<< "current curve value: " << tempValue <<std::endl;
        
        //multiply constant white noise times our curve every frame
		float s = 0;

        //left channel = right channel = s
		io.out(0) = io.out(1) = s;
	}
}

int main(){
    AudioIO io(256, 44100, audioCB, NULL, 2);
    Sync::master().spu(io.framesPerSecond());
    io.start();
    printf("Press 'enter' to quit...\n"); getchar();
    return 0;
}

