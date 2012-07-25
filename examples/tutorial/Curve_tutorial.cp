/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Example:		Generator / curve
	Description:	Curve used to interpolate between two numbers
 */

#include "../examples.h"

//length, curve, end, start
Curve<> curve1(10000,1, 50.0, 100.0);//,1.0,100.0);

void audioCB(AudioIOData& io){
    
    while(io()){
        //every frame we increment the internal clock of the attack timer.  It only returns true every 2.0 seconds. 
		if(curve1.done()){
            //curve1.set(10000,0,1.0,100.0);
            curve1.reset();
        }
        curve1();
        
        //print the current value of the curve every frame
        std::cout<< "current curve value: " << /*curve1.value()*/curve1() <<std::endl;
        
        //multiply constant white noise times our curve every frame
		float s = 0;

        //left channel = right channel = s
		io.out(0) = io.out(1) = s;
	}
}

RUN(audioCB);



