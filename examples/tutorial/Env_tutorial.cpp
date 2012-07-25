/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Example:		Generator / Envelope
	Description:	Using an Env to control white noise. 
*/

#include "../examples.h"

//Env is initialized with a number of segments, in this case 4.    
Env<4> envelope;

// Noise source
NoiseWhite<> src;           

//this is a timer with a period of 1/2.0 = 0.5, meaning it will fire 0.5 times a second or once every two seconds.
//in this example, this master timer is reseting our envelope() object, which causes an attack every two seconds.  
Accum<> attackTimer(1/2.0);

int frameCount = 0;

void audioCB(AudioIOData& io){


	while(io()){
	
        //every frame we increment the internal clock of the attack timer.  It only returns true every 2.0 seconds. 
		if(attackTimer()){
            envelope.reset();
            
            //change the breakpoint levels of the envelope, so in this case the attack will ramp up from 0 to 1, decay from 1 to 0.8, systain from 0.8 to 0.8, then release from 0.8 to 0.0.  There should always be one more breakpoint level than there are segments. 
            envelope.levels(0, 1, 0.8, 0.8, 0);
            
            //set the envelope lengths, so in this case, the attack will last 0.01 seconds, the decay will last 0.1, the sustain 0.9, and the release 0.1
            envelope.lengths(0.01, 0.1, 0.9, 0.1);
            
            //set the totalLength of the envelope to 1.5 seconds by scaling all of the segments
            envelope.totalLength(1.5);
            
            //totalLength() can take a second argument of the segment to extend (starting from 0), so instead of scaling all of the segments, it will just expand or contract segment 2 (sustain)
            envelope.totalLength(1.5, 2);
		}
        
        //print the current value of the envelope every frame
        if(frameCount%1000 == 0)
            std::cout<< "current envelope amplitude: " << envelope()<<std::endl;
        
        //multiply constant white noise times our envelope every frame
		float s = src() * envelope() * 0.2;

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
