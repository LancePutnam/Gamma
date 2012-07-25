/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Example:		Generator / Envelope
	Description:	Using an ADSR to control white noise. 
*/

#include "../examples.h"

// Noise source
NoiseWhite<> src;           

//ADSR is initialized with 6 values: attack time, decay time, sustain level (not sustain time), release time, amplidue, and curve value.  
ADSR<> ADSRenvelope(0.01, 0.1, 0.21, 0.3, 1, -4);

//this is a timer with a period of 1/2.0 = 0.5, meaning it will fire 0.5 times a second or once every two seconds.
//in this example, this master timer is reseting our ADSRenvelope() object, which causes an attack every two seconds.  
Accum<> attackTimer(1/2.0);

//release timer is another timer that fires after 0.7 seconds, telling the ADSRenvelope object to release().  
//It's only called AFTER the ADSR goes into sustain mode, so it determines the length that the ADSR is sustained for.  
Accum<> releaseTimer(1/0.7);

int frameCount = 0;

void audioCB(AudioIOData& io){

	while(io()){
	
        //every frame we increment the internal clock of the attack timer.  It only returns true every 2.0 seconds. 
		if(attackTimer()){
            ADSRenvelope.reset();
            releaseTimer.reset();
            
            //this is how you can change the ADSRenvelope's attach, decay, and release times:
            //ADSRenvelope.lengths(0.1,0.1,0.1);
            
            //if you want to turn off the sustain portion of the envelope, you can call sustainDisable():
            //ADSRenvelope.sustainDisable();
            
            //and if sustain is disabled, you can call loop(true) to make the ADR reset every time it completes:
            //ADSRenvelope.loop(true);
		}
        
        //ADSRenvelope.sustained() returns true if the envelope has entered its sustain mode. 
        if(ADSRenvelope.sustained()){
            //once the envelope is sustaining, we accumulate our release timer every frame until it fires, telling the ADSRenvelope to release().
            if(releaseTimer()){
                ADSRenvelope.release();
            }
        }
        
        //print the current value of the envelope every 1000 frames
        if(frameCount % 1000 == 0){
            std::cout<< "current envelope amplitude: " << ADSRenvelope()<<std::endl;
        }
        
        //multiply constant white noise times our envelope every frame
		float s = src() * ADSRenvelope() * 0.2;

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
