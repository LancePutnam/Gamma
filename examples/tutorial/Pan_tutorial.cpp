/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Example:		Tutorial / Pan_tutorial by Josh Dickinson, 2012
	Description:	Panning noise back and forth using the Pan<> class
*/

#include "../examples.h"

//Initialize a noise source
NoisePink<> src;           

//We initialize a Pan class, which is an equal power panning object.  
Pan<> panner;

//Initialize a Sine class to control the position of the panner
Sine<> pannerControl;

int frameCount = 0;

void audioCB(AudioIOData& io){
    
	while(io()){
        pannerControl.freq(0.2);
        
        float tempPanPosition = pannerControl();
        
        //set the panner's position from -1.0 (completely left) to 1.0 (completely right)
        panner.pos(tempPanPosition);
        
        //print the current pan position every 1000 frames
        if(frameCount% 1000 == 0){
            std::cout<< "current pan position: " << tempPanPosition << std::endl;
        }
        
        //get pink noise every frame multiplied by some gain 
		float s = src() * 0.2;
        float rightChannel;
        float leftChannel;
        
        //calling the panner object like this takes an input (s) and two destinations(leftChannel and rightChannel).  It sets leftChannel and rightChannel to an appropriate volume based on the current panning position. You can also pass it two inputs and two outputs in order to pan a stereo source.   
        panner(s, leftChannel, rightChannel);
    
		io.out(0) = leftChannel;
        io.out(1) = rightChannel;
        
        //increment frame count
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