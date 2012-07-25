/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Example:		Tutorial / Sine_tutorial by Josh Dickinson, 2012
	Description:	Using Sine as a signal generator
*/

#include "../examples.h"

Accum<> tmr(3.0);// Timer to change the frequency

//Sine uses a polynomial approximation to compute sine values instead of reading from a table.  This saves memory but costs a little bit more CPU.  Lance says that the performance is almost as good as reading from a table, so you probably don't need to worry unless you're using a whole bunch of them simultaneously. TLDR: Sine<> generates a sine wave by computing the values every frame. 
Sine<> src;

float fundamentalFrequency = 200.0;

int mode = 0;

using namespace std;

void audioCB(AudioIOData& io){
    
	while(io()){
        if(tmr()){
            //good way to get random numbers in gamma.  Uniform distribution between two numbers
            int tempHarmonic = rnd::uni(1,6);
 
            //set the frequency by calling freq(freqYouWant)
            src.freq(tempHarmonic*fundamentalFrequency);
            
            //get the frequency by calling freq()
            cout<< "harmonic: " << tempHarmonic << endl << "frequency: " << src.freq() << endl;
        }          
        src.phaseAdd(rnd::uni(0,2000));
        
        float s = src()*0.2;
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
