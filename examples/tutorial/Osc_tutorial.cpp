/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Example:		Tutorial / Osc_tutorial by Josh Dickinson, 2012
	Description:	Using the Osc class as a signal generator
*/

#include "../examples.h"

float ff = 110; // Fundamental frequency

Accum<> tmr(0.5);// Timer for modifying wavetable contents

//Osc is initialized with a frequency, a phase offset (usually 0), and a table length
Osc<> osc1(ff, 0, 4096);// Oscillator owning a 4096-element wavetable

//it can also be initialized with a reference to another oscillator, meaning it will share this oscillator's table. 
Osc<> osc2(ff + 0.17, 0, osc1);// Detuned oscillator sharing osc1's table

ArrayPow2<float> myTable(2048);

//it can also be initialized with reference to a pre-existing table 
Osc<> osc3(ff*6, 0, myTable);// oscillator using myTable as its lookup. 

int mode = -1;

int frameCount = 0;

void audioCB(AudioIOData& io){
    
	while(io()){
        
        if(tmr()){
            mode = (mode+1)%8;
            
            switch(mode){
                case 0:
                    //clear the oscillator's table
                    osc1.zero();
                    //we use the addWave() function to fill the osc1's table with a SAW wave.
                    addWave(osc1.table(), SAW);
                    break;
                case 1: 
                    //clear the oscillator's table
                    osc1.zero();
                    //you can also use any of these pre-defined waveform types: TRIANGLE, PARABOLIC, SQUARE, SAW, IMPULSE
                    addWave(osc1.table(), SQUARE);
                    //addWave(osc1.table(), TRIANGLE);
                    //addWave(osc1.table(), PARABOLIC);
                    //addWave(osc1.table(), SAW);
                    //addWave(osc1.table(), IMPULSE);
                    break;
                case 2:
                    //clear the oscillator's table
                    osc1.zero();
                    //there's a special function to fill a table with a sine wave
                    osc1.addSine(1);
                    break;
                case 3:
                    //clear the oscillator's table
                    osc1.zero();
                    //this function takes arguments for the number of cycles (harmonic from fundamental), amplitude, and phase offset. 
                    osc1.addSine(4, 1, 0);
                    break;
                case 4:
                    //clear the oscillator's table
                    osc1.zero();
                    //and if you don't zero() out the table in between, they stack on top of each other, allowing for easy additive synthesis
                    osc1.addSine(1);
                    osc1.addSine(2);
                    osc1.addSine(3);
                    osc1.addSine(4);
                    osc1.addSine(5);
                    osc1.addSine(6);
                    break;
                case 5:
                {
                    osc1.zero();
                    float amplitudes[] = {1,0.3,1,0.6,0.7,0.5,0.3,0.1};
                    //there's an easier way to add a bunch of sine tones together using the addSines (with an 's' at the end) function.  It takes arguments for the table we want to fill, an array of amplitude values, and a number of harmonics.  It then goes up the harmonic series (1, 2, 3, 4...) adding sine waves of the specified amplitudes.   
                    addSines(osc1.table(), amplitudes, 8);
                    break;
                }
                case 6:
                {
                    osc1.zero();
                    float amplitudes[] = {1,0.3,1,0.6,0.7,0.5,0.3,0.1};
                    float harmonics[] = {10, 11, 12, 13, 14, 15, 16, 17, 18};
                    //and if we want to specify exactly what harmonics we want, we can supply a second array with our choices.    
                    addSines(osc1.table(), amplitudes, harmonics, 8);
                    break;
                }
                case 7:
                    frameCount = 0;
                    osc1.zero();
                    addWave(osc1.table(), SQUARE);
                    break;
            }
            std::cout<< "current mode: " << mode << std::endl;
        }        
            
        
        if(mode == 7){
            if((frameCount%1000) == 0){
                //no matter what the wavetable holds, you can still change the osc's frequency just like you would any other generator
                osc1.freq(ff++);
            }
            if(ff > 300.0){
                ff = 110.0;
            }
        }
        //get the output of osc1 every frame 
		float s = osc1() * 0.2;
        
		io.out(0) = io.out(1) = s; 
        
        frameCount = frameCount + 1;
	}
}

int main(){
    AudioIO io(256, 44100, audioCB, NULL, 2);
    Sync::master().spu(io.framesPerSecond());
    io.start();
    printf("Press 'enter' to quit...\n"); getchar();
    return 0;
}
