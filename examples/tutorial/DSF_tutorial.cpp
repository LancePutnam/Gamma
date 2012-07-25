/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Example:		Tutorial / DSF_tutorial by Josh Dickinson, 2012
	Description:	Using DSF as a signal generator. 
*/

#include "../examples.h"

Accum<> tmr(0.5);// Timer to change the filter's center frequency

float fundamentalFrequency = 200.0;

//DSF is a Discrete summation formula (DSF) oscillator. It takes optional arguments for the fundamental frequency, the frequency ratio, the amplitude ratio, and the number of harmonics. It produces a finite set of harmonics whose amplitudes follow a geometric series.
DSF<> src(fundamentalFrequency, 1, 0.5, 8);

int mode = -1;

using namespace std;

void audioCB(AudioIOData& io){
    
	while(io()){
        if(tmr()){
            mode = (mode+1)%4;
            switch(mode){
                case 0:                    
                    //as with most generators, we set the fundamental frequency with freq()
                    src.freq(fundamentalFrequency);
                    
                    //use harmonic() to set the number of harmonics (including the fundamental).  An argument of 1 gives us a simple sine wave. 
                    src.harmonics(1);
                    break;
                case 1:
                    //two harmonics with a frequency ratio of 1.0, forms an octave relationship
                    src.harmonics(2);
                    src.freqRatio(1.0);
                    break;
                case 2:
                    //two harmonics with a frequency ratio of 1.5, forms an octave relationship
                    src.harmonics(3);
                    //harmonic frequencies follow the formula ((i*fr) + 1) where 'fr' is the frequency ratio.  If the frequency ratio is 0.1, then the first harmonic will be ((1*0.1)+1) = 1.1 times the fundamental frequency.  The second harmonic will be ((2*0.1)+1) = 1.2 times the fundamental, etc.
                    src.freqRatio(0.1);
                    break;
                case 3:
                    //two harmonics with a frequency ratio of 1.5, forms an octave relationship
                    src.harmonics(40);
                    src.freqRatio(1.0);
                    //The amplitude of harmonic i is ar^i where 'ar' is called the amplitude ratio.  Anything below 1 gets exponentially smaller as it goes up the scale.  
                    src.ampRatio(0.3);
                    break;
            }
            std::cout<< "current mode: " << mode << std::endl;
        }          
        
        //get the next sample from our generator each frame
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
