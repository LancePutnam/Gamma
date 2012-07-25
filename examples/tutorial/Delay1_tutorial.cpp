/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Example:		Tutorial / Delay1_tutorial by Josh Dickinson, 2012
	Description:	Using Delay1 - a one element (sample) delay line. 
*/

#include "../examples.h"

Accum<> tmr(0.5);// Timer to change the filter's center frequency

float fundamentalFrequency = 200.0;

DSF<> src(fundamentalFrequency, 1, 0.5, 8);

Delay1<> delay(0);

using namespace std;

void audioCB(AudioIOData& io){
    
	while(io()){
        if(tmr()){
            for(int i=1; i<5; ++i) cout << delay(i)<<endl;

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
