/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Example:		Effect / Binaural Beat
	Description:	Demonstration of a binaural beating effect. This effect is produced
					by outputing two sine waves with slightly different frequencies
					to the left and right channels. With headphones on, the result sounds 
					like a single sine wave rotating around one's head at a rate that is
					the difference between the two sine wave frequencies.
*/

#include "../examples.h"

Sine<> src1(220);
Sine<> src2(221);

void audioCB(AudioIOData& io){
	while(io()){
		io.out(0) = src1() * 0.2;
		io.out(1) = src2() * 0.2;
	}
}

RUN_AUDIO_MAIN
