/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Example:		IO / AudioDevice
	Description:	How to start an audio stream and send line-input to output.
*/

#include "Gamma/Gamma.h"
#include "Gamma/AudioIO.h"

// A user defined class that can be accessed from the audio callback
struct UserData{
	float ampL, ampR;
};


// create a callback for generating a block of samples
void audioCB(gam::AudioIOData& io){

	UserData& user = *(UserData *)io.user();
	float ampL = user.ampL;
	float ampR = user.ampR;

	// loop through the number of samples in the block
	for(int i=0; i<io.framesPerBuffer(); ++i){
		
		float s = io.in(0,i);		// get the line-in or microphone sample
		
		io.out(0,i) = s * ampL;		// set left  output sample
		io.out(1,i) = s * ampR;		// set right output sample
	}

	// an alternative way to loop through the sample buffers
	while(io()){
		io.out(0) *= 0.5;			// scale left  output sample
		io.out(1) *= 0.5;			// scale right output sample
	}

	// if looping again, you must reset the frame iterator
	io.frame(0);
	while(io()){
		io.out(0) *= 2;				// scale left  output sample
		io.out(1) *= 2;				// scale right output sample
	}	
}


int main(){

	// check if we have any audio devices
	if(gam::AudioDevice::numDevices() == 0){
		printf("Error: No audio devices detected. Exiting...\n");
		exit(EXIT_FAILURE);
	}

	// list all detected audio devices
	printf("Audio devices found:\n");
	gam::AudioDevice::printAll();
	printf("\n");

	// set parameters of audio stream
	int blockSize = 64;				// how many samples per block?
	float sampleRate = 44100;		// sampling rate (samples/second)
	int outputChannels = 2;			// how many output channels to open
	int inputChannels = 1;			// how many input channels to open
	UserData user = {-0.5, 0.5};	// external data to be passed into callback
	
	// create an audio i/o object using default input and output devices
	gam::AudioIO io(blockSize, sampleRate, audioCB, &user, outputChannels, inputChannels);

	// we can also set the input and output devices explicitly
	// use device 0 for input and output
	//io.deviceIn (AudioDevice(0));
	//io.deviceOut(AudioDevice(0));

	// use devices matching keyword in name
	//io.deviceIn (AudioDevice("Microphone", AudioDevice::INPUT));
	//io.deviceOut(AudioDevice("Built-in", AudioDevice::OUTPUT));

	if(io.channelsOut() < 2){
		printf("This example needs at least 2 output channels to start streaming. Exiting...\n");
		return 0;
	}
	if(io.channelsIn() < 1){
		printf("This example needs at least 1 input channel to start streaming. Exiting...\n");
		return 0;
	}
	
	// set the global sample rate "subject"
	//gam::sampleRate(io.framesPerSecond());
	
	// start the audio stream
	io.start();
	
	// print some information about the i/o streams
	printf("Audio stream info:\n");
	io.print();
	
	printf("\nPress 'enter' to quit...\n"); getchar();
	return 0;
}
