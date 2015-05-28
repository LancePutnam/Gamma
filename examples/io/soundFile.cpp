/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information
	
	Example:		Audio / IO
	Description:	Writing and reading a sound file.
*/

#include "Gamma/Gamma.h"
#include "Gamma/SoundFile.h"

int main(){
	using namespace gam;

	const char * path = "test.aiff";
	SoundFile sf(path);

	
	// Write a beating sine wave to file
	{
		const float sampleRate = 44100;
		const float lenSec = 2;
		const int numFrames = sampleRate * lenSec;
		
		float freq = 440;
		float buf[numFrames*2];		// Buffer for storing samples.
									// The samples are stored interleaved as 
									// sequential stereo frames.
		sf.format(SoundFile::AIFF);
		sf.encoding(SoundFile::PCM_16);
		sf.channels(2);
		sf.frameRate(sampleRate);
		
		printf("Open new file for writing... ");
		if(sf.openWrite()){	printf("OK\n"); }
		else{				printf("fail\n"); exit(-1); }

		for(int i=0; i<numFrames; i++){
			
			float s = i/float(numFrames) * lenSec;	// compute time in seconds
			float p1 = s * freq * M_2PI;
			float p2 = s * (freq+4) * M_2PI;
			
			buf[i*2  ] = sin(p1)*0.5;		// store left channel
			buf[i*2+1] = sin(p2)*0.5;		// store right channel
		}
		
		
		// write entire buffer contents to file
		printf("Write samples to file... ");
		if(sf.write(buf, numFrames) == numFrames){	printf("OK\n"); }
		else{										printf("fail\n"); exit(-1); }
		
		sf.close();
	}


	// Read samples from the file we just wrote...
	{
		printf("Open file for reading... ");
		if(sf.openRead()){	printf("OK\n"); }
		else{				printf("fail\n"); exit(-1); }
	
		//sf.print();
		
		int numFrames = sf.frames();
		float buf[sf.samples()];
		
		//sf.readAll(buf);		// read all samples from sound file into buffer

		printf("Read samples from file... ");
		if(sf.readAll(buf) == numFrames){	printf("OK\n"); }
		else{								printf("fail\n"); exit(-1); }

		sf.close();
	

		// Rewrite sound file to make sure we read samples correctly
		printf("Open existing file for rewriting... ");
		if(sf.openWrite()){	printf("OK\n"); }
		else{				printf("fail\n"); exit(-1); }

		printf("Rewrite samples to file... ");
		if(sf.write(buf, numFrames) == numFrames){	printf("OK\n"); }
		else{										printf("fail\n"); exit(-1); }

		sf.close();
	}

	return 0;
}
