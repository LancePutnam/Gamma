/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include <string>
#include <stdio.h>
#include "mem.h"
#include "SoundFile.h"

#define ULONG unsigned long

namespace gam{
using std::string;

SoundFile::SoundFile(string pathA){
	fp = 0;
	path(pathA);
	memset(&inst, 0, sizeof(inst));
}

SoundFile::SoundFile(string pathA, const SoundFile & infoSrc){
	fp = 0;
	path(pathA);
	memset(&inst, 0, sizeof(inst));
	info(infoSrc);
}

SoundFile::~SoundFile(){
	close();
}

bool SoundFile::openRead(){
	mInfo.format = 0;	// unless opening RAW
	fp = sf_open(mPath.c_str(), SFM_READ, &mInfo);
	return 0 != fp;
}

bool SoundFile::openWrite(){

	// set major format based on path
	formatMajor(formatMajor(mPath));

	fp = sf_open(mPath.c_str(), SFM_WRITE, &mInfo);
	if(fp) sf_command(fp, SFC_SET_CLIPPING, NULL, SF_TRUE) ;
	return 0 != fp;
}

bool SoundFile::close(){
	bool didClose = true;
	if(0 != fp && 0 == sf_close(fp))	fp = 0;
	else								didClose = false;
	return didClose;
}

/*
typedef struct
{	int gain ;
	char basenote, detune ;
	char velocity_lo, velocity_hi ;
	char key_lo, key_hi ;
	int loop_count ;

	struct
	{	int mode ;
		unsigned int start ;
		unsigned int end ;
		unsigned int count ;
	} loops [16] ; // make variable in a sensible way
} SF_INSTRUMENT ;

enum
{
	//The loop mode field in SF_INSTRUMENT will be one of the following.
	SF_LOOP_NONE = 800,
	SF_LOOP_FORWARD,
	SF_LOOP_BACKWARD,
	SF_LOOP_ALTERNATING
} ;

*/
bool SoundFile::instrument(){
	bool hasInst = SF_TRUE == sf_command(fp, SFC_GET_INSTRUMENT, &inst, sizeof (inst));
	if(hasInst){
		printf("Gain:      %d\n", inst.gain);
		printf("Base Note: %d\n", inst.basenote);
		printf("Detune:    %d\n", inst.detune);
		printf("Velocity:  [%d, %d]\n", inst.velocity_lo, inst.velocity_hi);
		printf("Key:       [%d, %d]\n", inst.key_lo, inst.key_hi);
		printf("Loops:     %d\n", inst.loop_count);
		for(int i=0; i<inst.loop_count; i++){
			printf("\t");
			switch(inst.loops[i].mode){
				case SF_LOOP_FORWARD:     printf("fwd "); break;
				case SF_LOOP_BACKWARD:    printf("bwd "); break;
				case SF_LOOP_ALTERNATING: printf("alt "); break;
				default:                  printf("none");
			}
			printf(" [%d, %d] x %d\n", inst.loops[i].start, inst.loops[i].end, inst.loops[i].count);
		}
	}
	return hasInst;
}


void SoundFile::info(const SoundFile & src){
	memcpy(&mInfo, &src.mInfo, sizeof(mInfo));
}

const char * SoundFile::extension(){
	formatInfoMajor();
	return formatInfo.extension;
}

void SoundFile::formatInfoMajor(){
	formatInfo.format = format() & SF_FORMAT_TYPEMASK;
	sf_command(fp, SFC_GET_FORMAT_INFO, &formatInfo, sizeof(formatInfo));
}

void SoundFile::formatInfoSubtype(){
	formatInfo.format = format() & SF_FORMAT_SUBMASK;
	sf_command(fp, SFC_GET_FORMAT_INFO, &formatInfo, sizeof(formatInfo));
}

int SoundFile::formatMajor(string sfpath){
	unsigned int pos = sfpath.find_last_of(".");
	
	if(string::npos != pos){
		string ext = sfpath.substr(pos+1);
		//printf("%s\n", ext.c_str());
		
		if(ext == "wav")  return SF_FORMAT_WAV;
		if(ext == "aiff") return SF_FORMAT_AIFF;
	}
	
	return 0;
}

void SoundFile::print(){
	printf("%s\n", mPath.c_str());
	printf("Sample rate: %f\n",  frameRate());
	printf("Frames:      %lu\n", frames());
	printf("Channels:    %lu\n", channels());
	printf("Samples:     %lu\n", samples());
	
	formatInfoMajor();
	printf("Format:      %s, %s, ", formatInfo.name, formatInfo.extension) ;
	
	formatInfoSubtype();
	printf("%s\n", formatInfo.name) ;
}

} // end namespace gam

#undef ULONG

