/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include <cctype> // tolower
#include <string>
#include <stdio.h>
#include "sndfile.h"
#include "Gamma/SoundFile.h"

namespace gam{
using std::string;

#ifdef GAM_NO_SOUNDFILE
class SoundFile::Impl {
public:

	double frameRate() const { return 1; }
	int frames() const { return 0; }
	int channels() const { return 0; }

	void channels(int num){}
	void frameRate(double hz){}

	bool close(){ return true; }

	int format() const { return 0; }
	int formatMajor() const { return 0; }
	int formatMinor() const { return 0; }
	void format(int v){}

	bool opened() const { return false; }

	bool openRead(const std::string& path){ return false; }
	bool openWrite(const std::string& path){ return false; }

	const char * extension(){ return ""; }
	int formatMajor(const string& sfpath){ return 0; }
	void formatMajor(int v){}
	void formatMinor(int v){}

	void info(const Impl& src){}

	template <class T> int read(T * dst, int numFrames){ return 0; }
	template <class T> int write(const T * src, int numFrames){ return 0; }

	void seek(int pos, int mode){}

	//void printInfo(){}
};
#else
class SoundFile::Impl{
public:
	Impl()
	:	fp(0)
	{
		mInfo.format = 0;
		memset(&mInst, 0, sizeof(mInst));
	}

	double frameRate() const { return (double)mInfo.samplerate; }
	int frames() const { return (int)mInfo.frames; }
	int channels() const { return (int)mInfo.channels; }

	void channels(int num){ mInfo.channels = (int)num; }
	void frameRate(double hz){ mInfo.samplerate = (int)hz; }

	bool close(){
		bool didClose = true;
		if(0 != fp && 0 == sf_close(fp))	fp = 0;
		else								didClose = false;
		return didClose;
	}

	void formatInfoMajor(){
		mFormatInfo.format = format() & SF_FORMAT_TYPEMASK;
		sf_command(fp, SFC_GET_FORMAT_INFO, &mFormatInfo, sizeof(mFormatInfo));
	}

	void formatInfoSubtype(){
		mFormatInfo.format = format() & SF_FORMAT_SUBMASK;
		sf_command(fp, SFC_GET_FORMAT_INFO, &mFormatInfo, sizeof(mFormatInfo));
	}

	int format() const { return mInfo.format; }
	int formatMajor() const { return mInfo.format & SF_FORMAT_TYPEMASK; }
	int formatMinor() const { return mInfo.format & SF_FORMAT_SUBMASK; }
	void format(int v){ mInfo.format = v; }

/*
When opening a file for read, the format field should be set to zero before 
calling sf_open(). The only exception to this is the case of RAW files where 
the caller has to set the samplerate, channels and format fields to valid 
values. All other fields of the structure are filled in by the library.
*/

	bool opened() const { return 0 != fp; }

	bool openRead(const std::string& path){
		if(formatMinor() != SF_FORMAT_RAW) mInfo.format = 0;
		fp = sf_open(path.c_str(), SFM_READ, &mInfo);

		//printf("%s\n", path.c_str());
		//printInfo();
		return opened();
	}

	bool openWrite(const std::string& path){

		// set major format based on path
		formatMajor(formatMajor(path));

		// if no encoding type, use sensible default...
		if(!formatMinor()) formatMinor(SF_FORMAT_PCM_16);

//		printInfo();

		fp = sf_open(path.c_str(), SFM_WRITE, &mInfo);

		if(fp) sf_command(fp, SFC_SET_CLIPPING, NULL, SF_TRUE);
		return opened();
	}


	const char * extension(){
		formatInfoMajor();
		return mFormatInfo.extension;
	}

	// get file format type based on extension
	int formatMajor(const string& sfpath){
		size_t pos = sfpath.find_last_of(".");
		
		if(string::npos != pos){
			string ext = sfpath.substr(pos+1);
			for(unsigned i=0; i<ext.size(); ++i) ext[i] = std::tolower(ext[i]);
			
			//printf("%s\n", ext.c_str());
			if(ext == "wav")					return SF_FORMAT_WAV;
			if(ext == "aiff" || ext == "aif")	return SF_FORMAT_AIFF;
			if(ext == "au")						return SF_FORMAT_AU;
			//if(ext == "flac")					return SF_FORMAT_FLAC;
		}
		return 0;
	}

	void formatMajor(int v){
		mInfo.format = (mInfo.format & (~SF_FORMAT_TYPEMASK)) | (v & SF_FORMAT_TYPEMASK);
		//printf("formatMajor(int v): %x\n", mInfo.format & SF_FORMAT_TYPEMASK);
	}

	void formatMinor(int v){
		mInfo.format = (mInfo.format & (~SF_FORMAT_SUBMASK)) | (v & SF_FORMAT_SUBMASK);
		//printf("formatMinor(int v): %x\n", mInfo.format & SF_FORMAT_SUBMASK);
	}

	void info(const Impl& src){
		memcpy(&mInfo, &src.mInfo, sizeof(mInfo));
	}

	template <class T> int read(T * dst, int numFrames);
	template <class T> int write(const T * src, int numFrames);

	void seek(int pos, int mode){
		sf_seek(fp, pos, mode);
	}

	void printInfo(){
		printf("\n");
		printf("frames=     %d\n", (int)mInfo.frames);
		printf("samplerate= %d\n", mInfo.samplerate);
		printf("channels=   %d\n", mInfo.channels);
		printf("format=     %x\n", mInfo.format);
		printf("sections=   %d\n", mInfo.sections);
		printf("seekable=   %d\n", mInfo.seekable);
	}

	SNDFILE * fp;
	SF_INFO mInfo;
	SF_INSTRUMENT mInst;
	SF_FORMAT_INFO mFormatInfo;
};

#define DEF_SPECIAL(type)\
template<>\
int SoundFile::Impl::read<type>(type * dst, int numFrames){\
	return sf_readf_##type(fp, dst, numFrames);\
}\
template<>\
int SoundFile::Impl::write<type>(const type * src, int numFrames){\
	return sf_writef_##type(fp, src, numFrames);\
}
DEF_SPECIAL(float)
DEF_SPECIAL(short)
DEF_SPECIAL(int)
DEF_SPECIAL(double)
#undef DEF_SPECIAL

#endif


SoundFile::SoundFile(const std::string& path_)
:	mImpl(new Impl)
{
	path(path_);
}

SoundFile::SoundFile(const std::string& path_, const SoundFile& infoSrc)
:	mImpl(new Impl)
{
	path(path_);
	info(infoSrc);
}

SoundFile::~SoundFile(){
	close();
	delete mImpl;
}

SoundFile& SoundFile::path(const std::string& v){ mPath=v; return *this; }
int SoundFile::samples() const { return frames() * channels(); }

const std::string& SoundFile::path() const { return mPath; }

bool SoundFile::opened() const { return mImpl->opened(); }

bool SoundFile::close(){ return mImpl->close(); }

SoundFile& SoundFile::info(const SoundFile& sf){ mImpl->info(*sf.mImpl); return *this; }

const char * SoundFile::extension(){ return mImpl->extension(); }
double SoundFile::frameRate() const { return mImpl->frameRate(); }
int SoundFile::frames() const { return mImpl->frames(); }
int SoundFile::channels() const { return mImpl->channels(); }	

SoundFile& SoundFile::channels(int num){ mImpl->channels(num); return *this; }
SoundFile& SoundFile::frameRate(double hz){ mImpl->frameRate(hz); return *this; }

SoundFile::Format SoundFile::format() const {
	#define CS(x) case SF_FORMAT_##x: return x;
	switch(mImpl->formatMajor()){
	CS(WAV) CS(AIFF) CS(AU) CS(RAW)
	//CS(FLAC)
	default: return Format(0);
	}
	#undef CS
}

SoundFile& SoundFile::format(Format v){
	#define CS(x) case x: mImpl->formatMajor(SF_FORMAT_##x); break;
	switch(v){
	CS(WAV) CS(AIFF) CS(AU) CS(RAW)
	//CS(FLAC)
	default:;
	}
	#undef CS
	return *this;
}

SoundFile::EncodingType SoundFile::encoding() const {
	#define CS(x) case SF_FORMAT_##x: return x;
	switch(mImpl->formatMinor()){
	CS(PCM_S8) CS(PCM_16) CS(PCM_24) CS(PCM_32) CS(PCM_U8)
	CS(FLOAT) CS(DOUBLE) CS(ULAW) CS(ALAW)
	default: return EncodingType(0);
	}
	#undef CS
}

SoundFile& SoundFile::encoding(EncodingType v){
	#define CS(x) case x: mImpl->formatMinor(SF_FORMAT_##x); break;
	switch(v){
	CS(PCM_S8) CS(PCM_16) CS(PCM_24) CS(PCM_32) CS(PCM_U8)
	CS(FLOAT) CS(DOUBLE) CS(ULAW) CS(ALAW)
	default:;
	}
	#undef CS
	return *this;
}

const char * SoundFile::toString(Format v){
	#define CS(x) case x: return #x;
	switch(v){
	CS(WAV) CS(AIFF) CS(AU) CS(RAW)
	//CS(FLAC)
	default: return "";
	}
	#undef CS
}
const char * SoundFile::toString(EncodingType v){
	#define CS(x) case x: return #x;
	switch(v){
	CS(PCM_S8) CS(PCM_16) CS(PCM_24) CS(PCM_32) CS(PCM_U8)
	CS(FLOAT) CS(DOUBLE) CS(ULAW) CS(ALAW)
	default: return "";
	}
	#undef CS
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
//bool SoundFile::instrument(){
//	bool hasInst = SF_TRUE == sf_command(fp, SFC_GET_INSTRUMENT, &inst, sizeof (inst));
//	if(hasInst){
//		printf("Gain:      %d\n", inst.gain);
//		printf("Base Note: %d\n", inst.basenote);
//		printf("Detune:    %d\n", inst.detune);
//		printf("Velocity:  [%d, %d]\n", inst.velocity_lo, inst.velocity_hi);
//		printf("Key:       [%d, %d]\n", inst.key_lo, inst.key_hi);
//		printf("Loops:     %d\n", inst.loop_count);
//		for(int i=0; i<inst.loop_count; i++){
//			printf("\t");
//			switch(inst.loops[i].mode){
//				case SF_LOOP_FORWARD:     printf("fwd "); break;
//				case SF_LOOP_BACKWARD:    printf("bwd "); break;
//				case SF_LOOP_ALTERNATING: printf("alt "); break;
//				default:                  printf("none");
//			}
//			printf(" [%d, %d] x %d\n", inst.loops[i].start, inst.loops[i].end, inst.loops[i].count);
//		}
//	}
//	return hasInst;
//}


bool SoundFile::openRead(){ return mImpl->openRead(path()); }
bool SoundFile::openRead(const std::string& path_){ path(path_); return openRead(); }

bool SoundFile::openWrite(){ return mImpl->openWrite(path()); }
bool SoundFile::openWrite(const std::string& path_){ path(path_); return openWrite(); }


void SoundFile::print(){
//	printf("Path:       %s\n", mPath.c_str());
	printf("%s (%s): %g frames/sec, %d chan, %d frames, %f sec\n",
		toString(format()), toString(encoding()), frameRate(), channels(), frames(), frames()/frameRate());
}

#define DEF_SPECIAL(type)\
	template<>\
	int SoundFile::read<type>(type * dst, int numFrames){\
		return mImpl->read(dst, numFrames);\
	}\
	template<>\
	int SoundFile::write<type>(const type * src, int numFrames){\
		return mImpl->write(src, numFrames);\
	}
	DEF_SPECIAL(float)
	DEF_SPECIAL(short)
	DEF_SPECIAL(int)
	DEF_SPECIAL(double)
#undef DEF_SPECIAL

void SoundFile::seek(int pos, int mode){
	mImpl->seek(pos, mode);
}

} // gam::
