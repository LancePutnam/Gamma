#include "Gamma/File.h"

#define ULONG unsigned long
#define UI4 unsigned long
#define UI8 unsigned long long

namespace gam{

File::File(const char * path, const char * mode, bool open_)
:	mPath(path), mMode(mode), mContent(0), mSizeBytes(0), mFP(0)
{	if(open_) open(); }

File::~File(){ close(); freeContent(); }

void File::allocContent(int n){
	if(mContent) freeContent();
	mContent = new char[n+1];
	mContent[n] = '\0';
}

void File::close(){ if(opened()){ fclose(mFP); mFP=0; } }

void File::freeContent(){ delete[] mContent; }

void File::getSize(){
	int r=0;
	if(opened()){
		fseek(mFP, 0, SEEK_END);
		r = ftell(mFP);
		rewind(mFP);
	}
	mSizeBytes = r;
}

bool File::open(){
	if(0 == mFP){
		if((mFP = fopen(mPath, mMode))){
			getSize();
			return true;
		}
	}
	return false;
}

char * File::readAll(){
	if(opened() && mMode[0]=='r'){
		int n = size();
		allocContent(n);
		int numRead = fread(mContent, sizeof(char), n, mFP);
		if(numRead < n){}
	}
	return mContent;
}

uint32_t File::write(const char * path, const void * v, int size, int items){
	File f(path, "w");
	uint32_t r = 0;
	if(f.open()){
		r = f.write(v, size, items);
		f.close();
	}
	return r;
}


bool File::exists(const char * path){ File f(path, "r"); return f.open(); }




DataFile::DataFile(const char * path)
	: mPath(path), fp(0)
{
}

DataFile::~DataFile(){
	close();
}

void DataFile::close(){
	if(fp){
		fclose(fp);
		fp = 0;
	}
}

bool DataFile::openWrite(){
	if(0 == fp){
		fp = fopen(mPath, "wb");
		if(fp) return true;
	}
	return false;
}

bool DataFile::openRead(){
	if(0 == fp){
		fp = fopen(mPath, "rb");
		if(fp) return true;
	}
	return false;
}

int DataFile::read4(void * element){
	return fscanf(fp, formatString<float >(), element);
}

int DataFile::read8(void * element){
	return fscanf(fp, formatString<double>(), element);
}

void DataFile::read4(void * dst, ULONG len){
	UI4 * d = (UI4 *)dst;
	for(ULONG i=0; i<len; i++) read4(d + i);
}

void DataFile::read8(void * dst, ULONG len){
	UI8 * d = (UI8 *)dst;
	for(ULONG i=0; i<len; i++) read8(d + i);
}

void DataFile::write4(const void * element){
	fprintf(fp, formatString<float >(), *(UI4 *)(element));
}

void DataFile::write8(const void * element){
	fprintf(fp, formatString<double>(), *(UI8 *)(element));
}

void DataFile::write4(const void * src, ULONG len){
	UI4 * s = (UI4 *)src;
	for(ULONG i=0; i<len; i++) write4(s + i);
}

void DataFile::write8(const void * src, ULONG len){
	UI8 * s = (UI8 *)src;
	for(ULONG i=0; i<len; i++) write8(s + i);
}

} // gam::

#undef ULONG
#undef UI4
#undef UI8

