#include "DataFile.h"

#define ULONG unsigned long
#define UI4 unsigned long
#define UI8 unsigned long long

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

void DataFile::read4(void * element){
	fscanf(fp, formatString<float >(), element);
}

void DataFile::read8(void * element){
	fscanf(fp, formatString<double>(), element);
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

#undef ULONG
#undef UI4
#undef UI8

