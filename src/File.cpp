#include "Gamma/File.h"

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

} // gam::
