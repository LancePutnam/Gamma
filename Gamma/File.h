#ifndef GAMMA_FILE_H_INC
#define GAMMA_FILE_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information

	File Description:
	Basic file i/o classes.
*/

#include <stdio.h>
#include "Gamma/pstdint.h"

namespace gam{

class File{
public:

	/// \param[in] path		path of file
	/// \param[in] mode		i/o mode w, r, wb, rb
	File(const char * path, const char * mode, bool open=false);

	~File();

	void close();	///< Close file
	bool open();	///< Open file with specified i/o mode

	void mode(const char * v){ mMode=v; }

	/// Write memory elements to file
	uint32_t write(const void * v, int size, int items=1){ return fwrite(v, size, items, mFP); }

	/// Quick and dirty write memory to file
	static uint32_t write(const char * path, const void * v, int size, int items=1);

	/// Returns character string of file contents (read mode only)
	char * readAll();

	/// Returns whether file is open
	bool opened() const { return 0 != mFP; }
	
	/// Returns file i/o mode string
	const char * mode() const { return mMode; }
	
	/// Returns path string
	const char * path() const { return mPath; }
	
	/// Returns size (in bytes) of file contents
	int size() const { return mSizeBytes; }

	/// Returns whether file exists
	static bool exists(const char * path);

protected:
	const char * mPath;
	const char * mMode;
	char * mContent;
	int mSizeBytes;
	FILE * mFP;
	
	void freeContent();
	void allocContent(int n);
	void getSize();
};

} // gam::

#endif

