#ifndef GAMMA_FILE_H_INC
#define GAMMA_FILE_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information

	File Description:
	Basic file i/o classes.
*/

#include <stdio.h>
#include "Gamma/pstdint.h"

#define TEM template<class T>
#define ULONG unsigned long
#define UI4 unsigned long
#define UI8 unsigned long long

namespace gam{

class File{
public:

	/// @param[in] path		path of file
	/// @param[in] mode		i/o mode w, r, wb, rb
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



/// Generic cross-endian file format for numerical data.

/// The file format is a mix between binary and human-readable text.  Numbers are
/// written in ascii in hexidecimal form to avoid messy endian issues.  At the
/// time, only 4- and 8-byte data elements are supported.
///
/// Example:
///		struct Point3D{ float x, y, z; }
///
///		Point3D p3d;
///		DataFile file("positions.dat");
///		file.openWrite();
///		file.write4(&p3d, 3);
///		file.write4(p3d.x); file.write4(p3d.y); file.write4(p3d.z);
///		file.close();
class DataFile{
public:
	DataFile(const char * path);
	virtual ~DataFile();
	
	void close();		///< Close file
	bool openWrite();	///< Open file for writing
	bool openRead();	///< Open file for reading
	

	int read4(void * element);								///< Read 4-byte element from file
	int read8(void * element);								///< Read 8-byte element from file
	template<ULONG N> void read(void * element);			///< Read N-byte element from file
	void read4(void * dst, ULONG len);						///< Read 'len' 4-byte elements from file
	void read8(void * dst, ULONG len);						///< Read 'len' 8-byte elements from file
	template<ULONG N> void read(void * dst, ULONG len);		///< Read 'len' N-byte elements from file

	void write4(const void * element);						///< Write 4-byte element to file
	void write8(const void * element);						///< Write 8-byte element to file
	template<ULONG N> void write(const void * element);		///< Write N-byte element to file
	void write4(const void * src, ULONG len);				///< Write 'len' 4-byte elements to file
	void write8(const void * src, ULONG len);				///< Write 'len' 8-byte elements to file
	template<ULONG N> void write(const void * src, ULONG len); ///< Write 'len' N-byte elements to file

	const char * path() const { return mPath; }
	void path(const char * p){ mPath=p; }

	//TEM void read(T * dst, ULONG len);
	//TEM void write(const T * src, ULONG len);

private:
	const char * mPath;
	FILE * fp;
	
	TEM const char * typeString();
	TEM const char * formatString();
};


//TEM void DataFile::write(const T * src, ULONG len){
//	fprintf(fp, formatString<double>(), (UI8)len);
//	fwrite(typeString<T>(), 3, 1, fp);
//	fwrite("\r\n", 2, 1, fp);
//	write8(src, len);	
//	fwrite("\r\n", 2, 1, fp);
//}
//
//

#define DEFINE_SPECIAL(numBytes)\
	template<> \
	inline void DataFile::read<numBytes>(void * element){\
		read##numBytes(element);\
	}\
	template<> \
	inline void DataFile::read<numBytes>(void * dst, ULONG len){\
		read##numBytes(dst, len);\
	}\
	template<>\
	inline void DataFile::write<numBytes>(const void * element){\
		write##numBytes(element);\
	}\
	template<>\
	inline void DataFile::write<numBytes>(const void * src, ULONG len){\
		write##numBytes(src, len);\
	}

	DEFINE_SPECIAL(4)
	DEFINE_SPECIAL(8)
	
#undef DEFINE_SPECIAL

#define DEFINE_SPECIAL(type, strType, strFormat) \
	template<> \
	inline const char * DataFile::typeString<type>(){\
		return strType;\
	}\
	template<> \
	inline const char * DataFile::formatString<type>(){\
		return strFormat;\
	}

	DEFINE_SPECIAL(float,  "f04", "%08x "   )
	DEFINE_SPECIAL(double, "f08", "%016llx ")
	//DEFINE_SPECIAL(unsigned long, u04)
	//DEFINE_SPECIAL(long, s04)
#undef DEFINE_SPECIAL

} // gam::

#undef ULONG
#undef UI4
#undef UI8
#undef TEM

#endif

