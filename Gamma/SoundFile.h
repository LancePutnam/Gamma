#ifndef GAMMA_SOUNDFILE_H_INC
#define GAMMA_SOUNDFILE_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include <string>
#include <vector>
#include <stdio.h>
#include "Gamma/mem.h"

namespace gam{
	
	
/// Class for reading and writing sound file data
class SoundFile{
public:

	/// Sound file formats
	enum Format{
		WAV = 1,	/**< Microsoft WAV format (little endian default). */
		AIFF,		/**< Apple/SGI AIFF format (big endian). */
		AU,			/**< Sun/NeXT AU format (big endian). */
		RAW,		/**< RAW PCM data. */
		NO_FORMAT
	};

	/// Sound file sample encoding types
	enum EncodingType{
		PCM_S8 = 1,	/**< Signed 8 bit data */
		PCM_16,		/**< Signed 16 bit data */
		PCM_24,		/**< Signed 24 bit data */
		PCM_32,		/**< Signed 32 bit data */
		PCM_U8,		/**< Unsigned 8 bit data (WAV and RAW only) */

		FLOAT,		/**< 32 bit float data */
		DOUBLE,		/**< 64 bit float data */

		ULAW,		/**< U-Law encoded. */
		ALAW,		/**< A-Law encoded. */
	};

	static const char * toString(Format v);
	static const char * toString(EncodingType v);


	/// Creates object given a path.

	/// The sound file info structure will be zero.
	///
	SoundFile(const std::string& path="");

	/// Creates object given a path and an other object from which to get its header info.
	SoundFile(const std::string& path, const SoundFile& src);

	/// The destructor will automatically close the file if it's open
	~SoundFile();


	/// Open sound file for reading

	/// Returns true on success and false otherwise.
	///
	bool openRead();
	bool openRead(const std::string& path);

	/// Open sound file for writing

	/// Before calling this method it is necessary to set the number of 
	/// channels and frame rate. If no EncodingType is set, then PCM_16 is 
	/// used. The Format of the sound file is derived from the extension.
	/// \returns true on success and false otherwise.
	bool openWrite();
	bool openWrite(const std::string& path);

	/// Close sound file

	/// \returns whether the sound file was closed successfully.
	///
	bool close();

	/// Read next chunk of frames from file into array

	/// From the libsndfile docs:\n
	/// The file read frames functions fill the array pointed to by ptr with
	/// the requested number of frames of data. The array must be large enough
	/// to hold the product of frames and the number of channels.
	//ULONG read(float * dst, ULONG numFrames);
	template<class T>
	int read(T * dst, int numFrames);

	/// Copy all contents of file into array interleaved. Returns number of frames read.
	template<class T>
	int readAll(T * dst);

	/// Copy all contents of file into array interleaved. Returns number of frames read.
	template<class T>
	int readAll(std::vector<T>& dst);

	/// Copy all contents of file into array deinterleaved. Returns number of frames read.

	/// If the number of channels is > 1, memory will be dynamically allocated
	///	and freed for the deinterleaving.
	template<class T>
	int readAllD(T * dst);

	/// Write interleaved frames from array to file

	/// From the libsndfile docs:\n
	/// The file write frames functions write the data in the array pointed to
	/// by ptr to the file. The array must be large enough to hold the product
	/// of frames and the number of channels.
	template<class T>
	int write(const T * src, int numFrames);

	template<class T>
	int write(const std::vector<T>& src);

	// Sound file properties
	bool opened() const;						///< Returns whether the sound file is open
	EncodingType encoding() const;				///< Get encoding type
	Format format() const;						///< Get format
	double frameRate() const;					///< Get frames/second
	int frames() const;							///< Get number of frames
	int channels() const;						///< Get number of channels
	int samples() const;						///< Get number of samples ( = frames x channels)
	const char * extension();					///< Get file extension
	const std::string& path() const;			///< Get path of sound file

	SoundFile& encoding(EncodingType v);		///< Set encoding type
	SoundFile& format(Format v);				///< Set format
	SoundFile& channels(int num);				///< Set number of channels
	SoundFile& frameRate(double hz);			///< Set frames/second
	SoundFile& info(const SoundFile& src);		///< Copy file information from an other file
	SoundFile& path(const std::string& path);	///< Set path of sound file

	void seek(int pos, int seekMode);

	void print();			///< Print information about file to stdout.

private:
	class Impl; Impl * mImpl;
	std::string mPath;
};




// Implementation_______________________________________________________________
template<class T>
int SoundFile::readAll(T * dst){
	seek(0, SEEK_SET);
	return read(dst, frames());
}

template<class T>
int SoundFile::readAll(std::vector<T>& dst){
	dst.resize(samples());
	return readAll(&dst[0]);
}

template<class T>
int SoundFile::readAllD(T * dst){
	int numChannels = channels();

	if(1 == numChannels) return readAll(dst);
	
	// Allocate memory for deinterleaving.  Don't know of any in-place methods.
	T * temp = new T[samples()];
	int framesRead = 0;
	
	if(temp){
		framesRead = readAll(temp);
		if(framesRead == frames()){
			if(2 == numChannels)	mem::deinterleave2(dst, temp, frames());
			else					mem::deinterleave(dst, temp, frames(), numChannels);
		}
		delete[] temp;
	}
	return framesRead;
}

template<class T>
int SoundFile::write(const std::vector<T>& src){
	return write(&src[0], src.size()/channels());
}

} // gam::

#endif

