#ifndef GAMMA_SOUNDFILE_H_INC
#define GAMMA_SOUNDFILE_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include <string>
#include "sndfile.h"
#include "Gamma/mem.h"
//#include "Gamma/Types.h"

#define TEM template<class T>

namespace gam{
	
	
/// Class for reading and writing sound file data.

/// This is a wrapper around Erik de Castro Lopo's excellent libsndfile 
/// http://www.mega-nerd.com/libsndfile/ .
/// Supports WAV, AIFF, AU, RAW, PAF, SVX, NIST, VOC, IRCAM, W64, MAT4
/// MAT5, PVF, XI, HTK, SDS, AVR, WAVEX, SD2, FLAC, CAF.\n
/// The read() and write() methods are specialized template methods and will
/// result in a compilation error if the data type is not supported by
/// libsndfile.  At the moment, the following types of data are supported:
/// short, int, float, and double.
class SoundFile{
public:

	/// Creates object given a path.
	
	/// The sound file info structure will be zero.
	///
	SoundFile(const std::string& path="");
	
	/// Creates object given a path and an other object from which to get its header info.
	SoundFile(const std::string& path, const SoundFile& src);
	
	/// The destructor will automatically close any open files.
	~SoundFile();
	
	/// Opens sound file for reading.
	
	/// Returns true on success and false otherwise.
	///
	bool openRead();
	
	/// Opens sound file for writing.
	
	/// Before calling this method it is necessary to set the format of the
	/// file using format(), e.g. format(SF_FORMAT_AIFF | SF_FORMAT_PCM_16).
	/// Returns true on success and false otherwise.
	bool openWrite();
	
	bool close();		///< Closes sound file.  Files are closed in destructor.
	
	/// Reads next chunk of frames from file into array.
	
	/// From the libsndfile docs:\n
	/// The file read frames functions fill the array pointed to by ptr with
	/// the requested number of frames of data. The array must be large enough
	/// to hold the product of frames and the number of channels.
	//ULONG read(float * dst, ULONG numFrames);
	TEM int read(T * dst, int numFrames);
	
	/// Copy all contents of file into array interleaved. Returns number of frames read.
	TEM int readAll(T * dst);

	/// Copy all contents of file into array deinterleaved. Returns number of frames read.
	
	/// If the number of channels is > 1, memory will be dynamically allocated
	///	and freed for the deinterleaving.
	TEM int readAllD(T * dst);
	
	/// Writes interleaved frames from array to file.
	
	/// From the libsndfile docs:\n
	/// The file write frames functions write the data in the array pointed to
	/// by ptr to the file. The array must be large enough to hold the product
	/// of frames and the number of channels.
	TEM int write(const T * src, int numFrames);

	// Sound file properties
	double frameRate() const;			///< Returns frames/second
	int frames() const;					///< Returns number of frames
	int channels() const;				///< Returns number of channels
	int samples() const;				///< Returns number of samples ( = frames x channels)
	int format() const;					///< Returns format field
	int formatMajor() const;
	int formatMinor() const;
	const char * extension();			///< Returns file extension
	const std::string& path() const;	///< Returns path of sound file
	
	void channels(int num);				///< Set number of channels
	void frameRate(double hz);			///< Set frames/second
	void format(int newFormat);			///< Set format field
	void formatMajor(int major);		///< Set major format field
	void formatMinor(int minor);		///< Set minor format field
	void info(const SoundFile& src);	///< Copy file information from an other file
	void path(const std::string& path);	///< Set path of sound file
	
	/// Gets instrument data from file.
	/// Returns whether or not it found the instrument data.
	bool instrument();
	
	void print();			///< Print information about file to stdout.
	
	const char * errorString() const;
	
	static int formatMajor(const std::string& sfpath);
	
private:
	SNDFILE * fp;
	SF_INFO mInfo;
	SF_INSTRUMENT inst;
	SF_FORMAT_INFO formatInfo;
	
	std::string mPath;

	void formatInfoMajor();
	void formatInfoSubtype();
};




// Implementation_______________________________________________________________

inline double SoundFile::frameRate() const { return (double)mInfo.samplerate; }
inline int SoundFile::frames() const { return (int)mInfo.frames; }
inline int SoundFile::channels() const { return (int)mInfo.channels; }
inline int SoundFile::samples() const { return frames() * channels(); }
inline int SoundFile::format() const { return mInfo.format; }
inline int SoundFile::formatMajor() const { return mInfo.format & SF_FORMAT_TYPEMASK; }
inline int SoundFile::formatMinor() const { return mInfo.format & SF_FORMAT_SUBMASK; }

inline const char * SoundFile::errorString() const { return sf_strerror(fp); }
inline const std::string& SoundFile::path() const { return mPath; }

inline void SoundFile::channels(int num){ mInfo.channels = (int)num; }
inline void SoundFile::frameRate(double hz){ mInfo.samplerate = (int)hz; }
inline void SoundFile::format(int newFormat){ mInfo.format = newFormat; }
inline void SoundFile::path(const std::string& p){ mPath = p; }

// specialized templates to hook into libsndfile functions
#define DEFINE_SPECIAL(type) \
	template<> \
	inline int SoundFile::read<type>(type * dst, int numFrames){\
		return sf_readf_##type(fp, dst, numFrames);\
	}\
	template<> \
	inline int SoundFile::write<type>(const type * src, int numFrames){\
		return sf_writef_##type(fp, src, numFrames);\
	}
	DEFINE_SPECIAL(float)
	DEFINE_SPECIAL(short)
	DEFINE_SPECIAL(int)
	DEFINE_SPECIAL(double)
#undef DEFINE_SPECIAL

TEM inline int SoundFile::readAll(T * dst){
	sf_seek(fp, 0, SEEK_SET);
	return read(dst, frames());
}

TEM int SoundFile::readAllD(T * dst){
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

} // gam::

#undef TEM
#endif

