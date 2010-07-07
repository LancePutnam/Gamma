#ifndef GAMMA_SOUNDFILE_H_INC
#define GAMMA_SOUNDFILE_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include <string>
#include <sndfile.h>
#include "Gamma/mem.h"
#include "Gamma/Types.h"

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
	SoundFile(std::string path="");
	
	/// Creates object given a path and an other object from which to get its header info.
	SoundFile(std::string path, const SoundFile& src);
	
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
	TEM uint32_t read(T * dst, uint32_t numFrames);
	
	/// Copy all contents of file into array interleaved. Returns number of frames read.
	TEM uint32_t readAll(T * dst);

	/// Copy all contents of file into array deinterleaved. Returns number of frames read.
	
	/// If the number of channels is > 1, memory will be dynamically allocated
	///	and freed for the deinterleaving.
	TEM uint32_t readAllD(T * dst);
	
	/// Writes interleaved frames from array to file.
	
	/// From the libsndfile docs:\n
	/// The file write frames functions write the data in the array pointed to
	/// by ptr to the file. The array must be large enough to hold the product
	/// of frames and the number of channels.
	TEM uint32_t write(const T * src, uint32_t numFrames);

	// Sound file properties
	double frameRate() const;			///< Returns frames/second
	uint32_t frames() const;			///< Returns number of frames
	uint32_t channels() const;			///< Returns number of channels
	uint32_t samples() const;			///< Returns number of samples ( = frames x channels)
	int format() const;					///< Returns format field
	int formatMajor() const;
	int formatMinor() const;
	const char * extension();			///< Returns file extension
	std::string path() const;			///< Returns path of sound file
	
	void channels(uint32_t num);		///< Set number of channels
	void frameRate(double hz);			///< Set frames/second
	void format(int newFormat);			///< Set format field
	void formatMajor(int major);		///< Set major format field
	void formatMinor(int minor);		///< Set minor format field
	void info(const SoundFile& src);	///< Copy file information from an other file
	void path(const char * path);		///< Set path of sound file
	void path(std::string path);		///< Set path of sound file
	
	/// Gets instrument data from file.
	/// Returns whether or not it found the instrument data.
	bool instrument();
	
	void print();			///< Print information about file to stdout.
	
	const char * errorString() const;
	
	static int formatMajor(std::string sfpath);
	
private:
	SNDFILE * fp;
	SF_INFO mInfo;
	SF_INSTRUMENT inst;
	SF_FORMAT_INFO formatInfo;
	
	std::string mPath;

	void formatInfoMajor();
	void formatInfoSubtype();
	
//	TEM static uint32_t read(SoundFile * sf, T * dst, uint32_t numFrames);
//	TEM static uint32_t readAll(SoundFile * sf, T * dst);
//	TEM static uint32_t readAllD(SoundFile * sf, T * dst);
//	
//	TEM static uint32_t write(SoundFile * sf, const T * src, uint32_t numFrames); 
};




// Implementation_______________________________________________________________

inline double SoundFile::frameRate() const { return (double)mInfo.samplerate; }
inline uint32_t SoundFile::frames() const { return (uint32_t)mInfo.frames; }
inline uint32_t SoundFile::channels() const { return (uint32_t)mInfo.channels; }
inline uint32_t SoundFile::samples() const { return frames() * channels(); }
inline int SoundFile::format() const { return mInfo.format; }
inline int SoundFile::formatMajor() const { return mInfo.format & SF_FORMAT_TYPEMASK; }
inline int SoundFile::formatMinor() const { return mInfo.format & SF_FORMAT_SUBMASK; }

inline const char * SoundFile::errorString() const { return sf_strerror(fp); }
inline std::string SoundFile::path() const { return mPath; }

inline void SoundFile::channels(uint32_t num){ mInfo.channels = (int)num; }
inline void SoundFile::frameRate(double hz){ mInfo.samplerate = (int)hz; }
inline void SoundFile::format(int newFormat){ mInfo.format = newFormat; }

inline void SoundFile::formatMajor(int major){
	mInfo.format = (mInfo.format & (~SF_FORMAT_TYPEMASK)) | (major & SF_FORMAT_TYPEMASK);
}

inline void SoundFile::formatMinor(int minor){
	mInfo.format = (mInfo.format & (~SF_FORMAT_SUBMASK)) | (minor & SF_FORMAT_SUBMASK);
}

inline void SoundFile::path(const char * p){ mPath = p; }
inline void SoundFile::path(std::string p){ mPath = p; }

// specialized templates to hook into libsndfile functions
#define DEFINE_SPECIAL(type) \
	template<> \
	inline uint32_t SoundFile::read<type>(type * dst, uint32_t numFrames){\
		return sf_readf_##type(fp, dst, numFrames);\
	}\
	template<> \
	inline uint32_t SoundFile::write<type>(const type * src, uint32_t numFrames){\
		return sf_writef_##type(fp, src, numFrames);\
	}
	DEFINE_SPECIAL(float)
	DEFINE_SPECIAL(short)
	DEFINE_SPECIAL(int)
	DEFINE_SPECIAL(double)
#undef DEFINE_SPECIAL

TEM inline uint32_t SoundFile::readAll(T * dst){
	sf_seek(fp, 0, SEEK_SET);
	return read(dst, frames());
}

TEM uint32_t SoundFile::readAllD(T * dst){
	uint32_t numChannels = channels();

	if(1 == numChannels) return readAll(dst);
	
	// Allocate memory for deinterleaving.  Don't know of any in-place methods.
	T * temp = new T[samples()];
	uint32_t framesRead = 0;
	
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

