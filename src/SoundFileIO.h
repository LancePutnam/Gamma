/* Minimal sound file reader and writer.

Supported formats:
- Wave
- AIFF
- au/snd

Supported encodings:
- PCM 8-bit, 16-bit, 24-bit and 32-bit integer, 32-bit and 64-bit float
- u-law
- a-law

Copyright (c) 2019 Lance Putnam

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include <cmath>
#include <cstdint>
#include <fstream>
#include <vector>
#include <string>

#define DPRINTF(...)
//#define DPRINTF(...) printf(__VA_ARGS__);

struct int24_t{ int32_t v; };
struct uint24_t{ uint32_t v; };

template <class T> constexpr size_t numBytes(){ return sizeof(T); }
template <> constexpr size_t numBytes< int24_t>(){ return 3; }
template <> constexpr size_t numBytes<uint24_t>(){ return 3; }

constexpr uint32_t ID(const char * s){
	return ((uint32_t)s[0]<<24) | ((uint32_t)s[1]<<16) | ((uint32_t)s[2]<<8) | s[3];
};

template <class To, class From>
To pun(From v){ union{From f; To t;}; f=v; return t; }

template <class T> T toBE(const unsigned char * bytes);
template <class T> T toLE(const unsigned char * bytes);

template<> uint16_t toBE(const unsigned char * bytes){
	return	((uint16_t)bytes[0] <<  8) |
			((uint16_t)bytes[1] <<  0);
}
template<> uint16_t toLE(const unsigned char * bytes){
	return	((uint16_t)bytes[1] <<  8) |
			((uint16_t)bytes[0] <<  0);
}
template<> int16_t toBE(const unsigned char * bytes){
	return pun<int16_t>(toBE<uint16_t>(bytes));
}
template<> int16_t toLE(const unsigned char * bytes){
	return pun<int16_t>(toLE<uint16_t>(bytes));
}

template<> uint24_t toBE(const unsigned char * bytes){
	return uint24_t{
			((uint32_t)bytes[0] << 16) |
			((uint32_t)bytes[1] <<  8) |
			((uint32_t)bytes[2] <<  0)
	};
}
template<> uint24_t toLE(const unsigned char * bytes){
	return uint24_t{
			((uint32_t)bytes[2] << 16) |
			((uint32_t)bytes[1] <<  8) |
			((uint32_t)bytes[0] <<  0)
	};
}
template<> int24_t toBE(const unsigned char * bytes){
	return int24_t{pun<int24_t>(toBE<uint24_t>(bytes).v<<8).v>>8};
}
template<> int24_t toLE(const unsigned char * bytes){
	return int24_t{pun<int24_t>(toLE<uint24_t>(bytes).v<<8).v>>8};
}

template<> uint32_t toBE(const unsigned char * bytes){
	return	((uint32_t)bytes[0] << 24) |
			((uint32_t)bytes[1] << 16) |
			((uint32_t)bytes[2] <<  8) |
			((uint32_t)bytes[3] <<  0);
}
template<> uint32_t toLE(const unsigned char * bytes){
	return	((uint32_t)bytes[3] << 24) |
			((uint32_t)bytes[2] << 16) |
			((uint32_t)bytes[1] <<  8) |
			((uint32_t)bytes[0] <<  0);
}
template<> int32_t toBE(const unsigned char * bytes){
	return pun<int32_t>(toBE<uint32_t>(bytes));
}
template<> int32_t toLE(const unsigned char * bytes){
	return pun<int32_t>(toLE<uint32_t>(bytes));
}

template<> uint64_t toBE(const unsigned char * bytes){
	return	((uint64_t)bytes[0] << 56) |
			((uint64_t)bytes[1] << 48) |
			((uint64_t)bytes[2] << 40) |
			((uint64_t)bytes[3] << 32) |
			((uint64_t)bytes[4] << 24) |
			((uint64_t)bytes[5] << 16) |
			((uint64_t)bytes[6] <<  8) |
			((uint64_t)bytes[7] <<  0);
}
template<> uint64_t toLE(const unsigned char * bytes){
	return	((uint64_t)bytes[7] << 56) |
			((uint64_t)bytes[6] << 48) |
			((uint64_t)bytes[5] << 40) |
			((uint64_t)bytes[4] << 32) |
			((uint64_t)bytes[3] << 24) |
			((uint64_t)bytes[2] << 16) |
			((uint64_t)bytes[1] <<  8) |
			((uint64_t)bytes[0] <<  0);
}

template<> float toBE(const unsigned char * bytes){
	return pun<float>(toBE<uint32_t>(bytes));
}
template<> float toLE(const unsigned char * bytes){
	return pun<float>(toLE<uint32_t>(bytes));
}

template<> double toBE(const unsigned char * bytes){
	return pun<double>(toBE<uint64_t>(bytes));
}
template<> double toLE(const unsigned char * bytes){
	return pun<double>(toLE<uint64_t>(bytes));
}

template <class T>
T readBE(std::ifstream& f, char * buf){
	f.read(buf, sizeof(T));
	return toBE<T>((unsigned char *)buf);
}
template <class T> T readBE(std::ifstream& f){
	char buf[sizeof(T)];
	return readBE<T>(f,buf);
}
template <class T>
void readBE(std::ifstream& f, T * dst, unsigned len){
	f.read((char *)dst, sizeof(T)*len);
	for(unsigned i=0; i<len; ++i){
		dst[i] = toBE<T>(((unsigned char *)dst) + sizeof(T)*i);
	}
}

template <class T>
T readLE(std::ifstream& f){
	unsigned char buf[sizeof(T)];
	f.read((char *)buf, sizeof(T));
	return toLE<T>(buf);
}
template <class T>
void readLE(std::ifstream& f, T * dst, unsigned len){
	f.read((char *)dst, sizeof(T)*len);
	for(unsigned i=0; i<len; ++i){
		dst[i] = toLE<T>(((unsigned char *)dst) + sizeof(T)*i);
	}
}

template <class From, class To> To convert(From v);
template <class From, class To> To convert(From v){ return To(v); }
// Precautions:
// * Bit shift on signed integer may be undefined
// * Casting between signed/unsigned types can under/overflow
template<> float   convert<int8_t, float  >(int8_t  v){ return v / 128.f; }
template<> double  convert<int8_t, double >(int8_t  v){ return v / 128.; }
template<> int16_t convert<int8_t, int16_t>(int8_t  v){ return int16_t(v) * 256; }
template<> int24_t convert<int8_t, int24_t>(int8_t  v){ return{int32_t(v) * 65536}; }
template<> int32_t convert<int8_t, int32_t>(int8_t  v){ return int32_t(v) * 16777216; }
template<> float   convert<uint8_t,float  >(uint8_t v){ return v / 128.f - 1.f; }
template<> double  convert<uint8_t,double >(uint8_t v){ return v / 128.0 - 1.0; }
template<> int16_t convert<uint8_t,int16_t>(uint8_t v){ return (int16_t(v)-128) * 256; }
template<> int24_t convert<uint8_t,int24_t>(uint8_t v){ return{(int32_t(v)-128) * 65536}; }
template<> int32_t convert<uint8_t,int32_t>(uint8_t v){ return (int32_t(v)-128) * 16777216; }
template<> float   convert<int16_t,float  >(int16_t v){ return v / 32768.f; }
template<> double  convert<int16_t,double >(int16_t v){ return v / 32768.; }
template<> int8_t  convert<int16_t,int8_t >(int16_t v){ return v / 256; }
template<> uint8_t convert<int16_t,uint8_t>(int16_t v){ return v / 256 + 128; }
template<> int24_t convert<int16_t,int24_t>(int16_t v){ return{int32_t(v) * 256}; }
template<> int32_t convert<int16_t,int32_t>(int16_t v){ return int32_t(v) * 65536; }
template<> float   convert<int24_t,float  >(int24_t v){ return v.v / 8388608.f; }
template<> double  convert<int24_t,double >(int24_t v){ return v.v / 8388608.; }
template<> int8_t  convert<int24_t,int8_t >(int24_t v){ return v.v / 65536; }
template<> uint8_t convert<int24_t,uint8_t>(int24_t v){ return v.v / 65536 + 128; }
template<> int16_t convert<int24_t,int16_t>(int24_t v){ return v.v / 256; }
template<> int32_t convert<int24_t,int32_t>(int24_t v){ return v.v * 256; }
template<> float   convert<int32_t,float  >(int32_t v){ return v / 2147483648.; }
template<> double  convert<int32_t,double >(int32_t v){ return v / 2147483648.; }
template<> int8_t  convert<int32_t,int8_t >(int32_t v){ return v / 16777216; }
template<> uint8_t convert<int32_t,uint8_t>(int32_t v){ return v / 16777216 + 128; }
template<> int16_t convert<int32_t,int16_t>(int32_t v){ return v / 65536; }
template<> int24_t convert<int32_t,int24_t>(int32_t v){ return{int32_t(v) / 256}; }
template<> int8_t  convert<float,  int8_t >(float   v){ return v * 127.999f; }
template<> uint8_t convert<float,  uint8_t>(float   v){ return v * 127.999f + 128.f; }
template<> int16_t convert<float,  int16_t>(float   v){ return v * 32767.999f; }
template<> int24_t convert<float,  int24_t>(float   v){ return{int32_t(v * 8388607.999)}; }
template<> int32_t convert<float,  int32_t>(float   v){ return v * 2147483647.999; }
template<> int8_t  convert<double, int8_t >(double  v){ return v * 127.99999; }
template<> uint8_t convert<double, uint8_t>(double  v){ return v * 127.99999 + 128.; }
template<> int16_t convert<double, int16_t>(double  v){ return v * 32767.99999; }
template<> int24_t convert<double, int24_t>(double  v){ return{int32_t(v * 8388607.99999)}; }
template<> int32_t convert<double, int32_t>(double  v){ return v * 2147483647.99999; }

template <class From, class To>
To convertLE(const unsigned char * bytes){ return convert<From,To>(toLE<From>(bytes)); }

template <class From, class To>
To convertBE(const unsigned char * bytes){ return convert<From,To>(toBE<From>(bytes)); }


class SoundFileInfo{
public:

	/// Sound file formats
	enum Format{
		WAV = 0,	/**< Microsoft WAV format */
		AIFF,		/**< Apple/SGI AIFF format */
		AU,			/**< Sun/NeXT AU format */
		NO_FORMAT
	};

	/// Sound file sample encoding types
	enum Encoding{
		PCM_U8 = 0,	/**< Unsigned 8 bit data */
		PCM_S8,		/**< Signed 8 bit data */
		PCM_16,		/**< Signed 16 bit data */
		PCM_24,		/**< Signed 24 bit data */
		PCM_32,		/**< Signed 32 bit data */
		FLOAT,		/**< 32 bit float data */
		DOUBLE,		/**< 64 bit float data */
		ULAW,		/**< 8-bit ITU-T G.711 mu-law */
		ALAW,		/**< 8-bit ITU-T G.711 A-law */
		NO_ENCODING
	};

	static const std::string& extension(Format v){
		static std::string e[] = {".wav",".aiff",".au",""};
		static_assert(sizeof(e)/sizeof(e[0])==(NO_FORMAT+1), "");
		return e[v];
	}

	SoundFileInfo& channels(unsigned v){ mChannels=v; return *this; }
	SoundFileInfo& frameRate(unsigned v){ mFrameRate=v; return *this; }
	SoundFileInfo& encoding(Encoding v){ mEncoding=v; return *this; }
	SoundFileInfo& format(Format v){
		switch(v){
		case WAV: mDataBigEndian = false; break;
		case AIFF:
		case AU:  mDataBigEndian = true; break;
		default: return *this;
		}
		mFormat=v;
		return *this;
	}

	unsigned channels() const { return mChannels; }
	unsigned frameRate() const { return mFrameRate; }
	Format format() const { return mFormat; }
	const std::string& extension() const { return extension(format()); }
	Encoding encoding() const { return mEncoding; }
	unsigned frames() const { return mFrames; }
	unsigned samples() const { return mChannels*mFrames; }
	unsigned bytesPerSample() const {
		static unsigned char s[] = {1,1,2,3,4,4,8,1,1,1};
		static_assert(sizeof(s)/sizeof(s[0])==(NO_ENCODING+1), "");
		return s[mEncoding];
	}
	unsigned bytesPerFrame() const { return bytesPerSample() * channels(); }

	struct Chunk{
		std::string name;
		std::vector<unsigned char> data;
	};

	typedef std::vector<Chunk> Chunks;

	const Chunks& chunks() const { return mChunks; }

	struct Loop{
		double begin = 0.;
		double end = 0.;
		uint32_t count = 0; // 0 means infinite
		uint8_t type = 0; // 0:forward 1:pingpong 2:backward
		bool valid() const { return begin!=end; }
	};

	struct SamplerInfo{
		float note = 60; // MIDI note (middle C is 60)
		float gain = 1.;
		int8_t noteLo = 0;
		int8_t noteHi = 127;
		int8_t velLo = 1;
		int8_t velHi = 127;
		Loop sustainLoop;
		Loop releaseLoop;
	};

	typedef std::vector<SamplerInfo> SamplerInfos;

	const SamplerInfos& samplerInfos() const { return mSamplerInfos; }

	void print() const {
		static const char * e[]={"PCM_U8","PCM_S8","PCM_16","PCM_24","PCM_32","FLOAT","DOUBLE","ULAW","ALAW","NO_ENCODING"};
		static_assert(sizeof(e)/sizeof(e[0])==(NO_ENCODING+1), "");
		printf("%s, %u Hz, %u chan, %u frames\n", e[mEncoding], mFrameRate, mChannels, mFrames);
	}

protected:
	Format mFormat = NO_FORMAT;
	Encoding mEncoding = NO_ENCODING;
	unsigned mChannels=0, mFrameRate=1, mFrames=0;
	Chunks mChunks;
	SamplerInfos mSamplerInfos;
	bool mDataBigEndian = false;
};


class SoundFileReader : public SoundFileInfo {
public:

	SoundFileReader(){}

	SoundFileReader(const SoundFileInfo& info)
	:	SoundFileInfo(info)
	{}


	/// Open file for reading frames

	/// \returns true if the file opened successfully. The file pointer will
	/// be at the start of the first frame.
	bool open(const std::string& path){
		mFile.open(path, std::ios::in | std::ios::binary);
		if(!mFile.is_open()) return false;

		union{
			char bufs[32];
			unsigned char buf[32];
		};
		char chunkName[4];

		mReadFrame = 0;
		std::streampos filePosOfData = -1;
		const auto fileID = readBE<uint32_t>(mFile);

		// For chunked formats wav, aiff, etc.
		auto storeChunk = [&](uint32_t chunkSize, int copyMode){
			if(0==copyMode) return;
			mChunks.emplace_back();
			auto& chunk = mChunks.back();
			chunk.name.assign(chunkName, 4);
			chunk.data.resize(chunkSize);
			if(2==copyMode) mFile.seekg(-int(chunkSize), mFile.cur);
			mFile.read((char*)&chunk.data[0], chunk.data.size());
		};

		if(ID(".snd") == fileID){
			format(AU);
			filePosOfData= readBE<uint32_t>(mFile);
			auto dataSize= readBE<uint32_t>(mFile);
			auto enc     = readBE<uint32_t>(mFile);
			mFrameRate   = readBE<uint32_t>(mFile);
			mChannels    = readBE<uint32_t>(mFile);
			switch(enc){
				case  2: mEncoding = PCM_S8; break;
				case  3: mEncoding = PCM_16; break;
				case  4: mEncoding = PCM_24; break;
				case  5: mEncoding = PCM_32; break;
				case  6: mEncoding = FLOAT; break;
				case  7: mEncoding = DOUBLE; break;
				case  1: mEncoding = ULAW; break;
				case 27: mEncoding = ALAW; break;
				default: mEncoding = NO_ENCODING;
			}
			if(0xFFFFFFFF == dataSize){
				mFile.seekg(0, mFile.end);
				dataSize = mFile.tellg() - filePosOfData;
			}
			mFrames = dataSize/mChannels/bytesPerSample();

		} else if(ID("RIFF") == fileID){
			mFile.read(bufs, 4); // main chunk size
			DPRINTF("RIFF (%d bytes)\n", toLE<uint32_t>(buf));
			if(ID("WAVE") == readBE<uint32_t>(mFile)){
				DPRINTF("Found WAVE file\n");
				format(WAV);
				unsigned bps = 1; // bits/sample

				while(!mFile.eof()){
					auto chunkID = readBE<uint32_t>(mFile, chunkName);
					if(mFile.gcount() != 4) break;
					DPRINTF("%c%c%c%c ", chunkName[0], chunkName[1], chunkName[2], chunkName[3]);
					auto chunkSize = readLE<uint32_t>(mFile);
					DPRINTF("(%u bytes)\n", chunkSize);
					int chunkStoreFlag = 0;

					switch(chunkID){
					case ID("fmt "):{
						auto format    = readLE<uint16_t>(mFile);
						mChannels      = readLE<uint16_t>(mFile);
						mFrameRate     = readLE<uint32_t>(mFile);
						auto BR        = readLE<uint32_t>(mFile);
						auto blockAlign= readLE<uint16_t>(mFile);
						bps = readLE<uint16_t>(mFile);
						DPRINTF("  format:%u chans:%u SR:%u BR:%u blockAlign:%u bps:%u\n", format, mChannels, mFrameRate, BR, blockAlign, bps);
						if(1 == format){ // PCM
							switch(bps){
							case  8: mEncoding = PCM_U8; break;
							case 16: mEncoding = PCM_16; break;
							case 24: mEncoding = PCM_24; break;
							case 32: mEncoding = PCM_32; break;
							default: goto error;
							}
						} else if(3 == format){
							switch(bps){
							case 32: mEncoding = FLOAT; break;
							case 64: mEncoding = DOUBLE; break;
							default: goto error;
							}
						} else if(6 == format){
							mEncoding = ALAW;
						} else if(7 == format){
							mEncoding = ULAW;
						} else if(0xFFFE == format){ // wave format extensible
						} else {
							goto error;
						}
						if(chunkSize > 16){ // extension data in chunk
							auto extSize = readLE<uint16_t>(mFile);
							DPRINTF("  ext size: %d\n", extSize);
							// ignore ext data; use chunk size in case extSize=0
							mFile.seekg(chunkSize-16-2, mFile.cur);
						}
						chunkStoreFlag = 2;
					} break;

					case ID("data"):{
						mFrames = chunkSize/mChannels/(bps/8);
						filePosOfData = mFile.tellg();
						if(mStoreDataChunk) chunkStoreFlag = 1;
						else mFile.seekg(chunkSize, mFile.cur);
					} break;

					case ID("fact"):{
						auto numSamps = readLE<uint32_t>(mFile);
						DPRINTF("  fact (num samples): %u\n", numSamps);
						if(chunkSize > 4) mFile.seekg(chunkSize-4, mFile.cur);
						chunkStoreFlag = 2;
					} break;

					/*case ID("LIST"):{
						mFile.read(bufs, 4);
						DPRINTF("%c%c%c%c\n", bufs[0], bufs[1], bufs[2], bufs[3]);
						// list of subchunks follows, can just parse as normal chunks
					} break;*/

					case ID("smpl"):{
						uint32_t d[9];
						readLE(mFile, d,9); // 0:manf 1:prod 2:period 3:note 4:semitone up 5:smpteFmt 6:smpteOff 7:loops 8:sampler data size
						DPRINTF("  manf:%u prod:%u period:%u note:%u semi:%u loops:%u data:%u\n", d[0],d[1],d[2],d[3],d[4],d[7],d[8]);
						for(unsigned i=0; i<d[7]; ++i){
							uint32_t l[6]; // 0:cue ID 1:type 2:beg 3:end 4:frac 5:count
							readLE(mFile, l,6);
							DPRINTF("  ID:%u type:%u beg:%u end:%u frac:%u cnt:%u\n", l[0],l[1],l[2],l[3],l[4],l[5]);
							SamplerInfo I;
							I.note = d[3] + d[4]/4294967296.;
							I.sustainLoop.type  = l[1];
							I.sustainLoop.begin = l[2];
							I.sustainLoop.end   = l[3] + l[4]/4294967296.;
							I.sustainLoop.count = l[5];
							mSamplerInfos.push_back(I);
						}
						if(d[8]) mFile.seekg(chunkSize-(36+d[7]*24), mFile.cur);
						chunkStoreFlag = 2;
					} break;

					case ID("inst"):{
						char d[7]; // 0:note 1:cents 2:gain(dB) 3:low note 4:high note 5:low vel 6:high vel
						mFile.read(d,7);
						SamplerInfo I;
						I.note = d[0] + d[1]/100.;
						I.gain = std::pow(10., d[2]/20.);
						I.noteLo = d[3];
						I.noteHi = d[4];
						I.velLo = d[5];
						I.velHi = d[6];
						mSamplerInfos.push_back(I);
						chunkStoreFlag = 2;
					} break;

					default: // unhandled chunk
						chunkStoreFlag = 1;
					}

					storeChunk(chunkSize, chunkStoreFlag);
				}
			}

		} else if(ID("FORM") == fileID){
			mFile.read(bufs, 4); // main chunk size
			auto aiffID = readBE<uint32_t>(mFile);
			if(ID("AIFF") == aiffID || ID("AIFC") == aiffID){
				const bool isAIFC = (aiffID == ID("AIFC"));
				DPRINTF("Found %s file\n", isAIFC ? "AIFC" : "AIFF");
				format(AIFF);
				unsigned bps = 1; // bits/sample

				auto getFloat10BE = [](const unsigned char * bytes){
					double f;
					int32_t expo= toBE<uint16_t>(bytes + 0) & 0x7FFF;
					auto mantHi = toBE<uint32_t>(bytes + 2);
					auto mantLo = toBE<uint32_t>(bytes + 6);
					if(expo==0 && mantHi==0 && mantLo==0) f=0.;
					else{
						if(expo == 0x7FFF) f = HUGE_VAL;
						else{
							expo -= 16383;
							auto uintToFloat = [](uint32_t u){ return (((double)((int32_t)(u - 2147483647L - 1))) + 2147483648.0); };
							f  = ldexp(uintToFloat(mantHi), expo-=31);
							f += ldexp(uintToFloat(mantLo), expo-=32);
						}
					}
					return bytes[0] & 0x80 ? -f : f;
				};

				while(!mFile.eof()){
					auto chunkID = readBE<uint32_t>(mFile, chunkName);
					if(mFile.gcount() != 4) break;
					DPRINTF("%c%c%c%c ", chunkName[0], chunkName[1], chunkName[2], chunkName[3]);
					auto chunkSize = readBE<uint32_t>(mFile);
					DPRINTF("(%u bytes)\n", chunkSize);
					int chunkStoreFlag = 0;

					switch(chunkID){
					case ID("COMM"):{
						mChannels = readBE<uint16_t>(mFile);
						mFrames   = readBE<uint32_t>(mFile);
						bps       = readBE<uint16_t>(mFile);
						mFile.read(bufs,10); mFrameRate = getFloat10BE(buf);
						mEncoding = NO_ENCODING;
						if(isAIFC){
							auto compID = readBE<uint32_t>(mFile, bufs);
							DPRINTF("%c%c%c%c\n", bufs[0], bufs[1], bufs[2], bufs[3]);
							if(chunkSize > 22) mFile.seekg(chunkSize-22, mFile.cur); // remaining bytes are compression name string
							switch(compID){
							case ID("sowt"): mDataBigEndian = false; break;
							case ID("fl32"):
							case ID("FL32"): mEncoding = FLOAT; break;
							case ID("fl64"):
							case ID("FL64"): mEncoding = DOUBLE; break;
							case ID("alaw"): mEncoding = ALAW; break;
							case ID("ulaw"): mEncoding = ULAW; break;
							case ID("raw "): mEncoding = PCM_U8; break;
							default: goto error; // unsupported compression
							}
						}
						if(mEncoding==NO_ENCODING){
							switch(bps){
							case  8: mEncoding = PCM_S8; break;
							case 16: mEncoding = PCM_16; break;
							case 24: mEncoding = PCM_24; break;
							case 32: mEncoding = PCM_32; break;
							default: goto error; // unsupported PCM bit depth
							}
						}
						DPRINTF("  chans:%u SR:%u frames:%u bps:%u\n", mChannels, mFrameRate, mFrames, bps);
						chunkStoreFlag = 2;
					} break;

					case ID("SSND"):{
						auto offset    = readBE<uint32_t>(mFile);
						auto blockSize = readBE<uint32_t>(mFile);
						DPRINTF("  offset:%u blockSize:%u\n", offset, blockSize);
						mFile.read(bufs, offset);
						filePosOfData = mFile.tellg();
						mFile.seekg(chunkSize - (8+offset), mFile.cur);
						if(mStoreDataChunk) chunkStoreFlag = 2;
					} break;

					case ID("INST"):{
						// http://www-mmsp.ece.mcgill.ca/Documents/AudioFormats/AIFF/Docs/AIFF-1.3.pdf
						SamplerInfo I;
						{	char d[6]; // 0:note 1:cents 2:low note 3:high note 4:low vel 5:high vel
							mFile.read(d,6);
							I.note = d[0] + d[1]/100.;
							I.noteLo = d[2];
							I.noteHi = d[3];
							I.velLo = d[4];
							I.velHi = d[5];
						}
						{	int16_t d[7]; // 0:gain 1:sus dir 2:sus mark beg 3:sus mark end 4:rel dir 5:rel mark beg 6:rel mark end
							readBE(mFile, d,7);
							I.gain = std::pow(10., d[0]/20.);
							I.sustainLoop.type = d[1];
							I.sustainLoop.begin = d[2];
							I.sustainLoop.end = d[3];
							I.releaseLoop.type = d[4];
							I.releaseLoop.begin = d[5];
							I.releaseLoop.end = d[6];
							// TODO: get loop samples from markers (in MARK chunk)
						}
						mSamplerInfos.push_back(I);
						chunkStoreFlag = 2;
					} break;
						
					case ID("MARK"):{
						auto cnt   = readBE<uint16_t>(mFile);
						for(uint16_t i=0; i<cnt; ++i){
							auto id  = readBE<int16_t>(mFile);
							auto pos = readBE<uint32_t>(mFile);
							uint8_t strlen = mFile.get();
							// next strlen bytes are name, must pad read out to even number of bytes
						}
						chunkStoreFlag = 2;
					} break;

					default: // unhandled chunk
						chunkStoreFlag = 1;
					}

					storeChunk(chunkSize, chunkStoreFlag);
				}
			}
		}

		if(filePosOfData != -1){
			mFile.clear(); //clear error flags (eof, etc.) so can seek
			mFile.seekg(filePosOfData);
			return true;
		}

		error:
		close();
		return false;
	}


	void close(){ mFile.close(); }


	/// Read a specific number of frames of sample data

	/// Destination array should have at least numFrames * numChannels elements.
	///
	template <class T>
	int read(T * dst, unsigned numFrames){
		if(!mFile.is_open() || mReadFrame >= mFrames) return 0;
		numFrames = std::min(numFrames, mFrames - mReadFrame);
		mReadFrame += numFrames;
		const auto numSamples = numFrames * mChannels;
		mFrameData.resize(numSamples * bytesPerSample());
		mFile.read(&mFrameData[0], mFrameData.size());
		const auto * src = (const uint8_t *)&mFrameData[0];
		const auto stride = bytesPerSample();
		switch(mEncoding){
		#define CS(ENC, type)\
			case ENC:\
				if(mDataBigEndian){\
					for(unsigned i=0; i<numSamples; ++i){\
						dst[i] = convertBE<type,T>(src);\
						src += stride;\
					}\
				} else {\
					for(unsigned i=0; i<numSamples; ++i){\
						dst[i] = convertLE<type,T>(src);\
						src += stride;\
					}\
				} break;
		CS(PCM_16, int16_t) CS(PCM_24, int24_t) CS(PCM_32, int32_t)
		CS(FLOAT, float) CS(DOUBLE, double)
		#undef CS

		case PCM_U8:{
			for(unsigned i=0; i<numSamples; ++i){
				dst[i] = convert<uint8_t,T>(*src);
				src += stride;
			}
		} break;
		case PCM_S8:{
			for(unsigned i=0; i<numSamples; ++i){
				dst[i] = convert<int8_t,T>(*src);
				src += stride;
			}
		} break;
		case ULAW:{
			for(unsigned i=0; i<numSamples; ++i){
				auto ui8 = *src ^ 0xFF;
				double v = (ui8 & 0x7F)/127.;
				static const double u = 255., invu = 1./u;
				v = (ui8&0x80?-invu:invu) * (std::pow(1.+u, v) - 1.);
				dst[i] = convert<double,T>(v);
				src += stride;
			}
		} break;
		case ALAW:{
			for(unsigned i=0; i<numSamples; ++i){
				auto ui8 = *src ^ 0x55;
				double v = (ui8 & 0x7F)/127.;
				static const double A = 87.6, invA = 1./A;
				auto vL = v*(1. + std::log(A));
				v = (ui8&0x80?invA:-invA) * (vL<1. ? vL : std::exp(vL-1.));
				dst[i] = convert<double,T>(v);
				src += stride;
			}
		} break;
		default: return 0;
		}

		return numFrames;
	}

	/// Copy all contents of file into array interleaved. Returns number of frames read.
	template<class T>
    int readAll(T * dst){
		return read(dst, mFrames);
	}

	SoundFileInfo& storeDataChunk(bool v){ mStoreDataChunk=v; return *this; }

private:
	std::ifstream mFile;
	std::vector<char> mFrameData;
	unsigned mReadFrame = 0;
	bool mStoreDataChunk = false;
};


//template <class T> char * toBE(char * dst, T v);
inline char * toBE(char * dst, uint16_t v){
	((uint8_t *)dst)[0] = v>> 8;
	((uint8_t *)dst)[1] = v;
	return dst;
}
inline char * toLE(char * dst, uint16_t v){
	((uint8_t *)dst)[1] = v>> 8;
	((uint8_t *)dst)[0] = v;
	return dst;
}

inline char * toBE(char * dst, int16_t v){ return toBE(dst, pun<uint16_t>(v)); }
inline char * toLE(char * dst, int16_t v){ return toLE(dst, pun<uint16_t>(v)); }

inline char * toBE(char * dst, uint24_t v){
	((uint8_t *)dst)[0] = v.v>>16;
	((uint8_t *)dst)[1] = v.v>> 8;
	((uint8_t *)dst)[2] = v.v;
	return dst;
}
inline char * toLE(char * dst, uint24_t v){
	((uint8_t *)dst)[2] = v.v>>16;
	((uint8_t *)dst)[1] = v.v>> 8;
	((uint8_t *)dst)[0] = v.v;
	return dst;
}

inline char * toBE(char * dst, int24_t v){ return toBE(dst, pun<uint24_t>(v)); }
inline char * toLE(char * dst, int24_t v){ return toLE(dst, pun<uint24_t>(v)); }

inline char * toBE(char * dst, uint32_t v){
	((uint8_t *)dst)[0] = v>>24;
	((uint8_t *)dst)[1] = v>>16;
	((uint8_t *)dst)[2] = v>> 8;
	((uint8_t *)dst)[3] = v;
	return dst;
}
inline char * toLE(char * dst, uint32_t v){
	((uint8_t *)dst)[3] = v>>24;
	((uint8_t *)dst)[2] = v>>16;
	((uint8_t *)dst)[1] = v>> 8;
	((uint8_t *)dst)[0] = v;
	return dst;
}

inline char * toBE(char * dst, int32_t v){ return toBE(dst, pun<uint32_t>(v)); }
inline char * toLE(char * dst, int32_t v){ return toLE(dst, pun<uint32_t>(v)); }

inline char * toBE(char * dst, uint64_t v){
	((uint8_t *)dst)[0] = v>>56;
	((uint8_t *)dst)[1] = v>>48;
	((uint8_t *)dst)[2] = v>>40;
	((uint8_t *)dst)[3] = v>>32;
	((uint8_t *)dst)[4] = v>>24;
	((uint8_t *)dst)[5] = v>>16;
	((uint8_t *)dst)[6] = v>> 8;
	((uint8_t *)dst)[7] = v;
	return dst;
}
inline char * toLE(char * dst, uint64_t v){
	((uint8_t *)dst)[7] = v>>56;
	((uint8_t *)dst)[6] = v>>48;
	((uint8_t *)dst)[5] = v>>40;
	((uint8_t *)dst)[4] = v>>32;
	((uint8_t *)dst)[3] = v>>24;
	((uint8_t *)dst)[2] = v>>16;
	((uint8_t *)dst)[1] = v>> 8;
	((uint8_t *)dst)[0] = v;
	return dst;
}

inline char * toBE(char * dst, float v){ return toBE(dst, pun<uint32_t>(v)); }
inline char * toLE(char * dst, float v){ return toLE(dst, pun<uint32_t>(v)); }

inline char * toBE(char * dst, double v){ return toBE(dst, pun<uint64_t>(v)); }
inline char * toLE(char * dst, double v){ return toLE(dst, pun<uint64_t>(v)); }

template <class T>
void writeBE(std::ofstream& f, const T& v){
	char buf[sizeof(T)];
	f.write(toBE(buf, v), sizeof(T));
}
template <class T>
void writeLE(std::ofstream& f, const T& v){
	char buf[sizeof(T)];
	f.write(toLE(buf, v), sizeof(T));
}

class SoundFileWriter : public SoundFileInfo {
public:

	SoundFileWriter(){}

	SoundFileWriter(const SoundFileInfo& info)
	:	SoundFileInfo(info)
	{}


	bool save(const std::string& path){
		if(mFrameData.empty()) return false;
		/*auto ext = path.substr(path.find_last_of('.'));
		if(ext.empty()) return false;
		ext = ext.substr(1);
		for(auto& c : ext) c = std::tolower(c);
		//printf("%s\n", ext.c_str());

		if("au"==ext || "snd"==ext) mFormat = AU;*/

		mFile.open(path, std::ios::out | std::ios::binary);
		if(!mFile.is_open()) return false;

		if(AU == mFormat){
			mFile.write(".snd", 4);
			writeBE(mFile, uint32_t(24));
			writeBE(mFile, uint32_t(mFrameData.size()));
			uint32_t enc = 0;
			switch(mEncoding){
			case PCM_S8: enc= 2; break;
			case PCM_16: enc= 3; break;
			case PCM_24: enc= 4; break;
			case PCM_32: enc= 5; break;
			case FLOAT:  enc= 6; break;
			case DOUBLE: enc= 7; break;
			case ULAW:   enc= 1; break;
			case ALAW:   enc=27; break;
			default: goto error;
			}
			writeBE(mFile, enc);
			writeBE(mFile, uint32_t(frameRate()));
			writeBE(mFile, uint32_t(channels()));
			mFile.write(&mFrameData[0], mFrameData.size());
		} else if (WAV == mFormat){
			uint16_t enc;
			uint32_t fmtExtSize=0, factSize=0;
			switch(mEncoding){
			case PCM_U8: case PCM_16: case PCM_24: case PCM_32:
				enc=1; break;
			case FLOAT: case DOUBLE:
				factSize=4; enc=3; break;
			case ALAW:
				fmtExtSize=2; factSize=4; enc=6; break;
			case ULAW:
				fmtExtSize=2; factSize=4; enc=7; break;
			default: goto error;
			}
			uint32_t fmtSize = 16 + fmtExtSize;
			int pad = mFrameData.size()&1;
			uint32_t dataSize = mFrameData.size() + pad;
			uint32_t mainSize = 4 + 8+fmtSize + (factSize?8+factSize:0) + 8+dataSize;
			mFile.write("RIFF",4); //printf("RIFF size: %u\n", mainSize);
			writeLE(mFile, mainSize); // size of file - 8
			mFile.write("WAVE",4);
			mFile.write("fmt ",4);
			writeLE(mFile, fmtSize);
			writeLE(mFile, enc);
			writeLE(mFile, uint16_t(channels()));
			writeLE(mFile, uint32_t(frameRate()));
			uint16_t bps = bytesPerSample()*8; // bits/sample
			writeLE(mFile, uint32_t(frameRate()*channels()*bps/8)); // byte rate
			writeLE(mFile, uint16_t(channels()*bps/8)); // block align
			writeLE(mFile, bps);
			if(fmtExtSize) writeLE(mFile, uint16_t(0)); // req'd for non-PCM
			if(factSize){
				mFile.write("fact",4);
				writeLE(mFile, uint32_t(4));
				writeLE(mFile, uint32_t(channels()*frames()));
			}
			mFile.write("data",4);
			writeLE(mFile, dataSize);
			mFile.write(&mFrameData[0], mFrameData.size());
			if(pad) mFile.write(0,1);
		} else if(AIFF == mFormat){
			auto setFloat10BE = [](double val, char * bytes){
				int32_t sign = (val<0)*0x8000;
				val *= (int(val>=0))*2-1;
				int32_t expon=0;
				uint32_t hiMant=0, loMant=0;
				if(val!=0.){
					double fMant = frexp(val, &expon);
					if((expon > 16384) || !(fMant < 1)){ // infinity or NaN
						expon = sign|0x7FFF; hiMant=0; loMant=0; // infinity
					} else { // finite
						expon += 16382;
						if(expon < 0){ // denormalized
							fMant = ldexp(fMant, expon);
							expon = 0;
						}
						auto floatToUnsigned = [](double v){ return ((uint32_t)(((int32_t)(v - 2147483648.0)) + 2147483647L) + 1); };
						expon |= sign;
						fMant = ldexp(fMant, 32);          
						double fsMant = floor(fMant);
						hiMant = floatToUnsigned(fsMant);
						fMant = ldexp(fMant - fsMant, 32); 
						fsMant = floor(fMant); 
						loMant = floatToUnsigned(fsMant);
					}
				}
				toBE(bytes, int16_t(expon));
				toBE(bytes+2, hiMant);
				toBE(bytes+6, loMant);
			};
			const char * encID = nullptr;
			switch(mEncoding){
			case PCM_U8: encID = "raw "; break;
			case FLOAT: encID = "fl32"; break;
			case DOUBLE: encID = "fl64"; break;
			case ALAW: encID = "alaw"; break;
			case ULAW: encID = "ulaw"; break;
			default:;
			}
			bool isAIFC = encID!=nullptr;
			uint32_t commSize = 18 + 4*isAIFC;
			uint32_t ssndSize = 8 + mFrameData.size();
			uint32_t mainSize = 4 + (8+4)*isAIFC + 8+commSize + 8+ssndSize;
			mFile.write("FORM",4);
			writeBE(mFile, mainSize);
			mFile.write(isAIFC?"AIFC":"AIFF",4);
			if(isAIFC){
				mFile.write("FVER",4);
				writeBE(mFile, uint32_t(4));
				writeBE(mFile, uint32_t(0xA2805140)); // May 23, 1990, 2:40 p.m.
			}
			mFile.write("COMM",4);
			writeBE(mFile, commSize);
			writeBE(mFile, uint16_t(channels()));
			writeBE(mFile, uint32_t(frames()));
			writeBE(mFile, uint16_t(bytesPerSample()*8));
			{char b[10]; setFloat10BE(frameRate(),b); mFile.write(b,10);}
			if(encID) mFile.write(encID,4);
			mFile.write("SSND",4);
			writeBE(mFile, ssndSize);
			writeBE(mFile, uint32_t(0)); // offset
			writeBE(mFile, uint32_t(0)); // blockAlign
			mFile.write(&mFrameData[0], mFrameData.size());
		}

		close();
		return true;

		error:
		close();
		return false;
	}

	void close(){ mFile.close(); }

	template <class T>
	int write(const T * src, unsigned numFrames){
		const auto numSamples = numFrames * mChannels;
		switch(mEncoding){
		#define CS(ENC, type)\
			case ENC:\
			if(mDataBigEndian){\
				for(unsigned i=0; i<numSamples; ++i){\
					auto v = convert<T,type>(src[i]);\
					char buf[numBytes<type>()]; toBE(buf, v);\
					for(auto c:buf) mFrameData.push_back(c);\
				}\
			} else {\
				for(unsigned i=0; i<numSamples; ++i){\
					auto v = convert<T,type>(src[i]);\
					char buf[numBytes<type>()]; toLE(buf, v);\
					for(auto c:buf) mFrameData.push_back(c);\
				}\
			} break;
		CS(PCM_16, int16_t)
		CS(PCM_24, int24_t)
		CS(PCM_32, int32_t)
		CS(FLOAT, float)
		CS(DOUBLE, double)
		#undef CS
		case PCM_S8:
			for(unsigned i=0; i<numSamples; ++i)
				mFrameData.push_back(convert<T,int8_t>(src[i]));
			break;
		case PCM_U8:
			for(unsigned i=0; i<numSamples; ++i)
				mFrameData.push_back(convert<T,uint8_t>(src[i]));
			break;
		case ULAW:
			for(unsigned i=0; i<numSamples; ++i){
				auto v = convert<T,double>(src[i]);
				auto ui8 = (uint8_t(v<0.f))<<7;
				static const double u = 255., c = 1./std::log(1. + u);
				v = std::log(1. + u*std::abs(v))*c;
				ui8 |= convert<double,int8_t>(v);
				mFrameData.push_back(~ui8);
			} break;
		default: return 0;
		}
		return numFrames;
	}

private:
	std::ofstream mFile;
	std::vector<char> mFrameData;
	unsigned mWriteFrame = 0;
};
