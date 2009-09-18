#ifndef INC_SERIALIZE_H
#define INC_SERIALIZE_H

/*
General serialization and deserialization for file and network i/o.

Either single elements or arrays of elements can be read or written to a
buffer. Each group of elements is preceeded by a header that tells the type of
elements and the size of the array.

A string is a C-style null-terminated array of characters.
A boolean is a single byte where a value of 0 is false and non-zero value is true.

Byte ordering is little endian. If the target architecture is big endian,
then pass in the preprocessor flag -DSER_IS_BIG_ENDIAN.

*/

#include <stddef.h>
#include <string.h>
#include <vector>
#include <string>



//#define SER_IS_BIG_ENDIAN

/* Serialization data types */
enum{
	SER_FLOAT32	= 'f',	/* 32-bit IEEE float */
	SER_FLOAT64	= 'd',	/* 64-bit IEEE float */
	SER_INT8	= 'h',	/*  8-bit 2's-complement signed integer */
	SER_INT16	= 'H',	/* 16-bit 2's-complement signed integer */
	SER_INT32	= 'i',	/* 32-bit 2's-complement signed integer */
	SER_INT64	= 'I',	/* 64-bit 2's-complement signed integer */
	SER_UINT8	= 't',	/*  8-bit unsigned integer */
	SER_UINT16	= 'T',	/* 16-bit unsigned integer */
	SER_UINT32	= 'u',	/* 32-bit unsigned integer */
	SER_UINT64	= 'U',	/* 64-bit unsigned integer */
	SER_SUB		= '_'	/* substructure */
};


/* Serialized data header */
struct SerHeader{
	uint8_t type;		/* type of data */
	uint32_t num;		/* number of data elements */
};


/* Serialized data header. Don't think this is needed... */
//struct SerStruct{
//	SerHeader header;
//	void * data;
//};


/* Copy 1-byte elements */
uint32_t serCopy1(void * dst, const void * src, uint32_t num);

/* Copy 2-byte elements in little endian byte order */
uint32_t serCopy2(void * dst, const void * src, uint32_t num);

/* Copy 4-byte elements in little endian byte order */
uint32_t serCopy4(void * dst, const void * src, uint32_t num);

/* Copy 8-byte elements in little endian byte order */
uint32_t serCopy8(void * dst, const void * src, uint32_t num);

/* Decode serialized data. Returns number of bytes parsed */
uint32_t serDecode(const char * b, void * data);

/*  */
SerHeader serGetHeader(const char * buf);

/* Returns header size in bytes */
inline int serHeaderSize(){ return 5; }

/* Returns size of data type in bytes */
int serTypeSize(uint8_t t);

uint32_t serElementsSize(const SerHeader * h);

void serSwap(char * a, char * b);
void serSwapBytes2(void * v);
void serSwapBytes4(void * v);
void serSwapBytes8(void * v);

/* Returns human-readable string of header */
const char * serStringifyHeader(const SerHeader * h);

/* Returns human-readable string of type */
const char * serStringifyType(uint8_t t);



	





#define SOH serHeaderSize()

inline void serHeaderWrite(char * b, uint8_t type, uint32_t num){
	b[0] = type;
	serCopy4(b+1, &num, 1);
}

#define DO(B, N)\
	serHeaderWrite(b, type, N);\
	return SOH + serCopy##B(b+SOH, v, N)
inline uint32_t serEncode1(uint8_t type, char * b, const void * v, uint32_t n){ DO(1,n); }
inline uint32_t serEncode2(uint8_t type, char * b, const void * v, uint32_t n){ DO(2,n); }
inline uint32_t serEncode4(uint8_t type, char * b, const void * v, uint32_t n){ DO(4,n); }
inline uint32_t serEncode8(uint8_t type, char * b, const void * v, uint32_t n){ DO(8,n); }
#undef DO

inline uint32_t serEncodeFloat32(char * b, const float * v   , uint32_t n){ return serEncode4(SER_FLOAT32, b,v,n); }
inline uint32_t serEncodeFloat64(char * b, const double * v  , uint32_t n){ return serEncode8(SER_FLOAT64, b,v,n); }
inline uint32_t serEncodeInt8   (char * b, const int8_t * v  , uint32_t n){ return serEncode1(SER_INT8   , b,v,n); }
inline uint32_t serEncodeInt16  (char * b, const int16_t * v , uint32_t n){ return serEncode2(SER_INT16  , b,v,n); }
inline uint32_t serEncodeInt32  (char * b, const int32_t * v , uint32_t n){ return serEncode4(SER_INT32  , b,v,n); }
inline uint32_t serEncodeInt64  (char * b, const int64_t * v , uint32_t n){ return serEncode8(SER_INT64  , b,v,n); }
inline uint32_t serEncodeBool   (char * b, const bool * v    , uint32_t n){ return serEncode1(SER_UINT8  , b,v,n); }
inline uint32_t serEncodeUInt8  (char * b, const uint8_t * v , uint32_t n){ return serEncode1(SER_UINT8  , b,v,n); }
inline uint32_t serEncodeUInt16 (char * b, const uint16_t * v, uint32_t n){ return serEncode2(SER_UINT16 , b,v,n); }
inline uint32_t serEncodeUInt32 (char * b, const uint32_t * v, uint32_t n){ return serEncode4(SER_UINT32 , b,v,n); }
inline uint32_t serEncodeUInt64 (char * b, const uint64_t * v, uint32_t n){ return serEncode8(SER_UINT64 , b,v,n); }






namespace ser{

template<class T> inline uint32_t encode(char * b, const T * v, uint32_t n){ return 0; }

#define DEF(t,T,c)\
template<> inline uint32_t encode(char * b, const t * v, uint32_t n){ return serEncode##T(b, c v,n); }
DEF(float,   Float32,)
DEF(double,  Float64,)
DEF(char,    Int8, (const int8_t *))
DEF(int8_t,  Int8,)
DEF(int16_t, Int16,)
DEF(int32_t, Int32,)
DEF(int64_t, Int64,)
DEF(uint8_t, UInt8,)
DEF(bool,    Bool,)
DEF(uint16_t,UInt16,)
DEF(uint32_t,UInt32,)
DEF(uint64_t,UInt64,)
#undef DEF

template<class T> uint8_t getType();
template<> inline uint8_t getType<float   >(){ return 'f'; }
template<> inline uint8_t getType<double  >(){ return 'd'; }
template<> inline uint8_t getType<bool    >(){ return 't'; }
template<> inline uint8_t getType<uint8_t >(){ return 't'; }
template<> inline uint8_t getType<uint16_t>(){ return 'T'; }
template<> inline uint8_t getType<uint32_t>(){ return 'u'; }
template<> inline uint8_t getType<uint64_t>(){ return 'U'; }
template<> inline uint8_t getType<int8_t  >(){ return 'h'; }
template<> inline uint8_t getType<int16_t >(){ return 'H'; }
template<> inline uint8_t getType<int32_t >(){ return 'i'; }
template<> inline uint8_t getType<int64_t >(){ return 'I'; }



///
struct Serializer{

	template <class T>
	Serializer& operator<< (T v);

	Serializer& operator<< (const char * v);
	Serializer& operator<< (const std::string& v);
	
	template <class T>
	Serializer& add(const T * v, uint32_t num);

	const std::vector<char>& buf() const;

private:
	std::vector<char> mBuf, mTemp;
	
	template <class T> void checkSize(uint32_t n=1);
};


///
struct Deserializer{
	
	Deserializer(const std::vector<char>& b);
	
	Deserializer(const char * b, uint32_t n);
	
	template <class T>
	Deserializer& operator>> (T& v);
	Deserializer& operator>> (char * v);
	Deserializer& operator>> (std::string& v);

	const std::vector<char>& buf() const;

private:
	int mStart;
	std::vector<char> mBuf;
	char * bufDec();
};




struct SyncedMemory{

	SyncedMemory(void * data, uint8_t type, uint32_t numElems=1)
	:	curr((char *)data), prev(0), type(type), numElems(numElems)
	{
		prev = new char[size()];
		bzero(prev, size());
	}

	~SyncedMemory(){ delete[] prev; }
	
	bool changed(){ return 0 != memcmp(curr, prev, size()); }
	void update(){ memcpy(prev, curr, size()); }
	
	/// Returns number of elements copied
	int copyTo(const SyncedMemory& m){
		if(numElems != m.numElems) return 0;

		if(type == m.type){
			memcpy(m.curr, curr, size());
		}
		
		else{
			#define CP(d,s)\
				for(unsigned i=0; i<numElems; ++i){\
					((d *)m.curr)[i] = ((s *)curr)[i];\
				}
			
			switch(type){
			case 'f':
				switch(m.type){
				case 'd': CP( double, float)
				case 't': CP(uint8_t, float)
				}
				break;

			case 'd':	
				switch(m.type){
				case 'f': CP(  float, double)
				case 't': CP(uint8_t, double)
				}
				break;
				
			case 't':	
				switch(m.type){
				case 'f': CP(  float, uint8_t)
				case 'd': CP( double, uint8_t)
				}
				break;
			}
			#undef CP
		}
		
		return numElems;
	}

	char * curr;	// pointer to memory
	char * prev;	// local copy
	uint8_t type;
	uint32_t numElems;	// number of elements
	uint32_t size(){ return serTypeSize(type)*numElems; }
};



// Watches two memory locations for changes and synchronizes them accordingly.

// If either memory location changes between calls to synchronize(), then
// the memory that changed will be copied to the other memory.

// If the two data types differ, then casting will be performed.
struct MemorySynchronizer{



	SyncedMemory data1;
	SyncedMemory data2;
};


template <class T> Serializer& Serializer::operator<< (T v){
	return add(&v, 1);
}

template <class T> Serializer& Serializer::add(const T * v, uint32_t num){
	checkSize<T>(num);
	uint32_t n = ser::encode(&mTemp[0], v, num);
	mBuf.insert(mBuf.end(), mTemp.begin(), mTemp.begin()+n);
	return *this;		
}

template <class T> void Serializer::checkSize(uint32_t n){
	uint32_t need = n*sizeof(T) + serHeaderSize();
	if(mTemp.size() < need) mTemp.resize(need);
}


template <class T> Deserializer& Deserializer::operator>> (T& v){
	uint32_t n = serDecode(bufDec(), &v);
	mStart += n;
	return *this;
}

} // ser::


inline void serSwap(char * a, char * b){ char t=*a; *a=*b; *b=t; }

inline void serSwapBytes2(void * v){	
	char * b = (char *)v;
	serSwap(b  , b+1);
}

inline void serSwapBytes4(void * v){	
	char * b = (char *)v;
	serSwap(b  , b+3);
	serSwap(b+1, b+2);
}

inline void serSwapBytes8(void * v){	
	char * b = (char *)v;
	serSwap(b  , b+7);
	serSwap(b+1, b+6);
	serSwap(b+2, b+5);
	serSwap(b+3, b+4);
}


inline uint32_t serCopy1(void * d, const void * s, uint32_t n){
	memcpy(d,s,n); return n;
}

#define DEF_LE(B, S)\
inline uint32_t serCopy##B(void * d, const void * s, uint32_t n){\
	n = n<<S;\
	memcpy(d,s,n);\
	return n;\
}

#define DEF_BE(p, B, S)\
inline uint32_t serCopy##B(void * d, const void * s, uint32_t n){\
	n = n<<S;\
	memcpy(d,s,n);\
	char * t = (char *)d;\
	for(uint32_t i=0; i<n; i+=B) swapBytes##B(t+i);\
	return n;\
}

#ifdef SER_IS_BIG_ENDIAN
DEF_BE(2,1) DEF_BE(4,2) DEF_BE(8,3)
#else
DEF_LE(2,1) DEF_LE(4,2) DEF_LE(8,3)
#endif

#undef DEF_BE
#undef DEF_LE

#undef SOH







#endif
