#ifndef GAMMA_DIMS_H_INC
#define GAMMA_DIMS_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include "Gamma/Access.h"
#include "Gamma/Conversion.h"
#include "Gamma/scl.h"

namespace gam{


/// Stores dimensions of a lattice (up to 4 dimensions) and provides conversion \
functions between real positions and indices.

/// The positions must lie within [-1,1) relative to the smallest dimension
/// other than those with a size of 1. The position interval of larger dimensions
/// is scaled up linearly in proportion to [-1,1).
/// Conversions to and from positions and indices do not perform any bounds
/// checking. Potential out-of-bounds positions and indices should be validated 
/// with the various contains() methods, otherwise undefined behavior will
/// result.
class Dims{
public:

	/// @param[in] s1	size of dimension 1
	/// @param[in] s2	size of dimension 2
	/// @param[in] s3	size of dimension 3
	/// @param[in] s4	size of dimension 4
	Dims(int s1, int s2=1, int s3=1, int s4=1){ resize(s1,s2,s3,s4); }

	/// Returns whether the point (p1, p2, -1, -1) is within bounds.
	bool contains(float p1, float p2) const;

	/// Returns whether the point (p1, p2, p3, -1) is within bounds.
	bool contains(float p1, float p2, float p3) const;

	/// Returns whether the point (p1, p2, p3, p4) is within bounds.
	bool contains(float p1, float p2, float p3, float p4) const;
	
	/// Returns whether flat index is within space.
	bool contains(uint64_t index) const;

	void indExpand(int i, int& i1, int& i2, int& i3) const;

	/// Convert dimension indices to flat lattice index.
	
	/// The indices lie in [0, N) where N is the size of the dimension.
	///
	int indFlatten(int i1, int i2) const;
	int indFlatten(int i1, int i2, int i3) const;

	/// Convert flat index to dimension positions.
	void indToPos(int i, float& p1, float& p2, float& p3){
		int i1,i2,i3; indExpand(i, i1,i2,i3);
		indToPos(i1,i2,i3, p1,p2,p3);
	}

	/// Convert dimension indices to Cartesian space coordinates.
	void indToPos(int i1, int i2, int i3, float& p1, float& p2, float& p3) const;
	
	float indToPos1(int i) const;	///< Converts index to position of dimension 1
	float indToPos2(int i) const;	///< Converts index to position of dimension 2
	float indToPos3(int i) const;	///< Converts index to position of dimension 3
	float indToPos4(int i) const;	///< Converts index to position of dimension 4

	/// Get neighboring indices. No bounds checking.
	void indNeighbors(int i1, int i2, int i3, int (&inds)[6], bool (&valid)[6]) const;
	void indNeighbors1(int i, int * inds, bool * valid) const;
	void indNeighbors2(int i, int * inds, bool * valid) const;
	void indNeighbors3(int i, int * inds, bool * valid) const;
	void indNeighbors4(int i, int * inds, bool * valid) const;

	/// Converts position (p1, p2, -1, -1) to flat index.
	int posToInd(float p1, float p2) const;
	
	/// Converts position (p1, p2, p3, -1) to flat index.
	int posToInd(float p1, float p2, float p3) const;

	/// Converts position (p1, p2, p3, p4) to flat index.
	int posToInd(float p1, float p2, float p3, float p4) const;

	int posToInd1(float p) const;	///< Converts position to index of dimension 1
	int posToInd2(float p) const;	///< Converts position to index of dimension 2
	int posToInd3(float p) const;	///< Converts position to index of dimension 3
	int posToInd4(float p) const;	///< Converts position to index of dimension 4

	void posToIndFrac(float p1, float p2, float p3, float& f1, float& f2, float& f3) const {
		f1 = posToIndFrac(p1, 0);
		f2 = posToIndFrac(p2, 1);
		f3 = posToIndFrac(p3, 2);
	}
	
	float max1() const { return mMax[0]; }	///< Returns maximum position of dimension 1
	float max2() const { return mMax[1]; }	///< Returns maximum position of dimension 2
	float max3() const { return mMax[2]; }	///< Returns maximum position of dimension 3
	float max4() const { return mMax[3]; }	///< Returns maximum position of dimension 4

	uint64_t size() const;					///< Returns total size of lattice
	int size1() const { return n[0]; }		///< Returns size of dimension 1
	int size2() const { return n[1]; }		///< Returns size of dimension 2
	int size3() const { return n[2]; }		///< Returns size of dimension 3
	int size4() const { return n[3]; }		///< Returns size of dimension 4
	uint64_t size12() const;				///< Returns size of plane 12
	uint64_t size123() const;				///< Returns size of cube 123

	/// Returns smallest dimension size (excluding size 1)
	int sizeMin() const { return mSizeMin; }
	
protected:

	int n[4];				// size of dimensions
	int mSizeMin;			// size of minimum >1 sized dimension
	float m1_m;				// 1/min
	float mMax[4];			// maximum positions
	float mMul[4], mAdd[4];	// values for converting to indices

	void resize(int s1, int s2, int s3, int s4){
	
		n[0]=s1; n[1]=s2; n[2]=s3; n[3]=s4;

		// find smallest >1 sized dimension
		// find max
		int sizeMax=0;
		for(int i=0; i<4; ++i){ if(n[i] > sizeMax) sizeMax=n[i]; }
		
		// find min >1
		mSizeMin=sizeMax;
		for(int i=0; i<4; ++i){ if(n[i] < sizeMax && n[i]>1) mSizeMin=n[i]; }
		
		m1_m = 1./mSizeMin;
		
		for(int i=0; i<4; ++i){
			mMax[i] = float(n[i])*m1_m;
			mAdd[i] = n[i]*0.5;
			mMul[i] = mAdd[i]/mMax[i];
		}
	}

	// Does position lie in dimension interval?
	bool hasPos(float p, int d) const { return scl::withinIE(p, -mMax[d], mMax[d]); }

	void indNeighbors(int i, int d, int * inds, bool * valid) const {
		inds[0] = i-1; inds[1] = i+1;
		valid[0] = inds[0] >= 0;
		valid[1] = inds[1] <  n[d];
	}

	// Convert a dimension index to position
	float indToPos(int i, int d) const { return ((i<<1) - n[d]) * m1_m; }
	
	// Direct conversion of position to dimension index
	//int posToInd(float p, int d) const { return scl::wrap((p/mMax[d])*0.5f + 0.5f, 0.9999f) * n[d]; }
	//int posToInd(float p, int d) const { return (p/mMax[d]*0.5f + 0.5f) * n[d]; }
	int posToInd(float p, int d) const { return castIntTrunc(posToIndFrac(p,d)); }
	float posToIndFrac(float p, int d) const { return p*mMul[d] + mAdd[d]; }
};





inline bool Dims::contains(float p1, float p2) const {
	return hasPos(p1,0) && hasPos(p2,1);
}

inline bool Dims::contains(float p1, float p2, float p3) const {
	return hasPos(p1,0) && hasPos(p2,1) && hasPos(p3,2);
}

inline bool Dims::contains(float p1, float p2, float p3, float p4) const {
	return hasPos(p1,0) && hasPos(p2,1) && hasPos(p3,2) && hasPos(p4,3);
}

inline bool Dims::contains(uint64_t index) const { return scl::withinIE(index, uint64_t(0), size()); }

inline void Dims::indExpand(int i, int& i1, int& i2, int& i3) const {
	i1 = i % size1();
	i2 = (i / size1()) % size2();
	i3 = i / size12();
 }

inline int Dims::indFlatten(int i1, int i2) const { return i1 + i2*size1(); }

inline int Dims::indFlatten(int i1, int i2, int i3) const {
	return index3to1(i1,i2,i3, size1(),size2());
}

inline void Dims::indNeighbors(int i1, int i2, int i3, int (&inds)[6], bool (&valid)[6]) const {
//	inds[0] = i1-1; inds[1] = i1+1;
//	inds[2] = i2-1;	inds[3] = i2+1;
//	inds[4] = i3-1;	inds[5] = i3+1;
//	
//	valid[0] = inds[0] >= 0;
//	valid[1] = inds[1] <  size1();
//	valid[2] = inds[2] >= 0;
//	valid[3] = inds[3] <  size2();
//	valid[4] = inds[4] >= 0;
//	valid[5] = inds[5] <  size3();
	
	indNeighbors1(i1, inds  , valid  );
	indNeighbors2(i2, inds+2, valid+2);
	indNeighbors3(i3, inds+4, valid+4);
}

inline void Dims::indNeighbors1(int i, int * inds, bool * v) const { indNeighbors(i,0,inds,v); }
inline void Dims::indNeighbors2(int i, int * inds, bool * v) const { indNeighbors(i,1,inds,v); }
inline void Dims::indNeighbors3(int i, int * inds, bool * v) const { indNeighbors(i,2,inds,v); }
inline void Dims::indNeighbors4(int i, int * inds, bool * v) const { indNeighbors(i,3,inds,v); }

inline void Dims::indToPos(int i1, int i2, int i3, float& p1, float& p2, float& p3) const {
	p1=indToPos1(i1); p2=indToPos2(i2); p3=indToPos3(i3);
}

inline float Dims::indToPos1(int i) const { return indToPos(i,0); }
inline float Dims::indToPos2(int i) const { return indToPos(i,1); }
inline float Dims::indToPos3(int i) const { return indToPos(i,2); }
inline float Dims::indToPos4(int i) const { return indToPos(i,3); }

inline int Dims::posToInd1(float v) const { return posToInd(v,0); }
inline int Dims::posToInd2(float v) const { return posToInd(v,1); }
inline int Dims::posToInd3(float v) const { return posToInd(v,2); }
inline int Dims::posToInd4(float v) const { return posToInd(v,3); }


inline int Dims::posToInd(float p1, float p2) const {
	return posToInd1(p1) + size1() * posToInd2(p2);
}

inline int Dims::posToInd(float p1, float p2, float p3) const {
	return posToInd1(p1) + size1() * (posToInd2(p2) + size2() * posToInd3(p3));
}

inline int Dims::posToInd(float p1, float p2, float p3, float p4) const {
	// i1 + s1*i2 + s1*s2*i3 + s1*s2*s3*i4
	return posToInd1(p1) + size1() * (posToInd2(p2) + size2() * (posToInd3(p3) + size3() * posToInd4(p4)));
}

//inline int Dims::posToIndFlat(float p1, float p2, float p3) const {
//
//	// map and wrap all positions into [0,1) interval
//	p1 = posToInd(p1, 0);
//	p2 = posToInd(p2, 1);
//	p3 = posToInd(p3, 2);
//
//	return int(p1*size1()) + size1()*(int(p2*size2()) + size2()*int(p3*size3()));
//}

inline uint64_t Dims::size() const { return size123()*uint64_t(size4()); }
inline uint64_t Dims::size12() const { return uint64_t(size1())*uint64_t(size2()); }
inline uint64_t Dims::size123() const { return size12()*uint64_t(size3()); }

} // gam::

#endif
