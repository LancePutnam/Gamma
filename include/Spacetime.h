#ifndef GAMMA_SPACETIME_H_INC
#define GAMMA_SPACETIME_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include "scl.h"

namespace gam{


/// Dimensions of a spacetime grid.

/// Conversions to and from positions and indices do not perform any bounds
/// checking. Potential out-of-bounds coordinates and indices should be validated 
/// with the various contains() methods, otherwise undefined behavior will
/// result.
class Dims{
public:

	Dims(int s1, int s2, int s3=1, int s4=1){
		resize(s1,s2,s3,s4);
	}

	/// Returns whether the point (p1, p2, -1, -1) is within bounds.
	bool contains(float p1, float p2) const;

	/// Returns whether the point (p1, p2, p3, -1) is within bounds.
	bool contains(float p1, float p2, float p3) const;

	/// Returns whether the point (p1, p2, p3, p4) is within bounds.
	bool contains(float p1, float p2, float p3, float p4) const;
	
	/// Returns whether flat index is within space.
	bool contains(int index) const;

	void indExpand(int i, int& i1, int& i2, int& i3) const;

	// Convert dimension indices to flat lattice index.
	// The indices lie in [0, N) where N is the size of the dimension.
	int indFlatten(int i1, int i2, int i3) const;

	// Convert flat index to dimension positions.
	void indToPos(int i, float& p1, float& p2, float& p3){
		int i1,i2,i3; indExpand(i, i1,i2,i3);
		indToPos(i1,i2,i3, p1,p2,p3);
	}

	// Convert dimension indices to Cartesian space coordinates.
	void indToPos(int i1, int i2, int i3, float& p1, float& p2, float& p3) const;

	// Get neighboring indices. No bounds checking.
	void indNeighbors(int i1, int i2, int i3, int (&inds)[6], bool (&valid)[6]) const;
	void indNeighbors1(int i, int * inds, bool * valid) const;
	void indNeighbors2(int i, int * inds, bool * valid) const;
	void indNeighbors3(int i, int * inds, bool * valid) const;
	void indNeighbors4(int i, int * inds, bool * valid) const;

	int posToInd1(float p) const;
	int posToInd2(float p) const;
	int posToInd3(float p) const;
	int posToInd4(float p) const;

	// Converts coordinate (p1, p2, -1, -1) to flat index.
	int posToInd(float p1, float p2) const;
	
	// Converts coordinate (p1, p2, p3, -1) to flat index.
	int posToInd(float p1, float p2, float p3) const;

	// Converts coordinate (p1, p2, p3, p4) to flat index.
	int posToInd(float p1, float p2, float p3, float p4) const;
	
	void posToIndFrac(float p1, float p2, float p3, float& f1, float& f2, float& f3) const {
		f1 = posToIndFrac(p1, 0);
		f2 = posToIndFrac(p2, 1);
		f3 = posToIndFrac(p3, 2);
	}
	
	float max1() const { return mMax[0]; }
	float max2() const { return mMax[1]; }
	float max3() const { return mMax[2]; }
	float max4() const { return mMax[3]; }

	int size() const { return size123()*size4(); }
	int size1() const { return n[0]; }
	int size2() const { return n[1]; }
	int size3() const { return n[2]; }
	int size4() const { return n[3]; }
	int size12() const { return size1()*size2(); }
	int size123() const { return size12()*size3(); }
	
	// Returns smallest dimension size (excluding size 1)
	int sizeMin() const { return mSizeMin; }
	
protected:

	int n[4];		// size of dimensions
	int mSizeMin;		// size of minimum >1 sized dimension
	float m1_m;		// 1/min
	float mMax[4];	// maximum positions
	float mMul[4];
	float mAdd[4];

	void resize(int s1, int s2, int s3, int s4){
	
		n[0]=s1; n[1]=s2; n[2]=s3; n[3]=s4;
	
		// find smallest >1 sized dimension
		int sz[4];
		for(int i=0; i<4; ++i) sz[i] = (n[i]!=1 ? n[i] : size());
		mSizeMin = scl::min(scl::min(sz[0], sz[1]), scl::min(sz[2], sz[3]));
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
	int posToInd(float p, int d) const { return scl::castIntTrunc(posToIndFrac(p,d)); }
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

inline bool Dims::contains(int index) const { return scl::withinIE(index, 0, size()); }

inline void Dims::indExpand(int i, int& i1, int& i2, int& i3) const {
	i1 = i % size1();
	i2 = (i / size1()) % size2();
	i3 = i / size12();
 }

inline int Dims::indFlatten(int i1, int i2, int i3) const {
	return scl::index3to1(i1,i2,i3, n[0],n[1]);
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

inline void Dims::indNeighbors1(int i, int * inds, bool * valid) const {
	indNeighbors(i,0,inds,valid);
}

inline void Dims::indNeighbors2(int i, int * inds, bool * valid) const {
	indNeighbors(i,1,inds,valid);
}

inline void Dims::indNeighbors3(int i, int * inds, bool * valid) const {
	indNeighbors(i,2,inds,valid);
}

inline void Dims::indNeighbors4(int i, int * inds, bool * valid) const {
	indNeighbors(i,3,inds,valid);
}

inline void Dims::indToPos(int i1, int i2, int i3, float& p1, float& p2, float& p3) const {
	p1 = indToPos(i1, 0);
	p2 = indToPos(i2, 1);
	p3 = indToPos(i3, 2);
}

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





template <class T>
class Spacetime : public Dims{

public:

	//Dims dims;

	/// 
	Spacetime(int size1, int size2, int size3=1, int size4=1);
	
	virtual ~Spacetime();
	
	/// Resize dimensions of space
	void resize(int size1, int size2, int size3=1, int size4=1);

	/// Access elements
	T  operator[](int i) const { return mElems[i]; }
	T& operator[](int i)       { return mElems[i]; }

	bool updated() const { return mUpdated; }

	static Spacetime<T> singularity;	// the default do-nothing space

protected:
	T * mElems;			// shortest memory distance is along x
	bool mUpdated;		// whether space values were just updated (for avatars)
};




// Implementation --------------------------------------------------------------

#define TEM template <class T>


TEM Spacetime<T> Spacetime<T>::singularity(0,0,0,0);

TEM Spacetime<T>::Spacetime(int s1, int s2, int s3, int s4)
:	Dims(0,0,0,0), mElems(0), mUpdated(false)
{
	resize(s1, s2, s3, s4);
}

TEM Spacetime<T>::~Spacetime(){
	delete [] mElems;
}

//TEM void Spacetime<T>::spaceRange(float m1, float m2, float m3, float m4){
//	mMax[0]=m1; mMax[1]=m2; mMax[2]=m3; mMax[3]=m4;
//
//	for(int i=0; i<4; ++i){
//		mAdd[i] = mSize[i]*0.5;
//		mMul[i] = mAdd[i]/mMax[i];
//	}
//}

TEM void Spacetime<T>::resize(int s1, int s2, int s3, int s4){

	int sizeOld = size();
	int sizeNew = s1*s2*s3*s4;

	if(sizeNew != sizeOld){
		delete [] mElems;
		mElems = new T[sizeNew + 1]; // add one to size for out-of-bounds positions
		
		Dims::resize(s1,s2,s3,s4);
	}
}

#undef TEM

} // gam::

#endif
