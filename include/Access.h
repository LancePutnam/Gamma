/*
 *  Access.h
 *  waveInteractions
 *
 *  Created by Lance on 10/2/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

namespace gam{


inline void neighborsWrap(int x1, int N, int& x0, int& x2){
	x0 = x1-1; if(x1< 0) x1+=N;
	x2 = x1+1; if(x2>=N) x2-=N;
}


class AccessStream1{
public:

	AccessStream1(int size, int begin=0)
	:	i0(0), i1(0), i2(0), mSizeM1(size-1)
	{
		i1=checkM1(begin-1);
		i2=begin;
	}

	void operator()(){
		i0=i1; i1=i2; i2 = checkP1(i2+1);
	}
	
	void operator()(int& a0, int& a1, int& a2){
		(*this)(); a0=i0; a1=i1; a2=i2;
	}

	int i0, i1, i2;

private:
	int mSizeM1;

	// wrap
	int checkM1(int v){ return v<0 ? mSizeM1 : v; }
	int checkP1(int v){ return v>mSizeM1 ? 0: v; }

	// clip
//	int checkM1(int v){ return v>0 ? v : 0; }
//	int checkP1(int v){ return v<mSizeM1 ? v : mSizeM1; }

	// fold
//	int checkM1(int v){ return v>=0 ? v : 1; }
//	int checkP1(int v){ return v<=mSizeM1 ? v : mSizeM1-1; }
};


template <class SBounds>
struct NeighborsCross2D{

	void operator()(int ix, int iy, int sizeX, int sizeY, const SBounds& sBounds){
		l = sBounds(ix-1,iy, sizeX, sizeY);
		r = sBounds(ix+1,iy, sizeX, sizeY);
		t = sBounds(ix,iy-1, sizeX, sizeY);
		b = sBounds(ix,iy+1, sizeX, sizeY);
	}

	union{
		struct{ int l, r, t, b; }; 
		int indices[4];
	};
};


template <class SBounds>
struct NeighborsDiag2D{

	void operator()(int ix, int iy, int sizeX, int sizeY, const SBounds& sBounds){
		tl = sBounds(ix-1,iy-1, sizeX, sizeY);
		br = sBounds(ix+1,iy+1, sizeX, sizeY);
		tr = sBounds(ix+1,iy-1, sizeX, sizeY);
		bl = sBounds(ix-1,iy+1, sizeX, sizeY);
	}

	union{
		struct{ int tl, br, tr, bl; }; 
		int indices[4];
	};
};



struct BoundaryClip{
	int operator()(int ix, int sizeX) const{
		return ix < 0 ? 0 : ix >= sizeX ? sizeX-1 : ix;
	}

	int operator()(int ix, int iy, int sizeX, int sizeY) const{		
		return (*this)(ix, sizeX) + (*this)(iy, sizeY)*sizeX;
	}
};


struct BoundaryWrap{
	int operator()(int ix, int sizeX) const{
		return ix < 0 ? sizeX+ix : ix >= sizeX ? ix-sizeX : ix;
	}

	int operator()(int ix, int iy, int sizeX, int sizeY) const{
		return (*this)(ix, sizeX) + (*this)(iy, sizeY)*sizeX;
	}
};




template <class Sneigh, class Sbounds=BoundaryWrap>
struct Neighbors2D : public Sneigh{

	Neighbors2D(int sizeX, int sizeY): sizeX(sizeX), sizeY(sizeY){}

	void operator()(int ix, int iy){
		(*this)(ix, iy, sizeX, sizeY, mBounds);
	}

private:
	Sbounds bounds;
	int sizeX, sizeY;
};




template <int N=1, class Sbounds=BoundaryWrap>
struct Neighbors1D {

	Neighbors1D(int sizeX): sizeX(sizeX){}

	void operator()(int ix){
		for(int i=0; i<N; ++i){
			l[i] = bounds(ix-(i+1), sizeX);
			r[i] = bounds(ix+(i+1), sizeX);
		}
	}
	
	template <class T>
	T access(const T * src, int ix){ return src[bounds(ix, sizeX)]; }
	
	int l[N], r[N];

private:
	Sbounds bounds;
	int sizeX;
};




/*

int sx = 64;
int sy = 32;
float data[sx*sy];

Neighbors2D<NeighborsCross2D, BoundaryWrap> n2d(sx, sy);

for(int j=0; j<sy; ++j){
	for(int i=0; i<sx; ++i){
		n2d(i, j);
		float d1xc = (data[n2d.r] - data[n2d.l])*0.5;
		float d1yc = (data[n2d.b] - data[n2d.t])*0.5;
	}
}

*/


} // gam::
