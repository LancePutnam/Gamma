#ifndef GAMMA_PARTICLE_H_INC
#define GAMMA_PARTICLE_H_INC

/*	Gamma - Generic processing library
	See COPYRIGHT file for authors and license information */

#include "Gamma/arr.h"

#define TEM template<class T>

namespace gam{

// Util class for doing operations on particle properties
namespace ParticleOp{

	TEM void boundNone(T * pos, T * dpos, long n, T lo, T hi){}
	
	TEM void boundWrap(T * pos, T * dpos, long n, T lo, T hi){
		for(long i=0; i<n; ++i) scl::wrap(pos[i], hi, lo);
	}
	
	TEM void boundReflect(T * pos, T * dpos, long n, T lo, T hi){
		T diff = hi - lo;

		for(long i=0; i<n; i++){
			T val = *pos;
			T dval = *dpos;
			
			if(val < lo){
				val = (T)2 * lo - val;
				dval = -dval; 
			}
			else if(val >= hi){
				val = (T)2 * hi - val;
				dval = -dval;
			}
			
			*pos++ = val;
			*dpos++ = dval;
		}	
	}
	
	TEM void boundClip(T * pos, T * dpos, long n, T lo, T hi){
		for(long i=0; i<n; ++i) scl::clip(pos[i], hi, lo);
	}
	
} // ParticleOp::



/// Group of space-time objects.
template <class T=gam::real>
class PointParticles {
public:

	/// @param[in] size			Number of particles
	/// @param[in] dimensions	Number of dimensions
	/// @param[in] rates		Number of rate derivatives (2 = pos & vel)
	PointParticles(uint32_t size, uint32_t dimensions=3, uint32_t rates=2);
	~PointParticles();

	void operator()();				///< Update space state of all particles.
	//void prev(uint32_t index);	///< Set indexed particle to previous state.

	T * at(uint32_t dim, uint32_t rate=0);	///< Returns pointer to particles with dimension and rate
	T * pos(uint32_t dim);					///< Returns pointer to positions in dimension
	T * vel(uint32_t dim);					///< Returns pointer to velocities in dimension
	T * acc(uint32_t dim);					///< Returns pointer to accelerations in dimension
	T * jrk(uint32_t dim);					///< Returns pointer to jerks in dimension
	
	uint32_t dims();	///< Returns number of dimensions.
	uint32_t rates();	///< Returns number of rates.
	uint32_t size();	///< Returns number of particles.

	// Setters for single particles.
	void pos(uint32_t dim, uint32_t index, T value);
	void vel(uint32_t dim, uint32_t index, T value);
	void acc(uint32_t dim, uint32_t index, T value);
	void jrk(uint32_t dim, uint32_t index, T value);

	void zero();	///< Zeros all space-time derivatives.

	void print();

protected:
	T * dt;				// n-dimensional lattice of space-time derivatives.
	uint32_t mSize;		// # particles
	uint32_t mDims;		// # spatial dimensions
	uint32_t mRates;	// # space-time derivatives
};



// Implementation_______________________________________________________________

// PointParticles

TEM PointParticles<T>::PointParticles(uint32_t size, uint32_t dims, uint32_t rates):
	dt(0), mSize(size), mDims(dims), mRates(rates)
{
	mem::resize(dt, 0, size * dims * rates);
	zero();
}


TEM PointParticles<T>::~PointParticles(){ mem::free(dt); }

TEM void PointParticles<T>::operator()(){

	const uint32_t dtSize = size() * dims();

	for(uint32_t r = dtSize * (rates() - 1); r > 0; r -= dtSize){
	
		T * d = dt + r;
		T * dm1 = d - dtSize;
	
		for(uint32_t p = 0; p < dtSize; p++){
			*dm1++ += *d++;
		}
	}
}

TEM void PointParticles<T>::zero(){ mem::deepZero(dt, size() * dims() * rates()); }

TEM void PointParticles<T>::print(){
	T * val = dt;
	for(uint32_t r=0; r<rates(); r++){
		for(uint32_t d=0; d<dims(); d++){
			printf(" [");
			for(uint32_t p=0; p<size(); p++){
				printf("% 4.2f ", (float)*val++);
			}
			printf("]");
		}
		printf("\n");
	}
}

TEM inline T * PointParticles<T>::at(uint32_t dim, uint32_t rate){
	return dt + size() * (dims() * rate + dim);
}
TEM inline T * PointParticles<T>::pos(uint32_t dim){ return dt + size() * dim; }
TEM inline T * PointParticles<T>::vel(uint32_t dim){ return dt + size() * (dims() + dim); }
TEM inline T * PointParticles<T>::acc(uint32_t dim){ return dt + size() * ((dims()<<1) + dim); }
TEM inline T * PointParticles<T>::jrk(uint32_t dim){ return dt + size() * (dims()*3 + dim); }
TEM inline void PointParticles<T>::pos(uint32_t dim, uint32_t i, T v){ pos(dim)[i] = v; }
TEM inline void PointParticles<T>::vel(uint32_t dim, uint32_t i, T v){ vel(dim)[i] = v; }
TEM inline void PointParticles<T>::acc(uint32_t dim, uint32_t i, T v){ acc(dim)[i] = v; }
TEM inline void PointParticles<T>::jrk(uint32_t dim, uint32_t i, T v){ jrk(dim)[i] = v; }
TEM inline uint32_t PointParticles<T>::dims(){ return mDims; }
TEM inline uint32_t PointParticles<T>::rates(){ return mRates; }
TEM inline uint32_t PointParticles<T>::size(){ return mSize; }

} // gam::

#undef TEM

#endif

