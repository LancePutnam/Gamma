/* 
 * This file is based largely on the following software distribution:
 * 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * 
 *                              FFTPACK
 * 
 * Reference                                                                                                                        
 *    P.N. Swarztrauber, Vectorizing the FFTs, in Parallel Computations
 *    (G. Rodrigue, ed.), Academic Press, 1982, pp. 51--83.                                                                                                                   
 * 
 *     http://www.netlib.org/fftpack/
 * 
 * Updated to single, double, and extended precision,
 * and translated to ISO-Standard C/C++ (without aliasing)
 * on 10 October 2005 by Andrew Fernandes <andrew_AT_fernandes.org>
 * 
 *                   Version 4  April 1985
 * 
 *      A Package of Fortran Subprograms for the Fast Fourier
 *       Transform of Periodic and other Symmetric Sequences
 * 
 *                          by
 * 
 *                   Paul N Swarztrauber
 * 
 *   National Center for Atmospheric Research, Boulder, Colorado 80307,
 * 
 *    which is sponsored by the National Science Foundation
 * 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * 
 * There appears to be no explicit license for FFTPACK. However, the
 * package has been incorporated verbatim into a large number of software
 * systems over the years with numerous types of license without complaint
 * from the original author; therefore it would appear
 * that the code is effectively public domain. If you are in doubt,
 * however, you will need to contact the author or the  National Center
 * for Atmospheric Research to be sure.
 * 
 * All the changes from the original FFTPACK to the current file
 * fall under the following BSD-style open-source license:
 * 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * 
 * Copyright (c) 2005, Andrew Fernandes (andrew@fernandes.org);
 * All rights reserved.
 *  
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *  
 * - Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 * 
 * - Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 * 
 * - Neither the name of the North Carolina State University nor the
 * names of its contributors may be used to endorse or promote products
 * derived from this software without specific prior written permission.
 *  
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 * ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 */
#ifndef FFTPACK_H
#define FFTPACK_H

#ifdef __cplusplus
extern "C" {
#endif

/*

      The 'ifac' parameter needs to be as long as the number of
	  bits in a signed integer, minus one. Therefore if you are
	  on a 64-bit machine, the longest fft you can compute is
	  2^63-1, and 'ifac' should have length 63.

	  
*/

/* Single-precision */

/* Backward transform (synthesis) of a complex Fourier coefficient array */
void cfftb1(int *n, float *c, float *wsave, int *ifac);
/* Forward transform of a complex periodic sequence */
void cfftf1(int *n, float *c, float *wsave, int *ifac);
/* Initialization routine for (fftpack/cfftf) and (fftpack/cfftb) */
void cffti1(int *n, float *wsave, int *ifac);

/* Backward transform (synthesis) of a real Fourier coefficient array */
void rfftb1(int *n, float *r, float *wsave, int *ifac);
/* Forward transform of a real periodic sequence */
void rfftf1(int *n, float *r, float *wsave, int *ifac);
/* Initialization routine for (fftpack/rfftf) and (fftpack/rfftb) */
void rffti1(int *n, float *wsave, int *ifac);

/* Backward transform (synthesis) of a real Fourier coefficient array, simplified and slower version of (fftpack/rfftb) */
void ezfftb1(int *n, float *r, float *azero, float *a, float *b, float *wsave, int *ifac);
/* Forward transform of a real periodic sequence, simplified and slower version of (fftpack/rfftf) */
void ezfftf1(int *n, float *r, float *azero, float *a, float *b, float *wsave, int *ifac);
/* Initialization routine for (fftpack/ezfftf) and (fftpack/ezfftb) */
void ezffti1(int *n, float *wsave, int *ifac);

/* Backward cosine transform (synthesis) with odd wave numbers */
void cosqb1(int *n, float *x, float *wsave, int *ifac);
/* Forward cosine transform with odd wave numbers */
void cosqf1(int *n, float *x, float *wsave, int *ifac);
/* Initialization routine for (fftpack/cosqf) and (fftpack/cosqb) */
void cosqi1(int *n, float *wsave, int *ifac);

/* Discrete cosine transform of a real even sequence */
void cost1(int *n, float *x, float *wsave, int *ifac);
/* Initialization routine for (fftpack/cost) */
void costi1(int *n, float *wsave, int *ifac);

/* Backward sine transform (synthesis) with odd wave numbers */
void sinqb1(int *n, float *x, float *wsave, int *ifac);
/* Forward sine transform with odd wave numbers */
void sinqf1(int *n, float *x, float *wsave, int *ifac);
/* Initialization routine for (fftpack/sinqf) and (fftpack/sinqb) */
void sinqi1(int *n, float *wsave, int *ifac);

/* Discrete sine transform of a real odd sequence */
void sint1(int *n, float *x, float *wsave, int *ifac);
/* Initialization routine for (fftpack/sint) */
void sinti1(int *n, float *wsave, int *ifac);


/* Double-precision */

/* Backward transform (synthesis) of a complex Fourier coefficient array */
void cfftb2(int *n, double *c, double *wsave, int *ifac);
/* Forward transform of a complex periodic sequence */
void cfftf2(int *n, double *c, double *wsave, int *ifac);
/* Initialization routine for (fftpack/cfftf) and (fftpack/cfftb) */
void cffti2(int *n, double *wsave, int *ifac);

/* Backward transform (synthesis) of a real Fourier coefficient array */
void rfftb2(int *n, double *r, double *wsave, int *ifac);
/* Forward transform of a real periodic sequence */
void rfftf2(int *n, double *r, double *wsave, int *ifac);
/* Initialization routine for (fftpack/rfftf) and (fftpack/rfftb) */
void rffti2(int *n, double *wsave, int *ifac);

/* Backward transform (synthesis) of a real Fourier coefficient array, simplified and slower version of (fftpack/rfftb) */
void ezfftb2(int *n, double *r, double *azero, double *a, double *b, double *wsave, int *ifac);
/* Forward transform of a real periodic sequence, simplified and slower version of (fftpack/rfftf) */
void ezfftf2(int *n, double *r, double *azero, double *a, double *b, double *wsave, int *ifac);
/* Initialization routine for (fftpack/ezfftf) and (fftpack/ezfftb) */
void ezffti2(int *n, double *wsave, int *ifac);

/* Backward cosine transform (synthesis) with odd wave numbers */
void cosqb2(int *n, double *x, double *wsave, int *ifac);
/* Forward cosine transform with odd wave numbers */
void cosqf2(int *n, double *x, double *wsave, int *ifac);
/* Initialization routine for (fftpack/cosqf) and (fftpack/cosqb) */
void cosqi2(int *n, double *wsave, int *ifac);

/* Discrete cosine transform of a real even sequence */
void cost2(int *n, double *x, double *wsave, int *ifac);
/* Initialization routine for (fftpack/cost) */
void costi2(int *n, double *wsave, int *ifac);

/* Backward sine transform (synthesis) with odd wave numbers */
void sinqb2(int *n, double *x, double *wsave, int *ifac);
/* Forward sine transform with odd wave numbers */
void sinqf2(int *n, double *x, double *wsave, int *ifac);
/* Initialization routine for (fftpack/sinqf) and (fftpack/sinqb) */
void sinqi2(int *n, double *wsave, int *ifac);

/* Discrete sine transform of a real odd sequence */
void sint2(int *n, double *x, double *wsave, int *ifac);
/* Initialization routine for (fftpack/sint) */
void sinti2(int *n, double *wsave, int *ifac);


#ifdef __cplusplus
}
#endif



#ifdef __cplusplus

namespace fftpack{

/* Single-precision */

/* Backward transform (synthesis) of a complex Fourier coefficient array */
inline void cfftb(int *n, float *c, float *wsave, int *ifac){ cfftb1(n,c,wsave,ifac); }
/* Forward transform of a complex periodic sequence */
inline void cfftf(int *n, float *c, float *wsave, int *ifac){ cfftf1(n,c,wsave,ifac); }
/* Initialization routine for fftpack::cfftf and fftpack::cfftb */
inline void cffti(int *n, float *wsave, int *ifac){ cffti1(n,wsave,ifac); }

/* Backward transform (synthesis) of a real Fourier coefficient array */
inline void rfftb(int *n, float *r, float *wsave, int *ifac){ rfftb1(n,r,wsave,ifac); }
/* Forward transform of a real periodic sequence */
inline void rfftf(int *n, float *r, float *wsave, int *ifac){ rfftf1(n,r,wsave,ifac); }
/* Initialization routine for fftpack::rfftf and fftpack::rfftb */
inline void rffti(int *n, float *wsave, int *ifac){ rffti1(n,wsave,ifac); }

/* Backward transform (synthesis) of a real Fourier coefficient array, simplified and slower version of fftpack::rfftb */
inline void ezfftb(int *n, float *r, float *azero, float *a, float *b, float *wsave, int *ifac){ ezfftb1(n,r,azero,a,b,wsave,ifac); }
/* Forward transform of a real periodic sequence, simplified and slower version of fftpack::rfftf */
inline void ezfftf(int *n, float *r, float *azero, float *a, float *b, float *wsave, int *ifac){ ezfftf1(n,r,azero,a,b,wsave,ifac); }
/* Initialization routine for fftpack::ezfftf and fftpack::ezfftb */
inline void ezffti(int *n, float *wsave, int *ifac){ ezffti1(n,wsave,ifac); }

/* Backward cosine transform (synthesis) with odd wave numbers */
inline void cosqb(int *n, float *x, float *wsave, int *ifac){ cosqb1(n,x,wsave,ifac); }
/* Forward cosine transform with odd wave numbers */
inline void cosqf(int *n, float *x, float *wsave, int *ifac){ cosqf1(n,x,wsave,ifac); }
/* Initialization routine for fftpack::cosqf and fftpack::cosqb */
inline void cosqi(int *n, float *wsave, int *ifac){ cosqi1(n,wsave,ifac); }

/* Discrete cosine transform of a real even sequence */
inline void cost(int *n, float *x, float *wsave, int *ifac){ cost1(n,x,wsave,ifac); }
/* Initialization routine for fftpack::cost */
inline void costi(int *n, float *wsave, int *ifac){ costi1(n,wsave,ifac); }

/* Backward sine transform (synthesis) with odd wave numbers */
inline void sinqb(int *n, float *x, float *wsave, int *ifac){ sinqb1(n,x,wsave,ifac); }
/* Forward sine transform with odd wave numbers */
inline void sinqf(int *n, float *x, float *wsave, int *ifac){ sinqf1(n,x,wsave,ifac); }
/* Initialization routine for fftpack::sinqf and fftpack::sinqb */
inline void sinqi(int *n, float *wsave, int *ifac){ sinqi1(n,wsave,ifac); }

/* Discrete sine transform of a real odd sequence */
inline void sint(int *n, float *x, float *wsave, int *ifac){ sint1(n,x,wsave,ifac); }
/* Initialization routine for fftpack::sint */
inline void sinti(int *n, float *wsave, int *ifac){ sinti1(n,wsave,ifac); }


/* Double-precision */

/* Backward transform (synthesis) of a complex Fourier coefficient array */
inline void cfftb(int *n, double *c, double *wsave, int *ifac){ cfftb2(n,c,wsave,ifac); }
/* Forward transform of a complex periodic sequence */
inline void cfftf(int *n, double *c, double *wsave, int *ifac){ cfftf2(n,c,wsave,ifac); }
/* Initialization routine for fftpack::cfftf and fftpack::cfftb */
inline void cffti(int *n, double *wsave, int *ifac){ cffti2(n,wsave,ifac); }

/* Backward transform (synthesis) of a real Fourier coefficient array */
inline void rfftb(int *n, double *r, double *wsave, int *ifac){ rfftb2(n,r,wsave,ifac); }
/* Forward transform of a real periodic sequence */
inline void rfftf(int *n, double *r, double *wsave, int *ifac){ rfftf2(n,r,wsave,ifac); }
/* Initialization routine for fftpack::rfftf and fftpack::rfftb */
inline void rffti(int *n, double *wsave, int *ifac){ rffti2(n,wsave,ifac); }

/* Backward transform (synthesis) of a real Fourier coefficient array, simplified and slower version of fftpack::rfftb */
inline void ezfftb(int *n, double *r, double *azero, double *a, double *b, double *wsave, int *ifac){ ezfftb2(n,r,azero,a,b,wsave,ifac); }
/* Forward transform of a real periodic sequence, simplified and slower version of fftpack::rfftf */
inline void ezfftf(int *n, double *r, double *azero, double *a, double *b, double *wsave, int *ifac){ ezfftf2(n,r,azero,a,b,wsave,ifac); }
/* Initialization routine for fftpack::ezfftf and fftpack::ezfftb */
inline void ezffti(int *n, double *wsave, int *ifac){ ezffti2(n,wsave,ifac); }

/* Backward cosine transform (synthesis) with odd wave numbers */
inline void cosqb(int *n, double *x, double *wsave, int *ifac){ cosqb2(n,x,wsave,ifac); }
/* Forward cosine transform with odd wave numbers */
inline void cosqf(int *n, double *x, double *wsave, int *ifac){ cosqf2(n,x,wsave,ifac); }
/* Initialization routine for fftpack::cosqf and fftpack::cosqb */
inline void cosqi(int *n, double *wsave, int *ifac){ cosqi2(n,wsave,ifac); }

/* Discrete cosine transform of a real even sequence */
inline void cost(int *n, double *x, double *wsave, int *ifac){ cost2(n,x,wsave,ifac); }
/* Initialization routine for fftpack::cost */
inline void costi(int *n, double *wsave, int *ifac){ costi2(n,wsave,ifac); }

/* Backward sine transform (synthesis) with odd wave numbers */
inline void sinqb(int *n, double *x, double *wsave, int *ifac){ sinqb2(n,x,wsave,ifac); }
/* Forward sine transform with odd wave numbers */
inline void sinqf(int *n, double *x, double *wsave, int *ifac){ sinqf2(n,x,wsave,ifac); }
/* Initialization routine for fftpack::sinqf and fftpack::sinqb */
inline void sinqi(int *n, double *wsave, int *ifac){ sinqi2(n,wsave,ifac); }

/* Discrete sine transform of a real odd sequence */
inline void sint(int *n, double *x, double *wsave, int *ifac){ sint2(n,x,wsave,ifac); }
/* Initialization routine for fftpack::sint */
inline void sinti(int *n, double *wsave, int *ifac){ sinti2(n,wsave,ifac); }


}; // fftpack::


#endif


#endif /* ! FFTPACK_H */



// // -*- C++ -*-
// /***************************************************************************
//  *
//  * The IPPL Framework
//  * 
//  *
//  * Visit http://people.web.psi.ch/adelmann/ for more details
//  *
//  ***************************************************************************/
// 
// #ifndef IPPL_FFT_FFTPACK_FFT_H
// #define IPPL_FFT_FFTPACK_FFT_H
// 
// // include files
// #include "Utility/PAssert.h"
// #include "Utility/IpplInfo.h"
// 
// 
// /**************************************************************************
//  * fftpack_FFT.h:  Prototypes for accessing Fortran 1D FFT routines from
//  * Netlib, and definitions for templated class FFTPACK, which acts as an
//  * FFT engine for the FFT class, providing storage for trigonometric
//  * information and performing the 1D FFTs as needed.
//  **************************************************************************/
// 
// // For platforms that do Fortran symbols in all caps.
// #if defined(IPPL_T3E)
// 
// #define cffti_ CFFTI
// #define cfftf_ CFFTF
// #define cfftb_ CFFTB
// #define rffti_ RFFTI
// #define rfftf_ RFFTF
// #define rfftb_ RFFTB
// #define sinti_ SINTI
// #define sint_ SINT
// #define fcffti_ FCFFTI
// #define fcfftf_ FCFFTF
// #define fcfftb_ FCFFTB
// #define frffti_ FRFFTI
// #define frfftf_ FRFFTF
// #define frfftb_ FRFFTB
// #define fsinti_ FSINTI
// #define fsint_ FSINT
// 
// #endif
// 
// // For platforms that do Fortran symbols just like C symbols.
// #if defined(IPPL_SP2)
// 
// #define cffti_ cffti
// #define cfftf_ cfftf
// #define cfftb_ cfftb
// #define rffti_ rffti
// #define rfftf_ rfftf
// #define rfftb_ rfftb
// #define sinti_ sinti
// #define sint_ sint
// #define fcffti_ fcffti
// #define fcfftf_ fcfftf
// #define fcfftb_ fcfftb
// #define frffti_ frffti
// #define frfftf_ frfftf
// #define frfftb_ frfftb
// #define fsinti_ fsinti
// #define fsint_ fsint
// 
// #endif
// 
// // FFTPACK function prototypes for Fortran routines
// extern "C" {
//   // double-precision CC FFT
//   void cffti_(int& n, double& wsave);
//   void cfftf_(int& n, double& r, double& wsave);
//   void cfftb_(int& n, double& r, double& wsave);
//   // double-precision RC FFT
//   void rffti_(int& n, double& wsave);
//   void rfftf_(int& n, double& r, double& wsave);
//   void rfftb_(int& n, double& r, double& wsave);
//   // double-precision sine transform
//   void sinti_(int& n, double& wsave);
//   void sint_(int& n, double& r, double& wsave);
//   // single-precision CC FFT
//   void fcffti_(int& n, float& wsave);
//   void fcfftf_(int& n, float& r, float& wsave);
//   void fcfftb_(int& n, float& r, float& wsave);
//   // single-precision RC FFT
//   void frffti_(int& n, float& wsave);
//   void frfftf_(int& n, float& r, float& wsave);
//   void frfftb_(int& n, float& r, float& wsave);
//   // single-precision sine transform
//   void fsinti_(int& n, float& wsave);
//   void fsint_(int& n, float& r, float& wsave);
// }
// 
// 
// // FFTPACK_wrap provides static functions that wrap the Fortran functions
// // in a common interface.  We specialize this class on precision type.
// template <class T>
// class FFTPACK_wrap {};
// 
// // Specialization for float
// template <>
// class FFTPACK_wrap<float> {
// 
// public:
//   // interface functions used by class FFTPACK
// 
//   // initialization functions for CC FFT, RC FFT, and sine transform
//   static void ccffti(int n, float* wsave) { fcffti_(n, *wsave); }
//   static void rcffti(int n, float* wsave) { frffti_(n, *wsave); }
//   static void rrffti(int n, float* wsave) { fsinti_(n, *wsave); }
//   // forward and backward CC FFT
//   static void ccfftf(int n, float* r, float* wsave) { fcfftf_(n, *r, *wsave); }
//   static void ccfftb(int n, float* r, float* wsave) { fcfftb_(n, *r, *wsave); }
//   // forward and backward RC FFT
//   static void rcfftf(int n, float* r, float* wsave) { frfftf_(n, *r, *wsave); }
//   static void rcfftb(int n, float* r, float* wsave) { frfftb_(n, *r, *wsave); }
//   // sine transform
//   static void rrfft(int n, float* r, float* wsave) { fsint_(n, *r, *wsave); }
// 
// };
// 
// // Specialization for double
// template <>
// class FFTPACK_wrap<double> {
// 
// public:
//   // interface functions used by class FFTPACK
// 
//   // initialization functions for CC FFT, RC FFT, and sine transform
//   static void ccffti(int n, double* wsave) { cffti_(n, *wsave); }
//   static void rcffti(int n, double* wsave) { rffti_(n, *wsave); }
//   static void rrffti(int n, double* wsave) { sinti_(n, *wsave); }
//   // forward and backward CC FFT
//   static void ccfftf(int n, double* r, double* wsave) {cfftf_(n, *r, *wsave);}
//   static void ccfftb(int n, double* r, double* wsave) {cfftb_(n, *r, *wsave);}
//   // forward and backward RC FFT
//   static void rcfftf(int n, double* r, double* wsave) {rfftf_(n, *r, *wsave);}
//   static void rcfftb(int n, double* r, double* wsave) {rfftb_(n, *r, *wsave);}
//   // sine transform
//   static void rrfft(int n, double* r, double* wsave) { sint_(n, *r, *wsave); }
// 
// };
// 
// 
// // Definition of FFT engine class FFTPACK
// template <class T>
// class FFTPACK {
// 
// public:
// 
//   // definition of complex type
// #ifdef IPPL_HAS_TEMPLATED_COMPLEX
//   typedef complex<T> Complex_t;
// #else
//   typedef complex Complex_t;
// #endif
// 
//   // Trivial constructor.  Do the real work in setup function.
//   FFTPACK(void) {}
// 
//   // destructor
//   ~FFTPACK(void);
// 
//   // setup internal storage and prepare to perform FFTs
//   // inputs are number of dimensions to transform, the transform types,
//   // and the lengths of these dimensions.
//   void setup(unsigned numTransformDims, const int* transformTypes,
//              const int* axisLengths);
// 
//   // invoke FFT on complex data for given dimension and direction
//   void callFFT(unsigned transformDim, int direction, Complex_t* data);
// 
//   // invoke FFT on real data for given dimension and direction
//   void callFFT(unsigned transformDim, int direction, T* data);
// 
// private:
// 
//   unsigned numTransformDims_m;  // number of dimensions to transform
//   int* transformType_m;         // transform type for each dimension
//   int* axisLength_m;            // length of each transform dimension
//   T** trig_m;                   // trigonometric tables
// 
// };
// 
// 
// // Inline member function definitions
// 
// // destructor
// template <class T>
// inline
// FFTPACK<T>::~FFTPACK(void) {
//   // delete storage
//   for (unsigned d=0; d<numTransformDims_m; ++d)
//     delete [] trig_m[d];
//   delete [] trig_m;
//   delete [] transformType_m;
//   delete [] axisLength_m;
// }
// 
// // setup internal storage to prepare for FFTs
// template <class T>
// inline void
// FFTPACK<T>::setup(unsigned numTransformDims, const int* transformTypes,
//                   const int* axisLengths) {
// 
//   // store transform types and lengths for each transform dim
//   numTransformDims_m = numTransformDims;
//   transformType_m = new int[numTransformDims_m];
//   axisLength_m = new int[numTransformDims_m];
//   unsigned d;
//   for (d=0; d<numTransformDims_m; ++d) {
//     transformType_m[d] = transformTypes[d];
//     axisLength_m[d] = axisLengths[d];
//   }
// 
//   // allocate and initialize trig table
//   trig_m = new T*[numTransformDims_m];
//   for (d=0; d<numTransformDims_m; ++d) {
//     switch (transformType_m[d]) {
//     case 0:  // CC FFT
//       trig_m[d] = new T[4 * axisLength_m[d] + 15];
//       FFTPACK_wrap<T>::ccffti(axisLength_m[d], trig_m[d]);
//       break;
//     case 1:  // RC FFT
//       trig_m[d] = new T[2 * axisLength_m[d] + 15];
//       FFTPACK_wrap<T>::rcffti(axisLength_m[d], trig_m[d]);
//       break;
//     case 2:  // Sine transform
//       trig_m[d] = new T[static_cast<int>(2.5 * axisLength_m[d] + 0.5) + 15];
//       FFTPACK_wrap<T>::rrffti(axisLength_m[d], trig_m[d]);
//       break;
//     default:
//       ERRORMSG("Unknown transform type requested!!" << endl);
//       break;
//     }
//   }
// 
//   return;
// }
// 
// // invoke FFT on complex data for given dimension and direction
// template <class T>
// inline void
// FFTPACK<T>::callFFT(unsigned transformDim, int direction,
//                     FFTPACK<T>::Complex_t* data) {
// 
//   // check transform dimension and direction arguments
//   PAssert(transformDim<numTransformDims_m);
//   PAssert(direction==+1 || direction==-1);
// 
//   // cast complex number pointer to T* for calling Fortran routines
//   T* rdata = reinterpret_cast<T*>(data);
// 
//   // branch on transform type for this dimension
//   switch (transformType_m[transformDim]) {
//   case 0:  // CC FFT
//     if (direction == +1) {
//       // call forward complex-to-complex FFT
//       FFTPACK_wrap<T>::ccfftf(axisLength_m[transformDim], rdata,
//                               trig_m[transformDim]);
//     }
//     else {
//       // call backward complex-to-complex FFT
//       FFTPACK_wrap<T>::ccfftb(axisLength_m[transformDim], rdata,
//                               trig_m[transformDim]);
//     }
//     break;
//   case 1:  // RC FFT
//     if (direction == +1) {
//       // call forward real-to-complex FFT
//       FFTPACK_wrap<T>::rcfftf(axisLength_m[transformDim], rdata,
//                               trig_m[transformDim]);
//       // rearrange output to conform with SGI format for complex result
//       int clen = axisLength_m[transformDim]/2 + 1;
//       data[clen-1] = Complex_t(imag(data[clen-2]),0.0);
//       for (int i = clen-2; i > 0; --i)
//         data[i] = Complex_t(imag(data[i-1]),real(data[i]));
//       data[0] = Complex_t(real(data[0]),0.0);
//     }
//     else {                
//       // rearrange input to conform with Netlib format for complex modes
//       int clen = axisLength_m[transformDim]/2 + 1;
//       data[0] = Complex_t(real(data[0]),real(data[1]));
//       for (int i = 1; i < clen-1; ++i)
//         data[i] = Complex_t(imag(data[i]),real(data[i+1]));
//       // call backward complex-to-real FFT
//       FFTPACK_wrap<T>::rcfftb(axisLength_m[transformDim], rdata,
//                               trig_m[transformDim]);
//     }
//     break;
//   case 2:  // Sine transform
//     ERRORMSG("Input for real-to-real FFT should be real!!" << endl);
//     break;
//   default:
//     ERRORMSG("Unknown transform type requested!!" << endl);
//     break;
//   }
// 
//   return;
// }
// 
// // invoke FFT on real data for given dimension and direction
// template <class T>
// inline void
// FFTPACK<T>::callFFT(unsigned transformDim, int direction, T* data) {
// 
//   // check transform dimension and direction arguments
//   PAssert(transformDim<numTransformDims_m);
//   PAssert(direction==+1 || direction==-1);
// 
//   // branch on transform type for this dimension
//   switch (transformType_m[transformDim]) {
//   case 0:  // CC FFT
//     ERRORMSG("Input for complex-to-complex FFT should be complex!!" << endl);
//     break;
//   case 1:  // RC FFT
//     ERRORMSG("real-to-complex FFT uses complex input!!" << endl);
//     break;
//   case 2:  // Sine transform
//     // invoke the real-to-real transform on the data
//     FFTPACK_wrap<T>::rrfft(axisLength_m[transformDim], data,
//                            trig_m[transformDim]);
//     break;
//   default:
//     ERRORMSG("Unknown transform type requested!!" << endl);
//     break;
//   }
// 
//   return;
// }
// 
// 
// #endif // IPPL_FFT_FFTPACK_FFT_H
// 
// /***************************************************************************
//  * $RCSfile: fftpack_FFT.h,v $   $Author: adelmann $
//  * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:26 $
//  * IPPL_VERSION_ID: $Id: fftpack_FFT.h,v 1.1.1.1 2003/01/23 07:40:26 adelmann Exp $ 
//  ***************************************************************************/
// 

