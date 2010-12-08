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
