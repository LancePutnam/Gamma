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

typedef int integer_t;            /* must be a signed integral type         */
typedef float real_t;            /* can be a float, double, or long-double */

/* Depending on the 'real_t' typedef, use the appropriate definition below  */
#define REAL_CONSTANT(x) x##f   // real_t <==> float                     */
/* #define REAL_CONSTANT(x) x      // real_t <==> double                    */
/* #define REAL_CONSTANT(x) x##l   // real_t <==> long double               */

/* Lastly, you will need to declare/define/alias the appropriate trig     */
/* functions 'cos/cosf/cosl' and 'sin/sinf/sinl' in fftpack.{c,cpp}       */


#ifdef __cplusplus
extern "C" {
#endif

/* Backward transform (synthesis) of a complex Fourier coefficient array */
void cfftb( integer_t *n, real_t *c, real_t *wsave, integer_t *ifac );
/* Forward transform of a complex periodic sequence */
void cfftf( integer_t *n, real_t *c, real_t *wsave, integer_t *ifac );
/* Initialization routine for (fftpack/cfftf) and (fftpack/cfftb) */
void cffti( integer_t *n, real_t *wsave, integer_t *ifac );


/* Backward cosine transform (synthesis) with odd wave numbers */
void cosqb( integer_t *n, real_t *x, real_t *wsave, integer_t *ifac );
/* Forward cosine transform with odd wave numbers */
void cosqf( integer_t *n, real_t *x, real_t *wsave, integer_t *ifac );
/* Initialization routine for (fftpack/cosqf) and (fftpack/cosqb) */
void cosqi( integer_t *n, real_t *wsave, integer_t *ifac );


/* Discrete cosine transform of a real even sequence */
void cost( integer_t *n, real_t *x, real_t *wsave, integer_t *ifac );
/* Initialization routine for (fftpack/cost) */
void costi( integer_t *n, real_t *wsave, integer_t *ifac );


/* Backward transform (synthesis) of a real Fourier coefficient array, simplified and slower version of (fftpack/rfftb) */
void ezfftb( integer_t *n, real_t *r, real_t *azero, real_t *a, real_t *b, real_t *wsave, integer_t *ifac );
/* Forward transform of a real periodic sequence, simplified and slower version of (fftpack/rfftf) */
void ezfftf( integer_t *n, real_t *r, real_t *azero, real_t *a, real_t *b, real_t *wsave, integer_t *ifac );
/* Initialization routine for (fftpack/ezfftf) and (fftpack/ezfftb) */
void ezffti( integer_t *n, real_t *wsave, integer_t *ifac );


/* Backward transform (synthesis) of a real Fourier coefficient array */
void rfftb( integer_t *n, real_t *r, real_t *wsave, integer_t *ifac );
/* Forward transform of a real periodic sequence */
void rfftf( integer_t *n, real_t *r, real_t *wsave, integer_t *ifac );
/* Initialization routine for (fftpack/rfftf) and (ftpack/rfftb) */
void rffti( integer_t *n, real_t *wsave, integer_t *ifac );


/* Backward sine transform (synthesis) with odd wave numbers */
void sinqb( integer_t *n, real_t *x, real_t *wsave, integer_t *ifac );
/* Forward sine transform with odd wave numbers */
void sinqf( integer_t *n, real_t *x, real_t *wsave, integer_t *ifac );
/* Initialization routtine for (fftpack/sinqf) and (fftpack/sinqb) */
void sinqi( integer_t *n, real_t *wsave, integer_t *ifac );


/* Discrete sine transform of a real odd sequence */
void sint( integer_t *n, real_t *x, real_t *wsave, integer_t *ifac );
/* Initialization routine for (fftpack/sint) */
void sinti( integer_t *n, real_t *wsave, integer_t *ifac );

#ifdef __cplusplus
}
#endif

#endif /* ! FFTPACK_H */
