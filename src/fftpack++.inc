/* 
 * This file is based largely on the following software distribution:
 * 
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * 
 *                              FFTPACK
 * 
 * Reference                                                                                                                        
 *	P.N. Swarztrauber, Vectorizing the FFTs, in Parallel Computations
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

static void s_cfftb1(int *n, real_t *c__, real_t *ch, real_t *wa, int *ifac);
static void s_cfftf1(int *n, real_t *c__, real_t *ch, real_t *wa, int *ifac);
static void s_cffti1(int *n, real_t *wa, int *ifac);
static void s_cosqb1(int *n, real_t *x, real_t *w, real_t *xh, int *ifac);
static void s_cosqf1(int *n, real_t *x, real_t *w, real_t *xh, int *ifac);
static void s_ezfft1(int *n, real_t *wa, int *ifac);
static void s_passb(int *nac, int *ido, int *ip, int *l1, int *idl1, real_t *cc, real_t *c1, real_t *c2, real_t *ch, real_t *ch2, real_t *wa);
static void s_passb2(int *ido, int *l1, real_t *cc, real_t *ch, real_t *wa1);
static void s_passb3(int *ido, int *l1, real_t *cc, real_t *ch, real_t *wa1, real_t *wa2);
static void s_passb4(int *ido, int *l1, real_t *cc, real_t *ch, real_t *wa1, real_t *wa2, real_t *wa3);
static void s_passb5(int *ido, int *l1, real_t *cc, real_t *ch, real_t *wa1, real_t *wa2, real_t *wa3, real_t *wa4);
static void s_passf(int *nac, int *ido, int *ip, int *l1, int *idl1, real_t *cc, real_t *c1, real_t *c2, real_t *ch, real_t *ch2, real_t *wa);
static void s_passf2(int *ido, int *l1, real_t *cc, real_t *ch, real_t *wa1);
static void s_passf3(int *ido, int *l1, real_t *cc, real_t *ch, real_t *wa1, real_t *wa2);
static void s_passf4(int *ido, int *l1, real_t *cc, real_t *ch, real_t *wa1, real_t *wa2, real_t *wa3);
static void s_passf5(int *ido, int *l1, real_t *cc, real_t *ch, real_t *wa1, real_t *wa2, real_t *wa3, real_t *wa4);
static void s_radb2(int *ido, int *l1, real_t *cc, real_t *ch, real_t *wa1);
static void s_radb3(int *ido, int *l1, real_t *cc, real_t *ch, real_t *wa1, real_t *wa2);
static void s_radb4(int *ido, int *l1, real_t *cc, real_t *ch, real_t *wa1, real_t *wa2, real_t *wa3);
static void s_radb5(int *ido, int *l1, real_t *cc, real_t *ch, real_t *wa1, real_t *wa2, real_t *wa3, real_t *wa4);
static void s_radbg(int *ido, int *ip, int *l1, int *idl1, real_t *cc, real_t *c1, real_t *c2, real_t *ch, real_t *ch2, real_t *wa);
static void s_radf2(int *ido, int *l1, real_t *cc, real_t *ch, real_t *wa1);
static void s_radf3(int *ido, int *l1, real_t *cc, real_t *ch, real_t *wa1, real_t *wa2);
static void s_radf4(int *ido, int *l1, real_t *cc, real_t *ch, real_t *wa1, real_t *wa2, real_t *wa3);
static void s_radf5(int *ido, int *l1, real_t *cc, real_t *ch, real_t *wa1, real_t *wa2, real_t *wa3, real_t *wa4);
static void s_radfg(int *ido, int *ip, int *l1, int *idl1, real_t *cc, real_t *c1, real_t *c2, real_t *ch, real_t *ch2, real_t *wa);
static void s_rfftb1(int *n, real_t *c__, real_t *ch, real_t *wa, int *ifac);
static void s_rfftf1(int *n, real_t *c__, real_t *ch, real_t *wa, int *ifac);
static void s_rffti1(int *n, real_t *wa, int *ifac);
static void s_sint1(int *n, real_t *war, real_t *was, real_t *xh, real_t *x, int *ifac);

/* Subroutine */ void FUNC(cfftb)(int *n, real_t *c__, real_t *wsave, int *ifac)
{
	int iw1;

	/* Parameter adjustments */
	--ifac;
	--wsave;
	--c__;

	/* Function Body */
	if (*n == 1) {
	return;
	}
	iw1 = *n + *n + 1;
	s_cfftb1(n, &c__[1], &wsave[1], &wsave[iw1], &ifac[1]);
	return;
} /* cfftb_ */

/* Subroutine */ static void s_cfftb1(int *n, real_t *c__, real_t *ch, real_t *wa, int *ifac)
{
	/* System generated locals */
	int i__1;

	/* Local variables */
	int i__, k1, l1, l2, n2, na, nf, ip, iw, ix2, ix3, ix4, nac, ido, idl1, idot;

	/* Parameter adjustments */
	--ifac;
	--wa;
	--ch;
	--c__;

	/* Function Body */
	nf = ifac[2];
	na = 0;
	l1 = 1;
	iw = 1;
	i__1 = nf;
	for (k1 = 1; k1 <= i__1; ++k1) {
	ip = ifac[k1 + 2];
	l2 = ip * l1;
	ido = *n / l2;
	idot = ido + ido;
	idl1 = idot * l1;
	if (ip != 4) {
	    goto L103;
	}
	ix2 = iw + idot;
	ix3 = ix2 + idot;
	if (na != 0) {
	    goto L101;
	}
	s_passb4(&idot, &l1, &c__[1], &ch[1], &wa[iw], &wa[ix2], &wa[ix3]);
	goto L102;
L101:
	s_passb4(&idot, &l1, &ch[1], &c__[1], &wa[iw], &wa[ix2], &wa[ix3]);
L102:
	na = 1 - na;
	goto L115;
L103:
	if (ip != 2) {
	    goto L106;
	}
	if (na != 0) {
	    goto L104;
	}
	s_passb2(&idot, &l1, &c__[1], &ch[1], &wa[iw]);
	goto L105;
L104:
	s_passb2(&idot, &l1, &ch[1], &c__[1], &wa[iw]);
L105:
	na = 1 - na;
	goto L115;
L106:
	if (ip != 3) {
	    goto L109;
	}
	ix2 = iw + idot;
	if (na != 0) {
	    goto L107;
	}
	s_passb3(&idot, &l1, &c__[1], &ch[1], &wa[iw], &wa[ix2]);
	goto L108;
L107:
	s_passb3(&idot, &l1, &ch[1], &c__[1], &wa[iw], &wa[ix2]);
L108:
	na = 1 - na;
	goto L115;
L109:
	if (ip != 5) {
	    goto L112;
	}
	ix2 = iw + idot;
	ix3 = ix2 + idot;
	ix4 = ix3 + idot;
	if (na != 0) {
	    goto L110;
	}
	s_passb5(&idot, &l1, &c__[1], &ch[1], &wa[iw], &wa[ix2], &wa[ix3], &wa[ix4]);
	goto L111;
L110:
	s_passb5(&idot, &l1, &ch[1], &c__[1], &wa[iw], &wa[ix2], &wa[ix3], &wa[ix4]);
L111:
	na = 1 - na;
	goto L115;
L112:
	if (na != 0) {
	    goto L113;
	}
	s_passb(&nac, &idot, &ip, &l1, &idl1, &c__[1], &c__[1], &c__[1], &ch[1], &ch[1], &wa[iw]);
	goto L114;
L113:
	s_passb(&nac, &idot, &ip, &l1, &idl1, &ch[1], &ch[1], &ch[1], &c__[1], &c__[1], &wa[iw]);
L114:
	if (nac != 0) {
	    na = 1 - na;
	}
L115:
	l1 = l2;
	iw += (ip - 1) * idot;
/* L116: */
	}
	if (na == 0) {
	return;
	}
	n2 = *n + *n;
	i__1 = n2;
	for (i__ = 1; i__ <= i__1; ++i__) {
	c__[i__] = ch[i__];
/* L117: */
	}
	return;
} /* s_cfftb1 */

/* Subroutine */ void FUNC(cfftf)(int *n, real_t *c__, real_t *wsave, int *ifac)
{
	int iw1;

	/* Parameter adjustments */
	--ifac;
	--wsave;
	--c__;

	/* Function Body */
	if (*n == 1) {
	return;
	}
	iw1 = *n + *n + 1;
	s_cfftf1(n, &c__[1], &wsave[1], &wsave[iw1], &ifac[1]);
	return;
} /* cfftf_ */

/* Subroutine */ static void s_cfftf1(int *n, real_t *c__, real_t *ch, real_t *wa, int *ifac)
{
	/* System generated locals */
	int i__1;

	/* Local variables */
	int i__, k1, l1, l2, n2, na, nf, ip, iw, ix2, ix3, ix4, nac, ido, idl1, idot;

	/* Parameter adjustments */
	--ifac;
	--wa;
	--ch;
	--c__;

	/* Function Body */
	nf = ifac[2];
	na = 0;
	l1 = 1;
	iw = 1;
	i__1 = nf;
	for (k1 = 1; k1 <= i__1; ++k1) {
	ip = ifac[k1 + 2];
	l2 = ip * l1;
	ido = *n / l2;
	idot = ido + ido;
	idl1 = idot * l1;
	if (ip != 4) {
	    goto L103;
	}
	ix2 = iw + idot;
	ix3 = ix2 + idot;
	if (na != 0) {
	    goto L101;
	}
	s_passf4(&idot, &l1, &c__[1], &ch[1], &wa[iw], &wa[ix2], &wa[ix3]);
	goto L102;
L101:
	s_passf4(&idot, &l1, &ch[1], &c__[1], &wa[iw], &wa[ix2], &wa[ix3]);
L102:
	na = 1 - na;
	goto L115;
L103:
	if (ip != 2) {
	    goto L106;
	}
	if (na != 0) {
	    goto L104;
	}
	s_passf2(&idot, &l1, &c__[1], &ch[1], &wa[iw]);
	goto L105;
L104:
	s_passf2(&idot, &l1, &ch[1], &c__[1], &wa[iw]);
L105:
	na = 1 - na;
	goto L115;
L106:
	if (ip != 3) {
	    goto L109;
	}
	ix2 = iw + idot;
	if (na != 0) {
	    goto L107;
	}
	s_passf3(&idot, &l1, &c__[1], &ch[1], &wa[iw], &wa[ix2]);
	goto L108;
L107:
	s_passf3(&idot, &l1, &ch[1], &c__[1], &wa[iw], &wa[ix2]);
L108:
	na = 1 - na;
	goto L115;
L109:
	if (ip != 5) {
	    goto L112;
	}
	ix2 = iw + idot;
	ix3 = ix2 + idot;
	ix4 = ix3 + idot;
	if (na != 0) {
	    goto L110;
	}
	s_passf5(&idot, &l1, &c__[1], &ch[1], &wa[iw], &wa[ix2], &wa[ix3], &wa[ix4]);
	goto L111;
L110:
	s_passf5(&idot, &l1, &ch[1], &c__[1], &wa[iw], &wa[ix2], &wa[ix3], &wa[ix4]);
L111:
	na = 1 - na;
	goto L115;
L112:
	if (na != 0) {
	    goto L113;
	}
	s_passf(&nac, &idot, &ip, &l1, &idl1, &c__[1], &c__[1], &c__[1], &ch[1], &ch[1], &wa[iw]);
	goto L114;
L113:
	s_passf(&nac, &idot, &ip, &l1, &idl1, &ch[1], &ch[1], &ch[1], &c__[1], &c__[1], &wa[iw]);
L114:
	if (nac != 0) {
	    na = 1 - na;
	}
L115:
	l1 = l2;
	iw += (ip - 1) * idot;
/* L116: */
	}
	if (na == 0) {
	return;
	}
	n2 = *n + *n;
	i__1 = n2;
	for (i__ = 1; i__ <= i__1; ++i__) {
	c__[i__] = ch[i__];
/* L117: */
	}
	return;
} /* cfftf1_ */

/* Subroutine */ void FUNC(cffti)(int *n, real_t *wsave, int *ifac)
{
	int iw1;

	/* Parameter adjustments */
	--ifac;
	--wsave;

	/* Function Body */
	if (*n == 1) {
	return;
	}
	iw1 = *n + *n + 1;
	s_cffti1(n, &wsave[iw1], &ifac[1]);
	return;
} /* cffti_ */

/* Subroutine */ static void s_cffti1(int *n, real_t *wa, int *ifac)
{
	/* Initialized data */

	static int ntryh[4] = { 3,4,2,5 };

	/* System generated locals */
	int i__1, i__2, i__3;

	/* Local variables */
	int i__, j, i1, k1, l1, l2, ib;
	real_t fi;
	int ld, ii, nf, ip, nl, nq, nr;
	real_t arg;
	int ido, ipm;
	real_t tpi, argh;
	int idot, ntry=0;
	real_t argld;

	/* Parameter adjustments */
	--ifac;
	--wa;

	/* Function Body */
	nl = *n;
	nf = 0;
	j = 0;
L101:
	++j;
	if (j - 4 <= 0) {
	goto L102;
	} else {
	goto L103;
	}
L102:
	ntry = ntryh[j - 1];
	goto L104;
L103:
	ntry += 2;
L104:
	nq = nl / ntry;
	nr = nl - ntry * nq;
	if (nr != 0) {
	goto L101;
	} else {
	goto L105;
	}
L105:
	++nf;
	ifac[nf + 2] = ntry;
	nl = nq;
	if (ntry != 2) {
	goto L107;
	}
	if (nf == 1) {
	goto L107;
	}
	i__1 = nf;
	for (i__ = 2; i__ <= i__1; ++i__) {
	ib = nf - i__ + 2;
	ifac[ib + 2] = ifac[ib + 1];
/* L106: */
	}
	ifac[3] = 2;
L107:
	if (nl != 1) {
	goto L104;
	}
	ifac[1] = *n;
	ifac[2] = nf;
	tpi = POST(6.283185307179586476925286766559005768394338798750211619498891846);
	argh = tpi / (real_t) (*n);
	i__ = 2;
	l1 = 1;
	i__1 = nf;
	for (k1 = 1; k1 <= i__1; ++k1) {
	ip = ifac[k1 + 2];
	ld = 0;
	l2 = l1 * ip;
	ido = *n / l2;
	idot = ido + ido + 2;
	ipm = ip - 1;
	i__2 = ipm;
	for (j = 1; j <= i__2; ++j) {
	    i1 = i__;
	    wa[i__ - 1] = POST(1.0);
	    wa[i__] = POST(0.0);
	    ld += l1;
	    fi = POST(0.0);
	    argld = (real_t) ld * argh;
	    i__3 = idot;
	    for (ii = 4; ii <= i__3; ii += 2) {
		i__ += 2;
		fi += POST(1.0);
		arg = fi * argld;
		wa[i__ - 1] = POST(cos)(arg);
		wa[i__] = POST(sin)(arg);
/* L108: */
	    }
	    if (ip <= 5) {
		goto L109;
	    }
	    wa[i1 - 1] = wa[i__ - 1];
	    wa[i1] = wa[i__];
L109:
	    ;
	}
	l1 = l2;
/* L110: */
	}
	return;
} /* s_cffti1 */

/* Subroutine */ void FUNC(cosqb)(int *n, real_t *x, real_t *wsave, int *ifac)
{
	/* Initialized data */

	static real_t tsqrt2 = POST(2.82842712474619009760337744841939615713934375053896146353359476);

	/* System generated locals */
	int i__1;

	/* Local variables */
	real_t x1;

	/* Parameter adjustments */
	--ifac;
	--wsave;
	--x;

	/* Function Body */
	if ((i__1 = *n - 2) < 0) {
	goto L101;
	} else if (i__1 == 0) {
	goto L102;
	} else {
	goto L103;
	}
L101:
	x[1] *= POST(4.0);
	return;
L102:
	x1 = (x[1] + x[2]) * POST(4.0);
	x[2] = tsqrt2 * (x[1] - x[2]);
	x[1] = x1;
	return;
L103:
	s_cosqb1(n, &x[1], &wsave[1], &wsave[*n + 1], &ifac[1]);
	return;
} /* cosqb_ */

/* Subroutine */ static void s_cosqb1(int *n, real_t *x, real_t *w, real_t *xh, int *ifac)
{
	/* System generated locals */
	int i__1;

	/* Local variables */
	int i__, k, kc, np2, ns2;
	real_t xim1;
	int modn;

	/* Parameter adjustments */
	--ifac;
	--xh;
	--w;
	--x;

	/* Function Body */
	ns2 = (*n + 1) / 2;
	np2 = *n + 2;
	i__1 = *n;
	for (i__ = 3; i__ <= i__1; i__ += 2) {
	xim1 = x[i__ - 1] + x[i__];
	x[i__] -= x[i__ - 1];
	x[i__ - 1] = xim1;
/* L101: */
	}
	x[1] += x[1];
	modn = *n % 2;
	if (modn == 0) {
	x[*n] += x[*n];
	}
	FUNC(rfftb)(n, &x[1], &xh[1], &ifac[1]);
	i__1 = ns2;
	for (k = 2; k <= i__1; ++k) {
	kc = np2 - k;
	xh[k] = w[k - 1] * x[kc] + w[kc - 1] * x[k];
	xh[kc] = w[k - 1] * x[k] - w[kc - 1] * x[kc];
/* L102: */
	}
	if (modn == 0) {
	x[ns2 + 1] = w[ns2] * (x[ns2 + 1] + x[ns2 + 1]);
	}
	i__1 = ns2;
	for (k = 2; k <= i__1; ++k) {
	kc = np2 - k;
	x[k] = xh[k] + xh[kc];
	x[kc] = xh[k] - xh[kc];
/* L103: */
	}
	x[1] += x[1];
	return;
} /* cosqb1_ */

/* Subroutine */ void FUNC(cosqf)(int *n, real_t *x, real_t *wsave, int *ifac)
{
	/* Initialized data */

	static real_t sqrt2 = POST(1.41421356237309504880168872420969807856967187536948073176679738);

	/* System generated locals */
	int i__1;

	/* Local variables */
	real_t tsqx;

	/* Parameter adjustments */
	--ifac;
	--wsave;
	--x;

	/* Function Body */
	if ((i__1 = *n - 2) < 0) {
	goto L102;
	} else if (i__1 == 0) {
	goto L101;
	} else {
	goto L103;
	}
L101:
	tsqx = sqrt2 * x[2];
	x[2] = x[1] - tsqx;
	x[1] += tsqx;
L102:
	return;
L103:
	s_cosqf1(n, &x[1], &wsave[1], &wsave[*n + 1], &ifac[1]);
	return;
} /* cosqf_ */

/* Subroutine */ static void s_cosqf1(int *n, real_t *x, real_t *w, 
	real_t *xh, int *ifac)
{
	/* System generated locals */
	int i__1;

	/* Local variables */
	int i__, k, kc, np2, ns2;
	real_t xim1;
	int modn;

	/* Parameter adjustments */
	--ifac;
	--xh;
	--w;
	--x;

	/* Function Body */
	ns2 = (*n + 1) / 2;
	np2 = *n + 2;
	i__1 = ns2;
	for (k = 2; k <= i__1; ++k) {
	kc = np2 - k;
	xh[k] = x[k] + x[kc];
	xh[kc] = x[k] - x[kc];
/* L101: */
	}
	modn = *n % 2;
	if (modn == 0) {
	xh[ns2 + 1] = x[ns2 + 1] + x[ns2 + 1];
	}
	i__1 = ns2;
	for (k = 2; k <= i__1; ++k) {
	kc = np2 - k;
	x[k] = w[k - 1] * xh[kc] + w[kc - 1] * xh[k];
	x[kc] = w[k - 1] * xh[k] - w[kc - 1] * xh[kc];
/* L102: */
	}
	if (modn == 0) {
	x[ns2 + 1] = w[ns2] * xh[ns2 + 1];
	}
	FUNC(rfftf)(n, &x[1], &xh[1], &ifac[1]);
	i__1 = *n;
	for (i__ = 3; i__ <= i__1; i__ += 2) {
	xim1 = x[i__ - 1] - x[i__];
	x[i__] = x[i__ - 1] + x[i__];
	x[i__ - 1] = xim1;
/* L103: */
	}
	return;
} /* cosqf1_ */

/* Subroutine */ void FUNC(cosqi)(int *n, real_t *wsave, int *ifac)
{
	/* Initialized data */

	static real_t pih = POST(1.570796326794896619231321691639751442098584699687529104874722962);

	/* System generated locals */
	int i__1;

	/* Local variables */
	int k;
	real_t fk, dt;

	/* Parameter adjustments */
	--ifac;
	--wsave;

	/* Function Body */
	dt = pih / (real_t) (*n);
	fk = POST(0.0);
	i__1 = *n;
	for (k = 1; k <= i__1; ++k) {
	fk += POST(1.0);
	wsave[k] = POST(cos)(fk * dt);
/* L101: */
	}
	FUNC(rffti)(n, &wsave[*n + 1], &ifac[1]);
	return;
} /* cosqi_ */

/* Subroutine */ void FUNC(cost)(int *n, real_t *x, real_t *wsave, int *ifac)
{
	/* System generated locals */
	int i__1;

	/* Local variables */
	int i__, k;
	real_t c1, t1, t2;
	int kc;
	real_t xi;
	int nm1, np1;
	real_t x1h;
	int ns2;
	real_t tx2, x1p3, xim2;
	int modn;

	/* Parameter adjustments */
	--ifac;
	--wsave;
	--x;

	/* Function Body */
	nm1 = *n - 1;
	np1 = *n + 1;
	ns2 = *n / 2;
	if ((i__1 = *n - 2) < 0) {
	goto L106;
	} else if (i__1 == 0) {
	goto L101;
	} else {
	goto L102;
	}
L101:
	x1h = x[1] + x[2];
	x[2] = x[1] - x[2];
	x[1] = x1h;
	return;
L102:
	if (*n > 3) {
	goto L103;
	}
	x1p3 = x[1] + x[3];
	tx2 = x[2] + x[2];
	x[2] = x[1] - x[3];
	x[1] = x1p3 + tx2;
	x[3] = x1p3 - tx2;
	return;
L103:
	c1 = x[1] - x[*n];
	x[1] += x[*n];
	i__1 = ns2;
	for (k = 2; k <= i__1; ++k) {
	kc = np1 - k;
	t1 = x[k] + x[kc];
	t2 = x[k] - x[kc];
	c1 += wsave[kc] * t2;
	t2 = wsave[k] * t2;
	x[k] = t1 - t2;
	x[kc] = t1 + t2;
/* L104: */
	}
	modn = *n % 2;
	if (modn != 0) {
	x[ns2 + 1] += x[ns2 + 1];
	}
	FUNC(rfftf)(&nm1, &x[1], &wsave[*n + 1], &ifac[1]);
	xim2 = x[2];
	x[2] = c1;
	i__1 = *n;
	for (i__ = 4; i__ <= i__1; i__ += 2) {
	xi = x[i__];
	x[i__] = x[i__ - 2] - x[i__ - 1];
	x[i__ - 1] = xim2;
	xim2 = xi;
/* L105: */
	}
	if (modn != 0) {
	x[*n] = xim2;
	}
L106:
	return;
} /* cost_ */

/* Subroutine */ void FUNC(costi)(int *n, real_t *wsave, int *ifac)
{
	/* Initialized data */

	static real_t pi = POST(3.141592653589793238462643383279502884197169399375158209749445923);

	/* System generated locals */
	int i__1;

	/* Local variables */
	int k, kc;
	real_t fk, dt;
	int nm1, np1, ns2;

	/* Parameter adjustments */
	--ifac;
	--wsave;

	/* Function Body */
	if (*n <= 3) {
	return;
	}
	nm1 = *n - 1;
	np1 = *n + 1;
	ns2 = *n / 2;
	dt = pi / (real_t) nm1;
	fk = POST(0.0);
	i__1 = ns2;
	for (k = 2; k <= i__1; ++k) {
	kc = np1 - k;
	fk += POST(1.0);
	wsave[k] = POST(sin)(fk * dt) * POST(2.0);
	wsave[kc] = POST(cos)(fk * dt) * POST(2.0);
/* L101: */
	}
	FUNC(rffti)(&nm1, &wsave[*n + 1], &ifac[1]);
	return;
} /* costi_ */

/* Subroutine */ static void s_ezfft1(int *n, real_t *wa, int *ifac)
{
	/* Initialized data */

	static int ntryh[4] = { 4,2,3,5 };
	static real_t tpi = POST(6.283185307179586476925286766559005768394338798750211419498891846);

	/* System generated locals */
	int i__1, i__2, i__3;

	/* Local variables */
	int i__, j, k1, l1, l2, ib, ii, nf, ip, nl, is, nq, nr;
	real_t ch1, sh1;
	int ido, ipm;
	real_t dch1, ch1h, arg1, dsh1;
	int nfm1;
	real_t argh;
	int ntry=0;

	/* Parameter adjustments */
	--ifac;
	--wa;

	/* Function Body */
	nl = *n;
	nf = 0;
	j = 0;
L101:
	++j;
	if (j - 4 <= 0) {
	goto L102;
	} else {
	goto L103;
	}
L102:
	ntry = ntryh[j - 1];
	goto L104;
L103:
	ntry += 2;
L104:
	nq = nl / ntry;
	nr = nl - ntry * nq;
	if (nr != 0) {
	goto L101;
	} else {
	goto L105;
	}
L105:
	++nf;
	ifac[nf + 2] = ntry;
	nl = nq;
	if (ntry != 2) {
	goto L107;
	}
	if (nf == 1) {
	goto L107;
	}
	i__1 = nf;
	for (i__ = 2; i__ <= i__1; ++i__) {
	ib = nf - i__ + 2;
	ifac[ib + 2] = ifac[ib + 1];
/* L106: */
	}
	ifac[3] = 2;
L107:
	if (nl != 1) {
	goto L104;
	}
	ifac[1] = *n;
	ifac[2] = nf;
	argh = tpi / (real_t) (*n);
	is = 0;
	nfm1 = nf - 1;
	l1 = 1;
	if (nfm1 == 0) {
	return;
	}
	i__1 = nfm1;
	for (k1 = 1; k1 <= i__1; ++k1) {
	ip = ifac[k1 + 2];
	l2 = l1 * ip;
	ido = *n / l2;
	ipm = ip - 1;
	arg1 = (real_t) l1 * argh;
	ch1 = POST(1.0);
	sh1 = POST(0.0);
	dch1 = POST(cos)(arg1);
	dsh1 = POST(sin)(arg1);
	i__2 = ipm;
	for (j = 1; j <= i__2; ++j) {
	    ch1h = dch1 * ch1 - dsh1 * sh1;
	    sh1 = dch1 * sh1 + dsh1 * ch1;
	    ch1 = ch1h;
	    i__ = is + 2;
	    wa[i__ - 1] = ch1;
	    wa[i__] = sh1;
	    if (ido < 5) {
		goto L109;
	    }
	    i__3 = ido;
	    for (ii = 5; ii <= i__3; ii += 2) {
		i__ += 2;
		wa[i__ - 1] = ch1 * wa[i__ - 3] - sh1 * wa[i__ - 2];
		wa[i__] = ch1 * wa[i__ - 2] + sh1 * wa[i__ - 3];
/* L108: */
	    }
L109:
	    is += ido;
/* L110: */
	}
	l1 = l2;
/* L111: */
	}
	return;
} /* ezfft1_ */

/* Subroutine */ void FUNC(ezfftb)(int *n, real_t *r__, real_t *azero, 
	real_t *a, real_t *b, real_t *wsave, int *ifac)
{
	/* System generated locals */
	int i__1;

	/* Local variables */
	int i__, ns2;

	/* Parameter adjustments */
	--ifac;
	--wsave;
	--b;
	--a;
	--r__;

	/* Function Body */
	if ((i__1 = *n - 2) < 0) {
	goto L101;
	} else if (i__1 == 0) {
	goto L102;
	} else {
	goto L103;
	}
L101:
	r__[1] = *azero;
	return;
L102:
	r__[1] = *azero + a[1];
	r__[2] = *azero - a[1];
	return;
L103:
	ns2 = (*n - 1) / 2;
	i__1 = ns2;
	for (i__ = 1; i__ <= i__1; ++i__) {
	r__[i__ * 2] = a[i__] * POST(0.5);
	r__[(i__ << 1) + 1] = b[i__] * POST(-0.5);
/* L104: */
	}
	r__[1] = *azero;
	if (*n % 2 == 0) {
	r__[*n] = a[ns2 + 1];
	}
	FUNC(rfftb)(n, &r__[1], &wsave[*n + 1], &ifac[1]);
	return;
} /* ezfftb_ */

/* Subroutine */ void FUNC(ezfftf)(int *n, real_t *r__, real_t *azero, 
	real_t *a, real_t *b, real_t *wsave, int *ifac)
{
	/* System generated locals */
	int i__1;

	/* Local variables */
	int i__;
	real_t cf;
	int ns2;
	real_t cfm;
	int ns2m;

/*                       VERSION 3  JUNE 1979 */

	/* Parameter adjustments */
	--ifac;
	--wsave;
	--b;
	--a;
	--r__;

	/* Function Body */
	if ((i__1 = *n - 2) < 0) {
	goto L101;
	} else if (i__1 == 0) {
	goto L102;
	} else {
	goto L103;
	}
L101:
	*azero = r__[1];
	return;
L102:
	*azero = (r__[1] + r__[2]) * POST(0.5);
	a[1] = (r__[1] - r__[2]) * POST(0.5);
	return;
L103:
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	wsave[i__] = r__[i__];
/* L104: */
	}
	FUNC(rfftf)(n, &wsave[1], &wsave[*n + 1], &ifac[1]);
	cf = 2.0 / (real_t) (*n);
	cfm = -cf;
	*azero = cf * 0.5 * wsave[1];
	ns2 = (*n + 1) / 2;
	ns2m = ns2 - 1;
	i__1 = ns2m;
	for (i__ = 1; i__ <= i__1; ++i__) {
	a[i__] = cf * wsave[i__ * 2];
	b[i__] = cfm * wsave[(i__ << 1) + 1];
/* L105: */
	}
	if (*n % 2 == 1) {
	return;
	}
	a[ns2] = cf * 0.5 * wsave[*n];
	b[ns2] = POST(0.0);
	return;
} /* ezfftf_ */

/* Subroutine */ void FUNC(ezffti)(int *n, real_t *wsave, int *ifac)
{
	/* Parameter adjustments */
	--ifac;
	--wsave;

	/* Function Body */
	if (*n == 1) {
	return;
	}
	s_ezfft1(n, &wsave[(*n << 1) + 1], &ifac[1]);
	return;
} /* ezffti_ */

/* Subroutine */ static void s_passb(int *nac, int *ido, int *ip, int *
	l1, int *idl1, real_t *cc, real_t *c1, real_t *c2, 
	real_t *ch, real_t *ch2, real_t *wa)
{
	/* System generated locals */
	int ch_dim1, ch_dim2, ch_offset, cc_dim1, cc_dim2, cc_offset, c1_dim1,
		c1_dim2, c1_offset, c2_dim1, c2_offset, ch2_dim1, ch2_offset, 
		i__1, i__2, i__3;

	/* Local variables */
	int i__, j, k, l, jc, lc, ik, nt, idj, idl, inc, idp;
	real_t wai, war;
	int ipp2, idij, idlj, idot, ipph;

	/* Parameter adjustments */
	ch_dim1 = *ido;
	ch_dim2 = *l1;
	ch_offset = 1 + ch_dim1 * (1 + ch_dim2);
	ch -= ch_offset;
	c1_dim1 = *ido;
	c1_dim2 = *l1;
	c1_offset = 1 + c1_dim1 * (1 + c1_dim2);
	c1 -= c1_offset;
	cc_dim1 = *ido;
	cc_dim2 = *ip;
	cc_offset = 1 + cc_dim1 * (1 + cc_dim2);
	cc -= cc_offset;
	ch2_dim1 = *idl1;
	ch2_offset = 1 + ch2_dim1;
	ch2 -= ch2_offset;
	c2_dim1 = *idl1;
	c2_offset = 1 + c2_dim1;
	c2 -= c2_offset;
	--wa;

	/* Function Body */
	idot = *ido / 2;
	nt = *ip * *idl1;
	ipp2 = *ip + 2;
	ipph = (*ip + 1) / 2;
	idp = *ip * *ido;

	if (*ido < *l1) {
	goto L106;
	}
	i__1 = ipph;
	for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    i__3 = *ido;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		ch[i__ + (k + j * ch_dim2) * ch_dim1] = cc[i__ + (j + k * 
			cc_dim2) * cc_dim1] + cc[i__ + (jc + k * cc_dim2) * cc_dim1];
		ch[i__ + (k + jc * ch_dim2) * ch_dim1] = cc[i__ + (j + k * 
			cc_dim2) * cc_dim1] - cc[i__ + (jc + k * cc_dim2) * cc_dim1];
/* L101: */
	    }
/* L102: */
	}
/* L103: */
	}
	i__1 = *l1;
	for (k = 1; k <= i__1; ++k) {
	i__2 = *ido;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    ch[i__ + (k + ch_dim2) * ch_dim1] = cc[i__ + (k * cc_dim2 + 1) * cc_dim1];
/* L104: */
	}
/* L105: */
	}
	goto L112;
L106:
	i__1 = ipph;
	for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	i__2 = *ido;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = *l1;
	    for (k = 1; k <= i__3; ++k) {
		ch[i__ + (k + j * ch_dim2) * ch_dim1] = cc[i__ + (j + k * 
			cc_dim2) * cc_dim1] + cc[i__ + (jc + k * cc_dim2) * cc_dim1];
		ch[i__ + (k + jc * ch_dim2) * ch_dim1] = cc[i__ + (j + k * 
			cc_dim2) * cc_dim1] - cc[i__ + (jc + k * cc_dim2) * cc_dim1];
/* L107: */
	    }
/* L108: */
	}
/* L109: */
	}
	i__1 = *ido;
	for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    ch[i__ + (k + ch_dim2) * ch_dim1] = cc[i__ + (k * cc_dim2 + 1) * cc_dim1];
/* L110: */
	}
/* L111: */
	}
L112:
	idl = 2 - *ido;
	inc = 0;
	i__1 = ipph;
	for (l = 2; l <= i__1; ++l) {
	lc = ipp2 - l;
	idl += *ido;
	i__2 = *idl1;
	for (ik = 1; ik <= i__2; ++ik) {
	    c2[ik + l * c2_dim1] = ch2[ik + ch2_dim1] + wa[idl - 1] * ch2[ik + (ch2_dim1 << 1)];
	    c2[ik + lc * c2_dim1] = wa[idl] * ch2[ik + *ip * ch2_dim1];
/* L113: */
	}
	idlj = idl;
	inc += *ido;
	i__2 = ipph;
	for (j = 3; j <= i__2; ++j) {
	    jc = ipp2 - j;
	    idlj += inc;
	    if (idlj > idp) {
		idlj -= idp;
	    }
	    war = wa[idlj - 1];
	    wai = wa[idlj];
	    i__3 = *idl1;
	    for (ik = 1; ik <= i__3; ++ik) {
		c2[ik + l * c2_dim1] += war * ch2[ik + j * ch2_dim1];
		c2[ik + lc * c2_dim1] += wai * ch2[ik + jc * ch2_dim1];
/* L114: */
	    }
/* L115: */
	}
/* L116: */
	}
	i__1 = ipph;
	for (j = 2; j <= i__1; ++j) {
	i__2 = *idl1;
	for (ik = 1; ik <= i__2; ++ik) {
	    ch2[ik + ch2_dim1] += ch2[ik + j * ch2_dim1];
/* L117: */
	}
/* L118: */
	}
	i__1 = ipph;
	for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	i__2 = *idl1;
	for (ik = 2; ik <= i__2; ik += 2) {
	    ch2[ik - 1 + j * ch2_dim1] = c2[ik - 1 + j * c2_dim1] - c2[ik + jc * c2_dim1];
	    ch2[ik - 1 + jc * ch2_dim1] = c2[ik - 1 + j * c2_dim1] + c2[ik + jc * c2_dim1];
	    ch2[ik + j * ch2_dim1] = c2[ik + j * c2_dim1] + c2[ik - 1 + jc * c2_dim1];
	    ch2[ik + jc * ch2_dim1] = c2[ik + j * c2_dim1] - c2[ik - 1 + jc * c2_dim1];
/* L119: */
	}
/* L120: */
	}
	*nac = 1;
	if (*ido == 2) {
	return;
	}
	*nac = 0;
	i__1 = *idl1;
	for (ik = 1; ik <= i__1; ++ik) {
	c2[ik + c2_dim1] = ch2[ik + ch2_dim1];
/* L121: */
	}
	i__1 = *ip;
	for (j = 2; j <= i__1; ++j) {
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    c1[(k + j * c1_dim2) * c1_dim1 + 1] = ch[(k + j * ch_dim2) * ch_dim1 + 1];
	    c1[(k + j * c1_dim2) * c1_dim1 + 2] = ch[(k + j * ch_dim2) * ch_dim1 + 2];
/* L122: */
	}
/* L123: */
	}
	if (idot > *l1) {
	goto L127;
	}
	idij = 0;
	i__1 = *ip;
	for (j = 2; j <= i__1; ++j) {
	idij += 2;
	i__2 = *ido;
	for (i__ = 4; i__ <= i__2; i__ += 2) {
	    idij += 2;
	    i__3 = *l1;
	    for (k = 1; k <= i__3; ++k) {
		c1[i__ - 1 + (k + j * c1_dim2) * c1_dim1] = wa[idij - 1] * ch[
			i__ - 1 + (k + j * ch_dim2) * ch_dim1] - wa[idij] * 
			ch[i__ + (k + j * ch_dim2) * ch_dim1];
		c1[i__ + (k + j * c1_dim2) * c1_dim1] = wa[idij - 1] * ch[i__ 
			+ (k + j * ch_dim2) * ch_dim1] + wa[idij] * ch[i__ - 
			1 + (k + j * ch_dim2) * ch_dim1];
/* L124: */
	    }
/* L125: */
	}
/* L126: */
	}
	return;
L127:
	idj = 2 - *ido;
	i__1 = *ip;
	for (j = 2; j <= i__1; ++j) {
	idj += *ido;
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    idij = idj;
	    i__3 = *ido;
	    for (i__ = 4; i__ <= i__3; i__ += 2) {
		idij += 2;
		c1[i__ - 1 + (k + j * c1_dim2) * c1_dim1] = wa[idij - 1] * ch[
			i__ - 1 + (k + j * ch_dim2) * ch_dim1] - wa[idij] * 
			ch[i__ + (k + j * ch_dim2) * ch_dim1];
		c1[i__ + (k + j * c1_dim2) * c1_dim1] = wa[idij - 1] * ch[i__ 
			+ (k + j * ch_dim2) * ch_dim1] + wa[idij] * ch[i__ - 
			1 + (k + j * ch_dim2) * ch_dim1];
/* L128: */
	    }
/* L129: */
	}
/* L130: */
	}
	return;
} /* passb_ */

/* Subroutine */ static void s_passb2(int *ido, int *l1, real_t *cc, 
	real_t *ch, real_t *wa1)
{
	/* System generated locals */
	int cc_dim1, cc_offset, ch_dim1, ch_dim2, ch_offset, i__1, i__2;

	/* Local variables */
	int i__, k;
	real_t ti2, tr2;

	/* Parameter adjustments */
	ch_dim1 = *ido;
	ch_dim2 = *l1;
	ch_offset = 1 + ch_dim1 * (1 + ch_dim2);
	ch -= ch_offset;
	cc_dim1 = *ido;
	cc_offset = 1 + cc_dim1 * 3;
	cc -= cc_offset;
	--wa1;

	/* Function Body */
	if (*ido > 2) {
	goto L102;
	}
	i__1 = *l1;
	for (k = 1; k <= i__1; ++k) {
	ch[(k + ch_dim2) * ch_dim1 + 1] = cc[((k << 1) + 1) * cc_dim1 + 1] + cc[((k << 1) + 2) * cc_dim1 + 1];
	ch[(k + (ch_dim2 << 1)) * ch_dim1 + 1] = cc[((k << 1) + 1) * cc_dim1 + 1] - cc[((k << 1) + 2) * cc_dim1 + 1];
	ch[(k + ch_dim2) * ch_dim1 + 2] = cc[((k << 1) + 1) * cc_dim1 + 2] + cc[((k << 1) + 2) * cc_dim1 + 2];
	ch[(k + (ch_dim2 << 1)) * ch_dim1 + 2] = cc[((k << 1) + 1) * cc_dim1 + 2] - cc[((k << 1) + 2) * cc_dim1 + 2];
/* L101: */
	}
	return;
L102:
	i__1 = *l1;
	for (k = 1; k <= i__1; ++k) {
	i__2 = *ido;
	for (i__ = 2; i__ <= i__2; i__ += 2) {
	    ch[i__ - 1 + (k + ch_dim2) * ch_dim1] = cc[i__ - 1 + ((k << 1) + 1) * cc_dim1] + cc[i__ - 1 + ((k << 1) + 2) * cc_dim1];
	    tr2 = cc[i__ - 1 + ((k << 1) + 1) * cc_dim1] - cc[i__ - 1 + ((k <<1) + 2) * cc_dim1];
	    ch[i__ + (k + ch_dim2) * ch_dim1] = cc[i__ + ((k << 1) + 1) * cc_dim1] + cc[i__ + ((k << 1) + 2) * cc_dim1];
	    ti2 = cc[i__ + ((k << 1) + 1) * cc_dim1] - cc[i__ + ((k << 1) + 2)* cc_dim1];
	    ch[i__ + (k + (ch_dim2 << 1)) * ch_dim1] = wa1[i__ - 1] * ti2 + wa1[i__] * tr2;
	    ch[i__ - 1 + (k + (ch_dim2 << 1)) * ch_dim1] = wa1[i__ - 1] * tr2 - wa1[i__] * ti2;
/* L103: */
	}
/* L104: */
	}
	return;
} /* passb2_ */

/* Subroutine */ static void s_passb3(int *ido, int *l1, real_t *cc, 
	real_t *ch, real_t *wa1, real_t *wa2)
{
	/* Initialized data */

	static real_t taur = POST(-0.5);
	static real_t taui = POST(0.8660254037844386467637231707529361834710262690519031402790348975);

	/* System generated locals */
	int cc_dim1, cc_offset, ch_dim1, ch_dim2, ch_offset, i__1, i__2;

	/* Local variables */
	int i__, k;
	real_t ci2, ci3, di2, di3, cr2, cr3, dr2, dr3, ti2, tr2;

	/* Parameter adjustments */
	ch_dim1 = *ido;
	ch_dim2 = *l1;
	ch_offset = 1 + ch_dim1 * (1 + ch_dim2);
	ch -= ch_offset;
	cc_dim1 = *ido;
	cc_offset = 1 + (cc_dim1 << 2);
	cc -= cc_offset;
	--wa1;
	--wa2;

	/* Function Body */
	if (*ido != 2) {
	goto L102;
	}
	i__1 = *l1;
	for (k = 1; k <= i__1; ++k) {
	tr2 = cc[(k * 3 + 2) * cc_dim1 + 1] + cc[(k * 3 + 3) * cc_dim1 + 1];
	cr2 = cc[(k * 3 + 1) * cc_dim1 + 1] + taur * tr2;
	ch[(k + ch_dim2) * ch_dim1 + 1] = cc[(k * 3 + 1) * cc_dim1 + 1] + tr2;
	ti2 = cc[(k * 3 + 2) * cc_dim1 + 2] + cc[(k * 3 + 3) * cc_dim1 + 2];
	ci2 = cc[(k * 3 + 1) * cc_dim1 + 2] + taur * ti2;
	ch[(k + ch_dim2) * ch_dim1 + 2] = cc[(k * 3 + 1) * cc_dim1 + 2] + ti2;
	cr3 = taui * (cc[(k * 3 + 2) * cc_dim1 + 1] - cc[(k * 3 + 3) * cc_dim1 + 1]);
	ci3 = taui * (cc[(k * 3 + 2) * cc_dim1 + 2] - cc[(k * 3 + 3) * cc_dim1 + 2]);
	ch[(k + (ch_dim2 << 1)) * ch_dim1 + 1] = cr2 - ci3;
	ch[(k + ch_dim2 * 3) * ch_dim1 + 1] = cr2 + ci3;
	ch[(k + (ch_dim2 << 1)) * ch_dim1 + 2] = ci2 + cr3;
	ch[(k + ch_dim2 * 3) * ch_dim1 + 2] = ci2 - cr3;
/* L101: */
	}
	return;
L102:
	i__1 = *l1;
	for (k = 1; k <= i__1; ++k) {
	i__2 = *ido;
	for (i__ = 2; i__ <= i__2; i__ += 2) {
	    tr2 = cc[i__ - 1 + (k * 3 + 2) * cc_dim1] + cc[i__ - 1 + (k * 3 + 3) * cc_dim1];
	    cr2 = cc[i__ - 1 + (k * 3 + 1) * cc_dim1] + taur * tr2;
	    ch[i__ - 1 + (k + ch_dim2) * ch_dim1] = cc[i__ - 1 + (k * 3 + 1) *cc_dim1] + tr2;
	    ti2 = cc[i__ + (k * 3 + 2) * cc_dim1] + cc[i__ + (k * 3 + 3) * cc_dim1];
	    ci2 = cc[i__ + (k * 3 + 1) * cc_dim1] + taur * ti2;
	    ch[i__ + (k + ch_dim2) * ch_dim1] = cc[i__ + (k * 3 + 1) * cc_dim1] + ti2;
	    cr3 = taui * (cc[i__ - 1 + (k * 3 + 2) * cc_dim1] - cc[i__ - 1 + (k * 3 + 3) * cc_dim1]);
	    ci3 = taui * (cc[i__ + (k * 3 + 2) * cc_dim1] - cc[i__ + (k * 3 + 3) * cc_dim1]);
	    dr2 = cr2 - ci3;
	    dr3 = cr2 + ci3;
	    di2 = ci2 + cr3;
	    di3 = ci2 - cr3;
	    ch[i__ + (k + (ch_dim2 << 1)) * ch_dim1] = wa1[i__ - 1] * di2 + wa1[i__] * dr2;
	    ch[i__ - 1 + (k + (ch_dim2 << 1)) * ch_dim1] = wa1[i__ - 1] * dr2 - wa1[i__] * di2;
	    ch[i__ + (k + ch_dim2 * 3) * ch_dim1] = wa2[i__ - 1] * di3 + wa2[i__] * dr3;
	    ch[i__ - 1 + (k + ch_dim2 * 3) * ch_dim1] = wa2[i__ - 1] * dr3 - wa2[i__] * di3;
/* L103: */
	}
/* L104: */
	}
	return;
} /* passb3_ */

/* Subroutine */ static void s_passb4(int *ido, int *l1, real_t *cc, 
	real_t *ch, real_t *wa1, real_t *wa2, real_t *wa3)
{
	/* System generated locals */
	int cc_dim1, cc_offset, ch_dim1, ch_dim2, ch_offset, i__1, i__2;

	/* Local variables */
	int i__, k;
	real_t ci2, ci3, ci4, cr2, cr3, cr4, ti1, ti2, ti3, ti4, tr1, tr2, tr3, tr4;

	/* Parameter adjustments */
	ch_dim1 = *ido;
	ch_dim2 = *l1;
	ch_offset = 1 + ch_dim1 * (1 + ch_dim2);
	ch -= ch_offset;
	cc_dim1 = *ido;
	cc_offset = 1 + cc_dim1 * 5;
	cc -= cc_offset;
	--wa1;
	--wa2;
	--wa3;

	/* Function Body */
	if (*ido != 2) {
	goto L102;
	}
	i__1 = *l1;
	for (k = 1; k <= i__1; ++k) {
	ti1 = cc[((k << 2) + 1) * cc_dim1 + 2] - cc[((k << 2) + 3) * cc_dim1 + 2];
	ti2 = cc[((k << 2) + 1) * cc_dim1 + 2] + cc[((k << 2) + 3) * cc_dim1 + 2];
	tr4 = cc[((k << 2) + 4) * cc_dim1 + 2] - cc[((k << 2) + 2) * cc_dim1 + 2];
	ti3 = cc[((k << 2) + 2) * cc_dim1 + 2] + cc[((k << 2) + 4) * cc_dim1 + 2];
	tr1 = cc[((k << 2) + 1) * cc_dim1 + 1] - cc[((k << 2) + 3) * cc_dim1 + 1];
	tr2 = cc[((k << 2) + 1) * cc_dim1 + 1] + cc[((k << 2) + 3) * cc_dim1 + 1];
	ti4 = cc[((k << 2) + 2) * cc_dim1 + 1] - cc[((k << 2) + 4) * cc_dim1 + 1];
	tr3 = cc[((k << 2) + 2) * cc_dim1 + 1] + cc[((k << 2) + 4) * cc_dim1 + 1];
	ch[(k + ch_dim2) * ch_dim1 + 1] = tr2 + tr3;
	ch[(k + ch_dim2 * 3) * ch_dim1 + 1] = tr2 - tr3;
	ch[(k + ch_dim2) * ch_dim1 + 2] = ti2 + ti3;
	ch[(k + ch_dim2 * 3) * ch_dim1 + 2] = ti2 - ti3;
	ch[(k + (ch_dim2 << 1)) * ch_dim1 + 1] = tr1 + tr4;
	ch[(k + (ch_dim2 << 2)) * ch_dim1 + 1] = tr1 - tr4;
	ch[(k + (ch_dim2 << 1)) * ch_dim1 + 2] = ti1 + ti4;
	ch[(k + (ch_dim2 << 2)) * ch_dim1 + 2] = ti1 - ti4;
/* L101: */
	}
	return;
L102:
	i__1 = *l1;
	for (k = 1; k <= i__1; ++k) {
	i__2 = *ido;
	for (i__ = 2; i__ <= i__2; i__ += 2) {
	    ti1 = cc[i__ + ((k << 2) + 1) * cc_dim1] - cc[i__ + ((k << 2) + 3)* cc_dim1];
	    ti2 = cc[i__ + ((k << 2) + 1) * cc_dim1] + cc[i__ + ((k << 2) + 3)* cc_dim1];
	    ti3 = cc[i__ + ((k << 2) + 2) * cc_dim1] + cc[i__ + ((k << 2) + 4)* cc_dim1];
	    tr4 = cc[i__ + ((k << 2) + 4) * cc_dim1] - cc[i__ + ((k << 2) + 2)* cc_dim1];
	    tr1 = cc[i__ - 1 + ((k << 2) + 1) * cc_dim1] - cc[i__ - 1 + ((k <<2) + 3) * cc_dim1];
	    tr2 = cc[i__ - 1 + ((k << 2) + 1) * cc_dim1] + cc[i__ - 1 + ((k <<2) + 3) * cc_dim1];
	    ti4 = cc[i__ - 1 + ((k << 2) + 2) * cc_dim1] - cc[i__ - 1 + ((k <<2) + 4) * cc_dim1];
	    tr3 = cc[i__ - 1 + ((k << 2) + 2) * cc_dim1] + cc[i__ - 1 + ((k <<2) + 4) * cc_dim1];
	    ch[i__ - 1 + (k + ch_dim2) * ch_dim1] = tr2 + tr3;
	    cr3 = tr2 - tr3;
	    ch[i__ + (k + ch_dim2) * ch_dim1] = ti2 + ti3;
	    ci3 = ti2 - ti3;
	    cr2 = tr1 + tr4;
	    cr4 = tr1 - tr4;
	    ci2 = ti1 + ti4;
	    ci4 = ti1 - ti4;
	    ch[i__ - 1 + (k + (ch_dim2 << 1)) * ch_dim1] = wa1[i__ - 1] * cr2 - wa1[i__] * ci2;
	    ch[i__ + (k + (ch_dim2 << 1)) * ch_dim1] = wa1[i__ - 1] * ci2 + wa1[i__] * cr2;
	    ch[i__ - 1 + (k + ch_dim2 * 3) * ch_dim1] = wa2[i__ - 1] * cr3 - wa2[i__] * ci3;
	    ch[i__ + (k + ch_dim2 * 3) * ch_dim1] = wa2[i__ - 1] * ci3 + wa2[i__] * cr3;
	    ch[i__ - 1 + (k + (ch_dim2 << 2)) * ch_dim1] = wa3[i__ - 1] * cr4 - wa3[i__] * ci4;
	    ch[i__ + (k + (ch_dim2 << 2)) * ch_dim1] = wa3[i__ - 1] * ci4 + wa3[i__] * cr4;
/* L103: */
	}
/* L104: */
	}
	return;
} /* passb4_ */

/* Subroutine */ static void s_passb5(int *ido, int *l1, real_t *cc, 
	real_t *ch, real_t *wa1, real_t *wa2, real_t *wa3, 
	real_t *wa4)
{
	/* Initialized data */

	static real_t tr11 = POST(0.3090169943749474241022934171828195886015458990288143106772431137);
	static real_t ti11 = POST(0.9510565162951535721164393337938214340569863412575022244730564442);
	static real_t tr12 = POST(-0.8090169943749474241022934171828190588601545899028814310677431135);
	static real_t ti12 = POST(0.5877852522924731291687059546390727685976524376431459107227248076);

	/* System generated locals */
	int cc_dim1, cc_offset, ch_dim1, ch_dim2, ch_offset, i__1, i__2;

	/* Local variables */
	int i__, k;
	real_t ci2, ci3, ci4, ci5, di3, di4, di5, di2, cr2, cr3, cr5, cr4, 
	    ti2, ti3, ti4, ti5, dr3, dr4, dr5, dr2, tr2, tr3, tr4, tr5;

	/* Parameter adjustments */
	ch_dim1 = *ido;
	ch_dim2 = *l1;
	ch_offset = 1 + ch_dim1 * (1 + ch_dim2);
	ch -= ch_offset;
	cc_dim1 = *ido;
	cc_offset = 1 + cc_dim1 * 6;
	cc -= cc_offset;
	--wa1;
	--wa2;
	--wa3;
	--wa4;

	/* Function Body */
	if (*ido != 2) {
	goto L102;
	}
	i__1 = *l1;
	for (k = 1; k <= i__1; ++k) {
	ti5 = cc[(k * 5 + 2) * cc_dim1 + 2] - cc[(k * 5 + 5) * cc_dim1 + 2];
	ti2 = cc[(k * 5 + 2) * cc_dim1 + 2] + cc[(k * 5 + 5) * cc_dim1 + 2];
	ti4 = cc[(k * 5 + 3) * cc_dim1 + 2] - cc[(k * 5 + 4) * cc_dim1 + 2];
	ti3 = cc[(k * 5 + 3) * cc_dim1 + 2] + cc[(k * 5 + 4) * cc_dim1 + 2];
	tr5 = cc[(k * 5 + 2) * cc_dim1 + 1] - cc[(k * 5 + 5) * cc_dim1 + 1];
	tr2 = cc[(k * 5 + 2) * cc_dim1 + 1] + cc[(k * 5 + 5) * cc_dim1 + 1];
	tr4 = cc[(k * 5 + 3) * cc_dim1 + 1] - cc[(k * 5 + 4) * cc_dim1 + 1];
	tr3 = cc[(k * 5 + 3) * cc_dim1 + 1] + cc[(k * 5 + 4) * cc_dim1 + 1];
	ch[(k + ch_dim2) * ch_dim1 + 1] = cc[(k * 5 + 1) * cc_dim1 + 1] + tr2 
		+ tr3;
	ch[(k + ch_dim2) * ch_dim1 + 2] = cc[(k * 5 + 1) * cc_dim1 + 2] + ti2 
		+ ti3;
	cr2 = cc[(k * 5 + 1) * cc_dim1 + 1] + tr11 * tr2 + tr12 * tr3;
	ci2 = cc[(k * 5 + 1) * cc_dim1 + 2] + tr11 * ti2 + tr12 * ti3;
	cr3 = cc[(k * 5 + 1) * cc_dim1 + 1] + tr12 * tr2 + tr11 * tr3;
	ci3 = cc[(k * 5 + 1) * cc_dim1 + 2] + tr12 * ti2 + tr11 * ti3;
	cr5 = ti11 * tr5 + ti12 * tr4;
	ci5 = ti11 * ti5 + ti12 * ti4;
	cr4 = ti12 * tr5 - ti11 * tr4;
	ci4 = ti12 * ti5 - ti11 * ti4;
	ch[(k + (ch_dim2 << 1)) * ch_dim1 + 1] = cr2 - ci5;
	ch[(k + ch_dim2 * 5) * ch_dim1 + 1] = cr2 + ci5;
	ch[(k + (ch_dim2 << 1)) * ch_dim1 + 2] = ci2 + cr5;
	ch[(k + ch_dim2 * 3) * ch_dim1 + 2] = ci3 + cr4;
	ch[(k + ch_dim2 * 3) * ch_dim1 + 1] = cr3 - ci4;
	ch[(k + (ch_dim2 << 2)) * ch_dim1 + 1] = cr3 + ci4;
	ch[(k + (ch_dim2 << 2)) * ch_dim1 + 2] = ci3 - cr4;
	ch[(k + ch_dim2 * 5) * ch_dim1 + 2] = ci2 - cr5;
/* L101: */
	}
	return;
L102:
	i__1 = *l1;
	for (k = 1; k <= i__1; ++k) {
	i__2 = *ido;
	for (i__ = 2; i__ <= i__2; i__ += 2) {
	    ti5 = cc[i__ + (k * 5 + 2) * cc_dim1] - cc[i__ + (k * 5 + 5) * cc_dim1];
	    ti2 = cc[i__ + (k * 5 + 2) * cc_dim1] + cc[i__ + (k * 5 + 5) * cc_dim1];
	    ti4 = cc[i__ + (k * 5 + 3) * cc_dim1] - cc[i__ + (k * 5 + 4) * cc_dim1];
	    ti3 = cc[i__ + (k * 5 + 3) * cc_dim1] + cc[i__ + (k * 5 + 4) * cc_dim1];
	    tr5 = cc[i__ - 1 + (k * 5 + 2) * cc_dim1] - cc[i__ - 1 + (k * 5 + 5) * cc_dim1];
	    tr2 = cc[i__ - 1 + (k * 5 + 2) * cc_dim1] + cc[i__ - 1 + (k * 5 + 5) * cc_dim1];
	    tr4 = cc[i__ - 1 + (k * 5 + 3) * cc_dim1] - cc[i__ - 1 + (k * 5 + 4) * cc_dim1];
	    tr3 = cc[i__ - 1 + (k * 5 + 3) * cc_dim1] + cc[i__ - 1 + (k * 5 + 4) * cc_dim1];
	    ch[i__ - 1 + (k + ch_dim2) * ch_dim1] = cc[i__ - 1 + (k * 5 + 1) *cc_dim1] + tr2 + tr3;
	    ch[i__ + (k + ch_dim2) * ch_dim1] = cc[i__ + (k * 5 + 1) * cc_dim1] + ti2 + ti3;
	    cr2 = cc[i__ - 1 + (k * 5 + 1) * cc_dim1] + tr11 * tr2 + tr12 * tr3;
	    ci2 = cc[i__ + (k * 5 + 1) * cc_dim1] + tr11 * ti2 + tr12 * ti3;
	    cr3 = cc[i__ - 1 + (k * 5 + 1) * cc_dim1] + tr12 * tr2 + tr11 * tr3;
	    ci3 = cc[i__ + (k * 5 + 1) * cc_dim1] + tr12 * ti2 + tr11 * ti3;
	    cr5 = ti11 * tr5 + ti12 * tr4;
	    ci5 = ti11 * ti5 + ti12 * ti4;
	    cr4 = ti12 * tr5 - ti11 * tr4;
	    ci4 = ti12 * ti5 - ti11 * ti4;
	    dr3 = cr3 - ci4;
	    dr4 = cr3 + ci4;
	    di3 = ci3 + cr4;
	    di4 = ci3 - cr4;
	    dr5 = cr2 + ci5;
	    dr2 = cr2 - ci5;
	    di5 = ci2 - cr5;
	    di2 = ci2 + cr5;
	    ch[i__ - 1 + (k + (ch_dim2 << 1)) * ch_dim1] = wa1[i__ - 1] * dr2 - wa1[i__] * di2;
	    ch[i__ + (k + (ch_dim2 << 1)) * ch_dim1] = wa1[i__ - 1] * di2 + wa1[i__] * dr2;
	    ch[i__ - 1 + (k + ch_dim2 * 3) * ch_dim1] = wa2[i__ - 1] * dr3 - wa2[i__] * di3;
	    ch[i__ + (k + ch_dim2 * 3) * ch_dim1] = wa2[i__ - 1] * di3 + wa2[i__] * dr3;
	    ch[i__ - 1 + (k + (ch_dim2 << 2)) * ch_dim1] = wa3[i__ - 1] * dr4 - wa3[i__] * di4;
	    ch[i__ + (k + (ch_dim2 << 2)) * ch_dim1] = wa3[i__ - 1] * di4 + wa3[i__] * dr4;
	    ch[i__ - 1 + (k + ch_dim2 * 5) * ch_dim1] = wa4[i__ - 1] * dr5 - wa4[i__] * di5;
	    ch[i__ + (k + ch_dim2 * 5) * ch_dim1] = wa4[i__ - 1] * di5 + wa4[i__] * dr5;
/* L103: */
	}
/* L104: */
	}
	return;
} /* passb5_ */

/* Subroutine */ static void s_passf(int *nac, int *ido, int *ip, int *
	l1, int *idl1, real_t *cc, real_t *c1, real_t *c2, 
	real_t *ch, real_t *ch2, real_t *wa)
{
	/* System generated locals */
	int ch_dim1, ch_dim2, ch_offset, cc_dim1, cc_dim2, cc_offset, c1_dim1,
		c1_dim2, c1_offset, c2_dim1, c2_offset, ch2_dim1, ch2_offset, 
		i__1, i__2, i__3;

	/* Local variables */
	int i__, j, k, l, jc, lc, ik, nt, idj, idl, inc, idp;
	real_t wai, war;
	int ipp2, idij, idlj, idot, ipph;

	/* Parameter adjustments */
	ch_dim1 = *ido;
	ch_dim2 = *l1;
	ch_offset = 1 + ch_dim1 * (1 + ch_dim2);
	ch -= ch_offset;
	c1_dim1 = *ido;
	c1_dim2 = *l1;
	c1_offset = 1 + c1_dim1 * (1 + c1_dim2);
	c1 -= c1_offset;
	cc_dim1 = *ido;
	cc_dim2 = *ip;
	cc_offset = 1 + cc_dim1 * (1 + cc_dim2);
	cc -= cc_offset;
	ch2_dim1 = *idl1;
	ch2_offset = 1 + ch2_dim1;
	ch2 -= ch2_offset;
	c2_dim1 = *idl1;
	c2_offset = 1 + c2_dim1;
	c2 -= c2_offset;
	--wa;

	/* Function Body */
	idot = *ido / 2;
	nt = *ip * *idl1;
	ipp2 = *ip + 2;
	ipph = (*ip + 1) / 2;
	idp = *ip * *ido;

	if (*ido < *l1) {
	goto L106;
	}
	i__1 = ipph;
	for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    i__3 = *ido;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		ch[i__ + (k + j * ch_dim2) * ch_dim1] = cc[i__ + (j + k * cc_dim2) * cc_dim1] + cc[i__ + (jc + k * cc_dim2) * cc_dim1];
		ch[i__ + (k + jc * ch_dim2) * ch_dim1] = cc[i__ + (j + k * cc_dim2) * cc_dim1] - cc[i__ + (jc + k * cc_dim2) * cc_dim1];
/* L101: */
	    }
/* L102: */
	}
/* L103: */
	}
	i__1 = *l1;
	for (k = 1; k <= i__1; ++k) {
	i__2 = *ido;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    ch[i__ + (k + ch_dim2) * ch_dim1] = cc[i__ + (k * cc_dim2 + 1) * cc_dim1];
/* L104: */
	}
/* L105: */
	}
	goto L112;
L106:
	i__1 = ipph;
	for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	i__2 = *ido;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = *l1;
	    for (k = 1; k <= i__3; ++k) {
		ch[i__ + (k + j * ch_dim2) * ch_dim1] = cc[i__ + (j + k * cc_dim2) * cc_dim1] + cc[i__ + (jc + k * cc_dim2) * cc_dim1];
		ch[i__ + (k + jc * ch_dim2) * ch_dim1] = cc[i__ + (j + k * cc_dim2) * cc_dim1] - cc[i__ + (jc + k * cc_dim2) * cc_dim1];
/* L107: */
	    }
/* L108: */
	}
/* L109: */
	}
	i__1 = *ido;
	for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    ch[i__ + (k + ch_dim2) * ch_dim1] = cc[i__ + (k * cc_dim2 + 1) * cc_dim1];
/* L110: */
	}
/* L111: */
	}
L112:
	idl = 2 - *ido;
	inc = 0;
	i__1 = ipph;
	for (l = 2; l <= i__1; ++l) {
	lc = ipp2 - l;
	idl += *ido;
	i__2 = *idl1;
	for (ik = 1; ik <= i__2; ++ik) {
	    c2[ik + l * c2_dim1] = ch2[ik + ch2_dim1] + wa[idl - 1] * ch2[ik + (ch2_dim1 << 1)];
	    c2[ik + lc * c2_dim1] = -wa[idl] * ch2[ik + *ip * ch2_dim1];
/* L113: */
	}
	idlj = idl;
	inc += *ido;
	i__2 = ipph;
	for (j = 3; j <= i__2; ++j) {
	    jc = ipp2 - j;
	    idlj += inc;
	    if (idlj > idp) {
		idlj -= idp;
	    }
	    war = wa[idlj - 1];
	    wai = wa[idlj];
	    i__3 = *idl1;
	    for (ik = 1; ik <= i__3; ++ik) {
		c2[ik + l * c2_dim1] += war * ch2[ik + j * ch2_dim1];
		c2[ik + lc * c2_dim1] -= wai * ch2[ik + jc * ch2_dim1];
/* L114: */
	    }
/* L115: */
	}
/* L116: */
	}
	i__1 = ipph;
	for (j = 2; j <= i__1; ++j) {
	i__2 = *idl1;
	for (ik = 1; ik <= i__2; ++ik) {
	    ch2[ik + ch2_dim1] += ch2[ik + j * ch2_dim1];
/* L117: */
	}
/* L118: */
	}
	i__1 = ipph;
	for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	i__2 = *idl1;
	for (ik = 2; ik <= i__2; ik += 2) {
	    ch2[ik - 1 + j * ch2_dim1] = c2[ik - 1 + j * c2_dim1] - c2[ik + jc * c2_dim1];
	    ch2[ik - 1 + jc * ch2_dim1] = c2[ik - 1 + j * c2_dim1] + c2[ik + jc * c2_dim1];
	    ch2[ik + j * ch2_dim1] = c2[ik + j * c2_dim1] + c2[ik - 1 + jc * c2_dim1];
	    ch2[ik + jc * ch2_dim1] = c2[ik + j * c2_dim1] - c2[ik - 1 + jc * c2_dim1];
/* L119: */
	}
/* L120: */
	}
	*nac = 1;
	if (*ido == 2) {
	return;
	}
	*nac = 0;
	i__1 = *idl1;
	for (ik = 1; ik <= i__1; ++ik) {
	c2[ik + c2_dim1] = ch2[ik + ch2_dim1];
/* L121: */
	}
	i__1 = *ip;
	for (j = 2; j <= i__1; ++j) {
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    c1[(k + j * c1_dim2) * c1_dim1 + 1] = ch[(k + j * ch_dim2) * ch_dim1 + 1];
	    c1[(k + j * c1_dim2) * c1_dim1 + 2] = ch[(k + j * ch_dim2) * ch_dim1 + 2];
/* L122: */
	}
/* L123: */
	}
	if (idot > *l1) {
	goto L127;
	}
	idij = 0;
	i__1 = *ip;
	for (j = 2; j <= i__1; ++j) {
	idij += 2;
	i__2 = *ido;
	for (i__ = 4; i__ <= i__2; i__ += 2) {
	    idij += 2;
	    i__3 = *l1;
	    for (k = 1; k <= i__3; ++k) {
		c1[	i__ - 1 + (k + j * c1_dim2) * c1_dim1] = wa[idij - 1] * ch[
			i__ - 1 + (k + j * ch_dim2) * ch_dim1] + wa[idij] * 
			ch[i__ + (k + j * ch_dim2) * ch_dim1];
		c1[	i__ + (k + j * c1_dim2) * c1_dim1] = wa[idij - 1] * ch[i__ 
			+ (k + j * ch_dim2) * ch_dim1] - wa[idij] * ch[i__ - 
			1 + (k + j * ch_dim2) * ch_dim1];
/* L124: */
	    }
/* L125: */
	}
/* L126: */
	}
	return;
L127:
	idj = 2 - *ido;
	i__1 = *ip;
	for (j = 2; j <= i__1; ++j) {
	idj += *ido;
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    idij = idj;
	    i__3 = *ido;
	    for (i__ = 4; i__ <= i__3; i__ += 2) {
		idij += 2;
		c1[i__ - 1 + (k + j * c1_dim2) * c1_dim1] = wa[idij - 1] * ch[
			i__ - 1 + (k + j * ch_dim2) * ch_dim1] + wa[idij] * 
			ch[i__ + (k + j * ch_dim2) * ch_dim1];
		c1[i__ + (k + j * c1_dim2) * c1_dim1] = wa[idij - 1] * ch[i__ 
			+ (k + j * ch_dim2) * ch_dim1] - wa[idij] * ch[i__ - 
			1 + (k + j * ch_dim2) * ch_dim1];
/* L128: */
	    }
/* L129: */
	}
/* L130: */
	}
	return;
} /* passf_ */

/* Subroutine */ static void s_passf2(int *ido, int *l1, real_t *cc, 
	real_t *ch, real_t *wa1)
{
	/* System generated locals */
	int cc_dim1, cc_offset, ch_dim1, ch_dim2, ch_offset, i__1, i__2;

	/* Local variables */
	int i__, k;
	real_t ti2, tr2;

	/* Parameter adjustments */
	ch_dim1 = *ido;
	ch_dim2 = *l1;
	ch_offset = 1 + ch_dim1 * (1 + ch_dim2);
	ch -= ch_offset;
	cc_dim1 = *ido;
	cc_offset = 1 + cc_dim1 * 3;
	cc -= cc_offset;
	--wa1;

	/* Function Body */
	if (*ido > 2) {
	goto L102;
	}
	i__1 = *l1;
	for (k = 1; k <= i__1; ++k) {
	ch[(k + ch_dim2) * ch_dim1 + 1] = cc[((k << 1) + 1) * cc_dim1 + 1] + cc[((k << 1) + 2) * cc_dim1 + 1];
	ch[(k + (ch_dim2 << 1)) * ch_dim1 + 1] = cc[((k << 1) + 1) * cc_dim1 + 1] - cc[((k << 1) + 2) * cc_dim1 + 1];
	ch[(k + ch_dim2) * ch_dim1 + 2] = cc[((k << 1) + 1) * cc_dim1 + 2] + cc[((k << 1) + 2) * cc_dim1 + 2];
	ch[(k + (ch_dim2 << 1)) * ch_dim1 + 2] = cc[((k << 1) + 1) * cc_dim1 + 2] - cc[((k << 1) + 2) * cc_dim1 + 2];
/* L101: */
	}
	return;
L102:
	i__1 = *l1;
	for (k = 1; k <= i__1; ++k) {
	i__2 = *ido;
	for (i__ = 2; i__ <= i__2; i__ += 2) {
	    ch[i__ - 1 + (k + ch_dim2) * ch_dim1] = cc[i__ - 1 + ((k << 1) + 1) * cc_dim1] + cc[i__ - 1 + ((k << 1) + 2) * cc_dim1];
	    tr2 = cc[i__ - 1 + ((k << 1) + 1) * cc_dim1] - cc[i__ - 1 + ((k <<1) + 2) * cc_dim1];
	    ch[i__ + (k + ch_dim2) * ch_dim1] = cc[i__ + ((k << 1) + 1) * cc_dim1] + cc[i__ + ((k << 1) + 2) * cc_dim1];
	    ti2 = cc[i__ + ((k << 1) + 1) * cc_dim1] - cc[i__ + ((k << 1) + 2)* cc_dim1];
	    ch[i__ + (k + (ch_dim2 << 1)) * ch_dim1] = wa1[i__ - 1] * ti2 - wa1[i__] * tr2;
	    ch[i__ - 1 + (k + (ch_dim2 << 1)) * ch_dim1] = wa1[i__ - 1] * tr2 + wa1[i__] * ti2;
/* L103: */
	}
/* L104: */
	}
	return;
} /* passf2_ */

/* Subroutine */ static void s_passf3(int *ido, int *l1, real_t *cc, 
	real_t *ch, real_t *wa1, real_t *wa2)
{
	/* Initialized data */

	static real_t taur = POST(-0.5);
	static real_t taui = POST(-0.8660254037844386467637231707529361834740262690519031402790348975);

	/* System generated locals */
	int cc_dim1, cc_offset, ch_dim1, ch_dim2, ch_offset, i__1, i__2;

	/* Local variables */
	int i__, k;
	real_t ci2, ci3, di2, di3, cr2, cr3, dr2, dr3, ti2, tr2;

	/* Parameter adjustments */
	ch_dim1 = *ido;
	ch_dim2 = *l1;
	ch_offset = 1 + ch_dim1 * (1 + ch_dim2);
	ch -= ch_offset;
	cc_dim1 = *ido;
	cc_offset = 1 + (cc_dim1 << 2);
	cc -= cc_offset;
	--wa1;
	--wa2;

	/* Function Body */
	if (*ido != 2) {
	goto L102;
	}
	i__1 = *l1;
	for (k = 1; k <= i__1; ++k) {
	tr2 = cc[(k * 3 + 2) * cc_dim1 + 1] + cc[(k * 3 + 3) * cc_dim1 + 1];
	cr2 = cc[(k * 3 + 1) * cc_dim1 + 1] + taur * tr2;
	ch[(k + ch_dim2) * ch_dim1 + 1] = cc[(k * 3 + 1) * cc_dim1 + 1] + tr2;
	ti2 = cc[(k * 3 + 2) * cc_dim1 + 2] + cc[(k * 3 + 3) * cc_dim1 + 2];
	ci2 = cc[(k * 3 + 1) * cc_dim1 + 2] + taur * ti2;
	ch[(k + ch_dim2) * ch_dim1 + 2] = cc[(k * 3 + 1) * cc_dim1 + 2] + ti2;
	cr3 = taui * (cc[(k * 3 + 2) * cc_dim1 + 1] - cc[(k * 3 + 3) * cc_dim1 + 1]);
	ci3 = taui * (cc[(k * 3 + 2) * cc_dim1 + 2] - cc[(k * 3 + 3) * cc_dim1 + 2]);
	ch[(k + (ch_dim2 << 1)) * ch_dim1 + 1] = cr2 - ci3;
	ch[(k + ch_dim2 * 3) * ch_dim1 + 1] = cr2 + ci3;
	ch[(k + (ch_dim2 << 1)) * ch_dim1 + 2] = ci2 + cr3;
	ch[(k + ch_dim2 * 3) * ch_dim1 + 2] = ci2 - cr3;
/* L101: */
	}
	return;
L102:
	i__1 = *l1;
	for (k = 1; k <= i__1; ++k) {
	i__2 = *ido;
	for (i__ = 2; i__ <= i__2; i__ += 2) {
	    tr2 = cc[i__ - 1 + (k * 3 + 2) * cc_dim1] + cc[i__ - 1 + (k * 3 + 3) * cc_dim1];
	    cr2 = cc[i__ - 1 + (k * 3 + 1) * cc_dim1] + taur * tr2;
	    ch[i__ - 1 + (k + ch_dim2) * ch_dim1] = cc[i__ - 1 + (k * 3 + 1) *cc_dim1] + tr2;
	    ti2 = cc[i__ + (k * 3 + 2) * cc_dim1] + cc[i__ + (k * 3 + 3) * cc_dim1];
	    ci2 = cc[i__ + (k * 3 + 1) * cc_dim1] + taur * ti2;
	    ch[i__ + (k + ch_dim2) * ch_dim1] = cc[i__ + (k * 3 + 1) * cc_dim1] + ti2;
	    cr3 = taui * (cc[i__ - 1 + (k * 3 + 2) * cc_dim1] - cc[i__ - 1 + (k * 3 + 3) * cc_dim1]);
	    ci3 = taui * (cc[i__ + (k * 3 + 2) * cc_dim1] - cc[i__ + (k * 3 + 3) * cc_dim1]);
	    dr2 = cr2 - ci3;
	    dr3 = cr2 + ci3;
	    di2 = ci2 + cr3;
	    di3 = ci2 - cr3;
	    ch[i__ + (k + (ch_dim2 << 1)) * ch_dim1] = wa1[i__ - 1] * di2 - wa1[i__] * dr2;
	    ch[i__ - 1 + (k + (ch_dim2 << 1)) * ch_dim1] = wa1[i__ - 1] * dr2 + wa1[i__] * di2;
	    ch[i__ + (k + ch_dim2 * 3) * ch_dim1] = wa2[i__ - 1] * di3 - wa2[i__] * dr3;
	    ch[i__ - 1 + (k + ch_dim2 * 3) * ch_dim1] = wa2[i__ - 1] * dr3 + wa2[i__] * di3;
/* L103: */
	}
/* L104: */
	}
	return;
} /* passf3_ */

/* Subroutine */ static void s_passf4(int *ido, int *l1, real_t *cc, 
	real_t *ch, real_t *wa1, real_t *wa2, real_t *wa3)
{
	/* System generated locals */
	int cc_dim1, cc_offset, ch_dim1, ch_dim2, ch_offset, i__1, i__2;

	/* Local variables */
	int i__, k;
	real_t ci2, ci3, ci4, cr2, cr3, cr4, ti1, ti2, ti3, ti4, tr1, tr2, 
	    tr3, tr4;

	/* Parameter adjustments */
	ch_dim1 = *ido;
	ch_dim2 = *l1;
	ch_offset = 1 + ch_dim1 * (1 + ch_dim2);
	ch -= ch_offset;
	cc_dim1 = *ido;
	cc_offset = 1 + cc_dim1 * 5;
	cc -= cc_offset;
	--wa1;
	--wa2;
	--wa3;

	/* Function Body */
	if (*ido != 2) {
	goto L102;
	}
	i__1 = *l1;
	for (k = 1; k <= i__1; ++k) {
	ti1 = cc[((k << 2) + 1) * cc_dim1 + 2] - cc[((k << 2) + 3) * cc_dim1 + 2];
	ti2 = cc[((k << 2) + 1) * cc_dim1 + 2] + cc[((k << 2) + 3) * cc_dim1 + 2];
	tr4 = cc[((k << 2) + 2) * cc_dim1 + 2] - cc[((k << 2) + 4) * cc_dim1 + 2];
	ti3 = cc[((k << 2) + 2) * cc_dim1 + 2] + cc[((k << 2) + 4) * cc_dim1 + 2];
	tr1 = cc[((k << 2) + 1) * cc_dim1 + 1] - cc[((k << 2) + 3) * cc_dim1 + 1];
	tr2 = cc[((k << 2) + 1) * cc_dim1 + 1] + cc[((k << 2) + 3) * cc_dim1 + 1];
	ti4 = cc[((k << 2) + 4) * cc_dim1 + 1] - cc[((k << 2) + 2) * cc_dim1 + 1];
	tr3 = cc[((k << 2) + 2) * cc_dim1 + 1] + cc[((k << 2) + 4) * cc_dim1 + 1];
	ch[(k + ch_dim2) * ch_dim1 + 1] = tr2 + tr3;
	ch[(k + ch_dim2 * 3) * ch_dim1 + 1] = tr2 - tr3;
	ch[(k + ch_dim2) * ch_dim1 + 2] = ti2 + ti3;
	ch[(k + ch_dim2 * 3) * ch_dim1 + 2] = ti2 - ti3;
	ch[(k + (ch_dim2 << 1)) * ch_dim1 + 1] = tr1 + tr4;
	ch[(k + (ch_dim2 << 2)) * ch_dim1 + 1] = tr1 - tr4;
	ch[(k + (ch_dim2 << 1)) * ch_dim1 + 2] = ti1 + ti4;
	ch[(k + (ch_dim2 << 2)) * ch_dim1 + 2] = ti1 - ti4;
/* L101: */
	}
	return;
L102:
	i__1 = *l1;
	for (k = 1; k <= i__1; ++k) {
	i__2 = *ido;
	for (i__ = 2; i__ <= i__2; i__ += 2) {
	    ti1 = cc[i__ + ((k << 2) + 1) * cc_dim1] - cc[i__ + ((k << 2) + 3)* cc_dim1];
	    ti2 = cc[i__ + ((k << 2) + 1) * cc_dim1] + cc[i__ + ((k << 2) + 3)* cc_dim1];
	    ti3 = cc[i__ + ((k << 2) + 2) * cc_dim1] + cc[i__ + ((k << 2) + 4)* cc_dim1];
	    tr4 = cc[i__ + ((k << 2) + 2) * cc_dim1] - cc[i__ + ((k << 2) + 4)* cc_dim1];
	    tr1 = cc[i__ - 1 + ((k << 2) + 1) * cc_dim1] - cc[i__ - 1 + ((k <<2) + 3) * cc_dim1];
	    tr2 = cc[i__ - 1 + ((k << 2) + 1) * cc_dim1] + cc[i__ - 1 + ((k <<2) + 3) * cc_dim1];
	    ti4 = cc[i__ - 1 + ((k << 2) + 4) * cc_dim1] - cc[i__ - 1 + ((k <<2) + 2) * cc_dim1];
	    tr3 = cc[i__ - 1 + ((k << 2) + 2) * cc_dim1] + cc[i__ - 1 + ((k <<2) + 4) * cc_dim1];
	    ch[i__ - 1 + (k + ch_dim2) * ch_dim1] = tr2 + tr3;
	    cr3 = tr2 - tr3;
	    ch[i__ + (k + ch_dim2) * ch_dim1] = ti2 + ti3;
	    ci3 = ti2 - ti3;
	    cr2 = tr1 + tr4;
	    cr4 = tr1 - tr4;
	    ci2 = ti1 + ti4;
	    ci4 = ti1 - ti4;
	    ch[i__ - 1 + (k + (ch_dim2 << 1)) * ch_dim1] = wa1[i__ - 1] * cr2 + wa1[i__] * ci2;
	    ch[i__ + (k + (ch_dim2 << 1)) * ch_dim1] = wa1[i__ - 1] * ci2 - wa1[i__] * cr2;
	    ch[i__ - 1 + (k + ch_dim2 * 3) * ch_dim1] = wa2[i__ - 1] * cr3 + wa2[i__] * ci3;
	    ch[i__ + (k + ch_dim2 * 3) * ch_dim1] = wa2[i__ - 1] * ci3 - wa2[i__] * cr3;
	    ch[i__ - 1 + (k + (ch_dim2 << 2)) * ch_dim1] = wa3[i__ - 1] * cr4 + wa3[i__] * ci4;
	    ch[i__ + (k + (ch_dim2 << 2)) * ch_dim1] = wa3[i__ - 1] * ci4 - wa3[i__] * cr4;
/* L103: */
	}
/* L104: */
	}
	return;
} /* passf4_ */

/* Subroutine */ static void s_passf5(int *ido, int *l1, real_t *cc, 
	real_t *ch, real_t *wa1, real_t *wa2, real_t *wa3, 
	real_t *wa4)
{
	/* Initialized data */

	static real_t tr11 = POST(0.3090169943749474241022934171828195886015458990288143106772431137);
	static real_t ti11 = POST(-0.9510565162951535721164393337938214340569863412575022244730564442);
	static real_t tr12 = POST(-0.8090169943749474241022934171828190588601545899028814310677431135);
	static real_t ti12 = POST(-0.5877852522924731291687059546390727685976524376431459107227248076);

	/* System generated locals */
	int cc_dim1, cc_offset, ch_dim1, ch_dim2, ch_offset, i__1, i__2;

	/* Local variables */
	int i__, k;
	real_t ci2, ci3, ci4, ci5, di3, di4, di5, di2, cr2, cr3, cr5, cr4, 
		ti2, ti3, ti4, ti5, dr3, dr4, dr5, dr2, tr2, tr3, tr4, tr5;

	/* Parameter adjustments */
	ch_dim1 = *ido;
	ch_dim2 = *l1;
	ch_offset = 1 + ch_dim1 * (1 + ch_dim2);
	ch -= ch_offset;
	cc_dim1 = *ido;
	cc_offset = 1 + cc_dim1 * 6;
	cc -= cc_offset;
	--wa1;
	--wa2;
	--wa3;
	--wa4;

	/* Function Body */
	if (*ido != 2) {
	goto L102;
	}
	i__1 = *l1;
	for (k = 1; k <= i__1; ++k) {
	ti5 = cc[(k * 5 + 2) * cc_dim1 + 2] - cc[(k * 5 + 5) * cc_dim1 + 2];
	ti2 = cc[(k * 5 + 2) * cc_dim1 + 2] + cc[(k * 5 + 5) * cc_dim1 + 2];
	ti4 = cc[(k * 5 + 3) * cc_dim1 + 2] - cc[(k * 5 + 4) * cc_dim1 + 2];
	ti3 = cc[(k * 5 + 3) * cc_dim1 + 2] + cc[(k * 5 + 4) * cc_dim1 + 2];
	tr5 = cc[(k * 5 + 2) * cc_dim1 + 1] - cc[(k * 5 + 5) * cc_dim1 + 1];
	tr2 = cc[(k * 5 + 2) * cc_dim1 + 1] + cc[(k * 5 + 5) * cc_dim1 + 1];
	tr4 = cc[(k * 5 + 3) * cc_dim1 + 1] - cc[(k * 5 + 4) * cc_dim1 + 1];
	tr3 = cc[(k * 5 + 3) * cc_dim1 + 1] + cc[(k * 5 + 4) * cc_dim1 + 1];
	ch[(k + ch_dim2) * ch_dim1 + 1] = cc[(k * 5 + 1) * cc_dim1 + 1] + tr2 + tr3;
	ch[(k + ch_dim2) * ch_dim1 + 2] = cc[(k * 5 + 1) * cc_dim1 + 2] + ti2 + ti3;
	cr2 = cc[(k * 5 + 1) * cc_dim1 + 1] + tr11 * tr2 + tr12 * tr3;
	ci2 = cc[(k * 5 + 1) * cc_dim1 + 2] + tr11 * ti2 + tr12 * ti3;
	cr3 = cc[(k * 5 + 1) * cc_dim1 + 1] + tr12 * tr2 + tr11 * tr3;
	ci3 = cc[(k * 5 + 1) * cc_dim1 + 2] + tr12 * ti2 + tr11 * ti3;
	cr5 = ti11 * tr5 + ti12 * tr4;
	ci5 = ti11 * ti5 + ti12 * ti4;
	cr4 = ti12 * tr5 - ti11 * tr4;
	ci4 = ti12 * ti5 - ti11 * ti4;
	ch[(k + (ch_dim2 << 1)) * ch_dim1 + 1] = cr2 - ci5;
	ch[(k + ch_dim2 * 5) * ch_dim1 + 1] = cr2 + ci5;
	ch[(k + (ch_dim2 << 1)) * ch_dim1 + 2] = ci2 + cr5;
	ch[(k + ch_dim2 * 3) * ch_dim1 + 2] = ci3 + cr4;
	ch[(k + ch_dim2 * 3) * ch_dim1 + 1] = cr3 - ci4;
	ch[(k + (ch_dim2 << 2)) * ch_dim1 + 1] = cr3 + ci4;
	ch[(k + (ch_dim2 << 2)) * ch_dim1 + 2] = ci3 - cr4;
	ch[(k + ch_dim2 * 5) * ch_dim1 + 2] = ci2 - cr5;
/* L101: */
	}
	return;
L102:
	i__1 = *l1;
	for (k = 1; k <= i__1; ++k) {
	i__2 = *ido;
	for (i__ = 2; i__ <= i__2; i__ += 2) {
	    ti5 = cc[i__ + (k * 5 + 2) * cc_dim1] - cc[i__ + (k * 5 + 5) * cc_dim1];
	    ti2 = cc[i__ + (k * 5 + 2) * cc_dim1] + cc[i__ + (k * 5 + 5) * cc_dim1];
	    ti4 = cc[i__ + (k * 5 + 3) * cc_dim1] - cc[i__ + (k * 5 + 4) * cc_dim1];
	    ti3 = cc[i__ + (k * 5 + 3) * cc_dim1] + cc[i__ + (k * 5 + 4) * cc_dim1];
	    tr5 = cc[i__ - 1 + (k * 5 + 2) * cc_dim1] - cc[i__ - 1 + (k * 5 + 5) * cc_dim1];
	    tr2 = cc[i__ - 1 + (k * 5 + 2) * cc_dim1] + cc[i__ - 1 + (k * 5 + 5) * cc_dim1];
	    tr4 = cc[i__ - 1 + (k * 5 + 3) * cc_dim1] - cc[i__ - 1 + (k * 5 + 4) * cc_dim1];
	    tr3 = cc[i__ - 1 + (k * 5 + 3) * cc_dim1] + cc[i__ - 1 + (k * 5 + 4) * cc_dim1];
	    ch[i__ - 1 + (k + ch_dim2) * ch_dim1] = cc[i__ - 1 + (k * 5 + 1) *cc_dim1] + tr2 + tr3;
	    ch[i__ + (k + ch_dim2) * ch_dim1] = cc[i__ + (k * 5 + 1) * cc_dim1] + ti2 + ti3;
	    cr2 = cc[i__ - 1 + (k * 5 + 1) * cc_dim1] + tr11 * tr2 + tr12 * tr3;
	    ci2 = cc[i__ + (k * 5 + 1) * cc_dim1] + tr11 * ti2 + tr12 * ti3;
	    cr3 = cc[i__ - 1 + (k * 5 + 1) * cc_dim1] + tr12 * tr2 + tr11 * tr3;
	    ci3 = cc[i__ + (k * 5 + 1) * cc_dim1] + tr12 * ti2 + tr11 * ti3;
	    cr5 = ti11 * tr5 + ti12 * tr4;
	    ci5 = ti11 * ti5 + ti12 * ti4;
	    cr4 = ti12 * tr5 - ti11 * tr4;
	    ci4 = ti12 * ti5 - ti11 * ti4;
	    dr3 = cr3 - ci4;
	    dr4 = cr3 + ci4;
	    di3 = ci3 + cr4;
	    di4 = ci3 - cr4;
	    dr5 = cr2 + ci5;
	    dr2 = cr2 - ci5;
	    di5 = ci2 - cr5;
	    di2 = ci2 + cr5;
	    ch[i__ - 1 + (k + (ch_dim2 << 1)) * ch_dim1] = wa1[i__ - 1] * dr2 + wa1[i__] * di2;
	    ch[i__ + (k + (ch_dim2 << 1)) * ch_dim1] = wa1[i__ - 1] * di2 - wa1[i__] * dr2;
	    ch[i__ - 1 + (k + ch_dim2 * 3) * ch_dim1] = wa2[i__ - 1] * dr3 + wa2[i__] * di3;
	    ch[i__ + (k + ch_dim2 * 3) * ch_dim1] = wa2[i__ - 1] * di3 - wa2[i__] * dr3;
	    ch[i__ - 1 + (k + (ch_dim2 << 2)) * ch_dim1] = wa3[i__ - 1] * dr4 + wa3[i__] * di4;
	    ch[i__ + (k + (ch_dim2 << 2)) * ch_dim1] = wa3[i__ - 1] * di4 - wa3[i__] * dr4;
	    ch[i__ - 1 + (k + ch_dim2 * 5) * ch_dim1] = wa4[i__ - 1] * dr5 + wa4[i__] * di5;
	    ch[i__ + (k + ch_dim2 * 5) * ch_dim1] = wa4[i__ - 1] * di5 - wa4[i__] * dr5;
/* L103: */
	}
/* L104: */
	}
	return;
} /* passf5_ */

/* Subroutine */ static void s_radb2(int *ido, int *l1, real_t *cc, 
	real_t *ch, real_t *wa1)
{
	/* System generated locals */
	int cc_dim1, cc_offset, ch_dim1, ch_dim2, ch_offset, i__1, i__2;

	/* Local variables */
	int i__, k, ic;
	real_t ti2, tr2;
	int idp2;

	/* Parameter adjustments */
	ch_dim1 = *ido;
	ch_dim2 = *l1;
	ch_offset = 1 + ch_dim1 * (1 + ch_dim2);
	ch -= ch_offset;
	cc_dim1 = *ido;
	cc_offset = 1 + cc_dim1 * 3;
	cc -= cc_offset;
	--wa1;

	/* Function Body */
	i__1 = *l1;
	for (k = 1; k <= i__1; ++k) {
	ch[(k + ch_dim2) * ch_dim1 + 1] = cc[((k << 1) + 1) * cc_dim1 + 1] + cc[*ido + ((k << 1) + 2) * cc_dim1];
	ch[(k + (ch_dim2 << 1)) * ch_dim1 + 1] = cc[((k << 1) + 1) * cc_dim1 + 1] - cc[*ido + ((k << 1) + 2) * cc_dim1];
/* L101: */
	}
	if ((i__1 = *ido - 2) < 0) {
	goto L107;
	} else if (i__1 == 0) {
	goto L105;
	} else {
	goto L102;
	}
L102:
	idp2 = *ido + 2;
	i__1 = *l1;
	for (k = 1; k <= i__1; ++k) {
	i__2 = *ido;
	for (i__ = 3; i__ <= i__2; i__ += 2) {
	    ic = idp2 - i__;
	    ch[i__ - 1 + (k + ch_dim2) * ch_dim1] = cc[i__ - 1 + ((k << 1) + 1) * cc_dim1] + cc[ic - 1 + ((k << 1) + 2) * cc_dim1];
	    tr2 = cc[i__ - 1 + ((k << 1) + 1) * cc_dim1] - cc[ic - 1 + ((k << 1) + 2) * cc_dim1];
	    ch[i__ + (k + ch_dim2) * ch_dim1] = cc[i__ + ((k << 1) + 1) * cc_dim1] - cc[ic + ((k << 1) + 2) * cc_dim1];
	    ti2 = cc[i__ + ((k << 1) + 1) * cc_dim1] + cc[ic + ((k << 1) + 2) * cc_dim1];
	    ch[i__ - 1 + (k + (ch_dim2 << 1)) * ch_dim1] = wa1[i__ - 2] * tr2 - wa1[i__ - 1] * ti2;
	    ch[i__ + (k + (ch_dim2 << 1)) * ch_dim1] = wa1[i__ - 2] * ti2 + wa1[i__ - 1] * tr2;
/* L103: */
	}
/* L104: */
	}
	if (*ido % 2 == 1) {
	return;
	}
L105:
	i__1 = *l1;
	for (k = 1; k <= i__1; ++k) {
	ch[*ido + (k + ch_dim2) * ch_dim1] = cc[*ido + ((k << 1) + 1) * cc_dim1] + cc[*ido + ((k << 1) + 1) * cc_dim1];
	ch[*ido + (k + (ch_dim2 << 1)) * ch_dim1] = -(cc[((k << 1) + 2) * cc_dim1 + 1] + cc[((k << 1) + 2) * cc_dim1 + 1]);
/* L106: */
	}
L107:
	return;
} /* radb2_ */

/* Subroutine */ static void s_radb3(int *ido, int *l1, real_t *cc, 
	real_t *ch, real_t *wa1, real_t *wa2)
{
	/* Initialized data */

	static real_t taur = POST(-0.5);
	static real_t taui = POST(0.8660254037844386467637231707529361834710262690519031402790348975);

	/* System generated locals */
	int cc_dim1, cc_offset, ch_dim1, ch_dim2, ch_offset, i__1, i__2;

	/* Local variables */
	int i__, k, ic;
	real_t ci2, ci3, di2, di3, cr2, cr3, dr2, dr3, ti2, tr2;
	int idp2;

	/* Parameter adjustments */
	ch_dim1 = *ido;
	ch_dim2 = *l1;
	ch_offset = 1 + ch_dim1 * (1 + ch_dim2);
	ch -= ch_offset;
	cc_dim1 = *ido;
	cc_offset = 1 + (cc_dim1 << 2);
	cc -= cc_offset;
	--wa1;
	--wa2;

	/* Function Body */
	i__1 = *l1;
	for (k = 1; k <= i__1; ++k) {
	tr2 = cc[*ido + (k * 3 + 2) * cc_dim1] + cc[*ido + (k * 3 + 2) * cc_dim1];
	cr2 = cc[(k * 3 + 1) * cc_dim1 + 1] + taur * tr2;
	ch[(k + ch_dim2) * ch_dim1 + 1] = cc[(k * 3 + 1) * cc_dim1 + 1] + tr2;
	ci3 = taui * (cc[(k * 3 + 3) * cc_dim1 + 1] + cc[(k * 3 + 3) * cc_dim1 + 1]);
	ch[(k + (ch_dim2 << 1)) * ch_dim1 + 1] = cr2 - ci3;
	ch[(k + ch_dim2 * 3) * ch_dim1 + 1] = cr2 + ci3;
/* L101: */
	}
	if (*ido == 1) {
	return;
	}
	idp2 = *ido + 2;
	i__1 = *l1;
	for (k = 1; k <= i__1; ++k) {
	i__2 = *ido;
	for (i__ = 3; i__ <= i__2; i__ += 2) {
	    ic = idp2 - i__;
	    tr2 = cc[i__ - 1 + (k * 3 + 3) * cc_dim1] + cc[ic - 1 + (k * 3 + 2) * cc_dim1];
	    cr2 = cc[i__ - 1 + (k * 3 + 1) * cc_dim1] + taur * tr2;
	    ch[i__ - 1 + (k + ch_dim2) * ch_dim1] = cc[i__ - 1 + (k * 3 + 1) *cc_dim1] + tr2;
	    ti2 = cc[i__ + (k * 3 + 3) * cc_dim1] - cc[ic + (k * 3 + 2) * cc_dim1];
	    ci2 = cc[i__ + (k * 3 + 1) * cc_dim1] + taur * ti2;
	    ch[i__ + (k + ch_dim2) * ch_dim1] = cc[i__ + (k * 3 + 1) * cc_dim1] + ti2;
	    cr3 = taui * (cc[i__ - 1 + (k * 3 + 3) * cc_dim1] - cc[ic - 1 + (k * 3 + 2) * cc_dim1]);
	    ci3 = taui * (cc[i__ + (k * 3 + 3) * cc_dim1] + cc[ic + (k * 3 + 2) * cc_dim1]);
	    dr2 = cr2 - ci3;
	    dr3 = cr2 + ci3;
	    di2 = ci2 + cr3;
	    di3 = ci2 - cr3;
	    ch[i__ - 1 + (k + (ch_dim2 << 1)) * ch_dim1] = wa1[i__ - 2] * dr2 - wa1[i__ - 1] * di2;
	    ch[i__ + (k + (ch_dim2 << 1)) * ch_dim1] = wa1[i__ - 2] * di2 + wa1[i__ - 1] * dr2;
	    ch[i__ - 1 + (k + ch_dim2 * 3) * ch_dim1] = wa2[i__ - 2] * dr3 - wa2[i__ - 1] * di3;
	    ch[i__ + (k + ch_dim2 * 3) * ch_dim1] = wa2[i__ - 2] * di3 + wa2[i__ - 1] * dr3;
/* L102: */
	}
/* L103: */
	}
	return;
} /* radb3_ */

/* Subroutine */ static void s_radb4(int *ido, int *l1, real_t *cc, 
	real_t *ch, real_t *wa1, real_t *wa2, real_t *wa3)
{
	/* Initialized data */

	static real_t sqrt2 = POST(1.41421356237309504880168872420969807856967187536948073176679738);

	/* System generated locals */
	int cc_dim1, cc_offset, ch_dim1, ch_dim2, ch_offset, i__1, i__2;

	/* Local variables */
	int i__, k, ic;
	real_t ci2, ci3, ci4, cr2, cr3, cr4, ti1, ti2, ti3, ti4, tr1, tr2, 
	    tr3, tr4;
	int idp2;

	/* Parameter adjustments */
	ch_dim1 = *ido;
	ch_dim2 = *l1;
	ch_offset = 1 + ch_dim1 * (1 + ch_dim2);
	ch -= ch_offset;
	cc_dim1 = *ido;
	cc_offset = 1 + cc_dim1 * 5;
	cc -= cc_offset;
	--wa1;
	--wa2;
	--wa3;

	/* Function Body */
	i__1 = *l1;
	for (k = 1; k <= i__1; ++k) {
	tr1 = cc[((k << 2) + 1) * cc_dim1 + 1] - cc[*ido + ((k << 2) + 4) * cc_dim1];
	tr2 = cc[((k << 2) + 1) * cc_dim1 + 1] + cc[*ido + ((k << 2) + 4) * cc_dim1];
	tr3 = cc[*ido + ((k << 2) + 2) * cc_dim1] + cc[*ido + ((k << 2) + 2) *cc_dim1];
	tr4 = cc[((k << 2) + 3) * cc_dim1 + 1] + cc[((k << 2) + 3) * cc_dim1 + 1];
	ch[(k + ch_dim2) * ch_dim1 + 1] = tr2 + tr3;
	ch[(k + (ch_dim2 << 1)) * ch_dim1 + 1] = tr1 - tr4;
	ch[(k + ch_dim2 * 3) * ch_dim1 + 1] = tr2 - tr3;
	ch[(k + (ch_dim2 << 2)) * ch_dim1 + 1] = tr1 + tr4;
/* L101: */
	}
	if ((i__1 = *ido - 2) < 0) {
	goto L107;
	} else if (i__1 == 0) {
	goto L105;
	} else {
	goto L102;
	}
L102:
	idp2 = *ido + 2;
	i__1 = *l1;
	for (k = 1; k <= i__1; ++k) {
	i__2 = *ido;
	for (i__ = 3; i__ <= i__2; i__ += 2) {
	    ic = idp2 - i__;
	    ti1 = cc[i__ + ((k << 2) + 1) * cc_dim1] + cc[ic + ((k << 2) + 4) * cc_dim1];
	    ti2 = cc[i__ + ((k << 2) + 1) * cc_dim1] - cc[ic + ((k << 2) + 4) * cc_dim1];
	    ti3 = cc[i__ + ((k << 2) + 3) * cc_dim1] - cc[ic + ((k << 2) + 2) * cc_dim1];
	    tr4 = cc[i__ + ((k << 2) + 3) * cc_dim1] + cc[ic + ((k << 2) + 2) * cc_dim1];
	    tr1 = cc[i__ - 1 + ((k << 2) + 1) * cc_dim1] - cc[ic - 1 + ((k << 2) + 4) * cc_dim1];
	    tr2 = cc[i__ - 1 + ((k << 2) + 1) * cc_dim1] + cc[ic - 1 + ((k << 2) + 4) * cc_dim1];
	    ti4 = cc[i__ - 1 + ((k << 2) + 3) * cc_dim1] - cc[ic - 1 + ((k << 2) + 2) * cc_dim1];
	    tr3 = cc[i__ - 1 + ((k << 2) + 3) * cc_dim1] + cc[ic - 1 + ((k << 2) + 2) * cc_dim1];
	    ch[i__ - 1 + (k + ch_dim2) * ch_dim1] = tr2 + tr3;
	    cr3 = tr2 - tr3;
	    ch[i__ + (k + ch_dim2) * ch_dim1] = ti2 + ti3;
	    ci3 = ti2 - ti3;
	    cr2 = tr1 - tr4;
	    cr4 = tr1 + tr4;
	    ci2 = ti1 + ti4;
	    ci4 = ti1 - ti4;
	    ch[i__ - 1 + (k + (ch_dim2 << 1)) * ch_dim1] = wa1[i__ - 2] * cr2 - wa1[i__ - 1] * ci2;
	    ch[i__ + (k + (ch_dim2 << 1)) * ch_dim1] = wa1[i__ - 2] * ci2 + wa1[i__ - 1] * cr2;
	    ch[i__ - 1 + (k + ch_dim2 * 3) * ch_dim1] = wa2[i__ - 2] * cr3 - wa2[i__ - 1] * ci3;
	    ch[i__ + (k + ch_dim2 * 3) * ch_dim1] = wa2[i__ - 2] * ci3 + wa2[i__ - 1] * cr3;
	    ch[i__ - 1 + (k + (ch_dim2 << 2)) * ch_dim1] = wa3[i__ - 2] * cr4 - wa3[i__ - 1] * ci4;
	    ch[i__ + (k + (ch_dim2 << 2)) * ch_dim1] = wa3[i__ - 2] * ci4 + wa3[i__ - 1] * cr4;
/* L103: */
	}
/* L104: */
	}
	if (*ido % 2 == 1) {
	return;
	}
L105:
	i__1 = *l1;
	for (k = 1; k <= i__1; ++k) {
	ti1 = cc[((k << 2) + 2) * cc_dim1 + 1] + cc[((k << 2) + 4) * cc_dim1 + 1];
	ti2 = cc[((k << 2) + 4) * cc_dim1 + 1] - cc[((k << 2) + 2) * cc_dim1 + 1];
	tr1 = cc[*ido + ((k << 2) + 1) * cc_dim1] - cc[*ido + ((k << 2) + 3) *cc_dim1];
	tr2 = cc[*ido + ((k << 2) + 1) * cc_dim1] + cc[*ido + ((k << 2) + 3) *cc_dim1];
	ch[*ido + (k + ch_dim2) * ch_dim1] = tr2 + tr2;
	ch[*ido + (k + (ch_dim2 << 1)) * ch_dim1] = sqrt2 * (tr1 - ti1);
	ch[*ido + (k + ch_dim2 * 3) * ch_dim1] = ti2 + ti2;
	ch[*ido + (k + (ch_dim2 << 2)) * ch_dim1] = -sqrt2 * (tr1 + ti1);
/* L106: */
	}
L107:
	return;
} /* radb4_ */

/* Subroutine */ static void s_radb5(int *ido, int *l1, real_t *cc, 
	real_t *ch, real_t *wa1, real_t *wa2, real_t *wa3, 
	real_t *wa4)
{
	/* Initialized data */

	static real_t tr11 = POST(0.3090169943749474241022934171828195886015458990288143106772431137);
	static real_t ti11 = POST(0.9510565162951535721164393337938214340569863412575022244730564442);
	static real_t tr12 = POST(-0.8090169943749474241022934171828190588601545899028814310677431135);
	static real_t ti12 = POST(0.5877852522924731291687059546390727685976524376431459107227248076);

	/* System generated locals */
	int cc_dim1, cc_offset, ch_dim1, ch_dim2, ch_offset, i__1, i__2;

	/* Local variables */
	int i__, k, ic;
	real_t ci2, ci3, ci4, ci5, di3, di4, di5, di2, cr2, cr3, cr5, cr4, 
		ti2, ti3, ti4, ti5, dr3, dr4, dr5, dr2, tr2, tr3, tr4, tr5;
	int idp2;

	/* Parameter adjustments */
	ch_dim1 = *ido;
	ch_dim2 = *l1;
	ch_offset = 1 + ch_dim1 * (1 + ch_dim2);
	ch -= ch_offset;
	cc_dim1 = *ido;
	cc_offset = 1 + cc_dim1 * 6;
	cc -= cc_offset;
	--wa1;
	--wa2;
	--wa3;
	--wa4;

	/* Function Body */
	i__1 = *l1;
	for (k = 1; k <= i__1; ++k) {
	ti5 = cc[(k * 5 + 3) * cc_dim1 + 1] + cc[(k * 5 + 3) * cc_dim1 + 1];
	ti4 = cc[(k * 5 + 5) * cc_dim1 + 1] + cc[(k * 5 + 5) * cc_dim1 + 1];
	tr2 = cc[*ido + (k * 5 + 2) * cc_dim1] + cc[*ido + (k * 5 + 2) * cc_dim1];
	tr3 = cc[*ido + (k * 5 + 4) * cc_dim1] + cc[*ido + (k * 5 + 4) * cc_dim1];
	ch[(k + ch_dim2) * ch_dim1 + 1] = cc[(k * 5 + 1) * cc_dim1 + 1] + tr2 + tr3;
	cr2 = cc[(k * 5 + 1) * cc_dim1 + 1] + tr11 * tr2 + tr12 * tr3;
	cr3 = cc[(k * 5 + 1) * cc_dim1 + 1] + tr12 * tr2 + tr11 * tr3;
	ci5 = ti11 * ti5 + ti12 * ti4;
	ci4 = ti12 * ti5 - ti11 * ti4;
	ch[(k + (ch_dim2 << 1)) * ch_dim1 + 1] = cr2 - ci5;
	ch[(k + ch_dim2 * 3) * ch_dim1 + 1] = cr3 - ci4;
	ch[(k + (ch_dim2 << 2)) * ch_dim1 + 1] = cr3 + ci4;
	ch[(k + ch_dim2 * 5) * ch_dim1 + 1] = cr2 + ci5;
/* L101: */
	}
	if (*ido == 1) {
	return;
	}
	idp2 = *ido + 2;
	i__1 = *l1;
	for (k = 1; k <= i__1; ++k) {
	i__2 = *ido;
	for (i__ = 3; i__ <= i__2; i__ += 2) {
	    ic = idp2 - i__;
	    ti5 = cc[i__ + (k * 5 + 3) * cc_dim1] + cc[ic + (k * 5 + 2) * cc_dim1];
	    ti2 = cc[i__ + (k * 5 + 3) * cc_dim1] - cc[ic + (k * 5 + 2) * cc_dim1];
	    ti4 = cc[i__ + (k * 5 + 5) * cc_dim1] + cc[ic + (k * 5 + 4) * cc_dim1];
	    ti3 = cc[i__ + (k * 5 + 5) * cc_dim1] - cc[ic + (k * 5 + 4) * cc_dim1];
	    tr5 = cc[i__ - 1 + (k * 5 + 3) * cc_dim1] - cc[ic - 1 + (k * 5 + 2) * cc_dim1];
	    tr2 = cc[i__ - 1 + (k * 5 + 3) * cc_dim1] + cc[ic - 1 + (k * 5 + 2) * cc_dim1];
	    tr4 = cc[i__ - 1 + (k * 5 + 5) * cc_dim1] - cc[ic - 1 + (k * 5 + 4) * cc_dim1];
	    tr3 = cc[i__ - 1 + (k * 5 + 5) * cc_dim1] + cc[ic - 1 + (k * 5 + 4) * cc_dim1];
	    ch[i__ - 1 + (k + ch_dim2) * ch_dim1] = cc[i__ - 1 + (k * 5 + 1) *cc_dim1] + tr2 + tr3;
	    ch[i__ + (k + ch_dim2) * ch_dim1] = cc[i__ + (k * 5 + 1) * cc_dim1] + ti2 + ti3;
	    cr2 = cc[i__ - 1 + (k * 5 + 1) * cc_dim1] + tr11 * tr2 + tr12 * tr3;
	    ci2 = cc[i__ + (k * 5 + 1) * cc_dim1] + tr11 * ti2 + tr12 * ti3;
	    cr3 = cc[i__ - 1 + (k * 5 + 1) * cc_dim1] + tr12 * tr2 + tr11 * tr3;
	    ci3 = cc[i__ + (k * 5 + 1) * cc_dim1] + tr12 * ti2 + tr11 * ti3;
	    cr5 = ti11 * tr5 + ti12 * tr4;
	    ci5 = ti11 * ti5 + ti12 * ti4;
	    cr4 = ti12 * tr5 - ti11 * tr4;
	    ci4 = ti12 * ti5 - ti11 * ti4;
	    dr3 = cr3 - ci4;
	    dr4 = cr3 + ci4;
	    di3 = ci3 + cr4;
	    di4 = ci3 - cr4;
	    dr5 = cr2 + ci5;
	    dr2 = cr2 - ci5;
	    di5 = ci2 - cr5;
	    di2 = ci2 + cr5;
	    ch[i__ - 1 + (k + (ch_dim2 << 1)) * ch_dim1] = wa1[i__ - 2] * dr2 - wa1[i__ - 1] * di2;
	    ch[i__ + (k + (ch_dim2 << 1)) * ch_dim1] = wa1[i__ - 2] * di2 + wa1[i__ - 1] * dr2;
	    ch[i__ - 1 + (k + ch_dim2 * 3) * ch_dim1] = wa2[i__ - 2] * dr3 - wa2[i__ - 1] * di3;
	    ch[i__ + (k + ch_dim2 * 3) * ch_dim1] = wa2[i__ - 2] * di3 + wa2[i__ - 1] * dr3;
	    ch[i__ - 1 + (k + (ch_dim2 << 2)) * ch_dim1] = wa3[i__ - 2] * dr4 - wa3[i__ - 1] * di4;
	    ch[i__ + (k + (ch_dim2 << 2)) * ch_dim1] = wa3[i__ - 2] * di4 + wa3[i__ - 1] * dr4;
	    ch[i__ - 1 + (k + ch_dim2 * 5) * ch_dim1] = wa4[i__ - 2] * dr5 - wa4[i__ - 1] * di5;
	    ch[i__ + (k + ch_dim2 * 5) * ch_dim1] = wa4[i__ - 2] * di5 + wa4[i__ - 1] * dr5;
/* L102: */
	}
/* L103: */
	}
	return;
} /* radb5_ */

/* Subroutine */ static void s_radbg(int *ido, int *ip, int *l1, int *
	idl1, real_t *cc, real_t *c1, real_t *c2, real_t *ch, 
	real_t *ch2, real_t *wa)
{
	/* Initialized data */

	static real_t tpi = POST(6.283185307179586476925286766559005768394338798750116419498891846);

	/* System generated locals */
	int ch_dim1, ch_dim2, ch_offset, cc_dim1, cc_dim2, cc_offset, c1_dim1,
		c1_dim2, c1_offset, c2_dim1, c2_offset, ch2_dim1, ch2_offset, 
		i__1, i__2, i__3;

	/* Local variables */
	int i__, j, k, l, j2, ic, jc, lc, ik, is;
	real_t dc2, ai1, ai2, ar1, ar2, ds2;
	int nbd;
	real_t dcp, arg, dsp, ar1h, ar2h;
	int idp2, ipp2, idij, ipph;

	/* Parameter adjustments */
	ch_dim1 = *ido;
	ch_dim2 = *l1;
	ch_offset = 1 + ch_dim1 * (1 + ch_dim2);
	ch -= ch_offset;
	c1_dim1 = *ido;
	c1_dim2 = *l1;
	c1_offset = 1 + c1_dim1 * (1 + c1_dim2);
	c1 -= c1_offset;
	cc_dim1 = *ido;
	cc_dim2 = *ip;
	cc_offset = 1 + cc_dim1 * (1 + cc_dim2);
	cc -= cc_offset;
	ch2_dim1 = *idl1;
	ch2_offset = 1 + ch2_dim1;
	ch2 -= ch2_offset;
	c2_dim1 = *idl1;
	c2_offset = 1 + c2_dim1;
	c2 -= c2_offset;
	--wa;

	/* Function Body */
	arg = tpi / (real_t) (*ip);
	dcp = POST(cos)(arg);
	dsp = POST(sin)(arg);
	idp2 = *ido + 2;
	nbd = (*ido - 1) / 2;
	ipp2 = *ip + 2;
	ipph = (*ip + 1) / 2;
	if (*ido < *l1) {
	goto L103;
	}
	i__1 = *l1;
	for (k = 1; k <= i__1; ++k) {
	i__2 = *ido;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    ch[i__ + (k + ch_dim2) * ch_dim1] = cc[i__ + (k * cc_dim2 + 1) * 
		    cc_dim1];
/* L101: */
	}
/* L102: */
	}
	goto L106;
L103:
	i__1 = *ido;
	for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    ch[i__ + (k + ch_dim2) * ch_dim1] = cc[i__ + (k * cc_dim2 + 1) * 
		    cc_dim1];
/* L104: */
	}
/* L105: */
	}
L106:
	i__1 = ipph;
	for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	j2 = j + j;
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    ch[(k + j * ch_dim2) * ch_dim1 + 1] = cc[*ido + (j2 - 2 + k * 
		    cc_dim2) * cc_dim1] + cc[*ido + (j2 - 2 + k * cc_dim2) * 
		    cc_dim1];
	    ch[(k + jc * ch_dim2) * ch_dim1 + 1] = cc[(j2 - 1 + k * cc_dim2) *
		     cc_dim1 + 1] + cc[(j2 - 1 + k * cc_dim2) * cc_dim1 + 1];
/* L107: */
	}
/* L108: */
	}
	if (*ido == 1) {
	goto L116;
	}
	if (nbd < *l1) {
	goto L112;
	}
	i__1 = ipph;
	for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    i__3 = *ido;
	    for (i__ = 3; i__ <= i__3; i__ += 2) {
		ic = idp2 - i__;
		ch[i__ - 1 + (k + j * ch_dim2) * ch_dim1] = cc[i__ - 1 + ((j 
			<< 1) - 1 + k * cc_dim2) * cc_dim1] + cc[ic - 1 + ((j 
			<< 1) - 2 + k * cc_dim2) * cc_dim1];
		ch[i__ - 1 + (k + jc * ch_dim2) * ch_dim1] = cc[i__ - 1 + ((j 
			<< 1) - 1 + k * cc_dim2) * cc_dim1] - cc[ic - 1 + ((j 
			<< 1) - 2 + k * cc_dim2) * cc_dim1];
		ch[i__ + (k + j * ch_dim2) * ch_dim1] = cc[i__ + ((j << 1) - 
			1 + k * cc_dim2) * cc_dim1] - cc[ic + ((j << 1) - 2 + 
			k * cc_dim2) * cc_dim1];
		ch[i__ + (k + jc * ch_dim2) * ch_dim1] = cc[i__ + ((j << 1) - 
			1 + k * cc_dim2) * cc_dim1] + cc[ic + ((j << 1) - 2 + 
			k * cc_dim2) * cc_dim1];
/* L109: */
	    }
/* L110: */
	}
/* L111: */
	}
	goto L116;
L112:
	i__1 = ipph;
	for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	i__2 = *ido;
	for (i__ = 3; i__ <= i__2; i__ += 2) {
	    ic = idp2 - i__;
	    i__3 = *l1;
	    for (k = 1; k <= i__3; ++k) {
		ch[i__ - 1 + (k + j * ch_dim2) * ch_dim1] = cc[i__ - 1 + ((j 
			<< 1) - 1 + k * cc_dim2) * cc_dim1] + cc[ic - 1 + ((j 
			<< 1) - 2 + k * cc_dim2) * cc_dim1];
		ch[i__ - 1 + (k + jc * ch_dim2) * ch_dim1] = cc[i__ - 1 + ((j 
			<< 1) - 1 + k * cc_dim2) * cc_dim1] - cc[ic - 1 + ((j 
			<< 1) - 2 + k * cc_dim2) * cc_dim1];
		ch[i__ + (k + j * ch_dim2) * ch_dim1] = cc[i__ + ((j << 1) - 
			1 + k * cc_dim2) * cc_dim1] - cc[ic + ((j << 1) - 2 + 
			k * cc_dim2) * cc_dim1];
		ch[i__ + (k + jc * ch_dim2) * ch_dim1] = cc[i__ + ((j << 1) - 
			1 + k * cc_dim2) * cc_dim1] + cc[ic + ((j << 1) - 2 + 
			k * cc_dim2) * cc_dim1];
/* L113: */
	    }
/* L114: */
	}
/* L115: */
	}
L116:
	ar1 = POST(1.0);
	ai1 = POST(0.0);
	i__1 = ipph;
	for (l = 2; l <= i__1; ++l) {
	lc = ipp2 - l;
	ar1h = dcp * ar1 - dsp * ai1;
	ai1 = dcp * ai1 + dsp * ar1;
	ar1 = ar1h;
	i__2 = *idl1;
	for (ik = 1; ik <= i__2; ++ik) {
	    c2[ik + l * c2_dim1] = ch2[ik + ch2_dim1] + ar1 * ch2[ik + (ch2_dim1 << 1)];
	    c2[ik + lc * c2_dim1] = ai1 * ch2[ik + *ip * ch2_dim1];
/* L117: */
	}
	dc2 = ar1;
	ds2 = ai1;
	ar2 = ar1;
	ai2 = ai1;
	i__2 = ipph;
	for (j = 3; j <= i__2; ++j) {
	    jc = ipp2 - j;
	    ar2h = dc2 * ar2 - ds2 * ai2;
	    ai2 = dc2 * ai2 + ds2 * ar2;
	    ar2 = ar2h;
	    i__3 = *idl1;
	    for (ik = 1; ik <= i__3; ++ik) {
		c2[ik + l * c2_dim1] += ar2 * ch2[ik + j * ch2_dim1];
		c2[ik + lc * c2_dim1] += ai2 * ch2[ik + jc * ch2_dim1];
/* L118: */
	    }
/* L119: */
	}
/* L120: */
	}
	i__1 = ipph;
	for (j = 2; j <= i__1; ++j) {
	i__2 = *idl1;
	for (ik = 1; ik <= i__2; ++ik) {
	    ch2[ik + ch2_dim1] += ch2[ik + j * ch2_dim1];
/* L121: */
	}
/* L122: */
	}
	i__1 = ipph;
	for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    ch[(k + j * ch_dim2) * ch_dim1 + 1] = c1[(k + j * c1_dim2) * c1_dim1 + 1] - c1[(k + jc * c1_dim2) * c1_dim1 + 1];
	    ch[(k + jc * ch_dim2) * ch_dim1 + 1] = c1[(k + j * c1_dim2) * c1_dim1 + 1] + c1[(k + jc * c1_dim2) * c1_dim1 + 1];
/* L123: */
	}
/* L124: */
	}
	if (*ido == 1) {
	goto L132;
	}
	if (nbd < *l1) {
	goto L128;
	}
	i__1 = ipph;
	for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    i__3 = *ido;
	    for (i__ = 3; i__ <= i__3; i__ += 2) {
		ch[i__ - 1 + (k + j * ch_dim2) * ch_dim1] = c1[i__ - 1 + (k + 
			j * c1_dim2) * c1_dim1] - c1[i__ + (k + jc * c1_dim2) * c1_dim1];
		ch[i__ - 1 + (k + jc * ch_dim2) * ch_dim1] = c1[i__ - 1 + (k 
			+ j * c1_dim2) * c1_dim1] + c1[i__ + (k + jc * c1_dim2) * c1_dim1];
		ch[i__ + (k + j * ch_dim2) * ch_dim1] = c1[i__ + (k + j * 
			c1_dim2) * c1_dim1] + c1[i__ - 1 + (k + jc * c1_dim2) * c1_dim1];
		ch[i__ + (k + jc * ch_dim2) * ch_dim1] = c1[i__ + (k + j * 
			c1_dim2) * c1_dim1] - c1[i__ - 1 + (k + jc * c1_dim2) * c1_dim1];
/* L125: */
	    }
/* L126: */
	}
/* L127: */
	}
	goto L132;
L128:
	i__1 = ipph;
	for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	i__2 = *ido;
	for (i__ = 3; i__ <= i__2; i__ += 2) {
	    i__3 = *l1;
	    for (k = 1; k <= i__3; ++k) {
		ch[i__ - 1 + (k + j * ch_dim2) * ch_dim1] = c1[i__ - 1 + (k + 
			j * c1_dim2) * c1_dim1] - c1[i__ + (k + jc * c1_dim2) * c1_dim1];
		ch[i__ - 1 + (k + jc * ch_dim2) * ch_dim1] = c1[i__ - 1 + (k 
			+ j * c1_dim2) * c1_dim1] + c1[i__ + (k + jc * c1_dim2) * c1_dim1];
		ch[i__ + (k + j * ch_dim2) * ch_dim1] = c1[i__ + (k + j * 
			c1_dim2) * c1_dim1] + c1[i__ - 1 + (k + jc * c1_dim2) * c1_dim1];
		ch[i__ + (k + jc * ch_dim2) * ch_dim1] = c1[i__ + (k + j * 
			c1_dim2) * c1_dim1] - c1[i__ - 1 + (k + jc * c1_dim2) * c1_dim1];
/* L129: */
	    }
/* L130: */
	}
/* L131: */
	}
L132:
	if (*ido == 1) {
	return;
	}
	i__1 = *idl1;
	for (ik = 1; ik <= i__1; ++ik) {
	c2[ik + c2_dim1] = ch2[ik + ch2_dim1];
/* L133: */
	}
	i__1 = *ip;
	for (j = 2; j <= i__1; ++j) {
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    c1[(k + j * c1_dim2) * c1_dim1 + 1] = ch[(k + j * ch_dim2) * ch_dim1 + 1];
/* L134: */
	}
/* L135: */
	}
	if (nbd > *l1) {
	goto L139;
	}
	is = -(*ido);
	i__1 = *ip;
	for (j = 2; j <= i__1; ++j) {
	is += *ido;
	idij = is;
	i__2 = *ido;
	for (i__ = 3; i__ <= i__2; i__ += 2) {
	    idij += 2;
	    i__3 = *l1;
	    for (k = 1; k <= i__3; ++k) {
		c1[i__ - 1 + (k + j * c1_dim2) * c1_dim1] = wa[idij - 1] * ch[
			i__ - 1 + (k + j * ch_dim2) * ch_dim1] - wa[idij] * 
			ch[i__ + (k + j * ch_dim2) * ch_dim1];
		c1[i__ + (k + j * c1_dim2) * c1_dim1] = wa[idij - 1] * ch[i__ 
			+ (k + j * ch_dim2) * ch_dim1] + wa[idij] * ch[i__ - 
			1 + (k + j * ch_dim2) * ch_dim1];
/* L136: */
	    }
/* L137: */
	}
/* L138: */
	}
	goto L143;
L139:
	is = -(*ido);
	i__1 = *ip;
	for (j = 2; j <= i__1; ++j) {
	is += *ido;
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    idij = is;
	    i__3 = *ido;
	    for (i__ = 3; i__ <= i__3; i__ += 2) {
		idij += 2;
		c1[i__ - 1 + (k + j * c1_dim2) * c1_dim1] = wa[idij - 1] * ch[
			i__ - 1 + (k + j * ch_dim2) * ch_dim1] - wa[idij] * 
			ch[i__ + (k + j * ch_dim2) * ch_dim1];
		c1[i__ + (k + j * c1_dim2) * c1_dim1] = wa[idij - 1] * ch[i__ 
			+ (k + j * ch_dim2) * ch_dim1] + wa[idij] * ch[i__ - 
			1 + (k + j * ch_dim2) * ch_dim1];
/* L140: */
	    }
/* L141: */
	}
/* L142: */
	}
L143:
	return;
} /* radbg_ */

/* Subroutine */ static void s_radf2(int *ido, int *l1, real_t *cc, 
	real_t *ch, real_t *wa1)
{
	/* System generated locals */
	int ch_dim1, ch_offset, cc_dim1, cc_dim2, cc_offset, i__1, i__2;

	/* Local variables */
	int i__, k, ic;
	real_t ti2, tr2;
	int idp2;

	/* Parameter adjustments */
	ch_dim1 = *ido;
	ch_offset = 1 + ch_dim1 * 3;
	ch -= ch_offset;
	cc_dim1 = *ido;
	cc_dim2 = *l1;
	cc_offset = 1 + cc_dim1 * (1 + cc_dim2);
	cc -= cc_offset;
	--wa1;

	/* Function Body */
	i__1 = *l1;
	for (k = 1; k <= i__1; ++k) {
	ch[((k << 1) + 1) * ch_dim1 + 1] = cc[(k + cc_dim2) * cc_dim1 + 1] + cc[(k + (cc_dim2 << 1)) * cc_dim1 + 1];
	ch[*ido + ((k << 1) + 2) * ch_dim1] = cc[(k + cc_dim2) * cc_dim1 + 1] - cc[(k + (cc_dim2 << 1)) * cc_dim1 + 1];
/* L101: */
	}
	if ((i__1 = *ido - 2) < 0) {
	goto L107;
	} else if (i__1 == 0) {
	goto L105;
	} else {
	goto L102;
	}
L102:
	idp2 = *ido + 2;
	i__1 = *l1;
	for (k = 1; k <= i__1; ++k) {
	i__2 = *ido;
	for (i__ = 3; i__ <= i__2; i__ += 2) {
	    ic = idp2 - i__;
	    tr2 = wa1[i__ - 2] * cc[i__ - 1 + (k + (cc_dim2 << 1)) * cc_dim1] 
		    + wa1[i__ - 1] * cc[i__ + (k + (cc_dim2 << 1)) * cc_dim1];
	    ti2 = wa1[i__ - 2] * cc[i__ + (k + (cc_dim2 << 1)) * cc_dim1] - 
		    wa1[i__ - 1] * cc[i__ - 1 + (k + (cc_dim2 << 1)) * cc_dim1];
	    ch[i__ + ((k << 1) + 1) * ch_dim1] = cc[i__ + (k + cc_dim2) * cc_dim1] + ti2;
	    ch[ic + ((k << 1) + 2) * ch_dim1] = ti2 - cc[i__ + (k + cc_dim2) *cc_dim1];
	    ch[i__ - 1 + ((k << 1) + 1) * ch_dim1] = cc[i__ - 1 + (k + cc_dim2) * cc_dim1] + tr2;
	    ch[ic - 1 + ((k << 1) + 2) * ch_dim1] = cc[i__ - 1 + (k + cc_dim2)* cc_dim1] - tr2;
/* L103: */
	}
/* L104: */
	}
	if (*ido % 2 == 1) {
	return;
	}
L105:
	i__1 = *l1;
	for (k = 1; k <= i__1; ++k) {
	ch[((k << 1) + 2) * ch_dim1 + 1] = -cc[*ido + (k + (cc_dim2 << 1)) * cc_dim1];
	ch[*ido + ((k << 1) + 1) * ch_dim1] = cc[*ido + (k + cc_dim2) * cc_dim1];
/* L106: */
	}
L107:
	return;
} /* radf2_ */

/* Subroutine */ static void s_radf3(int *ido, int *l1, real_t *cc, 
	real_t *ch, real_t *wa1, real_t *wa2)
{
	/* Initialized data */

	static real_t taur = POST(-0.5);
	static real_t taui = POST(0.8660254037844386467637231707529361834710262690519031402790348975);

	/* System generated locals */
	int ch_dim1, ch_offset, cc_dim1, cc_dim2, cc_offset, i__1, i__2;

	/* Local variables */
	int i__, k, ic;
	real_t ci2, di2, di3, cr2, dr2, dr3, ti2, ti3, tr2, tr3;
	int idp2;

	/* Parameter adjustments */
	ch_dim1 = *ido;
	ch_offset = 1 + (ch_dim1 << 2);
	ch -= ch_offset;
	cc_dim1 = *ido;
	cc_dim2 = *l1;
	cc_offset = 1 + cc_dim1 * (1 + cc_dim2);
	cc -= cc_offset;
	--wa1;
	--wa2;

	/* Function Body */
	i__1 = *l1;
	for (k = 1; k <= i__1; ++k) {
	cr2 = cc[(k + (cc_dim2 << 1)) * cc_dim1 + 1] + cc[(k + cc_dim2 * 3) * cc_dim1 + 1];
	ch[(k * 3 + 1) * ch_dim1 + 1] = cc[(k + cc_dim2) * cc_dim1 + 1] + cr2;
	ch[(k * 3 + 3) * ch_dim1 + 1] = taui * (cc[(k + cc_dim2 * 3) * cc_dim1 + 1] - cc[(k + (cc_dim2 << 1)) * cc_dim1 + 1]);
	ch[*ido + (k * 3 + 2) * ch_dim1] = cc[(k + cc_dim2) * cc_dim1 + 1] + taur * cr2;
/* L101: */
	}
	if (*ido == 1) {
	return;
	}
	idp2 = *ido + 2;
	i__1 = *l1;
	for (k = 1; k <= i__1; ++k) {
	i__2 = *ido;
	for (i__ = 3; i__ <= i__2; i__ += 2) {
	    ic = idp2 - i__;
	    dr2 = wa1[i__ - 2] * cc[i__ - 1 + (k + (cc_dim2 << 1)) * cc_dim1] + wa1[i__ - 1] * cc[i__ + (k + (cc_dim2 << 1)) * cc_dim1];
	    di2 = wa1[i__ - 2] * cc[i__ + (k + (cc_dim2 << 1)) * cc_dim1] - wa1[i__ - 1] * cc[i__ - 1 + (k + (cc_dim2 << 1)) * cc_dim1];
	    dr3 = wa2[i__ - 2] * cc[i__ - 1 + (k + cc_dim2 * 3) * cc_dim1] + wa2[i__ - 1] * cc[i__ + (k + cc_dim2 * 3) * cc_dim1];
	    di3 = wa2[i__ - 2] * cc[i__ + (k + cc_dim2 * 3) * cc_dim1] - wa2[i__ - 1] * cc[i__ - 1 + (k + cc_dim2 * 3) * cc_dim1];
	    cr2 = dr2 + dr3;
	    ci2 = di2 + di3;
	    ch[i__ - 1 + (k * 3 + 1) * ch_dim1] = cc[i__ - 1 + (k + cc_dim2) *cc_dim1] + cr2;
	    ch[i__ + (k * 3 + 1) * ch_dim1] = cc[i__ + (k + cc_dim2) * cc_dim1] + ci2;
	    tr2 = cc[i__ - 1 + (k + cc_dim2) * cc_dim1] + taur * cr2;
	    ti2 = cc[i__ + (k + cc_dim2) * cc_dim1] + taur * ci2;
	    tr3 = taui * (di2 - di3);
	    ti3 = taui * (dr3 - dr2);
	    ch[i__ - 1 + (k * 3 + 3) * ch_dim1] = tr2 + tr3;
	    ch[ic - 1 + (k * 3 + 2) * ch_dim1] = tr2 - tr3;
	    ch[i__ + (k * 3 + 3) * ch_dim1] = ti2 + ti3;
	    ch[ic + (k * 3 + 2) * ch_dim1] = ti3 - ti2;
/* L102: */
	}
/* L103: */
	}
	return;
} /* radf3_ */

/* Subroutine */ static void s_radf4(int *ido, int *l1, real_t *cc, 
	real_t *ch, real_t *wa1, real_t *wa2, real_t *wa3)
{
	/* Initialized data */

	static real_t hsqt2 = POST(0.70710678118654752440084436210484903928483593768474036588339869);

	/* System generated locals */
	int cc_dim1, cc_dim2, cc_offset, ch_dim1, ch_offset, i__1, i__2;

	/* Local variables */
	int i__, k, ic;
	real_t ci2, ci3, ci4, cr2, cr3, cr4, ti1, ti2, ti3, ti4, tr1, tr2, tr3, tr4;
	int idp2;

	/* Parameter adjustments */
	ch_dim1 = *ido;
	ch_offset = 1 + ch_dim1 * 5;
	ch -= ch_offset;
	cc_dim1 = *ido;
	cc_dim2 = *l1;
	cc_offset = 1 + cc_dim1 * (1 + cc_dim2);
	cc -= cc_offset;
	--wa1;
	--wa2;
	--wa3;

	/* Function Body */
	i__1 = *l1;
	for (k = 1; k <= i__1; ++k) {
	tr1 = cc[(k + (cc_dim2 << 1)) * cc_dim1 + 1] + cc[(k + (cc_dim2 << 2))* cc_dim1 + 1];
	tr2 = cc[(k + cc_dim2) * cc_dim1 + 1] + cc[(k + cc_dim2 * 3) * cc_dim1 + 1];
	ch[((k << 2) + 1) * ch_dim1 + 1] = tr1 + tr2;
	ch[*ido + ((k << 2) + 4) * ch_dim1] = tr2 - tr1;
	ch[*ido + ((k << 2) + 2) * ch_dim1] = cc[(k + cc_dim2) * cc_dim1 + 1] - cc[(k + cc_dim2 * 3) * cc_dim1 + 1];
	ch[((k << 2) + 3) * ch_dim1 + 1] = cc[(k + (cc_dim2 << 2)) * cc_dim1 + 1] - cc[(k + (cc_dim2 << 1)) * cc_dim1 + 1];
/* L101: */
	}
	if ((i__1 = *ido - 2) < 0) {
	goto L107;
	} else if (i__1 == 0) {
	goto L105;
	} else {
	goto L102;
	}
L102:
	idp2 = *ido + 2;
	i__1 = *l1;
	for (k = 1; k <= i__1; ++k) {
	i__2 = *ido;
	for (i__ = 3; i__ <= i__2; i__ += 2) {
	    ic = idp2 - i__;
	    cr2 = wa1[i__ - 2] * cc[i__ - 1 + (k + (cc_dim2 << 1)) * cc_dim1] + wa1[i__ - 1] * cc[i__ + (k + (cc_dim2 << 1)) * cc_dim1];
	    ci2 = wa1[i__ - 2] * cc[i__ + (k + (cc_dim2 << 1)) * cc_dim1] - wa1[i__ - 1] * cc[i__ - 1 + (k + (cc_dim2 << 1)) * cc_dim1];
	    cr3 = wa2[i__ - 2] * cc[i__ - 1 + (k + cc_dim2 * 3) * cc_dim1] + wa2[i__ - 1] * cc[i__ + (k + cc_dim2 * 3) * cc_dim1];
	    ci3 = wa2[i__ - 2] * cc[i__ + (k + cc_dim2 * 3) * cc_dim1] - wa2[i__ - 1] * cc[i__ - 1 + (k + cc_dim2 * 3) * cc_dim1];
	    cr4 = wa3[i__ - 2] * cc[i__ - 1 + (k + (cc_dim2 << 2)) * cc_dim1] + wa3[i__ - 1] * cc[i__ + (k + (cc_dim2 << 2)) * cc_dim1];
	    ci4 = wa3[i__ - 2] * cc[i__ + (k + (cc_dim2 << 2)) * cc_dim1] - wa3[i__ - 1] * cc[i__ - 1 + (k + (cc_dim2 << 2)) * cc_dim1];
	    tr1 = cr2 + cr4;
	    tr4 = cr4 - cr2;
	    ti1 = ci2 + ci4;
	    ti4 = ci2 - ci4;
	    ti2 = cc[i__ + (k + cc_dim2) * cc_dim1] + ci3;
	    ti3 = cc[i__ + (k + cc_dim2) * cc_dim1] - ci3;
	    tr2 = cc[i__ - 1 + (k + cc_dim2) * cc_dim1] + cr3;
	    tr3 = cc[i__ - 1 + (k + cc_dim2) * cc_dim1] - cr3;
	    ch[i__ - 1 + ((k << 2) + 1) * ch_dim1] = tr1 + tr2;
	    ch[ic - 1 + ((k << 2) + 4) * ch_dim1] = tr2 - tr1;
	    ch[i__ + ((k << 2) + 1) * ch_dim1] = ti1 + ti2;
	    ch[ic + ((k << 2) + 4) * ch_dim1] = ti1 - ti2;
	    ch[i__ - 1 + ((k << 2) + 3) * ch_dim1] = ti4 + tr3;
	    ch[ic - 1 + ((k << 2) + 2) * ch_dim1] = tr3 - ti4;
	    ch[i__ + ((k << 2) + 3) * ch_dim1] = tr4 + ti3;
	    ch[ic + ((k << 2) + 2) * ch_dim1] = tr4 - ti3;
/* L103: */
	}
/* L104: */
	}
	if (*ido % 2 == 1) {
	return;
	}
L105:
	i__1 = *l1;
	for (k = 1; k <= i__1; ++k) {
	ti1 = -hsqt2 * (cc[*ido + (k + (cc_dim2 << 1)) * cc_dim1] + cc[*ido + (k + (cc_dim2 << 2)) * cc_dim1]);
	tr1 = hsqt2 * (cc[*ido + (k + (cc_dim2 << 1)) * cc_dim1] - cc[*ido + (k + (cc_dim2 << 2)) * cc_dim1]);
	ch[*ido + ((k << 2) + 1) * ch_dim1] = tr1 + cc[*ido + (k + cc_dim2) * cc_dim1];
	ch[*ido + ((k << 2) + 3) * ch_dim1] = cc[*ido + (k + cc_dim2) * cc_dim1] - tr1;
	ch[((k << 2) + 2) * ch_dim1 + 1] = ti1 - cc[*ido + (k + cc_dim2 * 3) *cc_dim1];
	ch[((k << 2) + 4) * ch_dim1 + 1] = ti1 + cc[*ido + (k + cc_dim2 * 3) *cc_dim1];
/* L106: */
	}
L107:
	return;
} /* radf4_ */

/* Subroutine */ static void s_radf5(int *ido, int *l1, real_t *cc, 
	real_t *ch, real_t *wa1, real_t *wa2, real_t *wa3, 
	real_t *wa4)
{
	/* Initialized data */

	static real_t tr11 = POST(0.3090169943749474241022934171828195886015458990288143106772431137);
	static real_t ti11 = POST(0.9510565162951535721164393337938214340569863412575022244730564442);
	static real_t tr12 = POST(-0.8090169943749474241022934171828190588601545899028814310677431135);
	static real_t ti12 = POST(0.5877852522924731291687059546390727685976524376431459107227248076);

	/* System generated locals */
	int cc_dim1, cc_dim2, cc_offset, ch_dim1, ch_offset, i__1, i__2;

	/* Local variables */
	int i__, k, ic;
	real_t ci2, di2, ci4, ci5, di3, di4, di5, ci3, cr2, cr3, dr2, dr3, 
		dr4, dr5, cr5, cr4, ti2, ti3, ti5, ti4, tr2, tr3, tr4, tr5;
	int idp2;

	/* Parameter adjustments */
	ch_dim1 = *ido;
	ch_offset = 1 + ch_dim1 * 6;
	ch -= ch_offset;
	cc_dim1 = *ido;
	cc_dim2 = *l1;
	cc_offset = 1 + cc_dim1 * (1 + cc_dim2);
	cc -= cc_offset;
	--wa1;
	--wa2;
	--wa3;
	--wa4;

	/* Function Body */
	i__1 = *l1;
	for (k = 1; k <= i__1; ++k) {
	cr2 = cc[(k + cc_dim2 * 5) * cc_dim1 + 1] + cc[(k + (cc_dim2 << 1)) * cc_dim1 + 1];
	ci5 = cc[(k + cc_dim2 * 5) * cc_dim1 + 1] - cc[(k + (cc_dim2 << 1)) * cc_dim1 + 1];
	cr3 = cc[(k + (cc_dim2 << 2)) * cc_dim1 + 1] + cc[(k + cc_dim2 * 3) * cc_dim1 + 1];
	ci4 = cc[(k + (cc_dim2 << 2)) * cc_dim1 + 1] - cc[(k + cc_dim2 * 3) * cc_dim1 + 1];
	ch[(k * 5 + 1) * ch_dim1 + 1] = cc[(k + cc_dim2) * cc_dim1 + 1] + cr2 + cr3;
	ch[*ido + (k * 5 + 2) * ch_dim1] = cc[(k + cc_dim2) * cc_dim1 + 1] + tr11 * cr2 + tr12 * cr3;
	ch[(k * 5 + 3) * ch_dim1 + 1] = ti11 * ci5 + ti12 * ci4;
	ch[*ido + (k * 5 + 4) * ch_dim1] = cc[(k + cc_dim2) * cc_dim1 + 1] + tr12 * cr2 + tr11 * cr3;
	ch[(k * 5 + 5) * ch_dim1 + 1] = ti12 * ci5 - ti11 * ci4;
/* L101: */
	}
	if (*ido == 1) {
	return;
	}
	idp2 = *ido + 2;
	i__1 = *l1;
	for (k = 1; k <= i__1; ++k) {
	i__2 = *ido;
	for (i__ = 3; i__ <= i__2; i__ += 2) {
	    ic = idp2 - i__;
	    dr2 = wa1[i__ - 2] * cc[i__ - 1 + (k + (cc_dim2 << 1)) * cc_dim1] + wa1[i__ - 1] * cc[i__ + (k + (cc_dim2 << 1)) * cc_dim1];
	    di2 = wa1[i__ - 2] * cc[i__ + (k + (cc_dim2 << 1)) * cc_dim1] - wa1[i__ - 1] * cc[i__ - 1 + (k + (cc_dim2 << 1)) * cc_dim1];
	    dr3 = wa2[i__ - 2] * cc[i__ - 1 + (k + cc_dim2 * 3) * cc_dim1] + wa2[i__ - 1] * cc[i__ + (k + cc_dim2 * 3) * cc_dim1];
	    di3 = wa2[i__ - 2] * cc[i__ + (k + cc_dim2 * 3) * cc_dim1] - wa2[i__ - 1] * cc[i__ - 1 + (k + cc_dim2 * 3) * cc_dim1];
	    dr4 = wa3[i__ - 2] * cc[i__ - 1 + (k + (cc_dim2 << 2)) * cc_dim1] + wa3[i__ - 1] * cc[i__ + (k + (cc_dim2 << 2)) * cc_dim1];
	    di4 = wa3[i__ - 2] * cc[i__ + (k + (cc_dim2 << 2)) * cc_dim1] - wa3[i__ - 1] * cc[i__ - 1 + (k + (cc_dim2 << 2)) * cc_dim1];
	    dr5 = wa4[i__ - 2] * cc[i__ - 1 + (k + cc_dim2 * 5) * cc_dim1] + wa4[i__ - 1] * cc[i__ + (k + cc_dim2 * 5) * cc_dim1];
	    di5 = wa4[i__ - 2] * cc[i__ + (k + cc_dim2 * 5) * cc_dim1] - wa4[i__ - 1] * cc[i__ - 1 + (k + cc_dim2 * 5) * cc_dim1];
	    cr2 = dr2 + dr5;
	    ci5 = dr5 - dr2;
	    cr5 = di2 - di5;
	    ci2 = di2 + di5;
	    cr3 = dr3 + dr4;
	    ci4 = dr4 - dr3;
	    cr4 = di3 - di4;
	    ci3 = di3 + di4;
	    ch[i__ - 1 + (k * 5 + 1) * ch_dim1] = cc[i__ - 1 + (k + cc_dim2) *cc_dim1] + cr2 + cr3;
	    ch[i__ + (k * 5 + 1) * ch_dim1] = cc[i__ + (k + cc_dim2) * cc_dim1] + ci2 + ci3;
	    tr2 = cc[i__ - 1 + (k + cc_dim2) * cc_dim1] + tr11 * cr2 + tr12 * cr3;
	    ti2 = cc[i__ + (k + cc_dim2) * cc_dim1] + tr11 * ci2 + tr12 * ci3;
	    tr3 = cc[i__ - 1 + (k + cc_dim2) * cc_dim1] + tr12 * cr2 + tr11 * cr3;
	    ti3 = cc[i__ + (k + cc_dim2) * cc_dim1] + tr12 * ci2 + tr11 * ci3;
	    tr5 = ti11 * cr5 + ti12 * cr4;
	    ti5 = ti11 * ci5 + ti12 * ci4;
	    tr4 = ti12 * cr5 - ti11 * cr4;
	    ti4 = ti12 * ci5 - ti11 * ci4;
	    ch[i__ - 1 + (k * 5 + 3) * ch_dim1] = tr2 + tr5;
	    ch[ic - 1 + (k * 5 + 2) * ch_dim1] = tr2 - tr5;
	    ch[i__ + (k * 5 + 3) * ch_dim1] = ti2 + ti5;
	    ch[ic + (k * 5 + 2) * ch_dim1] = ti5 - ti2;
	    ch[i__ - 1 + (k * 5 + 5) * ch_dim1] = tr3 + tr4;
	    ch[ic - 1 + (k * 5 + 4) * ch_dim1] = tr3 - tr4;
	    ch[i__ + (k * 5 + 5) * ch_dim1] = ti3 + ti4;
	    ch[ic + (k * 5 + 4) * ch_dim1] = ti4 - ti3;
/* L102: */
	}
/* L103: */
	}
	return;
} /* radf5_ */

/* Subroutine */ static void s_radfg(int *ido, int *ip, int *l1, int *
	idl1, real_t *cc, real_t *c1, real_t *c2, real_t *ch, 
	real_t *ch2, real_t *wa)
{
	/* Initialized data */

	static real_t tpi = POST(6.283185307179586476925286766559005768394338798750116419498891846);

	/* System generated locals */
	int ch_dim1, ch_dim2, ch_offset, cc_dim1, cc_dim2, cc_offset, c1_dim1,
		c1_dim2, c1_offset, c2_dim1, c2_offset, ch2_dim1, ch2_offset, 
		i__1, i__2, i__3;

	/* Local variables */
	int i__, j, k, l, j2, ic, jc, lc, ik, is;
	real_t dc2, ai1, ai2, ar1, ar2, ds2;
	int nbd;
	real_t dcp, arg, dsp, ar1h, ar2h;
	int idp2, ipp2, idij, ipph;

	/* Parameter adjustments */
	ch_dim1 = *ido;
	ch_dim2 = *l1;
	ch_offset = 1 + ch_dim1 * (1 + ch_dim2);
	ch -= ch_offset;
	c1_dim1 = *ido;
	c1_dim2 = *l1;
	c1_offset = 1 + c1_dim1 * (1 + c1_dim2);
	c1 -= c1_offset;
	cc_dim1 = *ido;
	cc_dim2 = *ip;
	cc_offset = 1 + cc_dim1 * (1 + cc_dim2);
	cc -= cc_offset;
	ch2_dim1 = *idl1;
	ch2_offset = 1 + ch2_dim1;
	ch2 -= ch2_offset;
	c2_dim1 = *idl1;
	c2_offset = 1 + c2_dim1;
	c2 -= c2_offset;
	--wa;

	/* Function Body */
	arg = tpi / (real_t) (*ip);
	dcp = POST(cos)(arg);
	dsp = POST(sin)(arg);
	ipph = (*ip + 1) / 2;
	ipp2 = *ip + 2;
	idp2 = *ido + 2;
	nbd = (*ido - 1) / 2;
	if (*ido == 1) {
	goto L119;
	}
	i__1 = *idl1;
	for (ik = 1; ik <= i__1; ++ik) {
	ch2[ik + ch2_dim1] = c2[ik + c2_dim1];
/* L101: */
	}
	i__1 = *ip;
	for (j = 2; j <= i__1; ++j) {
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    ch[(k + j * ch_dim2) * ch_dim1 + 1] = c1[(k + j * c1_dim2) * c1_dim1 + 1];
/* L102: */
	}
/* L103: */
	}
	if (nbd > *l1) {
	goto L107;
	}
	is = -(*ido);
	i__1 = *ip;
	for (j = 2; j <= i__1; ++j) {
	is += *ido;
	idij = is;
	i__2 = *ido;
	for (i__ = 3; i__ <= i__2; i__ += 2) {
	    idij += 2;
	    i__3 = *l1;
	    for (k = 1; k <= i__3; ++k) {
		ch[i__ - 1 + (k + j * ch_dim2) * ch_dim1] = wa[idij - 1] * c1[
			i__ - 1 + (k + j * c1_dim2) * c1_dim1] + wa[idij] * 
			c1[i__ + (k + j * c1_dim2) * c1_dim1];
		ch[i__ + (k + j * ch_dim2) * ch_dim1] = wa[idij - 1] * c1[i__ 
			+ (k + j * c1_dim2) * c1_dim1] - wa[idij] * c1[i__ - 
			1 + (k + j * c1_dim2) * c1_dim1];
/* L104: */
	    }
/* L105: */
	}
/* L106: */
	}
	goto L111;
L107:
	is = -(*ido);
	i__1 = *ip;
	for (j = 2; j <= i__1; ++j) {
	is += *ido;
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    idij = is;
	    i__3 = *ido;
	    for (i__ = 3; i__ <= i__3; i__ += 2) {
		idij += 2;
		ch[i__ - 1 + (k + j * ch_dim2) * ch_dim1] = wa[idij - 1] * c1[
			i__ - 1 + (k + j * c1_dim2) * c1_dim1] + wa[idij] * 
			c1[i__ + (k + j * c1_dim2) * c1_dim1];
		ch[i__ + (k + j * ch_dim2) * ch_dim1] = wa[idij - 1] * c1[i__ 
			+ (k + j * c1_dim2) * c1_dim1] - wa[idij] * c1[i__ - 
			1 + (k + j * c1_dim2) * c1_dim1];
/* L108: */
	    }
/* L109: */
	}
/* L110: */
	}
L111:
	if (nbd < *l1) {
	goto L115;
	}
	i__1 = ipph;
	for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    i__3 = *ido;
	    for (i__ = 3; i__ <= i__3; i__ += 2) {
		c1[i__ - 1 + (k + j * c1_dim2) * c1_dim1] = ch[i__ - 1 + (k + j * ch_dim2) * ch_dim1] + ch[i__ - 1 + (k + jc * ch_dim2) * ch_dim1];
		c1[i__ - 1 + (k + jc * c1_dim2) * c1_dim1] = ch[i__ + (k + j *ch_dim2) * ch_dim1] - ch[i__ + (k + jc * ch_dim2) * ch_dim1];
		c1[i__ + (k + j * c1_dim2) * c1_dim1] = ch[i__ + (k + j * ch_dim2) * ch_dim1] + ch[i__ + (k + jc * ch_dim2) * ch_dim1];
		c1[i__ + (k + jc * c1_dim2) * c1_dim1] = ch[i__ - 1 + (k + jc * ch_dim2) * ch_dim1] - ch[i__ - 1 + (k + j * ch_dim2)* ch_dim1];
/* L112: */
	    }
/* L113: */
	}
/* L114: */
	}
	goto L121;
L115:
	i__1 = ipph;
	for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	i__2 = *ido;
	for (i__ = 3; i__ <= i__2; i__ += 2) {
	    i__3 = *l1;
	    for (k = 1; k <= i__3; ++k) {
		c1[i__ - 1 + (k + j * c1_dim2) * c1_dim1] = ch[i__ - 1 + (k + j * ch_dim2) * ch_dim1] + ch[i__ - 1 + (k + jc * ch_dim2) * ch_dim1];
		c1[i__ - 1 + (k + jc * c1_dim2) * c1_dim1] = ch[i__ + (k + j *ch_dim2) * ch_dim1] - ch[i__ + (k + jc * ch_dim2) * ch_dim1];
		c1[i__ + (k + j * c1_dim2) * c1_dim1] = ch[i__ + (k + j * ch_dim2) * ch_dim1] + ch[i__ + (k + jc * ch_dim2) * ch_dim1];
		c1[i__ + (k + jc * c1_dim2) * c1_dim1] = ch[i__ - 1 + (k + jc * ch_dim2) * ch_dim1] - ch[i__ - 1 + (k + j * ch_dim2)* ch_dim1];
/* L116: */
	    }
/* L117: */
	}
/* L118: */
	}
	goto L121;
L119:
	i__1 = *idl1;
	for (ik = 1; ik <= i__1; ++ik) {
	c2[ik + c2_dim1] = ch2[ik + ch2_dim1];
/* L120: */
	}
L121:
	i__1 = ipph;
	for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    c1[(k + j * c1_dim2) * c1_dim1 + 1] = ch[(k + j * ch_dim2) * ch_dim1 + 1] + ch[(k + jc * ch_dim2) * ch_dim1 + 1];
	    c1[(k + jc * c1_dim2) * c1_dim1 + 1] = ch[(k + jc * ch_dim2) * ch_dim1 + 1] - ch[(k + j * ch_dim2) * ch_dim1 + 1];
/* L122: */
	}
/* L123: */
	}

	ar1 = POST(1.0);
	ai1 = POST(0.0);
	i__1 = ipph;
	for (l = 2; l <= i__1; ++l) {
	lc = ipp2 - l;
	ar1h = dcp * ar1 - dsp * ai1;
	ai1 = dcp * ai1 + dsp * ar1;
	ar1 = ar1h;
	i__2 = *idl1;
	for (ik = 1; ik <= i__2; ++ik) {
	    ch2[ik + l * ch2_dim1] = c2[ik + c2_dim1] + ar1 * c2[ik + (c2_dim1 << 1)];
	    ch2[ik + lc * ch2_dim1] = ai1 * c2[ik + *ip * c2_dim1];
/* L124: */
	}
	dc2 = ar1;
	ds2 = ai1;
	ar2 = ar1;
	ai2 = ai1;
	i__2 = ipph;
	for (j = 3; j <= i__2; ++j) {
	    jc = ipp2 - j;
	    ar2h = dc2 * ar2 - ds2 * ai2;
	    ai2 = dc2 * ai2 + ds2 * ar2;
	    ar2 = ar2h;
	    i__3 = *idl1;
	    for (ik = 1; ik <= i__3; ++ik) {
		ch2[ik + l * ch2_dim1] += ar2 * c2[ik + j * c2_dim1];
		ch2[ik + lc * ch2_dim1] += ai2 * c2[ik + jc * c2_dim1];
/* L125: */
	    }
/* L126: */
	}
/* L127: */
	}
	i__1 = ipph;
	for (j = 2; j <= i__1; ++j) {
	i__2 = *idl1;
	for (ik = 1; ik <= i__2; ++ik) {
	    ch2[ik + ch2_dim1] += c2[ik + j * c2_dim1];
/* L128: */
	}
/* L129: */
	}

	if (*ido < *l1) {
	goto L132;
	}
	i__1 = *l1;
	for (k = 1; k <= i__1; ++k) {
	i__2 = *ido;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    cc[i__ + (k * cc_dim2 + 1) * cc_dim1] = ch[i__ + (k + ch_dim2) * ch_dim1];
/* L130: */
	}
/* L131: */
	}
	goto L135;
L132:
	i__1 = *ido;
	for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    cc[i__ + (k * cc_dim2 + 1) * cc_dim1] = ch[i__ + (k + ch_dim2) * ch_dim1];
/* L133: */
	}
/* L134: */
	}
L135:
	i__1 = ipph;
	for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	j2 = j + j;
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    cc[*ido + (j2 - 2 + k * cc_dim2) * cc_dim1] = ch[(k + j * ch_dim2)* ch_dim1 + 1];
	    cc[(j2 - 1 + k * cc_dim2) * cc_dim1 + 1] = ch[(k + jc * ch_dim2) *ch_dim1 + 1];
/* L136: */
	}
/* L137: */
	}
	if (*ido == 1) {
	return;
	}
	if (nbd < *l1) {
	goto L141;
	}
	i__1 = ipph;
	for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	j2 = j + j;
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    i__3 = *ido;
	    for (i__ = 3; i__ <= i__3; i__ += 2) {
		ic = idp2 - i__;
		cc[i__ - 1 + (j2 - 1 + k * cc_dim2) * cc_dim1]
			= ch[i__ - 1 + (k + j * ch_dim2) * ch_dim1] + ch[i__ - 1 + (k + jc * ch_dim2) * ch_dim1];
		cc[ic - 1 + (j2 - 2 + k * cc_dim2) * cc_dim1]
			= ch[i__ - 1 + (k + j * ch_dim2) * ch_dim1] - ch[i__ - 1 + (k + jc * ch_dim2) * ch_dim1];
		cc[i__ + (j2 - 1 + k * cc_dim2) * cc_dim1]
			= ch[i__ + (k + j *ch_dim2) * ch_dim1] + ch[i__ + (k + jc * ch_dim2) * ch_dim1];
		cc[ic + (j2 - 2 + k * cc_dim2) * cc_dim1]
			= ch[i__ + (k + jc *ch_dim2) * ch_dim1] - ch[i__ + (k + j * ch_dim2) * ch_dim1];
/* L138: */
	    }
/* L139: */
	}
/* L140: */
	}
	return;
L141:
	i__1 = ipph;
	for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	j2 = j + j;
	i__2 = *ido;
	for (i__ = 3; i__ <= i__2; i__ += 2) {
	    ic = idp2 - i__;
	    i__3 = *l1;
	    for (k = 1; k <= i__3; ++k) {
		cc[i__ - 1 + (j2 - 1 + k * cc_dim2) * cc_dim1]
			= ch[i__ - 1 + (k + j * ch_dim2) * ch_dim1] + ch[i__ - 1 + (k + jc * ch_dim2) * ch_dim1];
		cc[ic - 1 + (j2 - 2 + k * cc_dim2) * cc_dim1]
			= ch[i__ - 1 + (k + j * ch_dim2) * ch_dim1] - ch[i__ - 1 + (k + jc * ch_dim2) * ch_dim1];
		cc[i__ + (j2 - 1 + k * cc_dim2) * cc_dim1]
			= ch[i__ + (k + j *ch_dim2) * ch_dim1] + ch[i__ + (k + jc * ch_dim2) * ch_dim1];
		cc[ic + (j2 - 2 + k * cc_dim2) * cc_dim1]
			= ch[i__ + (k + jc *ch_dim2) * ch_dim1] - ch[i__ + (k + j * ch_dim2) * ch_dim1];
/* L142: */
	    }
/* L143: */
	}
/* L144: */
	}
	return;
} /* radfg_ */

/* Subroutine */ void FUNC(rfftb)(int *n, real_t *r__, real_t *wsave, int *ifac)
{
	/* Parameter adjustments */
	--ifac;
	--wsave;
	--r__;

	/* Function Body */
	if (*n == 1) {
	return;
	}
	s_rfftb1(n, &r__[1], &wsave[1], &wsave[*n + 1], &ifac[1]);
	return;
} /* rfftb_ */

/* Subroutine */ static void s_rfftb1(int *n, real_t *c__, real_t *ch, 
	real_t *wa, int *ifac)
{
	/* System generated locals */
	int i__1;

	/* Local variables */
	int i__, k1, l1, l2, na, nf, ip, iw, ix2, ix3, ix4, ido, idl1;

	/* Parameter adjustments */
	--ifac;
	--wa;
	--ch;
	--c__;

	/* Function Body */
	nf = ifac[2];
	na = 0;
	l1 = 1;
	iw = 1;
	i__1 = nf;
	for (k1 = 1; k1 <= i__1; ++k1) {
	ip = ifac[k1 + 2];
	l2 = ip * l1;
	ido = *n / l2;
	idl1 = ido * l1;
	if (ip != 4) {
	    goto L103;
	}
	ix2 = iw + ido;
	ix3 = ix2 + ido;
	if (na != 0) {
	    goto L101;
	}
	s_radb4(&ido, &l1, &c__[1], &ch[1], &wa[iw], &wa[ix2], &wa[ix3]);
	goto L102;
L101:
	s_radb4(&ido, &l1, &ch[1], &c__[1], &wa[iw], &wa[ix2], &wa[ix3]);
L102:
	na = 1 - na;
	goto L115;
L103:
	if (ip != 2) {
	    goto L106;
	}
	if (na != 0) {
	    goto L104;
	}
	s_radb2(&ido, &l1, &c__[1], &ch[1], &wa[iw]);
	goto L105;
L104:
	s_radb2(&ido, &l1, &ch[1], &c__[1], &wa[iw]);
L105:
	na = 1 - na;
	goto L115;
L106:
	if (ip != 3) {
	    goto L109;
	}
	ix2 = iw + ido;
	if (na != 0) {
	    goto L107;
	}
	s_radb3(&ido, &l1, &c__[1], &ch[1], &wa[iw], &wa[ix2]);
	goto L108;
L107:
	s_radb3(&ido, &l1, &ch[1], &c__[1], &wa[iw], &wa[ix2]);
L108:
	na = 1 - na;
	goto L115;
L109:
	if (ip != 5) {
	    goto L112;
	}
	ix2 = iw + ido;
	ix3 = ix2 + ido;
	ix4 = ix3 + ido;
	if (na != 0) {
	    goto L110;
	}
	s_radb5(&ido, &l1, &c__[1], &ch[1], &wa[iw], &wa[ix2], &wa[ix3], &wa[ix4]);
	goto L111;
L110:
	s_radb5(&ido, &l1, &ch[1], &c__[1], &wa[iw], &wa[ix2], &wa[ix3], &wa[ix4]);
L111:
	na = 1 - na;
	goto L115;
L112:
	if (na != 0) {
	    goto L113;
	}
	s_radbg(&ido, &ip, &l1, &idl1, &c__[1], &c__[1], &c__[1], &ch[1], &ch[1], &wa[iw]);
	goto L114;
L113:
	s_radbg(&ido, &ip, &l1, &idl1, &ch[1], &ch[1], &ch[1], &c__[1], &c__[1], &wa[iw]);
L114:
	if (ido == 1) {
	    na = 1 - na;
	}
L115:
	l1 = l2;
	iw += (ip - 1) * ido;
/* L116: */
	}
	if (na == 0) {
	return;
	}
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	c__[i__] = ch[i__];
/* L117: */
	}
	return;
} /* rfftb1_ */

/* Subroutine */ void FUNC(rfftf)(int *n, real_t *r__, real_t *wsave, 
	int *ifac)
{
	/* Parameter adjustments */
	--ifac;
	--wsave;
	--r__;

	/* Function Body */
	if (*n == 1) {
	return;
	}
	s_rfftf1(n, &r__[1], &wsave[1], &wsave[*n + 1], &ifac[1]);
	return;
} /* rfftf_ */

/* Subroutine */ static void s_rfftf1(int *n, real_t *c__, real_t *ch, 
	real_t *wa, int *ifac)
{
	/* System generated locals */
	int i__1;

	/* Local variables */
	int i__, k1, l1, l2, na, kh, nf, ip, iw, ix2, ix3, ix4, ido, idl1;

	/* Parameter adjustments */
	--ifac;
	--wa;
	--ch;
	--c__;

	/* Function Body */
	nf = ifac[2];
	na = 1;
	l2 = *n;
	iw = *n;
	i__1 = nf;
	for (k1 = 1; k1 <= i__1; ++k1) {
	kh = nf - k1;
	ip = ifac[kh + 3];
	l1 = l2 / ip;
	ido = *n / l2;
	idl1 = ido * l1;
	iw -= (ip - 1) * ido;
	na = 1 - na;
	if (ip != 4) {
	    goto L102;
	}
	ix2 = iw + ido;
	ix3 = ix2 + ido;
	if (na != 0) {
	    goto L101;
	}
	s_radf4(&ido, &l1, &c__[1], &ch[1], &wa[iw], &wa[ix2], &wa[ix3]);
	goto L110;
L101:
	s_radf4(&ido, &l1, &ch[1], &c__[1], &wa[iw], &wa[ix2], &wa[ix3]);
	goto L110;
L102:
	if (ip != 2) {
	    goto L104;
	}
	if (na != 0) {
	    goto L103;
	}
	s_radf2(&ido, &l1, &c__[1], &ch[1], &wa[iw]);
	goto L110;
L103:
	s_radf2(&ido, &l1, &ch[1], &c__[1], &wa[iw]);
	goto L110;
L104:
	if (ip != 3) {
	    goto L106;
	}
	ix2 = iw + ido;
	if (na != 0) {
	    goto L105;
	}
	s_radf3(&ido, &l1, &c__[1], &ch[1], &wa[iw], &wa[ix2]);
	goto L110;
L105:
	s_radf3(&ido, &l1, &ch[1], &c__[1], &wa[iw], &wa[ix2]);
	goto L110;
L106:
	if (ip != 5) {
	    goto L108;
	}
	ix2 = iw + ido;
	ix3 = ix2 + ido;
	ix4 = ix3 + ido;
	if (na != 0) {
	    goto L107;
	}
	s_radf5(&ido, &l1, &c__[1], &ch[1], &wa[iw], &wa[ix2], &wa[ix3], &wa[ix4]);
	goto L110;
L107:
	s_radf5(&ido, &l1, &ch[1], &c__[1], &wa[iw], &wa[ix2], &wa[ix3], &wa[ix4]);
	goto L110;
L108:
	if (ido == 1) {
	    na = 1 - na;
	}
	if (na != 0) {
	    goto L109;
	}
	s_radfg(&ido, &ip, &l1, &idl1, &c__[1], &c__[1], &c__[1], &ch[1], &ch[1], &wa[iw]);
	na = 1;
	goto L110;
L109:
	s_radfg(&ido, &ip, &l1, &idl1, &ch[1], &ch[1], &ch[1], &c__[1], &c__[1], &wa[iw]);
	na = 0;
L110:
	l2 = l1;
/* L111: */
	}
	if (na == 1) {
	return;
	}
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	c__[i__] = ch[i__];
/* L112: */
	}
	return;
} /* rfftf1_ */

/* Subroutine */ void FUNC(rffti)(int *n, real_t *wsave, int *ifac)
{
	/* Parameter adjustments */
	--ifac;
	--wsave;

	/* Function Body */
	if (*n == 1) {
	return;
	}
	s_rffti1(n, &wsave[*n + 1], &ifac[1]);
	return;
} /* rffti_ */

/* Subroutine */ static void s_rffti1(int *n, real_t *wa, int *ifac)
{
	/* Initialized data */

	static int ntryh[4] = { 4,2,3,5 };

	/* System generated locals */
	int i__1, i__2, i__3;

	/* Local variables */
	int i__, j, k1, l1, l2, ib;
	real_t fi;
	int ld, ii, nf, ip, nl, is, nq, nr;
	real_t arg;
	int ido, ipm;
	real_t tpi;
	int nfm1;
	real_t argh;
	int ntry=0;
	real_t argld;

	/* Parameter adjustments */
	--ifac;
	--wa;

	/* Function Body */
	nl = *n;
	nf = 0;
	j = 0;
L101:
	++j;
	if (j - 4 <= 0) {
	goto L102;
	} else {
	goto L103;
	}
L102:
	ntry = ntryh[j - 1];
	goto L104;
L103:
	ntry += 2;
L104:
	nq = nl / ntry;
	nr = nl - ntry * nq;
	if (nr != 0) {
	goto L101;
	} else {
	goto L105;
	}
L105:
	++nf;
	ifac[nf + 2] = ntry;
	nl = nq;
	if (ntry != 2) {
	goto L107;
	}
	if (nf == 1) {
	goto L107;
	}
	i__1 = nf;
	for (i__ = 2; i__ <= i__1; ++i__) {
	ib = nf - i__ + 2;
	ifac[ib + 2] = ifac[ib + 1];
/* L106: */
	}
	ifac[3] = 2;
L107:
	if (nl != 1) {
	goto L104;
	}
	ifac[1] = *n;
	ifac[2] = nf;
	tpi = POST(6.283185307179586476925286766559005768394338798750211619498891846);
	argh = tpi / (real_t) (*n);
	is = 0;
	nfm1 = nf - 1;
	l1 = 1;
	if (nfm1 == 0) {
	return;
	}
	i__1 = nfm1;
	for (k1 = 1; k1 <= i__1; ++k1) {
	ip = ifac[k1 + 2];
	ld = 0;
	l2 = l1 * ip;
	ido = *n / l2;
	ipm = ip - 1;
	i__2 = ipm;
	for (j = 1; j <= i__2; ++j) {
	    ld += l1;
	    i__ = is;
	    argld = (real_t) ld * argh;
	    fi = POST(0.0);
	    i__3 = ido;
	    for (ii = 3; ii <= i__3; ii += 2) {
		i__ += 2;
		fi += POST(1.0);
		arg = fi * argld;
		wa[i__ - 1] = POST(cos)(arg);
		wa[i__] = POST(sin)(arg);
/* L108: */
	    }
	    is += ido;
/* L109: */
	}
	l1 = l2;
/* L110: */
	}
	return;
} /* s_rffti1 */

/* Subroutine */ void FUNC(sinqb)(int *n, real_t *x, real_t *wsave, 
	int *ifac)
{
	/* System generated locals */
	int i__1;

	/* Local variables */
	int k, kc, ns2;
	real_t xhold;

	/* Parameter adjustments */
	--ifac;
	--wsave;
	--x;

	/* Function Body */
	if (*n > 1) {
	goto L101;
	}
	x[1] *= POST(4.0);
	return;
L101:
	ns2 = *n / 2;
	i__1 = *n;
	for (k = 2; k <= i__1; k += 2) {
	x[k] = -x[k];
/* L102: */
	}
	FUNC(cosqb)(n, &x[1], &wsave[1], &ifac[1]);
	i__1 = ns2;
	for (k = 1; k <= i__1; ++k) {
	kc = *n - k;
	xhold = x[k];
	x[k] = x[kc + 1];
	x[kc + 1] = xhold;
/* L103: */
	}
	return;
} /* sinqb_ */

/* Subroutine */ void FUNC(sinqf)(int *n, real_t *x, real_t *wsave, 
	int *ifac)
{
	/* System generated locals */
	int i__1;

	/* Local variables */
	int k, kc, ns2;
	real_t xhold;

	/* Parameter adjustments */
	--ifac;
	--wsave;
	--x;

	/* Function Body */
	if (*n == 1) {
	return;
	}
	ns2 = *n / 2;
	i__1 = ns2;
	for (k = 1; k <= i__1; ++k) {
	kc = *n - k;
	xhold = x[k];
	x[k] = x[kc + 1];
	x[kc + 1] = xhold;
/* L101: */
	}
	FUNC(cosqf)(n, &x[1], &wsave[1], &ifac[1]);
	i__1 = *n;
	for (k = 2; k <= i__1; k += 2) {
	x[k] = -x[k];
/* L102: */
	}
	return;
} /* sinqf_ */

/* Subroutine */ void FUNC(sinqi)(int *n, real_t *wsave, int *ifac)
{
	/* Parameter adjustments */
	--ifac;
	--wsave;

	/* Function Body */
	FUNC(cosqi)(n, &wsave[1], &ifac[1]);
	return;
} /* sinqi_ */

/* Subroutine */ void FUNC(sint)(int *n, real_t *x, real_t *wsave, 
	int *ifac)
{
	int np1, iw1, iw2;

	/* Parameter adjustments */
	--ifac;
	--wsave;
	--x;

	/* Function Body */
	np1 = *n + 1;
	iw1 = *n / 2 + 1;
	iw2 = iw1 + np1;
	s_sint1(n, &x[1], &wsave[1], &wsave[iw1], &wsave[iw2], &ifac[1]);
	return;
} /* sint_ */

/* Subroutine */ static void s_sint1(int *n, real_t *war, real_t *was, 
	real_t *xh, real_t *x, int *ifac)
{
	/* Initialized data */

	static real_t sqrt3 = POST(1.732050807568877293527446341505872366942805253803806280558069795);

	/* System generated locals */
	int i__1;

	/* Local variables */
	int i__, k;
	real_t t1, t2;
	int kc, np1, ns2, modn;
	real_t xhold;

	/* Parameter adjustments */
	--ifac;
	--x;
	--xh;
	--was;
	--war;

	/* Function Body */
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	xh[i__] = war[i__];
	war[i__] = x[i__];
/* L100: */
	}
	if ((i__1 = *n - 2) < 0) {
	goto L101;
	} else if (i__1 == 0) {
	goto L102;
	} else {
	goto L103;
	}
L101:
	xh[1] += xh[1];
	goto L106;
L102:
	xhold = sqrt3 * (xh[1] + xh[2]);
	xh[2] = sqrt3 * (xh[1] - xh[2]);
	xh[1] = xhold;
	goto L106;
L103:
	np1 = *n + 1;
	ns2 = *n / 2;
	x[1] = POST(0.0);
	i__1 = ns2;
	for (k = 1; k <= i__1; ++k) {
	kc = np1 - k;
	t1 = xh[k] - xh[kc];
	t2 = was[k] * (xh[k] + xh[kc]);
	x[k + 1] = t1 + t2;
	x[kc + 1] = t2 - t1;
/* L104: */
	}
	modn = *n % 2;
	if (modn != 0) {
	x[ns2 + 2] = xh[ns2 + 1] * POST(4.0);
	}
	s_rfftf1(&np1, &x[1], &xh[1], &war[1], &ifac[1]);
	xh[1] = x[1] * POST(0.5);
	i__1 = *n;
	for (i__ = 3; i__ <= i__1; i__ += 2) {
	xh[i__ - 1] = -x[i__];
	xh[i__] = xh[i__ - 2] + x[i__ - 1];
/* L105: */
	}
	if (modn != 0) {
	goto L106;
	}
	xh[*n] = -x[*n + 1];
L106:
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	x[i__] = war[i__];
	war[i__] = xh[i__];
/* L107: */
	}
	return;
} /* sint1_ */

/* Subroutine */ void FUNC(sinti)(int *n, real_t *wsave, int *ifac)
{
	/* Initialized data */

	static real_t pi = POST(3.141592653589793238462643383279502884197169399375158209749445923);

	/* System generated locals */
	int i__1;

	/* Local variables */
	int k;
	real_t dt;
	int np1, ns2;

	/* Parameter adjustments */
	--ifac;
	--wsave;

	/* Function Body */
	if (*n <= 1) {
	return;
	}
	ns2 = *n / 2;
	np1 = *n + 1;
	dt = pi / (real_t) np1;
	i__1 = ns2;
	for (k = 1; k <= i__1; ++k) {
	wsave[k] = POST(sin)(k * dt) * POST(2.0);
/* L101: */
	}
	FUNC(rffti)(&np1, &wsave[ns2 + 1], &ifac[1]);
	return;
} /* sinti_ */
