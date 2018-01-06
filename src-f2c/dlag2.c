#line 1 "dlag2.f"
/* dlag2.f -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

#line 1 "dlag2.f"
/* > \brief \b DLAG2 computes the eigenvalues of a 2-by-2 generalized eigenvalue problem, with scaling as nece
ssary to avoid over-/underflow. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLAG2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlag2.f
"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlag2.f
"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlag2.f
"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLAG2( A, LDA, B, LDB, SAFMIN, SCALE1, SCALE2, WR1, */
/*                         WR2, WI ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            LDA, LDB */
/*       DOUBLE PRECISION   SAFMIN, SCALE1, SCALE2, WI, WR1, WR2 */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLAG2 computes the eigenvalues of a 2 x 2 generalized eigenvalue */
/* > problem  A - w B, with scaling as necessary to avoid over-/underflow. */
/* > */
/* > The scaling factor "s" results in a modified eigenvalue equation */
/* > */
/* >     s A - w B */
/* > */
/* > where  s  is a non-negative scaling factor chosen so that  w,  w B, */
/* > and  s A  do not overflow and, if possible, do not underflow, either. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension (LDA, 2) */
/* >          On entry, the 2 x 2 matrix A.  It is assumed that its 1-norm */
/* >          is less than 1/SAFMIN.  Entries less than */
/* >          sqrt(SAFMIN)*norm(A) are subject to being treated as zero. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= 2. */
/* > \endverbatim */
/* > */
/* > \param[in] B */
/* > \verbatim */
/* >          B is DOUBLE PRECISION array, dimension (LDB, 2) */
/* >          On entry, the 2 x 2 upper triangular matrix B.  It is */
/* >          assumed that the one-norm of B is less than 1/SAFMIN.  The */
/* >          diagonals should be at least sqrt(SAFMIN) times the largest */
/* >          element of B (in absolute value); if a diagonal is smaller */
/* >          than that, then  +/- sqrt(SAFMIN) will be used instead of */
/* >          that diagonal. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >          The leading dimension of the array B.  LDB >= 2. */
/* > \endverbatim */
/* > */
/* > \param[in] SAFMIN */
/* > \verbatim */
/* >          SAFMIN is DOUBLE PRECISION */
/* >          The smallest positive number s.t. 1/SAFMIN does not */
/* >          overflow.  (This should always be DLAMCH('S') -- it is an */
/* >          argument in order to avoid having to call DLAMCH frequently.) */
/* > \endverbatim */
/* > */
/* > \param[out] SCALE1 */
/* > \verbatim */
/* >          SCALE1 is DOUBLE PRECISION */
/* >          A scaling factor used to avoid over-/underflow in the */
/* >          eigenvalue equation which defines the first eigenvalue.  If */
/* >          the eigenvalues are complex, then the eigenvalues are */
/* >          ( WR1  +/-  WI i ) / SCALE1  (which may lie outside the */
/* >          exponent range of the machine), SCALE1=SCALE2, and SCALE1 */
/* >          will always be positive.  If the eigenvalues are real, then */
/* >          the first (real) eigenvalue is  WR1 / SCALE1 , but this may */
/* >          overflow or underflow, and in fact, SCALE1 may be zero or */
/* >          less than the underflow threshold if the exact eigenvalue */
/* >          is sufficiently large. */
/* > \endverbatim */
/* > */
/* > \param[out] SCALE2 */
/* > \verbatim */
/* >          SCALE2 is DOUBLE PRECISION */
/* >          A scaling factor used to avoid over-/underflow in the */
/* >          eigenvalue equation which defines the second eigenvalue.  If */
/* >          the eigenvalues are complex, then SCALE2=SCALE1.  If the */
/* >          eigenvalues are real, then the second (real) eigenvalue is */
/* >          WR2 / SCALE2 , but this may overflow or underflow, and in */
/* >          fact, SCALE2 may be zero or less than the underflow */
/* >          threshold if the exact eigenvalue is sufficiently large. */
/* > \endverbatim */
/* > */
/* > \param[out] WR1 */
/* > \verbatim */
/* >          WR1 is DOUBLE PRECISION */
/* >          If the eigenvalue is real, then WR1 is SCALE1 times the */
/* >          eigenvalue closest to the (2,2) element of A B**(-1).  If the */
/* >          eigenvalue is complex, then WR1=WR2 is SCALE1 times the real */
/* >          part of the eigenvalues. */
/* > \endverbatim */
/* > */
/* > \param[out] WR2 */
/* > \verbatim */
/* >          WR2 is DOUBLE PRECISION */
/* >          If the eigenvalue is real, then WR2 is SCALE2 times the */
/* >          other eigenvalue.  If the eigenvalue is complex, then */
/* >          WR1=WR2 is SCALE1 times the real part of the eigenvalues. */
/* > \endverbatim */
/* > */
/* > \param[out] WI */
/* > \verbatim */
/* >          WI is DOUBLE PRECISION */
/* >          If the eigenvalue is real, then WI is zero.  If the */
/* >          eigenvalue is complex, then WI is SCALE1 times the imaginary */
/* >          part of the eigenvalues.  WI will always be non-negative. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date June 2016 */

/* > \ingroup doubleOTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int dlag2_(doublereal *a, integer *lda, doublereal *b, 
	integer *ldb, doublereal *safmin, doublereal *scale1, doublereal *
	scale2, doublereal *wr1, doublereal *wr2, doublereal *wi)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6;

    /* Builtin functions */
    double sqrt(doublereal), d_sign(doublereal *, doublereal *);

    /* Local variables */
    static doublereal r__, c1, c2, c3, c4, c5, s1, s2, a11, a12, a21, a22, 
	    b11, b12, b22, pp, qq, ss, as11, as12, as22, sum, abi22, diff, 
	    bmin, wbig, wabs, wdet, binv11, binv22, discr, anorm, bnorm, 
	    bsize, shift, rtmin, rtmax, wsize, ascale, bscale, wscale, safmax,
	     wsmall;


/*  -- LAPACK auxiliary routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     June 2016 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 195 "dlag2.f"
    /* Parameter adjustments */
#line 195 "dlag2.f"
    a_dim1 = *lda;
#line 195 "dlag2.f"
    a_offset = 1 + a_dim1;
#line 195 "dlag2.f"
    a -= a_offset;
#line 195 "dlag2.f"
    b_dim1 = *ldb;
#line 195 "dlag2.f"
    b_offset = 1 + b_dim1;
#line 195 "dlag2.f"
    b -= b_offset;
#line 195 "dlag2.f"

#line 195 "dlag2.f"
    /* Function Body */
#line 195 "dlag2.f"
    rtmin = sqrt(*safmin);
#line 196 "dlag2.f"
    rtmax = 1. / rtmin;
#line 197 "dlag2.f"
    safmax = 1. / *safmin;

/*     Scale A */

/* Computing MAX */
#line 201 "dlag2.f"
    d__5 = (d__1 = a[a_dim1 + 1], abs(d__1)) + (d__2 = a[a_dim1 + 2], abs(
	    d__2)), d__6 = (d__3 = a[(a_dim1 << 1) + 1], abs(d__3)) + (d__4 = 
	    a[(a_dim1 << 1) + 2], abs(d__4)), d__5 = max(d__5,d__6);
#line 201 "dlag2.f"
    anorm = max(d__5,*safmin);
#line 203 "dlag2.f"
    ascale = 1. / anorm;
#line 204 "dlag2.f"
    a11 = ascale * a[a_dim1 + 1];
#line 205 "dlag2.f"
    a21 = ascale * a[a_dim1 + 2];
#line 206 "dlag2.f"
    a12 = ascale * a[(a_dim1 << 1) + 1];
#line 207 "dlag2.f"
    a22 = ascale * a[(a_dim1 << 1) + 2];

/*     Perturb B if necessary to insure non-singularity */

#line 211 "dlag2.f"
    b11 = b[b_dim1 + 1];
#line 212 "dlag2.f"
    b12 = b[(b_dim1 << 1) + 1];
#line 213 "dlag2.f"
    b22 = b[(b_dim1 << 1) + 2];
/* Computing MAX */
#line 214 "dlag2.f"
    d__1 = abs(b11), d__2 = abs(b12), d__1 = max(d__1,d__2), d__2 = abs(b22), 
	    d__1 = max(d__1,d__2);
#line 214 "dlag2.f"
    bmin = rtmin * max(d__1,rtmin);
#line 215 "dlag2.f"
    if (abs(b11) < bmin) {
#line 215 "dlag2.f"
	b11 = d_sign(&bmin, &b11);
#line 215 "dlag2.f"
    }
#line 217 "dlag2.f"
    if (abs(b22) < bmin) {
#line 217 "dlag2.f"
	b22 = d_sign(&bmin, &b22);
#line 217 "dlag2.f"
    }

/*     Scale B */

/* Computing MAX */
#line 222 "dlag2.f"
    d__1 = abs(b11), d__2 = abs(b12) + abs(b22), d__1 = max(d__1,d__2);
#line 222 "dlag2.f"
    bnorm = max(d__1,*safmin);
/* Computing MAX */
#line 223 "dlag2.f"
    d__1 = abs(b11), d__2 = abs(b22);
#line 223 "dlag2.f"
    bsize = max(d__1,d__2);
#line 224 "dlag2.f"
    bscale = 1. / bsize;
#line 225 "dlag2.f"
    b11 *= bscale;
#line 226 "dlag2.f"
    b12 *= bscale;
#line 227 "dlag2.f"
    b22 *= bscale;

/*     Compute larger eigenvalue by method described by C. van Loan */

/*     ( AS is A shifted by -SHIFT*B ) */

#line 233 "dlag2.f"
    binv11 = 1. / b11;
#line 234 "dlag2.f"
    binv22 = 1. / b22;
#line 235 "dlag2.f"
    s1 = a11 * binv11;
#line 236 "dlag2.f"
    s2 = a22 * binv22;
#line 237 "dlag2.f"
    if (abs(s1) <= abs(s2)) {
#line 238 "dlag2.f"
	as12 = a12 - s1 * b12;
#line 239 "dlag2.f"
	as22 = a22 - s1 * b22;
#line 240 "dlag2.f"
	ss = a21 * (binv11 * binv22);
#line 241 "dlag2.f"
	abi22 = as22 * binv22 - ss * b12;
#line 242 "dlag2.f"
	pp = abi22 * .5;
#line 243 "dlag2.f"
	shift = s1;
#line 244 "dlag2.f"
    } else {
#line 245 "dlag2.f"
	as12 = a12 - s2 * b12;
#line 246 "dlag2.f"
	as11 = a11 - s2 * b11;
#line 247 "dlag2.f"
	ss = a21 * (binv11 * binv22);
#line 248 "dlag2.f"
	abi22 = -ss * b12;
#line 249 "dlag2.f"
	pp = (as11 * binv11 + abi22) * .5;
#line 250 "dlag2.f"
	shift = s2;
#line 251 "dlag2.f"
    }
#line 252 "dlag2.f"
    qq = ss * as12;
#line 253 "dlag2.f"
    if ((d__1 = pp * rtmin, abs(d__1)) >= 1.) {
/* Computing 2nd power */
#line 254 "dlag2.f"
	d__1 = rtmin * pp;
#line 254 "dlag2.f"
	discr = d__1 * d__1 + qq * *safmin;
#line 255 "dlag2.f"
	r__ = sqrt((abs(discr))) * rtmax;
#line 256 "dlag2.f"
    } else {
/* Computing 2nd power */
#line 257 "dlag2.f"
	d__1 = pp;
#line 257 "dlag2.f"
	if (d__1 * d__1 + abs(qq) <= *safmin) {
/* Computing 2nd power */
#line 258 "dlag2.f"
	    d__1 = rtmax * pp;
#line 258 "dlag2.f"
	    discr = d__1 * d__1 + qq * safmax;
#line 259 "dlag2.f"
	    r__ = sqrt((abs(discr))) * rtmin;
#line 260 "dlag2.f"
	} else {
/* Computing 2nd power */
#line 261 "dlag2.f"
	    d__1 = pp;
#line 261 "dlag2.f"
	    discr = d__1 * d__1 + qq;
#line 262 "dlag2.f"
	    r__ = sqrt((abs(discr)));
#line 263 "dlag2.f"
	}
#line 264 "dlag2.f"
    }

/*     Note: the test of R in the following IF is to cover the case when */
/*           DISCR is small and negative and is flushed to zero during */
/*           the calculation of R.  On machines which have a consistent */
/*           flush-to-zero threshold and handle numbers above that */
/*           threshold correctly, it would not be necessary. */

#line 272 "dlag2.f"
    if (discr >= 0. || r__ == 0.) {
#line 273 "dlag2.f"
	sum = pp + d_sign(&r__, &pp);
#line 274 "dlag2.f"
	diff = pp - d_sign(&r__, &pp);
#line 275 "dlag2.f"
	wbig = shift + sum;

/*        Compute smaller eigenvalue */

#line 279 "dlag2.f"
	wsmall = shift + diff;
/* Computing MAX */
#line 280 "dlag2.f"
	d__1 = abs(wsmall);
#line 280 "dlag2.f"
	if (abs(wbig) * .5 > max(d__1,*safmin)) {
#line 281 "dlag2.f"
	    wdet = (a11 * a22 - a12 * a21) * (binv11 * binv22);
#line 282 "dlag2.f"
	    wsmall = wdet / wbig;
#line 283 "dlag2.f"
	}

/*        Choose (real) eigenvalue closest to 2,2 element of A*B**(-1) */
/*        for WR1. */

#line 288 "dlag2.f"
	if (pp > abi22) {
#line 289 "dlag2.f"
	    *wr1 = min(wbig,wsmall);
#line 290 "dlag2.f"
	    *wr2 = max(wbig,wsmall);
#line 291 "dlag2.f"
	} else {
#line 292 "dlag2.f"
	    *wr1 = max(wbig,wsmall);
#line 293 "dlag2.f"
	    *wr2 = min(wbig,wsmall);
#line 294 "dlag2.f"
	}
#line 295 "dlag2.f"
	*wi = 0.;
#line 296 "dlag2.f"
    } else {

/*        Complex eigenvalues */

#line 300 "dlag2.f"
	*wr1 = shift + pp;
#line 301 "dlag2.f"
	*wr2 = *wr1;
#line 302 "dlag2.f"
	*wi = r__;
#line 303 "dlag2.f"
    }

/*     Further scaling to avoid underflow and overflow in computing */
/*     SCALE1 and overflow in computing w*B. */

/*     This scale factor (WSCALE) is bounded from above using C1 and C2, */
/*     and from below using C3 and C4. */
/*        C1 implements the condition  s A  must never overflow. */
/*        C2 implements the condition  w B  must never overflow. */
/*        C3, with C2, */
/*           implement the condition that s A - w B must never overflow. */
/*        C4 implements the condition  s    should not underflow. */
/*        C5 implements the condition  max(s,|w|) should be at least 2. */

#line 317 "dlag2.f"
    c1 = bsize * (*safmin * max(1.,ascale));
#line 318 "dlag2.f"
    c2 = *safmin * max(1.,bnorm);
#line 319 "dlag2.f"
    c3 = bsize * *safmin;
#line 320 "dlag2.f"
    if (ascale <= 1. && bsize <= 1.) {
/* Computing MIN */
#line 321 "dlag2.f"
	d__1 = 1., d__2 = ascale / *safmin * bsize;
#line 321 "dlag2.f"
	c4 = min(d__1,d__2);
#line 322 "dlag2.f"
    } else {
#line 323 "dlag2.f"
	c4 = 1.;
#line 324 "dlag2.f"
    }
#line 325 "dlag2.f"
    if (ascale <= 1. || bsize <= 1.) {
/* Computing MIN */
#line 326 "dlag2.f"
	d__1 = 1., d__2 = ascale * bsize;
#line 326 "dlag2.f"
	c5 = min(d__1,d__2);
#line 327 "dlag2.f"
    } else {
#line 328 "dlag2.f"
	c5 = 1.;
#line 329 "dlag2.f"
    }

/*     Scale first eigenvalue */

#line 333 "dlag2.f"
    wabs = abs(*wr1) + abs(*wi);
/* Computing MAX */
/* Computing MIN */
#line 334 "dlag2.f"
    d__3 = c4, d__4 = max(wabs,c5) * .5;
#line 334 "dlag2.f"
    d__1 = max(*safmin,c1), d__2 = (wabs * c2 + c3) * 1.0000100000000001, 
	    d__1 = max(d__1,d__2), d__2 = min(d__3,d__4);
#line 334 "dlag2.f"
    wsize = max(d__1,d__2);
#line 336 "dlag2.f"
    if (wsize != 1.) {
#line 337 "dlag2.f"
	wscale = 1. / wsize;
#line 338 "dlag2.f"
	if (wsize > 1.) {
#line 339 "dlag2.f"
	    *scale1 = max(ascale,bsize) * wscale * min(ascale,bsize);
#line 341 "dlag2.f"
	} else {
#line 342 "dlag2.f"
	    *scale1 = min(ascale,bsize) * wscale * max(ascale,bsize);
#line 344 "dlag2.f"
	}
#line 345 "dlag2.f"
	*wr1 *= wscale;
#line 346 "dlag2.f"
	if (*wi != 0.) {
#line 347 "dlag2.f"
	    *wi *= wscale;
#line 348 "dlag2.f"
	    *wr2 = *wr1;
#line 349 "dlag2.f"
	    *scale2 = *scale1;
#line 350 "dlag2.f"
	}
#line 351 "dlag2.f"
    } else {
#line 352 "dlag2.f"
	*scale1 = ascale * bsize;
#line 353 "dlag2.f"
	*scale2 = *scale1;
#line 354 "dlag2.f"
    }

/*     Scale second eigenvalue (if real) */

#line 358 "dlag2.f"
    if (*wi == 0.) {
/* Computing MAX */
/* Computing MIN */
/* Computing MAX */
#line 359 "dlag2.f"
	d__5 = abs(*wr2);
#line 359 "dlag2.f"
	d__3 = c4, d__4 = max(d__5,c5) * .5;
#line 359 "dlag2.f"
	d__1 = max(*safmin,c1), d__2 = (abs(*wr2) * c2 + c3) * 
		1.0000100000000001, d__1 = max(d__1,d__2), d__2 = min(d__3,
		d__4);
#line 359 "dlag2.f"
	wsize = max(d__1,d__2);
#line 361 "dlag2.f"
	if (wsize != 1.) {
#line 362 "dlag2.f"
	    wscale = 1. / wsize;
#line 363 "dlag2.f"
	    if (wsize > 1.) {
#line 364 "dlag2.f"
		*scale2 = max(ascale,bsize) * wscale * min(ascale,bsize);
#line 366 "dlag2.f"
	    } else {
#line 367 "dlag2.f"
		*scale2 = min(ascale,bsize) * wscale * max(ascale,bsize);
#line 369 "dlag2.f"
	    }
#line 370 "dlag2.f"
	    *wr2 *= wscale;
#line 371 "dlag2.f"
	} else {
#line 372 "dlag2.f"
	    *scale2 = ascale * bsize;
#line 373 "dlag2.f"
	}
#line 374 "dlag2.f"
    }

/*     End of DLAG2 */

#line 378 "dlag2.f"
    return 0;
} /* dlag2_ */

