#line 1 "dlasd3.f"
/* dlasd3.f -- translated by f2c (version 20100827).
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

#line 1 "dlasd3.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c__0 = 0;
static doublereal c_b13 = 1.;
static doublereal c_b26 = 0.;

/* > \brief \b DLASD3 finds all square roots of the roots of the secular equation, as defined by the values in
 D and Z, and then updates the singular vectors by matrix multiplication. Used by sbdsdc. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLASD3 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlasd3.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlasd3.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlasd3.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLASD3( NL, NR, SQRE, K, D, Q, LDQ, DSIGMA, U, LDU, U2, */
/*                          LDU2, VT, LDVT, VT2, LDVT2, IDXC, CTOT, Z, */
/*                          INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, K, LDQ, LDU, LDU2, LDVT, LDVT2, NL, NR, */
/*      $                   SQRE */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            CTOT( * ), IDXC( * ) */
/*       DOUBLE PRECISION   D( * ), DSIGMA( * ), Q( LDQ, * ), U( LDU, * ), */
/*      $                   U2( LDU2, * ), VT( LDVT, * ), VT2( LDVT2, * ), */
/*      $                   Z( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLASD3 finds all the square roots of the roots of the secular */
/* > equation, as defined by the values in D and Z.  It makes the */
/* > appropriate calls to DLASD4 and then updates the singular */
/* > vectors by matrix multiplication. */
/* > */
/* > This code makes very mild assumptions about floating point */
/* > arithmetic. It will work on machines with a guard digit in */
/* > add/subtract, or on those binary machines without guard digits */
/* > which subtract like the Cray XMP, Cray YMP, Cray C 90, or Cray 2. */
/* > It could conceivably fail on hexadecimal or decimal machines */
/* > without guard digits, but we know of none. */
/* > */
/* > DLASD3 is called from DLASD1. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] NL */
/* > \verbatim */
/* >          NL is INTEGER */
/* >         The row dimension of the upper block.  NL >= 1. */
/* > \endverbatim */
/* > */
/* > \param[in] NR */
/* > \verbatim */
/* >          NR is INTEGER */
/* >         The row dimension of the lower block.  NR >= 1. */
/* > \endverbatim */
/* > */
/* > \param[in] SQRE */
/* > \verbatim */
/* >          SQRE is INTEGER */
/* >         = 0: the lower block is an NR-by-NR square matrix. */
/* >         = 1: the lower block is an NR-by-(NR+1) rectangular matrix. */
/* > */
/* >         The bidiagonal matrix has N = NL + NR + 1 rows and */
/* >         M = N + SQRE >= N columns. */
/* > \endverbatim */
/* > */
/* > \param[in] K */
/* > \verbatim */
/* >          K is INTEGER */
/* >         The size of the secular equation, 1 =< K = < N. */
/* > \endverbatim */
/* > */
/* > \param[out] D */
/* > \verbatim */
/* >          D is DOUBLE PRECISION array, dimension(K) */
/* >         On exit the square roots of the roots of the secular equation, */
/* >         in ascending order. */
/* > \endverbatim */
/* > */
/* > \param[out] Q */
/* > \verbatim */
/* >          Q is DOUBLE PRECISION array, */
/* >                     dimension at least (LDQ,K). */
/* > \endverbatim */
/* > */
/* > \param[in] LDQ */
/* > \verbatim */
/* >          LDQ is INTEGER */
/* >         The leading dimension of the array Q.  LDQ >= K. */
/* > \endverbatim */
/* > */
/* > \param[in] DSIGMA */
/* > \verbatim */
/* >          DSIGMA is DOUBLE PRECISION array, dimension(K) */
/* >         The first K elements of this array contain the old roots */
/* >         of the deflated updating problem.  These are the poles */
/* >         of the secular equation. */
/* > \endverbatim */
/* > */
/* > \param[out] U */
/* > \verbatim */
/* >          U is DOUBLE PRECISION array, dimension (LDU, N) */
/* >         The last N - K columns of this matrix contain the deflated */
/* >         left singular vectors. */
/* > \endverbatim */
/* > */
/* > \param[in] LDU */
/* > \verbatim */
/* >          LDU is INTEGER */
/* >         The leading dimension of the array U.  LDU >= N. */
/* > \endverbatim */
/* > */
/* > \param[in,out] U2 */
/* > \verbatim */
/* >          U2 is DOUBLE PRECISION array, dimension (LDU2, N) */
/* >         The first K columns of this matrix contain the non-deflated */
/* >         left singular vectors for the split problem. */
/* > \endverbatim */
/* > */
/* > \param[in] LDU2 */
/* > \verbatim */
/* >          LDU2 is INTEGER */
/* >         The leading dimension of the array U2.  LDU2 >= N. */
/* > \endverbatim */
/* > */
/* > \param[out] VT */
/* > \verbatim */
/* >          VT is DOUBLE PRECISION array, dimension (LDVT, M) */
/* >         The last M - K columns of VT**T contain the deflated */
/* >         right singular vectors. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVT */
/* > \verbatim */
/* >          LDVT is INTEGER */
/* >         The leading dimension of the array VT.  LDVT >= N. */
/* > \endverbatim */
/* > */
/* > \param[in,out] VT2 */
/* > \verbatim */
/* >          VT2 is DOUBLE PRECISION array, dimension (LDVT2, N) */
/* >         The first K columns of VT2**T contain the non-deflated */
/* >         right singular vectors for the split problem. */
/* > \endverbatim */
/* > */
/* > \param[in] LDVT2 */
/* > \verbatim */
/* >          LDVT2 is INTEGER */
/* >         The leading dimension of the array VT2.  LDVT2 >= N. */
/* > \endverbatim */
/* > */
/* > \param[in] IDXC */
/* > \verbatim */
/* >          IDXC is INTEGER array, dimension ( N ) */
/* >         The permutation used to arrange the columns of U (and rows of */
/* >         VT) into three groups:  the first group contains non-zero */
/* >         entries only at and above (or before) NL +1; the second */
/* >         contains non-zero entries only at and below (or after) NL+2; */
/* >         and the third is dense. The first column of U and the row of */
/* >         VT are treated separately, however. */
/* > */
/* >         The rows of the singular vectors found by DLASD4 */
/* >         must be likewise permuted before the matrix multiplies can */
/* >         take place. */
/* > \endverbatim */
/* > */
/* > \param[in] CTOT */
/* > \verbatim */
/* >          CTOT is INTEGER array, dimension ( 4 ) */
/* >         A count of the total number of the various types of columns */
/* >         in U (or rows in VT), as described in IDXC. The fourth column */
/* >         type is any column which has been deflated. */
/* > \endverbatim */
/* > */
/* > \param[in] Z */
/* > \verbatim */
/* >          Z is DOUBLE PRECISION array, dimension (K) */
/* >         The first K elements of this array contain the components */
/* >         of the deflation-adjusted updating row vector. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >         = 0:  successful exit. */
/* >         < 0:  if INFO = -i, the i-th argument had an illegal value. */
/* >         > 0:  if INFO = 1, a singular value did not converge */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup auxOTHERauxiliary */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Ming Gu and Huan Ren, Computer Science Division, University of */
/* >     California at Berkeley, USA */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dlasd3_(integer *nl, integer *nr, integer *sqre, integer 
	*k, doublereal *d__, doublereal *q, integer *ldq, doublereal *dsigma, 
	doublereal *u, integer *ldu, doublereal *u2, integer *ldu2, 
	doublereal *vt, integer *ldvt, doublereal *vt2, integer *ldvt2, 
	integer *idxc, integer *ctot, doublereal *z__, integer *info)
{
    /* System generated locals */
    integer q_dim1, q_offset, u_dim1, u_offset, u2_dim1, u2_offset, vt_dim1, 
	    vt_offset, vt2_dim1, vt2_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal), d_sign(doublereal *, doublereal *);

    /* Local variables */
    static integer i__, j, m, n, jc;
    static doublereal rho;
    static integer nlp1, nlp2, nrp1;
    static doublereal temp;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    static integer ctemp;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer ktemp;
    extern doublereal dlamc3_(doublereal *, doublereal *);
    extern /* Subroutine */ int dlasd4_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *), dlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen), dlacpy_(char *, integer *, integer 
	    *, doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);


/*  -- LAPACK auxiliary routine (version 3.4.2) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     September 2012 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 269 "dlasd3.f"
    /* Parameter adjustments */
#line 269 "dlasd3.f"
    --d__;
#line 269 "dlasd3.f"
    q_dim1 = *ldq;
#line 269 "dlasd3.f"
    q_offset = 1 + q_dim1;
#line 269 "dlasd3.f"
    q -= q_offset;
#line 269 "dlasd3.f"
    --dsigma;
#line 269 "dlasd3.f"
    u_dim1 = *ldu;
#line 269 "dlasd3.f"
    u_offset = 1 + u_dim1;
#line 269 "dlasd3.f"
    u -= u_offset;
#line 269 "dlasd3.f"
    u2_dim1 = *ldu2;
#line 269 "dlasd3.f"
    u2_offset = 1 + u2_dim1;
#line 269 "dlasd3.f"
    u2 -= u2_offset;
#line 269 "dlasd3.f"
    vt_dim1 = *ldvt;
#line 269 "dlasd3.f"
    vt_offset = 1 + vt_dim1;
#line 269 "dlasd3.f"
    vt -= vt_offset;
#line 269 "dlasd3.f"
    vt2_dim1 = *ldvt2;
#line 269 "dlasd3.f"
    vt2_offset = 1 + vt2_dim1;
#line 269 "dlasd3.f"
    vt2 -= vt2_offset;
#line 269 "dlasd3.f"
    --idxc;
#line 269 "dlasd3.f"
    --ctot;
#line 269 "dlasd3.f"
    --z__;
#line 269 "dlasd3.f"

#line 269 "dlasd3.f"
    /* Function Body */
#line 269 "dlasd3.f"
    *info = 0;

#line 271 "dlasd3.f"
    if (*nl < 1) {
#line 272 "dlasd3.f"
	*info = -1;
#line 273 "dlasd3.f"
    } else if (*nr < 1) {
#line 274 "dlasd3.f"
	*info = -2;
#line 275 "dlasd3.f"
    } else if (*sqre != 1 && *sqre != 0) {
#line 276 "dlasd3.f"
	*info = -3;
#line 277 "dlasd3.f"
    }

#line 279 "dlasd3.f"
    n = *nl + *nr + 1;
#line 280 "dlasd3.f"
    m = n + *sqre;
#line 281 "dlasd3.f"
    nlp1 = *nl + 1;
#line 282 "dlasd3.f"
    nlp2 = *nl + 2;

#line 284 "dlasd3.f"
    if (*k < 1 || *k > n) {
#line 285 "dlasd3.f"
	*info = -4;
#line 286 "dlasd3.f"
    } else if (*ldq < *k) {
#line 287 "dlasd3.f"
	*info = -7;
#line 288 "dlasd3.f"
    } else if (*ldu < n) {
#line 289 "dlasd3.f"
	*info = -10;
#line 290 "dlasd3.f"
    } else if (*ldu2 < n) {
#line 291 "dlasd3.f"
	*info = -12;
#line 292 "dlasd3.f"
    } else if (*ldvt < m) {
#line 293 "dlasd3.f"
	*info = -14;
#line 294 "dlasd3.f"
    } else if (*ldvt2 < m) {
#line 295 "dlasd3.f"
	*info = -16;
#line 296 "dlasd3.f"
    }
#line 297 "dlasd3.f"
    if (*info != 0) {
#line 298 "dlasd3.f"
	i__1 = -(*info);
#line 298 "dlasd3.f"
	xerbla_("DLASD3", &i__1, (ftnlen)6);
#line 299 "dlasd3.f"
	return 0;
#line 300 "dlasd3.f"
    }

/*     Quick return if possible */

#line 304 "dlasd3.f"
    if (*k == 1) {
#line 305 "dlasd3.f"
	d__[1] = abs(z__[1]);
#line 306 "dlasd3.f"
	dcopy_(&m, &vt2[vt2_dim1 + 1], ldvt2, &vt[vt_dim1 + 1], ldvt);
#line 307 "dlasd3.f"
	if (z__[1] > 0.) {
#line 308 "dlasd3.f"
	    dcopy_(&n, &u2[u2_dim1 + 1], &c__1, &u[u_dim1 + 1], &c__1);
#line 309 "dlasd3.f"
	} else {
#line 310 "dlasd3.f"
	    i__1 = n;
#line 310 "dlasd3.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 311 "dlasd3.f"
		u[i__ + u_dim1] = -u2[i__ + u2_dim1];
#line 312 "dlasd3.f"
/* L10: */
#line 312 "dlasd3.f"
	    }
#line 313 "dlasd3.f"
	}
#line 314 "dlasd3.f"
	return 0;
#line 315 "dlasd3.f"
    }

/*     Modify values DSIGMA(i) to make sure all DSIGMA(i)-DSIGMA(j) can */
/*     be computed with high relative accuracy (barring over/underflow). */
/*     This is a problem on machines without a guard digit in */
/*     add/subtract (Cray XMP, Cray YMP, Cray C 90 and Cray 2). */
/*     The following code replaces DSIGMA(I) by 2*DSIGMA(I)-DSIGMA(I), */
/*     which on any of these machines zeros out the bottommost */
/*     bit of DSIGMA(I) if it is 1; this makes the subsequent */
/*     subtractions DSIGMA(I)-DSIGMA(J) unproblematic when cancellation */
/*     occurs. On binary machines with a guard digit (almost all */
/*     machines) it does not change DSIGMA(I) at all. On hexadecimal */
/*     and decimal machines with a guard digit, it slightly */
/*     changes the bottommost bits of DSIGMA(I). It does not account */
/*     for hexadecimal or decimal machines without guard digits */
/*     (we know of none). We use a subroutine call to compute */
/*     2*DSIGMA(I) to prevent optimizing compilers from eliminating */
/*     this code. */

#line 334 "dlasd3.f"
    i__1 = *k;
#line 334 "dlasd3.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 335 "dlasd3.f"
	dsigma[i__] = dlamc3_(&dsigma[i__], &dsigma[i__]) - dsigma[i__];
#line 336 "dlasd3.f"
/* L20: */
#line 336 "dlasd3.f"
    }

/*     Keep a copy of Z. */

#line 340 "dlasd3.f"
    dcopy_(k, &z__[1], &c__1, &q[q_offset], &c__1);

/*     Normalize Z. */

#line 344 "dlasd3.f"
    rho = dnrm2_(k, &z__[1], &c__1);
#line 345 "dlasd3.f"
    dlascl_("G", &c__0, &c__0, &rho, &c_b13, k, &c__1, &z__[1], k, info, (
	    ftnlen)1);
#line 346 "dlasd3.f"
    rho *= rho;

/*     Find the new singular values. */

#line 350 "dlasd3.f"
    i__1 = *k;
#line 350 "dlasd3.f"
    for (j = 1; j <= i__1; ++j) {
#line 351 "dlasd3.f"
	dlasd4_(k, &j, &dsigma[1], &z__[1], &u[j * u_dim1 + 1], &rho, &d__[j],
		 &vt[j * vt_dim1 + 1], info);

/*        If the zero finder fails, the computation is terminated. */

#line 356 "dlasd3.f"
	if (*info != 0) {
#line 357 "dlasd3.f"
	    return 0;
#line 358 "dlasd3.f"
	}
#line 359 "dlasd3.f"
/* L30: */
#line 359 "dlasd3.f"
    }

/*     Compute updated Z. */

#line 363 "dlasd3.f"
    i__1 = *k;
#line 363 "dlasd3.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 364 "dlasd3.f"
	z__[i__] = u[i__ + *k * u_dim1] * vt[i__ + *k * vt_dim1];
#line 365 "dlasd3.f"
	i__2 = i__ - 1;
#line 365 "dlasd3.f"
	for (j = 1; j <= i__2; ++j) {
#line 366 "dlasd3.f"
	    z__[i__] *= u[i__ + j * u_dim1] * vt[i__ + j * vt_dim1] / (dsigma[
		    i__] - dsigma[j]) / (dsigma[i__] + dsigma[j]);
#line 369 "dlasd3.f"
/* L40: */
#line 369 "dlasd3.f"
	}
#line 370 "dlasd3.f"
	i__2 = *k - 1;
#line 370 "dlasd3.f"
	for (j = i__; j <= i__2; ++j) {
#line 371 "dlasd3.f"
	    z__[i__] *= u[i__ + j * u_dim1] * vt[i__ + j * vt_dim1] / (dsigma[
		    i__] - dsigma[j + 1]) / (dsigma[i__] + dsigma[j + 1]);
#line 374 "dlasd3.f"
/* L50: */
#line 374 "dlasd3.f"
	}
#line 375 "dlasd3.f"
	d__2 = sqrt((d__1 = z__[i__], abs(d__1)));
#line 375 "dlasd3.f"
	z__[i__] = d_sign(&d__2, &q[i__ + q_dim1]);
#line 376 "dlasd3.f"
/* L60: */
#line 376 "dlasd3.f"
    }

/*     Compute left singular vectors of the modified diagonal matrix, */
/*     and store related information for the right singular vectors. */

#line 381 "dlasd3.f"
    i__1 = *k;
#line 381 "dlasd3.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 382 "dlasd3.f"
	vt[i__ * vt_dim1 + 1] = z__[1] / u[i__ * u_dim1 + 1] / vt[i__ * 
		vt_dim1 + 1];
#line 383 "dlasd3.f"
	u[i__ * u_dim1 + 1] = -1.;
#line 384 "dlasd3.f"
	i__2 = *k;
#line 384 "dlasd3.f"
	for (j = 2; j <= i__2; ++j) {
#line 385 "dlasd3.f"
	    vt[j + i__ * vt_dim1] = z__[j] / u[j + i__ * u_dim1] / vt[j + i__ 
		    * vt_dim1];
#line 386 "dlasd3.f"
	    u[j + i__ * u_dim1] = dsigma[j] * vt[j + i__ * vt_dim1];
#line 387 "dlasd3.f"
/* L70: */
#line 387 "dlasd3.f"
	}
#line 388 "dlasd3.f"
	temp = dnrm2_(k, &u[i__ * u_dim1 + 1], &c__1);
#line 389 "dlasd3.f"
	q[i__ * q_dim1 + 1] = u[i__ * u_dim1 + 1] / temp;
#line 390 "dlasd3.f"
	i__2 = *k;
#line 390 "dlasd3.f"
	for (j = 2; j <= i__2; ++j) {
#line 391 "dlasd3.f"
	    jc = idxc[j];
#line 392 "dlasd3.f"
	    q[j + i__ * q_dim1] = u[jc + i__ * u_dim1] / temp;
#line 393 "dlasd3.f"
/* L80: */
#line 393 "dlasd3.f"
	}
#line 394 "dlasd3.f"
/* L90: */
#line 394 "dlasd3.f"
    }

/*     Update the left singular vector matrix. */

#line 398 "dlasd3.f"
    if (*k == 2) {
#line 399 "dlasd3.f"
	dgemm_("N", "N", &n, k, k, &c_b13, &u2[u2_offset], ldu2, &q[q_offset],
		 ldq, &c_b26, &u[u_offset], ldu, (ftnlen)1, (ftnlen)1);
#line 401 "dlasd3.f"
	goto L100;
#line 402 "dlasd3.f"
    }
#line 403 "dlasd3.f"
    if (ctot[1] > 0) {
#line 404 "dlasd3.f"
	dgemm_("N", "N", nl, k, &ctot[1], &c_b13, &u2[(u2_dim1 << 1) + 1], 
		ldu2, &q[q_dim1 + 2], ldq, &c_b26, &u[u_dim1 + 1], ldu, (
		ftnlen)1, (ftnlen)1);
#line 406 "dlasd3.f"
	if (ctot[3] > 0) {
#line 407 "dlasd3.f"
	    ktemp = ctot[1] + 2 + ctot[2];
#line 408 "dlasd3.f"
	    dgemm_("N", "N", nl, k, &ctot[3], &c_b13, &u2[ktemp * u2_dim1 + 1]
		    , ldu2, &q[ktemp + q_dim1], ldq, &c_b13, &u[u_dim1 + 1], 
		    ldu, (ftnlen)1, (ftnlen)1);
#line 410 "dlasd3.f"
	}
#line 411 "dlasd3.f"
    } else if (ctot[3] > 0) {
#line 412 "dlasd3.f"
	ktemp = ctot[1] + 2 + ctot[2];
#line 413 "dlasd3.f"
	dgemm_("N", "N", nl, k, &ctot[3], &c_b13, &u2[ktemp * u2_dim1 + 1], 
		ldu2, &q[ktemp + q_dim1], ldq, &c_b26, &u[u_dim1 + 1], ldu, (
		ftnlen)1, (ftnlen)1);
#line 415 "dlasd3.f"
    } else {
#line 416 "dlasd3.f"
	dlacpy_("F", nl, k, &u2[u2_offset], ldu2, &u[u_offset], ldu, (ftnlen)
		1);
#line 417 "dlasd3.f"
    }
#line 418 "dlasd3.f"
    dcopy_(k, &q[q_dim1 + 1], ldq, &u[nlp1 + u_dim1], ldu);
#line 419 "dlasd3.f"
    ktemp = ctot[1] + 2;
#line 420 "dlasd3.f"
    ctemp = ctot[2] + ctot[3];
#line 421 "dlasd3.f"
    dgemm_("N", "N", nr, k, &ctemp, &c_b13, &u2[nlp2 + ktemp * u2_dim1], ldu2,
	     &q[ktemp + q_dim1], ldq, &c_b26, &u[nlp2 + u_dim1], ldu, (ftnlen)
	    1, (ftnlen)1);

/*     Generate the right singular vectors. */

#line 426 "dlasd3.f"
L100:
#line 427 "dlasd3.f"
    i__1 = *k;
#line 427 "dlasd3.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 428 "dlasd3.f"
	temp = dnrm2_(k, &vt[i__ * vt_dim1 + 1], &c__1);
#line 429 "dlasd3.f"
	q[i__ + q_dim1] = vt[i__ * vt_dim1 + 1] / temp;
#line 430 "dlasd3.f"
	i__2 = *k;
#line 430 "dlasd3.f"
	for (j = 2; j <= i__2; ++j) {
#line 431 "dlasd3.f"
	    jc = idxc[j];
#line 432 "dlasd3.f"
	    q[i__ + j * q_dim1] = vt[jc + i__ * vt_dim1] / temp;
#line 433 "dlasd3.f"
/* L110: */
#line 433 "dlasd3.f"
	}
#line 434 "dlasd3.f"
/* L120: */
#line 434 "dlasd3.f"
    }

/*     Update the right singular vector matrix. */

#line 438 "dlasd3.f"
    if (*k == 2) {
#line 439 "dlasd3.f"
	dgemm_("N", "N", k, &m, k, &c_b13, &q[q_offset], ldq, &vt2[vt2_offset]
		, ldvt2, &c_b26, &vt[vt_offset], ldvt, (ftnlen)1, (ftnlen)1);
#line 441 "dlasd3.f"
	return 0;
#line 442 "dlasd3.f"
    }
#line 443 "dlasd3.f"
    ktemp = ctot[1] + 1;
#line 444 "dlasd3.f"
    dgemm_("N", "N", k, &nlp1, &ktemp, &c_b13, &q[q_dim1 + 1], ldq, &vt2[
	    vt2_dim1 + 1], ldvt2, &c_b26, &vt[vt_dim1 + 1], ldvt, (ftnlen)1, (
	    ftnlen)1);
#line 446 "dlasd3.f"
    ktemp = ctot[1] + 2 + ctot[2];
#line 447 "dlasd3.f"
    if (ktemp <= *ldvt2) {
#line 447 "dlasd3.f"
	dgemm_("N", "N", k, &nlp1, &ctot[3], &c_b13, &q[ktemp * q_dim1 + 1], 
		ldq, &vt2[ktemp + vt2_dim1], ldvt2, &c_b13, &vt[vt_dim1 + 1], 
		ldvt, (ftnlen)1, (ftnlen)1);
#line 447 "dlasd3.f"
    }

#line 452 "dlasd3.f"
    ktemp = ctot[1] + 1;
#line 453 "dlasd3.f"
    nrp1 = *nr + *sqre;
#line 454 "dlasd3.f"
    if (ktemp > 1) {
#line 455 "dlasd3.f"
	i__1 = *k;
#line 455 "dlasd3.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 456 "dlasd3.f"
	    q[i__ + ktemp * q_dim1] = q[i__ + q_dim1];
#line 457 "dlasd3.f"
/* L130: */
#line 457 "dlasd3.f"
	}
#line 458 "dlasd3.f"
	i__1 = m;
#line 458 "dlasd3.f"
	for (i__ = nlp2; i__ <= i__1; ++i__) {
#line 459 "dlasd3.f"
	    vt2[ktemp + i__ * vt2_dim1] = vt2[i__ * vt2_dim1 + 1];
#line 460 "dlasd3.f"
/* L140: */
#line 460 "dlasd3.f"
	}
#line 461 "dlasd3.f"
    }
#line 462 "dlasd3.f"
    ctemp = ctot[2] + 1 + ctot[3];
#line 463 "dlasd3.f"
    dgemm_("N", "N", k, &nrp1, &ctemp, &c_b13, &q[ktemp * q_dim1 + 1], ldq, &
	    vt2[ktemp + nlp2 * vt2_dim1], ldvt2, &c_b26, &vt[nlp2 * vt_dim1 + 
	    1], ldvt, (ftnlen)1, (ftnlen)1);

#line 466 "dlasd3.f"
    return 0;

/*     End of DLASD3 */

} /* dlasd3_ */

