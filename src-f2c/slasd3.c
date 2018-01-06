#line 1 "slasd3.f"
/* slasd3.f -- translated by f2c (version 20100827).
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

#line 1 "slasd3.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c__0 = 0;
static doublereal c_b13 = 1.;
static doublereal c_b26 = 0.;

/* > \brief \b SLASD3 finds all square roots of the roots of the secular equation, as defined by the values in
 D and Z, and then updates the singular vectors by matrix multiplication. Used by sbdsdc. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLASD3 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slasd3.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slasd3.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slasd3.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLASD3( NL, NR, SQRE, K, D, Q, LDQ, DSIGMA, U, LDU, U2, */
/*                          LDU2, VT, LDVT, VT2, LDVT2, IDXC, CTOT, Z, */
/*                          INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, K, LDQ, LDU, LDU2, LDVT, LDVT2, NL, NR, */
/*      $                   SQRE */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            CTOT( * ), IDXC( * ) */
/*       REAL               D( * ), DSIGMA( * ), Q( LDQ, * ), U( LDU, * ), */
/*      $                   U2( LDU2, * ), VT( LDVT, * ), VT2( LDVT2, * ), */
/*      $                   Z( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLASD3 finds all the square roots of the roots of the secular */
/* > equation, as defined by the values in D and Z.  It makes the */
/* > appropriate calls to SLASD4 and then updates the singular */
/* > vectors by matrix multiplication. */
/* > */
/* > This code makes very mild assumptions about floating point */
/* > arithmetic. It will work on machines with a guard digit in */
/* > add/subtract, or on those binary machines without guard digits */
/* > which subtract like the Cray XMP, Cray YMP, Cray C 90, or Cray 2. */
/* > It could conceivably fail on hexadecimal or decimal machines */
/* > without guard digits, but we know of none. */
/* > */
/* > SLASD3 is called from SLASD1. */
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
/* >          D is REAL array, dimension(K) */
/* >         On exit the square roots of the roots of the secular equation, */
/* >         in ascending order. */
/* > \endverbatim */
/* > */
/* > \param[out] Q */
/* > \verbatim */
/* >          Q is REAL array, */
/* >                     dimension at least (LDQ,K). */
/* > \endverbatim */
/* > */
/* > \param[in] LDQ */
/* > \verbatim */
/* >          LDQ is INTEGER */
/* >         The leading dimension of the array Q.  LDQ >= K. */
/* > \endverbatim */
/* > */
/* > \param[in,out] DSIGMA */
/* > \verbatim */
/* >          DSIGMA is REAL array, dimension(K) */
/* >         The first K elements of this array contain the old roots */
/* >         of the deflated updating problem.  These are the poles */
/* >         of the secular equation. */
/* > \endverbatim */
/* > */
/* > \param[out] U */
/* > \verbatim */
/* >          U is REAL array, dimension (LDU, N) */
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
/* > \param[in] U2 */
/* > \verbatim */
/* >          U2 is REAL array, dimension (LDU2, N) */
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
/* >          VT is REAL array, dimension (LDVT, M) */
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
/* >          VT2 is REAL array, dimension (LDVT2, N) */
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
/* >          IDXC is INTEGER array, dimension (N) */
/* >         The permutation used to arrange the columns of U (and rows of */
/* >         VT) into three groups:  the first group contains non-zero */
/* >         entries only at and above (or before) NL +1; the second */
/* >         contains non-zero entries only at and below (or after) NL+2; */
/* >         and the third is dense. The first column of U and the row of */
/* >         VT are treated separately, however. */
/* > */
/* >         The rows of the singular vectors found by SLASD4 */
/* >         must be likewise permuted before the matrix multiplies can */
/* >         take place. */
/* > \endverbatim */
/* > */
/* > \param[in] CTOT */
/* > \verbatim */
/* >          CTOT is INTEGER array, dimension (4) */
/* >         A count of the total number of the various types of columns */
/* >         in U (or rows in VT), as described in IDXC. The fourth column */
/* >         type is any column which has been deflated. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Z */
/* > \verbatim */
/* >          Z is REAL array, dimension (K) */
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

/* > \date November 2015 */

/* > \ingroup auxOTHERauxiliary */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Ming Gu and Huan Ren, Computer Science Division, University of */
/* >     California at Berkeley, USA */
/* > */
/*  ===================================================================== */
/* Subroutine */ int slasd3_(integer *nl, integer *nr, integer *sqre, integer 
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
    extern doublereal snrm2_(integer *, doublereal *, integer *);
    static integer ctemp;
    extern /* Subroutine */ int sgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    static integer ktemp;
    extern /* Subroutine */ int scopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    extern doublereal slamc3_(doublereal *, doublereal *);
    extern /* Subroutine */ int slasd4_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *), xerbla_(char *, integer *, ftnlen), 
	    slascl_(char *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, integer *, doublereal *, integer *, integer *, ftnlen),
	     slacpy_(char *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen);


/*  -- LAPACK auxiliary routine (version 3.6.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2015 */

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

#line 269 "slasd3.f"
    /* Parameter adjustments */
#line 269 "slasd3.f"
    --d__;
#line 269 "slasd3.f"
    q_dim1 = *ldq;
#line 269 "slasd3.f"
    q_offset = 1 + q_dim1;
#line 269 "slasd3.f"
    q -= q_offset;
#line 269 "slasd3.f"
    --dsigma;
#line 269 "slasd3.f"
    u_dim1 = *ldu;
#line 269 "slasd3.f"
    u_offset = 1 + u_dim1;
#line 269 "slasd3.f"
    u -= u_offset;
#line 269 "slasd3.f"
    u2_dim1 = *ldu2;
#line 269 "slasd3.f"
    u2_offset = 1 + u2_dim1;
#line 269 "slasd3.f"
    u2 -= u2_offset;
#line 269 "slasd3.f"
    vt_dim1 = *ldvt;
#line 269 "slasd3.f"
    vt_offset = 1 + vt_dim1;
#line 269 "slasd3.f"
    vt -= vt_offset;
#line 269 "slasd3.f"
    vt2_dim1 = *ldvt2;
#line 269 "slasd3.f"
    vt2_offset = 1 + vt2_dim1;
#line 269 "slasd3.f"
    vt2 -= vt2_offset;
#line 269 "slasd3.f"
    --idxc;
#line 269 "slasd3.f"
    --ctot;
#line 269 "slasd3.f"
    --z__;
#line 269 "slasd3.f"

#line 269 "slasd3.f"
    /* Function Body */
#line 269 "slasd3.f"
    *info = 0;

#line 271 "slasd3.f"
    if (*nl < 1) {
#line 272 "slasd3.f"
	*info = -1;
#line 273 "slasd3.f"
    } else if (*nr < 1) {
#line 274 "slasd3.f"
	*info = -2;
#line 275 "slasd3.f"
    } else if (*sqre != 1 && *sqre != 0) {
#line 276 "slasd3.f"
	*info = -3;
#line 277 "slasd3.f"
    }

#line 279 "slasd3.f"
    n = *nl + *nr + 1;
#line 280 "slasd3.f"
    m = n + *sqre;
#line 281 "slasd3.f"
    nlp1 = *nl + 1;
#line 282 "slasd3.f"
    nlp2 = *nl + 2;

#line 284 "slasd3.f"
    if (*k < 1 || *k > n) {
#line 285 "slasd3.f"
	*info = -4;
#line 286 "slasd3.f"
    } else if (*ldq < *k) {
#line 287 "slasd3.f"
	*info = -7;
#line 288 "slasd3.f"
    } else if (*ldu < n) {
#line 289 "slasd3.f"
	*info = -10;
#line 290 "slasd3.f"
    } else if (*ldu2 < n) {
#line 291 "slasd3.f"
	*info = -12;
#line 292 "slasd3.f"
    } else if (*ldvt < m) {
#line 293 "slasd3.f"
	*info = -14;
#line 294 "slasd3.f"
    } else if (*ldvt2 < m) {
#line 295 "slasd3.f"
	*info = -16;
#line 296 "slasd3.f"
    }
#line 297 "slasd3.f"
    if (*info != 0) {
#line 298 "slasd3.f"
	i__1 = -(*info);
#line 298 "slasd3.f"
	xerbla_("SLASD3", &i__1, (ftnlen)6);
#line 299 "slasd3.f"
	return 0;
#line 300 "slasd3.f"
    }

/*     Quick return if possible */

#line 304 "slasd3.f"
    if (*k == 1) {
#line 305 "slasd3.f"
	d__[1] = abs(z__[1]);
#line 306 "slasd3.f"
	scopy_(&m, &vt2[vt2_dim1 + 1], ldvt2, &vt[vt_dim1 + 1], ldvt);
#line 307 "slasd3.f"
	if (z__[1] > 0.) {
#line 308 "slasd3.f"
	    scopy_(&n, &u2[u2_dim1 + 1], &c__1, &u[u_dim1 + 1], &c__1);
#line 309 "slasd3.f"
	} else {
#line 310 "slasd3.f"
	    i__1 = n;
#line 310 "slasd3.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 311 "slasd3.f"
		u[i__ + u_dim1] = -u2[i__ + u2_dim1];
#line 312 "slasd3.f"
/* L10: */
#line 312 "slasd3.f"
	    }
#line 313 "slasd3.f"
	}
#line 314 "slasd3.f"
	return 0;
#line 315 "slasd3.f"
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

#line 334 "slasd3.f"
    i__1 = *k;
#line 334 "slasd3.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 335 "slasd3.f"
	dsigma[i__] = slamc3_(&dsigma[i__], &dsigma[i__]) - dsigma[i__];
#line 336 "slasd3.f"
/* L20: */
#line 336 "slasd3.f"
    }

/*     Keep a copy of Z. */

#line 340 "slasd3.f"
    scopy_(k, &z__[1], &c__1, &q[q_offset], &c__1);

/*     Normalize Z. */

#line 344 "slasd3.f"
    rho = snrm2_(k, &z__[1], &c__1);
#line 345 "slasd3.f"
    slascl_("G", &c__0, &c__0, &rho, &c_b13, k, &c__1, &z__[1], k, info, (
	    ftnlen)1);
#line 346 "slasd3.f"
    rho *= rho;

/*     Find the new singular values. */

#line 350 "slasd3.f"
    i__1 = *k;
#line 350 "slasd3.f"
    for (j = 1; j <= i__1; ++j) {
#line 351 "slasd3.f"
	slasd4_(k, &j, &dsigma[1], &z__[1], &u[j * u_dim1 + 1], &rho, &d__[j],
		 &vt[j * vt_dim1 + 1], info);

/*        If the zero finder fails, report the convergence failure. */

#line 356 "slasd3.f"
	if (*info != 0) {
#line 357 "slasd3.f"
	    return 0;
#line 358 "slasd3.f"
	}
#line 359 "slasd3.f"
/* L30: */
#line 359 "slasd3.f"
    }

/*     Compute updated Z. */

#line 363 "slasd3.f"
    i__1 = *k;
#line 363 "slasd3.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 364 "slasd3.f"
	z__[i__] = u[i__ + *k * u_dim1] * vt[i__ + *k * vt_dim1];
#line 365 "slasd3.f"
	i__2 = i__ - 1;
#line 365 "slasd3.f"
	for (j = 1; j <= i__2; ++j) {
#line 366 "slasd3.f"
	    z__[i__] *= u[i__ + j * u_dim1] * vt[i__ + j * vt_dim1] / (dsigma[
		    i__] - dsigma[j]) / (dsigma[i__] + dsigma[j]);
#line 369 "slasd3.f"
/* L40: */
#line 369 "slasd3.f"
	}
#line 370 "slasd3.f"
	i__2 = *k - 1;
#line 370 "slasd3.f"
	for (j = i__; j <= i__2; ++j) {
#line 371 "slasd3.f"
	    z__[i__] *= u[i__ + j * u_dim1] * vt[i__ + j * vt_dim1] / (dsigma[
		    i__] - dsigma[j + 1]) / (dsigma[i__] + dsigma[j + 1]);
#line 374 "slasd3.f"
/* L50: */
#line 374 "slasd3.f"
	}
#line 375 "slasd3.f"
	d__2 = sqrt((d__1 = z__[i__], abs(d__1)));
#line 375 "slasd3.f"
	z__[i__] = d_sign(&d__2, &q[i__ + q_dim1]);
#line 376 "slasd3.f"
/* L60: */
#line 376 "slasd3.f"
    }

/*     Compute left singular vectors of the modified diagonal matrix, */
/*     and store related information for the right singular vectors. */

#line 381 "slasd3.f"
    i__1 = *k;
#line 381 "slasd3.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 382 "slasd3.f"
	vt[i__ * vt_dim1 + 1] = z__[1] / u[i__ * u_dim1 + 1] / vt[i__ * 
		vt_dim1 + 1];
#line 383 "slasd3.f"
	u[i__ * u_dim1 + 1] = -1.;
#line 384 "slasd3.f"
	i__2 = *k;
#line 384 "slasd3.f"
	for (j = 2; j <= i__2; ++j) {
#line 385 "slasd3.f"
	    vt[j + i__ * vt_dim1] = z__[j] / u[j + i__ * u_dim1] / vt[j + i__ 
		    * vt_dim1];
#line 386 "slasd3.f"
	    u[j + i__ * u_dim1] = dsigma[j] * vt[j + i__ * vt_dim1];
#line 387 "slasd3.f"
/* L70: */
#line 387 "slasd3.f"
	}
#line 388 "slasd3.f"
	temp = snrm2_(k, &u[i__ * u_dim1 + 1], &c__1);
#line 389 "slasd3.f"
	q[i__ * q_dim1 + 1] = u[i__ * u_dim1 + 1] / temp;
#line 390 "slasd3.f"
	i__2 = *k;
#line 390 "slasd3.f"
	for (j = 2; j <= i__2; ++j) {
#line 391 "slasd3.f"
	    jc = idxc[j];
#line 392 "slasd3.f"
	    q[j + i__ * q_dim1] = u[jc + i__ * u_dim1] / temp;
#line 393 "slasd3.f"
/* L80: */
#line 393 "slasd3.f"
	}
#line 394 "slasd3.f"
/* L90: */
#line 394 "slasd3.f"
    }

/*     Update the left singular vector matrix. */

#line 398 "slasd3.f"
    if (*k == 2) {
#line 399 "slasd3.f"
	sgemm_("N", "N", &n, k, k, &c_b13, &u2[u2_offset], ldu2, &q[q_offset],
		 ldq, &c_b26, &u[u_offset], ldu, (ftnlen)1, (ftnlen)1);
#line 401 "slasd3.f"
	goto L100;
#line 402 "slasd3.f"
    }
#line 403 "slasd3.f"
    if (ctot[1] > 0) {
#line 404 "slasd3.f"
	sgemm_("N", "N", nl, k, &ctot[1], &c_b13, &u2[(u2_dim1 << 1) + 1], 
		ldu2, &q[q_dim1 + 2], ldq, &c_b26, &u[u_dim1 + 1], ldu, (
		ftnlen)1, (ftnlen)1);
#line 406 "slasd3.f"
	if (ctot[3] > 0) {
#line 407 "slasd3.f"
	    ktemp = ctot[1] + 2 + ctot[2];
#line 408 "slasd3.f"
	    sgemm_("N", "N", nl, k, &ctot[3], &c_b13, &u2[ktemp * u2_dim1 + 1]
		    , ldu2, &q[ktemp + q_dim1], ldq, &c_b13, &u[u_dim1 + 1], 
		    ldu, (ftnlen)1, (ftnlen)1);
#line 410 "slasd3.f"
	}
#line 411 "slasd3.f"
    } else if (ctot[3] > 0) {
#line 412 "slasd3.f"
	ktemp = ctot[1] + 2 + ctot[2];
#line 413 "slasd3.f"
	sgemm_("N", "N", nl, k, &ctot[3], &c_b13, &u2[ktemp * u2_dim1 + 1], 
		ldu2, &q[ktemp + q_dim1], ldq, &c_b26, &u[u_dim1 + 1], ldu, (
		ftnlen)1, (ftnlen)1);
#line 415 "slasd3.f"
    } else {
#line 416 "slasd3.f"
	slacpy_("F", nl, k, &u2[u2_offset], ldu2, &u[u_offset], ldu, (ftnlen)
		1);
#line 417 "slasd3.f"
    }
#line 418 "slasd3.f"
    scopy_(k, &q[q_dim1 + 1], ldq, &u[nlp1 + u_dim1], ldu);
#line 419 "slasd3.f"
    ktemp = ctot[1] + 2;
#line 420 "slasd3.f"
    ctemp = ctot[2] + ctot[3];
#line 421 "slasd3.f"
    sgemm_("N", "N", nr, k, &ctemp, &c_b13, &u2[nlp2 + ktemp * u2_dim1], ldu2,
	     &q[ktemp + q_dim1], ldq, &c_b26, &u[nlp2 + u_dim1], ldu, (ftnlen)
	    1, (ftnlen)1);

/*     Generate the right singular vectors. */

#line 426 "slasd3.f"
L100:
#line 427 "slasd3.f"
    i__1 = *k;
#line 427 "slasd3.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 428 "slasd3.f"
	temp = snrm2_(k, &vt[i__ * vt_dim1 + 1], &c__1);
#line 429 "slasd3.f"
	q[i__ + q_dim1] = vt[i__ * vt_dim1 + 1] / temp;
#line 430 "slasd3.f"
	i__2 = *k;
#line 430 "slasd3.f"
	for (j = 2; j <= i__2; ++j) {
#line 431 "slasd3.f"
	    jc = idxc[j];
#line 432 "slasd3.f"
	    q[i__ + j * q_dim1] = vt[jc + i__ * vt_dim1] / temp;
#line 433 "slasd3.f"
/* L110: */
#line 433 "slasd3.f"
	}
#line 434 "slasd3.f"
/* L120: */
#line 434 "slasd3.f"
    }

/*     Update the right singular vector matrix. */

#line 438 "slasd3.f"
    if (*k == 2) {
#line 439 "slasd3.f"
	sgemm_("N", "N", k, &m, k, &c_b13, &q[q_offset], ldq, &vt2[vt2_offset]
		, ldvt2, &c_b26, &vt[vt_offset], ldvt, (ftnlen)1, (ftnlen)1);
#line 441 "slasd3.f"
	return 0;
#line 442 "slasd3.f"
    }
#line 443 "slasd3.f"
    ktemp = ctot[1] + 1;
#line 444 "slasd3.f"
    sgemm_("N", "N", k, &nlp1, &ktemp, &c_b13, &q[q_dim1 + 1], ldq, &vt2[
	    vt2_dim1 + 1], ldvt2, &c_b26, &vt[vt_dim1 + 1], ldvt, (ftnlen)1, (
	    ftnlen)1);
#line 446 "slasd3.f"
    ktemp = ctot[1] + 2 + ctot[2];
#line 447 "slasd3.f"
    if (ktemp <= *ldvt2) {
#line 447 "slasd3.f"
	sgemm_("N", "N", k, &nlp1, &ctot[3], &c_b13, &q[ktemp * q_dim1 + 1], 
		ldq, &vt2[ktemp + vt2_dim1], ldvt2, &c_b13, &vt[vt_dim1 + 1], 
		ldvt, (ftnlen)1, (ftnlen)1);
#line 447 "slasd3.f"
    }

#line 452 "slasd3.f"
    ktemp = ctot[1] + 1;
#line 453 "slasd3.f"
    nrp1 = *nr + *sqre;
#line 454 "slasd3.f"
    if (ktemp > 1) {
#line 455 "slasd3.f"
	i__1 = *k;
#line 455 "slasd3.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 456 "slasd3.f"
	    q[i__ + ktemp * q_dim1] = q[i__ + q_dim1];
#line 457 "slasd3.f"
/* L130: */
#line 457 "slasd3.f"
	}
#line 458 "slasd3.f"
	i__1 = m;
#line 458 "slasd3.f"
	for (i__ = nlp2; i__ <= i__1; ++i__) {
#line 459 "slasd3.f"
	    vt2[ktemp + i__ * vt2_dim1] = vt2[i__ * vt2_dim1 + 1];
#line 460 "slasd3.f"
/* L140: */
#line 460 "slasd3.f"
	}
#line 461 "slasd3.f"
    }
#line 462 "slasd3.f"
    ctemp = ctot[2] + 1 + ctot[3];
#line 463 "slasd3.f"
    sgemm_("N", "N", k, &nrp1, &ctemp, &c_b13, &q[ktemp * q_dim1 + 1], ldq, &
	    vt2[ktemp + nlp2 * vt2_dim1], ldvt2, &c_b26, &vt[nlp2 * vt_dim1 + 
	    1], ldvt, (ftnlen)1, (ftnlen)1);

#line 466 "slasd3.f"
    return 0;

/*     End of SLASD3 */

} /* slasd3_ */

