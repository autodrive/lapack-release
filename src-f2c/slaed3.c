#line 1 "slaed3.f"
/* slaed3.f -- translated by f2c (version 20100827).
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

#line 1 "slaed3.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b22 = 1.;
static doublereal c_b23 = 0.;

/* > \brief \b SLAED3 used by sstedc. Finds the roots of the secular equation and updates the eigenvectors. Us
ed when the original matrix is tridiagonal. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLAED3 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaed3.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaed3.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaed3.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLAED3( K, N, N1, D, Q, LDQ, RHO, DLAMDA, Q2, INDX, */
/*                          CTOT, W, S, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, K, LDQ, N, N1 */
/*       REAL               RHO */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            CTOT( * ), INDX( * ) */
/*       REAL               D( * ), DLAMDA( * ), Q( LDQ, * ), Q2( * ), */
/*      $                   S( * ), W( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLAED3 finds the roots of the secular equation, as defined by the */
/* > values in D, W, and RHO, between 1 and K.  It makes the */
/* > appropriate calls to SLAED4 and then updates the eigenvectors by */
/* > multiplying the matrix of eigenvectors of the pair of eigensystems */
/* > being combined by the matrix of eigenvectors of the K-by-K system */
/* > which is solved here. */
/* > */
/* > This code makes very mild assumptions about floating point */
/* > arithmetic. It will work on machines with a guard digit in */
/* > add/subtract, or on those binary machines without guard digits */
/* > which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or Cray-2. */
/* > It could conceivably fail on hexadecimal or decimal machines */
/* > without guard digits, but we know of none. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] K */
/* > \verbatim */
/* >          K is INTEGER */
/* >          The number of terms in the rational function to be solved by */
/* >          SLAED4.  K >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of rows and columns in the Q matrix. */
/* >          N >= K (deflation may result in N>K). */
/* > \endverbatim */
/* > */
/* > \param[in] N1 */
/* > \verbatim */
/* >          N1 is INTEGER */
/* >          The location of the last eigenvalue in the leading submatrix. */
/* >          min(1,N) <= N1 <= N/2. */
/* > \endverbatim */
/* > */
/* > \param[out] D */
/* > \verbatim */
/* >          D is REAL array, dimension (N) */
/* >          D(I) contains the updated eigenvalues for */
/* >          1 <= I <= K. */
/* > \endverbatim */
/* > */
/* > \param[out] Q */
/* > \verbatim */
/* >          Q is REAL array, dimension (LDQ,N) */
/* >          Initially the first K columns are used as workspace. */
/* >          On output the columns 1 to K contain */
/* >          the updated eigenvectors. */
/* > \endverbatim */
/* > */
/* > \param[in] LDQ */
/* > \verbatim */
/* >          LDQ is INTEGER */
/* >          The leading dimension of the array Q.  LDQ >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] RHO */
/* > \verbatim */
/* >          RHO is REAL */
/* >          The value of the parameter in the rank one update equation. */
/* >          RHO >= 0 required. */
/* > \endverbatim */
/* > */
/* > \param[in,out] DLAMDA */
/* > \verbatim */
/* >          DLAMDA is REAL array, dimension (K) */
/* >          The first K elements of this array contain the old roots */
/* >          of the deflated updating problem.  These are the poles */
/* >          of the secular equation. May be changed on output by */
/* >          having lowest order bit set to zero on Cray X-MP, Cray Y-MP, */
/* >          Cray-2, or Cray C-90, as described above. */
/* > \endverbatim */
/* > */
/* > \param[in] Q2 */
/* > \verbatim */
/* >          Q2 is REAL array, dimension (LDQ2, N) */
/* >          The first K columns of this matrix contain the non-deflated */
/* >          eigenvectors for the split problem. */
/* > \endverbatim */
/* > */
/* > \param[in] INDX */
/* > \verbatim */
/* >          INDX is INTEGER array, dimension (N) */
/* >          The permutation used to arrange the columns of the deflated */
/* >          Q matrix into three groups (see SLAED2). */
/* >          The rows of the eigenvectors found by SLAED4 must be likewise */
/* >          permuted before the matrix multiply can take place. */
/* > \endverbatim */
/* > */
/* > \param[in] CTOT */
/* > \verbatim */
/* >          CTOT is INTEGER array, dimension (4) */
/* >          A count of the total number of the various types of columns */
/* >          in Q, as described in INDX.  The fourth column type is any */
/* >          column which has been deflated. */
/* > \endverbatim */
/* > */
/* > \param[in,out] W */
/* > \verbatim */
/* >          W is REAL array, dimension (K) */
/* >          The first K elements of this array contain the components */
/* >          of the deflation-adjusted updating vector. Destroyed on */
/* >          output. */
/* > \endverbatim */
/* > */
/* > \param[out] S */
/* > \verbatim */
/* >          S is REAL array, dimension (N1 + 1)*K */
/* >          Will contain the eigenvectors of the repaired matrix which */
/* >          will be multiplied by the previously accumulated eigenvectors */
/* >          to update the system. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit. */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value. */
/* >          > 0:  if INFO = 1, an eigenvalue did not converge */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup auxOTHERcomputational */

/* > \par Contributors: */
/*  ================== */
/* > */
/* > Jeff Rutter, Computer Science Division, University of California */
/* > at Berkeley, USA \n */
/* >  Modified by Francoise Tisseur, University of Tennessee */
/* > */
/*  ===================================================================== */
/* Subroutine */ int slaed3_(integer *k, integer *n, integer *n1, doublereal *
	d__, doublereal *q, integer *ldq, doublereal *rho, doublereal *dlamda,
	 doublereal *q2, integer *indx, integer *ctot, doublereal *w, 
	doublereal *s, integer *info)
{
    /* System generated locals */
    integer q_dim1, q_offset, i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal), d_sign(doublereal *, doublereal *);

    /* Local variables */
    static integer i__, j, n2, n12, ii, n23, iq2;
    static doublereal temp;
    extern doublereal snrm2_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int sgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen),
	     scopy_(integer *, doublereal *, integer *, doublereal *, integer 
	    *), slaed4_(integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *);
    extern doublereal slamc3_(doublereal *, doublereal *);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), slacpy_(
	    char *, integer *, integer *, doublereal *, integer *, doublereal 
	    *, integer *, ftnlen), slaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen);


/*  -- LAPACK computational routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

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

#line 227 "slaed3.f"
    /* Parameter adjustments */
#line 227 "slaed3.f"
    --d__;
#line 227 "slaed3.f"
    q_dim1 = *ldq;
#line 227 "slaed3.f"
    q_offset = 1 + q_dim1;
#line 227 "slaed3.f"
    q -= q_offset;
#line 227 "slaed3.f"
    --dlamda;
#line 227 "slaed3.f"
    --q2;
#line 227 "slaed3.f"
    --indx;
#line 227 "slaed3.f"
    --ctot;
#line 227 "slaed3.f"
    --w;
#line 227 "slaed3.f"
    --s;
#line 227 "slaed3.f"

#line 227 "slaed3.f"
    /* Function Body */
#line 227 "slaed3.f"
    *info = 0;

#line 229 "slaed3.f"
    if (*k < 0) {
#line 230 "slaed3.f"
	*info = -1;
#line 231 "slaed3.f"
    } else if (*n < *k) {
#line 232 "slaed3.f"
	*info = -2;
#line 233 "slaed3.f"
    } else if (*ldq < max(1,*n)) {
#line 234 "slaed3.f"
	*info = -6;
#line 235 "slaed3.f"
    }
#line 236 "slaed3.f"
    if (*info != 0) {
#line 237 "slaed3.f"
	i__1 = -(*info);
#line 237 "slaed3.f"
	xerbla_("SLAED3", &i__1, (ftnlen)6);
#line 238 "slaed3.f"
	return 0;
#line 239 "slaed3.f"
    }

/*     Quick return if possible */

#line 243 "slaed3.f"
    if (*k == 0) {
#line 243 "slaed3.f"
	return 0;
#line 243 "slaed3.f"
    }

/*     Modify values DLAMDA(i) to make sure all DLAMDA(i)-DLAMDA(j) can */
/*     be computed with high relative accuracy (barring over/underflow). */
/*     This is a problem on machines without a guard digit in */
/*     add/subtract (Cray XMP, Cray YMP, Cray C 90 and Cray 2). */
/*     The following code replaces DLAMDA(I) by 2*DLAMDA(I)-DLAMDA(I), */
/*     which on any of these machines zeros out the bottommost */
/*     bit of DLAMDA(I) if it is 1; this makes the subsequent */
/*     subtractions DLAMDA(I)-DLAMDA(J) unproblematic when cancellation */
/*     occurs. On binary machines with a guard digit (almost all */
/*     machines) it does not change DLAMDA(I) at all. On hexadecimal */
/*     and decimal machines with a guard digit, it slightly */
/*     changes the bottommost bits of DLAMDA(I). It does not account */
/*     for hexadecimal or decimal machines without guard digits */
/*     (we know of none). We use a subroutine call to compute */
/*     2*DLAMBDA(I) to prevent optimizing compilers from eliminating */
/*     this code. */

#line 263 "slaed3.f"
    i__1 = *k;
#line 263 "slaed3.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 264 "slaed3.f"
	dlamda[i__] = slamc3_(&dlamda[i__], &dlamda[i__]) - dlamda[i__];
#line 265 "slaed3.f"
/* L10: */
#line 265 "slaed3.f"
    }

#line 267 "slaed3.f"
    i__1 = *k;
#line 267 "slaed3.f"
    for (j = 1; j <= i__1; ++j) {
#line 268 "slaed3.f"
	slaed4_(k, &j, &dlamda[1], &w[1], &q[j * q_dim1 + 1], rho, &d__[j], 
		info);

/*        If the zero finder fails, the computation is terminated. */

#line 272 "slaed3.f"
	if (*info != 0) {
#line 272 "slaed3.f"
	    goto L120;
#line 272 "slaed3.f"
	}
#line 274 "slaed3.f"
/* L20: */
#line 274 "slaed3.f"
    }

#line 276 "slaed3.f"
    if (*k == 1) {
#line 276 "slaed3.f"
	goto L110;
#line 276 "slaed3.f"
    }
#line 278 "slaed3.f"
    if (*k == 2) {
#line 279 "slaed3.f"
	i__1 = *k;
#line 279 "slaed3.f"
	for (j = 1; j <= i__1; ++j) {
#line 280 "slaed3.f"
	    w[1] = q[j * q_dim1 + 1];
#line 281 "slaed3.f"
	    w[2] = q[j * q_dim1 + 2];
#line 282 "slaed3.f"
	    ii = indx[1];
#line 283 "slaed3.f"
	    q[j * q_dim1 + 1] = w[ii];
#line 284 "slaed3.f"
	    ii = indx[2];
#line 285 "slaed3.f"
	    q[j * q_dim1 + 2] = w[ii];
#line 286 "slaed3.f"
/* L30: */
#line 286 "slaed3.f"
	}
#line 287 "slaed3.f"
	goto L110;
#line 288 "slaed3.f"
    }

/*     Compute updated W. */

#line 292 "slaed3.f"
    scopy_(k, &w[1], &c__1, &s[1], &c__1);

/*     Initialize W(I) = Q(I,I) */

#line 296 "slaed3.f"
    i__1 = *ldq + 1;
#line 296 "slaed3.f"
    scopy_(k, &q[q_offset], &i__1, &w[1], &c__1);
#line 297 "slaed3.f"
    i__1 = *k;
#line 297 "slaed3.f"
    for (j = 1; j <= i__1; ++j) {
#line 298 "slaed3.f"
	i__2 = j - 1;
#line 298 "slaed3.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 299 "slaed3.f"
	    w[i__] *= q[i__ + j * q_dim1] / (dlamda[i__] - dlamda[j]);
#line 300 "slaed3.f"
/* L40: */
#line 300 "slaed3.f"
	}
#line 301 "slaed3.f"
	i__2 = *k;
#line 301 "slaed3.f"
	for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 302 "slaed3.f"
	    w[i__] *= q[i__ + j * q_dim1] / (dlamda[i__] - dlamda[j]);
#line 303 "slaed3.f"
/* L50: */
#line 303 "slaed3.f"
	}
#line 304 "slaed3.f"
/* L60: */
#line 304 "slaed3.f"
    }
#line 305 "slaed3.f"
    i__1 = *k;
#line 305 "slaed3.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 306 "slaed3.f"
	d__1 = sqrt(-w[i__]);
#line 306 "slaed3.f"
	w[i__] = d_sign(&d__1, &s[i__]);
#line 307 "slaed3.f"
/* L70: */
#line 307 "slaed3.f"
    }

/*     Compute eigenvectors of the modified rank-1 modification. */

#line 311 "slaed3.f"
    i__1 = *k;
#line 311 "slaed3.f"
    for (j = 1; j <= i__1; ++j) {
#line 312 "slaed3.f"
	i__2 = *k;
#line 312 "slaed3.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 313 "slaed3.f"
	    s[i__] = w[i__] / q[i__ + j * q_dim1];
#line 314 "slaed3.f"
/* L80: */
#line 314 "slaed3.f"
	}
#line 315 "slaed3.f"
	temp = snrm2_(k, &s[1], &c__1);
#line 316 "slaed3.f"
	i__2 = *k;
#line 316 "slaed3.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 317 "slaed3.f"
	    ii = indx[i__];
#line 318 "slaed3.f"
	    q[i__ + j * q_dim1] = s[ii] / temp;
#line 319 "slaed3.f"
/* L90: */
#line 319 "slaed3.f"
	}
#line 320 "slaed3.f"
/* L100: */
#line 320 "slaed3.f"
    }

/*     Compute the updated eigenvectors. */

#line 324 "slaed3.f"
L110:

#line 326 "slaed3.f"
    n2 = *n - *n1;
#line 327 "slaed3.f"
    n12 = ctot[1] + ctot[2];
#line 328 "slaed3.f"
    n23 = ctot[2] + ctot[3];

#line 330 "slaed3.f"
    slacpy_("A", &n23, k, &q[ctot[1] + 1 + q_dim1], ldq, &s[1], &n23, (ftnlen)
	    1);
#line 331 "slaed3.f"
    iq2 = *n1 * n12 + 1;
#line 332 "slaed3.f"
    if (n23 != 0) {
#line 333 "slaed3.f"
	sgemm_("N", "N", &n2, k, &n23, &c_b22, &q2[iq2], &n2, &s[1], &n23, &
		c_b23, &q[*n1 + 1 + q_dim1], ldq, (ftnlen)1, (ftnlen)1);
#line 335 "slaed3.f"
    } else {
#line 336 "slaed3.f"
	slaset_("A", &n2, k, &c_b23, &c_b23, &q[*n1 + 1 + q_dim1], ldq, (
		ftnlen)1);
#line 337 "slaed3.f"
    }

#line 339 "slaed3.f"
    slacpy_("A", &n12, k, &q[q_offset], ldq, &s[1], &n12, (ftnlen)1);
#line 340 "slaed3.f"
    if (n12 != 0) {
#line 341 "slaed3.f"
	sgemm_("N", "N", n1, k, &n12, &c_b22, &q2[1], n1, &s[1], &n12, &c_b23,
		 &q[q_offset], ldq, (ftnlen)1, (ftnlen)1);
#line 343 "slaed3.f"
    } else {
#line 344 "slaed3.f"
	slaset_("A", n1, k, &c_b23, &c_b23, &q[q_dim1 + 1], ldq, (ftnlen)1);
#line 345 "slaed3.f"
    }


#line 348 "slaed3.f"
L120:
#line 349 "slaed3.f"
    return 0;

/*     End of SLAED3 */

} /* slaed3_ */

