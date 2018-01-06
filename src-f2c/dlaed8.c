#line 1 "dlaed8.f"
/* dlaed8.f -- translated by f2c (version 20100827).
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

#line 1 "dlaed8.f"
/* Table of constant values */

static doublereal c_b3 = -1.;
static integer c__1 = 1;

/* > \brief \b DLAED8 used by sstedc. Merges eigenvalues and deflates secular equation. Used when the original
 matrix is dense. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLAED8 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaed8.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaed8.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaed8.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLAED8( ICOMPQ, K, N, QSIZ, D, Q, LDQ, INDXQ, RHO, */
/*                          CUTPNT, Z, DLAMDA, Q2, LDQ2, W, PERM, GIVPTR, */
/*                          GIVCOL, GIVNUM, INDXP, INDX, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            CUTPNT, GIVPTR, ICOMPQ, INFO, K, LDQ, LDQ2, N, */
/*      $                   QSIZ */
/*       DOUBLE PRECISION   RHO */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            GIVCOL( 2, * ), INDX( * ), INDXP( * ), */
/*      $                   INDXQ( * ), PERM( * ) */
/*       DOUBLE PRECISION   D( * ), DLAMDA( * ), GIVNUM( 2, * ), */
/*      $                   Q( LDQ, * ), Q2( LDQ2, * ), W( * ), Z( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLAED8 merges the two sets of eigenvalues together into a single */
/* > sorted set.  Then it tries to deflate the size of the problem. */
/* > There are two ways in which deflation can occur:  when two or more */
/* > eigenvalues are close together or if there is a tiny element in the */
/* > Z vector.  For each such occurrence the order of the related secular */
/* > equation problem is reduced by one. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] ICOMPQ */
/* > \verbatim */
/* >          ICOMPQ is INTEGER */
/* >          = 0:  Compute eigenvalues only. */
/* >          = 1:  Compute eigenvectors of original dense symmetric matrix */
/* >                also.  On entry, Q contains the orthogonal matrix used */
/* >                to reduce the original matrix to tridiagonal form. */
/* > \endverbatim */
/* > */
/* > \param[out] K */
/* > \verbatim */
/* >          K is INTEGER */
/* >         The number of non-deflated eigenvalues, and the order of the */
/* >         related secular equation. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >         The dimension of the symmetric tridiagonal matrix.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] QSIZ */
/* > \verbatim */
/* >          QSIZ is INTEGER */
/* >         The dimension of the orthogonal matrix used to reduce */
/* >         the full matrix to tridiagonal form.  QSIZ >= N if ICOMPQ = 1. */
/* > \endverbatim */
/* > */
/* > \param[in,out] D */
/* > \verbatim */
/* >          D is DOUBLE PRECISION array, dimension (N) */
/* >         On entry, the eigenvalues of the two submatrices to be */
/* >         combined.  On exit, the trailing (N-K) updated eigenvalues */
/* >         (those which were deflated) sorted into increasing order. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Q */
/* > \verbatim */
/* >          Q is DOUBLE PRECISION array, dimension (LDQ,N) */
/* >         If ICOMPQ = 0, Q is not referenced.  Otherwise, */
/* >         on entry, Q contains the eigenvectors of the partially solved */
/* >         system which has been previously updated in matrix */
/* >         multiplies with other partially solved eigensystems. */
/* >         On exit, Q contains the trailing (N-K) updated eigenvectors */
/* >         (those which were deflated) in its last N-K columns. */
/* > \endverbatim */
/* > */
/* > \param[in] LDQ */
/* > \verbatim */
/* >          LDQ is INTEGER */
/* >         The leading dimension of the array Q.  LDQ >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] INDXQ */
/* > \verbatim */
/* >          INDXQ is INTEGER array, dimension (N) */
/* >         The permutation which separately sorts the two sub-problems */
/* >         in D into ascending order.  Note that elements in the second */
/* >         half of this permutation must first have CUTPNT added to */
/* >         their values in order to be accurate. */
/* > \endverbatim */
/* > */
/* > \param[in,out] RHO */
/* > \verbatim */
/* >          RHO is DOUBLE PRECISION */
/* >         On entry, the off-diagonal element associated with the rank-1 */
/* >         cut which originally split the two submatrices which are now */
/* >         being recombined. */
/* >         On exit, RHO has been modified to the value required by */
/* >         DLAED3. */
/* > \endverbatim */
/* > */
/* > \param[in] CUTPNT */
/* > \verbatim */
/* >          CUTPNT is INTEGER */
/* >         The location of the last eigenvalue in the leading */
/* >         sub-matrix.  min(1,N) <= CUTPNT <= N. */
/* > \endverbatim */
/* > */
/* > \param[in] Z */
/* > \verbatim */
/* >          Z is DOUBLE PRECISION array, dimension (N) */
/* >         On entry, Z contains the updating vector (the last row of */
/* >         the first sub-eigenvector matrix and the first row of the */
/* >         second sub-eigenvector matrix). */
/* >         On exit, the contents of Z are destroyed by the updating */
/* >         process. */
/* > \endverbatim */
/* > */
/* > \param[out] DLAMDA */
/* > \verbatim */
/* >          DLAMDA is DOUBLE PRECISION array, dimension (N) */
/* >         A copy of the first K eigenvalues which will be used by */
/* >         DLAED3 to form the secular equation. */
/* > \endverbatim */
/* > */
/* > \param[out] Q2 */
/* > \verbatim */
/* >          Q2 is DOUBLE PRECISION array, dimension (LDQ2,N) */
/* >         If ICOMPQ = 0, Q2 is not referenced.  Otherwise, */
/* >         a copy of the first K eigenvectors which will be used by */
/* >         DLAED7 in a matrix multiply (DGEMM) to update the new */
/* >         eigenvectors. */
/* > \endverbatim */
/* > */
/* > \param[in] LDQ2 */
/* > \verbatim */
/* >          LDQ2 is INTEGER */
/* >         The leading dimension of the array Q2.  LDQ2 >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] W */
/* > \verbatim */
/* >          W is DOUBLE PRECISION array, dimension (N) */
/* >         The first k values of the final deflation-altered z-vector and */
/* >         will be passed to DLAED3. */
/* > \endverbatim */
/* > */
/* > \param[out] PERM */
/* > \verbatim */
/* >          PERM is INTEGER array, dimension (N) */
/* >         The permutations (from deflation and sorting) to be applied */
/* >         to each eigenblock. */
/* > \endverbatim */
/* > */
/* > \param[out] GIVPTR */
/* > \verbatim */
/* >          GIVPTR is INTEGER */
/* >         The number of Givens rotations which took place in this */
/* >         subproblem. */
/* > \endverbatim */
/* > */
/* > \param[out] GIVCOL */
/* > \verbatim */
/* >          GIVCOL is INTEGER array, dimension (2, N) */
/* >         Each pair of numbers indicates a pair of columns to take place */
/* >         in a Givens rotation. */
/* > \endverbatim */
/* > */
/* > \param[out] GIVNUM */
/* > \verbatim */
/* >          GIVNUM is DOUBLE PRECISION array, dimension (2, N) */
/* >         Each number indicates the S value to be used in the */
/* >         corresponding Givens rotation. */
/* > \endverbatim */
/* > */
/* > \param[out] INDXP */
/* > \verbatim */
/* >          INDXP is INTEGER array, dimension (N) */
/* >         The permutation used to place deflated values of D at the end */
/* >         of the array.  INDXP(1:K) points to the nondeflated D-values */
/* >         and INDXP(K+1:N) points to the deflated eigenvalues. */
/* > \endverbatim */
/* > */
/* > \param[out] INDX */
/* > \verbatim */
/* >          INDX is INTEGER array, dimension (N) */
/* >         The permutation used to sort the contents of D into ascending */
/* >         order. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit. */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup auxOTHERcomputational */

/* > \par Contributors: */
/*  ================== */
/* > */
/* > Jeff Rutter, Computer Science Division, University of California */
/* > at Berkeley, USA */

/*  ===================================================================== */
/* Subroutine */ int dlaed8_(integer *icompq, integer *k, integer *n, integer 
	*qsiz, doublereal *d__, doublereal *q, integer *ldq, integer *indxq, 
	doublereal *rho, integer *cutpnt, doublereal *z__, doublereal *dlamda,
	 doublereal *q2, integer *ldq2, doublereal *w, integer *perm, integer 
	*givptr, integer *givcol, doublereal *givnum, integer *indxp, integer 
	*indx, integer *info)
{
    /* System generated locals */
    integer q_dim1, q_offset, q2_dim1, q2_offset, i__1;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal c__;
    static integer i__, j;
    static doublereal s, t;
    static integer k2, n1, n2, jp, n1p1;
    static doublereal eps, tau, tol;
    static integer jlam, imax, jmax;
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *), dscal_(
	    integer *, doublereal *, doublereal *, integer *), dcopy_(integer 
	    *, doublereal *, integer *, doublereal *, integer *);
    extern doublereal dlapy2_(doublereal *, doublereal *), dlamch_(char *, 
	    ftnlen);
    extern integer idamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int dlamrg_(integer *, integer *, doublereal *, 
	    integer *, integer *, integer *), dlacpy_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen), xerbla_(char *, integer *, ftnlen);


/*  -- LAPACK computational routine (version 3.4.2) -- */
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

#line 290 "dlaed8.f"
    /* Parameter adjustments */
#line 290 "dlaed8.f"
    --d__;
#line 290 "dlaed8.f"
    q_dim1 = *ldq;
#line 290 "dlaed8.f"
    q_offset = 1 + q_dim1;
#line 290 "dlaed8.f"
    q -= q_offset;
#line 290 "dlaed8.f"
    --indxq;
#line 290 "dlaed8.f"
    --z__;
#line 290 "dlaed8.f"
    --dlamda;
#line 290 "dlaed8.f"
    q2_dim1 = *ldq2;
#line 290 "dlaed8.f"
    q2_offset = 1 + q2_dim1;
#line 290 "dlaed8.f"
    q2 -= q2_offset;
#line 290 "dlaed8.f"
    --w;
#line 290 "dlaed8.f"
    --perm;
#line 290 "dlaed8.f"
    givcol -= 3;
#line 290 "dlaed8.f"
    givnum -= 3;
#line 290 "dlaed8.f"
    --indxp;
#line 290 "dlaed8.f"
    --indx;
#line 290 "dlaed8.f"

#line 290 "dlaed8.f"
    /* Function Body */
#line 290 "dlaed8.f"
    *info = 0;

#line 292 "dlaed8.f"
    if (*icompq < 0 || *icompq > 1) {
#line 293 "dlaed8.f"
	*info = -1;
#line 294 "dlaed8.f"
    } else if (*n < 0) {
#line 295 "dlaed8.f"
	*info = -3;
#line 296 "dlaed8.f"
    } else if (*icompq == 1 && *qsiz < *n) {
#line 297 "dlaed8.f"
	*info = -4;
#line 298 "dlaed8.f"
    } else if (*ldq < max(1,*n)) {
#line 299 "dlaed8.f"
	*info = -7;
#line 300 "dlaed8.f"
    } else if (*cutpnt < min(1,*n) || *cutpnt > *n) {
#line 301 "dlaed8.f"
	*info = -10;
#line 302 "dlaed8.f"
    } else if (*ldq2 < max(1,*n)) {
#line 303 "dlaed8.f"
	*info = -14;
#line 304 "dlaed8.f"
    }
#line 305 "dlaed8.f"
    if (*info != 0) {
#line 306 "dlaed8.f"
	i__1 = -(*info);
#line 306 "dlaed8.f"
	xerbla_("DLAED8", &i__1, (ftnlen)6);
#line 307 "dlaed8.f"
	return 0;
#line 308 "dlaed8.f"
    }

/*     Need to initialize GIVPTR to O here in case of quick exit */
/*     to prevent an unspecified code behavior (usually sigfault) */
/*     when IWORK array on entry to *stedc is not zeroed */
/*     (or at least some IWORK entries which used in *laed7 for GIVPTR). */

#line 315 "dlaed8.f"
    *givptr = 0;

/*     Quick return if possible */

#line 319 "dlaed8.f"
    if (*n == 0) {
#line 319 "dlaed8.f"
	return 0;
#line 319 "dlaed8.f"
    }

#line 322 "dlaed8.f"
    n1 = *cutpnt;
#line 323 "dlaed8.f"
    n2 = *n - n1;
#line 324 "dlaed8.f"
    n1p1 = n1 + 1;

#line 326 "dlaed8.f"
    if (*rho < 0.) {
#line 327 "dlaed8.f"
	dscal_(&n2, &c_b3, &z__[n1p1], &c__1);
#line 328 "dlaed8.f"
    }

/*     Normalize z so that norm(z) = 1 */

#line 332 "dlaed8.f"
    t = 1. / sqrt(2.);
#line 333 "dlaed8.f"
    i__1 = *n;
#line 333 "dlaed8.f"
    for (j = 1; j <= i__1; ++j) {
#line 334 "dlaed8.f"
	indx[j] = j;
#line 335 "dlaed8.f"
/* L10: */
#line 335 "dlaed8.f"
    }
#line 336 "dlaed8.f"
    dscal_(n, &t, &z__[1], &c__1);
#line 337 "dlaed8.f"
    *rho = (d__1 = *rho * 2., abs(d__1));

/*     Sort the eigenvalues into increasing order */

#line 341 "dlaed8.f"
    i__1 = *n;
#line 341 "dlaed8.f"
    for (i__ = *cutpnt + 1; i__ <= i__1; ++i__) {
#line 342 "dlaed8.f"
	indxq[i__] += *cutpnt;
#line 343 "dlaed8.f"
/* L20: */
#line 343 "dlaed8.f"
    }
#line 344 "dlaed8.f"
    i__1 = *n;
#line 344 "dlaed8.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 345 "dlaed8.f"
	dlamda[i__] = d__[indxq[i__]];
#line 346 "dlaed8.f"
	w[i__] = z__[indxq[i__]];
#line 347 "dlaed8.f"
/* L30: */
#line 347 "dlaed8.f"
    }
#line 348 "dlaed8.f"
    i__ = 1;
#line 349 "dlaed8.f"
    j = *cutpnt + 1;
#line 350 "dlaed8.f"
    dlamrg_(&n1, &n2, &dlamda[1], &c__1, &c__1, &indx[1]);
#line 351 "dlaed8.f"
    i__1 = *n;
#line 351 "dlaed8.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 352 "dlaed8.f"
	d__[i__] = dlamda[indx[i__]];
#line 353 "dlaed8.f"
	z__[i__] = w[indx[i__]];
#line 354 "dlaed8.f"
/* L40: */
#line 354 "dlaed8.f"
    }

/*     Calculate the allowable deflation tolerence */

#line 358 "dlaed8.f"
    imax = idamax_(n, &z__[1], &c__1);
#line 359 "dlaed8.f"
    jmax = idamax_(n, &d__[1], &c__1);
#line 360 "dlaed8.f"
    eps = dlamch_("Epsilon", (ftnlen)7);
#line 361 "dlaed8.f"
    tol = eps * 8. * (d__1 = d__[jmax], abs(d__1));

/*     If the rank-1 modifier is small enough, no more needs to be done */
/*     except to reorganize Q so that its columns correspond with the */
/*     elements in D. */

#line 367 "dlaed8.f"
    if (*rho * (d__1 = z__[imax], abs(d__1)) <= tol) {
#line 368 "dlaed8.f"
	*k = 0;
#line 369 "dlaed8.f"
	if (*icompq == 0) {
#line 370 "dlaed8.f"
	    i__1 = *n;
#line 370 "dlaed8.f"
	    for (j = 1; j <= i__1; ++j) {
#line 371 "dlaed8.f"
		perm[j] = indxq[indx[j]];
#line 372 "dlaed8.f"
/* L50: */
#line 372 "dlaed8.f"
	    }
#line 373 "dlaed8.f"
	} else {
#line 374 "dlaed8.f"
	    i__1 = *n;
#line 374 "dlaed8.f"
	    for (j = 1; j <= i__1; ++j) {
#line 375 "dlaed8.f"
		perm[j] = indxq[indx[j]];
#line 376 "dlaed8.f"
		dcopy_(qsiz, &q[perm[j] * q_dim1 + 1], &c__1, &q2[j * q2_dim1 
			+ 1], &c__1);
#line 377 "dlaed8.f"
/* L60: */
#line 377 "dlaed8.f"
	    }
#line 378 "dlaed8.f"
	    dlacpy_("A", qsiz, n, &q2[q2_dim1 + 1], ldq2, &q[q_dim1 + 1], ldq,
		     (ftnlen)1);
#line 380 "dlaed8.f"
	}
#line 381 "dlaed8.f"
	return 0;
#line 382 "dlaed8.f"
    }

/*     If there are multiple eigenvalues then the problem deflates.  Here */
/*     the number of equal eigenvalues are found.  As each equal */
/*     eigenvalue is found, an elementary reflector is computed to rotate */
/*     the corresponding eigensubspace so that the corresponding */
/*     components of Z are zero in this new basis. */

#line 390 "dlaed8.f"
    *k = 0;
#line 391 "dlaed8.f"
    k2 = *n + 1;
#line 392 "dlaed8.f"
    i__1 = *n;
#line 392 "dlaed8.f"
    for (j = 1; j <= i__1; ++j) {
#line 393 "dlaed8.f"
	if (*rho * (d__1 = z__[j], abs(d__1)) <= tol) {

/*           Deflate due to small z component. */

#line 397 "dlaed8.f"
	    --k2;
#line 398 "dlaed8.f"
	    indxp[k2] = j;
#line 399 "dlaed8.f"
	    if (j == *n) {
#line 399 "dlaed8.f"
		goto L110;
#line 399 "dlaed8.f"
	    }
#line 401 "dlaed8.f"
	} else {
#line 402 "dlaed8.f"
	    jlam = j;
#line 403 "dlaed8.f"
	    goto L80;
#line 404 "dlaed8.f"
	}
#line 405 "dlaed8.f"
/* L70: */
#line 405 "dlaed8.f"
    }
#line 406 "dlaed8.f"
L80:
#line 407 "dlaed8.f"
    ++j;
#line 408 "dlaed8.f"
    if (j > *n) {
#line 408 "dlaed8.f"
	goto L100;
#line 408 "dlaed8.f"
    }
#line 410 "dlaed8.f"
    if (*rho * (d__1 = z__[j], abs(d__1)) <= tol) {

/*        Deflate due to small z component. */

#line 414 "dlaed8.f"
	--k2;
#line 415 "dlaed8.f"
	indxp[k2] = j;
#line 416 "dlaed8.f"
    } else {

/*        Check if eigenvalues are close enough to allow deflation. */

#line 420 "dlaed8.f"
	s = z__[jlam];
#line 421 "dlaed8.f"
	c__ = z__[j];

/*        Find sqrt(a**2+b**2) without overflow or */
/*        destructive underflow. */

#line 426 "dlaed8.f"
	tau = dlapy2_(&c__, &s);
#line 427 "dlaed8.f"
	t = d__[j] - d__[jlam];
#line 428 "dlaed8.f"
	c__ /= tau;
#line 429 "dlaed8.f"
	s = -s / tau;
#line 430 "dlaed8.f"
	if ((d__1 = t * c__ * s, abs(d__1)) <= tol) {

/*           Deflation is possible. */

#line 434 "dlaed8.f"
	    z__[j] = tau;
#line 435 "dlaed8.f"
	    z__[jlam] = 0.;

/*           Record the appropriate Givens rotation */

#line 439 "dlaed8.f"
	    ++(*givptr);
#line 440 "dlaed8.f"
	    givcol[(*givptr << 1) + 1] = indxq[indx[jlam]];
#line 441 "dlaed8.f"
	    givcol[(*givptr << 1) + 2] = indxq[indx[j]];
#line 442 "dlaed8.f"
	    givnum[(*givptr << 1) + 1] = c__;
#line 443 "dlaed8.f"
	    givnum[(*givptr << 1) + 2] = s;
#line 444 "dlaed8.f"
	    if (*icompq == 1) {
#line 445 "dlaed8.f"
		drot_(qsiz, &q[indxq[indx[jlam]] * q_dim1 + 1], &c__1, &q[
			indxq[indx[j]] * q_dim1 + 1], &c__1, &c__, &s);
#line 447 "dlaed8.f"
	    }
#line 448 "dlaed8.f"
	    t = d__[jlam] * c__ * c__ + d__[j] * s * s;
#line 449 "dlaed8.f"
	    d__[j] = d__[jlam] * s * s + d__[j] * c__ * c__;
#line 450 "dlaed8.f"
	    d__[jlam] = t;
#line 451 "dlaed8.f"
	    --k2;
#line 452 "dlaed8.f"
	    i__ = 1;
#line 453 "dlaed8.f"
L90:
#line 454 "dlaed8.f"
	    if (k2 + i__ <= *n) {
#line 455 "dlaed8.f"
		if (d__[jlam] < d__[indxp[k2 + i__]]) {
#line 456 "dlaed8.f"
		    indxp[k2 + i__ - 1] = indxp[k2 + i__];
#line 457 "dlaed8.f"
		    indxp[k2 + i__] = jlam;
#line 458 "dlaed8.f"
		    ++i__;
#line 459 "dlaed8.f"
		    goto L90;
#line 460 "dlaed8.f"
		} else {
#line 461 "dlaed8.f"
		    indxp[k2 + i__ - 1] = jlam;
#line 462 "dlaed8.f"
		}
#line 463 "dlaed8.f"
	    } else {
#line 464 "dlaed8.f"
		indxp[k2 + i__ - 1] = jlam;
#line 465 "dlaed8.f"
	    }
#line 466 "dlaed8.f"
	    jlam = j;
#line 467 "dlaed8.f"
	} else {
#line 468 "dlaed8.f"
	    ++(*k);
#line 469 "dlaed8.f"
	    w[*k] = z__[jlam];
#line 470 "dlaed8.f"
	    dlamda[*k] = d__[jlam];
#line 471 "dlaed8.f"
	    indxp[*k] = jlam;
#line 472 "dlaed8.f"
	    jlam = j;
#line 473 "dlaed8.f"
	}
#line 474 "dlaed8.f"
    }
#line 475 "dlaed8.f"
    goto L80;
#line 476 "dlaed8.f"
L100:

/*     Record the last eigenvalue. */

#line 480 "dlaed8.f"
    ++(*k);
#line 481 "dlaed8.f"
    w[*k] = z__[jlam];
#line 482 "dlaed8.f"
    dlamda[*k] = d__[jlam];
#line 483 "dlaed8.f"
    indxp[*k] = jlam;

#line 485 "dlaed8.f"
L110:

/*     Sort the eigenvalues and corresponding eigenvectors into DLAMDA */
/*     and Q2 respectively.  The eigenvalues/vectors which were not */
/*     deflated go into the first K slots of DLAMDA and Q2 respectively, */
/*     while those which were deflated go into the last N - K slots. */

#line 492 "dlaed8.f"
    if (*icompq == 0) {
#line 493 "dlaed8.f"
	i__1 = *n;
#line 493 "dlaed8.f"
	for (j = 1; j <= i__1; ++j) {
#line 494 "dlaed8.f"
	    jp = indxp[j];
#line 495 "dlaed8.f"
	    dlamda[j] = d__[jp];
#line 496 "dlaed8.f"
	    perm[j] = indxq[indx[jp]];
#line 497 "dlaed8.f"
/* L120: */
#line 497 "dlaed8.f"
	}
#line 498 "dlaed8.f"
    } else {
#line 499 "dlaed8.f"
	i__1 = *n;
#line 499 "dlaed8.f"
	for (j = 1; j <= i__1; ++j) {
#line 500 "dlaed8.f"
	    jp = indxp[j];
#line 501 "dlaed8.f"
	    dlamda[j] = d__[jp];
#line 502 "dlaed8.f"
	    perm[j] = indxq[indx[jp]];
#line 503 "dlaed8.f"
	    dcopy_(qsiz, &q[perm[j] * q_dim1 + 1], &c__1, &q2[j * q2_dim1 + 1]
		    , &c__1);
#line 504 "dlaed8.f"
/* L130: */
#line 504 "dlaed8.f"
	}
#line 505 "dlaed8.f"
    }

/*     The deflated eigenvalues and their corresponding vectors go back */
/*     into the last N - K slots of D and Q respectively. */

#line 510 "dlaed8.f"
    if (*k < *n) {
#line 511 "dlaed8.f"
	if (*icompq == 0) {
#line 512 "dlaed8.f"
	    i__1 = *n - *k;
#line 512 "dlaed8.f"
	    dcopy_(&i__1, &dlamda[*k + 1], &c__1, &d__[*k + 1], &c__1);
#line 513 "dlaed8.f"
	} else {
#line 514 "dlaed8.f"
	    i__1 = *n - *k;
#line 514 "dlaed8.f"
	    dcopy_(&i__1, &dlamda[*k + 1], &c__1, &d__[*k + 1], &c__1);
#line 515 "dlaed8.f"
	    i__1 = *n - *k;
#line 515 "dlaed8.f"
	    dlacpy_("A", qsiz, &i__1, &q2[(*k + 1) * q2_dim1 + 1], ldq2, &q[(*
		    k + 1) * q_dim1 + 1], ldq, (ftnlen)1);
#line 517 "dlaed8.f"
	}
#line 518 "dlaed8.f"
    }

#line 520 "dlaed8.f"
    return 0;

/*     End of DLAED8 */

} /* dlaed8_ */

