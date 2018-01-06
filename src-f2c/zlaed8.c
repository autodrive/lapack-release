#line 1 "zlaed8.f"
/* zlaed8.f -- translated by f2c (version 20100827).
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

#line 1 "zlaed8.f"
/* Table of constant values */

static doublereal c_b3 = -1.;
static integer c__1 = 1;

/* > \brief \b ZLAED8 used by sstedc. Merges eigenvalues and deflates secular equation. Used when the original
 matrix is dense. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZLAED8 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlaed8.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlaed8.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlaed8.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZLAED8( K, N, QSIZ, Q, LDQ, D, RHO, CUTPNT, Z, DLAMDA, */
/*                          Q2, LDQ2, W, INDXP, INDX, INDXQ, PERM, GIVPTR, */
/*                          GIVCOL, GIVNUM, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            CUTPNT, GIVPTR, INFO, K, LDQ, LDQ2, N, QSIZ */
/*       DOUBLE PRECISION   RHO */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            GIVCOL( 2, * ), INDX( * ), INDXP( * ), */
/*      $                   INDXQ( * ), PERM( * ) */
/*       DOUBLE PRECISION   D( * ), DLAMDA( * ), GIVNUM( 2, * ), W( * ), */
/*      $                   Z( * ) */
/*       COMPLEX*16         Q( LDQ, * ), Q2( LDQ2, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLAED8 merges the two sets of eigenvalues together into a single */
/* > sorted set.  Then it tries to deflate the size of the problem. */
/* > There are two ways in which deflation can occur:  when two or more */
/* > eigenvalues are close together or if there is a tiny element in the */
/* > Z vector.  For each such occurrence the order of the related secular */
/* > equation problem is reduced by one. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[out] K */
/* > \verbatim */
/* >          K is INTEGER */
/* >         Contains the number of non-deflated eigenvalues. */
/* >         This is the order of the related secular equation. */
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
/* >         The dimension of the unitary matrix used to reduce */
/* >         the dense or band matrix to tridiagonal form. */
/* >         QSIZ >= N if ICOMPQ = 1. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Q */
/* > \verbatim */
/* >          Q is COMPLEX*16 array, dimension (LDQ,N) */
/* >         On entry, Q contains the eigenvectors of the partially solved */
/* >         system which has been previously updated in matrix */
/* >         multiplies with other partially solved eigensystems. */
/* >         On exit, Q contains the trailing (N-K) updated eigenvectors */
/* >         (those which were deflated) in its last N-K columns. */
/* > \endverbatim */
/* > */
/* > \param[in] LDQ */
/* > \verbatim */
/* >          LDQ is INTEGER */
/* >         The leading dimension of the array Q.  LDQ >= max( 1, N ). */
/* > \endverbatim */
/* > */
/* > \param[in,out] D */
/* > \verbatim */
/* >          D is DOUBLE PRECISION array, dimension (N) */
/* >         On entry, D contains the eigenvalues of the two submatrices to */
/* >         be combined.  On exit, D contains the trailing (N-K) updated */
/* >         eigenvalues (those which were deflated) sorted into increasing */
/* >         order. */
/* > \endverbatim */
/* > */
/* > \param[in,out] RHO */
/* > \verbatim */
/* >          RHO is DOUBLE PRECISION */
/* >         Contains the off diagonal element associated with the rank-1 */
/* >         cut which originally split the two submatrices which are now */
/* >         being recombined. RHO is modified during the computation to */
/* >         the value required by DLAED3. */
/* > \endverbatim */
/* > */
/* > \param[in] CUTPNT */
/* > \verbatim */
/* >          CUTPNT is INTEGER */
/* >         Contains the location of the last eigenvalue in the leading */
/* >         sub-matrix.  MIN(1,N) <= CUTPNT <= N. */
/* > \endverbatim */
/* > */
/* > \param[in] Z */
/* > \verbatim */
/* >          Z is DOUBLE PRECISION array, dimension (N) */
/* >         On input this vector contains the updating vector (the last */
/* >         row of the first sub-eigenvector matrix and the first row of */
/* >         the second sub-eigenvector matrix).  The contents of Z are */
/* >         destroyed during the updating process. */
/* > \endverbatim */
/* > */
/* > \param[out] DLAMDA */
/* > \verbatim */
/* >          DLAMDA is DOUBLE PRECISION array, dimension (N) */
/* >         Contains a copy of the first K eigenvalues which will be used */
/* >         by DLAED3 to form the secular equation. */
/* > \endverbatim */
/* > */
/* > \param[out] Q2 */
/* > \verbatim */
/* >          Q2 is COMPLEX*16 array, dimension (LDQ2,N) */
/* >         If ICOMPQ = 0, Q2 is not referenced.  Otherwise, */
/* >         Contains a copy of the first K eigenvectors which will be used */
/* >         by DLAED7 in a matrix multiply (DGEMM) to update the new */
/* >         eigenvectors. */
/* > \endverbatim */
/* > */
/* > \param[in] LDQ2 */
/* > \verbatim */
/* >          LDQ2 is INTEGER */
/* >         The leading dimension of the array Q2.  LDQ2 >= max( 1, N ). */
/* > \endverbatim */
/* > */
/* > \param[out] W */
/* > \verbatim */
/* >          W is DOUBLE PRECISION array, dimension (N) */
/* >         This will hold the first k values of the final */
/* >         deflation-altered z-vector and will be passed to DLAED3. */
/* > \endverbatim */
/* > */
/* > \param[out] INDXP */
/* > \verbatim */
/* >          INDXP is INTEGER array, dimension (N) */
/* >         This will contain the permutation used to place deflated */
/* >         values of D at the end of the array. On output INDXP(1:K) */
/* >         points to the nondeflated D-values and INDXP(K+1:N) */
/* >         points to the deflated eigenvalues. */
/* > \endverbatim */
/* > */
/* > \param[out] INDX */
/* > \verbatim */
/* >          INDX is INTEGER array, dimension (N) */
/* >         This will contain the permutation used to sort the contents of */
/* >         D into ascending order. */
/* > \endverbatim */
/* > */
/* > \param[in] INDXQ */
/* > \verbatim */
/* >          INDXQ is INTEGER array, dimension (N) */
/* >         This contains the permutation which separately sorts the two */
/* >         sub-problems in D into ascending order.  Note that elements in */
/* >         the second half of this permutation must first have CUTPNT */
/* >         added to their values in order to be accurate. */
/* > \endverbatim */
/* > */
/* > \param[out] PERM */
/* > \verbatim */
/* >          PERM is INTEGER array, dimension (N) */
/* >         Contains the permutations (from deflation and sorting) to be */
/* >         applied to each eigenblock. */
/* > \endverbatim */
/* > */
/* > \param[out] GIVPTR */
/* > \verbatim */
/* >          GIVPTR is INTEGER */
/* >         Contains the number of Givens rotations which took place in */
/* >         this subproblem. */
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

/* > \date December 2016 */

/* > \ingroup complex16OTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int zlaed8_(integer *k, integer *n, integer *qsiz, 
	doublecomplex *q, integer *ldq, doublereal *d__, doublereal *rho, 
	integer *cutpnt, doublereal *z__, doublereal *dlamda, doublecomplex *
	q2, integer *ldq2, doublereal *w, integer *indxp, integer *indx, 
	integer *indxq, integer *perm, integer *givptr, integer *givcol, 
	doublereal *givnum, integer *info)
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
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), dcopy_(integer *, doublereal *, integer *, doublereal 
	    *, integer *), zdrot_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublereal *, doublereal *), zcopy_(
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *)
	    ;
    extern doublereal dlapy2_(doublereal *, doublereal *), dlamch_(char *, 
	    ftnlen);
    extern integer idamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int dlamrg_(integer *, integer *, doublereal *, 
	    integer *, integer *, integer *), xerbla_(char *, integer *, 
	    ftnlen), zlacpy_(char *, integer *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, ftnlen);


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

#line 275 "zlaed8.f"
    /* Parameter adjustments */
#line 275 "zlaed8.f"
    q_dim1 = *ldq;
#line 275 "zlaed8.f"
    q_offset = 1 + q_dim1;
#line 275 "zlaed8.f"
    q -= q_offset;
#line 275 "zlaed8.f"
    --d__;
#line 275 "zlaed8.f"
    --z__;
#line 275 "zlaed8.f"
    --dlamda;
#line 275 "zlaed8.f"
    q2_dim1 = *ldq2;
#line 275 "zlaed8.f"
    q2_offset = 1 + q2_dim1;
#line 275 "zlaed8.f"
    q2 -= q2_offset;
#line 275 "zlaed8.f"
    --w;
#line 275 "zlaed8.f"
    --indxp;
#line 275 "zlaed8.f"
    --indx;
#line 275 "zlaed8.f"
    --indxq;
#line 275 "zlaed8.f"
    --perm;
#line 275 "zlaed8.f"
    givcol -= 3;
#line 275 "zlaed8.f"
    givnum -= 3;
#line 275 "zlaed8.f"

#line 275 "zlaed8.f"
    /* Function Body */
#line 275 "zlaed8.f"
    *info = 0;

#line 277 "zlaed8.f"
    if (*n < 0) {
#line 278 "zlaed8.f"
	*info = -2;
#line 279 "zlaed8.f"
    } else if (*qsiz < *n) {
#line 280 "zlaed8.f"
	*info = -3;
#line 281 "zlaed8.f"
    } else if (*ldq < max(1,*n)) {
#line 282 "zlaed8.f"
	*info = -5;
#line 283 "zlaed8.f"
    } else if (*cutpnt < min(1,*n) || *cutpnt > *n) {
#line 284 "zlaed8.f"
	*info = -8;
#line 285 "zlaed8.f"
    } else if (*ldq2 < max(1,*n)) {
#line 286 "zlaed8.f"
	*info = -12;
#line 287 "zlaed8.f"
    }
#line 288 "zlaed8.f"
    if (*info != 0) {
#line 289 "zlaed8.f"
	i__1 = -(*info);
#line 289 "zlaed8.f"
	xerbla_("ZLAED8", &i__1, (ftnlen)6);
#line 290 "zlaed8.f"
	return 0;
#line 291 "zlaed8.f"
    }

/*     Need to initialize GIVPTR to O here in case of quick exit */
/*     to prevent an unspecified code behavior (usually sigfault) */
/*     when IWORK array on entry to *stedc is not zeroed */
/*     (or at least some IWORK entries which used in *laed7 for GIVPTR). */

#line 298 "zlaed8.f"
    *givptr = 0;

/*     Quick return if possible */

#line 302 "zlaed8.f"
    if (*n == 0) {
#line 302 "zlaed8.f"
	return 0;
#line 302 "zlaed8.f"
    }

#line 305 "zlaed8.f"
    n1 = *cutpnt;
#line 306 "zlaed8.f"
    n2 = *n - n1;
#line 307 "zlaed8.f"
    n1p1 = n1 + 1;

#line 309 "zlaed8.f"
    if (*rho < 0.) {
#line 310 "zlaed8.f"
	dscal_(&n2, &c_b3, &z__[n1p1], &c__1);
#line 311 "zlaed8.f"
    }

/*     Normalize z so that norm(z) = 1 */

#line 315 "zlaed8.f"
    t = 1. / sqrt(2.);
#line 316 "zlaed8.f"
    i__1 = *n;
#line 316 "zlaed8.f"
    for (j = 1; j <= i__1; ++j) {
#line 317 "zlaed8.f"
	indx[j] = j;
#line 318 "zlaed8.f"
/* L10: */
#line 318 "zlaed8.f"
    }
#line 319 "zlaed8.f"
    dscal_(n, &t, &z__[1], &c__1);
#line 320 "zlaed8.f"
    *rho = (d__1 = *rho * 2., abs(d__1));

/*     Sort the eigenvalues into increasing order */

#line 324 "zlaed8.f"
    i__1 = *n;
#line 324 "zlaed8.f"
    for (i__ = *cutpnt + 1; i__ <= i__1; ++i__) {
#line 325 "zlaed8.f"
	indxq[i__] += *cutpnt;
#line 326 "zlaed8.f"
/* L20: */
#line 326 "zlaed8.f"
    }
#line 327 "zlaed8.f"
    i__1 = *n;
#line 327 "zlaed8.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 328 "zlaed8.f"
	dlamda[i__] = d__[indxq[i__]];
#line 329 "zlaed8.f"
	w[i__] = z__[indxq[i__]];
#line 330 "zlaed8.f"
/* L30: */
#line 330 "zlaed8.f"
    }
#line 331 "zlaed8.f"
    i__ = 1;
#line 332 "zlaed8.f"
    j = *cutpnt + 1;
#line 333 "zlaed8.f"
    dlamrg_(&n1, &n2, &dlamda[1], &c__1, &c__1, &indx[1]);
#line 334 "zlaed8.f"
    i__1 = *n;
#line 334 "zlaed8.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 335 "zlaed8.f"
	d__[i__] = dlamda[indx[i__]];
#line 336 "zlaed8.f"
	z__[i__] = w[indx[i__]];
#line 337 "zlaed8.f"
/* L40: */
#line 337 "zlaed8.f"
    }

/*     Calculate the allowable deflation tolerance */

#line 341 "zlaed8.f"
    imax = idamax_(n, &z__[1], &c__1);
#line 342 "zlaed8.f"
    jmax = idamax_(n, &d__[1], &c__1);
#line 343 "zlaed8.f"
    eps = dlamch_("Epsilon", (ftnlen)7);
#line 344 "zlaed8.f"
    tol = eps * 8. * (d__1 = d__[jmax], abs(d__1));

/*     If the rank-1 modifier is small enough, no more needs to be done */
/*     -- except to reorganize Q so that its columns correspond with the */
/*     elements in D. */

#line 350 "zlaed8.f"
    if (*rho * (d__1 = z__[imax], abs(d__1)) <= tol) {
#line 351 "zlaed8.f"
	*k = 0;
#line 352 "zlaed8.f"
	i__1 = *n;
#line 352 "zlaed8.f"
	for (j = 1; j <= i__1; ++j) {
#line 353 "zlaed8.f"
	    perm[j] = indxq[indx[j]];
#line 354 "zlaed8.f"
	    zcopy_(qsiz, &q[perm[j] * q_dim1 + 1], &c__1, &q2[j * q2_dim1 + 1]
		    , &c__1);
#line 355 "zlaed8.f"
/* L50: */
#line 355 "zlaed8.f"
	}
#line 356 "zlaed8.f"
	zlacpy_("A", qsiz, n, &q2[q2_dim1 + 1], ldq2, &q[q_dim1 + 1], ldq, (
		ftnlen)1);
#line 357 "zlaed8.f"
	return 0;
#line 358 "zlaed8.f"
    }

/*     If there are multiple eigenvalues then the problem deflates.  Here */
/*     the number of equal eigenvalues are found.  As each equal */
/*     eigenvalue is found, an elementary reflector is computed to rotate */
/*     the corresponding eigensubspace so that the corresponding */
/*     components of Z are zero in this new basis. */

#line 366 "zlaed8.f"
    *k = 0;
#line 367 "zlaed8.f"
    k2 = *n + 1;
#line 368 "zlaed8.f"
    i__1 = *n;
#line 368 "zlaed8.f"
    for (j = 1; j <= i__1; ++j) {
#line 369 "zlaed8.f"
	if (*rho * (d__1 = z__[j], abs(d__1)) <= tol) {

/*           Deflate due to small z component. */

#line 373 "zlaed8.f"
	    --k2;
#line 374 "zlaed8.f"
	    indxp[k2] = j;
#line 375 "zlaed8.f"
	    if (j == *n) {
#line 375 "zlaed8.f"
		goto L100;
#line 375 "zlaed8.f"
	    }
#line 377 "zlaed8.f"
	} else {
#line 378 "zlaed8.f"
	    jlam = j;
#line 379 "zlaed8.f"
	    goto L70;
#line 380 "zlaed8.f"
	}
#line 381 "zlaed8.f"
/* L60: */
#line 381 "zlaed8.f"
    }
#line 382 "zlaed8.f"
L70:
#line 383 "zlaed8.f"
    ++j;
#line 384 "zlaed8.f"
    if (j > *n) {
#line 384 "zlaed8.f"
	goto L90;
#line 384 "zlaed8.f"
    }
#line 386 "zlaed8.f"
    if (*rho * (d__1 = z__[j], abs(d__1)) <= tol) {

/*        Deflate due to small z component. */

#line 390 "zlaed8.f"
	--k2;
#line 391 "zlaed8.f"
	indxp[k2] = j;
#line 392 "zlaed8.f"
    } else {

/*        Check if eigenvalues are close enough to allow deflation. */

#line 396 "zlaed8.f"
	s = z__[jlam];
#line 397 "zlaed8.f"
	c__ = z__[j];

/*        Find sqrt(a**2+b**2) without overflow or */
/*        destructive underflow. */

#line 402 "zlaed8.f"
	tau = dlapy2_(&c__, &s);
#line 403 "zlaed8.f"
	t = d__[j] - d__[jlam];
#line 404 "zlaed8.f"
	c__ /= tau;
#line 405 "zlaed8.f"
	s = -s / tau;
#line 406 "zlaed8.f"
	if ((d__1 = t * c__ * s, abs(d__1)) <= tol) {

/*           Deflation is possible. */

#line 410 "zlaed8.f"
	    z__[j] = tau;
#line 411 "zlaed8.f"
	    z__[jlam] = 0.;

/*           Record the appropriate Givens rotation */

#line 415 "zlaed8.f"
	    ++(*givptr);
#line 416 "zlaed8.f"
	    givcol[(*givptr << 1) + 1] = indxq[indx[jlam]];
#line 417 "zlaed8.f"
	    givcol[(*givptr << 1) + 2] = indxq[indx[j]];
#line 418 "zlaed8.f"
	    givnum[(*givptr << 1) + 1] = c__;
#line 419 "zlaed8.f"
	    givnum[(*givptr << 1) + 2] = s;
#line 420 "zlaed8.f"
	    zdrot_(qsiz, &q[indxq[indx[jlam]] * q_dim1 + 1], &c__1, &q[indxq[
		    indx[j]] * q_dim1 + 1], &c__1, &c__, &s);
#line 422 "zlaed8.f"
	    t = d__[jlam] * c__ * c__ + d__[j] * s * s;
#line 423 "zlaed8.f"
	    d__[j] = d__[jlam] * s * s + d__[j] * c__ * c__;
#line 424 "zlaed8.f"
	    d__[jlam] = t;
#line 425 "zlaed8.f"
	    --k2;
#line 426 "zlaed8.f"
	    i__ = 1;
#line 427 "zlaed8.f"
L80:
#line 428 "zlaed8.f"
	    if (k2 + i__ <= *n) {
#line 429 "zlaed8.f"
		if (d__[jlam] < d__[indxp[k2 + i__]]) {
#line 430 "zlaed8.f"
		    indxp[k2 + i__ - 1] = indxp[k2 + i__];
#line 431 "zlaed8.f"
		    indxp[k2 + i__] = jlam;
#line 432 "zlaed8.f"
		    ++i__;
#line 433 "zlaed8.f"
		    goto L80;
#line 434 "zlaed8.f"
		} else {
#line 435 "zlaed8.f"
		    indxp[k2 + i__ - 1] = jlam;
#line 436 "zlaed8.f"
		}
#line 437 "zlaed8.f"
	    } else {
#line 438 "zlaed8.f"
		indxp[k2 + i__ - 1] = jlam;
#line 439 "zlaed8.f"
	    }
#line 440 "zlaed8.f"
	    jlam = j;
#line 441 "zlaed8.f"
	} else {
#line 442 "zlaed8.f"
	    ++(*k);
#line 443 "zlaed8.f"
	    w[*k] = z__[jlam];
#line 444 "zlaed8.f"
	    dlamda[*k] = d__[jlam];
#line 445 "zlaed8.f"
	    indxp[*k] = jlam;
#line 446 "zlaed8.f"
	    jlam = j;
#line 447 "zlaed8.f"
	}
#line 448 "zlaed8.f"
    }
#line 449 "zlaed8.f"
    goto L70;
#line 450 "zlaed8.f"
L90:

/*     Record the last eigenvalue. */

#line 454 "zlaed8.f"
    ++(*k);
#line 455 "zlaed8.f"
    w[*k] = z__[jlam];
#line 456 "zlaed8.f"
    dlamda[*k] = d__[jlam];
#line 457 "zlaed8.f"
    indxp[*k] = jlam;

#line 459 "zlaed8.f"
L100:

/*     Sort the eigenvalues and corresponding eigenvectors into DLAMDA */
/*     and Q2 respectively.  The eigenvalues/vectors which were not */
/*     deflated go into the first K slots of DLAMDA and Q2 respectively, */
/*     while those which were deflated go into the last N - K slots. */

#line 466 "zlaed8.f"
    i__1 = *n;
#line 466 "zlaed8.f"
    for (j = 1; j <= i__1; ++j) {
#line 467 "zlaed8.f"
	jp = indxp[j];
#line 468 "zlaed8.f"
	dlamda[j] = d__[jp];
#line 469 "zlaed8.f"
	perm[j] = indxq[indx[jp]];
#line 470 "zlaed8.f"
	zcopy_(qsiz, &q[perm[j] * q_dim1 + 1], &c__1, &q2[j * q2_dim1 + 1], &
		c__1);
#line 471 "zlaed8.f"
/* L110: */
#line 471 "zlaed8.f"
    }

/*     The deflated eigenvalues and their corresponding vectors go back */
/*     into the last N - K slots of D and Q respectively. */

#line 476 "zlaed8.f"
    if (*k < *n) {
#line 477 "zlaed8.f"
	i__1 = *n - *k;
#line 477 "zlaed8.f"
	dcopy_(&i__1, &dlamda[*k + 1], &c__1, &d__[*k + 1], &c__1);
#line 478 "zlaed8.f"
	i__1 = *n - *k;
#line 478 "zlaed8.f"
	zlacpy_("A", qsiz, &i__1, &q2[(*k + 1) * q2_dim1 + 1], ldq2, &q[(*k + 
		1) * q_dim1 + 1], ldq, (ftnlen)1);
#line 480 "zlaed8.f"
    }

#line 482 "zlaed8.f"
    return 0;

/*     End of ZLAED8 */

} /* zlaed8_ */

