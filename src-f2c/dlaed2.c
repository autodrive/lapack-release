#line 1 "dlaed2.f"
/* dlaed2.f -- translated by f2c (version 20100827).
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

#line 1 "dlaed2.f"
/* Table of constant values */

static doublereal c_b3 = -1.;
static integer c__1 = 1;

/* > \brief \b DLAED2 used by sstedc. Merges eigenvalues and deflates secular equation. Used when the original
 matrix is tridiagonal. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLAED2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaed2.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaed2.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaed2.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLAED2( K, N, N1, D, Q, LDQ, INDXQ, RHO, Z, DLAMDA, W, */
/*                          Q2, INDX, INDXC, INDXP, COLTYP, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, K, LDQ, N, N1 */
/*       DOUBLE PRECISION   RHO */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            COLTYP( * ), INDX( * ), INDXC( * ), INDXP( * ), */
/*      $                   INDXQ( * ) */
/*       DOUBLE PRECISION   D( * ), DLAMDA( * ), Q( LDQ, * ), Q2( * ), */
/*      $                   W( * ), Z( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLAED2 merges the two sets of eigenvalues together into a single */
/* > sorted set.  Then it tries to deflate the size of the problem. */
/* > There are two ways in which deflation can occur:  when two or more */
/* > eigenvalues are close together or if there is a tiny entry in the */
/* > Z vector.  For each such occurrence the order of the related secular */
/* > equation problem is reduced by one. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[out] K */
/* > \verbatim */
/* >          K is INTEGER */
/* >         The number of non-deflated eigenvalues, and the order of the */
/* >         related secular equation. 0 <= K <=N. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >         The dimension of the symmetric tridiagonal matrix.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N1 */
/* > \verbatim */
/* >          N1 is INTEGER */
/* >         The location of the last eigenvalue in the leading sub-matrix. */
/* >         min(1,N) <= N1 <= N/2. */
/* > \endverbatim */
/* > */
/* > \param[in,out] D */
/* > \verbatim */
/* >          D is DOUBLE PRECISION array, dimension (N) */
/* >         On entry, D contains the eigenvalues of the two submatrices to */
/* >         be combined. */
/* >         On exit, D contains the trailing (N-K) updated eigenvalues */
/* >         (those which were deflated) sorted into increasing order. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Q */
/* > \verbatim */
/* >          Q is DOUBLE PRECISION array, dimension (LDQ, N) */
/* >         On entry, Q contains the eigenvectors of two submatrices in */
/* >         the two square blocks with corners at (1,1), (N1,N1) */
/* >         and (N1+1, N1+1), (N,N). */
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
/* > \param[in,out] INDXQ */
/* > \verbatim */
/* >          INDXQ is INTEGER array, dimension (N) */
/* >         The permutation which separately sorts the two sub-problems */
/* >         in D into ascending order.  Note that elements in the second */
/* >         half of this permutation must first have N1 added to their */
/* >         values. Destroyed on exit. */
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
/* > \param[in] Z */
/* > \verbatim */
/* >          Z is DOUBLE PRECISION array, dimension (N) */
/* >         On entry, Z contains the updating vector (the last */
/* >         row of the first sub-eigenvector matrix and the first row of */
/* >         the second sub-eigenvector matrix). */
/* >         On exit, the contents of Z have been destroyed by the updating */
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
/* > \param[out] W */
/* > \verbatim */
/* >          W is DOUBLE PRECISION array, dimension (N) */
/* >         The first k values of the final deflation-altered z-vector */
/* >         which will be passed to DLAED3. */
/* > \endverbatim */
/* > */
/* > \param[out] Q2 */
/* > \verbatim */
/* >          Q2 is DOUBLE PRECISION array, dimension (N1**2+(N-N1)**2) */
/* >         A copy of the first K eigenvectors which will be used by */
/* >         DLAED3 in a matrix multiply (DGEMM) to solve for the new */
/* >         eigenvectors. */
/* > \endverbatim */
/* > */
/* > \param[out] INDX */
/* > \verbatim */
/* >          INDX is INTEGER array, dimension (N) */
/* >         The permutation used to sort the contents of DLAMDA into */
/* >         ascending order. */
/* > \endverbatim */
/* > */
/* > \param[out] INDXC */
/* > \verbatim */
/* >          INDXC is INTEGER array, dimension (N) */
/* >         The permutation used to arrange the columns of the deflated */
/* >         Q matrix into three groups:  the first group contains non-zero */
/* >         elements only at and above N1, the second contains */
/* >         non-zero elements only below N1, and the third is dense. */
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
/* > \param[out] COLTYP */
/* > \verbatim */
/* >          COLTYP is INTEGER array, dimension (N) */
/* >         During execution, a label which will indicate which of the */
/* >         following types a column in the Q2 matrix is: */
/* >         1 : non-zero in the upper half only; */
/* >         2 : dense; */
/* >         3 : non-zero in the lower half only; */
/* >         4 : deflated. */
/* >         On exit, COLTYP(i) is the number of columns of type i, */
/* >         for i=1 to 4 only. */
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

/* > \ingroup auxOTHERcomputational */

/* > \par Contributors: */
/*  ================== */
/* > */
/* > Jeff Rutter, Computer Science Division, University of California */
/* > at Berkeley, USA \n */
/* >  Modified by Francoise Tisseur, University of Tennessee */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dlaed2_(integer *k, integer *n, integer *n1, doublereal *
	d__, doublereal *q, integer *ldq, integer *indxq, doublereal *rho, 
	doublereal *z__, doublereal *dlamda, doublereal *w, doublereal *q2, 
	integer *indx, integer *indxc, integer *indxp, integer *coltyp, 
	integer *info)
{
    /* System generated locals */
    integer q_dim1, q_offset, i__1, i__2;
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal c__;
    static integer i__, j;
    static doublereal s, t;
    static integer k2, n2, ct, nj, pj, js, iq1, iq2, n1p1;
    static doublereal eps, tau, tol;
    static integer psm[4], imax, jmax;
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *);
    static integer ctot[4];
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), dcopy_(integer *, doublereal *, integer *, doublereal 
	    *, integer *);
    extern doublereal dlapy2_(doublereal *, doublereal *), dlamch_(char *, 
	    ftnlen);
    extern integer idamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int dlamrg_(integer *, integer *, doublereal *, 
	    integer *, integer *, integer *), dlacpy_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen), xerbla_(char *, integer *, ftnlen);


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
/*     .. Local Arrays .. */
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

#line 261 "dlaed2.f"
    /* Parameter adjustments */
#line 261 "dlaed2.f"
    --d__;
#line 261 "dlaed2.f"
    q_dim1 = *ldq;
#line 261 "dlaed2.f"
    q_offset = 1 + q_dim1;
#line 261 "dlaed2.f"
    q -= q_offset;
#line 261 "dlaed2.f"
    --indxq;
#line 261 "dlaed2.f"
    --z__;
#line 261 "dlaed2.f"
    --dlamda;
#line 261 "dlaed2.f"
    --w;
#line 261 "dlaed2.f"
    --q2;
#line 261 "dlaed2.f"
    --indx;
#line 261 "dlaed2.f"
    --indxc;
#line 261 "dlaed2.f"
    --indxp;
#line 261 "dlaed2.f"
    --coltyp;
#line 261 "dlaed2.f"

#line 261 "dlaed2.f"
    /* Function Body */
#line 261 "dlaed2.f"
    *info = 0;

#line 263 "dlaed2.f"
    if (*n < 0) {
#line 264 "dlaed2.f"
	*info = -2;
#line 265 "dlaed2.f"
    } else if (*ldq < max(1,*n)) {
#line 266 "dlaed2.f"
	*info = -6;
#line 267 "dlaed2.f"
    } else /* if(complicated condition) */ {
/* Computing MIN */
#line 267 "dlaed2.f"
	i__1 = 1, i__2 = *n / 2;
#line 267 "dlaed2.f"
	if (min(i__1,i__2) > *n1 || *n / 2 < *n1) {
#line 268 "dlaed2.f"
	    *info = -3;
#line 269 "dlaed2.f"
	}
#line 269 "dlaed2.f"
    }
#line 270 "dlaed2.f"
    if (*info != 0) {
#line 271 "dlaed2.f"
	i__1 = -(*info);
#line 271 "dlaed2.f"
	xerbla_("DLAED2", &i__1, (ftnlen)6);
#line 272 "dlaed2.f"
	return 0;
#line 273 "dlaed2.f"
    }

/*     Quick return if possible */

#line 277 "dlaed2.f"
    if (*n == 0) {
#line 277 "dlaed2.f"
	return 0;
#line 277 "dlaed2.f"
    }

#line 280 "dlaed2.f"
    n2 = *n - *n1;
#line 281 "dlaed2.f"
    n1p1 = *n1 + 1;

#line 283 "dlaed2.f"
    if (*rho < 0.) {
#line 284 "dlaed2.f"
	dscal_(&n2, &c_b3, &z__[n1p1], &c__1);
#line 285 "dlaed2.f"
    }

/*     Normalize z so that norm(z) = 1.  Since z is the concatenation of */
/*     two normalized vectors, norm2(z) = sqrt(2). */

#line 290 "dlaed2.f"
    t = 1. / sqrt(2.);
#line 291 "dlaed2.f"
    dscal_(n, &t, &z__[1], &c__1);

/*     RHO = ABS( norm(z)**2 * RHO ) */

#line 295 "dlaed2.f"
    *rho = (d__1 = *rho * 2., abs(d__1));

/*     Sort the eigenvalues into increasing order */

#line 299 "dlaed2.f"
    i__1 = *n;
#line 299 "dlaed2.f"
    for (i__ = n1p1; i__ <= i__1; ++i__) {
#line 300 "dlaed2.f"
	indxq[i__] += *n1;
#line 301 "dlaed2.f"
/* L10: */
#line 301 "dlaed2.f"
    }

/*     re-integrate the deflated parts from the last pass */

#line 305 "dlaed2.f"
    i__1 = *n;
#line 305 "dlaed2.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 306 "dlaed2.f"
	dlamda[i__] = d__[indxq[i__]];
#line 307 "dlaed2.f"
/* L20: */
#line 307 "dlaed2.f"
    }
#line 308 "dlaed2.f"
    dlamrg_(n1, &n2, &dlamda[1], &c__1, &c__1, &indxc[1]);
#line 309 "dlaed2.f"
    i__1 = *n;
#line 309 "dlaed2.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 310 "dlaed2.f"
	indx[i__] = indxq[indxc[i__]];
#line 311 "dlaed2.f"
/* L30: */
#line 311 "dlaed2.f"
    }

/*     Calculate the allowable deflation tolerance */

#line 315 "dlaed2.f"
    imax = idamax_(n, &z__[1], &c__1);
#line 316 "dlaed2.f"
    jmax = idamax_(n, &d__[1], &c__1);
#line 317 "dlaed2.f"
    eps = dlamch_("Epsilon", (ftnlen)7);
/* Computing MAX */
#line 318 "dlaed2.f"
    d__3 = (d__1 = d__[jmax], abs(d__1)), d__4 = (d__2 = z__[imax], abs(d__2))
	    ;
#line 318 "dlaed2.f"
    tol = eps * 8. * max(d__3,d__4);

/*     If the rank-1 modifier is small enough, no more needs to be done */
/*     except to reorganize Q so that its columns correspond with the */
/*     elements in D. */

#line 324 "dlaed2.f"
    if (*rho * (d__1 = z__[imax], abs(d__1)) <= tol) {
#line 325 "dlaed2.f"
	*k = 0;
#line 326 "dlaed2.f"
	iq2 = 1;
#line 327 "dlaed2.f"
	i__1 = *n;
#line 327 "dlaed2.f"
	for (j = 1; j <= i__1; ++j) {
#line 328 "dlaed2.f"
	    i__ = indx[j];
#line 329 "dlaed2.f"
	    dcopy_(n, &q[i__ * q_dim1 + 1], &c__1, &q2[iq2], &c__1);
#line 330 "dlaed2.f"
	    dlamda[j] = d__[i__];
#line 331 "dlaed2.f"
	    iq2 += *n;
#line 332 "dlaed2.f"
/* L40: */
#line 332 "dlaed2.f"
	}
#line 333 "dlaed2.f"
	dlacpy_("A", n, n, &q2[1], n, &q[q_offset], ldq, (ftnlen)1);
#line 334 "dlaed2.f"
	dcopy_(n, &dlamda[1], &c__1, &d__[1], &c__1);
#line 335 "dlaed2.f"
	goto L190;
#line 336 "dlaed2.f"
    }

/*     If there are multiple eigenvalues then the problem deflates.  Here */
/*     the number of equal eigenvalues are found.  As each equal */
/*     eigenvalue is found, an elementary reflector is computed to rotate */
/*     the corresponding eigensubspace so that the corresponding */
/*     components of Z are zero in this new basis. */

#line 344 "dlaed2.f"
    i__1 = *n1;
#line 344 "dlaed2.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 345 "dlaed2.f"
	coltyp[i__] = 1;
#line 346 "dlaed2.f"
/* L50: */
#line 346 "dlaed2.f"
    }
#line 347 "dlaed2.f"
    i__1 = *n;
#line 347 "dlaed2.f"
    for (i__ = n1p1; i__ <= i__1; ++i__) {
#line 348 "dlaed2.f"
	coltyp[i__] = 3;
#line 349 "dlaed2.f"
/* L60: */
#line 349 "dlaed2.f"
    }


#line 352 "dlaed2.f"
    *k = 0;
#line 353 "dlaed2.f"
    k2 = *n + 1;
#line 354 "dlaed2.f"
    i__1 = *n;
#line 354 "dlaed2.f"
    for (j = 1; j <= i__1; ++j) {
#line 355 "dlaed2.f"
	nj = indx[j];
#line 356 "dlaed2.f"
	if (*rho * (d__1 = z__[nj], abs(d__1)) <= tol) {

/*           Deflate due to small z component. */

#line 360 "dlaed2.f"
	    --k2;
#line 361 "dlaed2.f"
	    coltyp[nj] = 4;
#line 362 "dlaed2.f"
	    indxp[k2] = nj;
#line 363 "dlaed2.f"
	    if (j == *n) {
#line 363 "dlaed2.f"
		goto L100;
#line 363 "dlaed2.f"
	    }
#line 365 "dlaed2.f"
	} else {
#line 366 "dlaed2.f"
	    pj = nj;
#line 367 "dlaed2.f"
	    goto L80;
#line 368 "dlaed2.f"
	}
#line 369 "dlaed2.f"
/* L70: */
#line 369 "dlaed2.f"
    }
#line 370 "dlaed2.f"
L80:
#line 371 "dlaed2.f"
    ++j;
#line 372 "dlaed2.f"
    nj = indx[j];
#line 373 "dlaed2.f"
    if (j > *n) {
#line 373 "dlaed2.f"
	goto L100;
#line 373 "dlaed2.f"
    }
#line 375 "dlaed2.f"
    if (*rho * (d__1 = z__[nj], abs(d__1)) <= tol) {

/*        Deflate due to small z component. */

#line 379 "dlaed2.f"
	--k2;
#line 380 "dlaed2.f"
	coltyp[nj] = 4;
#line 381 "dlaed2.f"
	indxp[k2] = nj;
#line 382 "dlaed2.f"
    } else {

/*        Check if eigenvalues are close enough to allow deflation. */

#line 386 "dlaed2.f"
	s = z__[pj];
#line 387 "dlaed2.f"
	c__ = z__[nj];

/*        Find sqrt(a**2+b**2) without overflow or */
/*        destructive underflow. */

#line 392 "dlaed2.f"
	tau = dlapy2_(&c__, &s);
#line 393 "dlaed2.f"
	t = d__[nj] - d__[pj];
#line 394 "dlaed2.f"
	c__ /= tau;
#line 395 "dlaed2.f"
	s = -s / tau;
#line 396 "dlaed2.f"
	if ((d__1 = t * c__ * s, abs(d__1)) <= tol) {

/*           Deflation is possible. */

#line 400 "dlaed2.f"
	    z__[nj] = tau;
#line 401 "dlaed2.f"
	    z__[pj] = 0.;
#line 402 "dlaed2.f"
	    if (coltyp[nj] != coltyp[pj]) {
#line 402 "dlaed2.f"
		coltyp[nj] = 2;
#line 402 "dlaed2.f"
	    }
#line 404 "dlaed2.f"
	    coltyp[pj] = 4;
#line 405 "dlaed2.f"
	    drot_(n, &q[pj * q_dim1 + 1], &c__1, &q[nj * q_dim1 + 1], &c__1, &
		    c__, &s);
/* Computing 2nd power */
#line 406 "dlaed2.f"
	    d__1 = c__;
/* Computing 2nd power */
#line 406 "dlaed2.f"
	    d__2 = s;
#line 406 "dlaed2.f"
	    t = d__[pj] * (d__1 * d__1) + d__[nj] * (d__2 * d__2);
/* Computing 2nd power */
#line 407 "dlaed2.f"
	    d__1 = s;
/* Computing 2nd power */
#line 407 "dlaed2.f"
	    d__2 = c__;
#line 407 "dlaed2.f"
	    d__[nj] = d__[pj] * (d__1 * d__1) + d__[nj] * (d__2 * d__2);
#line 408 "dlaed2.f"
	    d__[pj] = t;
#line 409 "dlaed2.f"
	    --k2;
#line 410 "dlaed2.f"
	    i__ = 1;
#line 411 "dlaed2.f"
L90:
#line 412 "dlaed2.f"
	    if (k2 + i__ <= *n) {
#line 413 "dlaed2.f"
		if (d__[pj] < d__[indxp[k2 + i__]]) {
#line 414 "dlaed2.f"
		    indxp[k2 + i__ - 1] = indxp[k2 + i__];
#line 415 "dlaed2.f"
		    indxp[k2 + i__] = pj;
#line 416 "dlaed2.f"
		    ++i__;
#line 417 "dlaed2.f"
		    goto L90;
#line 418 "dlaed2.f"
		} else {
#line 419 "dlaed2.f"
		    indxp[k2 + i__ - 1] = pj;
#line 420 "dlaed2.f"
		}
#line 421 "dlaed2.f"
	    } else {
#line 422 "dlaed2.f"
		indxp[k2 + i__ - 1] = pj;
#line 423 "dlaed2.f"
	    }
#line 424 "dlaed2.f"
	    pj = nj;
#line 425 "dlaed2.f"
	} else {
#line 426 "dlaed2.f"
	    ++(*k);
#line 427 "dlaed2.f"
	    dlamda[*k] = d__[pj];
#line 428 "dlaed2.f"
	    w[*k] = z__[pj];
#line 429 "dlaed2.f"
	    indxp[*k] = pj;
#line 430 "dlaed2.f"
	    pj = nj;
#line 431 "dlaed2.f"
	}
#line 432 "dlaed2.f"
    }
#line 433 "dlaed2.f"
    goto L80;
#line 434 "dlaed2.f"
L100:

/*     Record the last eigenvalue. */

#line 438 "dlaed2.f"
    ++(*k);
#line 439 "dlaed2.f"
    dlamda[*k] = d__[pj];
#line 440 "dlaed2.f"
    w[*k] = z__[pj];
#line 441 "dlaed2.f"
    indxp[*k] = pj;

/*     Count up the total number of the various types of columns, then */
/*     form a permutation which positions the four column types into */
/*     four uniform groups (although one or more of these groups may be */
/*     empty). */

#line 448 "dlaed2.f"
    for (j = 1; j <= 4; ++j) {
#line 449 "dlaed2.f"
	ctot[j - 1] = 0;
#line 450 "dlaed2.f"
/* L110: */
#line 450 "dlaed2.f"
    }
#line 451 "dlaed2.f"
    i__1 = *n;
#line 451 "dlaed2.f"
    for (j = 1; j <= i__1; ++j) {
#line 452 "dlaed2.f"
	ct = coltyp[j];
#line 453 "dlaed2.f"
	++ctot[ct - 1];
#line 454 "dlaed2.f"
/* L120: */
#line 454 "dlaed2.f"
    }

/*     PSM(*) = Position in SubMatrix (of types 1 through 4) */

#line 458 "dlaed2.f"
    psm[0] = 1;
#line 459 "dlaed2.f"
    psm[1] = ctot[0] + 1;
#line 460 "dlaed2.f"
    psm[2] = psm[1] + ctot[1];
#line 461 "dlaed2.f"
    psm[3] = psm[2] + ctot[2];
#line 462 "dlaed2.f"
    *k = *n - ctot[3];

/*     Fill out the INDXC array so that the permutation which it induces */
/*     will place all type-1 columns first, all type-2 columns next, */
/*     then all type-3's, and finally all type-4's. */

#line 468 "dlaed2.f"
    i__1 = *n;
#line 468 "dlaed2.f"
    for (j = 1; j <= i__1; ++j) {
#line 469 "dlaed2.f"
	js = indxp[j];
#line 470 "dlaed2.f"
	ct = coltyp[js];
#line 471 "dlaed2.f"
	indx[psm[ct - 1]] = js;
#line 472 "dlaed2.f"
	indxc[psm[ct - 1]] = j;
#line 473 "dlaed2.f"
	++psm[ct - 1];
#line 474 "dlaed2.f"
/* L130: */
#line 474 "dlaed2.f"
    }

/*     Sort the eigenvalues and corresponding eigenvectors into DLAMDA */
/*     and Q2 respectively.  The eigenvalues/vectors which were not */
/*     deflated go into the first K slots of DLAMDA and Q2 respectively, */
/*     while those which were deflated go into the last N - K slots. */

#line 481 "dlaed2.f"
    i__ = 1;
#line 482 "dlaed2.f"
    iq1 = 1;
#line 483 "dlaed2.f"
    iq2 = (ctot[0] + ctot[1]) * *n1 + 1;
#line 484 "dlaed2.f"
    i__1 = ctot[0];
#line 484 "dlaed2.f"
    for (j = 1; j <= i__1; ++j) {
#line 485 "dlaed2.f"
	js = indx[i__];
#line 486 "dlaed2.f"
	dcopy_(n1, &q[js * q_dim1 + 1], &c__1, &q2[iq1], &c__1);
#line 487 "dlaed2.f"
	z__[i__] = d__[js];
#line 488 "dlaed2.f"
	++i__;
#line 489 "dlaed2.f"
	iq1 += *n1;
#line 490 "dlaed2.f"
/* L140: */
#line 490 "dlaed2.f"
    }

#line 492 "dlaed2.f"
    i__1 = ctot[1];
#line 492 "dlaed2.f"
    for (j = 1; j <= i__1; ++j) {
#line 493 "dlaed2.f"
	js = indx[i__];
#line 494 "dlaed2.f"
	dcopy_(n1, &q[js * q_dim1 + 1], &c__1, &q2[iq1], &c__1);
#line 495 "dlaed2.f"
	dcopy_(&n2, &q[*n1 + 1 + js * q_dim1], &c__1, &q2[iq2], &c__1);
#line 496 "dlaed2.f"
	z__[i__] = d__[js];
#line 497 "dlaed2.f"
	++i__;
#line 498 "dlaed2.f"
	iq1 += *n1;
#line 499 "dlaed2.f"
	iq2 += n2;
#line 500 "dlaed2.f"
/* L150: */
#line 500 "dlaed2.f"
    }

#line 502 "dlaed2.f"
    i__1 = ctot[2];
#line 502 "dlaed2.f"
    for (j = 1; j <= i__1; ++j) {
#line 503 "dlaed2.f"
	js = indx[i__];
#line 504 "dlaed2.f"
	dcopy_(&n2, &q[*n1 + 1 + js * q_dim1], &c__1, &q2[iq2], &c__1);
#line 505 "dlaed2.f"
	z__[i__] = d__[js];
#line 506 "dlaed2.f"
	++i__;
#line 507 "dlaed2.f"
	iq2 += n2;
#line 508 "dlaed2.f"
/* L160: */
#line 508 "dlaed2.f"
    }

#line 510 "dlaed2.f"
    iq1 = iq2;
#line 511 "dlaed2.f"
    i__1 = ctot[3];
#line 511 "dlaed2.f"
    for (j = 1; j <= i__1; ++j) {
#line 512 "dlaed2.f"
	js = indx[i__];
#line 513 "dlaed2.f"
	dcopy_(n, &q[js * q_dim1 + 1], &c__1, &q2[iq2], &c__1);
#line 514 "dlaed2.f"
	iq2 += *n;
#line 515 "dlaed2.f"
	z__[i__] = d__[js];
#line 516 "dlaed2.f"
	++i__;
#line 517 "dlaed2.f"
/* L170: */
#line 517 "dlaed2.f"
    }

/*     The deflated eigenvalues and their corresponding vectors go back */
/*     into the last N - K slots of D and Q respectively. */

#line 522 "dlaed2.f"
    if (*k < *n) {
#line 523 "dlaed2.f"
	dlacpy_("A", n, &ctot[3], &q2[iq1], n, &q[(*k + 1) * q_dim1 + 1], ldq,
		 (ftnlen)1);
#line 525 "dlaed2.f"
	i__1 = *n - *k;
#line 525 "dlaed2.f"
	dcopy_(&i__1, &z__[*k + 1], &c__1, &d__[*k + 1], &c__1);
#line 526 "dlaed2.f"
    }

/*     Copy CTOT into COLTYP for referencing in DLAED3. */

#line 530 "dlaed2.f"
    for (j = 1; j <= 4; ++j) {
#line 531 "dlaed2.f"
	coltyp[j] = ctot[j - 1];
#line 532 "dlaed2.f"
/* L180: */
#line 532 "dlaed2.f"
    }

#line 534 "dlaed2.f"
L190:
#line 535 "dlaed2.f"
    return 0;

/*     End of DLAED2 */

} /* dlaed2_ */

