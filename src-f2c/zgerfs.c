#line 1 "zgerfs.f"
/* zgerfs.f -- translated by f2c (version 20100827).
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

#line 1 "zgerfs.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static integer c__1 = 1;

/* > \brief \b ZGERFS */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZGERFS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgerfs.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgerfs.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgerfs.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZGERFS( TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV, B, LDB, */
/*                          X, LDX, FERR, BERR, WORK, RWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          TRANS */
/*       INTEGER            INFO, LDA, LDAF, LDB, LDX, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       DOUBLE PRECISION   BERR( * ), FERR( * ), RWORK( * ) */
/*       COMPLEX*16         A( LDA, * ), AF( LDAF, * ), B( LDB, * ), */
/*      $                   WORK( * ), X( LDX, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZGERFS improves the computed solution to a system of linear */
/* > equations and provides error bounds and backward error estimates for */
/* > the solution. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] TRANS */
/* > \verbatim */
/* >          TRANS is CHARACTER*1 */
/* >          Specifies the form of the system of equations: */
/* >          = 'N':  A * X = B     (No transpose) */
/* >          = 'T':  A**T * X = B  (Transpose) */
/* >          = 'C':  A**H * X = B  (Conjugate transpose) */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NRHS */
/* > \verbatim */
/* >          NRHS is INTEGER */
/* >          The number of right hand sides, i.e., the number of columns */
/* >          of the matrices B and X.  NRHS >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
/* >          The original N-by-N matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] AF */
/* > \verbatim */
/* >          AF is COMPLEX*16 array, dimension (LDAF,N) */
/* >          The factors L and U from the factorization A = P*L*U */
/* >          as computed by ZGETRF. */
/* > \endverbatim */
/* > */
/* > \param[in] LDAF */
/* > \verbatim */
/* >          LDAF is INTEGER */
/* >          The leading dimension of the array AF.  LDAF >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N) */
/* >          The pivot indices from ZGETRF; for 1<=i<=N, row i of the */
/* >          matrix was interchanged with row IPIV(i). */
/* > \endverbatim */
/* > */
/* > \param[in] B */
/* > \verbatim */
/* >          B is COMPLEX*16 array, dimension (LDB,NRHS) */
/* >          The right hand side matrix B. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >          The leading dimension of the array B.  LDB >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] X */
/* > \verbatim */
/* >          X is COMPLEX*16 array, dimension (LDX,NRHS) */
/* >          On entry, the solution matrix X, as computed by ZGETRS. */
/* >          On exit, the improved solution matrix X. */
/* > \endverbatim */
/* > */
/* > \param[in] LDX */
/* > \verbatim */
/* >          LDX is INTEGER */
/* >          The leading dimension of the array X.  LDX >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] FERR */
/* > \verbatim */
/* >          FERR is DOUBLE PRECISION array, dimension (NRHS) */
/* >          The estimated forward error bound for each solution vector */
/* >          X(j) (the j-th column of the solution matrix X). */
/* >          If XTRUE is the true solution corresponding to X(j), FERR(j) */
/* >          is an estimated upper bound for the magnitude of the largest */
/* >          element in (X(j) - XTRUE) divided by the magnitude of the */
/* >          largest element in X(j).  The estimate is as reliable as */
/* >          the estimate for RCOND, and is almost always a slight */
/* >          overestimate of the true error. */
/* > \endverbatim */
/* > */
/* > \param[out] BERR */
/* > \verbatim */
/* >          BERR is DOUBLE PRECISION array, dimension (NRHS) */
/* >          The componentwise relative backward error of each solution */
/* >          vector X(j) (i.e., the smallest relative change in */
/* >          any element of A or B that makes X(j) an exact solution). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX*16 array, dimension (2*N) */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is DOUBLE PRECISION array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */

/* > \par Internal Parameters: */
/*  ========================= */
/* > */
/* > \verbatim */
/* >  ITMAX is the maximum number of steps of iterative refinement. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup complex16GEcomputational */

/*  ===================================================================== */
/* Subroutine */ int zgerfs_(char *trans, integer *n, integer *nrhs, 
	doublecomplex *a, integer *lda, doublecomplex *af, integer *ldaf, 
	integer *ipiv, doublecomplex *b, integer *ldb, doublecomplex *x, 
	integer *ldx, doublereal *ferr, doublereal *berr, doublecomplex *work,
	 doublereal *rwork, integer *info, ftnlen trans_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, af_dim1, af_offset, b_dim1, b_offset, x_dim1, 
	    x_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2, d__3, d__4;
    doublecomplex z__1;

    /* Builtin functions */
    double d_imag(doublecomplex *);

    /* Local variables */
    static integer i__, j, k;
    static doublereal s, xk;
    static integer nz;
    static doublereal eps;
    static integer kase;
    static doublereal safe1, safe2;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer isave[3], count;
    extern /* Subroutine */ int zgemv_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, ftnlen), 
	    zcopy_(integer *, doublecomplex *, integer *, doublecomplex *, 
	    integer *), zaxpy_(integer *, doublecomplex *, doublecomplex *, 
	    integer *, doublecomplex *, integer *), zlacn2_(integer *, 
	    doublecomplex *, doublecomplex *, doublereal *, integer *, 
	    integer *);
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static logical notran;
    static char transn[1], transt[1];
    static doublereal lstres;
    extern /* Subroutine */ int zgetrs_(char *, integer *, integer *, 
	    doublecomplex *, integer *, integer *, doublecomplex *, integer *,
	     integer *, ftnlen);


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
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Statement Functions .. */
/*     .. */
/*     .. Statement Function definitions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 250 "zgerfs.f"
    /* Parameter adjustments */
#line 250 "zgerfs.f"
    a_dim1 = *lda;
#line 250 "zgerfs.f"
    a_offset = 1 + a_dim1;
#line 250 "zgerfs.f"
    a -= a_offset;
#line 250 "zgerfs.f"
    af_dim1 = *ldaf;
#line 250 "zgerfs.f"
    af_offset = 1 + af_dim1;
#line 250 "zgerfs.f"
    af -= af_offset;
#line 250 "zgerfs.f"
    --ipiv;
#line 250 "zgerfs.f"
    b_dim1 = *ldb;
#line 250 "zgerfs.f"
    b_offset = 1 + b_dim1;
#line 250 "zgerfs.f"
    b -= b_offset;
#line 250 "zgerfs.f"
    x_dim1 = *ldx;
#line 250 "zgerfs.f"
    x_offset = 1 + x_dim1;
#line 250 "zgerfs.f"
    x -= x_offset;
#line 250 "zgerfs.f"
    --ferr;
#line 250 "zgerfs.f"
    --berr;
#line 250 "zgerfs.f"
    --work;
#line 250 "zgerfs.f"
    --rwork;
#line 250 "zgerfs.f"

#line 250 "zgerfs.f"
    /* Function Body */
#line 250 "zgerfs.f"
    *info = 0;
#line 251 "zgerfs.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 252 "zgerfs.f"
    if (! notran && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1) && ! lsame_(
	    trans, "C", (ftnlen)1, (ftnlen)1)) {
#line 254 "zgerfs.f"
	*info = -1;
#line 255 "zgerfs.f"
    } else if (*n < 0) {
#line 256 "zgerfs.f"
	*info = -2;
#line 257 "zgerfs.f"
    } else if (*nrhs < 0) {
#line 258 "zgerfs.f"
	*info = -3;
#line 259 "zgerfs.f"
    } else if (*lda < max(1,*n)) {
#line 260 "zgerfs.f"
	*info = -5;
#line 261 "zgerfs.f"
    } else if (*ldaf < max(1,*n)) {
#line 262 "zgerfs.f"
	*info = -7;
#line 263 "zgerfs.f"
    } else if (*ldb < max(1,*n)) {
#line 264 "zgerfs.f"
	*info = -10;
#line 265 "zgerfs.f"
    } else if (*ldx < max(1,*n)) {
#line 266 "zgerfs.f"
	*info = -12;
#line 267 "zgerfs.f"
    }
#line 268 "zgerfs.f"
    if (*info != 0) {
#line 269 "zgerfs.f"
	i__1 = -(*info);
#line 269 "zgerfs.f"
	xerbla_("ZGERFS", &i__1, (ftnlen)6);
#line 270 "zgerfs.f"
	return 0;
#line 271 "zgerfs.f"
    }

/*     Quick return if possible */

#line 275 "zgerfs.f"
    if (*n == 0 || *nrhs == 0) {
#line 276 "zgerfs.f"
	i__1 = *nrhs;
#line 276 "zgerfs.f"
	for (j = 1; j <= i__1; ++j) {
#line 277 "zgerfs.f"
	    ferr[j] = 0.;
#line 278 "zgerfs.f"
	    berr[j] = 0.;
#line 279 "zgerfs.f"
/* L10: */
#line 279 "zgerfs.f"
	}
#line 280 "zgerfs.f"
	return 0;
#line 281 "zgerfs.f"
    }

#line 283 "zgerfs.f"
    if (notran) {
#line 284 "zgerfs.f"
	*(unsigned char *)transn = 'N';
#line 285 "zgerfs.f"
	*(unsigned char *)transt = 'C';
#line 286 "zgerfs.f"
    } else {
#line 287 "zgerfs.f"
	*(unsigned char *)transn = 'C';
#line 288 "zgerfs.f"
	*(unsigned char *)transt = 'N';
#line 289 "zgerfs.f"
    }

/*     NZ = maximum number of nonzero elements in each row of A, plus 1 */

#line 293 "zgerfs.f"
    nz = *n + 1;
#line 294 "zgerfs.f"
    eps = dlamch_("Epsilon", (ftnlen)7);
#line 295 "zgerfs.f"
    safmin = dlamch_("Safe minimum", (ftnlen)12);
#line 296 "zgerfs.f"
    safe1 = nz * safmin;
#line 297 "zgerfs.f"
    safe2 = safe1 / eps;

/*     Do for each right hand side */

#line 301 "zgerfs.f"
    i__1 = *nrhs;
#line 301 "zgerfs.f"
    for (j = 1; j <= i__1; ++j) {

#line 303 "zgerfs.f"
	count = 1;
#line 304 "zgerfs.f"
	lstres = 3.;
#line 305 "zgerfs.f"
L20:

/*        Loop until stopping criterion is satisfied. */

/*        Compute residual R = B - op(A) * X, */
/*        where op(A) = A, A**T, or A**H, depending on TRANS. */

#line 312 "zgerfs.f"
	zcopy_(n, &b[j * b_dim1 + 1], &c__1, &work[1], &c__1);
#line 313 "zgerfs.f"
	z__1.r = -1., z__1.i = -0.;
#line 313 "zgerfs.f"
	zgemv_(trans, n, n, &z__1, &a[a_offset], lda, &x[j * x_dim1 + 1], &
		c__1, &c_b1, &work[1], &c__1, (ftnlen)1);

/*        Compute componentwise relative backward error from formula */

/*        max(i) ( abs(R(i)) / ( abs(op(A))*abs(X) + abs(B) )(i) ) */

/*        where abs(Z) is the componentwise absolute value of the matrix */
/*        or vector Z.  If the i-th component of the denominator is less */
/*        than SAFE2, then SAFE1 is added to the i-th components of the */
/*        numerator and denominator before dividing. */

#line 325 "zgerfs.f"
	i__2 = *n;
#line 325 "zgerfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 326 "zgerfs.f"
	    i__3 = i__ + j * b_dim1;
#line 326 "zgerfs.f"
	    rwork[i__] = (d__1 = b[i__3].r, abs(d__1)) + (d__2 = d_imag(&b[
		    i__ + j * b_dim1]), abs(d__2));
#line 327 "zgerfs.f"
/* L30: */
#line 327 "zgerfs.f"
	}

/*        Compute abs(op(A))*abs(X) + abs(B). */

#line 331 "zgerfs.f"
	if (notran) {
#line 332 "zgerfs.f"
	    i__2 = *n;
#line 332 "zgerfs.f"
	    for (k = 1; k <= i__2; ++k) {
#line 333 "zgerfs.f"
		i__3 = k + j * x_dim1;
#line 333 "zgerfs.f"
		xk = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = d_imag(&x[k + j *
			 x_dim1]), abs(d__2));
#line 334 "zgerfs.f"
		i__3 = *n;
#line 334 "zgerfs.f"
		for (i__ = 1; i__ <= i__3; ++i__) {
#line 335 "zgerfs.f"
		    i__4 = i__ + k * a_dim1;
#line 335 "zgerfs.f"
		    rwork[i__] += ((d__1 = a[i__4].r, abs(d__1)) + (d__2 = 
			    d_imag(&a[i__ + k * a_dim1]), abs(d__2))) * xk;
#line 336 "zgerfs.f"
/* L40: */
#line 336 "zgerfs.f"
		}
#line 337 "zgerfs.f"
/* L50: */
#line 337 "zgerfs.f"
	    }
#line 338 "zgerfs.f"
	} else {
#line 339 "zgerfs.f"
	    i__2 = *n;
#line 339 "zgerfs.f"
	    for (k = 1; k <= i__2; ++k) {
#line 340 "zgerfs.f"
		s = 0.;
#line 341 "zgerfs.f"
		i__3 = *n;
#line 341 "zgerfs.f"
		for (i__ = 1; i__ <= i__3; ++i__) {
#line 342 "zgerfs.f"
		    i__4 = i__ + k * a_dim1;
#line 342 "zgerfs.f"
		    i__5 = i__ + j * x_dim1;
#line 342 "zgerfs.f"
		    s += ((d__1 = a[i__4].r, abs(d__1)) + (d__2 = d_imag(&a[
			    i__ + k * a_dim1]), abs(d__2))) * ((d__3 = x[i__5]
			    .r, abs(d__3)) + (d__4 = d_imag(&x[i__ + j * 
			    x_dim1]), abs(d__4)));
#line 343 "zgerfs.f"
/* L60: */
#line 343 "zgerfs.f"
		}
#line 344 "zgerfs.f"
		rwork[k] += s;
#line 345 "zgerfs.f"
/* L70: */
#line 345 "zgerfs.f"
	    }
#line 346 "zgerfs.f"
	}
#line 347 "zgerfs.f"
	s = 0.;
#line 348 "zgerfs.f"
	i__2 = *n;
#line 348 "zgerfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 349 "zgerfs.f"
	    if (rwork[i__] > safe2) {
/* Computing MAX */
#line 350 "zgerfs.f"
		i__3 = i__;
#line 350 "zgerfs.f"
		d__3 = s, d__4 = ((d__1 = work[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&work[i__]), abs(d__2))) / rwork[i__];
#line 350 "zgerfs.f"
		s = max(d__3,d__4);
#line 351 "zgerfs.f"
	    } else {
/* Computing MAX */
#line 352 "zgerfs.f"
		i__3 = i__;
#line 352 "zgerfs.f"
		d__3 = s, d__4 = ((d__1 = work[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&work[i__]), abs(d__2)) + safe1) / (rwork[i__] 
			+ safe1);
#line 352 "zgerfs.f"
		s = max(d__3,d__4);
#line 354 "zgerfs.f"
	    }
#line 355 "zgerfs.f"
/* L80: */
#line 355 "zgerfs.f"
	}
#line 356 "zgerfs.f"
	berr[j] = s;

/*        Test stopping criterion. Continue iterating if */
/*           1) The residual BERR(J) is larger than machine epsilon, and */
/*           2) BERR(J) decreased by at least a factor of 2 during the */
/*              last iteration, and */
/*           3) At most ITMAX iterations tried. */

#line 364 "zgerfs.f"
	if (berr[j] > eps && berr[j] * 2. <= lstres && count <= 5) {

/*           Update solution and try again. */

#line 369 "zgerfs.f"
	    zgetrs_(trans, n, &c__1, &af[af_offset], ldaf, &ipiv[1], &work[1],
		     n, info, (ftnlen)1);
#line 370 "zgerfs.f"
	    zaxpy_(n, &c_b1, &work[1], &c__1, &x[j * x_dim1 + 1], &c__1);
#line 371 "zgerfs.f"
	    lstres = berr[j];
#line 372 "zgerfs.f"
	    ++count;
#line 373 "zgerfs.f"
	    goto L20;
#line 374 "zgerfs.f"
	}

/*        Bound error from formula */

/*        norm(X - XTRUE) / norm(X) .le. FERR = */
/*        norm( abs(inv(op(A)))* */
/*           ( abs(R) + NZ*EPS*( abs(op(A))*abs(X)+abs(B) ))) / norm(X) */

/*        where */
/*          norm(Z) is the magnitude of the largest component of Z */
/*          inv(op(A)) is the inverse of op(A) */
/*          abs(Z) is the componentwise absolute value of the matrix or */
/*             vector Z */
/*          NZ is the maximum number of nonzeros in any row of A, plus 1 */
/*          EPS is machine epsilon */

/*        The i-th component of abs(R)+NZ*EPS*(abs(op(A))*abs(X)+abs(B)) */
/*        is incremented by SAFE1 if the i-th component of */
/*        abs(op(A))*abs(X) + abs(B) is less than SAFE2. */

/*        Use ZLACN2 to estimate the infinity-norm of the matrix */
/*           inv(op(A)) * diag(W), */
/*        where W = abs(R) + NZ*EPS*( abs(op(A))*abs(X)+abs(B) ))) */

#line 398 "zgerfs.f"
	i__2 = *n;
#line 398 "zgerfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 399 "zgerfs.f"
	    if (rwork[i__] > safe2) {
#line 400 "zgerfs.f"
		i__3 = i__;
#line 400 "zgerfs.f"
		rwork[i__] = (d__1 = work[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&work[i__]), abs(d__2)) + nz * eps * rwork[i__]
			;
#line 401 "zgerfs.f"
	    } else {
#line 402 "zgerfs.f"
		i__3 = i__;
#line 402 "zgerfs.f"
		rwork[i__] = (d__1 = work[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&work[i__]), abs(d__2)) + nz * eps * rwork[i__]
			 + safe1;
#line 404 "zgerfs.f"
	    }
#line 405 "zgerfs.f"
/* L90: */
#line 405 "zgerfs.f"
	}

#line 407 "zgerfs.f"
	kase = 0;
#line 408 "zgerfs.f"
L100:
#line 409 "zgerfs.f"
	zlacn2_(n, &work[*n + 1], &work[1], &ferr[j], &kase, isave);
#line 410 "zgerfs.f"
	if (kase != 0) {
#line 411 "zgerfs.f"
	    if (kase == 1) {

/*              Multiply by diag(W)*inv(op(A)**H). */

#line 415 "zgerfs.f"
		zgetrs_(transt, n, &c__1, &af[af_offset], ldaf, &ipiv[1], &
			work[1], n, info, (ftnlen)1);
#line 417 "zgerfs.f"
		i__2 = *n;
#line 417 "zgerfs.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 418 "zgerfs.f"
		    i__3 = i__;
#line 418 "zgerfs.f"
		    i__4 = i__;
#line 418 "zgerfs.f"
		    i__5 = i__;
#line 418 "zgerfs.f"
		    z__1.r = rwork[i__4] * work[i__5].r, z__1.i = rwork[i__4] 
			    * work[i__5].i;
#line 418 "zgerfs.f"
		    work[i__3].r = z__1.r, work[i__3].i = z__1.i;
#line 419 "zgerfs.f"
/* L110: */
#line 419 "zgerfs.f"
		}
#line 420 "zgerfs.f"
	    } else {

/*              Multiply by inv(op(A))*diag(W). */

#line 424 "zgerfs.f"
		i__2 = *n;
#line 424 "zgerfs.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 425 "zgerfs.f"
		    i__3 = i__;
#line 425 "zgerfs.f"
		    i__4 = i__;
#line 425 "zgerfs.f"
		    i__5 = i__;
#line 425 "zgerfs.f"
		    z__1.r = rwork[i__4] * work[i__5].r, z__1.i = rwork[i__4] 
			    * work[i__5].i;
#line 425 "zgerfs.f"
		    work[i__3].r = z__1.r, work[i__3].i = z__1.i;
#line 426 "zgerfs.f"
/* L120: */
#line 426 "zgerfs.f"
		}
#line 427 "zgerfs.f"
		zgetrs_(transn, n, &c__1, &af[af_offset], ldaf, &ipiv[1], &
			work[1], n, info, (ftnlen)1);
#line 429 "zgerfs.f"
	    }
#line 430 "zgerfs.f"
	    goto L100;
#line 431 "zgerfs.f"
	}

/*        Normalize error. */

#line 435 "zgerfs.f"
	lstres = 0.;
#line 436 "zgerfs.f"
	i__2 = *n;
#line 436 "zgerfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 437 "zgerfs.f"
	    i__3 = i__ + j * x_dim1;
#line 437 "zgerfs.f"
	    d__3 = lstres, d__4 = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = 
		    d_imag(&x[i__ + j * x_dim1]), abs(d__2));
#line 437 "zgerfs.f"
	    lstres = max(d__3,d__4);
#line 438 "zgerfs.f"
/* L130: */
#line 438 "zgerfs.f"
	}
#line 439 "zgerfs.f"
	if (lstres != 0.) {
#line 439 "zgerfs.f"
	    ferr[j] /= lstres;
#line 439 "zgerfs.f"
	}

#line 442 "zgerfs.f"
/* L140: */
#line 442 "zgerfs.f"
    }

#line 444 "zgerfs.f"
    return 0;

/*     End of ZGERFS */

} /* zgerfs_ */

