#line 1 "dsyrfs.f"
/* dsyrfs.f -- translated by f2c (version 20100827).
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

#line 1 "dsyrfs.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b12 = -1.;
static doublereal c_b14 = 1.;

/* > \brief \b DSYRFS */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DSYRFS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsyrfs.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsyrfs.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsyrfs.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSYRFS( UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV, B, LDB, */
/*                          X, LDX, FERR, BERR, WORK, IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDA, LDAF, LDB, LDX, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ), IWORK( * ) */
/*       DOUBLE PRECISION   A( LDA, * ), AF( LDAF, * ), B( LDB, * ), */
/*      $                   BERR( * ), FERR( * ), WORK( * ), X( LDX, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DSYRFS improves the computed solution to a system of linear */
/* > equations when the coefficient matrix is symmetric indefinite, and */
/* > provides error bounds and backward error estimates for the solution. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          = 'U':  Upper triangle of A is stored; */
/* >          = 'L':  Lower triangle of A is stored. */
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
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* >          The symmetric matrix A.  If UPLO = 'U', the leading N-by-N */
/* >          upper triangular part of A contains the upper triangular part */
/* >          of the matrix A, and the strictly lower triangular part of A */
/* >          is not referenced.  If UPLO = 'L', the leading N-by-N lower */
/* >          triangular part of A contains the lower triangular part of */
/* >          the matrix A, and the strictly upper triangular part of A is */
/* >          not referenced. */
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
/* >          AF is DOUBLE PRECISION array, dimension (LDAF,N) */
/* >          The factored form of the matrix A.  AF contains the block */
/* >          diagonal matrix D and the multipliers used to obtain the */
/* >          factor U or L from the factorization A = U*D*U**T or */
/* >          A = L*D*L**T as computed by DSYTRF. */
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
/* >          Details of the interchanges and the block structure of D */
/* >          as determined by DSYTRF. */
/* > \endverbatim */
/* > */
/* > \param[in] B */
/* > \verbatim */
/* >          B is DOUBLE PRECISION array, dimension (LDB,NRHS) */
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
/* >          X is DOUBLE PRECISION array, dimension (LDX,NRHS) */
/* >          On entry, the solution matrix X, as computed by DSYTRS. */
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
/* >          WORK is DOUBLE PRECISION array, dimension (3*N) */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* >          IWORK is INTEGER array, dimension (N) */
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

/* > \date November 2011 */

/* > \ingroup doubleSYcomputational */

/*  ===================================================================== */
/* Subroutine */ int dsyrfs_(char *uplo, integer *n, integer *nrhs, 
	doublereal *a, integer *lda, doublereal *af, integer *ldaf, integer *
	ipiv, doublereal *b, integer *ldb, doublereal *x, integer *ldx, 
	doublereal *ferr, doublereal *berr, doublereal *work, integer *iwork, 
	integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, af_dim1, af_offset, b_dim1, b_offset, x_dim1, 
	    x_offset, i__1, i__2, i__3;
    doublereal d__1, d__2, d__3;

    /* Local variables */
    static integer i__, j, k;
    static doublereal s, xk;
    static integer nz;
    static doublereal eps;
    static integer kase;
    static doublereal safe1, safe2;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer isave[3];
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), daxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *);
    static integer count;
    static logical upper;
    extern /* Subroutine */ int dsymv_(char *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, ftnlen), dlacn2_(integer *, doublereal *,
	     doublereal *, integer *, doublereal *, integer *, integer *);
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal lstres;
    extern /* Subroutine */ int dsytrs_(char *, integer *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, ftnlen);


/*  -- LAPACK computational routine (version 3.4.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */

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
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 246 "dsyrfs.f"
    /* Parameter adjustments */
#line 246 "dsyrfs.f"
    a_dim1 = *lda;
#line 246 "dsyrfs.f"
    a_offset = 1 + a_dim1;
#line 246 "dsyrfs.f"
    a -= a_offset;
#line 246 "dsyrfs.f"
    af_dim1 = *ldaf;
#line 246 "dsyrfs.f"
    af_offset = 1 + af_dim1;
#line 246 "dsyrfs.f"
    af -= af_offset;
#line 246 "dsyrfs.f"
    --ipiv;
#line 246 "dsyrfs.f"
    b_dim1 = *ldb;
#line 246 "dsyrfs.f"
    b_offset = 1 + b_dim1;
#line 246 "dsyrfs.f"
    b -= b_offset;
#line 246 "dsyrfs.f"
    x_dim1 = *ldx;
#line 246 "dsyrfs.f"
    x_offset = 1 + x_dim1;
#line 246 "dsyrfs.f"
    x -= x_offset;
#line 246 "dsyrfs.f"
    --ferr;
#line 246 "dsyrfs.f"
    --berr;
#line 246 "dsyrfs.f"
    --work;
#line 246 "dsyrfs.f"
    --iwork;
#line 246 "dsyrfs.f"

#line 246 "dsyrfs.f"
    /* Function Body */
#line 246 "dsyrfs.f"
    *info = 0;
#line 247 "dsyrfs.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 248 "dsyrfs.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 249 "dsyrfs.f"
	*info = -1;
#line 250 "dsyrfs.f"
    } else if (*n < 0) {
#line 251 "dsyrfs.f"
	*info = -2;
#line 252 "dsyrfs.f"
    } else if (*nrhs < 0) {
#line 253 "dsyrfs.f"
	*info = -3;
#line 254 "dsyrfs.f"
    } else if (*lda < max(1,*n)) {
#line 255 "dsyrfs.f"
	*info = -5;
#line 256 "dsyrfs.f"
    } else if (*ldaf < max(1,*n)) {
#line 257 "dsyrfs.f"
	*info = -7;
#line 258 "dsyrfs.f"
    } else if (*ldb < max(1,*n)) {
#line 259 "dsyrfs.f"
	*info = -10;
#line 260 "dsyrfs.f"
    } else if (*ldx < max(1,*n)) {
#line 261 "dsyrfs.f"
	*info = -12;
#line 262 "dsyrfs.f"
    }
#line 263 "dsyrfs.f"
    if (*info != 0) {
#line 264 "dsyrfs.f"
	i__1 = -(*info);
#line 264 "dsyrfs.f"
	xerbla_("DSYRFS", &i__1, (ftnlen)6);
#line 265 "dsyrfs.f"
	return 0;
#line 266 "dsyrfs.f"
    }

/*     Quick return if possible */

#line 270 "dsyrfs.f"
    if (*n == 0 || *nrhs == 0) {
#line 271 "dsyrfs.f"
	i__1 = *nrhs;
#line 271 "dsyrfs.f"
	for (j = 1; j <= i__1; ++j) {
#line 272 "dsyrfs.f"
	    ferr[j] = 0.;
#line 273 "dsyrfs.f"
	    berr[j] = 0.;
#line 274 "dsyrfs.f"
/* L10: */
#line 274 "dsyrfs.f"
	}
#line 275 "dsyrfs.f"
	return 0;
#line 276 "dsyrfs.f"
    }

/*     NZ = maximum number of nonzero elements in each row of A, plus 1 */

#line 280 "dsyrfs.f"
    nz = *n + 1;
#line 281 "dsyrfs.f"
    eps = dlamch_("Epsilon", (ftnlen)7);
#line 282 "dsyrfs.f"
    safmin = dlamch_("Safe minimum", (ftnlen)12);
#line 283 "dsyrfs.f"
    safe1 = nz * safmin;
#line 284 "dsyrfs.f"
    safe2 = safe1 / eps;

/*     Do for each right hand side */

#line 288 "dsyrfs.f"
    i__1 = *nrhs;
#line 288 "dsyrfs.f"
    for (j = 1; j <= i__1; ++j) {

#line 290 "dsyrfs.f"
	count = 1;
#line 291 "dsyrfs.f"
	lstres = 3.;
#line 292 "dsyrfs.f"
L20:

/*        Loop until stopping criterion is satisfied. */

/*        Compute residual R = B - A * X */

#line 298 "dsyrfs.f"
	dcopy_(n, &b[j * b_dim1 + 1], &c__1, &work[*n + 1], &c__1);
#line 299 "dsyrfs.f"
	dsymv_(uplo, n, &c_b12, &a[a_offset], lda, &x[j * x_dim1 + 1], &c__1, 
		&c_b14, &work[*n + 1], &c__1, (ftnlen)1);

/*        Compute componentwise relative backward error from formula */

/*        max(i) ( abs(R(i)) / ( abs(A)*abs(X) + abs(B) )(i) ) */

/*        where abs(Z) is the componentwise absolute value of the matrix */
/*        or vector Z.  If the i-th component of the denominator is less */
/*        than SAFE2, then SAFE1 is added to the i-th components of the */
/*        numerator and denominator before dividing. */

#line 311 "dsyrfs.f"
	i__2 = *n;
#line 311 "dsyrfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 312 "dsyrfs.f"
	    work[i__] = (d__1 = b[i__ + j * b_dim1], abs(d__1));
#line 313 "dsyrfs.f"
/* L30: */
#line 313 "dsyrfs.f"
	}

/*        Compute abs(A)*abs(X) + abs(B). */

#line 317 "dsyrfs.f"
	if (upper) {
#line 318 "dsyrfs.f"
	    i__2 = *n;
#line 318 "dsyrfs.f"
	    for (k = 1; k <= i__2; ++k) {
#line 319 "dsyrfs.f"
		s = 0.;
#line 320 "dsyrfs.f"
		xk = (d__1 = x[k + j * x_dim1], abs(d__1));
#line 321 "dsyrfs.f"
		i__3 = k - 1;
#line 321 "dsyrfs.f"
		for (i__ = 1; i__ <= i__3; ++i__) {
#line 322 "dsyrfs.f"
		    work[i__] += (d__1 = a[i__ + k * a_dim1], abs(d__1)) * xk;
#line 323 "dsyrfs.f"
		    s += (d__1 = a[i__ + k * a_dim1], abs(d__1)) * (d__2 = x[
			    i__ + j * x_dim1], abs(d__2));
#line 324 "dsyrfs.f"
/* L40: */
#line 324 "dsyrfs.f"
		}
#line 325 "dsyrfs.f"
		work[k] = work[k] + (d__1 = a[k + k * a_dim1], abs(d__1)) * 
			xk + s;
#line 326 "dsyrfs.f"
/* L50: */
#line 326 "dsyrfs.f"
	    }
#line 327 "dsyrfs.f"
	} else {
#line 328 "dsyrfs.f"
	    i__2 = *n;
#line 328 "dsyrfs.f"
	    for (k = 1; k <= i__2; ++k) {
#line 329 "dsyrfs.f"
		s = 0.;
#line 330 "dsyrfs.f"
		xk = (d__1 = x[k + j * x_dim1], abs(d__1));
#line 331 "dsyrfs.f"
		work[k] += (d__1 = a[k + k * a_dim1], abs(d__1)) * xk;
#line 332 "dsyrfs.f"
		i__3 = *n;
#line 332 "dsyrfs.f"
		for (i__ = k + 1; i__ <= i__3; ++i__) {
#line 333 "dsyrfs.f"
		    work[i__] += (d__1 = a[i__ + k * a_dim1], abs(d__1)) * xk;
#line 334 "dsyrfs.f"
		    s += (d__1 = a[i__ + k * a_dim1], abs(d__1)) * (d__2 = x[
			    i__ + j * x_dim1], abs(d__2));
#line 335 "dsyrfs.f"
/* L60: */
#line 335 "dsyrfs.f"
		}
#line 336 "dsyrfs.f"
		work[k] += s;
#line 337 "dsyrfs.f"
/* L70: */
#line 337 "dsyrfs.f"
	    }
#line 338 "dsyrfs.f"
	}
#line 339 "dsyrfs.f"
	s = 0.;
#line 340 "dsyrfs.f"
	i__2 = *n;
#line 340 "dsyrfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 341 "dsyrfs.f"
	    if (work[i__] > safe2) {
/* Computing MAX */
#line 342 "dsyrfs.f"
		d__2 = s, d__3 = (d__1 = work[*n + i__], abs(d__1)) / work[
			i__];
#line 342 "dsyrfs.f"
		s = max(d__2,d__3);
#line 343 "dsyrfs.f"
	    } else {
/* Computing MAX */
#line 344 "dsyrfs.f"
		d__2 = s, d__3 = ((d__1 = work[*n + i__], abs(d__1)) + safe1) 
			/ (work[i__] + safe1);
#line 344 "dsyrfs.f"
		s = max(d__2,d__3);
#line 346 "dsyrfs.f"
	    }
#line 347 "dsyrfs.f"
/* L80: */
#line 347 "dsyrfs.f"
	}
#line 348 "dsyrfs.f"
	berr[j] = s;

/*        Test stopping criterion. Continue iterating if */
/*           1) The residual BERR(J) is larger than machine epsilon, and */
/*           2) BERR(J) decreased by at least a factor of 2 during the */
/*              last iteration, and */
/*           3) At most ITMAX iterations tried. */

#line 356 "dsyrfs.f"
	if (berr[j] > eps && berr[j] * 2. <= lstres && count <= 5) {

/*           Update solution and try again. */

#line 361 "dsyrfs.f"
	    dsytrs_(uplo, n, &c__1, &af[af_offset], ldaf, &ipiv[1], &work[*n 
		    + 1], n, info, (ftnlen)1);
#line 363 "dsyrfs.f"
	    daxpy_(n, &c_b14, &work[*n + 1], &c__1, &x[j * x_dim1 + 1], &c__1)
		    ;
#line 364 "dsyrfs.f"
	    lstres = berr[j];
#line 365 "dsyrfs.f"
	    ++count;
#line 366 "dsyrfs.f"
	    goto L20;
#line 367 "dsyrfs.f"
	}

/*        Bound error from formula */

/*        norm(X - XTRUE) / norm(X) .le. FERR = */
/*        norm( abs(inv(A))* */
/*           ( abs(R) + NZ*EPS*( abs(A)*abs(X)+abs(B) ))) / norm(X) */

/*        where */
/*          norm(Z) is the magnitude of the largest component of Z */
/*          inv(A) is the inverse of A */
/*          abs(Z) is the componentwise absolute value of the matrix or */
/*             vector Z */
/*          NZ is the maximum number of nonzeros in any row of A, plus 1 */
/*          EPS is machine epsilon */

/*        The i-th component of abs(R)+NZ*EPS*(abs(A)*abs(X)+abs(B)) */
/*        is incremented by SAFE1 if the i-th component of */
/*        abs(A)*abs(X) + abs(B) is less than SAFE2. */

/*        Use DLACN2 to estimate the infinity-norm of the matrix */
/*           inv(A) * diag(W), */
/*        where W = abs(R) + NZ*EPS*( abs(A)*abs(X)+abs(B) ))) */

#line 391 "dsyrfs.f"
	i__2 = *n;
#line 391 "dsyrfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 392 "dsyrfs.f"
	    if (work[i__] > safe2) {
#line 393 "dsyrfs.f"
		work[i__] = (d__1 = work[*n + i__], abs(d__1)) + nz * eps * 
			work[i__];
#line 394 "dsyrfs.f"
	    } else {
#line 395 "dsyrfs.f"
		work[i__] = (d__1 = work[*n + i__], abs(d__1)) + nz * eps * 
			work[i__] + safe1;
#line 396 "dsyrfs.f"
	    }
#line 397 "dsyrfs.f"
/* L90: */
#line 397 "dsyrfs.f"
	}

#line 399 "dsyrfs.f"
	kase = 0;
#line 400 "dsyrfs.f"
L100:
#line 401 "dsyrfs.f"
	dlacn2_(n, &work[(*n << 1) + 1], &work[*n + 1], &iwork[1], &ferr[j], &
		kase, isave);
#line 403 "dsyrfs.f"
	if (kase != 0) {
#line 404 "dsyrfs.f"
	    if (kase == 1) {

/*              Multiply by diag(W)*inv(A**T). */

#line 408 "dsyrfs.f"
		dsytrs_(uplo, n, &c__1, &af[af_offset], ldaf, &ipiv[1], &work[
			*n + 1], n, info, (ftnlen)1);
#line 410 "dsyrfs.f"
		i__2 = *n;
#line 410 "dsyrfs.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 411 "dsyrfs.f"
		    work[*n + i__] = work[i__] * work[*n + i__];
#line 412 "dsyrfs.f"
/* L110: */
#line 412 "dsyrfs.f"
		}
#line 413 "dsyrfs.f"
	    } else if (kase == 2) {

/*              Multiply by inv(A)*diag(W). */

#line 417 "dsyrfs.f"
		i__2 = *n;
#line 417 "dsyrfs.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 418 "dsyrfs.f"
		    work[*n + i__] = work[i__] * work[*n + i__];
#line 419 "dsyrfs.f"
/* L120: */
#line 419 "dsyrfs.f"
		}
#line 420 "dsyrfs.f"
		dsytrs_(uplo, n, &c__1, &af[af_offset], ldaf, &ipiv[1], &work[
			*n + 1], n, info, (ftnlen)1);
#line 422 "dsyrfs.f"
	    }
#line 423 "dsyrfs.f"
	    goto L100;
#line 424 "dsyrfs.f"
	}

/*        Normalize error. */

#line 428 "dsyrfs.f"
	lstres = 0.;
#line 429 "dsyrfs.f"
	i__2 = *n;
#line 429 "dsyrfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 430 "dsyrfs.f"
	    d__2 = lstres, d__3 = (d__1 = x[i__ + j * x_dim1], abs(d__1));
#line 430 "dsyrfs.f"
	    lstres = max(d__2,d__3);
#line 431 "dsyrfs.f"
/* L130: */
#line 431 "dsyrfs.f"
	}
#line 432 "dsyrfs.f"
	if (lstres != 0.) {
#line 432 "dsyrfs.f"
	    ferr[j] /= lstres;
#line 432 "dsyrfs.f"
	}

#line 435 "dsyrfs.f"
/* L140: */
#line 435 "dsyrfs.f"
    }

#line 437 "dsyrfs.f"
    return 0;

/*     End of DSYRFS */

} /* dsyrfs_ */

