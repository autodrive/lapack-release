#line 1 "sporfs.f"
/* sporfs.f -- translated by f2c (version 20100827).
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

#line 1 "sporfs.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b12 = -1.;
static doublereal c_b14 = 1.;

/* > \brief \b SPORFS */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SPORFS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sporfs.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sporfs.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sporfs.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SPORFS( UPLO, N, NRHS, A, LDA, AF, LDAF, B, LDB, X, */
/*                          LDX, FERR, BERR, WORK, IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDA, LDAF, LDB, LDX, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       REAL               A( LDA, * ), AF( LDAF, * ), B( LDB, * ), */
/*      $                   BERR( * ), FERR( * ), WORK( * ), X( LDX, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SPORFS improves the computed solution to a system of linear */
/* > equations when the coefficient matrix is symmetric positive definite, */
/* > and provides error bounds and backward error estimates for the */
/* > solution. */
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
/* >          A is REAL array, dimension (LDA,N) */
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
/* >          AF is REAL array, dimension (LDAF,N) */
/* >          The triangular factor U or L from the Cholesky factorization */
/* >          A = U**T*U or A = L*L**T, as computed by SPOTRF. */
/* > \endverbatim */
/* > */
/* > \param[in] LDAF */
/* > \verbatim */
/* >          LDAF is INTEGER */
/* >          The leading dimension of the array AF.  LDAF >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] B */
/* > \verbatim */
/* >          B is REAL array, dimension (LDB,NRHS) */
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
/* >          X is REAL array, dimension (LDX,NRHS) */
/* >          On entry, the solution matrix X, as computed by SPOTRS. */
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
/* >          FERR is REAL array, dimension (NRHS) */
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
/* >          BERR is REAL array, dimension (NRHS) */
/* >          The componentwise relative backward error of each solution */
/* >          vector X(j) (i.e., the smallest relative change in */
/* >          any element of A or B that makes X(j) an exact solution). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is REAL array, dimension (3*N) */
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

/* > \date December 2016 */

/* > \ingroup realPOcomputational */

/*  ===================================================================== */
/* Subroutine */ int sporfs_(char *uplo, integer *n, integer *nrhs, 
	doublereal *a, integer *lda, doublereal *af, integer *ldaf, 
	doublereal *b, integer *ldb, doublereal *x, integer *ldx, doublereal *
	ferr, doublereal *berr, doublereal *work, integer *iwork, integer *
	info, ftnlen uplo_len)
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
    static integer isave[3], count;
    static logical upper;
    extern /* Subroutine */ int scopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), saxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *), ssymv_(char *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen), 
	    slacn2_(integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, integer *);
    extern doublereal slamch_(char *, ftnlen);
    static doublereal safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal lstres;
    extern /* Subroutine */ int spotrs_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen);


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
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 238 "sporfs.f"
    /* Parameter adjustments */
#line 238 "sporfs.f"
    a_dim1 = *lda;
#line 238 "sporfs.f"
    a_offset = 1 + a_dim1;
#line 238 "sporfs.f"
    a -= a_offset;
#line 238 "sporfs.f"
    af_dim1 = *ldaf;
#line 238 "sporfs.f"
    af_offset = 1 + af_dim1;
#line 238 "sporfs.f"
    af -= af_offset;
#line 238 "sporfs.f"
    b_dim1 = *ldb;
#line 238 "sporfs.f"
    b_offset = 1 + b_dim1;
#line 238 "sporfs.f"
    b -= b_offset;
#line 238 "sporfs.f"
    x_dim1 = *ldx;
#line 238 "sporfs.f"
    x_offset = 1 + x_dim1;
#line 238 "sporfs.f"
    x -= x_offset;
#line 238 "sporfs.f"
    --ferr;
#line 238 "sporfs.f"
    --berr;
#line 238 "sporfs.f"
    --work;
#line 238 "sporfs.f"
    --iwork;
#line 238 "sporfs.f"

#line 238 "sporfs.f"
    /* Function Body */
#line 238 "sporfs.f"
    *info = 0;
#line 239 "sporfs.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 240 "sporfs.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 241 "sporfs.f"
	*info = -1;
#line 242 "sporfs.f"
    } else if (*n < 0) {
#line 243 "sporfs.f"
	*info = -2;
#line 244 "sporfs.f"
    } else if (*nrhs < 0) {
#line 245 "sporfs.f"
	*info = -3;
#line 246 "sporfs.f"
    } else if (*lda < max(1,*n)) {
#line 247 "sporfs.f"
	*info = -5;
#line 248 "sporfs.f"
    } else if (*ldaf < max(1,*n)) {
#line 249 "sporfs.f"
	*info = -7;
#line 250 "sporfs.f"
    } else if (*ldb < max(1,*n)) {
#line 251 "sporfs.f"
	*info = -9;
#line 252 "sporfs.f"
    } else if (*ldx < max(1,*n)) {
#line 253 "sporfs.f"
	*info = -11;
#line 254 "sporfs.f"
    }
#line 255 "sporfs.f"
    if (*info != 0) {
#line 256 "sporfs.f"
	i__1 = -(*info);
#line 256 "sporfs.f"
	xerbla_("SPORFS", &i__1, (ftnlen)6);
#line 257 "sporfs.f"
	return 0;
#line 258 "sporfs.f"
    }

/*     Quick return if possible */

#line 262 "sporfs.f"
    if (*n == 0 || *nrhs == 0) {
#line 263 "sporfs.f"
	i__1 = *nrhs;
#line 263 "sporfs.f"
	for (j = 1; j <= i__1; ++j) {
#line 264 "sporfs.f"
	    ferr[j] = 0.;
#line 265 "sporfs.f"
	    berr[j] = 0.;
#line 266 "sporfs.f"
/* L10: */
#line 266 "sporfs.f"
	}
#line 267 "sporfs.f"
	return 0;
#line 268 "sporfs.f"
    }

/*     NZ = maximum number of nonzero elements in each row of A, plus 1 */

#line 272 "sporfs.f"
    nz = *n + 1;
#line 273 "sporfs.f"
    eps = slamch_("Epsilon", (ftnlen)7);
#line 274 "sporfs.f"
    safmin = slamch_("Safe minimum", (ftnlen)12);
#line 275 "sporfs.f"
    safe1 = nz * safmin;
#line 276 "sporfs.f"
    safe2 = safe1 / eps;

/*     Do for each right hand side */

#line 280 "sporfs.f"
    i__1 = *nrhs;
#line 280 "sporfs.f"
    for (j = 1; j <= i__1; ++j) {

#line 282 "sporfs.f"
	count = 1;
#line 283 "sporfs.f"
	lstres = 3.;
#line 284 "sporfs.f"
L20:

/*        Loop until stopping criterion is satisfied. */

/*        Compute residual R = B - A * X */

#line 290 "sporfs.f"
	scopy_(n, &b[j * b_dim1 + 1], &c__1, &work[*n + 1], &c__1);
#line 291 "sporfs.f"
	ssymv_(uplo, n, &c_b12, &a[a_offset], lda, &x[j * x_dim1 + 1], &c__1, 
		&c_b14, &work[*n + 1], &c__1, (ftnlen)1);

/*        Compute componentwise relative backward error from formula */

/*        max(i) ( abs(R(i)) / ( abs(A)*abs(X) + abs(B) )(i) ) */

/*        where abs(Z) is the componentwise absolute value of the matrix */
/*        or vector Z.  If the i-th component of the denominator is less */
/*        than SAFE2, then SAFE1 is added to the i-th components of the */
/*        numerator and denominator before dividing. */

#line 303 "sporfs.f"
	i__2 = *n;
#line 303 "sporfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 304 "sporfs.f"
	    work[i__] = (d__1 = b[i__ + j * b_dim1], abs(d__1));
#line 305 "sporfs.f"
/* L30: */
#line 305 "sporfs.f"
	}

/*        Compute abs(A)*abs(X) + abs(B). */

#line 309 "sporfs.f"
	if (upper) {
#line 310 "sporfs.f"
	    i__2 = *n;
#line 310 "sporfs.f"
	    for (k = 1; k <= i__2; ++k) {
#line 311 "sporfs.f"
		s = 0.;
#line 312 "sporfs.f"
		xk = (d__1 = x[k + j * x_dim1], abs(d__1));
#line 313 "sporfs.f"
		i__3 = k - 1;
#line 313 "sporfs.f"
		for (i__ = 1; i__ <= i__3; ++i__) {
#line 314 "sporfs.f"
		    work[i__] += (d__1 = a[i__ + k * a_dim1], abs(d__1)) * xk;
#line 315 "sporfs.f"
		    s += (d__1 = a[i__ + k * a_dim1], abs(d__1)) * (d__2 = x[
			    i__ + j * x_dim1], abs(d__2));
#line 316 "sporfs.f"
/* L40: */
#line 316 "sporfs.f"
		}
#line 317 "sporfs.f"
		work[k] = work[k] + (d__1 = a[k + k * a_dim1], abs(d__1)) * 
			xk + s;
#line 318 "sporfs.f"
/* L50: */
#line 318 "sporfs.f"
	    }
#line 319 "sporfs.f"
	} else {
#line 320 "sporfs.f"
	    i__2 = *n;
#line 320 "sporfs.f"
	    for (k = 1; k <= i__2; ++k) {
#line 321 "sporfs.f"
		s = 0.;
#line 322 "sporfs.f"
		xk = (d__1 = x[k + j * x_dim1], abs(d__1));
#line 323 "sporfs.f"
		work[k] += (d__1 = a[k + k * a_dim1], abs(d__1)) * xk;
#line 324 "sporfs.f"
		i__3 = *n;
#line 324 "sporfs.f"
		for (i__ = k + 1; i__ <= i__3; ++i__) {
#line 325 "sporfs.f"
		    work[i__] += (d__1 = a[i__ + k * a_dim1], abs(d__1)) * xk;
#line 326 "sporfs.f"
		    s += (d__1 = a[i__ + k * a_dim1], abs(d__1)) * (d__2 = x[
			    i__ + j * x_dim1], abs(d__2));
#line 327 "sporfs.f"
/* L60: */
#line 327 "sporfs.f"
		}
#line 328 "sporfs.f"
		work[k] += s;
#line 329 "sporfs.f"
/* L70: */
#line 329 "sporfs.f"
	    }
#line 330 "sporfs.f"
	}
#line 331 "sporfs.f"
	s = 0.;
#line 332 "sporfs.f"
	i__2 = *n;
#line 332 "sporfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 333 "sporfs.f"
	    if (work[i__] > safe2) {
/* Computing MAX */
#line 334 "sporfs.f"
		d__2 = s, d__3 = (d__1 = work[*n + i__], abs(d__1)) / work[
			i__];
#line 334 "sporfs.f"
		s = max(d__2,d__3);
#line 335 "sporfs.f"
	    } else {
/* Computing MAX */
#line 336 "sporfs.f"
		d__2 = s, d__3 = ((d__1 = work[*n + i__], abs(d__1)) + safe1) 
			/ (work[i__] + safe1);
#line 336 "sporfs.f"
		s = max(d__2,d__3);
#line 338 "sporfs.f"
	    }
#line 339 "sporfs.f"
/* L80: */
#line 339 "sporfs.f"
	}
#line 340 "sporfs.f"
	berr[j] = s;

/*        Test stopping criterion. Continue iterating if */
/*           1) The residual BERR(J) is larger than machine epsilon, and */
/*           2) BERR(J) decreased by at least a factor of 2 during the */
/*              last iteration, and */
/*           3) At most ITMAX iterations tried. */

#line 348 "sporfs.f"
	if (berr[j] > eps && berr[j] * 2. <= lstres && count <= 5) {

/*           Update solution and try again. */

#line 353 "sporfs.f"
	    spotrs_(uplo, n, &c__1, &af[af_offset], ldaf, &work[*n + 1], n, 
		    info, (ftnlen)1);
#line 354 "sporfs.f"
	    saxpy_(n, &c_b14, &work[*n + 1], &c__1, &x[j * x_dim1 + 1], &c__1)
		    ;
#line 355 "sporfs.f"
	    lstres = berr[j];
#line 356 "sporfs.f"
	    ++count;
#line 357 "sporfs.f"
	    goto L20;
#line 358 "sporfs.f"
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

/*        Use SLACN2 to estimate the infinity-norm of the matrix */
/*           inv(A) * diag(W), */
/*        where W = abs(R) + NZ*EPS*( abs(A)*abs(X)+abs(B) ))) */

#line 382 "sporfs.f"
	i__2 = *n;
#line 382 "sporfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 383 "sporfs.f"
	    if (work[i__] > safe2) {
#line 384 "sporfs.f"
		work[i__] = (d__1 = work[*n + i__], abs(d__1)) + nz * eps * 
			work[i__];
#line 385 "sporfs.f"
	    } else {
#line 386 "sporfs.f"
		work[i__] = (d__1 = work[*n + i__], abs(d__1)) + nz * eps * 
			work[i__] + safe1;
#line 387 "sporfs.f"
	    }
#line 388 "sporfs.f"
/* L90: */
#line 388 "sporfs.f"
	}

#line 390 "sporfs.f"
	kase = 0;
#line 391 "sporfs.f"
L100:
#line 392 "sporfs.f"
	slacn2_(n, &work[(*n << 1) + 1], &work[*n + 1], &iwork[1], &ferr[j], &
		kase, isave);
#line 394 "sporfs.f"
	if (kase != 0) {
#line 395 "sporfs.f"
	    if (kase == 1) {

/*              Multiply by diag(W)*inv(A**T). */

#line 399 "sporfs.f"
		spotrs_(uplo, n, &c__1, &af[af_offset], ldaf, &work[*n + 1], 
			n, info, (ftnlen)1);
#line 400 "sporfs.f"
		i__2 = *n;
#line 400 "sporfs.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 401 "sporfs.f"
		    work[*n + i__] = work[i__] * work[*n + i__];
#line 402 "sporfs.f"
/* L110: */
#line 402 "sporfs.f"
		}
#line 403 "sporfs.f"
	    } else if (kase == 2) {

/*              Multiply by inv(A)*diag(W). */

#line 407 "sporfs.f"
		i__2 = *n;
#line 407 "sporfs.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 408 "sporfs.f"
		    work[*n + i__] = work[i__] * work[*n + i__];
#line 409 "sporfs.f"
/* L120: */
#line 409 "sporfs.f"
		}
#line 410 "sporfs.f"
		spotrs_(uplo, n, &c__1, &af[af_offset], ldaf, &work[*n + 1], 
			n, info, (ftnlen)1);
#line 411 "sporfs.f"
	    }
#line 412 "sporfs.f"
	    goto L100;
#line 413 "sporfs.f"
	}

/*        Normalize error. */

#line 417 "sporfs.f"
	lstres = 0.;
#line 418 "sporfs.f"
	i__2 = *n;
#line 418 "sporfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 419 "sporfs.f"
	    d__2 = lstres, d__3 = (d__1 = x[i__ + j * x_dim1], abs(d__1));
#line 419 "sporfs.f"
	    lstres = max(d__2,d__3);
#line 420 "sporfs.f"
/* L130: */
#line 420 "sporfs.f"
	}
#line 421 "sporfs.f"
	if (lstres != 0.) {
#line 421 "sporfs.f"
	    ferr[j] /= lstres;
#line 421 "sporfs.f"
	}

#line 424 "sporfs.f"
/* L140: */
#line 424 "sporfs.f"
    }

#line 426 "sporfs.f"
    return 0;

/*     End of SPORFS */

} /* sporfs_ */

