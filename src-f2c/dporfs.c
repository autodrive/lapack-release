#line 1 "dporfs.f"
/* dporfs.f -- translated by f2c (version 20100827).
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

#line 1 "dporfs.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b12 = -1.;
static doublereal c_b14 = 1.;

/* > \brief \b DPORFS */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DPORFS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dporfs.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dporfs.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dporfs.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DPORFS( UPLO, N, NRHS, A, LDA, AF, LDAF, B, LDB, X, */
/*                          LDX, FERR, BERR, WORK, IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDA, LDAF, LDB, LDX, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       DOUBLE PRECISION   A( LDA, * ), AF( LDAF, * ), B( LDB, * ), */
/*      $                   BERR( * ), FERR( * ), WORK( * ), X( LDX, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DPORFS improves the computed solution to a system of linear */
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
/* >          The triangular factor U or L from the Cholesky factorization */
/* >          A = U**T*U or A = L*L**T, as computed by DPOTRF. */
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
/* >          On entry, the solution matrix X, as computed by DPOTRS. */
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

/* > \ingroup doublePOcomputational */

/*  ===================================================================== */
/* Subroutine */ int dporfs_(char *uplo, integer *n, integer *nrhs, 
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
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), dpotrs_(
	    char *, integer *, integer *, doublereal *, integer *, doublereal 
	    *, integer *, integer *, ftnlen);
    static doublereal lstres;


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

#line 238 "dporfs.f"
    /* Parameter adjustments */
#line 238 "dporfs.f"
    a_dim1 = *lda;
#line 238 "dporfs.f"
    a_offset = 1 + a_dim1;
#line 238 "dporfs.f"
    a -= a_offset;
#line 238 "dporfs.f"
    af_dim1 = *ldaf;
#line 238 "dporfs.f"
    af_offset = 1 + af_dim1;
#line 238 "dporfs.f"
    af -= af_offset;
#line 238 "dporfs.f"
    b_dim1 = *ldb;
#line 238 "dporfs.f"
    b_offset = 1 + b_dim1;
#line 238 "dporfs.f"
    b -= b_offset;
#line 238 "dporfs.f"
    x_dim1 = *ldx;
#line 238 "dporfs.f"
    x_offset = 1 + x_dim1;
#line 238 "dporfs.f"
    x -= x_offset;
#line 238 "dporfs.f"
    --ferr;
#line 238 "dporfs.f"
    --berr;
#line 238 "dporfs.f"
    --work;
#line 238 "dporfs.f"
    --iwork;
#line 238 "dporfs.f"

#line 238 "dporfs.f"
    /* Function Body */
#line 238 "dporfs.f"
    *info = 0;
#line 239 "dporfs.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 240 "dporfs.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 241 "dporfs.f"
	*info = -1;
#line 242 "dporfs.f"
    } else if (*n < 0) {
#line 243 "dporfs.f"
	*info = -2;
#line 244 "dporfs.f"
    } else if (*nrhs < 0) {
#line 245 "dporfs.f"
	*info = -3;
#line 246 "dporfs.f"
    } else if (*lda < max(1,*n)) {
#line 247 "dporfs.f"
	*info = -5;
#line 248 "dporfs.f"
    } else if (*ldaf < max(1,*n)) {
#line 249 "dporfs.f"
	*info = -7;
#line 250 "dporfs.f"
    } else if (*ldb < max(1,*n)) {
#line 251 "dporfs.f"
	*info = -9;
#line 252 "dporfs.f"
    } else if (*ldx < max(1,*n)) {
#line 253 "dporfs.f"
	*info = -11;
#line 254 "dporfs.f"
    }
#line 255 "dporfs.f"
    if (*info != 0) {
#line 256 "dporfs.f"
	i__1 = -(*info);
#line 256 "dporfs.f"
	xerbla_("DPORFS", &i__1, (ftnlen)6);
#line 257 "dporfs.f"
	return 0;
#line 258 "dporfs.f"
    }

/*     Quick return if possible */

#line 262 "dporfs.f"
    if (*n == 0 || *nrhs == 0) {
#line 263 "dporfs.f"
	i__1 = *nrhs;
#line 263 "dporfs.f"
	for (j = 1; j <= i__1; ++j) {
#line 264 "dporfs.f"
	    ferr[j] = 0.;
#line 265 "dporfs.f"
	    berr[j] = 0.;
#line 266 "dporfs.f"
/* L10: */
#line 266 "dporfs.f"
	}
#line 267 "dporfs.f"
	return 0;
#line 268 "dporfs.f"
    }

/*     NZ = maximum number of nonzero elements in each row of A, plus 1 */

#line 272 "dporfs.f"
    nz = *n + 1;
#line 273 "dporfs.f"
    eps = dlamch_("Epsilon", (ftnlen)7);
#line 274 "dporfs.f"
    safmin = dlamch_("Safe minimum", (ftnlen)12);
#line 275 "dporfs.f"
    safe1 = nz * safmin;
#line 276 "dporfs.f"
    safe2 = safe1 / eps;

/*     Do for each right hand side */

#line 280 "dporfs.f"
    i__1 = *nrhs;
#line 280 "dporfs.f"
    for (j = 1; j <= i__1; ++j) {

#line 282 "dporfs.f"
	count = 1;
#line 283 "dporfs.f"
	lstres = 3.;
#line 284 "dporfs.f"
L20:

/*        Loop until stopping criterion is satisfied. */

/*        Compute residual R = B - A * X */

#line 290 "dporfs.f"
	dcopy_(n, &b[j * b_dim1 + 1], &c__1, &work[*n + 1], &c__1);
#line 291 "dporfs.f"
	dsymv_(uplo, n, &c_b12, &a[a_offset], lda, &x[j * x_dim1 + 1], &c__1, 
		&c_b14, &work[*n + 1], &c__1, (ftnlen)1);

/*        Compute componentwise relative backward error from formula */

/*        max(i) ( abs(R(i)) / ( abs(A)*abs(X) + abs(B) )(i) ) */

/*        where abs(Z) is the componentwise absolute value of the matrix */
/*        or vector Z.  If the i-th component of the denominator is less */
/*        than SAFE2, then SAFE1 is added to the i-th components of the */
/*        numerator and denominator before dividing. */

#line 303 "dporfs.f"
	i__2 = *n;
#line 303 "dporfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 304 "dporfs.f"
	    work[i__] = (d__1 = b[i__ + j * b_dim1], abs(d__1));
#line 305 "dporfs.f"
/* L30: */
#line 305 "dporfs.f"
	}

/*        Compute abs(A)*abs(X) + abs(B). */

#line 309 "dporfs.f"
	if (upper) {
#line 310 "dporfs.f"
	    i__2 = *n;
#line 310 "dporfs.f"
	    for (k = 1; k <= i__2; ++k) {
#line 311 "dporfs.f"
		s = 0.;
#line 312 "dporfs.f"
		xk = (d__1 = x[k + j * x_dim1], abs(d__1));
#line 313 "dporfs.f"
		i__3 = k - 1;
#line 313 "dporfs.f"
		for (i__ = 1; i__ <= i__3; ++i__) {
#line 314 "dporfs.f"
		    work[i__] += (d__1 = a[i__ + k * a_dim1], abs(d__1)) * xk;
#line 315 "dporfs.f"
		    s += (d__1 = a[i__ + k * a_dim1], abs(d__1)) * (d__2 = x[
			    i__ + j * x_dim1], abs(d__2));
#line 316 "dporfs.f"
/* L40: */
#line 316 "dporfs.f"
		}
#line 317 "dporfs.f"
		work[k] = work[k] + (d__1 = a[k + k * a_dim1], abs(d__1)) * 
			xk + s;
#line 318 "dporfs.f"
/* L50: */
#line 318 "dporfs.f"
	    }
#line 319 "dporfs.f"
	} else {
#line 320 "dporfs.f"
	    i__2 = *n;
#line 320 "dporfs.f"
	    for (k = 1; k <= i__2; ++k) {
#line 321 "dporfs.f"
		s = 0.;
#line 322 "dporfs.f"
		xk = (d__1 = x[k + j * x_dim1], abs(d__1));
#line 323 "dporfs.f"
		work[k] += (d__1 = a[k + k * a_dim1], abs(d__1)) * xk;
#line 324 "dporfs.f"
		i__3 = *n;
#line 324 "dporfs.f"
		for (i__ = k + 1; i__ <= i__3; ++i__) {
#line 325 "dporfs.f"
		    work[i__] += (d__1 = a[i__ + k * a_dim1], abs(d__1)) * xk;
#line 326 "dporfs.f"
		    s += (d__1 = a[i__ + k * a_dim1], abs(d__1)) * (d__2 = x[
			    i__ + j * x_dim1], abs(d__2));
#line 327 "dporfs.f"
/* L60: */
#line 327 "dporfs.f"
		}
#line 328 "dporfs.f"
		work[k] += s;
#line 329 "dporfs.f"
/* L70: */
#line 329 "dporfs.f"
	    }
#line 330 "dporfs.f"
	}
#line 331 "dporfs.f"
	s = 0.;
#line 332 "dporfs.f"
	i__2 = *n;
#line 332 "dporfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 333 "dporfs.f"
	    if (work[i__] > safe2) {
/* Computing MAX */
#line 334 "dporfs.f"
		d__2 = s, d__3 = (d__1 = work[*n + i__], abs(d__1)) / work[
			i__];
#line 334 "dporfs.f"
		s = max(d__2,d__3);
#line 335 "dporfs.f"
	    } else {
/* Computing MAX */
#line 336 "dporfs.f"
		d__2 = s, d__3 = ((d__1 = work[*n + i__], abs(d__1)) + safe1) 
			/ (work[i__] + safe1);
#line 336 "dporfs.f"
		s = max(d__2,d__3);
#line 338 "dporfs.f"
	    }
#line 339 "dporfs.f"
/* L80: */
#line 339 "dporfs.f"
	}
#line 340 "dporfs.f"
	berr[j] = s;

/*        Test stopping criterion. Continue iterating if */
/*           1) The residual BERR(J) is larger than machine epsilon, and */
/*           2) BERR(J) decreased by at least a factor of 2 during the */
/*              last iteration, and */
/*           3) At most ITMAX iterations tried. */

#line 348 "dporfs.f"
	if (berr[j] > eps && berr[j] * 2. <= lstres && count <= 5) {

/*           Update solution and try again. */

#line 353 "dporfs.f"
	    dpotrs_(uplo, n, &c__1, &af[af_offset], ldaf, &work[*n + 1], n, 
		    info, (ftnlen)1);
#line 354 "dporfs.f"
	    daxpy_(n, &c_b14, &work[*n + 1], &c__1, &x[j * x_dim1 + 1], &c__1)
		    ;
#line 355 "dporfs.f"
	    lstres = berr[j];
#line 356 "dporfs.f"
	    ++count;
#line 357 "dporfs.f"
	    goto L20;
#line 358 "dporfs.f"
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

#line 382 "dporfs.f"
	i__2 = *n;
#line 382 "dporfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 383 "dporfs.f"
	    if (work[i__] > safe2) {
#line 384 "dporfs.f"
		work[i__] = (d__1 = work[*n + i__], abs(d__1)) + nz * eps * 
			work[i__];
#line 385 "dporfs.f"
	    } else {
#line 386 "dporfs.f"
		work[i__] = (d__1 = work[*n + i__], abs(d__1)) + nz * eps * 
			work[i__] + safe1;
#line 387 "dporfs.f"
	    }
#line 388 "dporfs.f"
/* L90: */
#line 388 "dporfs.f"
	}

#line 390 "dporfs.f"
	kase = 0;
#line 391 "dporfs.f"
L100:
#line 392 "dporfs.f"
	dlacn2_(n, &work[(*n << 1) + 1], &work[*n + 1], &iwork[1], &ferr[j], &
		kase, isave);
#line 394 "dporfs.f"
	if (kase != 0) {
#line 395 "dporfs.f"
	    if (kase == 1) {

/*              Multiply by diag(W)*inv(A**T). */

#line 399 "dporfs.f"
		dpotrs_(uplo, n, &c__1, &af[af_offset], ldaf, &work[*n + 1], 
			n, info, (ftnlen)1);
#line 400 "dporfs.f"
		i__2 = *n;
#line 400 "dporfs.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 401 "dporfs.f"
		    work[*n + i__] = work[i__] * work[*n + i__];
#line 402 "dporfs.f"
/* L110: */
#line 402 "dporfs.f"
		}
#line 403 "dporfs.f"
	    } else if (kase == 2) {

/*              Multiply by inv(A)*diag(W). */

#line 407 "dporfs.f"
		i__2 = *n;
#line 407 "dporfs.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 408 "dporfs.f"
		    work[*n + i__] = work[i__] * work[*n + i__];
#line 409 "dporfs.f"
/* L120: */
#line 409 "dporfs.f"
		}
#line 410 "dporfs.f"
		dpotrs_(uplo, n, &c__1, &af[af_offset], ldaf, &work[*n + 1], 
			n, info, (ftnlen)1);
#line 411 "dporfs.f"
	    }
#line 412 "dporfs.f"
	    goto L100;
#line 413 "dporfs.f"
	}

/*        Normalize error. */

#line 417 "dporfs.f"
	lstres = 0.;
#line 418 "dporfs.f"
	i__2 = *n;
#line 418 "dporfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 419 "dporfs.f"
	    d__2 = lstres, d__3 = (d__1 = x[i__ + j * x_dim1], abs(d__1));
#line 419 "dporfs.f"
	    lstres = max(d__2,d__3);
#line 420 "dporfs.f"
/* L130: */
#line 420 "dporfs.f"
	}
#line 421 "dporfs.f"
	if (lstres != 0.) {
#line 421 "dporfs.f"
	    ferr[j] /= lstres;
#line 421 "dporfs.f"
	}

#line 424 "dporfs.f"
/* L140: */
#line 424 "dporfs.f"
    }

#line 426 "dporfs.f"
    return 0;

/*     End of DPORFS */

} /* dporfs_ */

