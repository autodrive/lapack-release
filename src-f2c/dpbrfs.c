#line 1 "dpbrfs.f"
/* dpbrfs.f -- translated by f2c (version 20100827).
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

#line 1 "dpbrfs.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b12 = -1.;
static doublereal c_b14 = 1.;

/* > \brief \b DPBRFS */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DPBRFS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dpbrfs.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dpbrfs.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dpbrfs.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DPBRFS( UPLO, N, KD, NRHS, AB, LDAB, AFB, LDAFB, B, */
/*                          LDB, X, LDX, FERR, BERR, WORK, IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, KD, LDAB, LDAFB, LDB, LDX, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       DOUBLE PRECISION   AB( LDAB, * ), AFB( LDAFB, * ), B( LDB, * ), */
/*      $                   BERR( * ), FERR( * ), WORK( * ), X( LDX, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DPBRFS improves the computed solution to a system of linear */
/* > equations when the coefficient matrix is symmetric positive definite */
/* > and banded, and provides error bounds and backward error estimates */
/* > for the solution. */
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
/* > \param[in] KD */
/* > \verbatim */
/* >          KD is INTEGER */
/* >          The number of superdiagonals of the matrix A if UPLO = 'U', */
/* >          or the number of subdiagonals if UPLO = 'L'.  KD >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NRHS */
/* > \verbatim */
/* >          NRHS is INTEGER */
/* >          The number of right hand sides, i.e., the number of columns */
/* >          of the matrices B and X.  NRHS >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] AB */
/* > \verbatim */
/* >          AB is DOUBLE PRECISION array, dimension (LDAB,N) */
/* >          The upper or lower triangle of the symmetric band matrix A, */
/* >          stored in the first KD+1 rows of the array.  The j-th column */
/* >          of A is stored in the j-th column of the array AB as follows: */
/* >          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j; */
/* >          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd). */
/* > \endverbatim */
/* > */
/* > \param[in] LDAB */
/* > \verbatim */
/* >          LDAB is INTEGER */
/* >          The leading dimension of the array AB.  LDAB >= KD+1. */
/* > \endverbatim */
/* > */
/* > \param[in] AFB */
/* > \verbatim */
/* >          AFB is DOUBLE PRECISION array, dimension (LDAFB,N) */
/* >          The triangular factor U or L from the Cholesky factorization */
/* >          A = U**T*U or A = L*L**T of the band matrix A as computed by */
/* >          DPBTRF, in the same storage format as A (see AB). */
/* > \endverbatim */
/* > */
/* > \param[in] LDAFB */
/* > \verbatim */
/* >          LDAFB is INTEGER */
/* >          The leading dimension of the array AFB.  LDAFB >= KD+1. */
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
/* >          On entry, the solution matrix X, as computed by DPBTRS. */
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

/* > \ingroup doubleOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int dpbrfs_(char *uplo, integer *n, integer *kd, integer *
	nrhs, doublereal *ab, integer *ldab, doublereal *afb, integer *ldafb, 
	doublereal *b, integer *ldb, doublereal *x, integer *ldx, doublereal *
	ferr, doublereal *berr, doublereal *work, integer *iwork, integer *
	info, ftnlen uplo_len)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, afb_dim1, afb_offset, b_dim1, b_offset, 
	    x_dim1, x_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2, d__3;

    /* Local variables */
    static integer i__, j, k, l;
    static doublereal s, xk;
    static integer nz;
    static doublereal eps;
    static integer kase;
    static doublereal safe1, safe2;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer isave[3];
    extern /* Subroutine */ int dsbmv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen), dcopy_(integer *, 
	    doublereal *, integer *, doublereal *, integer *), daxpy_(integer 
	    *, doublereal *, doublereal *, integer *, doublereal *, integer *)
	    ;
    static integer count;
    static logical upper;
    extern /* Subroutine */ int dlacn2_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, integer *);
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), dpbtrs_(
	    char *, integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, ftnlen);
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

#line 244 "dpbrfs.f"
    /* Parameter adjustments */
#line 244 "dpbrfs.f"
    ab_dim1 = *ldab;
#line 244 "dpbrfs.f"
    ab_offset = 1 + ab_dim1;
#line 244 "dpbrfs.f"
    ab -= ab_offset;
#line 244 "dpbrfs.f"
    afb_dim1 = *ldafb;
#line 244 "dpbrfs.f"
    afb_offset = 1 + afb_dim1;
#line 244 "dpbrfs.f"
    afb -= afb_offset;
#line 244 "dpbrfs.f"
    b_dim1 = *ldb;
#line 244 "dpbrfs.f"
    b_offset = 1 + b_dim1;
#line 244 "dpbrfs.f"
    b -= b_offset;
#line 244 "dpbrfs.f"
    x_dim1 = *ldx;
#line 244 "dpbrfs.f"
    x_offset = 1 + x_dim1;
#line 244 "dpbrfs.f"
    x -= x_offset;
#line 244 "dpbrfs.f"
    --ferr;
#line 244 "dpbrfs.f"
    --berr;
#line 244 "dpbrfs.f"
    --work;
#line 244 "dpbrfs.f"
    --iwork;
#line 244 "dpbrfs.f"

#line 244 "dpbrfs.f"
    /* Function Body */
#line 244 "dpbrfs.f"
    *info = 0;
#line 245 "dpbrfs.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 246 "dpbrfs.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 247 "dpbrfs.f"
	*info = -1;
#line 248 "dpbrfs.f"
    } else if (*n < 0) {
#line 249 "dpbrfs.f"
	*info = -2;
#line 250 "dpbrfs.f"
    } else if (*kd < 0) {
#line 251 "dpbrfs.f"
	*info = -3;
#line 252 "dpbrfs.f"
    } else if (*nrhs < 0) {
#line 253 "dpbrfs.f"
	*info = -4;
#line 254 "dpbrfs.f"
    } else if (*ldab < *kd + 1) {
#line 255 "dpbrfs.f"
	*info = -6;
#line 256 "dpbrfs.f"
    } else if (*ldafb < *kd + 1) {
#line 257 "dpbrfs.f"
	*info = -8;
#line 258 "dpbrfs.f"
    } else if (*ldb < max(1,*n)) {
#line 259 "dpbrfs.f"
	*info = -10;
#line 260 "dpbrfs.f"
    } else if (*ldx < max(1,*n)) {
#line 261 "dpbrfs.f"
	*info = -12;
#line 262 "dpbrfs.f"
    }
#line 263 "dpbrfs.f"
    if (*info != 0) {
#line 264 "dpbrfs.f"
	i__1 = -(*info);
#line 264 "dpbrfs.f"
	xerbla_("DPBRFS", &i__1, (ftnlen)6);
#line 265 "dpbrfs.f"
	return 0;
#line 266 "dpbrfs.f"
    }

/*     Quick return if possible */

#line 270 "dpbrfs.f"
    if (*n == 0 || *nrhs == 0) {
#line 271 "dpbrfs.f"
	i__1 = *nrhs;
#line 271 "dpbrfs.f"
	for (j = 1; j <= i__1; ++j) {
#line 272 "dpbrfs.f"
	    ferr[j] = 0.;
#line 273 "dpbrfs.f"
	    berr[j] = 0.;
#line 274 "dpbrfs.f"
/* L10: */
#line 274 "dpbrfs.f"
	}
#line 275 "dpbrfs.f"
	return 0;
#line 276 "dpbrfs.f"
    }

/*     NZ = maximum number of nonzero elements in each row of A, plus 1 */

/* Computing MIN */
#line 280 "dpbrfs.f"
    i__1 = *n + 1, i__2 = (*kd << 1) + 2;
#line 280 "dpbrfs.f"
    nz = min(i__1,i__2);
#line 281 "dpbrfs.f"
    eps = dlamch_("Epsilon", (ftnlen)7);
#line 282 "dpbrfs.f"
    safmin = dlamch_("Safe minimum", (ftnlen)12);
#line 283 "dpbrfs.f"
    safe1 = nz * safmin;
#line 284 "dpbrfs.f"
    safe2 = safe1 / eps;

/*     Do for each right hand side */

#line 288 "dpbrfs.f"
    i__1 = *nrhs;
#line 288 "dpbrfs.f"
    for (j = 1; j <= i__1; ++j) {

#line 290 "dpbrfs.f"
	count = 1;
#line 291 "dpbrfs.f"
	lstres = 3.;
#line 292 "dpbrfs.f"
L20:

/*        Loop until stopping criterion is satisfied. */

/*        Compute residual R = B - A * X */

#line 298 "dpbrfs.f"
	dcopy_(n, &b[j * b_dim1 + 1], &c__1, &work[*n + 1], &c__1);
#line 299 "dpbrfs.f"
	dsbmv_(uplo, n, kd, &c_b12, &ab[ab_offset], ldab, &x[j * x_dim1 + 1], 
		&c__1, &c_b14, &work[*n + 1], &c__1, (ftnlen)1);

/*        Compute componentwise relative backward error from formula */

/*        max(i) ( abs(R(i)) / ( abs(A)*abs(X) + abs(B) )(i) ) */

/*        where abs(Z) is the componentwise absolute value of the matrix */
/*        or vector Z.  If the i-th component of the denominator is less */
/*        than SAFE2, then SAFE1 is added to the i-th components of the */
/*        numerator and denominator before dividing. */

#line 311 "dpbrfs.f"
	i__2 = *n;
#line 311 "dpbrfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 312 "dpbrfs.f"
	    work[i__] = (d__1 = b[i__ + j * b_dim1], abs(d__1));
#line 313 "dpbrfs.f"
/* L30: */
#line 313 "dpbrfs.f"
	}

/*        Compute abs(A)*abs(X) + abs(B). */

#line 317 "dpbrfs.f"
	if (upper) {
#line 318 "dpbrfs.f"
	    i__2 = *n;
#line 318 "dpbrfs.f"
	    for (k = 1; k <= i__2; ++k) {
#line 319 "dpbrfs.f"
		s = 0.;
#line 320 "dpbrfs.f"
		xk = (d__1 = x[k + j * x_dim1], abs(d__1));
#line 321 "dpbrfs.f"
		l = *kd + 1 - k;
/* Computing MAX */
#line 322 "dpbrfs.f"
		i__3 = 1, i__4 = k - *kd;
#line 322 "dpbrfs.f"
		i__5 = k - 1;
#line 322 "dpbrfs.f"
		for (i__ = max(i__3,i__4); i__ <= i__5; ++i__) {
#line 323 "dpbrfs.f"
		    work[i__] += (d__1 = ab[l + i__ + k * ab_dim1], abs(d__1))
			     * xk;
#line 324 "dpbrfs.f"
		    s += (d__1 = ab[l + i__ + k * ab_dim1], abs(d__1)) * (
			    d__2 = x[i__ + j * x_dim1], abs(d__2));
#line 325 "dpbrfs.f"
/* L40: */
#line 325 "dpbrfs.f"
		}
#line 326 "dpbrfs.f"
		work[k] = work[k] + (d__1 = ab[*kd + 1 + k * ab_dim1], abs(
			d__1)) * xk + s;
#line 327 "dpbrfs.f"
/* L50: */
#line 327 "dpbrfs.f"
	    }
#line 328 "dpbrfs.f"
	} else {
#line 329 "dpbrfs.f"
	    i__2 = *n;
#line 329 "dpbrfs.f"
	    for (k = 1; k <= i__2; ++k) {
#line 330 "dpbrfs.f"
		s = 0.;
#line 331 "dpbrfs.f"
		xk = (d__1 = x[k + j * x_dim1], abs(d__1));
#line 332 "dpbrfs.f"
		work[k] += (d__1 = ab[k * ab_dim1 + 1], abs(d__1)) * xk;
#line 333 "dpbrfs.f"
		l = 1 - k;
/* Computing MIN */
#line 334 "dpbrfs.f"
		i__3 = *n, i__4 = k + *kd;
#line 334 "dpbrfs.f"
		i__5 = min(i__3,i__4);
#line 334 "dpbrfs.f"
		for (i__ = k + 1; i__ <= i__5; ++i__) {
#line 335 "dpbrfs.f"
		    work[i__] += (d__1 = ab[l + i__ + k * ab_dim1], abs(d__1))
			     * xk;
#line 336 "dpbrfs.f"
		    s += (d__1 = ab[l + i__ + k * ab_dim1], abs(d__1)) * (
			    d__2 = x[i__ + j * x_dim1], abs(d__2));
#line 337 "dpbrfs.f"
/* L60: */
#line 337 "dpbrfs.f"
		}
#line 338 "dpbrfs.f"
		work[k] += s;
#line 339 "dpbrfs.f"
/* L70: */
#line 339 "dpbrfs.f"
	    }
#line 340 "dpbrfs.f"
	}
#line 341 "dpbrfs.f"
	s = 0.;
#line 342 "dpbrfs.f"
	i__2 = *n;
#line 342 "dpbrfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 343 "dpbrfs.f"
	    if (work[i__] > safe2) {
/* Computing MAX */
#line 344 "dpbrfs.f"
		d__2 = s, d__3 = (d__1 = work[*n + i__], abs(d__1)) / work[
			i__];
#line 344 "dpbrfs.f"
		s = max(d__2,d__3);
#line 345 "dpbrfs.f"
	    } else {
/* Computing MAX */
#line 346 "dpbrfs.f"
		d__2 = s, d__3 = ((d__1 = work[*n + i__], abs(d__1)) + safe1) 
			/ (work[i__] + safe1);
#line 346 "dpbrfs.f"
		s = max(d__2,d__3);
#line 348 "dpbrfs.f"
	    }
#line 349 "dpbrfs.f"
/* L80: */
#line 349 "dpbrfs.f"
	}
#line 350 "dpbrfs.f"
	berr[j] = s;

/*        Test stopping criterion. Continue iterating if */
/*           1) The residual BERR(J) is larger than machine epsilon, and */
/*           2) BERR(J) decreased by at least a factor of 2 during the */
/*              last iteration, and */
/*           3) At most ITMAX iterations tried. */

#line 358 "dpbrfs.f"
	if (berr[j] > eps && berr[j] * 2. <= lstres && count <= 5) {

/*           Update solution and try again. */

#line 363 "dpbrfs.f"
	    dpbtrs_(uplo, n, kd, &c__1, &afb[afb_offset], ldafb, &work[*n + 1]
		    , n, info, (ftnlen)1);
#line 365 "dpbrfs.f"
	    daxpy_(n, &c_b14, &work[*n + 1], &c__1, &x[j * x_dim1 + 1], &c__1)
		    ;
#line 366 "dpbrfs.f"
	    lstres = berr[j];
#line 367 "dpbrfs.f"
	    ++count;
#line 368 "dpbrfs.f"
	    goto L20;
#line 369 "dpbrfs.f"
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

#line 393 "dpbrfs.f"
	i__2 = *n;
#line 393 "dpbrfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 394 "dpbrfs.f"
	    if (work[i__] > safe2) {
#line 395 "dpbrfs.f"
		work[i__] = (d__1 = work[*n + i__], abs(d__1)) + nz * eps * 
			work[i__];
#line 396 "dpbrfs.f"
	    } else {
#line 397 "dpbrfs.f"
		work[i__] = (d__1 = work[*n + i__], abs(d__1)) + nz * eps * 
			work[i__] + safe1;
#line 398 "dpbrfs.f"
	    }
#line 399 "dpbrfs.f"
/* L90: */
#line 399 "dpbrfs.f"
	}

#line 401 "dpbrfs.f"
	kase = 0;
#line 402 "dpbrfs.f"
L100:
#line 403 "dpbrfs.f"
	dlacn2_(n, &work[(*n << 1) + 1], &work[*n + 1], &iwork[1], &ferr[j], &
		kase, isave);
#line 405 "dpbrfs.f"
	if (kase != 0) {
#line 406 "dpbrfs.f"
	    if (kase == 1) {

/*              Multiply by diag(W)*inv(A**T). */

#line 410 "dpbrfs.f"
		dpbtrs_(uplo, n, kd, &c__1, &afb[afb_offset], ldafb, &work[*n 
			+ 1], n, info, (ftnlen)1);
#line 412 "dpbrfs.f"
		i__2 = *n;
#line 412 "dpbrfs.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 413 "dpbrfs.f"
		    work[*n + i__] *= work[i__];
#line 414 "dpbrfs.f"
/* L110: */
#line 414 "dpbrfs.f"
		}
#line 415 "dpbrfs.f"
	    } else if (kase == 2) {

/*              Multiply by inv(A)*diag(W). */

#line 419 "dpbrfs.f"
		i__2 = *n;
#line 419 "dpbrfs.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 420 "dpbrfs.f"
		    work[*n + i__] *= work[i__];
#line 421 "dpbrfs.f"
/* L120: */
#line 421 "dpbrfs.f"
		}
#line 422 "dpbrfs.f"
		dpbtrs_(uplo, n, kd, &c__1, &afb[afb_offset], ldafb, &work[*n 
			+ 1], n, info, (ftnlen)1);
#line 424 "dpbrfs.f"
	    }
#line 425 "dpbrfs.f"
	    goto L100;
#line 426 "dpbrfs.f"
	}

/*        Normalize error. */

#line 430 "dpbrfs.f"
	lstres = 0.;
#line 431 "dpbrfs.f"
	i__2 = *n;
#line 431 "dpbrfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 432 "dpbrfs.f"
	    d__2 = lstres, d__3 = (d__1 = x[i__ + j * x_dim1], abs(d__1));
#line 432 "dpbrfs.f"
	    lstres = max(d__2,d__3);
#line 433 "dpbrfs.f"
/* L130: */
#line 433 "dpbrfs.f"
	}
#line 434 "dpbrfs.f"
	if (lstres != 0.) {
#line 434 "dpbrfs.f"
	    ferr[j] /= lstres;
#line 434 "dpbrfs.f"
	}

#line 437 "dpbrfs.f"
/* L140: */
#line 437 "dpbrfs.f"
    }

#line 439 "dpbrfs.f"
    return 0;

/*     End of DPBRFS */

} /* dpbrfs_ */

