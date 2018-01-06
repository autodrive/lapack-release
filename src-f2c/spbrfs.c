#line 1 "spbrfs.f"
/* spbrfs.f -- translated by f2c (version 20100827).
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

#line 1 "spbrfs.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b12 = -1.;
static doublereal c_b14 = 1.;

/* > \brief \b SPBRFS */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SPBRFS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/spbrfs.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/spbrfs.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/spbrfs.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SPBRFS( UPLO, N, KD, NRHS, AB, LDAB, AFB, LDAFB, B, */
/*                          LDB, X, LDX, FERR, BERR, WORK, IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, KD, LDAB, LDAFB, LDB, LDX, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       REAL               AB( LDAB, * ), AFB( LDAFB, * ), B( LDB, * ), */
/*      $                   BERR( * ), FERR( * ), WORK( * ), X( LDX, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SPBRFS improves the computed solution to a system of linear */
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
/* >          AB is REAL array, dimension (LDAB,N) */
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
/* >          AFB is REAL array, dimension (LDAFB,N) */
/* >          The triangular factor U or L from the Cholesky factorization */
/* >          A = U**T*U or A = L*L**T of the band matrix A as computed by */
/* >          SPBTRF, in the same storage format as A (see AB). */
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
/* >          On entry, the solution matrix X, as computed by SPBTRS. */
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

/* > \date November 2011 */

/* > \ingroup realOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int spbrfs_(char *uplo, integer *n, integer *kd, integer *
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
    static integer isave[3], count;
    extern /* Subroutine */ int ssbmv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen);
    static logical upper;
    extern /* Subroutine */ int scopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), saxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *), slacn2_(
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, integer *);
    extern doublereal slamch_(char *, ftnlen);
    static doublereal safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal lstres;
    extern /* Subroutine */ int spbtrs_(char *, integer *, integer *, integer 
	    *, doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen);


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

#line 244 "spbrfs.f"
    /* Parameter adjustments */
#line 244 "spbrfs.f"
    ab_dim1 = *ldab;
#line 244 "spbrfs.f"
    ab_offset = 1 + ab_dim1;
#line 244 "spbrfs.f"
    ab -= ab_offset;
#line 244 "spbrfs.f"
    afb_dim1 = *ldafb;
#line 244 "spbrfs.f"
    afb_offset = 1 + afb_dim1;
#line 244 "spbrfs.f"
    afb -= afb_offset;
#line 244 "spbrfs.f"
    b_dim1 = *ldb;
#line 244 "spbrfs.f"
    b_offset = 1 + b_dim1;
#line 244 "spbrfs.f"
    b -= b_offset;
#line 244 "spbrfs.f"
    x_dim1 = *ldx;
#line 244 "spbrfs.f"
    x_offset = 1 + x_dim1;
#line 244 "spbrfs.f"
    x -= x_offset;
#line 244 "spbrfs.f"
    --ferr;
#line 244 "spbrfs.f"
    --berr;
#line 244 "spbrfs.f"
    --work;
#line 244 "spbrfs.f"
    --iwork;
#line 244 "spbrfs.f"

#line 244 "spbrfs.f"
    /* Function Body */
#line 244 "spbrfs.f"
    *info = 0;
#line 245 "spbrfs.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 246 "spbrfs.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 247 "spbrfs.f"
	*info = -1;
#line 248 "spbrfs.f"
    } else if (*n < 0) {
#line 249 "spbrfs.f"
	*info = -2;
#line 250 "spbrfs.f"
    } else if (*kd < 0) {
#line 251 "spbrfs.f"
	*info = -3;
#line 252 "spbrfs.f"
    } else if (*nrhs < 0) {
#line 253 "spbrfs.f"
	*info = -4;
#line 254 "spbrfs.f"
    } else if (*ldab < *kd + 1) {
#line 255 "spbrfs.f"
	*info = -6;
#line 256 "spbrfs.f"
    } else if (*ldafb < *kd + 1) {
#line 257 "spbrfs.f"
	*info = -8;
#line 258 "spbrfs.f"
    } else if (*ldb < max(1,*n)) {
#line 259 "spbrfs.f"
	*info = -10;
#line 260 "spbrfs.f"
    } else if (*ldx < max(1,*n)) {
#line 261 "spbrfs.f"
	*info = -12;
#line 262 "spbrfs.f"
    }
#line 263 "spbrfs.f"
    if (*info != 0) {
#line 264 "spbrfs.f"
	i__1 = -(*info);
#line 264 "spbrfs.f"
	xerbla_("SPBRFS", &i__1, (ftnlen)6);
#line 265 "spbrfs.f"
	return 0;
#line 266 "spbrfs.f"
    }

/*     Quick return if possible */

#line 270 "spbrfs.f"
    if (*n == 0 || *nrhs == 0) {
#line 271 "spbrfs.f"
	i__1 = *nrhs;
#line 271 "spbrfs.f"
	for (j = 1; j <= i__1; ++j) {
#line 272 "spbrfs.f"
	    ferr[j] = 0.;
#line 273 "spbrfs.f"
	    berr[j] = 0.;
#line 274 "spbrfs.f"
/* L10: */
#line 274 "spbrfs.f"
	}
#line 275 "spbrfs.f"
	return 0;
#line 276 "spbrfs.f"
    }

/*     NZ = maximum number of nonzero elements in each row of A, plus 1 */

/* Computing MIN */
#line 280 "spbrfs.f"
    i__1 = *n + 1, i__2 = (*kd << 1) + 2;
#line 280 "spbrfs.f"
    nz = min(i__1,i__2);
#line 281 "spbrfs.f"
    eps = slamch_("Epsilon", (ftnlen)7);
#line 282 "spbrfs.f"
    safmin = slamch_("Safe minimum", (ftnlen)12);
#line 283 "spbrfs.f"
    safe1 = nz * safmin;
#line 284 "spbrfs.f"
    safe2 = safe1 / eps;

/*     Do for each right hand side */

#line 288 "spbrfs.f"
    i__1 = *nrhs;
#line 288 "spbrfs.f"
    for (j = 1; j <= i__1; ++j) {

#line 290 "spbrfs.f"
	count = 1;
#line 291 "spbrfs.f"
	lstres = 3.;
#line 292 "spbrfs.f"
L20:

/*        Loop until stopping criterion is satisfied. */

/*        Compute residual R = B - A * X */

#line 298 "spbrfs.f"
	scopy_(n, &b[j * b_dim1 + 1], &c__1, &work[*n + 1], &c__1);
#line 299 "spbrfs.f"
	ssbmv_(uplo, n, kd, &c_b12, &ab[ab_offset], ldab, &x[j * x_dim1 + 1], 
		&c__1, &c_b14, &work[*n + 1], &c__1, (ftnlen)1);

/*        Compute componentwise relative backward error from formula */

/*        max(i) ( abs(R(i)) / ( abs(A)*abs(X) + abs(B) )(i) ) */

/*        where abs(Z) is the componentwise absolute value of the matrix */
/*        or vector Z.  If the i-th component of the denominator is less */
/*        than SAFE2, then SAFE1 is added to the i-th components of the */
/*        numerator and denominator before dividing. */

#line 311 "spbrfs.f"
	i__2 = *n;
#line 311 "spbrfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 312 "spbrfs.f"
	    work[i__] = (d__1 = b[i__ + j * b_dim1], abs(d__1));
#line 313 "spbrfs.f"
/* L30: */
#line 313 "spbrfs.f"
	}

/*        Compute abs(A)*abs(X) + abs(B). */

#line 317 "spbrfs.f"
	if (upper) {
#line 318 "spbrfs.f"
	    i__2 = *n;
#line 318 "spbrfs.f"
	    for (k = 1; k <= i__2; ++k) {
#line 319 "spbrfs.f"
		s = 0.;
#line 320 "spbrfs.f"
		xk = (d__1 = x[k + j * x_dim1], abs(d__1));
#line 321 "spbrfs.f"
		l = *kd + 1 - k;
/* Computing MAX */
#line 322 "spbrfs.f"
		i__3 = 1, i__4 = k - *kd;
#line 322 "spbrfs.f"
		i__5 = k - 1;
#line 322 "spbrfs.f"
		for (i__ = max(i__3,i__4); i__ <= i__5; ++i__) {
#line 323 "spbrfs.f"
		    work[i__] += (d__1 = ab[l + i__ + k * ab_dim1], abs(d__1))
			     * xk;
#line 324 "spbrfs.f"
		    s += (d__1 = ab[l + i__ + k * ab_dim1], abs(d__1)) * (
			    d__2 = x[i__ + j * x_dim1], abs(d__2));
#line 325 "spbrfs.f"
/* L40: */
#line 325 "spbrfs.f"
		}
#line 326 "spbrfs.f"
		work[k] = work[k] + (d__1 = ab[*kd + 1 + k * ab_dim1], abs(
			d__1)) * xk + s;
#line 327 "spbrfs.f"
/* L50: */
#line 327 "spbrfs.f"
	    }
#line 328 "spbrfs.f"
	} else {
#line 329 "spbrfs.f"
	    i__2 = *n;
#line 329 "spbrfs.f"
	    for (k = 1; k <= i__2; ++k) {
#line 330 "spbrfs.f"
		s = 0.;
#line 331 "spbrfs.f"
		xk = (d__1 = x[k + j * x_dim1], abs(d__1));
#line 332 "spbrfs.f"
		work[k] += (d__1 = ab[k * ab_dim1 + 1], abs(d__1)) * xk;
#line 333 "spbrfs.f"
		l = 1 - k;
/* Computing MIN */
#line 334 "spbrfs.f"
		i__3 = *n, i__4 = k + *kd;
#line 334 "spbrfs.f"
		i__5 = min(i__3,i__4);
#line 334 "spbrfs.f"
		for (i__ = k + 1; i__ <= i__5; ++i__) {
#line 335 "spbrfs.f"
		    work[i__] += (d__1 = ab[l + i__ + k * ab_dim1], abs(d__1))
			     * xk;
#line 336 "spbrfs.f"
		    s += (d__1 = ab[l + i__ + k * ab_dim1], abs(d__1)) * (
			    d__2 = x[i__ + j * x_dim1], abs(d__2));
#line 337 "spbrfs.f"
/* L60: */
#line 337 "spbrfs.f"
		}
#line 338 "spbrfs.f"
		work[k] += s;
#line 339 "spbrfs.f"
/* L70: */
#line 339 "spbrfs.f"
	    }
#line 340 "spbrfs.f"
	}
#line 341 "spbrfs.f"
	s = 0.;
#line 342 "spbrfs.f"
	i__2 = *n;
#line 342 "spbrfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 343 "spbrfs.f"
	    if (work[i__] > safe2) {
/* Computing MAX */
#line 344 "spbrfs.f"
		d__2 = s, d__3 = (d__1 = work[*n + i__], abs(d__1)) / work[
			i__];
#line 344 "spbrfs.f"
		s = max(d__2,d__3);
#line 345 "spbrfs.f"
	    } else {
/* Computing MAX */
#line 346 "spbrfs.f"
		d__2 = s, d__3 = ((d__1 = work[*n + i__], abs(d__1)) + safe1) 
			/ (work[i__] + safe1);
#line 346 "spbrfs.f"
		s = max(d__2,d__3);
#line 348 "spbrfs.f"
	    }
#line 349 "spbrfs.f"
/* L80: */
#line 349 "spbrfs.f"
	}
#line 350 "spbrfs.f"
	berr[j] = s;

/*        Test stopping criterion. Continue iterating if */
/*           1) The residual BERR(J) is larger than machine epsilon, and */
/*           2) BERR(J) decreased by at least a factor of 2 during the */
/*              last iteration, and */
/*           3) At most ITMAX iterations tried. */

#line 358 "spbrfs.f"
	if (berr[j] > eps && berr[j] * 2. <= lstres && count <= 5) {

/*           Update solution and try again. */

#line 363 "spbrfs.f"
	    spbtrs_(uplo, n, kd, &c__1, &afb[afb_offset], ldafb, &work[*n + 1]
		    , n, info, (ftnlen)1);
#line 365 "spbrfs.f"
	    saxpy_(n, &c_b14, &work[*n + 1], &c__1, &x[j * x_dim1 + 1], &c__1)
		    ;
#line 366 "spbrfs.f"
	    lstres = berr[j];
#line 367 "spbrfs.f"
	    ++count;
#line 368 "spbrfs.f"
	    goto L20;
#line 369 "spbrfs.f"
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

#line 393 "spbrfs.f"
	i__2 = *n;
#line 393 "spbrfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 394 "spbrfs.f"
	    if (work[i__] > safe2) {
#line 395 "spbrfs.f"
		work[i__] = (d__1 = work[*n + i__], abs(d__1)) + nz * eps * 
			work[i__];
#line 396 "spbrfs.f"
	    } else {
#line 397 "spbrfs.f"
		work[i__] = (d__1 = work[*n + i__], abs(d__1)) + nz * eps * 
			work[i__] + safe1;
#line 398 "spbrfs.f"
	    }
#line 399 "spbrfs.f"
/* L90: */
#line 399 "spbrfs.f"
	}

#line 401 "spbrfs.f"
	kase = 0;
#line 402 "spbrfs.f"
L100:
#line 403 "spbrfs.f"
	slacn2_(n, &work[(*n << 1) + 1], &work[*n + 1], &iwork[1], &ferr[j], &
		kase, isave);
#line 405 "spbrfs.f"
	if (kase != 0) {
#line 406 "spbrfs.f"
	    if (kase == 1) {

/*              Multiply by diag(W)*inv(A**T). */

#line 410 "spbrfs.f"
		spbtrs_(uplo, n, kd, &c__1, &afb[afb_offset], ldafb, &work[*n 
			+ 1], n, info, (ftnlen)1);
#line 412 "spbrfs.f"
		i__2 = *n;
#line 412 "spbrfs.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 413 "spbrfs.f"
		    work[*n + i__] *= work[i__];
#line 414 "spbrfs.f"
/* L110: */
#line 414 "spbrfs.f"
		}
#line 415 "spbrfs.f"
	    } else if (kase == 2) {

/*              Multiply by inv(A)*diag(W). */

#line 419 "spbrfs.f"
		i__2 = *n;
#line 419 "spbrfs.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 420 "spbrfs.f"
		    work[*n + i__] *= work[i__];
#line 421 "spbrfs.f"
/* L120: */
#line 421 "spbrfs.f"
		}
#line 422 "spbrfs.f"
		spbtrs_(uplo, n, kd, &c__1, &afb[afb_offset], ldafb, &work[*n 
			+ 1], n, info, (ftnlen)1);
#line 424 "spbrfs.f"
	    }
#line 425 "spbrfs.f"
	    goto L100;
#line 426 "spbrfs.f"
	}

/*        Normalize error. */

#line 430 "spbrfs.f"
	lstres = 0.;
#line 431 "spbrfs.f"
	i__2 = *n;
#line 431 "spbrfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 432 "spbrfs.f"
	    d__2 = lstres, d__3 = (d__1 = x[i__ + j * x_dim1], abs(d__1));
#line 432 "spbrfs.f"
	    lstres = max(d__2,d__3);
#line 433 "spbrfs.f"
/* L130: */
#line 433 "spbrfs.f"
	}
#line 434 "spbrfs.f"
	if (lstres != 0.) {
#line 434 "spbrfs.f"
	    ferr[j] /= lstres;
#line 434 "spbrfs.f"
	}

#line 437 "spbrfs.f"
/* L140: */
#line 437 "spbrfs.f"
    }

#line 439 "spbrfs.f"
    return 0;

/*     End of SPBRFS */

} /* spbrfs_ */

