#line 1 "dtbrfs.f"
/* dtbrfs.f -- translated by f2c (version 20100827).
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

#line 1 "dtbrfs.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b19 = -1.;

/* > \brief \b DTBRFS */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DTBRFS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dtbrfs.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dtbrfs.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dtbrfs.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DTBRFS( UPLO, TRANS, DIAG, N, KD, NRHS, AB, LDAB, B, */
/*                          LDB, X, LDX, FERR, BERR, WORK, IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          DIAG, TRANS, UPLO */
/*       INTEGER            INFO, KD, LDAB, LDB, LDX, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       DOUBLE PRECISION   AB( LDAB, * ), B( LDB, * ), BERR( * ), */
/*      $                   FERR( * ), WORK( * ), X( LDX, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DTBRFS provides error bounds and backward error estimates for the */
/* > solution to a system of linear equations with a triangular band */
/* > coefficient matrix. */
/* > */
/* > The solution matrix X must be computed by DTBTRS or some other */
/* > means before entering this routine.  DTBRFS does not do iterative */
/* > refinement because doing so cannot improve the backward error. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          = 'U':  A is upper triangular; */
/* >          = 'L':  A is lower triangular. */
/* > \endverbatim */
/* > */
/* > \param[in] TRANS */
/* > \verbatim */
/* >          TRANS is CHARACTER*1 */
/* >          Specifies the form of the system of equations: */
/* >          = 'N':  A * X = B  (No transpose) */
/* >          = 'T':  A**T * X = B  (Transpose) */
/* >          = 'C':  A**H * X = B  (Conjugate transpose = Transpose) */
/* > \endverbatim */
/* > */
/* > \param[in] DIAG */
/* > \verbatim */
/* >          DIAG is CHARACTER*1 */
/* >          = 'N':  A is non-unit triangular; */
/* >          = 'U':  A is unit triangular. */
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
/* >          The number of superdiagonals or subdiagonals of the */
/* >          triangular band matrix A.  KD >= 0. */
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
/* >          The upper or lower triangular band matrix A, stored in the */
/* >          first kd+1 rows of the array. The j-th column of A is stored */
/* >          in the j-th column of the array AB as follows: */
/* >          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j; */
/* >          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd). */
/* >          If DIAG = 'U', the diagonal elements of A are not referenced */
/* >          and are assumed to be 1. */
/* > \endverbatim */
/* > */
/* > \param[in] LDAB */
/* > \verbatim */
/* >          LDAB is INTEGER */
/* >          The leading dimension of the array AB.  LDAB >= KD+1. */
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
/* > \param[in] X */
/* > \verbatim */
/* >          X is DOUBLE PRECISION array, dimension (LDX,NRHS) */
/* >          The solution matrix X. */
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

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup doubleOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int dtbrfs_(char *uplo, char *trans, char *diag, integer *n, 
	integer *kd, integer *nrhs, doublereal *ab, integer *ldab, doublereal 
	*b, integer *ldb, doublereal *x, integer *ldx, doublereal *ferr, 
	doublereal *berr, doublereal *work, integer *iwork, integer *info, 
	ftnlen uplo_len, ftnlen trans_len, ftnlen diag_len)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, b_dim1, b_offset, x_dim1, x_offset, i__1, 
	    i__2, i__3, i__4, i__5;
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
    extern /* Subroutine */ int dtbmv_(char *, char *, char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen, ftnlen, ftnlen), dcopy_(integer *, doublereal *, integer *
	    , doublereal *, integer *), dtbsv_(char *, char *, char *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, ftnlen, ftnlen, ftnlen), daxpy_(integer *, doublereal *
	    , doublereal *, integer *, doublereal *, integer *);
    static logical upper;
    extern /* Subroutine */ int dlacn2_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, integer *);
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static logical notran;
    static char transt[1];
    static logical nounit;
    static doublereal lstres;


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

#line 238 "dtbrfs.f"
    /* Parameter adjustments */
#line 238 "dtbrfs.f"
    ab_dim1 = *ldab;
#line 238 "dtbrfs.f"
    ab_offset = 1 + ab_dim1;
#line 238 "dtbrfs.f"
    ab -= ab_offset;
#line 238 "dtbrfs.f"
    b_dim1 = *ldb;
#line 238 "dtbrfs.f"
    b_offset = 1 + b_dim1;
#line 238 "dtbrfs.f"
    b -= b_offset;
#line 238 "dtbrfs.f"
    x_dim1 = *ldx;
#line 238 "dtbrfs.f"
    x_offset = 1 + x_dim1;
#line 238 "dtbrfs.f"
    x -= x_offset;
#line 238 "dtbrfs.f"
    --ferr;
#line 238 "dtbrfs.f"
    --berr;
#line 238 "dtbrfs.f"
    --work;
#line 238 "dtbrfs.f"
    --iwork;
#line 238 "dtbrfs.f"

#line 238 "dtbrfs.f"
    /* Function Body */
#line 238 "dtbrfs.f"
    *info = 0;
#line 239 "dtbrfs.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 240 "dtbrfs.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 241 "dtbrfs.f"
    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);

#line 243 "dtbrfs.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 244 "dtbrfs.f"
	*info = -1;
#line 245 "dtbrfs.f"
    } else if (! notran && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1) && ! 
	    lsame_(trans, "C", (ftnlen)1, (ftnlen)1)) {
#line 247 "dtbrfs.f"
	*info = -2;
#line 248 "dtbrfs.f"
    } else if (! nounit && ! lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 249 "dtbrfs.f"
	*info = -3;
#line 250 "dtbrfs.f"
    } else if (*n < 0) {
#line 251 "dtbrfs.f"
	*info = -4;
#line 252 "dtbrfs.f"
    } else if (*kd < 0) {
#line 253 "dtbrfs.f"
	*info = -5;
#line 254 "dtbrfs.f"
    } else if (*nrhs < 0) {
#line 255 "dtbrfs.f"
	*info = -6;
#line 256 "dtbrfs.f"
    } else if (*ldab < *kd + 1) {
#line 257 "dtbrfs.f"
	*info = -8;
#line 258 "dtbrfs.f"
    } else if (*ldb < max(1,*n)) {
#line 259 "dtbrfs.f"
	*info = -10;
#line 260 "dtbrfs.f"
    } else if (*ldx < max(1,*n)) {
#line 261 "dtbrfs.f"
	*info = -12;
#line 262 "dtbrfs.f"
    }
#line 263 "dtbrfs.f"
    if (*info != 0) {
#line 264 "dtbrfs.f"
	i__1 = -(*info);
#line 264 "dtbrfs.f"
	xerbla_("DTBRFS", &i__1, (ftnlen)6);
#line 265 "dtbrfs.f"
	return 0;
#line 266 "dtbrfs.f"
    }

/*     Quick return if possible */

#line 270 "dtbrfs.f"
    if (*n == 0 || *nrhs == 0) {
#line 271 "dtbrfs.f"
	i__1 = *nrhs;
#line 271 "dtbrfs.f"
	for (j = 1; j <= i__1; ++j) {
#line 272 "dtbrfs.f"
	    ferr[j] = 0.;
#line 273 "dtbrfs.f"
	    berr[j] = 0.;
#line 274 "dtbrfs.f"
/* L10: */
#line 274 "dtbrfs.f"
	}
#line 275 "dtbrfs.f"
	return 0;
#line 276 "dtbrfs.f"
    }

#line 278 "dtbrfs.f"
    if (notran) {
#line 279 "dtbrfs.f"
	*(unsigned char *)transt = 'T';
#line 280 "dtbrfs.f"
    } else {
#line 281 "dtbrfs.f"
	*(unsigned char *)transt = 'N';
#line 282 "dtbrfs.f"
    }

/*     NZ = maximum number of nonzero elements in each row of A, plus 1 */

#line 286 "dtbrfs.f"
    nz = *kd + 2;
#line 287 "dtbrfs.f"
    eps = dlamch_("Epsilon", (ftnlen)7);
#line 288 "dtbrfs.f"
    safmin = dlamch_("Safe minimum", (ftnlen)12);
#line 289 "dtbrfs.f"
    safe1 = nz * safmin;
#line 290 "dtbrfs.f"
    safe2 = safe1 / eps;

/*     Do for each right hand side */

#line 294 "dtbrfs.f"
    i__1 = *nrhs;
#line 294 "dtbrfs.f"
    for (j = 1; j <= i__1; ++j) {

/*        Compute residual R = B - op(A) * X, */
/*        where op(A) = A or A**T, depending on TRANS. */

#line 299 "dtbrfs.f"
	dcopy_(n, &x[j * x_dim1 + 1], &c__1, &work[*n + 1], &c__1);
#line 300 "dtbrfs.f"
	dtbmv_(uplo, trans, diag, n, kd, &ab[ab_offset], ldab, &work[*n + 1], 
		&c__1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 302 "dtbrfs.f"
	daxpy_(n, &c_b19, &b[j * b_dim1 + 1], &c__1, &work[*n + 1], &c__1);

/*        Compute componentwise relative backward error from formula */

/*        max(i) ( abs(R(i)) / ( abs(op(A))*abs(X) + abs(B) )(i) ) */

/*        where abs(Z) is the componentwise absolute value of the matrix */
/*        or vector Z.  If the i-th component of the denominator is less */
/*        than SAFE2, then SAFE1 is added to the i-th components of the */
/*        numerator and denominator before dividing. */

#line 313 "dtbrfs.f"
	i__2 = *n;
#line 313 "dtbrfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 314 "dtbrfs.f"
	    work[i__] = (d__1 = b[i__ + j * b_dim1], abs(d__1));
#line 315 "dtbrfs.f"
/* L20: */
#line 315 "dtbrfs.f"
	}

#line 317 "dtbrfs.f"
	if (notran) {

/*           Compute abs(A)*abs(X) + abs(B). */

#line 321 "dtbrfs.f"
	    if (upper) {
#line 322 "dtbrfs.f"
		if (nounit) {
#line 323 "dtbrfs.f"
		    i__2 = *n;
#line 323 "dtbrfs.f"
		    for (k = 1; k <= i__2; ++k) {
#line 324 "dtbrfs.f"
			xk = (d__1 = x[k + j * x_dim1], abs(d__1));
/* Computing MAX */
#line 325 "dtbrfs.f"
			i__3 = 1, i__4 = k - *kd;
#line 325 "dtbrfs.f"
			i__5 = k;
#line 325 "dtbrfs.f"
			for (i__ = max(i__3,i__4); i__ <= i__5; ++i__) {
#line 326 "dtbrfs.f"
			    work[i__] += (d__1 = ab[*kd + 1 + i__ - k + k * 
				    ab_dim1], abs(d__1)) * xk;
#line 328 "dtbrfs.f"
/* L30: */
#line 328 "dtbrfs.f"
			}
#line 329 "dtbrfs.f"
/* L40: */
#line 329 "dtbrfs.f"
		    }
#line 330 "dtbrfs.f"
		} else {
#line 331 "dtbrfs.f"
		    i__2 = *n;
#line 331 "dtbrfs.f"
		    for (k = 1; k <= i__2; ++k) {
#line 332 "dtbrfs.f"
			xk = (d__1 = x[k + j * x_dim1], abs(d__1));
/* Computing MAX */
#line 333 "dtbrfs.f"
			i__5 = 1, i__3 = k - *kd;
#line 333 "dtbrfs.f"
			i__4 = k - 1;
#line 333 "dtbrfs.f"
			for (i__ = max(i__5,i__3); i__ <= i__4; ++i__) {
#line 334 "dtbrfs.f"
			    work[i__] += (d__1 = ab[*kd + 1 + i__ - k + k * 
				    ab_dim1], abs(d__1)) * xk;
#line 336 "dtbrfs.f"
/* L50: */
#line 336 "dtbrfs.f"
			}
#line 337 "dtbrfs.f"
			work[k] += xk;
#line 338 "dtbrfs.f"
/* L60: */
#line 338 "dtbrfs.f"
		    }
#line 339 "dtbrfs.f"
		}
#line 340 "dtbrfs.f"
	    } else {
#line 341 "dtbrfs.f"
		if (nounit) {
#line 342 "dtbrfs.f"
		    i__2 = *n;
#line 342 "dtbrfs.f"
		    for (k = 1; k <= i__2; ++k) {
#line 343 "dtbrfs.f"
			xk = (d__1 = x[k + j * x_dim1], abs(d__1));
/* Computing MIN */
#line 344 "dtbrfs.f"
			i__5 = *n, i__3 = k + *kd;
#line 344 "dtbrfs.f"
			i__4 = min(i__5,i__3);
#line 344 "dtbrfs.f"
			for (i__ = k; i__ <= i__4; ++i__) {
#line 345 "dtbrfs.f"
			    work[i__] += (d__1 = ab[i__ + 1 - k + k * ab_dim1]
				    , abs(d__1)) * xk;
#line 346 "dtbrfs.f"
/* L70: */
#line 346 "dtbrfs.f"
			}
#line 347 "dtbrfs.f"
/* L80: */
#line 347 "dtbrfs.f"
		    }
#line 348 "dtbrfs.f"
		} else {
#line 349 "dtbrfs.f"
		    i__2 = *n;
#line 349 "dtbrfs.f"
		    for (k = 1; k <= i__2; ++k) {
#line 350 "dtbrfs.f"
			xk = (d__1 = x[k + j * x_dim1], abs(d__1));
/* Computing MIN */
#line 351 "dtbrfs.f"
			i__5 = *n, i__3 = k + *kd;
#line 351 "dtbrfs.f"
			i__4 = min(i__5,i__3);
#line 351 "dtbrfs.f"
			for (i__ = k + 1; i__ <= i__4; ++i__) {
#line 352 "dtbrfs.f"
			    work[i__] += (d__1 = ab[i__ + 1 - k + k * ab_dim1]
				    , abs(d__1)) * xk;
#line 353 "dtbrfs.f"
/* L90: */
#line 353 "dtbrfs.f"
			}
#line 354 "dtbrfs.f"
			work[k] += xk;
#line 355 "dtbrfs.f"
/* L100: */
#line 355 "dtbrfs.f"
		    }
#line 356 "dtbrfs.f"
		}
#line 357 "dtbrfs.f"
	    }
#line 358 "dtbrfs.f"
	} else {

/*           Compute abs(A**T)*abs(X) + abs(B). */

#line 362 "dtbrfs.f"
	    if (upper) {
#line 363 "dtbrfs.f"
		if (nounit) {
#line 364 "dtbrfs.f"
		    i__2 = *n;
#line 364 "dtbrfs.f"
		    for (k = 1; k <= i__2; ++k) {
#line 365 "dtbrfs.f"
			s = 0.;
/* Computing MAX */
#line 366 "dtbrfs.f"
			i__4 = 1, i__5 = k - *kd;
#line 366 "dtbrfs.f"
			i__3 = k;
#line 366 "dtbrfs.f"
			for (i__ = max(i__4,i__5); i__ <= i__3; ++i__) {
#line 367 "dtbrfs.f"
			    s += (d__1 = ab[*kd + 1 + i__ - k + k * ab_dim1], 
				    abs(d__1)) * (d__2 = x[i__ + j * x_dim1], 
				    abs(d__2));
#line 369 "dtbrfs.f"
/* L110: */
#line 369 "dtbrfs.f"
			}
#line 370 "dtbrfs.f"
			work[k] += s;
#line 371 "dtbrfs.f"
/* L120: */
#line 371 "dtbrfs.f"
		    }
#line 372 "dtbrfs.f"
		} else {
#line 373 "dtbrfs.f"
		    i__2 = *n;
#line 373 "dtbrfs.f"
		    for (k = 1; k <= i__2; ++k) {
#line 374 "dtbrfs.f"
			s = (d__1 = x[k + j * x_dim1], abs(d__1));
/* Computing MAX */
#line 375 "dtbrfs.f"
			i__3 = 1, i__4 = k - *kd;
#line 375 "dtbrfs.f"
			i__5 = k - 1;
#line 375 "dtbrfs.f"
			for (i__ = max(i__3,i__4); i__ <= i__5; ++i__) {
#line 376 "dtbrfs.f"
			    s += (d__1 = ab[*kd + 1 + i__ - k + k * ab_dim1], 
				    abs(d__1)) * (d__2 = x[i__ + j * x_dim1], 
				    abs(d__2));
#line 378 "dtbrfs.f"
/* L130: */
#line 378 "dtbrfs.f"
			}
#line 379 "dtbrfs.f"
			work[k] += s;
#line 380 "dtbrfs.f"
/* L140: */
#line 380 "dtbrfs.f"
		    }
#line 381 "dtbrfs.f"
		}
#line 382 "dtbrfs.f"
	    } else {
#line 383 "dtbrfs.f"
		if (nounit) {
#line 384 "dtbrfs.f"
		    i__2 = *n;
#line 384 "dtbrfs.f"
		    for (k = 1; k <= i__2; ++k) {
#line 385 "dtbrfs.f"
			s = 0.;
/* Computing MIN */
#line 386 "dtbrfs.f"
			i__3 = *n, i__4 = k + *kd;
#line 386 "dtbrfs.f"
			i__5 = min(i__3,i__4);
#line 386 "dtbrfs.f"
			for (i__ = k; i__ <= i__5; ++i__) {
#line 387 "dtbrfs.f"
			    s += (d__1 = ab[i__ + 1 - k + k * ab_dim1], abs(
				    d__1)) * (d__2 = x[i__ + j * x_dim1], abs(
				    d__2));
#line 388 "dtbrfs.f"
/* L150: */
#line 388 "dtbrfs.f"
			}
#line 389 "dtbrfs.f"
			work[k] += s;
#line 390 "dtbrfs.f"
/* L160: */
#line 390 "dtbrfs.f"
		    }
#line 391 "dtbrfs.f"
		} else {
#line 392 "dtbrfs.f"
		    i__2 = *n;
#line 392 "dtbrfs.f"
		    for (k = 1; k <= i__2; ++k) {
#line 393 "dtbrfs.f"
			s = (d__1 = x[k + j * x_dim1], abs(d__1));
/* Computing MIN */
#line 394 "dtbrfs.f"
			i__3 = *n, i__4 = k + *kd;
#line 394 "dtbrfs.f"
			i__5 = min(i__3,i__4);
#line 394 "dtbrfs.f"
			for (i__ = k + 1; i__ <= i__5; ++i__) {
#line 395 "dtbrfs.f"
			    s += (d__1 = ab[i__ + 1 - k + k * ab_dim1], abs(
				    d__1)) * (d__2 = x[i__ + j * x_dim1], abs(
				    d__2));
#line 396 "dtbrfs.f"
/* L170: */
#line 396 "dtbrfs.f"
			}
#line 397 "dtbrfs.f"
			work[k] += s;
#line 398 "dtbrfs.f"
/* L180: */
#line 398 "dtbrfs.f"
		    }
#line 399 "dtbrfs.f"
		}
#line 400 "dtbrfs.f"
	    }
#line 401 "dtbrfs.f"
	}
#line 402 "dtbrfs.f"
	s = 0.;
#line 403 "dtbrfs.f"
	i__2 = *n;
#line 403 "dtbrfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 404 "dtbrfs.f"
	    if (work[i__] > safe2) {
/* Computing MAX */
#line 405 "dtbrfs.f"
		d__2 = s, d__3 = (d__1 = work[*n + i__], abs(d__1)) / work[
			i__];
#line 405 "dtbrfs.f"
		s = max(d__2,d__3);
#line 406 "dtbrfs.f"
	    } else {
/* Computing MAX */
#line 407 "dtbrfs.f"
		d__2 = s, d__3 = ((d__1 = work[*n + i__], abs(d__1)) + safe1) 
			/ (work[i__] + safe1);
#line 407 "dtbrfs.f"
		s = max(d__2,d__3);
#line 409 "dtbrfs.f"
	    }
#line 410 "dtbrfs.f"
/* L190: */
#line 410 "dtbrfs.f"
	}
#line 411 "dtbrfs.f"
	berr[j] = s;

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

/*        Use DLACN2 to estimate the infinity-norm of the matrix */
/*           inv(op(A)) * diag(W), */
/*        where W = abs(R) + NZ*EPS*( abs(op(A))*abs(X)+abs(B) ))) */

#line 435 "dtbrfs.f"
	i__2 = *n;
#line 435 "dtbrfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 436 "dtbrfs.f"
	    if (work[i__] > safe2) {
#line 437 "dtbrfs.f"
		work[i__] = (d__1 = work[*n + i__], abs(d__1)) + nz * eps * 
			work[i__];
#line 438 "dtbrfs.f"
	    } else {
#line 439 "dtbrfs.f"
		work[i__] = (d__1 = work[*n + i__], abs(d__1)) + nz * eps * 
			work[i__] + safe1;
#line 440 "dtbrfs.f"
	    }
#line 441 "dtbrfs.f"
/* L200: */
#line 441 "dtbrfs.f"
	}

#line 443 "dtbrfs.f"
	kase = 0;
#line 444 "dtbrfs.f"
L210:
#line 445 "dtbrfs.f"
	dlacn2_(n, &work[(*n << 1) + 1], &work[*n + 1], &iwork[1], &ferr[j], &
		kase, isave);
#line 447 "dtbrfs.f"
	if (kase != 0) {
#line 448 "dtbrfs.f"
	    if (kase == 1) {

/*              Multiply by diag(W)*inv(op(A)**T). */

#line 452 "dtbrfs.f"
		dtbsv_(uplo, transt, diag, n, kd, &ab[ab_offset], ldab, &work[
			*n + 1], &c__1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 454 "dtbrfs.f"
		i__2 = *n;
#line 454 "dtbrfs.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 455 "dtbrfs.f"
		    work[*n + i__] = work[i__] * work[*n + i__];
#line 456 "dtbrfs.f"
/* L220: */
#line 456 "dtbrfs.f"
		}
#line 457 "dtbrfs.f"
	    } else {

/*              Multiply by inv(op(A))*diag(W). */

#line 461 "dtbrfs.f"
		i__2 = *n;
#line 461 "dtbrfs.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 462 "dtbrfs.f"
		    work[*n + i__] = work[i__] * work[*n + i__];
#line 463 "dtbrfs.f"
/* L230: */
#line 463 "dtbrfs.f"
		}
#line 464 "dtbrfs.f"
		dtbsv_(uplo, trans, diag, n, kd, &ab[ab_offset], ldab, &work[*
			n + 1], &c__1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 466 "dtbrfs.f"
	    }
#line 467 "dtbrfs.f"
	    goto L210;
#line 468 "dtbrfs.f"
	}

/*        Normalize error. */

#line 472 "dtbrfs.f"
	lstres = 0.;
#line 473 "dtbrfs.f"
	i__2 = *n;
#line 473 "dtbrfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 474 "dtbrfs.f"
	    d__2 = lstres, d__3 = (d__1 = x[i__ + j * x_dim1], abs(d__1));
#line 474 "dtbrfs.f"
	    lstres = max(d__2,d__3);
#line 475 "dtbrfs.f"
/* L240: */
#line 475 "dtbrfs.f"
	}
#line 476 "dtbrfs.f"
	if (lstres != 0.) {
#line 476 "dtbrfs.f"
	    ferr[j] /= lstres;
#line 476 "dtbrfs.f"
	}

#line 479 "dtbrfs.f"
/* L250: */
#line 479 "dtbrfs.f"
    }

#line 481 "dtbrfs.f"
    return 0;

/*     End of DTBRFS */

} /* dtbrfs_ */

