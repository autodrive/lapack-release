#line 1 "dtrrfs.f"
/* dtrrfs.f -- translated by f2c (version 20100827).
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

#line 1 "dtrrfs.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b19 = -1.;

/* > \brief \b DTRRFS */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DTRRFS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dtrrfs.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dtrrfs.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dtrrfs.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DTRRFS( UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDB, X, */
/*                          LDX, FERR, BERR, WORK, IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          DIAG, TRANS, UPLO */
/*       INTEGER            INFO, LDA, LDB, LDX, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), BERR( * ), FERR( * ), */
/*      $                   WORK( * ), X( LDX, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DTRRFS provides error bounds and backward error estimates for the */
/* > solution to a system of linear equations with a triangular */
/* > coefficient matrix. */
/* > */
/* > The solution matrix X must be computed by DTRTRS or some other */
/* > means before entering this routine.  DTRRFS does not do iterative */
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
/* >          The triangular matrix A.  If UPLO = 'U', the leading N-by-N */
/* >          upper triangular part of the array A contains the upper */
/* >          triangular matrix, and the strictly lower triangular part of */
/* >          A is not referenced.  If UPLO = 'L', the leading N-by-N lower */
/* >          triangular part of the array A contains the lower triangular */
/* >          matrix, and the strictly upper triangular part of A is not */
/* >          referenced.  If DIAG = 'U', the diagonal elements of A are */
/* >          also not referenced and are assumed to be 1. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
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
/* Subroutine */ int dtrrfs_(char *uplo, char *trans, char *diag, integer *n, 
	integer *nrhs, doublereal *a, integer *lda, doublereal *b, integer *
	ldb, doublereal *x, integer *ldx, doublereal *ferr, doublereal *berr, 
	doublereal *work, integer *iwork, integer *info, ftnlen uplo_len, 
	ftnlen trans_len, ftnlen diag_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, x_dim1, x_offset, i__1, i__2, 
	    i__3;
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
    static logical upper;
    extern /* Subroutine */ int dtrmv_(char *, char *, char *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen, ftnlen, 
	    ftnlen), dtrsv_(char *, char *, char *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, ftnlen, ftnlen, ftnlen), 
	    dlacn2_(integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, integer *);
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

#line 232 "dtrrfs.f"
    /* Parameter adjustments */
#line 232 "dtrrfs.f"
    a_dim1 = *lda;
#line 232 "dtrrfs.f"
    a_offset = 1 + a_dim1;
#line 232 "dtrrfs.f"
    a -= a_offset;
#line 232 "dtrrfs.f"
    b_dim1 = *ldb;
#line 232 "dtrrfs.f"
    b_offset = 1 + b_dim1;
#line 232 "dtrrfs.f"
    b -= b_offset;
#line 232 "dtrrfs.f"
    x_dim1 = *ldx;
#line 232 "dtrrfs.f"
    x_offset = 1 + x_dim1;
#line 232 "dtrrfs.f"
    x -= x_offset;
#line 232 "dtrrfs.f"
    --ferr;
#line 232 "dtrrfs.f"
    --berr;
#line 232 "dtrrfs.f"
    --work;
#line 232 "dtrrfs.f"
    --iwork;
#line 232 "dtrrfs.f"

#line 232 "dtrrfs.f"
    /* Function Body */
#line 232 "dtrrfs.f"
    *info = 0;
#line 233 "dtrrfs.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 234 "dtrrfs.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 235 "dtrrfs.f"
    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);

#line 237 "dtrrfs.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 238 "dtrrfs.f"
	*info = -1;
#line 239 "dtrrfs.f"
    } else if (! notran && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1) && ! 
	    lsame_(trans, "C", (ftnlen)1, (ftnlen)1)) {
#line 241 "dtrrfs.f"
	*info = -2;
#line 242 "dtrrfs.f"
    } else if (! nounit && ! lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 243 "dtrrfs.f"
	*info = -3;
#line 244 "dtrrfs.f"
    } else if (*n < 0) {
#line 245 "dtrrfs.f"
	*info = -4;
#line 246 "dtrrfs.f"
    } else if (*nrhs < 0) {
#line 247 "dtrrfs.f"
	*info = -5;
#line 248 "dtrrfs.f"
    } else if (*lda < max(1,*n)) {
#line 249 "dtrrfs.f"
	*info = -7;
#line 250 "dtrrfs.f"
    } else if (*ldb < max(1,*n)) {
#line 251 "dtrrfs.f"
	*info = -9;
#line 252 "dtrrfs.f"
    } else if (*ldx < max(1,*n)) {
#line 253 "dtrrfs.f"
	*info = -11;
#line 254 "dtrrfs.f"
    }
#line 255 "dtrrfs.f"
    if (*info != 0) {
#line 256 "dtrrfs.f"
	i__1 = -(*info);
#line 256 "dtrrfs.f"
	xerbla_("DTRRFS", &i__1, (ftnlen)6);
#line 257 "dtrrfs.f"
	return 0;
#line 258 "dtrrfs.f"
    }

/*     Quick return if possible */

#line 262 "dtrrfs.f"
    if (*n == 0 || *nrhs == 0) {
#line 263 "dtrrfs.f"
	i__1 = *nrhs;
#line 263 "dtrrfs.f"
	for (j = 1; j <= i__1; ++j) {
#line 264 "dtrrfs.f"
	    ferr[j] = 0.;
#line 265 "dtrrfs.f"
	    berr[j] = 0.;
#line 266 "dtrrfs.f"
/* L10: */
#line 266 "dtrrfs.f"
	}
#line 267 "dtrrfs.f"
	return 0;
#line 268 "dtrrfs.f"
    }

#line 270 "dtrrfs.f"
    if (notran) {
#line 271 "dtrrfs.f"
	*(unsigned char *)transt = 'T';
#line 272 "dtrrfs.f"
    } else {
#line 273 "dtrrfs.f"
	*(unsigned char *)transt = 'N';
#line 274 "dtrrfs.f"
    }

/*     NZ = maximum number of nonzero elements in each row of A, plus 1 */

#line 278 "dtrrfs.f"
    nz = *n + 1;
#line 279 "dtrrfs.f"
    eps = dlamch_("Epsilon", (ftnlen)7);
#line 280 "dtrrfs.f"
    safmin = dlamch_("Safe minimum", (ftnlen)12);
#line 281 "dtrrfs.f"
    safe1 = nz * safmin;
#line 282 "dtrrfs.f"
    safe2 = safe1 / eps;

/*     Do for each right hand side */

#line 286 "dtrrfs.f"
    i__1 = *nrhs;
#line 286 "dtrrfs.f"
    for (j = 1; j <= i__1; ++j) {

/*        Compute residual R = B - op(A) * X, */
/*        where op(A) = A or A**T, depending on TRANS. */

#line 291 "dtrrfs.f"
	dcopy_(n, &x[j * x_dim1 + 1], &c__1, &work[*n + 1], &c__1);
#line 292 "dtrrfs.f"
	dtrmv_(uplo, trans, diag, n, &a[a_offset], lda, &work[*n + 1], &c__1, 
		(ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 293 "dtrrfs.f"
	daxpy_(n, &c_b19, &b[j * b_dim1 + 1], &c__1, &work[*n + 1], &c__1);

/*        Compute componentwise relative backward error from formula */

/*        max(i) ( abs(R(i)) / ( abs(op(A))*abs(X) + abs(B) )(i) ) */

/*        where abs(Z) is the componentwise absolute value of the matrix */
/*        or vector Z.  If the i-th component of the denominator is less */
/*        than SAFE2, then SAFE1 is added to the i-th components of the */
/*        numerator and denominator before dividing. */

#line 304 "dtrrfs.f"
	i__2 = *n;
#line 304 "dtrrfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 305 "dtrrfs.f"
	    work[i__] = (d__1 = b[i__ + j * b_dim1], abs(d__1));
#line 306 "dtrrfs.f"
/* L20: */
#line 306 "dtrrfs.f"
	}

#line 308 "dtrrfs.f"
	if (notran) {

/*           Compute abs(A)*abs(X) + abs(B). */

#line 312 "dtrrfs.f"
	    if (upper) {
#line 313 "dtrrfs.f"
		if (nounit) {
#line 314 "dtrrfs.f"
		    i__2 = *n;
#line 314 "dtrrfs.f"
		    for (k = 1; k <= i__2; ++k) {
#line 315 "dtrrfs.f"
			xk = (d__1 = x[k + j * x_dim1], abs(d__1));
#line 316 "dtrrfs.f"
			i__3 = k;
#line 316 "dtrrfs.f"
			for (i__ = 1; i__ <= i__3; ++i__) {
#line 317 "dtrrfs.f"
			    work[i__] += (d__1 = a[i__ + k * a_dim1], abs(
				    d__1)) * xk;
#line 318 "dtrrfs.f"
/* L30: */
#line 318 "dtrrfs.f"
			}
#line 319 "dtrrfs.f"
/* L40: */
#line 319 "dtrrfs.f"
		    }
#line 320 "dtrrfs.f"
		} else {
#line 321 "dtrrfs.f"
		    i__2 = *n;
#line 321 "dtrrfs.f"
		    for (k = 1; k <= i__2; ++k) {
#line 322 "dtrrfs.f"
			xk = (d__1 = x[k + j * x_dim1], abs(d__1));
#line 323 "dtrrfs.f"
			i__3 = k - 1;
#line 323 "dtrrfs.f"
			for (i__ = 1; i__ <= i__3; ++i__) {
#line 324 "dtrrfs.f"
			    work[i__] += (d__1 = a[i__ + k * a_dim1], abs(
				    d__1)) * xk;
#line 325 "dtrrfs.f"
/* L50: */
#line 325 "dtrrfs.f"
			}
#line 326 "dtrrfs.f"
			work[k] += xk;
#line 327 "dtrrfs.f"
/* L60: */
#line 327 "dtrrfs.f"
		    }
#line 328 "dtrrfs.f"
		}
#line 329 "dtrrfs.f"
	    } else {
#line 330 "dtrrfs.f"
		if (nounit) {
#line 331 "dtrrfs.f"
		    i__2 = *n;
#line 331 "dtrrfs.f"
		    for (k = 1; k <= i__2; ++k) {
#line 332 "dtrrfs.f"
			xk = (d__1 = x[k + j * x_dim1], abs(d__1));
#line 333 "dtrrfs.f"
			i__3 = *n;
#line 333 "dtrrfs.f"
			for (i__ = k; i__ <= i__3; ++i__) {
#line 334 "dtrrfs.f"
			    work[i__] += (d__1 = a[i__ + k * a_dim1], abs(
				    d__1)) * xk;
#line 335 "dtrrfs.f"
/* L70: */
#line 335 "dtrrfs.f"
			}
#line 336 "dtrrfs.f"
/* L80: */
#line 336 "dtrrfs.f"
		    }
#line 337 "dtrrfs.f"
		} else {
#line 338 "dtrrfs.f"
		    i__2 = *n;
#line 338 "dtrrfs.f"
		    for (k = 1; k <= i__2; ++k) {
#line 339 "dtrrfs.f"
			xk = (d__1 = x[k + j * x_dim1], abs(d__1));
#line 340 "dtrrfs.f"
			i__3 = *n;
#line 340 "dtrrfs.f"
			for (i__ = k + 1; i__ <= i__3; ++i__) {
#line 341 "dtrrfs.f"
			    work[i__] += (d__1 = a[i__ + k * a_dim1], abs(
				    d__1)) * xk;
#line 342 "dtrrfs.f"
/* L90: */
#line 342 "dtrrfs.f"
			}
#line 343 "dtrrfs.f"
			work[k] += xk;
#line 344 "dtrrfs.f"
/* L100: */
#line 344 "dtrrfs.f"
		    }
#line 345 "dtrrfs.f"
		}
#line 346 "dtrrfs.f"
	    }
#line 347 "dtrrfs.f"
	} else {

/*           Compute abs(A**T)*abs(X) + abs(B). */

#line 351 "dtrrfs.f"
	    if (upper) {
#line 352 "dtrrfs.f"
		if (nounit) {
#line 353 "dtrrfs.f"
		    i__2 = *n;
#line 353 "dtrrfs.f"
		    for (k = 1; k <= i__2; ++k) {
#line 354 "dtrrfs.f"
			s = 0.;
#line 355 "dtrrfs.f"
			i__3 = k;
#line 355 "dtrrfs.f"
			for (i__ = 1; i__ <= i__3; ++i__) {
#line 356 "dtrrfs.f"
			    s += (d__1 = a[i__ + k * a_dim1], abs(d__1)) * (
				    d__2 = x[i__ + j * x_dim1], abs(d__2));
#line 357 "dtrrfs.f"
/* L110: */
#line 357 "dtrrfs.f"
			}
#line 358 "dtrrfs.f"
			work[k] += s;
#line 359 "dtrrfs.f"
/* L120: */
#line 359 "dtrrfs.f"
		    }
#line 360 "dtrrfs.f"
		} else {
#line 361 "dtrrfs.f"
		    i__2 = *n;
#line 361 "dtrrfs.f"
		    for (k = 1; k <= i__2; ++k) {
#line 362 "dtrrfs.f"
			s = (d__1 = x[k + j * x_dim1], abs(d__1));
#line 363 "dtrrfs.f"
			i__3 = k - 1;
#line 363 "dtrrfs.f"
			for (i__ = 1; i__ <= i__3; ++i__) {
#line 364 "dtrrfs.f"
			    s += (d__1 = a[i__ + k * a_dim1], abs(d__1)) * (
				    d__2 = x[i__ + j * x_dim1], abs(d__2));
#line 365 "dtrrfs.f"
/* L130: */
#line 365 "dtrrfs.f"
			}
#line 366 "dtrrfs.f"
			work[k] += s;
#line 367 "dtrrfs.f"
/* L140: */
#line 367 "dtrrfs.f"
		    }
#line 368 "dtrrfs.f"
		}
#line 369 "dtrrfs.f"
	    } else {
#line 370 "dtrrfs.f"
		if (nounit) {
#line 371 "dtrrfs.f"
		    i__2 = *n;
#line 371 "dtrrfs.f"
		    for (k = 1; k <= i__2; ++k) {
#line 372 "dtrrfs.f"
			s = 0.;
#line 373 "dtrrfs.f"
			i__3 = *n;
#line 373 "dtrrfs.f"
			for (i__ = k; i__ <= i__3; ++i__) {
#line 374 "dtrrfs.f"
			    s += (d__1 = a[i__ + k * a_dim1], abs(d__1)) * (
				    d__2 = x[i__ + j * x_dim1], abs(d__2));
#line 375 "dtrrfs.f"
/* L150: */
#line 375 "dtrrfs.f"
			}
#line 376 "dtrrfs.f"
			work[k] += s;
#line 377 "dtrrfs.f"
/* L160: */
#line 377 "dtrrfs.f"
		    }
#line 378 "dtrrfs.f"
		} else {
#line 379 "dtrrfs.f"
		    i__2 = *n;
#line 379 "dtrrfs.f"
		    for (k = 1; k <= i__2; ++k) {
#line 380 "dtrrfs.f"
			s = (d__1 = x[k + j * x_dim1], abs(d__1));
#line 381 "dtrrfs.f"
			i__3 = *n;
#line 381 "dtrrfs.f"
			for (i__ = k + 1; i__ <= i__3; ++i__) {
#line 382 "dtrrfs.f"
			    s += (d__1 = a[i__ + k * a_dim1], abs(d__1)) * (
				    d__2 = x[i__ + j * x_dim1], abs(d__2));
#line 383 "dtrrfs.f"
/* L170: */
#line 383 "dtrrfs.f"
			}
#line 384 "dtrrfs.f"
			work[k] += s;
#line 385 "dtrrfs.f"
/* L180: */
#line 385 "dtrrfs.f"
		    }
#line 386 "dtrrfs.f"
		}
#line 387 "dtrrfs.f"
	    }
#line 388 "dtrrfs.f"
	}
#line 389 "dtrrfs.f"
	s = 0.;
#line 390 "dtrrfs.f"
	i__2 = *n;
#line 390 "dtrrfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 391 "dtrrfs.f"
	    if (work[i__] > safe2) {
/* Computing MAX */
#line 392 "dtrrfs.f"
		d__2 = s, d__3 = (d__1 = work[*n + i__], abs(d__1)) / work[
			i__];
#line 392 "dtrrfs.f"
		s = max(d__2,d__3);
#line 393 "dtrrfs.f"
	    } else {
/* Computing MAX */
#line 394 "dtrrfs.f"
		d__2 = s, d__3 = ((d__1 = work[*n + i__], abs(d__1)) + safe1) 
			/ (work[i__] + safe1);
#line 394 "dtrrfs.f"
		s = max(d__2,d__3);
#line 396 "dtrrfs.f"
	    }
#line 397 "dtrrfs.f"
/* L190: */
#line 397 "dtrrfs.f"
	}
#line 398 "dtrrfs.f"
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

#line 422 "dtrrfs.f"
	i__2 = *n;
#line 422 "dtrrfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 423 "dtrrfs.f"
	    if (work[i__] > safe2) {
#line 424 "dtrrfs.f"
		work[i__] = (d__1 = work[*n + i__], abs(d__1)) + nz * eps * 
			work[i__];
#line 425 "dtrrfs.f"
	    } else {
#line 426 "dtrrfs.f"
		work[i__] = (d__1 = work[*n + i__], abs(d__1)) + nz * eps * 
			work[i__] + safe1;
#line 427 "dtrrfs.f"
	    }
#line 428 "dtrrfs.f"
/* L200: */
#line 428 "dtrrfs.f"
	}

#line 430 "dtrrfs.f"
	kase = 0;
#line 431 "dtrrfs.f"
L210:
#line 432 "dtrrfs.f"
	dlacn2_(n, &work[(*n << 1) + 1], &work[*n + 1], &iwork[1], &ferr[j], &
		kase, isave);
#line 434 "dtrrfs.f"
	if (kase != 0) {
#line 435 "dtrrfs.f"
	    if (kase == 1) {

/*              Multiply by diag(W)*inv(op(A)**T). */

#line 439 "dtrrfs.f"
		dtrsv_(uplo, transt, diag, n, &a[a_offset], lda, &work[*n + 1]
			, &c__1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 441 "dtrrfs.f"
		i__2 = *n;
#line 441 "dtrrfs.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 442 "dtrrfs.f"
		    work[*n + i__] = work[i__] * work[*n + i__];
#line 443 "dtrrfs.f"
/* L220: */
#line 443 "dtrrfs.f"
		}
#line 444 "dtrrfs.f"
	    } else {

/*              Multiply by inv(op(A))*diag(W). */

#line 448 "dtrrfs.f"
		i__2 = *n;
#line 448 "dtrrfs.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 449 "dtrrfs.f"
		    work[*n + i__] = work[i__] * work[*n + i__];
#line 450 "dtrrfs.f"
/* L230: */
#line 450 "dtrrfs.f"
		}
#line 451 "dtrrfs.f"
		dtrsv_(uplo, trans, diag, n, &a[a_offset], lda, &work[*n + 1],
			 &c__1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 453 "dtrrfs.f"
	    }
#line 454 "dtrrfs.f"
	    goto L210;
#line 455 "dtrrfs.f"
	}

/*        Normalize error. */

#line 459 "dtrrfs.f"
	lstres = 0.;
#line 460 "dtrrfs.f"
	i__2 = *n;
#line 460 "dtrrfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 461 "dtrrfs.f"
	    d__2 = lstres, d__3 = (d__1 = x[i__ + j * x_dim1], abs(d__1));
#line 461 "dtrrfs.f"
	    lstres = max(d__2,d__3);
#line 462 "dtrrfs.f"
/* L240: */
#line 462 "dtrrfs.f"
	}
#line 463 "dtrrfs.f"
	if (lstres != 0.) {
#line 463 "dtrrfs.f"
	    ferr[j] /= lstres;
#line 463 "dtrrfs.f"
	}

#line 466 "dtrrfs.f"
/* L250: */
#line 466 "dtrrfs.f"
    }

#line 468 "dtrrfs.f"
    return 0;

/*     End of DTRRFS */

} /* dtrrfs_ */

