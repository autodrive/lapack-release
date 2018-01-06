#line 1 "strrfs.f"
/* strrfs.f -- translated by f2c (version 20100827).
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

#line 1 "strrfs.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b19 = -1.;

/* > \brief \b STRRFS */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download STRRFS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/strrfs.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/strrfs.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/strrfs.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE STRRFS( UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDB, X, */
/*                          LDX, FERR, BERR, WORK, IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          DIAG, TRANS, UPLO */
/*       INTEGER            INFO, LDA, LDB, LDX, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       REAL               A( LDA, * ), B( LDB, * ), BERR( * ), FERR( * ), */
/*      $                   WORK( * ), X( LDX, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > STRRFS provides error bounds and backward error estimates for the */
/* > solution to a system of linear equations with a triangular */
/* > coefficient matrix. */
/* > */
/* > The solution matrix X must be computed by STRTRS or some other */
/* > means before entering this routine.  STRRFS does not do iterative */
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
/* >          A is REAL array, dimension (LDA,N) */
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
/* > \param[in] X */
/* > \verbatim */
/* >          X is REAL array, dimension (LDX,NRHS) */
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

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup realOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int strrfs_(char *uplo, char *trans, char *diag, integer *n, 
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
    static logical upper;
    extern /* Subroutine */ int scopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), saxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *), strmv_(char *, 
	    char *, char *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, ftnlen, ftnlen, ftnlen), strsv_(char *, char *, char *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen, ftnlen, ftnlen), slacn2_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, integer *);
    extern doublereal slamch_(char *, ftnlen);
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

#line 232 "strrfs.f"
    /* Parameter adjustments */
#line 232 "strrfs.f"
    a_dim1 = *lda;
#line 232 "strrfs.f"
    a_offset = 1 + a_dim1;
#line 232 "strrfs.f"
    a -= a_offset;
#line 232 "strrfs.f"
    b_dim1 = *ldb;
#line 232 "strrfs.f"
    b_offset = 1 + b_dim1;
#line 232 "strrfs.f"
    b -= b_offset;
#line 232 "strrfs.f"
    x_dim1 = *ldx;
#line 232 "strrfs.f"
    x_offset = 1 + x_dim1;
#line 232 "strrfs.f"
    x -= x_offset;
#line 232 "strrfs.f"
    --ferr;
#line 232 "strrfs.f"
    --berr;
#line 232 "strrfs.f"
    --work;
#line 232 "strrfs.f"
    --iwork;
#line 232 "strrfs.f"

#line 232 "strrfs.f"
    /* Function Body */
#line 232 "strrfs.f"
    *info = 0;
#line 233 "strrfs.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 234 "strrfs.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 235 "strrfs.f"
    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);

#line 237 "strrfs.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 238 "strrfs.f"
	*info = -1;
#line 239 "strrfs.f"
    } else if (! notran && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1) && ! 
	    lsame_(trans, "C", (ftnlen)1, (ftnlen)1)) {
#line 241 "strrfs.f"
	*info = -2;
#line 242 "strrfs.f"
    } else if (! nounit && ! lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 243 "strrfs.f"
	*info = -3;
#line 244 "strrfs.f"
    } else if (*n < 0) {
#line 245 "strrfs.f"
	*info = -4;
#line 246 "strrfs.f"
    } else if (*nrhs < 0) {
#line 247 "strrfs.f"
	*info = -5;
#line 248 "strrfs.f"
    } else if (*lda < max(1,*n)) {
#line 249 "strrfs.f"
	*info = -7;
#line 250 "strrfs.f"
    } else if (*ldb < max(1,*n)) {
#line 251 "strrfs.f"
	*info = -9;
#line 252 "strrfs.f"
    } else if (*ldx < max(1,*n)) {
#line 253 "strrfs.f"
	*info = -11;
#line 254 "strrfs.f"
    }
#line 255 "strrfs.f"
    if (*info != 0) {
#line 256 "strrfs.f"
	i__1 = -(*info);
#line 256 "strrfs.f"
	xerbla_("STRRFS", &i__1, (ftnlen)6);
#line 257 "strrfs.f"
	return 0;
#line 258 "strrfs.f"
    }

/*     Quick return if possible */

#line 262 "strrfs.f"
    if (*n == 0 || *nrhs == 0) {
#line 263 "strrfs.f"
	i__1 = *nrhs;
#line 263 "strrfs.f"
	for (j = 1; j <= i__1; ++j) {
#line 264 "strrfs.f"
	    ferr[j] = 0.;
#line 265 "strrfs.f"
	    berr[j] = 0.;
#line 266 "strrfs.f"
/* L10: */
#line 266 "strrfs.f"
	}
#line 267 "strrfs.f"
	return 0;
#line 268 "strrfs.f"
    }

#line 270 "strrfs.f"
    if (notran) {
#line 271 "strrfs.f"
	*(unsigned char *)transt = 'T';
#line 272 "strrfs.f"
    } else {
#line 273 "strrfs.f"
	*(unsigned char *)transt = 'N';
#line 274 "strrfs.f"
    }

/*     NZ = maximum number of nonzero elements in each row of A, plus 1 */

#line 278 "strrfs.f"
    nz = *n + 1;
#line 279 "strrfs.f"
    eps = slamch_("Epsilon", (ftnlen)7);
#line 280 "strrfs.f"
    safmin = slamch_("Safe minimum", (ftnlen)12);
#line 281 "strrfs.f"
    safe1 = nz * safmin;
#line 282 "strrfs.f"
    safe2 = safe1 / eps;

/*     Do for each right hand side */

#line 286 "strrfs.f"
    i__1 = *nrhs;
#line 286 "strrfs.f"
    for (j = 1; j <= i__1; ++j) {

/*        Compute residual R = B - op(A) * X, */
/*        where op(A) = A or A**T, depending on TRANS. */

#line 291 "strrfs.f"
	scopy_(n, &x[j * x_dim1 + 1], &c__1, &work[*n + 1], &c__1);
#line 292 "strrfs.f"
	strmv_(uplo, trans, diag, n, &a[a_offset], lda, &work[*n + 1], &c__1, 
		(ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 293 "strrfs.f"
	saxpy_(n, &c_b19, &b[j * b_dim1 + 1], &c__1, &work[*n + 1], &c__1);

/*        Compute componentwise relative backward error from formula */

/*        max(i) ( abs(R(i)) / ( abs(op(A))*abs(X) + abs(B) )(i) ) */

/*        where abs(Z) is the componentwise absolute value of the matrix */
/*        or vector Z.  If the i-th component of the denominator is less */
/*        than SAFE2, then SAFE1 is added to the i-th components of the */
/*        numerator and denominator before dividing. */

#line 304 "strrfs.f"
	i__2 = *n;
#line 304 "strrfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 305 "strrfs.f"
	    work[i__] = (d__1 = b[i__ + j * b_dim1], abs(d__1));
#line 306 "strrfs.f"
/* L20: */
#line 306 "strrfs.f"
	}

#line 308 "strrfs.f"
	if (notran) {

/*           Compute abs(A)*abs(X) + abs(B). */

#line 312 "strrfs.f"
	    if (upper) {
#line 313 "strrfs.f"
		if (nounit) {
#line 314 "strrfs.f"
		    i__2 = *n;
#line 314 "strrfs.f"
		    for (k = 1; k <= i__2; ++k) {
#line 315 "strrfs.f"
			xk = (d__1 = x[k + j * x_dim1], abs(d__1));
#line 316 "strrfs.f"
			i__3 = k;
#line 316 "strrfs.f"
			for (i__ = 1; i__ <= i__3; ++i__) {
#line 317 "strrfs.f"
			    work[i__] += (d__1 = a[i__ + k * a_dim1], abs(
				    d__1)) * xk;
#line 318 "strrfs.f"
/* L30: */
#line 318 "strrfs.f"
			}
#line 319 "strrfs.f"
/* L40: */
#line 319 "strrfs.f"
		    }
#line 320 "strrfs.f"
		} else {
#line 321 "strrfs.f"
		    i__2 = *n;
#line 321 "strrfs.f"
		    for (k = 1; k <= i__2; ++k) {
#line 322 "strrfs.f"
			xk = (d__1 = x[k + j * x_dim1], abs(d__1));
#line 323 "strrfs.f"
			i__3 = k - 1;
#line 323 "strrfs.f"
			for (i__ = 1; i__ <= i__3; ++i__) {
#line 324 "strrfs.f"
			    work[i__] += (d__1 = a[i__ + k * a_dim1], abs(
				    d__1)) * xk;
#line 325 "strrfs.f"
/* L50: */
#line 325 "strrfs.f"
			}
#line 326 "strrfs.f"
			work[k] += xk;
#line 327 "strrfs.f"
/* L60: */
#line 327 "strrfs.f"
		    }
#line 328 "strrfs.f"
		}
#line 329 "strrfs.f"
	    } else {
#line 330 "strrfs.f"
		if (nounit) {
#line 331 "strrfs.f"
		    i__2 = *n;
#line 331 "strrfs.f"
		    for (k = 1; k <= i__2; ++k) {
#line 332 "strrfs.f"
			xk = (d__1 = x[k + j * x_dim1], abs(d__1));
#line 333 "strrfs.f"
			i__3 = *n;
#line 333 "strrfs.f"
			for (i__ = k; i__ <= i__3; ++i__) {
#line 334 "strrfs.f"
			    work[i__] += (d__1 = a[i__ + k * a_dim1], abs(
				    d__1)) * xk;
#line 335 "strrfs.f"
/* L70: */
#line 335 "strrfs.f"
			}
#line 336 "strrfs.f"
/* L80: */
#line 336 "strrfs.f"
		    }
#line 337 "strrfs.f"
		} else {
#line 338 "strrfs.f"
		    i__2 = *n;
#line 338 "strrfs.f"
		    for (k = 1; k <= i__2; ++k) {
#line 339 "strrfs.f"
			xk = (d__1 = x[k + j * x_dim1], abs(d__1));
#line 340 "strrfs.f"
			i__3 = *n;
#line 340 "strrfs.f"
			for (i__ = k + 1; i__ <= i__3; ++i__) {
#line 341 "strrfs.f"
			    work[i__] += (d__1 = a[i__ + k * a_dim1], abs(
				    d__1)) * xk;
#line 342 "strrfs.f"
/* L90: */
#line 342 "strrfs.f"
			}
#line 343 "strrfs.f"
			work[k] += xk;
#line 344 "strrfs.f"
/* L100: */
#line 344 "strrfs.f"
		    }
#line 345 "strrfs.f"
		}
#line 346 "strrfs.f"
	    }
#line 347 "strrfs.f"
	} else {

/*           Compute abs(A**T)*abs(X) + abs(B). */

#line 351 "strrfs.f"
	    if (upper) {
#line 352 "strrfs.f"
		if (nounit) {
#line 353 "strrfs.f"
		    i__2 = *n;
#line 353 "strrfs.f"
		    for (k = 1; k <= i__2; ++k) {
#line 354 "strrfs.f"
			s = 0.;
#line 355 "strrfs.f"
			i__3 = k;
#line 355 "strrfs.f"
			for (i__ = 1; i__ <= i__3; ++i__) {
#line 356 "strrfs.f"
			    s += (d__1 = a[i__ + k * a_dim1], abs(d__1)) * (
				    d__2 = x[i__ + j * x_dim1], abs(d__2));
#line 357 "strrfs.f"
/* L110: */
#line 357 "strrfs.f"
			}
#line 358 "strrfs.f"
			work[k] += s;
#line 359 "strrfs.f"
/* L120: */
#line 359 "strrfs.f"
		    }
#line 360 "strrfs.f"
		} else {
#line 361 "strrfs.f"
		    i__2 = *n;
#line 361 "strrfs.f"
		    for (k = 1; k <= i__2; ++k) {
#line 362 "strrfs.f"
			s = (d__1 = x[k + j * x_dim1], abs(d__1));
#line 363 "strrfs.f"
			i__3 = k - 1;
#line 363 "strrfs.f"
			for (i__ = 1; i__ <= i__3; ++i__) {
#line 364 "strrfs.f"
			    s += (d__1 = a[i__ + k * a_dim1], abs(d__1)) * (
				    d__2 = x[i__ + j * x_dim1], abs(d__2));
#line 365 "strrfs.f"
/* L130: */
#line 365 "strrfs.f"
			}
#line 366 "strrfs.f"
			work[k] += s;
#line 367 "strrfs.f"
/* L140: */
#line 367 "strrfs.f"
		    }
#line 368 "strrfs.f"
		}
#line 369 "strrfs.f"
	    } else {
#line 370 "strrfs.f"
		if (nounit) {
#line 371 "strrfs.f"
		    i__2 = *n;
#line 371 "strrfs.f"
		    for (k = 1; k <= i__2; ++k) {
#line 372 "strrfs.f"
			s = 0.;
#line 373 "strrfs.f"
			i__3 = *n;
#line 373 "strrfs.f"
			for (i__ = k; i__ <= i__3; ++i__) {
#line 374 "strrfs.f"
			    s += (d__1 = a[i__ + k * a_dim1], abs(d__1)) * (
				    d__2 = x[i__ + j * x_dim1], abs(d__2));
#line 375 "strrfs.f"
/* L150: */
#line 375 "strrfs.f"
			}
#line 376 "strrfs.f"
			work[k] += s;
#line 377 "strrfs.f"
/* L160: */
#line 377 "strrfs.f"
		    }
#line 378 "strrfs.f"
		} else {
#line 379 "strrfs.f"
		    i__2 = *n;
#line 379 "strrfs.f"
		    for (k = 1; k <= i__2; ++k) {
#line 380 "strrfs.f"
			s = (d__1 = x[k + j * x_dim1], abs(d__1));
#line 381 "strrfs.f"
			i__3 = *n;
#line 381 "strrfs.f"
			for (i__ = k + 1; i__ <= i__3; ++i__) {
#line 382 "strrfs.f"
			    s += (d__1 = a[i__ + k * a_dim1], abs(d__1)) * (
				    d__2 = x[i__ + j * x_dim1], abs(d__2));
#line 383 "strrfs.f"
/* L170: */
#line 383 "strrfs.f"
			}
#line 384 "strrfs.f"
			work[k] += s;
#line 385 "strrfs.f"
/* L180: */
#line 385 "strrfs.f"
		    }
#line 386 "strrfs.f"
		}
#line 387 "strrfs.f"
	    }
#line 388 "strrfs.f"
	}
#line 389 "strrfs.f"
	s = 0.;
#line 390 "strrfs.f"
	i__2 = *n;
#line 390 "strrfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 391 "strrfs.f"
	    if (work[i__] > safe2) {
/* Computing MAX */
#line 392 "strrfs.f"
		d__2 = s, d__3 = (d__1 = work[*n + i__], abs(d__1)) / work[
			i__];
#line 392 "strrfs.f"
		s = max(d__2,d__3);
#line 393 "strrfs.f"
	    } else {
/* Computing MAX */
#line 394 "strrfs.f"
		d__2 = s, d__3 = ((d__1 = work[*n + i__], abs(d__1)) + safe1) 
			/ (work[i__] + safe1);
#line 394 "strrfs.f"
		s = max(d__2,d__3);
#line 396 "strrfs.f"
	    }
#line 397 "strrfs.f"
/* L190: */
#line 397 "strrfs.f"
	}
#line 398 "strrfs.f"
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

/*        Use SLACN2 to estimate the infinity-norm of the matrix */
/*           inv(op(A)) * diag(W), */
/*        where W = abs(R) + NZ*EPS*( abs(op(A))*abs(X)+abs(B) ))) */

#line 422 "strrfs.f"
	i__2 = *n;
#line 422 "strrfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 423 "strrfs.f"
	    if (work[i__] > safe2) {
#line 424 "strrfs.f"
		work[i__] = (d__1 = work[*n + i__], abs(d__1)) + nz * eps * 
			work[i__];
#line 425 "strrfs.f"
	    } else {
#line 426 "strrfs.f"
		work[i__] = (d__1 = work[*n + i__], abs(d__1)) + nz * eps * 
			work[i__] + safe1;
#line 427 "strrfs.f"
	    }
#line 428 "strrfs.f"
/* L200: */
#line 428 "strrfs.f"
	}

#line 430 "strrfs.f"
	kase = 0;
#line 431 "strrfs.f"
L210:
#line 432 "strrfs.f"
	slacn2_(n, &work[(*n << 1) + 1], &work[*n + 1], &iwork[1], &ferr[j], &
		kase, isave);
#line 434 "strrfs.f"
	if (kase != 0) {
#line 435 "strrfs.f"
	    if (kase == 1) {

/*              Multiply by diag(W)*inv(op(A)**T). */

#line 439 "strrfs.f"
		strsv_(uplo, transt, diag, n, &a[a_offset], lda, &work[*n + 1]
			, &c__1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 441 "strrfs.f"
		i__2 = *n;
#line 441 "strrfs.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 442 "strrfs.f"
		    work[*n + i__] = work[i__] * work[*n + i__];
#line 443 "strrfs.f"
/* L220: */
#line 443 "strrfs.f"
		}
#line 444 "strrfs.f"
	    } else {

/*              Multiply by inv(op(A))*diag(W). */

#line 448 "strrfs.f"
		i__2 = *n;
#line 448 "strrfs.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 449 "strrfs.f"
		    work[*n + i__] = work[i__] * work[*n + i__];
#line 450 "strrfs.f"
/* L230: */
#line 450 "strrfs.f"
		}
#line 451 "strrfs.f"
		strsv_(uplo, trans, diag, n, &a[a_offset], lda, &work[*n + 1],
			 &c__1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 453 "strrfs.f"
	    }
#line 454 "strrfs.f"
	    goto L210;
#line 455 "strrfs.f"
	}

/*        Normalize error. */

#line 459 "strrfs.f"
	lstres = 0.;
#line 460 "strrfs.f"
	i__2 = *n;
#line 460 "strrfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 461 "strrfs.f"
	    d__2 = lstres, d__3 = (d__1 = x[i__ + j * x_dim1], abs(d__1));
#line 461 "strrfs.f"
	    lstres = max(d__2,d__3);
#line 462 "strrfs.f"
/* L240: */
#line 462 "strrfs.f"
	}
#line 463 "strrfs.f"
	if (lstres != 0.) {
#line 463 "strrfs.f"
	    ferr[j] /= lstres;
#line 463 "strrfs.f"
	}

#line 466 "strrfs.f"
/* L250: */
#line 466 "strrfs.f"
    }

#line 468 "strrfs.f"
    return 0;

/*     End of STRRFS */

} /* strrfs_ */

