#line 1 "dtprfs.f"
/* dtprfs.f -- translated by f2c (version 20100827).
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

#line 1 "dtprfs.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b19 = -1.;

/* > \brief \b DTPRFS */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DTPRFS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dtprfs.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dtprfs.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dtprfs.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DTPRFS( UPLO, TRANS, DIAG, N, NRHS, AP, B, LDB, X, LDX, */
/*                          FERR, BERR, WORK, IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          DIAG, TRANS, UPLO */
/*       INTEGER            INFO, LDB, LDX, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       DOUBLE PRECISION   AP( * ), B( LDB, * ), BERR( * ), FERR( * ), */
/*      $                   WORK( * ), X( LDX, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DTPRFS provides error bounds and backward error estimates for the */
/* > solution to a system of linear equations with a triangular packed */
/* > coefficient matrix. */
/* > */
/* > The solution matrix X must be computed by DTPTRS or some other */
/* > means before entering this routine.  DTPRFS does not do iterative */
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
/* > \param[in] AP */
/* > \verbatim */
/* >          AP is DOUBLE PRECISION array, dimension (N*(N+1)/2) */
/* >          The upper or lower triangular matrix A, packed columnwise in */
/* >          a linear array.  The j-th column of A is stored in the array */
/* >          AP as follows: */
/* >          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j; */
/* >          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n. */
/* >          If DIAG = 'U', the diagonal elements of A are not referenced */
/* >          and are assumed to be 1. */
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
/* Subroutine */ int dtprfs_(char *uplo, char *trans, char *diag, integer *n, 
	integer *nrhs, doublereal *ap, doublereal *b, integer *ldb, 
	doublereal *x, integer *ldx, doublereal *ferr, doublereal *berr, 
	doublereal *work, integer *iwork, integer *info, ftnlen uplo_len, 
	ftnlen trans_len, ftnlen diag_len)
{
    /* System generated locals */
    integer b_dim1, b_offset, x_dim1, x_offset, i__1, i__2, i__3;
    doublereal d__1, d__2, d__3;

    /* Local variables */
    static integer i__, j, k;
    static doublereal s;
    static integer kc;
    static doublereal xk;
    static integer nz;
    static doublereal eps;
    static integer kase;
    static doublereal safe1, safe2;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer isave[3];
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), daxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *), dtpmv_(char *, 
	    char *, char *, integer *, doublereal *, doublereal *, integer *, 
	    ftnlen, ftnlen, ftnlen);
    static logical upper;
    extern /* Subroutine */ int dtpsv_(char *, char *, char *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen, ftnlen, ftnlen), 
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

#line 225 "dtprfs.f"
    /* Parameter adjustments */
#line 225 "dtprfs.f"
    --ap;
#line 225 "dtprfs.f"
    b_dim1 = *ldb;
#line 225 "dtprfs.f"
    b_offset = 1 + b_dim1;
#line 225 "dtprfs.f"
    b -= b_offset;
#line 225 "dtprfs.f"
    x_dim1 = *ldx;
#line 225 "dtprfs.f"
    x_offset = 1 + x_dim1;
#line 225 "dtprfs.f"
    x -= x_offset;
#line 225 "dtprfs.f"
    --ferr;
#line 225 "dtprfs.f"
    --berr;
#line 225 "dtprfs.f"
    --work;
#line 225 "dtprfs.f"
    --iwork;
#line 225 "dtprfs.f"

#line 225 "dtprfs.f"
    /* Function Body */
#line 225 "dtprfs.f"
    *info = 0;
#line 226 "dtprfs.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 227 "dtprfs.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 228 "dtprfs.f"
    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);

#line 230 "dtprfs.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 231 "dtprfs.f"
	*info = -1;
#line 232 "dtprfs.f"
    } else if (! notran && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1) && ! 
	    lsame_(trans, "C", (ftnlen)1, (ftnlen)1)) {
#line 234 "dtprfs.f"
	*info = -2;
#line 235 "dtprfs.f"
    } else if (! nounit && ! lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 236 "dtprfs.f"
	*info = -3;
#line 237 "dtprfs.f"
    } else if (*n < 0) {
#line 238 "dtprfs.f"
	*info = -4;
#line 239 "dtprfs.f"
    } else if (*nrhs < 0) {
#line 240 "dtprfs.f"
	*info = -5;
#line 241 "dtprfs.f"
    } else if (*ldb < max(1,*n)) {
#line 242 "dtprfs.f"
	*info = -8;
#line 243 "dtprfs.f"
    } else if (*ldx < max(1,*n)) {
#line 244 "dtprfs.f"
	*info = -10;
#line 245 "dtprfs.f"
    }
#line 246 "dtprfs.f"
    if (*info != 0) {
#line 247 "dtprfs.f"
	i__1 = -(*info);
#line 247 "dtprfs.f"
	xerbla_("DTPRFS", &i__1, (ftnlen)6);
#line 248 "dtprfs.f"
	return 0;
#line 249 "dtprfs.f"
    }

/*     Quick return if possible */

#line 253 "dtprfs.f"
    if (*n == 0 || *nrhs == 0) {
#line 254 "dtprfs.f"
	i__1 = *nrhs;
#line 254 "dtprfs.f"
	for (j = 1; j <= i__1; ++j) {
#line 255 "dtprfs.f"
	    ferr[j] = 0.;
#line 256 "dtprfs.f"
	    berr[j] = 0.;
#line 257 "dtprfs.f"
/* L10: */
#line 257 "dtprfs.f"
	}
#line 258 "dtprfs.f"
	return 0;
#line 259 "dtprfs.f"
    }

#line 261 "dtprfs.f"
    if (notran) {
#line 262 "dtprfs.f"
	*(unsigned char *)transt = 'T';
#line 263 "dtprfs.f"
    } else {
#line 264 "dtprfs.f"
	*(unsigned char *)transt = 'N';
#line 265 "dtprfs.f"
    }

/*     NZ = maximum number of nonzero elements in each row of A, plus 1 */

#line 269 "dtprfs.f"
    nz = *n + 1;
#line 270 "dtprfs.f"
    eps = dlamch_("Epsilon", (ftnlen)7);
#line 271 "dtprfs.f"
    safmin = dlamch_("Safe minimum", (ftnlen)12);
#line 272 "dtprfs.f"
    safe1 = nz * safmin;
#line 273 "dtprfs.f"
    safe2 = safe1 / eps;

/*     Do for each right hand side */

#line 277 "dtprfs.f"
    i__1 = *nrhs;
#line 277 "dtprfs.f"
    for (j = 1; j <= i__1; ++j) {

/*        Compute residual R = B - op(A) * X, */
/*        where op(A) = A or A**T, depending on TRANS. */

#line 282 "dtprfs.f"
	dcopy_(n, &x[j * x_dim1 + 1], &c__1, &work[*n + 1], &c__1);
#line 283 "dtprfs.f"
	dtpmv_(uplo, trans, diag, n, &ap[1], &work[*n + 1], &c__1, (ftnlen)1, 
		(ftnlen)1, (ftnlen)1);
#line 284 "dtprfs.f"
	daxpy_(n, &c_b19, &b[j * b_dim1 + 1], &c__1, &work[*n + 1], &c__1);

/*        Compute componentwise relative backward error from formula */

/*        max(i) ( abs(R(i)) / ( abs(op(A))*abs(X) + abs(B) )(i) ) */

/*        where abs(Z) is the componentwise absolute value of the matrix */
/*        or vector Z.  If the i-th component of the denominator is less */
/*        than SAFE2, then SAFE1 is added to the i-th components of the */
/*        numerator and denominator before dividing. */

#line 295 "dtprfs.f"
	i__2 = *n;
#line 295 "dtprfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 296 "dtprfs.f"
	    work[i__] = (d__1 = b[i__ + j * b_dim1], abs(d__1));
#line 297 "dtprfs.f"
/* L20: */
#line 297 "dtprfs.f"
	}

#line 299 "dtprfs.f"
	if (notran) {

/*           Compute abs(A)*abs(X) + abs(B). */

#line 303 "dtprfs.f"
	    if (upper) {
#line 304 "dtprfs.f"
		kc = 1;
#line 305 "dtprfs.f"
		if (nounit) {
#line 306 "dtprfs.f"
		    i__2 = *n;
#line 306 "dtprfs.f"
		    for (k = 1; k <= i__2; ++k) {
#line 307 "dtprfs.f"
			xk = (d__1 = x[k + j * x_dim1], abs(d__1));
#line 308 "dtprfs.f"
			i__3 = k;
#line 308 "dtprfs.f"
			for (i__ = 1; i__ <= i__3; ++i__) {
#line 309 "dtprfs.f"
			    work[i__] += (d__1 = ap[kc + i__ - 1], abs(d__1)) 
				    * xk;
#line 310 "dtprfs.f"
/* L30: */
#line 310 "dtprfs.f"
			}
#line 311 "dtprfs.f"
			kc += k;
#line 312 "dtprfs.f"
/* L40: */
#line 312 "dtprfs.f"
		    }
#line 313 "dtprfs.f"
		} else {
#line 314 "dtprfs.f"
		    i__2 = *n;
#line 314 "dtprfs.f"
		    for (k = 1; k <= i__2; ++k) {
#line 315 "dtprfs.f"
			xk = (d__1 = x[k + j * x_dim1], abs(d__1));
#line 316 "dtprfs.f"
			i__3 = k - 1;
#line 316 "dtprfs.f"
			for (i__ = 1; i__ <= i__3; ++i__) {
#line 317 "dtprfs.f"
			    work[i__] += (d__1 = ap[kc + i__ - 1], abs(d__1)) 
				    * xk;
#line 318 "dtprfs.f"
/* L50: */
#line 318 "dtprfs.f"
			}
#line 319 "dtprfs.f"
			work[k] += xk;
#line 320 "dtprfs.f"
			kc += k;
#line 321 "dtprfs.f"
/* L60: */
#line 321 "dtprfs.f"
		    }
#line 322 "dtprfs.f"
		}
#line 323 "dtprfs.f"
	    } else {
#line 324 "dtprfs.f"
		kc = 1;
#line 325 "dtprfs.f"
		if (nounit) {
#line 326 "dtprfs.f"
		    i__2 = *n;
#line 326 "dtprfs.f"
		    for (k = 1; k <= i__2; ++k) {
#line 327 "dtprfs.f"
			xk = (d__1 = x[k + j * x_dim1], abs(d__1));
#line 328 "dtprfs.f"
			i__3 = *n;
#line 328 "dtprfs.f"
			for (i__ = k; i__ <= i__3; ++i__) {
#line 329 "dtprfs.f"
			    work[i__] += (d__1 = ap[kc + i__ - k], abs(d__1)) 
				    * xk;
#line 330 "dtprfs.f"
/* L70: */
#line 330 "dtprfs.f"
			}
#line 331 "dtprfs.f"
			kc = kc + *n - k + 1;
#line 332 "dtprfs.f"
/* L80: */
#line 332 "dtprfs.f"
		    }
#line 333 "dtprfs.f"
		} else {
#line 334 "dtprfs.f"
		    i__2 = *n;
#line 334 "dtprfs.f"
		    for (k = 1; k <= i__2; ++k) {
#line 335 "dtprfs.f"
			xk = (d__1 = x[k + j * x_dim1], abs(d__1));
#line 336 "dtprfs.f"
			i__3 = *n;
#line 336 "dtprfs.f"
			for (i__ = k + 1; i__ <= i__3; ++i__) {
#line 337 "dtprfs.f"
			    work[i__] += (d__1 = ap[kc + i__ - k], abs(d__1)) 
				    * xk;
#line 338 "dtprfs.f"
/* L90: */
#line 338 "dtprfs.f"
			}
#line 339 "dtprfs.f"
			work[k] += xk;
#line 340 "dtprfs.f"
			kc = kc + *n - k + 1;
#line 341 "dtprfs.f"
/* L100: */
#line 341 "dtprfs.f"
		    }
#line 342 "dtprfs.f"
		}
#line 343 "dtprfs.f"
	    }
#line 344 "dtprfs.f"
	} else {

/*           Compute abs(A**T)*abs(X) + abs(B). */

#line 348 "dtprfs.f"
	    if (upper) {
#line 349 "dtprfs.f"
		kc = 1;
#line 350 "dtprfs.f"
		if (nounit) {
#line 351 "dtprfs.f"
		    i__2 = *n;
#line 351 "dtprfs.f"
		    for (k = 1; k <= i__2; ++k) {
#line 352 "dtprfs.f"
			s = 0.;
#line 353 "dtprfs.f"
			i__3 = k;
#line 353 "dtprfs.f"
			for (i__ = 1; i__ <= i__3; ++i__) {
#line 354 "dtprfs.f"
			    s += (d__1 = ap[kc + i__ - 1], abs(d__1)) * (d__2 
				    = x[i__ + j * x_dim1], abs(d__2));
#line 355 "dtprfs.f"
/* L110: */
#line 355 "dtprfs.f"
			}
#line 356 "dtprfs.f"
			work[k] += s;
#line 357 "dtprfs.f"
			kc += k;
#line 358 "dtprfs.f"
/* L120: */
#line 358 "dtprfs.f"
		    }
#line 359 "dtprfs.f"
		} else {
#line 360 "dtprfs.f"
		    i__2 = *n;
#line 360 "dtprfs.f"
		    for (k = 1; k <= i__2; ++k) {
#line 361 "dtprfs.f"
			s = (d__1 = x[k + j * x_dim1], abs(d__1));
#line 362 "dtprfs.f"
			i__3 = k - 1;
#line 362 "dtprfs.f"
			for (i__ = 1; i__ <= i__3; ++i__) {
#line 363 "dtprfs.f"
			    s += (d__1 = ap[kc + i__ - 1], abs(d__1)) * (d__2 
				    = x[i__ + j * x_dim1], abs(d__2));
#line 364 "dtprfs.f"
/* L130: */
#line 364 "dtprfs.f"
			}
#line 365 "dtprfs.f"
			work[k] += s;
#line 366 "dtprfs.f"
			kc += k;
#line 367 "dtprfs.f"
/* L140: */
#line 367 "dtprfs.f"
		    }
#line 368 "dtprfs.f"
		}
#line 369 "dtprfs.f"
	    } else {
#line 370 "dtprfs.f"
		kc = 1;
#line 371 "dtprfs.f"
		if (nounit) {
#line 372 "dtprfs.f"
		    i__2 = *n;
#line 372 "dtprfs.f"
		    for (k = 1; k <= i__2; ++k) {
#line 373 "dtprfs.f"
			s = 0.;
#line 374 "dtprfs.f"
			i__3 = *n;
#line 374 "dtprfs.f"
			for (i__ = k; i__ <= i__3; ++i__) {
#line 375 "dtprfs.f"
			    s += (d__1 = ap[kc + i__ - k], abs(d__1)) * (d__2 
				    = x[i__ + j * x_dim1], abs(d__2));
#line 376 "dtprfs.f"
/* L150: */
#line 376 "dtprfs.f"
			}
#line 377 "dtprfs.f"
			work[k] += s;
#line 378 "dtprfs.f"
			kc = kc + *n - k + 1;
#line 379 "dtprfs.f"
/* L160: */
#line 379 "dtprfs.f"
		    }
#line 380 "dtprfs.f"
		} else {
#line 381 "dtprfs.f"
		    i__2 = *n;
#line 381 "dtprfs.f"
		    for (k = 1; k <= i__2; ++k) {
#line 382 "dtprfs.f"
			s = (d__1 = x[k + j * x_dim1], abs(d__1));
#line 383 "dtprfs.f"
			i__3 = *n;
#line 383 "dtprfs.f"
			for (i__ = k + 1; i__ <= i__3; ++i__) {
#line 384 "dtprfs.f"
			    s += (d__1 = ap[kc + i__ - k], abs(d__1)) * (d__2 
				    = x[i__ + j * x_dim1], abs(d__2));
#line 385 "dtprfs.f"
/* L170: */
#line 385 "dtprfs.f"
			}
#line 386 "dtprfs.f"
			work[k] += s;
#line 387 "dtprfs.f"
			kc = kc + *n - k + 1;
#line 388 "dtprfs.f"
/* L180: */
#line 388 "dtprfs.f"
		    }
#line 389 "dtprfs.f"
		}
#line 390 "dtprfs.f"
	    }
#line 391 "dtprfs.f"
	}
#line 392 "dtprfs.f"
	s = 0.;
#line 393 "dtprfs.f"
	i__2 = *n;
#line 393 "dtprfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 394 "dtprfs.f"
	    if (work[i__] > safe2) {
/* Computing MAX */
#line 395 "dtprfs.f"
		d__2 = s, d__3 = (d__1 = work[*n + i__], abs(d__1)) / work[
			i__];
#line 395 "dtprfs.f"
		s = max(d__2,d__3);
#line 396 "dtprfs.f"
	    } else {
/* Computing MAX */
#line 397 "dtprfs.f"
		d__2 = s, d__3 = ((d__1 = work[*n + i__], abs(d__1)) + safe1) 
			/ (work[i__] + safe1);
#line 397 "dtprfs.f"
		s = max(d__2,d__3);
#line 399 "dtprfs.f"
	    }
#line 400 "dtprfs.f"
/* L190: */
#line 400 "dtprfs.f"
	}
#line 401 "dtprfs.f"
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

#line 425 "dtprfs.f"
	i__2 = *n;
#line 425 "dtprfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 426 "dtprfs.f"
	    if (work[i__] > safe2) {
#line 427 "dtprfs.f"
		work[i__] = (d__1 = work[*n + i__], abs(d__1)) + nz * eps * 
			work[i__];
#line 428 "dtprfs.f"
	    } else {
#line 429 "dtprfs.f"
		work[i__] = (d__1 = work[*n + i__], abs(d__1)) + nz * eps * 
			work[i__] + safe1;
#line 430 "dtprfs.f"
	    }
#line 431 "dtprfs.f"
/* L200: */
#line 431 "dtprfs.f"
	}

#line 433 "dtprfs.f"
	kase = 0;
#line 434 "dtprfs.f"
L210:
#line 435 "dtprfs.f"
	dlacn2_(n, &work[(*n << 1) + 1], &work[*n + 1], &iwork[1], &ferr[j], &
		kase, isave);
#line 437 "dtprfs.f"
	if (kase != 0) {
#line 438 "dtprfs.f"
	    if (kase == 1) {

/*              Multiply by diag(W)*inv(op(A)**T). */

#line 442 "dtprfs.f"
		dtpsv_(uplo, transt, diag, n, &ap[1], &work[*n + 1], &c__1, (
			ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 443 "dtprfs.f"
		i__2 = *n;
#line 443 "dtprfs.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 444 "dtprfs.f"
		    work[*n + i__] = work[i__] * work[*n + i__];
#line 445 "dtprfs.f"
/* L220: */
#line 445 "dtprfs.f"
		}
#line 446 "dtprfs.f"
	    } else {

/*              Multiply by inv(op(A))*diag(W). */

#line 450 "dtprfs.f"
		i__2 = *n;
#line 450 "dtprfs.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 451 "dtprfs.f"
		    work[*n + i__] = work[i__] * work[*n + i__];
#line 452 "dtprfs.f"
/* L230: */
#line 452 "dtprfs.f"
		}
#line 453 "dtprfs.f"
		dtpsv_(uplo, trans, diag, n, &ap[1], &work[*n + 1], &c__1, (
			ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 454 "dtprfs.f"
	    }
#line 455 "dtprfs.f"
	    goto L210;
#line 456 "dtprfs.f"
	}

/*        Normalize error. */

#line 460 "dtprfs.f"
	lstres = 0.;
#line 461 "dtprfs.f"
	i__2 = *n;
#line 461 "dtprfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 462 "dtprfs.f"
	    d__2 = lstres, d__3 = (d__1 = x[i__ + j * x_dim1], abs(d__1));
#line 462 "dtprfs.f"
	    lstres = max(d__2,d__3);
#line 463 "dtprfs.f"
/* L240: */
#line 463 "dtprfs.f"
	}
#line 464 "dtprfs.f"
	if (lstres != 0.) {
#line 464 "dtprfs.f"
	    ferr[j] /= lstres;
#line 464 "dtprfs.f"
	}

#line 467 "dtprfs.f"
/* L250: */
#line 467 "dtprfs.f"
    }

#line 469 "dtprfs.f"
    return 0;

/*     End of DTPRFS */

} /* dtprfs_ */

