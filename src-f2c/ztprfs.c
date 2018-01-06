#line 1 "ztprfs.f"
/* ztprfs.f -- translated by f2c (version 20100827).
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

#line 1 "ztprfs.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b ZTPRFS */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZTPRFS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztprfs.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztprfs.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztprfs.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZTPRFS( UPLO, TRANS, DIAG, N, NRHS, AP, B, LDB, X, LDX, */
/*                          FERR, BERR, WORK, RWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          DIAG, TRANS, UPLO */
/*       INTEGER            INFO, LDB, LDX, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   BERR( * ), FERR( * ), RWORK( * ) */
/*       COMPLEX*16         AP( * ), B( LDB, * ), WORK( * ), X( LDX, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZTPRFS provides error bounds and backward error estimates for the */
/* > solution to a system of linear equations with a triangular packed */
/* > coefficient matrix. */
/* > */
/* > The solution matrix X must be computed by ZTPTRS or some other */
/* > means before entering this routine.  ZTPRFS does not do iterative */
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
/* >          = 'N':  A * X = B     (No transpose) */
/* >          = 'T':  A**T * X = B  (Transpose) */
/* >          = 'C':  A**H * X = B  (Conjugate transpose) */
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
/* >          AP is COMPLEX*16 array, dimension (N*(N+1)/2) */
/* >          The upper or lower triangular matrix A, packed columnwise in */
/* >          a linear array.  The j-th column of A is stored in the array */
/* >          AP as follows: */
/* >          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j; */
/* >          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n. */
/* >          If DIAG = 'U', the diagonal elements of A are not referenced */
/* >          and are assumed to be 1. */
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
/* > \param[in] X */
/* > \verbatim */
/* >          X is COMPLEX*16 array, dimension (LDX,NRHS) */
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

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup complex16OTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int ztprfs_(char *uplo, char *trans, char *diag, integer *n, 
	integer *nrhs, doublecomplex *ap, doublecomplex *b, integer *ldb, 
	doublecomplex *x, integer *ldx, doublereal *ferr, doublereal *berr, 
	doublecomplex *work, doublereal *rwork, integer *info, ftnlen 
	uplo_len, ftnlen trans_len, ftnlen diag_len)
{
    /* System generated locals */
    integer b_dim1, b_offset, x_dim1, x_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2, d__3, d__4;
    doublecomplex z__1;

    /* Builtin functions */
    double d_imag(doublecomplex *);

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
    static logical upper;
    extern /* Subroutine */ int zcopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), zaxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *), ztpmv_(
	    char *, char *, char *, integer *, doublecomplex *, doublecomplex 
	    *, integer *, ftnlen, ftnlen, ftnlen), ztpsv_(char *, char *, 
	    char *, integer *, doublecomplex *, doublecomplex *, integer *, 
	    ftnlen, ftnlen, ftnlen), zlacn2_(integer *, doublecomplex *, 
	    doublecomplex *, doublereal *, integer *, integer *);
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static logical notran;
    static char transn[1], transt[1];
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
/*     .. Statement Functions .. */
/*     .. */
/*     .. Statement Function definitions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 230 "ztprfs.f"
    /* Parameter adjustments */
#line 230 "ztprfs.f"
    --ap;
#line 230 "ztprfs.f"
    b_dim1 = *ldb;
#line 230 "ztprfs.f"
    b_offset = 1 + b_dim1;
#line 230 "ztprfs.f"
    b -= b_offset;
#line 230 "ztprfs.f"
    x_dim1 = *ldx;
#line 230 "ztprfs.f"
    x_offset = 1 + x_dim1;
#line 230 "ztprfs.f"
    x -= x_offset;
#line 230 "ztprfs.f"
    --ferr;
#line 230 "ztprfs.f"
    --berr;
#line 230 "ztprfs.f"
    --work;
#line 230 "ztprfs.f"
    --rwork;
#line 230 "ztprfs.f"

#line 230 "ztprfs.f"
    /* Function Body */
#line 230 "ztprfs.f"
    *info = 0;
#line 231 "ztprfs.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 232 "ztprfs.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 233 "ztprfs.f"
    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);

#line 235 "ztprfs.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 236 "ztprfs.f"
	*info = -1;
#line 237 "ztprfs.f"
    } else if (! notran && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1) && ! 
	    lsame_(trans, "C", (ftnlen)1, (ftnlen)1)) {
#line 239 "ztprfs.f"
	*info = -2;
#line 240 "ztprfs.f"
    } else if (! nounit && ! lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 241 "ztprfs.f"
	*info = -3;
#line 242 "ztprfs.f"
    } else if (*n < 0) {
#line 243 "ztprfs.f"
	*info = -4;
#line 244 "ztprfs.f"
    } else if (*nrhs < 0) {
#line 245 "ztprfs.f"
	*info = -5;
#line 246 "ztprfs.f"
    } else if (*ldb < max(1,*n)) {
#line 247 "ztprfs.f"
	*info = -8;
#line 248 "ztprfs.f"
    } else if (*ldx < max(1,*n)) {
#line 249 "ztprfs.f"
	*info = -10;
#line 250 "ztprfs.f"
    }
#line 251 "ztprfs.f"
    if (*info != 0) {
#line 252 "ztprfs.f"
	i__1 = -(*info);
#line 252 "ztprfs.f"
	xerbla_("ZTPRFS", &i__1, (ftnlen)6);
#line 253 "ztprfs.f"
	return 0;
#line 254 "ztprfs.f"
    }

/*     Quick return if possible */

#line 258 "ztprfs.f"
    if (*n == 0 || *nrhs == 0) {
#line 259 "ztprfs.f"
	i__1 = *nrhs;
#line 259 "ztprfs.f"
	for (j = 1; j <= i__1; ++j) {
#line 260 "ztprfs.f"
	    ferr[j] = 0.;
#line 261 "ztprfs.f"
	    berr[j] = 0.;
#line 262 "ztprfs.f"
/* L10: */
#line 262 "ztprfs.f"
	}
#line 263 "ztprfs.f"
	return 0;
#line 264 "ztprfs.f"
    }

#line 266 "ztprfs.f"
    if (notran) {
#line 267 "ztprfs.f"
	*(unsigned char *)transn = 'N';
#line 268 "ztprfs.f"
	*(unsigned char *)transt = 'C';
#line 269 "ztprfs.f"
    } else {
#line 270 "ztprfs.f"
	*(unsigned char *)transn = 'C';
#line 271 "ztprfs.f"
	*(unsigned char *)transt = 'N';
#line 272 "ztprfs.f"
    }

/*     NZ = maximum number of nonzero elements in each row of A, plus 1 */

#line 276 "ztprfs.f"
    nz = *n + 1;
#line 277 "ztprfs.f"
    eps = dlamch_("Epsilon", (ftnlen)7);
#line 278 "ztprfs.f"
    safmin = dlamch_("Safe minimum", (ftnlen)12);
#line 279 "ztprfs.f"
    safe1 = nz * safmin;
#line 280 "ztprfs.f"
    safe2 = safe1 / eps;

/*     Do for each right hand side */

#line 284 "ztprfs.f"
    i__1 = *nrhs;
#line 284 "ztprfs.f"
    for (j = 1; j <= i__1; ++j) {

/*        Compute residual R = B - op(A) * X, */
/*        where op(A) = A, A**T, or A**H, depending on TRANS. */

#line 289 "ztprfs.f"
	zcopy_(n, &x[j * x_dim1 + 1], &c__1, &work[1], &c__1);
#line 290 "ztprfs.f"
	ztpmv_(uplo, trans, diag, n, &ap[1], &work[1], &c__1, (ftnlen)1, (
		ftnlen)1, (ftnlen)1);
#line 291 "ztprfs.f"
	z__1.r = -1., z__1.i = -0.;
#line 291 "ztprfs.f"
	zaxpy_(n, &z__1, &b[j * b_dim1 + 1], &c__1, &work[1], &c__1);

/*        Compute componentwise relative backward error from formula */

/*        max(i) ( abs(R(i)) / ( abs(op(A))*abs(X) + abs(B) )(i) ) */

/*        where abs(Z) is the componentwise absolute value of the matrix */
/*        or vector Z.  If the i-th component of the denominator is less */
/*        than SAFE2, then SAFE1 is added to the i-th components of the */
/*        numerator and denominator before dividing. */

#line 302 "ztprfs.f"
	i__2 = *n;
#line 302 "ztprfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 303 "ztprfs.f"
	    i__3 = i__ + j * b_dim1;
#line 303 "ztprfs.f"
	    rwork[i__] = (d__1 = b[i__3].r, abs(d__1)) + (d__2 = d_imag(&b[
		    i__ + j * b_dim1]), abs(d__2));
#line 304 "ztprfs.f"
/* L20: */
#line 304 "ztprfs.f"
	}

#line 306 "ztprfs.f"
	if (notran) {

/*           Compute abs(A)*abs(X) + abs(B). */

#line 310 "ztprfs.f"
	    if (upper) {
#line 311 "ztprfs.f"
		kc = 1;
#line 312 "ztprfs.f"
		if (nounit) {
#line 313 "ztprfs.f"
		    i__2 = *n;
#line 313 "ztprfs.f"
		    for (k = 1; k <= i__2; ++k) {
#line 314 "ztprfs.f"
			i__3 = k + j * x_dim1;
#line 314 "ztprfs.f"
			xk = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = d_imag(&
				x[k + j * x_dim1]), abs(d__2));
#line 315 "ztprfs.f"
			i__3 = k;
#line 315 "ztprfs.f"
			for (i__ = 1; i__ <= i__3; ++i__) {
#line 316 "ztprfs.f"
			    i__4 = kc + i__ - 1;
#line 316 "ztprfs.f"
			    rwork[i__] += ((d__1 = ap[i__4].r, abs(d__1)) + (
				    d__2 = d_imag(&ap[kc + i__ - 1]), abs(
				    d__2))) * xk;
#line 318 "ztprfs.f"
/* L30: */
#line 318 "ztprfs.f"
			}
#line 319 "ztprfs.f"
			kc += k;
#line 320 "ztprfs.f"
/* L40: */
#line 320 "ztprfs.f"
		    }
#line 321 "ztprfs.f"
		} else {
#line 322 "ztprfs.f"
		    i__2 = *n;
#line 322 "ztprfs.f"
		    for (k = 1; k <= i__2; ++k) {
#line 323 "ztprfs.f"
			i__3 = k + j * x_dim1;
#line 323 "ztprfs.f"
			xk = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = d_imag(&
				x[k + j * x_dim1]), abs(d__2));
#line 324 "ztprfs.f"
			i__3 = k - 1;
#line 324 "ztprfs.f"
			for (i__ = 1; i__ <= i__3; ++i__) {
#line 325 "ztprfs.f"
			    i__4 = kc + i__ - 1;
#line 325 "ztprfs.f"
			    rwork[i__] += ((d__1 = ap[i__4].r, abs(d__1)) + (
				    d__2 = d_imag(&ap[kc + i__ - 1]), abs(
				    d__2))) * xk;
#line 327 "ztprfs.f"
/* L50: */
#line 327 "ztprfs.f"
			}
#line 328 "ztprfs.f"
			rwork[k] += xk;
#line 329 "ztprfs.f"
			kc += k;
#line 330 "ztprfs.f"
/* L60: */
#line 330 "ztprfs.f"
		    }
#line 331 "ztprfs.f"
		}
#line 332 "ztprfs.f"
	    } else {
#line 333 "ztprfs.f"
		kc = 1;
#line 334 "ztprfs.f"
		if (nounit) {
#line 335 "ztprfs.f"
		    i__2 = *n;
#line 335 "ztprfs.f"
		    for (k = 1; k <= i__2; ++k) {
#line 336 "ztprfs.f"
			i__3 = k + j * x_dim1;
#line 336 "ztprfs.f"
			xk = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = d_imag(&
				x[k + j * x_dim1]), abs(d__2));
#line 337 "ztprfs.f"
			i__3 = *n;
#line 337 "ztprfs.f"
			for (i__ = k; i__ <= i__3; ++i__) {
#line 338 "ztprfs.f"
			    i__4 = kc + i__ - k;
#line 338 "ztprfs.f"
			    rwork[i__] += ((d__1 = ap[i__4].r, abs(d__1)) + (
				    d__2 = d_imag(&ap[kc + i__ - k]), abs(
				    d__2))) * xk;
#line 340 "ztprfs.f"
/* L70: */
#line 340 "ztprfs.f"
			}
#line 341 "ztprfs.f"
			kc = kc + *n - k + 1;
#line 342 "ztprfs.f"
/* L80: */
#line 342 "ztprfs.f"
		    }
#line 343 "ztprfs.f"
		} else {
#line 344 "ztprfs.f"
		    i__2 = *n;
#line 344 "ztprfs.f"
		    for (k = 1; k <= i__2; ++k) {
#line 345 "ztprfs.f"
			i__3 = k + j * x_dim1;
#line 345 "ztprfs.f"
			xk = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = d_imag(&
				x[k + j * x_dim1]), abs(d__2));
#line 346 "ztprfs.f"
			i__3 = *n;
#line 346 "ztprfs.f"
			for (i__ = k + 1; i__ <= i__3; ++i__) {
#line 347 "ztprfs.f"
			    i__4 = kc + i__ - k;
#line 347 "ztprfs.f"
			    rwork[i__] += ((d__1 = ap[i__4].r, abs(d__1)) + (
				    d__2 = d_imag(&ap[kc + i__ - k]), abs(
				    d__2))) * xk;
#line 349 "ztprfs.f"
/* L90: */
#line 349 "ztprfs.f"
			}
#line 350 "ztprfs.f"
			rwork[k] += xk;
#line 351 "ztprfs.f"
			kc = kc + *n - k + 1;
#line 352 "ztprfs.f"
/* L100: */
#line 352 "ztprfs.f"
		    }
#line 353 "ztprfs.f"
		}
#line 354 "ztprfs.f"
	    }
#line 355 "ztprfs.f"
	} else {

/*           Compute abs(A**H)*abs(X) + abs(B). */

#line 359 "ztprfs.f"
	    if (upper) {
#line 360 "ztprfs.f"
		kc = 1;
#line 361 "ztprfs.f"
		if (nounit) {
#line 362 "ztprfs.f"
		    i__2 = *n;
#line 362 "ztprfs.f"
		    for (k = 1; k <= i__2; ++k) {
#line 363 "ztprfs.f"
			s = 0.;
#line 364 "ztprfs.f"
			i__3 = k;
#line 364 "ztprfs.f"
			for (i__ = 1; i__ <= i__3; ++i__) {
#line 365 "ztprfs.f"
			    i__4 = kc + i__ - 1;
#line 365 "ztprfs.f"
			    i__5 = i__ + j * x_dim1;
#line 365 "ztprfs.f"
			    s += ((d__1 = ap[i__4].r, abs(d__1)) + (d__2 = 
				    d_imag(&ap[kc + i__ - 1]), abs(d__2))) * (
				    (d__3 = x[i__5].r, abs(d__3)) + (d__4 = 
				    d_imag(&x[i__ + j * x_dim1]), abs(d__4)));
#line 366 "ztprfs.f"
/* L110: */
#line 366 "ztprfs.f"
			}
#line 367 "ztprfs.f"
			rwork[k] += s;
#line 368 "ztprfs.f"
			kc += k;
#line 369 "ztprfs.f"
/* L120: */
#line 369 "ztprfs.f"
		    }
#line 370 "ztprfs.f"
		} else {
#line 371 "ztprfs.f"
		    i__2 = *n;
#line 371 "ztprfs.f"
		    for (k = 1; k <= i__2; ++k) {
#line 372 "ztprfs.f"
			i__3 = k + j * x_dim1;
#line 372 "ztprfs.f"
			s = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = d_imag(&x[
				k + j * x_dim1]), abs(d__2));
#line 373 "ztprfs.f"
			i__3 = k - 1;
#line 373 "ztprfs.f"
			for (i__ = 1; i__ <= i__3; ++i__) {
#line 374 "ztprfs.f"
			    i__4 = kc + i__ - 1;
#line 374 "ztprfs.f"
			    i__5 = i__ + j * x_dim1;
#line 374 "ztprfs.f"
			    s += ((d__1 = ap[i__4].r, abs(d__1)) + (d__2 = 
				    d_imag(&ap[kc + i__ - 1]), abs(d__2))) * (
				    (d__3 = x[i__5].r, abs(d__3)) + (d__4 = 
				    d_imag(&x[i__ + j * x_dim1]), abs(d__4)));
#line 375 "ztprfs.f"
/* L130: */
#line 375 "ztprfs.f"
			}
#line 376 "ztprfs.f"
			rwork[k] += s;
#line 377 "ztprfs.f"
			kc += k;
#line 378 "ztprfs.f"
/* L140: */
#line 378 "ztprfs.f"
		    }
#line 379 "ztprfs.f"
		}
#line 380 "ztprfs.f"
	    } else {
#line 381 "ztprfs.f"
		kc = 1;
#line 382 "ztprfs.f"
		if (nounit) {
#line 383 "ztprfs.f"
		    i__2 = *n;
#line 383 "ztprfs.f"
		    for (k = 1; k <= i__2; ++k) {
#line 384 "ztprfs.f"
			s = 0.;
#line 385 "ztprfs.f"
			i__3 = *n;
#line 385 "ztprfs.f"
			for (i__ = k; i__ <= i__3; ++i__) {
#line 386 "ztprfs.f"
			    i__4 = kc + i__ - k;
#line 386 "ztprfs.f"
			    i__5 = i__ + j * x_dim1;
#line 386 "ztprfs.f"
			    s += ((d__1 = ap[i__4].r, abs(d__1)) + (d__2 = 
				    d_imag(&ap[kc + i__ - k]), abs(d__2))) * (
				    (d__3 = x[i__5].r, abs(d__3)) + (d__4 = 
				    d_imag(&x[i__ + j * x_dim1]), abs(d__4)));
#line 387 "ztprfs.f"
/* L150: */
#line 387 "ztprfs.f"
			}
#line 388 "ztprfs.f"
			rwork[k] += s;
#line 389 "ztprfs.f"
			kc = kc + *n - k + 1;
#line 390 "ztprfs.f"
/* L160: */
#line 390 "ztprfs.f"
		    }
#line 391 "ztprfs.f"
		} else {
#line 392 "ztprfs.f"
		    i__2 = *n;
#line 392 "ztprfs.f"
		    for (k = 1; k <= i__2; ++k) {
#line 393 "ztprfs.f"
			i__3 = k + j * x_dim1;
#line 393 "ztprfs.f"
			s = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = d_imag(&x[
				k + j * x_dim1]), abs(d__2));
#line 394 "ztprfs.f"
			i__3 = *n;
#line 394 "ztprfs.f"
			for (i__ = k + 1; i__ <= i__3; ++i__) {
#line 395 "ztprfs.f"
			    i__4 = kc + i__ - k;
#line 395 "ztprfs.f"
			    i__5 = i__ + j * x_dim1;
#line 395 "ztprfs.f"
			    s += ((d__1 = ap[i__4].r, abs(d__1)) + (d__2 = 
				    d_imag(&ap[kc + i__ - k]), abs(d__2))) * (
				    (d__3 = x[i__5].r, abs(d__3)) + (d__4 = 
				    d_imag(&x[i__ + j * x_dim1]), abs(d__4)));
#line 396 "ztprfs.f"
/* L170: */
#line 396 "ztprfs.f"
			}
#line 397 "ztprfs.f"
			rwork[k] += s;
#line 398 "ztprfs.f"
			kc = kc + *n - k + 1;
#line 399 "ztprfs.f"
/* L180: */
#line 399 "ztprfs.f"
		    }
#line 400 "ztprfs.f"
		}
#line 401 "ztprfs.f"
	    }
#line 402 "ztprfs.f"
	}
#line 403 "ztprfs.f"
	s = 0.;
#line 404 "ztprfs.f"
	i__2 = *n;
#line 404 "ztprfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 405 "ztprfs.f"
	    if (rwork[i__] > safe2) {
/* Computing MAX */
#line 406 "ztprfs.f"
		i__3 = i__;
#line 406 "ztprfs.f"
		d__3 = s, d__4 = ((d__1 = work[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&work[i__]), abs(d__2))) / rwork[i__];
#line 406 "ztprfs.f"
		s = max(d__3,d__4);
#line 407 "ztprfs.f"
	    } else {
/* Computing MAX */
#line 408 "ztprfs.f"
		i__3 = i__;
#line 408 "ztprfs.f"
		d__3 = s, d__4 = ((d__1 = work[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&work[i__]), abs(d__2)) + safe1) / (rwork[i__] 
			+ safe1);
#line 408 "ztprfs.f"
		s = max(d__3,d__4);
#line 410 "ztprfs.f"
	    }
#line 411 "ztprfs.f"
/* L190: */
#line 411 "ztprfs.f"
	}
#line 412 "ztprfs.f"
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

/*        Use ZLACN2 to estimate the infinity-norm of the matrix */
/*           inv(op(A)) * diag(W), */
/*        where W = abs(R) + NZ*EPS*( abs(op(A))*abs(X)+abs(B) ))) */

#line 436 "ztprfs.f"
	i__2 = *n;
#line 436 "ztprfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 437 "ztprfs.f"
	    if (rwork[i__] > safe2) {
#line 438 "ztprfs.f"
		i__3 = i__;
#line 438 "ztprfs.f"
		rwork[i__] = (d__1 = work[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&work[i__]), abs(d__2)) + nz * eps * rwork[i__]
			;
#line 439 "ztprfs.f"
	    } else {
#line 440 "ztprfs.f"
		i__3 = i__;
#line 440 "ztprfs.f"
		rwork[i__] = (d__1 = work[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&work[i__]), abs(d__2)) + nz * eps * rwork[i__]
			 + safe1;
#line 442 "ztprfs.f"
	    }
#line 443 "ztprfs.f"
/* L200: */
#line 443 "ztprfs.f"
	}

#line 445 "ztprfs.f"
	kase = 0;
#line 446 "ztprfs.f"
L210:
#line 447 "ztprfs.f"
	zlacn2_(n, &work[*n + 1], &work[1], &ferr[j], &kase, isave);
#line 448 "ztprfs.f"
	if (kase != 0) {
#line 449 "ztprfs.f"
	    if (kase == 1) {

/*              Multiply by diag(W)*inv(op(A)**H). */

#line 453 "ztprfs.f"
		ztpsv_(uplo, transt, diag, n, &ap[1], &work[1], &c__1, (
			ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 454 "ztprfs.f"
		i__2 = *n;
#line 454 "ztprfs.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 455 "ztprfs.f"
		    i__3 = i__;
#line 455 "ztprfs.f"
		    i__4 = i__;
#line 455 "ztprfs.f"
		    i__5 = i__;
#line 455 "ztprfs.f"
		    z__1.r = rwork[i__4] * work[i__5].r, z__1.i = rwork[i__4] 
			    * work[i__5].i;
#line 455 "ztprfs.f"
		    work[i__3].r = z__1.r, work[i__3].i = z__1.i;
#line 456 "ztprfs.f"
/* L220: */
#line 456 "ztprfs.f"
		}
#line 457 "ztprfs.f"
	    } else {

/*              Multiply by inv(op(A))*diag(W). */

#line 461 "ztprfs.f"
		i__2 = *n;
#line 461 "ztprfs.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 462 "ztprfs.f"
		    i__3 = i__;
#line 462 "ztprfs.f"
		    i__4 = i__;
#line 462 "ztprfs.f"
		    i__5 = i__;
#line 462 "ztprfs.f"
		    z__1.r = rwork[i__4] * work[i__5].r, z__1.i = rwork[i__4] 
			    * work[i__5].i;
#line 462 "ztprfs.f"
		    work[i__3].r = z__1.r, work[i__3].i = z__1.i;
#line 463 "ztprfs.f"
/* L230: */
#line 463 "ztprfs.f"
		}
#line 464 "ztprfs.f"
		ztpsv_(uplo, transn, diag, n, &ap[1], &work[1], &c__1, (
			ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 465 "ztprfs.f"
	    }
#line 466 "ztprfs.f"
	    goto L210;
#line 467 "ztprfs.f"
	}

/*        Normalize error. */

#line 471 "ztprfs.f"
	lstres = 0.;
#line 472 "ztprfs.f"
	i__2 = *n;
#line 472 "ztprfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 473 "ztprfs.f"
	    i__3 = i__ + j * x_dim1;
#line 473 "ztprfs.f"
	    d__3 = lstres, d__4 = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = 
		    d_imag(&x[i__ + j * x_dim1]), abs(d__2));
#line 473 "ztprfs.f"
	    lstres = max(d__3,d__4);
#line 474 "ztprfs.f"
/* L240: */
#line 474 "ztprfs.f"
	}
#line 475 "ztprfs.f"
	if (lstres != 0.) {
#line 475 "ztprfs.f"
	    ferr[j] /= lstres;
#line 475 "ztprfs.f"
	}

#line 478 "ztprfs.f"
/* L250: */
#line 478 "ztprfs.f"
    }

#line 480 "ztprfs.f"
    return 0;

/*     End of ZTPRFS */

} /* ztprfs_ */

