#line 1 "ztbrfs.f"
/* ztbrfs.f -- translated by f2c (version 20100827).
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

#line 1 "ztbrfs.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b ZTBRFS */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZTBRFS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztbrfs.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztbrfs.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztbrfs.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZTBRFS( UPLO, TRANS, DIAG, N, KD, NRHS, AB, LDAB, B, */
/*                          LDB, X, LDX, FERR, BERR, WORK, RWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          DIAG, TRANS, UPLO */
/*       INTEGER            INFO, KD, LDAB, LDB, LDX, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   BERR( * ), FERR( * ), RWORK( * ) */
/*       COMPLEX*16         AB( LDAB, * ), B( LDB, * ), WORK( * ), */
/*      $                   X( LDX, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZTBRFS provides error bounds and backward error estimates for the */
/* > solution to a system of linear equations with a triangular band */
/* > coefficient matrix. */
/* > */
/* > The solution matrix X must be computed by ZTBTRS or some other */
/* > means before entering this routine.  ZTBRFS does not do iterative */
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
/* >          AB is COMPLEX*16 array, dimension (LDAB,N) */
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

/* > \date November 2011 */

/* > \ingroup complex16OTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int ztbrfs_(char *uplo, char *trans, char *diag, integer *n, 
	integer *kd, integer *nrhs, doublecomplex *ab, integer *ldab, 
	doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx, 
	doublereal *ferr, doublereal *berr, doublecomplex *work, doublereal *
	rwork, integer *info, ftnlen uplo_len, ftnlen trans_len, ftnlen 
	diag_len)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, b_dim1, b_offset, x_dim1, x_offset, i__1, 
	    i__2, i__3, i__4, i__5;
    doublereal d__1, d__2, d__3, d__4;
    doublecomplex z__1;

    /* Builtin functions */
    double d_imag(doublecomplex *);

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
    extern /* Subroutine */ int ztbmv_(char *, char *, char *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *,
	     ftnlen, ftnlen, ftnlen), zcopy_(integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *), ztbsv_(char *, char *, 
	    char *, integer *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, ftnlen, ftnlen, ftnlen), zaxpy_(
	    integer *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), zlacn2_(integer *, doublecomplex *, 
	    doublecomplex *, doublereal *, integer *, integer *);
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static logical notran;
    static char transn[1], transt[1];
    static logical nounit;
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
/*     .. Statement Functions .. */
/*     .. */
/*     .. Statement Function definitions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 245 "ztbrfs.f"
    /* Parameter adjustments */
#line 245 "ztbrfs.f"
    ab_dim1 = *ldab;
#line 245 "ztbrfs.f"
    ab_offset = 1 + ab_dim1;
#line 245 "ztbrfs.f"
    ab -= ab_offset;
#line 245 "ztbrfs.f"
    b_dim1 = *ldb;
#line 245 "ztbrfs.f"
    b_offset = 1 + b_dim1;
#line 245 "ztbrfs.f"
    b -= b_offset;
#line 245 "ztbrfs.f"
    x_dim1 = *ldx;
#line 245 "ztbrfs.f"
    x_offset = 1 + x_dim1;
#line 245 "ztbrfs.f"
    x -= x_offset;
#line 245 "ztbrfs.f"
    --ferr;
#line 245 "ztbrfs.f"
    --berr;
#line 245 "ztbrfs.f"
    --work;
#line 245 "ztbrfs.f"
    --rwork;
#line 245 "ztbrfs.f"

#line 245 "ztbrfs.f"
    /* Function Body */
#line 245 "ztbrfs.f"
    *info = 0;
#line 246 "ztbrfs.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 247 "ztbrfs.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 248 "ztbrfs.f"
    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);

#line 250 "ztbrfs.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 251 "ztbrfs.f"
	*info = -1;
#line 252 "ztbrfs.f"
    } else if (! notran && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1) && ! 
	    lsame_(trans, "C", (ftnlen)1, (ftnlen)1)) {
#line 254 "ztbrfs.f"
	*info = -2;
#line 255 "ztbrfs.f"
    } else if (! nounit && ! lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 256 "ztbrfs.f"
	*info = -3;
#line 257 "ztbrfs.f"
    } else if (*n < 0) {
#line 258 "ztbrfs.f"
	*info = -4;
#line 259 "ztbrfs.f"
    } else if (*kd < 0) {
#line 260 "ztbrfs.f"
	*info = -5;
#line 261 "ztbrfs.f"
    } else if (*nrhs < 0) {
#line 262 "ztbrfs.f"
	*info = -6;
#line 263 "ztbrfs.f"
    } else if (*ldab < *kd + 1) {
#line 264 "ztbrfs.f"
	*info = -8;
#line 265 "ztbrfs.f"
    } else if (*ldb < max(1,*n)) {
#line 266 "ztbrfs.f"
	*info = -10;
#line 267 "ztbrfs.f"
    } else if (*ldx < max(1,*n)) {
#line 268 "ztbrfs.f"
	*info = -12;
#line 269 "ztbrfs.f"
    }
#line 270 "ztbrfs.f"
    if (*info != 0) {
#line 271 "ztbrfs.f"
	i__1 = -(*info);
#line 271 "ztbrfs.f"
	xerbla_("ZTBRFS", &i__1, (ftnlen)6);
#line 272 "ztbrfs.f"
	return 0;
#line 273 "ztbrfs.f"
    }

/*     Quick return if possible */

#line 277 "ztbrfs.f"
    if (*n == 0 || *nrhs == 0) {
#line 278 "ztbrfs.f"
	i__1 = *nrhs;
#line 278 "ztbrfs.f"
	for (j = 1; j <= i__1; ++j) {
#line 279 "ztbrfs.f"
	    ferr[j] = 0.;
#line 280 "ztbrfs.f"
	    berr[j] = 0.;
#line 281 "ztbrfs.f"
/* L10: */
#line 281 "ztbrfs.f"
	}
#line 282 "ztbrfs.f"
	return 0;
#line 283 "ztbrfs.f"
    }

#line 285 "ztbrfs.f"
    if (notran) {
#line 286 "ztbrfs.f"
	*(unsigned char *)transn = 'N';
#line 287 "ztbrfs.f"
	*(unsigned char *)transt = 'C';
#line 288 "ztbrfs.f"
    } else {
#line 289 "ztbrfs.f"
	*(unsigned char *)transn = 'C';
#line 290 "ztbrfs.f"
	*(unsigned char *)transt = 'N';
#line 291 "ztbrfs.f"
    }

/*     NZ = maximum number of nonzero elements in each row of A, plus 1 */

#line 295 "ztbrfs.f"
    nz = *kd + 2;
#line 296 "ztbrfs.f"
    eps = dlamch_("Epsilon", (ftnlen)7);
#line 297 "ztbrfs.f"
    safmin = dlamch_("Safe minimum", (ftnlen)12);
#line 298 "ztbrfs.f"
    safe1 = nz * safmin;
#line 299 "ztbrfs.f"
    safe2 = safe1 / eps;

/*     Do for each right hand side */

#line 303 "ztbrfs.f"
    i__1 = *nrhs;
#line 303 "ztbrfs.f"
    for (j = 1; j <= i__1; ++j) {

/*        Compute residual R = B - op(A) * X, */
/*        where op(A) = A, A**T, or A**H, depending on TRANS. */

#line 308 "ztbrfs.f"
	zcopy_(n, &x[j * x_dim1 + 1], &c__1, &work[1], &c__1);
#line 309 "ztbrfs.f"
	ztbmv_(uplo, trans, diag, n, kd, &ab[ab_offset], ldab, &work[1], &
		c__1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 310 "ztbrfs.f"
	z__1.r = -1., z__1.i = -0.;
#line 310 "ztbrfs.f"
	zaxpy_(n, &z__1, &b[j * b_dim1 + 1], &c__1, &work[1], &c__1);

/*        Compute componentwise relative backward error from formula */

/*        max(i) ( abs(R(i)) / ( abs(op(A))*abs(X) + abs(B) )(i) ) */

/*        where abs(Z) is the componentwise absolute value of the matrix */
/*        or vector Z.  If the i-th component of the denominator is less */
/*        than SAFE2, then SAFE1 is added to the i-th components of the */
/*        numerator and denominator before dividing. */

#line 321 "ztbrfs.f"
	i__2 = *n;
#line 321 "ztbrfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 322 "ztbrfs.f"
	    i__3 = i__ + j * b_dim1;
#line 322 "ztbrfs.f"
	    rwork[i__] = (d__1 = b[i__3].r, abs(d__1)) + (d__2 = d_imag(&b[
		    i__ + j * b_dim1]), abs(d__2));
#line 323 "ztbrfs.f"
/* L20: */
#line 323 "ztbrfs.f"
	}

#line 325 "ztbrfs.f"
	if (notran) {

/*           Compute abs(A)*abs(X) + abs(B). */

#line 329 "ztbrfs.f"
	    if (upper) {
#line 330 "ztbrfs.f"
		if (nounit) {
#line 331 "ztbrfs.f"
		    i__2 = *n;
#line 331 "ztbrfs.f"
		    for (k = 1; k <= i__2; ++k) {
#line 332 "ztbrfs.f"
			i__3 = k + j * x_dim1;
#line 332 "ztbrfs.f"
			xk = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = d_imag(&
				x[k + j * x_dim1]), abs(d__2));
/* Computing MAX */
#line 333 "ztbrfs.f"
			i__3 = 1, i__4 = k - *kd;
#line 333 "ztbrfs.f"
			i__5 = k;
#line 333 "ztbrfs.f"
			for (i__ = max(i__3,i__4); i__ <= i__5; ++i__) {
#line 334 "ztbrfs.f"
			    i__3 = *kd + 1 + i__ - k + k * ab_dim1;
#line 334 "ztbrfs.f"
			    rwork[i__] += ((d__1 = ab[i__3].r, abs(d__1)) + (
				    d__2 = d_imag(&ab[*kd + 1 + i__ - k + k * 
				    ab_dim1]), abs(d__2))) * xk;
#line 336 "ztbrfs.f"
/* L30: */
#line 336 "ztbrfs.f"
			}
#line 337 "ztbrfs.f"
/* L40: */
#line 337 "ztbrfs.f"
		    }
#line 338 "ztbrfs.f"
		} else {
#line 339 "ztbrfs.f"
		    i__2 = *n;
#line 339 "ztbrfs.f"
		    for (k = 1; k <= i__2; ++k) {
#line 340 "ztbrfs.f"
			i__5 = k + j * x_dim1;
#line 340 "ztbrfs.f"
			xk = (d__1 = x[i__5].r, abs(d__1)) + (d__2 = d_imag(&
				x[k + j * x_dim1]), abs(d__2));
/* Computing MAX */
#line 341 "ztbrfs.f"
			i__5 = 1, i__3 = k - *kd;
#line 341 "ztbrfs.f"
			i__4 = k - 1;
#line 341 "ztbrfs.f"
			for (i__ = max(i__5,i__3); i__ <= i__4; ++i__) {
#line 342 "ztbrfs.f"
			    i__5 = *kd + 1 + i__ - k + k * ab_dim1;
#line 342 "ztbrfs.f"
			    rwork[i__] += ((d__1 = ab[i__5].r, abs(d__1)) + (
				    d__2 = d_imag(&ab[*kd + 1 + i__ - k + k * 
				    ab_dim1]), abs(d__2))) * xk;
#line 344 "ztbrfs.f"
/* L50: */
#line 344 "ztbrfs.f"
			}
#line 345 "ztbrfs.f"
			rwork[k] += xk;
#line 346 "ztbrfs.f"
/* L60: */
#line 346 "ztbrfs.f"
		    }
#line 347 "ztbrfs.f"
		}
#line 348 "ztbrfs.f"
	    } else {
#line 349 "ztbrfs.f"
		if (nounit) {
#line 350 "ztbrfs.f"
		    i__2 = *n;
#line 350 "ztbrfs.f"
		    for (k = 1; k <= i__2; ++k) {
#line 351 "ztbrfs.f"
			i__4 = k + j * x_dim1;
#line 351 "ztbrfs.f"
			xk = (d__1 = x[i__4].r, abs(d__1)) + (d__2 = d_imag(&
				x[k + j * x_dim1]), abs(d__2));
/* Computing MIN */
#line 352 "ztbrfs.f"
			i__5 = *n, i__3 = k + *kd;
#line 352 "ztbrfs.f"
			i__4 = min(i__5,i__3);
#line 352 "ztbrfs.f"
			for (i__ = k; i__ <= i__4; ++i__) {
#line 353 "ztbrfs.f"
			    i__5 = i__ + 1 - k + k * ab_dim1;
#line 353 "ztbrfs.f"
			    rwork[i__] += ((d__1 = ab[i__5].r, abs(d__1)) + (
				    d__2 = d_imag(&ab[i__ + 1 - k + k * 
				    ab_dim1]), abs(d__2))) * xk;
#line 355 "ztbrfs.f"
/* L70: */
#line 355 "ztbrfs.f"
			}
#line 356 "ztbrfs.f"
/* L80: */
#line 356 "ztbrfs.f"
		    }
#line 357 "ztbrfs.f"
		} else {
#line 358 "ztbrfs.f"
		    i__2 = *n;
#line 358 "ztbrfs.f"
		    for (k = 1; k <= i__2; ++k) {
#line 359 "ztbrfs.f"
			i__4 = k + j * x_dim1;
#line 359 "ztbrfs.f"
			xk = (d__1 = x[i__4].r, abs(d__1)) + (d__2 = d_imag(&
				x[k + j * x_dim1]), abs(d__2));
/* Computing MIN */
#line 360 "ztbrfs.f"
			i__5 = *n, i__3 = k + *kd;
#line 360 "ztbrfs.f"
			i__4 = min(i__5,i__3);
#line 360 "ztbrfs.f"
			for (i__ = k + 1; i__ <= i__4; ++i__) {
#line 361 "ztbrfs.f"
			    i__5 = i__ + 1 - k + k * ab_dim1;
#line 361 "ztbrfs.f"
			    rwork[i__] += ((d__1 = ab[i__5].r, abs(d__1)) + (
				    d__2 = d_imag(&ab[i__ + 1 - k + k * 
				    ab_dim1]), abs(d__2))) * xk;
#line 363 "ztbrfs.f"
/* L90: */
#line 363 "ztbrfs.f"
			}
#line 364 "ztbrfs.f"
			rwork[k] += xk;
#line 365 "ztbrfs.f"
/* L100: */
#line 365 "ztbrfs.f"
		    }
#line 366 "ztbrfs.f"
		}
#line 367 "ztbrfs.f"
	    }
#line 368 "ztbrfs.f"
	} else {

/*           Compute abs(A**H)*abs(X) + abs(B). */

#line 372 "ztbrfs.f"
	    if (upper) {
#line 373 "ztbrfs.f"
		if (nounit) {
#line 374 "ztbrfs.f"
		    i__2 = *n;
#line 374 "ztbrfs.f"
		    for (k = 1; k <= i__2; ++k) {
#line 375 "ztbrfs.f"
			s = 0.;
/* Computing MAX */
#line 376 "ztbrfs.f"
			i__4 = 1, i__5 = k - *kd;
#line 376 "ztbrfs.f"
			i__3 = k;
#line 376 "ztbrfs.f"
			for (i__ = max(i__4,i__5); i__ <= i__3; ++i__) {
#line 377 "ztbrfs.f"
			    i__4 = *kd + 1 + i__ - k + k * ab_dim1;
#line 377 "ztbrfs.f"
			    i__5 = i__ + j * x_dim1;
#line 377 "ztbrfs.f"
			    s += ((d__1 = ab[i__4].r, abs(d__1)) + (d__2 = 
				    d_imag(&ab[*kd + 1 + i__ - k + k * 
				    ab_dim1]), abs(d__2))) * ((d__3 = x[i__5]
				    .r, abs(d__3)) + (d__4 = d_imag(&x[i__ + 
				    j * x_dim1]), abs(d__4)));
#line 379 "ztbrfs.f"
/* L110: */
#line 379 "ztbrfs.f"
			}
#line 380 "ztbrfs.f"
			rwork[k] += s;
#line 381 "ztbrfs.f"
/* L120: */
#line 381 "ztbrfs.f"
		    }
#line 382 "ztbrfs.f"
		} else {
#line 383 "ztbrfs.f"
		    i__2 = *n;
#line 383 "ztbrfs.f"
		    for (k = 1; k <= i__2; ++k) {
#line 384 "ztbrfs.f"
			i__3 = k + j * x_dim1;
#line 384 "ztbrfs.f"
			s = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = d_imag(&x[
				k + j * x_dim1]), abs(d__2));
/* Computing MAX */
#line 385 "ztbrfs.f"
			i__3 = 1, i__4 = k - *kd;
#line 385 "ztbrfs.f"
			i__5 = k - 1;
#line 385 "ztbrfs.f"
			for (i__ = max(i__3,i__4); i__ <= i__5; ++i__) {
#line 386 "ztbrfs.f"
			    i__3 = *kd + 1 + i__ - k + k * ab_dim1;
#line 386 "ztbrfs.f"
			    i__4 = i__ + j * x_dim1;
#line 386 "ztbrfs.f"
			    s += ((d__1 = ab[i__3].r, abs(d__1)) + (d__2 = 
				    d_imag(&ab[*kd + 1 + i__ - k + k * 
				    ab_dim1]), abs(d__2))) * ((d__3 = x[i__4]
				    .r, abs(d__3)) + (d__4 = d_imag(&x[i__ + 
				    j * x_dim1]), abs(d__4)));
#line 388 "ztbrfs.f"
/* L130: */
#line 388 "ztbrfs.f"
			}
#line 389 "ztbrfs.f"
			rwork[k] += s;
#line 390 "ztbrfs.f"
/* L140: */
#line 390 "ztbrfs.f"
		    }
#line 391 "ztbrfs.f"
		}
#line 392 "ztbrfs.f"
	    } else {
#line 393 "ztbrfs.f"
		if (nounit) {
#line 394 "ztbrfs.f"
		    i__2 = *n;
#line 394 "ztbrfs.f"
		    for (k = 1; k <= i__2; ++k) {
#line 395 "ztbrfs.f"
			s = 0.;
/* Computing MIN */
#line 396 "ztbrfs.f"
			i__3 = *n, i__4 = k + *kd;
#line 396 "ztbrfs.f"
			i__5 = min(i__3,i__4);
#line 396 "ztbrfs.f"
			for (i__ = k; i__ <= i__5; ++i__) {
#line 397 "ztbrfs.f"
			    i__3 = i__ + 1 - k + k * ab_dim1;
#line 397 "ztbrfs.f"
			    i__4 = i__ + j * x_dim1;
#line 397 "ztbrfs.f"
			    s += ((d__1 = ab[i__3].r, abs(d__1)) + (d__2 = 
				    d_imag(&ab[i__ + 1 - k + k * ab_dim1]), 
				    abs(d__2))) * ((d__3 = x[i__4].r, abs(
				    d__3)) + (d__4 = d_imag(&x[i__ + j * 
				    x_dim1]), abs(d__4)));
#line 399 "ztbrfs.f"
/* L150: */
#line 399 "ztbrfs.f"
			}
#line 400 "ztbrfs.f"
			rwork[k] += s;
#line 401 "ztbrfs.f"
/* L160: */
#line 401 "ztbrfs.f"
		    }
#line 402 "ztbrfs.f"
		} else {
#line 403 "ztbrfs.f"
		    i__2 = *n;
#line 403 "ztbrfs.f"
		    for (k = 1; k <= i__2; ++k) {
#line 404 "ztbrfs.f"
			i__5 = k + j * x_dim1;
#line 404 "ztbrfs.f"
			s = (d__1 = x[i__5].r, abs(d__1)) + (d__2 = d_imag(&x[
				k + j * x_dim1]), abs(d__2));
/* Computing MIN */
#line 405 "ztbrfs.f"
			i__3 = *n, i__4 = k + *kd;
#line 405 "ztbrfs.f"
			i__5 = min(i__3,i__4);
#line 405 "ztbrfs.f"
			for (i__ = k + 1; i__ <= i__5; ++i__) {
#line 406 "ztbrfs.f"
			    i__3 = i__ + 1 - k + k * ab_dim1;
#line 406 "ztbrfs.f"
			    i__4 = i__ + j * x_dim1;
#line 406 "ztbrfs.f"
			    s += ((d__1 = ab[i__3].r, abs(d__1)) + (d__2 = 
				    d_imag(&ab[i__ + 1 - k + k * ab_dim1]), 
				    abs(d__2))) * ((d__3 = x[i__4].r, abs(
				    d__3)) + (d__4 = d_imag(&x[i__ + j * 
				    x_dim1]), abs(d__4)));
#line 408 "ztbrfs.f"
/* L170: */
#line 408 "ztbrfs.f"
			}
#line 409 "ztbrfs.f"
			rwork[k] += s;
#line 410 "ztbrfs.f"
/* L180: */
#line 410 "ztbrfs.f"
		    }
#line 411 "ztbrfs.f"
		}
#line 412 "ztbrfs.f"
	    }
#line 413 "ztbrfs.f"
	}
#line 414 "ztbrfs.f"
	s = 0.;
#line 415 "ztbrfs.f"
	i__2 = *n;
#line 415 "ztbrfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 416 "ztbrfs.f"
	    if (rwork[i__] > safe2) {
/* Computing MAX */
#line 417 "ztbrfs.f"
		i__5 = i__;
#line 417 "ztbrfs.f"
		d__3 = s, d__4 = ((d__1 = work[i__5].r, abs(d__1)) + (d__2 = 
			d_imag(&work[i__]), abs(d__2))) / rwork[i__];
#line 417 "ztbrfs.f"
		s = max(d__3,d__4);
#line 418 "ztbrfs.f"
	    } else {
/* Computing MAX */
#line 419 "ztbrfs.f"
		i__5 = i__;
#line 419 "ztbrfs.f"
		d__3 = s, d__4 = ((d__1 = work[i__5].r, abs(d__1)) + (d__2 = 
			d_imag(&work[i__]), abs(d__2)) + safe1) / (rwork[i__] 
			+ safe1);
#line 419 "ztbrfs.f"
		s = max(d__3,d__4);
#line 421 "ztbrfs.f"
	    }
#line 422 "ztbrfs.f"
/* L190: */
#line 422 "ztbrfs.f"
	}
#line 423 "ztbrfs.f"
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

#line 447 "ztbrfs.f"
	i__2 = *n;
#line 447 "ztbrfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 448 "ztbrfs.f"
	    if (rwork[i__] > safe2) {
#line 449 "ztbrfs.f"
		i__5 = i__;
#line 449 "ztbrfs.f"
		rwork[i__] = (d__1 = work[i__5].r, abs(d__1)) + (d__2 = 
			d_imag(&work[i__]), abs(d__2)) + nz * eps * rwork[i__]
			;
#line 450 "ztbrfs.f"
	    } else {
#line 451 "ztbrfs.f"
		i__5 = i__;
#line 451 "ztbrfs.f"
		rwork[i__] = (d__1 = work[i__5].r, abs(d__1)) + (d__2 = 
			d_imag(&work[i__]), abs(d__2)) + nz * eps * rwork[i__]
			 + safe1;
#line 453 "ztbrfs.f"
	    }
#line 454 "ztbrfs.f"
/* L200: */
#line 454 "ztbrfs.f"
	}

#line 456 "ztbrfs.f"
	kase = 0;
#line 457 "ztbrfs.f"
L210:
#line 458 "ztbrfs.f"
	zlacn2_(n, &work[*n + 1], &work[1], &ferr[j], &kase, isave);
#line 459 "ztbrfs.f"
	if (kase != 0) {
#line 460 "ztbrfs.f"
	    if (kase == 1) {

/*              Multiply by diag(W)*inv(op(A)**H). */

#line 464 "ztbrfs.f"
		ztbsv_(uplo, transt, diag, n, kd, &ab[ab_offset], ldab, &work[
			1], &c__1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 466 "ztbrfs.f"
		i__2 = *n;
#line 466 "ztbrfs.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 467 "ztbrfs.f"
		    i__5 = i__;
#line 467 "ztbrfs.f"
		    i__3 = i__;
#line 467 "ztbrfs.f"
		    i__4 = i__;
#line 467 "ztbrfs.f"
		    z__1.r = rwork[i__3] * work[i__4].r, z__1.i = rwork[i__3] 
			    * work[i__4].i;
#line 467 "ztbrfs.f"
		    work[i__5].r = z__1.r, work[i__5].i = z__1.i;
#line 468 "ztbrfs.f"
/* L220: */
#line 468 "ztbrfs.f"
		}
#line 469 "ztbrfs.f"
	    } else {

/*              Multiply by inv(op(A))*diag(W). */

#line 473 "ztbrfs.f"
		i__2 = *n;
#line 473 "ztbrfs.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 474 "ztbrfs.f"
		    i__5 = i__;
#line 474 "ztbrfs.f"
		    i__3 = i__;
#line 474 "ztbrfs.f"
		    i__4 = i__;
#line 474 "ztbrfs.f"
		    z__1.r = rwork[i__3] * work[i__4].r, z__1.i = rwork[i__3] 
			    * work[i__4].i;
#line 474 "ztbrfs.f"
		    work[i__5].r = z__1.r, work[i__5].i = z__1.i;
#line 475 "ztbrfs.f"
/* L230: */
#line 475 "ztbrfs.f"
		}
#line 476 "ztbrfs.f"
		ztbsv_(uplo, transn, diag, n, kd, &ab[ab_offset], ldab, &work[
			1], &c__1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 478 "ztbrfs.f"
	    }
#line 479 "ztbrfs.f"
	    goto L210;
#line 480 "ztbrfs.f"
	}

/*        Normalize error. */

#line 484 "ztbrfs.f"
	lstres = 0.;
#line 485 "ztbrfs.f"
	i__2 = *n;
#line 485 "ztbrfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 486 "ztbrfs.f"
	    i__5 = i__ + j * x_dim1;
#line 486 "ztbrfs.f"
	    d__3 = lstres, d__4 = (d__1 = x[i__5].r, abs(d__1)) + (d__2 = 
		    d_imag(&x[i__ + j * x_dim1]), abs(d__2));
#line 486 "ztbrfs.f"
	    lstres = max(d__3,d__4);
#line 487 "ztbrfs.f"
/* L240: */
#line 487 "ztbrfs.f"
	}
#line 488 "ztbrfs.f"
	if (lstres != 0.) {
#line 488 "ztbrfs.f"
	    ferr[j] /= lstres;
#line 488 "ztbrfs.f"
	}

#line 491 "ztbrfs.f"
/* L250: */
#line 491 "ztbrfs.f"
    }

#line 493 "ztbrfs.f"
    return 0;

/*     End of ZTBRFS */

} /* ztbrfs_ */

