#line 1 "ctbrfs.f"
/* ctbrfs.f -- translated by f2c (version 20100827).
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

#line 1 "ctbrfs.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b CTBRFS */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CTBRFS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ctbrfs.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ctbrfs.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ctbrfs.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CTBRFS( UPLO, TRANS, DIAG, N, KD, NRHS, AB, LDAB, B, */
/*                          LDB, X, LDX, FERR, BERR, WORK, RWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          DIAG, TRANS, UPLO */
/*       INTEGER            INFO, KD, LDAB, LDB, LDX, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               BERR( * ), FERR( * ), RWORK( * ) */
/*       COMPLEX            AB( LDAB, * ), B( LDB, * ), WORK( * ), */
/*      $                   X( LDX, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CTBRFS provides error bounds and backward error estimates for the */
/* > solution to a system of linear equations with a triangular band */
/* > coefficient matrix. */
/* > */
/* > The solution matrix X must be computed by CTBTRS or some other */
/* > means before entering this routine.  CTBRFS does not do iterative */
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
/* >          AB is COMPLEX array, dimension (LDAB,N) */
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
/* >          B is COMPLEX array, dimension (LDB,NRHS) */
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
/* >          X is COMPLEX array, dimension (LDX,NRHS) */
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
/* >          WORK is COMPLEX array, dimension (2*N) */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is REAL array, dimension (N) */
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

/* > \ingroup complexOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int ctbrfs_(char *uplo, char *trans, char *diag, integer *n, 
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
    extern /* Subroutine */ int ctbmv_(char *, char *, char *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *,
	     ftnlen, ftnlen, ftnlen), ccopy_(integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *), ctbsv_(char *, char *, 
	    char *, integer *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, ftnlen, ftnlen, ftnlen), caxpy_(
	    integer *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    static logical upper;
    extern /* Subroutine */ int clacn2_(integer *, doublecomplex *, 
	    doublecomplex *, doublereal *, integer *, integer *);
    extern doublereal slamch_(char *, ftnlen);
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

#line 245 "ctbrfs.f"
    /* Parameter adjustments */
#line 245 "ctbrfs.f"
    ab_dim1 = *ldab;
#line 245 "ctbrfs.f"
    ab_offset = 1 + ab_dim1;
#line 245 "ctbrfs.f"
    ab -= ab_offset;
#line 245 "ctbrfs.f"
    b_dim1 = *ldb;
#line 245 "ctbrfs.f"
    b_offset = 1 + b_dim1;
#line 245 "ctbrfs.f"
    b -= b_offset;
#line 245 "ctbrfs.f"
    x_dim1 = *ldx;
#line 245 "ctbrfs.f"
    x_offset = 1 + x_dim1;
#line 245 "ctbrfs.f"
    x -= x_offset;
#line 245 "ctbrfs.f"
    --ferr;
#line 245 "ctbrfs.f"
    --berr;
#line 245 "ctbrfs.f"
    --work;
#line 245 "ctbrfs.f"
    --rwork;
#line 245 "ctbrfs.f"

#line 245 "ctbrfs.f"
    /* Function Body */
#line 245 "ctbrfs.f"
    *info = 0;
#line 246 "ctbrfs.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 247 "ctbrfs.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 248 "ctbrfs.f"
    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);

#line 250 "ctbrfs.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 251 "ctbrfs.f"
	*info = -1;
#line 252 "ctbrfs.f"
    } else if (! notran && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1) && ! 
	    lsame_(trans, "C", (ftnlen)1, (ftnlen)1)) {
#line 254 "ctbrfs.f"
	*info = -2;
#line 255 "ctbrfs.f"
    } else if (! nounit && ! lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 256 "ctbrfs.f"
	*info = -3;
#line 257 "ctbrfs.f"
    } else if (*n < 0) {
#line 258 "ctbrfs.f"
	*info = -4;
#line 259 "ctbrfs.f"
    } else if (*kd < 0) {
#line 260 "ctbrfs.f"
	*info = -5;
#line 261 "ctbrfs.f"
    } else if (*nrhs < 0) {
#line 262 "ctbrfs.f"
	*info = -6;
#line 263 "ctbrfs.f"
    } else if (*ldab < *kd + 1) {
#line 264 "ctbrfs.f"
	*info = -8;
#line 265 "ctbrfs.f"
    } else if (*ldb < max(1,*n)) {
#line 266 "ctbrfs.f"
	*info = -10;
#line 267 "ctbrfs.f"
    } else if (*ldx < max(1,*n)) {
#line 268 "ctbrfs.f"
	*info = -12;
#line 269 "ctbrfs.f"
    }
#line 270 "ctbrfs.f"
    if (*info != 0) {
#line 271 "ctbrfs.f"
	i__1 = -(*info);
#line 271 "ctbrfs.f"
	xerbla_("CTBRFS", &i__1, (ftnlen)6);
#line 272 "ctbrfs.f"
	return 0;
#line 273 "ctbrfs.f"
    }

/*     Quick return if possible */

#line 277 "ctbrfs.f"
    if (*n == 0 || *nrhs == 0) {
#line 278 "ctbrfs.f"
	i__1 = *nrhs;
#line 278 "ctbrfs.f"
	for (j = 1; j <= i__1; ++j) {
#line 279 "ctbrfs.f"
	    ferr[j] = 0.;
#line 280 "ctbrfs.f"
	    berr[j] = 0.;
#line 281 "ctbrfs.f"
/* L10: */
#line 281 "ctbrfs.f"
	}
#line 282 "ctbrfs.f"
	return 0;
#line 283 "ctbrfs.f"
    }

#line 285 "ctbrfs.f"
    if (notran) {
#line 286 "ctbrfs.f"
	*(unsigned char *)transn = 'N';
#line 287 "ctbrfs.f"
	*(unsigned char *)transt = 'C';
#line 288 "ctbrfs.f"
    } else {
#line 289 "ctbrfs.f"
	*(unsigned char *)transn = 'C';
#line 290 "ctbrfs.f"
	*(unsigned char *)transt = 'N';
#line 291 "ctbrfs.f"
    }

/*     NZ = maximum number of nonzero elements in each row of A, plus 1 */

#line 295 "ctbrfs.f"
    nz = *kd + 2;
#line 296 "ctbrfs.f"
    eps = slamch_("Epsilon", (ftnlen)7);
#line 297 "ctbrfs.f"
    safmin = slamch_("Safe minimum", (ftnlen)12);
#line 298 "ctbrfs.f"
    safe1 = nz * safmin;
#line 299 "ctbrfs.f"
    safe2 = safe1 / eps;

/*     Do for each right hand side */

#line 303 "ctbrfs.f"
    i__1 = *nrhs;
#line 303 "ctbrfs.f"
    for (j = 1; j <= i__1; ++j) {

/*        Compute residual R = B - op(A) * X, */
/*        where op(A) = A, A**T, or A**H, depending on TRANS. */

#line 308 "ctbrfs.f"
	ccopy_(n, &x[j * x_dim1 + 1], &c__1, &work[1], &c__1);
#line 309 "ctbrfs.f"
	ctbmv_(uplo, trans, diag, n, kd, &ab[ab_offset], ldab, &work[1], &
		c__1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 310 "ctbrfs.f"
	z__1.r = -1., z__1.i = -0.;
#line 310 "ctbrfs.f"
	caxpy_(n, &z__1, &b[j * b_dim1 + 1], &c__1, &work[1], &c__1);

/*        Compute componentwise relative backward error from formula */

/*        max(i) ( abs(R(i)) / ( abs(op(A))*abs(X) + abs(B) )(i) ) */

/*        where abs(Z) is the componentwise absolute value of the matrix */
/*        or vector Z.  If the i-th component of the denominator is less */
/*        than SAFE2, then SAFE1 is added to the i-th components of the */
/*        numerator and denominator before dividing. */

#line 321 "ctbrfs.f"
	i__2 = *n;
#line 321 "ctbrfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 322 "ctbrfs.f"
	    i__3 = i__ + j * b_dim1;
#line 322 "ctbrfs.f"
	    rwork[i__] = (d__1 = b[i__3].r, abs(d__1)) + (d__2 = d_imag(&b[
		    i__ + j * b_dim1]), abs(d__2));
#line 323 "ctbrfs.f"
/* L20: */
#line 323 "ctbrfs.f"
	}

#line 325 "ctbrfs.f"
	if (notran) {

/*           Compute abs(A)*abs(X) + abs(B). */

#line 329 "ctbrfs.f"
	    if (upper) {
#line 330 "ctbrfs.f"
		if (nounit) {
#line 331 "ctbrfs.f"
		    i__2 = *n;
#line 331 "ctbrfs.f"
		    for (k = 1; k <= i__2; ++k) {
#line 332 "ctbrfs.f"
			i__3 = k + j * x_dim1;
#line 332 "ctbrfs.f"
			xk = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = d_imag(&
				x[k + j * x_dim1]), abs(d__2));
/* Computing MAX */
#line 333 "ctbrfs.f"
			i__3 = 1, i__4 = k - *kd;
#line 333 "ctbrfs.f"
			i__5 = k;
#line 333 "ctbrfs.f"
			for (i__ = max(i__3,i__4); i__ <= i__5; ++i__) {
#line 334 "ctbrfs.f"
			    i__3 = *kd + 1 + i__ - k + k * ab_dim1;
#line 334 "ctbrfs.f"
			    rwork[i__] += ((d__1 = ab[i__3].r, abs(d__1)) + (
				    d__2 = d_imag(&ab[*kd + 1 + i__ - k + k * 
				    ab_dim1]), abs(d__2))) * xk;
#line 336 "ctbrfs.f"
/* L30: */
#line 336 "ctbrfs.f"
			}
#line 337 "ctbrfs.f"
/* L40: */
#line 337 "ctbrfs.f"
		    }
#line 338 "ctbrfs.f"
		} else {
#line 339 "ctbrfs.f"
		    i__2 = *n;
#line 339 "ctbrfs.f"
		    for (k = 1; k <= i__2; ++k) {
#line 340 "ctbrfs.f"
			i__5 = k + j * x_dim1;
#line 340 "ctbrfs.f"
			xk = (d__1 = x[i__5].r, abs(d__1)) + (d__2 = d_imag(&
				x[k + j * x_dim1]), abs(d__2));
/* Computing MAX */
#line 341 "ctbrfs.f"
			i__5 = 1, i__3 = k - *kd;
#line 341 "ctbrfs.f"
			i__4 = k - 1;
#line 341 "ctbrfs.f"
			for (i__ = max(i__5,i__3); i__ <= i__4; ++i__) {
#line 342 "ctbrfs.f"
			    i__5 = *kd + 1 + i__ - k + k * ab_dim1;
#line 342 "ctbrfs.f"
			    rwork[i__] += ((d__1 = ab[i__5].r, abs(d__1)) + (
				    d__2 = d_imag(&ab[*kd + 1 + i__ - k + k * 
				    ab_dim1]), abs(d__2))) * xk;
#line 344 "ctbrfs.f"
/* L50: */
#line 344 "ctbrfs.f"
			}
#line 345 "ctbrfs.f"
			rwork[k] += xk;
#line 346 "ctbrfs.f"
/* L60: */
#line 346 "ctbrfs.f"
		    }
#line 347 "ctbrfs.f"
		}
#line 348 "ctbrfs.f"
	    } else {
#line 349 "ctbrfs.f"
		if (nounit) {
#line 350 "ctbrfs.f"
		    i__2 = *n;
#line 350 "ctbrfs.f"
		    for (k = 1; k <= i__2; ++k) {
#line 351 "ctbrfs.f"
			i__4 = k + j * x_dim1;
#line 351 "ctbrfs.f"
			xk = (d__1 = x[i__4].r, abs(d__1)) + (d__2 = d_imag(&
				x[k + j * x_dim1]), abs(d__2));
/* Computing MIN */
#line 352 "ctbrfs.f"
			i__5 = *n, i__3 = k + *kd;
#line 352 "ctbrfs.f"
			i__4 = min(i__5,i__3);
#line 352 "ctbrfs.f"
			for (i__ = k; i__ <= i__4; ++i__) {
#line 353 "ctbrfs.f"
			    i__5 = i__ + 1 - k + k * ab_dim1;
#line 353 "ctbrfs.f"
			    rwork[i__] += ((d__1 = ab[i__5].r, abs(d__1)) + (
				    d__2 = d_imag(&ab[i__ + 1 - k + k * 
				    ab_dim1]), abs(d__2))) * xk;
#line 355 "ctbrfs.f"
/* L70: */
#line 355 "ctbrfs.f"
			}
#line 356 "ctbrfs.f"
/* L80: */
#line 356 "ctbrfs.f"
		    }
#line 357 "ctbrfs.f"
		} else {
#line 358 "ctbrfs.f"
		    i__2 = *n;
#line 358 "ctbrfs.f"
		    for (k = 1; k <= i__2; ++k) {
#line 359 "ctbrfs.f"
			i__4 = k + j * x_dim1;
#line 359 "ctbrfs.f"
			xk = (d__1 = x[i__4].r, abs(d__1)) + (d__2 = d_imag(&
				x[k + j * x_dim1]), abs(d__2));
/* Computing MIN */
#line 360 "ctbrfs.f"
			i__5 = *n, i__3 = k + *kd;
#line 360 "ctbrfs.f"
			i__4 = min(i__5,i__3);
#line 360 "ctbrfs.f"
			for (i__ = k + 1; i__ <= i__4; ++i__) {
#line 361 "ctbrfs.f"
			    i__5 = i__ + 1 - k + k * ab_dim1;
#line 361 "ctbrfs.f"
			    rwork[i__] += ((d__1 = ab[i__5].r, abs(d__1)) + (
				    d__2 = d_imag(&ab[i__ + 1 - k + k * 
				    ab_dim1]), abs(d__2))) * xk;
#line 363 "ctbrfs.f"
/* L90: */
#line 363 "ctbrfs.f"
			}
#line 364 "ctbrfs.f"
			rwork[k] += xk;
#line 365 "ctbrfs.f"
/* L100: */
#line 365 "ctbrfs.f"
		    }
#line 366 "ctbrfs.f"
		}
#line 367 "ctbrfs.f"
	    }
#line 368 "ctbrfs.f"
	} else {

/*           Compute abs(A**H)*abs(X) + abs(B). */

#line 372 "ctbrfs.f"
	    if (upper) {
#line 373 "ctbrfs.f"
		if (nounit) {
#line 374 "ctbrfs.f"
		    i__2 = *n;
#line 374 "ctbrfs.f"
		    for (k = 1; k <= i__2; ++k) {
#line 375 "ctbrfs.f"
			s = 0.;
/* Computing MAX */
#line 376 "ctbrfs.f"
			i__4 = 1, i__5 = k - *kd;
#line 376 "ctbrfs.f"
			i__3 = k;
#line 376 "ctbrfs.f"
			for (i__ = max(i__4,i__5); i__ <= i__3; ++i__) {
#line 377 "ctbrfs.f"
			    i__4 = *kd + 1 + i__ - k + k * ab_dim1;
#line 377 "ctbrfs.f"
			    i__5 = i__ + j * x_dim1;
#line 377 "ctbrfs.f"
			    s += ((d__1 = ab[i__4].r, abs(d__1)) + (d__2 = 
				    d_imag(&ab[*kd + 1 + i__ - k + k * 
				    ab_dim1]), abs(d__2))) * ((d__3 = x[i__5]
				    .r, abs(d__3)) + (d__4 = d_imag(&x[i__ + 
				    j * x_dim1]), abs(d__4)));
#line 379 "ctbrfs.f"
/* L110: */
#line 379 "ctbrfs.f"
			}
#line 380 "ctbrfs.f"
			rwork[k] += s;
#line 381 "ctbrfs.f"
/* L120: */
#line 381 "ctbrfs.f"
		    }
#line 382 "ctbrfs.f"
		} else {
#line 383 "ctbrfs.f"
		    i__2 = *n;
#line 383 "ctbrfs.f"
		    for (k = 1; k <= i__2; ++k) {
#line 384 "ctbrfs.f"
			i__3 = k + j * x_dim1;
#line 384 "ctbrfs.f"
			s = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = d_imag(&x[
				k + j * x_dim1]), abs(d__2));
/* Computing MAX */
#line 385 "ctbrfs.f"
			i__3 = 1, i__4 = k - *kd;
#line 385 "ctbrfs.f"
			i__5 = k - 1;
#line 385 "ctbrfs.f"
			for (i__ = max(i__3,i__4); i__ <= i__5; ++i__) {
#line 386 "ctbrfs.f"
			    i__3 = *kd + 1 + i__ - k + k * ab_dim1;
#line 386 "ctbrfs.f"
			    i__4 = i__ + j * x_dim1;
#line 386 "ctbrfs.f"
			    s += ((d__1 = ab[i__3].r, abs(d__1)) + (d__2 = 
				    d_imag(&ab[*kd + 1 + i__ - k + k * 
				    ab_dim1]), abs(d__2))) * ((d__3 = x[i__4]
				    .r, abs(d__3)) + (d__4 = d_imag(&x[i__ + 
				    j * x_dim1]), abs(d__4)));
#line 388 "ctbrfs.f"
/* L130: */
#line 388 "ctbrfs.f"
			}
#line 389 "ctbrfs.f"
			rwork[k] += s;
#line 390 "ctbrfs.f"
/* L140: */
#line 390 "ctbrfs.f"
		    }
#line 391 "ctbrfs.f"
		}
#line 392 "ctbrfs.f"
	    } else {
#line 393 "ctbrfs.f"
		if (nounit) {
#line 394 "ctbrfs.f"
		    i__2 = *n;
#line 394 "ctbrfs.f"
		    for (k = 1; k <= i__2; ++k) {
#line 395 "ctbrfs.f"
			s = 0.;
/* Computing MIN */
#line 396 "ctbrfs.f"
			i__3 = *n, i__4 = k + *kd;
#line 396 "ctbrfs.f"
			i__5 = min(i__3,i__4);
#line 396 "ctbrfs.f"
			for (i__ = k; i__ <= i__5; ++i__) {
#line 397 "ctbrfs.f"
			    i__3 = i__ + 1 - k + k * ab_dim1;
#line 397 "ctbrfs.f"
			    i__4 = i__ + j * x_dim1;
#line 397 "ctbrfs.f"
			    s += ((d__1 = ab[i__3].r, abs(d__1)) + (d__2 = 
				    d_imag(&ab[i__ + 1 - k + k * ab_dim1]), 
				    abs(d__2))) * ((d__3 = x[i__4].r, abs(
				    d__3)) + (d__4 = d_imag(&x[i__ + j * 
				    x_dim1]), abs(d__4)));
#line 399 "ctbrfs.f"
/* L150: */
#line 399 "ctbrfs.f"
			}
#line 400 "ctbrfs.f"
			rwork[k] += s;
#line 401 "ctbrfs.f"
/* L160: */
#line 401 "ctbrfs.f"
		    }
#line 402 "ctbrfs.f"
		} else {
#line 403 "ctbrfs.f"
		    i__2 = *n;
#line 403 "ctbrfs.f"
		    for (k = 1; k <= i__2; ++k) {
#line 404 "ctbrfs.f"
			i__5 = k + j * x_dim1;
#line 404 "ctbrfs.f"
			s = (d__1 = x[i__5].r, abs(d__1)) + (d__2 = d_imag(&x[
				k + j * x_dim1]), abs(d__2));
/* Computing MIN */
#line 405 "ctbrfs.f"
			i__3 = *n, i__4 = k + *kd;
#line 405 "ctbrfs.f"
			i__5 = min(i__3,i__4);
#line 405 "ctbrfs.f"
			for (i__ = k + 1; i__ <= i__5; ++i__) {
#line 406 "ctbrfs.f"
			    i__3 = i__ + 1 - k + k * ab_dim1;
#line 406 "ctbrfs.f"
			    i__4 = i__ + j * x_dim1;
#line 406 "ctbrfs.f"
			    s += ((d__1 = ab[i__3].r, abs(d__1)) + (d__2 = 
				    d_imag(&ab[i__ + 1 - k + k * ab_dim1]), 
				    abs(d__2))) * ((d__3 = x[i__4].r, abs(
				    d__3)) + (d__4 = d_imag(&x[i__ + j * 
				    x_dim1]), abs(d__4)));
#line 408 "ctbrfs.f"
/* L170: */
#line 408 "ctbrfs.f"
			}
#line 409 "ctbrfs.f"
			rwork[k] += s;
#line 410 "ctbrfs.f"
/* L180: */
#line 410 "ctbrfs.f"
		    }
#line 411 "ctbrfs.f"
		}
#line 412 "ctbrfs.f"
	    }
#line 413 "ctbrfs.f"
	}
#line 414 "ctbrfs.f"
	s = 0.;
#line 415 "ctbrfs.f"
	i__2 = *n;
#line 415 "ctbrfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 416 "ctbrfs.f"
	    if (rwork[i__] > safe2) {
/* Computing MAX */
#line 417 "ctbrfs.f"
		i__5 = i__;
#line 417 "ctbrfs.f"
		d__3 = s, d__4 = ((d__1 = work[i__5].r, abs(d__1)) + (d__2 = 
			d_imag(&work[i__]), abs(d__2))) / rwork[i__];
#line 417 "ctbrfs.f"
		s = max(d__3,d__4);
#line 418 "ctbrfs.f"
	    } else {
/* Computing MAX */
#line 419 "ctbrfs.f"
		i__5 = i__;
#line 419 "ctbrfs.f"
		d__3 = s, d__4 = ((d__1 = work[i__5].r, abs(d__1)) + (d__2 = 
			d_imag(&work[i__]), abs(d__2)) + safe1) / (rwork[i__] 
			+ safe1);
#line 419 "ctbrfs.f"
		s = max(d__3,d__4);
#line 421 "ctbrfs.f"
	    }
#line 422 "ctbrfs.f"
/* L190: */
#line 422 "ctbrfs.f"
	}
#line 423 "ctbrfs.f"
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

/*        Use CLACN2 to estimate the infinity-norm of the matrix */
/*           inv(op(A)) * diag(W), */
/*        where W = abs(R) + NZ*EPS*( abs(op(A))*abs(X)+abs(B) ))) */

#line 447 "ctbrfs.f"
	i__2 = *n;
#line 447 "ctbrfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 448 "ctbrfs.f"
	    if (rwork[i__] > safe2) {
#line 449 "ctbrfs.f"
		i__5 = i__;
#line 449 "ctbrfs.f"
		rwork[i__] = (d__1 = work[i__5].r, abs(d__1)) + (d__2 = 
			d_imag(&work[i__]), abs(d__2)) + nz * eps * rwork[i__]
			;
#line 450 "ctbrfs.f"
	    } else {
#line 451 "ctbrfs.f"
		i__5 = i__;
#line 451 "ctbrfs.f"
		rwork[i__] = (d__1 = work[i__5].r, abs(d__1)) + (d__2 = 
			d_imag(&work[i__]), abs(d__2)) + nz * eps * rwork[i__]
			 + safe1;
#line 453 "ctbrfs.f"
	    }
#line 454 "ctbrfs.f"
/* L200: */
#line 454 "ctbrfs.f"
	}

#line 456 "ctbrfs.f"
	kase = 0;
#line 457 "ctbrfs.f"
L210:
#line 458 "ctbrfs.f"
	clacn2_(n, &work[*n + 1], &work[1], &ferr[j], &kase, isave);
#line 459 "ctbrfs.f"
	if (kase != 0) {
#line 460 "ctbrfs.f"
	    if (kase == 1) {

/*              Multiply by diag(W)*inv(op(A)**H). */

#line 464 "ctbrfs.f"
		ctbsv_(uplo, transt, diag, n, kd, &ab[ab_offset], ldab, &work[
			1], &c__1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 466 "ctbrfs.f"
		i__2 = *n;
#line 466 "ctbrfs.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 467 "ctbrfs.f"
		    i__5 = i__;
#line 467 "ctbrfs.f"
		    i__3 = i__;
#line 467 "ctbrfs.f"
		    i__4 = i__;
#line 467 "ctbrfs.f"
		    z__1.r = rwork[i__3] * work[i__4].r, z__1.i = rwork[i__3] 
			    * work[i__4].i;
#line 467 "ctbrfs.f"
		    work[i__5].r = z__1.r, work[i__5].i = z__1.i;
#line 468 "ctbrfs.f"
/* L220: */
#line 468 "ctbrfs.f"
		}
#line 469 "ctbrfs.f"
	    } else {

/*              Multiply by inv(op(A))*diag(W). */

#line 473 "ctbrfs.f"
		i__2 = *n;
#line 473 "ctbrfs.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 474 "ctbrfs.f"
		    i__5 = i__;
#line 474 "ctbrfs.f"
		    i__3 = i__;
#line 474 "ctbrfs.f"
		    i__4 = i__;
#line 474 "ctbrfs.f"
		    z__1.r = rwork[i__3] * work[i__4].r, z__1.i = rwork[i__3] 
			    * work[i__4].i;
#line 474 "ctbrfs.f"
		    work[i__5].r = z__1.r, work[i__5].i = z__1.i;
#line 475 "ctbrfs.f"
/* L230: */
#line 475 "ctbrfs.f"
		}
#line 476 "ctbrfs.f"
		ctbsv_(uplo, transn, diag, n, kd, &ab[ab_offset], ldab, &work[
			1], &c__1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 478 "ctbrfs.f"
	    }
#line 479 "ctbrfs.f"
	    goto L210;
#line 480 "ctbrfs.f"
	}

/*        Normalize error. */

#line 484 "ctbrfs.f"
	lstres = 0.;
#line 485 "ctbrfs.f"
	i__2 = *n;
#line 485 "ctbrfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 486 "ctbrfs.f"
	    i__5 = i__ + j * x_dim1;
#line 486 "ctbrfs.f"
	    d__3 = lstres, d__4 = (d__1 = x[i__5].r, abs(d__1)) + (d__2 = 
		    d_imag(&x[i__ + j * x_dim1]), abs(d__2));
#line 486 "ctbrfs.f"
	    lstres = max(d__3,d__4);
#line 487 "ctbrfs.f"
/* L240: */
#line 487 "ctbrfs.f"
	}
#line 488 "ctbrfs.f"
	if (lstres != 0.) {
#line 488 "ctbrfs.f"
	    ferr[j] /= lstres;
#line 488 "ctbrfs.f"
	}

#line 491 "ctbrfs.f"
/* L250: */
#line 491 "ctbrfs.f"
    }

#line 493 "ctbrfs.f"
    return 0;

/*     End of CTBRFS */

} /* ctbrfs_ */

