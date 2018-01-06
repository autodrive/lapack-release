#line 1 "ctrrfs.f"
/* ctrrfs.f -- translated by f2c (version 20100827).
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

#line 1 "ctrrfs.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b CTRRFS */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CTRRFS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ctrrfs.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ctrrfs.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ctrrfs.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CTRRFS( UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDB, X, */
/*                          LDX, FERR, BERR, WORK, RWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          DIAG, TRANS, UPLO */
/*       INTEGER            INFO, LDA, LDB, LDX, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               BERR( * ), FERR( * ), RWORK( * ) */
/*       COMPLEX            A( LDA, * ), B( LDB, * ), WORK( * ), */
/*      $                   X( LDX, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CTRRFS provides error bounds and backward error estimates for the */
/* > solution to a system of linear equations with a triangular */
/* > coefficient matrix. */
/* > */
/* > The solution matrix X must be computed by CTRTRS or some other */
/* > means before entering this routine.  CTRRFS does not do iterative */
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
/* > \param[in] A */
/* > \verbatim */
/* >          A is COMPLEX array, dimension (LDA,N) */
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
/* Subroutine */ int ctrrfs_(char *uplo, char *trans, char *diag, integer *n, 
	integer *nrhs, doublecomplex *a, integer *lda, doublecomplex *b, 
	integer *ldb, doublecomplex *x, integer *ldx, doublereal *ferr, 
	doublereal *berr, doublecomplex *work, doublereal *rwork, integer *
	info, ftnlen uplo_len, ftnlen trans_len, ftnlen diag_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, x_dim1, x_offset, i__1, i__2, 
	    i__3, i__4, i__5;
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
    extern /* Subroutine */ int ccopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), caxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    static logical upper;
    extern /* Subroutine */ int ctrmv_(char *, char *, char *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen, 
	    ftnlen, ftnlen), ctrsv_(char *, char *, char *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen, 
	    ftnlen, ftnlen), clacn2_(integer *, doublecomplex *, 
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

#line 239 "ctrrfs.f"
    /* Parameter adjustments */
#line 239 "ctrrfs.f"
    a_dim1 = *lda;
#line 239 "ctrrfs.f"
    a_offset = 1 + a_dim1;
#line 239 "ctrrfs.f"
    a -= a_offset;
#line 239 "ctrrfs.f"
    b_dim1 = *ldb;
#line 239 "ctrrfs.f"
    b_offset = 1 + b_dim1;
#line 239 "ctrrfs.f"
    b -= b_offset;
#line 239 "ctrrfs.f"
    x_dim1 = *ldx;
#line 239 "ctrrfs.f"
    x_offset = 1 + x_dim1;
#line 239 "ctrrfs.f"
    x -= x_offset;
#line 239 "ctrrfs.f"
    --ferr;
#line 239 "ctrrfs.f"
    --berr;
#line 239 "ctrrfs.f"
    --work;
#line 239 "ctrrfs.f"
    --rwork;
#line 239 "ctrrfs.f"

#line 239 "ctrrfs.f"
    /* Function Body */
#line 239 "ctrrfs.f"
    *info = 0;
#line 240 "ctrrfs.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 241 "ctrrfs.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 242 "ctrrfs.f"
    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);

#line 244 "ctrrfs.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 245 "ctrrfs.f"
	*info = -1;
#line 246 "ctrrfs.f"
    } else if (! notran && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1) && ! 
	    lsame_(trans, "C", (ftnlen)1, (ftnlen)1)) {
#line 248 "ctrrfs.f"
	*info = -2;
#line 249 "ctrrfs.f"
    } else if (! nounit && ! lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
#line 250 "ctrrfs.f"
	*info = -3;
#line 251 "ctrrfs.f"
    } else if (*n < 0) {
#line 252 "ctrrfs.f"
	*info = -4;
#line 253 "ctrrfs.f"
    } else if (*nrhs < 0) {
#line 254 "ctrrfs.f"
	*info = -5;
#line 255 "ctrrfs.f"
    } else if (*lda < max(1,*n)) {
#line 256 "ctrrfs.f"
	*info = -7;
#line 257 "ctrrfs.f"
    } else if (*ldb < max(1,*n)) {
#line 258 "ctrrfs.f"
	*info = -9;
#line 259 "ctrrfs.f"
    } else if (*ldx < max(1,*n)) {
#line 260 "ctrrfs.f"
	*info = -11;
#line 261 "ctrrfs.f"
    }
#line 262 "ctrrfs.f"
    if (*info != 0) {
#line 263 "ctrrfs.f"
	i__1 = -(*info);
#line 263 "ctrrfs.f"
	xerbla_("CTRRFS", &i__1, (ftnlen)6);
#line 264 "ctrrfs.f"
	return 0;
#line 265 "ctrrfs.f"
    }

/*     Quick return if possible */

#line 269 "ctrrfs.f"
    if (*n == 0 || *nrhs == 0) {
#line 270 "ctrrfs.f"
	i__1 = *nrhs;
#line 270 "ctrrfs.f"
	for (j = 1; j <= i__1; ++j) {
#line 271 "ctrrfs.f"
	    ferr[j] = 0.;
#line 272 "ctrrfs.f"
	    berr[j] = 0.;
#line 273 "ctrrfs.f"
/* L10: */
#line 273 "ctrrfs.f"
	}
#line 274 "ctrrfs.f"
	return 0;
#line 275 "ctrrfs.f"
    }

#line 277 "ctrrfs.f"
    if (notran) {
#line 278 "ctrrfs.f"
	*(unsigned char *)transn = 'N';
#line 279 "ctrrfs.f"
	*(unsigned char *)transt = 'C';
#line 280 "ctrrfs.f"
    } else {
#line 281 "ctrrfs.f"
	*(unsigned char *)transn = 'C';
#line 282 "ctrrfs.f"
	*(unsigned char *)transt = 'N';
#line 283 "ctrrfs.f"
    }

/*     NZ = maximum number of nonzero elements in each row of A, plus 1 */

#line 287 "ctrrfs.f"
    nz = *n + 1;
#line 288 "ctrrfs.f"
    eps = slamch_("Epsilon", (ftnlen)7);
#line 289 "ctrrfs.f"
    safmin = slamch_("Safe minimum", (ftnlen)12);
#line 290 "ctrrfs.f"
    safe1 = nz * safmin;
#line 291 "ctrrfs.f"
    safe2 = safe1 / eps;

/*     Do for each right hand side */

#line 295 "ctrrfs.f"
    i__1 = *nrhs;
#line 295 "ctrrfs.f"
    for (j = 1; j <= i__1; ++j) {

/*        Compute residual R = B - op(A) * X, */
/*        where op(A) = A, A**T, or A**H, depending on TRANS. */

#line 300 "ctrrfs.f"
	ccopy_(n, &x[j * x_dim1 + 1], &c__1, &work[1], &c__1);
#line 301 "ctrrfs.f"
	ctrmv_(uplo, trans, diag, n, &a[a_offset], lda, &work[1], &c__1, (
		ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 302 "ctrrfs.f"
	z__1.r = -1., z__1.i = -0.;
#line 302 "ctrrfs.f"
	caxpy_(n, &z__1, &b[j * b_dim1 + 1], &c__1, &work[1], &c__1);

/*        Compute componentwise relative backward error from formula */

/*        max(i) ( abs(R(i)) / ( abs(op(A))*abs(X) + abs(B) )(i) ) */

/*        where abs(Z) is the componentwise absolute value of the matrix */
/*        or vector Z.  If the i-th component of the denominator is less */
/*        than SAFE2, then SAFE1 is added to the i-th components of the */
/*        numerator and denominator before dividing. */

#line 313 "ctrrfs.f"
	i__2 = *n;
#line 313 "ctrrfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 314 "ctrrfs.f"
	    i__3 = i__ + j * b_dim1;
#line 314 "ctrrfs.f"
	    rwork[i__] = (d__1 = b[i__3].r, abs(d__1)) + (d__2 = d_imag(&b[
		    i__ + j * b_dim1]), abs(d__2));
#line 315 "ctrrfs.f"
/* L20: */
#line 315 "ctrrfs.f"
	}

#line 317 "ctrrfs.f"
	if (notran) {

/*           Compute abs(A)*abs(X) + abs(B). */

#line 321 "ctrrfs.f"
	    if (upper) {
#line 322 "ctrrfs.f"
		if (nounit) {
#line 323 "ctrrfs.f"
		    i__2 = *n;
#line 323 "ctrrfs.f"
		    for (k = 1; k <= i__2; ++k) {
#line 324 "ctrrfs.f"
			i__3 = k + j * x_dim1;
#line 324 "ctrrfs.f"
			xk = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = d_imag(&
				x[k + j * x_dim1]), abs(d__2));
#line 325 "ctrrfs.f"
			i__3 = k;
#line 325 "ctrrfs.f"
			for (i__ = 1; i__ <= i__3; ++i__) {
#line 326 "ctrrfs.f"
			    i__4 = i__ + k * a_dim1;
#line 326 "ctrrfs.f"
			    rwork[i__] += ((d__1 = a[i__4].r, abs(d__1)) + (
				    d__2 = d_imag(&a[i__ + k * a_dim1]), abs(
				    d__2))) * xk;
#line 327 "ctrrfs.f"
/* L30: */
#line 327 "ctrrfs.f"
			}
#line 328 "ctrrfs.f"
/* L40: */
#line 328 "ctrrfs.f"
		    }
#line 329 "ctrrfs.f"
		} else {
#line 330 "ctrrfs.f"
		    i__2 = *n;
#line 330 "ctrrfs.f"
		    for (k = 1; k <= i__2; ++k) {
#line 331 "ctrrfs.f"
			i__3 = k + j * x_dim1;
#line 331 "ctrrfs.f"
			xk = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = d_imag(&
				x[k + j * x_dim1]), abs(d__2));
#line 332 "ctrrfs.f"
			i__3 = k - 1;
#line 332 "ctrrfs.f"
			for (i__ = 1; i__ <= i__3; ++i__) {
#line 333 "ctrrfs.f"
			    i__4 = i__ + k * a_dim1;
#line 333 "ctrrfs.f"
			    rwork[i__] += ((d__1 = a[i__4].r, abs(d__1)) + (
				    d__2 = d_imag(&a[i__ + k * a_dim1]), abs(
				    d__2))) * xk;
#line 334 "ctrrfs.f"
/* L50: */
#line 334 "ctrrfs.f"
			}
#line 335 "ctrrfs.f"
			rwork[k] += xk;
#line 336 "ctrrfs.f"
/* L60: */
#line 336 "ctrrfs.f"
		    }
#line 337 "ctrrfs.f"
		}
#line 338 "ctrrfs.f"
	    } else {
#line 339 "ctrrfs.f"
		if (nounit) {
#line 340 "ctrrfs.f"
		    i__2 = *n;
#line 340 "ctrrfs.f"
		    for (k = 1; k <= i__2; ++k) {
#line 341 "ctrrfs.f"
			i__3 = k + j * x_dim1;
#line 341 "ctrrfs.f"
			xk = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = d_imag(&
				x[k + j * x_dim1]), abs(d__2));
#line 342 "ctrrfs.f"
			i__3 = *n;
#line 342 "ctrrfs.f"
			for (i__ = k; i__ <= i__3; ++i__) {
#line 343 "ctrrfs.f"
			    i__4 = i__ + k * a_dim1;
#line 343 "ctrrfs.f"
			    rwork[i__] += ((d__1 = a[i__4].r, abs(d__1)) + (
				    d__2 = d_imag(&a[i__ + k * a_dim1]), abs(
				    d__2))) * xk;
#line 344 "ctrrfs.f"
/* L70: */
#line 344 "ctrrfs.f"
			}
#line 345 "ctrrfs.f"
/* L80: */
#line 345 "ctrrfs.f"
		    }
#line 346 "ctrrfs.f"
		} else {
#line 347 "ctrrfs.f"
		    i__2 = *n;
#line 347 "ctrrfs.f"
		    for (k = 1; k <= i__2; ++k) {
#line 348 "ctrrfs.f"
			i__3 = k + j * x_dim1;
#line 348 "ctrrfs.f"
			xk = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = d_imag(&
				x[k + j * x_dim1]), abs(d__2));
#line 349 "ctrrfs.f"
			i__3 = *n;
#line 349 "ctrrfs.f"
			for (i__ = k + 1; i__ <= i__3; ++i__) {
#line 350 "ctrrfs.f"
			    i__4 = i__ + k * a_dim1;
#line 350 "ctrrfs.f"
			    rwork[i__] += ((d__1 = a[i__4].r, abs(d__1)) + (
				    d__2 = d_imag(&a[i__ + k * a_dim1]), abs(
				    d__2))) * xk;
#line 351 "ctrrfs.f"
/* L90: */
#line 351 "ctrrfs.f"
			}
#line 352 "ctrrfs.f"
			rwork[k] += xk;
#line 353 "ctrrfs.f"
/* L100: */
#line 353 "ctrrfs.f"
		    }
#line 354 "ctrrfs.f"
		}
#line 355 "ctrrfs.f"
	    }
#line 356 "ctrrfs.f"
	} else {

/*           Compute abs(A**H)*abs(X) + abs(B). */

#line 360 "ctrrfs.f"
	    if (upper) {
#line 361 "ctrrfs.f"
		if (nounit) {
#line 362 "ctrrfs.f"
		    i__2 = *n;
#line 362 "ctrrfs.f"
		    for (k = 1; k <= i__2; ++k) {
#line 363 "ctrrfs.f"
			s = 0.;
#line 364 "ctrrfs.f"
			i__3 = k;
#line 364 "ctrrfs.f"
			for (i__ = 1; i__ <= i__3; ++i__) {
#line 365 "ctrrfs.f"
			    i__4 = i__ + k * a_dim1;
#line 365 "ctrrfs.f"
			    i__5 = i__ + j * x_dim1;
#line 365 "ctrrfs.f"
			    s += ((d__1 = a[i__4].r, abs(d__1)) + (d__2 = 
				    d_imag(&a[i__ + k * a_dim1]), abs(d__2))) 
				    * ((d__3 = x[i__5].r, abs(d__3)) + (d__4 =
				     d_imag(&x[i__ + j * x_dim1]), abs(d__4)))
				    ;
#line 366 "ctrrfs.f"
/* L110: */
#line 366 "ctrrfs.f"
			}
#line 367 "ctrrfs.f"
			rwork[k] += s;
#line 368 "ctrrfs.f"
/* L120: */
#line 368 "ctrrfs.f"
		    }
#line 369 "ctrrfs.f"
		} else {
#line 370 "ctrrfs.f"
		    i__2 = *n;
#line 370 "ctrrfs.f"
		    for (k = 1; k <= i__2; ++k) {
#line 371 "ctrrfs.f"
			i__3 = k + j * x_dim1;
#line 371 "ctrrfs.f"
			s = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = d_imag(&x[
				k + j * x_dim1]), abs(d__2));
#line 372 "ctrrfs.f"
			i__3 = k - 1;
#line 372 "ctrrfs.f"
			for (i__ = 1; i__ <= i__3; ++i__) {
#line 373 "ctrrfs.f"
			    i__4 = i__ + k * a_dim1;
#line 373 "ctrrfs.f"
			    i__5 = i__ + j * x_dim1;
#line 373 "ctrrfs.f"
			    s += ((d__1 = a[i__4].r, abs(d__1)) + (d__2 = 
				    d_imag(&a[i__ + k * a_dim1]), abs(d__2))) 
				    * ((d__3 = x[i__5].r, abs(d__3)) + (d__4 =
				     d_imag(&x[i__ + j * x_dim1]), abs(d__4)))
				    ;
#line 374 "ctrrfs.f"
/* L130: */
#line 374 "ctrrfs.f"
			}
#line 375 "ctrrfs.f"
			rwork[k] += s;
#line 376 "ctrrfs.f"
/* L140: */
#line 376 "ctrrfs.f"
		    }
#line 377 "ctrrfs.f"
		}
#line 378 "ctrrfs.f"
	    } else {
#line 379 "ctrrfs.f"
		if (nounit) {
#line 380 "ctrrfs.f"
		    i__2 = *n;
#line 380 "ctrrfs.f"
		    for (k = 1; k <= i__2; ++k) {
#line 381 "ctrrfs.f"
			s = 0.;
#line 382 "ctrrfs.f"
			i__3 = *n;
#line 382 "ctrrfs.f"
			for (i__ = k; i__ <= i__3; ++i__) {
#line 383 "ctrrfs.f"
			    i__4 = i__ + k * a_dim1;
#line 383 "ctrrfs.f"
			    i__5 = i__ + j * x_dim1;
#line 383 "ctrrfs.f"
			    s += ((d__1 = a[i__4].r, abs(d__1)) + (d__2 = 
				    d_imag(&a[i__ + k * a_dim1]), abs(d__2))) 
				    * ((d__3 = x[i__5].r, abs(d__3)) + (d__4 =
				     d_imag(&x[i__ + j * x_dim1]), abs(d__4)))
				    ;
#line 384 "ctrrfs.f"
/* L150: */
#line 384 "ctrrfs.f"
			}
#line 385 "ctrrfs.f"
			rwork[k] += s;
#line 386 "ctrrfs.f"
/* L160: */
#line 386 "ctrrfs.f"
		    }
#line 387 "ctrrfs.f"
		} else {
#line 388 "ctrrfs.f"
		    i__2 = *n;
#line 388 "ctrrfs.f"
		    for (k = 1; k <= i__2; ++k) {
#line 389 "ctrrfs.f"
			i__3 = k + j * x_dim1;
#line 389 "ctrrfs.f"
			s = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = d_imag(&x[
				k + j * x_dim1]), abs(d__2));
#line 390 "ctrrfs.f"
			i__3 = *n;
#line 390 "ctrrfs.f"
			for (i__ = k + 1; i__ <= i__3; ++i__) {
#line 391 "ctrrfs.f"
			    i__4 = i__ + k * a_dim1;
#line 391 "ctrrfs.f"
			    i__5 = i__ + j * x_dim1;
#line 391 "ctrrfs.f"
			    s += ((d__1 = a[i__4].r, abs(d__1)) + (d__2 = 
				    d_imag(&a[i__ + k * a_dim1]), abs(d__2))) 
				    * ((d__3 = x[i__5].r, abs(d__3)) + (d__4 =
				     d_imag(&x[i__ + j * x_dim1]), abs(d__4)))
				    ;
#line 392 "ctrrfs.f"
/* L170: */
#line 392 "ctrrfs.f"
			}
#line 393 "ctrrfs.f"
			rwork[k] += s;
#line 394 "ctrrfs.f"
/* L180: */
#line 394 "ctrrfs.f"
		    }
#line 395 "ctrrfs.f"
		}
#line 396 "ctrrfs.f"
	    }
#line 397 "ctrrfs.f"
	}
#line 398 "ctrrfs.f"
	s = 0.;
#line 399 "ctrrfs.f"
	i__2 = *n;
#line 399 "ctrrfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 400 "ctrrfs.f"
	    if (rwork[i__] > safe2) {
/* Computing MAX */
#line 401 "ctrrfs.f"
		i__3 = i__;
#line 401 "ctrrfs.f"
		d__3 = s, d__4 = ((d__1 = work[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&work[i__]), abs(d__2))) / rwork[i__];
#line 401 "ctrrfs.f"
		s = max(d__3,d__4);
#line 402 "ctrrfs.f"
	    } else {
/* Computing MAX */
#line 403 "ctrrfs.f"
		i__3 = i__;
#line 403 "ctrrfs.f"
		d__3 = s, d__4 = ((d__1 = work[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&work[i__]), abs(d__2)) + safe1) / (rwork[i__] 
			+ safe1);
#line 403 "ctrrfs.f"
		s = max(d__3,d__4);
#line 405 "ctrrfs.f"
	    }
#line 406 "ctrrfs.f"
/* L190: */
#line 406 "ctrrfs.f"
	}
#line 407 "ctrrfs.f"
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

#line 431 "ctrrfs.f"
	i__2 = *n;
#line 431 "ctrrfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 432 "ctrrfs.f"
	    if (rwork[i__] > safe2) {
#line 433 "ctrrfs.f"
		i__3 = i__;
#line 433 "ctrrfs.f"
		rwork[i__] = (d__1 = work[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&work[i__]), abs(d__2)) + nz * eps * rwork[i__]
			;
#line 434 "ctrrfs.f"
	    } else {
#line 435 "ctrrfs.f"
		i__3 = i__;
#line 435 "ctrrfs.f"
		rwork[i__] = (d__1 = work[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&work[i__]), abs(d__2)) + nz * eps * rwork[i__]
			 + safe1;
#line 437 "ctrrfs.f"
	    }
#line 438 "ctrrfs.f"
/* L200: */
#line 438 "ctrrfs.f"
	}

#line 440 "ctrrfs.f"
	kase = 0;
#line 441 "ctrrfs.f"
L210:
#line 442 "ctrrfs.f"
	clacn2_(n, &work[*n + 1], &work[1], &ferr[j], &kase, isave);
#line 443 "ctrrfs.f"
	if (kase != 0) {
#line 444 "ctrrfs.f"
	    if (kase == 1) {

/*              Multiply by diag(W)*inv(op(A)**H). */

#line 448 "ctrrfs.f"
		ctrsv_(uplo, transt, diag, n, &a[a_offset], lda, &work[1], &
			c__1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 449 "ctrrfs.f"
		i__2 = *n;
#line 449 "ctrrfs.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 450 "ctrrfs.f"
		    i__3 = i__;
#line 450 "ctrrfs.f"
		    i__4 = i__;
#line 450 "ctrrfs.f"
		    i__5 = i__;
#line 450 "ctrrfs.f"
		    z__1.r = rwork[i__4] * work[i__5].r, z__1.i = rwork[i__4] 
			    * work[i__5].i;
#line 450 "ctrrfs.f"
		    work[i__3].r = z__1.r, work[i__3].i = z__1.i;
#line 451 "ctrrfs.f"
/* L220: */
#line 451 "ctrrfs.f"
		}
#line 452 "ctrrfs.f"
	    } else {

/*              Multiply by inv(op(A))*diag(W). */

#line 456 "ctrrfs.f"
		i__2 = *n;
#line 456 "ctrrfs.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 457 "ctrrfs.f"
		    i__3 = i__;
#line 457 "ctrrfs.f"
		    i__4 = i__;
#line 457 "ctrrfs.f"
		    i__5 = i__;
#line 457 "ctrrfs.f"
		    z__1.r = rwork[i__4] * work[i__5].r, z__1.i = rwork[i__4] 
			    * work[i__5].i;
#line 457 "ctrrfs.f"
		    work[i__3].r = z__1.r, work[i__3].i = z__1.i;
#line 458 "ctrrfs.f"
/* L230: */
#line 458 "ctrrfs.f"
		}
#line 459 "ctrrfs.f"
		ctrsv_(uplo, transn, diag, n, &a[a_offset], lda, &work[1], &
			c__1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 460 "ctrrfs.f"
	    }
#line 461 "ctrrfs.f"
	    goto L210;
#line 462 "ctrrfs.f"
	}

/*        Normalize error. */

#line 466 "ctrrfs.f"
	lstres = 0.;
#line 467 "ctrrfs.f"
	i__2 = *n;
#line 467 "ctrrfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 468 "ctrrfs.f"
	    i__3 = i__ + j * x_dim1;
#line 468 "ctrrfs.f"
	    d__3 = lstres, d__4 = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = 
		    d_imag(&x[i__ + j * x_dim1]), abs(d__2));
#line 468 "ctrrfs.f"
	    lstres = max(d__3,d__4);
#line 469 "ctrrfs.f"
/* L240: */
#line 469 "ctrrfs.f"
	}
#line 470 "ctrrfs.f"
	if (lstres != 0.) {
#line 470 "ctrrfs.f"
	    ferr[j] /= lstres;
#line 470 "ctrrfs.f"
	}

#line 473 "ctrrfs.f"
/* L250: */
#line 473 "ctrrfs.f"
    }

#line 475 "ctrrfs.f"
    return 0;

/*     End of CTRRFS */

} /* ctrrfs_ */

