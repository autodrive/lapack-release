#line 1 "zherfs.f"
/* zherfs.f -- translated by f2c (version 20100827).
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

#line 1 "zherfs.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static integer c__1 = 1;

/* > \brief \b ZHERFS */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZHERFS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zherfs.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zherfs.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zherfs.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZHERFS( UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV, B, LDB, */
/*                          X, LDX, FERR, BERR, WORK, RWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDA, LDAF, LDB, LDX, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       DOUBLE PRECISION   BERR( * ), FERR( * ), RWORK( * ) */
/*       COMPLEX*16         A( LDA, * ), AF( LDAF, * ), B( LDB, * ), */
/*      $                   WORK( * ), X( LDX, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZHERFS improves the computed solution to a system of linear */
/* > equations when the coefficient matrix is Hermitian indefinite, and */
/* > provides error bounds and backward error estimates for the solution. */
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
/* > \param[in] NRHS */
/* > \verbatim */
/* >          NRHS is INTEGER */
/* >          The number of right hand sides, i.e., the number of columns */
/* >          of the matrices B and X.  NRHS >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
/* >          The Hermitian matrix A.  If UPLO = 'U', the leading N-by-N */
/* >          upper triangular part of A contains the upper triangular part */
/* >          of the matrix A, and the strictly lower triangular part of A */
/* >          is not referenced.  If UPLO = 'L', the leading N-by-N lower */
/* >          triangular part of A contains the lower triangular part of */
/* >          the matrix A, and the strictly upper triangular part of A is */
/* >          not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] AF */
/* > \verbatim */
/* >          AF is COMPLEX*16 array, dimension (LDAF,N) */
/* >          The factored form of the matrix A.  AF contains the block */
/* >          diagonal matrix D and the multipliers used to obtain the */
/* >          factor U or L from the factorization A = U*D*U**H or */
/* >          A = L*D*L**H as computed by ZHETRF. */
/* > \endverbatim */
/* > */
/* > \param[in] LDAF */
/* > \verbatim */
/* >          LDAF is INTEGER */
/* >          The leading dimension of the array AF.  LDAF >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N) */
/* >          Details of the interchanges and the block structure of D */
/* >          as determined by ZHETRF. */
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
/* > \param[in,out] X */
/* > \verbatim */
/* >          X is COMPLEX*16 array, dimension (LDX,NRHS) */
/* >          On entry, the solution matrix X, as computed by ZHETRS. */
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

/* > \ingroup complex16HEcomputational */

/*  ===================================================================== */
/* Subroutine */ int zherfs_(char *uplo, integer *n, integer *nrhs, 
	doublecomplex *a, integer *lda, doublecomplex *af, integer *ldaf, 
	integer *ipiv, doublecomplex *b, integer *ldb, doublecomplex *x, 
	integer *ldx, doublereal *ferr, doublereal *berr, doublecomplex *work,
	 doublereal *rwork, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, af_dim1, af_offset, b_dim1, b_offset, x_dim1, 
	    x_offset, i__1, i__2, i__3, i__4, i__5;
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
    static integer isave[3], count;
    extern /* Subroutine */ int zhemv_(char *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, ftnlen);
    static logical upper;
    extern /* Subroutine */ int zcopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), zaxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *), zlacn2_(
	    integer *, doublecomplex *, doublecomplex *, doublereal *, 
	    integer *, integer *);
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal lstres;
    extern /* Subroutine */ int zhetrs_(char *, integer *, integer *, 
	    doublecomplex *, integer *, integer *, doublecomplex *, integer *,
	     integer *, ftnlen);


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

#line 255 "zherfs.f"
    /* Parameter adjustments */
#line 255 "zherfs.f"
    a_dim1 = *lda;
#line 255 "zherfs.f"
    a_offset = 1 + a_dim1;
#line 255 "zherfs.f"
    a -= a_offset;
#line 255 "zherfs.f"
    af_dim1 = *ldaf;
#line 255 "zherfs.f"
    af_offset = 1 + af_dim1;
#line 255 "zherfs.f"
    af -= af_offset;
#line 255 "zherfs.f"
    --ipiv;
#line 255 "zherfs.f"
    b_dim1 = *ldb;
#line 255 "zherfs.f"
    b_offset = 1 + b_dim1;
#line 255 "zherfs.f"
    b -= b_offset;
#line 255 "zherfs.f"
    x_dim1 = *ldx;
#line 255 "zherfs.f"
    x_offset = 1 + x_dim1;
#line 255 "zherfs.f"
    x -= x_offset;
#line 255 "zherfs.f"
    --ferr;
#line 255 "zherfs.f"
    --berr;
#line 255 "zherfs.f"
    --work;
#line 255 "zherfs.f"
    --rwork;
#line 255 "zherfs.f"

#line 255 "zherfs.f"
    /* Function Body */
#line 255 "zherfs.f"
    *info = 0;
#line 256 "zherfs.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 257 "zherfs.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 258 "zherfs.f"
	*info = -1;
#line 259 "zherfs.f"
    } else if (*n < 0) {
#line 260 "zherfs.f"
	*info = -2;
#line 261 "zherfs.f"
    } else if (*nrhs < 0) {
#line 262 "zherfs.f"
	*info = -3;
#line 263 "zherfs.f"
    } else if (*lda < max(1,*n)) {
#line 264 "zherfs.f"
	*info = -5;
#line 265 "zherfs.f"
    } else if (*ldaf < max(1,*n)) {
#line 266 "zherfs.f"
	*info = -7;
#line 267 "zherfs.f"
    } else if (*ldb < max(1,*n)) {
#line 268 "zherfs.f"
	*info = -10;
#line 269 "zherfs.f"
    } else if (*ldx < max(1,*n)) {
#line 270 "zherfs.f"
	*info = -12;
#line 271 "zherfs.f"
    }
#line 272 "zherfs.f"
    if (*info != 0) {
#line 273 "zherfs.f"
	i__1 = -(*info);
#line 273 "zherfs.f"
	xerbla_("ZHERFS", &i__1, (ftnlen)6);
#line 274 "zherfs.f"
	return 0;
#line 275 "zherfs.f"
    }

/*     Quick return if possible */

#line 279 "zherfs.f"
    if (*n == 0 || *nrhs == 0) {
#line 280 "zherfs.f"
	i__1 = *nrhs;
#line 280 "zherfs.f"
	for (j = 1; j <= i__1; ++j) {
#line 281 "zherfs.f"
	    ferr[j] = 0.;
#line 282 "zherfs.f"
	    berr[j] = 0.;
#line 283 "zherfs.f"
/* L10: */
#line 283 "zherfs.f"
	}
#line 284 "zherfs.f"
	return 0;
#line 285 "zherfs.f"
    }

/*     NZ = maximum number of nonzero elements in each row of A, plus 1 */

#line 289 "zherfs.f"
    nz = *n + 1;
#line 290 "zherfs.f"
    eps = dlamch_("Epsilon", (ftnlen)7);
#line 291 "zherfs.f"
    safmin = dlamch_("Safe minimum", (ftnlen)12);
#line 292 "zherfs.f"
    safe1 = nz * safmin;
#line 293 "zherfs.f"
    safe2 = safe1 / eps;

/*     Do for each right hand side */

#line 297 "zherfs.f"
    i__1 = *nrhs;
#line 297 "zherfs.f"
    for (j = 1; j <= i__1; ++j) {

#line 299 "zherfs.f"
	count = 1;
#line 300 "zherfs.f"
	lstres = 3.;
#line 301 "zherfs.f"
L20:

/*        Loop until stopping criterion is satisfied. */

/*        Compute residual R = B - A * X */

#line 307 "zherfs.f"
	zcopy_(n, &b[j * b_dim1 + 1], &c__1, &work[1], &c__1);
#line 308 "zherfs.f"
	z__1.r = -1., z__1.i = -0.;
#line 308 "zherfs.f"
	zhemv_(uplo, n, &z__1, &a[a_offset], lda, &x[j * x_dim1 + 1], &c__1, &
		c_b1, &work[1], &c__1, (ftnlen)1);

/*        Compute componentwise relative backward error from formula */

/*        max(i) ( abs(R(i)) / ( abs(A)*abs(X) + abs(B) )(i) ) */

/*        where abs(Z) is the componentwise absolute value of the matrix */
/*        or vector Z.  If the i-th component of the denominator is less */
/*        than SAFE2, then SAFE1 is added to the i-th components of the */
/*        numerator and denominator before dividing. */

#line 319 "zherfs.f"
	i__2 = *n;
#line 319 "zherfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 320 "zherfs.f"
	    i__3 = i__ + j * b_dim1;
#line 320 "zherfs.f"
	    rwork[i__] = (d__1 = b[i__3].r, abs(d__1)) + (d__2 = d_imag(&b[
		    i__ + j * b_dim1]), abs(d__2));
#line 321 "zherfs.f"
/* L30: */
#line 321 "zherfs.f"
	}

/*        Compute abs(A)*abs(X) + abs(B). */

#line 325 "zherfs.f"
	if (upper) {
#line 326 "zherfs.f"
	    i__2 = *n;
#line 326 "zherfs.f"
	    for (k = 1; k <= i__2; ++k) {
#line 327 "zherfs.f"
		s = 0.;
#line 328 "zherfs.f"
		i__3 = k + j * x_dim1;
#line 328 "zherfs.f"
		xk = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = d_imag(&x[k + j *
			 x_dim1]), abs(d__2));
#line 329 "zherfs.f"
		i__3 = k - 1;
#line 329 "zherfs.f"
		for (i__ = 1; i__ <= i__3; ++i__) {
#line 330 "zherfs.f"
		    i__4 = i__ + k * a_dim1;
#line 330 "zherfs.f"
		    rwork[i__] += ((d__1 = a[i__4].r, abs(d__1)) + (d__2 = 
			    d_imag(&a[i__ + k * a_dim1]), abs(d__2))) * xk;
#line 331 "zherfs.f"
		    i__4 = i__ + k * a_dim1;
#line 331 "zherfs.f"
		    i__5 = i__ + j * x_dim1;
#line 331 "zherfs.f"
		    s += ((d__1 = a[i__4].r, abs(d__1)) + (d__2 = d_imag(&a[
			    i__ + k * a_dim1]), abs(d__2))) * ((d__3 = x[i__5]
			    .r, abs(d__3)) + (d__4 = d_imag(&x[i__ + j * 
			    x_dim1]), abs(d__4)));
#line 332 "zherfs.f"
/* L40: */
#line 332 "zherfs.f"
		}
#line 333 "zherfs.f"
		i__3 = k + k * a_dim1;
#line 333 "zherfs.f"
		rwork[k] = rwork[k] + (d__1 = a[i__3].r, abs(d__1)) * xk + s;
#line 334 "zherfs.f"
/* L50: */
#line 334 "zherfs.f"
	    }
#line 335 "zherfs.f"
	} else {
#line 336 "zherfs.f"
	    i__2 = *n;
#line 336 "zherfs.f"
	    for (k = 1; k <= i__2; ++k) {
#line 337 "zherfs.f"
		s = 0.;
#line 338 "zherfs.f"
		i__3 = k + j * x_dim1;
#line 338 "zherfs.f"
		xk = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = d_imag(&x[k + j *
			 x_dim1]), abs(d__2));
#line 339 "zherfs.f"
		i__3 = k + k * a_dim1;
#line 339 "zherfs.f"
		rwork[k] += (d__1 = a[i__3].r, abs(d__1)) * xk;
#line 340 "zherfs.f"
		i__3 = *n;
#line 340 "zherfs.f"
		for (i__ = k + 1; i__ <= i__3; ++i__) {
#line 341 "zherfs.f"
		    i__4 = i__ + k * a_dim1;
#line 341 "zherfs.f"
		    rwork[i__] += ((d__1 = a[i__4].r, abs(d__1)) + (d__2 = 
			    d_imag(&a[i__ + k * a_dim1]), abs(d__2))) * xk;
#line 342 "zherfs.f"
		    i__4 = i__ + k * a_dim1;
#line 342 "zherfs.f"
		    i__5 = i__ + j * x_dim1;
#line 342 "zherfs.f"
		    s += ((d__1 = a[i__4].r, abs(d__1)) + (d__2 = d_imag(&a[
			    i__ + k * a_dim1]), abs(d__2))) * ((d__3 = x[i__5]
			    .r, abs(d__3)) + (d__4 = d_imag(&x[i__ + j * 
			    x_dim1]), abs(d__4)));
#line 343 "zherfs.f"
/* L60: */
#line 343 "zherfs.f"
		}
#line 344 "zherfs.f"
		rwork[k] += s;
#line 345 "zherfs.f"
/* L70: */
#line 345 "zherfs.f"
	    }
#line 346 "zherfs.f"
	}
#line 347 "zherfs.f"
	s = 0.;
#line 348 "zherfs.f"
	i__2 = *n;
#line 348 "zherfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 349 "zherfs.f"
	    if (rwork[i__] > safe2) {
/* Computing MAX */
#line 350 "zherfs.f"
		i__3 = i__;
#line 350 "zherfs.f"
		d__3 = s, d__4 = ((d__1 = work[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&work[i__]), abs(d__2))) / rwork[i__];
#line 350 "zherfs.f"
		s = max(d__3,d__4);
#line 351 "zherfs.f"
	    } else {
/* Computing MAX */
#line 352 "zherfs.f"
		i__3 = i__;
#line 352 "zherfs.f"
		d__3 = s, d__4 = ((d__1 = work[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&work[i__]), abs(d__2)) + safe1) / (rwork[i__] 
			+ safe1);
#line 352 "zherfs.f"
		s = max(d__3,d__4);
#line 354 "zherfs.f"
	    }
#line 355 "zherfs.f"
/* L80: */
#line 355 "zherfs.f"
	}
#line 356 "zherfs.f"
	berr[j] = s;

/*        Test stopping criterion. Continue iterating if */
/*           1) The residual BERR(J) is larger than machine epsilon, and */
/*           2) BERR(J) decreased by at least a factor of 2 during the */
/*              last iteration, and */
/*           3) At most ITMAX iterations tried. */

#line 364 "zherfs.f"
	if (berr[j] > eps && berr[j] * 2. <= lstres && count <= 5) {

/*           Update solution and try again. */

#line 369 "zherfs.f"
	    zhetrs_(uplo, n, &c__1, &af[af_offset], ldaf, &ipiv[1], &work[1], 
		    n, info, (ftnlen)1);
#line 370 "zherfs.f"
	    zaxpy_(n, &c_b1, &work[1], &c__1, &x[j * x_dim1 + 1], &c__1);
#line 371 "zherfs.f"
	    lstres = berr[j];
#line 372 "zherfs.f"
	    ++count;
#line 373 "zherfs.f"
	    goto L20;
#line 374 "zherfs.f"
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

/*        Use ZLACN2 to estimate the infinity-norm of the matrix */
/*           inv(A) * diag(W), */
/*        where W = abs(R) + NZ*EPS*( abs(A)*abs(X)+abs(B) ))) */

#line 398 "zherfs.f"
	i__2 = *n;
#line 398 "zherfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 399 "zherfs.f"
	    if (rwork[i__] > safe2) {
#line 400 "zherfs.f"
		i__3 = i__;
#line 400 "zherfs.f"
		rwork[i__] = (d__1 = work[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&work[i__]), abs(d__2)) + nz * eps * rwork[i__]
			;
#line 401 "zherfs.f"
	    } else {
#line 402 "zherfs.f"
		i__3 = i__;
#line 402 "zherfs.f"
		rwork[i__] = (d__1 = work[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&work[i__]), abs(d__2)) + nz * eps * rwork[i__]
			 + safe1;
#line 404 "zherfs.f"
	    }
#line 405 "zherfs.f"
/* L90: */
#line 405 "zherfs.f"
	}

#line 407 "zherfs.f"
	kase = 0;
#line 408 "zherfs.f"
L100:
#line 409 "zherfs.f"
	zlacn2_(n, &work[*n + 1], &work[1], &ferr[j], &kase, isave);
#line 410 "zherfs.f"
	if (kase != 0) {
#line 411 "zherfs.f"
	    if (kase == 1) {

/*              Multiply by diag(W)*inv(A**H). */

#line 415 "zherfs.f"
		zhetrs_(uplo, n, &c__1, &af[af_offset], ldaf, &ipiv[1], &work[
			1], n, info, (ftnlen)1);
#line 416 "zherfs.f"
		i__2 = *n;
#line 416 "zherfs.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 417 "zherfs.f"
		    i__3 = i__;
#line 417 "zherfs.f"
		    i__4 = i__;
#line 417 "zherfs.f"
		    i__5 = i__;
#line 417 "zherfs.f"
		    z__1.r = rwork[i__4] * work[i__5].r, z__1.i = rwork[i__4] 
			    * work[i__5].i;
#line 417 "zherfs.f"
		    work[i__3].r = z__1.r, work[i__3].i = z__1.i;
#line 418 "zherfs.f"
/* L110: */
#line 418 "zherfs.f"
		}
#line 419 "zherfs.f"
	    } else if (kase == 2) {

/*              Multiply by inv(A)*diag(W). */

#line 423 "zherfs.f"
		i__2 = *n;
#line 423 "zherfs.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 424 "zherfs.f"
		    i__3 = i__;
#line 424 "zherfs.f"
		    i__4 = i__;
#line 424 "zherfs.f"
		    i__5 = i__;
#line 424 "zherfs.f"
		    z__1.r = rwork[i__4] * work[i__5].r, z__1.i = rwork[i__4] 
			    * work[i__5].i;
#line 424 "zherfs.f"
		    work[i__3].r = z__1.r, work[i__3].i = z__1.i;
#line 425 "zherfs.f"
/* L120: */
#line 425 "zherfs.f"
		}
#line 426 "zherfs.f"
		zhetrs_(uplo, n, &c__1, &af[af_offset], ldaf, &ipiv[1], &work[
			1], n, info, (ftnlen)1);
#line 427 "zherfs.f"
	    }
#line 428 "zherfs.f"
	    goto L100;
#line 429 "zherfs.f"
	}

/*        Normalize error. */

#line 433 "zherfs.f"
	lstres = 0.;
#line 434 "zherfs.f"
	i__2 = *n;
#line 434 "zherfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 435 "zherfs.f"
	    i__3 = i__ + j * x_dim1;
#line 435 "zherfs.f"
	    d__3 = lstres, d__4 = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = 
		    d_imag(&x[i__ + j * x_dim1]), abs(d__2));
#line 435 "zherfs.f"
	    lstres = max(d__3,d__4);
#line 436 "zherfs.f"
/* L130: */
#line 436 "zherfs.f"
	}
#line 437 "zherfs.f"
	if (lstres != 0.) {
#line 437 "zherfs.f"
	    ferr[j] /= lstres;
#line 437 "zherfs.f"
	}

#line 440 "zherfs.f"
/* L140: */
#line 440 "zherfs.f"
    }

#line 442 "zherfs.f"
    return 0;

/*     End of ZHERFS */

} /* zherfs_ */

