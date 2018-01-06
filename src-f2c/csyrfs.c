#line 1 "csyrfs.f"
/* csyrfs.f -- translated by f2c (version 20100827).
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

#line 1 "csyrfs.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static integer c__1 = 1;

/* > \brief \b CSYRFS */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CSYRFS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/csyrfs.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/csyrfs.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/csyrfs.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CSYRFS( UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV, B, LDB, */
/*                          X, LDX, FERR, BERR, WORK, RWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDA, LDAF, LDB, LDX, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       REAL               BERR( * ), FERR( * ), RWORK( * ) */
/*       COMPLEX            A( LDA, * ), AF( LDAF, * ), B( LDB, * ), */
/*      $                   WORK( * ), X( LDX, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CSYRFS improves the computed solution to a system of linear */
/* > equations when the coefficient matrix is symmetric indefinite, and */
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
/* >          A is COMPLEX array, dimension (LDA,N) */
/* >          The symmetric matrix A.  If UPLO = 'U', the leading N-by-N */
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
/* >          AF is COMPLEX array, dimension (LDAF,N) */
/* >          The factored form of the matrix A.  AF contains the block */
/* >          diagonal matrix D and the multipliers used to obtain the */
/* >          factor U or L from the factorization A = U*D*U**T or */
/* >          A = L*D*L**T as computed by CSYTRF. */
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
/* >          as determined by CSYTRF. */
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
/* > \param[in,out] X */
/* > \verbatim */
/* >          X is COMPLEX array, dimension (LDX,NRHS) */
/* >          On entry, the solution matrix X, as computed by CSYTRS. */
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

/* > \date December 2016 */

/* > \ingroup complexSYcomputational */

/*  ===================================================================== */
/* Subroutine */ int csyrfs_(char *uplo, integer *n, integer *nrhs, 
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
    static integer isave[3];
    extern /* Subroutine */ int ccopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), caxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    static integer count;
    static logical upper;
    extern /* Subroutine */ int csymv_(char *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, ftnlen), clacn2_(
	    integer *, doublecomplex *, doublecomplex *, doublereal *, 
	    integer *, integer *);
    extern doublereal slamch_(char *, ftnlen);
    static doublereal safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal lstres;
    extern /* Subroutine */ int csytrs_(char *, integer *, integer *, 
	    doublecomplex *, integer *, integer *, doublecomplex *, integer *,
	     integer *, ftnlen);


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

#line 255 "csyrfs.f"
    /* Parameter adjustments */
#line 255 "csyrfs.f"
    a_dim1 = *lda;
#line 255 "csyrfs.f"
    a_offset = 1 + a_dim1;
#line 255 "csyrfs.f"
    a -= a_offset;
#line 255 "csyrfs.f"
    af_dim1 = *ldaf;
#line 255 "csyrfs.f"
    af_offset = 1 + af_dim1;
#line 255 "csyrfs.f"
    af -= af_offset;
#line 255 "csyrfs.f"
    --ipiv;
#line 255 "csyrfs.f"
    b_dim1 = *ldb;
#line 255 "csyrfs.f"
    b_offset = 1 + b_dim1;
#line 255 "csyrfs.f"
    b -= b_offset;
#line 255 "csyrfs.f"
    x_dim1 = *ldx;
#line 255 "csyrfs.f"
    x_offset = 1 + x_dim1;
#line 255 "csyrfs.f"
    x -= x_offset;
#line 255 "csyrfs.f"
    --ferr;
#line 255 "csyrfs.f"
    --berr;
#line 255 "csyrfs.f"
    --work;
#line 255 "csyrfs.f"
    --rwork;
#line 255 "csyrfs.f"

#line 255 "csyrfs.f"
    /* Function Body */
#line 255 "csyrfs.f"
    *info = 0;
#line 256 "csyrfs.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 257 "csyrfs.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 258 "csyrfs.f"
	*info = -1;
#line 259 "csyrfs.f"
    } else if (*n < 0) {
#line 260 "csyrfs.f"
	*info = -2;
#line 261 "csyrfs.f"
    } else if (*nrhs < 0) {
#line 262 "csyrfs.f"
	*info = -3;
#line 263 "csyrfs.f"
    } else if (*lda < max(1,*n)) {
#line 264 "csyrfs.f"
	*info = -5;
#line 265 "csyrfs.f"
    } else if (*ldaf < max(1,*n)) {
#line 266 "csyrfs.f"
	*info = -7;
#line 267 "csyrfs.f"
    } else if (*ldb < max(1,*n)) {
#line 268 "csyrfs.f"
	*info = -10;
#line 269 "csyrfs.f"
    } else if (*ldx < max(1,*n)) {
#line 270 "csyrfs.f"
	*info = -12;
#line 271 "csyrfs.f"
    }
#line 272 "csyrfs.f"
    if (*info != 0) {
#line 273 "csyrfs.f"
	i__1 = -(*info);
#line 273 "csyrfs.f"
	xerbla_("CSYRFS", &i__1, (ftnlen)6);
#line 274 "csyrfs.f"
	return 0;
#line 275 "csyrfs.f"
    }

/*     Quick return if possible */

#line 279 "csyrfs.f"
    if (*n == 0 || *nrhs == 0) {
#line 280 "csyrfs.f"
	i__1 = *nrhs;
#line 280 "csyrfs.f"
	for (j = 1; j <= i__1; ++j) {
#line 281 "csyrfs.f"
	    ferr[j] = 0.;
#line 282 "csyrfs.f"
	    berr[j] = 0.;
#line 283 "csyrfs.f"
/* L10: */
#line 283 "csyrfs.f"
	}
#line 284 "csyrfs.f"
	return 0;
#line 285 "csyrfs.f"
    }

/*     NZ = maximum number of nonzero elements in each row of A, plus 1 */

#line 289 "csyrfs.f"
    nz = *n + 1;
#line 290 "csyrfs.f"
    eps = slamch_("Epsilon", (ftnlen)7);
#line 291 "csyrfs.f"
    safmin = slamch_("Safe minimum", (ftnlen)12);
#line 292 "csyrfs.f"
    safe1 = nz * safmin;
#line 293 "csyrfs.f"
    safe2 = safe1 / eps;

/*     Do for each right hand side */

#line 297 "csyrfs.f"
    i__1 = *nrhs;
#line 297 "csyrfs.f"
    for (j = 1; j <= i__1; ++j) {

#line 299 "csyrfs.f"
	count = 1;
#line 300 "csyrfs.f"
	lstres = 3.;
#line 301 "csyrfs.f"
L20:

/*        Loop until stopping criterion is satisfied. */

/*        Compute residual R = B - A * X */

#line 307 "csyrfs.f"
	ccopy_(n, &b[j * b_dim1 + 1], &c__1, &work[1], &c__1);
#line 308 "csyrfs.f"
	z__1.r = -1., z__1.i = -0.;
#line 308 "csyrfs.f"
	csymv_(uplo, n, &z__1, &a[a_offset], lda, &x[j * x_dim1 + 1], &c__1, &
		c_b1, &work[1], &c__1, (ftnlen)1);

/*        Compute componentwise relative backward error from formula */

/*        max(i) ( abs(R(i)) / ( abs(A)*abs(X) + abs(B) )(i) ) */

/*        where abs(Z) is the componentwise absolute value of the matrix */
/*        or vector Z.  If the i-th component of the denominator is less */
/*        than SAFE2, then SAFE1 is added to the i-th components of the */
/*        numerator and denominator before dividing. */

#line 319 "csyrfs.f"
	i__2 = *n;
#line 319 "csyrfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 320 "csyrfs.f"
	    i__3 = i__ + j * b_dim1;
#line 320 "csyrfs.f"
	    rwork[i__] = (d__1 = b[i__3].r, abs(d__1)) + (d__2 = d_imag(&b[
		    i__ + j * b_dim1]), abs(d__2));
#line 321 "csyrfs.f"
/* L30: */
#line 321 "csyrfs.f"
	}

/*        Compute abs(A)*abs(X) + abs(B). */

#line 325 "csyrfs.f"
	if (upper) {
#line 326 "csyrfs.f"
	    i__2 = *n;
#line 326 "csyrfs.f"
	    for (k = 1; k <= i__2; ++k) {
#line 327 "csyrfs.f"
		s = 0.;
#line 328 "csyrfs.f"
		i__3 = k + j * x_dim1;
#line 328 "csyrfs.f"
		xk = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = d_imag(&x[k + j *
			 x_dim1]), abs(d__2));
#line 329 "csyrfs.f"
		i__3 = k - 1;
#line 329 "csyrfs.f"
		for (i__ = 1; i__ <= i__3; ++i__) {
#line 330 "csyrfs.f"
		    i__4 = i__ + k * a_dim1;
#line 330 "csyrfs.f"
		    rwork[i__] += ((d__1 = a[i__4].r, abs(d__1)) + (d__2 = 
			    d_imag(&a[i__ + k * a_dim1]), abs(d__2))) * xk;
#line 331 "csyrfs.f"
		    i__4 = i__ + k * a_dim1;
#line 331 "csyrfs.f"
		    i__5 = i__ + j * x_dim1;
#line 331 "csyrfs.f"
		    s += ((d__1 = a[i__4].r, abs(d__1)) + (d__2 = d_imag(&a[
			    i__ + k * a_dim1]), abs(d__2))) * ((d__3 = x[i__5]
			    .r, abs(d__3)) + (d__4 = d_imag(&x[i__ + j * 
			    x_dim1]), abs(d__4)));
#line 332 "csyrfs.f"
/* L40: */
#line 332 "csyrfs.f"
		}
#line 333 "csyrfs.f"
		i__3 = k + k * a_dim1;
#line 333 "csyrfs.f"
		rwork[k] = rwork[k] + ((d__1 = a[i__3].r, abs(d__1)) + (d__2 =
			 d_imag(&a[k + k * a_dim1]), abs(d__2))) * xk + s;
#line 334 "csyrfs.f"
/* L50: */
#line 334 "csyrfs.f"
	    }
#line 335 "csyrfs.f"
	} else {
#line 336 "csyrfs.f"
	    i__2 = *n;
#line 336 "csyrfs.f"
	    for (k = 1; k <= i__2; ++k) {
#line 337 "csyrfs.f"
		s = 0.;
#line 338 "csyrfs.f"
		i__3 = k + j * x_dim1;
#line 338 "csyrfs.f"
		xk = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = d_imag(&x[k + j *
			 x_dim1]), abs(d__2));
#line 339 "csyrfs.f"
		i__3 = k + k * a_dim1;
#line 339 "csyrfs.f"
		rwork[k] += ((d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&
			a[k + k * a_dim1]), abs(d__2))) * xk;
#line 340 "csyrfs.f"
		i__3 = *n;
#line 340 "csyrfs.f"
		for (i__ = k + 1; i__ <= i__3; ++i__) {
#line 341 "csyrfs.f"
		    i__4 = i__ + k * a_dim1;
#line 341 "csyrfs.f"
		    rwork[i__] += ((d__1 = a[i__4].r, abs(d__1)) + (d__2 = 
			    d_imag(&a[i__ + k * a_dim1]), abs(d__2))) * xk;
#line 342 "csyrfs.f"
		    i__4 = i__ + k * a_dim1;
#line 342 "csyrfs.f"
		    i__5 = i__ + j * x_dim1;
#line 342 "csyrfs.f"
		    s += ((d__1 = a[i__4].r, abs(d__1)) + (d__2 = d_imag(&a[
			    i__ + k * a_dim1]), abs(d__2))) * ((d__3 = x[i__5]
			    .r, abs(d__3)) + (d__4 = d_imag(&x[i__ + j * 
			    x_dim1]), abs(d__4)));
#line 343 "csyrfs.f"
/* L60: */
#line 343 "csyrfs.f"
		}
#line 344 "csyrfs.f"
		rwork[k] += s;
#line 345 "csyrfs.f"
/* L70: */
#line 345 "csyrfs.f"
	    }
#line 346 "csyrfs.f"
	}
#line 347 "csyrfs.f"
	s = 0.;
#line 348 "csyrfs.f"
	i__2 = *n;
#line 348 "csyrfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 349 "csyrfs.f"
	    if (rwork[i__] > safe2) {
/* Computing MAX */
#line 350 "csyrfs.f"
		i__3 = i__;
#line 350 "csyrfs.f"
		d__3 = s, d__4 = ((d__1 = work[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&work[i__]), abs(d__2))) / rwork[i__];
#line 350 "csyrfs.f"
		s = max(d__3,d__4);
#line 351 "csyrfs.f"
	    } else {
/* Computing MAX */
#line 352 "csyrfs.f"
		i__3 = i__;
#line 352 "csyrfs.f"
		d__3 = s, d__4 = ((d__1 = work[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&work[i__]), abs(d__2)) + safe1) / (rwork[i__] 
			+ safe1);
#line 352 "csyrfs.f"
		s = max(d__3,d__4);
#line 354 "csyrfs.f"
	    }
#line 355 "csyrfs.f"
/* L80: */
#line 355 "csyrfs.f"
	}
#line 356 "csyrfs.f"
	berr[j] = s;

/*        Test stopping criterion. Continue iterating if */
/*           1) The residual BERR(J) is larger than machine epsilon, and */
/*           2) BERR(J) decreased by at least a factor of 2 during the */
/*              last iteration, and */
/*           3) At most ITMAX iterations tried. */

#line 364 "csyrfs.f"
	if (berr[j] > eps && berr[j] * 2. <= lstres && count <= 5) {

/*           Update solution and try again. */

#line 369 "csyrfs.f"
	    csytrs_(uplo, n, &c__1, &af[af_offset], ldaf, &ipiv[1], &work[1], 
		    n, info, (ftnlen)1);
#line 370 "csyrfs.f"
	    caxpy_(n, &c_b1, &work[1], &c__1, &x[j * x_dim1 + 1], &c__1);
#line 371 "csyrfs.f"
	    lstres = berr[j];
#line 372 "csyrfs.f"
	    ++count;
#line 373 "csyrfs.f"
	    goto L20;
#line 374 "csyrfs.f"
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

/*        Use CLACN2 to estimate the infinity-norm of the matrix */
/*           inv(A) * diag(W), */
/*        where W = abs(R) + NZ*EPS*( abs(A)*abs(X)+abs(B) ))) */

#line 398 "csyrfs.f"
	i__2 = *n;
#line 398 "csyrfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 399 "csyrfs.f"
	    if (rwork[i__] > safe2) {
#line 400 "csyrfs.f"
		i__3 = i__;
#line 400 "csyrfs.f"
		rwork[i__] = (d__1 = work[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&work[i__]), abs(d__2)) + nz * eps * rwork[i__]
			;
#line 401 "csyrfs.f"
	    } else {
#line 402 "csyrfs.f"
		i__3 = i__;
#line 402 "csyrfs.f"
		rwork[i__] = (d__1 = work[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&work[i__]), abs(d__2)) + nz * eps * rwork[i__]
			 + safe1;
#line 404 "csyrfs.f"
	    }
#line 405 "csyrfs.f"
/* L90: */
#line 405 "csyrfs.f"
	}

#line 407 "csyrfs.f"
	kase = 0;
#line 408 "csyrfs.f"
L100:
#line 409 "csyrfs.f"
	clacn2_(n, &work[*n + 1], &work[1], &ferr[j], &kase, isave);
#line 410 "csyrfs.f"
	if (kase != 0) {
#line 411 "csyrfs.f"
	    if (kase == 1) {

/*              Multiply by diag(W)*inv(A**T). */

#line 415 "csyrfs.f"
		csytrs_(uplo, n, &c__1, &af[af_offset], ldaf, &ipiv[1], &work[
			1], n, info, (ftnlen)1);
#line 416 "csyrfs.f"
		i__2 = *n;
#line 416 "csyrfs.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 417 "csyrfs.f"
		    i__3 = i__;
#line 417 "csyrfs.f"
		    i__4 = i__;
#line 417 "csyrfs.f"
		    i__5 = i__;
#line 417 "csyrfs.f"
		    z__1.r = rwork[i__4] * work[i__5].r, z__1.i = rwork[i__4] 
			    * work[i__5].i;
#line 417 "csyrfs.f"
		    work[i__3].r = z__1.r, work[i__3].i = z__1.i;
#line 418 "csyrfs.f"
/* L110: */
#line 418 "csyrfs.f"
		}
#line 419 "csyrfs.f"
	    } else if (kase == 2) {

/*              Multiply by inv(A)*diag(W). */

#line 423 "csyrfs.f"
		i__2 = *n;
#line 423 "csyrfs.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 424 "csyrfs.f"
		    i__3 = i__;
#line 424 "csyrfs.f"
		    i__4 = i__;
#line 424 "csyrfs.f"
		    i__5 = i__;
#line 424 "csyrfs.f"
		    z__1.r = rwork[i__4] * work[i__5].r, z__1.i = rwork[i__4] 
			    * work[i__5].i;
#line 424 "csyrfs.f"
		    work[i__3].r = z__1.r, work[i__3].i = z__1.i;
#line 425 "csyrfs.f"
/* L120: */
#line 425 "csyrfs.f"
		}
#line 426 "csyrfs.f"
		csytrs_(uplo, n, &c__1, &af[af_offset], ldaf, &ipiv[1], &work[
			1], n, info, (ftnlen)1);
#line 427 "csyrfs.f"
	    }
#line 428 "csyrfs.f"
	    goto L100;
#line 429 "csyrfs.f"
	}

/*        Normalize error. */

#line 433 "csyrfs.f"
	lstres = 0.;
#line 434 "csyrfs.f"
	i__2 = *n;
#line 434 "csyrfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 435 "csyrfs.f"
	    i__3 = i__ + j * x_dim1;
#line 435 "csyrfs.f"
	    d__3 = lstres, d__4 = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = 
		    d_imag(&x[i__ + j * x_dim1]), abs(d__2));
#line 435 "csyrfs.f"
	    lstres = max(d__3,d__4);
#line 436 "csyrfs.f"
/* L130: */
#line 436 "csyrfs.f"
	}
#line 437 "csyrfs.f"
	if (lstres != 0.) {
#line 437 "csyrfs.f"
	    ferr[j] /= lstres;
#line 437 "csyrfs.f"
	}

#line 440 "csyrfs.f"
/* L140: */
#line 440 "csyrfs.f"
    }

#line 442 "csyrfs.f"
    return 0;

/*     End of CSYRFS */

} /* csyrfs_ */

