#line 1 "cporfs.f"
/* cporfs.f -- translated by f2c (version 20100827).
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

#line 1 "cporfs.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static integer c__1 = 1;

/* > \brief \b CPORFS */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CPORFS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cporfs.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cporfs.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cporfs.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CPORFS( UPLO, N, NRHS, A, LDA, AF, LDAF, B, LDB, X, */
/*                          LDX, FERR, BERR, WORK, RWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDA, LDAF, LDB, LDX, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               BERR( * ), FERR( * ), RWORK( * ) */
/*       COMPLEX            A( LDA, * ), AF( LDAF, * ), B( LDB, * ), */
/*      $                   WORK( * ), X( LDX, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CPORFS improves the computed solution to a system of linear */
/* > equations when the coefficient matrix is Hermitian positive definite, */
/* > and provides error bounds and backward error estimates for the */
/* > solution. */
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
/* >          AF is COMPLEX array, dimension (LDAF,N) */
/* >          The triangular factor U or L from the Cholesky factorization */
/* >          A = U**H*U or A = L*L**H, as computed by CPOTRF. */
/* > \endverbatim */
/* > */
/* > \param[in] LDAF */
/* > \verbatim */
/* >          LDAF is INTEGER */
/* >          The leading dimension of the array AF.  LDAF >= max(1,N). */
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
/* >          On entry, the solution matrix X, as computed by CPOTRS. */
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

/* > \date November 2011 */

/* > \ingroup complexPOcomputational */

/*  ===================================================================== */
/* Subroutine */ int cporfs_(char *uplo, integer *n, integer *nrhs, 
	doublecomplex *a, integer *lda, doublecomplex *af, integer *ldaf, 
	doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx, 
	doublereal *ferr, doublereal *berr, doublecomplex *work, doublereal *
	rwork, integer *info, ftnlen uplo_len)
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
    extern /* Subroutine */ int chemv_(char *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, ftnlen);
    static integer isave[3];
    extern /* Subroutine */ int ccopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), caxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    static integer count;
    static logical upper;
    extern /* Subroutine */ int clacn2_(integer *, doublecomplex *, 
	    doublecomplex *, doublereal *, integer *, integer *);
    extern doublereal slamch_(char *, ftnlen);
    static doublereal safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), cpotrs_(
	    char *, integer *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, integer *, ftnlen);
    static doublereal lstres;


/*  -- LAPACK computational routine (version 3.4.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ==================================================================== */

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

#line 245 "cporfs.f"
    /* Parameter adjustments */
#line 245 "cporfs.f"
    a_dim1 = *lda;
#line 245 "cporfs.f"
    a_offset = 1 + a_dim1;
#line 245 "cporfs.f"
    a -= a_offset;
#line 245 "cporfs.f"
    af_dim1 = *ldaf;
#line 245 "cporfs.f"
    af_offset = 1 + af_dim1;
#line 245 "cporfs.f"
    af -= af_offset;
#line 245 "cporfs.f"
    b_dim1 = *ldb;
#line 245 "cporfs.f"
    b_offset = 1 + b_dim1;
#line 245 "cporfs.f"
    b -= b_offset;
#line 245 "cporfs.f"
    x_dim1 = *ldx;
#line 245 "cporfs.f"
    x_offset = 1 + x_dim1;
#line 245 "cporfs.f"
    x -= x_offset;
#line 245 "cporfs.f"
    --ferr;
#line 245 "cporfs.f"
    --berr;
#line 245 "cporfs.f"
    --work;
#line 245 "cporfs.f"
    --rwork;
#line 245 "cporfs.f"

#line 245 "cporfs.f"
    /* Function Body */
#line 245 "cporfs.f"
    *info = 0;
#line 246 "cporfs.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 247 "cporfs.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 248 "cporfs.f"
	*info = -1;
#line 249 "cporfs.f"
    } else if (*n < 0) {
#line 250 "cporfs.f"
	*info = -2;
#line 251 "cporfs.f"
    } else if (*nrhs < 0) {
#line 252 "cporfs.f"
	*info = -3;
#line 253 "cporfs.f"
    } else if (*lda < max(1,*n)) {
#line 254 "cporfs.f"
	*info = -5;
#line 255 "cporfs.f"
    } else if (*ldaf < max(1,*n)) {
#line 256 "cporfs.f"
	*info = -7;
#line 257 "cporfs.f"
    } else if (*ldb < max(1,*n)) {
#line 258 "cporfs.f"
	*info = -9;
#line 259 "cporfs.f"
    } else if (*ldx < max(1,*n)) {
#line 260 "cporfs.f"
	*info = -11;
#line 261 "cporfs.f"
    }
#line 262 "cporfs.f"
    if (*info != 0) {
#line 263 "cporfs.f"
	i__1 = -(*info);
#line 263 "cporfs.f"
	xerbla_("CPORFS", &i__1, (ftnlen)6);
#line 264 "cporfs.f"
	return 0;
#line 265 "cporfs.f"
    }

/*     Quick return if possible */

#line 269 "cporfs.f"
    if (*n == 0 || *nrhs == 0) {
#line 270 "cporfs.f"
	i__1 = *nrhs;
#line 270 "cporfs.f"
	for (j = 1; j <= i__1; ++j) {
#line 271 "cporfs.f"
	    ferr[j] = 0.;
#line 272 "cporfs.f"
	    berr[j] = 0.;
#line 273 "cporfs.f"
/* L10: */
#line 273 "cporfs.f"
	}
#line 274 "cporfs.f"
	return 0;
#line 275 "cporfs.f"
    }

/*     NZ = maximum number of nonzero elements in each row of A, plus 1 */

#line 279 "cporfs.f"
    nz = *n + 1;
#line 280 "cporfs.f"
    eps = slamch_("Epsilon", (ftnlen)7);
#line 281 "cporfs.f"
    safmin = slamch_("Safe minimum", (ftnlen)12);
#line 282 "cporfs.f"
    safe1 = nz * safmin;
#line 283 "cporfs.f"
    safe2 = safe1 / eps;

/*     Do for each right hand side */

#line 287 "cporfs.f"
    i__1 = *nrhs;
#line 287 "cporfs.f"
    for (j = 1; j <= i__1; ++j) {

#line 289 "cporfs.f"
	count = 1;
#line 290 "cporfs.f"
	lstres = 3.;
#line 291 "cporfs.f"
L20:

/*        Loop until stopping criterion is satisfied. */

/*        Compute residual R = B - A * X */

#line 297 "cporfs.f"
	ccopy_(n, &b[j * b_dim1 + 1], &c__1, &work[1], &c__1);
#line 298 "cporfs.f"
	z__1.r = -1., z__1.i = -0.;
#line 298 "cporfs.f"
	chemv_(uplo, n, &z__1, &a[a_offset], lda, &x[j * x_dim1 + 1], &c__1, &
		c_b1, &work[1], &c__1, (ftnlen)1);

/*        Compute componentwise relative backward error from formula */

/*        max(i) ( abs(R(i)) / ( abs(A)*abs(X) + abs(B) )(i) ) */

/*        where abs(Z) is the componentwise absolute value of the matrix */
/*        or vector Z.  If the i-th component of the denominator is less */
/*        than SAFE2, then SAFE1 is added to the i-th components of the */
/*        numerator and denominator before dividing. */

#line 309 "cporfs.f"
	i__2 = *n;
#line 309 "cporfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 310 "cporfs.f"
	    i__3 = i__ + j * b_dim1;
#line 310 "cporfs.f"
	    rwork[i__] = (d__1 = b[i__3].r, abs(d__1)) + (d__2 = d_imag(&b[
		    i__ + j * b_dim1]), abs(d__2));
#line 311 "cporfs.f"
/* L30: */
#line 311 "cporfs.f"
	}

/*        Compute abs(A)*abs(X) + abs(B). */

#line 315 "cporfs.f"
	if (upper) {
#line 316 "cporfs.f"
	    i__2 = *n;
#line 316 "cporfs.f"
	    for (k = 1; k <= i__2; ++k) {
#line 317 "cporfs.f"
		s = 0.;
#line 318 "cporfs.f"
		i__3 = k + j * x_dim1;
#line 318 "cporfs.f"
		xk = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = d_imag(&x[k + j *
			 x_dim1]), abs(d__2));
#line 319 "cporfs.f"
		i__3 = k - 1;
#line 319 "cporfs.f"
		for (i__ = 1; i__ <= i__3; ++i__) {
#line 320 "cporfs.f"
		    i__4 = i__ + k * a_dim1;
#line 320 "cporfs.f"
		    rwork[i__] += ((d__1 = a[i__4].r, abs(d__1)) + (d__2 = 
			    d_imag(&a[i__ + k * a_dim1]), abs(d__2))) * xk;
#line 321 "cporfs.f"
		    i__4 = i__ + k * a_dim1;
#line 321 "cporfs.f"
		    i__5 = i__ + j * x_dim1;
#line 321 "cporfs.f"
		    s += ((d__1 = a[i__4].r, abs(d__1)) + (d__2 = d_imag(&a[
			    i__ + k * a_dim1]), abs(d__2))) * ((d__3 = x[i__5]
			    .r, abs(d__3)) + (d__4 = d_imag(&x[i__ + j * 
			    x_dim1]), abs(d__4)));
#line 322 "cporfs.f"
/* L40: */
#line 322 "cporfs.f"
		}
#line 323 "cporfs.f"
		i__3 = k + k * a_dim1;
#line 323 "cporfs.f"
		rwork[k] = rwork[k] + (d__1 = a[i__3].r, abs(d__1)) * xk + s;
#line 324 "cporfs.f"
/* L50: */
#line 324 "cporfs.f"
	    }
#line 325 "cporfs.f"
	} else {
#line 326 "cporfs.f"
	    i__2 = *n;
#line 326 "cporfs.f"
	    for (k = 1; k <= i__2; ++k) {
#line 327 "cporfs.f"
		s = 0.;
#line 328 "cporfs.f"
		i__3 = k + j * x_dim1;
#line 328 "cporfs.f"
		xk = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = d_imag(&x[k + j *
			 x_dim1]), abs(d__2));
#line 329 "cporfs.f"
		i__3 = k + k * a_dim1;
#line 329 "cporfs.f"
		rwork[k] += (d__1 = a[i__3].r, abs(d__1)) * xk;
#line 330 "cporfs.f"
		i__3 = *n;
#line 330 "cporfs.f"
		for (i__ = k + 1; i__ <= i__3; ++i__) {
#line 331 "cporfs.f"
		    i__4 = i__ + k * a_dim1;
#line 331 "cporfs.f"
		    rwork[i__] += ((d__1 = a[i__4].r, abs(d__1)) + (d__2 = 
			    d_imag(&a[i__ + k * a_dim1]), abs(d__2))) * xk;
#line 332 "cporfs.f"
		    i__4 = i__ + k * a_dim1;
#line 332 "cporfs.f"
		    i__5 = i__ + j * x_dim1;
#line 332 "cporfs.f"
		    s += ((d__1 = a[i__4].r, abs(d__1)) + (d__2 = d_imag(&a[
			    i__ + k * a_dim1]), abs(d__2))) * ((d__3 = x[i__5]
			    .r, abs(d__3)) + (d__4 = d_imag(&x[i__ + j * 
			    x_dim1]), abs(d__4)));
#line 333 "cporfs.f"
/* L60: */
#line 333 "cporfs.f"
		}
#line 334 "cporfs.f"
		rwork[k] += s;
#line 335 "cporfs.f"
/* L70: */
#line 335 "cporfs.f"
	    }
#line 336 "cporfs.f"
	}
#line 337 "cporfs.f"
	s = 0.;
#line 338 "cporfs.f"
	i__2 = *n;
#line 338 "cporfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 339 "cporfs.f"
	    if (rwork[i__] > safe2) {
/* Computing MAX */
#line 340 "cporfs.f"
		i__3 = i__;
#line 340 "cporfs.f"
		d__3 = s, d__4 = ((d__1 = work[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&work[i__]), abs(d__2))) / rwork[i__];
#line 340 "cporfs.f"
		s = max(d__3,d__4);
#line 341 "cporfs.f"
	    } else {
/* Computing MAX */
#line 342 "cporfs.f"
		i__3 = i__;
#line 342 "cporfs.f"
		d__3 = s, d__4 = ((d__1 = work[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&work[i__]), abs(d__2)) + safe1) / (rwork[i__] 
			+ safe1);
#line 342 "cporfs.f"
		s = max(d__3,d__4);
#line 344 "cporfs.f"
	    }
#line 345 "cporfs.f"
/* L80: */
#line 345 "cporfs.f"
	}
#line 346 "cporfs.f"
	berr[j] = s;

/*        Test stopping criterion. Continue iterating if */
/*           1) The residual BERR(J) is larger than machine epsilon, and */
/*           2) BERR(J) decreased by at least a factor of 2 during the */
/*              last iteration, and */
/*           3) At most ITMAX iterations tried. */

#line 354 "cporfs.f"
	if (berr[j] > eps && berr[j] * 2. <= lstres && count <= 5) {

/*           Update solution and try again. */

#line 359 "cporfs.f"
	    cpotrs_(uplo, n, &c__1, &af[af_offset], ldaf, &work[1], n, info, (
		    ftnlen)1);
#line 360 "cporfs.f"
	    caxpy_(n, &c_b1, &work[1], &c__1, &x[j * x_dim1 + 1], &c__1);
#line 361 "cporfs.f"
	    lstres = berr[j];
#line 362 "cporfs.f"
	    ++count;
#line 363 "cporfs.f"
	    goto L20;
#line 364 "cporfs.f"
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

#line 388 "cporfs.f"
	i__2 = *n;
#line 388 "cporfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 389 "cporfs.f"
	    if (rwork[i__] > safe2) {
#line 390 "cporfs.f"
		i__3 = i__;
#line 390 "cporfs.f"
		rwork[i__] = (d__1 = work[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&work[i__]), abs(d__2)) + nz * eps * rwork[i__]
			;
#line 391 "cporfs.f"
	    } else {
#line 392 "cporfs.f"
		i__3 = i__;
#line 392 "cporfs.f"
		rwork[i__] = (d__1 = work[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&work[i__]), abs(d__2)) + nz * eps * rwork[i__]
			 + safe1;
#line 394 "cporfs.f"
	    }
#line 395 "cporfs.f"
/* L90: */
#line 395 "cporfs.f"
	}

#line 397 "cporfs.f"
	kase = 0;
#line 398 "cporfs.f"
L100:
#line 399 "cporfs.f"
	clacn2_(n, &work[*n + 1], &work[1], &ferr[j], &kase, isave);
#line 400 "cporfs.f"
	if (kase != 0) {
#line 401 "cporfs.f"
	    if (kase == 1) {

/*              Multiply by diag(W)*inv(A**H). */

#line 405 "cporfs.f"
		cpotrs_(uplo, n, &c__1, &af[af_offset], ldaf, &work[1], n, 
			info, (ftnlen)1);
#line 406 "cporfs.f"
		i__2 = *n;
#line 406 "cporfs.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 407 "cporfs.f"
		    i__3 = i__;
#line 407 "cporfs.f"
		    i__4 = i__;
#line 407 "cporfs.f"
		    i__5 = i__;
#line 407 "cporfs.f"
		    z__1.r = rwork[i__4] * work[i__5].r, z__1.i = rwork[i__4] 
			    * work[i__5].i;
#line 407 "cporfs.f"
		    work[i__3].r = z__1.r, work[i__3].i = z__1.i;
#line 408 "cporfs.f"
/* L110: */
#line 408 "cporfs.f"
		}
#line 409 "cporfs.f"
	    } else if (kase == 2) {

/*              Multiply by inv(A)*diag(W). */

#line 413 "cporfs.f"
		i__2 = *n;
#line 413 "cporfs.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 414 "cporfs.f"
		    i__3 = i__;
#line 414 "cporfs.f"
		    i__4 = i__;
#line 414 "cporfs.f"
		    i__5 = i__;
#line 414 "cporfs.f"
		    z__1.r = rwork[i__4] * work[i__5].r, z__1.i = rwork[i__4] 
			    * work[i__5].i;
#line 414 "cporfs.f"
		    work[i__3].r = z__1.r, work[i__3].i = z__1.i;
#line 415 "cporfs.f"
/* L120: */
#line 415 "cporfs.f"
		}
#line 416 "cporfs.f"
		cpotrs_(uplo, n, &c__1, &af[af_offset], ldaf, &work[1], n, 
			info, (ftnlen)1);
#line 417 "cporfs.f"
	    }
#line 418 "cporfs.f"
	    goto L100;
#line 419 "cporfs.f"
	}

/*        Normalize error. */

#line 423 "cporfs.f"
	lstres = 0.;
#line 424 "cporfs.f"
	i__2 = *n;
#line 424 "cporfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 425 "cporfs.f"
	    i__3 = i__ + j * x_dim1;
#line 425 "cporfs.f"
	    d__3 = lstres, d__4 = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = 
		    d_imag(&x[i__ + j * x_dim1]), abs(d__2));
#line 425 "cporfs.f"
	    lstres = max(d__3,d__4);
#line 426 "cporfs.f"
/* L130: */
#line 426 "cporfs.f"
	}
#line 427 "cporfs.f"
	if (lstres != 0.) {
#line 427 "cporfs.f"
	    ferr[j] /= lstres;
#line 427 "cporfs.f"
	}

#line 430 "cporfs.f"
/* L140: */
#line 430 "cporfs.f"
    }

#line 432 "cporfs.f"
    return 0;

/*     End of CPORFS */

} /* cporfs_ */

