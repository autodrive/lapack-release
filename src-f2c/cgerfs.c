#line 1 "cgerfs.f"
/* cgerfs.f -- translated by f2c (version 20100827).
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

#line 1 "cgerfs.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static integer c__1 = 1;

/* > \brief \b CGERFS */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CGERFS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgerfs.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgerfs.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgerfs.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CGERFS( TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV, B, LDB, */
/*                          X, LDX, FERR, BERR, WORK, RWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          TRANS */
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
/* > CGERFS improves the computed solution to a system of linear */
/* > equations and provides error bounds and backward error estimates for */
/* > the solution. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] TRANS */
/* > \verbatim */
/* >          TRANS is CHARACTER*1 */
/* >          Specifies the form of the system of equations: */
/* >          = 'N':  A * X = B     (No transpose) */
/* >          = 'T':  A**T * X = B  (Transpose) */
/* >          = 'C':  A**H * X = B  (Conjugate transpose) */
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
/* >          The original N-by-N matrix A. */
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
/* >          The factors L and U from the factorization A = P*L*U */
/* >          as computed by CGETRF. */
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
/* >          The pivot indices from CGETRF; for 1<=i<=N, row i of the */
/* >          matrix was interchanged with row IPIV(i). */
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
/* >          On entry, the solution matrix X, as computed by CGETRS. */
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

/* > \ingroup complexGEcomputational */

/*  ===================================================================== */
/* Subroutine */ int cgerfs_(char *trans, integer *n, integer *nrhs, 
	doublecomplex *a, integer *lda, doublecomplex *af, integer *ldaf, 
	integer *ipiv, doublecomplex *b, integer *ldb, doublecomplex *x, 
	integer *ldx, doublereal *ferr, doublereal *berr, doublecomplex *work,
	 doublereal *rwork, integer *info, ftnlen trans_len)
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
    extern /* Subroutine */ int cgemv_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, ftnlen);
    static integer isave[3];
    extern /* Subroutine */ int ccopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), caxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    static integer count;
    extern /* Subroutine */ int clacn2_(integer *, doublecomplex *, 
	    doublecomplex *, doublereal *, integer *, integer *);
    extern doublereal slamch_(char *, ftnlen);
    static doublereal safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), cgetrs_(
	    char *, integer *, integer *, doublecomplex *, integer *, integer 
	    *, doublecomplex *, integer *, integer *, ftnlen);
    static logical notran;
    static char transn[1], transt[1];
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
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Statement Functions .. */
/*     .. */
/*     .. Statement Function definitions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 250 "cgerfs.f"
    /* Parameter adjustments */
#line 250 "cgerfs.f"
    a_dim1 = *lda;
#line 250 "cgerfs.f"
    a_offset = 1 + a_dim1;
#line 250 "cgerfs.f"
    a -= a_offset;
#line 250 "cgerfs.f"
    af_dim1 = *ldaf;
#line 250 "cgerfs.f"
    af_offset = 1 + af_dim1;
#line 250 "cgerfs.f"
    af -= af_offset;
#line 250 "cgerfs.f"
    --ipiv;
#line 250 "cgerfs.f"
    b_dim1 = *ldb;
#line 250 "cgerfs.f"
    b_offset = 1 + b_dim1;
#line 250 "cgerfs.f"
    b -= b_offset;
#line 250 "cgerfs.f"
    x_dim1 = *ldx;
#line 250 "cgerfs.f"
    x_offset = 1 + x_dim1;
#line 250 "cgerfs.f"
    x -= x_offset;
#line 250 "cgerfs.f"
    --ferr;
#line 250 "cgerfs.f"
    --berr;
#line 250 "cgerfs.f"
    --work;
#line 250 "cgerfs.f"
    --rwork;
#line 250 "cgerfs.f"

#line 250 "cgerfs.f"
    /* Function Body */
#line 250 "cgerfs.f"
    *info = 0;
#line 251 "cgerfs.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 252 "cgerfs.f"
    if (! notran && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1) && ! lsame_(
	    trans, "C", (ftnlen)1, (ftnlen)1)) {
#line 254 "cgerfs.f"
	*info = -1;
#line 255 "cgerfs.f"
    } else if (*n < 0) {
#line 256 "cgerfs.f"
	*info = -2;
#line 257 "cgerfs.f"
    } else if (*nrhs < 0) {
#line 258 "cgerfs.f"
	*info = -3;
#line 259 "cgerfs.f"
    } else if (*lda < max(1,*n)) {
#line 260 "cgerfs.f"
	*info = -5;
#line 261 "cgerfs.f"
    } else if (*ldaf < max(1,*n)) {
#line 262 "cgerfs.f"
	*info = -7;
#line 263 "cgerfs.f"
    } else if (*ldb < max(1,*n)) {
#line 264 "cgerfs.f"
	*info = -10;
#line 265 "cgerfs.f"
    } else if (*ldx < max(1,*n)) {
#line 266 "cgerfs.f"
	*info = -12;
#line 267 "cgerfs.f"
    }
#line 268 "cgerfs.f"
    if (*info != 0) {
#line 269 "cgerfs.f"
	i__1 = -(*info);
#line 269 "cgerfs.f"
	xerbla_("CGERFS", &i__1, (ftnlen)6);
#line 270 "cgerfs.f"
	return 0;
#line 271 "cgerfs.f"
    }

/*     Quick return if possible */

#line 275 "cgerfs.f"
    if (*n == 0 || *nrhs == 0) {
#line 276 "cgerfs.f"
	i__1 = *nrhs;
#line 276 "cgerfs.f"
	for (j = 1; j <= i__1; ++j) {
#line 277 "cgerfs.f"
	    ferr[j] = 0.;
#line 278 "cgerfs.f"
	    berr[j] = 0.;
#line 279 "cgerfs.f"
/* L10: */
#line 279 "cgerfs.f"
	}
#line 280 "cgerfs.f"
	return 0;
#line 281 "cgerfs.f"
    }

#line 283 "cgerfs.f"
    if (notran) {
#line 284 "cgerfs.f"
	*(unsigned char *)transn = 'N';
#line 285 "cgerfs.f"
	*(unsigned char *)transt = 'C';
#line 286 "cgerfs.f"
    } else {
#line 287 "cgerfs.f"
	*(unsigned char *)transn = 'C';
#line 288 "cgerfs.f"
	*(unsigned char *)transt = 'N';
#line 289 "cgerfs.f"
    }

/*     NZ = maximum number of nonzero elements in each row of A, plus 1 */

#line 293 "cgerfs.f"
    nz = *n + 1;
#line 294 "cgerfs.f"
    eps = slamch_("Epsilon", (ftnlen)7);
#line 295 "cgerfs.f"
    safmin = slamch_("Safe minimum", (ftnlen)12);
#line 296 "cgerfs.f"
    safe1 = nz * safmin;
#line 297 "cgerfs.f"
    safe2 = safe1 / eps;

/*     Do for each right hand side */

#line 301 "cgerfs.f"
    i__1 = *nrhs;
#line 301 "cgerfs.f"
    for (j = 1; j <= i__1; ++j) {

#line 303 "cgerfs.f"
	count = 1;
#line 304 "cgerfs.f"
	lstres = 3.;
#line 305 "cgerfs.f"
L20:

/*        Loop until stopping criterion is satisfied. */

/*        Compute residual R = B - op(A) * X, */
/*        where op(A) = A, A**T, or A**H, depending on TRANS. */

#line 312 "cgerfs.f"
	ccopy_(n, &b[j * b_dim1 + 1], &c__1, &work[1], &c__1);
#line 313 "cgerfs.f"
	z__1.r = -1., z__1.i = -0.;
#line 313 "cgerfs.f"
	cgemv_(trans, n, n, &z__1, &a[a_offset], lda, &x[j * x_dim1 + 1], &
		c__1, &c_b1, &work[1], &c__1, (ftnlen)1);

/*        Compute componentwise relative backward error from formula */

/*        max(i) ( abs(R(i)) / ( abs(op(A))*abs(X) + abs(B) )(i) ) */

/*        where abs(Z) is the componentwise absolute value of the matrix */
/*        or vector Z.  If the i-th component of the denominator is less */
/*        than SAFE2, then SAFE1 is added to the i-th components of the */
/*        numerator and denominator before dividing. */

#line 325 "cgerfs.f"
	i__2 = *n;
#line 325 "cgerfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 326 "cgerfs.f"
	    i__3 = i__ + j * b_dim1;
#line 326 "cgerfs.f"
	    rwork[i__] = (d__1 = b[i__3].r, abs(d__1)) + (d__2 = d_imag(&b[
		    i__ + j * b_dim1]), abs(d__2));
#line 327 "cgerfs.f"
/* L30: */
#line 327 "cgerfs.f"
	}

/*        Compute abs(op(A))*abs(X) + abs(B). */

#line 331 "cgerfs.f"
	if (notran) {
#line 332 "cgerfs.f"
	    i__2 = *n;
#line 332 "cgerfs.f"
	    for (k = 1; k <= i__2; ++k) {
#line 333 "cgerfs.f"
		i__3 = k + j * x_dim1;
#line 333 "cgerfs.f"
		xk = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = d_imag(&x[k + j *
			 x_dim1]), abs(d__2));
#line 334 "cgerfs.f"
		i__3 = *n;
#line 334 "cgerfs.f"
		for (i__ = 1; i__ <= i__3; ++i__) {
#line 335 "cgerfs.f"
		    i__4 = i__ + k * a_dim1;
#line 335 "cgerfs.f"
		    rwork[i__] += ((d__1 = a[i__4].r, abs(d__1)) + (d__2 = 
			    d_imag(&a[i__ + k * a_dim1]), abs(d__2))) * xk;
#line 336 "cgerfs.f"
/* L40: */
#line 336 "cgerfs.f"
		}
#line 337 "cgerfs.f"
/* L50: */
#line 337 "cgerfs.f"
	    }
#line 338 "cgerfs.f"
	} else {
#line 339 "cgerfs.f"
	    i__2 = *n;
#line 339 "cgerfs.f"
	    for (k = 1; k <= i__2; ++k) {
#line 340 "cgerfs.f"
		s = 0.;
#line 341 "cgerfs.f"
		i__3 = *n;
#line 341 "cgerfs.f"
		for (i__ = 1; i__ <= i__3; ++i__) {
#line 342 "cgerfs.f"
		    i__4 = i__ + k * a_dim1;
#line 342 "cgerfs.f"
		    i__5 = i__ + j * x_dim1;
#line 342 "cgerfs.f"
		    s += ((d__1 = a[i__4].r, abs(d__1)) + (d__2 = d_imag(&a[
			    i__ + k * a_dim1]), abs(d__2))) * ((d__3 = x[i__5]
			    .r, abs(d__3)) + (d__4 = d_imag(&x[i__ + j * 
			    x_dim1]), abs(d__4)));
#line 343 "cgerfs.f"
/* L60: */
#line 343 "cgerfs.f"
		}
#line 344 "cgerfs.f"
		rwork[k] += s;
#line 345 "cgerfs.f"
/* L70: */
#line 345 "cgerfs.f"
	    }
#line 346 "cgerfs.f"
	}
#line 347 "cgerfs.f"
	s = 0.;
#line 348 "cgerfs.f"
	i__2 = *n;
#line 348 "cgerfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 349 "cgerfs.f"
	    if (rwork[i__] > safe2) {
/* Computing MAX */
#line 350 "cgerfs.f"
		i__3 = i__;
#line 350 "cgerfs.f"
		d__3 = s, d__4 = ((d__1 = work[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&work[i__]), abs(d__2))) / rwork[i__];
#line 350 "cgerfs.f"
		s = max(d__3,d__4);
#line 351 "cgerfs.f"
	    } else {
/* Computing MAX */
#line 352 "cgerfs.f"
		i__3 = i__;
#line 352 "cgerfs.f"
		d__3 = s, d__4 = ((d__1 = work[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&work[i__]), abs(d__2)) + safe1) / (rwork[i__] 
			+ safe1);
#line 352 "cgerfs.f"
		s = max(d__3,d__4);
#line 354 "cgerfs.f"
	    }
#line 355 "cgerfs.f"
/* L80: */
#line 355 "cgerfs.f"
	}
#line 356 "cgerfs.f"
	berr[j] = s;

/*        Test stopping criterion. Continue iterating if */
/*           1) The residual BERR(J) is larger than machine epsilon, and */
/*           2) BERR(J) decreased by at least a factor of 2 during the */
/*              last iteration, and */
/*           3) At most ITMAX iterations tried. */

#line 364 "cgerfs.f"
	if (berr[j] > eps && berr[j] * 2. <= lstres && count <= 5) {

/*           Update solution and try again. */

#line 369 "cgerfs.f"
	    cgetrs_(trans, n, &c__1, &af[af_offset], ldaf, &ipiv[1], &work[1],
		     n, info, (ftnlen)1);
#line 370 "cgerfs.f"
	    caxpy_(n, &c_b1, &work[1], &c__1, &x[j * x_dim1 + 1], &c__1);
#line 371 "cgerfs.f"
	    lstres = berr[j];
#line 372 "cgerfs.f"
	    ++count;
#line 373 "cgerfs.f"
	    goto L20;
#line 374 "cgerfs.f"
	}

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

#line 398 "cgerfs.f"
	i__2 = *n;
#line 398 "cgerfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 399 "cgerfs.f"
	    if (rwork[i__] > safe2) {
#line 400 "cgerfs.f"
		i__3 = i__;
#line 400 "cgerfs.f"
		rwork[i__] = (d__1 = work[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&work[i__]), abs(d__2)) + nz * eps * rwork[i__]
			;
#line 401 "cgerfs.f"
	    } else {
#line 402 "cgerfs.f"
		i__3 = i__;
#line 402 "cgerfs.f"
		rwork[i__] = (d__1 = work[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&work[i__]), abs(d__2)) + nz * eps * rwork[i__]
			 + safe1;
#line 404 "cgerfs.f"
	    }
#line 405 "cgerfs.f"
/* L90: */
#line 405 "cgerfs.f"
	}

#line 407 "cgerfs.f"
	kase = 0;
#line 408 "cgerfs.f"
L100:
#line 409 "cgerfs.f"
	clacn2_(n, &work[*n + 1], &work[1], &ferr[j], &kase, isave);
#line 410 "cgerfs.f"
	if (kase != 0) {
#line 411 "cgerfs.f"
	    if (kase == 1) {

/*              Multiply by diag(W)*inv(op(A)**H). */

#line 415 "cgerfs.f"
		cgetrs_(transt, n, &c__1, &af[af_offset], ldaf, &ipiv[1], &
			work[1], n, info, (ftnlen)1);
#line 417 "cgerfs.f"
		i__2 = *n;
#line 417 "cgerfs.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 418 "cgerfs.f"
		    i__3 = i__;
#line 418 "cgerfs.f"
		    i__4 = i__;
#line 418 "cgerfs.f"
		    i__5 = i__;
#line 418 "cgerfs.f"
		    z__1.r = rwork[i__4] * work[i__5].r, z__1.i = rwork[i__4] 
			    * work[i__5].i;
#line 418 "cgerfs.f"
		    work[i__3].r = z__1.r, work[i__3].i = z__1.i;
#line 419 "cgerfs.f"
/* L110: */
#line 419 "cgerfs.f"
		}
#line 420 "cgerfs.f"
	    } else {

/*              Multiply by inv(op(A))*diag(W). */

#line 424 "cgerfs.f"
		i__2 = *n;
#line 424 "cgerfs.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 425 "cgerfs.f"
		    i__3 = i__;
#line 425 "cgerfs.f"
		    i__4 = i__;
#line 425 "cgerfs.f"
		    i__5 = i__;
#line 425 "cgerfs.f"
		    z__1.r = rwork[i__4] * work[i__5].r, z__1.i = rwork[i__4] 
			    * work[i__5].i;
#line 425 "cgerfs.f"
		    work[i__3].r = z__1.r, work[i__3].i = z__1.i;
#line 426 "cgerfs.f"
/* L120: */
#line 426 "cgerfs.f"
		}
#line 427 "cgerfs.f"
		cgetrs_(transn, n, &c__1, &af[af_offset], ldaf, &ipiv[1], &
			work[1], n, info, (ftnlen)1);
#line 429 "cgerfs.f"
	    }
#line 430 "cgerfs.f"
	    goto L100;
#line 431 "cgerfs.f"
	}

/*        Normalize error. */

#line 435 "cgerfs.f"
	lstres = 0.;
#line 436 "cgerfs.f"
	i__2 = *n;
#line 436 "cgerfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 437 "cgerfs.f"
	    i__3 = i__ + j * x_dim1;
#line 437 "cgerfs.f"
	    d__3 = lstres, d__4 = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = 
		    d_imag(&x[i__ + j * x_dim1]), abs(d__2));
#line 437 "cgerfs.f"
	    lstres = max(d__3,d__4);
#line 438 "cgerfs.f"
/* L130: */
#line 438 "cgerfs.f"
	}
#line 439 "cgerfs.f"
	if (lstres != 0.) {
#line 439 "cgerfs.f"
	    ferr[j] /= lstres;
#line 439 "cgerfs.f"
	}

#line 442 "cgerfs.f"
/* L140: */
#line 442 "cgerfs.f"
    }

#line 444 "cgerfs.f"
    return 0;

/*     End of CGERFS */

} /* cgerfs_ */

