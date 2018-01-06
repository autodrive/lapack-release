#line 1 "dgtrfs.f"
/* dgtrfs.f -- translated by f2c (version 20100827).
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

#line 1 "dgtrfs.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b18 = -1.;
static doublereal c_b19 = 1.;

/* > \brief \b DGTRFS */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DGTRFS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgtrfs.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgtrfs.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgtrfs.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DGTRFS( TRANS, N, NRHS, DL, D, DU, DLF, DF, DUF, DU2, */
/*                          IPIV, B, LDB, X, LDX, FERR, BERR, WORK, IWORK, */
/*                          INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          TRANS */
/*       INTEGER            INFO, LDB, LDX, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ), IWORK( * ) */
/*       DOUBLE PRECISION   B( LDB, * ), BERR( * ), D( * ), DF( * ), */
/*      $                   DL( * ), DLF( * ), DU( * ), DU2( * ), DUF( * ), */
/*      $                   FERR( * ), WORK( * ), X( LDX, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DGTRFS improves the computed solution to a system of linear */
/* > equations when the coefficient matrix is tridiagonal, and provides */
/* > error bounds and backward error estimates for the solution. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] TRANS */
/* > \verbatim */
/* >          TRANS is CHARACTER*1 */
/* >          Specifies the form of the system of equations: */
/* >          = 'N':  A * X = B     (No transpose) */
/* >          = 'T':  A**T * X = B  (Transpose) */
/* >          = 'C':  A**H * X = B  (Conjugate transpose = Transpose) */
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
/* >          of the matrix B.  NRHS >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] DL */
/* > \verbatim */
/* >          DL is DOUBLE PRECISION array, dimension (N-1) */
/* >          The (n-1) subdiagonal elements of A. */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* >          D is DOUBLE PRECISION array, dimension (N) */
/* >          The diagonal elements of A. */
/* > \endverbatim */
/* > */
/* > \param[in] DU */
/* > \verbatim */
/* >          DU is DOUBLE PRECISION array, dimension (N-1) */
/* >          The (n-1) superdiagonal elements of A. */
/* > \endverbatim */
/* > */
/* > \param[in] DLF */
/* > \verbatim */
/* >          DLF is DOUBLE PRECISION array, dimension (N-1) */
/* >          The (n-1) multipliers that define the matrix L from the */
/* >          LU factorization of A as computed by DGTTRF. */
/* > \endverbatim */
/* > */
/* > \param[in] DF */
/* > \verbatim */
/* >          DF is DOUBLE PRECISION array, dimension (N) */
/* >          The n diagonal elements of the upper triangular matrix U from */
/* >          the LU factorization of A. */
/* > \endverbatim */
/* > */
/* > \param[in] DUF */
/* > \verbatim */
/* >          DUF is DOUBLE PRECISION array, dimension (N-1) */
/* >          The (n-1) elements of the first superdiagonal of U. */
/* > \endverbatim */
/* > */
/* > \param[in] DU2 */
/* > \verbatim */
/* >          DU2 is DOUBLE PRECISION array, dimension (N-2) */
/* >          The (n-2) elements of the second superdiagonal of U. */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N) */
/* >          The pivot indices; for 1 <= i <= n, row i of the matrix was */
/* >          interchanged with row IPIV(i).  IPIV(i) will always be either */
/* >          i or i+1; IPIV(i) = i indicates a row interchange was not */
/* >          required. */
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
/* > \param[in,out] X */
/* > \verbatim */
/* >          X is DOUBLE PRECISION array, dimension (LDX,NRHS) */
/* >          On entry, the solution matrix X, as computed by DGTTRS. */
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

/* > \ingroup doubleGTcomputational */

/*  ===================================================================== */
/* Subroutine */ int dgtrfs_(char *trans, integer *n, integer *nrhs, 
	doublereal *dl, doublereal *d__, doublereal *du, doublereal *dlf, 
	doublereal *df, doublereal *duf, doublereal *du2, integer *ipiv, 
	doublereal *b, integer *ldb, doublereal *x, integer *ldx, doublereal *
	ferr, doublereal *berr, doublereal *work, integer *iwork, integer *
	info, ftnlen trans_len)
{
    /* System generated locals */
    integer b_dim1, b_offset, x_dim1, x_offset, i__1, i__2;
    doublereal d__1, d__2, d__3, d__4;

    /* Local variables */
    static integer i__, j;
    static doublereal s;
    static integer nz;
    static doublereal eps;
    static integer kase;
    static doublereal safe1, safe2;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer isave[3];
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), daxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *);
    static integer count;
    extern /* Subroutine */ int dlacn2_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, integer *);
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int dlagtm_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    ftnlen);
    static doublereal safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static logical notran;
    static char transn[1];
    extern /* Subroutine */ int dgttrs_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     doublereal *, integer *, integer *, ftnlen);
    static char transt[1];
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

#line 264 "dgtrfs.f"
    /* Parameter adjustments */
#line 264 "dgtrfs.f"
    --dl;
#line 264 "dgtrfs.f"
    --d__;
#line 264 "dgtrfs.f"
    --du;
#line 264 "dgtrfs.f"
    --dlf;
#line 264 "dgtrfs.f"
    --df;
#line 264 "dgtrfs.f"
    --duf;
#line 264 "dgtrfs.f"
    --du2;
#line 264 "dgtrfs.f"
    --ipiv;
#line 264 "dgtrfs.f"
    b_dim1 = *ldb;
#line 264 "dgtrfs.f"
    b_offset = 1 + b_dim1;
#line 264 "dgtrfs.f"
    b -= b_offset;
#line 264 "dgtrfs.f"
    x_dim1 = *ldx;
#line 264 "dgtrfs.f"
    x_offset = 1 + x_dim1;
#line 264 "dgtrfs.f"
    x -= x_offset;
#line 264 "dgtrfs.f"
    --ferr;
#line 264 "dgtrfs.f"
    --berr;
#line 264 "dgtrfs.f"
    --work;
#line 264 "dgtrfs.f"
    --iwork;
#line 264 "dgtrfs.f"

#line 264 "dgtrfs.f"
    /* Function Body */
#line 264 "dgtrfs.f"
    *info = 0;
#line 265 "dgtrfs.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 266 "dgtrfs.f"
    if (! notran && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1) && ! lsame_(
	    trans, "C", (ftnlen)1, (ftnlen)1)) {
#line 268 "dgtrfs.f"
	*info = -1;
#line 269 "dgtrfs.f"
    } else if (*n < 0) {
#line 270 "dgtrfs.f"
	*info = -2;
#line 271 "dgtrfs.f"
    } else if (*nrhs < 0) {
#line 272 "dgtrfs.f"
	*info = -3;
#line 273 "dgtrfs.f"
    } else if (*ldb < max(1,*n)) {
#line 274 "dgtrfs.f"
	*info = -13;
#line 275 "dgtrfs.f"
    } else if (*ldx < max(1,*n)) {
#line 276 "dgtrfs.f"
	*info = -15;
#line 277 "dgtrfs.f"
    }
#line 278 "dgtrfs.f"
    if (*info != 0) {
#line 279 "dgtrfs.f"
	i__1 = -(*info);
#line 279 "dgtrfs.f"
	xerbla_("DGTRFS", &i__1, (ftnlen)6);
#line 280 "dgtrfs.f"
	return 0;
#line 281 "dgtrfs.f"
    }

/*     Quick return if possible */

#line 285 "dgtrfs.f"
    if (*n == 0 || *nrhs == 0) {
#line 286 "dgtrfs.f"
	i__1 = *nrhs;
#line 286 "dgtrfs.f"
	for (j = 1; j <= i__1; ++j) {
#line 287 "dgtrfs.f"
	    ferr[j] = 0.;
#line 288 "dgtrfs.f"
	    berr[j] = 0.;
#line 289 "dgtrfs.f"
/* L10: */
#line 289 "dgtrfs.f"
	}
#line 290 "dgtrfs.f"
	return 0;
#line 291 "dgtrfs.f"
    }

#line 293 "dgtrfs.f"
    if (notran) {
#line 294 "dgtrfs.f"
	*(unsigned char *)transn = 'N';
#line 295 "dgtrfs.f"
	*(unsigned char *)transt = 'T';
#line 296 "dgtrfs.f"
    } else {
#line 297 "dgtrfs.f"
	*(unsigned char *)transn = 'T';
#line 298 "dgtrfs.f"
	*(unsigned char *)transt = 'N';
#line 299 "dgtrfs.f"
    }

/*     NZ = maximum number of nonzero elements in each row of A, plus 1 */

#line 303 "dgtrfs.f"
    nz = 4;
#line 304 "dgtrfs.f"
    eps = dlamch_("Epsilon", (ftnlen)7);
#line 305 "dgtrfs.f"
    safmin = dlamch_("Safe minimum", (ftnlen)12);
#line 306 "dgtrfs.f"
    safe1 = nz * safmin;
#line 307 "dgtrfs.f"
    safe2 = safe1 / eps;

/*     Do for each right hand side */

#line 311 "dgtrfs.f"
    i__1 = *nrhs;
#line 311 "dgtrfs.f"
    for (j = 1; j <= i__1; ++j) {

#line 313 "dgtrfs.f"
	count = 1;
#line 314 "dgtrfs.f"
	lstres = 3.;
#line 315 "dgtrfs.f"
L20:

/*        Loop until stopping criterion is satisfied. */

/*        Compute residual R = B - op(A) * X, */
/*        where op(A) = A, A**T, or A**H, depending on TRANS. */

#line 322 "dgtrfs.f"
	dcopy_(n, &b[j * b_dim1 + 1], &c__1, &work[*n + 1], &c__1);
#line 323 "dgtrfs.f"
	dlagtm_(trans, n, &c__1, &c_b18, &dl[1], &d__[1], &du[1], &x[j * 
		x_dim1 + 1], ldx, &c_b19, &work[*n + 1], n, (ftnlen)1);

/*        Compute abs(op(A))*abs(x) + abs(b) for use in the backward */
/*        error bound. */

#line 329 "dgtrfs.f"
	if (notran) {
#line 330 "dgtrfs.f"
	    if (*n == 1) {
#line 331 "dgtrfs.f"
		work[1] = (d__1 = b[j * b_dim1 + 1], abs(d__1)) + (d__2 = d__[
			1] * x[j * x_dim1 + 1], abs(d__2));
#line 332 "dgtrfs.f"
	    } else {
#line 333 "dgtrfs.f"
		work[1] = (d__1 = b[j * b_dim1 + 1], abs(d__1)) + (d__2 = d__[
			1] * x[j * x_dim1 + 1], abs(d__2)) + (d__3 = du[1] * 
			x[j * x_dim1 + 2], abs(d__3));
#line 335 "dgtrfs.f"
		i__2 = *n - 1;
#line 335 "dgtrfs.f"
		for (i__ = 2; i__ <= i__2; ++i__) {
#line 336 "dgtrfs.f"
		    work[i__] = (d__1 = b[i__ + j * b_dim1], abs(d__1)) + (
			    d__2 = dl[i__ - 1] * x[i__ - 1 + j * x_dim1], abs(
			    d__2)) + (d__3 = d__[i__] * x[i__ + j * x_dim1], 
			    abs(d__3)) + (d__4 = du[i__] * x[i__ + 1 + j * 
			    x_dim1], abs(d__4));
#line 340 "dgtrfs.f"
/* L30: */
#line 340 "dgtrfs.f"
		}
#line 341 "dgtrfs.f"
		work[*n] = (d__1 = b[*n + j * b_dim1], abs(d__1)) + (d__2 = 
			dl[*n - 1] * x[*n - 1 + j * x_dim1], abs(d__2)) + (
			d__3 = d__[*n] * x[*n + j * x_dim1], abs(d__3));
#line 344 "dgtrfs.f"
	    }
#line 345 "dgtrfs.f"
	} else {
#line 346 "dgtrfs.f"
	    if (*n == 1) {
#line 347 "dgtrfs.f"
		work[1] = (d__1 = b[j * b_dim1 + 1], abs(d__1)) + (d__2 = d__[
			1] * x[j * x_dim1 + 1], abs(d__2));
#line 348 "dgtrfs.f"
	    } else {
#line 349 "dgtrfs.f"
		work[1] = (d__1 = b[j * b_dim1 + 1], abs(d__1)) + (d__2 = d__[
			1] * x[j * x_dim1 + 1], abs(d__2)) + (d__3 = dl[1] * 
			x[j * x_dim1 + 2], abs(d__3));
#line 351 "dgtrfs.f"
		i__2 = *n - 1;
#line 351 "dgtrfs.f"
		for (i__ = 2; i__ <= i__2; ++i__) {
#line 352 "dgtrfs.f"
		    work[i__] = (d__1 = b[i__ + j * b_dim1], abs(d__1)) + (
			    d__2 = du[i__ - 1] * x[i__ - 1 + j * x_dim1], abs(
			    d__2)) + (d__3 = d__[i__] * x[i__ + j * x_dim1], 
			    abs(d__3)) + (d__4 = dl[i__] * x[i__ + 1 + j * 
			    x_dim1], abs(d__4));
#line 356 "dgtrfs.f"
/* L40: */
#line 356 "dgtrfs.f"
		}
#line 357 "dgtrfs.f"
		work[*n] = (d__1 = b[*n + j * b_dim1], abs(d__1)) + (d__2 = 
			du[*n - 1] * x[*n - 1 + j * x_dim1], abs(d__2)) + (
			d__3 = d__[*n] * x[*n + j * x_dim1], abs(d__3));
#line 360 "dgtrfs.f"
	    }
#line 361 "dgtrfs.f"
	}

/*        Compute componentwise relative backward error from formula */

/*        max(i) ( abs(R(i)) / ( abs(op(A))*abs(X) + abs(B) )(i) ) */

/*        where abs(Z) is the componentwise absolute value of the matrix */
/*        or vector Z.  If the i-th component of the denominator is less */
/*        than SAFE2, then SAFE1 is added to the i-th components of the */
/*        numerator and denominator before dividing. */

#line 372 "dgtrfs.f"
	s = 0.;
#line 373 "dgtrfs.f"
	i__2 = *n;
#line 373 "dgtrfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 374 "dgtrfs.f"
	    if (work[i__] > safe2) {
/* Computing MAX */
#line 375 "dgtrfs.f"
		d__2 = s, d__3 = (d__1 = work[*n + i__], abs(d__1)) / work[
			i__];
#line 375 "dgtrfs.f"
		s = max(d__2,d__3);
#line 376 "dgtrfs.f"
	    } else {
/* Computing MAX */
#line 377 "dgtrfs.f"
		d__2 = s, d__3 = ((d__1 = work[*n + i__], abs(d__1)) + safe1) 
			/ (work[i__] + safe1);
#line 377 "dgtrfs.f"
		s = max(d__2,d__3);
#line 379 "dgtrfs.f"
	    }
#line 380 "dgtrfs.f"
/* L50: */
#line 380 "dgtrfs.f"
	}
#line 381 "dgtrfs.f"
	berr[j] = s;

/*        Test stopping criterion. Continue iterating if */
/*           1) The residual BERR(J) is larger than machine epsilon, and */
/*           2) BERR(J) decreased by at least a factor of 2 during the */
/*              last iteration, and */
/*           3) At most ITMAX iterations tried. */

#line 389 "dgtrfs.f"
	if (berr[j] > eps && berr[j] * 2. <= lstres && count <= 5) {

/*           Update solution and try again. */

#line 394 "dgtrfs.f"
	    dgttrs_(trans, n, &c__1, &dlf[1], &df[1], &duf[1], &du2[1], &ipiv[
		    1], &work[*n + 1], n, info, (ftnlen)1);
#line 396 "dgtrfs.f"
	    daxpy_(n, &c_b19, &work[*n + 1], &c__1, &x[j * x_dim1 + 1], &c__1)
		    ;
#line 397 "dgtrfs.f"
	    lstres = berr[j];
#line 398 "dgtrfs.f"
	    ++count;
#line 399 "dgtrfs.f"
	    goto L20;
#line 400 "dgtrfs.f"
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

/*        Use DLACN2 to estimate the infinity-norm of the matrix */
/*           inv(op(A)) * diag(W), */
/*        where W = abs(R) + NZ*EPS*( abs(op(A))*abs(X)+abs(B) ))) */

#line 424 "dgtrfs.f"
	i__2 = *n;
#line 424 "dgtrfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 425 "dgtrfs.f"
	    if (work[i__] > safe2) {
#line 426 "dgtrfs.f"
		work[i__] = (d__1 = work[*n + i__], abs(d__1)) + nz * eps * 
			work[i__];
#line 427 "dgtrfs.f"
	    } else {
#line 428 "dgtrfs.f"
		work[i__] = (d__1 = work[*n + i__], abs(d__1)) + nz * eps * 
			work[i__] + safe1;
#line 429 "dgtrfs.f"
	    }
#line 430 "dgtrfs.f"
/* L60: */
#line 430 "dgtrfs.f"
	}

#line 432 "dgtrfs.f"
	kase = 0;
#line 433 "dgtrfs.f"
L70:
#line 434 "dgtrfs.f"
	dlacn2_(n, &work[(*n << 1) + 1], &work[*n + 1], &iwork[1], &ferr[j], &
		kase, isave);
#line 436 "dgtrfs.f"
	if (kase != 0) {
#line 437 "dgtrfs.f"
	    if (kase == 1) {

/*              Multiply by diag(W)*inv(op(A)**T). */

#line 441 "dgtrfs.f"
		dgttrs_(transt, n, &c__1, &dlf[1], &df[1], &duf[1], &du2[1], &
			ipiv[1], &work[*n + 1], n, info, (ftnlen)1);
#line 443 "dgtrfs.f"
		i__2 = *n;
#line 443 "dgtrfs.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 444 "dgtrfs.f"
		    work[*n + i__] = work[i__] * work[*n + i__];
#line 445 "dgtrfs.f"
/* L80: */
#line 445 "dgtrfs.f"
		}
#line 446 "dgtrfs.f"
	    } else {

/*              Multiply by inv(op(A))*diag(W). */

#line 450 "dgtrfs.f"
		i__2 = *n;
#line 450 "dgtrfs.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 451 "dgtrfs.f"
		    work[*n + i__] = work[i__] * work[*n + i__];
#line 452 "dgtrfs.f"
/* L90: */
#line 452 "dgtrfs.f"
		}
#line 453 "dgtrfs.f"
		dgttrs_(transn, n, &c__1, &dlf[1], &df[1], &duf[1], &du2[1], &
			ipiv[1], &work[*n + 1], n, info, (ftnlen)1);
#line 455 "dgtrfs.f"
	    }
#line 456 "dgtrfs.f"
	    goto L70;
#line 457 "dgtrfs.f"
	}

/*        Normalize error. */

#line 461 "dgtrfs.f"
	lstres = 0.;
#line 462 "dgtrfs.f"
	i__2 = *n;
#line 462 "dgtrfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 463 "dgtrfs.f"
	    d__2 = lstres, d__3 = (d__1 = x[i__ + j * x_dim1], abs(d__1));
#line 463 "dgtrfs.f"
	    lstres = max(d__2,d__3);
#line 464 "dgtrfs.f"
/* L100: */
#line 464 "dgtrfs.f"
	}
#line 465 "dgtrfs.f"
	if (lstres != 0.) {
#line 465 "dgtrfs.f"
	    ferr[j] /= lstres;
#line 465 "dgtrfs.f"
	}

#line 468 "dgtrfs.f"
/* L110: */
#line 468 "dgtrfs.f"
    }

#line 470 "dgtrfs.f"
    return 0;

/*     End of DGTRFS */

} /* dgtrfs_ */

