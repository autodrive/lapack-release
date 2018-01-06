#line 1 "sgtrfs.f"
/* sgtrfs.f -- translated by f2c (version 20100827).
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

#line 1 "sgtrfs.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b18 = -1.;
static doublereal c_b19 = 1.;

/* > \brief \b SGTRFS */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SGTRFS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgtrfs.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgtrfs.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgtrfs.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SGTRFS( TRANS, N, NRHS, DL, D, DU, DLF, DF, DUF, DU2, */
/*                          IPIV, B, LDB, X, LDX, FERR, BERR, WORK, IWORK, */
/*                          INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          TRANS */
/*       INTEGER            INFO, LDB, LDX, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ), IWORK( * ) */
/*       REAL               B( LDB, * ), BERR( * ), D( * ), DF( * ), */
/*      $                   DL( * ), DLF( * ), DU( * ), DU2( * ), DUF( * ), */
/*      $                   FERR( * ), WORK( * ), X( LDX, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SGTRFS improves the computed solution to a system of linear */
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
/* >          DL is REAL array, dimension (N-1) */
/* >          The (n-1) subdiagonal elements of A. */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* >          D is REAL array, dimension (N) */
/* >          The diagonal elements of A. */
/* > \endverbatim */
/* > */
/* > \param[in] DU */
/* > \verbatim */
/* >          DU is REAL array, dimension (N-1) */
/* >          The (n-1) superdiagonal elements of A. */
/* > \endverbatim */
/* > */
/* > \param[in] DLF */
/* > \verbatim */
/* >          DLF is REAL array, dimension (N-1) */
/* >          The (n-1) multipliers that define the matrix L from the */
/* >          LU factorization of A as computed by SGTTRF. */
/* > \endverbatim */
/* > */
/* > \param[in] DF */
/* > \verbatim */
/* >          DF is REAL array, dimension (N) */
/* >          The n diagonal elements of the upper triangular matrix U from */
/* >          the LU factorization of A. */
/* > \endverbatim */
/* > */
/* > \param[in] DUF */
/* > \verbatim */
/* >          DUF is REAL array, dimension (N-1) */
/* >          The (n-1) elements of the first superdiagonal of U. */
/* > \endverbatim */
/* > */
/* > \param[in] DU2 */
/* > \verbatim */
/* >          DU2 is REAL array, dimension (N-2) */
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
/* >          B is REAL array, dimension (LDB,NRHS) */
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
/* >          X is REAL array, dimension (LDX,NRHS) */
/* >          On entry, the solution matrix X, as computed by SGTTRS. */
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
/* >          WORK is REAL array, dimension (3*N) */
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

/* > \ingroup realGTcomputational */

/*  ===================================================================== */
/* Subroutine */ int sgtrfs_(char *trans, integer *n, integer *nrhs, 
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
    static integer isave[3], count;
    extern /* Subroutine */ int scopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), saxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *), slacn2_(
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, integer *);
    extern doublereal slamch_(char *, ftnlen);
    static doublereal safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), slagtm_(
	    char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     doublereal *, integer *, ftnlen);
    static logical notran;
    static char transn[1], transt[1];
    static doublereal lstres;
    extern /* Subroutine */ int sgttrs_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     doublereal *, integer *, integer *, ftnlen);


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

#line 264 "sgtrfs.f"
    /* Parameter adjustments */
#line 264 "sgtrfs.f"
    --dl;
#line 264 "sgtrfs.f"
    --d__;
#line 264 "sgtrfs.f"
    --du;
#line 264 "sgtrfs.f"
    --dlf;
#line 264 "sgtrfs.f"
    --df;
#line 264 "sgtrfs.f"
    --duf;
#line 264 "sgtrfs.f"
    --du2;
#line 264 "sgtrfs.f"
    --ipiv;
#line 264 "sgtrfs.f"
    b_dim1 = *ldb;
#line 264 "sgtrfs.f"
    b_offset = 1 + b_dim1;
#line 264 "sgtrfs.f"
    b -= b_offset;
#line 264 "sgtrfs.f"
    x_dim1 = *ldx;
#line 264 "sgtrfs.f"
    x_offset = 1 + x_dim1;
#line 264 "sgtrfs.f"
    x -= x_offset;
#line 264 "sgtrfs.f"
    --ferr;
#line 264 "sgtrfs.f"
    --berr;
#line 264 "sgtrfs.f"
    --work;
#line 264 "sgtrfs.f"
    --iwork;
#line 264 "sgtrfs.f"

#line 264 "sgtrfs.f"
    /* Function Body */
#line 264 "sgtrfs.f"
    *info = 0;
#line 265 "sgtrfs.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 266 "sgtrfs.f"
    if (! notran && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1) && ! lsame_(
	    trans, "C", (ftnlen)1, (ftnlen)1)) {
#line 268 "sgtrfs.f"
	*info = -1;
#line 269 "sgtrfs.f"
    } else if (*n < 0) {
#line 270 "sgtrfs.f"
	*info = -2;
#line 271 "sgtrfs.f"
    } else if (*nrhs < 0) {
#line 272 "sgtrfs.f"
	*info = -3;
#line 273 "sgtrfs.f"
    } else if (*ldb < max(1,*n)) {
#line 274 "sgtrfs.f"
	*info = -13;
#line 275 "sgtrfs.f"
    } else if (*ldx < max(1,*n)) {
#line 276 "sgtrfs.f"
	*info = -15;
#line 277 "sgtrfs.f"
    }
#line 278 "sgtrfs.f"
    if (*info != 0) {
#line 279 "sgtrfs.f"
	i__1 = -(*info);
#line 279 "sgtrfs.f"
	xerbla_("SGTRFS", &i__1, (ftnlen)6);
#line 280 "sgtrfs.f"
	return 0;
#line 281 "sgtrfs.f"
    }

/*     Quick return if possible */

#line 285 "sgtrfs.f"
    if (*n == 0 || *nrhs == 0) {
#line 286 "sgtrfs.f"
	i__1 = *nrhs;
#line 286 "sgtrfs.f"
	for (j = 1; j <= i__1; ++j) {
#line 287 "sgtrfs.f"
	    ferr[j] = 0.;
#line 288 "sgtrfs.f"
	    berr[j] = 0.;
#line 289 "sgtrfs.f"
/* L10: */
#line 289 "sgtrfs.f"
	}
#line 290 "sgtrfs.f"
	return 0;
#line 291 "sgtrfs.f"
    }

#line 293 "sgtrfs.f"
    if (notran) {
#line 294 "sgtrfs.f"
	*(unsigned char *)transn = 'N';
#line 295 "sgtrfs.f"
	*(unsigned char *)transt = 'T';
#line 296 "sgtrfs.f"
    } else {
#line 297 "sgtrfs.f"
	*(unsigned char *)transn = 'T';
#line 298 "sgtrfs.f"
	*(unsigned char *)transt = 'N';
#line 299 "sgtrfs.f"
    }

/*     NZ = maximum number of nonzero elements in each row of A, plus 1 */

#line 303 "sgtrfs.f"
    nz = 4;
#line 304 "sgtrfs.f"
    eps = slamch_("Epsilon", (ftnlen)7);
#line 305 "sgtrfs.f"
    safmin = slamch_("Safe minimum", (ftnlen)12);
#line 306 "sgtrfs.f"
    safe1 = nz * safmin;
#line 307 "sgtrfs.f"
    safe2 = safe1 / eps;

/*     Do for each right hand side */

#line 311 "sgtrfs.f"
    i__1 = *nrhs;
#line 311 "sgtrfs.f"
    for (j = 1; j <= i__1; ++j) {

#line 313 "sgtrfs.f"
	count = 1;
#line 314 "sgtrfs.f"
	lstres = 3.;
#line 315 "sgtrfs.f"
L20:

/*        Loop until stopping criterion is satisfied. */

/*        Compute residual R = B - op(A) * X, */
/*        where op(A) = A, A**T, or A**H, depending on TRANS. */

#line 322 "sgtrfs.f"
	scopy_(n, &b[j * b_dim1 + 1], &c__1, &work[*n + 1], &c__1);
#line 323 "sgtrfs.f"
	slagtm_(trans, n, &c__1, &c_b18, &dl[1], &d__[1], &du[1], &x[j * 
		x_dim1 + 1], ldx, &c_b19, &work[*n + 1], n, (ftnlen)1);

/*        Compute abs(op(A))*abs(x) + abs(b) for use in the backward */
/*        error bound. */

#line 329 "sgtrfs.f"
	if (notran) {
#line 330 "sgtrfs.f"
	    if (*n == 1) {
#line 331 "sgtrfs.f"
		work[1] = (d__1 = b[j * b_dim1 + 1], abs(d__1)) + (d__2 = d__[
			1] * x[j * x_dim1 + 1], abs(d__2));
#line 332 "sgtrfs.f"
	    } else {
#line 333 "sgtrfs.f"
		work[1] = (d__1 = b[j * b_dim1 + 1], abs(d__1)) + (d__2 = d__[
			1] * x[j * x_dim1 + 1], abs(d__2)) + (d__3 = du[1] * 
			x[j * x_dim1 + 2], abs(d__3));
#line 335 "sgtrfs.f"
		i__2 = *n - 1;
#line 335 "sgtrfs.f"
		for (i__ = 2; i__ <= i__2; ++i__) {
#line 336 "sgtrfs.f"
		    work[i__] = (d__1 = b[i__ + j * b_dim1], abs(d__1)) + (
			    d__2 = dl[i__ - 1] * x[i__ - 1 + j * x_dim1], abs(
			    d__2)) + (d__3 = d__[i__] * x[i__ + j * x_dim1], 
			    abs(d__3)) + (d__4 = du[i__] * x[i__ + 1 + j * 
			    x_dim1], abs(d__4));
#line 340 "sgtrfs.f"
/* L30: */
#line 340 "sgtrfs.f"
		}
#line 341 "sgtrfs.f"
		work[*n] = (d__1 = b[*n + j * b_dim1], abs(d__1)) + (d__2 = 
			dl[*n - 1] * x[*n - 1 + j * x_dim1], abs(d__2)) + (
			d__3 = d__[*n] * x[*n + j * x_dim1], abs(d__3));
#line 344 "sgtrfs.f"
	    }
#line 345 "sgtrfs.f"
	} else {
#line 346 "sgtrfs.f"
	    if (*n == 1) {
#line 347 "sgtrfs.f"
		work[1] = (d__1 = b[j * b_dim1 + 1], abs(d__1)) + (d__2 = d__[
			1] * x[j * x_dim1 + 1], abs(d__2));
#line 348 "sgtrfs.f"
	    } else {
#line 349 "sgtrfs.f"
		work[1] = (d__1 = b[j * b_dim1 + 1], abs(d__1)) + (d__2 = d__[
			1] * x[j * x_dim1 + 1], abs(d__2)) + (d__3 = dl[1] * 
			x[j * x_dim1 + 2], abs(d__3));
#line 351 "sgtrfs.f"
		i__2 = *n - 1;
#line 351 "sgtrfs.f"
		for (i__ = 2; i__ <= i__2; ++i__) {
#line 352 "sgtrfs.f"
		    work[i__] = (d__1 = b[i__ + j * b_dim1], abs(d__1)) + (
			    d__2 = du[i__ - 1] * x[i__ - 1 + j * x_dim1], abs(
			    d__2)) + (d__3 = d__[i__] * x[i__ + j * x_dim1], 
			    abs(d__3)) + (d__4 = dl[i__] * x[i__ + 1 + j * 
			    x_dim1], abs(d__4));
#line 356 "sgtrfs.f"
/* L40: */
#line 356 "sgtrfs.f"
		}
#line 357 "sgtrfs.f"
		work[*n] = (d__1 = b[*n + j * b_dim1], abs(d__1)) + (d__2 = 
			du[*n - 1] * x[*n - 1 + j * x_dim1], abs(d__2)) + (
			d__3 = d__[*n] * x[*n + j * x_dim1], abs(d__3));
#line 360 "sgtrfs.f"
	    }
#line 361 "sgtrfs.f"
	}

/*        Compute componentwise relative backward error from formula */

/*        max(i) ( abs(R(i)) / ( abs(op(A))*abs(X) + abs(B) )(i) ) */

/*        where abs(Z) is the componentwise absolute value of the matrix */
/*        or vector Z.  If the i-th component of the denominator is less */
/*        than SAFE2, then SAFE1 is added to the i-th components of the */
/*        numerator and denominator before dividing. */

#line 372 "sgtrfs.f"
	s = 0.;
#line 373 "sgtrfs.f"
	i__2 = *n;
#line 373 "sgtrfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 374 "sgtrfs.f"
	    if (work[i__] > safe2) {
/* Computing MAX */
#line 375 "sgtrfs.f"
		d__2 = s, d__3 = (d__1 = work[*n + i__], abs(d__1)) / work[
			i__];
#line 375 "sgtrfs.f"
		s = max(d__2,d__3);
#line 376 "sgtrfs.f"
	    } else {
/* Computing MAX */
#line 377 "sgtrfs.f"
		d__2 = s, d__3 = ((d__1 = work[*n + i__], abs(d__1)) + safe1) 
			/ (work[i__] + safe1);
#line 377 "sgtrfs.f"
		s = max(d__2,d__3);
#line 379 "sgtrfs.f"
	    }
#line 380 "sgtrfs.f"
/* L50: */
#line 380 "sgtrfs.f"
	}
#line 381 "sgtrfs.f"
	berr[j] = s;

/*        Test stopping criterion. Continue iterating if */
/*           1) The residual BERR(J) is larger than machine epsilon, and */
/*           2) BERR(J) decreased by at least a factor of 2 during the */
/*              last iteration, and */
/*           3) At most ITMAX iterations tried. */

#line 389 "sgtrfs.f"
	if (berr[j] > eps && berr[j] * 2. <= lstres && count <= 5) {

/*           Update solution and try again. */

#line 394 "sgtrfs.f"
	    sgttrs_(trans, n, &c__1, &dlf[1], &df[1], &duf[1], &du2[1], &ipiv[
		    1], &work[*n + 1], n, info, (ftnlen)1);
#line 396 "sgtrfs.f"
	    saxpy_(n, &c_b19, &work[*n + 1], &c__1, &x[j * x_dim1 + 1], &c__1)
		    ;
#line 397 "sgtrfs.f"
	    lstres = berr[j];
#line 398 "sgtrfs.f"
	    ++count;
#line 399 "sgtrfs.f"
	    goto L20;
#line 400 "sgtrfs.f"
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

/*        Use SLACN2 to estimate the infinity-norm of the matrix */
/*           inv(op(A)) * diag(W), */
/*        where W = abs(R) + NZ*EPS*( abs(op(A))*abs(X)+abs(B) ))) */

#line 424 "sgtrfs.f"
	i__2 = *n;
#line 424 "sgtrfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 425 "sgtrfs.f"
	    if (work[i__] > safe2) {
#line 426 "sgtrfs.f"
		work[i__] = (d__1 = work[*n + i__], abs(d__1)) + nz * eps * 
			work[i__];
#line 427 "sgtrfs.f"
	    } else {
#line 428 "sgtrfs.f"
		work[i__] = (d__1 = work[*n + i__], abs(d__1)) + nz * eps * 
			work[i__] + safe1;
#line 429 "sgtrfs.f"
	    }
#line 430 "sgtrfs.f"
/* L60: */
#line 430 "sgtrfs.f"
	}

#line 432 "sgtrfs.f"
	kase = 0;
#line 433 "sgtrfs.f"
L70:
#line 434 "sgtrfs.f"
	slacn2_(n, &work[(*n << 1) + 1], &work[*n + 1], &iwork[1], &ferr[j], &
		kase, isave);
#line 436 "sgtrfs.f"
	if (kase != 0) {
#line 437 "sgtrfs.f"
	    if (kase == 1) {

/*              Multiply by diag(W)*inv(op(A)**T). */

#line 441 "sgtrfs.f"
		sgttrs_(transt, n, &c__1, &dlf[1], &df[1], &duf[1], &du2[1], &
			ipiv[1], &work[*n + 1], n, info, (ftnlen)1);
#line 443 "sgtrfs.f"
		i__2 = *n;
#line 443 "sgtrfs.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 444 "sgtrfs.f"
		    work[*n + i__] = work[i__] * work[*n + i__];
#line 445 "sgtrfs.f"
/* L80: */
#line 445 "sgtrfs.f"
		}
#line 446 "sgtrfs.f"
	    } else {

/*              Multiply by inv(op(A))*diag(W). */

#line 450 "sgtrfs.f"
		i__2 = *n;
#line 450 "sgtrfs.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 451 "sgtrfs.f"
		    work[*n + i__] = work[i__] * work[*n + i__];
#line 452 "sgtrfs.f"
/* L90: */
#line 452 "sgtrfs.f"
		}
#line 453 "sgtrfs.f"
		sgttrs_(transn, n, &c__1, &dlf[1], &df[1], &duf[1], &du2[1], &
			ipiv[1], &work[*n + 1], n, info, (ftnlen)1);
#line 455 "sgtrfs.f"
	    }
#line 456 "sgtrfs.f"
	    goto L70;
#line 457 "sgtrfs.f"
	}

/*        Normalize error. */

#line 461 "sgtrfs.f"
	lstres = 0.;
#line 462 "sgtrfs.f"
	i__2 = *n;
#line 462 "sgtrfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 463 "sgtrfs.f"
	    d__2 = lstres, d__3 = (d__1 = x[i__ + j * x_dim1], abs(d__1));
#line 463 "sgtrfs.f"
	    lstres = max(d__2,d__3);
#line 464 "sgtrfs.f"
/* L100: */
#line 464 "sgtrfs.f"
	}
#line 465 "sgtrfs.f"
	if (lstres != 0.) {
#line 465 "sgtrfs.f"
	    ferr[j] /= lstres;
#line 465 "sgtrfs.f"
	}

#line 468 "sgtrfs.f"
/* L110: */
#line 468 "sgtrfs.f"
    }

#line 470 "sgtrfs.f"
    return 0;

/*     End of SGTRFS */

} /* sgtrfs_ */

