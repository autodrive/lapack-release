#line 1 "dgbrfs.f"
/* dgbrfs.f -- translated by f2c (version 20100827).
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

#line 1 "dgbrfs.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b15 = -1.;
static doublereal c_b17 = 1.;

/* > \brief \b DGBRFS */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DGBRFS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgbrfs.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgbrfs.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgbrfs.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DGBRFS( TRANS, N, KL, KU, NRHS, AB, LDAB, AFB, LDAFB, */
/*                          IPIV, B, LDB, X, LDX, FERR, BERR, WORK, IWORK, */
/*                          INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          TRANS */
/*       INTEGER            INFO, KL, KU, LDAB, LDAFB, LDB, LDX, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ), IWORK( * ) */
/*       DOUBLE PRECISION   AB( LDAB, * ), AFB( LDAFB, * ), B( LDB, * ), */
/*      $                   BERR( * ), FERR( * ), WORK( * ), X( LDX, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DGBRFS improves the computed solution to a system of linear */
/* > equations when the coefficient matrix is banded, and provides */
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
/* > \param[in] KL */
/* > \verbatim */
/* >          KL is INTEGER */
/* >          The number of subdiagonals within the band of A.  KL >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] KU */
/* > \verbatim */
/* >          KU is INTEGER */
/* >          The number of superdiagonals within the band of A.  KU >= 0. */
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
/* >          AB is DOUBLE PRECISION array, dimension (LDAB,N) */
/* >          The original band matrix A, stored in rows 1 to KL+KU+1. */
/* >          The j-th column of A is stored in the j-th column of the */
/* >          array AB as follows: */
/* >          AB(ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(n,j+kl). */
/* > \endverbatim */
/* > */
/* > \param[in] LDAB */
/* > \verbatim */
/* >          LDAB is INTEGER */
/* >          The leading dimension of the array AB.  LDAB >= KL+KU+1. */
/* > \endverbatim */
/* > */
/* > \param[in] AFB */
/* > \verbatim */
/* >          AFB is DOUBLE PRECISION array, dimension (LDAFB,N) */
/* >          Details of the LU factorization of the band matrix A, as */
/* >          computed by DGBTRF.  U is stored as an upper triangular band */
/* >          matrix with KL+KU superdiagonals in rows 1 to KL+KU+1, and */
/* >          the multipliers used during the factorization are stored in */
/* >          rows KL+KU+2 to 2*KL+KU+1. */
/* > \endverbatim */
/* > */
/* > \param[in] LDAFB */
/* > \verbatim */
/* >          LDAFB is INTEGER */
/* >          The leading dimension of the array AFB.  LDAFB >= 2*KL*KU+1. */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N) */
/* >          The pivot indices from DGBTRF; for 1<=i<=N, row i of the */
/* >          matrix was interchanged with row IPIV(i). */
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
/* >          On entry, the solution matrix X, as computed by DGBTRS. */
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

/* > \date November 2011 */

/* > \ingroup doubleGBcomputational */

/*  ===================================================================== */
/* Subroutine */ int dgbrfs_(char *trans, integer *n, integer *kl, integer *
	ku, integer *nrhs, doublereal *ab, integer *ldab, doublereal *afb, 
	integer *ldafb, integer *ipiv, doublereal *b, integer *ldb, 
	doublereal *x, integer *ldx, doublereal *ferr, doublereal *berr, 
	doublereal *work, integer *iwork, integer *info, ftnlen trans_len)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, afb_dim1, afb_offset, b_dim1, b_offset, 
	    x_dim1, x_offset, i__1, i__2, i__3, i__4, i__5, i__6, i__7;
    doublereal d__1, d__2, d__3;

    /* Local variables */
    static integer i__, j, k;
    static doublereal s;
    static integer kk;
    static doublereal xk;
    static integer nz;
    static doublereal eps;
    static integer kase;
    static doublereal safe1, safe2;
    extern /* Subroutine */ int dgbmv_(char *, integer *, integer *, integer *
	    , integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer isave[3];
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), daxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *);
    static integer count;
    extern /* Subroutine */ int dlacn2_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, integer *);
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), dgbtrs_(
	    char *, integer *, integer *, integer *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, integer *, integer *, ftnlen);
    static logical notran;
    static char transt[1];
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
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 261 "dgbrfs.f"
    /* Parameter adjustments */
#line 261 "dgbrfs.f"
    ab_dim1 = *ldab;
#line 261 "dgbrfs.f"
    ab_offset = 1 + ab_dim1;
#line 261 "dgbrfs.f"
    ab -= ab_offset;
#line 261 "dgbrfs.f"
    afb_dim1 = *ldafb;
#line 261 "dgbrfs.f"
    afb_offset = 1 + afb_dim1;
#line 261 "dgbrfs.f"
    afb -= afb_offset;
#line 261 "dgbrfs.f"
    --ipiv;
#line 261 "dgbrfs.f"
    b_dim1 = *ldb;
#line 261 "dgbrfs.f"
    b_offset = 1 + b_dim1;
#line 261 "dgbrfs.f"
    b -= b_offset;
#line 261 "dgbrfs.f"
    x_dim1 = *ldx;
#line 261 "dgbrfs.f"
    x_offset = 1 + x_dim1;
#line 261 "dgbrfs.f"
    x -= x_offset;
#line 261 "dgbrfs.f"
    --ferr;
#line 261 "dgbrfs.f"
    --berr;
#line 261 "dgbrfs.f"
    --work;
#line 261 "dgbrfs.f"
    --iwork;
#line 261 "dgbrfs.f"

#line 261 "dgbrfs.f"
    /* Function Body */
#line 261 "dgbrfs.f"
    *info = 0;
#line 262 "dgbrfs.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 263 "dgbrfs.f"
    if (! notran && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1) && ! lsame_(
	    trans, "C", (ftnlen)1, (ftnlen)1)) {
#line 265 "dgbrfs.f"
	*info = -1;
#line 266 "dgbrfs.f"
    } else if (*n < 0) {
#line 267 "dgbrfs.f"
	*info = -2;
#line 268 "dgbrfs.f"
    } else if (*kl < 0) {
#line 269 "dgbrfs.f"
	*info = -3;
#line 270 "dgbrfs.f"
    } else if (*ku < 0) {
#line 271 "dgbrfs.f"
	*info = -4;
#line 272 "dgbrfs.f"
    } else if (*nrhs < 0) {
#line 273 "dgbrfs.f"
	*info = -5;
#line 274 "dgbrfs.f"
    } else if (*ldab < *kl + *ku + 1) {
#line 275 "dgbrfs.f"
	*info = -7;
#line 276 "dgbrfs.f"
    } else if (*ldafb < (*kl << 1) + *ku + 1) {
#line 277 "dgbrfs.f"
	*info = -9;
#line 278 "dgbrfs.f"
    } else if (*ldb < max(1,*n)) {
#line 279 "dgbrfs.f"
	*info = -12;
#line 280 "dgbrfs.f"
    } else if (*ldx < max(1,*n)) {
#line 281 "dgbrfs.f"
	*info = -14;
#line 282 "dgbrfs.f"
    }
#line 283 "dgbrfs.f"
    if (*info != 0) {
#line 284 "dgbrfs.f"
	i__1 = -(*info);
#line 284 "dgbrfs.f"
	xerbla_("DGBRFS", &i__1, (ftnlen)6);
#line 285 "dgbrfs.f"
	return 0;
#line 286 "dgbrfs.f"
    }

/*     Quick return if possible */

#line 290 "dgbrfs.f"
    if (*n == 0 || *nrhs == 0) {
#line 291 "dgbrfs.f"
	i__1 = *nrhs;
#line 291 "dgbrfs.f"
	for (j = 1; j <= i__1; ++j) {
#line 292 "dgbrfs.f"
	    ferr[j] = 0.;
#line 293 "dgbrfs.f"
	    berr[j] = 0.;
#line 294 "dgbrfs.f"
/* L10: */
#line 294 "dgbrfs.f"
	}
#line 295 "dgbrfs.f"
	return 0;
#line 296 "dgbrfs.f"
    }

#line 298 "dgbrfs.f"
    if (notran) {
#line 299 "dgbrfs.f"
	*(unsigned char *)transt = 'T';
#line 300 "dgbrfs.f"
    } else {
#line 301 "dgbrfs.f"
	*(unsigned char *)transt = 'N';
#line 302 "dgbrfs.f"
    }

/*     NZ = maximum number of nonzero elements in each row of A, plus 1 */

/* Computing MIN */
#line 306 "dgbrfs.f"
    i__1 = *kl + *ku + 2, i__2 = *n + 1;
#line 306 "dgbrfs.f"
    nz = min(i__1,i__2);
#line 307 "dgbrfs.f"
    eps = dlamch_("Epsilon", (ftnlen)7);
#line 308 "dgbrfs.f"
    safmin = dlamch_("Safe minimum", (ftnlen)12);
#line 309 "dgbrfs.f"
    safe1 = nz * safmin;
#line 310 "dgbrfs.f"
    safe2 = safe1 / eps;

/*     Do for each right hand side */

#line 314 "dgbrfs.f"
    i__1 = *nrhs;
#line 314 "dgbrfs.f"
    for (j = 1; j <= i__1; ++j) {

#line 316 "dgbrfs.f"
	count = 1;
#line 317 "dgbrfs.f"
	lstres = 3.;
#line 318 "dgbrfs.f"
L20:

/*        Loop until stopping criterion is satisfied. */

/*        Compute residual R = B - op(A) * X, */
/*        where op(A) = A, A**T, or A**H, depending on TRANS. */

#line 325 "dgbrfs.f"
	dcopy_(n, &b[j * b_dim1 + 1], &c__1, &work[*n + 1], &c__1);
#line 326 "dgbrfs.f"
	dgbmv_(trans, n, n, kl, ku, &c_b15, &ab[ab_offset], ldab, &x[j * 
		x_dim1 + 1], &c__1, &c_b17, &work[*n + 1], &c__1, (ftnlen)1);

/*        Compute componentwise relative backward error from formula */

/*        max(i) ( abs(R(i)) / ( abs(op(A))*abs(X) + abs(B) )(i) ) */

/*        where abs(Z) is the componentwise absolute value of the matrix */
/*        or vector Z.  If the i-th component of the denominator is less */
/*        than SAFE2, then SAFE1 is added to the i-th components of the */
/*        numerator and denominator before dividing. */

#line 338 "dgbrfs.f"
	i__2 = *n;
#line 338 "dgbrfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 339 "dgbrfs.f"
	    work[i__] = (d__1 = b[i__ + j * b_dim1], abs(d__1));
#line 340 "dgbrfs.f"
/* L30: */
#line 340 "dgbrfs.f"
	}

/*        Compute abs(op(A))*abs(X) + abs(B). */

#line 344 "dgbrfs.f"
	if (notran) {
#line 345 "dgbrfs.f"
	    i__2 = *n;
#line 345 "dgbrfs.f"
	    for (k = 1; k <= i__2; ++k) {
#line 346 "dgbrfs.f"
		kk = *ku + 1 - k;
#line 347 "dgbrfs.f"
		xk = (d__1 = x[k + j * x_dim1], abs(d__1));
/* Computing MAX */
#line 348 "dgbrfs.f"
		i__3 = 1, i__4 = k - *ku;
/* Computing MIN */
#line 348 "dgbrfs.f"
		i__6 = *n, i__7 = k + *kl;
#line 348 "dgbrfs.f"
		i__5 = min(i__6,i__7);
#line 348 "dgbrfs.f"
		for (i__ = max(i__3,i__4); i__ <= i__5; ++i__) {
#line 349 "dgbrfs.f"
		    work[i__] += (d__1 = ab[kk + i__ + k * ab_dim1], abs(d__1)
			    ) * xk;
#line 350 "dgbrfs.f"
/* L40: */
#line 350 "dgbrfs.f"
		}
#line 351 "dgbrfs.f"
/* L50: */
#line 351 "dgbrfs.f"
	    }
#line 352 "dgbrfs.f"
	} else {
#line 353 "dgbrfs.f"
	    i__2 = *n;
#line 353 "dgbrfs.f"
	    for (k = 1; k <= i__2; ++k) {
#line 354 "dgbrfs.f"
		s = 0.;
#line 355 "dgbrfs.f"
		kk = *ku + 1 - k;
/* Computing MAX */
#line 356 "dgbrfs.f"
		i__5 = 1, i__3 = k - *ku;
/* Computing MIN */
#line 356 "dgbrfs.f"
		i__6 = *n, i__7 = k + *kl;
#line 356 "dgbrfs.f"
		i__4 = min(i__6,i__7);
#line 356 "dgbrfs.f"
		for (i__ = max(i__5,i__3); i__ <= i__4; ++i__) {
#line 357 "dgbrfs.f"
		    s += (d__1 = ab[kk + i__ + k * ab_dim1], abs(d__1)) * (
			    d__2 = x[i__ + j * x_dim1], abs(d__2));
#line 358 "dgbrfs.f"
/* L60: */
#line 358 "dgbrfs.f"
		}
#line 359 "dgbrfs.f"
		work[k] += s;
#line 360 "dgbrfs.f"
/* L70: */
#line 360 "dgbrfs.f"
	    }
#line 361 "dgbrfs.f"
	}
#line 362 "dgbrfs.f"
	s = 0.;
#line 363 "dgbrfs.f"
	i__2 = *n;
#line 363 "dgbrfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 364 "dgbrfs.f"
	    if (work[i__] > safe2) {
/* Computing MAX */
#line 365 "dgbrfs.f"
		d__2 = s, d__3 = (d__1 = work[*n + i__], abs(d__1)) / work[
			i__];
#line 365 "dgbrfs.f"
		s = max(d__2,d__3);
#line 366 "dgbrfs.f"
	    } else {
/* Computing MAX */
#line 367 "dgbrfs.f"
		d__2 = s, d__3 = ((d__1 = work[*n + i__], abs(d__1)) + safe1) 
			/ (work[i__] + safe1);
#line 367 "dgbrfs.f"
		s = max(d__2,d__3);
#line 369 "dgbrfs.f"
	    }
#line 370 "dgbrfs.f"
/* L80: */
#line 370 "dgbrfs.f"
	}
#line 371 "dgbrfs.f"
	berr[j] = s;

/*        Test stopping criterion. Continue iterating if */
/*           1) The residual BERR(J) is larger than machine epsilon, and */
/*           2) BERR(J) decreased by at least a factor of 2 during the */
/*              last iteration, and */
/*           3) At most ITMAX iterations tried. */

#line 379 "dgbrfs.f"
	if (berr[j] > eps && berr[j] * 2. <= lstres && count <= 5) {

/*           Update solution and try again. */

#line 384 "dgbrfs.f"
	    dgbtrs_(trans, n, kl, ku, &c__1, &afb[afb_offset], ldafb, &ipiv[1]
		    , &work[*n + 1], n, info, (ftnlen)1);
#line 386 "dgbrfs.f"
	    daxpy_(n, &c_b17, &work[*n + 1], &c__1, &x[j * x_dim1 + 1], &c__1)
		    ;
#line 387 "dgbrfs.f"
	    lstres = berr[j];
#line 388 "dgbrfs.f"
	    ++count;
#line 389 "dgbrfs.f"
	    goto L20;
#line 390 "dgbrfs.f"
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

#line 414 "dgbrfs.f"
	i__2 = *n;
#line 414 "dgbrfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 415 "dgbrfs.f"
	    if (work[i__] > safe2) {
#line 416 "dgbrfs.f"
		work[i__] = (d__1 = work[*n + i__], abs(d__1)) + nz * eps * 
			work[i__];
#line 417 "dgbrfs.f"
	    } else {
#line 418 "dgbrfs.f"
		work[i__] = (d__1 = work[*n + i__], abs(d__1)) + nz * eps * 
			work[i__] + safe1;
#line 419 "dgbrfs.f"
	    }
#line 420 "dgbrfs.f"
/* L90: */
#line 420 "dgbrfs.f"
	}

#line 422 "dgbrfs.f"
	kase = 0;
#line 423 "dgbrfs.f"
L100:
#line 424 "dgbrfs.f"
	dlacn2_(n, &work[(*n << 1) + 1], &work[*n + 1], &iwork[1], &ferr[j], &
		kase, isave);
#line 426 "dgbrfs.f"
	if (kase != 0) {
#line 427 "dgbrfs.f"
	    if (kase == 1) {

/*              Multiply by diag(W)*inv(op(A)**T). */

#line 431 "dgbrfs.f"
		dgbtrs_(transt, n, kl, ku, &c__1, &afb[afb_offset], ldafb, &
			ipiv[1], &work[*n + 1], n, info, (ftnlen)1);
#line 433 "dgbrfs.f"
		i__2 = *n;
#line 433 "dgbrfs.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 434 "dgbrfs.f"
		    work[*n + i__] *= work[i__];
#line 435 "dgbrfs.f"
/* L110: */
#line 435 "dgbrfs.f"
		}
#line 436 "dgbrfs.f"
	    } else {

/*              Multiply by inv(op(A))*diag(W). */

#line 440 "dgbrfs.f"
		i__2 = *n;
#line 440 "dgbrfs.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 441 "dgbrfs.f"
		    work[*n + i__] *= work[i__];
#line 442 "dgbrfs.f"
/* L120: */
#line 442 "dgbrfs.f"
		}
#line 443 "dgbrfs.f"
		dgbtrs_(trans, n, kl, ku, &c__1, &afb[afb_offset], ldafb, &
			ipiv[1], &work[*n + 1], n, info, (ftnlen)1);
#line 445 "dgbrfs.f"
	    }
#line 446 "dgbrfs.f"
	    goto L100;
#line 447 "dgbrfs.f"
	}

/*        Normalize error. */

#line 451 "dgbrfs.f"
	lstres = 0.;
#line 452 "dgbrfs.f"
	i__2 = *n;
#line 452 "dgbrfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 453 "dgbrfs.f"
	    d__2 = lstres, d__3 = (d__1 = x[i__ + j * x_dim1], abs(d__1));
#line 453 "dgbrfs.f"
	    lstres = max(d__2,d__3);
#line 454 "dgbrfs.f"
/* L130: */
#line 454 "dgbrfs.f"
	}
#line 455 "dgbrfs.f"
	if (lstres != 0.) {
#line 455 "dgbrfs.f"
	    ferr[j] /= lstres;
#line 455 "dgbrfs.f"
	}

#line 458 "dgbrfs.f"
/* L140: */
#line 458 "dgbrfs.f"
    }

#line 460 "dgbrfs.f"
    return 0;

/*     End of DGBRFS */

} /* dgbrfs_ */

