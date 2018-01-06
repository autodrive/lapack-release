#line 1 "dgerfs.f"
/* dgerfs.f -- translated by f2c (version 20100827).
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

#line 1 "dgerfs.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b15 = -1.;
static doublereal c_b17 = 1.;

/* > \brief \b DGERFS */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DGERFS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgerfs.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgerfs.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgerfs.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DGERFS( TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV, B, LDB, */
/*                          X, LDX, FERR, BERR, WORK, IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          TRANS */
/*       INTEGER            INFO, LDA, LDAF, LDB, LDX, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ), IWORK( * ) */
/*       DOUBLE PRECISION   A( LDA, * ), AF( LDAF, * ), B( LDB, * ), */
/*      $                   BERR( * ), FERR( * ), WORK( * ), X( LDX, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DGERFS improves the computed solution to a system of linear */
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
/* >          of the matrices B and X.  NRHS >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
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
/* >          AF is DOUBLE PRECISION array, dimension (LDAF,N) */
/* >          The factors L and U from the factorization A = P*L*U */
/* >          as computed by DGETRF. */
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
/* >          The pivot indices from DGETRF; for 1<=i<=N, row i of the */
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
/* >          On entry, the solution matrix X, as computed by DGETRS. */
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

/* > \ingroup doubleGEcomputational */

/*  ===================================================================== */
/* Subroutine */ int dgerfs_(char *trans, integer *n, integer *nrhs, 
	doublereal *a, integer *lda, doublereal *af, integer *ldaf, integer *
	ipiv, doublereal *b, integer *ldb, doublereal *x, integer *ldx, 
	doublereal *ferr, doublereal *berr, doublereal *work, integer *iwork, 
	integer *info, ftnlen trans_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, af_dim1, af_offset, b_dim1, b_offset, x_dim1, 
	    x_offset, i__1, i__2, i__3;
    doublereal d__1, d__2, d__3;

    /* Local variables */
    static integer i__, j, k;
    static doublereal s, xk;
    static integer nz;
    static doublereal eps;
    static integer kase;
    static doublereal safe1, safe2;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen);
    static integer isave[3];
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), daxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *);
    static integer count;
    extern /* Subroutine */ int dlacn2_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, integer *);
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), dgetrs_(
	    char *, integer *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, integer *, integer *, ftnlen);
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

#line 241 "dgerfs.f"
    /* Parameter adjustments */
#line 241 "dgerfs.f"
    a_dim1 = *lda;
#line 241 "dgerfs.f"
    a_offset = 1 + a_dim1;
#line 241 "dgerfs.f"
    a -= a_offset;
#line 241 "dgerfs.f"
    af_dim1 = *ldaf;
#line 241 "dgerfs.f"
    af_offset = 1 + af_dim1;
#line 241 "dgerfs.f"
    af -= af_offset;
#line 241 "dgerfs.f"
    --ipiv;
#line 241 "dgerfs.f"
    b_dim1 = *ldb;
#line 241 "dgerfs.f"
    b_offset = 1 + b_dim1;
#line 241 "dgerfs.f"
    b -= b_offset;
#line 241 "dgerfs.f"
    x_dim1 = *ldx;
#line 241 "dgerfs.f"
    x_offset = 1 + x_dim1;
#line 241 "dgerfs.f"
    x -= x_offset;
#line 241 "dgerfs.f"
    --ferr;
#line 241 "dgerfs.f"
    --berr;
#line 241 "dgerfs.f"
    --work;
#line 241 "dgerfs.f"
    --iwork;
#line 241 "dgerfs.f"

#line 241 "dgerfs.f"
    /* Function Body */
#line 241 "dgerfs.f"
    *info = 0;
#line 242 "dgerfs.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 243 "dgerfs.f"
    if (! notran && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1) && ! lsame_(
	    trans, "C", (ftnlen)1, (ftnlen)1)) {
#line 245 "dgerfs.f"
	*info = -1;
#line 246 "dgerfs.f"
    } else if (*n < 0) {
#line 247 "dgerfs.f"
	*info = -2;
#line 248 "dgerfs.f"
    } else if (*nrhs < 0) {
#line 249 "dgerfs.f"
	*info = -3;
#line 250 "dgerfs.f"
    } else if (*lda < max(1,*n)) {
#line 251 "dgerfs.f"
	*info = -5;
#line 252 "dgerfs.f"
    } else if (*ldaf < max(1,*n)) {
#line 253 "dgerfs.f"
	*info = -7;
#line 254 "dgerfs.f"
    } else if (*ldb < max(1,*n)) {
#line 255 "dgerfs.f"
	*info = -10;
#line 256 "dgerfs.f"
    } else if (*ldx < max(1,*n)) {
#line 257 "dgerfs.f"
	*info = -12;
#line 258 "dgerfs.f"
    }
#line 259 "dgerfs.f"
    if (*info != 0) {
#line 260 "dgerfs.f"
	i__1 = -(*info);
#line 260 "dgerfs.f"
	xerbla_("DGERFS", &i__1, (ftnlen)6);
#line 261 "dgerfs.f"
	return 0;
#line 262 "dgerfs.f"
    }

/*     Quick return if possible */

#line 266 "dgerfs.f"
    if (*n == 0 || *nrhs == 0) {
#line 267 "dgerfs.f"
	i__1 = *nrhs;
#line 267 "dgerfs.f"
	for (j = 1; j <= i__1; ++j) {
#line 268 "dgerfs.f"
	    ferr[j] = 0.;
#line 269 "dgerfs.f"
	    berr[j] = 0.;
#line 270 "dgerfs.f"
/* L10: */
#line 270 "dgerfs.f"
	}
#line 271 "dgerfs.f"
	return 0;
#line 272 "dgerfs.f"
    }

#line 274 "dgerfs.f"
    if (notran) {
#line 275 "dgerfs.f"
	*(unsigned char *)transt = 'T';
#line 276 "dgerfs.f"
    } else {
#line 277 "dgerfs.f"
	*(unsigned char *)transt = 'N';
#line 278 "dgerfs.f"
    }

/*     NZ = maximum number of nonzero elements in each row of A, plus 1 */

#line 282 "dgerfs.f"
    nz = *n + 1;
#line 283 "dgerfs.f"
    eps = dlamch_("Epsilon", (ftnlen)7);
#line 284 "dgerfs.f"
    safmin = dlamch_("Safe minimum", (ftnlen)12);
#line 285 "dgerfs.f"
    safe1 = nz * safmin;
#line 286 "dgerfs.f"
    safe2 = safe1 / eps;

/*     Do for each right hand side */

#line 290 "dgerfs.f"
    i__1 = *nrhs;
#line 290 "dgerfs.f"
    for (j = 1; j <= i__1; ++j) {

#line 292 "dgerfs.f"
	count = 1;
#line 293 "dgerfs.f"
	lstres = 3.;
#line 294 "dgerfs.f"
L20:

/*        Loop until stopping criterion is satisfied. */

/*        Compute residual R = B - op(A) * X, */
/*        where op(A) = A, A**T, or A**H, depending on TRANS. */

#line 301 "dgerfs.f"
	dcopy_(n, &b[j * b_dim1 + 1], &c__1, &work[*n + 1], &c__1);
#line 302 "dgerfs.f"
	dgemv_(trans, n, n, &c_b15, &a[a_offset], lda, &x[j * x_dim1 + 1], &
		c__1, &c_b17, &work[*n + 1], &c__1, (ftnlen)1);

/*        Compute componentwise relative backward error from formula */

/*        max(i) ( abs(R(i)) / ( abs(op(A))*abs(X) + abs(B) )(i) ) */

/*        where abs(Z) is the componentwise absolute value of the matrix */
/*        or vector Z.  If the i-th component of the denominator is less */
/*        than SAFE2, then SAFE1 is added to the i-th components of the */
/*        numerator and denominator before dividing. */

#line 314 "dgerfs.f"
	i__2 = *n;
#line 314 "dgerfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 315 "dgerfs.f"
	    work[i__] = (d__1 = b[i__ + j * b_dim1], abs(d__1));
#line 316 "dgerfs.f"
/* L30: */
#line 316 "dgerfs.f"
	}

/*        Compute abs(op(A))*abs(X) + abs(B). */

#line 320 "dgerfs.f"
	if (notran) {
#line 321 "dgerfs.f"
	    i__2 = *n;
#line 321 "dgerfs.f"
	    for (k = 1; k <= i__2; ++k) {
#line 322 "dgerfs.f"
		xk = (d__1 = x[k + j * x_dim1], abs(d__1));
#line 323 "dgerfs.f"
		i__3 = *n;
#line 323 "dgerfs.f"
		for (i__ = 1; i__ <= i__3; ++i__) {
#line 324 "dgerfs.f"
		    work[i__] += (d__1 = a[i__ + k * a_dim1], abs(d__1)) * xk;
#line 325 "dgerfs.f"
/* L40: */
#line 325 "dgerfs.f"
		}
#line 326 "dgerfs.f"
/* L50: */
#line 326 "dgerfs.f"
	    }
#line 327 "dgerfs.f"
	} else {
#line 328 "dgerfs.f"
	    i__2 = *n;
#line 328 "dgerfs.f"
	    for (k = 1; k <= i__2; ++k) {
#line 329 "dgerfs.f"
		s = 0.;
#line 330 "dgerfs.f"
		i__3 = *n;
#line 330 "dgerfs.f"
		for (i__ = 1; i__ <= i__3; ++i__) {
#line 331 "dgerfs.f"
		    s += (d__1 = a[i__ + k * a_dim1], abs(d__1)) * (d__2 = x[
			    i__ + j * x_dim1], abs(d__2));
#line 332 "dgerfs.f"
/* L60: */
#line 332 "dgerfs.f"
		}
#line 333 "dgerfs.f"
		work[k] += s;
#line 334 "dgerfs.f"
/* L70: */
#line 334 "dgerfs.f"
	    }
#line 335 "dgerfs.f"
	}
#line 336 "dgerfs.f"
	s = 0.;
#line 337 "dgerfs.f"
	i__2 = *n;
#line 337 "dgerfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 338 "dgerfs.f"
	    if (work[i__] > safe2) {
/* Computing MAX */
#line 339 "dgerfs.f"
		d__2 = s, d__3 = (d__1 = work[*n + i__], abs(d__1)) / work[
			i__];
#line 339 "dgerfs.f"
		s = max(d__2,d__3);
#line 340 "dgerfs.f"
	    } else {
/* Computing MAX */
#line 341 "dgerfs.f"
		d__2 = s, d__3 = ((d__1 = work[*n + i__], abs(d__1)) + safe1) 
			/ (work[i__] + safe1);
#line 341 "dgerfs.f"
		s = max(d__2,d__3);
#line 343 "dgerfs.f"
	    }
#line 344 "dgerfs.f"
/* L80: */
#line 344 "dgerfs.f"
	}
#line 345 "dgerfs.f"
	berr[j] = s;

/*        Test stopping criterion. Continue iterating if */
/*           1) The residual BERR(J) is larger than machine epsilon, and */
/*           2) BERR(J) decreased by at least a factor of 2 during the */
/*              last iteration, and */
/*           3) At most ITMAX iterations tried. */

#line 353 "dgerfs.f"
	if (berr[j] > eps && berr[j] * 2. <= lstres && count <= 5) {

/*           Update solution and try again. */

#line 358 "dgerfs.f"
	    dgetrs_(trans, n, &c__1, &af[af_offset], ldaf, &ipiv[1], &work[*n 
		    + 1], n, info, (ftnlen)1);
#line 360 "dgerfs.f"
	    daxpy_(n, &c_b17, &work[*n + 1], &c__1, &x[j * x_dim1 + 1], &c__1)
		    ;
#line 361 "dgerfs.f"
	    lstres = berr[j];
#line 362 "dgerfs.f"
	    ++count;
#line 363 "dgerfs.f"
	    goto L20;
#line 364 "dgerfs.f"
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

#line 388 "dgerfs.f"
	i__2 = *n;
#line 388 "dgerfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 389 "dgerfs.f"
	    if (work[i__] > safe2) {
#line 390 "dgerfs.f"
		work[i__] = (d__1 = work[*n + i__], abs(d__1)) + nz * eps * 
			work[i__];
#line 391 "dgerfs.f"
	    } else {
#line 392 "dgerfs.f"
		work[i__] = (d__1 = work[*n + i__], abs(d__1)) + nz * eps * 
			work[i__] + safe1;
#line 393 "dgerfs.f"
	    }
#line 394 "dgerfs.f"
/* L90: */
#line 394 "dgerfs.f"
	}

#line 396 "dgerfs.f"
	kase = 0;
#line 397 "dgerfs.f"
L100:
#line 398 "dgerfs.f"
	dlacn2_(n, &work[(*n << 1) + 1], &work[*n + 1], &iwork[1], &ferr[j], &
		kase, isave);
#line 400 "dgerfs.f"
	if (kase != 0) {
#line 401 "dgerfs.f"
	    if (kase == 1) {

/*              Multiply by diag(W)*inv(op(A)**T). */

#line 405 "dgerfs.f"
		dgetrs_(transt, n, &c__1, &af[af_offset], ldaf, &ipiv[1], &
			work[*n + 1], n, info, (ftnlen)1);
#line 407 "dgerfs.f"
		i__2 = *n;
#line 407 "dgerfs.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 408 "dgerfs.f"
		    work[*n + i__] = work[i__] * work[*n + i__];
#line 409 "dgerfs.f"
/* L110: */
#line 409 "dgerfs.f"
		}
#line 410 "dgerfs.f"
	    } else {

/*              Multiply by inv(op(A))*diag(W). */

#line 414 "dgerfs.f"
		i__2 = *n;
#line 414 "dgerfs.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 415 "dgerfs.f"
		    work[*n + i__] = work[i__] * work[*n + i__];
#line 416 "dgerfs.f"
/* L120: */
#line 416 "dgerfs.f"
		}
#line 417 "dgerfs.f"
		dgetrs_(trans, n, &c__1, &af[af_offset], ldaf, &ipiv[1], &
			work[*n + 1], n, info, (ftnlen)1);
#line 419 "dgerfs.f"
	    }
#line 420 "dgerfs.f"
	    goto L100;
#line 421 "dgerfs.f"
	}

/*        Normalize error. */

#line 425 "dgerfs.f"
	lstres = 0.;
#line 426 "dgerfs.f"
	i__2 = *n;
#line 426 "dgerfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 427 "dgerfs.f"
	    d__2 = lstres, d__3 = (d__1 = x[i__ + j * x_dim1], abs(d__1));
#line 427 "dgerfs.f"
	    lstres = max(d__2,d__3);
#line 428 "dgerfs.f"
/* L130: */
#line 428 "dgerfs.f"
	}
#line 429 "dgerfs.f"
	if (lstres != 0.) {
#line 429 "dgerfs.f"
	    ferr[j] /= lstres;
#line 429 "dgerfs.f"
	}

#line 432 "dgerfs.f"
/* L140: */
#line 432 "dgerfs.f"
    }

#line 434 "dgerfs.f"
    return 0;

/*     End of DGERFS */

} /* dgerfs_ */

