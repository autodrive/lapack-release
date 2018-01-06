#line 1 "sgerfs.f"
/* sgerfs.f -- translated by f2c (version 20100827).
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

#line 1 "sgerfs.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b15 = -1.;
static doublereal c_b17 = 1.;

/* > \brief \b SGERFS */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SGERFS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgerfs.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgerfs.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgerfs.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SGERFS( TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV, B, LDB, */
/*                          X, LDX, FERR, BERR, WORK, IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          TRANS */
/*       INTEGER            INFO, LDA, LDAF, LDB, LDX, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ), IWORK( * ) */
/*       REAL               A( LDA, * ), AF( LDAF, * ), B( LDB, * ), */
/*      $                   BERR( * ), FERR( * ), WORK( * ), X( LDX, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SGERFS improves the computed solution to a system of linear */
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
/* >          A is REAL array, dimension (LDA,N) */
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
/* >          AF is REAL array, dimension (LDAF,N) */
/* >          The factors L and U from the factorization A = P*L*U */
/* >          as computed by SGETRF. */
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
/* >          The pivot indices from SGETRF; for 1<=i<=N, row i of the */
/* >          matrix was interchanged with row IPIV(i). */
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
/* >          On entry, the solution matrix X, as computed by SGETRS. */
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

/* > \ingroup realGEcomputational */

/*  ===================================================================== */
/* Subroutine */ int sgerfs_(char *trans, integer *n, integer *nrhs, 
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
    static integer isave[3];
    extern /* Subroutine */ int sgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen);
    static integer count;
    extern /* Subroutine */ int scopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), saxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *), slacn2_(
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, integer *);
    extern doublereal slamch_(char *, ftnlen);
    static doublereal safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static logical notran;
    extern /* Subroutine */ int sgetrs_(char *, integer *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, ftnlen);
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

#line 241 "sgerfs.f"
    /* Parameter adjustments */
#line 241 "sgerfs.f"
    a_dim1 = *lda;
#line 241 "sgerfs.f"
    a_offset = 1 + a_dim1;
#line 241 "sgerfs.f"
    a -= a_offset;
#line 241 "sgerfs.f"
    af_dim1 = *ldaf;
#line 241 "sgerfs.f"
    af_offset = 1 + af_dim1;
#line 241 "sgerfs.f"
    af -= af_offset;
#line 241 "sgerfs.f"
    --ipiv;
#line 241 "sgerfs.f"
    b_dim1 = *ldb;
#line 241 "sgerfs.f"
    b_offset = 1 + b_dim1;
#line 241 "sgerfs.f"
    b -= b_offset;
#line 241 "sgerfs.f"
    x_dim1 = *ldx;
#line 241 "sgerfs.f"
    x_offset = 1 + x_dim1;
#line 241 "sgerfs.f"
    x -= x_offset;
#line 241 "sgerfs.f"
    --ferr;
#line 241 "sgerfs.f"
    --berr;
#line 241 "sgerfs.f"
    --work;
#line 241 "sgerfs.f"
    --iwork;
#line 241 "sgerfs.f"

#line 241 "sgerfs.f"
    /* Function Body */
#line 241 "sgerfs.f"
    *info = 0;
#line 242 "sgerfs.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 243 "sgerfs.f"
    if (! notran && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1) && ! lsame_(
	    trans, "C", (ftnlen)1, (ftnlen)1)) {
#line 245 "sgerfs.f"
	*info = -1;
#line 246 "sgerfs.f"
    } else if (*n < 0) {
#line 247 "sgerfs.f"
	*info = -2;
#line 248 "sgerfs.f"
    } else if (*nrhs < 0) {
#line 249 "sgerfs.f"
	*info = -3;
#line 250 "sgerfs.f"
    } else if (*lda < max(1,*n)) {
#line 251 "sgerfs.f"
	*info = -5;
#line 252 "sgerfs.f"
    } else if (*ldaf < max(1,*n)) {
#line 253 "sgerfs.f"
	*info = -7;
#line 254 "sgerfs.f"
    } else if (*ldb < max(1,*n)) {
#line 255 "sgerfs.f"
	*info = -10;
#line 256 "sgerfs.f"
    } else if (*ldx < max(1,*n)) {
#line 257 "sgerfs.f"
	*info = -12;
#line 258 "sgerfs.f"
    }
#line 259 "sgerfs.f"
    if (*info != 0) {
#line 260 "sgerfs.f"
	i__1 = -(*info);
#line 260 "sgerfs.f"
	xerbla_("SGERFS", &i__1, (ftnlen)6);
#line 261 "sgerfs.f"
	return 0;
#line 262 "sgerfs.f"
    }

/*     Quick return if possible */

#line 266 "sgerfs.f"
    if (*n == 0 || *nrhs == 0) {
#line 267 "sgerfs.f"
	i__1 = *nrhs;
#line 267 "sgerfs.f"
	for (j = 1; j <= i__1; ++j) {
#line 268 "sgerfs.f"
	    ferr[j] = 0.;
#line 269 "sgerfs.f"
	    berr[j] = 0.;
#line 270 "sgerfs.f"
/* L10: */
#line 270 "sgerfs.f"
	}
#line 271 "sgerfs.f"
	return 0;
#line 272 "sgerfs.f"
    }

#line 274 "sgerfs.f"
    if (notran) {
#line 275 "sgerfs.f"
	*(unsigned char *)transt = 'T';
#line 276 "sgerfs.f"
    } else {
#line 277 "sgerfs.f"
	*(unsigned char *)transt = 'N';
#line 278 "sgerfs.f"
    }

/*     NZ = maximum number of nonzero elements in each row of A, plus 1 */

#line 282 "sgerfs.f"
    nz = *n + 1;
#line 283 "sgerfs.f"
    eps = slamch_("Epsilon", (ftnlen)7);
#line 284 "sgerfs.f"
    safmin = slamch_("Safe minimum", (ftnlen)12);
#line 285 "sgerfs.f"
    safe1 = nz * safmin;
#line 286 "sgerfs.f"
    safe2 = safe1 / eps;

/*     Do for each right hand side */

#line 290 "sgerfs.f"
    i__1 = *nrhs;
#line 290 "sgerfs.f"
    for (j = 1; j <= i__1; ++j) {

#line 292 "sgerfs.f"
	count = 1;
#line 293 "sgerfs.f"
	lstres = 3.;
#line 294 "sgerfs.f"
L20:

/*        Loop until stopping criterion is satisfied. */

/*        Compute residual R = B - op(A) * X, */
/*        where op(A) = A, A**T, or A**H, depending on TRANS. */

#line 301 "sgerfs.f"
	scopy_(n, &b[j * b_dim1 + 1], &c__1, &work[*n + 1], &c__1);
#line 302 "sgerfs.f"
	sgemv_(trans, n, n, &c_b15, &a[a_offset], lda, &x[j * x_dim1 + 1], &
		c__1, &c_b17, &work[*n + 1], &c__1, (ftnlen)1);

/*        Compute componentwise relative backward error from formula */

/*        max(i) ( abs(R(i)) / ( abs(op(A))*abs(X) + abs(B) )(i) ) */

/*        where abs(Z) is the componentwise absolute value of the matrix */
/*        or vector Z.  If the i-th component of the denominator is less */
/*        than SAFE2, then SAFE1 is added to the i-th components of the */
/*        numerator and denominator before dividing. */

#line 314 "sgerfs.f"
	i__2 = *n;
#line 314 "sgerfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 315 "sgerfs.f"
	    work[i__] = (d__1 = b[i__ + j * b_dim1], abs(d__1));
#line 316 "sgerfs.f"
/* L30: */
#line 316 "sgerfs.f"
	}

/*        Compute abs(op(A))*abs(X) + abs(B). */

#line 320 "sgerfs.f"
	if (notran) {
#line 321 "sgerfs.f"
	    i__2 = *n;
#line 321 "sgerfs.f"
	    for (k = 1; k <= i__2; ++k) {
#line 322 "sgerfs.f"
		xk = (d__1 = x[k + j * x_dim1], abs(d__1));
#line 323 "sgerfs.f"
		i__3 = *n;
#line 323 "sgerfs.f"
		for (i__ = 1; i__ <= i__3; ++i__) {
#line 324 "sgerfs.f"
		    work[i__] += (d__1 = a[i__ + k * a_dim1], abs(d__1)) * xk;
#line 325 "sgerfs.f"
/* L40: */
#line 325 "sgerfs.f"
		}
#line 326 "sgerfs.f"
/* L50: */
#line 326 "sgerfs.f"
	    }
#line 327 "sgerfs.f"
	} else {
#line 328 "sgerfs.f"
	    i__2 = *n;
#line 328 "sgerfs.f"
	    for (k = 1; k <= i__2; ++k) {
#line 329 "sgerfs.f"
		s = 0.;
#line 330 "sgerfs.f"
		i__3 = *n;
#line 330 "sgerfs.f"
		for (i__ = 1; i__ <= i__3; ++i__) {
#line 331 "sgerfs.f"
		    s += (d__1 = a[i__ + k * a_dim1], abs(d__1)) * (d__2 = x[
			    i__ + j * x_dim1], abs(d__2));
#line 332 "sgerfs.f"
/* L60: */
#line 332 "sgerfs.f"
		}
#line 333 "sgerfs.f"
		work[k] += s;
#line 334 "sgerfs.f"
/* L70: */
#line 334 "sgerfs.f"
	    }
#line 335 "sgerfs.f"
	}
#line 336 "sgerfs.f"
	s = 0.;
#line 337 "sgerfs.f"
	i__2 = *n;
#line 337 "sgerfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 338 "sgerfs.f"
	    if (work[i__] > safe2) {
/* Computing MAX */
#line 339 "sgerfs.f"
		d__2 = s, d__3 = (d__1 = work[*n + i__], abs(d__1)) / work[
			i__];
#line 339 "sgerfs.f"
		s = max(d__2,d__3);
#line 340 "sgerfs.f"
	    } else {
/* Computing MAX */
#line 341 "sgerfs.f"
		d__2 = s, d__3 = ((d__1 = work[*n + i__], abs(d__1)) + safe1) 
			/ (work[i__] + safe1);
#line 341 "sgerfs.f"
		s = max(d__2,d__3);
#line 343 "sgerfs.f"
	    }
#line 344 "sgerfs.f"
/* L80: */
#line 344 "sgerfs.f"
	}
#line 345 "sgerfs.f"
	berr[j] = s;

/*        Test stopping criterion. Continue iterating if */
/*           1) The residual BERR(J) is larger than machine epsilon, and */
/*           2) BERR(J) decreased by at least a factor of 2 during the */
/*              last iteration, and */
/*           3) At most ITMAX iterations tried. */

#line 353 "sgerfs.f"
	if (berr[j] > eps && berr[j] * 2. <= lstres && count <= 5) {

/*           Update solution and try again. */

#line 358 "sgerfs.f"
	    sgetrs_(trans, n, &c__1, &af[af_offset], ldaf, &ipiv[1], &work[*n 
		    + 1], n, info, (ftnlen)1);
#line 360 "sgerfs.f"
	    saxpy_(n, &c_b17, &work[*n + 1], &c__1, &x[j * x_dim1 + 1], &c__1)
		    ;
#line 361 "sgerfs.f"
	    lstres = berr[j];
#line 362 "sgerfs.f"
	    ++count;
#line 363 "sgerfs.f"
	    goto L20;
#line 364 "sgerfs.f"
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

#line 388 "sgerfs.f"
	i__2 = *n;
#line 388 "sgerfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 389 "sgerfs.f"
	    if (work[i__] > safe2) {
#line 390 "sgerfs.f"
		work[i__] = (d__1 = work[*n + i__], abs(d__1)) + nz * eps * 
			work[i__];
#line 391 "sgerfs.f"
	    } else {
#line 392 "sgerfs.f"
		work[i__] = (d__1 = work[*n + i__], abs(d__1)) + nz * eps * 
			work[i__] + safe1;
#line 393 "sgerfs.f"
	    }
#line 394 "sgerfs.f"
/* L90: */
#line 394 "sgerfs.f"
	}

#line 396 "sgerfs.f"
	kase = 0;
#line 397 "sgerfs.f"
L100:
#line 398 "sgerfs.f"
	slacn2_(n, &work[(*n << 1) + 1], &work[*n + 1], &iwork[1], &ferr[j], &
		kase, isave);
#line 400 "sgerfs.f"
	if (kase != 0) {
#line 401 "sgerfs.f"
	    if (kase == 1) {

/*              Multiply by diag(W)*inv(op(A)**T). */

#line 405 "sgerfs.f"
		sgetrs_(transt, n, &c__1, &af[af_offset], ldaf, &ipiv[1], &
			work[*n + 1], n, info, (ftnlen)1);
#line 407 "sgerfs.f"
		i__2 = *n;
#line 407 "sgerfs.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 408 "sgerfs.f"
		    work[*n + i__] = work[i__] * work[*n + i__];
#line 409 "sgerfs.f"
/* L110: */
#line 409 "sgerfs.f"
		}
#line 410 "sgerfs.f"
	    } else {

/*              Multiply by inv(op(A))*diag(W). */

#line 414 "sgerfs.f"
		i__2 = *n;
#line 414 "sgerfs.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 415 "sgerfs.f"
		    work[*n + i__] = work[i__] * work[*n + i__];
#line 416 "sgerfs.f"
/* L120: */
#line 416 "sgerfs.f"
		}
#line 417 "sgerfs.f"
		sgetrs_(trans, n, &c__1, &af[af_offset], ldaf, &ipiv[1], &
			work[*n + 1], n, info, (ftnlen)1);
#line 419 "sgerfs.f"
	    }
#line 420 "sgerfs.f"
	    goto L100;
#line 421 "sgerfs.f"
	}

/*        Normalize error. */

#line 425 "sgerfs.f"
	lstres = 0.;
#line 426 "sgerfs.f"
	i__2 = *n;
#line 426 "sgerfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 427 "sgerfs.f"
	    d__2 = lstres, d__3 = (d__1 = x[i__ + j * x_dim1], abs(d__1));
#line 427 "sgerfs.f"
	    lstres = max(d__2,d__3);
#line 428 "sgerfs.f"
/* L130: */
#line 428 "sgerfs.f"
	}
#line 429 "sgerfs.f"
	if (lstres != 0.) {
#line 429 "sgerfs.f"
	    ferr[j] /= lstres;
#line 429 "sgerfs.f"
	}

#line 432 "sgerfs.f"
/* L140: */
#line 432 "sgerfs.f"
    }

#line 434 "sgerfs.f"
    return 0;

/*     End of SGERFS */

} /* sgerfs_ */

