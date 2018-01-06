#line 1 "cgbrfs.f"
/* cgbrfs.f -- translated by f2c (version 20100827).
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

#line 1 "cgbrfs.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static integer c__1 = 1;

/* > \brief \b CGBRFS */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CGBRFS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgbrfs.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgbrfs.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgbrfs.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CGBRFS( TRANS, N, KL, KU, NRHS, AB, LDAB, AFB, LDAFB, */
/*                          IPIV, B, LDB, X, LDX, FERR, BERR, WORK, RWORK, */
/*                          INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          TRANS */
/*       INTEGER            INFO, KL, KU, LDAB, LDAFB, LDB, LDX, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       REAL               BERR( * ), FERR( * ), RWORK( * ) */
/*       COMPLEX            AB( LDAB, * ), AFB( LDAFB, * ), B( LDB, * ), */
/*      $                   WORK( * ), X( LDX, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CGBRFS improves the computed solution to a system of linear */
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
/* >          = 'C':  A**H * X = B  (Conjugate transpose) */
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
/* >          AB is COMPLEX array, dimension (LDAB,N) */
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
/* >          AFB is COMPLEX array, dimension (LDAFB,N) */
/* >          Details of the LU factorization of the band matrix A, as */
/* >          computed by CGBTRF.  U is stored as an upper triangular band */
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
/* >          The pivot indices from CGBTRF; for 1<=i<=N, row i of the */
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
/* >          On entry, the solution matrix X, as computed by CGBTRS. */
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

/* > \ingroup complexGBcomputational */

/*  ===================================================================== */
/* Subroutine */ int cgbrfs_(char *trans, integer *n, integer *kl, integer *
	ku, integer *nrhs, doublecomplex *ab, integer *ldab, doublecomplex *
	afb, integer *ldafb, integer *ipiv, doublecomplex *b, integer *ldb, 
	doublecomplex *x, integer *ldx, doublereal *ferr, doublereal *berr, 
	doublecomplex *work, doublereal *rwork, integer *info, ftnlen 
	trans_len)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, afb_dim1, afb_offset, b_dim1, b_offset, 
	    x_dim1, x_offset, i__1, i__2, i__3, i__4, i__5, i__6, i__7;
    doublereal d__1, d__2, d__3, d__4;
    doublecomplex z__1;

    /* Builtin functions */
    double d_imag(doublecomplex *);

    /* Local variables */
    static integer i__, j, k;
    static doublereal s;
    static integer kk;
    static doublereal xk;
    static integer nz;
    static doublereal eps;
    static integer kase;
    static doublereal safe1, safe2;
    extern /* Subroutine */ int cgbmv_(char *, integer *, integer *, integer *
	    , integer *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer isave[3];
    extern /* Subroutine */ int ccopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), caxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    static integer count;
    extern /* Subroutine */ int clacn2_(integer *, doublecomplex *, 
	    doublecomplex *, doublereal *, integer *, integer *);
    extern doublereal slamch_(char *, ftnlen);
    static doublereal safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), cgbtrs_(
	    char *, integer *, integer *, integer *, integer *, doublecomplex 
	    *, integer *, integer *, doublecomplex *, integer *, integer *, 
	    ftnlen);
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

#line 270 "cgbrfs.f"
    /* Parameter adjustments */
#line 270 "cgbrfs.f"
    ab_dim1 = *ldab;
#line 270 "cgbrfs.f"
    ab_offset = 1 + ab_dim1;
#line 270 "cgbrfs.f"
    ab -= ab_offset;
#line 270 "cgbrfs.f"
    afb_dim1 = *ldafb;
#line 270 "cgbrfs.f"
    afb_offset = 1 + afb_dim1;
#line 270 "cgbrfs.f"
    afb -= afb_offset;
#line 270 "cgbrfs.f"
    --ipiv;
#line 270 "cgbrfs.f"
    b_dim1 = *ldb;
#line 270 "cgbrfs.f"
    b_offset = 1 + b_dim1;
#line 270 "cgbrfs.f"
    b -= b_offset;
#line 270 "cgbrfs.f"
    x_dim1 = *ldx;
#line 270 "cgbrfs.f"
    x_offset = 1 + x_dim1;
#line 270 "cgbrfs.f"
    x -= x_offset;
#line 270 "cgbrfs.f"
    --ferr;
#line 270 "cgbrfs.f"
    --berr;
#line 270 "cgbrfs.f"
    --work;
#line 270 "cgbrfs.f"
    --rwork;
#line 270 "cgbrfs.f"

#line 270 "cgbrfs.f"
    /* Function Body */
#line 270 "cgbrfs.f"
    *info = 0;
#line 271 "cgbrfs.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 272 "cgbrfs.f"
    if (! notran && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1) && ! lsame_(
	    trans, "C", (ftnlen)1, (ftnlen)1)) {
#line 274 "cgbrfs.f"
	*info = -1;
#line 275 "cgbrfs.f"
    } else if (*n < 0) {
#line 276 "cgbrfs.f"
	*info = -2;
#line 277 "cgbrfs.f"
    } else if (*kl < 0) {
#line 278 "cgbrfs.f"
	*info = -3;
#line 279 "cgbrfs.f"
    } else if (*ku < 0) {
#line 280 "cgbrfs.f"
	*info = -4;
#line 281 "cgbrfs.f"
    } else if (*nrhs < 0) {
#line 282 "cgbrfs.f"
	*info = -5;
#line 283 "cgbrfs.f"
    } else if (*ldab < *kl + *ku + 1) {
#line 284 "cgbrfs.f"
	*info = -7;
#line 285 "cgbrfs.f"
    } else if (*ldafb < (*kl << 1) + *ku + 1) {
#line 286 "cgbrfs.f"
	*info = -9;
#line 287 "cgbrfs.f"
    } else if (*ldb < max(1,*n)) {
#line 288 "cgbrfs.f"
	*info = -12;
#line 289 "cgbrfs.f"
    } else if (*ldx < max(1,*n)) {
#line 290 "cgbrfs.f"
	*info = -14;
#line 291 "cgbrfs.f"
    }
#line 292 "cgbrfs.f"
    if (*info != 0) {
#line 293 "cgbrfs.f"
	i__1 = -(*info);
#line 293 "cgbrfs.f"
	xerbla_("CGBRFS", &i__1, (ftnlen)6);
#line 294 "cgbrfs.f"
	return 0;
#line 295 "cgbrfs.f"
    }

/*     Quick return if possible */

#line 299 "cgbrfs.f"
    if (*n == 0 || *nrhs == 0) {
#line 300 "cgbrfs.f"
	i__1 = *nrhs;
#line 300 "cgbrfs.f"
	for (j = 1; j <= i__1; ++j) {
#line 301 "cgbrfs.f"
	    ferr[j] = 0.;
#line 302 "cgbrfs.f"
	    berr[j] = 0.;
#line 303 "cgbrfs.f"
/* L10: */
#line 303 "cgbrfs.f"
	}
#line 304 "cgbrfs.f"
	return 0;
#line 305 "cgbrfs.f"
    }

#line 307 "cgbrfs.f"
    if (notran) {
#line 308 "cgbrfs.f"
	*(unsigned char *)transn = 'N';
#line 309 "cgbrfs.f"
	*(unsigned char *)transt = 'C';
#line 310 "cgbrfs.f"
    } else {
#line 311 "cgbrfs.f"
	*(unsigned char *)transn = 'C';
#line 312 "cgbrfs.f"
	*(unsigned char *)transt = 'N';
#line 313 "cgbrfs.f"
    }

/*     NZ = maximum number of nonzero elements in each row of A, plus 1 */

/* Computing MIN */
#line 317 "cgbrfs.f"
    i__1 = *kl + *ku + 2, i__2 = *n + 1;
#line 317 "cgbrfs.f"
    nz = min(i__1,i__2);
#line 318 "cgbrfs.f"
    eps = slamch_("Epsilon", (ftnlen)7);
#line 319 "cgbrfs.f"
    safmin = slamch_("Safe minimum", (ftnlen)12);
#line 320 "cgbrfs.f"
    safe1 = nz * safmin;
#line 321 "cgbrfs.f"
    safe2 = safe1 / eps;

/*     Do for each right hand side */

#line 325 "cgbrfs.f"
    i__1 = *nrhs;
#line 325 "cgbrfs.f"
    for (j = 1; j <= i__1; ++j) {

#line 327 "cgbrfs.f"
	count = 1;
#line 328 "cgbrfs.f"
	lstres = 3.;
#line 329 "cgbrfs.f"
L20:

/*        Loop until stopping criterion is satisfied. */

/*        Compute residual R = B - op(A) * X, */
/*        where op(A) = A, A**T, or A**H, depending on TRANS. */

#line 336 "cgbrfs.f"
	ccopy_(n, &b[j * b_dim1 + 1], &c__1, &work[1], &c__1);
#line 337 "cgbrfs.f"
	z__1.r = -1., z__1.i = -0.;
#line 337 "cgbrfs.f"
	cgbmv_(trans, n, n, kl, ku, &z__1, &ab[ab_offset], ldab, &x[j * 
		x_dim1 + 1], &c__1, &c_b1, &work[1], &c__1, (ftnlen)1);

/*        Compute componentwise relative backward error from formula */

/*        max(i) ( abs(R(i)) / ( abs(op(A))*abs(X) + abs(B) )(i) ) */

/*        where abs(Z) is the componentwise absolute value of the matrix */
/*        or vector Z.  If the i-th component of the denominator is less */
/*        than SAFE2, then SAFE1 is added to the i-th components of the */
/*        numerator and denominator before dividing. */

#line 349 "cgbrfs.f"
	i__2 = *n;
#line 349 "cgbrfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 350 "cgbrfs.f"
	    i__3 = i__ + j * b_dim1;
#line 350 "cgbrfs.f"
	    rwork[i__] = (d__1 = b[i__3].r, abs(d__1)) + (d__2 = d_imag(&b[
		    i__ + j * b_dim1]), abs(d__2));
#line 351 "cgbrfs.f"
/* L30: */
#line 351 "cgbrfs.f"
	}

/*        Compute abs(op(A))*abs(X) + abs(B). */

#line 355 "cgbrfs.f"
	if (notran) {
#line 356 "cgbrfs.f"
	    i__2 = *n;
#line 356 "cgbrfs.f"
	    for (k = 1; k <= i__2; ++k) {
#line 357 "cgbrfs.f"
		kk = *ku + 1 - k;
#line 358 "cgbrfs.f"
		i__3 = k + j * x_dim1;
#line 358 "cgbrfs.f"
		xk = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = d_imag(&x[k + j *
			 x_dim1]), abs(d__2));
/* Computing MAX */
#line 359 "cgbrfs.f"
		i__3 = 1, i__4 = k - *ku;
/* Computing MIN */
#line 359 "cgbrfs.f"
		i__6 = *n, i__7 = k + *kl;
#line 359 "cgbrfs.f"
		i__5 = min(i__6,i__7);
#line 359 "cgbrfs.f"
		for (i__ = max(i__3,i__4); i__ <= i__5; ++i__) {
#line 360 "cgbrfs.f"
		    i__3 = kk + i__ + k * ab_dim1;
#line 360 "cgbrfs.f"
		    rwork[i__] += ((d__1 = ab[i__3].r, abs(d__1)) + (d__2 = 
			    d_imag(&ab[kk + i__ + k * ab_dim1]), abs(d__2))) *
			     xk;
#line 361 "cgbrfs.f"
/* L40: */
#line 361 "cgbrfs.f"
		}
#line 362 "cgbrfs.f"
/* L50: */
#line 362 "cgbrfs.f"
	    }
#line 363 "cgbrfs.f"
	} else {
#line 364 "cgbrfs.f"
	    i__2 = *n;
#line 364 "cgbrfs.f"
	    for (k = 1; k <= i__2; ++k) {
#line 365 "cgbrfs.f"
		s = 0.;
#line 366 "cgbrfs.f"
		kk = *ku + 1 - k;
/* Computing MAX */
#line 367 "cgbrfs.f"
		i__5 = 1, i__3 = k - *ku;
/* Computing MIN */
#line 367 "cgbrfs.f"
		i__6 = *n, i__7 = k + *kl;
#line 367 "cgbrfs.f"
		i__4 = min(i__6,i__7);
#line 367 "cgbrfs.f"
		for (i__ = max(i__5,i__3); i__ <= i__4; ++i__) {
#line 368 "cgbrfs.f"
		    i__5 = kk + i__ + k * ab_dim1;
#line 368 "cgbrfs.f"
		    i__3 = i__ + j * x_dim1;
#line 368 "cgbrfs.f"
		    s += ((d__1 = ab[i__5].r, abs(d__1)) + (d__2 = d_imag(&ab[
			    kk + i__ + k * ab_dim1]), abs(d__2))) * ((d__3 = 
			    x[i__3].r, abs(d__3)) + (d__4 = d_imag(&x[i__ + j 
			    * x_dim1]), abs(d__4)));
#line 369 "cgbrfs.f"
/* L60: */
#line 369 "cgbrfs.f"
		}
#line 370 "cgbrfs.f"
		rwork[k] += s;
#line 371 "cgbrfs.f"
/* L70: */
#line 371 "cgbrfs.f"
	    }
#line 372 "cgbrfs.f"
	}
#line 373 "cgbrfs.f"
	s = 0.;
#line 374 "cgbrfs.f"
	i__2 = *n;
#line 374 "cgbrfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 375 "cgbrfs.f"
	    if (rwork[i__] > safe2) {
/* Computing MAX */
#line 376 "cgbrfs.f"
		i__4 = i__;
#line 376 "cgbrfs.f"
		d__3 = s, d__4 = ((d__1 = work[i__4].r, abs(d__1)) + (d__2 = 
			d_imag(&work[i__]), abs(d__2))) / rwork[i__];
#line 376 "cgbrfs.f"
		s = max(d__3,d__4);
#line 377 "cgbrfs.f"
	    } else {
/* Computing MAX */
#line 378 "cgbrfs.f"
		i__4 = i__;
#line 378 "cgbrfs.f"
		d__3 = s, d__4 = ((d__1 = work[i__4].r, abs(d__1)) + (d__2 = 
			d_imag(&work[i__]), abs(d__2)) + safe1) / (rwork[i__] 
			+ safe1);
#line 378 "cgbrfs.f"
		s = max(d__3,d__4);
#line 380 "cgbrfs.f"
	    }
#line 381 "cgbrfs.f"
/* L80: */
#line 381 "cgbrfs.f"
	}
#line 382 "cgbrfs.f"
	berr[j] = s;

/*        Test stopping criterion. Continue iterating if */
/*           1) The residual BERR(J) is larger than machine epsilon, and */
/*           2) BERR(J) decreased by at least a factor of 2 during the */
/*              last iteration, and */
/*           3) At most ITMAX iterations tried. */

#line 390 "cgbrfs.f"
	if (berr[j] > eps && berr[j] * 2. <= lstres && count <= 5) {

/*           Update solution and try again. */

#line 395 "cgbrfs.f"
	    cgbtrs_(trans, n, kl, ku, &c__1, &afb[afb_offset], ldafb, &ipiv[1]
		    , &work[1], n, info, (ftnlen)1);
#line 397 "cgbrfs.f"
	    caxpy_(n, &c_b1, &work[1], &c__1, &x[j * x_dim1 + 1], &c__1);
#line 398 "cgbrfs.f"
	    lstres = berr[j];
#line 399 "cgbrfs.f"
	    ++count;
#line 400 "cgbrfs.f"
	    goto L20;
#line 401 "cgbrfs.f"
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

#line 425 "cgbrfs.f"
	i__2 = *n;
#line 425 "cgbrfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 426 "cgbrfs.f"
	    if (rwork[i__] > safe2) {
#line 427 "cgbrfs.f"
		i__4 = i__;
#line 427 "cgbrfs.f"
		rwork[i__] = (d__1 = work[i__4].r, abs(d__1)) + (d__2 = 
			d_imag(&work[i__]), abs(d__2)) + nz * eps * rwork[i__]
			;
#line 428 "cgbrfs.f"
	    } else {
#line 429 "cgbrfs.f"
		i__4 = i__;
#line 429 "cgbrfs.f"
		rwork[i__] = (d__1 = work[i__4].r, abs(d__1)) + (d__2 = 
			d_imag(&work[i__]), abs(d__2)) + nz * eps * rwork[i__]
			 + safe1;
#line 431 "cgbrfs.f"
	    }
#line 432 "cgbrfs.f"
/* L90: */
#line 432 "cgbrfs.f"
	}

#line 434 "cgbrfs.f"
	kase = 0;
#line 435 "cgbrfs.f"
L100:
#line 436 "cgbrfs.f"
	clacn2_(n, &work[*n + 1], &work[1], &ferr[j], &kase, isave);
#line 437 "cgbrfs.f"
	if (kase != 0) {
#line 438 "cgbrfs.f"
	    if (kase == 1) {

/*              Multiply by diag(W)*inv(op(A)**H). */

#line 442 "cgbrfs.f"
		cgbtrs_(transt, n, kl, ku, &c__1, &afb[afb_offset], ldafb, &
			ipiv[1], &work[1], n, info, (ftnlen)1);
#line 444 "cgbrfs.f"
		i__2 = *n;
#line 444 "cgbrfs.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 445 "cgbrfs.f"
		    i__4 = i__;
#line 445 "cgbrfs.f"
		    i__5 = i__;
#line 445 "cgbrfs.f"
		    i__3 = i__;
#line 445 "cgbrfs.f"
		    z__1.r = rwork[i__5] * work[i__3].r, z__1.i = rwork[i__5] 
			    * work[i__3].i;
#line 445 "cgbrfs.f"
		    work[i__4].r = z__1.r, work[i__4].i = z__1.i;
#line 446 "cgbrfs.f"
/* L110: */
#line 446 "cgbrfs.f"
		}
#line 447 "cgbrfs.f"
	    } else {

/*              Multiply by inv(op(A))*diag(W). */

#line 451 "cgbrfs.f"
		i__2 = *n;
#line 451 "cgbrfs.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 452 "cgbrfs.f"
		    i__4 = i__;
#line 452 "cgbrfs.f"
		    i__5 = i__;
#line 452 "cgbrfs.f"
		    i__3 = i__;
#line 452 "cgbrfs.f"
		    z__1.r = rwork[i__5] * work[i__3].r, z__1.i = rwork[i__5] 
			    * work[i__3].i;
#line 452 "cgbrfs.f"
		    work[i__4].r = z__1.r, work[i__4].i = z__1.i;
#line 453 "cgbrfs.f"
/* L120: */
#line 453 "cgbrfs.f"
		}
#line 454 "cgbrfs.f"
		cgbtrs_(transn, n, kl, ku, &c__1, &afb[afb_offset], ldafb, &
			ipiv[1], &work[1], n, info, (ftnlen)1);
#line 456 "cgbrfs.f"
	    }
#line 457 "cgbrfs.f"
	    goto L100;
#line 458 "cgbrfs.f"
	}

/*        Normalize error. */

#line 462 "cgbrfs.f"
	lstres = 0.;
#line 463 "cgbrfs.f"
	i__2 = *n;
#line 463 "cgbrfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 464 "cgbrfs.f"
	    i__4 = i__ + j * x_dim1;
#line 464 "cgbrfs.f"
	    d__3 = lstres, d__4 = (d__1 = x[i__4].r, abs(d__1)) + (d__2 = 
		    d_imag(&x[i__ + j * x_dim1]), abs(d__2));
#line 464 "cgbrfs.f"
	    lstres = max(d__3,d__4);
#line 465 "cgbrfs.f"
/* L130: */
#line 465 "cgbrfs.f"
	}
#line 466 "cgbrfs.f"
	if (lstres != 0.) {
#line 466 "cgbrfs.f"
	    ferr[j] /= lstres;
#line 466 "cgbrfs.f"
	}

#line 469 "cgbrfs.f"
/* L140: */
#line 469 "cgbrfs.f"
    }

#line 471 "cgbrfs.f"
    return 0;

/*     End of CGBRFS */

} /* cgbrfs_ */

