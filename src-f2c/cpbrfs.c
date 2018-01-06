#line 1 "cpbrfs.f"
/* cpbrfs.f -- translated by f2c (version 20100827).
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

#line 1 "cpbrfs.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static integer c__1 = 1;

/* > \brief \b CPBRFS */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CPBRFS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cpbrfs.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cpbrfs.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cpbrfs.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CPBRFS( UPLO, N, KD, NRHS, AB, LDAB, AFB, LDAFB, B, */
/*                          LDB, X, LDX, FERR, BERR, WORK, RWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, KD, LDAB, LDAFB, LDB, LDX, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               BERR( * ), FERR( * ), RWORK( * ) */
/*       COMPLEX            AB( LDAB, * ), AFB( LDAFB, * ), B( LDB, * ), */
/*      $                   WORK( * ), X( LDX, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CPBRFS improves the computed solution to a system of linear */
/* > equations when the coefficient matrix is Hermitian positive definite */
/* > and banded, and provides error bounds and backward error estimates */
/* > for the solution. */
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
/* > \param[in] KD */
/* > \verbatim */
/* >          KD is INTEGER */
/* >          The number of superdiagonals of the matrix A if UPLO = 'U', */
/* >          or the number of subdiagonals if UPLO = 'L'.  KD >= 0. */
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
/* >          The upper or lower triangle of the Hermitian band matrix A, */
/* >          stored in the first KD+1 rows of the array.  The j-th column */
/* >          of A is stored in the j-th column of the array AB as follows: */
/* >          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j; */
/* >          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd). */
/* > \endverbatim */
/* > */
/* > \param[in] LDAB */
/* > \verbatim */
/* >          LDAB is INTEGER */
/* >          The leading dimension of the array AB.  LDAB >= KD+1. */
/* > \endverbatim */
/* > */
/* > \param[in] AFB */
/* > \verbatim */
/* >          AFB is COMPLEX array, dimension (LDAFB,N) */
/* >          The triangular factor U or L from the Cholesky factorization */
/* >          A = U**H*U or A = L*L**H of the band matrix A as computed by */
/* >          CPBTRF, in the same storage format as A (see AB). */
/* > \endverbatim */
/* > */
/* > \param[in] LDAFB */
/* > \verbatim */
/* >          LDAFB is INTEGER */
/* >          The leading dimension of the array AFB.  LDAFB >= KD+1. */
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
/* >          On entry, the solution matrix X, as computed by CPBTRS. */
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

/* > \ingroup complexOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int cpbrfs_(char *uplo, integer *n, integer *kd, integer *
	nrhs, doublecomplex *ab, integer *ldab, doublecomplex *afb, integer *
	ldafb, doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx,
	 doublereal *ferr, doublereal *berr, doublecomplex *work, doublereal *
	rwork, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, afb_dim1, afb_offset, b_dim1, b_offset, 
	    x_dim1, x_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2, d__3, d__4;
    doublecomplex z__1;

    /* Builtin functions */
    double d_imag(doublecomplex *);

    /* Local variables */
    static integer i__, j, k, l;
    static doublereal s, xk;
    static integer nz;
    static doublereal eps;
    static integer kase;
    static doublereal safe1, safe2;
    extern /* Subroutine */ int chbmv_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
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
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), cpbtrs_(
	    char *, integer *, integer *, integer *, doublecomplex *, integer 
	    *, doublecomplex *, integer *, integer *, ftnlen);
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

#line 251 "cpbrfs.f"
    /* Parameter adjustments */
#line 251 "cpbrfs.f"
    ab_dim1 = *ldab;
#line 251 "cpbrfs.f"
    ab_offset = 1 + ab_dim1;
#line 251 "cpbrfs.f"
    ab -= ab_offset;
#line 251 "cpbrfs.f"
    afb_dim1 = *ldafb;
#line 251 "cpbrfs.f"
    afb_offset = 1 + afb_dim1;
#line 251 "cpbrfs.f"
    afb -= afb_offset;
#line 251 "cpbrfs.f"
    b_dim1 = *ldb;
#line 251 "cpbrfs.f"
    b_offset = 1 + b_dim1;
#line 251 "cpbrfs.f"
    b -= b_offset;
#line 251 "cpbrfs.f"
    x_dim1 = *ldx;
#line 251 "cpbrfs.f"
    x_offset = 1 + x_dim1;
#line 251 "cpbrfs.f"
    x -= x_offset;
#line 251 "cpbrfs.f"
    --ferr;
#line 251 "cpbrfs.f"
    --berr;
#line 251 "cpbrfs.f"
    --work;
#line 251 "cpbrfs.f"
    --rwork;
#line 251 "cpbrfs.f"

#line 251 "cpbrfs.f"
    /* Function Body */
#line 251 "cpbrfs.f"
    *info = 0;
#line 252 "cpbrfs.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 253 "cpbrfs.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 254 "cpbrfs.f"
	*info = -1;
#line 255 "cpbrfs.f"
    } else if (*n < 0) {
#line 256 "cpbrfs.f"
	*info = -2;
#line 257 "cpbrfs.f"
    } else if (*kd < 0) {
#line 258 "cpbrfs.f"
	*info = -3;
#line 259 "cpbrfs.f"
    } else if (*nrhs < 0) {
#line 260 "cpbrfs.f"
	*info = -4;
#line 261 "cpbrfs.f"
    } else if (*ldab < *kd + 1) {
#line 262 "cpbrfs.f"
	*info = -6;
#line 263 "cpbrfs.f"
    } else if (*ldafb < *kd + 1) {
#line 264 "cpbrfs.f"
	*info = -8;
#line 265 "cpbrfs.f"
    } else if (*ldb < max(1,*n)) {
#line 266 "cpbrfs.f"
	*info = -10;
#line 267 "cpbrfs.f"
    } else if (*ldx < max(1,*n)) {
#line 268 "cpbrfs.f"
	*info = -12;
#line 269 "cpbrfs.f"
    }
#line 270 "cpbrfs.f"
    if (*info != 0) {
#line 271 "cpbrfs.f"
	i__1 = -(*info);
#line 271 "cpbrfs.f"
	xerbla_("CPBRFS", &i__1, (ftnlen)6);
#line 272 "cpbrfs.f"
	return 0;
#line 273 "cpbrfs.f"
    }

/*     Quick return if possible */

#line 277 "cpbrfs.f"
    if (*n == 0 || *nrhs == 0) {
#line 278 "cpbrfs.f"
	i__1 = *nrhs;
#line 278 "cpbrfs.f"
	for (j = 1; j <= i__1; ++j) {
#line 279 "cpbrfs.f"
	    ferr[j] = 0.;
#line 280 "cpbrfs.f"
	    berr[j] = 0.;
#line 281 "cpbrfs.f"
/* L10: */
#line 281 "cpbrfs.f"
	}
#line 282 "cpbrfs.f"
	return 0;
#line 283 "cpbrfs.f"
    }

/*     NZ = maximum number of nonzero elements in each row of A, plus 1 */

/* Computing MIN */
#line 287 "cpbrfs.f"
    i__1 = *n + 1, i__2 = (*kd << 1) + 2;
#line 287 "cpbrfs.f"
    nz = min(i__1,i__2);
#line 288 "cpbrfs.f"
    eps = slamch_("Epsilon", (ftnlen)7);
#line 289 "cpbrfs.f"
    safmin = slamch_("Safe minimum", (ftnlen)12);
#line 290 "cpbrfs.f"
    safe1 = nz * safmin;
#line 291 "cpbrfs.f"
    safe2 = safe1 / eps;

/*     Do for each right hand side */

#line 295 "cpbrfs.f"
    i__1 = *nrhs;
#line 295 "cpbrfs.f"
    for (j = 1; j <= i__1; ++j) {

#line 297 "cpbrfs.f"
	count = 1;
#line 298 "cpbrfs.f"
	lstres = 3.;
#line 299 "cpbrfs.f"
L20:

/*        Loop until stopping criterion is satisfied. */

/*        Compute residual R = B - A * X */

#line 305 "cpbrfs.f"
	ccopy_(n, &b[j * b_dim1 + 1], &c__1, &work[1], &c__1);
#line 306 "cpbrfs.f"
	z__1.r = -1., z__1.i = -0.;
#line 306 "cpbrfs.f"
	chbmv_(uplo, n, kd, &z__1, &ab[ab_offset], ldab, &x[j * x_dim1 + 1], &
		c__1, &c_b1, &work[1], &c__1, (ftnlen)1);

/*        Compute componentwise relative backward error from formula */

/*        max(i) ( abs(R(i)) / ( abs(A)*abs(X) + abs(B) )(i) ) */

/*        where abs(Z) is the componentwise absolute value of the matrix */
/*        or vector Z.  If the i-th component of the denominator is less */
/*        than SAFE2, then SAFE1 is added to the i-th components of the */
/*        numerator and denominator before dividing. */

#line 318 "cpbrfs.f"
	i__2 = *n;
#line 318 "cpbrfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 319 "cpbrfs.f"
	    i__3 = i__ + j * b_dim1;
#line 319 "cpbrfs.f"
	    rwork[i__] = (d__1 = b[i__3].r, abs(d__1)) + (d__2 = d_imag(&b[
		    i__ + j * b_dim1]), abs(d__2));
#line 320 "cpbrfs.f"
/* L30: */
#line 320 "cpbrfs.f"
	}

/*        Compute abs(A)*abs(X) + abs(B). */

#line 324 "cpbrfs.f"
	if (upper) {
#line 325 "cpbrfs.f"
	    i__2 = *n;
#line 325 "cpbrfs.f"
	    for (k = 1; k <= i__2; ++k) {
#line 326 "cpbrfs.f"
		s = 0.;
#line 327 "cpbrfs.f"
		i__3 = k + j * x_dim1;
#line 327 "cpbrfs.f"
		xk = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = d_imag(&x[k + j *
			 x_dim1]), abs(d__2));
#line 328 "cpbrfs.f"
		l = *kd + 1 - k;
/* Computing MAX */
#line 329 "cpbrfs.f"
		i__3 = 1, i__4 = k - *kd;
#line 329 "cpbrfs.f"
		i__5 = k - 1;
#line 329 "cpbrfs.f"
		for (i__ = max(i__3,i__4); i__ <= i__5; ++i__) {
#line 330 "cpbrfs.f"
		    i__3 = l + i__ + k * ab_dim1;
#line 330 "cpbrfs.f"
		    rwork[i__] += ((d__1 = ab[i__3].r, abs(d__1)) + (d__2 = 
			    d_imag(&ab[l + i__ + k * ab_dim1]), abs(d__2))) * 
			    xk;
#line 331 "cpbrfs.f"
		    i__3 = l + i__ + k * ab_dim1;
#line 331 "cpbrfs.f"
		    i__4 = i__ + j * x_dim1;
#line 331 "cpbrfs.f"
		    s += ((d__1 = ab[i__3].r, abs(d__1)) + (d__2 = d_imag(&ab[
			    l + i__ + k * ab_dim1]), abs(d__2))) * ((d__3 = x[
			    i__4].r, abs(d__3)) + (d__4 = d_imag(&x[i__ + j * 
			    x_dim1]), abs(d__4)));
#line 332 "cpbrfs.f"
/* L40: */
#line 332 "cpbrfs.f"
		}
#line 333 "cpbrfs.f"
		i__5 = *kd + 1 + k * ab_dim1;
#line 333 "cpbrfs.f"
		rwork[k] = rwork[k] + (d__1 = ab[i__5].r, abs(d__1)) * xk + s;
#line 335 "cpbrfs.f"
/* L50: */
#line 335 "cpbrfs.f"
	    }
#line 336 "cpbrfs.f"
	} else {
#line 337 "cpbrfs.f"
	    i__2 = *n;
#line 337 "cpbrfs.f"
	    for (k = 1; k <= i__2; ++k) {
#line 338 "cpbrfs.f"
		s = 0.;
#line 339 "cpbrfs.f"
		i__5 = k + j * x_dim1;
#line 339 "cpbrfs.f"
		xk = (d__1 = x[i__5].r, abs(d__1)) + (d__2 = d_imag(&x[k + j *
			 x_dim1]), abs(d__2));
#line 340 "cpbrfs.f"
		i__5 = k * ab_dim1 + 1;
#line 340 "cpbrfs.f"
		rwork[k] += (d__1 = ab[i__5].r, abs(d__1)) * xk;
#line 341 "cpbrfs.f"
		l = 1 - k;
/* Computing MIN */
#line 342 "cpbrfs.f"
		i__3 = *n, i__4 = k + *kd;
#line 342 "cpbrfs.f"
		i__5 = min(i__3,i__4);
#line 342 "cpbrfs.f"
		for (i__ = k + 1; i__ <= i__5; ++i__) {
#line 343 "cpbrfs.f"
		    i__3 = l + i__ + k * ab_dim1;
#line 343 "cpbrfs.f"
		    rwork[i__] += ((d__1 = ab[i__3].r, abs(d__1)) + (d__2 = 
			    d_imag(&ab[l + i__ + k * ab_dim1]), abs(d__2))) * 
			    xk;
#line 344 "cpbrfs.f"
		    i__3 = l + i__ + k * ab_dim1;
#line 344 "cpbrfs.f"
		    i__4 = i__ + j * x_dim1;
#line 344 "cpbrfs.f"
		    s += ((d__1 = ab[i__3].r, abs(d__1)) + (d__2 = d_imag(&ab[
			    l + i__ + k * ab_dim1]), abs(d__2))) * ((d__3 = x[
			    i__4].r, abs(d__3)) + (d__4 = d_imag(&x[i__ + j * 
			    x_dim1]), abs(d__4)));
#line 345 "cpbrfs.f"
/* L60: */
#line 345 "cpbrfs.f"
		}
#line 346 "cpbrfs.f"
		rwork[k] += s;
#line 347 "cpbrfs.f"
/* L70: */
#line 347 "cpbrfs.f"
	    }
#line 348 "cpbrfs.f"
	}
#line 349 "cpbrfs.f"
	s = 0.;
#line 350 "cpbrfs.f"
	i__2 = *n;
#line 350 "cpbrfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 351 "cpbrfs.f"
	    if (rwork[i__] > safe2) {
/* Computing MAX */
#line 352 "cpbrfs.f"
		i__5 = i__;
#line 352 "cpbrfs.f"
		d__3 = s, d__4 = ((d__1 = work[i__5].r, abs(d__1)) + (d__2 = 
			d_imag(&work[i__]), abs(d__2))) / rwork[i__];
#line 352 "cpbrfs.f"
		s = max(d__3,d__4);
#line 353 "cpbrfs.f"
	    } else {
/* Computing MAX */
#line 354 "cpbrfs.f"
		i__5 = i__;
#line 354 "cpbrfs.f"
		d__3 = s, d__4 = ((d__1 = work[i__5].r, abs(d__1)) + (d__2 = 
			d_imag(&work[i__]), abs(d__2)) + safe1) / (rwork[i__] 
			+ safe1);
#line 354 "cpbrfs.f"
		s = max(d__3,d__4);
#line 356 "cpbrfs.f"
	    }
#line 357 "cpbrfs.f"
/* L80: */
#line 357 "cpbrfs.f"
	}
#line 358 "cpbrfs.f"
	berr[j] = s;

/*        Test stopping criterion. Continue iterating if */
/*           1) The residual BERR(J) is larger than machine epsilon, and */
/*           2) BERR(J) decreased by at least a factor of 2 during the */
/*              last iteration, and */
/*           3) At most ITMAX iterations tried. */

#line 366 "cpbrfs.f"
	if (berr[j] > eps && berr[j] * 2. <= lstres && count <= 5) {

/*           Update solution and try again. */

#line 371 "cpbrfs.f"
	    cpbtrs_(uplo, n, kd, &c__1, &afb[afb_offset], ldafb, &work[1], n, 
		    info, (ftnlen)1);
#line 372 "cpbrfs.f"
	    caxpy_(n, &c_b1, &work[1], &c__1, &x[j * x_dim1 + 1], &c__1);
#line 373 "cpbrfs.f"
	    lstres = berr[j];
#line 374 "cpbrfs.f"
	    ++count;
#line 375 "cpbrfs.f"
	    goto L20;
#line 376 "cpbrfs.f"
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

#line 400 "cpbrfs.f"
	i__2 = *n;
#line 400 "cpbrfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 401 "cpbrfs.f"
	    if (rwork[i__] > safe2) {
#line 402 "cpbrfs.f"
		i__5 = i__;
#line 402 "cpbrfs.f"
		rwork[i__] = (d__1 = work[i__5].r, abs(d__1)) + (d__2 = 
			d_imag(&work[i__]), abs(d__2)) + nz * eps * rwork[i__]
			;
#line 403 "cpbrfs.f"
	    } else {
#line 404 "cpbrfs.f"
		i__5 = i__;
#line 404 "cpbrfs.f"
		rwork[i__] = (d__1 = work[i__5].r, abs(d__1)) + (d__2 = 
			d_imag(&work[i__]), abs(d__2)) + nz * eps * rwork[i__]
			 + safe1;
#line 406 "cpbrfs.f"
	    }
#line 407 "cpbrfs.f"
/* L90: */
#line 407 "cpbrfs.f"
	}

#line 409 "cpbrfs.f"
	kase = 0;
#line 410 "cpbrfs.f"
L100:
#line 411 "cpbrfs.f"
	clacn2_(n, &work[*n + 1], &work[1], &ferr[j], &kase, isave);
#line 412 "cpbrfs.f"
	if (kase != 0) {
#line 413 "cpbrfs.f"
	    if (kase == 1) {

/*              Multiply by diag(W)*inv(A**H). */

#line 417 "cpbrfs.f"
		cpbtrs_(uplo, n, kd, &c__1, &afb[afb_offset], ldafb, &work[1],
			 n, info, (ftnlen)1);
#line 418 "cpbrfs.f"
		i__2 = *n;
#line 418 "cpbrfs.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 419 "cpbrfs.f"
		    i__5 = i__;
#line 419 "cpbrfs.f"
		    i__3 = i__;
#line 419 "cpbrfs.f"
		    i__4 = i__;
#line 419 "cpbrfs.f"
		    z__1.r = rwork[i__3] * work[i__4].r, z__1.i = rwork[i__3] 
			    * work[i__4].i;
#line 419 "cpbrfs.f"
		    work[i__5].r = z__1.r, work[i__5].i = z__1.i;
#line 420 "cpbrfs.f"
/* L110: */
#line 420 "cpbrfs.f"
		}
#line 421 "cpbrfs.f"
	    } else if (kase == 2) {

/*              Multiply by inv(A)*diag(W). */

#line 425 "cpbrfs.f"
		i__2 = *n;
#line 425 "cpbrfs.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 426 "cpbrfs.f"
		    i__5 = i__;
#line 426 "cpbrfs.f"
		    i__3 = i__;
#line 426 "cpbrfs.f"
		    i__4 = i__;
#line 426 "cpbrfs.f"
		    z__1.r = rwork[i__3] * work[i__4].r, z__1.i = rwork[i__3] 
			    * work[i__4].i;
#line 426 "cpbrfs.f"
		    work[i__5].r = z__1.r, work[i__5].i = z__1.i;
#line 427 "cpbrfs.f"
/* L120: */
#line 427 "cpbrfs.f"
		}
#line 428 "cpbrfs.f"
		cpbtrs_(uplo, n, kd, &c__1, &afb[afb_offset], ldafb, &work[1],
			 n, info, (ftnlen)1);
#line 429 "cpbrfs.f"
	    }
#line 430 "cpbrfs.f"
	    goto L100;
#line 431 "cpbrfs.f"
	}

/*        Normalize error. */

#line 435 "cpbrfs.f"
	lstres = 0.;
#line 436 "cpbrfs.f"
	i__2 = *n;
#line 436 "cpbrfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 437 "cpbrfs.f"
	    i__5 = i__ + j * x_dim1;
#line 437 "cpbrfs.f"
	    d__3 = lstres, d__4 = (d__1 = x[i__5].r, abs(d__1)) + (d__2 = 
		    d_imag(&x[i__ + j * x_dim1]), abs(d__2));
#line 437 "cpbrfs.f"
	    lstres = max(d__3,d__4);
#line 438 "cpbrfs.f"
/* L130: */
#line 438 "cpbrfs.f"
	}
#line 439 "cpbrfs.f"
	if (lstres != 0.) {
#line 439 "cpbrfs.f"
	    ferr[j] /= lstres;
#line 439 "cpbrfs.f"
	}

#line 442 "cpbrfs.f"
/* L140: */
#line 442 "cpbrfs.f"
    }

#line 444 "cpbrfs.f"
    return 0;

/*     End of CPBRFS */

} /* cpbrfs_ */

