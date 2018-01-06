#line 1 "cpprfs.f"
/* cpprfs.f -- translated by f2c (version 20100827).
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

#line 1 "cpprfs.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static integer c__1 = 1;

/* > \brief \b CPPRFS */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CPPRFS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cpprfs.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cpprfs.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cpprfs.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CPPRFS( UPLO, N, NRHS, AP, AFP, B, LDB, X, LDX, FERR, */
/*                          BERR, WORK, RWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDB, LDX, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               BERR( * ), FERR( * ), RWORK( * ) */
/*       COMPLEX            AFP( * ), AP( * ), B( LDB, * ), WORK( * ), */
/*      $                   X( LDX, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CPPRFS improves the computed solution to a system of linear */
/* > equations when the coefficient matrix is Hermitian positive definite */
/* > and packed, and provides error bounds and backward error estimates */
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
/* > \param[in] NRHS */
/* > \verbatim */
/* >          NRHS is INTEGER */
/* >          The number of right hand sides, i.e., the number of columns */
/* >          of the matrices B and X.  NRHS >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] AP */
/* > \verbatim */
/* >          AP is COMPLEX array, dimension (N*(N+1)/2) */
/* >          The upper or lower triangle of the Hermitian matrix A, packed */
/* >          columnwise in a linear array.  The j-th column of A is stored */
/* >          in the array AP as follows: */
/* >          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j; */
/* >          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n. */
/* > \endverbatim */
/* > */
/* > \param[in] AFP */
/* > \verbatim */
/* >          AFP is COMPLEX array, dimension (N*(N+1)/2) */
/* >          The triangular factor U or L from the Cholesky factorization */
/* >          A = U**H*U or A = L*L**H, as computed by SPPTRF/CPPTRF, */
/* >          packed columnwise in a linear array in the same format as A */
/* >          (see AP). */
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
/* >          On entry, the solution matrix X, as computed by CPPTRS. */
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

/* > \ingroup complexOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int cpprfs_(char *uplo, integer *n, integer *nrhs, 
	doublecomplex *ap, doublecomplex *afp, doublecomplex *b, integer *ldb,
	 doublecomplex *x, integer *ldx, doublereal *ferr, doublereal *berr, 
	doublecomplex *work, doublereal *rwork, integer *info, ftnlen 
	uplo_len)
{
    /* System generated locals */
    integer b_dim1, b_offset, x_dim1, x_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2, d__3, d__4;
    doublecomplex z__1;

    /* Builtin functions */
    double d_imag(doublecomplex *);

    /* Local variables */
    static integer i__, j, k;
    static doublereal s;
    static integer ik, kk;
    static doublereal xk;
    static integer nz;
    static doublereal eps;
    static integer kase;
    static doublereal safe1, safe2;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer isave[3];
    extern /* Subroutine */ int ccopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), chpmv_(char *, integer *, 
	    doublecomplex *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, ftnlen), caxpy_(
	    integer *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    static integer count;
    static logical upper;
    extern /* Subroutine */ int clacn2_(integer *, doublecomplex *, 
	    doublecomplex *, doublereal *, integer *, integer *);
    extern doublereal slamch_(char *, ftnlen);
    static doublereal safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), cpptrs_(
	    char *, integer *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, integer *, ftnlen);
    static doublereal lstres;


/*  -- LAPACK computational routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

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

#line 233 "cpprfs.f"
    /* Parameter adjustments */
#line 233 "cpprfs.f"
    --ap;
#line 233 "cpprfs.f"
    --afp;
#line 233 "cpprfs.f"
    b_dim1 = *ldb;
#line 233 "cpprfs.f"
    b_offset = 1 + b_dim1;
#line 233 "cpprfs.f"
    b -= b_offset;
#line 233 "cpprfs.f"
    x_dim1 = *ldx;
#line 233 "cpprfs.f"
    x_offset = 1 + x_dim1;
#line 233 "cpprfs.f"
    x -= x_offset;
#line 233 "cpprfs.f"
    --ferr;
#line 233 "cpprfs.f"
    --berr;
#line 233 "cpprfs.f"
    --work;
#line 233 "cpprfs.f"
    --rwork;
#line 233 "cpprfs.f"

#line 233 "cpprfs.f"
    /* Function Body */
#line 233 "cpprfs.f"
    *info = 0;
#line 234 "cpprfs.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 235 "cpprfs.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 236 "cpprfs.f"
	*info = -1;
#line 237 "cpprfs.f"
    } else if (*n < 0) {
#line 238 "cpprfs.f"
	*info = -2;
#line 239 "cpprfs.f"
    } else if (*nrhs < 0) {
#line 240 "cpprfs.f"
	*info = -3;
#line 241 "cpprfs.f"
    } else if (*ldb < max(1,*n)) {
#line 242 "cpprfs.f"
	*info = -7;
#line 243 "cpprfs.f"
    } else if (*ldx < max(1,*n)) {
#line 244 "cpprfs.f"
	*info = -9;
#line 245 "cpprfs.f"
    }
#line 246 "cpprfs.f"
    if (*info != 0) {
#line 247 "cpprfs.f"
	i__1 = -(*info);
#line 247 "cpprfs.f"
	xerbla_("CPPRFS", &i__1, (ftnlen)6);
#line 248 "cpprfs.f"
	return 0;
#line 249 "cpprfs.f"
    }

/*     Quick return if possible */

#line 253 "cpprfs.f"
    if (*n == 0 || *nrhs == 0) {
#line 254 "cpprfs.f"
	i__1 = *nrhs;
#line 254 "cpprfs.f"
	for (j = 1; j <= i__1; ++j) {
#line 255 "cpprfs.f"
	    ferr[j] = 0.;
#line 256 "cpprfs.f"
	    berr[j] = 0.;
#line 257 "cpprfs.f"
/* L10: */
#line 257 "cpprfs.f"
	}
#line 258 "cpprfs.f"
	return 0;
#line 259 "cpprfs.f"
    }

/*     NZ = maximum number of nonzero elements in each row of A, plus 1 */

#line 263 "cpprfs.f"
    nz = *n + 1;
#line 264 "cpprfs.f"
    eps = slamch_("Epsilon", (ftnlen)7);
#line 265 "cpprfs.f"
    safmin = slamch_("Safe minimum", (ftnlen)12);
#line 266 "cpprfs.f"
    safe1 = nz * safmin;
#line 267 "cpprfs.f"
    safe2 = safe1 / eps;

/*     Do for each right hand side */

#line 271 "cpprfs.f"
    i__1 = *nrhs;
#line 271 "cpprfs.f"
    for (j = 1; j <= i__1; ++j) {

#line 273 "cpprfs.f"
	count = 1;
#line 274 "cpprfs.f"
	lstres = 3.;
#line 275 "cpprfs.f"
L20:

/*        Loop until stopping criterion is satisfied. */

/*        Compute residual R = B - A * X */

#line 281 "cpprfs.f"
	ccopy_(n, &b[j * b_dim1 + 1], &c__1, &work[1], &c__1);
#line 282 "cpprfs.f"
	z__1.r = -1., z__1.i = -0.;
#line 282 "cpprfs.f"
	chpmv_(uplo, n, &z__1, &ap[1], &x[j * x_dim1 + 1], &c__1, &c_b1, &
		work[1], &c__1, (ftnlen)1);

/*        Compute componentwise relative backward error from formula */

/*        max(i) ( abs(R(i)) / ( abs(A)*abs(X) + abs(B) )(i) ) */

/*        where abs(Z) is the componentwise absolute value of the matrix */
/*        or vector Z.  If the i-th component of the denominator is less */
/*        than SAFE2, then SAFE1 is added to the i-th components of the */
/*        numerator and denominator before dividing. */

#line 293 "cpprfs.f"
	i__2 = *n;
#line 293 "cpprfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 294 "cpprfs.f"
	    i__3 = i__ + j * b_dim1;
#line 294 "cpprfs.f"
	    rwork[i__] = (d__1 = b[i__3].r, abs(d__1)) + (d__2 = d_imag(&b[
		    i__ + j * b_dim1]), abs(d__2));
#line 295 "cpprfs.f"
/* L30: */
#line 295 "cpprfs.f"
	}

/*        Compute abs(A)*abs(X) + abs(B). */

#line 299 "cpprfs.f"
	kk = 1;
#line 300 "cpprfs.f"
	if (upper) {
#line 301 "cpprfs.f"
	    i__2 = *n;
#line 301 "cpprfs.f"
	    for (k = 1; k <= i__2; ++k) {
#line 302 "cpprfs.f"
		s = 0.;
#line 303 "cpprfs.f"
		i__3 = k + j * x_dim1;
#line 303 "cpprfs.f"
		xk = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = d_imag(&x[k + j *
			 x_dim1]), abs(d__2));
#line 304 "cpprfs.f"
		ik = kk;
#line 305 "cpprfs.f"
		i__3 = k - 1;
#line 305 "cpprfs.f"
		for (i__ = 1; i__ <= i__3; ++i__) {
#line 306 "cpprfs.f"
		    i__4 = ik;
#line 306 "cpprfs.f"
		    rwork[i__] += ((d__1 = ap[i__4].r, abs(d__1)) + (d__2 = 
			    d_imag(&ap[ik]), abs(d__2))) * xk;
#line 307 "cpprfs.f"
		    i__4 = ik;
#line 307 "cpprfs.f"
		    i__5 = i__ + j * x_dim1;
#line 307 "cpprfs.f"
		    s += ((d__1 = ap[i__4].r, abs(d__1)) + (d__2 = d_imag(&ap[
			    ik]), abs(d__2))) * ((d__3 = x[i__5].r, abs(d__3))
			     + (d__4 = d_imag(&x[i__ + j * x_dim1]), abs(d__4)
			    ));
#line 308 "cpprfs.f"
		    ++ik;
#line 309 "cpprfs.f"
/* L40: */
#line 309 "cpprfs.f"
		}
#line 310 "cpprfs.f"
		i__3 = kk + k - 1;
#line 310 "cpprfs.f"
		rwork[k] = rwork[k] + (d__1 = ap[i__3].r, abs(d__1)) * xk + s;
#line 312 "cpprfs.f"
		kk += k;
#line 313 "cpprfs.f"
/* L50: */
#line 313 "cpprfs.f"
	    }
#line 314 "cpprfs.f"
	} else {
#line 315 "cpprfs.f"
	    i__2 = *n;
#line 315 "cpprfs.f"
	    for (k = 1; k <= i__2; ++k) {
#line 316 "cpprfs.f"
		s = 0.;
#line 317 "cpprfs.f"
		i__3 = k + j * x_dim1;
#line 317 "cpprfs.f"
		xk = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = d_imag(&x[k + j *
			 x_dim1]), abs(d__2));
#line 318 "cpprfs.f"
		i__3 = kk;
#line 318 "cpprfs.f"
		rwork[k] += (d__1 = ap[i__3].r, abs(d__1)) * xk;
#line 319 "cpprfs.f"
		ik = kk + 1;
#line 320 "cpprfs.f"
		i__3 = *n;
#line 320 "cpprfs.f"
		for (i__ = k + 1; i__ <= i__3; ++i__) {
#line 321 "cpprfs.f"
		    i__4 = ik;
#line 321 "cpprfs.f"
		    rwork[i__] += ((d__1 = ap[i__4].r, abs(d__1)) + (d__2 = 
			    d_imag(&ap[ik]), abs(d__2))) * xk;
#line 322 "cpprfs.f"
		    i__4 = ik;
#line 322 "cpprfs.f"
		    i__5 = i__ + j * x_dim1;
#line 322 "cpprfs.f"
		    s += ((d__1 = ap[i__4].r, abs(d__1)) + (d__2 = d_imag(&ap[
			    ik]), abs(d__2))) * ((d__3 = x[i__5].r, abs(d__3))
			     + (d__4 = d_imag(&x[i__ + j * x_dim1]), abs(d__4)
			    ));
#line 323 "cpprfs.f"
		    ++ik;
#line 324 "cpprfs.f"
/* L60: */
#line 324 "cpprfs.f"
		}
#line 325 "cpprfs.f"
		rwork[k] += s;
#line 326 "cpprfs.f"
		kk += *n - k + 1;
#line 327 "cpprfs.f"
/* L70: */
#line 327 "cpprfs.f"
	    }
#line 328 "cpprfs.f"
	}
#line 329 "cpprfs.f"
	s = 0.;
#line 330 "cpprfs.f"
	i__2 = *n;
#line 330 "cpprfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 331 "cpprfs.f"
	    if (rwork[i__] > safe2) {
/* Computing MAX */
#line 332 "cpprfs.f"
		i__3 = i__;
#line 332 "cpprfs.f"
		d__3 = s, d__4 = ((d__1 = work[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&work[i__]), abs(d__2))) / rwork[i__];
#line 332 "cpprfs.f"
		s = max(d__3,d__4);
#line 333 "cpprfs.f"
	    } else {
/* Computing MAX */
#line 334 "cpprfs.f"
		i__3 = i__;
#line 334 "cpprfs.f"
		d__3 = s, d__4 = ((d__1 = work[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&work[i__]), abs(d__2)) + safe1) / (rwork[i__] 
			+ safe1);
#line 334 "cpprfs.f"
		s = max(d__3,d__4);
#line 336 "cpprfs.f"
	    }
#line 337 "cpprfs.f"
/* L80: */
#line 337 "cpprfs.f"
	}
#line 338 "cpprfs.f"
	berr[j] = s;

/*        Test stopping criterion. Continue iterating if */
/*           1) The residual BERR(J) is larger than machine epsilon, and */
/*           2) BERR(J) decreased by at least a factor of 2 during the */
/*              last iteration, and */
/*           3) At most ITMAX iterations tried. */

#line 346 "cpprfs.f"
	if (berr[j] > eps && berr[j] * 2. <= lstres && count <= 5) {

/*           Update solution and try again. */

#line 351 "cpprfs.f"
	    cpptrs_(uplo, n, &c__1, &afp[1], &work[1], n, info, (ftnlen)1);
#line 352 "cpprfs.f"
	    caxpy_(n, &c_b1, &work[1], &c__1, &x[j * x_dim1 + 1], &c__1);
#line 353 "cpprfs.f"
	    lstres = berr[j];
#line 354 "cpprfs.f"
	    ++count;
#line 355 "cpprfs.f"
	    goto L20;
#line 356 "cpprfs.f"
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

#line 380 "cpprfs.f"
	i__2 = *n;
#line 380 "cpprfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 381 "cpprfs.f"
	    if (rwork[i__] > safe2) {
#line 382 "cpprfs.f"
		i__3 = i__;
#line 382 "cpprfs.f"
		rwork[i__] = (d__1 = work[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&work[i__]), abs(d__2)) + nz * eps * rwork[i__]
			;
#line 383 "cpprfs.f"
	    } else {
#line 384 "cpprfs.f"
		i__3 = i__;
#line 384 "cpprfs.f"
		rwork[i__] = (d__1 = work[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&work[i__]), abs(d__2)) + nz * eps * rwork[i__]
			 + safe1;
#line 386 "cpprfs.f"
	    }
#line 387 "cpprfs.f"
/* L90: */
#line 387 "cpprfs.f"
	}

#line 389 "cpprfs.f"
	kase = 0;
#line 390 "cpprfs.f"
L100:
#line 391 "cpprfs.f"
	clacn2_(n, &work[*n + 1], &work[1], &ferr[j], &kase, isave);
#line 392 "cpprfs.f"
	if (kase != 0) {
#line 393 "cpprfs.f"
	    if (kase == 1) {

/*              Multiply by diag(W)*inv(A**H). */

#line 397 "cpprfs.f"
		cpptrs_(uplo, n, &c__1, &afp[1], &work[1], n, info, (ftnlen)1)
			;
#line 398 "cpprfs.f"
		i__2 = *n;
#line 398 "cpprfs.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 399 "cpprfs.f"
		    i__3 = i__;
#line 399 "cpprfs.f"
		    i__4 = i__;
#line 399 "cpprfs.f"
		    i__5 = i__;
#line 399 "cpprfs.f"
		    z__1.r = rwork[i__4] * work[i__5].r, z__1.i = rwork[i__4] 
			    * work[i__5].i;
#line 399 "cpprfs.f"
		    work[i__3].r = z__1.r, work[i__3].i = z__1.i;
#line 400 "cpprfs.f"
/* L110: */
#line 400 "cpprfs.f"
		}
#line 401 "cpprfs.f"
	    } else if (kase == 2) {

/*              Multiply by inv(A)*diag(W). */

#line 405 "cpprfs.f"
		i__2 = *n;
#line 405 "cpprfs.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 406 "cpprfs.f"
		    i__3 = i__;
#line 406 "cpprfs.f"
		    i__4 = i__;
#line 406 "cpprfs.f"
		    i__5 = i__;
#line 406 "cpprfs.f"
		    z__1.r = rwork[i__4] * work[i__5].r, z__1.i = rwork[i__4] 
			    * work[i__5].i;
#line 406 "cpprfs.f"
		    work[i__3].r = z__1.r, work[i__3].i = z__1.i;
#line 407 "cpprfs.f"
/* L120: */
#line 407 "cpprfs.f"
		}
#line 408 "cpprfs.f"
		cpptrs_(uplo, n, &c__1, &afp[1], &work[1], n, info, (ftnlen)1)
			;
#line 409 "cpprfs.f"
	    }
#line 410 "cpprfs.f"
	    goto L100;
#line 411 "cpprfs.f"
	}

/*        Normalize error. */

#line 415 "cpprfs.f"
	lstres = 0.;
#line 416 "cpprfs.f"
	i__2 = *n;
#line 416 "cpprfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 417 "cpprfs.f"
	    i__3 = i__ + j * x_dim1;
#line 417 "cpprfs.f"
	    d__3 = lstres, d__4 = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = 
		    d_imag(&x[i__ + j * x_dim1]), abs(d__2));
#line 417 "cpprfs.f"
	    lstres = max(d__3,d__4);
#line 418 "cpprfs.f"
/* L130: */
#line 418 "cpprfs.f"
	}
#line 419 "cpprfs.f"
	if (lstres != 0.) {
#line 419 "cpprfs.f"
	    ferr[j] /= lstres;
#line 419 "cpprfs.f"
	}

#line 422 "cpprfs.f"
/* L140: */
#line 422 "cpprfs.f"
    }

#line 424 "cpprfs.f"
    return 0;

/*     End of CPPRFS */

} /* cpprfs_ */

