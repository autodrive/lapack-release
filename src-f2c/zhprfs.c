#line 1 "zhprfs.f"
/* zhprfs.f -- translated by f2c (version 20100827).
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

#line 1 "zhprfs.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static integer c__1 = 1;

/* > \brief \b ZHPRFS */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZHPRFS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhprfs.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhprfs.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhprfs.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZHPRFS( UPLO, N, NRHS, AP, AFP, IPIV, B, LDB, X, LDX, */
/*                          FERR, BERR, WORK, RWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDB, LDX, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       DOUBLE PRECISION   BERR( * ), FERR( * ), RWORK( * ) */
/*       COMPLEX*16         AFP( * ), AP( * ), B( LDB, * ), WORK( * ), */
/*      $                   X( LDX, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZHPRFS improves the computed solution to a system of linear */
/* > equations when the coefficient matrix is Hermitian indefinite */
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
/* >          AP is COMPLEX*16 array, dimension (N*(N+1)/2) */
/* >          The upper or lower triangle of the Hermitian matrix A, packed */
/* >          columnwise in a linear array.  The j-th column of A is stored */
/* >          in the array AP as follows: */
/* >          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j; */
/* >          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n. */
/* > \endverbatim */
/* > */
/* > \param[in] AFP */
/* > \verbatim */
/* >          AFP is COMPLEX*16 array, dimension (N*(N+1)/2) */
/* >          The factored form of the matrix A.  AFP contains the block */
/* >          diagonal matrix D and the multipliers used to obtain the */
/* >          factor U or L from the factorization A = U*D*U**H or */
/* >          A = L*D*L**H as computed by ZHPTRF, stored as a packed */
/* >          triangular matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N) */
/* >          Details of the interchanges and the block structure of D */
/* >          as determined by ZHPTRF. */
/* > \endverbatim */
/* > */
/* > \param[in] B */
/* > \verbatim */
/* >          B is COMPLEX*16 array, dimension (LDB,NRHS) */
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
/* >          X is COMPLEX*16 array, dimension (LDX,NRHS) */
/* >          On entry, the solution matrix X, as computed by ZHPTRS. */
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
/* >          WORK is COMPLEX*16 array, dimension (2*N) */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is DOUBLE PRECISION array, dimension (N) */
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

/* > \ingroup complex16OTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int zhprfs_(char *uplo, integer *n, integer *nrhs, 
	doublecomplex *ap, doublecomplex *afp, integer *ipiv, doublecomplex *
	b, integer *ldb, doublecomplex *x, integer *ldx, doublereal *ferr, 
	doublereal *berr, doublecomplex *work, doublereal *rwork, integer *
	info, ftnlen uplo_len)
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
    static integer isave[3], count;
    static logical upper;
    extern /* Subroutine */ int zcopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), zhpmv_(char *, integer *, 
	    doublecomplex *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, ftnlen), zaxpy_(
	    integer *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), zlacn2_(integer *, doublecomplex *, 
	    doublecomplex *, doublereal *, integer *, integer *);
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal lstres;
    extern /* Subroutine */ int zhptrs_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     ftnlen);


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

#line 243 "zhprfs.f"
    /* Parameter adjustments */
#line 243 "zhprfs.f"
    --ap;
#line 243 "zhprfs.f"
    --afp;
#line 243 "zhprfs.f"
    --ipiv;
#line 243 "zhprfs.f"
    b_dim1 = *ldb;
#line 243 "zhprfs.f"
    b_offset = 1 + b_dim1;
#line 243 "zhprfs.f"
    b -= b_offset;
#line 243 "zhprfs.f"
    x_dim1 = *ldx;
#line 243 "zhprfs.f"
    x_offset = 1 + x_dim1;
#line 243 "zhprfs.f"
    x -= x_offset;
#line 243 "zhprfs.f"
    --ferr;
#line 243 "zhprfs.f"
    --berr;
#line 243 "zhprfs.f"
    --work;
#line 243 "zhprfs.f"
    --rwork;
#line 243 "zhprfs.f"

#line 243 "zhprfs.f"
    /* Function Body */
#line 243 "zhprfs.f"
    *info = 0;
#line 244 "zhprfs.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 245 "zhprfs.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 246 "zhprfs.f"
	*info = -1;
#line 247 "zhprfs.f"
    } else if (*n < 0) {
#line 248 "zhprfs.f"
	*info = -2;
#line 249 "zhprfs.f"
    } else if (*nrhs < 0) {
#line 250 "zhprfs.f"
	*info = -3;
#line 251 "zhprfs.f"
    } else if (*ldb < max(1,*n)) {
#line 252 "zhprfs.f"
	*info = -8;
#line 253 "zhprfs.f"
    } else if (*ldx < max(1,*n)) {
#line 254 "zhprfs.f"
	*info = -10;
#line 255 "zhprfs.f"
    }
#line 256 "zhprfs.f"
    if (*info != 0) {
#line 257 "zhprfs.f"
	i__1 = -(*info);
#line 257 "zhprfs.f"
	xerbla_("ZHPRFS", &i__1, (ftnlen)6);
#line 258 "zhprfs.f"
	return 0;
#line 259 "zhprfs.f"
    }

/*     Quick return if possible */

#line 263 "zhprfs.f"
    if (*n == 0 || *nrhs == 0) {
#line 264 "zhprfs.f"
	i__1 = *nrhs;
#line 264 "zhprfs.f"
	for (j = 1; j <= i__1; ++j) {
#line 265 "zhprfs.f"
	    ferr[j] = 0.;
#line 266 "zhprfs.f"
	    berr[j] = 0.;
#line 267 "zhprfs.f"
/* L10: */
#line 267 "zhprfs.f"
	}
#line 268 "zhprfs.f"
	return 0;
#line 269 "zhprfs.f"
    }

/*     NZ = maximum number of nonzero elements in each row of A, plus 1 */

#line 273 "zhprfs.f"
    nz = *n + 1;
#line 274 "zhprfs.f"
    eps = dlamch_("Epsilon", (ftnlen)7);
#line 275 "zhprfs.f"
    safmin = dlamch_("Safe minimum", (ftnlen)12);
#line 276 "zhprfs.f"
    safe1 = nz * safmin;
#line 277 "zhprfs.f"
    safe2 = safe1 / eps;

/*     Do for each right hand side */

#line 281 "zhprfs.f"
    i__1 = *nrhs;
#line 281 "zhprfs.f"
    for (j = 1; j <= i__1; ++j) {

#line 283 "zhprfs.f"
	count = 1;
#line 284 "zhprfs.f"
	lstres = 3.;
#line 285 "zhprfs.f"
L20:

/*        Loop until stopping criterion is satisfied. */

/*        Compute residual R = B - A * X */

#line 291 "zhprfs.f"
	zcopy_(n, &b[j * b_dim1 + 1], &c__1, &work[1], &c__1);
#line 292 "zhprfs.f"
	z__1.r = -1., z__1.i = -0.;
#line 292 "zhprfs.f"
	zhpmv_(uplo, n, &z__1, &ap[1], &x[j * x_dim1 + 1], &c__1, &c_b1, &
		work[1], &c__1, (ftnlen)1);

/*        Compute componentwise relative backward error from formula */

/*        max(i) ( abs(R(i)) / ( abs(A)*abs(X) + abs(B) )(i) ) */

/*        where abs(Z) is the componentwise absolute value of the matrix */
/*        or vector Z.  If the i-th component of the denominator is less */
/*        than SAFE2, then SAFE1 is added to the i-th components of the */
/*        numerator and denominator before dividing. */

#line 303 "zhprfs.f"
	i__2 = *n;
#line 303 "zhprfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 304 "zhprfs.f"
	    i__3 = i__ + j * b_dim1;
#line 304 "zhprfs.f"
	    rwork[i__] = (d__1 = b[i__3].r, abs(d__1)) + (d__2 = d_imag(&b[
		    i__ + j * b_dim1]), abs(d__2));
#line 305 "zhprfs.f"
/* L30: */
#line 305 "zhprfs.f"
	}

/*        Compute abs(A)*abs(X) + abs(B). */

#line 309 "zhprfs.f"
	kk = 1;
#line 310 "zhprfs.f"
	if (upper) {
#line 311 "zhprfs.f"
	    i__2 = *n;
#line 311 "zhprfs.f"
	    for (k = 1; k <= i__2; ++k) {
#line 312 "zhprfs.f"
		s = 0.;
#line 313 "zhprfs.f"
		i__3 = k + j * x_dim1;
#line 313 "zhprfs.f"
		xk = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = d_imag(&x[k + j *
			 x_dim1]), abs(d__2));
#line 314 "zhprfs.f"
		ik = kk;
#line 315 "zhprfs.f"
		i__3 = k - 1;
#line 315 "zhprfs.f"
		for (i__ = 1; i__ <= i__3; ++i__) {
#line 316 "zhprfs.f"
		    i__4 = ik;
#line 316 "zhprfs.f"
		    rwork[i__] += ((d__1 = ap[i__4].r, abs(d__1)) + (d__2 = 
			    d_imag(&ap[ik]), abs(d__2))) * xk;
#line 317 "zhprfs.f"
		    i__4 = ik;
#line 317 "zhprfs.f"
		    i__5 = i__ + j * x_dim1;
#line 317 "zhprfs.f"
		    s += ((d__1 = ap[i__4].r, abs(d__1)) + (d__2 = d_imag(&ap[
			    ik]), abs(d__2))) * ((d__3 = x[i__5].r, abs(d__3))
			     + (d__4 = d_imag(&x[i__ + j * x_dim1]), abs(d__4)
			    ));
#line 318 "zhprfs.f"
		    ++ik;
#line 319 "zhprfs.f"
/* L40: */
#line 319 "zhprfs.f"
		}
#line 320 "zhprfs.f"
		i__3 = kk + k - 1;
#line 320 "zhprfs.f"
		rwork[k] = rwork[k] + (d__1 = ap[i__3].r, abs(d__1)) * xk + s;
#line 322 "zhprfs.f"
		kk += k;
#line 323 "zhprfs.f"
/* L50: */
#line 323 "zhprfs.f"
	    }
#line 324 "zhprfs.f"
	} else {
#line 325 "zhprfs.f"
	    i__2 = *n;
#line 325 "zhprfs.f"
	    for (k = 1; k <= i__2; ++k) {
#line 326 "zhprfs.f"
		s = 0.;
#line 327 "zhprfs.f"
		i__3 = k + j * x_dim1;
#line 327 "zhprfs.f"
		xk = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = d_imag(&x[k + j *
			 x_dim1]), abs(d__2));
#line 328 "zhprfs.f"
		i__3 = kk;
#line 328 "zhprfs.f"
		rwork[k] += (d__1 = ap[i__3].r, abs(d__1)) * xk;
#line 329 "zhprfs.f"
		ik = kk + 1;
#line 330 "zhprfs.f"
		i__3 = *n;
#line 330 "zhprfs.f"
		for (i__ = k + 1; i__ <= i__3; ++i__) {
#line 331 "zhprfs.f"
		    i__4 = ik;
#line 331 "zhprfs.f"
		    rwork[i__] += ((d__1 = ap[i__4].r, abs(d__1)) + (d__2 = 
			    d_imag(&ap[ik]), abs(d__2))) * xk;
#line 332 "zhprfs.f"
		    i__4 = ik;
#line 332 "zhprfs.f"
		    i__5 = i__ + j * x_dim1;
#line 332 "zhprfs.f"
		    s += ((d__1 = ap[i__4].r, abs(d__1)) + (d__2 = d_imag(&ap[
			    ik]), abs(d__2))) * ((d__3 = x[i__5].r, abs(d__3))
			     + (d__4 = d_imag(&x[i__ + j * x_dim1]), abs(d__4)
			    ));
#line 333 "zhprfs.f"
		    ++ik;
#line 334 "zhprfs.f"
/* L60: */
#line 334 "zhprfs.f"
		}
#line 335 "zhprfs.f"
		rwork[k] += s;
#line 336 "zhprfs.f"
		kk += *n - k + 1;
#line 337 "zhprfs.f"
/* L70: */
#line 337 "zhprfs.f"
	    }
#line 338 "zhprfs.f"
	}
#line 339 "zhprfs.f"
	s = 0.;
#line 340 "zhprfs.f"
	i__2 = *n;
#line 340 "zhprfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 341 "zhprfs.f"
	    if (rwork[i__] > safe2) {
/* Computing MAX */
#line 342 "zhprfs.f"
		i__3 = i__;
#line 342 "zhprfs.f"
		d__3 = s, d__4 = ((d__1 = work[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&work[i__]), abs(d__2))) / rwork[i__];
#line 342 "zhprfs.f"
		s = max(d__3,d__4);
#line 343 "zhprfs.f"
	    } else {
/* Computing MAX */
#line 344 "zhprfs.f"
		i__3 = i__;
#line 344 "zhprfs.f"
		d__3 = s, d__4 = ((d__1 = work[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&work[i__]), abs(d__2)) + safe1) / (rwork[i__] 
			+ safe1);
#line 344 "zhprfs.f"
		s = max(d__3,d__4);
#line 346 "zhprfs.f"
	    }
#line 347 "zhprfs.f"
/* L80: */
#line 347 "zhprfs.f"
	}
#line 348 "zhprfs.f"
	berr[j] = s;

/*        Test stopping criterion. Continue iterating if */
/*           1) The residual BERR(J) is larger than machine epsilon, and */
/*           2) BERR(J) decreased by at least a factor of 2 during the */
/*              last iteration, and */
/*           3) At most ITMAX iterations tried. */

#line 356 "zhprfs.f"
	if (berr[j] > eps && berr[j] * 2. <= lstres && count <= 5) {

/*           Update solution and try again. */

#line 361 "zhprfs.f"
	    zhptrs_(uplo, n, &c__1, &afp[1], &ipiv[1], &work[1], n, info, (
		    ftnlen)1);
#line 362 "zhprfs.f"
	    zaxpy_(n, &c_b1, &work[1], &c__1, &x[j * x_dim1 + 1], &c__1);
#line 363 "zhprfs.f"
	    lstres = berr[j];
#line 364 "zhprfs.f"
	    ++count;
#line 365 "zhprfs.f"
	    goto L20;
#line 366 "zhprfs.f"
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

/*        Use ZLACN2 to estimate the infinity-norm of the matrix */
/*           inv(A) * diag(W), */
/*        where W = abs(R) + NZ*EPS*( abs(A)*abs(X)+abs(B) ))) */

#line 390 "zhprfs.f"
	i__2 = *n;
#line 390 "zhprfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 391 "zhprfs.f"
	    if (rwork[i__] > safe2) {
#line 392 "zhprfs.f"
		i__3 = i__;
#line 392 "zhprfs.f"
		rwork[i__] = (d__1 = work[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&work[i__]), abs(d__2)) + nz * eps * rwork[i__]
			;
#line 393 "zhprfs.f"
	    } else {
#line 394 "zhprfs.f"
		i__3 = i__;
#line 394 "zhprfs.f"
		rwork[i__] = (d__1 = work[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&work[i__]), abs(d__2)) + nz * eps * rwork[i__]
			 + safe1;
#line 396 "zhprfs.f"
	    }
#line 397 "zhprfs.f"
/* L90: */
#line 397 "zhprfs.f"
	}

#line 399 "zhprfs.f"
	kase = 0;
#line 400 "zhprfs.f"
L100:
#line 401 "zhprfs.f"
	zlacn2_(n, &work[*n + 1], &work[1], &ferr[j], &kase, isave);
#line 402 "zhprfs.f"
	if (kase != 0) {
#line 403 "zhprfs.f"
	    if (kase == 1) {

/*              Multiply by diag(W)*inv(A**H). */

#line 407 "zhprfs.f"
		zhptrs_(uplo, n, &c__1, &afp[1], &ipiv[1], &work[1], n, info, 
			(ftnlen)1);
#line 408 "zhprfs.f"
		i__2 = *n;
#line 408 "zhprfs.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 409 "zhprfs.f"
		    i__3 = i__;
#line 409 "zhprfs.f"
		    i__4 = i__;
#line 409 "zhprfs.f"
		    i__5 = i__;
#line 409 "zhprfs.f"
		    z__1.r = rwork[i__4] * work[i__5].r, z__1.i = rwork[i__4] 
			    * work[i__5].i;
#line 409 "zhprfs.f"
		    work[i__3].r = z__1.r, work[i__3].i = z__1.i;
#line 410 "zhprfs.f"
/* L110: */
#line 410 "zhprfs.f"
		}
#line 411 "zhprfs.f"
	    } else if (kase == 2) {

/*              Multiply by inv(A)*diag(W). */

#line 415 "zhprfs.f"
		i__2 = *n;
#line 415 "zhprfs.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 416 "zhprfs.f"
		    i__3 = i__;
#line 416 "zhprfs.f"
		    i__4 = i__;
#line 416 "zhprfs.f"
		    i__5 = i__;
#line 416 "zhprfs.f"
		    z__1.r = rwork[i__4] * work[i__5].r, z__1.i = rwork[i__4] 
			    * work[i__5].i;
#line 416 "zhprfs.f"
		    work[i__3].r = z__1.r, work[i__3].i = z__1.i;
#line 417 "zhprfs.f"
/* L120: */
#line 417 "zhprfs.f"
		}
#line 418 "zhprfs.f"
		zhptrs_(uplo, n, &c__1, &afp[1], &ipiv[1], &work[1], n, info, 
			(ftnlen)1);
#line 419 "zhprfs.f"
	    }
#line 420 "zhprfs.f"
	    goto L100;
#line 421 "zhprfs.f"
	}

/*        Normalize error. */

#line 425 "zhprfs.f"
	lstres = 0.;
#line 426 "zhprfs.f"
	i__2 = *n;
#line 426 "zhprfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 427 "zhprfs.f"
	    i__3 = i__ + j * x_dim1;
#line 427 "zhprfs.f"
	    d__3 = lstres, d__4 = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = 
		    d_imag(&x[i__ + j * x_dim1]), abs(d__2));
#line 427 "zhprfs.f"
	    lstres = max(d__3,d__4);
#line 428 "zhprfs.f"
/* L130: */
#line 428 "zhprfs.f"
	}
#line 429 "zhprfs.f"
	if (lstres != 0.) {
#line 429 "zhprfs.f"
	    ferr[j] /= lstres;
#line 429 "zhprfs.f"
	}

#line 432 "zhprfs.f"
/* L140: */
#line 432 "zhprfs.f"
    }

#line 434 "zhprfs.f"
    return 0;

/*     End of ZHPRFS */

} /* zhprfs_ */

