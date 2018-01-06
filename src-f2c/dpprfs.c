#line 1 "dpprfs.f"
/* dpprfs.f -- translated by f2c (version 20100827).
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

#line 1 "dpprfs.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b12 = -1.;
static doublereal c_b14 = 1.;

/* > \brief \b DPPRFS */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DPPRFS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dpprfs.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dpprfs.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dpprfs.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DPPRFS( UPLO, N, NRHS, AP, AFP, B, LDB, X, LDX, FERR, */
/*                          BERR, WORK, IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDB, LDX, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       DOUBLE PRECISION   AFP( * ), AP( * ), B( LDB, * ), BERR( * ), */
/*      $                   FERR( * ), WORK( * ), X( LDX, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DPPRFS improves the computed solution to a system of linear */
/* > equations when the coefficient matrix is symmetric positive definite */
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
/* >          AP is DOUBLE PRECISION array, dimension (N*(N+1)/2) */
/* >          The upper or lower triangle of the symmetric matrix A, packed */
/* >          columnwise in a linear array.  The j-th column of A is stored */
/* >          in the array AP as follows: */
/* >          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j; */
/* >          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n. */
/* > \endverbatim */
/* > */
/* > \param[in] AFP */
/* > \verbatim */
/* >          AFP is DOUBLE PRECISION array, dimension (N*(N+1)/2) */
/* >          The triangular factor U or L from the Cholesky factorization */
/* >          A = U**T*U or A = L*L**T, as computed by DPPTRF/ZPPTRF, */
/* >          packed columnwise in a linear array in the same format as A */
/* >          (see AP). */
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
/* >          On entry, the solution matrix X, as computed by DPPTRS. */
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

/* > \ingroup doubleOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int dpprfs_(char *uplo, integer *n, integer *nrhs, 
	doublereal *ap, doublereal *afp, doublereal *b, integer *ldb, 
	doublereal *x, integer *ldx, doublereal *ferr, doublereal *berr, 
	doublereal *work, integer *iwork, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer b_dim1, b_offset, x_dim1, x_offset, i__1, i__2, i__3;
    doublereal d__1, d__2, d__3;

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
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), daxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *);
    static integer count;
    extern /* Subroutine */ int dspmv_(char *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     integer *, ftnlen);
    static logical upper;
    extern /* Subroutine */ int dlacn2_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, integer *);
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal lstres;
    extern /* Subroutine */ int dpptrs_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, ftnlen);


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

#line 226 "dpprfs.f"
    /* Parameter adjustments */
#line 226 "dpprfs.f"
    --ap;
#line 226 "dpprfs.f"
    --afp;
#line 226 "dpprfs.f"
    b_dim1 = *ldb;
#line 226 "dpprfs.f"
    b_offset = 1 + b_dim1;
#line 226 "dpprfs.f"
    b -= b_offset;
#line 226 "dpprfs.f"
    x_dim1 = *ldx;
#line 226 "dpprfs.f"
    x_offset = 1 + x_dim1;
#line 226 "dpprfs.f"
    x -= x_offset;
#line 226 "dpprfs.f"
    --ferr;
#line 226 "dpprfs.f"
    --berr;
#line 226 "dpprfs.f"
    --work;
#line 226 "dpprfs.f"
    --iwork;
#line 226 "dpprfs.f"

#line 226 "dpprfs.f"
    /* Function Body */
#line 226 "dpprfs.f"
    *info = 0;
#line 227 "dpprfs.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 228 "dpprfs.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 229 "dpprfs.f"
	*info = -1;
#line 230 "dpprfs.f"
    } else if (*n < 0) {
#line 231 "dpprfs.f"
	*info = -2;
#line 232 "dpprfs.f"
    } else if (*nrhs < 0) {
#line 233 "dpprfs.f"
	*info = -3;
#line 234 "dpprfs.f"
    } else if (*ldb < max(1,*n)) {
#line 235 "dpprfs.f"
	*info = -7;
#line 236 "dpprfs.f"
    } else if (*ldx < max(1,*n)) {
#line 237 "dpprfs.f"
	*info = -9;
#line 238 "dpprfs.f"
    }
#line 239 "dpprfs.f"
    if (*info != 0) {
#line 240 "dpprfs.f"
	i__1 = -(*info);
#line 240 "dpprfs.f"
	xerbla_("DPPRFS", &i__1, (ftnlen)6);
#line 241 "dpprfs.f"
	return 0;
#line 242 "dpprfs.f"
    }

/*     Quick return if possible */

#line 246 "dpprfs.f"
    if (*n == 0 || *nrhs == 0) {
#line 247 "dpprfs.f"
	i__1 = *nrhs;
#line 247 "dpprfs.f"
	for (j = 1; j <= i__1; ++j) {
#line 248 "dpprfs.f"
	    ferr[j] = 0.;
#line 249 "dpprfs.f"
	    berr[j] = 0.;
#line 250 "dpprfs.f"
/* L10: */
#line 250 "dpprfs.f"
	}
#line 251 "dpprfs.f"
	return 0;
#line 252 "dpprfs.f"
    }

/*     NZ = maximum number of nonzero elements in each row of A, plus 1 */

#line 256 "dpprfs.f"
    nz = *n + 1;
#line 257 "dpprfs.f"
    eps = dlamch_("Epsilon", (ftnlen)7);
#line 258 "dpprfs.f"
    safmin = dlamch_("Safe minimum", (ftnlen)12);
#line 259 "dpprfs.f"
    safe1 = nz * safmin;
#line 260 "dpprfs.f"
    safe2 = safe1 / eps;

/*     Do for each right hand side */

#line 264 "dpprfs.f"
    i__1 = *nrhs;
#line 264 "dpprfs.f"
    for (j = 1; j <= i__1; ++j) {

#line 266 "dpprfs.f"
	count = 1;
#line 267 "dpprfs.f"
	lstres = 3.;
#line 268 "dpprfs.f"
L20:

/*        Loop until stopping criterion is satisfied. */

/*        Compute residual R = B - A * X */

#line 274 "dpprfs.f"
	dcopy_(n, &b[j * b_dim1 + 1], &c__1, &work[*n + 1], &c__1);
#line 275 "dpprfs.f"
	dspmv_(uplo, n, &c_b12, &ap[1], &x[j * x_dim1 + 1], &c__1, &c_b14, &
		work[*n + 1], &c__1, (ftnlen)1);

/*        Compute componentwise relative backward error from formula */

/*        max(i) ( abs(R(i)) / ( abs(A)*abs(X) + abs(B) )(i) ) */

/*        where abs(Z) is the componentwise absolute value of the matrix */
/*        or vector Z.  If the i-th component of the denominator is less */
/*        than SAFE2, then SAFE1 is added to the i-th components of the */
/*        numerator and denominator before dividing. */

#line 287 "dpprfs.f"
	i__2 = *n;
#line 287 "dpprfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 288 "dpprfs.f"
	    work[i__] = (d__1 = b[i__ + j * b_dim1], abs(d__1));
#line 289 "dpprfs.f"
/* L30: */
#line 289 "dpprfs.f"
	}

/*        Compute abs(A)*abs(X) + abs(B). */

#line 293 "dpprfs.f"
	kk = 1;
#line 294 "dpprfs.f"
	if (upper) {
#line 295 "dpprfs.f"
	    i__2 = *n;
#line 295 "dpprfs.f"
	    for (k = 1; k <= i__2; ++k) {
#line 296 "dpprfs.f"
		s = 0.;
#line 297 "dpprfs.f"
		xk = (d__1 = x[k + j * x_dim1], abs(d__1));
#line 298 "dpprfs.f"
		ik = kk;
#line 299 "dpprfs.f"
		i__3 = k - 1;
#line 299 "dpprfs.f"
		for (i__ = 1; i__ <= i__3; ++i__) {
#line 300 "dpprfs.f"
		    work[i__] += (d__1 = ap[ik], abs(d__1)) * xk;
#line 301 "dpprfs.f"
		    s += (d__1 = ap[ik], abs(d__1)) * (d__2 = x[i__ + j * 
			    x_dim1], abs(d__2));
#line 302 "dpprfs.f"
		    ++ik;
#line 303 "dpprfs.f"
/* L40: */
#line 303 "dpprfs.f"
		}
#line 304 "dpprfs.f"
		work[k] = work[k] + (d__1 = ap[kk + k - 1], abs(d__1)) * xk + 
			s;
#line 305 "dpprfs.f"
		kk += k;
#line 306 "dpprfs.f"
/* L50: */
#line 306 "dpprfs.f"
	    }
#line 307 "dpprfs.f"
	} else {
#line 308 "dpprfs.f"
	    i__2 = *n;
#line 308 "dpprfs.f"
	    for (k = 1; k <= i__2; ++k) {
#line 309 "dpprfs.f"
		s = 0.;
#line 310 "dpprfs.f"
		xk = (d__1 = x[k + j * x_dim1], abs(d__1));
#line 311 "dpprfs.f"
		work[k] += (d__1 = ap[kk], abs(d__1)) * xk;
#line 312 "dpprfs.f"
		ik = kk + 1;
#line 313 "dpprfs.f"
		i__3 = *n;
#line 313 "dpprfs.f"
		for (i__ = k + 1; i__ <= i__3; ++i__) {
#line 314 "dpprfs.f"
		    work[i__] += (d__1 = ap[ik], abs(d__1)) * xk;
#line 315 "dpprfs.f"
		    s += (d__1 = ap[ik], abs(d__1)) * (d__2 = x[i__ + j * 
			    x_dim1], abs(d__2));
#line 316 "dpprfs.f"
		    ++ik;
#line 317 "dpprfs.f"
/* L60: */
#line 317 "dpprfs.f"
		}
#line 318 "dpprfs.f"
		work[k] += s;
#line 319 "dpprfs.f"
		kk += *n - k + 1;
#line 320 "dpprfs.f"
/* L70: */
#line 320 "dpprfs.f"
	    }
#line 321 "dpprfs.f"
	}
#line 322 "dpprfs.f"
	s = 0.;
#line 323 "dpprfs.f"
	i__2 = *n;
#line 323 "dpprfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 324 "dpprfs.f"
	    if (work[i__] > safe2) {
/* Computing MAX */
#line 325 "dpprfs.f"
		d__2 = s, d__3 = (d__1 = work[*n + i__], abs(d__1)) / work[
			i__];
#line 325 "dpprfs.f"
		s = max(d__2,d__3);
#line 326 "dpprfs.f"
	    } else {
/* Computing MAX */
#line 327 "dpprfs.f"
		d__2 = s, d__3 = ((d__1 = work[*n + i__], abs(d__1)) + safe1) 
			/ (work[i__] + safe1);
#line 327 "dpprfs.f"
		s = max(d__2,d__3);
#line 329 "dpprfs.f"
	    }
#line 330 "dpprfs.f"
/* L80: */
#line 330 "dpprfs.f"
	}
#line 331 "dpprfs.f"
	berr[j] = s;

/*        Test stopping criterion. Continue iterating if */
/*           1) The residual BERR(J) is larger than machine epsilon, and */
/*           2) BERR(J) decreased by at least a factor of 2 during the */
/*              last iteration, and */
/*           3) At most ITMAX iterations tried. */

#line 339 "dpprfs.f"
	if (berr[j] > eps && berr[j] * 2. <= lstres && count <= 5) {

/*           Update solution and try again. */

#line 344 "dpprfs.f"
	    dpptrs_(uplo, n, &c__1, &afp[1], &work[*n + 1], n, info, (ftnlen)
		    1);
#line 345 "dpprfs.f"
	    daxpy_(n, &c_b14, &work[*n + 1], &c__1, &x[j * x_dim1 + 1], &c__1)
		    ;
#line 346 "dpprfs.f"
	    lstres = berr[j];
#line 347 "dpprfs.f"
	    ++count;
#line 348 "dpprfs.f"
	    goto L20;
#line 349 "dpprfs.f"
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

/*        Use DLACN2 to estimate the infinity-norm of the matrix */
/*           inv(A) * diag(W), */
/*        where W = abs(R) + NZ*EPS*( abs(A)*abs(X)+abs(B) ))) */

#line 373 "dpprfs.f"
	i__2 = *n;
#line 373 "dpprfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 374 "dpprfs.f"
	    if (work[i__] > safe2) {
#line 375 "dpprfs.f"
		work[i__] = (d__1 = work[*n + i__], abs(d__1)) + nz * eps * 
			work[i__];
#line 376 "dpprfs.f"
	    } else {
#line 377 "dpprfs.f"
		work[i__] = (d__1 = work[*n + i__], abs(d__1)) + nz * eps * 
			work[i__] + safe1;
#line 378 "dpprfs.f"
	    }
#line 379 "dpprfs.f"
/* L90: */
#line 379 "dpprfs.f"
	}

#line 381 "dpprfs.f"
	kase = 0;
#line 382 "dpprfs.f"
L100:
#line 383 "dpprfs.f"
	dlacn2_(n, &work[(*n << 1) + 1], &work[*n + 1], &iwork[1], &ferr[j], &
		kase, isave);
#line 385 "dpprfs.f"
	if (kase != 0) {
#line 386 "dpprfs.f"
	    if (kase == 1) {

/*              Multiply by diag(W)*inv(A**T). */

#line 390 "dpprfs.f"
		dpptrs_(uplo, n, &c__1, &afp[1], &work[*n + 1], n, info, (
			ftnlen)1);
#line 391 "dpprfs.f"
		i__2 = *n;
#line 391 "dpprfs.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 392 "dpprfs.f"
		    work[*n + i__] = work[i__] * work[*n + i__];
#line 393 "dpprfs.f"
/* L110: */
#line 393 "dpprfs.f"
		}
#line 394 "dpprfs.f"
	    } else if (kase == 2) {

/*              Multiply by inv(A)*diag(W). */

#line 398 "dpprfs.f"
		i__2 = *n;
#line 398 "dpprfs.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 399 "dpprfs.f"
		    work[*n + i__] = work[i__] * work[*n + i__];
#line 400 "dpprfs.f"
/* L120: */
#line 400 "dpprfs.f"
		}
#line 401 "dpprfs.f"
		dpptrs_(uplo, n, &c__1, &afp[1], &work[*n + 1], n, info, (
			ftnlen)1);
#line 402 "dpprfs.f"
	    }
#line 403 "dpprfs.f"
	    goto L100;
#line 404 "dpprfs.f"
	}

/*        Normalize error. */

#line 408 "dpprfs.f"
	lstres = 0.;
#line 409 "dpprfs.f"
	i__2 = *n;
#line 409 "dpprfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 410 "dpprfs.f"
	    d__2 = lstres, d__3 = (d__1 = x[i__ + j * x_dim1], abs(d__1));
#line 410 "dpprfs.f"
	    lstres = max(d__2,d__3);
#line 411 "dpprfs.f"
/* L130: */
#line 411 "dpprfs.f"
	}
#line 412 "dpprfs.f"
	if (lstres != 0.) {
#line 412 "dpprfs.f"
	    ferr[j] /= lstres;
#line 412 "dpprfs.f"
	}

#line 415 "dpprfs.f"
/* L140: */
#line 415 "dpprfs.f"
    }

#line 417 "dpprfs.f"
    return 0;

/*     End of DPPRFS */

} /* dpprfs_ */

