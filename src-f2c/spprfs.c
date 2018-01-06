#line 1 "spprfs.f"
/* spprfs.f -- translated by f2c (version 20100827).
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

#line 1 "spprfs.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b12 = -1.;
static doublereal c_b14 = 1.;

/* > \brief \b SPPRFS */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SPPRFS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/spprfs.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/spprfs.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/spprfs.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SPPRFS( UPLO, N, NRHS, AP, AFP, B, LDB, X, LDX, FERR, */
/*                          BERR, WORK, IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDB, LDX, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       REAL               AFP( * ), AP( * ), B( LDB, * ), BERR( * ), */
/*      $                   FERR( * ), WORK( * ), X( LDX, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SPPRFS improves the computed solution to a system of linear */
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
/* >          AP is REAL array, dimension (N*(N+1)/2) */
/* >          The upper or lower triangle of the symmetric matrix A, packed */
/* >          columnwise in a linear array.  The j-th column of A is stored */
/* >          in the array AP as follows: */
/* >          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j; */
/* >          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n. */
/* > \endverbatim */
/* > */
/* > \param[in] AFP */
/* > \verbatim */
/* >          AFP is REAL array, dimension (N*(N+1)/2) */
/* >          The triangular factor U or L from the Cholesky factorization */
/* >          A = U**T*U or A = L*L**T, as computed by SPPTRF/CPPTRF, */
/* >          packed columnwise in a linear array in the same format as A */
/* >          (see AP). */
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
/* >          On entry, the solution matrix X, as computed by SPPTRS. */
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

/* > \date November 2011 */

/* > \ingroup realOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int spprfs_(char *uplo, integer *n, integer *nrhs, 
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
    static integer isave[3], count;
    static logical upper;
    extern /* Subroutine */ int scopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), saxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *), sspmv_(char *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen), slacn2_(integer *,
	     doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    integer *);
    extern doublereal slamch_(char *, ftnlen);
    static doublereal safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal lstres;
    extern /* Subroutine */ int spptrs_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, ftnlen);


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

#line 226 "spprfs.f"
    /* Parameter adjustments */
#line 226 "spprfs.f"
    --ap;
#line 226 "spprfs.f"
    --afp;
#line 226 "spprfs.f"
    b_dim1 = *ldb;
#line 226 "spprfs.f"
    b_offset = 1 + b_dim1;
#line 226 "spprfs.f"
    b -= b_offset;
#line 226 "spprfs.f"
    x_dim1 = *ldx;
#line 226 "spprfs.f"
    x_offset = 1 + x_dim1;
#line 226 "spprfs.f"
    x -= x_offset;
#line 226 "spprfs.f"
    --ferr;
#line 226 "spprfs.f"
    --berr;
#line 226 "spprfs.f"
    --work;
#line 226 "spprfs.f"
    --iwork;
#line 226 "spprfs.f"

#line 226 "spprfs.f"
    /* Function Body */
#line 226 "spprfs.f"
    *info = 0;
#line 227 "spprfs.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 228 "spprfs.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 229 "spprfs.f"
	*info = -1;
#line 230 "spprfs.f"
    } else if (*n < 0) {
#line 231 "spprfs.f"
	*info = -2;
#line 232 "spprfs.f"
    } else if (*nrhs < 0) {
#line 233 "spprfs.f"
	*info = -3;
#line 234 "spprfs.f"
    } else if (*ldb < max(1,*n)) {
#line 235 "spprfs.f"
	*info = -7;
#line 236 "spprfs.f"
    } else if (*ldx < max(1,*n)) {
#line 237 "spprfs.f"
	*info = -9;
#line 238 "spprfs.f"
    }
#line 239 "spprfs.f"
    if (*info != 0) {
#line 240 "spprfs.f"
	i__1 = -(*info);
#line 240 "spprfs.f"
	xerbla_("SPPRFS", &i__1, (ftnlen)6);
#line 241 "spprfs.f"
	return 0;
#line 242 "spprfs.f"
    }

/*     Quick return if possible */

#line 246 "spprfs.f"
    if (*n == 0 || *nrhs == 0) {
#line 247 "spprfs.f"
	i__1 = *nrhs;
#line 247 "spprfs.f"
	for (j = 1; j <= i__1; ++j) {
#line 248 "spprfs.f"
	    ferr[j] = 0.;
#line 249 "spprfs.f"
	    berr[j] = 0.;
#line 250 "spprfs.f"
/* L10: */
#line 250 "spprfs.f"
	}
#line 251 "spprfs.f"
	return 0;
#line 252 "spprfs.f"
    }

/*     NZ = maximum number of nonzero elements in each row of A, plus 1 */

#line 256 "spprfs.f"
    nz = *n + 1;
#line 257 "spprfs.f"
    eps = slamch_("Epsilon", (ftnlen)7);
#line 258 "spprfs.f"
    safmin = slamch_("Safe minimum", (ftnlen)12);
#line 259 "spprfs.f"
    safe1 = nz * safmin;
#line 260 "spprfs.f"
    safe2 = safe1 / eps;

/*     Do for each right hand side */

#line 264 "spprfs.f"
    i__1 = *nrhs;
#line 264 "spprfs.f"
    for (j = 1; j <= i__1; ++j) {

#line 266 "spprfs.f"
	count = 1;
#line 267 "spprfs.f"
	lstres = 3.;
#line 268 "spprfs.f"
L20:

/*        Loop until stopping criterion is satisfied. */

/*        Compute residual R = B - A * X */

#line 274 "spprfs.f"
	scopy_(n, &b[j * b_dim1 + 1], &c__1, &work[*n + 1], &c__1);
#line 275 "spprfs.f"
	sspmv_(uplo, n, &c_b12, &ap[1], &x[j * x_dim1 + 1], &c__1, &c_b14, &
		work[*n + 1], &c__1, (ftnlen)1);

/*        Compute componentwise relative backward error from formula */

/*        max(i) ( abs(R(i)) / ( abs(A)*abs(X) + abs(B) )(i) ) */

/*        where abs(Z) is the componentwise absolute value of the matrix */
/*        or vector Z.  If the i-th component of the denominator is less */
/*        than SAFE2, then SAFE1 is added to the i-th components of the */
/*        numerator and denominator before dividing. */

#line 287 "spprfs.f"
	i__2 = *n;
#line 287 "spprfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 288 "spprfs.f"
	    work[i__] = (d__1 = b[i__ + j * b_dim1], abs(d__1));
#line 289 "spprfs.f"
/* L30: */
#line 289 "spprfs.f"
	}

/*        Compute abs(A)*abs(X) + abs(B). */

#line 293 "spprfs.f"
	kk = 1;
#line 294 "spprfs.f"
	if (upper) {
#line 295 "spprfs.f"
	    i__2 = *n;
#line 295 "spprfs.f"
	    for (k = 1; k <= i__2; ++k) {
#line 296 "spprfs.f"
		s = 0.;
#line 297 "spprfs.f"
		xk = (d__1 = x[k + j * x_dim1], abs(d__1));
#line 298 "spprfs.f"
		ik = kk;
#line 299 "spprfs.f"
		i__3 = k - 1;
#line 299 "spprfs.f"
		for (i__ = 1; i__ <= i__3; ++i__) {
#line 300 "spprfs.f"
		    work[i__] += (d__1 = ap[ik], abs(d__1)) * xk;
#line 301 "spprfs.f"
		    s += (d__1 = ap[ik], abs(d__1)) * (d__2 = x[i__ + j * 
			    x_dim1], abs(d__2));
#line 302 "spprfs.f"
		    ++ik;
#line 303 "spprfs.f"
/* L40: */
#line 303 "spprfs.f"
		}
#line 304 "spprfs.f"
		work[k] = work[k] + (d__1 = ap[kk + k - 1], abs(d__1)) * xk + 
			s;
#line 305 "spprfs.f"
		kk += k;
#line 306 "spprfs.f"
/* L50: */
#line 306 "spprfs.f"
	    }
#line 307 "spprfs.f"
	} else {
#line 308 "spprfs.f"
	    i__2 = *n;
#line 308 "spprfs.f"
	    for (k = 1; k <= i__2; ++k) {
#line 309 "spprfs.f"
		s = 0.;
#line 310 "spprfs.f"
		xk = (d__1 = x[k + j * x_dim1], abs(d__1));
#line 311 "spprfs.f"
		work[k] += (d__1 = ap[kk], abs(d__1)) * xk;
#line 312 "spprfs.f"
		ik = kk + 1;
#line 313 "spprfs.f"
		i__3 = *n;
#line 313 "spprfs.f"
		for (i__ = k + 1; i__ <= i__3; ++i__) {
#line 314 "spprfs.f"
		    work[i__] += (d__1 = ap[ik], abs(d__1)) * xk;
#line 315 "spprfs.f"
		    s += (d__1 = ap[ik], abs(d__1)) * (d__2 = x[i__ + j * 
			    x_dim1], abs(d__2));
#line 316 "spprfs.f"
		    ++ik;
#line 317 "spprfs.f"
/* L60: */
#line 317 "spprfs.f"
		}
#line 318 "spprfs.f"
		work[k] += s;
#line 319 "spprfs.f"
		kk += *n - k + 1;
#line 320 "spprfs.f"
/* L70: */
#line 320 "spprfs.f"
	    }
#line 321 "spprfs.f"
	}
#line 322 "spprfs.f"
	s = 0.;
#line 323 "spprfs.f"
	i__2 = *n;
#line 323 "spprfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 324 "spprfs.f"
	    if (work[i__] > safe2) {
/* Computing MAX */
#line 325 "spprfs.f"
		d__2 = s, d__3 = (d__1 = work[*n + i__], abs(d__1)) / work[
			i__];
#line 325 "spprfs.f"
		s = max(d__2,d__3);
#line 326 "spprfs.f"
	    } else {
/* Computing MAX */
#line 327 "spprfs.f"
		d__2 = s, d__3 = ((d__1 = work[*n + i__], abs(d__1)) + safe1) 
			/ (work[i__] + safe1);
#line 327 "spprfs.f"
		s = max(d__2,d__3);
#line 329 "spprfs.f"
	    }
#line 330 "spprfs.f"
/* L80: */
#line 330 "spprfs.f"
	}
#line 331 "spprfs.f"
	berr[j] = s;

/*        Test stopping criterion. Continue iterating if */
/*           1) The residual BERR(J) is larger than machine epsilon, and */
/*           2) BERR(J) decreased by at least a factor of 2 during the */
/*              last iteration, and */
/*           3) At most ITMAX iterations tried. */

#line 339 "spprfs.f"
	if (berr[j] > eps && berr[j] * 2. <= lstres && count <= 5) {

/*           Update solution and try again. */

#line 344 "spprfs.f"
	    spptrs_(uplo, n, &c__1, &afp[1], &work[*n + 1], n, info, (ftnlen)
		    1);
#line 345 "spprfs.f"
	    saxpy_(n, &c_b14, &work[*n + 1], &c__1, &x[j * x_dim1 + 1], &c__1)
		    ;
#line 346 "spprfs.f"
	    lstres = berr[j];
#line 347 "spprfs.f"
	    ++count;
#line 348 "spprfs.f"
	    goto L20;
#line 349 "spprfs.f"
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

/*        Use SLACN2 to estimate the infinity-norm of the matrix */
/*           inv(A) * diag(W), */
/*        where W = abs(R) + NZ*EPS*( abs(A)*abs(X)+abs(B) ))) */

#line 373 "spprfs.f"
	i__2 = *n;
#line 373 "spprfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 374 "spprfs.f"
	    if (work[i__] > safe2) {
#line 375 "spprfs.f"
		work[i__] = (d__1 = work[*n + i__], abs(d__1)) + nz * eps * 
			work[i__];
#line 376 "spprfs.f"
	    } else {
#line 377 "spprfs.f"
		work[i__] = (d__1 = work[*n + i__], abs(d__1)) + nz * eps * 
			work[i__] + safe1;
#line 378 "spprfs.f"
	    }
#line 379 "spprfs.f"
/* L90: */
#line 379 "spprfs.f"
	}

#line 381 "spprfs.f"
	kase = 0;
#line 382 "spprfs.f"
L100:
#line 383 "spprfs.f"
	slacn2_(n, &work[(*n << 1) + 1], &work[*n + 1], &iwork[1], &ferr[j], &
		kase, isave);
#line 385 "spprfs.f"
	if (kase != 0) {
#line 386 "spprfs.f"
	    if (kase == 1) {

/*              Multiply by diag(W)*inv(A**T). */

#line 390 "spprfs.f"
		spptrs_(uplo, n, &c__1, &afp[1], &work[*n + 1], n, info, (
			ftnlen)1);
#line 391 "spprfs.f"
		i__2 = *n;
#line 391 "spprfs.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 392 "spprfs.f"
		    work[*n + i__] = work[i__] * work[*n + i__];
#line 393 "spprfs.f"
/* L110: */
#line 393 "spprfs.f"
		}
#line 394 "spprfs.f"
	    } else if (kase == 2) {

/*              Multiply by inv(A)*diag(W). */

#line 398 "spprfs.f"
		i__2 = *n;
#line 398 "spprfs.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 399 "spprfs.f"
		    work[*n + i__] = work[i__] * work[*n + i__];
#line 400 "spprfs.f"
/* L120: */
#line 400 "spprfs.f"
		}
#line 401 "spprfs.f"
		spptrs_(uplo, n, &c__1, &afp[1], &work[*n + 1], n, info, (
			ftnlen)1);
#line 402 "spprfs.f"
	    }
#line 403 "spprfs.f"
	    goto L100;
#line 404 "spprfs.f"
	}

/*        Normalize error. */

#line 408 "spprfs.f"
	lstres = 0.;
#line 409 "spprfs.f"
	i__2 = *n;
#line 409 "spprfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 410 "spprfs.f"
	    d__2 = lstres, d__3 = (d__1 = x[i__ + j * x_dim1], abs(d__1));
#line 410 "spprfs.f"
	    lstres = max(d__2,d__3);
#line 411 "spprfs.f"
/* L130: */
#line 411 "spprfs.f"
	}
#line 412 "spprfs.f"
	if (lstres != 0.) {
#line 412 "spprfs.f"
	    ferr[j] /= lstres;
#line 412 "spprfs.f"
	}

#line 415 "spprfs.f"
/* L140: */
#line 415 "spprfs.f"
    }

#line 417 "spprfs.f"
    return 0;

/*     End of SPPRFS */

} /* spprfs_ */

