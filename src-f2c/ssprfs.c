#line 1 "ssprfs.f"
/* ssprfs.f -- translated by f2c (version 20100827).
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

#line 1 "ssprfs.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b12 = -1.;
static doublereal c_b14 = 1.;

/* > \brief \b SSPRFS */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SSPRFS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssprfs.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssprfs.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssprfs.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SSPRFS( UPLO, N, NRHS, AP, AFP, IPIV, B, LDB, X, LDX, */
/*                          FERR, BERR, WORK, IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDB, LDX, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ), IWORK( * ) */
/*       REAL               AFP( * ), AP( * ), B( LDB, * ), BERR( * ), */
/*      $                   FERR( * ), WORK( * ), X( LDX, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SSPRFS improves the computed solution to a system of linear */
/* > equations when the coefficient matrix is symmetric indefinite */
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
/* >          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n. */
/* > \endverbatim */
/* > */
/* > \param[in] AFP */
/* > \verbatim */
/* >          AFP is REAL array, dimension (N*(N+1)/2) */
/* >          The factored form of the matrix A.  AFP contains the block */
/* >          diagonal matrix D and the multipliers used to obtain the */
/* >          factor U or L from the factorization A = U*D*U**T or */
/* >          A = L*D*L**T as computed by SSPTRF, stored as a packed */
/* >          triangular matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N) */
/* >          Details of the interchanges and the block structure of D */
/* >          as determined by SSPTRF. */
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
/* >          On entry, the solution matrix X, as computed by SSPTRS. */
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

/* > \ingroup realOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int ssprfs_(char *uplo, integer *n, integer *nrhs, 
	doublereal *ap, doublereal *afp, integer *ipiv, doublereal *b, 
	integer *ldb, doublereal *x, integer *ldx, doublereal *ferr, 
	doublereal *berr, doublereal *work, integer *iwork, integer *info, 
	ftnlen uplo_len)
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
    extern /* Subroutine */ int ssptrs_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen);


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

#line 234 "ssprfs.f"
    /* Parameter adjustments */
#line 234 "ssprfs.f"
    --ap;
#line 234 "ssprfs.f"
    --afp;
#line 234 "ssprfs.f"
    --ipiv;
#line 234 "ssprfs.f"
    b_dim1 = *ldb;
#line 234 "ssprfs.f"
    b_offset = 1 + b_dim1;
#line 234 "ssprfs.f"
    b -= b_offset;
#line 234 "ssprfs.f"
    x_dim1 = *ldx;
#line 234 "ssprfs.f"
    x_offset = 1 + x_dim1;
#line 234 "ssprfs.f"
    x -= x_offset;
#line 234 "ssprfs.f"
    --ferr;
#line 234 "ssprfs.f"
    --berr;
#line 234 "ssprfs.f"
    --work;
#line 234 "ssprfs.f"
    --iwork;
#line 234 "ssprfs.f"

#line 234 "ssprfs.f"
    /* Function Body */
#line 234 "ssprfs.f"
    *info = 0;
#line 235 "ssprfs.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 236 "ssprfs.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 237 "ssprfs.f"
	*info = -1;
#line 238 "ssprfs.f"
    } else if (*n < 0) {
#line 239 "ssprfs.f"
	*info = -2;
#line 240 "ssprfs.f"
    } else if (*nrhs < 0) {
#line 241 "ssprfs.f"
	*info = -3;
#line 242 "ssprfs.f"
    } else if (*ldb < max(1,*n)) {
#line 243 "ssprfs.f"
	*info = -8;
#line 244 "ssprfs.f"
    } else if (*ldx < max(1,*n)) {
#line 245 "ssprfs.f"
	*info = -10;
#line 246 "ssprfs.f"
    }
#line 247 "ssprfs.f"
    if (*info != 0) {
#line 248 "ssprfs.f"
	i__1 = -(*info);
#line 248 "ssprfs.f"
	xerbla_("SSPRFS", &i__1, (ftnlen)6);
#line 249 "ssprfs.f"
	return 0;
#line 250 "ssprfs.f"
    }

/*     Quick return if possible */

#line 254 "ssprfs.f"
    if (*n == 0 || *nrhs == 0) {
#line 255 "ssprfs.f"
	i__1 = *nrhs;
#line 255 "ssprfs.f"
	for (j = 1; j <= i__1; ++j) {
#line 256 "ssprfs.f"
	    ferr[j] = 0.;
#line 257 "ssprfs.f"
	    berr[j] = 0.;
#line 258 "ssprfs.f"
/* L10: */
#line 258 "ssprfs.f"
	}
#line 259 "ssprfs.f"
	return 0;
#line 260 "ssprfs.f"
    }

/*     NZ = maximum number of nonzero elements in each row of A, plus 1 */

#line 264 "ssprfs.f"
    nz = *n + 1;
#line 265 "ssprfs.f"
    eps = slamch_("Epsilon", (ftnlen)7);
#line 266 "ssprfs.f"
    safmin = slamch_("Safe minimum", (ftnlen)12);
#line 267 "ssprfs.f"
    safe1 = nz * safmin;
#line 268 "ssprfs.f"
    safe2 = safe1 / eps;

/*     Do for each right hand side */

#line 272 "ssprfs.f"
    i__1 = *nrhs;
#line 272 "ssprfs.f"
    for (j = 1; j <= i__1; ++j) {

#line 274 "ssprfs.f"
	count = 1;
#line 275 "ssprfs.f"
	lstres = 3.;
#line 276 "ssprfs.f"
L20:

/*        Loop until stopping criterion is satisfied. */

/*        Compute residual R = B - A * X */

#line 282 "ssprfs.f"
	scopy_(n, &b[j * b_dim1 + 1], &c__1, &work[*n + 1], &c__1);
#line 283 "ssprfs.f"
	sspmv_(uplo, n, &c_b12, &ap[1], &x[j * x_dim1 + 1], &c__1, &c_b14, &
		work[*n + 1], &c__1, (ftnlen)1);

/*        Compute componentwise relative backward error from formula */

/*        max(i) ( abs(R(i)) / ( abs(A)*abs(X) + abs(B) )(i) ) */

/*        where abs(Z) is the componentwise absolute value of the matrix */
/*        or vector Z.  If the i-th component of the denominator is less */
/*        than SAFE2, then SAFE1 is added to the i-th components of the */
/*        numerator and denominator before dividing. */

#line 295 "ssprfs.f"
	i__2 = *n;
#line 295 "ssprfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 296 "ssprfs.f"
	    work[i__] = (d__1 = b[i__ + j * b_dim1], abs(d__1));
#line 297 "ssprfs.f"
/* L30: */
#line 297 "ssprfs.f"
	}

/*        Compute abs(A)*abs(X) + abs(B). */

#line 301 "ssprfs.f"
	kk = 1;
#line 302 "ssprfs.f"
	if (upper) {
#line 303 "ssprfs.f"
	    i__2 = *n;
#line 303 "ssprfs.f"
	    for (k = 1; k <= i__2; ++k) {
#line 304 "ssprfs.f"
		s = 0.;
#line 305 "ssprfs.f"
		xk = (d__1 = x[k + j * x_dim1], abs(d__1));
#line 306 "ssprfs.f"
		ik = kk;
#line 307 "ssprfs.f"
		i__3 = k - 1;
#line 307 "ssprfs.f"
		for (i__ = 1; i__ <= i__3; ++i__) {
#line 308 "ssprfs.f"
		    work[i__] += (d__1 = ap[ik], abs(d__1)) * xk;
#line 309 "ssprfs.f"
		    s += (d__1 = ap[ik], abs(d__1)) * (d__2 = x[i__ + j * 
			    x_dim1], abs(d__2));
#line 310 "ssprfs.f"
		    ++ik;
#line 311 "ssprfs.f"
/* L40: */
#line 311 "ssprfs.f"
		}
#line 312 "ssprfs.f"
		work[k] = work[k] + (d__1 = ap[kk + k - 1], abs(d__1)) * xk + 
			s;
#line 313 "ssprfs.f"
		kk += k;
#line 314 "ssprfs.f"
/* L50: */
#line 314 "ssprfs.f"
	    }
#line 315 "ssprfs.f"
	} else {
#line 316 "ssprfs.f"
	    i__2 = *n;
#line 316 "ssprfs.f"
	    for (k = 1; k <= i__2; ++k) {
#line 317 "ssprfs.f"
		s = 0.;
#line 318 "ssprfs.f"
		xk = (d__1 = x[k + j * x_dim1], abs(d__1));
#line 319 "ssprfs.f"
		work[k] += (d__1 = ap[kk], abs(d__1)) * xk;
#line 320 "ssprfs.f"
		ik = kk + 1;
#line 321 "ssprfs.f"
		i__3 = *n;
#line 321 "ssprfs.f"
		for (i__ = k + 1; i__ <= i__3; ++i__) {
#line 322 "ssprfs.f"
		    work[i__] += (d__1 = ap[ik], abs(d__1)) * xk;
#line 323 "ssprfs.f"
		    s += (d__1 = ap[ik], abs(d__1)) * (d__2 = x[i__ + j * 
			    x_dim1], abs(d__2));
#line 324 "ssprfs.f"
		    ++ik;
#line 325 "ssprfs.f"
/* L60: */
#line 325 "ssprfs.f"
		}
#line 326 "ssprfs.f"
		work[k] += s;
#line 327 "ssprfs.f"
		kk += *n - k + 1;
#line 328 "ssprfs.f"
/* L70: */
#line 328 "ssprfs.f"
	    }
#line 329 "ssprfs.f"
	}
#line 330 "ssprfs.f"
	s = 0.;
#line 331 "ssprfs.f"
	i__2 = *n;
#line 331 "ssprfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 332 "ssprfs.f"
	    if (work[i__] > safe2) {
/* Computing MAX */
#line 333 "ssprfs.f"
		d__2 = s, d__3 = (d__1 = work[*n + i__], abs(d__1)) / work[
			i__];
#line 333 "ssprfs.f"
		s = max(d__2,d__3);
#line 334 "ssprfs.f"
	    } else {
/* Computing MAX */
#line 335 "ssprfs.f"
		d__2 = s, d__3 = ((d__1 = work[*n + i__], abs(d__1)) + safe1) 
			/ (work[i__] + safe1);
#line 335 "ssprfs.f"
		s = max(d__2,d__3);
#line 337 "ssprfs.f"
	    }
#line 338 "ssprfs.f"
/* L80: */
#line 338 "ssprfs.f"
	}
#line 339 "ssprfs.f"
	berr[j] = s;

/*        Test stopping criterion. Continue iterating if */
/*           1) The residual BERR(J) is larger than machine epsilon, and */
/*           2) BERR(J) decreased by at least a factor of 2 during the */
/*              last iteration, and */
/*           3) At most ITMAX iterations tried. */

#line 347 "ssprfs.f"
	if (berr[j] > eps && berr[j] * 2. <= lstres && count <= 5) {

/*           Update solution and try again. */

#line 352 "ssprfs.f"
	    ssptrs_(uplo, n, &c__1, &afp[1], &ipiv[1], &work[*n + 1], n, info,
		     (ftnlen)1);
#line 353 "ssprfs.f"
	    saxpy_(n, &c_b14, &work[*n + 1], &c__1, &x[j * x_dim1 + 1], &c__1)
		    ;
#line 354 "ssprfs.f"
	    lstres = berr[j];
#line 355 "ssprfs.f"
	    ++count;
#line 356 "ssprfs.f"
	    goto L20;
#line 357 "ssprfs.f"
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

#line 381 "ssprfs.f"
	i__2 = *n;
#line 381 "ssprfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 382 "ssprfs.f"
	    if (work[i__] > safe2) {
#line 383 "ssprfs.f"
		work[i__] = (d__1 = work[*n + i__], abs(d__1)) + nz * eps * 
			work[i__];
#line 384 "ssprfs.f"
	    } else {
#line 385 "ssprfs.f"
		work[i__] = (d__1 = work[*n + i__], abs(d__1)) + nz * eps * 
			work[i__] + safe1;
#line 386 "ssprfs.f"
	    }
#line 387 "ssprfs.f"
/* L90: */
#line 387 "ssprfs.f"
	}

#line 389 "ssprfs.f"
	kase = 0;
#line 390 "ssprfs.f"
L100:
#line 391 "ssprfs.f"
	slacn2_(n, &work[(*n << 1) + 1], &work[*n + 1], &iwork[1], &ferr[j], &
		kase, isave);
#line 393 "ssprfs.f"
	if (kase != 0) {
#line 394 "ssprfs.f"
	    if (kase == 1) {

/*              Multiply by diag(W)*inv(A**T). */

#line 398 "ssprfs.f"
		ssptrs_(uplo, n, &c__1, &afp[1], &ipiv[1], &work[*n + 1], n, 
			info, (ftnlen)1);
#line 400 "ssprfs.f"
		i__2 = *n;
#line 400 "ssprfs.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 401 "ssprfs.f"
		    work[*n + i__] = work[i__] * work[*n + i__];
#line 402 "ssprfs.f"
/* L110: */
#line 402 "ssprfs.f"
		}
#line 403 "ssprfs.f"
	    } else if (kase == 2) {

/*              Multiply by inv(A)*diag(W). */

#line 407 "ssprfs.f"
		i__2 = *n;
#line 407 "ssprfs.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 408 "ssprfs.f"
		    work[*n + i__] = work[i__] * work[*n + i__];
#line 409 "ssprfs.f"
/* L120: */
#line 409 "ssprfs.f"
		}
#line 410 "ssprfs.f"
		ssptrs_(uplo, n, &c__1, &afp[1], &ipiv[1], &work[*n + 1], n, 
			info, (ftnlen)1);
#line 412 "ssprfs.f"
	    }
#line 413 "ssprfs.f"
	    goto L100;
#line 414 "ssprfs.f"
	}

/*        Normalize error. */

#line 418 "ssprfs.f"
	lstres = 0.;
#line 419 "ssprfs.f"
	i__2 = *n;
#line 419 "ssprfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 420 "ssprfs.f"
	    d__2 = lstres, d__3 = (d__1 = x[i__ + j * x_dim1], abs(d__1));
#line 420 "ssprfs.f"
	    lstres = max(d__2,d__3);
#line 421 "ssprfs.f"
/* L130: */
#line 421 "ssprfs.f"
	}
#line 422 "ssprfs.f"
	if (lstres != 0.) {
#line 422 "ssprfs.f"
	    ferr[j] /= lstres;
#line 422 "ssprfs.f"
	}

#line 425 "ssprfs.f"
/* L140: */
#line 425 "ssprfs.f"
    }

#line 427 "ssprfs.f"
    return 0;

/*     End of SSPRFS */

} /* ssprfs_ */

