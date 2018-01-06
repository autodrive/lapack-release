#line 1 "dsprfs.f"
/* dsprfs.f -- translated by f2c (version 20100827).
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

#line 1 "dsprfs.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b12 = -1.;
static doublereal c_b14 = 1.;

/* > \brief \b DSPRFS */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DSPRFS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsprfs.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsprfs.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsprfs.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSPRFS( UPLO, N, NRHS, AP, AFP, IPIV, B, LDB, X, LDX, */
/*                          FERR, BERR, WORK, IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDB, LDX, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ), IWORK( * ) */
/*       DOUBLE PRECISION   AFP( * ), AP( * ), B( LDB, * ), BERR( * ), */
/*      $                   FERR( * ), WORK( * ), X( LDX, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DSPRFS improves the computed solution to a system of linear */
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
/* >          AP is DOUBLE PRECISION array, dimension (N*(N+1)/2) */
/* >          The upper or lower triangle of the symmetric matrix A, packed */
/* >          columnwise in a linear array.  The j-th column of A is stored */
/* >          in the array AP as follows: */
/* >          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j; */
/* >          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n. */
/* > \endverbatim */
/* > */
/* > \param[in] AFP */
/* > \verbatim */
/* >          AFP is DOUBLE PRECISION array, dimension (N*(N+1)/2) */
/* >          The factored form of the matrix A.  AFP contains the block */
/* >          diagonal matrix D and the multipliers used to obtain the */
/* >          factor U or L from the factorization A = U*D*U**T or */
/* >          A = L*D*L**T as computed by DSPTRF, stored as a packed */
/* >          triangular matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N) */
/* >          Details of the interchanges and the block structure of D */
/* >          as determined by DSPTRF. */
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
/* >          On entry, the solution matrix X, as computed by DSPTRS. */
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

/* > \ingroup doubleOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int dsprfs_(char *uplo, integer *n, integer *nrhs, 
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
    extern /* Subroutine */ int dsptrs_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
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
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 234 "dsprfs.f"
    /* Parameter adjustments */
#line 234 "dsprfs.f"
    --ap;
#line 234 "dsprfs.f"
    --afp;
#line 234 "dsprfs.f"
    --ipiv;
#line 234 "dsprfs.f"
    b_dim1 = *ldb;
#line 234 "dsprfs.f"
    b_offset = 1 + b_dim1;
#line 234 "dsprfs.f"
    b -= b_offset;
#line 234 "dsprfs.f"
    x_dim1 = *ldx;
#line 234 "dsprfs.f"
    x_offset = 1 + x_dim1;
#line 234 "dsprfs.f"
    x -= x_offset;
#line 234 "dsprfs.f"
    --ferr;
#line 234 "dsprfs.f"
    --berr;
#line 234 "dsprfs.f"
    --work;
#line 234 "dsprfs.f"
    --iwork;
#line 234 "dsprfs.f"

#line 234 "dsprfs.f"
    /* Function Body */
#line 234 "dsprfs.f"
    *info = 0;
#line 235 "dsprfs.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 236 "dsprfs.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 237 "dsprfs.f"
	*info = -1;
#line 238 "dsprfs.f"
    } else if (*n < 0) {
#line 239 "dsprfs.f"
	*info = -2;
#line 240 "dsprfs.f"
    } else if (*nrhs < 0) {
#line 241 "dsprfs.f"
	*info = -3;
#line 242 "dsprfs.f"
    } else if (*ldb < max(1,*n)) {
#line 243 "dsprfs.f"
	*info = -8;
#line 244 "dsprfs.f"
    } else if (*ldx < max(1,*n)) {
#line 245 "dsprfs.f"
	*info = -10;
#line 246 "dsprfs.f"
    }
#line 247 "dsprfs.f"
    if (*info != 0) {
#line 248 "dsprfs.f"
	i__1 = -(*info);
#line 248 "dsprfs.f"
	xerbla_("DSPRFS", &i__1, (ftnlen)6);
#line 249 "dsprfs.f"
	return 0;
#line 250 "dsprfs.f"
    }

/*     Quick return if possible */

#line 254 "dsprfs.f"
    if (*n == 0 || *nrhs == 0) {
#line 255 "dsprfs.f"
	i__1 = *nrhs;
#line 255 "dsprfs.f"
	for (j = 1; j <= i__1; ++j) {
#line 256 "dsprfs.f"
	    ferr[j] = 0.;
#line 257 "dsprfs.f"
	    berr[j] = 0.;
#line 258 "dsprfs.f"
/* L10: */
#line 258 "dsprfs.f"
	}
#line 259 "dsprfs.f"
	return 0;
#line 260 "dsprfs.f"
    }

/*     NZ = maximum number of nonzero elements in each row of A, plus 1 */

#line 264 "dsprfs.f"
    nz = *n + 1;
#line 265 "dsprfs.f"
    eps = dlamch_("Epsilon", (ftnlen)7);
#line 266 "dsprfs.f"
    safmin = dlamch_("Safe minimum", (ftnlen)12);
#line 267 "dsprfs.f"
    safe1 = nz * safmin;
#line 268 "dsprfs.f"
    safe2 = safe1 / eps;

/*     Do for each right hand side */

#line 272 "dsprfs.f"
    i__1 = *nrhs;
#line 272 "dsprfs.f"
    for (j = 1; j <= i__1; ++j) {

#line 274 "dsprfs.f"
	count = 1;
#line 275 "dsprfs.f"
	lstres = 3.;
#line 276 "dsprfs.f"
L20:

/*        Loop until stopping criterion is satisfied. */

/*        Compute residual R = B - A * X */

#line 282 "dsprfs.f"
	dcopy_(n, &b[j * b_dim1 + 1], &c__1, &work[*n + 1], &c__1);
#line 283 "dsprfs.f"
	dspmv_(uplo, n, &c_b12, &ap[1], &x[j * x_dim1 + 1], &c__1, &c_b14, &
		work[*n + 1], &c__1, (ftnlen)1);

/*        Compute componentwise relative backward error from formula */

/*        max(i) ( abs(R(i)) / ( abs(A)*abs(X) + abs(B) )(i) ) */

/*        where abs(Z) is the componentwise absolute value of the matrix */
/*        or vector Z.  If the i-th component of the denominator is less */
/*        than SAFE2, then SAFE1 is added to the i-th components of the */
/*        numerator and denominator before dividing. */

#line 295 "dsprfs.f"
	i__2 = *n;
#line 295 "dsprfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 296 "dsprfs.f"
	    work[i__] = (d__1 = b[i__ + j * b_dim1], abs(d__1));
#line 297 "dsprfs.f"
/* L30: */
#line 297 "dsprfs.f"
	}

/*        Compute abs(A)*abs(X) + abs(B). */

#line 301 "dsprfs.f"
	kk = 1;
#line 302 "dsprfs.f"
	if (upper) {
#line 303 "dsprfs.f"
	    i__2 = *n;
#line 303 "dsprfs.f"
	    for (k = 1; k <= i__2; ++k) {
#line 304 "dsprfs.f"
		s = 0.;
#line 305 "dsprfs.f"
		xk = (d__1 = x[k + j * x_dim1], abs(d__1));
#line 306 "dsprfs.f"
		ik = kk;
#line 307 "dsprfs.f"
		i__3 = k - 1;
#line 307 "dsprfs.f"
		for (i__ = 1; i__ <= i__3; ++i__) {
#line 308 "dsprfs.f"
		    work[i__] += (d__1 = ap[ik], abs(d__1)) * xk;
#line 309 "dsprfs.f"
		    s += (d__1 = ap[ik], abs(d__1)) * (d__2 = x[i__ + j * 
			    x_dim1], abs(d__2));
#line 310 "dsprfs.f"
		    ++ik;
#line 311 "dsprfs.f"
/* L40: */
#line 311 "dsprfs.f"
		}
#line 312 "dsprfs.f"
		work[k] = work[k] + (d__1 = ap[kk + k - 1], abs(d__1)) * xk + 
			s;
#line 313 "dsprfs.f"
		kk += k;
#line 314 "dsprfs.f"
/* L50: */
#line 314 "dsprfs.f"
	    }
#line 315 "dsprfs.f"
	} else {
#line 316 "dsprfs.f"
	    i__2 = *n;
#line 316 "dsprfs.f"
	    for (k = 1; k <= i__2; ++k) {
#line 317 "dsprfs.f"
		s = 0.;
#line 318 "dsprfs.f"
		xk = (d__1 = x[k + j * x_dim1], abs(d__1));
#line 319 "dsprfs.f"
		work[k] += (d__1 = ap[kk], abs(d__1)) * xk;
#line 320 "dsprfs.f"
		ik = kk + 1;
#line 321 "dsprfs.f"
		i__3 = *n;
#line 321 "dsprfs.f"
		for (i__ = k + 1; i__ <= i__3; ++i__) {
#line 322 "dsprfs.f"
		    work[i__] += (d__1 = ap[ik], abs(d__1)) * xk;
#line 323 "dsprfs.f"
		    s += (d__1 = ap[ik], abs(d__1)) * (d__2 = x[i__ + j * 
			    x_dim1], abs(d__2));
#line 324 "dsprfs.f"
		    ++ik;
#line 325 "dsprfs.f"
/* L60: */
#line 325 "dsprfs.f"
		}
#line 326 "dsprfs.f"
		work[k] += s;
#line 327 "dsprfs.f"
		kk += *n - k + 1;
#line 328 "dsprfs.f"
/* L70: */
#line 328 "dsprfs.f"
	    }
#line 329 "dsprfs.f"
	}
#line 330 "dsprfs.f"
	s = 0.;
#line 331 "dsprfs.f"
	i__2 = *n;
#line 331 "dsprfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 332 "dsprfs.f"
	    if (work[i__] > safe2) {
/* Computing MAX */
#line 333 "dsprfs.f"
		d__2 = s, d__3 = (d__1 = work[*n + i__], abs(d__1)) / work[
			i__];
#line 333 "dsprfs.f"
		s = max(d__2,d__3);
#line 334 "dsprfs.f"
	    } else {
/* Computing MAX */
#line 335 "dsprfs.f"
		d__2 = s, d__3 = ((d__1 = work[*n + i__], abs(d__1)) + safe1) 
			/ (work[i__] + safe1);
#line 335 "dsprfs.f"
		s = max(d__2,d__3);
#line 337 "dsprfs.f"
	    }
#line 338 "dsprfs.f"
/* L80: */
#line 338 "dsprfs.f"
	}
#line 339 "dsprfs.f"
	berr[j] = s;

/*        Test stopping criterion. Continue iterating if */
/*           1) The residual BERR(J) is larger than machine epsilon, and */
/*           2) BERR(J) decreased by at least a factor of 2 during the */
/*              last iteration, and */
/*           3) At most ITMAX iterations tried. */

#line 347 "dsprfs.f"
	if (berr[j] > eps && berr[j] * 2. <= lstres && count <= 5) {

/*           Update solution and try again. */

#line 352 "dsprfs.f"
	    dsptrs_(uplo, n, &c__1, &afp[1], &ipiv[1], &work[*n + 1], n, info,
		     (ftnlen)1);
#line 353 "dsprfs.f"
	    daxpy_(n, &c_b14, &work[*n + 1], &c__1, &x[j * x_dim1 + 1], &c__1)
		    ;
#line 354 "dsprfs.f"
	    lstres = berr[j];
#line 355 "dsprfs.f"
	    ++count;
#line 356 "dsprfs.f"
	    goto L20;
#line 357 "dsprfs.f"
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

#line 381 "dsprfs.f"
	i__2 = *n;
#line 381 "dsprfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 382 "dsprfs.f"
	    if (work[i__] > safe2) {
#line 383 "dsprfs.f"
		work[i__] = (d__1 = work[*n + i__], abs(d__1)) + nz * eps * 
			work[i__];
#line 384 "dsprfs.f"
	    } else {
#line 385 "dsprfs.f"
		work[i__] = (d__1 = work[*n + i__], abs(d__1)) + nz * eps * 
			work[i__] + safe1;
#line 386 "dsprfs.f"
	    }
#line 387 "dsprfs.f"
/* L90: */
#line 387 "dsprfs.f"
	}

#line 389 "dsprfs.f"
	kase = 0;
#line 390 "dsprfs.f"
L100:
#line 391 "dsprfs.f"
	dlacn2_(n, &work[(*n << 1) + 1], &work[*n + 1], &iwork[1], &ferr[j], &
		kase, isave);
#line 393 "dsprfs.f"
	if (kase != 0) {
#line 394 "dsprfs.f"
	    if (kase == 1) {

/*              Multiply by diag(W)*inv(A**T). */

#line 398 "dsprfs.f"
		dsptrs_(uplo, n, &c__1, &afp[1], &ipiv[1], &work[*n + 1], n, 
			info, (ftnlen)1);
#line 400 "dsprfs.f"
		i__2 = *n;
#line 400 "dsprfs.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 401 "dsprfs.f"
		    work[*n + i__] = work[i__] * work[*n + i__];
#line 402 "dsprfs.f"
/* L110: */
#line 402 "dsprfs.f"
		}
#line 403 "dsprfs.f"
	    } else if (kase == 2) {

/*              Multiply by inv(A)*diag(W). */

#line 407 "dsprfs.f"
		i__2 = *n;
#line 407 "dsprfs.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 408 "dsprfs.f"
		    work[*n + i__] = work[i__] * work[*n + i__];
#line 409 "dsprfs.f"
/* L120: */
#line 409 "dsprfs.f"
		}
#line 410 "dsprfs.f"
		dsptrs_(uplo, n, &c__1, &afp[1], &ipiv[1], &work[*n + 1], n, 
			info, (ftnlen)1);
#line 412 "dsprfs.f"
	    }
#line 413 "dsprfs.f"
	    goto L100;
#line 414 "dsprfs.f"
	}

/*        Normalize error. */

#line 418 "dsprfs.f"
	lstres = 0.;
#line 419 "dsprfs.f"
	i__2 = *n;
#line 419 "dsprfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 420 "dsprfs.f"
	    d__2 = lstres, d__3 = (d__1 = x[i__ + j * x_dim1], abs(d__1));
#line 420 "dsprfs.f"
	    lstres = max(d__2,d__3);
#line 421 "dsprfs.f"
/* L130: */
#line 421 "dsprfs.f"
	}
#line 422 "dsprfs.f"
	if (lstres != 0.) {
#line 422 "dsprfs.f"
	    ferr[j] /= lstres;
#line 422 "dsprfs.f"
	}

#line 425 "dsprfs.f"
/* L140: */
#line 425 "dsprfs.f"
    }

#line 427 "dsprfs.f"
    return 0;

/*     End of DSPRFS */

} /* dsprfs_ */

