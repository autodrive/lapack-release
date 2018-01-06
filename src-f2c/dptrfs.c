#line 1 "dptrfs.f"
/* dptrfs.f -- translated by f2c (version 20100827).
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

#line 1 "dptrfs.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b11 = 1.;

/* > \brief \b DPTRFS */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DPTRFS + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dptrfs.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dptrfs.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dptrfs.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DPTRFS( N, NRHS, D, E, DF, EF, B, LDB, X, LDX, FERR, */
/*                          BERR, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDB, LDX, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   B( LDB, * ), BERR( * ), D( * ), DF( * ), */
/*      $                   E( * ), EF( * ), FERR( * ), WORK( * ), */
/*      $                   X( LDX, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DPTRFS improves the computed solution to a system of linear */
/* > equations when the coefficient matrix is symmetric positive definite */
/* > and tridiagonal, and provides error bounds and backward error */
/* > estimates for the solution. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

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
/* >          of the matrix B.  NRHS >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* >          D is DOUBLE PRECISION array, dimension (N) */
/* >          The n diagonal elements of the tridiagonal matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] E */
/* > \verbatim */
/* >          E is DOUBLE PRECISION array, dimension (N-1) */
/* >          The (n-1) subdiagonal elements of the tridiagonal matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] DF */
/* > \verbatim */
/* >          DF is DOUBLE PRECISION array, dimension (N) */
/* >          The n diagonal elements of the diagonal matrix D from the */
/* >          factorization computed by DPTTRF. */
/* > \endverbatim */
/* > */
/* > \param[in] EF */
/* > \verbatim */
/* >          EF is DOUBLE PRECISION array, dimension (N-1) */
/* >          The (n-1) subdiagonal elements of the unit bidiagonal factor */
/* >          L from the factorization computed by DPTTRF. */
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
/* >          On entry, the solution matrix X, as computed by DPTTRS. */
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
/* >          The forward error bound for each solution vector */
/* >          X(j) (the j-th column of the solution matrix X). */
/* >          If XTRUE is the true solution corresponding to X(j), FERR(j) */
/* >          is an estimated upper bound for the magnitude of the largest */
/* >          element in (X(j) - XTRUE) divided by the magnitude of the */
/* >          largest element in X(j). */
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
/* >          WORK is DOUBLE PRECISION array, dimension (2*N) */
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

/* > \ingroup doublePTcomputational */

/*  ===================================================================== */
/* Subroutine */ int dptrfs_(integer *n, integer *nrhs, doublereal *d__, 
	doublereal *e, doublereal *df, doublereal *ef, doublereal *b, integer 
	*ldb, doublereal *x, integer *ldx, doublereal *ferr, doublereal *berr,
	 doublereal *work, integer *info)
{
    /* System generated locals */
    integer b_dim1, b_offset, x_dim1, x_offset, i__1, i__2;
    doublereal d__1, d__2, d__3;

    /* Local variables */
    static integer i__, j;
    static doublereal s, bi, cx, dx, ex;
    static integer ix, nz;
    static doublereal eps, safe1, safe2;
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    static integer count;
    extern doublereal dlamch_(char *, ftnlen);
    extern integer idamax_(integer *, doublereal *, integer *);
    static doublereal safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal lstres;
    extern /* Subroutine */ int dpttrs_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *);


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
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 214 "dptrfs.f"
    /* Parameter adjustments */
#line 214 "dptrfs.f"
    --d__;
#line 214 "dptrfs.f"
    --e;
#line 214 "dptrfs.f"
    --df;
#line 214 "dptrfs.f"
    --ef;
#line 214 "dptrfs.f"
    b_dim1 = *ldb;
#line 214 "dptrfs.f"
    b_offset = 1 + b_dim1;
#line 214 "dptrfs.f"
    b -= b_offset;
#line 214 "dptrfs.f"
    x_dim1 = *ldx;
#line 214 "dptrfs.f"
    x_offset = 1 + x_dim1;
#line 214 "dptrfs.f"
    x -= x_offset;
#line 214 "dptrfs.f"
    --ferr;
#line 214 "dptrfs.f"
    --berr;
#line 214 "dptrfs.f"
    --work;
#line 214 "dptrfs.f"

#line 214 "dptrfs.f"
    /* Function Body */
#line 214 "dptrfs.f"
    *info = 0;
#line 215 "dptrfs.f"
    if (*n < 0) {
#line 216 "dptrfs.f"
	*info = -1;
#line 217 "dptrfs.f"
    } else if (*nrhs < 0) {
#line 218 "dptrfs.f"
	*info = -2;
#line 219 "dptrfs.f"
    } else if (*ldb < max(1,*n)) {
#line 220 "dptrfs.f"
	*info = -8;
#line 221 "dptrfs.f"
    } else if (*ldx < max(1,*n)) {
#line 222 "dptrfs.f"
	*info = -10;
#line 223 "dptrfs.f"
    }
#line 224 "dptrfs.f"
    if (*info != 0) {
#line 225 "dptrfs.f"
	i__1 = -(*info);
#line 225 "dptrfs.f"
	xerbla_("DPTRFS", &i__1, (ftnlen)6);
#line 226 "dptrfs.f"
	return 0;
#line 227 "dptrfs.f"
    }

/*     Quick return if possible */

#line 231 "dptrfs.f"
    if (*n == 0 || *nrhs == 0) {
#line 232 "dptrfs.f"
	i__1 = *nrhs;
#line 232 "dptrfs.f"
	for (j = 1; j <= i__1; ++j) {
#line 233 "dptrfs.f"
	    ferr[j] = 0.;
#line 234 "dptrfs.f"
	    berr[j] = 0.;
#line 235 "dptrfs.f"
/* L10: */
#line 235 "dptrfs.f"
	}
#line 236 "dptrfs.f"
	return 0;
#line 237 "dptrfs.f"
    }

/*     NZ = maximum number of nonzero elements in each row of A, plus 1 */

#line 241 "dptrfs.f"
    nz = 4;
#line 242 "dptrfs.f"
    eps = dlamch_("Epsilon", (ftnlen)7);
#line 243 "dptrfs.f"
    safmin = dlamch_("Safe minimum", (ftnlen)12);
#line 244 "dptrfs.f"
    safe1 = nz * safmin;
#line 245 "dptrfs.f"
    safe2 = safe1 / eps;

/*     Do for each right hand side */

#line 249 "dptrfs.f"
    i__1 = *nrhs;
#line 249 "dptrfs.f"
    for (j = 1; j <= i__1; ++j) {

#line 251 "dptrfs.f"
	count = 1;
#line 252 "dptrfs.f"
	lstres = 3.;
#line 253 "dptrfs.f"
L20:

/*        Loop until stopping criterion is satisfied. */

/*        Compute residual R = B - A * X.  Also compute */
/*        abs(A)*abs(x) + abs(b) for use in the backward error bound. */

#line 260 "dptrfs.f"
	if (*n == 1) {
#line 261 "dptrfs.f"
	    bi = b[j * b_dim1 + 1];
#line 262 "dptrfs.f"
	    dx = d__[1] * x[j * x_dim1 + 1];
#line 263 "dptrfs.f"
	    work[*n + 1] = bi - dx;
#line 264 "dptrfs.f"
	    work[1] = abs(bi) + abs(dx);
#line 265 "dptrfs.f"
	} else {
#line 266 "dptrfs.f"
	    bi = b[j * b_dim1 + 1];
#line 267 "dptrfs.f"
	    dx = d__[1] * x[j * x_dim1 + 1];
#line 268 "dptrfs.f"
	    ex = e[1] * x[j * x_dim1 + 2];
#line 269 "dptrfs.f"
	    work[*n + 1] = bi - dx - ex;
#line 270 "dptrfs.f"
	    work[1] = abs(bi) + abs(dx) + abs(ex);
#line 271 "dptrfs.f"
	    i__2 = *n - 1;
#line 271 "dptrfs.f"
	    for (i__ = 2; i__ <= i__2; ++i__) {
#line 272 "dptrfs.f"
		bi = b[i__ + j * b_dim1];
#line 273 "dptrfs.f"
		cx = e[i__ - 1] * x[i__ - 1 + j * x_dim1];
#line 274 "dptrfs.f"
		dx = d__[i__] * x[i__ + j * x_dim1];
#line 275 "dptrfs.f"
		ex = e[i__] * x[i__ + 1 + j * x_dim1];
#line 276 "dptrfs.f"
		work[*n + i__] = bi - cx - dx - ex;
#line 277 "dptrfs.f"
		work[i__] = abs(bi) + abs(cx) + abs(dx) + abs(ex);
#line 278 "dptrfs.f"
/* L30: */
#line 278 "dptrfs.f"
	    }
#line 279 "dptrfs.f"
	    bi = b[*n + j * b_dim1];
#line 280 "dptrfs.f"
	    cx = e[*n - 1] * x[*n - 1 + j * x_dim1];
#line 281 "dptrfs.f"
	    dx = d__[*n] * x[*n + j * x_dim1];
#line 282 "dptrfs.f"
	    work[*n + *n] = bi - cx - dx;
#line 283 "dptrfs.f"
	    work[*n] = abs(bi) + abs(cx) + abs(dx);
#line 284 "dptrfs.f"
	}

/*        Compute componentwise relative backward error from formula */

/*        max(i) ( abs(R(i)) / ( abs(A)*abs(X) + abs(B) )(i) ) */

/*        where abs(Z) is the componentwise absolute value of the matrix */
/*        or vector Z.  If the i-th component of the denominator is less */
/*        than SAFE2, then SAFE1 is added to the i-th components of the */
/*        numerator and denominator before dividing. */

#line 295 "dptrfs.f"
	s = 0.;
#line 296 "dptrfs.f"
	i__2 = *n;
#line 296 "dptrfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 297 "dptrfs.f"
	    if (work[i__] > safe2) {
/* Computing MAX */
#line 298 "dptrfs.f"
		d__2 = s, d__3 = (d__1 = work[*n + i__], abs(d__1)) / work[
			i__];
#line 298 "dptrfs.f"
		s = max(d__2,d__3);
#line 299 "dptrfs.f"
	    } else {
/* Computing MAX */
#line 300 "dptrfs.f"
		d__2 = s, d__3 = ((d__1 = work[*n + i__], abs(d__1)) + safe1) 
			/ (work[i__] + safe1);
#line 300 "dptrfs.f"
		s = max(d__2,d__3);
#line 302 "dptrfs.f"
	    }
#line 303 "dptrfs.f"
/* L40: */
#line 303 "dptrfs.f"
	}
#line 304 "dptrfs.f"
	berr[j] = s;

/*        Test stopping criterion. Continue iterating if */
/*           1) The residual BERR(J) is larger than machine epsilon, and */
/*           2) BERR(J) decreased by at least a factor of 2 during the */
/*              last iteration, and */
/*           3) At most ITMAX iterations tried. */

#line 312 "dptrfs.f"
	if (berr[j] > eps && berr[j] * 2. <= lstres && count <= 5) {

/*           Update solution and try again. */

#line 317 "dptrfs.f"
	    dpttrs_(n, &c__1, &df[1], &ef[1], &work[*n + 1], n, info);
#line 318 "dptrfs.f"
	    daxpy_(n, &c_b11, &work[*n + 1], &c__1, &x[j * x_dim1 + 1], &c__1)
		    ;
#line 319 "dptrfs.f"
	    lstres = berr[j];
#line 320 "dptrfs.f"
	    ++count;
#line 321 "dptrfs.f"
	    goto L20;
#line 322 "dptrfs.f"
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

#line 342 "dptrfs.f"
	i__2 = *n;
#line 342 "dptrfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 343 "dptrfs.f"
	    if (work[i__] > safe2) {
#line 344 "dptrfs.f"
		work[i__] = (d__1 = work[*n + i__], abs(d__1)) + nz * eps * 
			work[i__];
#line 345 "dptrfs.f"
	    } else {
#line 346 "dptrfs.f"
		work[i__] = (d__1 = work[*n + i__], abs(d__1)) + nz * eps * 
			work[i__] + safe1;
#line 347 "dptrfs.f"
	    }
#line 348 "dptrfs.f"
/* L50: */
#line 348 "dptrfs.f"
	}
#line 349 "dptrfs.f"
	ix = idamax_(n, &work[1], &c__1);
#line 350 "dptrfs.f"
	ferr[j] = work[ix];

/*        Estimate the norm of inv(A). */

/*        Solve M(A) * x = e, where M(A) = (m(i,j)) is given by */

/*           m(i,j) =  abs(A(i,j)), i = j, */
/*           m(i,j) = -abs(A(i,j)), i .ne. j, */

/*        and e = [ 1, 1, ..., 1 ]**T.  Note M(A) = M(L)*D*M(L)**T. */

/*        Solve M(L) * x = e. */

#line 363 "dptrfs.f"
	work[1] = 1.;
#line 364 "dptrfs.f"
	i__2 = *n;
#line 364 "dptrfs.f"
	for (i__ = 2; i__ <= i__2; ++i__) {
#line 365 "dptrfs.f"
	    work[i__] = work[i__ - 1] * (d__1 = ef[i__ - 1], abs(d__1)) + 1.;
#line 366 "dptrfs.f"
/* L60: */
#line 366 "dptrfs.f"
	}

/*        Solve D * M(L)**T * x = b. */

#line 370 "dptrfs.f"
	work[*n] /= df[*n];
#line 371 "dptrfs.f"
	for (i__ = *n - 1; i__ >= 1; --i__) {
#line 372 "dptrfs.f"
	    work[i__] = work[i__] / df[i__] + work[i__ + 1] * (d__1 = ef[i__],
		     abs(d__1));
#line 373 "dptrfs.f"
/* L70: */
#line 373 "dptrfs.f"
	}

/*        Compute norm(inv(A)) = max(x(i)), 1<=i<=n. */

#line 377 "dptrfs.f"
	ix = idamax_(n, &work[1], &c__1);
#line 378 "dptrfs.f"
	ferr[j] *= (d__1 = work[ix], abs(d__1));

/*        Normalize error. */

#line 382 "dptrfs.f"
	lstres = 0.;
#line 383 "dptrfs.f"
	i__2 = *n;
#line 383 "dptrfs.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 384 "dptrfs.f"
	    d__2 = lstres, d__3 = (d__1 = x[i__ + j * x_dim1], abs(d__1));
#line 384 "dptrfs.f"
	    lstres = max(d__2,d__3);
#line 385 "dptrfs.f"
/* L80: */
#line 385 "dptrfs.f"
	}
#line 386 "dptrfs.f"
	if (lstres != 0.) {
#line 386 "dptrfs.f"
	    ferr[j] /= lstres;
#line 386 "dptrfs.f"
	}

#line 389 "dptrfs.f"
/* L90: */
#line 389 "dptrfs.f"
    }

#line 391 "dptrfs.f"
    return 0;

/*     End of DPTRFS */

} /* dptrfs_ */

