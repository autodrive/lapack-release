#line 1 "dsposv.f"
/* dsposv.f -- translated by f2c (version 20100827).
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

#line 1 "dsposv.f"
/* Table of constant values */

static doublereal c_b10 = -1.;
static doublereal c_b11 = 1.;
static integer c__1 = 1;

/* > \brief <b> DSPOSV computes the solution to system of linear equations A * X = B for PO matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DSPOSV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsposv.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsposv.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsposv.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSPOSV( UPLO, N, NRHS, A, LDA, B, LDB, X, LDX, WORK, */
/*                          SWORK, ITER, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, ITER, LDA, LDB, LDX, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               SWORK( * ) */
/*       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), WORK( N, * ), */
/*      $                   X( LDX, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DSPOSV computes the solution to a real system of linear equations */
/* >    A * X = B, */
/* > where A is an N-by-N symmetric positive definite matrix and X and B */
/* > are N-by-NRHS matrices. */
/* > */
/* > DSPOSV first attempts to factorize the matrix in SINGLE PRECISION */
/* > and use this factorization within an iterative refinement procedure */
/* > to produce a solution with DOUBLE PRECISION normwise backward error */
/* > quality (see below). If the approach fails the method switches to a */
/* > DOUBLE PRECISION factorization and solve. */
/* > */
/* > The iterative refinement is not going to be a winning strategy if */
/* > the ratio SINGLE PRECISION performance over DOUBLE PRECISION */
/* > performance is too small. A reasonable strategy should take the */
/* > number of right-hand sides and the size of the matrix into account. */
/* > This might be done with a call to ILAENV in the future. Up to now, we */
/* > always try iterative refinement. */
/* > */
/* > The iterative refinement process is stopped if */
/* >     ITER > ITERMAX */
/* > or for all the RHS we have: */
/* >     RNRM < SQRT(N)*XNRM*ANRM*EPS*BWDMAX */
/* > where */
/* >     o ITER is the number of the current iteration in the iterative */
/* >       refinement process */
/* >     o RNRM is the infinity-norm of the residual */
/* >     o XNRM is the infinity-norm of the solution */
/* >     o ANRM is the infinity-operator-norm of the matrix A */
/* >     o EPS is the machine epsilon returned by DLAMCH('Epsilon') */
/* > The value ITERMAX and BWDMAX are fixed to 30 and 1.0D+00 */
/* > respectively. */
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
/* >          The number of linear equations, i.e., the order of the */
/* >          matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NRHS */
/* > \verbatim */
/* >          NRHS is INTEGER */
/* >          The number of right hand sides, i.e., the number of columns */
/* >          of the matrix B.  NRHS >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, */
/* >          dimension (LDA,N) */
/* >          On entry, the symmetric matrix A.  If UPLO = 'U', the leading */
/* >          N-by-N upper triangular part of A contains the upper */
/* >          triangular part of the matrix A, and the strictly lower */
/* >          triangular part of A is not referenced.  If UPLO = 'L', the */
/* >          leading N-by-N lower triangular part of A contains the lower */
/* >          triangular part of the matrix A, and the strictly upper */
/* >          triangular part of A is not referenced. */
/* >          On exit, if iterative refinement has been successfully used */
/* >          (INFO.EQ.0 and ITER.GE.0, see description below), then A is */
/* >          unchanged, if double precision factorization has been used */
/* >          (INFO.EQ.0 and ITER.LT.0, see description below), then the */
/* >          array A contains the factor U or L from the Cholesky */
/* >          factorization A = U**T*U or A = L*L**T. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] B */
/* > \verbatim */
/* >          B is DOUBLE PRECISION array, dimension (LDB,NRHS) */
/* >          The N-by-NRHS right hand side matrix B. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >          The leading dimension of the array B.  LDB >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] X */
/* > \verbatim */
/* >          X is DOUBLE PRECISION array, dimension (LDX,NRHS) */
/* >          If INFO = 0, the N-by-NRHS solution matrix X. */
/* > \endverbatim */
/* > */
/* > \param[in] LDX */
/* > \verbatim */
/* >          LDX is INTEGER */
/* >          The leading dimension of the array X.  LDX >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension (N,NRHS) */
/* >          This array is used to hold the residual vectors. */
/* > \endverbatim */
/* > */
/* > \param[out] SWORK */
/* > \verbatim */
/* >          SWORK is REAL array, dimension (N*(N+NRHS)) */
/* >          This array is used to use the single precision matrix and the */
/* >          right-hand sides or solutions in single precision. */
/* > \endverbatim */
/* > */
/* > \param[out] ITER */
/* > \verbatim */
/* >          ITER is INTEGER */
/* >          < 0: iterative refinement has failed, double precision */
/* >               factorization has been performed */
/* >               -1 : the routine fell back to full precision for */
/* >                    implementation- or machine-specific reasons */
/* >               -2 : narrowing the precision induced an overflow, */
/* >                    the routine fell back to full precision */
/* >               -3 : failure of SPOTRF */
/* >               -31: stop the iterative refinement after the 30th */
/* >                    iterations */
/* >          > 0: iterative refinement has been successfully used. */
/* >               Returns the number of iterations */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* >          > 0:  if INFO = i, the leading minor of order i of (DOUBLE */
/* >                PRECISION) A is not positive definite, so the */
/* >                factorization could not be completed, and the solution */
/* >                has not been computed. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date June 2016 */

/* > \ingroup doublePOsolve */

/*  ===================================================================== */
/* Subroutine */ int dsposv_(char *uplo, integer *n, integer *nrhs, 
	doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *
	x, integer *ldx, doublereal *work, doublereal *swork, integer *iter, 
	integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, work_dim1, work_offset, 
	    x_dim1, x_offset, i__1;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__;
    static doublereal cte, eps, anrm;
    static integer ptsa;
    static doublereal rnrm, xnrm;
    static integer ptsx;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer iiter;
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *), dsymm_(char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    ftnlen, ftnlen), dlag2s_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *), slag2d_(integer *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *), dlat2s_(char *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, ftnlen);
    extern doublereal dlamch_(char *, ftnlen);
    extern integer idamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    extern doublereal dlansy_(char *, char *, integer *, doublereal *, 
	    integer *, doublereal *, ftnlen, ftnlen);
    extern /* Subroutine */ int dpotrf_(char *, integer *, doublereal *, 
	    integer *, integer *, ftnlen), dpotrs_(char *, integer *, integer 
	    *, doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen), spotrf_(char *, integer *, doublereal *, integer *, 
	    integer *, ftnlen), spotrs_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen);


/*  -- LAPACK driver routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     June 2016 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Parameters .. */




/*     .. Local Scalars .. */

/*     .. External Subroutines .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 251 "dsposv.f"
    /* Parameter adjustments */
#line 251 "dsposv.f"
    work_dim1 = *n;
#line 251 "dsposv.f"
    work_offset = 1 + work_dim1;
#line 251 "dsposv.f"
    work -= work_offset;
#line 251 "dsposv.f"
    a_dim1 = *lda;
#line 251 "dsposv.f"
    a_offset = 1 + a_dim1;
#line 251 "dsposv.f"
    a -= a_offset;
#line 251 "dsposv.f"
    b_dim1 = *ldb;
#line 251 "dsposv.f"
    b_offset = 1 + b_dim1;
#line 251 "dsposv.f"
    b -= b_offset;
#line 251 "dsposv.f"
    x_dim1 = *ldx;
#line 251 "dsposv.f"
    x_offset = 1 + x_dim1;
#line 251 "dsposv.f"
    x -= x_offset;
#line 251 "dsposv.f"
    --swork;
#line 251 "dsposv.f"

#line 251 "dsposv.f"
    /* Function Body */
#line 251 "dsposv.f"
    *info = 0;
#line 252 "dsposv.f"
    *iter = 0;

/*     Test the input parameters. */

#line 256 "dsposv.f"
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
#line 257 "dsposv.f"
	*info = -1;
#line 258 "dsposv.f"
    } else if (*n < 0) {
#line 259 "dsposv.f"
	*info = -2;
#line 260 "dsposv.f"
    } else if (*nrhs < 0) {
#line 261 "dsposv.f"
	*info = -3;
#line 262 "dsposv.f"
    } else if (*lda < max(1,*n)) {
#line 263 "dsposv.f"
	*info = -5;
#line 264 "dsposv.f"
    } else if (*ldb < max(1,*n)) {
#line 265 "dsposv.f"
	*info = -7;
#line 266 "dsposv.f"
    } else if (*ldx < max(1,*n)) {
#line 267 "dsposv.f"
	*info = -9;
#line 268 "dsposv.f"
    }
#line 269 "dsposv.f"
    if (*info != 0) {
#line 270 "dsposv.f"
	i__1 = -(*info);
#line 270 "dsposv.f"
	xerbla_("DSPOSV", &i__1, (ftnlen)6);
#line 271 "dsposv.f"
	return 0;
#line 272 "dsposv.f"
    }

/*     Quick return if (N.EQ.0). */

#line 276 "dsposv.f"
    if (*n == 0) {
#line 276 "dsposv.f"
	return 0;
#line 276 "dsposv.f"
    }

/*     Skip single precision iterative refinement if a priori slower */
/*     than double precision factorization. */

#line 282 "dsposv.f"
    if (FALSE_) {
#line 283 "dsposv.f"
	*iter = -1;
#line 284 "dsposv.f"
	goto L40;
#line 285 "dsposv.f"
    }

/*     Compute some constants. */

#line 289 "dsposv.f"
    anrm = dlansy_("I", uplo, n, &a[a_offset], lda, &work[work_offset], (
	    ftnlen)1, (ftnlen)1);
#line 290 "dsposv.f"
    eps = dlamch_("Epsilon", (ftnlen)7);
#line 291 "dsposv.f"
    cte = anrm * eps * sqrt((doublereal) (*n)) * 1.;

/*     Set the indices PTSA, PTSX for referencing SA and SX in SWORK. */

#line 295 "dsposv.f"
    ptsa = 1;
#line 296 "dsposv.f"
    ptsx = ptsa + *n * *n;

/*     Convert B from double precision to single precision and store the */
/*     result in SX. */

#line 301 "dsposv.f"
    dlag2s_(n, nrhs, &b[b_offset], ldb, &swork[ptsx], n, info);

#line 303 "dsposv.f"
    if (*info != 0) {
#line 304 "dsposv.f"
	*iter = -2;
#line 305 "dsposv.f"
	goto L40;
#line 306 "dsposv.f"
    }

/*     Convert A from double precision to single precision and store the */
/*     result in SA. */

#line 311 "dsposv.f"
    dlat2s_(uplo, n, &a[a_offset], lda, &swork[ptsa], n, info, (ftnlen)1);

#line 313 "dsposv.f"
    if (*info != 0) {
#line 314 "dsposv.f"
	*iter = -2;
#line 315 "dsposv.f"
	goto L40;
#line 316 "dsposv.f"
    }

/*     Compute the Cholesky factorization of SA. */

#line 320 "dsposv.f"
    spotrf_(uplo, n, &swork[ptsa], n, info, (ftnlen)1);

#line 322 "dsposv.f"
    if (*info != 0) {
#line 323 "dsposv.f"
	*iter = -3;
#line 324 "dsposv.f"
	goto L40;
#line 325 "dsposv.f"
    }

/*     Solve the system SA*SX = SB. */

#line 329 "dsposv.f"
    spotrs_(uplo, n, nrhs, &swork[ptsa], n, &swork[ptsx], n, info, (ftnlen)1);

/*     Convert SX back to double precision */

#line 334 "dsposv.f"
    slag2d_(n, nrhs, &swork[ptsx], n, &x[x_offset], ldx, info);

/*     Compute R = B - AX (R is WORK). */

#line 338 "dsposv.f"
    dlacpy_("All", n, nrhs, &b[b_offset], ldb, &work[work_offset], n, (ftnlen)
	    3);

#line 340 "dsposv.f"
    dsymm_("Left", uplo, n, nrhs, &c_b10, &a[a_offset], lda, &x[x_offset], 
	    ldx, &c_b11, &work[work_offset], n, (ftnlen)4, (ftnlen)1);

/*     Check whether the NRHS normwise backward errors satisfy the */
/*     stopping criterion. If yes, set ITER=0 and return. */

#line 346 "dsposv.f"
    i__1 = *nrhs;
#line 346 "dsposv.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 347 "dsposv.f"
	xnrm = (d__1 = x[idamax_(n, &x[i__ * x_dim1 + 1], &c__1) + i__ * 
		x_dim1], abs(d__1));
#line 348 "dsposv.f"
	rnrm = (d__1 = work[idamax_(n, &work[i__ * work_dim1 + 1], &c__1) + 
		i__ * work_dim1], abs(d__1));
#line 349 "dsposv.f"
	if (rnrm > xnrm * cte) {
#line 349 "dsposv.f"
	    goto L10;
#line 349 "dsposv.f"
	}
#line 351 "dsposv.f"
    }

/*     If we are here, the NRHS normwise backward errors satisfy the */
/*     stopping criterion. We are good to exit. */

#line 356 "dsposv.f"
    *iter = 0;
#line 357 "dsposv.f"
    return 0;

#line 359 "dsposv.f"
L10:

#line 361 "dsposv.f"
    for (iiter = 1; iiter <= 30; ++iiter) {

/*        Convert R (in WORK) from double precision to single precision */
/*        and store the result in SX. */

#line 366 "dsposv.f"
	dlag2s_(n, nrhs, &work[work_offset], n, &swork[ptsx], n, info);

#line 368 "dsposv.f"
	if (*info != 0) {
#line 369 "dsposv.f"
	    *iter = -2;
#line 370 "dsposv.f"
	    goto L40;
#line 371 "dsposv.f"
	}

/*        Solve the system SA*SX = SR. */

#line 375 "dsposv.f"
	spotrs_(uplo, n, nrhs, &swork[ptsa], n, &swork[ptsx], n, info, (
		ftnlen)1);

/*        Convert SX back to double precision and update the current */
/*        iterate. */

#line 381 "dsposv.f"
	slag2d_(n, nrhs, &swork[ptsx], n, &work[work_offset], n, info);

#line 383 "dsposv.f"
	i__1 = *nrhs;
#line 383 "dsposv.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 384 "dsposv.f"
	    daxpy_(n, &c_b11, &work[i__ * work_dim1 + 1], &c__1, &x[i__ * 
		    x_dim1 + 1], &c__1);
#line 385 "dsposv.f"
	}

/*        Compute R = B - AX (R is WORK). */

#line 389 "dsposv.f"
	dlacpy_("All", n, nrhs, &b[b_offset], ldb, &work[work_offset], n, (
		ftnlen)3);

#line 391 "dsposv.f"
	dsymm_("L", uplo, n, nrhs, &c_b10, &a[a_offset], lda, &x[x_offset], 
		ldx, &c_b11, &work[work_offset], n, (ftnlen)1, (ftnlen)1);

/*        Check whether the NRHS normwise backward errors satisfy the */
/*        stopping criterion. If yes, set ITER=IITER>0 and return. */

#line 397 "dsposv.f"
	i__1 = *nrhs;
#line 397 "dsposv.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 398 "dsposv.f"
	    xnrm = (d__1 = x[idamax_(n, &x[i__ * x_dim1 + 1], &c__1) + i__ * 
		    x_dim1], abs(d__1));
#line 399 "dsposv.f"
	    rnrm = (d__1 = work[idamax_(n, &work[i__ * work_dim1 + 1], &c__1) 
		    + i__ * work_dim1], abs(d__1));
#line 400 "dsposv.f"
	    if (rnrm > xnrm * cte) {
#line 400 "dsposv.f"
		goto L20;
#line 400 "dsposv.f"
	    }
#line 402 "dsposv.f"
	}

/*        If we are here, the NRHS normwise backward errors satisfy the */
/*        stopping criterion, we are good to exit. */

#line 407 "dsposv.f"
	*iter = iiter;

#line 409 "dsposv.f"
	return 0;

#line 411 "dsposv.f"
L20:

#line 413 "dsposv.f"
/* L30: */
#line 413 "dsposv.f"
	;
#line 413 "dsposv.f"
    }

/*     If we are at this place of the code, this is because we have */
/*     performed ITER=ITERMAX iterations and never satisified the */
/*     stopping criterion, set up the ITER flag accordingly and follow */
/*     up on double precision routine. */

#line 420 "dsposv.f"
    *iter = -31;

#line 422 "dsposv.f"
L40:

/*     Single-precision iterative refinement failed to converge to a */
/*     satisfactory solution, so we resort to double precision. */

#line 427 "dsposv.f"
    dpotrf_(uplo, n, &a[a_offset], lda, info, (ftnlen)1);

#line 429 "dsposv.f"
    if (*info != 0) {
#line 429 "dsposv.f"
	return 0;
#line 429 "dsposv.f"
    }

#line 432 "dsposv.f"
    dlacpy_("All", n, nrhs, &b[b_offset], ldb, &x[x_offset], ldx, (ftnlen)3);
#line 433 "dsposv.f"
    dpotrs_(uplo, n, nrhs, &a[a_offset], lda, &x[x_offset], ldx, info, (
	    ftnlen)1);

#line 435 "dsposv.f"
    return 0;

/*     End of DSPOSV. */

} /* dsposv_ */

