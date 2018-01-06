#line 1 "dsgesv.f"
/* dsgesv.f -- translated by f2c (version 20100827).
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

#line 1 "dsgesv.f"
/* Table of constant values */

static doublereal c_b10 = -1.;
static doublereal c_b11 = 1.;
static integer c__1 = 1;

/* > \brief <b> DSGESV computes the solution to system of linear equations A * X = B for GE matrices</b> (mixe
d precision with iterative refinement) */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DSGESV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsgesv.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsgesv.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsgesv.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSGESV( N, NRHS, A, LDA, IPIV, B, LDB, X, LDX, WORK, */
/*                          SWORK, ITER, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, ITER, LDA, LDB, LDX, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       REAL               SWORK( * ) */
/*       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), WORK( N, * ), */
/*      $                   X( LDX, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DSGESV computes the solution to a real system of linear equations */
/* >    A * X = B, */
/* > where A is an N-by-N matrix and X and B are N-by-NRHS matrices. */
/* > */
/* > DSGESV first attempts to factorize the matrix in SINGLE PRECISION */
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
/* >          On entry, the N-by-N coefficient matrix A. */
/* >          On exit, if iterative refinement has been successfully used */
/* >          (INFO.EQ.0 and ITER.GE.0, see description below), then A is */
/* >          unchanged, if double precision factorization has been used */
/* >          (INFO.EQ.0 and ITER.LT.0, see description below), then the */
/* >          array A contains the factors L and U from the factorization */
/* >          A = P*L*U; the unit diagonal elements of L are not stored. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N) */
/* >          The pivot indices that define the permutation matrix P; */
/* >          row i of the matrix was interchanged with row IPIV(i). */
/* >          Corresponds either to the single precision factorization */
/* >          (if INFO.EQ.0 and ITER.GE.0) or the double precision */
/* >          factorization (if INFO.EQ.0 and ITER.LT.0). */
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
/* >               -3 : failure of SGETRF */
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
/* >          > 0:  if INFO = i, U(i,i) computed in DOUBLE PRECISION is */
/* >                exactly zero.  The factorization has been completed, */
/* >                but the factor U is exactly singular, so the solution */
/* >                could not be computed. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date June 2016 */

/* > \ingroup doubleGEsolve */

/*  ===================================================================== */
/* Subroutine */ int dsgesv_(integer *n, integer *nrhs, doublereal *a, 
	integer *lda, integer *ipiv, doublereal *b, integer *ldb, doublereal *
	x, integer *ldx, doublereal *work, doublereal *swork, integer *iter, 
	integer *info)
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
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    static integer iiter;
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *), dlag2s_(integer *, integer *,
	     doublereal *, integer *, doublereal *, integer *, integer *), 
	    slag2d_(integer *, integer *, doublereal *, integer *, doublereal 
	    *, integer *, integer *);
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern integer idamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen), dgetrf_(integer *, integer *, 
	    doublereal *, integer *, integer *, integer *), dgetrs_(char *, 
	    integer *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, integer *, integer *, ftnlen), sgetrf_(integer *, 
	    integer *, doublereal *, integer *, integer *, integer *), 
	    sgetrs_(char *, integer *, integer *, doublereal *, integer *, 
	    integer *, doublereal *, integer *, integer *, ftnlen);


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

#line 246 "dsgesv.f"
    /* Parameter adjustments */
#line 246 "dsgesv.f"
    work_dim1 = *n;
#line 246 "dsgesv.f"
    work_offset = 1 + work_dim1;
#line 246 "dsgesv.f"
    work -= work_offset;
#line 246 "dsgesv.f"
    a_dim1 = *lda;
#line 246 "dsgesv.f"
    a_offset = 1 + a_dim1;
#line 246 "dsgesv.f"
    a -= a_offset;
#line 246 "dsgesv.f"
    --ipiv;
#line 246 "dsgesv.f"
    b_dim1 = *ldb;
#line 246 "dsgesv.f"
    b_offset = 1 + b_dim1;
#line 246 "dsgesv.f"
    b -= b_offset;
#line 246 "dsgesv.f"
    x_dim1 = *ldx;
#line 246 "dsgesv.f"
    x_offset = 1 + x_dim1;
#line 246 "dsgesv.f"
    x -= x_offset;
#line 246 "dsgesv.f"
    --swork;
#line 246 "dsgesv.f"

#line 246 "dsgesv.f"
    /* Function Body */
#line 246 "dsgesv.f"
    *info = 0;
#line 247 "dsgesv.f"
    *iter = 0;

/*     Test the input parameters. */

#line 251 "dsgesv.f"
    if (*n < 0) {
#line 252 "dsgesv.f"
	*info = -1;
#line 253 "dsgesv.f"
    } else if (*nrhs < 0) {
#line 254 "dsgesv.f"
	*info = -2;
#line 255 "dsgesv.f"
    } else if (*lda < max(1,*n)) {
#line 256 "dsgesv.f"
	*info = -4;
#line 257 "dsgesv.f"
    } else if (*ldb < max(1,*n)) {
#line 258 "dsgesv.f"
	*info = -7;
#line 259 "dsgesv.f"
    } else if (*ldx < max(1,*n)) {
#line 260 "dsgesv.f"
	*info = -9;
#line 261 "dsgesv.f"
    }
#line 262 "dsgesv.f"
    if (*info != 0) {
#line 263 "dsgesv.f"
	i__1 = -(*info);
#line 263 "dsgesv.f"
	xerbla_("DSGESV", &i__1, (ftnlen)6);
#line 264 "dsgesv.f"
	return 0;
#line 265 "dsgesv.f"
    }

/*     Quick return if (N.EQ.0). */

#line 269 "dsgesv.f"
    if (*n == 0) {
#line 269 "dsgesv.f"
	return 0;
#line 269 "dsgesv.f"
    }

/*     Skip single precision iterative refinement if a priori slower */
/*     than double precision factorization. */

#line 275 "dsgesv.f"
    if (FALSE_) {
#line 276 "dsgesv.f"
	*iter = -1;
#line 277 "dsgesv.f"
	goto L40;
#line 278 "dsgesv.f"
    }

/*     Compute some constants. */

#line 282 "dsgesv.f"
    anrm = dlange_("I", n, n, &a[a_offset], lda, &work[work_offset], (ftnlen)
	    1);
#line 283 "dsgesv.f"
    eps = dlamch_("Epsilon", (ftnlen)7);
#line 284 "dsgesv.f"
    cte = anrm * eps * sqrt((doublereal) (*n)) * 1.;

/*     Set the indices PTSA, PTSX for referencing SA and SX in SWORK. */

#line 288 "dsgesv.f"
    ptsa = 1;
#line 289 "dsgesv.f"
    ptsx = ptsa + *n * *n;

/*     Convert B from double precision to single precision and store the */
/*     result in SX. */

#line 294 "dsgesv.f"
    dlag2s_(n, nrhs, &b[b_offset], ldb, &swork[ptsx], n, info);

#line 296 "dsgesv.f"
    if (*info != 0) {
#line 297 "dsgesv.f"
	*iter = -2;
#line 298 "dsgesv.f"
	goto L40;
#line 299 "dsgesv.f"
    }

/*     Convert A from double precision to single precision and store the */
/*     result in SA. */

#line 304 "dsgesv.f"
    dlag2s_(n, n, &a[a_offset], lda, &swork[ptsa], n, info);

#line 306 "dsgesv.f"
    if (*info != 0) {
#line 307 "dsgesv.f"
	*iter = -2;
#line 308 "dsgesv.f"
	goto L40;
#line 309 "dsgesv.f"
    }

/*     Compute the LU factorization of SA. */

#line 313 "dsgesv.f"
    sgetrf_(n, n, &swork[ptsa], n, &ipiv[1], info);

#line 315 "dsgesv.f"
    if (*info != 0) {
#line 316 "dsgesv.f"
	*iter = -3;
#line 317 "dsgesv.f"
	goto L40;
#line 318 "dsgesv.f"
    }

/*     Solve the system SA*SX = SB. */

#line 322 "dsgesv.f"
    sgetrs_("No transpose", n, nrhs, &swork[ptsa], n, &ipiv[1], &swork[ptsx], 
	    n, info, (ftnlen)12);

/*     Convert SX back to double precision */

#line 327 "dsgesv.f"
    slag2d_(n, nrhs, &swork[ptsx], n, &x[x_offset], ldx, info);

/*     Compute R = B - AX (R is WORK). */

#line 331 "dsgesv.f"
    dlacpy_("All", n, nrhs, &b[b_offset], ldb, &work[work_offset], n, (ftnlen)
	    3);

#line 333 "dsgesv.f"
    dgemm_("No Transpose", "No Transpose", n, nrhs, n, &c_b10, &a[a_offset], 
	    lda, &x[x_offset], ldx, &c_b11, &work[work_offset], n, (ftnlen)12,
	     (ftnlen)12);

/*     Check whether the NRHS normwise backward errors satisfy the */
/*     stopping criterion. If yes, set ITER=0 and return. */

#line 339 "dsgesv.f"
    i__1 = *nrhs;
#line 339 "dsgesv.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 340 "dsgesv.f"
	xnrm = (d__1 = x[idamax_(n, &x[i__ * x_dim1 + 1], &c__1) + i__ * 
		x_dim1], abs(d__1));
#line 341 "dsgesv.f"
	rnrm = (d__1 = work[idamax_(n, &work[i__ * work_dim1 + 1], &c__1) + 
		i__ * work_dim1], abs(d__1));
#line 342 "dsgesv.f"
	if (rnrm > xnrm * cte) {
#line 342 "dsgesv.f"
	    goto L10;
#line 342 "dsgesv.f"
	}
#line 344 "dsgesv.f"
    }

/*     If we are here, the NRHS normwise backward errors satisfy the */
/*     stopping criterion. We are good to exit. */

#line 349 "dsgesv.f"
    *iter = 0;
#line 350 "dsgesv.f"
    return 0;

#line 352 "dsgesv.f"
L10:

#line 354 "dsgesv.f"
    for (iiter = 1; iiter <= 30; ++iiter) {

/*        Convert R (in WORK) from double precision to single precision */
/*        and store the result in SX. */

#line 359 "dsgesv.f"
	dlag2s_(n, nrhs, &work[work_offset], n, &swork[ptsx], n, info);

#line 361 "dsgesv.f"
	if (*info != 0) {
#line 362 "dsgesv.f"
	    *iter = -2;
#line 363 "dsgesv.f"
	    goto L40;
#line 364 "dsgesv.f"
	}

/*        Solve the system SA*SX = SR. */

#line 368 "dsgesv.f"
	sgetrs_("No transpose", n, nrhs, &swork[ptsa], n, &ipiv[1], &swork[
		ptsx], n, info, (ftnlen)12);

/*        Convert SX back to double precision and update the current */
/*        iterate. */

#line 374 "dsgesv.f"
	slag2d_(n, nrhs, &swork[ptsx], n, &work[work_offset], n, info);

#line 376 "dsgesv.f"
	i__1 = *nrhs;
#line 376 "dsgesv.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 377 "dsgesv.f"
	    daxpy_(n, &c_b11, &work[i__ * work_dim1 + 1], &c__1, &x[i__ * 
		    x_dim1 + 1], &c__1);
#line 378 "dsgesv.f"
	}

/*        Compute R = B - AX (R is WORK). */

#line 382 "dsgesv.f"
	dlacpy_("All", n, nrhs, &b[b_offset], ldb, &work[work_offset], n, (
		ftnlen)3);

#line 384 "dsgesv.f"
	dgemm_("No Transpose", "No Transpose", n, nrhs, n, &c_b10, &a[
		a_offset], lda, &x[x_offset], ldx, &c_b11, &work[work_offset],
		 n, (ftnlen)12, (ftnlen)12);

/*        Check whether the NRHS normwise backward errors satisfy the */
/*        stopping criterion. If yes, set ITER=IITER>0 and return. */

#line 390 "dsgesv.f"
	i__1 = *nrhs;
#line 390 "dsgesv.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 391 "dsgesv.f"
	    xnrm = (d__1 = x[idamax_(n, &x[i__ * x_dim1 + 1], &c__1) + i__ * 
		    x_dim1], abs(d__1));
#line 392 "dsgesv.f"
	    rnrm = (d__1 = work[idamax_(n, &work[i__ * work_dim1 + 1], &c__1) 
		    + i__ * work_dim1], abs(d__1));
#line 393 "dsgesv.f"
	    if (rnrm > xnrm * cte) {
#line 393 "dsgesv.f"
		goto L20;
#line 393 "dsgesv.f"
	    }
#line 395 "dsgesv.f"
	}

/*        If we are here, the NRHS normwise backward errors satisfy the */
/*        stopping criterion, we are good to exit. */

#line 400 "dsgesv.f"
	*iter = iiter;

#line 402 "dsgesv.f"
	return 0;

#line 404 "dsgesv.f"
L20:

#line 406 "dsgesv.f"
/* L30: */
#line 406 "dsgesv.f"
	;
#line 406 "dsgesv.f"
    }

/*     If we are at this place of the code, this is because we have */
/*     performed ITER=ITERMAX iterations and never satisified the */
/*     stopping criterion, set up the ITER flag accordingly and follow up */
/*     on double precision routine. */

#line 413 "dsgesv.f"
    *iter = -31;

#line 415 "dsgesv.f"
L40:

/*     Single-precision iterative refinement failed to converge to a */
/*     satisfactory solution, so we resort to double precision. */

#line 420 "dsgesv.f"
    dgetrf_(n, n, &a[a_offset], lda, &ipiv[1], info);

#line 422 "dsgesv.f"
    if (*info != 0) {
#line 422 "dsgesv.f"
	return 0;
#line 422 "dsgesv.f"
    }

#line 425 "dsgesv.f"
    dlacpy_("All", n, nrhs, &b[b_offset], ldb, &x[x_offset], ldx, (ftnlen)3);
#line 426 "dsgesv.f"
    dgetrs_("No transpose", n, nrhs, &a[a_offset], lda, &ipiv[1], &x[x_offset]
	    , ldx, info, (ftnlen)12);

#line 429 "dsgesv.f"
    return 0;

/*     End of DSGESV. */

} /* dsgesv_ */

