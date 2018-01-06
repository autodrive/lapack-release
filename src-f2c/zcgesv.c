#line 1 "zcgesv.f"
/* zcgesv.f -- translated by f2c (version 20100827).
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

#line 1 "zcgesv.f"
/* Table of constant values */

static doublecomplex c_b1 = {-1.,0.};
static doublecomplex c_b2 = {1.,0.};
static integer c__1 = 1;

/* > \brief <b> ZCGESV computes the solution to system of linear equations A * X = B for GE matrices</b> (mixe
d precision with iterative refinement) */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZCGESV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zcgesv.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zcgesv.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zcgesv.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZCGESV( N, NRHS, A, LDA, IPIV, B, LDB, X, LDX, WORK, */
/*                          SWORK, RWORK, ITER, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, ITER, LDA, LDB, LDX, N, NRHS */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       DOUBLE PRECISION   RWORK( * ) */
/*       COMPLEX            SWORK( * ) */
/*       COMPLEX*16         A( LDA, * ), B( LDB, * ), WORK( N, * ), */
/*      $                   X( LDX, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZCGESV computes the solution to a complex system of linear equations */
/* >    A * X = B, */
/* > where A is an N-by-N matrix and X and B are N-by-NRHS matrices. */
/* > */
/* > ZCGESV first attempts to factorize the matrix in COMPLEX and use this */
/* > factorization within an iterative refinement procedure to produce a */
/* > solution with COMPLEX*16 normwise backward error quality (see below). */
/* > If the approach fails the method switches to a COMPLEX*16 */
/* > factorization and solve. */
/* > */
/* > The iterative refinement is not going to be a winning strategy if */
/* > the ratio COMPLEX performance over COMPLEX*16 performance is too */
/* > small. A reasonable strategy should take the number of right-hand */
/* > sides and the size of the matrix into account. This might be done */
/* > with a call to ILAENV in the future. Up to now, we always try */
/* > iterative refinement. */
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
/* >          A is COMPLEX*16 array, */
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
/* >          B is COMPLEX*16 array, dimension (LDB,NRHS) */
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
/* >          X is COMPLEX*16 array, dimension (LDX,NRHS) */
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
/* >          WORK is COMPLEX*16 array, dimension (N*NRHS) */
/* >          This array is used to hold the residual vectors. */
/* > \endverbatim */
/* > */
/* > \param[out] SWORK */
/* > \verbatim */
/* >          SWORK is COMPLEX array, dimension (N*(N+NRHS)) */
/* >          This array is used to use the single precision matrix and the */
/* >          right-hand sides or solutions in single precision. */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is DOUBLE PRECISION array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] ITER */
/* > \verbatim */
/* >          ITER is INTEGER */
/* >          < 0: iterative refinement has failed, COMPLEX*16 */
/* >               factorization has been performed */
/* >               -1 : the routine fell back to full precision for */
/* >                    implementation- or machine-specific reasons */
/* >               -2 : narrowing the precision induced an overflow, */
/* >                    the routine fell back to full precision */
/* >               -3 : failure of CGETRF */
/* >               -31: stop the iterative refinement after the 30th */
/* >                    iterations */
/* >          > 0: iterative refinement has been sucessfully used. */
/* >               Returns the number of iterations */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* >          > 0:  if INFO = i, U(i,i) computed in COMPLEX*16 is exactly */
/* >                zero.  The factorization has been completed, but the */
/* >                factor U is exactly singular, so the solution */
/* >                could not be computed. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup complex16GEsolve */

/*  ===================================================================== */
/* Subroutine */ int zcgesv_(integer *n, integer *nrhs, doublecomplex *a, 
	integer *lda, integer *ipiv, doublecomplex *b, integer *ldb, 
	doublecomplex *x, integer *ldx, doublecomplex *work, doublecomplex *
	swork, doublereal *rwork, integer *iter, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, work_dim1, work_offset, 
	    x_dim1, x_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal), d_imag(doublecomplex *);

    /* Local variables */
    static integer i__;
    static doublereal cte, eps, anrm;
    static integer ptsa;
    static doublereal rnrm, xnrm;
    static integer ptsx, iiter;
    extern /* Subroutine */ int zgemm_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, ftnlen, ftnlen), zaxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *), clag2z_(
	    integer *, integer *, doublecomplex *, integer *, doublecomplex *,
	     integer *, integer *), zlag2c_(integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *)
	    ;
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int cgetrf_(integer *, integer *, doublecomplex *,
	     integer *, integer *, integer *), xerbla_(char *, integer *, 
	    ftnlen);
    extern doublereal zlange_(char *, integer *, integer *, doublecomplex *, 
	    integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int cgetrs_(char *, integer *, integer *, 
	    doublecomplex *, integer *, integer *, doublecomplex *, integer *,
	     integer *, ftnlen);
    extern integer izamax_(integer *, doublecomplex *, integer *);
    extern /* Subroutine */ int zlacpy_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen), 
	    zgetrf_(integer *, integer *, doublecomplex *, integer *, integer 
	    *, integer *), zgetrs_(char *, integer *, integer *, 
	    doublecomplex *, integer *, integer *, doublecomplex *, integer *,
	     integer *, ftnlen);


/*  -- LAPACK driver routine (version 3.4.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */

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
/*     .. Statement Functions .. */
/*     .. */
/*     .. Statement Function definitions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 261 "zcgesv.f"
    /* Parameter adjustments */
#line 261 "zcgesv.f"
    work_dim1 = *n;
#line 261 "zcgesv.f"
    work_offset = 1 + work_dim1;
#line 261 "zcgesv.f"
    work -= work_offset;
#line 261 "zcgesv.f"
    a_dim1 = *lda;
#line 261 "zcgesv.f"
    a_offset = 1 + a_dim1;
#line 261 "zcgesv.f"
    a -= a_offset;
#line 261 "zcgesv.f"
    --ipiv;
#line 261 "zcgesv.f"
    b_dim1 = *ldb;
#line 261 "zcgesv.f"
    b_offset = 1 + b_dim1;
#line 261 "zcgesv.f"
    b -= b_offset;
#line 261 "zcgesv.f"
    x_dim1 = *ldx;
#line 261 "zcgesv.f"
    x_offset = 1 + x_dim1;
#line 261 "zcgesv.f"
    x -= x_offset;
#line 261 "zcgesv.f"
    --swork;
#line 261 "zcgesv.f"
    --rwork;
#line 261 "zcgesv.f"

#line 261 "zcgesv.f"
    /* Function Body */
#line 261 "zcgesv.f"
    *info = 0;
#line 262 "zcgesv.f"
    *iter = 0;

/*     Test the input parameters. */

#line 266 "zcgesv.f"
    if (*n < 0) {
#line 267 "zcgesv.f"
	*info = -1;
#line 268 "zcgesv.f"
    } else if (*nrhs < 0) {
#line 269 "zcgesv.f"
	*info = -2;
#line 270 "zcgesv.f"
    } else if (*lda < max(1,*n)) {
#line 271 "zcgesv.f"
	*info = -4;
#line 272 "zcgesv.f"
    } else if (*ldb < max(1,*n)) {
#line 273 "zcgesv.f"
	*info = -7;
#line 274 "zcgesv.f"
    } else if (*ldx < max(1,*n)) {
#line 275 "zcgesv.f"
	*info = -9;
#line 276 "zcgesv.f"
    }
#line 277 "zcgesv.f"
    if (*info != 0) {
#line 278 "zcgesv.f"
	i__1 = -(*info);
#line 278 "zcgesv.f"
	xerbla_("ZCGESV", &i__1, (ftnlen)6);
#line 279 "zcgesv.f"
	return 0;
#line 280 "zcgesv.f"
    }

/*     Quick return if (N.EQ.0). */

#line 284 "zcgesv.f"
    if (*n == 0) {
#line 284 "zcgesv.f"
	return 0;
#line 284 "zcgesv.f"
    }

/*     Skip single precision iterative refinement if a priori slower */
/*     than double precision factorization. */

#line 290 "zcgesv.f"
    if (FALSE_) {
#line 291 "zcgesv.f"
	*iter = -1;
#line 292 "zcgesv.f"
	goto L40;
#line 293 "zcgesv.f"
    }

/*     Compute some constants. */

#line 297 "zcgesv.f"
    anrm = zlange_("I", n, n, &a[a_offset], lda, &rwork[1], (ftnlen)1);
#line 298 "zcgesv.f"
    eps = dlamch_("Epsilon", (ftnlen)7);
#line 299 "zcgesv.f"
    cte = anrm * eps * sqrt((doublereal) (*n)) * 1.;

/*     Set the indices PTSA, PTSX for referencing SA and SX in SWORK. */

#line 303 "zcgesv.f"
    ptsa = 1;
#line 304 "zcgesv.f"
    ptsx = ptsa + *n * *n;

/*     Convert B from double precision to single precision and store the */
/*     result in SX. */

#line 309 "zcgesv.f"
    zlag2c_(n, nrhs, &b[b_offset], ldb, &swork[ptsx], n, info);

#line 311 "zcgesv.f"
    if (*info != 0) {
#line 312 "zcgesv.f"
	*iter = -2;
#line 313 "zcgesv.f"
	goto L40;
#line 314 "zcgesv.f"
    }

/*     Convert A from double precision to single precision and store the */
/*     result in SA. */

#line 319 "zcgesv.f"
    zlag2c_(n, n, &a[a_offset], lda, &swork[ptsa], n, info);

#line 321 "zcgesv.f"
    if (*info != 0) {
#line 322 "zcgesv.f"
	*iter = -2;
#line 323 "zcgesv.f"
	goto L40;
#line 324 "zcgesv.f"
    }

/*     Compute the LU factorization of SA. */

#line 328 "zcgesv.f"
    cgetrf_(n, n, &swork[ptsa], n, &ipiv[1], info);

#line 330 "zcgesv.f"
    if (*info != 0) {
#line 331 "zcgesv.f"
	*iter = -3;
#line 332 "zcgesv.f"
	goto L40;
#line 333 "zcgesv.f"
    }

/*     Solve the system SA*SX = SB. */

#line 337 "zcgesv.f"
    cgetrs_("No transpose", n, nrhs, &swork[ptsa], n, &ipiv[1], &swork[ptsx], 
	    n, info, (ftnlen)12);

/*     Convert SX back to double precision */

#line 342 "zcgesv.f"
    clag2z_(n, nrhs, &swork[ptsx], n, &x[x_offset], ldx, info);

/*     Compute R = B - AX (R is WORK). */

#line 346 "zcgesv.f"
    zlacpy_("All", n, nrhs, &b[b_offset], ldb, &work[work_offset], n, (ftnlen)
	    3);

#line 348 "zcgesv.f"
    zgemm_("No Transpose", "No Transpose", n, nrhs, n, &c_b1, &a[a_offset], 
	    lda, &x[x_offset], ldx, &c_b2, &work[work_offset], n, (ftnlen)12, 
	    (ftnlen)12);

/*     Check whether the NRHS normwise backward errors satisfy the */
/*     stopping criterion. If yes, set ITER=0 and return. */

#line 354 "zcgesv.f"
    i__1 = *nrhs;
#line 354 "zcgesv.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 355 "zcgesv.f"
	i__2 = izamax_(n, &x[i__ * x_dim1 + 1], &c__1) + i__ * x_dim1;
#line 355 "zcgesv.f"
	xnrm = (d__1 = x[i__2].r, abs(d__1)) + (d__2 = d_imag(&x[izamax_(n, &
		x[i__ * x_dim1 + 1], &c__1) + i__ * x_dim1]), abs(d__2));
#line 356 "zcgesv.f"
	i__2 = izamax_(n, &work[i__ * work_dim1 + 1], &c__1) + i__ * 
		work_dim1;
#line 356 "zcgesv.f"
	rnrm = (d__1 = work[i__2].r, abs(d__1)) + (d__2 = d_imag(&work[
		izamax_(n, &work[i__ * work_dim1 + 1], &c__1) + i__ * 
		work_dim1]), abs(d__2));
#line 357 "zcgesv.f"
	if (rnrm > xnrm * cte) {
#line 357 "zcgesv.f"
	    goto L10;
#line 357 "zcgesv.f"
	}
#line 359 "zcgesv.f"
    }

/*     If we are here, the NRHS normwise backward errors satisfy the */
/*     stopping criterion. We are good to exit. */

#line 364 "zcgesv.f"
    *iter = 0;
#line 365 "zcgesv.f"
    return 0;

#line 367 "zcgesv.f"
L10:

#line 369 "zcgesv.f"
    for (iiter = 1; iiter <= 30; ++iiter) {

/*        Convert R (in WORK) from double precision to single precision */
/*        and store the result in SX. */

#line 374 "zcgesv.f"
	zlag2c_(n, nrhs, &work[work_offset], n, &swork[ptsx], n, info);

#line 376 "zcgesv.f"
	if (*info != 0) {
#line 377 "zcgesv.f"
	    *iter = -2;
#line 378 "zcgesv.f"
	    goto L40;
#line 379 "zcgesv.f"
	}

/*        Solve the system SA*SX = SR. */

#line 383 "zcgesv.f"
	cgetrs_("No transpose", n, nrhs, &swork[ptsa], n, &ipiv[1], &swork[
		ptsx], n, info, (ftnlen)12);

/*        Convert SX back to double precision and update the current */
/*        iterate. */

#line 389 "zcgesv.f"
	clag2z_(n, nrhs, &swork[ptsx], n, &work[work_offset], n, info);

#line 391 "zcgesv.f"
	i__1 = *nrhs;
#line 391 "zcgesv.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 392 "zcgesv.f"
	    zaxpy_(n, &c_b2, &work[i__ * work_dim1 + 1], &c__1, &x[i__ * 
		    x_dim1 + 1], &c__1);
#line 393 "zcgesv.f"
	}

/*        Compute R = B - AX (R is WORK). */

#line 397 "zcgesv.f"
	zlacpy_("All", n, nrhs, &b[b_offset], ldb, &work[work_offset], n, (
		ftnlen)3);

#line 399 "zcgesv.f"
	zgemm_("No Transpose", "No Transpose", n, nrhs, n, &c_b1, &a[a_offset]
		, lda, &x[x_offset], ldx, &c_b2, &work[work_offset], n, (
		ftnlen)12, (ftnlen)12);

/*        Check whether the NRHS normwise backward errors satisfy the */
/*        stopping criterion. If yes, set ITER=IITER>0 and return. */

#line 405 "zcgesv.f"
	i__1 = *nrhs;
#line 405 "zcgesv.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 406 "zcgesv.f"
	    i__2 = izamax_(n, &x[i__ * x_dim1 + 1], &c__1) + i__ * x_dim1;
#line 406 "zcgesv.f"
	    xnrm = (d__1 = x[i__2].r, abs(d__1)) + (d__2 = d_imag(&x[izamax_(
		    n, &x[i__ * x_dim1 + 1], &c__1) + i__ * x_dim1]), abs(
		    d__2));
#line 407 "zcgesv.f"
	    i__2 = izamax_(n, &work[i__ * work_dim1 + 1], &c__1) + i__ * 
		    work_dim1;
#line 407 "zcgesv.f"
	    rnrm = (d__1 = work[i__2].r, abs(d__1)) + (d__2 = d_imag(&work[
		    izamax_(n, &work[i__ * work_dim1 + 1], &c__1) + i__ * 
		    work_dim1]), abs(d__2));
#line 408 "zcgesv.f"
	    if (rnrm > xnrm * cte) {
#line 408 "zcgesv.f"
		goto L20;
#line 408 "zcgesv.f"
	    }
#line 410 "zcgesv.f"
	}

/*        If we are here, the NRHS normwise backward errors satisfy the */
/*        stopping criterion, we are good to exit. */

#line 415 "zcgesv.f"
	*iter = iiter;

#line 417 "zcgesv.f"
	return 0;

#line 419 "zcgesv.f"
L20:

#line 421 "zcgesv.f"
/* L30: */
#line 421 "zcgesv.f"
	;
#line 421 "zcgesv.f"
    }

/*     If we are at this place of the code, this is because we have */
/*     performed ITER=ITERMAX iterations and never satisified the stopping */
/*     criterion, set up the ITER flag accordingly and follow up on double */
/*     precision routine. */

#line 428 "zcgesv.f"
    *iter = -31;

#line 430 "zcgesv.f"
L40:

/*     Single-precision iterative refinement failed to converge to a */
/*     satisfactory solution, so we resort to double precision. */

#line 435 "zcgesv.f"
    zgetrf_(n, n, &a[a_offset], lda, &ipiv[1], info);

#line 437 "zcgesv.f"
    if (*info != 0) {
#line 437 "zcgesv.f"
	return 0;
#line 437 "zcgesv.f"
    }

#line 440 "zcgesv.f"
    zlacpy_("All", n, nrhs, &b[b_offset], ldb, &x[x_offset], ldx, (ftnlen)3);
#line 441 "zcgesv.f"
    zgetrs_("No transpose", n, nrhs, &a[a_offset], lda, &ipiv[1], &x[x_offset]
	    , ldx, info, (ftnlen)12);

#line 444 "zcgesv.f"
    return 0;

/*     End of ZCGESV. */

} /* zcgesv_ */

