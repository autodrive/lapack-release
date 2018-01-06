#line 1 "sgtsvx.f"
/* sgtsvx.f -- translated by f2c (version 20100827).
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

#line 1 "sgtsvx.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief <b> SGTSVX computes the solution to system of linear equations A * X = B for GT matrices <b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SGTSVX + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgtsvx.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgtsvx.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgtsvx.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SGTSVX( FACT, TRANS, N, NRHS, DL, D, DU, DLF, DF, DUF, */
/*                          DU2, IPIV, B, LDB, X, LDX, RCOND, FERR, BERR, */
/*                          WORK, IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          FACT, TRANS */
/*       INTEGER            INFO, LDB, LDX, N, NRHS */
/*       REAL               RCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ), IWORK( * ) */
/*       REAL               B( LDB, * ), BERR( * ), D( * ), DF( * ), */
/*      $                   DL( * ), DLF( * ), DU( * ), DU2( * ), DUF( * ), */
/*      $                   FERR( * ), WORK( * ), X( LDX, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SGTSVX uses the LU factorization to compute the solution to a real */
/* > system of linear equations A * X = B or A**T * X = B, */
/* > where A is a tridiagonal matrix of order N and X and B are N-by-NRHS */
/* > matrices. */
/* > */
/* > Error bounds on the solution and a condition estimate are also */
/* > provided. */
/* > \endverbatim */

/* > \par Description: */
/*  ================= */
/* > */
/* > \verbatim */
/* > */
/* > The following steps are performed: */
/* > */
/* > 1. If FACT = 'N', the LU decomposition is used to factor the matrix A */
/* >    as A = L * U, where L is a product of permutation and unit lower */
/* >    bidiagonal matrices and U is upper triangular with nonzeros in */
/* >    only the main diagonal and first two superdiagonals. */
/* > */
/* > 2. If some U(i,i)=0, so that U is exactly singular, then the routine */
/* >    returns with INFO = i. Otherwise, the factored form of A is used */
/* >    to estimate the condition number of the matrix A.  If the */
/* >    reciprocal of the condition number is less than machine precision, */
/* >    INFO = N+1 is returned as a warning, but the routine still goes on */
/* >    to solve for X and compute error bounds as described below. */
/* > */
/* > 3. The system of equations is solved for X using the factored form */
/* >    of A. */
/* > */
/* > 4. Iterative refinement is applied to improve the computed solution */
/* >    matrix and calculate error bounds and backward error estimates */
/* >    for it. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] FACT */
/* > \verbatim */
/* >          FACT is CHARACTER*1 */
/* >          Specifies whether or not the factored form of A has been */
/* >          supplied on entry. */
/* >          = 'F':  DLF, DF, DUF, DU2, and IPIV contain the factored */
/* >                  form of A; DL, D, DU, DLF, DF, DUF, DU2 and IPIV */
/* >                  will not be modified. */
/* >          = 'N':  The matrix will be copied to DLF, DF, and DUF */
/* >                  and factored. */
/* > \endverbatim */
/* > */
/* > \param[in] TRANS */
/* > \verbatim */
/* >          TRANS is CHARACTER*1 */
/* >          Specifies the form of the system of equations: */
/* >          = 'N':  A * X = B     (No transpose) */
/* >          = 'T':  A**T * X = B  (Transpose) */
/* >          = 'C':  A**H * X = B  (Conjugate transpose = Transpose) */
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
/* >          of the matrix B.  NRHS >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] DL */
/* > \verbatim */
/* >          DL is REAL array, dimension (N-1) */
/* >          The (n-1) subdiagonal elements of A. */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* >          D is REAL array, dimension (N) */
/* >          The n diagonal elements of A. */
/* > \endverbatim */
/* > */
/* > \param[in] DU */
/* > \verbatim */
/* >          DU is REAL array, dimension (N-1) */
/* >          The (n-1) superdiagonal elements of A. */
/* > \endverbatim */
/* > */
/* > \param[in,out] DLF */
/* > \verbatim */
/* >          DLF is REAL array, dimension (N-1) */
/* >          If FACT = 'F', then DLF is an input argument and on entry */
/* >          contains the (n-1) multipliers that define the matrix L from */
/* >          the LU factorization of A as computed by SGTTRF. */
/* > */
/* >          If FACT = 'N', then DLF is an output argument and on exit */
/* >          contains the (n-1) multipliers that define the matrix L from */
/* >          the LU factorization of A. */
/* > \endverbatim */
/* > */
/* > \param[in,out] DF */
/* > \verbatim */
/* >          DF is REAL array, dimension (N) */
/* >          If FACT = 'F', then DF is an input argument and on entry */
/* >          contains the n diagonal elements of the upper triangular */
/* >          matrix U from the LU factorization of A. */
/* > */
/* >          If FACT = 'N', then DF is an output argument and on exit */
/* >          contains the n diagonal elements of the upper triangular */
/* >          matrix U from the LU factorization of A. */
/* > \endverbatim */
/* > */
/* > \param[in,out] DUF */
/* > \verbatim */
/* >          DUF is REAL array, dimension (N-1) */
/* >          If FACT = 'F', then DUF is an input argument and on entry */
/* >          contains the (n-1) elements of the first superdiagonal of U. */
/* > */
/* >          If FACT = 'N', then DUF is an output argument and on exit */
/* >          contains the (n-1) elements of the first superdiagonal of U. */
/* > \endverbatim */
/* > */
/* > \param[in,out] DU2 */
/* > \verbatim */
/* >          DU2 is REAL array, dimension (N-2) */
/* >          If FACT = 'F', then DU2 is an input argument and on entry */
/* >          contains the (n-2) elements of the second superdiagonal of */
/* >          U. */
/* > */
/* >          If FACT = 'N', then DU2 is an output argument and on exit */
/* >          contains the (n-2) elements of the second superdiagonal of */
/* >          U. */
/* > \endverbatim */
/* > */
/* > \param[in,out] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N) */
/* >          If FACT = 'F', then IPIV is an input argument and on entry */
/* >          contains the pivot indices from the LU factorization of A as */
/* >          computed by SGTTRF. */
/* > */
/* >          If FACT = 'N', then IPIV is an output argument and on exit */
/* >          contains the pivot indices from the LU factorization of A; */
/* >          row i of the matrix was interchanged with row IPIV(i). */
/* >          IPIV(i) will always be either i or i+1; IPIV(i) = i indicates */
/* >          a row interchange was not required. */
/* > \endverbatim */
/* > */
/* > \param[in] B */
/* > \verbatim */
/* >          B is REAL array, dimension (LDB,NRHS) */
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
/* >          X is REAL array, dimension (LDX,NRHS) */
/* >          If INFO = 0 or INFO = N+1, the N-by-NRHS solution matrix X. */
/* > \endverbatim */
/* > */
/* > \param[in] LDX */
/* > \verbatim */
/* >          LDX is INTEGER */
/* >          The leading dimension of the array X.  LDX >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] RCOND */
/* > \verbatim */
/* >          RCOND is REAL */
/* >          The estimate of the reciprocal condition number of the matrix */
/* >          A.  If RCOND is less than the machine precision (in */
/* >          particular, if RCOND = 0), the matrix is singular to working */
/* >          precision.  This condition is indicated by a return code of */
/* >          INFO > 0. */
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
/* >          > 0:  if INFO = i, and i is */
/* >                <= N:  U(i,i) is exactly zero.  The factorization */
/* >                       has not been completed unless i = N, but the */
/* >                       factor U is exactly singular, so the solution */
/* >                       and error bounds could not be computed. */
/* >                       RCOND = 0 is returned. */
/* >                = N+1: U is nonsingular, but RCOND is less than machine */
/* >                       precision, meaning that the matrix is singular */
/* >                       to working precision.  Nevertheless, the */
/* >                       solution and error bounds are computed because */
/* >                       there are a number of situations where the */
/* >                       computed solution can be more accurate than the */
/* >                       value of RCOND would suggest. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup realGTsolve */

/*  ===================================================================== */
/* Subroutine */ int sgtsvx_(char *fact, char *trans, integer *n, integer *
	nrhs, doublereal *dl, doublereal *d__, doublereal *du, doublereal *
	dlf, doublereal *df, doublereal *duf, doublereal *du2, integer *ipiv, 
	doublereal *b, integer *ldb, doublereal *x, integer *ldx, doublereal *
	rcond, doublereal *ferr, doublereal *berr, doublereal *work, integer *
	iwork, integer *info, ftnlen fact_len, ftnlen trans_len)
{
    /* System generated locals */
    integer b_dim1, b_offset, x_dim1, x_offset, i__1;

    /* Local variables */
    static char norm[1];
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal anorm;
    extern /* Subroutine */ int scopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    extern doublereal slamch_(char *, ftnlen);
    static logical nofact;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern doublereal slangt_(char *, integer *, doublereal *, doublereal *, 
	    doublereal *, ftnlen);
    extern /* Subroutine */ int slacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    sgtcon_(char *, integer *, doublereal *, doublereal *, doublereal 
	    *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, ftnlen);
    static logical notran;
    extern /* Subroutine */ int sgtrfs_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, ftnlen), sgttrf_(integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *), sgttrs_(char *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     integer *, integer *, ftnlen);


/*  -- LAPACK driver routine (version 3.4.2) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     September 2012 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 338 "sgtsvx.f"
    /* Parameter adjustments */
#line 338 "sgtsvx.f"
    --dl;
#line 338 "sgtsvx.f"
    --d__;
#line 338 "sgtsvx.f"
    --du;
#line 338 "sgtsvx.f"
    --dlf;
#line 338 "sgtsvx.f"
    --df;
#line 338 "sgtsvx.f"
    --duf;
#line 338 "sgtsvx.f"
    --du2;
#line 338 "sgtsvx.f"
    --ipiv;
#line 338 "sgtsvx.f"
    b_dim1 = *ldb;
#line 338 "sgtsvx.f"
    b_offset = 1 + b_dim1;
#line 338 "sgtsvx.f"
    b -= b_offset;
#line 338 "sgtsvx.f"
    x_dim1 = *ldx;
#line 338 "sgtsvx.f"
    x_offset = 1 + x_dim1;
#line 338 "sgtsvx.f"
    x -= x_offset;
#line 338 "sgtsvx.f"
    --ferr;
#line 338 "sgtsvx.f"
    --berr;
#line 338 "sgtsvx.f"
    --work;
#line 338 "sgtsvx.f"
    --iwork;
#line 338 "sgtsvx.f"

#line 338 "sgtsvx.f"
    /* Function Body */
#line 338 "sgtsvx.f"
    *info = 0;
#line 339 "sgtsvx.f"
    nofact = lsame_(fact, "N", (ftnlen)1, (ftnlen)1);
#line 340 "sgtsvx.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 341 "sgtsvx.f"
    if (! nofact && ! lsame_(fact, "F", (ftnlen)1, (ftnlen)1)) {
#line 342 "sgtsvx.f"
	*info = -1;
#line 343 "sgtsvx.f"
    } else if (! notran && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1) && ! 
	    lsame_(trans, "C", (ftnlen)1, (ftnlen)1)) {
#line 345 "sgtsvx.f"
	*info = -2;
#line 346 "sgtsvx.f"
    } else if (*n < 0) {
#line 347 "sgtsvx.f"
	*info = -3;
#line 348 "sgtsvx.f"
    } else if (*nrhs < 0) {
#line 349 "sgtsvx.f"
	*info = -4;
#line 350 "sgtsvx.f"
    } else if (*ldb < max(1,*n)) {
#line 351 "sgtsvx.f"
	*info = -14;
#line 352 "sgtsvx.f"
    } else if (*ldx < max(1,*n)) {
#line 353 "sgtsvx.f"
	*info = -16;
#line 354 "sgtsvx.f"
    }
#line 355 "sgtsvx.f"
    if (*info != 0) {
#line 356 "sgtsvx.f"
	i__1 = -(*info);
#line 356 "sgtsvx.f"
	xerbla_("SGTSVX", &i__1, (ftnlen)6);
#line 357 "sgtsvx.f"
	return 0;
#line 358 "sgtsvx.f"
    }

#line 360 "sgtsvx.f"
    if (nofact) {

/*        Compute the LU factorization of A. */

#line 364 "sgtsvx.f"
	scopy_(n, &d__[1], &c__1, &df[1], &c__1);
#line 365 "sgtsvx.f"
	if (*n > 1) {
#line 366 "sgtsvx.f"
	    i__1 = *n - 1;
#line 366 "sgtsvx.f"
	    scopy_(&i__1, &dl[1], &c__1, &dlf[1], &c__1);
#line 367 "sgtsvx.f"
	    i__1 = *n - 1;
#line 367 "sgtsvx.f"
	    scopy_(&i__1, &du[1], &c__1, &duf[1], &c__1);
#line 368 "sgtsvx.f"
	}
#line 369 "sgtsvx.f"
	sgttrf_(n, &dlf[1], &df[1], &duf[1], &du2[1], &ipiv[1], info);

/*        Return if INFO is non-zero. */

#line 373 "sgtsvx.f"
	if (*info > 0) {
#line 374 "sgtsvx.f"
	    *rcond = 0.;
#line 375 "sgtsvx.f"
	    return 0;
#line 376 "sgtsvx.f"
	}
#line 377 "sgtsvx.f"
    }

/*     Compute the norm of the matrix A. */

#line 381 "sgtsvx.f"
    if (notran) {
#line 382 "sgtsvx.f"
	*(unsigned char *)norm = '1';
#line 383 "sgtsvx.f"
    } else {
#line 384 "sgtsvx.f"
	*(unsigned char *)norm = 'I';
#line 385 "sgtsvx.f"
    }
#line 386 "sgtsvx.f"
    anorm = slangt_(norm, n, &dl[1], &d__[1], &du[1], (ftnlen)1);

/*     Compute the reciprocal of the condition number of A. */

#line 390 "sgtsvx.f"
    sgtcon_(norm, n, &dlf[1], &df[1], &duf[1], &du2[1], &ipiv[1], &anorm, 
	    rcond, &work[1], &iwork[1], info, (ftnlen)1);

/*     Compute the solution vectors X. */

#line 395 "sgtsvx.f"
    slacpy_("Full", n, nrhs, &b[b_offset], ldb, &x[x_offset], ldx, (ftnlen)4);
#line 396 "sgtsvx.f"
    sgttrs_(trans, n, nrhs, &dlf[1], &df[1], &duf[1], &du2[1], &ipiv[1], &x[
	    x_offset], ldx, info, (ftnlen)1);

/*     Use iterative refinement to improve the computed solutions and */
/*     compute error bounds and backward error estimates for them. */

#line 402 "sgtsvx.f"
    sgtrfs_(trans, n, nrhs, &dl[1], &d__[1], &du[1], &dlf[1], &df[1], &duf[1],
	     &du2[1], &ipiv[1], &b[b_offset], ldb, &x[x_offset], ldx, &ferr[1]
	    , &berr[1], &work[1], &iwork[1], info, (ftnlen)1);

/*     Set INFO = N+1 if the matrix is singular to working precision. */

#line 407 "sgtsvx.f"
    if (*rcond < slamch_("Epsilon", (ftnlen)7)) {
#line 407 "sgtsvx.f"
	*info = *n + 1;
#line 407 "sgtsvx.f"
    }

#line 410 "sgtsvx.f"
    return 0;

/*     End of SGTSVX */

} /* sgtsvx_ */

