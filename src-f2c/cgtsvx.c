#line 1 "cgtsvx.f"
/* cgtsvx.f -- translated by f2c (version 20100827).
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

#line 1 "cgtsvx.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief <b> CGTSVX computes the solution to system of linear equations A * X = B for GT matrices <b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CGTSVX + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgtsvx.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgtsvx.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgtsvx.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CGTSVX( FACT, TRANS, N, NRHS, DL, D, DU, DLF, DF, DUF, */
/*                          DU2, IPIV, B, LDB, X, LDX, RCOND, FERR, BERR, */
/*                          WORK, RWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          FACT, TRANS */
/*       INTEGER            INFO, LDB, LDX, N, NRHS */
/*       REAL               RCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       REAL               BERR( * ), FERR( * ), RWORK( * ) */
/*       COMPLEX            B( LDB, * ), D( * ), DF( * ), DL( * ), */
/*      $                   DLF( * ), DU( * ), DU2( * ), DUF( * ), */
/*      $                   WORK( * ), X( LDX, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CGTSVX uses the LU factorization to compute the solution to a complex */
/* > system of linear equations A * X = B, A**T * X = B, or A**H * X = B, */
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
/* >          = 'F':  DLF, DF, DUF, DU2, and IPIV contain the factored form */
/* >                  of A; DL, D, DU, DLF, DF, DUF, DU2 and IPIV will not */
/* >                  be modified. */
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
/* >          = 'C':  A**H * X = B  (Conjugate transpose) */
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
/* >          DL is COMPLEX array, dimension (N-1) */
/* >          The (n-1) subdiagonal elements of A. */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* >          D is COMPLEX array, dimension (N) */
/* >          The n diagonal elements of A. */
/* > \endverbatim */
/* > */
/* > \param[in] DU */
/* > \verbatim */
/* >          DU is COMPLEX array, dimension (N-1) */
/* >          The (n-1) superdiagonal elements of A. */
/* > \endverbatim */
/* > */
/* > \param[in,out] DLF */
/* > \verbatim */
/* >          DLF is COMPLEX array, dimension (N-1) */
/* >          If FACT = 'F', then DLF is an input argument and on entry */
/* >          contains the (n-1) multipliers that define the matrix L from */
/* >          the LU factorization of A as computed by CGTTRF. */
/* > */
/* >          If FACT = 'N', then DLF is an output argument and on exit */
/* >          contains the (n-1) multipliers that define the matrix L from */
/* >          the LU factorization of A. */
/* > \endverbatim */
/* > */
/* > \param[in,out] DF */
/* > \verbatim */
/* >          DF is COMPLEX array, dimension (N) */
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
/* >          DUF is COMPLEX array, dimension (N-1) */
/* >          If FACT = 'F', then DUF is an input argument and on entry */
/* >          contains the (n-1) elements of the first superdiagonal of U. */
/* > */
/* >          If FACT = 'N', then DUF is an output argument and on exit */
/* >          contains the (n-1) elements of the first superdiagonal of U. */
/* > \endverbatim */
/* > */
/* > \param[in,out] DU2 */
/* > \verbatim */
/* >          DU2 is COMPLEX array, dimension (N-2) */
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
/* >          computed by CGTTRF. */
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
/* >          B is COMPLEX array, dimension (LDB,NRHS) */
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
/* >          X is COMPLEX array, dimension (LDX,NRHS) */
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
/* >          WORK is COMPLEX array, dimension (2*N) */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is REAL array, dimension (N) */
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

/* > \ingroup complexGTsolve */

/*  ===================================================================== */
/* Subroutine */ int cgtsvx_(char *fact, char *trans, integer *n, integer *
	nrhs, doublecomplex *dl, doublecomplex *d__, doublecomplex *du, 
	doublecomplex *dlf, doublecomplex *df, doublecomplex *duf, 
	doublecomplex *du2, integer *ipiv, doublecomplex *b, integer *ldb, 
	doublecomplex *x, integer *ldx, doublereal *rcond, doublereal *ferr, 
	doublereal *berr, doublecomplex *work, doublereal *rwork, integer *
	info, ftnlen fact_len, ftnlen trans_len)
{
    /* System generated locals */
    integer b_dim1, b_offset, x_dim1, x_offset, i__1;

    /* Local variables */
    static char norm[1];
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal anorm;
    extern /* Subroutine */ int ccopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    extern doublereal slamch_(char *, ftnlen), clangt_(char *, integer *, 
	    doublecomplex *, doublecomplex *, doublecomplex *, ftnlen);
    static logical nofact;
    extern /* Subroutine */ int clacpy_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen), 
	    cgtcon_(char *, integer *, doublecomplex *, doublecomplex *, 
	    doublecomplex *, doublecomplex *, integer *, doublereal *, 
	    doublereal *, doublecomplex *, integer *, ftnlen), xerbla_(char *,
	     integer *, ftnlen), cgtrfs_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, doublecomplex *, doublecomplex *
	    , doublecomplex *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublereal *, doublereal *, doublecomplex *, doublereal *, 
	    integer *, ftnlen), cgttrf_(integer *, doublecomplex *, 
	    doublecomplex *, doublecomplex *, doublecomplex *, integer *, 
	    integer *);
    static logical notran;
    extern /* Subroutine */ int cgttrs_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, doublecomplex *, doublecomplex *
	    , integer *, doublecomplex *, integer *, integer *, ftnlen);


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

#line 340 "cgtsvx.f"
    /* Parameter adjustments */
#line 340 "cgtsvx.f"
    --dl;
#line 340 "cgtsvx.f"
    --d__;
#line 340 "cgtsvx.f"
    --du;
#line 340 "cgtsvx.f"
    --dlf;
#line 340 "cgtsvx.f"
    --df;
#line 340 "cgtsvx.f"
    --duf;
#line 340 "cgtsvx.f"
    --du2;
#line 340 "cgtsvx.f"
    --ipiv;
#line 340 "cgtsvx.f"
    b_dim1 = *ldb;
#line 340 "cgtsvx.f"
    b_offset = 1 + b_dim1;
#line 340 "cgtsvx.f"
    b -= b_offset;
#line 340 "cgtsvx.f"
    x_dim1 = *ldx;
#line 340 "cgtsvx.f"
    x_offset = 1 + x_dim1;
#line 340 "cgtsvx.f"
    x -= x_offset;
#line 340 "cgtsvx.f"
    --ferr;
#line 340 "cgtsvx.f"
    --berr;
#line 340 "cgtsvx.f"
    --work;
#line 340 "cgtsvx.f"
    --rwork;
#line 340 "cgtsvx.f"

#line 340 "cgtsvx.f"
    /* Function Body */
#line 340 "cgtsvx.f"
    *info = 0;
#line 341 "cgtsvx.f"
    nofact = lsame_(fact, "N", (ftnlen)1, (ftnlen)1);
#line 342 "cgtsvx.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 343 "cgtsvx.f"
    if (! nofact && ! lsame_(fact, "F", (ftnlen)1, (ftnlen)1)) {
#line 344 "cgtsvx.f"
	*info = -1;
#line 345 "cgtsvx.f"
    } else if (! notran && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1) && ! 
	    lsame_(trans, "C", (ftnlen)1, (ftnlen)1)) {
#line 347 "cgtsvx.f"
	*info = -2;
#line 348 "cgtsvx.f"
    } else if (*n < 0) {
#line 349 "cgtsvx.f"
	*info = -3;
#line 350 "cgtsvx.f"
    } else if (*nrhs < 0) {
#line 351 "cgtsvx.f"
	*info = -4;
#line 352 "cgtsvx.f"
    } else if (*ldb < max(1,*n)) {
#line 353 "cgtsvx.f"
	*info = -14;
#line 354 "cgtsvx.f"
    } else if (*ldx < max(1,*n)) {
#line 355 "cgtsvx.f"
	*info = -16;
#line 356 "cgtsvx.f"
    }
#line 357 "cgtsvx.f"
    if (*info != 0) {
#line 358 "cgtsvx.f"
	i__1 = -(*info);
#line 358 "cgtsvx.f"
	xerbla_("CGTSVX", &i__1, (ftnlen)6);
#line 359 "cgtsvx.f"
	return 0;
#line 360 "cgtsvx.f"
    }

#line 362 "cgtsvx.f"
    if (nofact) {

/*        Compute the LU factorization of A. */

#line 366 "cgtsvx.f"
	ccopy_(n, &d__[1], &c__1, &df[1], &c__1);
#line 367 "cgtsvx.f"
	if (*n > 1) {
#line 368 "cgtsvx.f"
	    i__1 = *n - 1;
#line 368 "cgtsvx.f"
	    ccopy_(&i__1, &dl[1], &c__1, &dlf[1], &c__1);
#line 369 "cgtsvx.f"
	    i__1 = *n - 1;
#line 369 "cgtsvx.f"
	    ccopy_(&i__1, &du[1], &c__1, &duf[1], &c__1);
#line 370 "cgtsvx.f"
	}
#line 371 "cgtsvx.f"
	cgttrf_(n, &dlf[1], &df[1], &duf[1], &du2[1], &ipiv[1], info);

/*        Return if INFO is non-zero. */

#line 375 "cgtsvx.f"
	if (*info > 0) {
#line 376 "cgtsvx.f"
	    *rcond = 0.;
#line 377 "cgtsvx.f"
	    return 0;
#line 378 "cgtsvx.f"
	}
#line 379 "cgtsvx.f"
    }

/*     Compute the norm of the matrix A. */

#line 383 "cgtsvx.f"
    if (notran) {
#line 384 "cgtsvx.f"
	*(unsigned char *)norm = '1';
#line 385 "cgtsvx.f"
    } else {
#line 386 "cgtsvx.f"
	*(unsigned char *)norm = 'I';
#line 387 "cgtsvx.f"
    }
#line 388 "cgtsvx.f"
    anorm = clangt_(norm, n, &dl[1], &d__[1], &du[1], (ftnlen)1);

/*     Compute the reciprocal of the condition number of A. */

#line 392 "cgtsvx.f"
    cgtcon_(norm, n, &dlf[1], &df[1], &duf[1], &du2[1], &ipiv[1], &anorm, 
	    rcond, &work[1], info, (ftnlen)1);

/*     Compute the solution vectors X. */

#line 397 "cgtsvx.f"
    clacpy_("Full", n, nrhs, &b[b_offset], ldb, &x[x_offset], ldx, (ftnlen)4);
#line 398 "cgtsvx.f"
    cgttrs_(trans, n, nrhs, &dlf[1], &df[1], &duf[1], &du2[1], &ipiv[1], &x[
	    x_offset], ldx, info, (ftnlen)1);

/*     Use iterative refinement to improve the computed solutions and */
/*     compute error bounds and backward error estimates for them. */

#line 404 "cgtsvx.f"
    cgtrfs_(trans, n, nrhs, &dl[1], &d__[1], &du[1], &dlf[1], &df[1], &duf[1],
	     &du2[1], &ipiv[1], &b[b_offset], ldb, &x[x_offset], ldx, &ferr[1]
	    , &berr[1], &work[1], &rwork[1], info, (ftnlen)1);

/*     Set INFO = N+1 if the matrix is singular to working precision. */

#line 409 "cgtsvx.f"
    if (*rcond < slamch_("Epsilon", (ftnlen)7)) {
#line 409 "cgtsvx.f"
	*info = *n + 1;
#line 409 "cgtsvx.f"
    }

#line 412 "cgtsvx.f"
    return 0;

/*     End of CGTSVX */

} /* cgtsvx_ */

