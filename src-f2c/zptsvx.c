#line 1 "zptsvx.f"
/* zptsvx.f -- translated by f2c (version 20100827).
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

#line 1 "zptsvx.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief <b> ZPTSVX computes the solution to system of linear equations A * X = B for PT matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZPTSVX + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zptsvx.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zptsvx.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zptsvx.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZPTSVX( FACT, N, NRHS, D, E, DF, EF, B, LDB, X, LDX, */
/*                          RCOND, FERR, BERR, WORK, RWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          FACT */
/*       INTEGER            INFO, LDB, LDX, N, NRHS */
/*       DOUBLE PRECISION   RCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   BERR( * ), D( * ), DF( * ), FERR( * ), */
/*      $                   RWORK( * ) */
/*       COMPLEX*16         B( LDB, * ), E( * ), EF( * ), WORK( * ), */
/*      $                   X( LDX, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZPTSVX uses the factorization A = L*D*L**H to compute the solution */
/* > to a complex system of linear equations A*X = B, where A is an */
/* > N-by-N Hermitian positive definite tridiagonal matrix and X and B */
/* > are N-by-NRHS matrices. */
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
/* > 1. If FACT = 'N', the matrix A is factored as A = L*D*L**H, where L */
/* >    is a unit lower bidiagonal matrix and D is diagonal.  The */
/* >    factorization can also be regarded as having the form */
/* >    A = U**H*D*U. */
/* > */
/* > 2. If the leading i-by-i principal minor is not positive definite, */
/* >    then the routine returns with INFO = i. Otherwise, the factored */
/* >    form of A is used to estimate the condition number of the matrix */
/* >    A.  If the reciprocal of the condition number is less than machine */
/* >    precision, INFO = N+1 is returned as a warning, but the routine */
/* >    still goes on to solve for X and compute error bounds as */
/* >    described below. */
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
/* >          Specifies whether or not the factored form of the matrix */
/* >          A is supplied on entry. */
/* >          = 'F':  On entry, DF and EF contain the factored form of A. */
/* >                  D, E, DF, and EF will not be modified. */
/* >          = 'N':  The matrix A will be copied to DF and EF and */
/* >                  factored. */
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
/* > \param[in] D */
/* > \verbatim */
/* >          D is DOUBLE PRECISION array, dimension (N) */
/* >          The n diagonal elements of the tridiagonal matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] E */
/* > \verbatim */
/* >          E is COMPLEX*16 array, dimension (N-1) */
/* >          The (n-1) subdiagonal elements of the tridiagonal matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in,out] DF */
/* > \verbatim */
/* >          DF is DOUBLE PRECISION array, dimension (N) */
/* >          If FACT = 'F', then DF is an input argument and on entry */
/* >          contains the n diagonal elements of the diagonal matrix D */
/* >          from the L*D*L**H factorization of A. */
/* >          If FACT = 'N', then DF is an output argument and on exit */
/* >          contains the n diagonal elements of the diagonal matrix D */
/* >          from the L*D*L**H factorization of A. */
/* > \endverbatim */
/* > */
/* > \param[in,out] EF */
/* > \verbatim */
/* >          EF is COMPLEX*16 array, dimension (N-1) */
/* >          If FACT = 'F', then EF is an input argument and on entry */
/* >          contains the (n-1) subdiagonal elements of the unit */
/* >          bidiagonal factor L from the L*D*L**H factorization of A. */
/* >          If FACT = 'N', then EF is an output argument and on exit */
/* >          contains the (n-1) subdiagonal elements of the unit */
/* >          bidiagonal factor L from the L*D*L**H factorization of A. */
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
/* >          RCOND is DOUBLE PRECISION */
/* >          The reciprocal condition number of the matrix A.  If RCOND */
/* >          is less than the machine precision (in particular, if */
/* >          RCOND = 0), the matrix is singular to working precision. */
/* >          This condition is indicated by a return code of INFO > 0. */
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
/* >          vector X(j) (i.e., the smallest relative change in any */
/* >          element of A or B that makes X(j) an exact solution). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX*16 array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is DOUBLE PRECISION array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* >          > 0:  if INFO = i, and i is */
/* >                <= N:  the leading minor of order i of A is */
/* >                       not positive definite, so the factorization */
/* >                       could not be completed, and the solution has not */
/* >                       been computed. RCOND = 0 is returned. */
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

/* > \ingroup complex16PTsolve */

/*  ===================================================================== */
/* Subroutine */ int zptsvx_(char *fact, integer *n, integer *nrhs, 
	doublereal *d__, doublecomplex *e, doublereal *df, doublecomplex *ef, 
	doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx, 
	doublereal *rcond, doublereal *ferr, doublereal *berr, doublecomplex *
	work, doublereal *rwork, integer *info, ftnlen fact_len)
{
    /* System generated locals */
    integer b_dim1, b_offset, x_dim1, x_offset, i__1;

    /* Local variables */
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal anorm;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), zcopy_(integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *);
    extern doublereal dlamch_(char *, ftnlen);
    static logical nofact;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern doublereal zlanht_(char *, integer *, doublereal *, doublecomplex *
	    , ftnlen);
    extern /* Subroutine */ int zlacpy_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen), 
	    zptcon_(integer *, doublereal *, doublecomplex *, doublereal *, 
	    doublereal *, doublereal *, integer *), zptrfs_(char *, integer *,
	     integer *, doublereal *, doublecomplex *, doublereal *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublereal *, doublereal *, doublecomplex *, 
	    doublereal *, integer *, ftnlen), zpttrf_(integer *, doublereal *,
	     doublecomplex *, integer *), zpttrs_(char *, integer *, integer *
	    , doublereal *, doublecomplex *, doublecomplex *, integer *, 
	    integer *, ftnlen);


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

/*     Test the input parameters. */

#line 280 "zptsvx.f"
    /* Parameter adjustments */
#line 280 "zptsvx.f"
    --d__;
#line 280 "zptsvx.f"
    --e;
#line 280 "zptsvx.f"
    --df;
#line 280 "zptsvx.f"
    --ef;
#line 280 "zptsvx.f"
    b_dim1 = *ldb;
#line 280 "zptsvx.f"
    b_offset = 1 + b_dim1;
#line 280 "zptsvx.f"
    b -= b_offset;
#line 280 "zptsvx.f"
    x_dim1 = *ldx;
#line 280 "zptsvx.f"
    x_offset = 1 + x_dim1;
#line 280 "zptsvx.f"
    x -= x_offset;
#line 280 "zptsvx.f"
    --ferr;
#line 280 "zptsvx.f"
    --berr;
#line 280 "zptsvx.f"
    --work;
#line 280 "zptsvx.f"
    --rwork;
#line 280 "zptsvx.f"

#line 280 "zptsvx.f"
    /* Function Body */
#line 280 "zptsvx.f"
    *info = 0;
#line 281 "zptsvx.f"
    nofact = lsame_(fact, "N", (ftnlen)1, (ftnlen)1);
#line 282 "zptsvx.f"
    if (! nofact && ! lsame_(fact, "F", (ftnlen)1, (ftnlen)1)) {
#line 283 "zptsvx.f"
	*info = -1;
#line 284 "zptsvx.f"
    } else if (*n < 0) {
#line 285 "zptsvx.f"
	*info = -2;
#line 286 "zptsvx.f"
    } else if (*nrhs < 0) {
#line 287 "zptsvx.f"
	*info = -3;
#line 288 "zptsvx.f"
    } else if (*ldb < max(1,*n)) {
#line 289 "zptsvx.f"
	*info = -9;
#line 290 "zptsvx.f"
    } else if (*ldx < max(1,*n)) {
#line 291 "zptsvx.f"
	*info = -11;
#line 292 "zptsvx.f"
    }
#line 293 "zptsvx.f"
    if (*info != 0) {
#line 294 "zptsvx.f"
	i__1 = -(*info);
#line 294 "zptsvx.f"
	xerbla_("ZPTSVX", &i__1, (ftnlen)6);
#line 295 "zptsvx.f"
	return 0;
#line 296 "zptsvx.f"
    }

#line 298 "zptsvx.f"
    if (nofact) {

/*        Compute the L*D*L**H (or U**H*D*U) factorization of A. */

#line 302 "zptsvx.f"
	dcopy_(n, &d__[1], &c__1, &df[1], &c__1);
#line 303 "zptsvx.f"
	if (*n > 1) {
#line 303 "zptsvx.f"
	    i__1 = *n - 1;
#line 303 "zptsvx.f"
	    zcopy_(&i__1, &e[1], &c__1, &ef[1], &c__1);
#line 303 "zptsvx.f"
	}
#line 305 "zptsvx.f"
	zpttrf_(n, &df[1], &ef[1], info);

/*        Return if INFO is non-zero. */

#line 309 "zptsvx.f"
	if (*info > 0) {
#line 310 "zptsvx.f"
	    *rcond = 0.;
#line 311 "zptsvx.f"
	    return 0;
#line 312 "zptsvx.f"
	}
#line 313 "zptsvx.f"
    }

/*     Compute the norm of the matrix A. */

#line 317 "zptsvx.f"
    anorm = zlanht_("1", n, &d__[1], &e[1], (ftnlen)1);

/*     Compute the reciprocal of the condition number of A. */

#line 321 "zptsvx.f"
    zptcon_(n, &df[1], &ef[1], &anorm, rcond, &rwork[1], info);

/*     Compute the solution vectors X. */

#line 325 "zptsvx.f"
    zlacpy_("Full", n, nrhs, &b[b_offset], ldb, &x[x_offset], ldx, (ftnlen)4);
#line 326 "zptsvx.f"
    zpttrs_("Lower", n, nrhs, &df[1], &ef[1], &x[x_offset], ldx, info, (
	    ftnlen)5);

/*     Use iterative refinement to improve the computed solutions and */
/*     compute error bounds and backward error estimates for them. */

#line 331 "zptsvx.f"
    zptrfs_("Lower", n, nrhs, &d__[1], &e[1], &df[1], &ef[1], &b[b_offset], 
	    ldb, &x[x_offset], ldx, &ferr[1], &berr[1], &work[1], &rwork[1], 
	    info, (ftnlen)5);

/*     Set INFO = N+1 if the matrix is singular to working precision. */

#line 336 "zptsvx.f"
    if (*rcond < dlamch_("Epsilon", (ftnlen)7)) {
#line 336 "zptsvx.f"
	*info = *n + 1;
#line 336 "zptsvx.f"
    }

#line 339 "zptsvx.f"
    return 0;

/*     End of ZPTSVX */

} /* zptsvx_ */

