#line 1 "cptsvx.f"
/* cptsvx.f -- translated by f2c (version 20100827).
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

#line 1 "cptsvx.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief <b> CPTSVX computes the solution to system of linear equations A * X = B for PT matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CPTSVX + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cptsvx.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cptsvx.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cptsvx.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CPTSVX( FACT, N, NRHS, D, E, DF, EF, B, LDB, X, LDX, */
/*                          RCOND, FERR, BERR, WORK, RWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          FACT */
/*       INTEGER            INFO, LDB, LDX, N, NRHS */
/*       REAL               RCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               BERR( * ), D( * ), DF( * ), FERR( * ), */
/*      $                   RWORK( * ) */
/*       COMPLEX            B( LDB, * ), E( * ), EF( * ), WORK( * ), */
/*      $                   X( LDX, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CPTSVX uses the factorization A = L*D*L**H to compute the solution */
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
/* >          D is REAL array, dimension (N) */
/* >          The n diagonal elements of the tridiagonal matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] E */
/* > \verbatim */
/* >          E is COMPLEX array, dimension (N-1) */
/* >          The (n-1) subdiagonal elements of the tridiagonal matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in,out] DF */
/* > \verbatim */
/* >          DF is REAL array, dimension (N) */
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
/* >          EF is COMPLEX array, dimension (N-1) */
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
/* >          The reciprocal condition number of the matrix A.  If RCOND */
/* >          is less than the machine precision (in particular, if */
/* >          RCOND = 0), the matrix is singular to working precision. */
/* >          This condition is indicated by a return code of INFO > 0. */
/* > \endverbatim */
/* > */
/* > \param[out] FERR */
/* > \verbatim */
/* >          FERR is REAL array, dimension (NRHS) */
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
/* >          BERR is REAL array, dimension (NRHS) */
/* >          The componentwise relative backward error of each solution */
/* >          vector X(j) (i.e., the smallest relative change in any */
/* >          element of A or B that makes X(j) an exact solution). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX array, dimension (N) */
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

/* > \ingroup complexPTsolve */

/*  ===================================================================== */
/* Subroutine */ int cptsvx_(char *fact, integer *n, integer *nrhs, 
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
    extern /* Subroutine */ int ccopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), scopy_(integer *, doublereal *, 
	    integer *, doublereal *, integer *);
    extern doublereal slamch_(char *, ftnlen), clanht_(char *, integer *, 
	    doublereal *, doublecomplex *, ftnlen);
    static logical nofact;
    extern /* Subroutine */ int clacpy_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen), cptcon_(integer *, doublereal 
	    *, doublecomplex *, doublereal *, doublereal *, doublereal *, 
	    integer *), cptrfs_(char *, integer *, integer *, doublereal *, 
	    doublecomplex *, doublereal *, doublecomplex *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, doublereal *, doublereal *,
	     doublecomplex *, doublereal *, integer *, ftnlen), cpttrf_(
	    integer *, doublereal *, doublecomplex *, integer *), cpttrs_(
	    char *, integer *, integer *, doublereal *, doublecomplex *, 
	    doublecomplex *, integer *, integer *, ftnlen);


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

#line 280 "cptsvx.f"
    /* Parameter adjustments */
#line 280 "cptsvx.f"
    --d__;
#line 280 "cptsvx.f"
    --e;
#line 280 "cptsvx.f"
    --df;
#line 280 "cptsvx.f"
    --ef;
#line 280 "cptsvx.f"
    b_dim1 = *ldb;
#line 280 "cptsvx.f"
    b_offset = 1 + b_dim1;
#line 280 "cptsvx.f"
    b -= b_offset;
#line 280 "cptsvx.f"
    x_dim1 = *ldx;
#line 280 "cptsvx.f"
    x_offset = 1 + x_dim1;
#line 280 "cptsvx.f"
    x -= x_offset;
#line 280 "cptsvx.f"
    --ferr;
#line 280 "cptsvx.f"
    --berr;
#line 280 "cptsvx.f"
    --work;
#line 280 "cptsvx.f"
    --rwork;
#line 280 "cptsvx.f"

#line 280 "cptsvx.f"
    /* Function Body */
#line 280 "cptsvx.f"
    *info = 0;
#line 281 "cptsvx.f"
    nofact = lsame_(fact, "N", (ftnlen)1, (ftnlen)1);
#line 282 "cptsvx.f"
    if (! nofact && ! lsame_(fact, "F", (ftnlen)1, (ftnlen)1)) {
#line 283 "cptsvx.f"
	*info = -1;
#line 284 "cptsvx.f"
    } else if (*n < 0) {
#line 285 "cptsvx.f"
	*info = -2;
#line 286 "cptsvx.f"
    } else if (*nrhs < 0) {
#line 287 "cptsvx.f"
	*info = -3;
#line 288 "cptsvx.f"
    } else if (*ldb < max(1,*n)) {
#line 289 "cptsvx.f"
	*info = -9;
#line 290 "cptsvx.f"
    } else if (*ldx < max(1,*n)) {
#line 291 "cptsvx.f"
	*info = -11;
#line 292 "cptsvx.f"
    }
#line 293 "cptsvx.f"
    if (*info != 0) {
#line 294 "cptsvx.f"
	i__1 = -(*info);
#line 294 "cptsvx.f"
	xerbla_("CPTSVX", &i__1, (ftnlen)6);
#line 295 "cptsvx.f"
	return 0;
#line 296 "cptsvx.f"
    }

#line 298 "cptsvx.f"
    if (nofact) {

/*        Compute the L*D*L**H (or U**H*D*U) factorization of A. */

#line 302 "cptsvx.f"
	scopy_(n, &d__[1], &c__1, &df[1], &c__1);
#line 303 "cptsvx.f"
	if (*n > 1) {
#line 303 "cptsvx.f"
	    i__1 = *n - 1;
#line 303 "cptsvx.f"
	    ccopy_(&i__1, &e[1], &c__1, &ef[1], &c__1);
#line 303 "cptsvx.f"
	}
#line 305 "cptsvx.f"
	cpttrf_(n, &df[1], &ef[1], info);

/*        Return if INFO is non-zero. */

#line 309 "cptsvx.f"
	if (*info > 0) {
#line 310 "cptsvx.f"
	    *rcond = 0.;
#line 311 "cptsvx.f"
	    return 0;
#line 312 "cptsvx.f"
	}
#line 313 "cptsvx.f"
    }

/*     Compute the norm of the matrix A. */

#line 317 "cptsvx.f"
    anorm = clanht_("1", n, &d__[1], &e[1], (ftnlen)1);

/*     Compute the reciprocal of the condition number of A. */

#line 321 "cptsvx.f"
    cptcon_(n, &df[1], &ef[1], &anorm, rcond, &rwork[1], info);

/*     Compute the solution vectors X. */

#line 325 "cptsvx.f"
    clacpy_("Full", n, nrhs, &b[b_offset], ldb, &x[x_offset], ldx, (ftnlen)4);
#line 326 "cptsvx.f"
    cpttrs_("Lower", n, nrhs, &df[1], &ef[1], &x[x_offset], ldx, info, (
	    ftnlen)5);

/*     Use iterative refinement to improve the computed solutions and */
/*     compute error bounds and backward error estimates for them. */

#line 331 "cptsvx.f"
    cptrfs_("Lower", n, nrhs, &d__[1], &e[1], &df[1], &ef[1], &b[b_offset], 
	    ldb, &x[x_offset], ldx, &ferr[1], &berr[1], &work[1], &rwork[1], 
	    info, (ftnlen)5);

/*     Set INFO = N+1 if the matrix is singular to working precision. */

#line 336 "cptsvx.f"
    if (*rcond < slamch_("Epsilon", (ftnlen)7)) {
#line 336 "cptsvx.f"
	*info = *n + 1;
#line 336 "cptsvx.f"
    }

#line 339 "cptsvx.f"
    return 0;

/*     End of CPTSVX */

} /* cptsvx_ */

