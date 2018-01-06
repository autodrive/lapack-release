#line 1 "cspsvx.f"
/* cspsvx.f -- translated by f2c (version 20100827).
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

#line 1 "cspsvx.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief <b> CSPSVX computes the solution to system of linear equations A * X = B for OTHER matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CSPSVX + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cspsvx.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cspsvx.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cspsvx.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CSPSVX( FACT, UPLO, N, NRHS, AP, AFP, IPIV, B, LDB, X, */
/*                          LDX, RCOND, FERR, BERR, WORK, RWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          FACT, UPLO */
/*       INTEGER            INFO, LDB, LDX, N, NRHS */
/*       REAL               RCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       REAL               BERR( * ), FERR( * ), RWORK( * ) */
/*       COMPLEX            AFP( * ), AP( * ), B( LDB, * ), WORK( * ), */
/*      $                   X( LDX, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CSPSVX uses the diagonal pivoting factorization A = U*D*U**T or */
/* > A = L*D*L**T to compute the solution to a complex system of linear */
/* > equations A * X = B, where A is an N-by-N symmetric matrix stored */
/* > in packed format and X and B are N-by-NRHS matrices. */
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
/* > 1. If FACT = 'N', the diagonal pivoting method is used to factor A as */
/* >       A = U * D * U**T,  if UPLO = 'U', or */
/* >       A = L * D * L**T,  if UPLO = 'L', */
/* >    where U (or L) is a product of permutation and unit upper (lower) */
/* >    triangular matrices and D is symmetric and block diagonal with */
/* >    1-by-1 and 2-by-2 diagonal blocks. */
/* > */
/* > 2. If some D(i,i)=0, so that D is exactly singular, then the routine */
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
/* >          = 'F':  On entry, AFP and IPIV contain the factored form */
/* >                  of A.  AP, AFP and IPIV will not be modified. */
/* >          = 'N':  The matrix A will be copied to AFP and factored. */
/* > \endverbatim */
/* > */
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
/* >          of the matrices B and X.  NRHS >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] AP */
/* > \verbatim */
/* >          AP is COMPLEX array, dimension (N*(N+1)/2) */
/* >          The upper or lower triangle of the symmetric matrix A, packed */
/* >          columnwise in a linear array.  The j-th column of A is stored */
/* >          in the array AP as follows: */
/* >          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j; */
/* >          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n. */
/* >          See below for further details. */
/* > \endverbatim */
/* > */
/* > \param[in,out] AFP */
/* > \verbatim */
/* >          AFP is COMPLEX array, dimension (N*(N+1)/2) */
/* >          If FACT = 'F', then AFP is an input argument and on entry */
/* >          contains the block diagonal matrix D and the multipliers used */
/* >          to obtain the factor U or L from the factorization */
/* >          A = U*D*U**T or A = L*D*L**T as computed by CSPTRF, stored as */
/* >          a packed triangular matrix in the same storage format as A. */
/* > */
/* >          If FACT = 'N', then AFP is an output argument and on exit */
/* >          contains the block diagonal matrix D and the multipliers used */
/* >          to obtain the factor U or L from the factorization */
/* >          A = U*D*U**T or A = L*D*L**T as computed by CSPTRF, stored as */
/* >          a packed triangular matrix in the same storage format as A. */
/* > \endverbatim */
/* > */
/* > \param[in,out] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N) */
/* >          If FACT = 'F', then IPIV is an input argument and on entry */
/* >          contains details of the interchanges and the block structure */
/* >          of D, as determined by CSPTRF. */
/* >          If IPIV(k) > 0, then rows and columns k and IPIV(k) were */
/* >          interchanged and D(k,k) is a 1-by-1 diagonal block. */
/* >          If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0, then rows and */
/* >          columns k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k) */
/* >          is a 2-by-2 diagonal block.  If UPLO = 'L' and IPIV(k) = */
/* >          IPIV(k+1) < 0, then rows and columns k+1 and -IPIV(k) were */
/* >          interchanged and D(k:k+1,k:k+1) is a 2-by-2 diagonal block. */
/* > */
/* >          If FACT = 'N', then IPIV is an output argument and on exit */
/* >          contains details of the interchanges and the block structure */
/* >          of D, as determined by CSPTRF. */
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
/* >          = 0: successful exit */
/* >          < 0: if INFO = -i, the i-th argument had an illegal value */
/* >          > 0:  if INFO = i, and i is */
/* >                <= N:  D(i,i) is exactly zero.  The factorization */
/* >                       has been completed but the factor D is exactly */
/* >                       singular, so the solution and error bounds could */
/* >                       not be computed. RCOND = 0 is returned. */
/* >                = N+1: D is nonsingular, but RCOND is less than machine */
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

/* > \date April 2012 */

/* > \ingroup complexOTHERsolve */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  The packed storage scheme is illustrated by the following example */
/* >  when N = 4, UPLO = 'U': */
/* > */
/* >  Two-dimensional storage of the symmetric matrix A: */
/* > */
/* >     a11 a12 a13 a14 */
/* >         a22 a23 a24 */
/* >             a33 a34     (aij = aji) */
/* >                 a44 */
/* > */
/* >  Packed storage of the upper triangle of A: */
/* > */
/* >  AP = [ a11, a12, a22, a13, a23, a33, a14, a24, a34, a44 ] */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int cspsvx_(char *fact, char *uplo, integer *n, integer *
	nrhs, doublecomplex *ap, doublecomplex *afp, integer *ipiv, 
	doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx, 
	doublereal *rcond, doublereal *ferr, doublereal *berr, doublecomplex *
	work, doublereal *rwork, integer *info, ftnlen fact_len, ftnlen 
	uplo_len)
{
    /* System generated locals */
    integer b_dim1, b_offset, x_dim1, x_offset, i__1;

    /* Local variables */
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal anorm;
    extern /* Subroutine */ int ccopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    extern doublereal slamch_(char *, ftnlen);
    static logical nofact;
    extern /* Subroutine */ int clacpy_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    extern doublereal clansp_(char *, char *, integer *, doublecomplex *, 
	    doublereal *, ftnlen, ftnlen);
    extern /* Subroutine */ int cspcon_(char *, integer *, doublecomplex *, 
	    integer *, doublereal *, doublereal *, doublecomplex *, integer *,
	     ftnlen), csprfs_(char *, integer *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublereal *, doublereal *, 
	    doublecomplex *, doublereal *, integer *, ftnlen), csptrf_(char *,
	     integer *, doublecomplex *, integer *, integer *, ftnlen), 
	    csptrs_(char *, integer *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, integer *, ftnlen);


/*  -- LAPACK driver routine (version 3.4.1) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     April 2012 */

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

#line 323 "cspsvx.f"
    /* Parameter adjustments */
#line 323 "cspsvx.f"
    --ap;
#line 323 "cspsvx.f"
    --afp;
#line 323 "cspsvx.f"
    --ipiv;
#line 323 "cspsvx.f"
    b_dim1 = *ldb;
#line 323 "cspsvx.f"
    b_offset = 1 + b_dim1;
#line 323 "cspsvx.f"
    b -= b_offset;
#line 323 "cspsvx.f"
    x_dim1 = *ldx;
#line 323 "cspsvx.f"
    x_offset = 1 + x_dim1;
#line 323 "cspsvx.f"
    x -= x_offset;
#line 323 "cspsvx.f"
    --ferr;
#line 323 "cspsvx.f"
    --berr;
#line 323 "cspsvx.f"
    --work;
#line 323 "cspsvx.f"
    --rwork;
#line 323 "cspsvx.f"

#line 323 "cspsvx.f"
    /* Function Body */
#line 323 "cspsvx.f"
    *info = 0;
#line 324 "cspsvx.f"
    nofact = lsame_(fact, "N", (ftnlen)1, (ftnlen)1);
#line 325 "cspsvx.f"
    if (! nofact && ! lsame_(fact, "F", (ftnlen)1, (ftnlen)1)) {
#line 326 "cspsvx.f"
	*info = -1;
#line 327 "cspsvx.f"
    } else if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, 
	    "L", (ftnlen)1, (ftnlen)1)) {
#line 329 "cspsvx.f"
	*info = -2;
#line 330 "cspsvx.f"
    } else if (*n < 0) {
#line 331 "cspsvx.f"
	*info = -3;
#line 332 "cspsvx.f"
    } else if (*nrhs < 0) {
#line 333 "cspsvx.f"
	*info = -4;
#line 334 "cspsvx.f"
    } else if (*ldb < max(1,*n)) {
#line 335 "cspsvx.f"
	*info = -9;
#line 336 "cspsvx.f"
    } else if (*ldx < max(1,*n)) {
#line 337 "cspsvx.f"
	*info = -11;
#line 338 "cspsvx.f"
    }
#line 339 "cspsvx.f"
    if (*info != 0) {
#line 340 "cspsvx.f"
	i__1 = -(*info);
#line 340 "cspsvx.f"
	xerbla_("CSPSVX", &i__1, (ftnlen)6);
#line 341 "cspsvx.f"
	return 0;
#line 342 "cspsvx.f"
    }

#line 344 "cspsvx.f"
    if (nofact) {

/*        Compute the factorization A = U*D*U**T or A = L*D*L**T. */

#line 348 "cspsvx.f"
	i__1 = *n * (*n + 1) / 2;
#line 348 "cspsvx.f"
	ccopy_(&i__1, &ap[1], &c__1, &afp[1], &c__1);
#line 349 "cspsvx.f"
	csptrf_(uplo, n, &afp[1], &ipiv[1], info, (ftnlen)1);

/*        Return if INFO is non-zero. */

#line 353 "cspsvx.f"
	if (*info > 0) {
#line 354 "cspsvx.f"
	    *rcond = 0.;
#line 355 "cspsvx.f"
	    return 0;
#line 356 "cspsvx.f"
	}
#line 357 "cspsvx.f"
    }

/*     Compute the norm of the matrix A. */

#line 361 "cspsvx.f"
    anorm = clansp_("I", uplo, n, &ap[1], &rwork[1], (ftnlen)1, (ftnlen)1);

/*     Compute the reciprocal of the condition number of A. */

#line 365 "cspsvx.f"
    cspcon_(uplo, n, &afp[1], &ipiv[1], &anorm, rcond, &work[1], info, (
	    ftnlen)1);

/*     Compute the solution vectors X. */

#line 369 "cspsvx.f"
    clacpy_("Full", n, nrhs, &b[b_offset], ldb, &x[x_offset], ldx, (ftnlen)4);
#line 370 "cspsvx.f"
    csptrs_(uplo, n, nrhs, &afp[1], &ipiv[1], &x[x_offset], ldx, info, (
	    ftnlen)1);

/*     Use iterative refinement to improve the computed solutions and */
/*     compute error bounds and backward error estimates for them. */

#line 375 "cspsvx.f"
    csprfs_(uplo, n, nrhs, &ap[1], &afp[1], &ipiv[1], &b[b_offset], ldb, &x[
	    x_offset], ldx, &ferr[1], &berr[1], &work[1], &rwork[1], info, (
	    ftnlen)1);

/*     Set INFO = N+1 if the matrix is singular to working precision. */

#line 380 "cspsvx.f"
    if (*rcond < slamch_("Epsilon", (ftnlen)7)) {
#line 380 "cspsvx.f"
	*info = *n + 1;
#line 380 "cspsvx.f"
    }

#line 383 "cspsvx.f"
    return 0;

/*     End of CSPSVX */

} /* cspsvx_ */

