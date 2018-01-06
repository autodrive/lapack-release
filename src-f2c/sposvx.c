#line 1 "sposvx.f"
/* sposvx.f -- translated by f2c (version 20100827).
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

#line 1 "sposvx.f"
/* > \brief <b> SPOSVX computes the solution to system of linear equations A * X = B for PO matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SPOSVX + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sposvx.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sposvx.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sposvx.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SPOSVX( FACT, UPLO, N, NRHS, A, LDA, AF, LDAF, EQUED, */
/*                          S, B, LDB, X, LDX, RCOND, FERR, BERR, WORK, */
/*                          IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          EQUED, FACT, UPLO */
/*       INTEGER            INFO, LDA, LDAF, LDB, LDX, N, NRHS */
/*       REAL               RCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       REAL               A( LDA, * ), AF( LDAF, * ), B( LDB, * ), */
/*      $                   BERR( * ), FERR( * ), S( * ), WORK( * ), */
/*      $                   X( LDX, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SPOSVX uses the Cholesky factorization A = U**T*U or A = L*L**T to */
/* > compute the solution to a real system of linear equations */
/* >    A * X = B, */
/* > where A is an N-by-N symmetric positive definite matrix and X and B */
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
/* > 1. If FACT = 'E', real scaling factors are computed to equilibrate */
/* >    the system: */
/* >       diag(S) * A * diag(S) * inv(diag(S)) * X = diag(S) * B */
/* >    Whether or not the system will be equilibrated depends on the */
/* >    scaling of the matrix A, but if equilibration is used, A is */
/* >    overwritten by diag(S)*A*diag(S) and B by diag(S)*B. */
/* > */
/* > 2. If FACT = 'N' or 'E', the Cholesky decomposition is used to */
/* >    factor the matrix A (after equilibration if FACT = 'E') as */
/* >       A = U**T* U,  if UPLO = 'U', or */
/* >       A = L * L**T,  if UPLO = 'L', */
/* >    where U is an upper triangular matrix and L is a lower triangular */
/* >    matrix. */
/* > */
/* > 3. If the leading i-by-i principal minor is not positive definite, */
/* >    then the routine returns with INFO = i. Otherwise, the factored */
/* >    form of A is used to estimate the condition number of the matrix */
/* >    A.  If the reciprocal of the condition number is less than machine */
/* >    precision, INFO = N+1 is returned as a warning, but the routine */
/* >    still goes on to solve for X and compute error bounds as */
/* >    described below. */
/* > */
/* > 4. The system of equations is solved for X using the factored form */
/* >    of A. */
/* > */
/* > 5. Iterative refinement is applied to improve the computed solution */
/* >    matrix and calculate error bounds and backward error estimates */
/* >    for it. */
/* > */
/* > 6. If equilibration was used, the matrix X is premultiplied by */
/* >    diag(S) so that it solves the original system before */
/* >    equilibration. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] FACT */
/* > \verbatim */
/* >          FACT is CHARACTER*1 */
/* >          Specifies whether or not the factored form of the matrix A is */
/* >          supplied on entry, and if not, whether the matrix A should be */
/* >          equilibrated before it is factored. */
/* >          = 'F':  On entry, AF contains the factored form of A. */
/* >                  If EQUED = 'Y', the matrix A has been equilibrated */
/* >                  with scaling factors given by S.  A and AF will not */
/* >                  be modified. */
/* >          = 'N':  The matrix A will be copied to AF and factored. */
/* >          = 'E':  The matrix A will be equilibrated if necessary, then */
/* >                  copied to AF and factored. */
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
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is REAL array, dimension (LDA,N) */
/* >          On entry, the symmetric matrix A, except if FACT = 'F' and */
/* >          EQUED = 'Y', then A must contain the equilibrated matrix */
/* >          diag(S)*A*diag(S).  If UPLO = 'U', the leading */
/* >          N-by-N upper triangular part of A contains the upper */
/* >          triangular part of the matrix A, and the strictly lower */
/* >          triangular part of A is not referenced.  If UPLO = 'L', the */
/* >          leading N-by-N lower triangular part of A contains the lower */
/* >          triangular part of the matrix A, and the strictly upper */
/* >          triangular part of A is not referenced.  A is not modified if */
/* >          FACT = 'F' or 'N', or if FACT = 'E' and EQUED = 'N' on exit. */
/* > */
/* >          On exit, if FACT = 'E' and EQUED = 'Y', A is overwritten by */
/* >          diag(S)*A*diag(S). */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] AF */
/* > \verbatim */
/* >          AF is REAL array, dimension (LDAF,N) */
/* >          If FACT = 'F', then AF is an input argument and on entry */
/* >          contains the triangular factor U or L from the Cholesky */
/* >          factorization A = U**T*U or A = L*L**T, in the same storage */
/* >          format as A.  If EQUED .ne. 'N', then AF is the factored form */
/* >          of the equilibrated matrix diag(S)*A*diag(S). */
/* > */
/* >          If FACT = 'N', then AF is an output argument and on exit */
/* >          returns the triangular factor U or L from the Cholesky */
/* >          factorization A = U**T*U or A = L*L**T of the original */
/* >          matrix A. */
/* > */
/* >          If FACT = 'E', then AF is an output argument and on exit */
/* >          returns the triangular factor U or L from the Cholesky */
/* >          factorization A = U**T*U or A = L*L**T of the equilibrated */
/* >          matrix A (see the description of A for the form of the */
/* >          equilibrated matrix). */
/* > \endverbatim */
/* > */
/* > \param[in] LDAF */
/* > \verbatim */
/* >          LDAF is INTEGER */
/* >          The leading dimension of the array AF.  LDAF >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] EQUED */
/* > \verbatim */
/* >          EQUED is CHARACTER*1 */
/* >          Specifies the form of equilibration that was done. */
/* >          = 'N':  No equilibration (always true if FACT = 'N'). */
/* >          = 'Y':  Equilibration was done, i.e., A has been replaced by */
/* >                  diag(S) * A * diag(S). */
/* >          EQUED is an input argument if FACT = 'F'; otherwise, it is an */
/* >          output argument. */
/* > \endverbatim */
/* > */
/* > \param[in,out] S */
/* > \verbatim */
/* >          S is REAL array, dimension (N) */
/* >          The scale factors for A; not accessed if EQUED = 'N'.  S is */
/* >          an input argument if FACT = 'F'; otherwise, S is an output */
/* >          argument.  If FACT = 'F' and EQUED = 'Y', each element of S */
/* >          must be positive. */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* >          B is REAL array, dimension (LDB,NRHS) */
/* >          On entry, the N-by-NRHS right hand side matrix B. */
/* >          On exit, if EQUED = 'N', B is not modified; if EQUED = 'Y', */
/* >          B is overwritten by diag(S) * B. */
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
/* >          If INFO = 0 or INFO = N+1, the N-by-NRHS solution matrix X to */
/* >          the original system of equations.  Note that if EQUED = 'Y', */
/* >          A and B are modified on exit, and the solution to the */
/* >          equilibrated system is inv(diag(S))*X. */
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
/* >          A after equilibration (if done).  If RCOND is less than the */
/* >          machine precision (in particular, if RCOND = 0), the matrix */
/* >          is singular to working precision.  This condition is */
/* >          indicated by a return code of INFO > 0. */
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
/* >          = 0: successful exit */
/* >          < 0: if INFO = -i, the i-th argument had an illegal value */
/* >          > 0: if INFO = i, and i is */
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

/* > \date April 2012 */

/* > \ingroup realPOsolve */

/*  ===================================================================== */
/* Subroutine */ int sposvx_(char *fact, char *uplo, integer *n, integer *
	nrhs, doublereal *a, integer *lda, doublereal *af, integer *ldaf, 
	char *equed, doublereal *s, doublereal *b, integer *ldb, doublereal *
	x, integer *ldx, doublereal *rcond, doublereal *ferr, doublereal *
	berr, doublereal *work, integer *iwork, integer *info, ftnlen 
	fact_len, ftnlen uplo_len, ftnlen equed_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, af_dim1, af_offset, b_dim1, b_offset, x_dim1, 
	    x_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, j;
    static doublereal amax, smin, smax;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal scond, anorm;
    static logical equil, rcequ;
    extern doublereal slamch_(char *, ftnlen);
    static logical nofact;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal bignum;
    static integer infequ;
    extern /* Subroutine */ int slacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    spocon_(char *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, ftnlen);
    extern doublereal slansy_(char *, char *, integer *, doublereal *, 
	    integer *, doublereal *, ftnlen, ftnlen);
    static doublereal smlnum;
    extern /* Subroutine */ int slaqsy_(char *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, char *, 
	    ftnlen, ftnlen), spoequ_(integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *), sporfs_(
	    char *, integer *, integer *, doublereal *, integer *, doublereal 
	    *, integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    ftnlen), spotrf_(char *, integer *, doublereal *, integer *, 
	    integer *, ftnlen), spotrs_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen);


/*  -- LAPACK driver routine (version 3.7.0) -- */
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

#line 352 "sposvx.f"
    /* Parameter adjustments */
#line 352 "sposvx.f"
    a_dim1 = *lda;
#line 352 "sposvx.f"
    a_offset = 1 + a_dim1;
#line 352 "sposvx.f"
    a -= a_offset;
#line 352 "sposvx.f"
    af_dim1 = *ldaf;
#line 352 "sposvx.f"
    af_offset = 1 + af_dim1;
#line 352 "sposvx.f"
    af -= af_offset;
#line 352 "sposvx.f"
    --s;
#line 352 "sposvx.f"
    b_dim1 = *ldb;
#line 352 "sposvx.f"
    b_offset = 1 + b_dim1;
#line 352 "sposvx.f"
    b -= b_offset;
#line 352 "sposvx.f"
    x_dim1 = *ldx;
#line 352 "sposvx.f"
    x_offset = 1 + x_dim1;
#line 352 "sposvx.f"
    x -= x_offset;
#line 352 "sposvx.f"
    --ferr;
#line 352 "sposvx.f"
    --berr;
#line 352 "sposvx.f"
    --work;
#line 352 "sposvx.f"
    --iwork;
#line 352 "sposvx.f"

#line 352 "sposvx.f"
    /* Function Body */
#line 352 "sposvx.f"
    *info = 0;
#line 353 "sposvx.f"
    nofact = lsame_(fact, "N", (ftnlen)1, (ftnlen)1);
#line 354 "sposvx.f"
    equil = lsame_(fact, "E", (ftnlen)1, (ftnlen)1);
#line 355 "sposvx.f"
    if (nofact || equil) {
#line 356 "sposvx.f"
	*(unsigned char *)equed = 'N';
#line 357 "sposvx.f"
	rcequ = FALSE_;
#line 358 "sposvx.f"
    } else {
#line 359 "sposvx.f"
	rcequ = lsame_(equed, "Y", (ftnlen)1, (ftnlen)1);
#line 360 "sposvx.f"
	smlnum = slamch_("Safe minimum", (ftnlen)12);
#line 361 "sposvx.f"
	bignum = 1. / smlnum;
#line 362 "sposvx.f"
    }

/*     Test the input parameters. */

#line 366 "sposvx.f"
    if (! nofact && ! equil && ! lsame_(fact, "F", (ftnlen)1, (ftnlen)1)) {
#line 368 "sposvx.f"
	*info = -1;
#line 369 "sposvx.f"
    } else if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, 
	    "L", (ftnlen)1, (ftnlen)1)) {
#line 371 "sposvx.f"
	*info = -2;
#line 372 "sposvx.f"
    } else if (*n < 0) {
#line 373 "sposvx.f"
	*info = -3;
#line 374 "sposvx.f"
    } else if (*nrhs < 0) {
#line 375 "sposvx.f"
	*info = -4;
#line 376 "sposvx.f"
    } else if (*lda < max(1,*n)) {
#line 377 "sposvx.f"
	*info = -6;
#line 378 "sposvx.f"
    } else if (*ldaf < max(1,*n)) {
#line 379 "sposvx.f"
	*info = -8;
#line 380 "sposvx.f"
    } else if (lsame_(fact, "F", (ftnlen)1, (ftnlen)1) && ! (rcequ || lsame_(
	    equed, "N", (ftnlen)1, (ftnlen)1))) {
#line 382 "sposvx.f"
	*info = -9;
#line 383 "sposvx.f"
    } else {
#line 384 "sposvx.f"
	if (rcequ) {
#line 385 "sposvx.f"
	    smin = bignum;
#line 386 "sposvx.f"
	    smax = 0.;
#line 387 "sposvx.f"
	    i__1 = *n;
#line 387 "sposvx.f"
	    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 388 "sposvx.f"
		d__1 = smin, d__2 = s[j];
#line 388 "sposvx.f"
		smin = min(d__1,d__2);
/* Computing MAX */
#line 389 "sposvx.f"
		d__1 = smax, d__2 = s[j];
#line 389 "sposvx.f"
		smax = max(d__1,d__2);
#line 390 "sposvx.f"
/* L10: */
#line 390 "sposvx.f"
	    }
#line 391 "sposvx.f"
	    if (smin <= 0.) {
#line 392 "sposvx.f"
		*info = -10;
#line 393 "sposvx.f"
	    } else if (*n > 0) {
#line 394 "sposvx.f"
		scond = max(smin,smlnum) / min(smax,bignum);
#line 395 "sposvx.f"
	    } else {
#line 396 "sposvx.f"
		scond = 1.;
#line 397 "sposvx.f"
	    }
#line 398 "sposvx.f"
	}
#line 399 "sposvx.f"
	if (*info == 0) {
#line 400 "sposvx.f"
	    if (*ldb < max(1,*n)) {
#line 401 "sposvx.f"
		*info = -12;
#line 402 "sposvx.f"
	    } else if (*ldx < max(1,*n)) {
#line 403 "sposvx.f"
		*info = -14;
#line 404 "sposvx.f"
	    }
#line 405 "sposvx.f"
	}
#line 406 "sposvx.f"
    }

#line 408 "sposvx.f"
    if (*info != 0) {
#line 409 "sposvx.f"
	i__1 = -(*info);
#line 409 "sposvx.f"
	xerbla_("SPOSVX", &i__1, (ftnlen)6);
#line 410 "sposvx.f"
	return 0;
#line 411 "sposvx.f"
    }

#line 413 "sposvx.f"
    if (equil) {

/*        Compute row and column scalings to equilibrate the matrix A. */

#line 417 "sposvx.f"
	spoequ_(n, &a[a_offset], lda, &s[1], &scond, &amax, &infequ);
#line 418 "sposvx.f"
	if (infequ == 0) {

/*           Equilibrate the matrix. */

#line 422 "sposvx.f"
	    slaqsy_(uplo, n, &a[a_offset], lda, &s[1], &scond, &amax, equed, (
		    ftnlen)1, (ftnlen)1);
#line 423 "sposvx.f"
	    rcequ = lsame_(equed, "Y", (ftnlen)1, (ftnlen)1);
#line 424 "sposvx.f"
	}
#line 425 "sposvx.f"
    }

/*     Scale the right hand side. */

#line 429 "sposvx.f"
    if (rcequ) {
#line 430 "sposvx.f"
	i__1 = *nrhs;
#line 430 "sposvx.f"
	for (j = 1; j <= i__1; ++j) {
#line 431 "sposvx.f"
	    i__2 = *n;
#line 431 "sposvx.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 432 "sposvx.f"
		b[i__ + j * b_dim1] = s[i__] * b[i__ + j * b_dim1];
#line 433 "sposvx.f"
/* L20: */
#line 433 "sposvx.f"
	    }
#line 434 "sposvx.f"
/* L30: */
#line 434 "sposvx.f"
	}
#line 435 "sposvx.f"
    }

#line 437 "sposvx.f"
    if (nofact || equil) {

/*        Compute the Cholesky factorization A = U**T *U or A = L*L**T. */

#line 441 "sposvx.f"
	slacpy_(uplo, n, n, &a[a_offset], lda, &af[af_offset], ldaf, (ftnlen)
		1);
#line 442 "sposvx.f"
	spotrf_(uplo, n, &af[af_offset], ldaf, info, (ftnlen)1);

/*        Return if INFO is non-zero. */

#line 446 "sposvx.f"
	if (*info > 0) {
#line 447 "sposvx.f"
	    *rcond = 0.;
#line 448 "sposvx.f"
	    return 0;
#line 449 "sposvx.f"
	}
#line 450 "sposvx.f"
    }

/*     Compute the norm of the matrix A. */

#line 454 "sposvx.f"
    anorm = slansy_("1", uplo, n, &a[a_offset], lda, &work[1], (ftnlen)1, (
	    ftnlen)1);

/*     Compute the reciprocal of the condition number of A. */

#line 458 "sposvx.f"
    spocon_(uplo, n, &af[af_offset], ldaf, &anorm, rcond, &work[1], &iwork[1],
	     info, (ftnlen)1);

/*     Compute the solution matrix X. */

#line 462 "sposvx.f"
    slacpy_("Full", n, nrhs, &b[b_offset], ldb, &x[x_offset], ldx, (ftnlen)4);
#line 463 "sposvx.f"
    spotrs_(uplo, n, nrhs, &af[af_offset], ldaf, &x[x_offset], ldx, info, (
	    ftnlen)1);

/*     Use iterative refinement to improve the computed solution and */
/*     compute error bounds and backward error estimates for it. */

#line 468 "sposvx.f"
    sporfs_(uplo, n, nrhs, &a[a_offset], lda, &af[af_offset], ldaf, &b[
	    b_offset], ldb, &x[x_offset], ldx, &ferr[1], &berr[1], &work[1], &
	    iwork[1], info, (ftnlen)1);

/*     Transform the solution matrix X to a solution of the original */
/*     system. */

#line 474 "sposvx.f"
    if (rcequ) {
#line 475 "sposvx.f"
	i__1 = *nrhs;
#line 475 "sposvx.f"
	for (j = 1; j <= i__1; ++j) {
#line 476 "sposvx.f"
	    i__2 = *n;
#line 476 "sposvx.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 477 "sposvx.f"
		x[i__ + j * x_dim1] = s[i__] * x[i__ + j * x_dim1];
#line 478 "sposvx.f"
/* L40: */
#line 478 "sposvx.f"
	    }
#line 479 "sposvx.f"
/* L50: */
#line 479 "sposvx.f"
	}
#line 480 "sposvx.f"
	i__1 = *nrhs;
#line 480 "sposvx.f"
	for (j = 1; j <= i__1; ++j) {
#line 481 "sposvx.f"
	    ferr[j] /= scond;
#line 482 "sposvx.f"
/* L60: */
#line 482 "sposvx.f"
	}
#line 483 "sposvx.f"
    }

/*     Set INFO = N+1 if the matrix is singular to working precision. */

#line 487 "sposvx.f"
    if (*rcond < slamch_("Epsilon", (ftnlen)7)) {
#line 487 "sposvx.f"
	*info = *n + 1;
#line 487 "sposvx.f"
    }

#line 490 "sposvx.f"
    return 0;

/*     End of SPOSVX */

} /* sposvx_ */

