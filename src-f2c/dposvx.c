#line 1 "dposvx.f"
/* dposvx.f -- translated by f2c (version 20100827).
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

#line 1 "dposvx.f"
/* > \brief <b> DPOSVX computes the solution to system of linear equations A * X = B for PO matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DPOSVX + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dposvx.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dposvx.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dposvx.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DPOSVX( FACT, UPLO, N, NRHS, A, LDA, AF, LDAF, EQUED, */
/*                          S, B, LDB, X, LDX, RCOND, FERR, BERR, WORK, */
/*                          IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          EQUED, FACT, UPLO */
/*       INTEGER            INFO, LDA, LDAF, LDB, LDX, N, NRHS */
/*       DOUBLE PRECISION   RCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       DOUBLE PRECISION   A( LDA, * ), AF( LDAF, * ), B( LDB, * ), */
/*      $                   BERR( * ), FERR( * ), S( * ), WORK( * ), */
/*      $                   X( LDX, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DPOSVX uses the Cholesky factorization A = U**T*U or A = L*L**T to */
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
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
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
/* >          AF is DOUBLE PRECISION array, dimension (LDAF,N) */
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
/* >          S is DOUBLE PRECISION array, dimension (N) */
/* >          The scale factors for A; not accessed if EQUED = 'N'.  S is */
/* >          an input argument if FACT = 'F'; otherwise, S is an output */
/* >          argument.  If FACT = 'F' and EQUED = 'Y', each element of S */
/* >          must be positive. */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* >          B is DOUBLE PRECISION array, dimension (LDB,NRHS) */
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
/* >          X is DOUBLE PRECISION array, dimension (LDX,NRHS) */
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
/* >          RCOND is DOUBLE PRECISION */
/* >          The estimate of the reciprocal condition number of the matrix */
/* >          A after equilibration (if done).  If RCOND is less than the */
/* >          machine precision (in particular, if RCOND = 0), the matrix */
/* >          is singular to working precision.  This condition is */
/* >          indicated by a return code of INFO > 0. */
/* > \endverbatim */
/* > */
/* > \param[out] FERR */
/* > \verbatim */
/* >          FERR is DOUBLE PRECISION array, dimension (NRHS) */
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
/* >          BERR is DOUBLE PRECISION array, dimension (NRHS) */
/* >          The componentwise relative backward error of each solution */
/* >          vector X(j) (i.e., the smallest relative change in */
/* >          any element of A or B that makes X(j) an exact solution). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension (3*N) */
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

/* > \ingroup doublePOsolve */

/*  ===================================================================== */
/* Subroutine */ int dposvx_(char *fact, char *uplo, integer *n, integer *
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
    extern doublereal dlamch_(char *, ftnlen);
    static logical nofact;
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    static doublereal bignum;
    extern /* Subroutine */ int dpocon_(char *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    integer *, ftnlen);
    static integer infequ;
    extern doublereal dlansy_(char *, char *, integer *, doublereal *, 
	    integer *, doublereal *, ftnlen, ftnlen);
    extern /* Subroutine */ int dlaqsy_(char *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, char *, 
	    ftnlen, ftnlen), dpoequ_(integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *), dporfs_(
	    char *, integer *, integer *, doublereal *, integer *, doublereal 
	    *, integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    ftnlen), dpotrf_(char *, integer *, doublereal *, integer *, 
	    integer *, ftnlen);
    static doublereal smlnum;
    extern /* Subroutine */ int dpotrs_(char *, integer *, integer *, 
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

#line 352 "dposvx.f"
    /* Parameter adjustments */
#line 352 "dposvx.f"
    a_dim1 = *lda;
#line 352 "dposvx.f"
    a_offset = 1 + a_dim1;
#line 352 "dposvx.f"
    a -= a_offset;
#line 352 "dposvx.f"
    af_dim1 = *ldaf;
#line 352 "dposvx.f"
    af_offset = 1 + af_dim1;
#line 352 "dposvx.f"
    af -= af_offset;
#line 352 "dposvx.f"
    --s;
#line 352 "dposvx.f"
    b_dim1 = *ldb;
#line 352 "dposvx.f"
    b_offset = 1 + b_dim1;
#line 352 "dposvx.f"
    b -= b_offset;
#line 352 "dposvx.f"
    x_dim1 = *ldx;
#line 352 "dposvx.f"
    x_offset = 1 + x_dim1;
#line 352 "dposvx.f"
    x -= x_offset;
#line 352 "dposvx.f"
    --ferr;
#line 352 "dposvx.f"
    --berr;
#line 352 "dposvx.f"
    --work;
#line 352 "dposvx.f"
    --iwork;
#line 352 "dposvx.f"

#line 352 "dposvx.f"
    /* Function Body */
#line 352 "dposvx.f"
    *info = 0;
#line 353 "dposvx.f"
    nofact = lsame_(fact, "N", (ftnlen)1, (ftnlen)1);
#line 354 "dposvx.f"
    equil = lsame_(fact, "E", (ftnlen)1, (ftnlen)1);
#line 355 "dposvx.f"
    if (nofact || equil) {
#line 356 "dposvx.f"
	*(unsigned char *)equed = 'N';
#line 357 "dposvx.f"
	rcequ = FALSE_;
#line 358 "dposvx.f"
    } else {
#line 359 "dposvx.f"
	rcequ = lsame_(equed, "Y", (ftnlen)1, (ftnlen)1);
#line 360 "dposvx.f"
	smlnum = dlamch_("Safe minimum", (ftnlen)12);
#line 361 "dposvx.f"
	bignum = 1. / smlnum;
#line 362 "dposvx.f"
    }

/*     Test the input parameters. */

#line 366 "dposvx.f"
    if (! nofact && ! equil && ! lsame_(fact, "F", (ftnlen)1, (ftnlen)1)) {
#line 368 "dposvx.f"
	*info = -1;
#line 369 "dposvx.f"
    } else if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, 
	    "L", (ftnlen)1, (ftnlen)1)) {
#line 371 "dposvx.f"
	*info = -2;
#line 372 "dposvx.f"
    } else if (*n < 0) {
#line 373 "dposvx.f"
	*info = -3;
#line 374 "dposvx.f"
    } else if (*nrhs < 0) {
#line 375 "dposvx.f"
	*info = -4;
#line 376 "dposvx.f"
    } else if (*lda < max(1,*n)) {
#line 377 "dposvx.f"
	*info = -6;
#line 378 "dposvx.f"
    } else if (*ldaf < max(1,*n)) {
#line 379 "dposvx.f"
	*info = -8;
#line 380 "dposvx.f"
    } else if (lsame_(fact, "F", (ftnlen)1, (ftnlen)1) && ! (rcequ || lsame_(
	    equed, "N", (ftnlen)1, (ftnlen)1))) {
#line 382 "dposvx.f"
	*info = -9;
#line 383 "dposvx.f"
    } else {
#line 384 "dposvx.f"
	if (rcequ) {
#line 385 "dposvx.f"
	    smin = bignum;
#line 386 "dposvx.f"
	    smax = 0.;
#line 387 "dposvx.f"
	    i__1 = *n;
#line 387 "dposvx.f"
	    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 388 "dposvx.f"
		d__1 = smin, d__2 = s[j];
#line 388 "dposvx.f"
		smin = min(d__1,d__2);
/* Computing MAX */
#line 389 "dposvx.f"
		d__1 = smax, d__2 = s[j];
#line 389 "dposvx.f"
		smax = max(d__1,d__2);
#line 390 "dposvx.f"
/* L10: */
#line 390 "dposvx.f"
	    }
#line 391 "dposvx.f"
	    if (smin <= 0.) {
#line 392 "dposvx.f"
		*info = -10;
#line 393 "dposvx.f"
	    } else if (*n > 0) {
#line 394 "dposvx.f"
		scond = max(smin,smlnum) / min(smax,bignum);
#line 395 "dposvx.f"
	    } else {
#line 396 "dposvx.f"
		scond = 1.;
#line 397 "dposvx.f"
	    }
#line 398 "dposvx.f"
	}
#line 399 "dposvx.f"
	if (*info == 0) {
#line 400 "dposvx.f"
	    if (*ldb < max(1,*n)) {
#line 401 "dposvx.f"
		*info = -12;
#line 402 "dposvx.f"
	    } else if (*ldx < max(1,*n)) {
#line 403 "dposvx.f"
		*info = -14;
#line 404 "dposvx.f"
	    }
#line 405 "dposvx.f"
	}
#line 406 "dposvx.f"
    }

#line 408 "dposvx.f"
    if (*info != 0) {
#line 409 "dposvx.f"
	i__1 = -(*info);
#line 409 "dposvx.f"
	xerbla_("DPOSVX", &i__1, (ftnlen)6);
#line 410 "dposvx.f"
	return 0;
#line 411 "dposvx.f"
    }

#line 413 "dposvx.f"
    if (equil) {

/*        Compute row and column scalings to equilibrate the matrix A. */

#line 417 "dposvx.f"
	dpoequ_(n, &a[a_offset], lda, &s[1], &scond, &amax, &infequ);
#line 418 "dposvx.f"
	if (infequ == 0) {

/*           Equilibrate the matrix. */

#line 422 "dposvx.f"
	    dlaqsy_(uplo, n, &a[a_offset], lda, &s[1], &scond, &amax, equed, (
		    ftnlen)1, (ftnlen)1);
#line 423 "dposvx.f"
	    rcequ = lsame_(equed, "Y", (ftnlen)1, (ftnlen)1);
#line 424 "dposvx.f"
	}
#line 425 "dposvx.f"
    }

/*     Scale the right hand side. */

#line 429 "dposvx.f"
    if (rcequ) {
#line 430 "dposvx.f"
	i__1 = *nrhs;
#line 430 "dposvx.f"
	for (j = 1; j <= i__1; ++j) {
#line 431 "dposvx.f"
	    i__2 = *n;
#line 431 "dposvx.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 432 "dposvx.f"
		b[i__ + j * b_dim1] = s[i__] * b[i__ + j * b_dim1];
#line 433 "dposvx.f"
/* L20: */
#line 433 "dposvx.f"
	    }
#line 434 "dposvx.f"
/* L30: */
#line 434 "dposvx.f"
	}
#line 435 "dposvx.f"
    }

#line 437 "dposvx.f"
    if (nofact || equil) {

/*        Compute the Cholesky factorization A = U**T *U or A = L*L**T. */

#line 441 "dposvx.f"
	dlacpy_(uplo, n, n, &a[a_offset], lda, &af[af_offset], ldaf, (ftnlen)
		1);
#line 442 "dposvx.f"
	dpotrf_(uplo, n, &af[af_offset], ldaf, info, (ftnlen)1);

/*        Return if INFO is non-zero. */

#line 446 "dposvx.f"
	if (*info > 0) {
#line 447 "dposvx.f"
	    *rcond = 0.;
#line 448 "dposvx.f"
	    return 0;
#line 449 "dposvx.f"
	}
#line 450 "dposvx.f"
    }

/*     Compute the norm of the matrix A. */

#line 454 "dposvx.f"
    anorm = dlansy_("1", uplo, n, &a[a_offset], lda, &work[1], (ftnlen)1, (
	    ftnlen)1);

/*     Compute the reciprocal of the condition number of A. */

#line 458 "dposvx.f"
    dpocon_(uplo, n, &af[af_offset], ldaf, &anorm, rcond, &work[1], &iwork[1],
	     info, (ftnlen)1);

/*     Compute the solution matrix X. */

#line 462 "dposvx.f"
    dlacpy_("Full", n, nrhs, &b[b_offset], ldb, &x[x_offset], ldx, (ftnlen)4);
#line 463 "dposvx.f"
    dpotrs_(uplo, n, nrhs, &af[af_offset], ldaf, &x[x_offset], ldx, info, (
	    ftnlen)1);

/*     Use iterative refinement to improve the computed solution and */
/*     compute error bounds and backward error estimates for it. */

#line 468 "dposvx.f"
    dporfs_(uplo, n, nrhs, &a[a_offset], lda, &af[af_offset], ldaf, &b[
	    b_offset], ldb, &x[x_offset], ldx, &ferr[1], &berr[1], &work[1], &
	    iwork[1], info, (ftnlen)1);

/*     Transform the solution matrix X to a solution of the original */
/*     system. */

#line 474 "dposvx.f"
    if (rcequ) {
#line 475 "dposvx.f"
	i__1 = *nrhs;
#line 475 "dposvx.f"
	for (j = 1; j <= i__1; ++j) {
#line 476 "dposvx.f"
	    i__2 = *n;
#line 476 "dposvx.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 477 "dposvx.f"
		x[i__ + j * x_dim1] = s[i__] * x[i__ + j * x_dim1];
#line 478 "dposvx.f"
/* L40: */
#line 478 "dposvx.f"
	    }
#line 479 "dposvx.f"
/* L50: */
#line 479 "dposvx.f"
	}
#line 480 "dposvx.f"
	i__1 = *nrhs;
#line 480 "dposvx.f"
	for (j = 1; j <= i__1; ++j) {
#line 481 "dposvx.f"
	    ferr[j] /= scond;
#line 482 "dposvx.f"
/* L60: */
#line 482 "dposvx.f"
	}
#line 483 "dposvx.f"
    }

/*     Set INFO = N+1 if the matrix is singular to working precision. */

#line 487 "dposvx.f"
    if (*rcond < dlamch_("Epsilon", (ftnlen)7)) {
#line 487 "dposvx.f"
	*info = *n + 1;
#line 487 "dposvx.f"
    }

#line 490 "dposvx.f"
    return 0;

/*     End of DPOSVX */

} /* dposvx_ */

