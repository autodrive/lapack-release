#line 1 "zposvx.f"
/* zposvx.f -- translated by f2c (version 20100827).
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

#line 1 "zposvx.f"
/* > \brief <b> ZPOSVX computes the solution to system of linear equations A * X = B for PO matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZPOSVX + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zposvx.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zposvx.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zposvx.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZPOSVX( FACT, UPLO, N, NRHS, A, LDA, AF, LDAF, EQUED, */
/*                          S, B, LDB, X, LDX, RCOND, FERR, BERR, WORK, */
/*                          RWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          EQUED, FACT, UPLO */
/*       INTEGER            INFO, LDA, LDAF, LDB, LDX, N, NRHS */
/*       DOUBLE PRECISION   RCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   BERR( * ), FERR( * ), RWORK( * ), S( * ) */
/*       COMPLEX*16         A( LDA, * ), AF( LDAF, * ), B( LDB, * ), */
/*      $                   WORK( * ), X( LDX, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZPOSVX uses the Cholesky factorization A = U**H*U or A = L*L**H to */
/* > compute the solution to a complex system of linear equations */
/* >    A * X = B, */
/* > where A is an N-by-N Hermitian positive definite matrix and X and B */
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
/* >       A = U**H* U,  if UPLO = 'U', or */
/* >       A = L * L**H,  if UPLO = 'L', */
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
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
/* >          On entry, the Hermitian matrix A, except if FACT = 'F' and */
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
/* >          AF is COMPLEX*16 array, dimension (LDAF,N) */
/* >          If FACT = 'F', then AF is an input argument and on entry */
/* >          contains the triangular factor U or L from the Cholesky */
/* >          factorization A = U**H *U or A = L*L**H, in the same storage */
/* >          format as A.  If EQUED .ne. 'N', then AF is the factored form */
/* >          of the equilibrated matrix diag(S)*A*diag(S). */
/* > */
/* >          If FACT = 'N', then AF is an output argument and on exit */
/* >          returns the triangular factor U or L from the Cholesky */
/* >          factorization A = U**H *U or A = L*L**H of the original */
/* >          matrix A. */
/* > */
/* >          If FACT = 'E', then AF is an output argument and on exit */
/* >          returns the triangular factor U or L from the Cholesky */
/* >          factorization A = U**H *U or A = L*L**H of the equilibrated */
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
/* >          B is COMPLEX*16 array, dimension (LDB,NRHS) */
/* >          On entry, the N-by-NRHS righthand side matrix B. */
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
/* >          X is COMPLEX*16 array, dimension (LDX,NRHS) */
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
/* >          WORK is COMPLEX*16 array, dimension (2*N) */
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

/* > \ingroup complex16POsolve */

/*  ===================================================================== */
/* Subroutine */ int zposvx_(char *fact, char *uplo, integer *n, integer *
	nrhs, doublecomplex *a, integer *lda, doublecomplex *af, integer *
	ldaf, char *equed, doublereal *s, doublecomplex *b, integer *ldb, 
	doublecomplex *x, integer *ldx, doublereal *rcond, doublereal *ferr, 
	doublereal *berr, doublecomplex *work, doublereal *rwork, integer *
	info, ftnlen fact_len, ftnlen uplo_len, ftnlen equed_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, af_dim1, af_offset, b_dim1, b_offset, x_dim1, 
	    x_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2;
    doublecomplex z__1;

    /* Local variables */
    static integer i__, j;
    static doublereal amax, smin, smax;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal scond, anorm;
    static logical equil, rcequ;
    extern doublereal dlamch_(char *, ftnlen);
    static logical nofact;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal bignum;
    extern doublereal zlanhe_(char *, char *, integer *, doublecomplex *, 
	    integer *, doublereal *, ftnlen, ftnlen);
    extern /* Subroutine */ int zlaqhe_(char *, integer *, doublecomplex *, 
	    integer *, doublereal *, doublereal *, doublereal *, char *, 
	    ftnlen, ftnlen);
    static integer infequ;
    extern /* Subroutine */ int zlacpy_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen), 
	    zpocon_(char *, integer *, doublecomplex *, integer *, doublereal 
	    *, doublereal *, doublecomplex *, doublereal *, integer *, ftnlen)
	    ;
    static doublereal smlnum;
    extern /* Subroutine */ int zpoequ_(integer *, doublecomplex *, integer *,
	     doublereal *, doublereal *, doublereal *, integer *), zporfs_(
	    char *, integer *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublereal *, doublereal *, 
	    doublecomplex *, doublereal *, integer *, ftnlen), zpotrf_(char *,
	     integer *, doublecomplex *, integer *, integer *, ftnlen), 
	    zpotrs_(char *, integer *, integer *, doublecomplex *, integer *, 
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

#line 350 "zposvx.f"
    /* Parameter adjustments */
#line 350 "zposvx.f"
    a_dim1 = *lda;
#line 350 "zposvx.f"
    a_offset = 1 + a_dim1;
#line 350 "zposvx.f"
    a -= a_offset;
#line 350 "zposvx.f"
    af_dim1 = *ldaf;
#line 350 "zposvx.f"
    af_offset = 1 + af_dim1;
#line 350 "zposvx.f"
    af -= af_offset;
#line 350 "zposvx.f"
    --s;
#line 350 "zposvx.f"
    b_dim1 = *ldb;
#line 350 "zposvx.f"
    b_offset = 1 + b_dim1;
#line 350 "zposvx.f"
    b -= b_offset;
#line 350 "zposvx.f"
    x_dim1 = *ldx;
#line 350 "zposvx.f"
    x_offset = 1 + x_dim1;
#line 350 "zposvx.f"
    x -= x_offset;
#line 350 "zposvx.f"
    --ferr;
#line 350 "zposvx.f"
    --berr;
#line 350 "zposvx.f"
    --work;
#line 350 "zposvx.f"
    --rwork;
#line 350 "zposvx.f"

#line 350 "zposvx.f"
    /* Function Body */
#line 350 "zposvx.f"
    *info = 0;
#line 351 "zposvx.f"
    nofact = lsame_(fact, "N", (ftnlen)1, (ftnlen)1);
#line 352 "zposvx.f"
    equil = lsame_(fact, "E", (ftnlen)1, (ftnlen)1);
#line 353 "zposvx.f"
    if (nofact || equil) {
#line 354 "zposvx.f"
	*(unsigned char *)equed = 'N';
#line 355 "zposvx.f"
	rcequ = FALSE_;
#line 356 "zposvx.f"
    } else {
#line 357 "zposvx.f"
	rcequ = lsame_(equed, "Y", (ftnlen)1, (ftnlen)1);
#line 358 "zposvx.f"
	smlnum = dlamch_("Safe minimum", (ftnlen)12);
#line 359 "zposvx.f"
	bignum = 1. / smlnum;
#line 360 "zposvx.f"
    }

/*     Test the input parameters. */

#line 364 "zposvx.f"
    if (! nofact && ! equil && ! lsame_(fact, "F", (ftnlen)1, (ftnlen)1)) {
#line 366 "zposvx.f"
	*info = -1;
#line 367 "zposvx.f"
    } else if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, 
	    "L", (ftnlen)1, (ftnlen)1)) {
#line 369 "zposvx.f"
	*info = -2;
#line 370 "zposvx.f"
    } else if (*n < 0) {
#line 371 "zposvx.f"
	*info = -3;
#line 372 "zposvx.f"
    } else if (*nrhs < 0) {
#line 373 "zposvx.f"
	*info = -4;
#line 374 "zposvx.f"
    } else if (*lda < max(1,*n)) {
#line 375 "zposvx.f"
	*info = -6;
#line 376 "zposvx.f"
    } else if (*ldaf < max(1,*n)) {
#line 377 "zposvx.f"
	*info = -8;
#line 378 "zposvx.f"
    } else if (lsame_(fact, "F", (ftnlen)1, (ftnlen)1) && ! (rcequ || lsame_(
	    equed, "N", (ftnlen)1, (ftnlen)1))) {
#line 380 "zposvx.f"
	*info = -9;
#line 381 "zposvx.f"
    } else {
#line 382 "zposvx.f"
	if (rcequ) {
#line 383 "zposvx.f"
	    smin = bignum;
#line 384 "zposvx.f"
	    smax = 0.;
#line 385 "zposvx.f"
	    i__1 = *n;
#line 385 "zposvx.f"
	    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 386 "zposvx.f"
		d__1 = smin, d__2 = s[j];
#line 386 "zposvx.f"
		smin = min(d__1,d__2);
/* Computing MAX */
#line 387 "zposvx.f"
		d__1 = smax, d__2 = s[j];
#line 387 "zposvx.f"
		smax = max(d__1,d__2);
#line 388 "zposvx.f"
/* L10: */
#line 388 "zposvx.f"
	    }
#line 389 "zposvx.f"
	    if (smin <= 0.) {
#line 390 "zposvx.f"
		*info = -10;
#line 391 "zposvx.f"
	    } else if (*n > 0) {
#line 392 "zposvx.f"
		scond = max(smin,smlnum) / min(smax,bignum);
#line 393 "zposvx.f"
	    } else {
#line 394 "zposvx.f"
		scond = 1.;
#line 395 "zposvx.f"
	    }
#line 396 "zposvx.f"
	}
#line 397 "zposvx.f"
	if (*info == 0) {
#line 398 "zposvx.f"
	    if (*ldb < max(1,*n)) {
#line 399 "zposvx.f"
		*info = -12;
#line 400 "zposvx.f"
	    } else if (*ldx < max(1,*n)) {
#line 401 "zposvx.f"
		*info = -14;
#line 402 "zposvx.f"
	    }
#line 403 "zposvx.f"
	}
#line 404 "zposvx.f"
    }

#line 406 "zposvx.f"
    if (*info != 0) {
#line 407 "zposvx.f"
	i__1 = -(*info);
#line 407 "zposvx.f"
	xerbla_("ZPOSVX", &i__1, (ftnlen)6);
#line 408 "zposvx.f"
	return 0;
#line 409 "zposvx.f"
    }

#line 411 "zposvx.f"
    if (equil) {

/*        Compute row and column scalings to equilibrate the matrix A. */

#line 415 "zposvx.f"
	zpoequ_(n, &a[a_offset], lda, &s[1], &scond, &amax, &infequ);
#line 416 "zposvx.f"
	if (infequ == 0) {

/*           Equilibrate the matrix. */

#line 420 "zposvx.f"
	    zlaqhe_(uplo, n, &a[a_offset], lda, &s[1], &scond, &amax, equed, (
		    ftnlen)1, (ftnlen)1);
#line 421 "zposvx.f"
	    rcequ = lsame_(equed, "Y", (ftnlen)1, (ftnlen)1);
#line 422 "zposvx.f"
	}
#line 423 "zposvx.f"
    }

/*     Scale the right hand side. */

#line 427 "zposvx.f"
    if (rcequ) {
#line 428 "zposvx.f"
	i__1 = *nrhs;
#line 428 "zposvx.f"
	for (j = 1; j <= i__1; ++j) {
#line 429 "zposvx.f"
	    i__2 = *n;
#line 429 "zposvx.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 430 "zposvx.f"
		i__3 = i__ + j * b_dim1;
#line 430 "zposvx.f"
		i__4 = i__;
#line 430 "zposvx.f"
		i__5 = i__ + j * b_dim1;
#line 430 "zposvx.f"
		z__1.r = s[i__4] * b[i__5].r, z__1.i = s[i__4] * b[i__5].i;
#line 430 "zposvx.f"
		b[i__3].r = z__1.r, b[i__3].i = z__1.i;
#line 431 "zposvx.f"
/* L20: */
#line 431 "zposvx.f"
	    }
#line 432 "zposvx.f"
/* L30: */
#line 432 "zposvx.f"
	}
#line 433 "zposvx.f"
    }

#line 435 "zposvx.f"
    if (nofact || equil) {

/*        Compute the Cholesky factorization A = U**H *U or A = L*L**H. */

#line 439 "zposvx.f"
	zlacpy_(uplo, n, n, &a[a_offset], lda, &af[af_offset], ldaf, (ftnlen)
		1);
#line 440 "zposvx.f"
	zpotrf_(uplo, n, &af[af_offset], ldaf, info, (ftnlen)1);

/*        Return if INFO is non-zero. */

#line 444 "zposvx.f"
	if (*info > 0) {
#line 445 "zposvx.f"
	    *rcond = 0.;
#line 446 "zposvx.f"
	    return 0;
#line 447 "zposvx.f"
	}
#line 448 "zposvx.f"
    }

/*     Compute the norm of the matrix A. */

#line 452 "zposvx.f"
    anorm = zlanhe_("1", uplo, n, &a[a_offset], lda, &rwork[1], (ftnlen)1, (
	    ftnlen)1);

/*     Compute the reciprocal of the condition number of A. */

#line 456 "zposvx.f"
    zpocon_(uplo, n, &af[af_offset], ldaf, &anorm, rcond, &work[1], &rwork[1],
	     info, (ftnlen)1);

/*     Compute the solution matrix X. */

#line 460 "zposvx.f"
    zlacpy_("Full", n, nrhs, &b[b_offset], ldb, &x[x_offset], ldx, (ftnlen)4);
#line 461 "zposvx.f"
    zpotrs_(uplo, n, nrhs, &af[af_offset], ldaf, &x[x_offset], ldx, info, (
	    ftnlen)1);

/*     Use iterative refinement to improve the computed solution and */
/*     compute error bounds and backward error estimates for it. */

#line 466 "zposvx.f"
    zporfs_(uplo, n, nrhs, &a[a_offset], lda, &af[af_offset], ldaf, &b[
	    b_offset], ldb, &x[x_offset], ldx, &ferr[1], &berr[1], &work[1], &
	    rwork[1], info, (ftnlen)1);

/*     Transform the solution matrix X to a solution of the original */
/*     system. */

#line 472 "zposvx.f"
    if (rcequ) {
#line 473 "zposvx.f"
	i__1 = *nrhs;
#line 473 "zposvx.f"
	for (j = 1; j <= i__1; ++j) {
#line 474 "zposvx.f"
	    i__2 = *n;
#line 474 "zposvx.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 475 "zposvx.f"
		i__3 = i__ + j * x_dim1;
#line 475 "zposvx.f"
		i__4 = i__;
#line 475 "zposvx.f"
		i__5 = i__ + j * x_dim1;
#line 475 "zposvx.f"
		z__1.r = s[i__4] * x[i__5].r, z__1.i = s[i__4] * x[i__5].i;
#line 475 "zposvx.f"
		x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 476 "zposvx.f"
/* L40: */
#line 476 "zposvx.f"
	    }
#line 477 "zposvx.f"
/* L50: */
#line 477 "zposvx.f"
	}
#line 478 "zposvx.f"
	i__1 = *nrhs;
#line 478 "zposvx.f"
	for (j = 1; j <= i__1; ++j) {
#line 479 "zposvx.f"
	    ferr[j] /= scond;
#line 480 "zposvx.f"
/* L60: */
#line 480 "zposvx.f"
	}
#line 481 "zposvx.f"
    }

/*     Set INFO = N+1 if the matrix is singular to working precision. */

#line 485 "zposvx.f"
    if (*rcond < dlamch_("Epsilon", (ftnlen)7)) {
#line 485 "zposvx.f"
	*info = *n + 1;
#line 485 "zposvx.f"
    }

#line 488 "zposvx.f"
    return 0;

/*     End of ZPOSVX */

} /* zposvx_ */

