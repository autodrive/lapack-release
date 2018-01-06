#line 1 "cppsvx.f"
/* cppsvx.f -- translated by f2c (version 20100827).
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

#line 1 "cppsvx.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief <b> CPPSVX computes the solution to system of linear equations A * X = B for OTHER matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CPPSVX + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cppsvx.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cppsvx.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cppsvx.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CPPSVX( FACT, UPLO, N, NRHS, AP, AFP, EQUED, S, B, LDB, */
/*                          X, LDX, RCOND, FERR, BERR, WORK, RWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          EQUED, FACT, UPLO */
/*       INTEGER            INFO, LDB, LDX, N, NRHS */
/*       REAL               RCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               BERR( * ), FERR( * ), RWORK( * ), S( * ) */
/*       COMPLEX            AFP( * ), AP( * ), B( LDB, * ), WORK( * ), */
/*      $                   X( LDX, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CPPSVX uses the Cholesky factorization A = U**H*U or A = L*L**H to */
/* > compute the solution to a complex system of linear equations */
/* >    A * X = B, */
/* > where A is an N-by-N Hermitian positive definite matrix stored in */
/* > packed format and X and B are N-by-NRHS matrices. */
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
/* >       A = U**H * U ,  if UPLO = 'U', or */
/* >       A = L * L**H,  if UPLO = 'L', */
/* >    where U is an upper triangular matrix, L is a lower triangular */
/* >    matrix, and **H indicates conjugate transpose. */
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
/* >          = 'F':  On entry, AFP contains the factored form of A. */
/* >                  If EQUED = 'Y', the matrix A has been equilibrated */
/* >                  with scaling factors given by S.  AP and AFP will not */
/* >                  be modified. */
/* >          = 'N':  The matrix A will be copied to AFP and factored. */
/* >          = 'E':  The matrix A will be equilibrated if necessary, then */
/* >                  copied to AFP and factored. */
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
/* > \param[in,out] AP */
/* > \verbatim */
/* >          AP is COMPLEX array, dimension (N*(N+1)/2) */
/* >          On entry, the upper or lower triangle of the Hermitian matrix */
/* >          A, packed columnwise in a linear array, except if FACT = 'F' */
/* >          and EQUED = 'Y', then A must contain the equilibrated matrix */
/* >          diag(S)*A*diag(S).  The j-th column of A is stored in the */
/* >          array AP as follows: */
/* >          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j; */
/* >          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n. */
/* >          See below for further details.  A is not modified if */
/* >          FACT = 'F' or 'N', or if FACT = 'E' and EQUED = 'N' on exit. */
/* > */
/* >          On exit, if FACT = 'E' and EQUED = 'Y', A is overwritten by */
/* >          diag(S)*A*diag(S). */
/* > \endverbatim */
/* > */
/* > \param[in,out] AFP */
/* > \verbatim */
/* >          AFP is COMPLEX array, dimension (N*(N+1)/2) */
/* >          If FACT = 'F', then AFP is an input argument and on entry */
/* >          contains the triangular factor U or L from the Cholesky */
/* >          factorization A = U**H*U or A = L*L**H, in the same storage */
/* >          format as A.  If EQUED .ne. 'N', then AFP is the factored */
/* >          form of the equilibrated matrix A. */
/* > */
/* >          If FACT = 'N', then AFP is an output argument and on exit */
/* >          returns the triangular factor U or L from the Cholesky */
/* >          factorization A = U**H * U or A = L * L**H of the original */
/* >          matrix A. */
/* > */
/* >          If FACT = 'E', then AFP is an output argument and on exit */
/* >          returns the triangular factor U or L from the Cholesky */
/* >          factorization A = U**H*U or A = L*L**H of the equilibrated */
/* >          matrix A (see the description of AP for the form of the */
/* >          equilibrated matrix). */
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
/* >          B is COMPLEX array, dimension (LDB,NRHS) */
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
/* >          X is COMPLEX array, dimension (LDX,NRHS) */
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

/* > \ingroup complexOTHERsolve */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  The packed storage scheme is illustrated by the following example */
/* >  when N = 4, UPLO = 'U': */
/* > */
/* >  Two-dimensional storage of the Hermitian matrix A: */
/* > */
/* >     a11 a12 a13 a14 */
/* >         a22 a23 a24 */
/* >             a33 a34     (aij = conjg(aji)) */
/* >                 a44 */
/* > */
/* >  Packed storage of the upper triangle of A: */
/* > */
/* >  AP = [ a11, a12, a22, a13, a23, a33, a14, a24, a34, a44 ] */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int cppsvx_(char *fact, char *uplo, integer *n, integer *
	nrhs, doublecomplex *ap, doublecomplex *afp, char *equed, doublereal *
	s, doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx, 
	doublereal *rcond, doublereal *ferr, doublereal *berr, doublecomplex *
	work, doublereal *rwork, integer *info, ftnlen fact_len, ftnlen 
	uplo_len, ftnlen equed_len)
{
    /* System generated locals */
    integer b_dim1, b_offset, x_dim1, x_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2;
    doublecomplex z__1;

    /* Local variables */
    static integer i__, j;
    static doublereal amax, smin, smax;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal scond, anorm;
    extern /* Subroutine */ int ccopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    static logical equil, rcequ;
    extern doublereal clanhp_(char *, char *, integer *, doublecomplex *, 
	    doublereal *, ftnlen, ftnlen), slamch_(char *, ftnlen);
    extern /* Subroutine */ int claqhp_(char *, integer *, doublecomplex *, 
	    doublereal *, doublereal *, doublereal *, char *, ftnlen, ftnlen);
    static logical nofact;
    extern /* Subroutine */ int clacpy_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    static doublereal bignum;
    extern /* Subroutine */ int cppcon_(char *, integer *, doublecomplex *, 
	    doublereal *, doublereal *, doublecomplex *, doublereal *, 
	    integer *, ftnlen);
    static integer infequ;
    extern /* Subroutine */ int cppequ_(char *, integer *, doublecomplex *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    cpprfs_(char *, integer *, integer *, doublecomplex *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublereal *, doublereal *, doublecomplex *, 
	    doublereal *, integer *, ftnlen), cpptrf_(char *, integer *, 
	    doublecomplex *, integer *, ftnlen);
    static doublereal smlnum;
    extern /* Subroutine */ int cpptrs_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, integer *, ftnlen);


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

#line 355 "cppsvx.f"
    /* Parameter adjustments */
#line 355 "cppsvx.f"
    --ap;
#line 355 "cppsvx.f"
    --afp;
#line 355 "cppsvx.f"
    --s;
#line 355 "cppsvx.f"
    b_dim1 = *ldb;
#line 355 "cppsvx.f"
    b_offset = 1 + b_dim1;
#line 355 "cppsvx.f"
    b -= b_offset;
#line 355 "cppsvx.f"
    x_dim1 = *ldx;
#line 355 "cppsvx.f"
    x_offset = 1 + x_dim1;
#line 355 "cppsvx.f"
    x -= x_offset;
#line 355 "cppsvx.f"
    --ferr;
#line 355 "cppsvx.f"
    --berr;
#line 355 "cppsvx.f"
    --work;
#line 355 "cppsvx.f"
    --rwork;
#line 355 "cppsvx.f"

#line 355 "cppsvx.f"
    /* Function Body */
#line 355 "cppsvx.f"
    *info = 0;
#line 356 "cppsvx.f"
    nofact = lsame_(fact, "N", (ftnlen)1, (ftnlen)1);
#line 357 "cppsvx.f"
    equil = lsame_(fact, "E", (ftnlen)1, (ftnlen)1);
#line 358 "cppsvx.f"
    if (nofact || equil) {
#line 359 "cppsvx.f"
	*(unsigned char *)equed = 'N';
#line 360 "cppsvx.f"
	rcequ = FALSE_;
#line 361 "cppsvx.f"
    } else {
#line 362 "cppsvx.f"
	rcequ = lsame_(equed, "Y", (ftnlen)1, (ftnlen)1);
#line 363 "cppsvx.f"
	smlnum = slamch_("Safe minimum", (ftnlen)12);
#line 364 "cppsvx.f"
	bignum = 1. / smlnum;
#line 365 "cppsvx.f"
    }

/*     Test the input parameters. */

#line 369 "cppsvx.f"
    if (! nofact && ! equil && ! lsame_(fact, "F", (ftnlen)1, (ftnlen)1)) {
#line 371 "cppsvx.f"
	*info = -1;
#line 372 "cppsvx.f"
    } else if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, 
	    "L", (ftnlen)1, (ftnlen)1)) {
#line 374 "cppsvx.f"
	*info = -2;
#line 375 "cppsvx.f"
    } else if (*n < 0) {
#line 376 "cppsvx.f"
	*info = -3;
#line 377 "cppsvx.f"
    } else if (*nrhs < 0) {
#line 378 "cppsvx.f"
	*info = -4;
#line 379 "cppsvx.f"
    } else if (lsame_(fact, "F", (ftnlen)1, (ftnlen)1) && ! (rcequ || lsame_(
	    equed, "N", (ftnlen)1, (ftnlen)1))) {
#line 381 "cppsvx.f"
	*info = -7;
#line 382 "cppsvx.f"
    } else {
#line 383 "cppsvx.f"
	if (rcequ) {
#line 384 "cppsvx.f"
	    smin = bignum;
#line 385 "cppsvx.f"
	    smax = 0.;
#line 386 "cppsvx.f"
	    i__1 = *n;
#line 386 "cppsvx.f"
	    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 387 "cppsvx.f"
		d__1 = smin, d__2 = s[j];
#line 387 "cppsvx.f"
		smin = min(d__1,d__2);
/* Computing MAX */
#line 388 "cppsvx.f"
		d__1 = smax, d__2 = s[j];
#line 388 "cppsvx.f"
		smax = max(d__1,d__2);
#line 389 "cppsvx.f"
/* L10: */
#line 389 "cppsvx.f"
	    }
#line 390 "cppsvx.f"
	    if (smin <= 0.) {
#line 391 "cppsvx.f"
		*info = -8;
#line 392 "cppsvx.f"
	    } else if (*n > 0) {
#line 393 "cppsvx.f"
		scond = max(smin,smlnum) / min(smax,bignum);
#line 394 "cppsvx.f"
	    } else {
#line 395 "cppsvx.f"
		scond = 1.;
#line 396 "cppsvx.f"
	    }
#line 397 "cppsvx.f"
	}
#line 398 "cppsvx.f"
	if (*info == 0) {
#line 399 "cppsvx.f"
	    if (*ldb < max(1,*n)) {
#line 400 "cppsvx.f"
		*info = -10;
#line 401 "cppsvx.f"
	    } else if (*ldx < max(1,*n)) {
#line 402 "cppsvx.f"
		*info = -12;
#line 403 "cppsvx.f"
	    }
#line 404 "cppsvx.f"
	}
#line 405 "cppsvx.f"
    }

#line 407 "cppsvx.f"
    if (*info != 0) {
#line 408 "cppsvx.f"
	i__1 = -(*info);
#line 408 "cppsvx.f"
	xerbla_("CPPSVX", &i__1, (ftnlen)6);
#line 409 "cppsvx.f"
	return 0;
#line 410 "cppsvx.f"
    }

#line 412 "cppsvx.f"
    if (equil) {

/*        Compute row and column scalings to equilibrate the matrix A. */

#line 416 "cppsvx.f"
	cppequ_(uplo, n, &ap[1], &s[1], &scond, &amax, &infequ, (ftnlen)1);
#line 417 "cppsvx.f"
	if (infequ == 0) {

/*           Equilibrate the matrix. */

#line 421 "cppsvx.f"
	    claqhp_(uplo, n, &ap[1], &s[1], &scond, &amax, equed, (ftnlen)1, (
		    ftnlen)1);
#line 422 "cppsvx.f"
	    rcequ = lsame_(equed, "Y", (ftnlen)1, (ftnlen)1);
#line 423 "cppsvx.f"
	}
#line 424 "cppsvx.f"
    }

/*     Scale the right-hand side. */

#line 428 "cppsvx.f"
    if (rcequ) {
#line 429 "cppsvx.f"
	i__1 = *nrhs;
#line 429 "cppsvx.f"
	for (j = 1; j <= i__1; ++j) {
#line 430 "cppsvx.f"
	    i__2 = *n;
#line 430 "cppsvx.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 431 "cppsvx.f"
		i__3 = i__ + j * b_dim1;
#line 431 "cppsvx.f"
		i__4 = i__;
#line 431 "cppsvx.f"
		i__5 = i__ + j * b_dim1;
#line 431 "cppsvx.f"
		z__1.r = s[i__4] * b[i__5].r, z__1.i = s[i__4] * b[i__5].i;
#line 431 "cppsvx.f"
		b[i__3].r = z__1.r, b[i__3].i = z__1.i;
#line 432 "cppsvx.f"
/* L20: */
#line 432 "cppsvx.f"
	    }
#line 433 "cppsvx.f"
/* L30: */
#line 433 "cppsvx.f"
	}
#line 434 "cppsvx.f"
    }

#line 436 "cppsvx.f"
    if (nofact || equil) {

/*        Compute the Cholesky factorization A = U**H * U or A = L * L**H. */

#line 440 "cppsvx.f"
	i__1 = *n * (*n + 1) / 2;
#line 440 "cppsvx.f"
	ccopy_(&i__1, &ap[1], &c__1, &afp[1], &c__1);
#line 441 "cppsvx.f"
	cpptrf_(uplo, n, &afp[1], info, (ftnlen)1);

/*        Return if INFO is non-zero. */

#line 445 "cppsvx.f"
	if (*info > 0) {
#line 446 "cppsvx.f"
	    *rcond = 0.;
#line 447 "cppsvx.f"
	    return 0;
#line 448 "cppsvx.f"
	}
#line 449 "cppsvx.f"
    }

/*     Compute the norm of the matrix A. */

#line 453 "cppsvx.f"
    anorm = clanhp_("I", uplo, n, &ap[1], &rwork[1], (ftnlen)1, (ftnlen)1);

/*     Compute the reciprocal of the condition number of A. */

#line 457 "cppsvx.f"
    cppcon_(uplo, n, &afp[1], &anorm, rcond, &work[1], &rwork[1], info, (
	    ftnlen)1);

/*     Compute the solution matrix X. */

#line 461 "cppsvx.f"
    clacpy_("Full", n, nrhs, &b[b_offset], ldb, &x[x_offset], ldx, (ftnlen)4);
#line 462 "cppsvx.f"
    cpptrs_(uplo, n, nrhs, &afp[1], &x[x_offset], ldx, info, (ftnlen)1);

/*     Use iterative refinement to improve the computed solution and */
/*     compute error bounds and backward error estimates for it. */

#line 467 "cppsvx.f"
    cpprfs_(uplo, n, nrhs, &ap[1], &afp[1], &b[b_offset], ldb, &x[x_offset], 
	    ldx, &ferr[1], &berr[1], &work[1], &rwork[1], info, (ftnlen)1);

/*     Transform the solution matrix X to a solution of the original */
/*     system. */

#line 473 "cppsvx.f"
    if (rcequ) {
#line 474 "cppsvx.f"
	i__1 = *nrhs;
#line 474 "cppsvx.f"
	for (j = 1; j <= i__1; ++j) {
#line 475 "cppsvx.f"
	    i__2 = *n;
#line 475 "cppsvx.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 476 "cppsvx.f"
		i__3 = i__ + j * x_dim1;
#line 476 "cppsvx.f"
		i__4 = i__;
#line 476 "cppsvx.f"
		i__5 = i__ + j * x_dim1;
#line 476 "cppsvx.f"
		z__1.r = s[i__4] * x[i__5].r, z__1.i = s[i__4] * x[i__5].i;
#line 476 "cppsvx.f"
		x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 477 "cppsvx.f"
/* L40: */
#line 477 "cppsvx.f"
	    }
#line 478 "cppsvx.f"
/* L50: */
#line 478 "cppsvx.f"
	}
#line 479 "cppsvx.f"
	i__1 = *nrhs;
#line 479 "cppsvx.f"
	for (j = 1; j <= i__1; ++j) {
#line 480 "cppsvx.f"
	    ferr[j] /= scond;
#line 481 "cppsvx.f"
/* L60: */
#line 481 "cppsvx.f"
	}
#line 482 "cppsvx.f"
    }

/*     Set INFO = N+1 if the matrix is singular to working precision. */

#line 486 "cppsvx.f"
    if (*rcond < slamch_("Epsilon", (ftnlen)7)) {
#line 486 "cppsvx.f"
	*info = *n + 1;
#line 486 "cppsvx.f"
    }

#line 489 "cppsvx.f"
    return 0;

/*     End of CPPSVX */

} /* cppsvx_ */

