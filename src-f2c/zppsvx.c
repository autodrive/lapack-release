#line 1 "zppsvx.f"
/* zppsvx.f -- translated by f2c (version 20100827).
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

#line 1 "zppsvx.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief <b> ZPPSVX computes the solution to system of linear equations A * X = B for OTHER matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZPPSVX + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zppsvx.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zppsvx.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zppsvx.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZPPSVX( FACT, UPLO, N, NRHS, AP, AFP, EQUED, S, B, LDB, */
/*                          X, LDX, RCOND, FERR, BERR, WORK, RWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          EQUED, FACT, UPLO */
/*       INTEGER            INFO, LDB, LDX, N, NRHS */
/*       DOUBLE PRECISION   RCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   BERR( * ), FERR( * ), RWORK( * ), S( * ) */
/*       COMPLEX*16         AFP( * ), AP( * ), B( LDB, * ), WORK( * ), */
/*      $                   X( LDX, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZPPSVX uses the Cholesky factorization A = U**H * U or A = L * L**H to */
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
/* >          AP is COMPLEX*16 array, dimension (N*(N+1)/2) */
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
/* >          AFP is COMPLEX*16 array, dimension (N*(N+1)/2) */
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
/* >          factorization A = U**H * U or A = L * L**H of the equilibrated */
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

/* > \ingroup complex16OTHERsolve */

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
/* Subroutine */ int zppsvx_(char *fact, char *uplo, integer *n, integer *
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
    static logical equil, rcequ;
    extern /* Subroutine */ int zcopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    extern doublereal dlamch_(char *, ftnlen);
    static logical nofact;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal bignum;
    static integer infequ;
    extern doublereal zlanhp_(char *, char *, integer *, doublecomplex *, 
	    doublereal *, ftnlen, ftnlen);
    extern /* Subroutine */ int zlaqhp_(char *, integer *, doublecomplex *, 
	    doublereal *, doublereal *, doublereal *, char *, ftnlen, ftnlen),
	     zlacpy_(char *, integer *, integer *, doublecomplex *, integer *,
	     doublecomplex *, integer *, ftnlen), zppcon_(char *, integer *, 
	    doublecomplex *, doublereal *, doublereal *, doublecomplex *, 
	    doublereal *, integer *, ftnlen);
    static doublereal smlnum;
    extern /* Subroutine */ int zppequ_(char *, integer *, doublecomplex *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    zpprfs_(char *, integer *, integer *, doublecomplex *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublereal *, doublereal *, doublecomplex *, 
	    doublereal *, integer *, ftnlen), zpptrf_(char *, integer *, 
	    doublecomplex *, integer *, ftnlen), zpptrs_(char *, integer *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, integer *,
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

#line 355 "zppsvx.f"
    /* Parameter adjustments */
#line 355 "zppsvx.f"
    --ap;
#line 355 "zppsvx.f"
    --afp;
#line 355 "zppsvx.f"
    --s;
#line 355 "zppsvx.f"
    b_dim1 = *ldb;
#line 355 "zppsvx.f"
    b_offset = 1 + b_dim1;
#line 355 "zppsvx.f"
    b -= b_offset;
#line 355 "zppsvx.f"
    x_dim1 = *ldx;
#line 355 "zppsvx.f"
    x_offset = 1 + x_dim1;
#line 355 "zppsvx.f"
    x -= x_offset;
#line 355 "zppsvx.f"
    --ferr;
#line 355 "zppsvx.f"
    --berr;
#line 355 "zppsvx.f"
    --work;
#line 355 "zppsvx.f"
    --rwork;
#line 355 "zppsvx.f"

#line 355 "zppsvx.f"
    /* Function Body */
#line 355 "zppsvx.f"
    *info = 0;
#line 356 "zppsvx.f"
    nofact = lsame_(fact, "N", (ftnlen)1, (ftnlen)1);
#line 357 "zppsvx.f"
    equil = lsame_(fact, "E", (ftnlen)1, (ftnlen)1);
#line 358 "zppsvx.f"
    if (nofact || equil) {
#line 359 "zppsvx.f"
	*(unsigned char *)equed = 'N';
#line 360 "zppsvx.f"
	rcequ = FALSE_;
#line 361 "zppsvx.f"
    } else {
#line 362 "zppsvx.f"
	rcequ = lsame_(equed, "Y", (ftnlen)1, (ftnlen)1);
#line 363 "zppsvx.f"
	smlnum = dlamch_("Safe minimum", (ftnlen)12);
#line 364 "zppsvx.f"
	bignum = 1. / smlnum;
#line 365 "zppsvx.f"
    }

/*     Test the input parameters. */

#line 369 "zppsvx.f"
    if (! nofact && ! equil && ! lsame_(fact, "F", (ftnlen)1, (ftnlen)1)) {
#line 371 "zppsvx.f"
	*info = -1;
#line 372 "zppsvx.f"
    } else if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, 
	    "L", (ftnlen)1, (ftnlen)1)) {
#line 374 "zppsvx.f"
	*info = -2;
#line 375 "zppsvx.f"
    } else if (*n < 0) {
#line 376 "zppsvx.f"
	*info = -3;
#line 377 "zppsvx.f"
    } else if (*nrhs < 0) {
#line 378 "zppsvx.f"
	*info = -4;
#line 379 "zppsvx.f"
    } else if (lsame_(fact, "F", (ftnlen)1, (ftnlen)1) && ! (rcequ || lsame_(
	    equed, "N", (ftnlen)1, (ftnlen)1))) {
#line 381 "zppsvx.f"
	*info = -7;
#line 382 "zppsvx.f"
    } else {
#line 383 "zppsvx.f"
	if (rcequ) {
#line 384 "zppsvx.f"
	    smin = bignum;
#line 385 "zppsvx.f"
	    smax = 0.;
#line 386 "zppsvx.f"
	    i__1 = *n;
#line 386 "zppsvx.f"
	    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 387 "zppsvx.f"
		d__1 = smin, d__2 = s[j];
#line 387 "zppsvx.f"
		smin = min(d__1,d__2);
/* Computing MAX */
#line 388 "zppsvx.f"
		d__1 = smax, d__2 = s[j];
#line 388 "zppsvx.f"
		smax = max(d__1,d__2);
#line 389 "zppsvx.f"
/* L10: */
#line 389 "zppsvx.f"
	    }
#line 390 "zppsvx.f"
	    if (smin <= 0.) {
#line 391 "zppsvx.f"
		*info = -8;
#line 392 "zppsvx.f"
	    } else if (*n > 0) {
#line 393 "zppsvx.f"
		scond = max(smin,smlnum) / min(smax,bignum);
#line 394 "zppsvx.f"
	    } else {
#line 395 "zppsvx.f"
		scond = 1.;
#line 396 "zppsvx.f"
	    }
#line 397 "zppsvx.f"
	}
#line 398 "zppsvx.f"
	if (*info == 0) {
#line 399 "zppsvx.f"
	    if (*ldb < max(1,*n)) {
#line 400 "zppsvx.f"
		*info = -10;
#line 401 "zppsvx.f"
	    } else if (*ldx < max(1,*n)) {
#line 402 "zppsvx.f"
		*info = -12;
#line 403 "zppsvx.f"
	    }
#line 404 "zppsvx.f"
	}
#line 405 "zppsvx.f"
    }

#line 407 "zppsvx.f"
    if (*info != 0) {
#line 408 "zppsvx.f"
	i__1 = -(*info);
#line 408 "zppsvx.f"
	xerbla_("ZPPSVX", &i__1, (ftnlen)6);
#line 409 "zppsvx.f"
	return 0;
#line 410 "zppsvx.f"
    }

#line 412 "zppsvx.f"
    if (equil) {

/*        Compute row and column scalings to equilibrate the matrix A. */

#line 416 "zppsvx.f"
	zppequ_(uplo, n, &ap[1], &s[1], &scond, &amax, &infequ, (ftnlen)1);
#line 417 "zppsvx.f"
	if (infequ == 0) {

/*           Equilibrate the matrix. */

#line 421 "zppsvx.f"
	    zlaqhp_(uplo, n, &ap[1], &s[1], &scond, &amax, equed, (ftnlen)1, (
		    ftnlen)1);
#line 422 "zppsvx.f"
	    rcequ = lsame_(equed, "Y", (ftnlen)1, (ftnlen)1);
#line 423 "zppsvx.f"
	}
#line 424 "zppsvx.f"
    }

/*     Scale the right-hand side. */

#line 428 "zppsvx.f"
    if (rcequ) {
#line 429 "zppsvx.f"
	i__1 = *nrhs;
#line 429 "zppsvx.f"
	for (j = 1; j <= i__1; ++j) {
#line 430 "zppsvx.f"
	    i__2 = *n;
#line 430 "zppsvx.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 431 "zppsvx.f"
		i__3 = i__ + j * b_dim1;
#line 431 "zppsvx.f"
		i__4 = i__;
#line 431 "zppsvx.f"
		i__5 = i__ + j * b_dim1;
#line 431 "zppsvx.f"
		z__1.r = s[i__4] * b[i__5].r, z__1.i = s[i__4] * b[i__5].i;
#line 431 "zppsvx.f"
		b[i__3].r = z__1.r, b[i__3].i = z__1.i;
#line 432 "zppsvx.f"
/* L20: */
#line 432 "zppsvx.f"
	    }
#line 433 "zppsvx.f"
/* L30: */
#line 433 "zppsvx.f"
	}
#line 434 "zppsvx.f"
    }

#line 436 "zppsvx.f"
    if (nofact || equil) {

/*        Compute the Cholesky factorization A = U**H * U or A = L * L**H. */

#line 440 "zppsvx.f"
	i__1 = *n * (*n + 1) / 2;
#line 440 "zppsvx.f"
	zcopy_(&i__1, &ap[1], &c__1, &afp[1], &c__1);
#line 441 "zppsvx.f"
	zpptrf_(uplo, n, &afp[1], info, (ftnlen)1);

/*        Return if INFO is non-zero. */

#line 445 "zppsvx.f"
	if (*info > 0) {
#line 446 "zppsvx.f"
	    *rcond = 0.;
#line 447 "zppsvx.f"
	    return 0;
#line 448 "zppsvx.f"
	}
#line 449 "zppsvx.f"
    }

/*     Compute the norm of the matrix A. */

#line 453 "zppsvx.f"
    anorm = zlanhp_("I", uplo, n, &ap[1], &rwork[1], (ftnlen)1, (ftnlen)1);

/*     Compute the reciprocal of the condition number of A. */

#line 457 "zppsvx.f"
    zppcon_(uplo, n, &afp[1], &anorm, rcond, &work[1], &rwork[1], info, (
	    ftnlen)1);

/*     Compute the solution matrix X. */

#line 461 "zppsvx.f"
    zlacpy_("Full", n, nrhs, &b[b_offset], ldb, &x[x_offset], ldx, (ftnlen)4);
#line 462 "zppsvx.f"
    zpptrs_(uplo, n, nrhs, &afp[1], &x[x_offset], ldx, info, (ftnlen)1);

/*     Use iterative refinement to improve the computed solution and */
/*     compute error bounds and backward error estimates for it. */

#line 467 "zppsvx.f"
    zpprfs_(uplo, n, nrhs, &ap[1], &afp[1], &b[b_offset], ldb, &x[x_offset], 
	    ldx, &ferr[1], &berr[1], &work[1], &rwork[1], info, (ftnlen)1);

/*     Transform the solution matrix X to a solution of the original */
/*     system. */

#line 473 "zppsvx.f"
    if (rcequ) {
#line 474 "zppsvx.f"
	i__1 = *nrhs;
#line 474 "zppsvx.f"
	for (j = 1; j <= i__1; ++j) {
#line 475 "zppsvx.f"
	    i__2 = *n;
#line 475 "zppsvx.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 476 "zppsvx.f"
		i__3 = i__ + j * x_dim1;
#line 476 "zppsvx.f"
		i__4 = i__;
#line 476 "zppsvx.f"
		i__5 = i__ + j * x_dim1;
#line 476 "zppsvx.f"
		z__1.r = s[i__4] * x[i__5].r, z__1.i = s[i__4] * x[i__5].i;
#line 476 "zppsvx.f"
		x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 477 "zppsvx.f"
/* L40: */
#line 477 "zppsvx.f"
	    }
#line 478 "zppsvx.f"
/* L50: */
#line 478 "zppsvx.f"
	}
#line 479 "zppsvx.f"
	i__1 = *nrhs;
#line 479 "zppsvx.f"
	for (j = 1; j <= i__1; ++j) {
#line 480 "zppsvx.f"
	    ferr[j] /= scond;
#line 481 "zppsvx.f"
/* L60: */
#line 481 "zppsvx.f"
	}
#line 482 "zppsvx.f"
    }

/*     Set INFO = N+1 if the matrix is singular to working precision. */

#line 486 "zppsvx.f"
    if (*rcond < dlamch_("Epsilon", (ftnlen)7)) {
#line 486 "zppsvx.f"
	*info = *n + 1;
#line 486 "zppsvx.f"
    }

#line 489 "zppsvx.f"
    return 0;

/*     End of ZPPSVX */

} /* zppsvx_ */

