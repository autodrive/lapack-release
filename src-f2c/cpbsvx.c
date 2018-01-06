#line 1 "cpbsvx.f"
/* cpbsvx.f -- translated by f2c (version 20100827).
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

#line 1 "cpbsvx.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief <b> CPBSVX computes the solution to system of linear equations A * X = B for OTHER matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CPBSVX + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cpbsvx.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cpbsvx.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cpbsvx.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CPBSVX( FACT, UPLO, N, KD, NRHS, AB, LDAB, AFB, LDAFB, */
/*                          EQUED, S, B, LDB, X, LDX, RCOND, FERR, BERR, */
/*                          WORK, RWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          EQUED, FACT, UPLO */
/*       INTEGER            INFO, KD, LDAB, LDAFB, LDB, LDX, N, NRHS */
/*       REAL               RCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               BERR( * ), FERR( * ), RWORK( * ), S( * ) */
/*       COMPLEX            AB( LDAB, * ), AFB( LDAFB, * ), B( LDB, * ), */
/*      $                   WORK( * ), X( LDX, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CPBSVX uses the Cholesky factorization A = U**H*U or A = L*L**H to */
/* > compute the solution to a complex system of linear equations */
/* >    A * X = B, */
/* > where A is an N-by-N Hermitian positive definite band matrix and X */
/* > and B are N-by-NRHS matrices. */
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
/* >       A = U**H * U,  if UPLO = 'U', or */
/* >       A = L * L**H,  if UPLO = 'L', */
/* >    where U is an upper triangular band matrix, and L is a lower */
/* >    triangular band matrix. */
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
/* >          = 'F':  On entry, AFB contains the factored form of A. */
/* >                  If EQUED = 'Y', the matrix A has been equilibrated */
/* >                  with scaling factors given by S.  AB and AFB will not */
/* >                  be modified. */
/* >          = 'N':  The matrix A will be copied to AFB and factored. */
/* >          = 'E':  The matrix A will be equilibrated if necessary, then */
/* >                  copied to AFB and factored. */
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
/* > \param[in] KD */
/* > \verbatim */
/* >          KD is INTEGER */
/* >          The number of superdiagonals of the matrix A if UPLO = 'U', */
/* >          or the number of subdiagonals if UPLO = 'L'.  KD >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NRHS */
/* > \verbatim */
/* >          NRHS is INTEGER */
/* >          The number of right-hand sides, i.e., the number of columns */
/* >          of the matrices B and X.  NRHS >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] AB */
/* > \verbatim */
/* >          AB is COMPLEX array, dimension (LDAB,N) */
/* >          On entry, the upper or lower triangle of the Hermitian band */
/* >          matrix A, stored in the first KD+1 rows of the array, except */
/* >          if FACT = 'F' and EQUED = 'Y', then A must contain the */
/* >          equilibrated matrix diag(S)*A*diag(S).  The j-th column of A */
/* >          is stored in the j-th column of the array AB as follows: */
/* >          if UPLO = 'U', AB(KD+1+i-j,j) = A(i,j) for max(1,j-KD)<=i<=j; */
/* >          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(N,j+KD). */
/* >          See below for further details. */
/* > */
/* >          On exit, if FACT = 'E' and EQUED = 'Y', A is overwritten by */
/* >          diag(S)*A*diag(S). */
/* > \endverbatim */
/* > */
/* > \param[in] LDAB */
/* > \verbatim */
/* >          LDAB is INTEGER */
/* >          The leading dimension of the array A.  LDAB >= KD+1. */
/* > \endverbatim */
/* > */
/* > \param[in,out] AFB */
/* > \verbatim */
/* >          AFB is COMPLEX array, dimension (LDAFB,N) */
/* >          If FACT = 'F', then AFB is an input argument and on entry */
/* >          contains the triangular factor U or L from the Cholesky */
/* >          factorization A = U**H*U or A = L*L**H of the band matrix */
/* >          A, in the same storage format as A (see AB).  If EQUED = 'Y', */
/* >          then AFB is the factored form of the equilibrated matrix A. */
/* > */
/* >          If FACT = 'N', then AFB is an output argument and on exit */
/* >          returns the triangular factor U or L from the Cholesky */
/* >          factorization A = U**H*U or A = L*L**H. */
/* > */
/* >          If FACT = 'E', then AFB is an output argument and on exit */
/* >          returns the triangular factor U or L from the Cholesky */
/* >          factorization A = U**H*U or A = L*L**H of the equilibrated */
/* >          matrix A (see the description of A for the form of the */
/* >          equilibrated matrix). */
/* > \endverbatim */
/* > */
/* > \param[in] LDAFB */
/* > \verbatim */
/* >          LDAFB is INTEGER */
/* >          The leading dimension of the array AFB.  LDAFB >= KD+1. */
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

/* > \ingroup complexOTHERsolve */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  The band storage scheme is illustrated by the following example, when */
/* >  N = 6, KD = 2, and UPLO = 'U': */
/* > */
/* >  Two-dimensional storage of the Hermitian matrix A: */
/* > */
/* >     a11  a12  a13 */
/* >          a22  a23  a24 */
/* >               a33  a34  a35 */
/* >                    a44  a45  a46 */
/* >                         a55  a56 */
/* >     (aij=conjg(aji))         a66 */
/* > */
/* >  Band storage of the upper triangle of A: */
/* > */
/* >      *    *   a13  a24  a35  a46 */
/* >      *   a12  a23  a34  a45  a56 */
/* >     a11  a22  a33  a44  a55  a66 */
/* > */
/* >  Similarly, if UPLO = 'L' the format of A is as follows: */
/* > */
/* >     a11  a22  a33  a44  a55  a66 */
/* >     a21  a32  a43  a54  a65   * */
/* >     a31  a42  a53  a64   *    * */
/* > */
/* >  Array elements marked * are not used by the routine. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int cpbsvx_(char *fact, char *uplo, integer *n, integer *kd, 
	integer *nrhs, doublecomplex *ab, integer *ldab, doublecomplex *afb, 
	integer *ldafb, char *equed, doublereal *s, doublecomplex *b, integer 
	*ldb, doublecomplex *x, integer *ldx, doublereal *rcond, doublereal *
	ferr, doublereal *berr, doublecomplex *work, doublereal *rwork, 
	integer *info, ftnlen fact_len, ftnlen uplo_len, ftnlen equed_len)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, afb_dim1, afb_offset, b_dim1, b_offset, 
	    x_dim1, x_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2;
    doublecomplex z__1;

    /* Local variables */
    static integer i__, j, j1, j2;
    static doublereal amax, smin, smax;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal scond, anorm;
    extern /* Subroutine */ int ccopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    static logical equil, rcequ, upper;
    extern doublereal clanhb_(char *, char *, integer *, integer *, 
	    doublecomplex *, integer *, doublereal *, ftnlen, ftnlen);
    extern /* Subroutine */ int claqhb_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublereal *, doublereal *, 
	    doublereal *, char *, ftnlen, ftnlen), cpbcon_(char *, integer *, 
	    integer *, doublecomplex *, integer *, doublereal *, doublereal *,
	     doublecomplex *, doublereal *, integer *, ftnlen);
    extern doublereal slamch_(char *, ftnlen);
    static logical nofact;
    extern /* Subroutine */ int clacpy_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen), cpbequ_(char *, integer *, 
	    integer *, doublecomplex *, integer *, doublereal *, doublereal *,
	     doublereal *, integer *, ftnlen), cpbrfs_(char *, integer *, 
	    integer *, integer *, doublecomplex *, integer *, doublecomplex *,
	     integer *, doublecomplex *, integer *, doublecomplex *, integer *
	    , doublereal *, doublereal *, doublecomplex *, doublereal *, 
	    integer *, ftnlen);
    static doublereal bignum;
    extern /* Subroutine */ int cpbtrf_(char *, integer *, integer *, 
	    doublecomplex *, integer *, integer *, ftnlen);
    static integer infequ;
    extern /* Subroutine */ int cpbtrs_(char *, integer *, integer *, integer 
	    *, doublecomplex *, integer *, doublecomplex *, integer *, 
	    integer *, ftnlen);
    static doublereal smlnum;


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

#line 386 "cpbsvx.f"
    /* Parameter adjustments */
#line 386 "cpbsvx.f"
    ab_dim1 = *ldab;
#line 386 "cpbsvx.f"
    ab_offset = 1 + ab_dim1;
#line 386 "cpbsvx.f"
    ab -= ab_offset;
#line 386 "cpbsvx.f"
    afb_dim1 = *ldafb;
#line 386 "cpbsvx.f"
    afb_offset = 1 + afb_dim1;
#line 386 "cpbsvx.f"
    afb -= afb_offset;
#line 386 "cpbsvx.f"
    --s;
#line 386 "cpbsvx.f"
    b_dim1 = *ldb;
#line 386 "cpbsvx.f"
    b_offset = 1 + b_dim1;
#line 386 "cpbsvx.f"
    b -= b_offset;
#line 386 "cpbsvx.f"
    x_dim1 = *ldx;
#line 386 "cpbsvx.f"
    x_offset = 1 + x_dim1;
#line 386 "cpbsvx.f"
    x -= x_offset;
#line 386 "cpbsvx.f"
    --ferr;
#line 386 "cpbsvx.f"
    --berr;
#line 386 "cpbsvx.f"
    --work;
#line 386 "cpbsvx.f"
    --rwork;
#line 386 "cpbsvx.f"

#line 386 "cpbsvx.f"
    /* Function Body */
#line 386 "cpbsvx.f"
    *info = 0;
#line 387 "cpbsvx.f"
    nofact = lsame_(fact, "N", (ftnlen)1, (ftnlen)1);
#line 388 "cpbsvx.f"
    equil = lsame_(fact, "E", (ftnlen)1, (ftnlen)1);
#line 389 "cpbsvx.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 390 "cpbsvx.f"
    if (nofact || equil) {
#line 391 "cpbsvx.f"
	*(unsigned char *)equed = 'N';
#line 392 "cpbsvx.f"
	rcequ = FALSE_;
#line 393 "cpbsvx.f"
    } else {
#line 394 "cpbsvx.f"
	rcequ = lsame_(equed, "Y", (ftnlen)1, (ftnlen)1);
#line 395 "cpbsvx.f"
	smlnum = slamch_("Safe minimum", (ftnlen)12);
#line 396 "cpbsvx.f"
	bignum = 1. / smlnum;
#line 397 "cpbsvx.f"
    }

/*     Test the input parameters. */

#line 401 "cpbsvx.f"
    if (! nofact && ! equil && ! lsame_(fact, "F", (ftnlen)1, (ftnlen)1)) {
#line 403 "cpbsvx.f"
	*info = -1;
#line 404 "cpbsvx.f"
    } else if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 405 "cpbsvx.f"
	*info = -2;
#line 406 "cpbsvx.f"
    } else if (*n < 0) {
#line 407 "cpbsvx.f"
	*info = -3;
#line 408 "cpbsvx.f"
    } else if (*kd < 0) {
#line 409 "cpbsvx.f"
	*info = -4;
#line 410 "cpbsvx.f"
    } else if (*nrhs < 0) {
#line 411 "cpbsvx.f"
	*info = -5;
#line 412 "cpbsvx.f"
    } else if (*ldab < *kd + 1) {
#line 413 "cpbsvx.f"
	*info = -7;
#line 414 "cpbsvx.f"
    } else if (*ldafb < *kd + 1) {
#line 415 "cpbsvx.f"
	*info = -9;
#line 416 "cpbsvx.f"
    } else if (lsame_(fact, "F", (ftnlen)1, (ftnlen)1) && ! (rcequ || lsame_(
	    equed, "N", (ftnlen)1, (ftnlen)1))) {
#line 418 "cpbsvx.f"
	*info = -10;
#line 419 "cpbsvx.f"
    } else {
#line 420 "cpbsvx.f"
	if (rcequ) {
#line 421 "cpbsvx.f"
	    smin = bignum;
#line 422 "cpbsvx.f"
	    smax = 0.;
#line 423 "cpbsvx.f"
	    i__1 = *n;
#line 423 "cpbsvx.f"
	    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 424 "cpbsvx.f"
		d__1 = smin, d__2 = s[j];
#line 424 "cpbsvx.f"
		smin = min(d__1,d__2);
/* Computing MAX */
#line 425 "cpbsvx.f"
		d__1 = smax, d__2 = s[j];
#line 425 "cpbsvx.f"
		smax = max(d__1,d__2);
#line 426 "cpbsvx.f"
/* L10: */
#line 426 "cpbsvx.f"
	    }
#line 427 "cpbsvx.f"
	    if (smin <= 0.) {
#line 428 "cpbsvx.f"
		*info = -11;
#line 429 "cpbsvx.f"
	    } else if (*n > 0) {
#line 430 "cpbsvx.f"
		scond = max(smin,smlnum) / min(smax,bignum);
#line 431 "cpbsvx.f"
	    } else {
#line 432 "cpbsvx.f"
		scond = 1.;
#line 433 "cpbsvx.f"
	    }
#line 434 "cpbsvx.f"
	}
#line 435 "cpbsvx.f"
	if (*info == 0) {
#line 436 "cpbsvx.f"
	    if (*ldb < max(1,*n)) {
#line 437 "cpbsvx.f"
		*info = -13;
#line 438 "cpbsvx.f"
	    } else if (*ldx < max(1,*n)) {
#line 439 "cpbsvx.f"
		*info = -15;
#line 440 "cpbsvx.f"
	    }
#line 441 "cpbsvx.f"
	}
#line 442 "cpbsvx.f"
    }

#line 444 "cpbsvx.f"
    if (*info != 0) {
#line 445 "cpbsvx.f"
	i__1 = -(*info);
#line 445 "cpbsvx.f"
	xerbla_("CPBSVX", &i__1, (ftnlen)6);
#line 446 "cpbsvx.f"
	return 0;
#line 447 "cpbsvx.f"
    }

#line 449 "cpbsvx.f"
    if (equil) {

/*        Compute row and column scalings to equilibrate the matrix A. */

#line 453 "cpbsvx.f"
	cpbequ_(uplo, n, kd, &ab[ab_offset], ldab, &s[1], &scond, &amax, &
		infequ, (ftnlen)1);
#line 454 "cpbsvx.f"
	if (infequ == 0) {

/*           Equilibrate the matrix. */

#line 458 "cpbsvx.f"
	    claqhb_(uplo, n, kd, &ab[ab_offset], ldab, &s[1], &scond, &amax, 
		    equed, (ftnlen)1, (ftnlen)1);
#line 459 "cpbsvx.f"
	    rcequ = lsame_(equed, "Y", (ftnlen)1, (ftnlen)1);
#line 460 "cpbsvx.f"
	}
#line 461 "cpbsvx.f"
    }

/*     Scale the right-hand side. */

#line 465 "cpbsvx.f"
    if (rcequ) {
#line 466 "cpbsvx.f"
	i__1 = *nrhs;
#line 466 "cpbsvx.f"
	for (j = 1; j <= i__1; ++j) {
#line 467 "cpbsvx.f"
	    i__2 = *n;
#line 467 "cpbsvx.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 468 "cpbsvx.f"
		i__3 = i__ + j * b_dim1;
#line 468 "cpbsvx.f"
		i__4 = i__;
#line 468 "cpbsvx.f"
		i__5 = i__ + j * b_dim1;
#line 468 "cpbsvx.f"
		z__1.r = s[i__4] * b[i__5].r, z__1.i = s[i__4] * b[i__5].i;
#line 468 "cpbsvx.f"
		b[i__3].r = z__1.r, b[i__3].i = z__1.i;
#line 469 "cpbsvx.f"
/* L20: */
#line 469 "cpbsvx.f"
	    }
#line 470 "cpbsvx.f"
/* L30: */
#line 470 "cpbsvx.f"
	}
#line 471 "cpbsvx.f"
    }

#line 473 "cpbsvx.f"
    if (nofact || equil) {

/*        Compute the Cholesky factorization A = U**H *U or A = L*L**H. */

#line 477 "cpbsvx.f"
	if (upper) {
#line 478 "cpbsvx.f"
	    i__1 = *n;
#line 478 "cpbsvx.f"
	    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
#line 479 "cpbsvx.f"
		i__2 = j - *kd;
#line 479 "cpbsvx.f"
		j1 = max(i__2,1);
#line 480 "cpbsvx.f"
		i__2 = j - j1 + 1;
#line 480 "cpbsvx.f"
		ccopy_(&i__2, &ab[*kd + 1 - j + j1 + j * ab_dim1], &c__1, &
			afb[*kd + 1 - j + j1 + j * afb_dim1], &c__1);
#line 482 "cpbsvx.f"
/* L40: */
#line 482 "cpbsvx.f"
	    }
#line 483 "cpbsvx.f"
	} else {
#line 484 "cpbsvx.f"
	    i__1 = *n;
#line 484 "cpbsvx.f"
	    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 485 "cpbsvx.f"
		i__2 = j + *kd;
#line 485 "cpbsvx.f"
		j2 = min(i__2,*n);
#line 486 "cpbsvx.f"
		i__2 = j2 - j + 1;
#line 486 "cpbsvx.f"
		ccopy_(&i__2, &ab[j * ab_dim1 + 1], &c__1, &afb[j * afb_dim1 
			+ 1], &c__1);
#line 487 "cpbsvx.f"
/* L50: */
#line 487 "cpbsvx.f"
	    }
#line 488 "cpbsvx.f"
	}

#line 490 "cpbsvx.f"
	cpbtrf_(uplo, n, kd, &afb[afb_offset], ldafb, info, (ftnlen)1);

/*        Return if INFO is non-zero. */

#line 494 "cpbsvx.f"
	if (*info > 0) {
#line 495 "cpbsvx.f"
	    *rcond = 0.;
#line 496 "cpbsvx.f"
	    return 0;
#line 497 "cpbsvx.f"
	}
#line 498 "cpbsvx.f"
    }

/*     Compute the norm of the matrix A. */

#line 502 "cpbsvx.f"
    anorm = clanhb_("1", uplo, n, kd, &ab[ab_offset], ldab, &rwork[1], (
	    ftnlen)1, (ftnlen)1);

/*     Compute the reciprocal of the condition number of A. */

#line 506 "cpbsvx.f"
    cpbcon_(uplo, n, kd, &afb[afb_offset], ldafb, &anorm, rcond, &work[1], &
	    rwork[1], info, (ftnlen)1);

/*     Compute the solution matrix X. */

#line 511 "cpbsvx.f"
    clacpy_("Full", n, nrhs, &b[b_offset], ldb, &x[x_offset], ldx, (ftnlen)4);
#line 512 "cpbsvx.f"
    cpbtrs_(uplo, n, kd, nrhs, &afb[afb_offset], ldafb, &x[x_offset], ldx, 
	    info, (ftnlen)1);

/*     Use iterative refinement to improve the computed solution and */
/*     compute error bounds and backward error estimates for it. */

#line 517 "cpbsvx.f"
    cpbrfs_(uplo, n, kd, nrhs, &ab[ab_offset], ldab, &afb[afb_offset], ldafb, 
	    &b[b_offset], ldb, &x[x_offset], ldx, &ferr[1], &berr[1], &work[1]
	    , &rwork[1], info, (ftnlen)1);

/*     Transform the solution matrix X to a solution of the original */
/*     system. */

#line 523 "cpbsvx.f"
    if (rcequ) {
#line 524 "cpbsvx.f"
	i__1 = *nrhs;
#line 524 "cpbsvx.f"
	for (j = 1; j <= i__1; ++j) {
#line 525 "cpbsvx.f"
	    i__2 = *n;
#line 525 "cpbsvx.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 526 "cpbsvx.f"
		i__3 = i__ + j * x_dim1;
#line 526 "cpbsvx.f"
		i__4 = i__;
#line 526 "cpbsvx.f"
		i__5 = i__ + j * x_dim1;
#line 526 "cpbsvx.f"
		z__1.r = s[i__4] * x[i__5].r, z__1.i = s[i__4] * x[i__5].i;
#line 526 "cpbsvx.f"
		x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 527 "cpbsvx.f"
/* L60: */
#line 527 "cpbsvx.f"
	    }
#line 528 "cpbsvx.f"
/* L70: */
#line 528 "cpbsvx.f"
	}
#line 529 "cpbsvx.f"
	i__1 = *nrhs;
#line 529 "cpbsvx.f"
	for (j = 1; j <= i__1; ++j) {
#line 530 "cpbsvx.f"
	    ferr[j] /= scond;
#line 531 "cpbsvx.f"
/* L80: */
#line 531 "cpbsvx.f"
	}
#line 532 "cpbsvx.f"
    }

/*     Set INFO = N+1 if the matrix is singular to working precision. */

#line 536 "cpbsvx.f"
    if (*rcond < slamch_("Epsilon", (ftnlen)7)) {
#line 536 "cpbsvx.f"
	*info = *n + 1;
#line 536 "cpbsvx.f"
    }

#line 539 "cpbsvx.f"
    return 0;

/*     End of CPBSVX */

} /* cpbsvx_ */

