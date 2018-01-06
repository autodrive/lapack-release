#line 1 "dpbsvx.f"
/* dpbsvx.f -- translated by f2c (version 20100827).
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

#line 1 "dpbsvx.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief <b> DPBSVX computes the solution to system of linear equations A * X = B for OTHER matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DPBSVX + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dpbsvx.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dpbsvx.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dpbsvx.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DPBSVX( FACT, UPLO, N, KD, NRHS, AB, LDAB, AFB, LDAFB, */
/*                          EQUED, S, B, LDB, X, LDX, RCOND, FERR, BERR, */
/*                          WORK, IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          EQUED, FACT, UPLO */
/*       INTEGER            INFO, KD, LDAB, LDAFB, LDB, LDX, N, NRHS */
/*       DOUBLE PRECISION   RCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       DOUBLE PRECISION   AB( LDAB, * ), AFB( LDAFB, * ), B( LDB, * ), */
/*      $                   BERR( * ), FERR( * ), S( * ), WORK( * ), */
/*      $                   X( LDX, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DPBSVX uses the Cholesky factorization A = U**T*U or A = L*L**T to */
/* > compute the solution to a real system of linear equations */
/* >    A * X = B, */
/* > where A is an N-by-N symmetric positive definite band matrix and X */
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
/* >       A = U**T * U,  if UPLO = 'U', or */
/* >       A = L * L**T,  if UPLO = 'L', */
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
/* >          AB is DOUBLE PRECISION array, dimension (LDAB,N) */
/* >          On entry, the upper or lower triangle of the symmetric band */
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
/* >          AFB is DOUBLE PRECISION array, dimension (LDAFB,N) */
/* >          If FACT = 'F', then AFB is an input argument and on entry */
/* >          contains the triangular factor U or L from the Cholesky */
/* >          factorization A = U**T*U or A = L*L**T of the band matrix */
/* >          A, in the same storage format as A (see AB).  If EQUED = 'Y', */
/* >          then AFB is the factored form of the equilibrated matrix A. */
/* > */
/* >          If FACT = 'N', then AFB is an output argument and on exit */
/* >          returns the triangular factor U or L from the Cholesky */
/* >          factorization A = U**T*U or A = L*L**T. */
/* > */
/* >          If FACT = 'E', then AFB is an output argument and on exit */
/* >          returns the triangular factor U or L from the Cholesky */
/* >          factorization A = U**T*U or A = L*L**T of the equilibrated */
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

/* > \ingroup doubleOTHERsolve */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  The band storage scheme is illustrated by the following example, when */
/* >  N = 6, KD = 2, and UPLO = 'U': */
/* > */
/* >  Two-dimensional storage of the symmetric matrix A: */
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
/* Subroutine */ int dpbsvx_(char *fact, char *uplo, integer *n, integer *kd, 
	integer *nrhs, doublereal *ab, integer *ldab, doublereal *afb, 
	integer *ldafb, char *equed, doublereal *s, doublereal *b, integer *
	ldb, doublereal *x, integer *ldx, doublereal *rcond, doublereal *ferr,
	 doublereal *berr, doublereal *work, integer *iwork, integer *info, 
	ftnlen fact_len, ftnlen uplo_len, ftnlen equed_len)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, afb_dim1, afb_offset, b_dim1, b_offset, 
	    x_dim1, x_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, j, j1, j2;
    static doublereal amax, smin, smax;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal scond, anorm;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static logical equil, rcequ, upper;
    extern doublereal dlamch_(char *, ftnlen), dlansb_(char *, char *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    ftnlen, ftnlen);
    extern /* Subroutine */ int dpbcon_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     integer *, integer *, ftnlen), dlaqsb_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, char *, ftnlen, ftnlen);
    static logical nofact;
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen), dpbequ_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen);
    static doublereal bignum;
    extern /* Subroutine */ int dpbrfs_(char *, integer *, integer *, integer 
	    *, doublereal *, integer *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, ftnlen), dpbtrf_(char *, 
	    integer *, integer *, doublereal *, integer *, integer *, ftnlen);
    static integer infequ;
    extern /* Subroutine */ int dpbtrs_(char *, integer *, integer *, integer 
	    *, doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen);
    static doublereal smlnum;


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

#line 388 "dpbsvx.f"
    /* Parameter adjustments */
#line 388 "dpbsvx.f"
    ab_dim1 = *ldab;
#line 388 "dpbsvx.f"
    ab_offset = 1 + ab_dim1;
#line 388 "dpbsvx.f"
    ab -= ab_offset;
#line 388 "dpbsvx.f"
    afb_dim1 = *ldafb;
#line 388 "dpbsvx.f"
    afb_offset = 1 + afb_dim1;
#line 388 "dpbsvx.f"
    afb -= afb_offset;
#line 388 "dpbsvx.f"
    --s;
#line 388 "dpbsvx.f"
    b_dim1 = *ldb;
#line 388 "dpbsvx.f"
    b_offset = 1 + b_dim1;
#line 388 "dpbsvx.f"
    b -= b_offset;
#line 388 "dpbsvx.f"
    x_dim1 = *ldx;
#line 388 "dpbsvx.f"
    x_offset = 1 + x_dim1;
#line 388 "dpbsvx.f"
    x -= x_offset;
#line 388 "dpbsvx.f"
    --ferr;
#line 388 "dpbsvx.f"
    --berr;
#line 388 "dpbsvx.f"
    --work;
#line 388 "dpbsvx.f"
    --iwork;
#line 388 "dpbsvx.f"

#line 388 "dpbsvx.f"
    /* Function Body */
#line 388 "dpbsvx.f"
    *info = 0;
#line 389 "dpbsvx.f"
    nofact = lsame_(fact, "N", (ftnlen)1, (ftnlen)1);
#line 390 "dpbsvx.f"
    equil = lsame_(fact, "E", (ftnlen)1, (ftnlen)1);
#line 391 "dpbsvx.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 392 "dpbsvx.f"
    if (nofact || equil) {
#line 393 "dpbsvx.f"
	*(unsigned char *)equed = 'N';
#line 394 "dpbsvx.f"
	rcequ = FALSE_;
#line 395 "dpbsvx.f"
    } else {
#line 396 "dpbsvx.f"
	rcequ = lsame_(equed, "Y", (ftnlen)1, (ftnlen)1);
#line 397 "dpbsvx.f"
	smlnum = dlamch_("Safe minimum", (ftnlen)12);
#line 398 "dpbsvx.f"
	bignum = 1. / smlnum;
#line 399 "dpbsvx.f"
    }

/*     Test the input parameters. */

#line 403 "dpbsvx.f"
    if (! nofact && ! equil && ! lsame_(fact, "F", (ftnlen)1, (ftnlen)1)) {
#line 405 "dpbsvx.f"
	*info = -1;
#line 406 "dpbsvx.f"
    } else if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 407 "dpbsvx.f"
	*info = -2;
#line 408 "dpbsvx.f"
    } else if (*n < 0) {
#line 409 "dpbsvx.f"
	*info = -3;
#line 410 "dpbsvx.f"
    } else if (*kd < 0) {
#line 411 "dpbsvx.f"
	*info = -4;
#line 412 "dpbsvx.f"
    } else if (*nrhs < 0) {
#line 413 "dpbsvx.f"
	*info = -5;
#line 414 "dpbsvx.f"
    } else if (*ldab < *kd + 1) {
#line 415 "dpbsvx.f"
	*info = -7;
#line 416 "dpbsvx.f"
    } else if (*ldafb < *kd + 1) {
#line 417 "dpbsvx.f"
	*info = -9;
#line 418 "dpbsvx.f"
    } else if (lsame_(fact, "F", (ftnlen)1, (ftnlen)1) && ! (rcequ || lsame_(
	    equed, "N", (ftnlen)1, (ftnlen)1))) {
#line 420 "dpbsvx.f"
	*info = -10;
#line 421 "dpbsvx.f"
    } else {
#line 422 "dpbsvx.f"
	if (rcequ) {
#line 423 "dpbsvx.f"
	    smin = bignum;
#line 424 "dpbsvx.f"
	    smax = 0.;
#line 425 "dpbsvx.f"
	    i__1 = *n;
#line 425 "dpbsvx.f"
	    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 426 "dpbsvx.f"
		d__1 = smin, d__2 = s[j];
#line 426 "dpbsvx.f"
		smin = min(d__1,d__2);
/* Computing MAX */
#line 427 "dpbsvx.f"
		d__1 = smax, d__2 = s[j];
#line 427 "dpbsvx.f"
		smax = max(d__1,d__2);
#line 428 "dpbsvx.f"
/* L10: */
#line 428 "dpbsvx.f"
	    }
#line 429 "dpbsvx.f"
	    if (smin <= 0.) {
#line 430 "dpbsvx.f"
		*info = -11;
#line 431 "dpbsvx.f"
	    } else if (*n > 0) {
#line 432 "dpbsvx.f"
		scond = max(smin,smlnum) / min(smax,bignum);
#line 433 "dpbsvx.f"
	    } else {
#line 434 "dpbsvx.f"
		scond = 1.;
#line 435 "dpbsvx.f"
	    }
#line 436 "dpbsvx.f"
	}
#line 437 "dpbsvx.f"
	if (*info == 0) {
#line 438 "dpbsvx.f"
	    if (*ldb < max(1,*n)) {
#line 439 "dpbsvx.f"
		*info = -13;
#line 440 "dpbsvx.f"
	    } else if (*ldx < max(1,*n)) {
#line 441 "dpbsvx.f"
		*info = -15;
#line 442 "dpbsvx.f"
	    }
#line 443 "dpbsvx.f"
	}
#line 444 "dpbsvx.f"
    }

#line 446 "dpbsvx.f"
    if (*info != 0) {
#line 447 "dpbsvx.f"
	i__1 = -(*info);
#line 447 "dpbsvx.f"
	xerbla_("DPBSVX", &i__1, (ftnlen)6);
#line 448 "dpbsvx.f"
	return 0;
#line 449 "dpbsvx.f"
    }

#line 451 "dpbsvx.f"
    if (equil) {

/*        Compute row and column scalings to equilibrate the matrix A. */

#line 455 "dpbsvx.f"
	dpbequ_(uplo, n, kd, &ab[ab_offset], ldab, &s[1], &scond, &amax, &
		infequ, (ftnlen)1);
#line 456 "dpbsvx.f"
	if (infequ == 0) {

/*           Equilibrate the matrix. */

#line 460 "dpbsvx.f"
	    dlaqsb_(uplo, n, kd, &ab[ab_offset], ldab, &s[1], &scond, &amax, 
		    equed, (ftnlen)1, (ftnlen)1);
#line 461 "dpbsvx.f"
	    rcequ = lsame_(equed, "Y", (ftnlen)1, (ftnlen)1);
#line 462 "dpbsvx.f"
	}
#line 463 "dpbsvx.f"
    }

/*     Scale the right-hand side. */

#line 467 "dpbsvx.f"
    if (rcequ) {
#line 468 "dpbsvx.f"
	i__1 = *nrhs;
#line 468 "dpbsvx.f"
	for (j = 1; j <= i__1; ++j) {
#line 469 "dpbsvx.f"
	    i__2 = *n;
#line 469 "dpbsvx.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 470 "dpbsvx.f"
		b[i__ + j * b_dim1] = s[i__] * b[i__ + j * b_dim1];
#line 471 "dpbsvx.f"
/* L20: */
#line 471 "dpbsvx.f"
	    }
#line 472 "dpbsvx.f"
/* L30: */
#line 472 "dpbsvx.f"
	}
#line 473 "dpbsvx.f"
    }

#line 475 "dpbsvx.f"
    if (nofact || equil) {

/*        Compute the Cholesky factorization A = U**T *U or A = L*L**T. */

#line 479 "dpbsvx.f"
	if (upper) {
#line 480 "dpbsvx.f"
	    i__1 = *n;
#line 480 "dpbsvx.f"
	    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
#line 481 "dpbsvx.f"
		i__2 = j - *kd;
#line 481 "dpbsvx.f"
		j1 = max(i__2,1);
#line 482 "dpbsvx.f"
		i__2 = j - j1 + 1;
#line 482 "dpbsvx.f"
		dcopy_(&i__2, &ab[*kd + 1 - j + j1 + j * ab_dim1], &c__1, &
			afb[*kd + 1 - j + j1 + j * afb_dim1], &c__1);
#line 484 "dpbsvx.f"
/* L40: */
#line 484 "dpbsvx.f"
	    }
#line 485 "dpbsvx.f"
	} else {
#line 486 "dpbsvx.f"
	    i__1 = *n;
#line 486 "dpbsvx.f"
	    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 487 "dpbsvx.f"
		i__2 = j + *kd;
#line 487 "dpbsvx.f"
		j2 = min(i__2,*n);
#line 488 "dpbsvx.f"
		i__2 = j2 - j + 1;
#line 488 "dpbsvx.f"
		dcopy_(&i__2, &ab[j * ab_dim1 + 1], &c__1, &afb[j * afb_dim1 
			+ 1], &c__1);
#line 489 "dpbsvx.f"
/* L50: */
#line 489 "dpbsvx.f"
	    }
#line 490 "dpbsvx.f"
	}

#line 492 "dpbsvx.f"
	dpbtrf_(uplo, n, kd, &afb[afb_offset], ldafb, info, (ftnlen)1);

/*        Return if INFO is non-zero. */

#line 496 "dpbsvx.f"
	if (*info > 0) {
#line 497 "dpbsvx.f"
	    *rcond = 0.;
#line 498 "dpbsvx.f"
	    return 0;
#line 499 "dpbsvx.f"
	}
#line 500 "dpbsvx.f"
    }

/*     Compute the norm of the matrix A. */

#line 504 "dpbsvx.f"
    anorm = dlansb_("1", uplo, n, kd, &ab[ab_offset], ldab, &work[1], (ftnlen)
	    1, (ftnlen)1);

/*     Compute the reciprocal of the condition number of A. */

#line 508 "dpbsvx.f"
    dpbcon_(uplo, n, kd, &afb[afb_offset], ldafb, &anorm, rcond, &work[1], &
	    iwork[1], info, (ftnlen)1);

/*     Compute the solution matrix X. */

#line 513 "dpbsvx.f"
    dlacpy_("Full", n, nrhs, &b[b_offset], ldb, &x[x_offset], ldx, (ftnlen)4);
#line 514 "dpbsvx.f"
    dpbtrs_(uplo, n, kd, nrhs, &afb[afb_offset], ldafb, &x[x_offset], ldx, 
	    info, (ftnlen)1);

/*     Use iterative refinement to improve the computed solution and */
/*     compute error bounds and backward error estimates for it. */

#line 519 "dpbsvx.f"
    dpbrfs_(uplo, n, kd, nrhs, &ab[ab_offset], ldab, &afb[afb_offset], ldafb, 
	    &b[b_offset], ldb, &x[x_offset], ldx, &ferr[1], &berr[1], &work[1]
	    , &iwork[1], info, (ftnlen)1);

/*     Transform the solution matrix X to a solution of the original */
/*     system. */

#line 525 "dpbsvx.f"
    if (rcequ) {
#line 526 "dpbsvx.f"
	i__1 = *nrhs;
#line 526 "dpbsvx.f"
	for (j = 1; j <= i__1; ++j) {
#line 527 "dpbsvx.f"
	    i__2 = *n;
#line 527 "dpbsvx.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 528 "dpbsvx.f"
		x[i__ + j * x_dim1] = s[i__] * x[i__ + j * x_dim1];
#line 529 "dpbsvx.f"
/* L60: */
#line 529 "dpbsvx.f"
	    }
#line 530 "dpbsvx.f"
/* L70: */
#line 530 "dpbsvx.f"
	}
#line 531 "dpbsvx.f"
	i__1 = *nrhs;
#line 531 "dpbsvx.f"
	for (j = 1; j <= i__1; ++j) {
#line 532 "dpbsvx.f"
	    ferr[j] /= scond;
#line 533 "dpbsvx.f"
/* L80: */
#line 533 "dpbsvx.f"
	}
#line 534 "dpbsvx.f"
    }

/*     Set INFO = N+1 if the matrix is singular to working precision. */

#line 538 "dpbsvx.f"
    if (*rcond < dlamch_("Epsilon", (ftnlen)7)) {
#line 538 "dpbsvx.f"
	*info = *n + 1;
#line 538 "dpbsvx.f"
    }

#line 541 "dpbsvx.f"
    return 0;

/*     End of DPBSVX */

} /* dpbsvx_ */

