#line 1 "dgesvx.f"
/* dgesvx.f -- translated by f2c (version 20100827).
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

#line 1 "dgesvx.f"
/* > \brief <b> DGESVX computes the solution to system of linear equations A * X = B for GE matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DGESVX + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgesvx.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgesvx.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgesvx.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DGESVX( FACT, TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV, */
/*                          EQUED, R, C, B, LDB, X, LDX, RCOND, FERR, BERR, */
/*                          WORK, IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          EQUED, FACT, TRANS */
/*       INTEGER            INFO, LDA, LDAF, LDB, LDX, N, NRHS */
/*       DOUBLE PRECISION   RCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ), IWORK( * ) */
/*       DOUBLE PRECISION   A( LDA, * ), AF( LDAF, * ), B( LDB, * ), */
/*      $                   BERR( * ), C( * ), FERR( * ), R( * ), */
/*      $                   WORK( * ), X( LDX, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DGESVX uses the LU factorization to compute the solution to a real */
/* > system of linear equations */
/* >    A * X = B, */
/* > where A is an N-by-N matrix and X and B are N-by-NRHS matrices. */
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
/* >       TRANS = 'N':  diag(R)*A*diag(C)     *inv(diag(C))*X = diag(R)*B */
/* >       TRANS = 'T': (diag(R)*A*diag(C))**T *inv(diag(R))*X = diag(C)*B */
/* >       TRANS = 'C': (diag(R)*A*diag(C))**H *inv(diag(R))*X = diag(C)*B */
/* >    Whether or not the system will be equilibrated depends on the */
/* >    scaling of the matrix A, but if equilibration is used, A is */
/* >    overwritten by diag(R)*A*diag(C) and B by diag(R)*B (if TRANS='N') */
/* >    or diag(C)*B (if TRANS = 'T' or 'C'). */
/* > */
/* > 2. If FACT = 'N' or 'E', the LU decomposition is used to factor the */
/* >    matrix A (after equilibration if FACT = 'E') as */
/* >       A = P * L * U, */
/* >    where P is a permutation matrix, L is a unit lower triangular */
/* >    matrix, and U is upper triangular. */
/* > */
/* > 3. If some U(i,i)=0, so that U is exactly singular, then the routine */
/* >    returns with INFO = i. Otherwise, the factored form of A is used */
/* >    to estimate the condition number of the matrix A.  If the */
/* >    reciprocal of the condition number is less than machine precision, */
/* >    INFO = N+1 is returned as a warning, but the routine still goes on */
/* >    to solve for X and compute error bounds as described below. */
/* > */
/* > 4. The system of equations is solved for X using the factored form */
/* >    of A. */
/* > */
/* > 5. Iterative refinement is applied to improve the computed solution */
/* >    matrix and calculate error bounds and backward error estimates */
/* >    for it. */
/* > */
/* > 6. If equilibration was used, the matrix X is premultiplied by */
/* >    diag(C) (if TRANS = 'N') or diag(R) (if TRANS = 'T' or 'C') so */
/* >    that it solves the original system before equilibration. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] FACT */
/* > \verbatim */
/* >          FACT is CHARACTER*1 */
/* >          Specifies whether or not the factored form of the matrix A is */
/* >          supplied on entry, and if not, whether the matrix A should be */
/* >          equilibrated before it is factored. */
/* >          = 'F':  On entry, AF and IPIV contain the factored form of A. */
/* >                  If EQUED is not 'N', the matrix A has been */
/* >                  equilibrated with scaling factors given by R and C. */
/* >                  A, AF, and IPIV are not modified. */
/* >          = 'N':  The matrix A will be copied to AF and factored. */
/* >          = 'E':  The matrix A will be equilibrated if necessary, then */
/* >                  copied to AF and factored. */
/* > \endverbatim */
/* > */
/* > \param[in] TRANS */
/* > \verbatim */
/* >          TRANS is CHARACTER*1 */
/* >          Specifies the form of the system of equations: */
/* >          = 'N':  A * X = B     (No transpose) */
/* >          = 'T':  A**T * X = B  (Transpose) */
/* >          = 'C':  A**H * X = B  (Transpose) */
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
/* >          On entry, the N-by-N matrix A.  If FACT = 'F' and EQUED is */
/* >          not 'N', then A must have been equilibrated by the scaling */
/* >          factors in R and/or C.  A is not modified if FACT = 'F' or */
/* >          'N', or if FACT = 'E' and EQUED = 'N' on exit. */
/* > */
/* >          On exit, if EQUED .ne. 'N', A is scaled as follows: */
/* >          EQUED = 'R':  A := diag(R) * A */
/* >          EQUED = 'C':  A := A * diag(C) */
/* >          EQUED = 'B':  A := diag(R) * A * diag(C). */
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
/* >          contains the factors L and U from the factorization */
/* >          A = P*L*U as computed by DGETRF.  If EQUED .ne. 'N', then */
/* >          AF is the factored form of the equilibrated matrix A. */
/* > */
/* >          If FACT = 'N', then AF is an output argument and on exit */
/* >          returns the factors L and U from the factorization A = P*L*U */
/* >          of the original matrix A. */
/* > */
/* >          If FACT = 'E', then AF is an output argument and on exit */
/* >          returns the factors L and U from the factorization A = P*L*U */
/* >          of the equilibrated matrix A (see the description of A for */
/* >          the form of the equilibrated matrix). */
/* > \endverbatim */
/* > */
/* > \param[in] LDAF */
/* > \verbatim */
/* >          LDAF is INTEGER */
/* >          The leading dimension of the array AF.  LDAF >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N) */
/* >          If FACT = 'F', then IPIV is an input argument and on entry */
/* >          contains the pivot indices from the factorization A = P*L*U */
/* >          as computed by DGETRF; row i of the matrix was interchanged */
/* >          with row IPIV(i). */
/* > */
/* >          If FACT = 'N', then IPIV is an output argument and on exit */
/* >          contains the pivot indices from the factorization A = P*L*U */
/* >          of the original matrix A. */
/* > */
/* >          If FACT = 'E', then IPIV is an output argument and on exit */
/* >          contains the pivot indices from the factorization A = P*L*U */
/* >          of the equilibrated matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in,out] EQUED */
/* > \verbatim */
/* >          EQUED is CHARACTER*1 */
/* >          Specifies the form of equilibration that was done. */
/* >          = 'N':  No equilibration (always true if FACT = 'N'). */
/* >          = 'R':  Row equilibration, i.e., A has been premultiplied by */
/* >                  diag(R). */
/* >          = 'C':  Column equilibration, i.e., A has been postmultiplied */
/* >                  by diag(C). */
/* >          = 'B':  Both row and column equilibration, i.e., A has been */
/* >                  replaced by diag(R) * A * diag(C). */
/* >          EQUED is an input argument if FACT = 'F'; otherwise, it is an */
/* >          output argument. */
/* > \endverbatim */
/* > */
/* > \param[in,out] R */
/* > \verbatim */
/* >          R is DOUBLE PRECISION array, dimension (N) */
/* >          The row scale factors for A.  If EQUED = 'R' or 'B', A is */
/* >          multiplied on the left by diag(R); if EQUED = 'N' or 'C', R */
/* >          is not accessed.  R is an input argument if FACT = 'F'; */
/* >          otherwise, R is an output argument.  If FACT = 'F' and */
/* >          EQUED = 'R' or 'B', each element of R must be positive. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* >          C is DOUBLE PRECISION array, dimension (N) */
/* >          The column scale factors for A.  If EQUED = 'C' or 'B', A is */
/* >          multiplied on the right by diag(C); if EQUED = 'N' or 'R', C */
/* >          is not accessed.  C is an input argument if FACT = 'F'; */
/* >          otherwise, C is an output argument.  If FACT = 'F' and */
/* >          EQUED = 'C' or 'B', each element of C must be positive. */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* >          B is DOUBLE PRECISION array, dimension (LDB,NRHS) */
/* >          On entry, the N-by-NRHS right hand side matrix B. */
/* >          On exit, */
/* >          if EQUED = 'N', B is not modified; */
/* >          if TRANS = 'N' and EQUED = 'R' or 'B', B is overwritten by */
/* >          diag(R)*B; */
/* >          if TRANS = 'T' or 'C' and EQUED = 'C' or 'B', B is */
/* >          overwritten by diag(C)*B. */
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
/* >          If INFO = 0 or INFO = N+1, the N-by-NRHS solution matrix X */
/* >          to the original system of equations.  Note that A and B are */
/* >          modified on exit if EQUED .ne. 'N', and the solution to the */
/* >          equilibrated system is inv(diag(C))*X if TRANS = 'N' and */
/* >          EQUED = 'C' or 'B', or inv(diag(R))*X if TRANS = 'T' or 'C' */
/* >          and EQUED = 'R' or 'B'. */
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
/* >          WORK is DOUBLE PRECISION array, dimension (4*N) */
/* >          On exit, WORK(1) contains the reciprocal pivot growth */
/* >          factor norm(A)/norm(U). The "max absolute element" norm is */
/* >          used. If WORK(1) is much less than 1, then the stability */
/* >          of the LU factorization of the (equilibrated) matrix A */
/* >          could be poor. This also means that the solution X, condition */
/* >          estimator RCOND, and forward error bound FERR could be */
/* >          unreliable. If factorization fails with 0<INFO<=N, then */
/* >          WORK(1) contains the reciprocal pivot growth factor for the */
/* >          leading INFO columns of A. */
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
/* >                <= N:  U(i,i) is exactly zero.  The factorization has */
/* >                       been completed, but the factor U is exactly */
/* >                       singular, so the solution and error bounds */
/* >                       could not be computed. RCOND = 0 is returned. */
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

/* > \ingroup doubleGEsolve */

/*  ===================================================================== */
/* Subroutine */ int dgesvx_(char *fact, char *trans, integer *n, integer *
	nrhs, doublereal *a, integer *lda, doublereal *af, integer *ldaf, 
	integer *ipiv, char *equed, doublereal *r__, doublereal *c__, 
	doublereal *b, integer *ldb, doublereal *x, integer *ldx, doublereal *
	rcond, doublereal *ferr, doublereal *berr, doublereal *work, integer *
	iwork, integer *info, ftnlen fact_len, ftnlen trans_len, ftnlen 
	equed_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, af_dim1, af_offset, b_dim1, b_offset, x_dim1, 
	    x_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, j;
    static doublereal amax;
    static char norm[1];
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal rcmin, rcmax, anorm;
    static logical equil;
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int dlaqge_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, char *, ftnlen), dgecon_(char *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     integer *, integer *, ftnlen);
    static doublereal colcnd;
    static logical nofact;
    extern /* Subroutine */ int dgeequ_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, integer *), dgerfs_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, ftnlen), 
	    dgetrf_(integer *, integer *, doublereal *, integer *, integer *, 
	    integer *), dlacpy_(char *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, ftnlen), xerbla_(char *, 
	    integer *, ftnlen);
    static doublereal bignum;
    extern doublereal dlantr_(char *, char *, char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, ftnlen, ftnlen, ftnlen);
    static integer infequ;
    static logical colequ;
    extern /* Subroutine */ int dgetrs_(char *, integer *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, ftnlen);
    static doublereal rowcnd;
    static logical notran;
    static doublereal smlnum;
    static logical rowequ;
    static doublereal rpvgrw;


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

#line 396 "dgesvx.f"
    /* Parameter adjustments */
#line 396 "dgesvx.f"
    a_dim1 = *lda;
#line 396 "dgesvx.f"
    a_offset = 1 + a_dim1;
#line 396 "dgesvx.f"
    a -= a_offset;
#line 396 "dgesvx.f"
    af_dim1 = *ldaf;
#line 396 "dgesvx.f"
    af_offset = 1 + af_dim1;
#line 396 "dgesvx.f"
    af -= af_offset;
#line 396 "dgesvx.f"
    --ipiv;
#line 396 "dgesvx.f"
    --r__;
#line 396 "dgesvx.f"
    --c__;
#line 396 "dgesvx.f"
    b_dim1 = *ldb;
#line 396 "dgesvx.f"
    b_offset = 1 + b_dim1;
#line 396 "dgesvx.f"
    b -= b_offset;
#line 396 "dgesvx.f"
    x_dim1 = *ldx;
#line 396 "dgesvx.f"
    x_offset = 1 + x_dim1;
#line 396 "dgesvx.f"
    x -= x_offset;
#line 396 "dgesvx.f"
    --ferr;
#line 396 "dgesvx.f"
    --berr;
#line 396 "dgesvx.f"
    --work;
#line 396 "dgesvx.f"
    --iwork;
#line 396 "dgesvx.f"

#line 396 "dgesvx.f"
    /* Function Body */
#line 396 "dgesvx.f"
    *info = 0;
#line 397 "dgesvx.f"
    nofact = lsame_(fact, "N", (ftnlen)1, (ftnlen)1);
#line 398 "dgesvx.f"
    equil = lsame_(fact, "E", (ftnlen)1, (ftnlen)1);
#line 399 "dgesvx.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 400 "dgesvx.f"
    if (nofact || equil) {
#line 401 "dgesvx.f"
	*(unsigned char *)equed = 'N';
#line 402 "dgesvx.f"
	rowequ = FALSE_;
#line 403 "dgesvx.f"
	colequ = FALSE_;
#line 404 "dgesvx.f"
    } else {
#line 405 "dgesvx.f"
	rowequ = lsame_(equed, "R", (ftnlen)1, (ftnlen)1) || lsame_(equed, 
		"B", (ftnlen)1, (ftnlen)1);
#line 406 "dgesvx.f"
	colequ = lsame_(equed, "C", (ftnlen)1, (ftnlen)1) || lsame_(equed, 
		"B", (ftnlen)1, (ftnlen)1);
#line 407 "dgesvx.f"
	smlnum = dlamch_("Safe minimum", (ftnlen)12);
#line 408 "dgesvx.f"
	bignum = 1. / smlnum;
#line 409 "dgesvx.f"
    }

/*     Test the input parameters. */

#line 413 "dgesvx.f"
    if (! nofact && ! equil && ! lsame_(fact, "F", (ftnlen)1, (ftnlen)1)) {
#line 415 "dgesvx.f"
	*info = -1;
#line 416 "dgesvx.f"
    } else if (! notran && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1) && ! 
	    lsame_(trans, "C", (ftnlen)1, (ftnlen)1)) {
#line 418 "dgesvx.f"
	*info = -2;
#line 419 "dgesvx.f"
    } else if (*n < 0) {
#line 420 "dgesvx.f"
	*info = -3;
#line 421 "dgesvx.f"
    } else if (*nrhs < 0) {
#line 422 "dgesvx.f"
	*info = -4;
#line 423 "dgesvx.f"
    } else if (*lda < max(1,*n)) {
#line 424 "dgesvx.f"
	*info = -6;
#line 425 "dgesvx.f"
    } else if (*ldaf < max(1,*n)) {
#line 426 "dgesvx.f"
	*info = -8;
#line 427 "dgesvx.f"
    } else if (lsame_(fact, "F", (ftnlen)1, (ftnlen)1) && ! (rowequ || colequ 
	    || lsame_(equed, "N", (ftnlen)1, (ftnlen)1))) {
#line 429 "dgesvx.f"
	*info = -10;
#line 430 "dgesvx.f"
    } else {
#line 431 "dgesvx.f"
	if (rowequ) {
#line 432 "dgesvx.f"
	    rcmin = bignum;
#line 433 "dgesvx.f"
	    rcmax = 0.;
#line 434 "dgesvx.f"
	    i__1 = *n;
#line 434 "dgesvx.f"
	    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 435 "dgesvx.f"
		d__1 = rcmin, d__2 = r__[j];
#line 435 "dgesvx.f"
		rcmin = min(d__1,d__2);
/* Computing MAX */
#line 436 "dgesvx.f"
		d__1 = rcmax, d__2 = r__[j];
#line 436 "dgesvx.f"
		rcmax = max(d__1,d__2);
#line 437 "dgesvx.f"
/* L10: */
#line 437 "dgesvx.f"
	    }
#line 438 "dgesvx.f"
	    if (rcmin <= 0.) {
#line 439 "dgesvx.f"
		*info = -11;
#line 440 "dgesvx.f"
	    } else if (*n > 0) {
#line 441 "dgesvx.f"
		rowcnd = max(rcmin,smlnum) / min(rcmax,bignum);
#line 442 "dgesvx.f"
	    } else {
#line 443 "dgesvx.f"
		rowcnd = 1.;
#line 444 "dgesvx.f"
	    }
#line 445 "dgesvx.f"
	}
#line 446 "dgesvx.f"
	if (colequ && *info == 0) {
#line 447 "dgesvx.f"
	    rcmin = bignum;
#line 448 "dgesvx.f"
	    rcmax = 0.;
#line 449 "dgesvx.f"
	    i__1 = *n;
#line 449 "dgesvx.f"
	    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 450 "dgesvx.f"
		d__1 = rcmin, d__2 = c__[j];
#line 450 "dgesvx.f"
		rcmin = min(d__1,d__2);
/* Computing MAX */
#line 451 "dgesvx.f"
		d__1 = rcmax, d__2 = c__[j];
#line 451 "dgesvx.f"
		rcmax = max(d__1,d__2);
#line 452 "dgesvx.f"
/* L20: */
#line 452 "dgesvx.f"
	    }
#line 453 "dgesvx.f"
	    if (rcmin <= 0.) {
#line 454 "dgesvx.f"
		*info = -12;
#line 455 "dgesvx.f"
	    } else if (*n > 0) {
#line 456 "dgesvx.f"
		colcnd = max(rcmin,smlnum) / min(rcmax,bignum);
#line 457 "dgesvx.f"
	    } else {
#line 458 "dgesvx.f"
		colcnd = 1.;
#line 459 "dgesvx.f"
	    }
#line 460 "dgesvx.f"
	}
#line 461 "dgesvx.f"
	if (*info == 0) {
#line 462 "dgesvx.f"
	    if (*ldb < max(1,*n)) {
#line 463 "dgesvx.f"
		*info = -14;
#line 464 "dgesvx.f"
	    } else if (*ldx < max(1,*n)) {
#line 465 "dgesvx.f"
		*info = -16;
#line 466 "dgesvx.f"
	    }
#line 467 "dgesvx.f"
	}
#line 468 "dgesvx.f"
    }

#line 470 "dgesvx.f"
    if (*info != 0) {
#line 471 "dgesvx.f"
	i__1 = -(*info);
#line 471 "dgesvx.f"
	xerbla_("DGESVX", &i__1, (ftnlen)6);
#line 472 "dgesvx.f"
	return 0;
#line 473 "dgesvx.f"
    }

#line 475 "dgesvx.f"
    if (equil) {

/*        Compute row and column scalings to equilibrate the matrix A. */

#line 479 "dgesvx.f"
	dgeequ_(n, n, &a[a_offset], lda, &r__[1], &c__[1], &rowcnd, &colcnd, &
		amax, &infequ);
#line 480 "dgesvx.f"
	if (infequ == 0) {

/*           Equilibrate the matrix. */

#line 484 "dgesvx.f"
	    dlaqge_(n, n, &a[a_offset], lda, &r__[1], &c__[1], &rowcnd, &
		    colcnd, &amax, equed, (ftnlen)1);
#line 486 "dgesvx.f"
	    rowequ = lsame_(equed, "R", (ftnlen)1, (ftnlen)1) || lsame_(equed,
		     "B", (ftnlen)1, (ftnlen)1);
#line 487 "dgesvx.f"
	    colequ = lsame_(equed, "C", (ftnlen)1, (ftnlen)1) || lsame_(equed,
		     "B", (ftnlen)1, (ftnlen)1);
#line 488 "dgesvx.f"
	}
#line 489 "dgesvx.f"
    }

/*     Scale the right hand side. */

#line 493 "dgesvx.f"
    if (notran) {
#line 494 "dgesvx.f"
	if (rowequ) {
#line 495 "dgesvx.f"
	    i__1 = *nrhs;
#line 495 "dgesvx.f"
	    for (j = 1; j <= i__1; ++j) {
#line 496 "dgesvx.f"
		i__2 = *n;
#line 496 "dgesvx.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 497 "dgesvx.f"
		    b[i__ + j * b_dim1] = r__[i__] * b[i__ + j * b_dim1];
#line 498 "dgesvx.f"
/* L30: */
#line 498 "dgesvx.f"
		}
#line 499 "dgesvx.f"
/* L40: */
#line 499 "dgesvx.f"
	    }
#line 500 "dgesvx.f"
	}
#line 501 "dgesvx.f"
    } else if (colequ) {
#line 502 "dgesvx.f"
	i__1 = *nrhs;
#line 502 "dgesvx.f"
	for (j = 1; j <= i__1; ++j) {
#line 503 "dgesvx.f"
	    i__2 = *n;
#line 503 "dgesvx.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 504 "dgesvx.f"
		b[i__ + j * b_dim1] = c__[i__] * b[i__ + j * b_dim1];
#line 505 "dgesvx.f"
/* L50: */
#line 505 "dgesvx.f"
	    }
#line 506 "dgesvx.f"
/* L60: */
#line 506 "dgesvx.f"
	}
#line 507 "dgesvx.f"
    }

#line 509 "dgesvx.f"
    if (nofact || equil) {

/*        Compute the LU factorization of A. */

#line 513 "dgesvx.f"
	dlacpy_("Full", n, n, &a[a_offset], lda, &af[af_offset], ldaf, (
		ftnlen)4);
#line 514 "dgesvx.f"
	dgetrf_(n, n, &af[af_offset], ldaf, &ipiv[1], info);

/*        Return if INFO is non-zero. */

#line 518 "dgesvx.f"
	if (*info > 0) {

/*           Compute the reciprocal pivot growth factor of the */
/*           leading rank-deficient INFO columns of A. */

#line 523 "dgesvx.f"
	    rpvgrw = dlantr_("M", "U", "N", info, info, &af[af_offset], ldaf, 
		    &work[1], (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 525 "dgesvx.f"
	    if (rpvgrw == 0.) {
#line 526 "dgesvx.f"
		rpvgrw = 1.;
#line 527 "dgesvx.f"
	    } else {
#line 528 "dgesvx.f"
		rpvgrw = dlange_("M", n, info, &a[a_offset], lda, &work[1], (
			ftnlen)1) / rpvgrw;
#line 529 "dgesvx.f"
	    }
#line 530 "dgesvx.f"
	    work[1] = rpvgrw;
#line 531 "dgesvx.f"
	    *rcond = 0.;
#line 532 "dgesvx.f"
	    return 0;
#line 533 "dgesvx.f"
	}
#line 534 "dgesvx.f"
    }

/*     Compute the norm of the matrix A and the */
/*     reciprocal pivot growth factor RPVGRW. */

#line 539 "dgesvx.f"
    if (notran) {
#line 540 "dgesvx.f"
	*(unsigned char *)norm = '1';
#line 541 "dgesvx.f"
    } else {
#line 542 "dgesvx.f"
	*(unsigned char *)norm = 'I';
#line 543 "dgesvx.f"
    }
#line 544 "dgesvx.f"
    anorm = dlange_(norm, n, n, &a[a_offset], lda, &work[1], (ftnlen)1);
#line 545 "dgesvx.f"
    rpvgrw = dlantr_("M", "U", "N", n, n, &af[af_offset], ldaf, &work[1], (
	    ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 546 "dgesvx.f"
    if (rpvgrw == 0.) {
#line 547 "dgesvx.f"
	rpvgrw = 1.;
#line 548 "dgesvx.f"
    } else {
#line 549 "dgesvx.f"
	rpvgrw = dlange_("M", n, n, &a[a_offset], lda, &work[1], (ftnlen)1) / 
		rpvgrw;
#line 550 "dgesvx.f"
    }

/*     Compute the reciprocal of the condition number of A. */

#line 554 "dgesvx.f"
    dgecon_(norm, n, &af[af_offset], ldaf, &anorm, rcond, &work[1], &iwork[1],
	     info, (ftnlen)1);

/*     Compute the solution matrix X. */

#line 558 "dgesvx.f"
    dlacpy_("Full", n, nrhs, &b[b_offset], ldb, &x[x_offset], ldx, (ftnlen)4);
#line 559 "dgesvx.f"
    dgetrs_(trans, n, nrhs, &af[af_offset], ldaf, &ipiv[1], &x[x_offset], ldx,
	     info, (ftnlen)1);

/*     Use iterative refinement to improve the computed solution and */
/*     compute error bounds and backward error estimates for it. */

#line 564 "dgesvx.f"
    dgerfs_(trans, n, nrhs, &a[a_offset], lda, &af[af_offset], ldaf, &ipiv[1],
	     &b[b_offset], ldb, &x[x_offset], ldx, &ferr[1], &berr[1], &work[
	    1], &iwork[1], info, (ftnlen)1);

/*     Transform the solution matrix X to a solution of the original */
/*     system. */

#line 570 "dgesvx.f"
    if (notran) {
#line 571 "dgesvx.f"
	if (colequ) {
#line 572 "dgesvx.f"
	    i__1 = *nrhs;
#line 572 "dgesvx.f"
	    for (j = 1; j <= i__1; ++j) {
#line 573 "dgesvx.f"
		i__2 = *n;
#line 573 "dgesvx.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 574 "dgesvx.f"
		    x[i__ + j * x_dim1] = c__[i__] * x[i__ + j * x_dim1];
#line 575 "dgesvx.f"
/* L70: */
#line 575 "dgesvx.f"
		}
#line 576 "dgesvx.f"
/* L80: */
#line 576 "dgesvx.f"
	    }
#line 577 "dgesvx.f"
	    i__1 = *nrhs;
#line 577 "dgesvx.f"
	    for (j = 1; j <= i__1; ++j) {
#line 578 "dgesvx.f"
		ferr[j] /= colcnd;
#line 579 "dgesvx.f"
/* L90: */
#line 579 "dgesvx.f"
	    }
#line 580 "dgesvx.f"
	}
#line 581 "dgesvx.f"
    } else if (rowequ) {
#line 582 "dgesvx.f"
	i__1 = *nrhs;
#line 582 "dgesvx.f"
	for (j = 1; j <= i__1; ++j) {
#line 583 "dgesvx.f"
	    i__2 = *n;
#line 583 "dgesvx.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 584 "dgesvx.f"
		x[i__ + j * x_dim1] = r__[i__] * x[i__ + j * x_dim1];
#line 585 "dgesvx.f"
/* L100: */
#line 585 "dgesvx.f"
	    }
#line 586 "dgesvx.f"
/* L110: */
#line 586 "dgesvx.f"
	}
#line 587 "dgesvx.f"
	i__1 = *nrhs;
#line 587 "dgesvx.f"
	for (j = 1; j <= i__1; ++j) {
#line 588 "dgesvx.f"
	    ferr[j] /= rowcnd;
#line 589 "dgesvx.f"
/* L120: */
#line 589 "dgesvx.f"
	}
#line 590 "dgesvx.f"
    }

#line 592 "dgesvx.f"
    work[1] = rpvgrw;

/*     Set INFO = N+1 if the matrix is singular to working precision. */

#line 596 "dgesvx.f"
    if (*rcond < dlamch_("Epsilon", (ftnlen)7)) {
#line 596 "dgesvx.f"
	*info = *n + 1;
#line 596 "dgesvx.f"
    }
#line 598 "dgesvx.f"
    return 0;

/*     End of DGESVX */

} /* dgesvx_ */

