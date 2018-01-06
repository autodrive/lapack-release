#line 1 "sgesvx.f"
/* sgesvx.f -- translated by f2c (version 20100827).
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

#line 1 "sgesvx.f"
/* > \brief <b> SGESVX computes the solution to system of linear equations A * X = B for GE matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SGESVX + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgesvx.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgesvx.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgesvx.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SGESVX( FACT, TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV, */
/*                          EQUED, R, C, B, LDB, X, LDX, RCOND, FERR, BERR, */
/*                          WORK, IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          EQUED, FACT, TRANS */
/*       INTEGER            INFO, LDA, LDAF, LDB, LDX, N, NRHS */
/*       REAL               RCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ), IWORK( * ) */
/*       REAL               A( LDA, * ), AF( LDAF, * ), B( LDB, * ), */
/*      $                   BERR( * ), C( * ), FERR( * ), R( * ), */
/*      $                   WORK( * ), X( LDX, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SGESVX uses the LU factorization to compute the solution to a real */
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
/* >          A is REAL array, dimension (LDA,N) */
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
/* >          AF is REAL array, dimension (LDAF,N) */
/* >          If FACT = 'F', then AF is an input argument and on entry */
/* >          contains the factors L and U from the factorization */
/* >          A = P*L*U as computed by SGETRF.  If EQUED .ne. 'N', then */
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
/* >          as computed by SGETRF; row i of the matrix was interchanged */
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
/* >          R is REAL array, dimension (N) */
/* >          The row scale factors for A.  If EQUED = 'R' or 'B', A is */
/* >          multiplied on the left by diag(R); if EQUED = 'N' or 'C', R */
/* >          is not accessed.  R is an input argument if FACT = 'F'; */
/* >          otherwise, R is an output argument.  If FACT = 'F' and */
/* >          EQUED = 'R' or 'B', each element of R must be positive. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* >          C is REAL array, dimension (N) */
/* >          The column scale factors for A.  If EQUED = 'C' or 'B', A is */
/* >          multiplied on the right by diag(C); if EQUED = 'N' or 'R', C */
/* >          is not accessed.  C is an input argument if FACT = 'F'; */
/* >          otherwise, C is an output argument.  If FACT = 'F' and */
/* >          EQUED = 'C' or 'B', each element of C must be positive. */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* >          B is REAL array, dimension (LDB,NRHS) */
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
/* >          X is REAL array, dimension (LDX,NRHS) */
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
/* >          WORK is REAL array, dimension (4*N) */
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

/* > \ingroup realGEsolve */

/*  ===================================================================== */
/* Subroutine */ int sgesvx_(char *fact, char *trans, integer *n, integer *
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
    static doublereal colcnd;
    extern doublereal slamch_(char *, ftnlen), slange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    static logical nofact;
    extern /* Subroutine */ int slaqge_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, char *, ftnlen), xerbla_(char *, integer *, ftnlen)
	    , sgecon_(char *, integer *, doublereal *, integer *, doublereal *
	    , doublereal *, doublereal *, integer *, integer *, ftnlen);
    static doublereal bignum;
    static integer infequ;
    static logical colequ;
    extern /* Subroutine */ int sgeequ_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, integer *), sgerfs_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, ftnlen), 
	    sgetrf_(integer *, integer *, doublereal *, integer *, integer *, 
	    integer *);
    static doublereal rowcnd;
    extern /* Subroutine */ int slacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen);
    static logical notran;
    extern doublereal slantr_(char *, char *, char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, ftnlen, ftnlen, ftnlen);
    extern /* Subroutine */ int sgetrs_(char *, integer *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, ftnlen);
    static doublereal smlnum;
    static logical rowequ;
    static doublereal rpvgrw;


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

#line 396 "sgesvx.f"
    /* Parameter adjustments */
#line 396 "sgesvx.f"
    a_dim1 = *lda;
#line 396 "sgesvx.f"
    a_offset = 1 + a_dim1;
#line 396 "sgesvx.f"
    a -= a_offset;
#line 396 "sgesvx.f"
    af_dim1 = *ldaf;
#line 396 "sgesvx.f"
    af_offset = 1 + af_dim1;
#line 396 "sgesvx.f"
    af -= af_offset;
#line 396 "sgesvx.f"
    --ipiv;
#line 396 "sgesvx.f"
    --r__;
#line 396 "sgesvx.f"
    --c__;
#line 396 "sgesvx.f"
    b_dim1 = *ldb;
#line 396 "sgesvx.f"
    b_offset = 1 + b_dim1;
#line 396 "sgesvx.f"
    b -= b_offset;
#line 396 "sgesvx.f"
    x_dim1 = *ldx;
#line 396 "sgesvx.f"
    x_offset = 1 + x_dim1;
#line 396 "sgesvx.f"
    x -= x_offset;
#line 396 "sgesvx.f"
    --ferr;
#line 396 "sgesvx.f"
    --berr;
#line 396 "sgesvx.f"
    --work;
#line 396 "sgesvx.f"
    --iwork;
#line 396 "sgesvx.f"

#line 396 "sgesvx.f"
    /* Function Body */
#line 396 "sgesvx.f"
    *info = 0;
#line 397 "sgesvx.f"
    nofact = lsame_(fact, "N", (ftnlen)1, (ftnlen)1);
#line 398 "sgesvx.f"
    equil = lsame_(fact, "E", (ftnlen)1, (ftnlen)1);
#line 399 "sgesvx.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 400 "sgesvx.f"
    if (nofact || equil) {
#line 401 "sgesvx.f"
	*(unsigned char *)equed = 'N';
#line 402 "sgesvx.f"
	rowequ = FALSE_;
#line 403 "sgesvx.f"
	colequ = FALSE_;
#line 404 "sgesvx.f"
    } else {
#line 405 "sgesvx.f"
	rowequ = lsame_(equed, "R", (ftnlen)1, (ftnlen)1) || lsame_(equed, 
		"B", (ftnlen)1, (ftnlen)1);
#line 406 "sgesvx.f"
	colequ = lsame_(equed, "C", (ftnlen)1, (ftnlen)1) || lsame_(equed, 
		"B", (ftnlen)1, (ftnlen)1);
#line 407 "sgesvx.f"
	smlnum = slamch_("Safe minimum", (ftnlen)12);
#line 408 "sgesvx.f"
	bignum = 1. / smlnum;
#line 409 "sgesvx.f"
    }

/*     Test the input parameters. */

#line 413 "sgesvx.f"
    if (! nofact && ! equil && ! lsame_(fact, "F", (ftnlen)1, (ftnlen)1)) {
#line 415 "sgesvx.f"
	*info = -1;
#line 416 "sgesvx.f"
    } else if (! notran && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1) && ! 
	    lsame_(trans, "C", (ftnlen)1, (ftnlen)1)) {
#line 418 "sgesvx.f"
	*info = -2;
#line 419 "sgesvx.f"
    } else if (*n < 0) {
#line 420 "sgesvx.f"
	*info = -3;
#line 421 "sgesvx.f"
    } else if (*nrhs < 0) {
#line 422 "sgesvx.f"
	*info = -4;
#line 423 "sgesvx.f"
    } else if (*lda < max(1,*n)) {
#line 424 "sgesvx.f"
	*info = -6;
#line 425 "sgesvx.f"
    } else if (*ldaf < max(1,*n)) {
#line 426 "sgesvx.f"
	*info = -8;
#line 427 "sgesvx.f"
    } else if (lsame_(fact, "F", (ftnlen)1, (ftnlen)1) && ! (rowequ || colequ 
	    || lsame_(equed, "N", (ftnlen)1, (ftnlen)1))) {
#line 429 "sgesvx.f"
	*info = -10;
#line 430 "sgesvx.f"
    } else {
#line 431 "sgesvx.f"
	if (rowequ) {
#line 432 "sgesvx.f"
	    rcmin = bignum;
#line 433 "sgesvx.f"
	    rcmax = 0.;
#line 434 "sgesvx.f"
	    i__1 = *n;
#line 434 "sgesvx.f"
	    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 435 "sgesvx.f"
		d__1 = rcmin, d__2 = r__[j];
#line 435 "sgesvx.f"
		rcmin = min(d__1,d__2);
/* Computing MAX */
#line 436 "sgesvx.f"
		d__1 = rcmax, d__2 = r__[j];
#line 436 "sgesvx.f"
		rcmax = max(d__1,d__2);
#line 437 "sgesvx.f"
/* L10: */
#line 437 "sgesvx.f"
	    }
#line 438 "sgesvx.f"
	    if (rcmin <= 0.) {
#line 439 "sgesvx.f"
		*info = -11;
#line 440 "sgesvx.f"
	    } else if (*n > 0) {
#line 441 "sgesvx.f"
		rowcnd = max(rcmin,smlnum) / min(rcmax,bignum);
#line 442 "sgesvx.f"
	    } else {
#line 443 "sgesvx.f"
		rowcnd = 1.;
#line 444 "sgesvx.f"
	    }
#line 445 "sgesvx.f"
	}
#line 446 "sgesvx.f"
	if (colequ && *info == 0) {
#line 447 "sgesvx.f"
	    rcmin = bignum;
#line 448 "sgesvx.f"
	    rcmax = 0.;
#line 449 "sgesvx.f"
	    i__1 = *n;
#line 449 "sgesvx.f"
	    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 450 "sgesvx.f"
		d__1 = rcmin, d__2 = c__[j];
#line 450 "sgesvx.f"
		rcmin = min(d__1,d__2);
/* Computing MAX */
#line 451 "sgesvx.f"
		d__1 = rcmax, d__2 = c__[j];
#line 451 "sgesvx.f"
		rcmax = max(d__1,d__2);
#line 452 "sgesvx.f"
/* L20: */
#line 452 "sgesvx.f"
	    }
#line 453 "sgesvx.f"
	    if (rcmin <= 0.) {
#line 454 "sgesvx.f"
		*info = -12;
#line 455 "sgesvx.f"
	    } else if (*n > 0) {
#line 456 "sgesvx.f"
		colcnd = max(rcmin,smlnum) / min(rcmax,bignum);
#line 457 "sgesvx.f"
	    } else {
#line 458 "sgesvx.f"
		colcnd = 1.;
#line 459 "sgesvx.f"
	    }
#line 460 "sgesvx.f"
	}
#line 461 "sgesvx.f"
	if (*info == 0) {
#line 462 "sgesvx.f"
	    if (*ldb < max(1,*n)) {
#line 463 "sgesvx.f"
		*info = -14;
#line 464 "sgesvx.f"
	    } else if (*ldx < max(1,*n)) {
#line 465 "sgesvx.f"
		*info = -16;
#line 466 "sgesvx.f"
	    }
#line 467 "sgesvx.f"
	}
#line 468 "sgesvx.f"
    }

#line 470 "sgesvx.f"
    if (*info != 0) {
#line 471 "sgesvx.f"
	i__1 = -(*info);
#line 471 "sgesvx.f"
	xerbla_("SGESVX", &i__1, (ftnlen)6);
#line 472 "sgesvx.f"
	return 0;
#line 473 "sgesvx.f"
    }

#line 475 "sgesvx.f"
    if (equil) {

/*        Compute row and column scalings to equilibrate the matrix A. */

#line 479 "sgesvx.f"
	sgeequ_(n, n, &a[a_offset], lda, &r__[1], &c__[1], &rowcnd, &colcnd, &
		amax, &infequ);
#line 480 "sgesvx.f"
	if (infequ == 0) {

/*           Equilibrate the matrix. */

#line 484 "sgesvx.f"
	    slaqge_(n, n, &a[a_offset], lda, &r__[1], &c__[1], &rowcnd, &
		    colcnd, &amax, equed, (ftnlen)1);
#line 486 "sgesvx.f"
	    rowequ = lsame_(equed, "R", (ftnlen)1, (ftnlen)1) || lsame_(equed,
		     "B", (ftnlen)1, (ftnlen)1);
#line 487 "sgesvx.f"
	    colequ = lsame_(equed, "C", (ftnlen)1, (ftnlen)1) || lsame_(equed,
		     "B", (ftnlen)1, (ftnlen)1);
#line 488 "sgesvx.f"
	}
#line 489 "sgesvx.f"
    }

/*     Scale the right hand side. */

#line 493 "sgesvx.f"
    if (notran) {
#line 494 "sgesvx.f"
	if (rowequ) {
#line 495 "sgesvx.f"
	    i__1 = *nrhs;
#line 495 "sgesvx.f"
	    for (j = 1; j <= i__1; ++j) {
#line 496 "sgesvx.f"
		i__2 = *n;
#line 496 "sgesvx.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 497 "sgesvx.f"
		    b[i__ + j * b_dim1] = r__[i__] * b[i__ + j * b_dim1];
#line 498 "sgesvx.f"
/* L30: */
#line 498 "sgesvx.f"
		}
#line 499 "sgesvx.f"
/* L40: */
#line 499 "sgesvx.f"
	    }
#line 500 "sgesvx.f"
	}
#line 501 "sgesvx.f"
    } else if (colequ) {
#line 502 "sgesvx.f"
	i__1 = *nrhs;
#line 502 "sgesvx.f"
	for (j = 1; j <= i__1; ++j) {
#line 503 "sgesvx.f"
	    i__2 = *n;
#line 503 "sgesvx.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 504 "sgesvx.f"
		b[i__ + j * b_dim1] = c__[i__] * b[i__ + j * b_dim1];
#line 505 "sgesvx.f"
/* L50: */
#line 505 "sgesvx.f"
	    }
#line 506 "sgesvx.f"
/* L60: */
#line 506 "sgesvx.f"
	}
#line 507 "sgesvx.f"
    }

#line 509 "sgesvx.f"
    if (nofact || equil) {

/*        Compute the LU factorization of A. */

#line 513 "sgesvx.f"
	slacpy_("Full", n, n, &a[a_offset], lda, &af[af_offset], ldaf, (
		ftnlen)4);
#line 514 "sgesvx.f"
	sgetrf_(n, n, &af[af_offset], ldaf, &ipiv[1], info);

/*        Return if INFO is non-zero. */

#line 518 "sgesvx.f"
	if (*info > 0) {

/*           Compute the reciprocal pivot growth factor of the */
/*           leading rank-deficient INFO columns of A. */

#line 523 "sgesvx.f"
	    rpvgrw = slantr_("M", "U", "N", info, info, &af[af_offset], ldaf, 
		    &work[1], (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 525 "sgesvx.f"
	    if (rpvgrw == 0.) {
#line 526 "sgesvx.f"
		rpvgrw = 1.;
#line 527 "sgesvx.f"
	    } else {
#line 528 "sgesvx.f"
		rpvgrw = slange_("M", n, info, &a[a_offset], lda, &work[1], (
			ftnlen)1) / rpvgrw;
#line 529 "sgesvx.f"
	    }
#line 530 "sgesvx.f"
	    work[1] = rpvgrw;
#line 531 "sgesvx.f"
	    *rcond = 0.;
#line 532 "sgesvx.f"
	    return 0;
#line 533 "sgesvx.f"
	}
#line 534 "sgesvx.f"
    }

/*     Compute the norm of the matrix A and the */
/*     reciprocal pivot growth factor RPVGRW. */

#line 539 "sgesvx.f"
    if (notran) {
#line 540 "sgesvx.f"
	*(unsigned char *)norm = '1';
#line 541 "sgesvx.f"
    } else {
#line 542 "sgesvx.f"
	*(unsigned char *)norm = 'I';
#line 543 "sgesvx.f"
    }
#line 544 "sgesvx.f"
    anorm = slange_(norm, n, n, &a[a_offset], lda, &work[1], (ftnlen)1);
#line 545 "sgesvx.f"
    rpvgrw = slantr_("M", "U", "N", n, n, &af[af_offset], ldaf, &work[1], (
	    ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 546 "sgesvx.f"
    if (rpvgrw == 0.) {
#line 547 "sgesvx.f"
	rpvgrw = 1.;
#line 548 "sgesvx.f"
    } else {
#line 549 "sgesvx.f"
	rpvgrw = slange_("M", n, n, &a[a_offset], lda, &work[1], (ftnlen)1) / 
		rpvgrw;
#line 550 "sgesvx.f"
    }

/*     Compute the reciprocal of the condition number of A. */

#line 554 "sgesvx.f"
    sgecon_(norm, n, &af[af_offset], ldaf, &anorm, rcond, &work[1], &iwork[1],
	     info, (ftnlen)1);

/*     Compute the solution matrix X. */

#line 558 "sgesvx.f"
    slacpy_("Full", n, nrhs, &b[b_offset], ldb, &x[x_offset], ldx, (ftnlen)4);
#line 559 "sgesvx.f"
    sgetrs_(trans, n, nrhs, &af[af_offset], ldaf, &ipiv[1], &x[x_offset], ldx,
	     info, (ftnlen)1);

/*     Use iterative refinement to improve the computed solution and */
/*     compute error bounds and backward error estimates for it. */

#line 564 "sgesvx.f"
    sgerfs_(trans, n, nrhs, &a[a_offset], lda, &af[af_offset], ldaf, &ipiv[1],
	     &b[b_offset], ldb, &x[x_offset], ldx, &ferr[1], &berr[1], &work[
	    1], &iwork[1], info, (ftnlen)1);

/*     Transform the solution matrix X to a solution of the original */
/*     system. */

#line 570 "sgesvx.f"
    if (notran) {
#line 571 "sgesvx.f"
	if (colequ) {
#line 572 "sgesvx.f"
	    i__1 = *nrhs;
#line 572 "sgesvx.f"
	    for (j = 1; j <= i__1; ++j) {
#line 573 "sgesvx.f"
		i__2 = *n;
#line 573 "sgesvx.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 574 "sgesvx.f"
		    x[i__ + j * x_dim1] = c__[i__] * x[i__ + j * x_dim1];
#line 575 "sgesvx.f"
/* L70: */
#line 575 "sgesvx.f"
		}
#line 576 "sgesvx.f"
/* L80: */
#line 576 "sgesvx.f"
	    }
#line 577 "sgesvx.f"
	    i__1 = *nrhs;
#line 577 "sgesvx.f"
	    for (j = 1; j <= i__1; ++j) {
#line 578 "sgesvx.f"
		ferr[j] /= colcnd;
#line 579 "sgesvx.f"
/* L90: */
#line 579 "sgesvx.f"
	    }
#line 580 "sgesvx.f"
	}
#line 581 "sgesvx.f"
    } else if (rowequ) {
#line 582 "sgesvx.f"
	i__1 = *nrhs;
#line 582 "sgesvx.f"
	for (j = 1; j <= i__1; ++j) {
#line 583 "sgesvx.f"
	    i__2 = *n;
#line 583 "sgesvx.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 584 "sgesvx.f"
		x[i__ + j * x_dim1] = r__[i__] * x[i__ + j * x_dim1];
#line 585 "sgesvx.f"
/* L100: */
#line 585 "sgesvx.f"
	    }
#line 586 "sgesvx.f"
/* L110: */
#line 586 "sgesvx.f"
	}
#line 587 "sgesvx.f"
	i__1 = *nrhs;
#line 587 "sgesvx.f"
	for (j = 1; j <= i__1; ++j) {
#line 588 "sgesvx.f"
	    ferr[j] /= rowcnd;
#line 589 "sgesvx.f"
/* L120: */
#line 589 "sgesvx.f"
	}
#line 590 "sgesvx.f"
    }

/*     Set INFO = N+1 if the matrix is singular to working precision. */

#line 594 "sgesvx.f"
    if (*rcond < slamch_("Epsilon", (ftnlen)7)) {
#line 594 "sgesvx.f"
	*info = *n + 1;
#line 594 "sgesvx.f"
    }

#line 597 "sgesvx.f"
    work[1] = rpvgrw;
#line 598 "sgesvx.f"
    return 0;

/*     End of SGESVX */

} /* sgesvx_ */

