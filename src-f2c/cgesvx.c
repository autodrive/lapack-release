#line 1 "cgesvx.f"
/* cgesvx.f -- translated by f2c (version 20100827).
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

#line 1 "cgesvx.f"
/* > \brief <b> CGESVX computes the solution to system of linear equations A * X = B for GE matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CGESVX + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgesvx.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgesvx.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgesvx.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CGESVX( FACT, TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV, */
/*                          EQUED, R, C, B, LDB, X, LDX, RCOND, FERR, BERR, */
/*                          WORK, RWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          EQUED, FACT, TRANS */
/*       INTEGER            INFO, LDA, LDAF, LDB, LDX, N, NRHS */
/*       REAL               RCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       REAL               BERR( * ), C( * ), FERR( * ), R( * ), */
/*      $                   RWORK( * ) */
/*       COMPLEX            A( LDA, * ), AF( LDAF, * ), B( LDB, * ), */
/*      $                   WORK( * ), X( LDX, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CGESVX uses the LU factorization to compute the solution to a complex */
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
/* >          = 'C':  A**H * X = B  (Conjugate transpose) */
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
/* >          A is COMPLEX array, dimension (LDA,N) */
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
/* >          AF is COMPLEX array, dimension (LDAF,N) */
/* >          If FACT = 'F', then AF is an input argument and on entry */
/* >          contains the factors L and U from the factorization */
/* >          A = P*L*U as computed by CGETRF.  If EQUED .ne. 'N', then */
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
/* >          as computed by CGETRF; row i of the matrix was interchanged */
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
/* >          B is COMPLEX array, dimension (LDB,NRHS) */
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
/* >          X is COMPLEX array, dimension (LDX,NRHS) */
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
/* >          WORK is COMPLEX array, dimension (2*N) */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is REAL array, dimension (2*N) */
/* >          On exit, RWORK(1) contains the reciprocal pivot growth */
/* >          factor norm(A)/norm(U). The "max absolute element" norm is */
/* >          used. If RWORK(1) is much less than 1, then the stability */
/* >          of the LU factorization of the (equilibrated) matrix A */
/* >          could be poor. This also means that the solution X, condition */
/* >          estimator RCOND, and forward error bound FERR could be */
/* >          unreliable. If factorization fails with 0<INFO<=N, then */
/* >          RWORK(1) contains the reciprocal pivot growth factor for the */
/* >          leading INFO columns of A. */
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

/* > \ingroup complexGEsolve */

/*  ===================================================================== */
/* Subroutine */ int cgesvx_(char *fact, char *trans, integer *n, integer *
	nrhs, doublecomplex *a, integer *lda, doublecomplex *af, integer *
	ldaf, integer *ipiv, char *equed, doublereal *r__, doublereal *c__, 
	doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx, 
	doublereal *rcond, doublereal *ferr, doublereal *berr, doublecomplex *
	work, doublereal *rwork, integer *info, ftnlen fact_len, ftnlen 
	trans_len, ftnlen equed_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, af_dim1, af_offset, b_dim1, b_offset, x_dim1, 
	    x_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2;
    doublecomplex z__1;

    /* Local variables */
    static integer i__, j;
    static doublereal amax;
    static char norm[1];
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal rcmin, rcmax, anorm;
    static logical equil;
    extern doublereal clange_(char *, integer *, integer *, doublecomplex *, 
	    integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int claqge_(integer *, integer *, doublecomplex *,
	     integer *, doublereal *, doublereal *, doublereal *, doublereal *
	    , doublereal *, char *, ftnlen), cgecon_(char *, integer *, 
	    doublecomplex *, integer *, doublereal *, doublereal *, 
	    doublecomplex *, doublereal *, integer *, ftnlen);
    static doublereal colcnd;
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int cgeequ_(integer *, integer *, doublecomplex *,
	     integer *, doublereal *, doublereal *, doublereal *, doublereal *
	    , doublereal *, integer *);
    static logical nofact;
    extern /* Subroutine */ int cgerfs_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublereal *, doublereal *, doublecomplex *, doublereal *, 
	    integer *, ftnlen), cgetrf_(integer *, integer *, doublecomplex *,
	     integer *, integer *, integer *), clacpy_(char *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *,
	     ftnlen), xerbla_(char *, integer *, ftnlen);
    static doublereal bignum;
    extern doublereal clantr_(char *, char *, char *, integer *, integer *, 
	    doublecomplex *, integer *, doublereal *, ftnlen, ftnlen, ftnlen);
    static integer infequ;
    static logical colequ;
    extern /* Subroutine */ int cgetrs_(char *, integer *, integer *, 
	    doublecomplex *, integer *, integer *, doublecomplex *, integer *,
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

#line 398 "cgesvx.f"
    /* Parameter adjustments */
#line 398 "cgesvx.f"
    a_dim1 = *lda;
#line 398 "cgesvx.f"
    a_offset = 1 + a_dim1;
#line 398 "cgesvx.f"
    a -= a_offset;
#line 398 "cgesvx.f"
    af_dim1 = *ldaf;
#line 398 "cgesvx.f"
    af_offset = 1 + af_dim1;
#line 398 "cgesvx.f"
    af -= af_offset;
#line 398 "cgesvx.f"
    --ipiv;
#line 398 "cgesvx.f"
    --r__;
#line 398 "cgesvx.f"
    --c__;
#line 398 "cgesvx.f"
    b_dim1 = *ldb;
#line 398 "cgesvx.f"
    b_offset = 1 + b_dim1;
#line 398 "cgesvx.f"
    b -= b_offset;
#line 398 "cgesvx.f"
    x_dim1 = *ldx;
#line 398 "cgesvx.f"
    x_offset = 1 + x_dim1;
#line 398 "cgesvx.f"
    x -= x_offset;
#line 398 "cgesvx.f"
    --ferr;
#line 398 "cgesvx.f"
    --berr;
#line 398 "cgesvx.f"
    --work;
#line 398 "cgesvx.f"
    --rwork;
#line 398 "cgesvx.f"

#line 398 "cgesvx.f"
    /* Function Body */
#line 398 "cgesvx.f"
    *info = 0;
#line 399 "cgesvx.f"
    nofact = lsame_(fact, "N", (ftnlen)1, (ftnlen)1);
#line 400 "cgesvx.f"
    equil = lsame_(fact, "E", (ftnlen)1, (ftnlen)1);
#line 401 "cgesvx.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 402 "cgesvx.f"
    if (nofact || equil) {
#line 403 "cgesvx.f"
	*(unsigned char *)equed = 'N';
#line 404 "cgesvx.f"
	rowequ = FALSE_;
#line 405 "cgesvx.f"
	colequ = FALSE_;
#line 406 "cgesvx.f"
    } else {
#line 407 "cgesvx.f"
	rowequ = lsame_(equed, "R", (ftnlen)1, (ftnlen)1) || lsame_(equed, 
		"B", (ftnlen)1, (ftnlen)1);
#line 408 "cgesvx.f"
	colequ = lsame_(equed, "C", (ftnlen)1, (ftnlen)1) || lsame_(equed, 
		"B", (ftnlen)1, (ftnlen)1);
#line 409 "cgesvx.f"
	smlnum = slamch_("Safe minimum", (ftnlen)12);
#line 410 "cgesvx.f"
	bignum = 1. / smlnum;
#line 411 "cgesvx.f"
    }

/*     Test the input parameters. */

#line 415 "cgesvx.f"
    if (! nofact && ! equil && ! lsame_(fact, "F", (ftnlen)1, (ftnlen)1)) {
#line 417 "cgesvx.f"
	*info = -1;
#line 418 "cgesvx.f"
    } else if (! notran && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1) && ! 
	    lsame_(trans, "C", (ftnlen)1, (ftnlen)1)) {
#line 420 "cgesvx.f"
	*info = -2;
#line 421 "cgesvx.f"
    } else if (*n < 0) {
#line 422 "cgesvx.f"
	*info = -3;
#line 423 "cgesvx.f"
    } else if (*nrhs < 0) {
#line 424 "cgesvx.f"
	*info = -4;
#line 425 "cgesvx.f"
    } else if (*lda < max(1,*n)) {
#line 426 "cgesvx.f"
	*info = -6;
#line 427 "cgesvx.f"
    } else if (*ldaf < max(1,*n)) {
#line 428 "cgesvx.f"
	*info = -8;
#line 429 "cgesvx.f"
    } else if (lsame_(fact, "F", (ftnlen)1, (ftnlen)1) && ! (rowequ || colequ 
	    || lsame_(equed, "N", (ftnlen)1, (ftnlen)1))) {
#line 431 "cgesvx.f"
	*info = -10;
#line 432 "cgesvx.f"
    } else {
#line 433 "cgesvx.f"
	if (rowequ) {
#line 434 "cgesvx.f"
	    rcmin = bignum;
#line 435 "cgesvx.f"
	    rcmax = 0.;
#line 436 "cgesvx.f"
	    i__1 = *n;
#line 436 "cgesvx.f"
	    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 437 "cgesvx.f"
		d__1 = rcmin, d__2 = r__[j];
#line 437 "cgesvx.f"
		rcmin = min(d__1,d__2);
/* Computing MAX */
#line 438 "cgesvx.f"
		d__1 = rcmax, d__2 = r__[j];
#line 438 "cgesvx.f"
		rcmax = max(d__1,d__2);
#line 439 "cgesvx.f"
/* L10: */
#line 439 "cgesvx.f"
	    }
#line 440 "cgesvx.f"
	    if (rcmin <= 0.) {
#line 441 "cgesvx.f"
		*info = -11;
#line 442 "cgesvx.f"
	    } else if (*n > 0) {
#line 443 "cgesvx.f"
		rowcnd = max(rcmin,smlnum) / min(rcmax,bignum);
#line 444 "cgesvx.f"
	    } else {
#line 445 "cgesvx.f"
		rowcnd = 1.;
#line 446 "cgesvx.f"
	    }
#line 447 "cgesvx.f"
	}
#line 448 "cgesvx.f"
	if (colequ && *info == 0) {
#line 449 "cgesvx.f"
	    rcmin = bignum;
#line 450 "cgesvx.f"
	    rcmax = 0.;
#line 451 "cgesvx.f"
	    i__1 = *n;
#line 451 "cgesvx.f"
	    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 452 "cgesvx.f"
		d__1 = rcmin, d__2 = c__[j];
#line 452 "cgesvx.f"
		rcmin = min(d__1,d__2);
/* Computing MAX */
#line 453 "cgesvx.f"
		d__1 = rcmax, d__2 = c__[j];
#line 453 "cgesvx.f"
		rcmax = max(d__1,d__2);
#line 454 "cgesvx.f"
/* L20: */
#line 454 "cgesvx.f"
	    }
#line 455 "cgesvx.f"
	    if (rcmin <= 0.) {
#line 456 "cgesvx.f"
		*info = -12;
#line 457 "cgesvx.f"
	    } else if (*n > 0) {
#line 458 "cgesvx.f"
		colcnd = max(rcmin,smlnum) / min(rcmax,bignum);
#line 459 "cgesvx.f"
	    } else {
#line 460 "cgesvx.f"
		colcnd = 1.;
#line 461 "cgesvx.f"
	    }
#line 462 "cgesvx.f"
	}
#line 463 "cgesvx.f"
	if (*info == 0) {
#line 464 "cgesvx.f"
	    if (*ldb < max(1,*n)) {
#line 465 "cgesvx.f"
		*info = -14;
#line 466 "cgesvx.f"
	    } else if (*ldx < max(1,*n)) {
#line 467 "cgesvx.f"
		*info = -16;
#line 468 "cgesvx.f"
	    }
#line 469 "cgesvx.f"
	}
#line 470 "cgesvx.f"
    }

#line 472 "cgesvx.f"
    if (*info != 0) {
#line 473 "cgesvx.f"
	i__1 = -(*info);
#line 473 "cgesvx.f"
	xerbla_("CGESVX", &i__1, (ftnlen)6);
#line 474 "cgesvx.f"
	return 0;
#line 475 "cgesvx.f"
    }

#line 477 "cgesvx.f"
    if (equil) {

/*        Compute row and column scalings to equilibrate the matrix A. */

#line 481 "cgesvx.f"
	cgeequ_(n, n, &a[a_offset], lda, &r__[1], &c__[1], &rowcnd, &colcnd, &
		amax, &infequ);
#line 482 "cgesvx.f"
	if (infequ == 0) {

/*           Equilibrate the matrix. */

#line 486 "cgesvx.f"
	    claqge_(n, n, &a[a_offset], lda, &r__[1], &c__[1], &rowcnd, &
		    colcnd, &amax, equed, (ftnlen)1);
#line 488 "cgesvx.f"
	    rowequ = lsame_(equed, "R", (ftnlen)1, (ftnlen)1) || lsame_(equed,
		     "B", (ftnlen)1, (ftnlen)1);
#line 489 "cgesvx.f"
	    colequ = lsame_(equed, "C", (ftnlen)1, (ftnlen)1) || lsame_(equed,
		     "B", (ftnlen)1, (ftnlen)1);
#line 490 "cgesvx.f"
	}
#line 491 "cgesvx.f"
    }

/*     Scale the right hand side. */

#line 495 "cgesvx.f"
    if (notran) {
#line 496 "cgesvx.f"
	if (rowequ) {
#line 497 "cgesvx.f"
	    i__1 = *nrhs;
#line 497 "cgesvx.f"
	    for (j = 1; j <= i__1; ++j) {
#line 498 "cgesvx.f"
		i__2 = *n;
#line 498 "cgesvx.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 499 "cgesvx.f"
		    i__3 = i__ + j * b_dim1;
#line 499 "cgesvx.f"
		    i__4 = i__;
#line 499 "cgesvx.f"
		    i__5 = i__ + j * b_dim1;
#line 499 "cgesvx.f"
		    z__1.r = r__[i__4] * b[i__5].r, z__1.i = r__[i__4] * b[
			    i__5].i;
#line 499 "cgesvx.f"
		    b[i__3].r = z__1.r, b[i__3].i = z__1.i;
#line 500 "cgesvx.f"
/* L30: */
#line 500 "cgesvx.f"
		}
#line 501 "cgesvx.f"
/* L40: */
#line 501 "cgesvx.f"
	    }
#line 502 "cgesvx.f"
	}
#line 503 "cgesvx.f"
    } else if (colequ) {
#line 504 "cgesvx.f"
	i__1 = *nrhs;
#line 504 "cgesvx.f"
	for (j = 1; j <= i__1; ++j) {
#line 505 "cgesvx.f"
	    i__2 = *n;
#line 505 "cgesvx.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 506 "cgesvx.f"
		i__3 = i__ + j * b_dim1;
#line 506 "cgesvx.f"
		i__4 = i__;
#line 506 "cgesvx.f"
		i__5 = i__ + j * b_dim1;
#line 506 "cgesvx.f"
		z__1.r = c__[i__4] * b[i__5].r, z__1.i = c__[i__4] * b[i__5]
			.i;
#line 506 "cgesvx.f"
		b[i__3].r = z__1.r, b[i__3].i = z__1.i;
#line 507 "cgesvx.f"
/* L50: */
#line 507 "cgesvx.f"
	    }
#line 508 "cgesvx.f"
/* L60: */
#line 508 "cgesvx.f"
	}
#line 509 "cgesvx.f"
    }

#line 511 "cgesvx.f"
    if (nofact || equil) {

/*        Compute the LU factorization of A. */

#line 515 "cgesvx.f"
	clacpy_("Full", n, n, &a[a_offset], lda, &af[af_offset], ldaf, (
		ftnlen)4);
#line 516 "cgesvx.f"
	cgetrf_(n, n, &af[af_offset], ldaf, &ipiv[1], info);

/*        Return if INFO is non-zero. */

#line 520 "cgesvx.f"
	if (*info > 0) {

/*           Compute the reciprocal pivot growth factor of the */
/*           leading rank-deficient INFO columns of A. */

#line 525 "cgesvx.f"
	    rpvgrw = clantr_("M", "U", "N", info, info, &af[af_offset], ldaf, 
		    &rwork[1], (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 527 "cgesvx.f"
	    if (rpvgrw == 0.) {
#line 528 "cgesvx.f"
		rpvgrw = 1.;
#line 529 "cgesvx.f"
	    } else {
#line 530 "cgesvx.f"
		rpvgrw = clange_("M", n, info, &a[a_offset], lda, &rwork[1], (
			ftnlen)1) / rpvgrw;
#line 532 "cgesvx.f"
	    }
#line 533 "cgesvx.f"
	    rwork[1] = rpvgrw;
#line 534 "cgesvx.f"
	    *rcond = 0.;
#line 535 "cgesvx.f"
	    return 0;
#line 536 "cgesvx.f"
	}
#line 537 "cgesvx.f"
    }

/*     Compute the norm of the matrix A and the */
/*     reciprocal pivot growth factor RPVGRW. */

#line 542 "cgesvx.f"
    if (notran) {
#line 543 "cgesvx.f"
	*(unsigned char *)norm = '1';
#line 544 "cgesvx.f"
    } else {
#line 545 "cgesvx.f"
	*(unsigned char *)norm = 'I';
#line 546 "cgesvx.f"
    }
#line 547 "cgesvx.f"
    anorm = clange_(norm, n, n, &a[a_offset], lda, &rwork[1], (ftnlen)1);
#line 548 "cgesvx.f"
    rpvgrw = clantr_("M", "U", "N", n, n, &af[af_offset], ldaf, &rwork[1], (
	    ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 549 "cgesvx.f"
    if (rpvgrw == 0.) {
#line 550 "cgesvx.f"
	rpvgrw = 1.;
#line 551 "cgesvx.f"
    } else {
#line 552 "cgesvx.f"
	rpvgrw = clange_("M", n, n, &a[a_offset], lda, &rwork[1], (ftnlen)1) /
		 rpvgrw;
#line 553 "cgesvx.f"
    }

/*     Compute the reciprocal of the condition number of A. */

#line 557 "cgesvx.f"
    cgecon_(norm, n, &af[af_offset], ldaf, &anorm, rcond, &work[1], &rwork[1],
	     info, (ftnlen)1);

/*     Compute the solution matrix X. */

#line 561 "cgesvx.f"
    clacpy_("Full", n, nrhs, &b[b_offset], ldb, &x[x_offset], ldx, (ftnlen)4);
#line 562 "cgesvx.f"
    cgetrs_(trans, n, nrhs, &af[af_offset], ldaf, &ipiv[1], &x[x_offset], ldx,
	     info, (ftnlen)1);

/*     Use iterative refinement to improve the computed solution and */
/*     compute error bounds and backward error estimates for it. */

#line 567 "cgesvx.f"
    cgerfs_(trans, n, nrhs, &a[a_offset], lda, &af[af_offset], ldaf, &ipiv[1],
	     &b[b_offset], ldb, &x[x_offset], ldx, &ferr[1], &berr[1], &work[
	    1], &rwork[1], info, (ftnlen)1);

/*     Transform the solution matrix X to a solution of the original */
/*     system. */

#line 573 "cgesvx.f"
    if (notran) {
#line 574 "cgesvx.f"
	if (colequ) {
#line 575 "cgesvx.f"
	    i__1 = *nrhs;
#line 575 "cgesvx.f"
	    for (j = 1; j <= i__1; ++j) {
#line 576 "cgesvx.f"
		i__2 = *n;
#line 576 "cgesvx.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 577 "cgesvx.f"
		    i__3 = i__ + j * x_dim1;
#line 577 "cgesvx.f"
		    i__4 = i__;
#line 577 "cgesvx.f"
		    i__5 = i__ + j * x_dim1;
#line 577 "cgesvx.f"
		    z__1.r = c__[i__4] * x[i__5].r, z__1.i = c__[i__4] * x[
			    i__5].i;
#line 577 "cgesvx.f"
		    x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 578 "cgesvx.f"
/* L70: */
#line 578 "cgesvx.f"
		}
#line 579 "cgesvx.f"
/* L80: */
#line 579 "cgesvx.f"
	    }
#line 580 "cgesvx.f"
	    i__1 = *nrhs;
#line 580 "cgesvx.f"
	    for (j = 1; j <= i__1; ++j) {
#line 581 "cgesvx.f"
		ferr[j] /= colcnd;
#line 582 "cgesvx.f"
/* L90: */
#line 582 "cgesvx.f"
	    }
#line 583 "cgesvx.f"
	}
#line 584 "cgesvx.f"
    } else if (rowequ) {
#line 585 "cgesvx.f"
	i__1 = *nrhs;
#line 585 "cgesvx.f"
	for (j = 1; j <= i__1; ++j) {
#line 586 "cgesvx.f"
	    i__2 = *n;
#line 586 "cgesvx.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 587 "cgesvx.f"
		i__3 = i__ + j * x_dim1;
#line 587 "cgesvx.f"
		i__4 = i__;
#line 587 "cgesvx.f"
		i__5 = i__ + j * x_dim1;
#line 587 "cgesvx.f"
		z__1.r = r__[i__4] * x[i__5].r, z__1.i = r__[i__4] * x[i__5]
			.i;
#line 587 "cgesvx.f"
		x[i__3].r = z__1.r, x[i__3].i = z__1.i;
#line 588 "cgesvx.f"
/* L100: */
#line 588 "cgesvx.f"
	    }
#line 589 "cgesvx.f"
/* L110: */
#line 589 "cgesvx.f"
	}
#line 590 "cgesvx.f"
	i__1 = *nrhs;
#line 590 "cgesvx.f"
	for (j = 1; j <= i__1; ++j) {
#line 591 "cgesvx.f"
	    ferr[j] /= rowcnd;
#line 592 "cgesvx.f"
/* L120: */
#line 592 "cgesvx.f"
	}
#line 593 "cgesvx.f"
    }

/*     Set INFO = N+1 if the matrix is singular to working precision. */

#line 597 "cgesvx.f"
    if (*rcond < slamch_("Epsilon", (ftnlen)7)) {
#line 597 "cgesvx.f"
	*info = *n + 1;
#line 597 "cgesvx.f"
    }

#line 600 "cgesvx.f"
    rwork[1] = rpvgrw;
#line 601 "cgesvx.f"
    return 0;

/*     End of CGESVX */

} /* cgesvx_ */

