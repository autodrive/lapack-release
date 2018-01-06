#line 1 "sgbsvx.f"
/* sgbsvx.f -- translated by f2c (version 20100827).
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

#line 1 "sgbsvx.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief <b> SGBSVX computes the solution to system of linear equations A * X = B for GB matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SGBSVX + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgbsvx.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgbsvx.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgbsvx.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SGBSVX( FACT, TRANS, N, KL, KU, NRHS, AB, LDAB, AFB, */
/*                          LDAFB, IPIV, EQUED, R, C, B, LDB, X, LDX, */
/*                          RCOND, FERR, BERR, WORK, IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          EQUED, FACT, TRANS */
/*       INTEGER            INFO, KL, KU, LDAB, LDAFB, LDB, LDX, N, NRHS */
/*       REAL               RCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ), IWORK( * ) */
/*       REAL               AB( LDAB, * ), AFB( LDAFB, * ), B( LDB, * ), */
/*      $                   BERR( * ), C( * ), FERR( * ), R( * ), */
/*      $                   WORK( * ), X( LDX, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SGBSVX uses the LU factorization to compute the solution to a real */
/* > system of linear equations A * X = B, A**T * X = B, or A**H * X = B, */
/* > where A is a band matrix of order N with KL subdiagonals and KU */
/* > superdiagonals, and X and B are N-by-NRHS matrices. */
/* > */
/* > Error bounds on the solution and a condition estimate are also */
/* > provided. */
/* > \endverbatim */

/* > \par Description: */
/*  ================= */
/* > */
/* > \verbatim */
/* > */
/* > The following steps are performed by this subroutine: */
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
/* >       A = L * U, */
/* >    where L is a product of permutation and unit lower triangular */
/* >    matrices with KL subdiagonals, and U is upper triangular with */
/* >    KL+KU superdiagonals. */
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
/* >          = 'F':  On entry, AFB and IPIV contain the factored form of */
/* >                  A.  If EQUED is not 'N', the matrix A has been */
/* >                  equilibrated with scaling factors given by R and C. */
/* >                  AB, AFB, and IPIV are not modified. */
/* >          = 'N':  The matrix A will be copied to AFB and factored. */
/* >          = 'E':  The matrix A will be equilibrated if necessary, then */
/* >                  copied to AFB and factored. */
/* > \endverbatim */
/* > */
/* > \param[in] TRANS */
/* > \verbatim */
/* >          TRANS is CHARACTER*1 */
/* >          Specifies the form of the system of equations. */
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
/* > \param[in] KL */
/* > \verbatim */
/* >          KL is INTEGER */
/* >          The number of subdiagonals within the band of A.  KL >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] KU */
/* > \verbatim */
/* >          KU is INTEGER */
/* >          The number of superdiagonals within the band of A.  KU >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NRHS */
/* > \verbatim */
/* >          NRHS is INTEGER */
/* >          The number of right hand sides, i.e., the number of columns */
/* >          of the matrices B and X.  NRHS >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] AB */
/* > \verbatim */
/* >          AB is REAL array, dimension (LDAB,N) */
/* >          On entry, the matrix A in band storage, in rows 1 to KL+KU+1. */
/* >          The j-th column of A is stored in the j-th column of the */
/* >          array AB as follows: */
/* >          AB(KU+1+i-j,j) = A(i,j) for max(1,j-KU)<=i<=min(N,j+kl) */
/* > */
/* >          If FACT = 'F' and EQUED is not 'N', then A must have been */
/* >          equilibrated by the scaling factors in R and/or C.  AB is not */
/* >          modified if FACT = 'F' or 'N', or if FACT = 'E' and */
/* >          EQUED = 'N' on exit. */
/* > */
/* >          On exit, if EQUED .ne. 'N', A is scaled as follows: */
/* >          EQUED = 'R':  A := diag(R) * A */
/* >          EQUED = 'C':  A := A * diag(C) */
/* >          EQUED = 'B':  A := diag(R) * A * diag(C). */
/* > \endverbatim */
/* > */
/* > \param[in] LDAB */
/* > \verbatim */
/* >          LDAB is INTEGER */
/* >          The leading dimension of the array AB.  LDAB >= KL+KU+1. */
/* > \endverbatim */
/* > */
/* > \param[in,out] AFB */
/* > \verbatim */
/* >          AFB is REAL array, dimension (LDAFB,N) */
/* >          If FACT = 'F', then AFB is an input argument and on entry */
/* >          contains details of the LU factorization of the band matrix */
/* >          A, as computed by SGBTRF.  U is stored as an upper triangular */
/* >          band matrix with KL+KU superdiagonals in rows 1 to KL+KU+1, */
/* >          and the multipliers used during the factorization are stored */
/* >          in rows KL+KU+2 to 2*KL+KU+1.  If EQUED .ne. 'N', then AFB is */
/* >          the factored form of the equilibrated matrix A. */
/* > */
/* >          If FACT = 'N', then AFB is an output argument and on exit */
/* >          returns details of the LU factorization of A. */
/* > */
/* >          If FACT = 'E', then AFB is an output argument and on exit */
/* >          returns details of the LU factorization of the equilibrated */
/* >          matrix A (see the description of AB for the form of the */
/* >          equilibrated matrix). */
/* > \endverbatim */
/* > */
/* > \param[in] LDAFB */
/* > \verbatim */
/* >          LDAFB is INTEGER */
/* >          The leading dimension of the array AFB.  LDAFB >= 2*KL+KU+1. */
/* > \endverbatim */
/* > */
/* > \param[in,out] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N) */
/* >          If FACT = 'F', then IPIV is an input argument and on entry */
/* >          contains the pivot indices from the factorization A = L*U */
/* >          as computed by SGBTRF; row i of the matrix was interchanged */
/* >          with row IPIV(i). */
/* > */
/* >          If FACT = 'N', then IPIV is an output argument and on exit */
/* >          contains the pivot indices from the factorization A = L*U */
/* >          of the original matrix A. */
/* > */
/* >          If FACT = 'E', then IPIV is an output argument and on exit */
/* >          contains the pivot indices from the factorization A = L*U */
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
/* >          On entry, the right hand side matrix B. */
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
/* >          WORK is REAL array, dimension (3*N) */
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
/* >                <= N:  U(i,i) is exactly zero.  The factorization */
/* >                       has been completed, but the factor U is exactly */
/* >                       singular, so the solution and error bounds */
/* >                       could not be computed. RCOND = 0 is returned. */
/* >                = N+1: U is nonsingular, but RCOND is less than machine */
/* >                       precision, meaning that the matrix is singular */
/* >                       to working precision.  Nevertheless, the */
/* >                       solution and error bounds are computed because */
/* >                       there are a number of situations where the */
/* >                       computed solution can be more accurate than the */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date April 2012 */

/* > \ingroup realGBsolve */

/*  ===================================================================== */
/* Subroutine */ int sgbsvx_(char *fact, char *trans, integer *n, integer *kl,
	 integer *ku, integer *nrhs, doublereal *ab, integer *ldab, 
	doublereal *afb, integer *ldafb, integer *ipiv, char *equed, 
	doublereal *r__, doublereal *c__, doublereal *b, integer *ldb, 
	doublereal *x, integer *ldx, doublereal *rcond, doublereal *ferr, 
	doublereal *berr, doublereal *work, integer *iwork, integer *info, 
	ftnlen fact_len, ftnlen trans_len, ftnlen equed_len)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, afb_dim1, afb_offset, b_dim1, b_offset, 
	    x_dim1, x_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2, d__3;

    /* Local variables */
    static integer i__, j, j1, j2;
    static doublereal amax;
    static char norm[1];
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal rcmin, rcmax, anorm;
    static logical equil;
    extern /* Subroutine */ int scopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static doublereal colcnd;
    extern doublereal slangb_(char *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, ftnlen), slamch_(char *, 
	    ftnlen);
    extern /* Subroutine */ int slaqgb_(integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, char *, ftnlen);
    static logical nofact;
    extern /* Subroutine */ int sgbcon_(char *, integer *, integer *, integer 
	    *, doublereal *, integer *, integer *, doublereal *, doublereal *,
	     doublereal *, integer *, integer *, ftnlen), xerbla_(char *, 
	    integer *, ftnlen);
    static doublereal bignum;
    extern doublereal slantb_(char *, char *, char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, ftnlen, ftnlen, ftnlen);
    extern /* Subroutine */ int sgbequ_(integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *);
    static integer infequ;
    static logical colequ;
    extern /* Subroutine */ int sgbrfs_(char *, integer *, integer *, integer 
	    *, integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    ftnlen), sgbtrf_(integer *, integer *, integer *, integer *, 
	    doublereal *, integer *, integer *, integer *), slacpy_(char *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, ftnlen);
    static doublereal rowcnd;
    static logical notran;
    extern /* Subroutine */ int sgbtrs_(char *, integer *, integer *, integer 
	    *, integer *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen);
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
/*  Moved setting of INFO = N+1 so INFO does not subsequently get */
/*  overwritten.  Sven, 17 Mar 05. */
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

#line 418 "sgbsvx.f"
    /* Parameter adjustments */
#line 418 "sgbsvx.f"
    ab_dim1 = *ldab;
#line 418 "sgbsvx.f"
    ab_offset = 1 + ab_dim1;
#line 418 "sgbsvx.f"
    ab -= ab_offset;
#line 418 "sgbsvx.f"
    afb_dim1 = *ldafb;
#line 418 "sgbsvx.f"
    afb_offset = 1 + afb_dim1;
#line 418 "sgbsvx.f"
    afb -= afb_offset;
#line 418 "sgbsvx.f"
    --ipiv;
#line 418 "sgbsvx.f"
    --r__;
#line 418 "sgbsvx.f"
    --c__;
#line 418 "sgbsvx.f"
    b_dim1 = *ldb;
#line 418 "sgbsvx.f"
    b_offset = 1 + b_dim1;
#line 418 "sgbsvx.f"
    b -= b_offset;
#line 418 "sgbsvx.f"
    x_dim1 = *ldx;
#line 418 "sgbsvx.f"
    x_offset = 1 + x_dim1;
#line 418 "sgbsvx.f"
    x -= x_offset;
#line 418 "sgbsvx.f"
    --ferr;
#line 418 "sgbsvx.f"
    --berr;
#line 418 "sgbsvx.f"
    --work;
#line 418 "sgbsvx.f"
    --iwork;
#line 418 "sgbsvx.f"

#line 418 "sgbsvx.f"
    /* Function Body */
#line 418 "sgbsvx.f"
    *info = 0;
#line 419 "sgbsvx.f"
    nofact = lsame_(fact, "N", (ftnlen)1, (ftnlen)1);
#line 420 "sgbsvx.f"
    equil = lsame_(fact, "E", (ftnlen)1, (ftnlen)1);
#line 421 "sgbsvx.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 422 "sgbsvx.f"
    if (nofact || equil) {
#line 423 "sgbsvx.f"
	*(unsigned char *)equed = 'N';
#line 424 "sgbsvx.f"
	rowequ = FALSE_;
#line 425 "sgbsvx.f"
	colequ = FALSE_;
#line 426 "sgbsvx.f"
    } else {
#line 427 "sgbsvx.f"
	rowequ = lsame_(equed, "R", (ftnlen)1, (ftnlen)1) || lsame_(equed, 
		"B", (ftnlen)1, (ftnlen)1);
#line 428 "sgbsvx.f"
	colequ = lsame_(equed, "C", (ftnlen)1, (ftnlen)1) || lsame_(equed, 
		"B", (ftnlen)1, (ftnlen)1);
#line 429 "sgbsvx.f"
	smlnum = slamch_("Safe minimum", (ftnlen)12);
#line 430 "sgbsvx.f"
	bignum = 1. / smlnum;
#line 431 "sgbsvx.f"
    }

/*     Test the input parameters. */

#line 435 "sgbsvx.f"
    if (! nofact && ! equil && ! lsame_(fact, "F", (ftnlen)1, (ftnlen)1)) {
#line 437 "sgbsvx.f"
	*info = -1;
#line 438 "sgbsvx.f"
    } else if (! notran && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1) && ! 
	    lsame_(trans, "C", (ftnlen)1, (ftnlen)1)) {
#line 440 "sgbsvx.f"
	*info = -2;
#line 441 "sgbsvx.f"
    } else if (*n < 0) {
#line 442 "sgbsvx.f"
	*info = -3;
#line 443 "sgbsvx.f"
    } else if (*kl < 0) {
#line 444 "sgbsvx.f"
	*info = -4;
#line 445 "sgbsvx.f"
    } else if (*ku < 0) {
#line 446 "sgbsvx.f"
	*info = -5;
#line 447 "sgbsvx.f"
    } else if (*nrhs < 0) {
#line 448 "sgbsvx.f"
	*info = -6;
#line 449 "sgbsvx.f"
    } else if (*ldab < *kl + *ku + 1) {
#line 450 "sgbsvx.f"
	*info = -8;
#line 451 "sgbsvx.f"
    } else if (*ldafb < (*kl << 1) + *ku + 1) {
#line 452 "sgbsvx.f"
	*info = -10;
#line 453 "sgbsvx.f"
    } else if (lsame_(fact, "F", (ftnlen)1, (ftnlen)1) && ! (rowequ || colequ 
	    || lsame_(equed, "N", (ftnlen)1, (ftnlen)1))) {
#line 455 "sgbsvx.f"
	*info = -12;
#line 456 "sgbsvx.f"
    } else {
#line 457 "sgbsvx.f"
	if (rowequ) {
#line 458 "sgbsvx.f"
	    rcmin = bignum;
#line 459 "sgbsvx.f"
	    rcmax = 0.;
#line 460 "sgbsvx.f"
	    i__1 = *n;
#line 460 "sgbsvx.f"
	    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 461 "sgbsvx.f"
		d__1 = rcmin, d__2 = r__[j];
#line 461 "sgbsvx.f"
		rcmin = min(d__1,d__2);
/* Computing MAX */
#line 462 "sgbsvx.f"
		d__1 = rcmax, d__2 = r__[j];
#line 462 "sgbsvx.f"
		rcmax = max(d__1,d__2);
#line 463 "sgbsvx.f"
/* L10: */
#line 463 "sgbsvx.f"
	    }
#line 464 "sgbsvx.f"
	    if (rcmin <= 0.) {
#line 465 "sgbsvx.f"
		*info = -13;
#line 466 "sgbsvx.f"
	    } else if (*n > 0) {
#line 467 "sgbsvx.f"
		rowcnd = max(rcmin,smlnum) / min(rcmax,bignum);
#line 468 "sgbsvx.f"
	    } else {
#line 469 "sgbsvx.f"
		rowcnd = 1.;
#line 470 "sgbsvx.f"
	    }
#line 471 "sgbsvx.f"
	}
#line 472 "sgbsvx.f"
	if (colequ && *info == 0) {
#line 473 "sgbsvx.f"
	    rcmin = bignum;
#line 474 "sgbsvx.f"
	    rcmax = 0.;
#line 475 "sgbsvx.f"
	    i__1 = *n;
#line 475 "sgbsvx.f"
	    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 476 "sgbsvx.f"
		d__1 = rcmin, d__2 = c__[j];
#line 476 "sgbsvx.f"
		rcmin = min(d__1,d__2);
/* Computing MAX */
#line 477 "sgbsvx.f"
		d__1 = rcmax, d__2 = c__[j];
#line 477 "sgbsvx.f"
		rcmax = max(d__1,d__2);
#line 478 "sgbsvx.f"
/* L20: */
#line 478 "sgbsvx.f"
	    }
#line 479 "sgbsvx.f"
	    if (rcmin <= 0.) {
#line 480 "sgbsvx.f"
		*info = -14;
#line 481 "sgbsvx.f"
	    } else if (*n > 0) {
#line 482 "sgbsvx.f"
		colcnd = max(rcmin,smlnum) / min(rcmax,bignum);
#line 483 "sgbsvx.f"
	    } else {
#line 484 "sgbsvx.f"
		colcnd = 1.;
#line 485 "sgbsvx.f"
	    }
#line 486 "sgbsvx.f"
	}
#line 487 "sgbsvx.f"
	if (*info == 0) {
#line 488 "sgbsvx.f"
	    if (*ldb < max(1,*n)) {
#line 489 "sgbsvx.f"
		*info = -16;
#line 490 "sgbsvx.f"
	    } else if (*ldx < max(1,*n)) {
#line 491 "sgbsvx.f"
		*info = -18;
#line 492 "sgbsvx.f"
	    }
#line 493 "sgbsvx.f"
	}
#line 494 "sgbsvx.f"
    }

#line 496 "sgbsvx.f"
    if (*info != 0) {
#line 497 "sgbsvx.f"
	i__1 = -(*info);
#line 497 "sgbsvx.f"
	xerbla_("SGBSVX", &i__1, (ftnlen)6);
#line 498 "sgbsvx.f"
	return 0;
#line 499 "sgbsvx.f"
    }

#line 501 "sgbsvx.f"
    if (equil) {

/*        Compute row and column scalings to equilibrate the matrix A. */

#line 505 "sgbsvx.f"
	sgbequ_(n, n, kl, ku, &ab[ab_offset], ldab, &r__[1], &c__[1], &rowcnd,
		 &colcnd, &amax, &infequ);
#line 507 "sgbsvx.f"
	if (infequ == 0) {

/*           Equilibrate the matrix. */

#line 511 "sgbsvx.f"
	    slaqgb_(n, n, kl, ku, &ab[ab_offset], ldab, &r__[1], &c__[1], &
		    rowcnd, &colcnd, &amax, equed, (ftnlen)1);
#line 513 "sgbsvx.f"
	    rowequ = lsame_(equed, "R", (ftnlen)1, (ftnlen)1) || lsame_(equed,
		     "B", (ftnlen)1, (ftnlen)1);
#line 514 "sgbsvx.f"
	    colequ = lsame_(equed, "C", (ftnlen)1, (ftnlen)1) || lsame_(equed,
		     "B", (ftnlen)1, (ftnlen)1);
#line 515 "sgbsvx.f"
	}
#line 516 "sgbsvx.f"
    }

/*     Scale the right hand side. */

#line 520 "sgbsvx.f"
    if (notran) {
#line 521 "sgbsvx.f"
	if (rowequ) {
#line 522 "sgbsvx.f"
	    i__1 = *nrhs;
#line 522 "sgbsvx.f"
	    for (j = 1; j <= i__1; ++j) {
#line 523 "sgbsvx.f"
		i__2 = *n;
#line 523 "sgbsvx.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 524 "sgbsvx.f"
		    b[i__ + j * b_dim1] = r__[i__] * b[i__ + j * b_dim1];
#line 525 "sgbsvx.f"
/* L30: */
#line 525 "sgbsvx.f"
		}
#line 526 "sgbsvx.f"
/* L40: */
#line 526 "sgbsvx.f"
	    }
#line 527 "sgbsvx.f"
	}
#line 528 "sgbsvx.f"
    } else if (colequ) {
#line 529 "sgbsvx.f"
	i__1 = *nrhs;
#line 529 "sgbsvx.f"
	for (j = 1; j <= i__1; ++j) {
#line 530 "sgbsvx.f"
	    i__2 = *n;
#line 530 "sgbsvx.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 531 "sgbsvx.f"
		b[i__ + j * b_dim1] = c__[i__] * b[i__ + j * b_dim1];
#line 532 "sgbsvx.f"
/* L50: */
#line 532 "sgbsvx.f"
	    }
#line 533 "sgbsvx.f"
/* L60: */
#line 533 "sgbsvx.f"
	}
#line 534 "sgbsvx.f"
    }

#line 536 "sgbsvx.f"
    if (nofact || equil) {

/*        Compute the LU factorization of the band matrix A. */

#line 540 "sgbsvx.f"
	i__1 = *n;
#line 540 "sgbsvx.f"
	for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
#line 541 "sgbsvx.f"
	    i__2 = j - *ku;
#line 541 "sgbsvx.f"
	    j1 = max(i__2,1);
/* Computing MIN */
#line 542 "sgbsvx.f"
	    i__2 = j + *kl;
#line 542 "sgbsvx.f"
	    j2 = min(i__2,*n);
#line 543 "sgbsvx.f"
	    i__2 = j2 - j1 + 1;
#line 543 "sgbsvx.f"
	    scopy_(&i__2, &ab[*ku + 1 - j + j1 + j * ab_dim1], &c__1, &afb[*
		    kl + *ku + 1 - j + j1 + j * afb_dim1], &c__1);
#line 545 "sgbsvx.f"
/* L70: */
#line 545 "sgbsvx.f"
	}

#line 547 "sgbsvx.f"
	sgbtrf_(n, n, kl, ku, &afb[afb_offset], ldafb, &ipiv[1], info);

/*        Return if INFO is non-zero. */

#line 551 "sgbsvx.f"
	if (*info > 0) {

/*           Compute the reciprocal pivot growth factor of the */
/*           leading rank-deficient INFO columns of A. */

#line 556 "sgbsvx.f"
	    anorm = 0.;
#line 557 "sgbsvx.f"
	    i__1 = *info;
#line 557 "sgbsvx.f"
	    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
#line 558 "sgbsvx.f"
		i__2 = *ku + 2 - j;
/* Computing MIN */
#line 558 "sgbsvx.f"
		i__4 = *n + *ku + 1 - j, i__5 = *kl + *ku + 1;
#line 558 "sgbsvx.f"
		i__3 = min(i__4,i__5);
#line 558 "sgbsvx.f"
		for (i__ = max(i__2,1); i__ <= i__3; ++i__) {
/* Computing MAX */
#line 559 "sgbsvx.f"
		    d__2 = anorm, d__3 = (d__1 = ab[i__ + j * ab_dim1], abs(
			    d__1));
#line 559 "sgbsvx.f"
		    anorm = max(d__2,d__3);
#line 560 "sgbsvx.f"
/* L80: */
#line 560 "sgbsvx.f"
		}
#line 561 "sgbsvx.f"
/* L90: */
#line 561 "sgbsvx.f"
	    }
/* Computing MIN */
#line 562 "sgbsvx.f"
	    i__3 = *info - 1, i__2 = *kl + *ku;
#line 562 "sgbsvx.f"
	    i__1 = min(i__3,i__2);
/* Computing MAX */
#line 562 "sgbsvx.f"
	    i__4 = 1, i__5 = *kl + *ku + 2 - *info;
#line 562 "sgbsvx.f"
	    rpvgrw = slantb_("M", "U", "N", info, &i__1, &afb[max(i__4,i__5) 
		    + afb_dim1], ldafb, &work[1], (ftnlen)1, (ftnlen)1, (
		    ftnlen)1);
#line 565 "sgbsvx.f"
	    if (rpvgrw == 0.) {
#line 566 "sgbsvx.f"
		rpvgrw = 1.;
#line 567 "sgbsvx.f"
	    } else {
#line 568 "sgbsvx.f"
		rpvgrw = anorm / rpvgrw;
#line 569 "sgbsvx.f"
	    }
#line 570 "sgbsvx.f"
	    work[1] = rpvgrw;
#line 571 "sgbsvx.f"
	    *rcond = 0.;
#line 572 "sgbsvx.f"
	    return 0;
#line 573 "sgbsvx.f"
	}
#line 574 "sgbsvx.f"
    }

/*     Compute the norm of the matrix A and the */
/*     reciprocal pivot growth factor RPVGRW. */

#line 579 "sgbsvx.f"
    if (notran) {
#line 580 "sgbsvx.f"
	*(unsigned char *)norm = '1';
#line 581 "sgbsvx.f"
    } else {
#line 582 "sgbsvx.f"
	*(unsigned char *)norm = 'I';
#line 583 "sgbsvx.f"
    }
#line 584 "sgbsvx.f"
    anorm = slangb_(norm, n, kl, ku, &ab[ab_offset], ldab, &work[1], (ftnlen)
	    1);
#line 585 "sgbsvx.f"
    i__1 = *kl + *ku;
#line 585 "sgbsvx.f"
    rpvgrw = slantb_("M", "U", "N", n, &i__1, &afb[afb_offset], ldafb, &work[
	    1], (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 586 "sgbsvx.f"
    if (rpvgrw == 0.) {
#line 587 "sgbsvx.f"
	rpvgrw = 1.;
#line 588 "sgbsvx.f"
    } else {
#line 589 "sgbsvx.f"
	rpvgrw = slangb_("M", n, kl, ku, &ab[ab_offset], ldab, &work[1], (
		ftnlen)1) / rpvgrw;
#line 590 "sgbsvx.f"
    }

/*     Compute the reciprocal of the condition number of A. */

#line 594 "sgbsvx.f"
    sgbcon_(norm, n, kl, ku, &afb[afb_offset], ldafb, &ipiv[1], &anorm, rcond,
	     &work[1], &iwork[1], info, (ftnlen)1);

/*     Compute the solution matrix X. */

#line 599 "sgbsvx.f"
    slacpy_("Full", n, nrhs, &b[b_offset], ldb, &x[x_offset], ldx, (ftnlen)4);
#line 600 "sgbsvx.f"
    sgbtrs_(trans, n, kl, ku, nrhs, &afb[afb_offset], ldafb, &ipiv[1], &x[
	    x_offset], ldx, info, (ftnlen)1);

/*     Use iterative refinement to improve the computed solution and */
/*     compute error bounds and backward error estimates for it. */

#line 606 "sgbsvx.f"
    sgbrfs_(trans, n, kl, ku, nrhs, &ab[ab_offset], ldab, &afb[afb_offset], 
	    ldafb, &ipiv[1], &b[b_offset], ldb, &x[x_offset], ldx, &ferr[1], &
	    berr[1], &work[1], &iwork[1], info, (ftnlen)1);

/*     Transform the solution matrix X to a solution of the original */
/*     system. */

#line 612 "sgbsvx.f"
    if (notran) {
#line 613 "sgbsvx.f"
	if (colequ) {
#line 614 "sgbsvx.f"
	    i__1 = *nrhs;
#line 614 "sgbsvx.f"
	    for (j = 1; j <= i__1; ++j) {
#line 615 "sgbsvx.f"
		i__3 = *n;
#line 615 "sgbsvx.f"
		for (i__ = 1; i__ <= i__3; ++i__) {
#line 616 "sgbsvx.f"
		    x[i__ + j * x_dim1] = c__[i__] * x[i__ + j * x_dim1];
#line 617 "sgbsvx.f"
/* L100: */
#line 617 "sgbsvx.f"
		}
#line 618 "sgbsvx.f"
/* L110: */
#line 618 "sgbsvx.f"
	    }
#line 619 "sgbsvx.f"
	    i__1 = *nrhs;
#line 619 "sgbsvx.f"
	    for (j = 1; j <= i__1; ++j) {
#line 620 "sgbsvx.f"
		ferr[j] /= colcnd;
#line 621 "sgbsvx.f"
/* L120: */
#line 621 "sgbsvx.f"
	    }
#line 622 "sgbsvx.f"
	}
#line 623 "sgbsvx.f"
    } else if (rowequ) {
#line 624 "sgbsvx.f"
	i__1 = *nrhs;
#line 624 "sgbsvx.f"
	for (j = 1; j <= i__1; ++j) {
#line 625 "sgbsvx.f"
	    i__3 = *n;
#line 625 "sgbsvx.f"
	    for (i__ = 1; i__ <= i__3; ++i__) {
#line 626 "sgbsvx.f"
		x[i__ + j * x_dim1] = r__[i__] * x[i__ + j * x_dim1];
#line 627 "sgbsvx.f"
/* L130: */
#line 627 "sgbsvx.f"
	    }
#line 628 "sgbsvx.f"
/* L140: */
#line 628 "sgbsvx.f"
	}
#line 629 "sgbsvx.f"
	i__1 = *nrhs;
#line 629 "sgbsvx.f"
	for (j = 1; j <= i__1; ++j) {
#line 630 "sgbsvx.f"
	    ferr[j] /= rowcnd;
#line 631 "sgbsvx.f"
/* L150: */
#line 631 "sgbsvx.f"
	}
#line 632 "sgbsvx.f"
    }

/*     Set INFO = N+1 if the matrix is singular to working precision. */

#line 636 "sgbsvx.f"
    if (*rcond < slamch_("Epsilon", (ftnlen)7)) {
#line 636 "sgbsvx.f"
	*info = *n + 1;
#line 636 "sgbsvx.f"
    }

#line 639 "sgbsvx.f"
    work[1] = rpvgrw;
#line 640 "sgbsvx.f"
    return 0;

/*     End of SGBSVX */

} /* sgbsvx_ */

