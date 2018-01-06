#line 1 "cgbsvx.f"
/* cgbsvx.f -- translated by f2c (version 20100827).
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

#line 1 "cgbsvx.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief <b> CGBSVX computes the solution to system of linear equations A * X = B for GB matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CGBSVX + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgbsvx.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgbsvx.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgbsvx.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CGBSVX( FACT, TRANS, N, KL, KU, NRHS, AB, LDAB, AFB, */
/*                          LDAFB, IPIV, EQUED, R, C, B, LDB, X, LDX, */
/*                          RCOND, FERR, BERR, WORK, RWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          EQUED, FACT, TRANS */
/*       INTEGER            INFO, KL, KU, LDAB, LDAFB, LDB, LDX, N, NRHS */
/*       REAL               RCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       REAL               BERR( * ), C( * ), FERR( * ), R( * ), */
/*      $                   RWORK( * ) */
/*       COMPLEX            AB( LDAB, * ), AFB( LDAFB, * ), B( LDB, * ), */
/*      $                   WORK( * ), X( LDX, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CGBSVX uses the LU factorization to compute the solution to a complex */
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
/* >          AB is COMPLEX array, dimension (LDAB,N) */
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
/* >          AFB is COMPLEX array, dimension (LDAFB,N) */
/* >          If FACT = 'F', then AFB is an input argument and on entry */
/* >          contains details of the LU factorization of the band matrix */
/* >          A, as computed by CGBTRF.  U is stored as an upper triangular */
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
/* >          as computed by CGBTRF; row i of the matrix was interchanged */
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
/* >          B is COMPLEX array, dimension (LDB,NRHS) */
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
/* >          RWORK is REAL array, dimension (N) */
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
/* >                       value of RCOND would suggest. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date April 2012 */

/* > \ingroup complexGBsolve */

/*  ===================================================================== */
/* Subroutine */ int cgbsvx_(char *fact, char *trans, integer *n, integer *kl,
	 integer *ku, integer *nrhs, doublecomplex *ab, integer *ldab, 
	doublecomplex *afb, integer *ldafb, integer *ipiv, char *equed, 
	doublereal *r__, doublereal *c__, doublecomplex *b, integer *ldb, 
	doublecomplex *x, integer *ldx, doublereal *rcond, doublereal *ferr, 
	doublereal *berr, doublecomplex *work, doublereal *rwork, integer *
	info, ftnlen fact_len, ftnlen trans_len, ftnlen equed_len)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, afb_dim1, afb_offset, b_dim1, b_offset, 
	    x_dim1, x_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2;
    doublecomplex z__1;

    /* Builtin functions */
    double z_abs(doublecomplex *);

    /* Local variables */
    static integer i__, j, j1, j2;
    static doublereal amax;
    static char norm[1];
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal rcmin, rcmax, anorm;
    extern /* Subroutine */ int ccopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    static logical equil;
    extern doublereal clangb_(char *, integer *, integer *, integer *, 
	    doublecomplex *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int claqgb_(integer *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, char *, ftnlen), 
	    cgbcon_(char *, integer *, integer *, integer *, doublecomplex *, 
	    integer *, integer *, doublereal *, doublereal *, doublecomplex *,
	     doublereal *, integer *, ftnlen);
    static doublereal colcnd;
    extern doublereal clantb_(char *, char *, char *, integer *, integer *, 
	    doublecomplex *, integer *, doublereal *, ftnlen, ftnlen, ftnlen);
    extern /* Subroutine */ int cgbequ_(integer *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, integer *);
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int cgbrfs_(char *, integer *, integer *, integer 
	    *, integer *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, integer *, doublecomplex *, integer *, doublecomplex *,
	     integer *, doublereal *, doublereal *, doublecomplex *, 
	    doublereal *, integer *, ftnlen), cgbtrf_(integer *, integer *, 
	    integer *, integer *, doublecomplex *, integer *, integer *, 
	    integer *);
    static logical nofact;
    extern /* Subroutine */ int clacpy_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    static doublereal bignum;
    extern /* Subroutine */ int cgbtrs_(char *, integer *, integer *, integer 
	    *, integer *, doublecomplex *, integer *, integer *, 
	    doublecomplex *, integer *, integer *, ftnlen);
    static integer infequ;
    static logical colequ;
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

#line 421 "cgbsvx.f"
    /* Parameter adjustments */
#line 421 "cgbsvx.f"
    ab_dim1 = *ldab;
#line 421 "cgbsvx.f"
    ab_offset = 1 + ab_dim1;
#line 421 "cgbsvx.f"
    ab -= ab_offset;
#line 421 "cgbsvx.f"
    afb_dim1 = *ldafb;
#line 421 "cgbsvx.f"
    afb_offset = 1 + afb_dim1;
#line 421 "cgbsvx.f"
    afb -= afb_offset;
#line 421 "cgbsvx.f"
    --ipiv;
#line 421 "cgbsvx.f"
    --r__;
#line 421 "cgbsvx.f"
    --c__;
#line 421 "cgbsvx.f"
    b_dim1 = *ldb;
#line 421 "cgbsvx.f"
    b_offset = 1 + b_dim1;
#line 421 "cgbsvx.f"
    b -= b_offset;
#line 421 "cgbsvx.f"
    x_dim1 = *ldx;
#line 421 "cgbsvx.f"
    x_offset = 1 + x_dim1;
#line 421 "cgbsvx.f"
    x -= x_offset;
#line 421 "cgbsvx.f"
    --ferr;
#line 421 "cgbsvx.f"
    --berr;
#line 421 "cgbsvx.f"
    --work;
#line 421 "cgbsvx.f"
    --rwork;
#line 421 "cgbsvx.f"

#line 421 "cgbsvx.f"
    /* Function Body */
#line 421 "cgbsvx.f"
    *info = 0;
#line 422 "cgbsvx.f"
    nofact = lsame_(fact, "N", (ftnlen)1, (ftnlen)1);
#line 423 "cgbsvx.f"
    equil = lsame_(fact, "E", (ftnlen)1, (ftnlen)1);
#line 424 "cgbsvx.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 425 "cgbsvx.f"
    if (nofact || equil) {
#line 426 "cgbsvx.f"
	*(unsigned char *)equed = 'N';
#line 427 "cgbsvx.f"
	rowequ = FALSE_;
#line 428 "cgbsvx.f"
	colequ = FALSE_;
#line 429 "cgbsvx.f"
    } else {
#line 430 "cgbsvx.f"
	rowequ = lsame_(equed, "R", (ftnlen)1, (ftnlen)1) || lsame_(equed, 
		"B", (ftnlen)1, (ftnlen)1);
#line 431 "cgbsvx.f"
	colequ = lsame_(equed, "C", (ftnlen)1, (ftnlen)1) || lsame_(equed, 
		"B", (ftnlen)1, (ftnlen)1);
#line 432 "cgbsvx.f"
	smlnum = slamch_("Safe minimum", (ftnlen)12);
#line 433 "cgbsvx.f"
	bignum = 1. / smlnum;
#line 434 "cgbsvx.f"
    }

/*     Test the input parameters. */

#line 438 "cgbsvx.f"
    if (! nofact && ! equil && ! lsame_(fact, "F", (ftnlen)1, (ftnlen)1)) {
#line 440 "cgbsvx.f"
	*info = -1;
#line 441 "cgbsvx.f"
    } else if (! notran && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1) && ! 
	    lsame_(trans, "C", (ftnlen)1, (ftnlen)1)) {
#line 443 "cgbsvx.f"
	*info = -2;
#line 444 "cgbsvx.f"
    } else if (*n < 0) {
#line 445 "cgbsvx.f"
	*info = -3;
#line 446 "cgbsvx.f"
    } else if (*kl < 0) {
#line 447 "cgbsvx.f"
	*info = -4;
#line 448 "cgbsvx.f"
    } else if (*ku < 0) {
#line 449 "cgbsvx.f"
	*info = -5;
#line 450 "cgbsvx.f"
    } else if (*nrhs < 0) {
#line 451 "cgbsvx.f"
	*info = -6;
#line 452 "cgbsvx.f"
    } else if (*ldab < *kl + *ku + 1) {
#line 453 "cgbsvx.f"
	*info = -8;
#line 454 "cgbsvx.f"
    } else if (*ldafb < (*kl << 1) + *ku + 1) {
#line 455 "cgbsvx.f"
	*info = -10;
#line 456 "cgbsvx.f"
    } else if (lsame_(fact, "F", (ftnlen)1, (ftnlen)1) && ! (rowequ || colequ 
	    || lsame_(equed, "N", (ftnlen)1, (ftnlen)1))) {
#line 458 "cgbsvx.f"
	*info = -12;
#line 459 "cgbsvx.f"
    } else {
#line 460 "cgbsvx.f"
	if (rowequ) {
#line 461 "cgbsvx.f"
	    rcmin = bignum;
#line 462 "cgbsvx.f"
	    rcmax = 0.;
#line 463 "cgbsvx.f"
	    i__1 = *n;
#line 463 "cgbsvx.f"
	    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 464 "cgbsvx.f"
		d__1 = rcmin, d__2 = r__[j];
#line 464 "cgbsvx.f"
		rcmin = min(d__1,d__2);
/* Computing MAX */
#line 465 "cgbsvx.f"
		d__1 = rcmax, d__2 = r__[j];
#line 465 "cgbsvx.f"
		rcmax = max(d__1,d__2);
#line 466 "cgbsvx.f"
/* L10: */
#line 466 "cgbsvx.f"
	    }
#line 467 "cgbsvx.f"
	    if (rcmin <= 0.) {
#line 468 "cgbsvx.f"
		*info = -13;
#line 469 "cgbsvx.f"
	    } else if (*n > 0) {
#line 470 "cgbsvx.f"
		rowcnd = max(rcmin,smlnum) / min(rcmax,bignum);
#line 471 "cgbsvx.f"
	    } else {
#line 472 "cgbsvx.f"
		rowcnd = 1.;
#line 473 "cgbsvx.f"
	    }
#line 474 "cgbsvx.f"
	}
#line 475 "cgbsvx.f"
	if (colequ && *info == 0) {
#line 476 "cgbsvx.f"
	    rcmin = bignum;
#line 477 "cgbsvx.f"
	    rcmax = 0.;
#line 478 "cgbsvx.f"
	    i__1 = *n;
#line 478 "cgbsvx.f"
	    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 479 "cgbsvx.f"
		d__1 = rcmin, d__2 = c__[j];
#line 479 "cgbsvx.f"
		rcmin = min(d__1,d__2);
/* Computing MAX */
#line 480 "cgbsvx.f"
		d__1 = rcmax, d__2 = c__[j];
#line 480 "cgbsvx.f"
		rcmax = max(d__1,d__2);
#line 481 "cgbsvx.f"
/* L20: */
#line 481 "cgbsvx.f"
	    }
#line 482 "cgbsvx.f"
	    if (rcmin <= 0.) {
#line 483 "cgbsvx.f"
		*info = -14;
#line 484 "cgbsvx.f"
	    } else if (*n > 0) {
#line 485 "cgbsvx.f"
		colcnd = max(rcmin,smlnum) / min(rcmax,bignum);
#line 486 "cgbsvx.f"
	    } else {
#line 487 "cgbsvx.f"
		colcnd = 1.;
#line 488 "cgbsvx.f"
	    }
#line 489 "cgbsvx.f"
	}
#line 490 "cgbsvx.f"
	if (*info == 0) {
#line 491 "cgbsvx.f"
	    if (*ldb < max(1,*n)) {
#line 492 "cgbsvx.f"
		*info = -16;
#line 493 "cgbsvx.f"
	    } else if (*ldx < max(1,*n)) {
#line 494 "cgbsvx.f"
		*info = -18;
#line 495 "cgbsvx.f"
	    }
#line 496 "cgbsvx.f"
	}
#line 497 "cgbsvx.f"
    }

#line 499 "cgbsvx.f"
    if (*info != 0) {
#line 500 "cgbsvx.f"
	i__1 = -(*info);
#line 500 "cgbsvx.f"
	xerbla_("CGBSVX", &i__1, (ftnlen)6);
#line 501 "cgbsvx.f"
	return 0;
#line 502 "cgbsvx.f"
    }

#line 504 "cgbsvx.f"
    if (equil) {

/*        Compute row and column scalings to equilibrate the matrix A. */

#line 508 "cgbsvx.f"
	cgbequ_(n, n, kl, ku, &ab[ab_offset], ldab, &r__[1], &c__[1], &rowcnd,
		 &colcnd, &amax, &infequ);
#line 510 "cgbsvx.f"
	if (infequ == 0) {

/*           Equilibrate the matrix. */

#line 514 "cgbsvx.f"
	    claqgb_(n, n, kl, ku, &ab[ab_offset], ldab, &r__[1], &c__[1], &
		    rowcnd, &colcnd, &amax, equed, (ftnlen)1);
#line 516 "cgbsvx.f"
	    rowequ = lsame_(equed, "R", (ftnlen)1, (ftnlen)1) || lsame_(equed,
		     "B", (ftnlen)1, (ftnlen)1);
#line 517 "cgbsvx.f"
	    colequ = lsame_(equed, "C", (ftnlen)1, (ftnlen)1) || lsame_(equed,
		     "B", (ftnlen)1, (ftnlen)1);
#line 518 "cgbsvx.f"
	}
#line 519 "cgbsvx.f"
    }

/*     Scale the right hand side. */

#line 523 "cgbsvx.f"
    if (notran) {
#line 524 "cgbsvx.f"
	if (rowequ) {
#line 525 "cgbsvx.f"
	    i__1 = *nrhs;
#line 525 "cgbsvx.f"
	    for (j = 1; j <= i__1; ++j) {
#line 526 "cgbsvx.f"
		i__2 = *n;
#line 526 "cgbsvx.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 527 "cgbsvx.f"
		    i__3 = i__ + j * b_dim1;
#line 527 "cgbsvx.f"
		    i__4 = i__;
#line 527 "cgbsvx.f"
		    i__5 = i__ + j * b_dim1;
#line 527 "cgbsvx.f"
		    z__1.r = r__[i__4] * b[i__5].r, z__1.i = r__[i__4] * b[
			    i__5].i;
#line 527 "cgbsvx.f"
		    b[i__3].r = z__1.r, b[i__3].i = z__1.i;
#line 528 "cgbsvx.f"
/* L30: */
#line 528 "cgbsvx.f"
		}
#line 529 "cgbsvx.f"
/* L40: */
#line 529 "cgbsvx.f"
	    }
#line 530 "cgbsvx.f"
	}
#line 531 "cgbsvx.f"
    } else if (colequ) {
#line 532 "cgbsvx.f"
	i__1 = *nrhs;
#line 532 "cgbsvx.f"
	for (j = 1; j <= i__1; ++j) {
#line 533 "cgbsvx.f"
	    i__2 = *n;
#line 533 "cgbsvx.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 534 "cgbsvx.f"
		i__3 = i__ + j * b_dim1;
#line 534 "cgbsvx.f"
		i__4 = i__;
#line 534 "cgbsvx.f"
		i__5 = i__ + j * b_dim1;
#line 534 "cgbsvx.f"
		z__1.r = c__[i__4] * b[i__5].r, z__1.i = c__[i__4] * b[i__5]
			.i;
#line 534 "cgbsvx.f"
		b[i__3].r = z__1.r, b[i__3].i = z__1.i;
#line 535 "cgbsvx.f"
/* L50: */
#line 535 "cgbsvx.f"
	    }
#line 536 "cgbsvx.f"
/* L60: */
#line 536 "cgbsvx.f"
	}
#line 537 "cgbsvx.f"
    }

#line 539 "cgbsvx.f"
    if (nofact || equil) {

/*        Compute the LU factorization of the band matrix A. */

#line 543 "cgbsvx.f"
	i__1 = *n;
#line 543 "cgbsvx.f"
	for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
#line 544 "cgbsvx.f"
	    i__2 = j - *ku;
#line 544 "cgbsvx.f"
	    j1 = max(i__2,1);
/* Computing MIN */
#line 545 "cgbsvx.f"
	    i__2 = j + *kl;
#line 545 "cgbsvx.f"
	    j2 = min(i__2,*n);
#line 546 "cgbsvx.f"
	    i__2 = j2 - j1 + 1;
#line 546 "cgbsvx.f"
	    ccopy_(&i__2, &ab[*ku + 1 - j + j1 + j * ab_dim1], &c__1, &afb[*
		    kl + *ku + 1 - j + j1 + j * afb_dim1], &c__1);
#line 548 "cgbsvx.f"
/* L70: */
#line 548 "cgbsvx.f"
	}

#line 550 "cgbsvx.f"
	cgbtrf_(n, n, kl, ku, &afb[afb_offset], ldafb, &ipiv[1], info);

/*        Return if INFO is non-zero. */

#line 554 "cgbsvx.f"
	if (*info > 0) {

/*           Compute the reciprocal pivot growth factor of the */
/*           leading rank-deficient INFO columns of A. */

#line 559 "cgbsvx.f"
	    anorm = 0.;
#line 560 "cgbsvx.f"
	    i__1 = *info;
#line 560 "cgbsvx.f"
	    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
#line 561 "cgbsvx.f"
		i__2 = *ku + 2 - j;
/* Computing MIN */
#line 561 "cgbsvx.f"
		i__4 = *n + *ku + 1 - j, i__5 = *kl + *ku + 1;
#line 561 "cgbsvx.f"
		i__3 = min(i__4,i__5);
#line 561 "cgbsvx.f"
		for (i__ = max(i__2,1); i__ <= i__3; ++i__) {
/* Computing MAX */
#line 562 "cgbsvx.f"
		    d__1 = anorm, d__2 = z_abs(&ab[i__ + j * ab_dim1]);
#line 562 "cgbsvx.f"
		    anorm = max(d__1,d__2);
#line 563 "cgbsvx.f"
/* L80: */
#line 563 "cgbsvx.f"
		}
#line 564 "cgbsvx.f"
/* L90: */
#line 564 "cgbsvx.f"
	    }
/* Computing MIN */
#line 565 "cgbsvx.f"
	    i__3 = *info - 1, i__2 = *kl + *ku;
#line 565 "cgbsvx.f"
	    i__1 = min(i__3,i__2);
/* Computing MAX */
#line 565 "cgbsvx.f"
	    i__4 = 1, i__5 = *kl + *ku + 2 - *info;
#line 565 "cgbsvx.f"
	    rpvgrw = clantb_("M", "U", "N", info, &i__1, &afb[max(i__4,i__5) 
		    + afb_dim1], ldafb, &rwork[1], (ftnlen)1, (ftnlen)1, (
		    ftnlen)1);
#line 568 "cgbsvx.f"
	    if (rpvgrw == 0.) {
#line 569 "cgbsvx.f"
		rpvgrw = 1.;
#line 570 "cgbsvx.f"
	    } else {
#line 571 "cgbsvx.f"
		rpvgrw = anorm / rpvgrw;
#line 572 "cgbsvx.f"
	    }
#line 573 "cgbsvx.f"
	    rwork[1] = rpvgrw;
#line 574 "cgbsvx.f"
	    *rcond = 0.;
#line 575 "cgbsvx.f"
	    return 0;
#line 576 "cgbsvx.f"
	}
#line 577 "cgbsvx.f"
    }

/*     Compute the norm of the matrix A and the */
/*     reciprocal pivot growth factor RPVGRW. */

#line 582 "cgbsvx.f"
    if (notran) {
#line 583 "cgbsvx.f"
	*(unsigned char *)norm = '1';
#line 584 "cgbsvx.f"
    } else {
#line 585 "cgbsvx.f"
	*(unsigned char *)norm = 'I';
#line 586 "cgbsvx.f"
    }
#line 587 "cgbsvx.f"
    anorm = clangb_(norm, n, kl, ku, &ab[ab_offset], ldab, &rwork[1], (ftnlen)
	    1);
#line 588 "cgbsvx.f"
    i__1 = *kl + *ku;
#line 588 "cgbsvx.f"
    rpvgrw = clantb_("M", "U", "N", n, &i__1, &afb[afb_offset], ldafb, &rwork[
	    1], (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 589 "cgbsvx.f"
    if (rpvgrw == 0.) {
#line 590 "cgbsvx.f"
	rpvgrw = 1.;
#line 591 "cgbsvx.f"
    } else {
#line 592 "cgbsvx.f"
	rpvgrw = clangb_("M", n, kl, ku, &ab[ab_offset], ldab, &rwork[1], (
		ftnlen)1) / rpvgrw;
#line 593 "cgbsvx.f"
    }

/*     Compute the reciprocal of the condition number of A. */

#line 597 "cgbsvx.f"
    cgbcon_(norm, n, kl, ku, &afb[afb_offset], ldafb, &ipiv[1], &anorm, rcond,
	     &work[1], &rwork[1], info, (ftnlen)1);

/*     Compute the solution matrix X. */

#line 602 "cgbsvx.f"
    clacpy_("Full", n, nrhs, &b[b_offset], ldb, &x[x_offset], ldx, (ftnlen)4);
#line 603 "cgbsvx.f"
    cgbtrs_(trans, n, kl, ku, nrhs, &afb[afb_offset], ldafb, &ipiv[1], &x[
	    x_offset], ldx, info, (ftnlen)1);

/*     Use iterative refinement to improve the computed solution and */
/*     compute error bounds and backward error estimates for it. */

#line 609 "cgbsvx.f"
    cgbrfs_(trans, n, kl, ku, nrhs, &ab[ab_offset], ldab, &afb[afb_offset], 
	    ldafb, &ipiv[1], &b[b_offset], ldb, &x[x_offset], ldx, &ferr[1], &
	    berr[1], &work[1], &rwork[1], info, (ftnlen)1);

/*     Transform the solution matrix X to a solution of the original */
/*     system. */

#line 615 "cgbsvx.f"
    if (notran) {
#line 616 "cgbsvx.f"
	if (colequ) {
#line 617 "cgbsvx.f"
	    i__1 = *nrhs;
#line 617 "cgbsvx.f"
	    for (j = 1; j <= i__1; ++j) {
#line 618 "cgbsvx.f"
		i__3 = *n;
#line 618 "cgbsvx.f"
		for (i__ = 1; i__ <= i__3; ++i__) {
#line 619 "cgbsvx.f"
		    i__2 = i__ + j * x_dim1;
#line 619 "cgbsvx.f"
		    i__4 = i__;
#line 619 "cgbsvx.f"
		    i__5 = i__ + j * x_dim1;
#line 619 "cgbsvx.f"
		    z__1.r = c__[i__4] * x[i__5].r, z__1.i = c__[i__4] * x[
			    i__5].i;
#line 619 "cgbsvx.f"
		    x[i__2].r = z__1.r, x[i__2].i = z__1.i;
#line 620 "cgbsvx.f"
/* L100: */
#line 620 "cgbsvx.f"
		}
#line 621 "cgbsvx.f"
/* L110: */
#line 621 "cgbsvx.f"
	    }
#line 622 "cgbsvx.f"
	    i__1 = *nrhs;
#line 622 "cgbsvx.f"
	    for (j = 1; j <= i__1; ++j) {
#line 623 "cgbsvx.f"
		ferr[j] /= colcnd;
#line 624 "cgbsvx.f"
/* L120: */
#line 624 "cgbsvx.f"
	    }
#line 625 "cgbsvx.f"
	}
#line 626 "cgbsvx.f"
    } else if (rowequ) {
#line 627 "cgbsvx.f"
	i__1 = *nrhs;
#line 627 "cgbsvx.f"
	for (j = 1; j <= i__1; ++j) {
#line 628 "cgbsvx.f"
	    i__3 = *n;
#line 628 "cgbsvx.f"
	    for (i__ = 1; i__ <= i__3; ++i__) {
#line 629 "cgbsvx.f"
		i__2 = i__ + j * x_dim1;
#line 629 "cgbsvx.f"
		i__4 = i__;
#line 629 "cgbsvx.f"
		i__5 = i__ + j * x_dim1;
#line 629 "cgbsvx.f"
		z__1.r = r__[i__4] * x[i__5].r, z__1.i = r__[i__4] * x[i__5]
			.i;
#line 629 "cgbsvx.f"
		x[i__2].r = z__1.r, x[i__2].i = z__1.i;
#line 630 "cgbsvx.f"
/* L130: */
#line 630 "cgbsvx.f"
	    }
#line 631 "cgbsvx.f"
/* L140: */
#line 631 "cgbsvx.f"
	}
#line 632 "cgbsvx.f"
	i__1 = *nrhs;
#line 632 "cgbsvx.f"
	for (j = 1; j <= i__1; ++j) {
#line 633 "cgbsvx.f"
	    ferr[j] /= rowcnd;
#line 634 "cgbsvx.f"
/* L150: */
#line 634 "cgbsvx.f"
	}
#line 635 "cgbsvx.f"
    }

/*     Set INFO = N+1 if the matrix is singular to working precision. */

#line 639 "cgbsvx.f"
    if (*rcond < slamch_("Epsilon", (ftnlen)7)) {
#line 639 "cgbsvx.f"
	*info = *n + 1;
#line 639 "cgbsvx.f"
    }

#line 642 "cgbsvx.f"
    rwork[1] = rpvgrw;
#line 643 "cgbsvx.f"
    return 0;

/*     End of CGBSVX */

} /* cgbsvx_ */

