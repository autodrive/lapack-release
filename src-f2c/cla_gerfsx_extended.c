#line 1 "cla_gerfsx_extended.f"
/* cla_gerfsx_extended.f -- translated by f2c (version 20100827).
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

#line 1 "cla_gerfsx_extended.f"
/* Table of constant values */

static integer c__1 = 1;
static doublecomplex c_b6 = {-1.,0.};
static doublecomplex c_b8 = {1.,0.};
static doublereal c_b31 = 1.;

/* > \brief \b CLA_GERFSX_EXTENDED */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLA_GERFSX_EXTENDED + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cla_ger
fsx_extended.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cla_ger
fsx_extended.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cla_ger
fsx_extended.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CLA_GERFSX_EXTENDED( PREC_TYPE, TRANS_TYPE, N, NRHS, A, */
/*                                       LDA, AF, LDAF, IPIV, COLEQU, C, B, */
/*                                       LDB, Y, LDY, BERR_OUT, N_NORMS, */
/*                                       ERRS_N, ERRS_C, RES, AYB, DY, */
/*                                       Y_TAIL, RCOND, ITHRESH, RTHRESH, */
/*                                       DZ_UB, IGNORE_CWISE, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, LDAF, LDB, LDY, N, NRHS, PREC_TYPE, */
/*      $                   TRANS_TYPE, N_NORMS */
/*       LOGICAL            COLEQU, IGNORE_CWISE */
/*       INTEGER            ITHRESH */
/*       REAL               RTHRESH, DZ_UB */
/*       .. */
/*       .. Array Arguments */
/*       INTEGER            IPIV( * ) */
/*       COMPLEX            A( LDA, * ), AF( LDAF, * ), B( LDB, * ), */
/*      $                   Y( LDY, * ), RES( * ), DY( * ), Y_TAIL( * ) */
/*       REAL               C( * ), AYB( * ), RCOND, BERR_OUT( * ), */
/*      $                   ERRS_N( NRHS, * ), ERRS_C( NRHS, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > */
/* > CLA_GERFSX_EXTENDED improves the computed solution to a system of */
/* > linear equations by performing extra-precise iterative refinement */
/* > and provides error bounds and backward error estimates for the solution. */
/* > This subroutine is called by CGERFSX to perform iterative refinement. */
/* > In addition to normwise error bound, the code provides maximum */
/* > componentwise error bound if possible. See comments for ERRS_N */
/* > and ERRS_C for details of the error bounds. Note that this */
/* > subroutine is only resonsible for setting the second fields of */
/* > ERRS_N and ERRS_C. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] PREC_TYPE */
/* > \verbatim */
/* >          PREC_TYPE is INTEGER */
/* >     Specifies the intermediate precision to be used in refinement. */
/* >     The value is defined by ILAPREC(P) where P is a CHARACTER and */
/* >     P    = 'S':  Single */
/* >          = 'D':  Double */
/* >          = 'I':  Indigenous */
/* >          = 'X', 'E':  Extra */
/* > \endverbatim */
/* > */
/* > \param[in] TRANS_TYPE */
/* > \verbatim */
/* >          TRANS_TYPE is INTEGER */
/* >     Specifies the transposition operation on A. */
/* >     The value is defined by ILATRANS(T) where T is a CHARACTER and */
/* >     T    = 'N':  No transpose */
/* >          = 'T':  Transpose */
/* >          = 'C':  Conjugate transpose */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >     The number of linear equations, i.e., the order of the */
/* >     matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NRHS */
/* > \verbatim */
/* >          NRHS is INTEGER */
/* >     The number of right-hand-sides, i.e., the number of columns of the */
/* >     matrix B. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is COMPLEX array, dimension (LDA,N) */
/* >     On entry, the N-by-N matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >     The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] AF */
/* > \verbatim */
/* >          AF is COMPLEX array, dimension (LDAF,N) */
/* >     The factors L and U from the factorization */
/* >     A = P*L*U as computed by CGETRF. */
/* > \endverbatim */
/* > */
/* > \param[in] LDAF */
/* > \verbatim */
/* >          LDAF is INTEGER */
/* >     The leading dimension of the array AF.  LDAF >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N) */
/* >     The pivot indices from the factorization A = P*L*U */
/* >     as computed by CGETRF; row i of the matrix was interchanged */
/* >     with row IPIV(i). */
/* > \endverbatim */
/* > */
/* > \param[in] COLEQU */
/* > \verbatim */
/* >          COLEQU is LOGICAL */
/* >     If .TRUE. then column equilibration was done to A before calling */
/* >     this routine. This is needed to compute the solution and error */
/* >     bounds correctly. */
/* > \endverbatim */
/* > */
/* > \param[in] C */
/* > \verbatim */
/* >          C is REAL array, dimension (N) */
/* >     The column scale factors for A. If COLEQU = .FALSE., C */
/* >     is not accessed. If C is input, each element of C should be a power */
/* >     of the radix to ensure a reliable solution and error estimates. */
/* >     Scaling by powers of the radix does not cause rounding errors unless */
/* >     the result underflows or overflows. Rounding errors during scaling */
/* >     lead to refining with a matrix that is not equivalent to the */
/* >     input matrix, producing error estimates that may not be */
/* >     reliable. */
/* > \endverbatim */
/* > */
/* > \param[in] B */
/* > \verbatim */
/* >          B is COMPLEX array, dimension (LDB,NRHS) */
/* >     The right-hand-side matrix B. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >     The leading dimension of the array B.  LDB >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] Y */
/* > \verbatim */
/* >          Y is COMPLEX array, dimension (LDY,NRHS) */
/* >     On entry, the solution matrix X, as computed by CGETRS. */
/* >     On exit, the improved solution matrix Y. */
/* > \endverbatim */
/* > */
/* > \param[in] LDY */
/* > \verbatim */
/* >          LDY is INTEGER */
/* >     The leading dimension of the array Y.  LDY >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] BERR_OUT */
/* > \verbatim */
/* >          BERR_OUT is REAL array, dimension (NRHS) */
/* >     On exit, BERR_OUT(j) contains the componentwise relative backward */
/* >     error for right-hand-side j from the formula */
/* >         max(i) ( abs(RES(i)) / ( abs(op(A_s))*abs(Y) + abs(B_s) )(i) ) */
/* >     where abs(Z) is the componentwise absolute value of the matrix */
/* >     or vector Z. This is computed by CLA_LIN_BERR. */
/* > \endverbatim */
/* > */
/* > \param[in] N_NORMS */
/* > \verbatim */
/* >          N_NORMS is INTEGER */
/* >     Determines which error bounds to return (see ERRS_N */
/* >     and ERRS_C). */
/* >     If N_NORMS >= 1 return normwise error bounds. */
/* >     If N_NORMS >= 2 return componentwise error bounds. */
/* > \endverbatim */
/* > */
/* > \param[in,out] ERRS_N */
/* > \verbatim */
/* >          ERRS_N is REAL array, dimension (NRHS, N_ERR_BNDS) */
/* >     For each right-hand side, this array contains information about */
/* >     various error bounds and condition numbers corresponding to the */
/* >     normwise relative error, which is defined as follows: */
/* > */
/* >     Normwise relative error in the ith solution vector: */
/* >             max_j (abs(XTRUE(j,i) - X(j,i))) */
/* >            ------------------------------ */
/* >                  max_j abs(X(j,i)) */
/* > */
/* >     The array is indexed by the type of error information as described */
/* >     below. There currently are up to three pieces of information */
/* >     returned. */
/* > */
/* >     The first index in ERRS_N(i,:) corresponds to the ith */
/* >     right-hand side. */
/* > */
/* >     The second index in ERRS_N(:,err) contains the following */
/* >     three fields: */
/* >     err = 1 "Trust/don't trust" boolean. Trust the answer if the */
/* >              reciprocal condition number is less than the threshold */
/* >              sqrt(n) * slamch('Epsilon'). */
/* > */
/* >     err = 2 "Guaranteed" error bound: The estimated forward error, */
/* >              almost certainly within a factor of 10 of the true error */
/* >              so long as the next entry is greater than the threshold */
/* >              sqrt(n) * slamch('Epsilon'). This error bound should only */
/* >              be trusted if the previous boolean is true. */
/* > */
/* >     err = 3  Reciprocal condition number: Estimated normwise */
/* >              reciprocal condition number.  Compared with the threshold */
/* >              sqrt(n) * slamch('Epsilon') to determine if the error */
/* >              estimate is "guaranteed". These reciprocal condition */
/* >              numbers are 1 / (norm(Z^{-1},inf) * norm(Z,inf)) for some */
/* >              appropriately scaled matrix Z. */
/* >              Let Z = S*A, where S scales each row by a power of the */
/* >              radix so all absolute row sums of Z are approximately 1. */
/* > */
/* >     This subroutine is only responsible for setting the second field */
/* >     above. */
/* >     See Lapack Working Note 165 for further details and extra */
/* >     cautions. */
/* > \endverbatim */
/* > */
/* > \param[in,out] ERRS_C */
/* > \verbatim */
/* >          ERRS_C is REAL array, dimension (NRHS, N_ERR_BNDS) */
/* >     For each right-hand side, this array contains information about */
/* >     various error bounds and condition numbers corresponding to the */
/* >     componentwise relative error, which is defined as follows: */
/* > */
/* >     Componentwise relative error in the ith solution vector: */
/* >                    abs(XTRUE(j,i) - X(j,i)) */
/* >             max_j ---------------------- */
/* >                         abs(X(j,i)) */
/* > */
/* >     The array is indexed by the right-hand side i (on which the */
/* >     componentwise relative error depends), and the type of error */
/* >     information as described below. There currently are up to three */
/* >     pieces of information returned for each right-hand side. If */
/* >     componentwise accuracy is not requested (PARAMS(3) = 0.0), then */
/* >     ERRS_C is not accessed.  If N_ERR_BNDS .LT. 3, then at most */
/* >     the first (:,N_ERR_BNDS) entries are returned. */
/* > */
/* >     The first index in ERRS_C(i,:) corresponds to the ith */
/* >     right-hand side. */
/* > */
/* >     The second index in ERRS_C(:,err) contains the following */
/* >     three fields: */
/* >     err = 1 "Trust/don't trust" boolean. Trust the answer if the */
/* >              reciprocal condition number is less than the threshold */
/* >              sqrt(n) * slamch('Epsilon'). */
/* > */
/* >     err = 2 "Guaranteed" error bound: The estimated forward error, */
/* >              almost certainly within a factor of 10 of the true error */
/* >              so long as the next entry is greater than the threshold */
/* >              sqrt(n) * slamch('Epsilon'). This error bound should only */
/* >              be trusted if the previous boolean is true. */
/* > */
/* >     err = 3  Reciprocal condition number: Estimated componentwise */
/* >              reciprocal condition number.  Compared with the threshold */
/* >              sqrt(n) * slamch('Epsilon') to determine if the error */
/* >              estimate is "guaranteed". These reciprocal condition */
/* >              numbers are 1 / (norm(Z^{-1},inf) * norm(Z,inf)) for some */
/* >              appropriately scaled matrix Z. */
/* >              Let Z = S*(A*diag(x)), where x is the solution for the */
/* >              current right-hand side and S scales each row of */
/* >              A*diag(x) by a power of the radix so all absolute row */
/* >              sums of Z are approximately 1. */
/* > */
/* >     This subroutine is only responsible for setting the second field */
/* >     above. */
/* >     See Lapack Working Note 165 for further details and extra */
/* >     cautions. */
/* > \endverbatim */
/* > */
/* > \param[in] RES */
/* > \verbatim */
/* >          RES is COMPLEX array, dimension (N) */
/* >     Workspace to hold the intermediate residual. */
/* > \endverbatim */
/* > */
/* > \param[in] AYB */
/* > \verbatim */
/* >          AYB is REAL array, dimension (N) */
/* >     Workspace. */
/* > \endverbatim */
/* > */
/* > \param[in] DY */
/* > \verbatim */
/* >          DY is COMPLEX array, dimension (N) */
/* >     Workspace to hold the intermediate solution. */
/* > \endverbatim */
/* > */
/* > \param[in] Y_TAIL */
/* > \verbatim */
/* >          Y_TAIL is COMPLEX array, dimension (N) */
/* >     Workspace to hold the trailing bits of the intermediate solution. */
/* > \endverbatim */
/* > */
/* > \param[in] RCOND */
/* > \verbatim */
/* >          RCOND is REAL */
/* >     Reciprocal scaled condition number.  This is an estimate of the */
/* >     reciprocal Skeel condition number of the matrix A after */
/* >     equilibration (if done).  If this is less than the machine */
/* >     precision (in particular, if it is zero), the matrix is singular */
/* >     to working precision.  Note that the error may still be small even */
/* >     if this number is very small and the matrix appears ill- */
/* >     conditioned. */
/* > \endverbatim */
/* > */
/* > \param[in] ITHRESH */
/* > \verbatim */
/* >          ITHRESH is INTEGER */
/* >     The maximum number of residual computations allowed for */
/* >     refinement. The default is 10. For 'aggressive' set to 100 to */
/* >     permit convergence using approximate factorizations or */
/* >     factorizations other than LU. If the factorization uses a */
/* >     technique other than Gaussian elimination, the guarantees in */
/* >     ERRS_N and ERRS_C may no longer be trustworthy. */
/* > \endverbatim */
/* > */
/* > \param[in] RTHRESH */
/* > \verbatim */
/* >          RTHRESH is REAL */
/* >     Determines when to stop refinement if the error estimate stops */
/* >     decreasing. Refinement will stop when the next solution no longer */
/* >     satisfies norm(dx_{i+1}) < RTHRESH * norm(dx_i) where norm(Z) is */
/* >     the infinity norm of Z. RTHRESH satisfies 0 < RTHRESH <= 1. The */
/* >     default value is 0.5. For 'aggressive' set to 0.9 to permit */
/* >     convergence on extremely ill-conditioned matrices. See LAWN 165 */
/* >     for more details. */
/* > \endverbatim */
/* > */
/* > \param[in] DZ_UB */
/* > \verbatim */
/* >          DZ_UB is REAL */
/* >     Determines when to start considering componentwise convergence. */
/* >     Componentwise convergence is only considered after each component */
/* >     of the solution Y is stable, which we definte as the relative */
/* >     change in each component being less than DZ_UB. The default value */
/* >     is 0.25, requiring the first bit to be stable. See LAWN 165 for */
/* >     more details. */
/* > \endverbatim */
/* > */
/* > \param[in] IGNORE_CWISE */
/* > \verbatim */
/* >          IGNORE_CWISE is LOGICAL */
/* >     If .TRUE. then ignore componentwise convergence. Default value */
/* >     is .FALSE.. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >       = 0:  Successful exit. */
/* >       < 0:  if INFO = -i, the ith argument to CGETRS had an illegal */
/* >             value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup complexGEcomputational */

/*  ===================================================================== */
/* Subroutine */ int cla_gerfsx_extended__(integer *prec_type__, integer *
	trans_type__, integer *n, integer *nrhs, doublecomplex *a, integer *
	lda, doublecomplex *af, integer *ldaf, integer *ipiv, logical *colequ,
	 doublereal *c__, doublecomplex *b, integer *ldb, doublecomplex *y, 
	integer *ldy, doublereal *berr_out__, integer *n_norms__, doublereal *
	errs_n__, doublereal *errs_c__, doublecomplex *res, doublereal *ayb, 
	doublecomplex *dy, doublecomplex *y_tail__, doublereal *rcond, 
	integer *ithresh, doublereal *rthresh, doublereal *dz_ub__, logical *
	ignore_cwise__, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, af_dim1, af_offset, b_dim1, b_offset, y_dim1, 
	    y_offset, errs_n_dim1, errs_n_offset, errs_c_dim1, errs_c_offset, 
	    i__1, i__2, i__3, i__4;
    doublereal d__1, d__2;
    char ch__1[1];

    /* Builtin functions */
    double d_imag(doublecomplex *);

    /* Local variables */
    static doublereal dxratmax, dzratmax;
    static integer i__, j;
    extern /* Subroutine */ int cla_geamv__(integer *, integer *, integer *, 
	    doublereal *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublereal *, doublereal *, integer *);
    static logical incr_prec__;
    static doublereal prev_dz_z__, yk, final_dx_x__;
    extern /* Subroutine */ int cla_wwaddw__(integer *, doublecomplex *, 
	    doublecomplex *, doublecomplex *);
    static doublereal final_dz_z__, prevnormdx;
    static integer cnt;
    static doublereal dyk, eps, incr_thresh__, dx_x__, dz_z__;
    extern /* Subroutine */ int cla_lin_berr__(integer *, integer *, integer *
	    , doublecomplex *, doublereal *, doublereal *);
    static doublereal ymin;
    extern /* Subroutine */ int blas_cgemv_x__(integer *, integer *, integer *
	    , doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, integer *)
	    ;
    static integer y_prec_state__;
    extern /* Subroutine */ int blas_cgemv2_x__(integer *, integer *, integer 
	    *, doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, integer *), cgemv_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, ftnlen), 
	    ccopy_(integer *, doublecomplex *, integer *, doublecomplex *, 
	    integer *);
    static doublereal dxrat, dzrat;
    extern /* Subroutine */ int caxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    static char trans[1];
    static doublereal normx, normy;
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int cgetrs_(char *, integer *, integer *, 
	    doublecomplex *, integer *, integer *, doublecomplex *, integer *,
	     integer *, ftnlen);
    static doublereal normdx;
    extern /* Character */ VOID chla_transtype__(char *, ftnlen, integer *);
    static doublereal hugeval;
    static integer x_state__, z_state__;


/*  -- LAPACK computational routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments */
/*     .. */

/*  ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. Parameters .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Statement Functions .. */
/*     .. */
/*     .. Statement Function Definitions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 477 "cla_gerfsx_extended.f"
    /* Parameter adjustments */
#line 477 "cla_gerfsx_extended.f"
    errs_c_dim1 = *nrhs;
#line 477 "cla_gerfsx_extended.f"
    errs_c_offset = 1 + errs_c_dim1;
#line 477 "cla_gerfsx_extended.f"
    errs_c__ -= errs_c_offset;
#line 477 "cla_gerfsx_extended.f"
    errs_n_dim1 = *nrhs;
#line 477 "cla_gerfsx_extended.f"
    errs_n_offset = 1 + errs_n_dim1;
#line 477 "cla_gerfsx_extended.f"
    errs_n__ -= errs_n_offset;
#line 477 "cla_gerfsx_extended.f"
    a_dim1 = *lda;
#line 477 "cla_gerfsx_extended.f"
    a_offset = 1 + a_dim1;
#line 477 "cla_gerfsx_extended.f"
    a -= a_offset;
#line 477 "cla_gerfsx_extended.f"
    af_dim1 = *ldaf;
#line 477 "cla_gerfsx_extended.f"
    af_offset = 1 + af_dim1;
#line 477 "cla_gerfsx_extended.f"
    af -= af_offset;
#line 477 "cla_gerfsx_extended.f"
    --ipiv;
#line 477 "cla_gerfsx_extended.f"
    --c__;
#line 477 "cla_gerfsx_extended.f"
    b_dim1 = *ldb;
#line 477 "cla_gerfsx_extended.f"
    b_offset = 1 + b_dim1;
#line 477 "cla_gerfsx_extended.f"
    b -= b_offset;
#line 477 "cla_gerfsx_extended.f"
    y_dim1 = *ldy;
#line 477 "cla_gerfsx_extended.f"
    y_offset = 1 + y_dim1;
#line 477 "cla_gerfsx_extended.f"
    y -= y_offset;
#line 477 "cla_gerfsx_extended.f"
    --berr_out__;
#line 477 "cla_gerfsx_extended.f"
    --res;
#line 477 "cla_gerfsx_extended.f"
    --ayb;
#line 477 "cla_gerfsx_extended.f"
    --dy;
#line 477 "cla_gerfsx_extended.f"
    --y_tail__;
#line 477 "cla_gerfsx_extended.f"

#line 477 "cla_gerfsx_extended.f"
    /* Function Body */
#line 477 "cla_gerfsx_extended.f"
    if (*info != 0) {
#line 477 "cla_gerfsx_extended.f"
	return 0;
#line 477 "cla_gerfsx_extended.f"
    }
#line 478 "cla_gerfsx_extended.f"
    chla_transtype__(ch__1, (ftnlen)1, trans_type__);
#line 478 "cla_gerfsx_extended.f"
    *(unsigned char *)trans = *(unsigned char *)&ch__1[0];
#line 479 "cla_gerfsx_extended.f"
    eps = slamch_("Epsilon", (ftnlen)7);
#line 480 "cla_gerfsx_extended.f"
    hugeval = slamch_("Overflow", (ftnlen)8);
/*     Force HUGEVAL to Inf */
#line 482 "cla_gerfsx_extended.f"
    hugeval *= hugeval;
/*     Using HUGEVAL may lead to spurious underflows. */
#line 484 "cla_gerfsx_extended.f"
    incr_thresh__ = (doublereal) (*n) * eps;

#line 486 "cla_gerfsx_extended.f"
    i__1 = *nrhs;
#line 486 "cla_gerfsx_extended.f"
    for (j = 1; j <= i__1; ++j) {
#line 487 "cla_gerfsx_extended.f"
	y_prec_state__ = 1;
#line 488 "cla_gerfsx_extended.f"
	if (y_prec_state__ == 2) {
#line 489 "cla_gerfsx_extended.f"
	    i__2 = *n;
#line 489 "cla_gerfsx_extended.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 490 "cla_gerfsx_extended.f"
		i__3 = i__;
#line 490 "cla_gerfsx_extended.f"
		y_tail__[i__3].r = 0., y_tail__[i__3].i = 0.;
#line 491 "cla_gerfsx_extended.f"
	    }
#line 492 "cla_gerfsx_extended.f"
	}
#line 494 "cla_gerfsx_extended.f"
	dxrat = 0.;
#line 495 "cla_gerfsx_extended.f"
	dxratmax = 0.;
#line 496 "cla_gerfsx_extended.f"
	dzrat = 0.;
#line 497 "cla_gerfsx_extended.f"
	dzratmax = 0.;
#line 498 "cla_gerfsx_extended.f"
	final_dx_x__ = hugeval;
#line 499 "cla_gerfsx_extended.f"
	final_dz_z__ = hugeval;
#line 500 "cla_gerfsx_extended.f"
	prevnormdx = hugeval;
#line 501 "cla_gerfsx_extended.f"
	prev_dz_z__ = hugeval;
#line 502 "cla_gerfsx_extended.f"
	dz_z__ = hugeval;
#line 503 "cla_gerfsx_extended.f"
	dx_x__ = hugeval;
#line 505 "cla_gerfsx_extended.f"
	x_state__ = 1;
#line 506 "cla_gerfsx_extended.f"
	z_state__ = 0;
#line 507 "cla_gerfsx_extended.f"
	incr_prec__ = FALSE_;
#line 509 "cla_gerfsx_extended.f"
	i__2 = *ithresh;
#line 509 "cla_gerfsx_extended.f"
	for (cnt = 1; cnt <= i__2; ++cnt) {

/*         Compute residual RES = B_s - op(A_s) * Y, */
/*             op(A) = A, A**T, or A**H depending on TRANS (and type). */

#line 514 "cla_gerfsx_extended.f"
	    ccopy_(n, &b[j * b_dim1 + 1], &c__1, &res[1], &c__1);
#line 515 "cla_gerfsx_extended.f"
	    if (y_prec_state__ == 0) {
#line 516 "cla_gerfsx_extended.f"
		cgemv_(trans, n, n, &c_b6, &a[a_offset], lda, &y[j * y_dim1 + 
			1], &c__1, &c_b8, &res[1], &c__1, (ftnlen)1);
#line 518 "cla_gerfsx_extended.f"
	    } else if (y_prec_state__ == 1) {
#line 519 "cla_gerfsx_extended.f"
		blas_cgemv_x__(trans_type__, n, n, &c_b6, &a[a_offset], lda, &
			y[j * y_dim1 + 1], &c__1, &c_b8, &res[1], &c__1, 
			prec_type__);
#line 522 "cla_gerfsx_extended.f"
	    } else {
#line 523 "cla_gerfsx_extended.f"
		blas_cgemv2_x__(trans_type__, n, n, &c_b6, &a[a_offset], lda, 
			&y[j * y_dim1 + 1], &y_tail__[1], &c__1, &c_b8, &res[
			1], &c__1, prec_type__);
#line 526 "cla_gerfsx_extended.f"
	    }
/*         XXX: RES is no longer needed. */
#line 529 "cla_gerfsx_extended.f"
	    ccopy_(n, &res[1], &c__1, &dy[1], &c__1);
#line 530 "cla_gerfsx_extended.f"
	    cgetrs_(trans, n, &c__1, &af[af_offset], ldaf, &ipiv[1], &dy[1], 
		    n, info, (ftnlen)1);

/*         Calculate relative changes DX_X, DZ_Z and ratios DXRAT, DZRAT. */

#line 534 "cla_gerfsx_extended.f"
	    normx = 0.;
#line 535 "cla_gerfsx_extended.f"
	    normy = 0.;
#line 536 "cla_gerfsx_extended.f"
	    normdx = 0.;
#line 537 "cla_gerfsx_extended.f"
	    dz_z__ = 0.;
#line 538 "cla_gerfsx_extended.f"
	    ymin = hugeval;

#line 540 "cla_gerfsx_extended.f"
	    i__3 = *n;
#line 540 "cla_gerfsx_extended.f"
	    for (i__ = 1; i__ <= i__3; ++i__) {
#line 541 "cla_gerfsx_extended.f"
		i__4 = i__ + j * y_dim1;
#line 541 "cla_gerfsx_extended.f"
		yk = (d__1 = y[i__4].r, abs(d__1)) + (d__2 = d_imag(&y[i__ + 
			j * y_dim1]), abs(d__2));
#line 542 "cla_gerfsx_extended.f"
		i__4 = i__;
#line 542 "cla_gerfsx_extended.f"
		dyk = (d__1 = dy[i__4].r, abs(d__1)) + (d__2 = d_imag(&dy[i__]
			), abs(d__2));
#line 544 "cla_gerfsx_extended.f"
		if (yk != 0.) {
/* Computing MAX */
#line 545 "cla_gerfsx_extended.f"
		    d__1 = dz_z__, d__2 = dyk / yk;
#line 545 "cla_gerfsx_extended.f"
		    dz_z__ = max(d__1,d__2);
#line 546 "cla_gerfsx_extended.f"
		} else if (dyk != 0.) {
#line 547 "cla_gerfsx_extended.f"
		    dz_z__ = hugeval;
#line 548 "cla_gerfsx_extended.f"
		}
#line 550 "cla_gerfsx_extended.f"
		ymin = min(ymin,yk);
#line 552 "cla_gerfsx_extended.f"
		normy = max(normy,yk);
#line 554 "cla_gerfsx_extended.f"
		if (*colequ) {
/* Computing MAX */
#line 555 "cla_gerfsx_extended.f"
		    d__1 = normx, d__2 = yk * c__[i__];
#line 555 "cla_gerfsx_extended.f"
		    normx = max(d__1,d__2);
/* Computing MAX */
#line 556 "cla_gerfsx_extended.f"
		    d__1 = normdx, d__2 = dyk * c__[i__];
#line 556 "cla_gerfsx_extended.f"
		    normdx = max(d__1,d__2);
#line 557 "cla_gerfsx_extended.f"
		} else {
#line 558 "cla_gerfsx_extended.f"
		    normx = normy;
#line 559 "cla_gerfsx_extended.f"
		    normdx = max(normdx,dyk);
#line 560 "cla_gerfsx_extended.f"
		}
#line 561 "cla_gerfsx_extended.f"
	    }
#line 563 "cla_gerfsx_extended.f"
	    if (normx != 0.) {
#line 564 "cla_gerfsx_extended.f"
		dx_x__ = normdx / normx;
#line 565 "cla_gerfsx_extended.f"
	    } else if (normdx == 0.) {
#line 566 "cla_gerfsx_extended.f"
		dx_x__ = 0.;
#line 567 "cla_gerfsx_extended.f"
	    } else {
#line 568 "cla_gerfsx_extended.f"
		dx_x__ = hugeval;
#line 569 "cla_gerfsx_extended.f"
	    }
#line 571 "cla_gerfsx_extended.f"
	    dxrat = normdx / prevnormdx;
#line 572 "cla_gerfsx_extended.f"
	    dzrat = dz_z__ / prev_dz_z__;

/*         Check termination criteria */

#line 576 "cla_gerfsx_extended.f"
	    if (! (*ignore_cwise__) && ymin * *rcond < incr_thresh__ * normy 
		    && y_prec_state__ < 2) {
#line 576 "cla_gerfsx_extended.f"
		incr_prec__ = TRUE_;
#line 576 "cla_gerfsx_extended.f"
	    }
#line 581 "cla_gerfsx_extended.f"
	    if (x_state__ == 3 && dxrat <= *rthresh) {
#line 581 "cla_gerfsx_extended.f"
		x_state__ = 1;
#line 581 "cla_gerfsx_extended.f"
	    }
#line 583 "cla_gerfsx_extended.f"
	    if (x_state__ == 1) {
#line 584 "cla_gerfsx_extended.f"
		if (dx_x__ <= eps) {
#line 585 "cla_gerfsx_extended.f"
		    x_state__ = 2;
#line 586 "cla_gerfsx_extended.f"
		} else if (dxrat > *rthresh) {
#line 587 "cla_gerfsx_extended.f"
		    if (y_prec_state__ != 2) {
#line 588 "cla_gerfsx_extended.f"
			incr_prec__ = TRUE_;
#line 589 "cla_gerfsx_extended.f"
		    } else {
#line 590 "cla_gerfsx_extended.f"
			x_state__ = 3;
#line 591 "cla_gerfsx_extended.f"
		    }
#line 592 "cla_gerfsx_extended.f"
		} else {
#line 593 "cla_gerfsx_extended.f"
		    if (dxrat > dxratmax) {
#line 593 "cla_gerfsx_extended.f"
			dxratmax = dxrat;
#line 593 "cla_gerfsx_extended.f"
		    }
#line 594 "cla_gerfsx_extended.f"
		}
#line 595 "cla_gerfsx_extended.f"
		if (x_state__ > 1) {
#line 595 "cla_gerfsx_extended.f"
		    final_dx_x__ = dx_x__;
#line 595 "cla_gerfsx_extended.f"
		}
#line 596 "cla_gerfsx_extended.f"
	    }
#line 598 "cla_gerfsx_extended.f"
	    if (z_state__ == 0 && dz_z__ <= *dz_ub__) {
#line 598 "cla_gerfsx_extended.f"
		z_state__ = 1;
#line 598 "cla_gerfsx_extended.f"
	    }
#line 600 "cla_gerfsx_extended.f"
	    if (z_state__ == 3 && dzrat <= *rthresh) {
#line 600 "cla_gerfsx_extended.f"
		z_state__ = 1;
#line 600 "cla_gerfsx_extended.f"
	    }
#line 602 "cla_gerfsx_extended.f"
	    if (z_state__ == 1) {
#line 603 "cla_gerfsx_extended.f"
		if (dz_z__ <= eps) {
#line 604 "cla_gerfsx_extended.f"
		    z_state__ = 2;
#line 605 "cla_gerfsx_extended.f"
		} else if (dz_z__ > *dz_ub__) {
#line 606 "cla_gerfsx_extended.f"
		    z_state__ = 0;
#line 607 "cla_gerfsx_extended.f"
		    dzratmax = 0.;
#line 608 "cla_gerfsx_extended.f"
		    final_dz_z__ = hugeval;
#line 609 "cla_gerfsx_extended.f"
		} else if (dzrat > *rthresh) {
#line 610 "cla_gerfsx_extended.f"
		    if (y_prec_state__ != 2) {
#line 611 "cla_gerfsx_extended.f"
			incr_prec__ = TRUE_;
#line 612 "cla_gerfsx_extended.f"
		    } else {
#line 613 "cla_gerfsx_extended.f"
			z_state__ = 3;
#line 614 "cla_gerfsx_extended.f"
		    }
#line 615 "cla_gerfsx_extended.f"
		} else {
#line 616 "cla_gerfsx_extended.f"
		    if (dzrat > dzratmax) {
#line 616 "cla_gerfsx_extended.f"
			dzratmax = dzrat;
#line 616 "cla_gerfsx_extended.f"
		    }
#line 617 "cla_gerfsx_extended.f"
		}
#line 618 "cla_gerfsx_extended.f"
		if (z_state__ > 1) {
#line 618 "cla_gerfsx_extended.f"
		    final_dz_z__ = dz_z__;
#line 618 "cla_gerfsx_extended.f"
		}
#line 619 "cla_gerfsx_extended.f"
	    }

/*           Exit if both normwise and componentwise stopped working, */
/*           but if componentwise is unstable, let it go at least two */
/*           iterations. */

#line 625 "cla_gerfsx_extended.f"
	    if (x_state__ != 1) {
#line 626 "cla_gerfsx_extended.f"
		if (*ignore_cwise__) {
#line 626 "cla_gerfsx_extended.f"
		    goto L666;
#line 626 "cla_gerfsx_extended.f"
		}
#line 627 "cla_gerfsx_extended.f"
		if (z_state__ == 3 || z_state__ == 2) {
#line 627 "cla_gerfsx_extended.f"
		    goto L666;
#line 627 "cla_gerfsx_extended.f"
		}
#line 629 "cla_gerfsx_extended.f"
		if (z_state__ == 0 && cnt > 1) {
#line 629 "cla_gerfsx_extended.f"
		    goto L666;
#line 629 "cla_gerfsx_extended.f"
		}
#line 630 "cla_gerfsx_extended.f"
	    }
#line 632 "cla_gerfsx_extended.f"
	    if (incr_prec__) {
#line 633 "cla_gerfsx_extended.f"
		incr_prec__ = FALSE_;
#line 634 "cla_gerfsx_extended.f"
		++y_prec_state__;
#line 635 "cla_gerfsx_extended.f"
		i__3 = *n;
#line 635 "cla_gerfsx_extended.f"
		for (i__ = 1; i__ <= i__3; ++i__) {
#line 636 "cla_gerfsx_extended.f"
		    i__4 = i__;
#line 636 "cla_gerfsx_extended.f"
		    y_tail__[i__4].r = 0., y_tail__[i__4].i = 0.;
#line 637 "cla_gerfsx_extended.f"
		}
#line 638 "cla_gerfsx_extended.f"
	    }
#line 640 "cla_gerfsx_extended.f"
	    prevnormdx = normdx;
#line 641 "cla_gerfsx_extended.f"
	    prev_dz_z__ = dz_z__;

/*           Update soluton. */

#line 645 "cla_gerfsx_extended.f"
	    if (y_prec_state__ < 2) {
#line 646 "cla_gerfsx_extended.f"
		caxpy_(n, &c_b8, &dy[1], &c__1, &y[j * y_dim1 + 1], &c__1);
#line 647 "cla_gerfsx_extended.f"
	    } else {
#line 648 "cla_gerfsx_extended.f"
		cla_wwaddw__(n, &y[j * y_dim1 + 1], &y_tail__[1], &dy[1]);
#line 649 "cla_gerfsx_extended.f"
	    }
#line 651 "cla_gerfsx_extended.f"
	}
/*        Target of "IF (Z_STOP .AND. X_STOP)".  Sun's f77 won't EXIT. */
#line 653 "cla_gerfsx_extended.f"
L666:

/*     Set final_* when cnt hits ithresh */

#line 657 "cla_gerfsx_extended.f"
	if (x_state__ == 1) {
#line 657 "cla_gerfsx_extended.f"
	    final_dx_x__ = dx_x__;
#line 657 "cla_gerfsx_extended.f"
	}
#line 658 "cla_gerfsx_extended.f"
	if (z_state__ == 1) {
#line 658 "cla_gerfsx_extended.f"
	    final_dz_z__ = dz_z__;
#line 658 "cla_gerfsx_extended.f"
	}

/*     Compute error bounds */

#line 662 "cla_gerfsx_extended.f"
	if (*n_norms__ >= 1) {
#line 663 "cla_gerfsx_extended.f"
	    errs_n__[j + (errs_n_dim1 << 1)] = final_dx_x__ / (1 - dxratmax);
#line 665 "cla_gerfsx_extended.f"
	}
#line 666 "cla_gerfsx_extended.f"
	if (*n_norms__ >= 2) {
#line 667 "cla_gerfsx_extended.f"
	    errs_c__[j + (errs_c_dim1 << 1)] = final_dz_z__ / (1 - dzratmax);
#line 668 "cla_gerfsx_extended.f"
	}

/*     Compute componentwise relative backward error from formula */
/*         max(i) ( abs(R(i)) / ( abs(op(A_s))*abs(Y) + abs(B_s) )(i) ) */
/*     where abs(Z) is the componentwise absolute value of the matrix */
/*     or vector Z. */

/*        Compute residual RES = B_s - op(A_s) * Y, */
/*            op(A) = A, A**T, or A**H depending on TRANS (and type). */

#line 678 "cla_gerfsx_extended.f"
	ccopy_(n, &b[j * b_dim1 + 1], &c__1, &res[1], &c__1);
#line 679 "cla_gerfsx_extended.f"
	cgemv_(trans, n, n, &c_b6, &a[a_offset], lda, &y[j * y_dim1 + 1], &
		c__1, &c_b8, &res[1], &c__1, (ftnlen)1);
#line 682 "cla_gerfsx_extended.f"
	i__2 = *n;
#line 682 "cla_gerfsx_extended.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 683 "cla_gerfsx_extended.f"
	    i__3 = i__ + j * b_dim1;
#line 683 "cla_gerfsx_extended.f"
	    ayb[i__] = (d__1 = b[i__3].r, abs(d__1)) + (d__2 = d_imag(&b[i__ 
		    + j * b_dim1]), abs(d__2));
#line 684 "cla_gerfsx_extended.f"
	}

/*     Compute abs(op(A_s))*abs(Y) + abs(B_s). */

#line 688 "cla_gerfsx_extended.f"
	cla_geamv__(trans_type__, n, n, &c_b31, &a[a_offset], lda, &y[j * 
		y_dim1 + 1], &c__1, &c_b31, &ayb[1], &c__1);
#line 691 "cla_gerfsx_extended.f"
	cla_lin_berr__(n, n, &c__1, &res[1], &ayb[1], &berr_out__[j]);

/*     End of loop for each RHS. */

#line 695 "cla_gerfsx_extended.f"
    }

#line 697 "cla_gerfsx_extended.f"
    return 0;
} /* cla_gerfsx_extended__ */

