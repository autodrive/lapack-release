#line 1 "sla_gerfsx_extended.f"
/* sla_gerfsx_extended.f -- translated by f2c (version 20100827).
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

#line 1 "sla_gerfsx_extended.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b6 = -1.;
static doublereal c_b8 = 1.;

/* > \brief \b SLA_GERFSX_EXTENDED improves the computed solution to a system of linear equations for general 
matrices by performing extra-precise iterative refinement and provides error bounds and backward error
 estimates for the solution. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLA_GERFSX_EXTENDED + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sla_ger
fsx_extended.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sla_ger
fsx_extended.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sla_ger
fsx_extended.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLA_GERFSX_EXTENDED( PREC_TYPE, TRANS_TYPE, N, NRHS, A, */
/*                                       LDA, AF, LDAF, IPIV, COLEQU, C, B, */
/*                                       LDB, Y, LDY, BERR_OUT, N_NORMS, */
/*                                       ERRS_N, ERRS_C, RES, */
/*                                       AYB, DY, Y_TAIL, RCOND, ITHRESH, */
/*                                       RTHRESH, DZ_UB, IGNORE_CWISE, */
/*                                       INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, LDAF, LDB, LDY, N, NRHS, PREC_TYPE, */
/*      $                   TRANS_TYPE, N_NORMS, ITHRESH */
/*       LOGICAL            COLEQU, IGNORE_CWISE */
/*       REAL               RTHRESH, DZ_UB */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       REAL               A( LDA, * ), AF( LDAF, * ), B( LDB, * ), */
/*      $                   Y( LDY, * ), RES( * ), DY( * ), Y_TAIL( * ) */
/*       REAL               C( * ), AYB( * ), RCOND, BERR_OUT( * ), */
/*      $                   ERRS_N( NRHS, * ), */
/*      $                   ERRS_C( NRHS, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLA_GERFSX_EXTENDED improves the computed solution to a system of */
/* > linear equations by performing extra-precise iterative refinement */
/* > and provides error bounds and backward error estimates for the solution. */
/* > This subroutine is called by SGERFSX to perform iterative refinement. */
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
/* >          A is REAL array, dimension (LDA,N) */
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
/* >          AF is REAL array, dimension (LDAF,N) */
/* >     The factors L and U from the factorization */
/* >     A = P*L*U as computed by SGETRF. */
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
/* >     as computed by SGETRF; row i of the matrix was interchanged */
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
/* >          B is REAL array, dimension (LDB,NRHS) */
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
/* >          Y is REAL array, dimension (LDY,NRHS) */
/* >     On entry, the solution matrix X, as computed by SGETRS. */
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
/* >     or vector Z. This is computed by SLA_LIN_BERR. */
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
/* >          RES is REAL array, dimension (N) */
/* >     Workspace to hold the intermediate residual. */
/* > \endverbatim */
/* > */
/* > \param[in] AYB */
/* > \verbatim */
/* >          AYB is REAL array, dimension (N) */
/* >     Workspace. This can be the same workspace passed for Y_TAIL. */
/* > \endverbatim */
/* > */
/* > \param[in] DY */
/* > \verbatim */
/* >          DY is REAL array, dimension (N) */
/* >     Workspace to hold the intermediate solution. */
/* > \endverbatim */
/* > */
/* > \param[in] Y_TAIL */
/* > \verbatim */
/* >          Y_TAIL is REAL array, dimension (N) */
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
/* >       < 0:  if INFO = -i, the ith argument to SGETRS had an illegal */
/* >             value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup realGEcomputational */

/*  ===================================================================== */
/* Subroutine */ int sla_gerfsx_extended__(integer *prec_type__, integer *
	trans_type__, integer *n, integer *nrhs, doublereal *a, integer *lda, 
	doublereal *af, integer *ldaf, integer *ipiv, logical *colequ, 
	doublereal *c__, doublereal *b, integer *ldb, doublereal *y, integer *
	ldy, doublereal *berr_out__, integer *n_norms__, doublereal *errs_n__,
	 doublereal *errs_c__, doublereal *res, doublereal *ayb, doublereal *
	dy, doublereal *y_tail__, doublereal *rcond, integer *ithresh, 
	doublereal *rthresh, doublereal *dz_ub__, logical *ignore_cwise__, 
	integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, af_dim1, af_offset, b_dim1, b_offset, y_dim1, 
	    y_offset, errs_n_dim1, errs_n_offset, errs_c_dim1, errs_c_offset, 
	    i__1, i__2, i__3;
    doublereal d__1, d__2;
    char ch__1[1];

    /* Local variables */
    static doublereal dxratmax, dzratmax;
    static integer i__, j;
    extern /* Subroutine */ int sla_geamv__(integer *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *);
    static logical incr_prec__;
    static doublereal prev_dz_z__, yk, final_dx_x__, final_dz_z__;
    extern /* Subroutine */ int sla_wwaddw__(integer *, doublereal *, 
	    doublereal *, doublereal *);
    static doublereal prevnormdx;
    static integer cnt;
    static doublereal dyk, eps, incr_thresh__, dx_x__, dz_z__, ymin;
    extern /* Subroutine */ int sla_lin_berr__(integer *, integer *, integer *
	    , doublereal *, doublereal *, doublereal *), blas_sgemv_x__(
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, integer *);
    static integer y_prec_state__;
    extern /* Subroutine */ int blas_sgemv2_x__(integer *, integer *, integer 
	    *, doublereal *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *), sgemv_(char *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, ftnlen);
    static doublereal dxrat, dzrat;
    static char trans[1];
    extern /* Subroutine */ int scopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static doublereal normx, normy;
    extern /* Subroutine */ int saxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    extern doublereal slamch_(char *, ftnlen);
    static doublereal normdx;
    extern /* Subroutine */ int sgetrs_(char *, integer *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, ftnlen);
    extern /* Character */ VOID chla_transtype__(char *, ftnlen, integer *);
    static doublereal hugeval;
    static integer x_state__, z_state__;


/*  -- LAPACK computational routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
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
/*     .. Executable Statements .. */

#line 470 "sla_gerfsx_extended.f"
    /* Parameter adjustments */
#line 470 "sla_gerfsx_extended.f"
    errs_c_dim1 = *nrhs;
#line 470 "sla_gerfsx_extended.f"
    errs_c_offset = 1 + errs_c_dim1;
#line 470 "sla_gerfsx_extended.f"
    errs_c__ -= errs_c_offset;
#line 470 "sla_gerfsx_extended.f"
    errs_n_dim1 = *nrhs;
#line 470 "sla_gerfsx_extended.f"
    errs_n_offset = 1 + errs_n_dim1;
#line 470 "sla_gerfsx_extended.f"
    errs_n__ -= errs_n_offset;
#line 470 "sla_gerfsx_extended.f"
    a_dim1 = *lda;
#line 470 "sla_gerfsx_extended.f"
    a_offset = 1 + a_dim1;
#line 470 "sla_gerfsx_extended.f"
    a -= a_offset;
#line 470 "sla_gerfsx_extended.f"
    af_dim1 = *ldaf;
#line 470 "sla_gerfsx_extended.f"
    af_offset = 1 + af_dim1;
#line 470 "sla_gerfsx_extended.f"
    af -= af_offset;
#line 470 "sla_gerfsx_extended.f"
    --ipiv;
#line 470 "sla_gerfsx_extended.f"
    --c__;
#line 470 "sla_gerfsx_extended.f"
    b_dim1 = *ldb;
#line 470 "sla_gerfsx_extended.f"
    b_offset = 1 + b_dim1;
#line 470 "sla_gerfsx_extended.f"
    b -= b_offset;
#line 470 "sla_gerfsx_extended.f"
    y_dim1 = *ldy;
#line 470 "sla_gerfsx_extended.f"
    y_offset = 1 + y_dim1;
#line 470 "sla_gerfsx_extended.f"
    y -= y_offset;
#line 470 "sla_gerfsx_extended.f"
    --berr_out__;
#line 470 "sla_gerfsx_extended.f"
    --res;
#line 470 "sla_gerfsx_extended.f"
    --ayb;
#line 470 "sla_gerfsx_extended.f"
    --dy;
#line 470 "sla_gerfsx_extended.f"
    --y_tail__;
#line 470 "sla_gerfsx_extended.f"

#line 470 "sla_gerfsx_extended.f"
    /* Function Body */
#line 470 "sla_gerfsx_extended.f"
    if (*info != 0) {
#line 470 "sla_gerfsx_extended.f"
	return 0;
#line 470 "sla_gerfsx_extended.f"
    }
#line 471 "sla_gerfsx_extended.f"
    chla_transtype__(ch__1, (ftnlen)1, trans_type__);
#line 471 "sla_gerfsx_extended.f"
    *(unsigned char *)trans = *(unsigned char *)&ch__1[0];
#line 472 "sla_gerfsx_extended.f"
    eps = slamch_("Epsilon", (ftnlen)7);
#line 473 "sla_gerfsx_extended.f"
    hugeval = slamch_("Overflow", (ftnlen)8);
/*     Force HUGEVAL to Inf */
#line 475 "sla_gerfsx_extended.f"
    hugeval *= hugeval;
/*     Using HUGEVAL may lead to spurious underflows. */
#line 477 "sla_gerfsx_extended.f"
    incr_thresh__ = (doublereal) (*n) * eps;

#line 479 "sla_gerfsx_extended.f"
    i__1 = *nrhs;
#line 479 "sla_gerfsx_extended.f"
    for (j = 1; j <= i__1; ++j) {
#line 480 "sla_gerfsx_extended.f"
	y_prec_state__ = 1;
#line 481 "sla_gerfsx_extended.f"
	if (y_prec_state__ == 2) {
#line 482 "sla_gerfsx_extended.f"
	    i__2 = *n;
#line 482 "sla_gerfsx_extended.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 483 "sla_gerfsx_extended.f"
		y_tail__[i__] = 0.;
#line 484 "sla_gerfsx_extended.f"
	    }
#line 485 "sla_gerfsx_extended.f"
	}
#line 487 "sla_gerfsx_extended.f"
	dxrat = 0.;
#line 488 "sla_gerfsx_extended.f"
	dxratmax = 0.;
#line 489 "sla_gerfsx_extended.f"
	dzrat = 0.;
#line 490 "sla_gerfsx_extended.f"
	dzratmax = 0.;
#line 491 "sla_gerfsx_extended.f"
	final_dx_x__ = hugeval;
#line 492 "sla_gerfsx_extended.f"
	final_dz_z__ = hugeval;
#line 493 "sla_gerfsx_extended.f"
	prevnormdx = hugeval;
#line 494 "sla_gerfsx_extended.f"
	prev_dz_z__ = hugeval;
#line 495 "sla_gerfsx_extended.f"
	dz_z__ = hugeval;
#line 496 "sla_gerfsx_extended.f"
	dx_x__ = hugeval;
#line 498 "sla_gerfsx_extended.f"
	x_state__ = 1;
#line 499 "sla_gerfsx_extended.f"
	z_state__ = 0;
#line 500 "sla_gerfsx_extended.f"
	incr_prec__ = FALSE_;
#line 502 "sla_gerfsx_extended.f"
	i__2 = *ithresh;
#line 502 "sla_gerfsx_extended.f"
	for (cnt = 1; cnt <= i__2; ++cnt) {

/*         Compute residual RES = B_s - op(A_s) * Y, */
/*             op(A) = A, A**T, or A**H depending on TRANS (and type). */

#line 507 "sla_gerfsx_extended.f"
	    scopy_(n, &b[j * b_dim1 + 1], &c__1, &res[1], &c__1);
#line 508 "sla_gerfsx_extended.f"
	    if (y_prec_state__ == 0) {
#line 509 "sla_gerfsx_extended.f"
		sgemv_(trans, n, n, &c_b6, &a[a_offset], lda, &y[j * y_dim1 + 
			1], &c__1, &c_b8, &res[1], &c__1, (ftnlen)1);
#line 511 "sla_gerfsx_extended.f"
	    } else if (y_prec_state__ == 1) {
#line 512 "sla_gerfsx_extended.f"
		blas_sgemv_x__(trans_type__, n, n, &c_b6, &a[a_offset], lda, &
			y[j * y_dim1 + 1], &c__1, &c_b8, &res[1], &c__1, 
			prec_type__);
#line 514 "sla_gerfsx_extended.f"
	    } else {
#line 515 "sla_gerfsx_extended.f"
		blas_sgemv2_x__(trans_type__, n, n, &c_b6, &a[a_offset], lda, 
			&y[j * y_dim1 + 1], &y_tail__[1], &c__1, &c_b8, &res[
			1], &c__1, prec_type__);
#line 517 "sla_gerfsx_extended.f"
	    }
/*        XXX: RES is no longer needed. */
#line 520 "sla_gerfsx_extended.f"
	    scopy_(n, &res[1], &c__1, &dy[1], &c__1);
#line 521 "sla_gerfsx_extended.f"
	    sgetrs_(trans, n, &c__1, &af[af_offset], ldaf, &ipiv[1], &dy[1], 
		    n, info, (ftnlen)1);

/*         Calculate relative changes DX_X, DZ_Z and ratios DXRAT, DZRAT. */

#line 525 "sla_gerfsx_extended.f"
	    normx = 0.;
#line 526 "sla_gerfsx_extended.f"
	    normy = 0.;
#line 527 "sla_gerfsx_extended.f"
	    normdx = 0.;
#line 528 "sla_gerfsx_extended.f"
	    dz_z__ = 0.;
#line 529 "sla_gerfsx_extended.f"
	    ymin = hugeval;

#line 531 "sla_gerfsx_extended.f"
	    i__3 = *n;
#line 531 "sla_gerfsx_extended.f"
	    for (i__ = 1; i__ <= i__3; ++i__) {
#line 532 "sla_gerfsx_extended.f"
		yk = (d__1 = y[i__ + j * y_dim1], abs(d__1));
#line 533 "sla_gerfsx_extended.f"
		dyk = (d__1 = dy[i__], abs(d__1));
#line 535 "sla_gerfsx_extended.f"
		if (yk != 0.) {
/* Computing MAX */
#line 536 "sla_gerfsx_extended.f"
		    d__1 = dz_z__, d__2 = dyk / yk;
#line 536 "sla_gerfsx_extended.f"
		    dz_z__ = max(d__1,d__2);
#line 537 "sla_gerfsx_extended.f"
		} else if (dyk != 0.) {
#line 538 "sla_gerfsx_extended.f"
		    dz_z__ = hugeval;
#line 539 "sla_gerfsx_extended.f"
		}
#line 541 "sla_gerfsx_extended.f"
		ymin = min(ymin,yk);
#line 543 "sla_gerfsx_extended.f"
		normy = max(normy,yk);
#line 545 "sla_gerfsx_extended.f"
		if (*colequ) {
/* Computing MAX */
#line 546 "sla_gerfsx_extended.f"
		    d__1 = normx, d__2 = yk * c__[i__];
#line 546 "sla_gerfsx_extended.f"
		    normx = max(d__1,d__2);
/* Computing MAX */
#line 547 "sla_gerfsx_extended.f"
		    d__1 = normdx, d__2 = dyk * c__[i__];
#line 547 "sla_gerfsx_extended.f"
		    normdx = max(d__1,d__2);
#line 548 "sla_gerfsx_extended.f"
		} else {
#line 549 "sla_gerfsx_extended.f"
		    normx = normy;
#line 550 "sla_gerfsx_extended.f"
		    normdx = max(normdx,dyk);
#line 551 "sla_gerfsx_extended.f"
		}
#line 552 "sla_gerfsx_extended.f"
	    }
#line 554 "sla_gerfsx_extended.f"
	    if (normx != 0.) {
#line 555 "sla_gerfsx_extended.f"
		dx_x__ = normdx / normx;
#line 556 "sla_gerfsx_extended.f"
	    } else if (normdx == 0.) {
#line 557 "sla_gerfsx_extended.f"
		dx_x__ = 0.;
#line 558 "sla_gerfsx_extended.f"
	    } else {
#line 559 "sla_gerfsx_extended.f"
		dx_x__ = hugeval;
#line 560 "sla_gerfsx_extended.f"
	    }
#line 562 "sla_gerfsx_extended.f"
	    dxrat = normdx / prevnormdx;
#line 563 "sla_gerfsx_extended.f"
	    dzrat = dz_z__ / prev_dz_z__;

/*         Check termination criteria */

#line 567 "sla_gerfsx_extended.f"
	    if (! (*ignore_cwise__) && ymin * *rcond < incr_thresh__ * normy 
		    && y_prec_state__ < 2) {
#line 567 "sla_gerfsx_extended.f"
		incr_prec__ = TRUE_;
#line 567 "sla_gerfsx_extended.f"
	    }
#line 572 "sla_gerfsx_extended.f"
	    if (x_state__ == 3 && dxrat <= *rthresh) {
#line 572 "sla_gerfsx_extended.f"
		x_state__ = 1;
#line 572 "sla_gerfsx_extended.f"
	    }
#line 574 "sla_gerfsx_extended.f"
	    if (x_state__ == 1) {
#line 575 "sla_gerfsx_extended.f"
		if (dx_x__ <= eps) {
#line 576 "sla_gerfsx_extended.f"
		    x_state__ = 2;
#line 577 "sla_gerfsx_extended.f"
		} else if (dxrat > *rthresh) {
#line 578 "sla_gerfsx_extended.f"
		    if (y_prec_state__ != 2) {
#line 579 "sla_gerfsx_extended.f"
			incr_prec__ = TRUE_;
#line 580 "sla_gerfsx_extended.f"
		    } else {
#line 581 "sla_gerfsx_extended.f"
			x_state__ = 3;
#line 582 "sla_gerfsx_extended.f"
		    }
#line 583 "sla_gerfsx_extended.f"
		} else {
#line 584 "sla_gerfsx_extended.f"
		    if (dxrat > dxratmax) {
#line 584 "sla_gerfsx_extended.f"
			dxratmax = dxrat;
#line 584 "sla_gerfsx_extended.f"
		    }
#line 585 "sla_gerfsx_extended.f"
		}
#line 586 "sla_gerfsx_extended.f"
		if (x_state__ > 1) {
#line 586 "sla_gerfsx_extended.f"
		    final_dx_x__ = dx_x__;
#line 586 "sla_gerfsx_extended.f"
		}
#line 587 "sla_gerfsx_extended.f"
	    }
#line 589 "sla_gerfsx_extended.f"
	    if (z_state__ == 0 && dz_z__ <= *dz_ub__) {
#line 589 "sla_gerfsx_extended.f"
		z_state__ = 1;
#line 589 "sla_gerfsx_extended.f"
	    }
#line 591 "sla_gerfsx_extended.f"
	    if (z_state__ == 3 && dzrat <= *rthresh) {
#line 591 "sla_gerfsx_extended.f"
		z_state__ = 1;
#line 591 "sla_gerfsx_extended.f"
	    }
#line 593 "sla_gerfsx_extended.f"
	    if (z_state__ == 1) {
#line 594 "sla_gerfsx_extended.f"
		if (dz_z__ <= eps) {
#line 595 "sla_gerfsx_extended.f"
		    z_state__ = 2;
#line 596 "sla_gerfsx_extended.f"
		} else if (dz_z__ > *dz_ub__) {
#line 597 "sla_gerfsx_extended.f"
		    z_state__ = 0;
#line 598 "sla_gerfsx_extended.f"
		    dzratmax = 0.;
#line 599 "sla_gerfsx_extended.f"
		    final_dz_z__ = hugeval;
#line 600 "sla_gerfsx_extended.f"
		} else if (dzrat > *rthresh) {
#line 601 "sla_gerfsx_extended.f"
		    if (y_prec_state__ != 2) {
#line 602 "sla_gerfsx_extended.f"
			incr_prec__ = TRUE_;
#line 603 "sla_gerfsx_extended.f"
		    } else {
#line 604 "sla_gerfsx_extended.f"
			z_state__ = 3;
#line 605 "sla_gerfsx_extended.f"
		    }
#line 606 "sla_gerfsx_extended.f"
		} else {
#line 607 "sla_gerfsx_extended.f"
		    if (dzrat > dzratmax) {
#line 607 "sla_gerfsx_extended.f"
			dzratmax = dzrat;
#line 607 "sla_gerfsx_extended.f"
		    }
#line 608 "sla_gerfsx_extended.f"
		}
#line 609 "sla_gerfsx_extended.f"
		if (z_state__ > 1) {
#line 609 "sla_gerfsx_extended.f"
		    final_dz_z__ = dz_z__;
#line 609 "sla_gerfsx_extended.f"
		}
#line 610 "sla_gerfsx_extended.f"
	    }

/*           Exit if both normwise and componentwise stopped working, */
/*           but if componentwise is unstable, let it go at least two */
/*           iterations. */

#line 616 "sla_gerfsx_extended.f"
	    if (x_state__ != 1) {
#line 617 "sla_gerfsx_extended.f"
		if (*ignore_cwise__) {
#line 617 "sla_gerfsx_extended.f"
		    goto L666;
#line 617 "sla_gerfsx_extended.f"
		}
#line 618 "sla_gerfsx_extended.f"
		if (z_state__ == 3 || z_state__ == 2) {
#line 618 "sla_gerfsx_extended.f"
		    goto L666;
#line 618 "sla_gerfsx_extended.f"
		}
#line 620 "sla_gerfsx_extended.f"
		if (z_state__ == 0 && cnt > 1) {
#line 620 "sla_gerfsx_extended.f"
		    goto L666;
#line 620 "sla_gerfsx_extended.f"
		}
#line 621 "sla_gerfsx_extended.f"
	    }
#line 623 "sla_gerfsx_extended.f"
	    if (incr_prec__) {
#line 624 "sla_gerfsx_extended.f"
		incr_prec__ = FALSE_;
#line 625 "sla_gerfsx_extended.f"
		++y_prec_state__;
#line 626 "sla_gerfsx_extended.f"
		i__3 = *n;
#line 626 "sla_gerfsx_extended.f"
		for (i__ = 1; i__ <= i__3; ++i__) {
#line 627 "sla_gerfsx_extended.f"
		    y_tail__[i__] = 0.;
#line 628 "sla_gerfsx_extended.f"
		}
#line 629 "sla_gerfsx_extended.f"
	    }
#line 631 "sla_gerfsx_extended.f"
	    prevnormdx = normdx;
#line 632 "sla_gerfsx_extended.f"
	    prev_dz_z__ = dz_z__;

/*           Update soluton. */

#line 636 "sla_gerfsx_extended.f"
	    if (y_prec_state__ < 2) {
#line 637 "sla_gerfsx_extended.f"
		saxpy_(n, &c_b8, &dy[1], &c__1, &y[j * y_dim1 + 1], &c__1);
#line 638 "sla_gerfsx_extended.f"
	    } else {
#line 639 "sla_gerfsx_extended.f"
		sla_wwaddw__(n, &y[j * y_dim1 + 1], &y_tail__[1], &dy[1]);
#line 640 "sla_gerfsx_extended.f"
	    }
#line 642 "sla_gerfsx_extended.f"
	}
/*        Target of "IF (Z_STOP .AND. X_STOP)".  Sun's f77 won't EXIT. */
#line 644 "sla_gerfsx_extended.f"
L666:

/*     Set final_* when cnt hits ithresh. */

#line 648 "sla_gerfsx_extended.f"
	if (x_state__ == 1) {
#line 648 "sla_gerfsx_extended.f"
	    final_dx_x__ = dx_x__;
#line 648 "sla_gerfsx_extended.f"
	}
#line 649 "sla_gerfsx_extended.f"
	if (z_state__ == 1) {
#line 649 "sla_gerfsx_extended.f"
	    final_dz_z__ = dz_z__;
#line 649 "sla_gerfsx_extended.f"
	}

/*     Compute error bounds */

#line 653 "sla_gerfsx_extended.f"
	if (*n_norms__ >= 1) {
#line 654 "sla_gerfsx_extended.f"
	    errs_n__[j + (errs_n_dim1 << 1)] = final_dx_x__ / (1 - dxratmax);
#line 656 "sla_gerfsx_extended.f"
	}
#line 657 "sla_gerfsx_extended.f"
	if (*n_norms__ >= 2) {
#line 658 "sla_gerfsx_extended.f"
	    errs_c__[j + (errs_c_dim1 << 1)] = final_dz_z__ / (1 - dzratmax);
#line 660 "sla_gerfsx_extended.f"
	}

/*     Compute componentwise relative backward error from formula */
/*         max(i) ( abs(R(i)) / ( abs(op(A_s))*abs(Y) + abs(B_s) )(i) ) */
/*     where abs(Z) is the componentwise absolute value of the matrix */
/*     or vector Z. */

/*         Compute residual RES = B_s - op(A_s) * Y, */
/*             op(A) = A, A**T, or A**H depending on TRANS (and type). */

#line 670 "sla_gerfsx_extended.f"
	scopy_(n, &b[j * b_dim1 + 1], &c__1, &res[1], &c__1);
#line 671 "sla_gerfsx_extended.f"
	sgemv_(trans, n, n, &c_b6, &a[a_offset], lda, &y[j * y_dim1 + 1], &
		c__1, &c_b8, &res[1], &c__1, (ftnlen)1);
#line 673 "sla_gerfsx_extended.f"
	i__2 = *n;
#line 673 "sla_gerfsx_extended.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 674 "sla_gerfsx_extended.f"
	    ayb[i__] = (d__1 = b[i__ + j * b_dim1], abs(d__1));
#line 675 "sla_gerfsx_extended.f"
	}

/*     Compute abs(op(A_s))*abs(Y) + abs(B_s). */

#line 679 "sla_gerfsx_extended.f"
	sla_geamv__(trans_type__, n, n, &c_b8, &a[a_offset], lda, &y[j * 
		y_dim1 + 1], &c__1, &c_b8, &ayb[1], &c__1);
#line 682 "sla_gerfsx_extended.f"
	sla_lin_berr__(n, n, &c__1, &res[1], &ayb[1], &berr_out__[j]);

/*     End of loop for each RHS. */

#line 686 "sla_gerfsx_extended.f"
    }

#line 688 "sla_gerfsx_extended.f"
    return 0;
} /* sla_gerfsx_extended__ */

