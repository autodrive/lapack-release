#line 1 "sla_gbrfsx_extended.f"
/* sla_gbrfsx_extended.f -- translated by f2c (version 20100827).
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

#line 1 "sla_gbrfsx_extended.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b6 = -1.;
static doublereal c_b8 = 1.;

/* > \brief \b SLA_GBRFSX_EXTENDED improves the computed solution to a system of linear equations for general 
banded matrices by performing extra-precise iterative refinement and provides error bounds and backwar
d error estimates for the solution. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLA_GBRFSX_EXTENDED + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sla_gbr
fsx_extended.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sla_gbr
fsx_extended.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sla_gbr
fsx_extended.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLA_GBRFSX_EXTENDED( PREC_TYPE, TRANS_TYPE, N, KL, KU, */
/*                                       NRHS, AB, LDAB, AFB, LDAFB, IPIV, */
/*                                       COLEQU, C, B, LDB, Y, LDY, */
/*                                       BERR_OUT, N_NORMS, ERR_BNDS_NORM, */
/*                                       ERR_BNDS_COMP, RES, AYB, DY, */
/*                                       Y_TAIL, RCOND, ITHRESH, RTHRESH, */
/*                                       DZ_UB, IGNORE_CWISE, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDAB, LDAFB, LDB, LDY, N, KL, KU, NRHS, */
/*      $                   PREC_TYPE, TRANS_TYPE, N_NORMS, ITHRESH */
/*       LOGICAL            COLEQU, IGNORE_CWISE */
/*       REAL               RTHRESH, DZ_UB */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       REAL               AB( LDAB, * ), AFB( LDAFB, * ), B( LDB, * ), */
/*      $                   Y( LDY, * ), RES(*), DY(*), Y_TAIL(*) */
/*       REAL               C( * ), AYB(*), RCOND, BERR_OUT(*), */
/*      $                   ERR_BNDS_NORM( NRHS, * ), */
/*      $                   ERR_BNDS_COMP( NRHS, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLA_GBRFSX_EXTENDED improves the computed solution to a system of */
/* > linear equations by performing extra-precise iterative refinement */
/* > and provides error bounds and backward error estimates for the solution. */
/* > This subroutine is called by SGBRFSX to perform iterative refinement. */
/* > In addition to normwise error bound, the code provides maximum */
/* > componentwise error bound if possible. See comments for ERR_BNDS_NORM */
/* > and ERR_BNDS_COMP for details of the error bounds. Note that this */
/* > subroutine is only resonsible for setting the second fields of */
/* > ERR_BNDS_NORM and ERR_BNDS_COMP. */
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
/* > \param[in] KL */
/* > \verbatim */
/* >          KL is INTEGER */
/* >     The number of subdiagonals within the band of A.  KL >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] KU */
/* > \verbatim */
/* >          KU is INTEGER */
/* >     The number of superdiagonals within the band of A.  KU >= 0 */
/* > \endverbatim */
/* > */
/* > \param[in] NRHS */
/* > \verbatim */
/* >          NRHS is INTEGER */
/* >     The number of right-hand-sides, i.e., the number of columns of the */
/* >     matrix B. */
/* > \endverbatim */
/* > */
/* > \param[in] AB */
/* > \verbatim */
/* >          AB is REAL array, dimension (LDAB,N) */
/* >     On entry, the N-by-N matrix AB. */
/* > \endverbatim */
/* > */
/* > \param[in] LDAB */
/* > \verbatim */
/* >          LDAB is INTEGER */
/* >     The leading dimension of the array AB.  LDAB >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] AFB */
/* > \verbatim */
/* >          AFB is REAL array, dimension (LDAFB,N) */
/* >     The factors L and U from the factorization */
/* >     A = P*L*U as computed by SGBTRF. */
/* > \endverbatim */
/* > */
/* > \param[in] LDAFB */
/* > \verbatim */
/* >          LDAFB is INTEGER */
/* >     The leading dimension of the array AF.  LDAFB >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N) */
/* >     The pivot indices from the factorization A = P*L*U */
/* >     as computed by SGBTRF; row i of the matrix was interchanged */
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
/* >     On entry, the solution matrix X, as computed by SGBTRS. */
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
/* >     Determines which error bounds to return (see ERR_BNDS_NORM */
/* >     and ERR_BNDS_COMP). */
/* >     If N_NORMS >= 1 return normwise error bounds. */
/* >     If N_NORMS >= 2 return componentwise error bounds. */
/* > \endverbatim */
/* > */
/* > \param[in,out] ERR_BNDS_NORM */
/* > \verbatim */
/* >          ERR_BNDS_NORM is REAL array, dimension */
/* >                    (NRHS, N_ERR_BNDS) */
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
/* >     The first index in ERR_BNDS_NORM(i,:) corresponds to the ith */
/* >     right-hand side. */
/* > */
/* >     The second index in ERR_BNDS_NORM(:,err) contains the following */
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
/* > \param[in,out] ERR_BNDS_COMP */
/* > \verbatim */
/* >          ERR_BNDS_COMP is REAL array, dimension */
/* >                    (NRHS, N_ERR_BNDS) */
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
/* >     ERR_BNDS_COMP is not accessed.  If N_ERR_BNDS .LT. 3, then at most */
/* >     the first (:,N_ERR_BNDS) entries are returned. */
/* > */
/* >     The first index in ERR_BNDS_COMP(i,:) corresponds to the ith */
/* >     right-hand side. */
/* > */
/* >     The second index in ERR_BNDS_COMP(:,err) contains the following */
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
/* >     ERR_BNDS_NORM and ERR_BNDS_COMP may no longer be trustworthy. */
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
/* >       < 0:  if INFO = -i, the ith argument to SGBTRS had an illegal */
/* >             value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup realGBcomputational */

/*  ===================================================================== */
/* Subroutine */ int sla_gbrfsx_extended__(integer *prec_type__, integer *
	trans_type__, integer *n, integer *kl, integer *ku, integer *nrhs, 
	doublereal *ab, integer *ldab, doublereal *afb, integer *ldafb, 
	integer *ipiv, logical *colequ, doublereal *c__, doublereal *b, 
	integer *ldb, doublereal *y, integer *ldy, doublereal *berr_out__, 
	integer *n_norms__, doublereal *err_bnds_norm__, doublereal *
	err_bnds_comp__, doublereal *res, doublereal *ayb, doublereal *dy, 
	doublereal *y_tail__, doublereal *rcond, integer *ithresh, doublereal 
	*rthresh, doublereal *dz_ub__, logical *ignore_cwise__, integer *info)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, afb_dim1, afb_offset, b_dim1, b_offset, 
	    y_dim1, y_offset, err_bnds_norm_dim1, err_bnds_norm_offset, 
	    err_bnds_comp_dim1, err_bnds_comp_offset, i__1, i__2, i__3;
    doublereal d__1, d__2;
    char ch__1[1];

    /* Local variables */
    static doublereal dxratmax, dzratmax;
    static integer i__, j, m;
    extern /* Subroutine */ int sla_gbamv__(integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *);
    static logical incr_prec__;
    static doublereal prev_dz_z__, yk, final_dx_x__, final_dz_z__;
    extern /* Subroutine */ int sla_wwaddw__(integer *, doublereal *, 
	    doublereal *, doublereal *);
    static doublereal prevnormdx;
    static integer cnt;
    static doublereal dyk, eps, incr_thresh__, dx_x__, dz_z__, ymin;
    extern /* Subroutine */ int sla_lin_berr__(integer *, integer *, integer *
	    , doublereal *, doublereal *, doublereal *), blas_sgbmv_x__(
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *, doublereal *, integer *, doublereal *, integer *, doublereal *,
	     doublereal *, integer *, integer *);
    static integer y_prec_state__;
    extern /* Subroutine */ int blas_sgbmv2_x__(integer *, integer *, integer 
	    *, integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     integer *, integer *), sgbmv_(char *, integer *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    ftnlen);
    static doublereal dxrat, dzrat;
    static char trans[1];
    extern /* Subroutine */ int scopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static doublereal normx, normy;
    extern /* Subroutine */ int saxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int sgbtrs_(char *, integer *, integer *, integer 
	    *, integer *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen);
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

#line 484 "sla_gbrfsx_extended.f"
    /* Parameter adjustments */
#line 484 "sla_gbrfsx_extended.f"
    err_bnds_comp_dim1 = *nrhs;
#line 484 "sla_gbrfsx_extended.f"
    err_bnds_comp_offset = 1 + err_bnds_comp_dim1;
#line 484 "sla_gbrfsx_extended.f"
    err_bnds_comp__ -= err_bnds_comp_offset;
#line 484 "sla_gbrfsx_extended.f"
    err_bnds_norm_dim1 = *nrhs;
#line 484 "sla_gbrfsx_extended.f"
    err_bnds_norm_offset = 1 + err_bnds_norm_dim1;
#line 484 "sla_gbrfsx_extended.f"
    err_bnds_norm__ -= err_bnds_norm_offset;
#line 484 "sla_gbrfsx_extended.f"
    ab_dim1 = *ldab;
#line 484 "sla_gbrfsx_extended.f"
    ab_offset = 1 + ab_dim1;
#line 484 "sla_gbrfsx_extended.f"
    ab -= ab_offset;
#line 484 "sla_gbrfsx_extended.f"
    afb_dim1 = *ldafb;
#line 484 "sla_gbrfsx_extended.f"
    afb_offset = 1 + afb_dim1;
#line 484 "sla_gbrfsx_extended.f"
    afb -= afb_offset;
#line 484 "sla_gbrfsx_extended.f"
    --ipiv;
#line 484 "sla_gbrfsx_extended.f"
    --c__;
#line 484 "sla_gbrfsx_extended.f"
    b_dim1 = *ldb;
#line 484 "sla_gbrfsx_extended.f"
    b_offset = 1 + b_dim1;
#line 484 "sla_gbrfsx_extended.f"
    b -= b_offset;
#line 484 "sla_gbrfsx_extended.f"
    y_dim1 = *ldy;
#line 484 "sla_gbrfsx_extended.f"
    y_offset = 1 + y_dim1;
#line 484 "sla_gbrfsx_extended.f"
    y -= y_offset;
#line 484 "sla_gbrfsx_extended.f"
    --berr_out__;
#line 484 "sla_gbrfsx_extended.f"
    --res;
#line 484 "sla_gbrfsx_extended.f"
    --ayb;
#line 484 "sla_gbrfsx_extended.f"
    --dy;
#line 484 "sla_gbrfsx_extended.f"
    --y_tail__;
#line 484 "sla_gbrfsx_extended.f"

#line 484 "sla_gbrfsx_extended.f"
    /* Function Body */
#line 484 "sla_gbrfsx_extended.f"
    if (*info != 0) {
#line 484 "sla_gbrfsx_extended.f"
	return 0;
#line 484 "sla_gbrfsx_extended.f"
    }
#line 485 "sla_gbrfsx_extended.f"
    chla_transtype__(ch__1, (ftnlen)1, trans_type__);
#line 485 "sla_gbrfsx_extended.f"
    *(unsigned char *)trans = *(unsigned char *)&ch__1[0];
#line 486 "sla_gbrfsx_extended.f"
    eps = slamch_("Epsilon", (ftnlen)7);
#line 487 "sla_gbrfsx_extended.f"
    hugeval = slamch_("Overflow", (ftnlen)8);
/*     Force HUGEVAL to Inf */
#line 489 "sla_gbrfsx_extended.f"
    hugeval *= hugeval;
/*     Using HUGEVAL may lead to spurious underflows. */
#line 491 "sla_gbrfsx_extended.f"
    incr_thresh__ = (doublereal) (*n) * eps;
#line 492 "sla_gbrfsx_extended.f"
    m = *kl + *ku + 1;
#line 494 "sla_gbrfsx_extended.f"
    i__1 = *nrhs;
#line 494 "sla_gbrfsx_extended.f"
    for (j = 1; j <= i__1; ++j) {
#line 495 "sla_gbrfsx_extended.f"
	y_prec_state__ = 1;
#line 496 "sla_gbrfsx_extended.f"
	if (y_prec_state__ == 2) {
#line 497 "sla_gbrfsx_extended.f"
	    i__2 = *n;
#line 497 "sla_gbrfsx_extended.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 498 "sla_gbrfsx_extended.f"
		y_tail__[i__] = 0.;
#line 499 "sla_gbrfsx_extended.f"
	    }
#line 500 "sla_gbrfsx_extended.f"
	}
#line 502 "sla_gbrfsx_extended.f"
	dxrat = 0.;
#line 503 "sla_gbrfsx_extended.f"
	dxratmax = 0.;
#line 504 "sla_gbrfsx_extended.f"
	dzrat = 0.;
#line 505 "sla_gbrfsx_extended.f"
	dzratmax = 0.;
#line 506 "sla_gbrfsx_extended.f"
	final_dx_x__ = hugeval;
#line 507 "sla_gbrfsx_extended.f"
	final_dz_z__ = hugeval;
#line 508 "sla_gbrfsx_extended.f"
	prevnormdx = hugeval;
#line 509 "sla_gbrfsx_extended.f"
	prev_dz_z__ = hugeval;
#line 510 "sla_gbrfsx_extended.f"
	dz_z__ = hugeval;
#line 511 "sla_gbrfsx_extended.f"
	dx_x__ = hugeval;
#line 513 "sla_gbrfsx_extended.f"
	x_state__ = 1;
#line 514 "sla_gbrfsx_extended.f"
	z_state__ = 0;
#line 515 "sla_gbrfsx_extended.f"
	incr_prec__ = FALSE_;
#line 517 "sla_gbrfsx_extended.f"
	i__2 = *ithresh;
#line 517 "sla_gbrfsx_extended.f"
	for (cnt = 1; cnt <= i__2; ++cnt) {

/*        Compute residual RES = B_s - op(A_s) * Y, */
/*            op(A) = A, A**T, or A**H depending on TRANS (and type). */

#line 522 "sla_gbrfsx_extended.f"
	    scopy_(n, &b[j * b_dim1 + 1], &c__1, &res[1], &c__1);
#line 523 "sla_gbrfsx_extended.f"
	    if (y_prec_state__ == 0) {
#line 524 "sla_gbrfsx_extended.f"
		sgbmv_(trans, &m, n, kl, ku, &c_b6, &ab[ab_offset], ldab, &y[
			j * y_dim1 + 1], &c__1, &c_b8, &res[1], &c__1, (
			ftnlen)1);
#line 526 "sla_gbrfsx_extended.f"
	    } else if (y_prec_state__ == 1) {
#line 527 "sla_gbrfsx_extended.f"
		blas_sgbmv_x__(trans_type__, n, n, kl, ku, &c_b6, &ab[
			ab_offset], ldab, &y[j * y_dim1 + 1], &c__1, &c_b8, &
			res[1], &c__1, prec_type__);
#line 530 "sla_gbrfsx_extended.f"
	    } else {
#line 531 "sla_gbrfsx_extended.f"
		blas_sgbmv2_x__(trans_type__, n, n, kl, ku, &c_b6, &ab[
			ab_offset], ldab, &y[j * y_dim1 + 1], &y_tail__[1], &
			c__1, &c_b8, &res[1], &c__1, prec_type__);
#line 534 "sla_gbrfsx_extended.f"
	    }
/*        XXX: RES is no longer needed. */
#line 537 "sla_gbrfsx_extended.f"
	    scopy_(n, &res[1], &c__1, &dy[1], &c__1);
#line 538 "sla_gbrfsx_extended.f"
	    sgbtrs_(trans, n, kl, ku, &c__1, &afb[afb_offset], ldafb, &ipiv[1]
		    , &dy[1], n, info, (ftnlen)1);

/*         Calculate relative changes DX_X, DZ_Z and ratios DXRAT, DZRAT. */

#line 543 "sla_gbrfsx_extended.f"
	    normx = 0.;
#line 544 "sla_gbrfsx_extended.f"
	    normy = 0.;
#line 545 "sla_gbrfsx_extended.f"
	    normdx = 0.;
#line 546 "sla_gbrfsx_extended.f"
	    dz_z__ = 0.;
#line 547 "sla_gbrfsx_extended.f"
	    ymin = hugeval;
#line 549 "sla_gbrfsx_extended.f"
	    i__3 = *n;
#line 549 "sla_gbrfsx_extended.f"
	    for (i__ = 1; i__ <= i__3; ++i__) {
#line 550 "sla_gbrfsx_extended.f"
		yk = (d__1 = y[i__ + j * y_dim1], abs(d__1));
#line 551 "sla_gbrfsx_extended.f"
		dyk = (d__1 = dy[i__], abs(d__1));
#line 553 "sla_gbrfsx_extended.f"
		if (yk != 0.) {
/* Computing MAX */
#line 554 "sla_gbrfsx_extended.f"
		    d__1 = dz_z__, d__2 = dyk / yk;
#line 554 "sla_gbrfsx_extended.f"
		    dz_z__ = max(d__1,d__2);
#line 555 "sla_gbrfsx_extended.f"
		} else if (dyk != 0.) {
#line 556 "sla_gbrfsx_extended.f"
		    dz_z__ = hugeval;
#line 557 "sla_gbrfsx_extended.f"
		}
#line 559 "sla_gbrfsx_extended.f"
		ymin = min(ymin,yk);
#line 561 "sla_gbrfsx_extended.f"
		normy = max(normy,yk);
#line 563 "sla_gbrfsx_extended.f"
		if (*colequ) {
/* Computing MAX */
#line 564 "sla_gbrfsx_extended.f"
		    d__1 = normx, d__2 = yk * c__[i__];
#line 564 "sla_gbrfsx_extended.f"
		    normx = max(d__1,d__2);
/* Computing MAX */
#line 565 "sla_gbrfsx_extended.f"
		    d__1 = normdx, d__2 = dyk * c__[i__];
#line 565 "sla_gbrfsx_extended.f"
		    normdx = max(d__1,d__2);
#line 566 "sla_gbrfsx_extended.f"
		} else {
#line 567 "sla_gbrfsx_extended.f"
		    normx = normy;
#line 568 "sla_gbrfsx_extended.f"
		    normdx = max(normdx,dyk);
#line 569 "sla_gbrfsx_extended.f"
		}
#line 570 "sla_gbrfsx_extended.f"
	    }
#line 572 "sla_gbrfsx_extended.f"
	    if (normx != 0.) {
#line 573 "sla_gbrfsx_extended.f"
		dx_x__ = normdx / normx;
#line 574 "sla_gbrfsx_extended.f"
	    } else if (normdx == 0.) {
#line 575 "sla_gbrfsx_extended.f"
		dx_x__ = 0.;
#line 576 "sla_gbrfsx_extended.f"
	    } else {
#line 577 "sla_gbrfsx_extended.f"
		dx_x__ = hugeval;
#line 578 "sla_gbrfsx_extended.f"
	    }
#line 580 "sla_gbrfsx_extended.f"
	    dxrat = normdx / prevnormdx;
#line 581 "sla_gbrfsx_extended.f"
	    dzrat = dz_z__ / prev_dz_z__;

/*         Check termination criteria. */

#line 585 "sla_gbrfsx_extended.f"
	    if (! (*ignore_cwise__) && ymin * *rcond < incr_thresh__ * normy 
		    && y_prec_state__ < 2) {
#line 585 "sla_gbrfsx_extended.f"
		incr_prec__ = TRUE_;
#line 585 "sla_gbrfsx_extended.f"
	    }
#line 590 "sla_gbrfsx_extended.f"
	    if (x_state__ == 3 && dxrat <= *rthresh) {
#line 590 "sla_gbrfsx_extended.f"
		x_state__ = 1;
#line 590 "sla_gbrfsx_extended.f"
	    }
#line 592 "sla_gbrfsx_extended.f"
	    if (x_state__ == 1) {
#line 593 "sla_gbrfsx_extended.f"
		if (dx_x__ <= eps) {
#line 594 "sla_gbrfsx_extended.f"
		    x_state__ = 2;
#line 595 "sla_gbrfsx_extended.f"
		} else if (dxrat > *rthresh) {
#line 596 "sla_gbrfsx_extended.f"
		    if (y_prec_state__ != 2) {
#line 597 "sla_gbrfsx_extended.f"
			incr_prec__ = TRUE_;
#line 598 "sla_gbrfsx_extended.f"
		    } else {
#line 599 "sla_gbrfsx_extended.f"
			x_state__ = 3;
#line 600 "sla_gbrfsx_extended.f"
		    }
#line 601 "sla_gbrfsx_extended.f"
		} else {
#line 602 "sla_gbrfsx_extended.f"
		    if (dxrat > dxratmax) {
#line 602 "sla_gbrfsx_extended.f"
			dxratmax = dxrat;
#line 602 "sla_gbrfsx_extended.f"
		    }
#line 603 "sla_gbrfsx_extended.f"
		}
#line 604 "sla_gbrfsx_extended.f"
		if (x_state__ > 1) {
#line 604 "sla_gbrfsx_extended.f"
		    final_dx_x__ = dx_x__;
#line 604 "sla_gbrfsx_extended.f"
		}
#line 605 "sla_gbrfsx_extended.f"
	    }
#line 607 "sla_gbrfsx_extended.f"
	    if (z_state__ == 0 && dz_z__ <= *dz_ub__) {
#line 607 "sla_gbrfsx_extended.f"
		z_state__ = 1;
#line 607 "sla_gbrfsx_extended.f"
	    }
#line 609 "sla_gbrfsx_extended.f"
	    if (z_state__ == 3 && dzrat <= *rthresh) {
#line 609 "sla_gbrfsx_extended.f"
		z_state__ = 1;
#line 609 "sla_gbrfsx_extended.f"
	    }
#line 611 "sla_gbrfsx_extended.f"
	    if (z_state__ == 1) {
#line 612 "sla_gbrfsx_extended.f"
		if (dz_z__ <= eps) {
#line 613 "sla_gbrfsx_extended.f"
		    z_state__ = 2;
#line 614 "sla_gbrfsx_extended.f"
		} else if (dz_z__ > *dz_ub__) {
#line 615 "sla_gbrfsx_extended.f"
		    z_state__ = 0;
#line 616 "sla_gbrfsx_extended.f"
		    dzratmax = 0.;
#line 617 "sla_gbrfsx_extended.f"
		    final_dz_z__ = hugeval;
#line 618 "sla_gbrfsx_extended.f"
		} else if (dzrat > *rthresh) {
#line 619 "sla_gbrfsx_extended.f"
		    if (y_prec_state__ != 2) {
#line 620 "sla_gbrfsx_extended.f"
			incr_prec__ = TRUE_;
#line 621 "sla_gbrfsx_extended.f"
		    } else {
#line 622 "sla_gbrfsx_extended.f"
			z_state__ = 3;
#line 623 "sla_gbrfsx_extended.f"
		    }
#line 624 "sla_gbrfsx_extended.f"
		} else {
#line 625 "sla_gbrfsx_extended.f"
		    if (dzrat > dzratmax) {
#line 625 "sla_gbrfsx_extended.f"
			dzratmax = dzrat;
#line 625 "sla_gbrfsx_extended.f"
		    }
#line 626 "sla_gbrfsx_extended.f"
		}
#line 627 "sla_gbrfsx_extended.f"
		if (z_state__ > 1) {
#line 627 "sla_gbrfsx_extended.f"
		    final_dz_z__ = dz_z__;
#line 627 "sla_gbrfsx_extended.f"
		}
#line 628 "sla_gbrfsx_extended.f"
	    }

/*           Exit if both normwise and componentwise stopped working, */
/*           but if componentwise is unstable, let it go at least two */
/*           iterations. */

#line 634 "sla_gbrfsx_extended.f"
	    if (x_state__ != 1) {
#line 635 "sla_gbrfsx_extended.f"
		if (*ignore_cwise__) {
#line 635 "sla_gbrfsx_extended.f"
		    goto L666;
#line 635 "sla_gbrfsx_extended.f"
		}
#line 636 "sla_gbrfsx_extended.f"
		if (z_state__ == 3 || z_state__ == 2) {
#line 636 "sla_gbrfsx_extended.f"
		    goto L666;
#line 636 "sla_gbrfsx_extended.f"
		}
#line 638 "sla_gbrfsx_extended.f"
		if (z_state__ == 0 && cnt > 1) {
#line 638 "sla_gbrfsx_extended.f"
		    goto L666;
#line 638 "sla_gbrfsx_extended.f"
		}
#line 639 "sla_gbrfsx_extended.f"
	    }
#line 641 "sla_gbrfsx_extended.f"
	    if (incr_prec__) {
#line 642 "sla_gbrfsx_extended.f"
		incr_prec__ = FALSE_;
#line 643 "sla_gbrfsx_extended.f"
		++y_prec_state__;
#line 644 "sla_gbrfsx_extended.f"
		i__3 = *n;
#line 644 "sla_gbrfsx_extended.f"
		for (i__ = 1; i__ <= i__3; ++i__) {
#line 645 "sla_gbrfsx_extended.f"
		    y_tail__[i__] = 0.;
#line 646 "sla_gbrfsx_extended.f"
		}
#line 647 "sla_gbrfsx_extended.f"
	    }
#line 649 "sla_gbrfsx_extended.f"
	    prevnormdx = normdx;
#line 650 "sla_gbrfsx_extended.f"
	    prev_dz_z__ = dz_z__;

/*           Update soluton. */

#line 654 "sla_gbrfsx_extended.f"
	    if (y_prec_state__ < 2) {
#line 655 "sla_gbrfsx_extended.f"
		saxpy_(n, &c_b8, &dy[1], &c__1, &y[j * y_dim1 + 1], &c__1);
#line 656 "sla_gbrfsx_extended.f"
	    } else {
#line 657 "sla_gbrfsx_extended.f"
		sla_wwaddw__(n, &y[j * y_dim1 + 1], &y_tail__[1], &dy[1]);
#line 658 "sla_gbrfsx_extended.f"
	    }
#line 660 "sla_gbrfsx_extended.f"
	}
/*        Target of "IF (Z_STOP .AND. X_STOP)".  Sun's f77 won't EXIT. */
#line 662 "sla_gbrfsx_extended.f"
L666:

/*     Set final_* when cnt hits ithresh. */

#line 666 "sla_gbrfsx_extended.f"
	if (x_state__ == 1) {
#line 666 "sla_gbrfsx_extended.f"
	    final_dx_x__ = dx_x__;
#line 666 "sla_gbrfsx_extended.f"
	}
#line 667 "sla_gbrfsx_extended.f"
	if (z_state__ == 1) {
#line 667 "sla_gbrfsx_extended.f"
	    final_dz_z__ = dz_z__;
#line 667 "sla_gbrfsx_extended.f"
	}

/*     Compute error bounds. */

#line 671 "sla_gbrfsx_extended.f"
	if (*n_norms__ >= 1) {
#line 672 "sla_gbrfsx_extended.f"
	    err_bnds_norm__[j + (err_bnds_norm_dim1 << 1)] = final_dx_x__ / (
		    1 - dxratmax);
#line 674 "sla_gbrfsx_extended.f"
	}
#line 675 "sla_gbrfsx_extended.f"
	if (*n_norms__ >= 2) {
#line 676 "sla_gbrfsx_extended.f"
	    err_bnds_comp__[j + (err_bnds_comp_dim1 << 1)] = final_dz_z__ / (
		    1 - dzratmax);
#line 678 "sla_gbrfsx_extended.f"
	}

/*     Compute componentwise relative backward error from formula */
/*         max(i) ( abs(R(i)) / ( abs(op(A_s))*abs(Y) + abs(B_s) )(i) ) */
/*     where abs(Z) is the componentwise absolute value of the matrix */
/*     or vector Z. */

/*        Compute residual RES = B_s - op(A_s) * Y, */
/*            op(A) = A, A**T, or A**H depending on TRANS (and type). */

#line 688 "sla_gbrfsx_extended.f"
	scopy_(n, &b[j * b_dim1 + 1], &c__1, &res[1], &c__1);
#line 689 "sla_gbrfsx_extended.f"
	sgbmv_(trans, n, n, kl, ku, &c_b6, &ab[ab_offset], ldab, &y[j * 
		y_dim1 + 1], &c__1, &c_b8, &res[1], &c__1, (ftnlen)1);
#line 692 "sla_gbrfsx_extended.f"
	i__2 = *n;
#line 692 "sla_gbrfsx_extended.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 693 "sla_gbrfsx_extended.f"
	    ayb[i__] = (d__1 = b[i__ + j * b_dim1], abs(d__1));
#line 694 "sla_gbrfsx_extended.f"
	}

/*     Compute abs(op(A_s))*abs(Y) + abs(B_s). */

#line 698 "sla_gbrfsx_extended.f"
	sla_gbamv__(trans_type__, n, n, kl, ku, &c_b8, &ab[ab_offset], ldab, &
		y[j * y_dim1 + 1], &c__1, &c_b8, &ayb[1], &c__1);
#line 701 "sla_gbrfsx_extended.f"
	sla_lin_berr__(n, n, &c__1, &res[1], &ayb[1], &berr_out__[j]);

/*     End of loop for each RHS */

#line 705 "sla_gbrfsx_extended.f"
    }

#line 707 "sla_gbrfsx_extended.f"
    return 0;
} /* sla_gbrfsx_extended__ */

