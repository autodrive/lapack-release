#line 1 "sla_syrfsx_extended.f"
/* sla_syrfsx_extended.f -- translated by f2c (version 20100827).
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

#line 1 "sla_syrfsx_extended.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b12 = -1.;
static doublereal c_b14 = 1.;

/* > \brief \b SLA_SYRFSX_EXTENDED improves the computed solution to a system of linear equations for symmetri
c indefinite matrices by performing extra-precise iterative refinement and provides error bounds and b
ackward error estimates for the solution. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLA_SYRFSX_EXTENDED + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sla_syr
fsx_extended.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sla_syr
fsx_extended.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sla_syr
fsx_extended.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLA_SYRFSX_EXTENDED( PREC_TYPE, UPLO, N, NRHS, A, LDA, */
/*                                       AF, LDAF, IPIV, COLEQU, C, B, LDB, */
/*                                       Y, LDY, BERR_OUT, N_NORMS, */
/*                                       ERR_BNDS_NORM, ERR_BNDS_COMP, RES, */
/*                                       AYB, DY, Y_TAIL, RCOND, ITHRESH, */
/*                                       RTHRESH, DZ_UB, IGNORE_CWISE, */
/*                                       INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, LDAF, LDB, LDY, N, NRHS, PREC_TYPE, */
/*      $                   N_NORMS, ITHRESH */
/*       CHARACTER          UPLO */
/*       LOGICAL            COLEQU, IGNORE_CWISE */
/*       REAL               RTHRESH, DZ_UB */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       REAL               A( LDA, * ), AF( LDAF, * ), B( LDB, * ), */
/*      $                   Y( LDY, * ), RES( * ), DY( * ), Y_TAIL( * ) */
/*       REAL               C( * ), AYB( * ), RCOND, BERR_OUT( * ), */
/*      $                   ERR_BNDS_NORM( NRHS, * ), */
/*      $                   ERR_BNDS_COMP( NRHS, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > */
/* > SLA_SYRFSX_EXTENDED improves the computed solution to a system of */
/* > linear equations by performing extra-precise iterative refinement */
/* > and provides error bounds and backward error estimates for the solution. */
/* > This subroutine is called by SSYRFSX to perform iterative refinement. */
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
/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >       = 'U':  Upper triangle of A is stored; */
/* >       = 'L':  Lower triangle of A is stored. */
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
/* >     The block diagonal matrix D and the multipliers used to */
/* >     obtain the factor U or L as computed by SSYTRF. */
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
/* >     Details of the interchanges and the block structure of D */
/* >     as determined by SSYTRF. */
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
/* >     On entry, the solution matrix X, as computed by SSYTRS. */
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
/* >          ERR_BNDS_NORM is REAL array, dimension (NRHS, N_ERR_BNDS) */
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
/* >          ERR_BNDS_COMP is REAL array, dimension (NRHS, N_ERR_BNDS) */
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
/* >       < 0:  if INFO = -i, the ith argument to SLA_SYRFSX_EXTENDED had an illegal */
/* >             value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup realSYcomputational */

/*  ===================================================================== */
/* Subroutine */ int sla_syrfsx_extended__(integer *prec_type__, char *uplo, 
	integer *n, integer *nrhs, doublereal *a, integer *lda, doublereal *
	af, integer *ldaf, integer *ipiv, logical *colequ, doublereal *c__, 
	doublereal *b, integer *ldb, doublereal *y, integer *ldy, doublereal *
	berr_out__, integer *n_norms__, doublereal *err_bnds_norm__, 
	doublereal *err_bnds_comp__, doublereal *res, doublereal *ayb, 
	doublereal *dy, doublereal *y_tail__, doublereal *rcond, integer *
	ithresh, doublereal *rthresh, doublereal *dz_ub__, logical *
	ignore_cwise__, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, af_dim1, af_offset, b_dim1, b_offset, y_dim1, 
	    y_offset, err_bnds_norm_dim1, err_bnds_norm_offset, 
	    err_bnds_comp_dim1, err_bnds_comp_offset, i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Local variables */
    static doublereal dxratmax, dzratmax;
    static integer i__, j;
    static logical incr_prec__;
    extern /* Subroutine */ int sla_syamv__(integer *, integer *, doublereal *
	    , doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *);
    static doublereal prev_dz_z__, yk, final_dx_x__, final_dz_z__;
    extern /* Subroutine */ int sla_wwaddw__(integer *, doublereal *, 
	    doublereal *, doublereal *);
    static doublereal prevnormdx;
    static integer cnt;
    static doublereal dyk, eps, incr_thresh__, dx_x__, dz_z__, ymin;
    extern /* Subroutine */ int sla_lin_berr__(integer *, integer *, integer *
	    , doublereal *, doublereal *, doublereal *);
    static integer y_prec_state__, uplo2;
    extern /* Subroutine */ int blas_ssymv_x__(integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal dxrat, dzrat;
    extern /* Subroutine */ int blas_ssymv2_x__(integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     integer *, doublereal *, doublereal *, integer *, integer *);
    static logical upper;
    extern /* Subroutine */ int scopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static doublereal normx, normy;
    extern /* Subroutine */ int saxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *), ssymv_(char *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen);
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal normdx;
    extern /* Subroutine */ int ssytrs_(char *, integer *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, ftnlen);
    static doublereal hugeval;
    extern integer ilauplo_(char *, ftnlen);
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
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 472 "sla_syrfsx_extended.f"
    /* Parameter adjustments */
#line 472 "sla_syrfsx_extended.f"
    err_bnds_comp_dim1 = *nrhs;
#line 472 "sla_syrfsx_extended.f"
    err_bnds_comp_offset = 1 + err_bnds_comp_dim1;
#line 472 "sla_syrfsx_extended.f"
    err_bnds_comp__ -= err_bnds_comp_offset;
#line 472 "sla_syrfsx_extended.f"
    err_bnds_norm_dim1 = *nrhs;
#line 472 "sla_syrfsx_extended.f"
    err_bnds_norm_offset = 1 + err_bnds_norm_dim1;
#line 472 "sla_syrfsx_extended.f"
    err_bnds_norm__ -= err_bnds_norm_offset;
#line 472 "sla_syrfsx_extended.f"
    a_dim1 = *lda;
#line 472 "sla_syrfsx_extended.f"
    a_offset = 1 + a_dim1;
#line 472 "sla_syrfsx_extended.f"
    a -= a_offset;
#line 472 "sla_syrfsx_extended.f"
    af_dim1 = *ldaf;
#line 472 "sla_syrfsx_extended.f"
    af_offset = 1 + af_dim1;
#line 472 "sla_syrfsx_extended.f"
    af -= af_offset;
#line 472 "sla_syrfsx_extended.f"
    --ipiv;
#line 472 "sla_syrfsx_extended.f"
    --c__;
#line 472 "sla_syrfsx_extended.f"
    b_dim1 = *ldb;
#line 472 "sla_syrfsx_extended.f"
    b_offset = 1 + b_dim1;
#line 472 "sla_syrfsx_extended.f"
    b -= b_offset;
#line 472 "sla_syrfsx_extended.f"
    y_dim1 = *ldy;
#line 472 "sla_syrfsx_extended.f"
    y_offset = 1 + y_dim1;
#line 472 "sla_syrfsx_extended.f"
    y -= y_offset;
#line 472 "sla_syrfsx_extended.f"
    --berr_out__;
#line 472 "sla_syrfsx_extended.f"
    --res;
#line 472 "sla_syrfsx_extended.f"
    --ayb;
#line 472 "sla_syrfsx_extended.f"
    --dy;
#line 472 "sla_syrfsx_extended.f"
    --y_tail__;
#line 472 "sla_syrfsx_extended.f"

#line 472 "sla_syrfsx_extended.f"
    /* Function Body */
#line 472 "sla_syrfsx_extended.f"
    *info = 0;
#line 473 "sla_syrfsx_extended.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 474 "sla_syrfsx_extended.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 475 "sla_syrfsx_extended.f"
	*info = -2;
#line 476 "sla_syrfsx_extended.f"
    } else if (*n < 0) {
#line 477 "sla_syrfsx_extended.f"
	*info = -3;
#line 478 "sla_syrfsx_extended.f"
    } else if (*nrhs < 0) {
#line 479 "sla_syrfsx_extended.f"
	*info = -4;
#line 480 "sla_syrfsx_extended.f"
    } else if (*lda < max(1,*n)) {
#line 481 "sla_syrfsx_extended.f"
	*info = -6;
#line 482 "sla_syrfsx_extended.f"
    } else if (*ldaf < max(1,*n)) {
#line 483 "sla_syrfsx_extended.f"
	*info = -8;
#line 484 "sla_syrfsx_extended.f"
    } else if (*ldb < max(1,*n)) {
#line 485 "sla_syrfsx_extended.f"
	*info = -13;
#line 486 "sla_syrfsx_extended.f"
    } else if (*ldy < max(1,*n)) {
#line 487 "sla_syrfsx_extended.f"
	*info = -15;
#line 488 "sla_syrfsx_extended.f"
    }
#line 489 "sla_syrfsx_extended.f"
    if (*info != 0) {
#line 490 "sla_syrfsx_extended.f"
	i__1 = -(*info);
#line 490 "sla_syrfsx_extended.f"
	xerbla_("SLA_SYRFSX_EXTENDED", &i__1, (ftnlen)19);
#line 491 "sla_syrfsx_extended.f"
	return 0;
#line 492 "sla_syrfsx_extended.f"
    }
#line 493 "sla_syrfsx_extended.f"
    eps = slamch_("Epsilon", (ftnlen)7);
#line 494 "sla_syrfsx_extended.f"
    hugeval = slamch_("Overflow", (ftnlen)8);
/*     Force HUGEVAL to Inf */
#line 496 "sla_syrfsx_extended.f"
    hugeval *= hugeval;
/*     Using HUGEVAL may lead to spurious underflows. */
#line 498 "sla_syrfsx_extended.f"
    incr_thresh__ = (doublereal) (*n) * eps;
#line 500 "sla_syrfsx_extended.f"
    if (lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 501 "sla_syrfsx_extended.f"
	uplo2 = ilauplo_("L", (ftnlen)1);
#line 502 "sla_syrfsx_extended.f"
    } else {
#line 503 "sla_syrfsx_extended.f"
	uplo2 = ilauplo_("U", (ftnlen)1);
#line 504 "sla_syrfsx_extended.f"
    }
#line 506 "sla_syrfsx_extended.f"
    i__1 = *nrhs;
#line 506 "sla_syrfsx_extended.f"
    for (j = 1; j <= i__1; ++j) {
#line 507 "sla_syrfsx_extended.f"
	y_prec_state__ = 1;
#line 508 "sla_syrfsx_extended.f"
	if (y_prec_state__ == 2) {
#line 509 "sla_syrfsx_extended.f"
	    i__2 = *n;
#line 509 "sla_syrfsx_extended.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 510 "sla_syrfsx_extended.f"
		y_tail__[i__] = 0.;
#line 511 "sla_syrfsx_extended.f"
	    }
#line 512 "sla_syrfsx_extended.f"
	}
#line 514 "sla_syrfsx_extended.f"
	dxrat = 0.;
#line 515 "sla_syrfsx_extended.f"
	dxratmax = 0.;
#line 516 "sla_syrfsx_extended.f"
	dzrat = 0.;
#line 517 "sla_syrfsx_extended.f"
	dzratmax = 0.;
#line 518 "sla_syrfsx_extended.f"
	final_dx_x__ = hugeval;
#line 519 "sla_syrfsx_extended.f"
	final_dz_z__ = hugeval;
#line 520 "sla_syrfsx_extended.f"
	prevnormdx = hugeval;
#line 521 "sla_syrfsx_extended.f"
	prev_dz_z__ = hugeval;
#line 522 "sla_syrfsx_extended.f"
	dz_z__ = hugeval;
#line 523 "sla_syrfsx_extended.f"
	dx_x__ = hugeval;
#line 525 "sla_syrfsx_extended.f"
	x_state__ = 1;
#line 526 "sla_syrfsx_extended.f"
	z_state__ = 0;
#line 527 "sla_syrfsx_extended.f"
	incr_prec__ = FALSE_;
#line 529 "sla_syrfsx_extended.f"
	i__2 = *ithresh;
#line 529 "sla_syrfsx_extended.f"
	for (cnt = 1; cnt <= i__2; ++cnt) {

/*        Compute residual RES = B_s - op(A_s) * Y, */
/*            op(A) = A, A**T, or A**H depending on TRANS (and type). */

#line 534 "sla_syrfsx_extended.f"
	    scopy_(n, &b[j * b_dim1 + 1], &c__1, &res[1], &c__1);
#line 535 "sla_syrfsx_extended.f"
	    if (y_prec_state__ == 0) {
#line 536 "sla_syrfsx_extended.f"
		ssymv_(uplo, n, &c_b12, &a[a_offset], lda, &y[j * y_dim1 + 1],
			 &c__1, &c_b14, &res[1], &c__1, (ftnlen)1);
#line 538 "sla_syrfsx_extended.f"
	    } else if (y_prec_state__ == 1) {
#line 539 "sla_syrfsx_extended.f"
		blas_ssymv_x__(&uplo2, n, &c_b12, &a[a_offset], lda, &y[j * 
			y_dim1 + 1], &c__1, &c_b14, &res[1], &c__1, 
			prec_type__);
#line 541 "sla_syrfsx_extended.f"
	    } else {
#line 542 "sla_syrfsx_extended.f"
		blas_ssymv2_x__(&uplo2, n, &c_b12, &a[a_offset], lda, &y[j * 
			y_dim1 + 1], &y_tail__[1], &c__1, &c_b14, &res[1], &
			c__1, prec_type__);
#line 544 "sla_syrfsx_extended.f"
	    }
/*         XXX: RES is no longer needed. */
#line 547 "sla_syrfsx_extended.f"
	    scopy_(n, &res[1], &c__1, &dy[1], &c__1);
#line 548 "sla_syrfsx_extended.f"
	    ssytrs_(uplo, n, &c__1, &af[af_offset], ldaf, &ipiv[1], &dy[1], n,
		     info, (ftnlen)1);

/*         Calculate relative changes DX_X, DZ_Z and ratios DXRAT, DZRAT. */

#line 552 "sla_syrfsx_extended.f"
	    normx = 0.;
#line 553 "sla_syrfsx_extended.f"
	    normy = 0.;
#line 554 "sla_syrfsx_extended.f"
	    normdx = 0.;
#line 555 "sla_syrfsx_extended.f"
	    dz_z__ = 0.;
#line 556 "sla_syrfsx_extended.f"
	    ymin = hugeval;
#line 558 "sla_syrfsx_extended.f"
	    i__3 = *n;
#line 558 "sla_syrfsx_extended.f"
	    for (i__ = 1; i__ <= i__3; ++i__) {
#line 559 "sla_syrfsx_extended.f"
		yk = (d__1 = y[i__ + j * y_dim1], abs(d__1));
#line 560 "sla_syrfsx_extended.f"
		dyk = (d__1 = dy[i__], abs(d__1));
#line 562 "sla_syrfsx_extended.f"
		if (yk != 0.) {
/* Computing MAX */
#line 563 "sla_syrfsx_extended.f"
		    d__1 = dz_z__, d__2 = dyk / yk;
#line 563 "sla_syrfsx_extended.f"
		    dz_z__ = max(d__1,d__2);
#line 564 "sla_syrfsx_extended.f"
		} else if (dyk != 0.) {
#line 565 "sla_syrfsx_extended.f"
		    dz_z__ = hugeval;
#line 566 "sla_syrfsx_extended.f"
		}
#line 568 "sla_syrfsx_extended.f"
		ymin = min(ymin,yk);
#line 570 "sla_syrfsx_extended.f"
		normy = max(normy,yk);
#line 572 "sla_syrfsx_extended.f"
		if (*colequ) {
/* Computing MAX */
#line 573 "sla_syrfsx_extended.f"
		    d__1 = normx, d__2 = yk * c__[i__];
#line 573 "sla_syrfsx_extended.f"
		    normx = max(d__1,d__2);
/* Computing MAX */
#line 574 "sla_syrfsx_extended.f"
		    d__1 = normdx, d__2 = dyk * c__[i__];
#line 574 "sla_syrfsx_extended.f"
		    normdx = max(d__1,d__2);
#line 575 "sla_syrfsx_extended.f"
		} else {
#line 576 "sla_syrfsx_extended.f"
		    normx = normy;
#line 577 "sla_syrfsx_extended.f"
		    normdx = max(normdx,dyk);
#line 578 "sla_syrfsx_extended.f"
		}
#line 579 "sla_syrfsx_extended.f"
	    }
#line 581 "sla_syrfsx_extended.f"
	    if (normx != 0.) {
#line 582 "sla_syrfsx_extended.f"
		dx_x__ = normdx / normx;
#line 583 "sla_syrfsx_extended.f"
	    } else if (normdx == 0.) {
#line 584 "sla_syrfsx_extended.f"
		dx_x__ = 0.;
#line 585 "sla_syrfsx_extended.f"
	    } else {
#line 586 "sla_syrfsx_extended.f"
		dx_x__ = hugeval;
#line 587 "sla_syrfsx_extended.f"
	    }
#line 589 "sla_syrfsx_extended.f"
	    dxrat = normdx / prevnormdx;
#line 590 "sla_syrfsx_extended.f"
	    dzrat = dz_z__ / prev_dz_z__;

/*         Check termination criteria. */

#line 594 "sla_syrfsx_extended.f"
	    if (ymin * *rcond < incr_thresh__ * normy && y_prec_state__ < 2) {
#line 594 "sla_syrfsx_extended.f"
		incr_prec__ = TRUE_;
#line 594 "sla_syrfsx_extended.f"
	    }
#line 598 "sla_syrfsx_extended.f"
	    if (x_state__ == 3 && dxrat <= *rthresh) {
#line 598 "sla_syrfsx_extended.f"
		x_state__ = 1;
#line 598 "sla_syrfsx_extended.f"
	    }
#line 600 "sla_syrfsx_extended.f"
	    if (x_state__ == 1) {
#line 601 "sla_syrfsx_extended.f"
		if (dx_x__ <= eps) {
#line 602 "sla_syrfsx_extended.f"
		    x_state__ = 2;
#line 603 "sla_syrfsx_extended.f"
		} else if (dxrat > *rthresh) {
#line 604 "sla_syrfsx_extended.f"
		    if (y_prec_state__ != 2) {
#line 605 "sla_syrfsx_extended.f"
			incr_prec__ = TRUE_;
#line 606 "sla_syrfsx_extended.f"
		    } else {
#line 607 "sla_syrfsx_extended.f"
			x_state__ = 3;
#line 608 "sla_syrfsx_extended.f"
		    }
#line 609 "sla_syrfsx_extended.f"
		} else {
#line 610 "sla_syrfsx_extended.f"
		    if (dxrat > dxratmax) {
#line 610 "sla_syrfsx_extended.f"
			dxratmax = dxrat;
#line 610 "sla_syrfsx_extended.f"
		    }
#line 611 "sla_syrfsx_extended.f"
		}
#line 612 "sla_syrfsx_extended.f"
		if (x_state__ > 1) {
#line 612 "sla_syrfsx_extended.f"
		    final_dx_x__ = dx_x__;
#line 612 "sla_syrfsx_extended.f"
		}
#line 613 "sla_syrfsx_extended.f"
	    }
#line 615 "sla_syrfsx_extended.f"
	    if (z_state__ == 0 && dz_z__ <= *dz_ub__) {
#line 615 "sla_syrfsx_extended.f"
		z_state__ = 1;
#line 615 "sla_syrfsx_extended.f"
	    }
#line 617 "sla_syrfsx_extended.f"
	    if (z_state__ == 3 && dzrat <= *rthresh) {
#line 617 "sla_syrfsx_extended.f"
		z_state__ = 1;
#line 617 "sla_syrfsx_extended.f"
	    }
#line 619 "sla_syrfsx_extended.f"
	    if (z_state__ == 1) {
#line 620 "sla_syrfsx_extended.f"
		if (dz_z__ <= eps) {
#line 621 "sla_syrfsx_extended.f"
		    z_state__ = 2;
#line 622 "sla_syrfsx_extended.f"
		} else if (dz_z__ > *dz_ub__) {
#line 623 "sla_syrfsx_extended.f"
		    z_state__ = 0;
#line 624 "sla_syrfsx_extended.f"
		    dzratmax = 0.;
#line 625 "sla_syrfsx_extended.f"
		    final_dz_z__ = hugeval;
#line 626 "sla_syrfsx_extended.f"
		} else if (dzrat > *rthresh) {
#line 627 "sla_syrfsx_extended.f"
		    if (y_prec_state__ != 2) {
#line 628 "sla_syrfsx_extended.f"
			incr_prec__ = TRUE_;
#line 629 "sla_syrfsx_extended.f"
		    } else {
#line 630 "sla_syrfsx_extended.f"
			z_state__ = 3;
#line 631 "sla_syrfsx_extended.f"
		    }
#line 632 "sla_syrfsx_extended.f"
		} else {
#line 633 "sla_syrfsx_extended.f"
		    if (dzrat > dzratmax) {
#line 633 "sla_syrfsx_extended.f"
			dzratmax = dzrat;
#line 633 "sla_syrfsx_extended.f"
		    }
#line 634 "sla_syrfsx_extended.f"
		}
#line 635 "sla_syrfsx_extended.f"
		if (z_state__ > 1) {
#line 635 "sla_syrfsx_extended.f"
		    final_dz_z__ = dz_z__;
#line 635 "sla_syrfsx_extended.f"
		}
#line 636 "sla_syrfsx_extended.f"
	    }
#line 638 "sla_syrfsx_extended.f"
	    if (x_state__ != 1 && (*ignore_cwise__ || z_state__ != 1)) {
#line 638 "sla_syrfsx_extended.f"
		goto L666;
#line 638 "sla_syrfsx_extended.f"
	    }
#line 642 "sla_syrfsx_extended.f"
	    if (incr_prec__) {
#line 643 "sla_syrfsx_extended.f"
		incr_prec__ = FALSE_;
#line 644 "sla_syrfsx_extended.f"
		++y_prec_state__;
#line 645 "sla_syrfsx_extended.f"
		i__3 = *n;
#line 645 "sla_syrfsx_extended.f"
		for (i__ = 1; i__ <= i__3; ++i__) {
#line 646 "sla_syrfsx_extended.f"
		    y_tail__[i__] = 0.;
#line 647 "sla_syrfsx_extended.f"
		}
#line 648 "sla_syrfsx_extended.f"
	    }
#line 650 "sla_syrfsx_extended.f"
	    prevnormdx = normdx;
#line 651 "sla_syrfsx_extended.f"
	    prev_dz_z__ = dz_z__;

/*           Update soluton. */

#line 655 "sla_syrfsx_extended.f"
	    if (y_prec_state__ < 2) {
#line 656 "sla_syrfsx_extended.f"
		saxpy_(n, &c_b14, &dy[1], &c__1, &y[j * y_dim1 + 1], &c__1);
#line 657 "sla_syrfsx_extended.f"
	    } else {
#line 658 "sla_syrfsx_extended.f"
		sla_wwaddw__(n, &y[j * y_dim1 + 1], &y_tail__[1], &dy[1]);
#line 659 "sla_syrfsx_extended.f"
	    }
#line 661 "sla_syrfsx_extended.f"
	}
/*        Target of "IF (Z_STOP .AND. X_STOP)".  Sun's f77 won't EXIT. */
#line 663 "sla_syrfsx_extended.f"
L666:

/*     Set final_* when cnt hits ithresh. */

#line 667 "sla_syrfsx_extended.f"
	if (x_state__ == 1) {
#line 667 "sla_syrfsx_extended.f"
	    final_dx_x__ = dx_x__;
#line 667 "sla_syrfsx_extended.f"
	}
#line 668 "sla_syrfsx_extended.f"
	if (z_state__ == 1) {
#line 668 "sla_syrfsx_extended.f"
	    final_dz_z__ = dz_z__;
#line 668 "sla_syrfsx_extended.f"
	}

/*     Compute error bounds. */

#line 672 "sla_syrfsx_extended.f"
	if (*n_norms__ >= 1) {
#line 673 "sla_syrfsx_extended.f"
	    err_bnds_norm__[j + (err_bnds_norm_dim1 << 1)] = final_dx_x__ / (
		    1 - dxratmax);
#line 675 "sla_syrfsx_extended.f"
	}
#line 676 "sla_syrfsx_extended.f"
	if (*n_norms__ >= 2) {
#line 677 "sla_syrfsx_extended.f"
	    err_bnds_comp__[j + (err_bnds_comp_dim1 << 1)] = final_dz_z__ / (
		    1 - dzratmax);
#line 679 "sla_syrfsx_extended.f"
	}

/*     Compute componentwise relative backward error from formula */
/*         max(i) ( abs(R(i)) / ( abs(op(A_s))*abs(Y) + abs(B_s) )(i) ) */
/*     where abs(Z) is the componentwise absolute value of the matrix */
/*     or vector Z. */

/*        Compute residual RES = B_s - op(A_s) * Y, */
/*            op(A) = A, A**T, or A**H depending on TRANS (and type). */
#line 688 "sla_syrfsx_extended.f"
	scopy_(n, &b[j * b_dim1 + 1], &c__1, &res[1], &c__1);
#line 689 "sla_syrfsx_extended.f"
	ssymv_(uplo, n, &c_b12, &a[a_offset], lda, &y[j * y_dim1 + 1], &c__1, 
		&c_b14, &res[1], &c__1, (ftnlen)1);
#line 691 "sla_syrfsx_extended.f"
	i__2 = *n;
#line 691 "sla_syrfsx_extended.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 692 "sla_syrfsx_extended.f"
	    ayb[i__] = (d__1 = b[i__ + j * b_dim1], abs(d__1));
#line 693 "sla_syrfsx_extended.f"
	}

/*     Compute abs(op(A_s))*abs(Y) + abs(B_s). */

#line 697 "sla_syrfsx_extended.f"
	sla_syamv__(&uplo2, n, &c_b14, &a[a_offset], lda, &y[j * y_dim1 + 1], 
		&c__1, &c_b14, &ayb[1], &c__1);
#line 700 "sla_syrfsx_extended.f"
	sla_lin_berr__(n, n, &c__1, &res[1], &ayb[1], &berr_out__[j]);

/*     End of loop for each RHS. */

#line 704 "sla_syrfsx_extended.f"
    }

#line 706 "sla_syrfsx_extended.f"
    return 0;
} /* sla_syrfsx_extended__ */

