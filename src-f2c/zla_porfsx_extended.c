#line 1 "zla_porfsx_extended.f"
/* zla_porfsx_extended.f -- translated by f2c (version 20100827).
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

#line 1 "zla_porfsx_extended.f"
/* Table of constant values */

static integer c__1 = 1;
static doublecomplex c_b11 = {-1.,0.};
static doublecomplex c_b12 = {1.,0.};
static doublereal c_b34 = 1.;

/* > \brief \b ZLA_PORFSX_EXTENDED improves the computed solution to a system of linear equations for symmetri
c or Hermitian positive-definite matrices by performing extra-precise iterative refinement and provide
s error bounds and backward error estimates fo */
/* r the solution. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZLA_PORFSX_EXTENDED + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zla_por
fsx_extended.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zla_por
fsx_extended.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zla_por
fsx_extended.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZLA_PORFSX_EXTENDED( PREC_TYPE, UPLO, N, NRHS, A, LDA, */
/*                                       AF, LDAF, COLEQU, C, B, LDB, Y, */
/*                                       LDY, BERR_OUT, N_NORMS, */
/*                                       ERR_BNDS_NORM, ERR_BNDS_COMP, RES, */
/*                                       AYB, DY, Y_TAIL, RCOND, ITHRESH, */
/*                                       RTHRESH, DZ_UB, IGNORE_CWISE, */
/*                                       INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, LDAF, LDB, LDY, N, NRHS, PREC_TYPE, */
/*      $                   N_NORMS, ITHRESH */
/*       CHARACTER          UPLO */
/*       LOGICAL            COLEQU, IGNORE_CWISE */
/*       DOUBLE PRECISION   RTHRESH, DZ_UB */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         A( LDA, * ), AF( LDAF, * ), B( LDB, * ), */
/*      $                   Y( LDY, * ), RES( * ), DY( * ), Y_TAIL( * ) */
/*       DOUBLE PRECISION   C( * ), AYB( * ), RCOND, BERR_OUT( * ), */
/*      $                   ERR_BNDS_NORM( NRHS, * ), */
/*      $                   ERR_BNDS_COMP( NRHS, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLA_PORFSX_EXTENDED improves the computed solution to a system of */
/* > linear equations by performing extra-precise iterative refinement */
/* > and provides error bounds and backward error estimates for the solution. */
/* > This subroutine is called by ZPORFSX to perform iterative refinement. */
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
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
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
/* >          AF is COMPLEX*16 array, dimension (LDAF,N) */
/* >     The triangular factor U or L from the Cholesky factorization */
/* >     A = U**T*U or A = L*L**T, as computed by ZPOTRF. */
/* > \endverbatim */
/* > */
/* > \param[in] LDAF */
/* > \verbatim */
/* >          LDAF is INTEGER */
/* >     The leading dimension of the array AF.  LDAF >= max(1,N). */
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
/* >          C is DOUBLE PRECISION array, dimension (N) */
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
/* >          B is COMPLEX*16 array, dimension (LDB,NRHS) */
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
/* >          Y is COMPLEX*16 array, dimension */
/* >                    (LDY,NRHS) */
/* >     On entry, the solution matrix X, as computed by ZPOTRS. */
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
/* >          BERR_OUT is DOUBLE PRECISION array, dimension (NRHS) */
/* >     On exit, BERR_OUT(j) contains the componentwise relative backward */
/* >     error for right-hand-side j from the formula */
/* >         max(i) ( abs(RES(i)) / ( abs(op(A_s))*abs(Y) + abs(B_s) )(i) ) */
/* >     where abs(Z) is the componentwise absolute value of the matrix */
/* >     or vector Z. This is computed by ZLA_LIN_BERR. */
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
/* >          ERR_BNDS_NORM is DOUBLE PRECISION array, dimension */
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
/* >          ERR_BNDS_COMP is DOUBLE PRECISION array, dimension */
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
/* >          RES is COMPLEX*16 array, dimension (N) */
/* >     Workspace to hold the intermediate residual. */
/* > \endverbatim */
/* > */
/* > \param[in] AYB */
/* > \verbatim */
/* >          AYB is DOUBLE PRECISION array, dimension (N) */
/* >     Workspace. */
/* > \endverbatim */
/* > */
/* > \param[in] DY */
/* > \verbatim */
/* >          DY is COMPLEX*16 PRECISION array, dimension (N) */
/* >     Workspace to hold the intermediate solution. */
/* > \endverbatim */
/* > */
/* > \param[in] Y_TAIL */
/* > \verbatim */
/* >          Y_TAIL is COMPLEX*16 array, dimension (N) */
/* >     Workspace to hold the trailing bits of the intermediate solution. */
/* > \endverbatim */
/* > */
/* > \param[in] RCOND */
/* > \verbatim */
/* >          RCOND is DOUBLE PRECISION */
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
/* >          RTHRESH is DOUBLE PRECISION */
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
/* >          DZ_UB is DOUBLE PRECISION */
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
/* >       < 0:  if INFO = -i, the ith argument to ZPOTRS had an illegal */
/* >             value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup complex16POcomputational */

/*  ===================================================================== */
/* Subroutine */ int zla_porfsx_extended__(integer *prec_type__, char *uplo, 
	integer *n, integer *nrhs, doublecomplex *a, integer *lda, 
	doublecomplex *af, integer *ldaf, logical *colequ, doublereal *c__, 
	doublecomplex *b, integer *ldb, doublecomplex *y, integer *ldy, 
	doublereal *berr_out__, integer *n_norms__, doublereal *
	err_bnds_norm__, doublereal *err_bnds_comp__, doublecomplex *res, 
	doublereal *ayb, doublecomplex *dy, doublecomplex *y_tail__, 
	doublereal *rcond, integer *ithresh, doublereal *rthresh, doublereal *
	dz_ub__, logical *ignore_cwise__, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, af_dim1, af_offset, b_dim1, b_offset, y_dim1, 
	    y_offset, err_bnds_norm_dim1, err_bnds_norm_offset, 
	    err_bnds_comp_dim1, err_bnds_comp_offset, i__1, i__2, i__3, i__4;
    doublereal d__1, d__2;

    /* Builtin functions */
    double d_imag(doublecomplex *);

    /* Local variables */
    static doublereal dxratmax, dzratmax;
    static integer i__, j;
    static logical incr_prec__;
    extern /* Subroutine */ int zla_heamv__(integer *, integer *, doublereal *
	    , doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublereal *, doublereal *, integer *);
    static doublereal prev_dz_z__, yk, final_dx_x__, final_dz_z__;
    extern /* Subroutine */ int zla_wwaddw__(integer *, doublecomplex *, 
	    doublecomplex *, doublecomplex *);
    static doublereal prevnormdx;
    static integer cnt;
    static doublereal dyk, eps, incr_thresh__, dx_x__, dz_z__, ymin;
    extern /* Subroutine */ int zla_lin_berr__(integer *, integer *, integer *
	    , doublecomplex *, doublereal *, doublereal *);
    static integer y_prec_state__;
    extern /* Subroutine */ int blas_zhemv_x__(integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, integer *)
	    ;
    static integer uplo2;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int blas_zhemv2_x__(integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, integer *);
    static doublereal dxrat, dzrat;
    extern /* Subroutine */ int zhemv_(char *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, ftnlen);
    static doublereal normx, normy;
    extern /* Subroutine */ int zcopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), zaxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal normdx;
    extern /* Subroutine */ int zpotrs_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     ftnlen);
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
/*     .. Statement Functions .. */
/*     .. */
/*     .. Statement Function Definitions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 473 "zla_porfsx_extended.f"
    /* Parameter adjustments */
#line 473 "zla_porfsx_extended.f"
    err_bnds_comp_dim1 = *nrhs;
#line 473 "zla_porfsx_extended.f"
    err_bnds_comp_offset = 1 + err_bnds_comp_dim1;
#line 473 "zla_porfsx_extended.f"
    err_bnds_comp__ -= err_bnds_comp_offset;
#line 473 "zla_porfsx_extended.f"
    err_bnds_norm_dim1 = *nrhs;
#line 473 "zla_porfsx_extended.f"
    err_bnds_norm_offset = 1 + err_bnds_norm_dim1;
#line 473 "zla_porfsx_extended.f"
    err_bnds_norm__ -= err_bnds_norm_offset;
#line 473 "zla_porfsx_extended.f"
    a_dim1 = *lda;
#line 473 "zla_porfsx_extended.f"
    a_offset = 1 + a_dim1;
#line 473 "zla_porfsx_extended.f"
    a -= a_offset;
#line 473 "zla_porfsx_extended.f"
    af_dim1 = *ldaf;
#line 473 "zla_porfsx_extended.f"
    af_offset = 1 + af_dim1;
#line 473 "zla_porfsx_extended.f"
    af -= af_offset;
#line 473 "zla_porfsx_extended.f"
    --c__;
#line 473 "zla_porfsx_extended.f"
    b_dim1 = *ldb;
#line 473 "zla_porfsx_extended.f"
    b_offset = 1 + b_dim1;
#line 473 "zla_porfsx_extended.f"
    b -= b_offset;
#line 473 "zla_porfsx_extended.f"
    y_dim1 = *ldy;
#line 473 "zla_porfsx_extended.f"
    y_offset = 1 + y_dim1;
#line 473 "zla_porfsx_extended.f"
    y -= y_offset;
#line 473 "zla_porfsx_extended.f"
    --berr_out__;
#line 473 "zla_porfsx_extended.f"
    --res;
#line 473 "zla_porfsx_extended.f"
    --ayb;
#line 473 "zla_porfsx_extended.f"
    --dy;
#line 473 "zla_porfsx_extended.f"
    --y_tail__;
#line 473 "zla_porfsx_extended.f"

#line 473 "zla_porfsx_extended.f"
    /* Function Body */
#line 473 "zla_porfsx_extended.f"
    if (*info != 0) {
#line 473 "zla_porfsx_extended.f"
	return 0;
#line 473 "zla_porfsx_extended.f"
    }
#line 474 "zla_porfsx_extended.f"
    eps = dlamch_("Epsilon", (ftnlen)7);
#line 475 "zla_porfsx_extended.f"
    hugeval = dlamch_("Overflow", (ftnlen)8);
/*     Force HUGEVAL to Inf */
#line 477 "zla_porfsx_extended.f"
    hugeval *= hugeval;
/*     Using HUGEVAL may lead to spurious underflows. */
#line 479 "zla_porfsx_extended.f"
    incr_thresh__ = (doublereal) (*n) * eps;
#line 481 "zla_porfsx_extended.f"
    if (lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 482 "zla_porfsx_extended.f"
	uplo2 = ilauplo_("L", (ftnlen)1);
#line 483 "zla_porfsx_extended.f"
    } else {
#line 484 "zla_porfsx_extended.f"
	uplo2 = ilauplo_("U", (ftnlen)1);
#line 485 "zla_porfsx_extended.f"
    }
#line 487 "zla_porfsx_extended.f"
    i__1 = *nrhs;
#line 487 "zla_porfsx_extended.f"
    for (j = 1; j <= i__1; ++j) {
#line 488 "zla_porfsx_extended.f"
	y_prec_state__ = 1;
#line 489 "zla_porfsx_extended.f"
	if (y_prec_state__ == 2) {
#line 490 "zla_porfsx_extended.f"
	    i__2 = *n;
#line 490 "zla_porfsx_extended.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 491 "zla_porfsx_extended.f"
		i__3 = i__;
#line 491 "zla_porfsx_extended.f"
		y_tail__[i__3].r = 0., y_tail__[i__3].i = 0.;
#line 492 "zla_porfsx_extended.f"
	    }
#line 493 "zla_porfsx_extended.f"
	}
#line 495 "zla_porfsx_extended.f"
	dxrat = 0.;
#line 496 "zla_porfsx_extended.f"
	dxratmax = 0.;
#line 497 "zla_porfsx_extended.f"
	dzrat = 0.;
#line 498 "zla_porfsx_extended.f"
	dzratmax = 0.;
#line 499 "zla_porfsx_extended.f"
	final_dx_x__ = hugeval;
#line 500 "zla_porfsx_extended.f"
	final_dz_z__ = hugeval;
#line 501 "zla_porfsx_extended.f"
	prevnormdx = hugeval;
#line 502 "zla_porfsx_extended.f"
	prev_dz_z__ = hugeval;
#line 503 "zla_porfsx_extended.f"
	dz_z__ = hugeval;
#line 504 "zla_porfsx_extended.f"
	dx_x__ = hugeval;
#line 506 "zla_porfsx_extended.f"
	x_state__ = 1;
#line 507 "zla_porfsx_extended.f"
	z_state__ = 0;
#line 508 "zla_porfsx_extended.f"
	incr_prec__ = FALSE_;
#line 510 "zla_porfsx_extended.f"
	i__2 = *ithresh;
#line 510 "zla_porfsx_extended.f"
	for (cnt = 1; cnt <= i__2; ++cnt) {

/*         Compute residual RES = B_s - op(A_s) * Y, */
/*             op(A) = A, A**T, or A**H depending on TRANS (and type). */

#line 515 "zla_porfsx_extended.f"
	    zcopy_(n, &b[j * b_dim1 + 1], &c__1, &res[1], &c__1);
#line 516 "zla_porfsx_extended.f"
	    if (y_prec_state__ == 0) {
#line 517 "zla_porfsx_extended.f"
		zhemv_(uplo, n, &c_b11, &a[a_offset], lda, &y[j * y_dim1 + 1],
			 &c__1, &c_b12, &res[1], &c__1, (ftnlen)1);
#line 519 "zla_porfsx_extended.f"
	    } else if (y_prec_state__ == 1) {
#line 520 "zla_porfsx_extended.f"
		blas_zhemv_x__(&uplo2, n, &c_b11, &a[a_offset], lda, &y[j * 
			y_dim1 + 1], &c__1, &c_b12, &res[1], &c__1, 
			prec_type__);
#line 522 "zla_porfsx_extended.f"
	    } else {
#line 523 "zla_porfsx_extended.f"
		blas_zhemv2_x__(&uplo2, n, &c_b11, &a[a_offset], lda, &y[j * 
			y_dim1 + 1], &y_tail__[1], &c__1, &c_b12, &res[1], &
			c__1, prec_type__);
#line 526 "zla_porfsx_extended.f"
	    }
/*         XXX: RES is no longer needed. */
#line 529 "zla_porfsx_extended.f"
	    zcopy_(n, &res[1], &c__1, &dy[1], &c__1);
#line 530 "zla_porfsx_extended.f"
	    zpotrs_(uplo, n, &c__1, &af[af_offset], ldaf, &dy[1], n, info, (
		    ftnlen)1);

/*         Calculate relative changes DX_X, DZ_Z and ratios DXRAT, DZRAT. */

#line 534 "zla_porfsx_extended.f"
	    normx = 0.;
#line 535 "zla_porfsx_extended.f"
	    normy = 0.;
#line 536 "zla_porfsx_extended.f"
	    normdx = 0.;
#line 537 "zla_porfsx_extended.f"
	    dz_z__ = 0.;
#line 538 "zla_porfsx_extended.f"
	    ymin = hugeval;
#line 540 "zla_porfsx_extended.f"
	    i__3 = *n;
#line 540 "zla_porfsx_extended.f"
	    for (i__ = 1; i__ <= i__3; ++i__) {
#line 541 "zla_porfsx_extended.f"
		i__4 = i__ + j * y_dim1;
#line 541 "zla_porfsx_extended.f"
		yk = (d__1 = y[i__4].r, abs(d__1)) + (d__2 = d_imag(&y[i__ + 
			j * y_dim1]), abs(d__2));
#line 542 "zla_porfsx_extended.f"
		i__4 = i__;
#line 542 "zla_porfsx_extended.f"
		dyk = (d__1 = dy[i__4].r, abs(d__1)) + (d__2 = d_imag(&dy[i__]
			), abs(d__2));
#line 544 "zla_porfsx_extended.f"
		if (yk != 0.) {
/* Computing MAX */
#line 545 "zla_porfsx_extended.f"
		    d__1 = dz_z__, d__2 = dyk / yk;
#line 545 "zla_porfsx_extended.f"
		    dz_z__ = max(d__1,d__2);
#line 546 "zla_porfsx_extended.f"
		} else if (dyk != 0.) {
#line 547 "zla_porfsx_extended.f"
		    dz_z__ = hugeval;
#line 548 "zla_porfsx_extended.f"
		}
#line 550 "zla_porfsx_extended.f"
		ymin = min(ymin,yk);
#line 552 "zla_porfsx_extended.f"
		normy = max(normy,yk);
#line 554 "zla_porfsx_extended.f"
		if (*colequ) {
/* Computing MAX */
#line 555 "zla_porfsx_extended.f"
		    d__1 = normx, d__2 = yk * c__[i__];
#line 555 "zla_porfsx_extended.f"
		    normx = max(d__1,d__2);
/* Computing MAX */
#line 556 "zla_porfsx_extended.f"
		    d__1 = normdx, d__2 = dyk * c__[i__];
#line 556 "zla_porfsx_extended.f"
		    normdx = max(d__1,d__2);
#line 557 "zla_porfsx_extended.f"
		} else {
#line 558 "zla_porfsx_extended.f"
		    normx = normy;
#line 559 "zla_porfsx_extended.f"
		    normdx = max(normdx,dyk);
#line 560 "zla_porfsx_extended.f"
		}
#line 561 "zla_porfsx_extended.f"
	    }
#line 563 "zla_porfsx_extended.f"
	    if (normx != 0.) {
#line 564 "zla_porfsx_extended.f"
		dx_x__ = normdx / normx;
#line 565 "zla_porfsx_extended.f"
	    } else if (normdx == 0.) {
#line 566 "zla_porfsx_extended.f"
		dx_x__ = 0.;
#line 567 "zla_porfsx_extended.f"
	    } else {
#line 568 "zla_porfsx_extended.f"
		dx_x__ = hugeval;
#line 569 "zla_porfsx_extended.f"
	    }
#line 571 "zla_porfsx_extended.f"
	    dxrat = normdx / prevnormdx;
#line 572 "zla_porfsx_extended.f"
	    dzrat = dz_z__ / prev_dz_z__;

/*         Check termination criteria. */

#line 576 "zla_porfsx_extended.f"
	    if (ymin * *rcond < incr_thresh__ * normy && y_prec_state__ < 2) {
#line 576 "zla_porfsx_extended.f"
		incr_prec__ = TRUE_;
#line 576 "zla_porfsx_extended.f"
	    }
#line 580 "zla_porfsx_extended.f"
	    if (x_state__ == 3 && dxrat <= *rthresh) {
#line 580 "zla_porfsx_extended.f"
		x_state__ = 1;
#line 580 "zla_porfsx_extended.f"
	    }
#line 582 "zla_porfsx_extended.f"
	    if (x_state__ == 1) {
#line 583 "zla_porfsx_extended.f"
		if (dx_x__ <= eps) {
#line 584 "zla_porfsx_extended.f"
		    x_state__ = 2;
#line 585 "zla_porfsx_extended.f"
		} else if (dxrat > *rthresh) {
#line 586 "zla_porfsx_extended.f"
		    if (y_prec_state__ != 2) {
#line 587 "zla_porfsx_extended.f"
			incr_prec__ = TRUE_;
#line 588 "zla_porfsx_extended.f"
		    } else {
#line 589 "zla_porfsx_extended.f"
			x_state__ = 3;
#line 590 "zla_porfsx_extended.f"
		    }
#line 591 "zla_porfsx_extended.f"
		} else {
#line 592 "zla_porfsx_extended.f"
		    if (dxrat > dxratmax) {
#line 592 "zla_porfsx_extended.f"
			dxratmax = dxrat;
#line 592 "zla_porfsx_extended.f"
		    }
#line 593 "zla_porfsx_extended.f"
		}
#line 594 "zla_porfsx_extended.f"
		if (x_state__ > 1) {
#line 594 "zla_porfsx_extended.f"
		    final_dx_x__ = dx_x__;
#line 594 "zla_porfsx_extended.f"
		}
#line 595 "zla_porfsx_extended.f"
	    }
#line 597 "zla_porfsx_extended.f"
	    if (z_state__ == 0 && dz_z__ <= *dz_ub__) {
#line 597 "zla_porfsx_extended.f"
		z_state__ = 1;
#line 597 "zla_porfsx_extended.f"
	    }
#line 599 "zla_porfsx_extended.f"
	    if (z_state__ == 3 && dzrat <= *rthresh) {
#line 599 "zla_porfsx_extended.f"
		z_state__ = 1;
#line 599 "zla_porfsx_extended.f"
	    }
#line 601 "zla_porfsx_extended.f"
	    if (z_state__ == 1) {
#line 602 "zla_porfsx_extended.f"
		if (dz_z__ <= eps) {
#line 603 "zla_porfsx_extended.f"
		    z_state__ = 2;
#line 604 "zla_porfsx_extended.f"
		} else if (dz_z__ > *dz_ub__) {
#line 605 "zla_porfsx_extended.f"
		    z_state__ = 0;
#line 606 "zla_porfsx_extended.f"
		    dzratmax = 0.;
#line 607 "zla_porfsx_extended.f"
		    final_dz_z__ = hugeval;
#line 608 "zla_porfsx_extended.f"
		} else if (dzrat > *rthresh) {
#line 609 "zla_porfsx_extended.f"
		    if (y_prec_state__ != 2) {
#line 610 "zla_porfsx_extended.f"
			incr_prec__ = TRUE_;
#line 611 "zla_porfsx_extended.f"
		    } else {
#line 612 "zla_porfsx_extended.f"
			z_state__ = 3;
#line 613 "zla_porfsx_extended.f"
		    }
#line 614 "zla_porfsx_extended.f"
		} else {
#line 615 "zla_porfsx_extended.f"
		    if (dzrat > dzratmax) {
#line 615 "zla_porfsx_extended.f"
			dzratmax = dzrat;
#line 615 "zla_porfsx_extended.f"
		    }
#line 616 "zla_porfsx_extended.f"
		}
#line 617 "zla_porfsx_extended.f"
		if (z_state__ > 1) {
#line 617 "zla_porfsx_extended.f"
		    final_dz_z__ = dz_z__;
#line 617 "zla_porfsx_extended.f"
		}
#line 618 "zla_porfsx_extended.f"
	    }
#line 620 "zla_porfsx_extended.f"
	    if (x_state__ != 1 && (*ignore_cwise__ || z_state__ != 1)) {
#line 620 "zla_porfsx_extended.f"
		goto L666;
#line 620 "zla_porfsx_extended.f"
	    }
#line 624 "zla_porfsx_extended.f"
	    if (incr_prec__) {
#line 625 "zla_porfsx_extended.f"
		incr_prec__ = FALSE_;
#line 626 "zla_porfsx_extended.f"
		++y_prec_state__;
#line 627 "zla_porfsx_extended.f"
		i__3 = *n;
#line 627 "zla_porfsx_extended.f"
		for (i__ = 1; i__ <= i__3; ++i__) {
#line 628 "zla_porfsx_extended.f"
		    i__4 = i__;
#line 628 "zla_porfsx_extended.f"
		    y_tail__[i__4].r = 0., y_tail__[i__4].i = 0.;
#line 629 "zla_porfsx_extended.f"
		}
#line 630 "zla_porfsx_extended.f"
	    }
#line 632 "zla_porfsx_extended.f"
	    prevnormdx = normdx;
#line 633 "zla_porfsx_extended.f"
	    prev_dz_z__ = dz_z__;

/*           Update soluton. */

#line 637 "zla_porfsx_extended.f"
	    if (y_prec_state__ < 2) {
#line 638 "zla_porfsx_extended.f"
		zaxpy_(n, &c_b12, &dy[1], &c__1, &y[j * y_dim1 + 1], &c__1);
#line 639 "zla_porfsx_extended.f"
	    } else {
#line 640 "zla_porfsx_extended.f"
		zla_wwaddw__(n, &y[j * y_dim1 + 1], &y_tail__[1], &dy[1]);
#line 641 "zla_porfsx_extended.f"
	    }
#line 643 "zla_porfsx_extended.f"
	}
/*        Target of "IF (Z_STOP .AND. X_STOP)".  Sun's f77 won't EXIT. */
#line 645 "zla_porfsx_extended.f"
L666:

/*     Set final_* when cnt hits ithresh. */

#line 649 "zla_porfsx_extended.f"
	if (x_state__ == 1) {
#line 649 "zla_porfsx_extended.f"
	    final_dx_x__ = dx_x__;
#line 649 "zla_porfsx_extended.f"
	}
#line 650 "zla_porfsx_extended.f"
	if (z_state__ == 1) {
#line 650 "zla_porfsx_extended.f"
	    final_dz_z__ = dz_z__;
#line 650 "zla_porfsx_extended.f"
	}

/*     Compute error bounds. */

#line 654 "zla_porfsx_extended.f"
	if (*n_norms__ >= 1) {
#line 655 "zla_porfsx_extended.f"
	    err_bnds_norm__[j + (err_bnds_norm_dim1 << 1)] = final_dx_x__ / (
		    1 - dxratmax);
#line 657 "zla_porfsx_extended.f"
	}
#line 658 "zla_porfsx_extended.f"
	if (*n_norms__ >= 2) {
#line 659 "zla_porfsx_extended.f"
	    err_bnds_comp__[j + (err_bnds_comp_dim1 << 1)] = final_dz_z__ / (
		    1 - dzratmax);
#line 661 "zla_porfsx_extended.f"
	}

/*     Compute componentwise relative backward error from formula */
/*         max(i) ( abs(R(i)) / ( abs(op(A_s))*abs(Y) + abs(B_s) )(i) ) */
/*     where abs(Z) is the componentwise absolute value of the matrix */
/*     or vector Z. */

/*        Compute residual RES = B_s - op(A_s) * Y, */
/*            op(A) = A, A**T, or A**H depending on TRANS (and type). */

#line 671 "zla_porfsx_extended.f"
	zcopy_(n, &b[j * b_dim1 + 1], &c__1, &res[1], &c__1);
#line 672 "zla_porfsx_extended.f"
	zhemv_(uplo, n, &c_b11, &a[a_offset], lda, &y[j * y_dim1 + 1], &c__1, 
		&c_b12, &res[1], &c__1, (ftnlen)1);
#line 675 "zla_porfsx_extended.f"
	i__2 = *n;
#line 675 "zla_porfsx_extended.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 676 "zla_porfsx_extended.f"
	    i__3 = i__ + j * b_dim1;
#line 676 "zla_porfsx_extended.f"
	    ayb[i__] = (d__1 = b[i__3].r, abs(d__1)) + (d__2 = d_imag(&b[i__ 
		    + j * b_dim1]), abs(d__2));
#line 677 "zla_porfsx_extended.f"
	}

/*     Compute abs(op(A_s))*abs(Y) + abs(B_s). */

#line 681 "zla_porfsx_extended.f"
	zla_heamv__(&uplo2, n, &c_b34, &a[a_offset], lda, &y[j * y_dim1 + 1], 
		&c__1, &c_b34, &ayb[1], &c__1);
#line 684 "zla_porfsx_extended.f"
	zla_lin_berr__(n, n, &c__1, &res[1], &ayb[1], &berr_out__[j]);

/*     End of loop for each RHS. */

#line 688 "zla_porfsx_extended.f"
    }

#line 690 "zla_porfsx_extended.f"
    return 0;
} /* zla_porfsx_extended__ */

