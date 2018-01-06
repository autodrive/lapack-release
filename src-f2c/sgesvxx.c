#line 1 "sgesvxx.f"
/* sgesvxx.f -- translated by f2c (version 20100827).
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

#line 1 "sgesvxx.f"
/* > \brief <b> SGESVXX computes the solution to system of linear equations A * X = B for GE matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SGESVXX + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgesvxx
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgesvxx
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgesvxx
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SGESVXX( FACT, TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV, */
/*                           EQUED, R, C, B, LDB, X, LDX, RCOND, RPVGRW, */
/*                           BERR, N_ERR_BNDS, ERR_BNDS_NORM, */
/*                           ERR_BNDS_COMP, NPARAMS, PARAMS, WORK, IWORK, */
/*                           INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          EQUED, FACT, TRANS */
/*       INTEGER            INFO, LDA, LDAF, LDB, LDX, N, NRHS, NPARAMS, */
/*      $                   N_ERR_BNDS */
/*       REAL               RCOND, RPVGRW */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ), IWORK( * ) */
/*       REAL               A( LDA, * ), AF( LDAF, * ), B( LDB, * ), */
/*      $                   X( LDX , * ),WORK( * ) */
/*       REAL               R( * ), C( * ), PARAMS( * ), BERR( * ), */
/*      $                   ERR_BNDS_NORM( NRHS, * ), */
/*      $                   ERR_BNDS_COMP( NRHS, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >    SGESVXX uses the LU factorization to compute the solution to a */
/* >    real system of linear equations  A * X = B,  where A is an */
/* >    N-by-N matrix and X and B are N-by-NRHS matrices. */
/* > */
/* >    If requested, both normwise and maximum componentwise error bounds */
/* >    are returned. SGESVXX will return a solution with a tiny */
/* >    guaranteed error (O(eps) where eps is the working machine */
/* >    precision) unless the matrix is very ill-conditioned, in which */
/* >    case a warning is returned. Relevant condition numbers also are */
/* >    calculated and returned. */
/* > */
/* >    SGESVXX accepts user-provided factorizations and equilibration */
/* >    factors; see the definitions of the FACT and EQUED options. */
/* >    Solving with refinement and using a factorization from a previous */
/* >    SGESVXX call will also produce a solution with either O(eps) */
/* >    errors or warnings, but we cannot make that claim for general */
/* >    user-provided factorizations and equilibration factors if they */
/* >    differ from what SGESVXX would itself produce. */
/* > \endverbatim */

/* > \par Description: */
/*  ================= */
/* > */
/* > \verbatim */
/* > */
/* >    The following steps are performed: */
/* > */
/* >    1. If FACT = 'E', real scaling factors are computed to equilibrate */
/* >    the system: */
/* > */
/* >      TRANS = 'N':  diag(R)*A*diag(C)     *inv(diag(C))*X = diag(R)*B */
/* >      TRANS = 'T': (diag(R)*A*diag(C))**T *inv(diag(R))*X = diag(C)*B */
/* >      TRANS = 'C': (diag(R)*A*diag(C))**H *inv(diag(R))*X = diag(C)*B */
/* > */
/* >    Whether or not the system will be equilibrated depends on the */
/* >    scaling of the matrix A, but if equilibration is used, A is */
/* >    overwritten by diag(R)*A*diag(C) and B by diag(R)*B (if TRANS='N') */
/* >    or diag(C)*B (if TRANS = 'T' or 'C'). */
/* > */
/* >    2. If FACT = 'N' or 'E', the LU decomposition is used to factor */
/* >    the matrix A (after equilibration if FACT = 'E') as */
/* > */
/* >      A = P * L * U, */
/* > */
/* >    where P is a permutation matrix, L is a unit lower triangular */
/* >    matrix, and U is upper triangular. */
/* > */
/* >    3. If some U(i,i)=0, so that U is exactly singular, then the */
/* >    routine returns with INFO = i. Otherwise, the factored form of A */
/* >    is used to estimate the condition number of the matrix A (see */
/* >    argument RCOND). If the reciprocal of the condition number is less */
/* >    than machine precision, the routine still goes on to solve for X */
/* >    and compute error bounds as described below. */
/* > */
/* >    4. The system of equations is solved for X using the factored form */
/* >    of A. */
/* > */
/* >    5. By default (unless PARAMS(LA_LINRX_ITREF_I) is set to zero), */
/* >    the routine will use iterative refinement to try to get a small */
/* >    error and error bounds.  Refinement calculates the residual to at */
/* >    least twice the working precision. */
/* > */
/* >    6. If equilibration was used, the matrix X is premultiplied by */
/* >    diag(C) (if TRANS = 'N') or diag(R) (if TRANS = 'T' or 'C') so */
/* >    that it solves the original system before equilibration. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \verbatim */
/* >     Some optional parameters are bundled in the PARAMS array.  These */
/* >     settings determine how refinement is performed, but often the */
/* >     defaults are acceptable.  If the defaults are acceptable, users */
/* >     can pass NPARAMS = 0 which prevents the source code from accessing */
/* >     the PARAMS argument. */
/* > \endverbatim */
/* > */
/* > \param[in] FACT */
/* > \verbatim */
/* >          FACT is CHARACTER*1 */
/* >     Specifies whether or not the factored form of the matrix A is */
/* >     supplied on entry, and if not, whether the matrix A should be */
/* >     equilibrated before it is factored. */
/* >       = 'F':  On entry, AF and IPIV contain the factored form of A. */
/* >               If EQUED is not 'N', the matrix A has been */
/* >               equilibrated with scaling factors given by R and C. */
/* >               A, AF, and IPIV are not modified. */
/* >       = 'N':  The matrix A will be copied to AF and factored. */
/* >       = 'E':  The matrix A will be equilibrated if necessary, then */
/* >               copied to AF and factored. */
/* > \endverbatim */
/* > */
/* > \param[in] TRANS */
/* > \verbatim */
/* >          TRANS is CHARACTER*1 */
/* >     Specifies the form of the system of equations: */
/* >       = 'N':  A * X = B     (No transpose) */
/* >       = 'T':  A**T * X = B  (Transpose) */
/* >       = 'C':  A**H * X = B  (Conjugate Transpose = Transpose) */
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
/* >     The number of right hand sides, i.e., the number of columns */
/* >     of the matrices B and X.  NRHS >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is REAL array, dimension (LDA,N) */
/* >     On entry, the N-by-N matrix A.  If FACT = 'F' and EQUED is */
/* >     not 'N', then A must have been equilibrated by the scaling */
/* >     factors in R and/or C.  A is not modified if FACT = 'F' or */
/* >     'N', or if FACT = 'E' and EQUED = 'N' on exit. */
/* > */
/* >     On exit, if EQUED .ne. 'N', A is scaled as follows: */
/* >     EQUED = 'R':  A := diag(R) * A */
/* >     EQUED = 'C':  A := A * diag(C) */
/* >     EQUED = 'B':  A := diag(R) * A * diag(C). */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >     The leading dimension of the array A.  LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] AF */
/* > \verbatim */
/* >          AF is REAL array, dimension (LDAF,N) */
/* >     If FACT = 'F', then AF is an input argument and on entry */
/* >     contains the factors L and U from the factorization */
/* >     A = P*L*U as computed by SGETRF.  If EQUED .ne. 'N', then */
/* >     AF is the factored form of the equilibrated matrix A. */
/* > */
/* >     If FACT = 'N', then AF is an output argument and on exit */
/* >     returns the factors L and U from the factorization A = P*L*U */
/* >     of the original matrix A. */
/* > */
/* >     If FACT = 'E', then AF is an output argument and on exit */
/* >     returns the factors L and U from the factorization A = P*L*U */
/* >     of the equilibrated matrix A (see the description of A for */
/* >     the form of the equilibrated matrix). */
/* > \endverbatim */
/* > */
/* > \param[in] LDAF */
/* > \verbatim */
/* >          LDAF is INTEGER */
/* >     The leading dimension of the array AF.  LDAF >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] IPIV */
/* > \verbatim */
/* >          IPIV is INTEGER array, dimension (N) */
/* >     If FACT = 'F', then IPIV is an input argument and on entry */
/* >     contains the pivot indices from the factorization A = P*L*U */
/* >     as computed by SGETRF; row i of the matrix was interchanged */
/* >     with row IPIV(i). */
/* > */
/* >     If FACT = 'N', then IPIV is an output argument and on exit */
/* >     contains the pivot indices from the factorization A = P*L*U */
/* >     of the original matrix A. */
/* > */
/* >     If FACT = 'E', then IPIV is an output argument and on exit */
/* >     contains the pivot indices from the factorization A = P*L*U */
/* >     of the equilibrated matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in,out] EQUED */
/* > \verbatim */
/* >          EQUED is CHARACTER*1 */
/* >     Specifies the form of equilibration that was done. */
/* >       = 'N':  No equilibration (always true if FACT = 'N'). */
/* >       = 'R':  Row equilibration, i.e., A has been premultiplied by */
/* >               diag(R). */
/* >       = 'C':  Column equilibration, i.e., A has been postmultiplied */
/* >               by diag(C). */
/* >       = 'B':  Both row and column equilibration, i.e., A has been */
/* >               replaced by diag(R) * A * diag(C). */
/* >     EQUED is an input argument if FACT = 'F'; otherwise, it is an */
/* >     output argument. */
/* > \endverbatim */
/* > */
/* > \param[in,out] R */
/* > \verbatim */
/* >          R is REAL array, dimension (N) */
/* >     The row scale factors for A.  If EQUED = 'R' or 'B', A is */
/* >     multiplied on the left by diag(R); if EQUED = 'N' or 'C', R */
/* >     is not accessed.  R is an input argument if FACT = 'F'; */
/* >     otherwise, R is an output argument.  If FACT = 'F' and */
/* >     EQUED = 'R' or 'B', each element of R must be positive. */
/* >     If R is output, each element of R is a power of the radix. */
/* >     If R is input, each element of R should be a power of the radix */
/* >     to ensure a reliable solution and error estimates. Scaling by */
/* >     powers of the radix does not cause rounding errors unless the */
/* >     result underflows or overflows. Rounding errors during scaling */
/* >     lead to refining with a matrix that is not equivalent to the */
/* >     input matrix, producing error estimates that may not be */
/* >     reliable. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* >          C is REAL array, dimension (N) */
/* >     The column scale factors for A.  If EQUED = 'C' or 'B', A is */
/* >     multiplied on the right by diag(C); if EQUED = 'N' or 'R', C */
/* >     is not accessed.  C is an input argument if FACT = 'F'; */
/* >     otherwise, C is an output argument.  If FACT = 'F' and */
/* >     EQUED = 'C' or 'B', each element of C must be positive. */
/* >     If C is output, each element of C is a power of the radix. */
/* >     If C is input, each element of C should be a power of the radix */
/* >     to ensure a reliable solution and error estimates. Scaling by */
/* >     powers of the radix does not cause rounding errors unless the */
/* >     result underflows or overflows. Rounding errors during scaling */
/* >     lead to refining with a matrix that is not equivalent to the */
/* >     input matrix, producing error estimates that may not be */
/* >     reliable. */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* >          B is REAL array, dimension (LDB,NRHS) */
/* >     On entry, the N-by-NRHS right hand side matrix B. */
/* >     On exit, */
/* >     if EQUED = 'N', B is not modified; */
/* >     if TRANS = 'N' and EQUED = 'R' or 'B', B is overwritten by */
/* >        diag(R)*B; */
/* >     if TRANS = 'T' or 'C' and EQUED = 'C' or 'B', B is */
/* >        overwritten by diag(C)*B. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >     The leading dimension of the array B.  LDB >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] X */
/* > \verbatim */
/* >          X is REAL array, dimension (LDX,NRHS) */
/* >     If INFO = 0, the N-by-NRHS solution matrix X to the original */
/* >     system of equations.  Note that A and B are modified on exit */
/* >     if EQUED .ne. 'N', and the solution to the equilibrated system is */
/* >     inv(diag(C))*X if TRANS = 'N' and EQUED = 'C' or 'B', or */
/* >     inv(diag(R))*X if TRANS = 'T' or 'C' and EQUED = 'R' or 'B'. */
/* > \endverbatim */
/* > */
/* > \param[in] LDX */
/* > \verbatim */
/* >          LDX is INTEGER */
/* >     The leading dimension of the array X.  LDX >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] RCOND */
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
/* > \param[out] RPVGRW */
/* > \verbatim */
/* >          RPVGRW is REAL */
/* >     Reciprocal pivot growth.  On exit, this contains the reciprocal */
/* >     pivot growth factor norm(A)/norm(U). The "max absolute element" */
/* >     norm is used.  If this is much less than 1, then the stability of */
/* >     the LU factorization of the (equilibrated) matrix A could be poor. */
/* >     This also means that the solution X, estimated condition numbers, */
/* >     and error bounds could be unreliable. If factorization fails with */
/* >     0<INFO<=N, then this contains the reciprocal pivot growth factor */
/* >     for the leading INFO columns of A.  In SGESVX, this quantity is */
/* >     returned in WORK(1). */
/* > \endverbatim */
/* > */
/* > \param[out] BERR */
/* > \verbatim */
/* >          BERR is REAL array, dimension (NRHS) */
/* >     Componentwise relative backward error.  This is the */
/* >     componentwise relative backward error of each solution vector X(j) */
/* >     (i.e., the smallest relative change in any element of A or B that */
/* >     makes X(j) an exact solution). */
/* > \endverbatim */
/* > */
/* > \param[in] N_ERR_BNDS */
/* > \verbatim */
/* >          N_ERR_BNDS is INTEGER */
/* >     Number of error bounds to return for each right hand side */
/* >     and each type (normwise or componentwise).  See ERR_BNDS_NORM and */
/* >     ERR_BNDS_COMP below. */
/* > \endverbatim */
/* > */
/* > \param[out] ERR_BNDS_NORM */
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
/* >     See Lapack Working Note 165 for further details and extra */
/* >     cautions. */
/* > \endverbatim */
/* > */
/* > \param[out] ERR_BNDS_COMP */
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
/* >     See Lapack Working Note 165 for further details and extra */
/* >     cautions. */
/* > \endverbatim */
/* > */
/* > \param[in] NPARAMS */
/* > \verbatim */
/* >          NPARAMS is INTEGER */
/* >     Specifies the number of parameters set in PARAMS.  If .LE. 0, the */
/* >     PARAMS array is never referenced and default values are used. */
/* > \endverbatim */
/* > */
/* > \param[in,out] PARAMS */
/* > \verbatim */
/* >          PARAMS is REAL array, dimension NPARAMS */
/* >     Specifies algorithm parameters.  If an entry is .LT. 0.0, then */
/* >     that entry will be filled with default value used for that */
/* >     parameter.  Only positions up to NPARAMS are accessed; defaults */
/* >     are used for higher-numbered parameters. */
/* > */
/* >       PARAMS(LA_LINRX_ITREF_I = 1) : Whether to perform iterative */
/* >            refinement or not. */
/* >         Default: 1.0 */
/* >            = 0.0 : No refinement is performed, and no error bounds are */
/* >                    computed. */
/* >            = 1.0 : Use the double-precision refinement algorithm, */
/* >                    possibly with doubled-single computations if the */
/* >                    compilation environment does not support DOUBLE */
/* >                    PRECISION. */
/* >              (other values are reserved for future use) */
/* > */
/* >       PARAMS(LA_LINRX_ITHRESH_I = 2) : Maximum number of residual */
/* >            computations allowed for refinement. */
/* >         Default: 10 */
/* >         Aggressive: Set to 100 to permit convergence using approximate */
/* >                     factorizations or factorizations other than LU. If */
/* >                     the factorization uses a technique other than */
/* >                     Gaussian elimination, the guarantees in */
/* >                     err_bnds_norm and err_bnds_comp may no longer be */
/* >                     trustworthy. */
/* > */
/* >       PARAMS(LA_LINRX_CWISE_I = 3) : Flag determining if the code */
/* >            will attempt to find a solution with small componentwise */
/* >            relative error in the double-precision algorithm.  Positive */
/* >            is true, 0.0 is false. */
/* >         Default: 1.0 (attempt componentwise convergence) */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is REAL array, dimension (4*N) */
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
/* >       = 0:  Successful exit. The solution to every right-hand side is */
/* >         guaranteed. */
/* >       < 0:  If INFO = -i, the i-th argument had an illegal value */
/* >       > 0 and <= N:  U(INFO,INFO) is exactly zero.  The factorization */
/* >         has been completed, but the factor U is exactly singular, so */
/* >         the solution and error bounds could not be computed. RCOND = 0 */
/* >         is returned. */
/* >       = N+J: The solution corresponding to the Jth right-hand side is */
/* >         not guaranteed. The solutions corresponding to other right- */
/* >         hand sides K with K > J may not be guaranteed as well, but */
/* >         only the first such right-hand side is reported. If a small */
/* >         componentwise error is not requested (PARAMS(3) = 0.0) then */
/* >         the Jth right-hand side is the first with a normwise error */
/* >         bound that is not guaranteed (the smallest J such */
/* >         that ERR_BNDS_NORM(J,1) = 0.0). By default (PARAMS(3) = 1.0) */
/* >         the Jth right-hand side is the first with either a normwise or */
/* >         componentwise error bound that is not guaranteed (the smallest */
/* >         J such that either ERR_BNDS_NORM(J,1) = 0.0 or */
/* >         ERR_BNDS_COMP(J,1) = 0.0). See the definition of */
/* >         ERR_BNDS_NORM(:,1) and ERR_BNDS_COMP(:,1). To get information */
/* >         about all of the right-hand sides check ERR_BNDS_NORM or */
/* >         ERR_BNDS_COMP. */
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
/* Subroutine */ int sgesvxx_(char *fact, char *trans, integer *n, integer *
	nrhs, doublereal *a, integer *lda, doublereal *af, integer *ldaf, 
	integer *ipiv, char *equed, doublereal *r__, doublereal *c__, 
	doublereal *b, integer *ldb, doublereal *x, integer *ldx, doublereal *
	rcond, doublereal *rpvgrw, doublereal *berr, integer *n_err_bnds__, 
	doublereal *err_bnds_norm__, doublereal *err_bnds_comp__, integer *
	nparams, doublereal *params, doublereal *work, integer *iwork, 
	integer *info, ftnlen fact_len, ftnlen trans_len, ftnlen equed_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, af_dim1, af_offset, b_dim1, b_offset, x_dim1, 
	    x_offset, err_bnds_norm_dim1, err_bnds_norm_offset, 
	    err_bnds_comp_dim1, err_bnds_comp_offset, i__1;
    doublereal d__1, d__2;

    /* Local variables */
    static integer j;
    static doublereal amax;
    extern doublereal sla_gerpvgrw__(integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal rcmin, rcmax;
    static logical equil;
    static doublereal colcnd;
    extern doublereal slamch_(char *, ftnlen);
    static logical nofact;
    extern /* Subroutine */ int slaqge_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, char *, ftnlen), xerbla_(char *, integer *, ftnlen)
	    ;
    static doublereal bignum;
    static integer infequ;
    static logical colequ;
    extern /* Subroutine */ int sgetrf_(integer *, integer *, doublereal *, 
	    integer *, integer *, integer *), slacpy_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen);
    static doublereal rowcnd;
    static logical notran;
    extern /* Subroutine */ int sgetrs_(char *, integer *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, ftnlen);
    static doublereal smlnum;
    static logical rowequ;
    extern /* Subroutine */ int slascl2_(integer *, integer *, doublereal *, 
	    doublereal *, integer *), sgeequb_(integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, integer *), sgerfsx_(char *, char *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, ftnlen, ftnlen);


/*  -- LAPACK driver routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     April 2012 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ================================================================== */

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

#line 600 "sgesvxx.f"
    /* Parameter adjustments */
#line 600 "sgesvxx.f"
    err_bnds_comp_dim1 = *nrhs;
#line 600 "sgesvxx.f"
    err_bnds_comp_offset = 1 + err_bnds_comp_dim1;
#line 600 "sgesvxx.f"
    err_bnds_comp__ -= err_bnds_comp_offset;
#line 600 "sgesvxx.f"
    err_bnds_norm_dim1 = *nrhs;
#line 600 "sgesvxx.f"
    err_bnds_norm_offset = 1 + err_bnds_norm_dim1;
#line 600 "sgesvxx.f"
    err_bnds_norm__ -= err_bnds_norm_offset;
#line 600 "sgesvxx.f"
    a_dim1 = *lda;
#line 600 "sgesvxx.f"
    a_offset = 1 + a_dim1;
#line 600 "sgesvxx.f"
    a -= a_offset;
#line 600 "sgesvxx.f"
    af_dim1 = *ldaf;
#line 600 "sgesvxx.f"
    af_offset = 1 + af_dim1;
#line 600 "sgesvxx.f"
    af -= af_offset;
#line 600 "sgesvxx.f"
    --ipiv;
#line 600 "sgesvxx.f"
    --r__;
#line 600 "sgesvxx.f"
    --c__;
#line 600 "sgesvxx.f"
    b_dim1 = *ldb;
#line 600 "sgesvxx.f"
    b_offset = 1 + b_dim1;
#line 600 "sgesvxx.f"
    b -= b_offset;
#line 600 "sgesvxx.f"
    x_dim1 = *ldx;
#line 600 "sgesvxx.f"
    x_offset = 1 + x_dim1;
#line 600 "sgesvxx.f"
    x -= x_offset;
#line 600 "sgesvxx.f"
    --berr;
#line 600 "sgesvxx.f"
    --params;
#line 600 "sgesvxx.f"
    --work;
#line 600 "sgesvxx.f"
    --iwork;
#line 600 "sgesvxx.f"

#line 600 "sgesvxx.f"
    /* Function Body */
#line 600 "sgesvxx.f"
    *info = 0;
#line 601 "sgesvxx.f"
    nofact = lsame_(fact, "N", (ftnlen)1, (ftnlen)1);
#line 602 "sgesvxx.f"
    equil = lsame_(fact, "E", (ftnlen)1, (ftnlen)1);
#line 603 "sgesvxx.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 604 "sgesvxx.f"
    smlnum = slamch_("Safe minimum", (ftnlen)12);
#line 605 "sgesvxx.f"
    bignum = 1. / smlnum;
#line 606 "sgesvxx.f"
    if (nofact || equil) {
#line 607 "sgesvxx.f"
	*(unsigned char *)equed = 'N';
#line 608 "sgesvxx.f"
	rowequ = FALSE_;
#line 609 "sgesvxx.f"
	colequ = FALSE_;
#line 610 "sgesvxx.f"
    } else {
#line 611 "sgesvxx.f"
	rowequ = lsame_(equed, "R", (ftnlen)1, (ftnlen)1) || lsame_(equed, 
		"B", (ftnlen)1, (ftnlen)1);
#line 612 "sgesvxx.f"
	colequ = lsame_(equed, "C", (ftnlen)1, (ftnlen)1) || lsame_(equed, 
		"B", (ftnlen)1, (ftnlen)1);
#line 613 "sgesvxx.f"
    }

/*     Default is failure.  If an input parameter is wrong or */
/*     factorization fails, make everything look horrible.  Only the */
/*     pivot growth is set here, the rest is initialized in SGERFSX. */

#line 619 "sgesvxx.f"
    *rpvgrw = 0.;

/*     Test the input parameters.  PARAMS is not tested until SGERFSX. */

#line 623 "sgesvxx.f"
    if (! nofact && ! equil && ! lsame_(fact, "F", (ftnlen)1, (ftnlen)1)) {
#line 625 "sgesvxx.f"
	*info = -1;
#line 626 "sgesvxx.f"
    } else if (! notran && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1) && ! 
	    lsame_(trans, "C", (ftnlen)1, (ftnlen)1)) {
#line 628 "sgesvxx.f"
	*info = -2;
#line 629 "sgesvxx.f"
    } else if (*n < 0) {
#line 630 "sgesvxx.f"
	*info = -3;
#line 631 "sgesvxx.f"
    } else if (*nrhs < 0) {
#line 632 "sgesvxx.f"
	*info = -4;
#line 633 "sgesvxx.f"
    } else if (*lda < max(1,*n)) {
#line 634 "sgesvxx.f"
	*info = -6;
#line 635 "sgesvxx.f"
    } else if (*ldaf < max(1,*n)) {
#line 636 "sgesvxx.f"
	*info = -8;
#line 637 "sgesvxx.f"
    } else if (lsame_(fact, "F", (ftnlen)1, (ftnlen)1) && ! (rowequ || colequ 
	    || lsame_(equed, "N", (ftnlen)1, (ftnlen)1))) {
#line 639 "sgesvxx.f"
	*info = -10;
#line 640 "sgesvxx.f"
    } else {
#line 641 "sgesvxx.f"
	if (rowequ) {
#line 642 "sgesvxx.f"
	    rcmin = bignum;
#line 643 "sgesvxx.f"
	    rcmax = 0.;
#line 644 "sgesvxx.f"
	    i__1 = *n;
#line 644 "sgesvxx.f"
	    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 645 "sgesvxx.f"
		d__1 = rcmin, d__2 = r__[j];
#line 645 "sgesvxx.f"
		rcmin = min(d__1,d__2);
/* Computing MAX */
#line 646 "sgesvxx.f"
		d__1 = rcmax, d__2 = r__[j];
#line 646 "sgesvxx.f"
		rcmax = max(d__1,d__2);
#line 647 "sgesvxx.f"
/* L10: */
#line 647 "sgesvxx.f"
	    }
#line 648 "sgesvxx.f"
	    if (rcmin <= 0.) {
#line 649 "sgesvxx.f"
		*info = -11;
#line 650 "sgesvxx.f"
	    } else if (*n > 0) {
#line 651 "sgesvxx.f"
		rowcnd = max(rcmin,smlnum) / min(rcmax,bignum);
#line 652 "sgesvxx.f"
	    } else {
#line 653 "sgesvxx.f"
		rowcnd = 1.;
#line 654 "sgesvxx.f"
	    }
#line 655 "sgesvxx.f"
	}
#line 656 "sgesvxx.f"
	if (colequ && *info == 0) {
#line 657 "sgesvxx.f"
	    rcmin = bignum;
#line 658 "sgesvxx.f"
	    rcmax = 0.;
#line 659 "sgesvxx.f"
	    i__1 = *n;
#line 659 "sgesvxx.f"
	    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 660 "sgesvxx.f"
		d__1 = rcmin, d__2 = c__[j];
#line 660 "sgesvxx.f"
		rcmin = min(d__1,d__2);
/* Computing MAX */
#line 661 "sgesvxx.f"
		d__1 = rcmax, d__2 = c__[j];
#line 661 "sgesvxx.f"
		rcmax = max(d__1,d__2);
#line 662 "sgesvxx.f"
/* L20: */
#line 662 "sgesvxx.f"
	    }
#line 663 "sgesvxx.f"
	    if (rcmin <= 0.) {
#line 664 "sgesvxx.f"
		*info = -12;
#line 665 "sgesvxx.f"
	    } else if (*n > 0) {
#line 666 "sgesvxx.f"
		colcnd = max(rcmin,smlnum) / min(rcmax,bignum);
#line 667 "sgesvxx.f"
	    } else {
#line 668 "sgesvxx.f"
		colcnd = 1.;
#line 669 "sgesvxx.f"
	    }
#line 670 "sgesvxx.f"
	}
#line 671 "sgesvxx.f"
	if (*info == 0) {
#line 672 "sgesvxx.f"
	    if (*ldb < max(1,*n)) {
#line 673 "sgesvxx.f"
		*info = -14;
#line 674 "sgesvxx.f"
	    } else if (*ldx < max(1,*n)) {
#line 675 "sgesvxx.f"
		*info = -16;
#line 676 "sgesvxx.f"
	    }
#line 677 "sgesvxx.f"
	}
#line 678 "sgesvxx.f"
    }

#line 680 "sgesvxx.f"
    if (*info != 0) {
#line 681 "sgesvxx.f"
	i__1 = -(*info);
#line 681 "sgesvxx.f"
	xerbla_("SGESVXX", &i__1, (ftnlen)7);
#line 682 "sgesvxx.f"
	return 0;
#line 683 "sgesvxx.f"
    }

#line 685 "sgesvxx.f"
    if (equil) {

/*     Compute row and column scalings to equilibrate the matrix A. */

#line 689 "sgesvxx.f"
	sgeequb_(n, n, &a[a_offset], lda, &r__[1], &c__[1], &rowcnd, &colcnd, 
		&amax, &infequ);
#line 691 "sgesvxx.f"
	if (infequ == 0) {

/*     Equilibrate the matrix. */

#line 695 "sgesvxx.f"
	    slaqge_(n, n, &a[a_offset], lda, &r__[1], &c__[1], &rowcnd, &
		    colcnd, &amax, equed, (ftnlen)1);
#line 697 "sgesvxx.f"
	    rowequ = lsame_(equed, "R", (ftnlen)1, (ftnlen)1) || lsame_(equed,
		     "B", (ftnlen)1, (ftnlen)1);
#line 698 "sgesvxx.f"
	    colequ = lsame_(equed, "C", (ftnlen)1, (ftnlen)1) || lsame_(equed,
		     "B", (ftnlen)1, (ftnlen)1);
#line 699 "sgesvxx.f"
	}

/*     If the scaling factors are not applied, set them to 1.0. */

#line 703 "sgesvxx.f"
	if (! rowequ) {
#line 704 "sgesvxx.f"
	    i__1 = *n;
#line 704 "sgesvxx.f"
	    for (j = 1; j <= i__1; ++j) {
#line 705 "sgesvxx.f"
		r__[j] = 1.;
#line 706 "sgesvxx.f"
	    }
#line 707 "sgesvxx.f"
	}
#line 708 "sgesvxx.f"
	if (! colequ) {
#line 709 "sgesvxx.f"
	    i__1 = *n;
#line 709 "sgesvxx.f"
	    for (j = 1; j <= i__1; ++j) {
#line 710 "sgesvxx.f"
		c__[j] = 1.;
#line 711 "sgesvxx.f"
	    }
#line 712 "sgesvxx.f"
	}
#line 713 "sgesvxx.f"
    }

/*     Scale the right-hand side. */

#line 717 "sgesvxx.f"
    if (notran) {
#line 718 "sgesvxx.f"
	if (rowequ) {
#line 718 "sgesvxx.f"
	    slascl2_(n, nrhs, &r__[1], &b[b_offset], ldb);
#line 718 "sgesvxx.f"
	}
#line 719 "sgesvxx.f"
    } else {
#line 720 "sgesvxx.f"
	if (colequ) {
#line 720 "sgesvxx.f"
	    slascl2_(n, nrhs, &c__[1], &b[b_offset], ldb);
#line 720 "sgesvxx.f"
	}
#line 721 "sgesvxx.f"
    }

#line 723 "sgesvxx.f"
    if (nofact || equil) {

/*        Compute the LU factorization of A. */

#line 727 "sgesvxx.f"
	slacpy_("Full", n, n, &a[a_offset], lda, &af[af_offset], ldaf, (
		ftnlen)4);
#line 728 "sgesvxx.f"
	sgetrf_(n, n, &af[af_offset], ldaf, &ipiv[1], info);

/*        Return if INFO is non-zero. */

#line 732 "sgesvxx.f"
	if (*info > 0) {

/*           Pivot in column INFO is exactly 0 */
/*           Compute the reciprocal pivot growth factor of the */
/*           leading rank-deficient INFO columns of A. */

#line 738 "sgesvxx.f"
	    *rpvgrw = sla_gerpvgrw__(n, info, &a[a_offset], lda, &af[
		    af_offset], ldaf);
#line 739 "sgesvxx.f"
	    return 0;
#line 740 "sgesvxx.f"
	}
#line 741 "sgesvxx.f"
    }

/*     Compute the reciprocal pivot growth factor RPVGRW. */

#line 745 "sgesvxx.f"
    *rpvgrw = sla_gerpvgrw__(n, n, &a[a_offset], lda, &af[af_offset], ldaf);

/*     Compute the solution matrix X. */

#line 749 "sgesvxx.f"
    slacpy_("Full", n, nrhs, &b[b_offset], ldb, &x[x_offset], ldx, (ftnlen)4);
#line 750 "sgesvxx.f"
    sgetrs_(trans, n, nrhs, &af[af_offset], ldaf, &ipiv[1], &x[x_offset], ldx,
	     info, (ftnlen)1);

/*     Use iterative refinement to improve the computed solution and */
/*     compute error bounds and backward error estimates for it. */

#line 755 "sgesvxx.f"
    sgerfsx_(trans, equed, n, nrhs, &a[a_offset], lda, &af[af_offset], ldaf, &
	    ipiv[1], &r__[1], &c__[1], &b[b_offset], ldb, &x[x_offset], ldx, 
	    rcond, &berr[1], n_err_bnds__, &err_bnds_norm__[
	    err_bnds_norm_offset], &err_bnds_comp__[err_bnds_comp_offset], 
	    nparams, &params[1], &work[1], &iwork[1], info, (ftnlen)1, (
	    ftnlen)1);

/*     Scale solutions. */

#line 762 "sgesvxx.f"
    if (colequ && notran) {
#line 763 "sgesvxx.f"
	slascl2_(n, nrhs, &c__[1], &x[x_offset], ldx);
#line 764 "sgesvxx.f"
    } else if (rowequ && ! notran) {
#line 765 "sgesvxx.f"
	slascl2_(n, nrhs, &r__[1], &x[x_offset], ldx);
#line 766 "sgesvxx.f"
    }

#line 768 "sgesvxx.f"
    return 0;

/*     End of SGESVXX */
} /* sgesvxx_ */

