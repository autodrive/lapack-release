#line 1 "sgbsvxx.f"
/* sgbsvxx.f -- translated by f2c (version 20100827).
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

#line 1 "sgbsvxx.f"
/* > \brief <b> SGBSVXX computes the solution to system of linear equations A * X = B for GB matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SGBSVXX + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgbsvxx
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgbsvxx
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgbsvxx
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SGBSVXX( FACT, TRANS, N, KL, KU, NRHS, AB, LDAB, AFB, */
/*                           LDAFB, IPIV, EQUED, R, C, B, LDB, X, LDX, */
/*                           RCOND, RPVGRW, BERR, N_ERR_BNDS, */
/*                           ERR_BNDS_NORM, ERR_BNDS_COMP, NPARAMS, PARAMS, */
/*                           WORK, IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          EQUED, FACT, TRANS */
/*       INTEGER            INFO, LDAB, LDAFB, LDB, LDX, N, NRHS, NPARAMS, */
/*      $                   N_ERR_BNDS */
/*       REAL               RCOND, RPVGRW */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ), IWORK( * ) */
/*       REAL               AB( LDAB, * ), AFB( LDAFB, * ), B( LDB, * ), */
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
/* >    SGBSVXX uses the LU factorization to compute the solution to a */
/* >    real system of linear equations  A * X = B,  where A is an */
/* >    N-by-N matrix and X and B are N-by-NRHS matrices. */
/* > */
/* >    If requested, both normwise and maximum componentwise error bounds */
/* >    are returned. SGBSVXX will return a solution with a tiny */
/* >    guaranteed error (O(eps) where eps is the working machine */
/* >    precision) unless the matrix is very ill-conditioned, in which */
/* >    case a warning is returned. Relevant condition numbers also are */
/* >    calculated and returned. */
/* > */
/* >    SGBSVXX accepts user-provided factorizations and equilibration */
/* >    factors; see the definitions of the FACT and EQUED options. */
/* >    Solving with refinement and using a factorization from a previous */
/* >    SGBSVXX call will also produce a solution with either O(eps) */
/* >    errors or warnings, but we cannot make that claim for general */
/* >    user-provided factorizations and equilibration factors if they */
/* >    differ from what SGBSVXX would itself produce. */
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
/* > \param[in] KL */
/* > \verbatim */
/* >          KL is INTEGER */
/* >     The number of subdiagonals within the band of A.  KL >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] KU */
/* > \verbatim */
/* >          KU is INTEGER */
/* >     The number of superdiagonals within the band of A.  KU >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] NRHS */
/* > \verbatim */
/* >          NRHS is INTEGER */
/* >     The number of right hand sides, i.e., the number of columns */
/* >     of the matrices B and X.  NRHS >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] AB */
/* > \verbatim */
/* >          AB is REAL array, dimension (LDAB,N) */
/* >     On entry, the matrix A in band storage, in rows 1 to KL+KU+1. */
/* >     The j-th column of A is stored in the j-th column of the */
/* >     array AB as follows: */
/* >     AB(KU+1+i-j,j) = A(i,j) for max(1,j-KU)<=i<=min(N,j+kl) */
/* > */
/* >     If FACT = 'F' and EQUED is not 'N', then AB must have been */
/* >     equilibrated by the scaling factors in R and/or C.  AB is not */
/* >     modified if FACT = 'F' or 'N', or if FACT = 'E' and */
/* >     EQUED = 'N' on exit. */
/* > */
/* >     On exit, if EQUED .ne. 'N', A is scaled as follows: */
/* >     EQUED = 'R':  A := diag(R) * A */
/* >     EQUED = 'C':  A := A * diag(C) */
/* >     EQUED = 'B':  A := diag(R) * A * diag(C). */
/* > \endverbatim */
/* > */
/* > \param[in] LDAB */
/* > \verbatim */
/* >          LDAB is INTEGER */
/* >     The leading dimension of the array AB.  LDAB >= KL+KU+1. */
/* > \endverbatim */
/* > */
/* > \param[in,out] AFB */
/* > \verbatim */
/* >          AFB is REAL array, dimension (LDAFB,N) */
/* >     If FACT = 'F', then AFB is an input argument and on entry */
/* >     contains details of the LU factorization of the band matrix */
/* >     A, as computed by SGBTRF.  U is stored as an upper triangular */
/* >     band matrix with KL+KU superdiagonals in rows 1 to KL+KU+1, */
/* >     and the multipliers used during the factorization are stored */
/* >     in rows KL+KU+2 to 2*KL+KU+1.  If EQUED .ne. 'N', then AFB is */
/* >     the factored form of the equilibrated matrix A. */
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
/* > \param[in] LDAFB */
/* > \verbatim */
/* >          LDAFB is INTEGER */
/* >     The leading dimension of the array AFB.  LDAFB >= 2*KL+KU+1. */
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

/* > \ingroup realGBsolve */

/*  ===================================================================== */
/* Subroutine */ int sgbsvxx_(char *fact, char *trans, integer *n, integer *
	kl, integer *ku, integer *nrhs, doublereal *ab, integer *ldab, 
	doublereal *afb, integer *ldafb, integer *ipiv, char *equed, 
	doublereal *r__, doublereal *c__, doublereal *b, integer *ldb, 
	doublereal *x, integer *ldx, doublereal *rcond, doublereal *rpvgrw, 
	doublereal *berr, integer *n_err_bnds__, doublereal *err_bnds_norm__, 
	doublereal *err_bnds_comp__, integer *nparams, doublereal *params, 
	doublereal *work, integer *iwork, integer *info, ftnlen fact_len, 
	ftnlen trans_len, ftnlen equed_len)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, afb_dim1, afb_offset, b_dim1, b_offset, 
	    x_dim1, x_offset, err_bnds_norm_dim1, err_bnds_norm_offset, 
	    err_bnds_comp_dim1, err_bnds_comp_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, j;
    static doublereal amax;
    extern doublereal sla_gbrpvgrw__(integer *, integer *, integer *, integer 
	    *, doublereal *, integer *, doublereal *, integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal rcmin, rcmax;
    static logical equil;
    static doublereal colcnd;
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int slaqgb_(integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, char *, ftnlen);
    static logical nofact;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal bignum;
    static integer infequ;
    static logical colequ;
    extern /* Subroutine */ int sgbtrf_(integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, integer *, integer *), 
	    slacpy_(char *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen);
    static doublereal rowcnd;
    static logical notran;
    extern /* Subroutine */ int sgbtrs_(char *, integer *, integer *, integer 
	    *, integer *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen);
    static doublereal smlnum;
    static logical rowequ;
    extern /* Subroutine */ int slascl2_(integer *, integer *, doublereal *, 
	    doublereal *, integer *), sgbequb_(integer *, integer *, integer *
	    , integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *), sgbrfsx_(
	    char *, char *, integer *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen);


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

#line 620 "sgbsvxx.f"
    /* Parameter adjustments */
#line 620 "sgbsvxx.f"
    err_bnds_comp_dim1 = *nrhs;
#line 620 "sgbsvxx.f"
    err_bnds_comp_offset = 1 + err_bnds_comp_dim1;
#line 620 "sgbsvxx.f"
    err_bnds_comp__ -= err_bnds_comp_offset;
#line 620 "sgbsvxx.f"
    err_bnds_norm_dim1 = *nrhs;
#line 620 "sgbsvxx.f"
    err_bnds_norm_offset = 1 + err_bnds_norm_dim1;
#line 620 "sgbsvxx.f"
    err_bnds_norm__ -= err_bnds_norm_offset;
#line 620 "sgbsvxx.f"
    ab_dim1 = *ldab;
#line 620 "sgbsvxx.f"
    ab_offset = 1 + ab_dim1;
#line 620 "sgbsvxx.f"
    ab -= ab_offset;
#line 620 "sgbsvxx.f"
    afb_dim1 = *ldafb;
#line 620 "sgbsvxx.f"
    afb_offset = 1 + afb_dim1;
#line 620 "sgbsvxx.f"
    afb -= afb_offset;
#line 620 "sgbsvxx.f"
    --ipiv;
#line 620 "sgbsvxx.f"
    --r__;
#line 620 "sgbsvxx.f"
    --c__;
#line 620 "sgbsvxx.f"
    b_dim1 = *ldb;
#line 620 "sgbsvxx.f"
    b_offset = 1 + b_dim1;
#line 620 "sgbsvxx.f"
    b -= b_offset;
#line 620 "sgbsvxx.f"
    x_dim1 = *ldx;
#line 620 "sgbsvxx.f"
    x_offset = 1 + x_dim1;
#line 620 "sgbsvxx.f"
    x -= x_offset;
#line 620 "sgbsvxx.f"
    --berr;
#line 620 "sgbsvxx.f"
    --params;
#line 620 "sgbsvxx.f"
    --work;
#line 620 "sgbsvxx.f"
    --iwork;
#line 620 "sgbsvxx.f"

#line 620 "sgbsvxx.f"
    /* Function Body */
#line 620 "sgbsvxx.f"
    *info = 0;
#line 621 "sgbsvxx.f"
    nofact = lsame_(fact, "N", (ftnlen)1, (ftnlen)1);
#line 622 "sgbsvxx.f"
    equil = lsame_(fact, "E", (ftnlen)1, (ftnlen)1);
#line 623 "sgbsvxx.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 624 "sgbsvxx.f"
    smlnum = slamch_("Safe minimum", (ftnlen)12);
#line 625 "sgbsvxx.f"
    bignum = 1. / smlnum;
#line 626 "sgbsvxx.f"
    if (nofact || equil) {
#line 627 "sgbsvxx.f"
	*(unsigned char *)equed = 'N';
#line 628 "sgbsvxx.f"
	rowequ = FALSE_;
#line 629 "sgbsvxx.f"
	colequ = FALSE_;
#line 630 "sgbsvxx.f"
    } else {
#line 631 "sgbsvxx.f"
	rowequ = lsame_(equed, "R", (ftnlen)1, (ftnlen)1) || lsame_(equed, 
		"B", (ftnlen)1, (ftnlen)1);
#line 632 "sgbsvxx.f"
	colequ = lsame_(equed, "C", (ftnlen)1, (ftnlen)1) || lsame_(equed, 
		"B", (ftnlen)1, (ftnlen)1);
#line 633 "sgbsvxx.f"
    }

/*     Default is failure.  If an input parameter is wrong or */
/*     factorization fails, make everything look horrible.  Only the */
/*     pivot growth is set here, the rest is initialized in SGBRFSX. */

#line 639 "sgbsvxx.f"
    *rpvgrw = 0.;

/*     Test the input parameters.  PARAMS is not tested until SGBRFSX. */

#line 643 "sgbsvxx.f"
    if (! nofact && ! equil && ! lsame_(fact, "F", (ftnlen)1, (ftnlen)1)) {
#line 645 "sgbsvxx.f"
	*info = -1;
#line 646 "sgbsvxx.f"
    } else if (! notran && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1) && ! 
	    lsame_(trans, "C", (ftnlen)1, (ftnlen)1)) {
#line 648 "sgbsvxx.f"
	*info = -2;
#line 649 "sgbsvxx.f"
    } else if (*n < 0) {
#line 650 "sgbsvxx.f"
	*info = -3;
#line 651 "sgbsvxx.f"
    } else if (*kl < 0) {
#line 652 "sgbsvxx.f"
	*info = -4;
#line 653 "sgbsvxx.f"
    } else if (*ku < 0) {
#line 654 "sgbsvxx.f"
	*info = -5;
#line 655 "sgbsvxx.f"
    } else if (*nrhs < 0) {
#line 656 "sgbsvxx.f"
	*info = -6;
#line 657 "sgbsvxx.f"
    } else if (*ldab < *kl + *ku + 1) {
#line 658 "sgbsvxx.f"
	*info = -8;
#line 659 "sgbsvxx.f"
    } else if (*ldafb < (*kl << 1) + *ku + 1) {
#line 660 "sgbsvxx.f"
	*info = -10;
#line 661 "sgbsvxx.f"
    } else if (lsame_(fact, "F", (ftnlen)1, (ftnlen)1) && ! (rowequ || colequ 
	    || lsame_(equed, "N", (ftnlen)1, (ftnlen)1))) {
#line 663 "sgbsvxx.f"
	*info = -12;
#line 664 "sgbsvxx.f"
    } else {
#line 665 "sgbsvxx.f"
	if (rowequ) {
#line 666 "sgbsvxx.f"
	    rcmin = bignum;
#line 667 "sgbsvxx.f"
	    rcmax = 0.;
#line 668 "sgbsvxx.f"
	    i__1 = *n;
#line 668 "sgbsvxx.f"
	    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 669 "sgbsvxx.f"
		d__1 = rcmin, d__2 = r__[j];
#line 669 "sgbsvxx.f"
		rcmin = min(d__1,d__2);
/* Computing MAX */
#line 670 "sgbsvxx.f"
		d__1 = rcmax, d__2 = r__[j];
#line 670 "sgbsvxx.f"
		rcmax = max(d__1,d__2);
#line 671 "sgbsvxx.f"
/* L10: */
#line 671 "sgbsvxx.f"
	    }
#line 672 "sgbsvxx.f"
	    if (rcmin <= 0.) {
#line 673 "sgbsvxx.f"
		*info = -13;
#line 674 "sgbsvxx.f"
	    } else if (*n > 0) {
#line 675 "sgbsvxx.f"
		rowcnd = max(rcmin,smlnum) / min(rcmax,bignum);
#line 676 "sgbsvxx.f"
	    } else {
#line 677 "sgbsvxx.f"
		rowcnd = 1.;
#line 678 "sgbsvxx.f"
	    }
#line 679 "sgbsvxx.f"
	}
#line 680 "sgbsvxx.f"
	if (colequ && *info == 0) {
#line 681 "sgbsvxx.f"
	    rcmin = bignum;
#line 682 "sgbsvxx.f"
	    rcmax = 0.;
#line 683 "sgbsvxx.f"
	    i__1 = *n;
#line 683 "sgbsvxx.f"
	    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 684 "sgbsvxx.f"
		d__1 = rcmin, d__2 = c__[j];
#line 684 "sgbsvxx.f"
		rcmin = min(d__1,d__2);
/* Computing MAX */
#line 685 "sgbsvxx.f"
		d__1 = rcmax, d__2 = c__[j];
#line 685 "sgbsvxx.f"
		rcmax = max(d__1,d__2);
#line 686 "sgbsvxx.f"
/* L20: */
#line 686 "sgbsvxx.f"
	    }
#line 687 "sgbsvxx.f"
	    if (rcmin <= 0.) {
#line 688 "sgbsvxx.f"
		*info = -14;
#line 689 "sgbsvxx.f"
	    } else if (*n > 0) {
#line 690 "sgbsvxx.f"
		colcnd = max(rcmin,smlnum) / min(rcmax,bignum);
#line 691 "sgbsvxx.f"
	    } else {
#line 692 "sgbsvxx.f"
		colcnd = 1.;
#line 693 "sgbsvxx.f"
	    }
#line 694 "sgbsvxx.f"
	}
#line 695 "sgbsvxx.f"
	if (*info == 0) {
#line 696 "sgbsvxx.f"
	    if (*ldb < max(1,*n)) {
#line 697 "sgbsvxx.f"
		*info = -15;
#line 698 "sgbsvxx.f"
	    } else if (*ldx < max(1,*n)) {
#line 699 "sgbsvxx.f"
		*info = -16;
#line 700 "sgbsvxx.f"
	    }
#line 701 "sgbsvxx.f"
	}
#line 702 "sgbsvxx.f"
    }

#line 704 "sgbsvxx.f"
    if (*info != 0) {
#line 705 "sgbsvxx.f"
	i__1 = -(*info);
#line 705 "sgbsvxx.f"
	xerbla_("SGBSVXX", &i__1, (ftnlen)7);
#line 706 "sgbsvxx.f"
	return 0;
#line 707 "sgbsvxx.f"
    }

#line 709 "sgbsvxx.f"
    if (equil) {

/*     Compute row and column scalings to equilibrate the matrix A. */

#line 713 "sgbsvxx.f"
	sgbequb_(n, n, kl, ku, &ab[ab_offset], ldab, &r__[1], &c__[1], &
		rowcnd, &colcnd, &amax, &infequ);
#line 715 "sgbsvxx.f"
	if (infequ == 0) {

/*     Equilibrate the matrix. */

#line 719 "sgbsvxx.f"
	    slaqgb_(n, n, kl, ku, &ab[ab_offset], ldab, &r__[1], &c__[1], &
		    rowcnd, &colcnd, &amax, equed, (ftnlen)1);
#line 721 "sgbsvxx.f"
	    rowequ = lsame_(equed, "R", (ftnlen)1, (ftnlen)1) || lsame_(equed,
		     "B", (ftnlen)1, (ftnlen)1);
#line 722 "sgbsvxx.f"
	    colequ = lsame_(equed, "C", (ftnlen)1, (ftnlen)1) || lsame_(equed,
		     "B", (ftnlen)1, (ftnlen)1);
#line 723 "sgbsvxx.f"
	}

/*     If the scaling factors are not applied, set them to 1.0. */

#line 727 "sgbsvxx.f"
	if (! rowequ) {
#line 728 "sgbsvxx.f"
	    i__1 = *n;
#line 728 "sgbsvxx.f"
	    for (j = 1; j <= i__1; ++j) {
#line 729 "sgbsvxx.f"
		r__[j] = 1.;
#line 730 "sgbsvxx.f"
	    }
#line 731 "sgbsvxx.f"
	}
#line 732 "sgbsvxx.f"
	if (! colequ) {
#line 733 "sgbsvxx.f"
	    i__1 = *n;
#line 733 "sgbsvxx.f"
	    for (j = 1; j <= i__1; ++j) {
#line 734 "sgbsvxx.f"
		c__[j] = 1.;
#line 735 "sgbsvxx.f"
	    }
#line 736 "sgbsvxx.f"
	}
#line 737 "sgbsvxx.f"
    }

/*     Scale the right hand side. */

#line 741 "sgbsvxx.f"
    if (notran) {
#line 742 "sgbsvxx.f"
	if (rowequ) {
#line 742 "sgbsvxx.f"
	    slascl2_(n, nrhs, &r__[1], &b[b_offset], ldb);
#line 742 "sgbsvxx.f"
	}
#line 743 "sgbsvxx.f"
    } else {
#line 744 "sgbsvxx.f"
	if (colequ) {
#line 744 "sgbsvxx.f"
	    slascl2_(n, nrhs, &c__[1], &b[b_offset], ldb);
#line 744 "sgbsvxx.f"
	}
#line 745 "sgbsvxx.f"
    }

#line 747 "sgbsvxx.f"
    if (nofact || equil) {

/*        Compute the LU factorization of A. */

#line 751 "sgbsvxx.f"
	i__1 = *n;
#line 751 "sgbsvxx.f"
	for (j = 1; j <= i__1; ++j) {
#line 752 "sgbsvxx.f"
	    i__2 = (*kl << 1) + *ku + 1;
#line 752 "sgbsvxx.f"
	    for (i__ = *kl + 1; i__ <= i__2; ++i__) {
#line 753 "sgbsvxx.f"
		afb[i__ + j * afb_dim1] = ab[i__ - *kl + j * ab_dim1];
#line 754 "sgbsvxx.f"
/* L30: */
#line 754 "sgbsvxx.f"
	    }
#line 755 "sgbsvxx.f"
/* L40: */
#line 755 "sgbsvxx.f"
	}
#line 756 "sgbsvxx.f"
	sgbtrf_(n, n, kl, ku, &afb[afb_offset], ldafb, &ipiv[1], info);

/*        Return if INFO is non-zero. */

#line 760 "sgbsvxx.f"
	if (*info > 0) {

/*           Pivot in column INFO is exactly 0 */
/*           Compute the reciprocal pivot growth factor of the */
/*           leading rank-deficient INFO columns of A. */

#line 766 "sgbsvxx.f"
	    *rpvgrw = sla_gbrpvgrw__(n, kl, ku, info, &ab[ab_offset], ldab, &
		    afb[afb_offset], ldafb);
#line 768 "sgbsvxx.f"
	    return 0;
#line 769 "sgbsvxx.f"
	}
#line 770 "sgbsvxx.f"
    }

/*     Compute the reciprocal pivot growth factor RPVGRW. */

#line 774 "sgbsvxx.f"
    *rpvgrw = sla_gbrpvgrw__(n, kl, ku, n, &ab[ab_offset], ldab, &afb[
	    afb_offset], ldafb);

/*     Compute the solution matrix X. */

#line 778 "sgbsvxx.f"
    slacpy_("Full", n, nrhs, &b[b_offset], ldb, &x[x_offset], ldx, (ftnlen)4);
#line 779 "sgbsvxx.f"
    sgbtrs_(trans, n, kl, ku, nrhs, &afb[afb_offset], ldafb, &ipiv[1], &x[
	    x_offset], ldx, info, (ftnlen)1);

/*     Use iterative refinement to improve the computed solution and */
/*     compute error bounds and backward error estimates for it. */

#line 785 "sgbsvxx.f"
    sgbrfsx_(trans, equed, n, kl, ku, nrhs, &ab[ab_offset], ldab, &afb[
	    afb_offset], ldafb, &ipiv[1], &r__[1], &c__[1], &b[b_offset], ldb,
	     &x[x_offset], ldx, rcond, &berr[1], n_err_bnds__, &
	    err_bnds_norm__[err_bnds_norm_offset], &err_bnds_comp__[
	    err_bnds_comp_offset], nparams, &params[1], &work[1], &iwork[1], 
	    info, (ftnlen)1, (ftnlen)1);

/*     Scale solutions. */

#line 792 "sgbsvxx.f"
    if (colequ && notran) {
#line 793 "sgbsvxx.f"
	slascl2_(n, nrhs, &c__[1], &x[x_offset], ldx);
#line 794 "sgbsvxx.f"
    } else if (rowequ && ! notran) {
#line 795 "sgbsvxx.f"
	slascl2_(n, nrhs, &r__[1], &x[x_offset], ldx);
#line 796 "sgbsvxx.f"
    }

#line 798 "sgbsvxx.f"
    return 0;

/*     End of SGBSVXX */

} /* sgbsvxx_ */

