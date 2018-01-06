#line 1 "zgbsvxx.f"
/* zgbsvxx.f -- translated by f2c (version 20100827).
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

#line 1 "zgbsvxx.f"
/* > \brief <b> ZGBSVXX computes the solution to system of linear equations A * X = B for GB matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZGBSVXX + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgbsvxx
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgbsvxx
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgbsvxx
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZGBSVXX( FACT, TRANS, N, KL, KU, NRHS, AB, LDAB, AFB, */
/*                           LDAFB, IPIV, EQUED, R, C, B, LDB, X, LDX, */
/*                           RCOND, RPVGRW, BERR, N_ERR_BNDS, */
/*                           ERR_BNDS_NORM, ERR_BNDS_COMP, NPARAMS, PARAMS, */
/*                           WORK, RWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          EQUED, FACT, TRANS */
/*       INTEGER            INFO, LDAB, LDAFB, LDB, LDX, N, NRHS, NPARAMS, */
/*      $                   N_ERR_BNDS */
/*       DOUBLE PRECISION   RCOND, RPVGRW */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IPIV( * ) */
/*       COMPLEX*16         AB( LDAB, * ), AFB( LDAFB, * ), B( LDB, * ), */
/*      $                   X( LDX , * ),WORK( * ) */
/*       DOUBLE PRECISION   R( * ), C( * ), PARAMS( * ), BERR( * ), */
/*      $                   ERR_BNDS_NORM( NRHS, * ), */
/*      $                   ERR_BNDS_COMP( NRHS, * ), RWORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* >    ZGBSVXX uses the LU factorization to compute the solution to a */
/* >    complex*16 system of linear equations  A * X = B,  where A is an */
/* >    N-by-N matrix and X and B are N-by-NRHS matrices. */
/* > */
/* >    If requested, both normwise and maximum componentwise error bounds */
/* >    are returned. ZGBSVXX will return a solution with a tiny */
/* >    guaranteed error (O(eps) where eps is the working machine */
/* >    precision) unless the matrix is very ill-conditioned, in which */
/* >    case a warning is returned. Relevant condition numbers also are */
/* >    calculated and returned. */
/* > */
/* >    ZGBSVXX accepts user-provided factorizations and equilibration */
/* >    factors; see the definitions of the FACT and EQUED options. */
/* >    Solving with refinement and using a factorization from a previous */
/* >    ZGBSVXX call will also produce a solution with either O(eps) */
/* >    errors or warnings, but we cannot make that claim for general */
/* >    user-provided factorizations and equilibration factors if they */
/* >    differ from what ZGBSVXX would itself produce. */
/* > \endverbatim */

/* > \par Description: */
/*  ================= */
/* > */
/* > \verbatim */
/* > */
/* >    The following steps are performed: */
/* > */
/* >    1. If FACT = 'E', double precision scaling factors are computed to equilibrate */
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
/* >          AB is COMPLEX*16 array, dimension (LDAB,N) */
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
/* >          AFB is COMPLEX*16 array, dimension (LDAFB,N) */
/* >     If FACT = 'F', then AFB is an input argument and on entry */
/* >     contains details of the LU factorization of the band matrix */
/* >     A, as computed by ZGBTRF.  U is stored as an upper triangular */
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
/* >     as computed by DGETRF; row i of the matrix was interchanged */
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
/* >          R is DOUBLE PRECISION array, dimension (N) */
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
/* >          C is DOUBLE PRECISION array, dimension (N) */
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
/* >          B is COMPLEX*16 array, dimension (LDB,NRHS) */
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
/* >          X is COMPLEX*16 array, dimension (LDX,NRHS) */
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
/* > \param[out] RPVGRW */
/* > \verbatim */
/* >          RPVGRW is DOUBLE PRECISION */
/* >     Reciprocal pivot growth.  On exit, this contains the reciprocal */
/* >     pivot growth factor norm(A)/norm(U). The "max absolute element" */
/* >     norm is used.  If this is much less than 1, then the stability of */
/* >     the LU factorization of the (equilibrated) matrix A could be poor. */
/* >     This also means that the solution X, estimated condition numbers, */
/* >     and error bounds could be unreliable. If factorization fails with */
/* >     0<INFO<=N, then this contains the reciprocal pivot growth factor */
/* >     for the leading INFO columns of A.  In DGESVX, this quantity is */
/* >     returned in WORK(1). */
/* > \endverbatim */
/* > */
/* > \param[out] BERR */
/* > \verbatim */
/* >          BERR is DOUBLE PRECISION array, dimension (NRHS) */
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
/* >          ERR_BNDS_NORM is DOUBLE PRECISION array, dimension (NRHS, N_ERR_BNDS) */
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
/* >              sqrt(n) * dlamch('Epsilon'). */
/* > */
/* >     err = 2 "Guaranteed" error bound: The estimated forward error, */
/* >              almost certainly within a factor of 10 of the true error */
/* >              so long as the next entry is greater than the threshold */
/* >              sqrt(n) * dlamch('Epsilon'). This error bound should only */
/* >              be trusted if the previous boolean is true. */
/* > */
/* >     err = 3  Reciprocal condition number: Estimated normwise */
/* >              reciprocal condition number.  Compared with the threshold */
/* >              sqrt(n) * dlamch('Epsilon') to determine if the error */
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
/* >          ERR_BNDS_COMP is DOUBLE PRECISION array, dimension (NRHS, N_ERR_BNDS) */
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
/* >              sqrt(n) * dlamch('Epsilon'). */
/* > */
/* >     err = 2 "Guaranteed" error bound: The estimated forward error, */
/* >              almost certainly within a factor of 10 of the true error */
/* >              so long as the next entry is greater than the threshold */
/* >              sqrt(n) * dlamch('Epsilon'). This error bound should only */
/* >              be trusted if the previous boolean is true. */
/* > */
/* >     err = 3  Reciprocal condition number: Estimated componentwise */
/* >              reciprocal condition number.  Compared with the threshold */
/* >              sqrt(n) * dlamch('Epsilon') to determine if the error */
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
/* >          PARAMS is DOUBLE PRECISION array, dimension NPARAMS */
/* >     Specifies algorithm parameters.  If an entry is .LT. 0.0, then */
/* >     that entry will be filled with default value used for that */
/* >     parameter.  Only positions up to NPARAMS are accessed; defaults */
/* >     are used for higher-numbered parameters. */
/* > */
/* >       PARAMS(LA_LINRX_ITREF_I = 1) : Whether to perform iterative */
/* >            refinement or not. */
/* >         Default: 1.0D+0 */
/* >            = 0.0 : No refinement is performed, and no error bounds are */
/* >                    computed. */
/* >            = 1.0 : Use the extra-precise refinement algorithm. */
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
/* >          WORK is COMPLEX*16 array, dimension (2*N) */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is DOUBLE PRECISION array, dimension (2*N) */
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

/* > \ingroup complex16GBsolve */

/*  ===================================================================== */
/* Subroutine */ int zgbsvxx_(char *fact, char *trans, integer *n, integer *
	kl, integer *ku, integer *nrhs, doublecomplex *ab, integer *ldab, 
	doublecomplex *afb, integer *ldafb, integer *ipiv, char *equed, 
	doublereal *r__, doublereal *c__, doublecomplex *b, integer *ldb, 
	doublecomplex *x, integer *ldx, doublereal *rcond, doublereal *rpvgrw,
	 doublereal *berr, integer *n_err_bnds__, doublereal *err_bnds_norm__,
	 doublereal *err_bnds_comp__, integer *nparams, doublereal *params, 
	doublecomplex *work, doublereal *rwork, integer *info, ftnlen 
	fact_len, ftnlen trans_len, ftnlen equed_len)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, afb_dim1, afb_offset, b_dim1, b_offset, 
	    x_dim1, x_offset, err_bnds_norm_dim1, err_bnds_norm_offset, 
	    err_bnds_comp_dim1, err_bnds_comp_offset, i__1, i__2, i__3, i__4;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, j;
    static doublereal amax;
    extern doublereal zla_gbrpvgrw__(integer *, integer *, integer *, integer 
	    *, doublecomplex *, integer *, doublecomplex *, integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal rcmin, rcmax;
    static logical equil;
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal colcnd;
    static logical nofact;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), zlaqgb_(
	    integer *, integer *, integer *, integer *, doublecomplex *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, char *, ftnlen);
    static doublereal bignum;
    static integer infequ;
    static logical colequ;
    static doublereal rowcnd;
    extern /* Subroutine */ int zgbtrf_(integer *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, integer *, integer *);
    static logical notran;
    extern /* Subroutine */ int zlacpy_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen);
    static doublereal smlnum;
    extern /* Subroutine */ int zgbtrs_(char *, integer *, integer *, integer 
	    *, integer *, doublecomplex *, integer *, integer *, 
	    doublecomplex *, integer *, integer *, ftnlen);
    static logical rowequ;
    extern /* Subroutine */ int zlascl2_(integer *, integer *, doublereal *, 
	    doublecomplex *, integer *), zgbequb_(integer *, integer *, 
	    integer *, integer *, doublecomplex *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *)
	    , zgbrfsx_(char *, char *, integer *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *,
	     integer *, doublereal *, doublereal *, doublecomplex *, integer *
	    , doublecomplex *, integer *, doublereal *, doublereal *, integer 
	    *, doublereal *, doublereal *, integer *, doublereal *, 
	    doublecomplex *, doublereal *, integer *, ftnlen, ftnlen);


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

#line 617 "zgbsvxx.f"
    /* Parameter adjustments */
#line 617 "zgbsvxx.f"
    err_bnds_comp_dim1 = *nrhs;
#line 617 "zgbsvxx.f"
    err_bnds_comp_offset = 1 + err_bnds_comp_dim1;
#line 617 "zgbsvxx.f"
    err_bnds_comp__ -= err_bnds_comp_offset;
#line 617 "zgbsvxx.f"
    err_bnds_norm_dim1 = *nrhs;
#line 617 "zgbsvxx.f"
    err_bnds_norm_offset = 1 + err_bnds_norm_dim1;
#line 617 "zgbsvxx.f"
    err_bnds_norm__ -= err_bnds_norm_offset;
#line 617 "zgbsvxx.f"
    ab_dim1 = *ldab;
#line 617 "zgbsvxx.f"
    ab_offset = 1 + ab_dim1;
#line 617 "zgbsvxx.f"
    ab -= ab_offset;
#line 617 "zgbsvxx.f"
    afb_dim1 = *ldafb;
#line 617 "zgbsvxx.f"
    afb_offset = 1 + afb_dim1;
#line 617 "zgbsvxx.f"
    afb -= afb_offset;
#line 617 "zgbsvxx.f"
    --ipiv;
#line 617 "zgbsvxx.f"
    --r__;
#line 617 "zgbsvxx.f"
    --c__;
#line 617 "zgbsvxx.f"
    b_dim1 = *ldb;
#line 617 "zgbsvxx.f"
    b_offset = 1 + b_dim1;
#line 617 "zgbsvxx.f"
    b -= b_offset;
#line 617 "zgbsvxx.f"
    x_dim1 = *ldx;
#line 617 "zgbsvxx.f"
    x_offset = 1 + x_dim1;
#line 617 "zgbsvxx.f"
    x -= x_offset;
#line 617 "zgbsvxx.f"
    --berr;
#line 617 "zgbsvxx.f"
    --params;
#line 617 "zgbsvxx.f"
    --work;
#line 617 "zgbsvxx.f"
    --rwork;
#line 617 "zgbsvxx.f"

#line 617 "zgbsvxx.f"
    /* Function Body */
#line 617 "zgbsvxx.f"
    *info = 0;
#line 618 "zgbsvxx.f"
    nofact = lsame_(fact, "N", (ftnlen)1, (ftnlen)1);
#line 619 "zgbsvxx.f"
    equil = lsame_(fact, "E", (ftnlen)1, (ftnlen)1);
#line 620 "zgbsvxx.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 621 "zgbsvxx.f"
    smlnum = dlamch_("Safe minimum", (ftnlen)12);
#line 622 "zgbsvxx.f"
    bignum = 1. / smlnum;
#line 623 "zgbsvxx.f"
    if (nofact || equil) {
#line 624 "zgbsvxx.f"
	*(unsigned char *)equed = 'N';
#line 625 "zgbsvxx.f"
	rowequ = FALSE_;
#line 626 "zgbsvxx.f"
	colequ = FALSE_;
#line 627 "zgbsvxx.f"
    } else {
#line 628 "zgbsvxx.f"
	rowequ = lsame_(equed, "R", (ftnlen)1, (ftnlen)1) || lsame_(equed, 
		"B", (ftnlen)1, (ftnlen)1);
#line 629 "zgbsvxx.f"
	colequ = lsame_(equed, "C", (ftnlen)1, (ftnlen)1) || lsame_(equed, 
		"B", (ftnlen)1, (ftnlen)1);
#line 630 "zgbsvxx.f"
    }

/*     Default is failure.  If an input parameter is wrong or */
/*     factorization fails, make everything look horrible.  Only the */
/*     pivot growth is set here, the rest is initialized in ZGBRFSX. */

#line 636 "zgbsvxx.f"
    *rpvgrw = 0.;

/*     Test the input parameters.  PARAMS is not tested until DGERFSX. */

#line 640 "zgbsvxx.f"
    if (! nofact && ! equil && ! lsame_(fact, "F", (ftnlen)1, (ftnlen)1)) {
#line 642 "zgbsvxx.f"
	*info = -1;
#line 643 "zgbsvxx.f"
    } else if (! notran && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1) && ! 
	    lsame_(trans, "C", (ftnlen)1, (ftnlen)1)) {
#line 645 "zgbsvxx.f"
	*info = -2;
#line 646 "zgbsvxx.f"
    } else if (*n < 0) {
#line 647 "zgbsvxx.f"
	*info = -3;
#line 648 "zgbsvxx.f"
    } else if (*kl < 0) {
#line 649 "zgbsvxx.f"
	*info = -4;
#line 650 "zgbsvxx.f"
    } else if (*ku < 0) {
#line 651 "zgbsvxx.f"
	*info = -5;
#line 652 "zgbsvxx.f"
    } else if (*nrhs < 0) {
#line 653 "zgbsvxx.f"
	*info = -6;
#line 654 "zgbsvxx.f"
    } else if (*ldab < *kl + *ku + 1) {
#line 655 "zgbsvxx.f"
	*info = -8;
#line 656 "zgbsvxx.f"
    } else if (*ldafb < (*kl << 1) + *ku + 1) {
#line 657 "zgbsvxx.f"
	*info = -10;
#line 658 "zgbsvxx.f"
    } else if (lsame_(fact, "F", (ftnlen)1, (ftnlen)1) && ! (rowequ || colequ 
	    || lsame_(equed, "N", (ftnlen)1, (ftnlen)1))) {
#line 660 "zgbsvxx.f"
	*info = -12;
#line 661 "zgbsvxx.f"
    } else {
#line 662 "zgbsvxx.f"
	if (rowequ) {
#line 663 "zgbsvxx.f"
	    rcmin = bignum;
#line 664 "zgbsvxx.f"
	    rcmax = 0.;
#line 665 "zgbsvxx.f"
	    i__1 = *n;
#line 665 "zgbsvxx.f"
	    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 666 "zgbsvxx.f"
		d__1 = rcmin, d__2 = r__[j];
#line 666 "zgbsvxx.f"
		rcmin = min(d__1,d__2);
/* Computing MAX */
#line 667 "zgbsvxx.f"
		d__1 = rcmax, d__2 = r__[j];
#line 667 "zgbsvxx.f"
		rcmax = max(d__1,d__2);
#line 668 "zgbsvxx.f"
/* L10: */
#line 668 "zgbsvxx.f"
	    }
#line 669 "zgbsvxx.f"
	    if (rcmin <= 0.) {
#line 670 "zgbsvxx.f"
		*info = -13;
#line 671 "zgbsvxx.f"
	    } else if (*n > 0) {
#line 672 "zgbsvxx.f"
		rowcnd = max(rcmin,smlnum) / min(rcmax,bignum);
#line 673 "zgbsvxx.f"
	    } else {
#line 674 "zgbsvxx.f"
		rowcnd = 1.;
#line 675 "zgbsvxx.f"
	    }
#line 676 "zgbsvxx.f"
	}
#line 677 "zgbsvxx.f"
	if (colequ && *info == 0) {
#line 678 "zgbsvxx.f"
	    rcmin = bignum;
#line 679 "zgbsvxx.f"
	    rcmax = 0.;
#line 680 "zgbsvxx.f"
	    i__1 = *n;
#line 680 "zgbsvxx.f"
	    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 681 "zgbsvxx.f"
		d__1 = rcmin, d__2 = c__[j];
#line 681 "zgbsvxx.f"
		rcmin = min(d__1,d__2);
/* Computing MAX */
#line 682 "zgbsvxx.f"
		d__1 = rcmax, d__2 = c__[j];
#line 682 "zgbsvxx.f"
		rcmax = max(d__1,d__2);
#line 683 "zgbsvxx.f"
/* L20: */
#line 683 "zgbsvxx.f"
	    }
#line 684 "zgbsvxx.f"
	    if (rcmin <= 0.) {
#line 685 "zgbsvxx.f"
		*info = -14;
#line 686 "zgbsvxx.f"
	    } else if (*n > 0) {
#line 687 "zgbsvxx.f"
		colcnd = max(rcmin,smlnum) / min(rcmax,bignum);
#line 688 "zgbsvxx.f"
	    } else {
#line 689 "zgbsvxx.f"
		colcnd = 1.;
#line 690 "zgbsvxx.f"
	    }
#line 691 "zgbsvxx.f"
	}
#line 692 "zgbsvxx.f"
	if (*info == 0) {
#line 693 "zgbsvxx.f"
	    if (*ldb < max(1,*n)) {
#line 694 "zgbsvxx.f"
		*info = -15;
#line 695 "zgbsvxx.f"
	    } else if (*ldx < max(1,*n)) {
#line 696 "zgbsvxx.f"
		*info = -16;
#line 697 "zgbsvxx.f"
	    }
#line 698 "zgbsvxx.f"
	}
#line 699 "zgbsvxx.f"
    }

#line 701 "zgbsvxx.f"
    if (*info != 0) {
#line 702 "zgbsvxx.f"
	i__1 = -(*info);
#line 702 "zgbsvxx.f"
	xerbla_("ZGBSVXX", &i__1, (ftnlen)7);
#line 703 "zgbsvxx.f"
	return 0;
#line 704 "zgbsvxx.f"
    }

#line 706 "zgbsvxx.f"
    if (equil) {

/*     Compute row and column scalings to equilibrate the matrix A. */

#line 710 "zgbsvxx.f"
	zgbequb_(n, n, kl, ku, &ab[ab_offset], ldab, &r__[1], &c__[1], &
		rowcnd, &colcnd, &amax, &infequ);
#line 712 "zgbsvxx.f"
	if (infequ == 0) {

/*     Equilibrate the matrix. */

#line 716 "zgbsvxx.f"
	    zlaqgb_(n, n, kl, ku, &ab[ab_offset], ldab, &r__[1], &c__[1], &
		    rowcnd, &colcnd, &amax, equed, (ftnlen)1);
#line 718 "zgbsvxx.f"
	    rowequ = lsame_(equed, "R", (ftnlen)1, (ftnlen)1) || lsame_(equed,
		     "B", (ftnlen)1, (ftnlen)1);
#line 719 "zgbsvxx.f"
	    colequ = lsame_(equed, "C", (ftnlen)1, (ftnlen)1) || lsame_(equed,
		     "B", (ftnlen)1, (ftnlen)1);
#line 720 "zgbsvxx.f"
	}

/*     If the scaling factors are not applied, set them to 1.0. */

#line 724 "zgbsvxx.f"
	if (! rowequ) {
#line 725 "zgbsvxx.f"
	    i__1 = *n;
#line 725 "zgbsvxx.f"
	    for (j = 1; j <= i__1; ++j) {
#line 726 "zgbsvxx.f"
		r__[j] = 1.;
#line 727 "zgbsvxx.f"
	    }
#line 728 "zgbsvxx.f"
	}
#line 729 "zgbsvxx.f"
	if (! colequ) {
#line 730 "zgbsvxx.f"
	    i__1 = *n;
#line 730 "zgbsvxx.f"
	    for (j = 1; j <= i__1; ++j) {
#line 731 "zgbsvxx.f"
		c__[j] = 1.;
#line 732 "zgbsvxx.f"
	    }
#line 733 "zgbsvxx.f"
	}
#line 734 "zgbsvxx.f"
    }

/*     Scale the right-hand side. */

#line 738 "zgbsvxx.f"
    if (notran) {
#line 739 "zgbsvxx.f"
	if (rowequ) {
#line 739 "zgbsvxx.f"
	    zlascl2_(n, nrhs, &r__[1], &b[b_offset], ldb);
#line 739 "zgbsvxx.f"
	}
#line 740 "zgbsvxx.f"
    } else {
#line 741 "zgbsvxx.f"
	if (colequ) {
#line 741 "zgbsvxx.f"
	    zlascl2_(n, nrhs, &c__[1], &b[b_offset], ldb);
#line 741 "zgbsvxx.f"
	}
#line 742 "zgbsvxx.f"
    }

#line 744 "zgbsvxx.f"
    if (nofact || equil) {

/*        Compute the LU factorization of A. */

#line 748 "zgbsvxx.f"
	i__1 = *n;
#line 748 "zgbsvxx.f"
	for (j = 1; j <= i__1; ++j) {
#line 749 "zgbsvxx.f"
	    i__2 = (*kl << 1) + *ku + 1;
#line 749 "zgbsvxx.f"
	    for (i__ = *kl + 1; i__ <= i__2; ++i__) {
#line 750 "zgbsvxx.f"
		i__3 = i__ + j * afb_dim1;
#line 750 "zgbsvxx.f"
		i__4 = i__ - *kl + j * ab_dim1;
#line 750 "zgbsvxx.f"
		afb[i__3].r = ab[i__4].r, afb[i__3].i = ab[i__4].i;
#line 751 "zgbsvxx.f"
/* L30: */
#line 751 "zgbsvxx.f"
	    }
#line 752 "zgbsvxx.f"
/* L40: */
#line 752 "zgbsvxx.f"
	}
#line 753 "zgbsvxx.f"
	zgbtrf_(n, n, kl, ku, &afb[afb_offset], ldafb, &ipiv[1], info);

/*        Return if INFO is non-zero. */

#line 757 "zgbsvxx.f"
	if (*info > 0) {

/*           Pivot in column INFO is exactly 0 */
/*           Compute the reciprocal pivot growth factor of the */
/*           leading rank-deficient INFO columns of A. */

#line 763 "zgbsvxx.f"
	    *rpvgrw = zla_gbrpvgrw__(n, kl, ku, info, &ab[ab_offset], ldab, &
		    afb[afb_offset], ldafb);
#line 765 "zgbsvxx.f"
	    return 0;
#line 766 "zgbsvxx.f"
	}
#line 767 "zgbsvxx.f"
    }

/*     Compute the reciprocal pivot growth factor RPVGRW. */

#line 771 "zgbsvxx.f"
    *rpvgrw = zla_gbrpvgrw__(n, kl, ku, n, &ab[ab_offset], ldab, &afb[
	    afb_offset], ldafb);

/*     Compute the solution matrix X. */

#line 775 "zgbsvxx.f"
    zlacpy_("Full", n, nrhs, &b[b_offset], ldb, &x[x_offset], ldx, (ftnlen)4);
#line 776 "zgbsvxx.f"
    zgbtrs_(trans, n, kl, ku, nrhs, &afb[afb_offset], ldafb, &ipiv[1], &x[
	    x_offset], ldx, info, (ftnlen)1);

/*     Use iterative refinement to improve the computed solution and */
/*     compute error bounds and backward error estimates for it. */

#line 782 "zgbsvxx.f"
    zgbrfsx_(trans, equed, n, kl, ku, nrhs, &ab[ab_offset], ldab, &afb[
	    afb_offset], ldafb, &ipiv[1], &r__[1], &c__[1], &b[b_offset], ldb,
	     &x[x_offset], ldx, rcond, &berr[1], n_err_bnds__, &
	    err_bnds_norm__[err_bnds_norm_offset], &err_bnds_comp__[
	    err_bnds_comp_offset], nparams, &params[1], &work[1], &rwork[1], 
	    info, (ftnlen)1, (ftnlen)1);

/*     Scale solutions. */

#line 790 "zgbsvxx.f"
    if (colequ && notran) {
#line 791 "zgbsvxx.f"
	zlascl2_(n, nrhs, &c__[1], &x[x_offset], ldx);
#line 792 "zgbsvxx.f"
    } else if (rowequ && ! notran) {
#line 793 "zgbsvxx.f"
	zlascl2_(n, nrhs, &r__[1], &x[x_offset], ldx);
#line 794 "zgbsvxx.f"
    }

#line 796 "zgbsvxx.f"
    return 0;

/*     End of ZGBSVXX */

} /* zgbsvxx_ */

