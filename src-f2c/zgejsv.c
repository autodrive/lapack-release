#line 1 "zgejsv.f"
/* zgejsv.f -- translated by f2c (version 20100827).
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

#line 1 "zgejsv.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};
static integer c__1 = 1;
static integer c__0 = 0;
static doublereal c_b80 = 1.;
static doublereal c_b120 = 0.;
static integer c_n1 = -1;

/* > \brief \b ZGEJSV */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZGEJSV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgejsv.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgejsv.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgejsv.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*     SUBROUTINE ZGEJSV( JOBA, JOBU, JOBV, JOBR, JOBT, JOBP, */
/*                         M, N, A, LDA, SVA, U, LDU, V, LDV, */
/*                         CWORK, LWORK, RWORK, LRWORK, IWORK, INFO ) */

/*     .. Scalar Arguments .. */
/*     IMPLICIT    NONE */
/*     INTEGER     INFO, LDA, LDU, LDV, LWORK, M, N */
/*     .. */
/*     .. Array Arguments .. */
/*     DOUBLE COMPLEX     A( LDA, * ),  U( LDU, * ), V( LDV, * ), CWORK( LWORK ) */
/*     DOUBLE PRECISION   SVA( N ), RWORK( LRWORK ) */
/*     INTEGER     IWORK( * ) */
/*     CHARACTER*1 JOBA, JOBP, JOBR, JOBT, JOBU, JOBV */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* ZGEJSV computes the singular value decomposition (SVD) of a real M-by-N */
/* matrix [A], where M >= N. The SVD of [A] is written as */

/*              [A] = [U] * [SIGMA] * [V]^*, */

/* where [SIGMA] is an N-by-N (M-by-N) matrix which is zero except for its N */
/* diagonal elements, [U] is an M-by-N (or M-by-M) orthonormal matrix, and */
/* [V] is an N-by-N orthogonal matrix. The diagonal elements of [SIGMA] are */
/* the singular values of [A]. The columns of [U] and [V] are the left and */
/* the right singular vectors of [A], respectively. The matrices [U] and [V] */
/* are computed and stored in the arrays U and V, respectively. The diagonal */
/* of [SIGMA] is computed and stored in the array SVA. */

/*  Arguments: */
/*  ========== */
/* > */
/* > \param[in] JOBA */
/* > \verbatim */
/* >          JOBA is CHARACTER*1 */
/* >         Specifies the level of accuracy: */
/* >       = 'C': This option works well (high relative accuracy) if A = B * D, */
/* >              with well-conditioned B and arbitrary diagonal matrix D. */
/* >              The accuracy cannot be spoiled by COLUMN scaling. The */
/* >              accuracy of the computed output depends on the condition of */
/* >              B, and the procedure aims at the best theoretical accuracy. */
/* >              The relative error max_{i=1:N}|d sigma_i| / sigma_i is */
/* >              bounded by f(M,N)*epsilon* cond(B), independent of D. */
/* >              The input matrix is preprocessed with the QRF with column */
/* >              pivoting. This initial preprocessing and preconditioning by */
/* >              a rank revealing QR factorization is common for all values of */
/* >              JOBA. Additional actions are specified as follows: */
/* >       = 'E': Computation as with 'C' with an additional estimate of the */
/* >              condition number of B. It provides a realistic error bound. */
/* >       = 'F': If A = D1 * C * D2 with ill-conditioned diagonal scalings */
/* >              D1, D2, and well-conditioned matrix C, this option gives */
/* >              higher accuracy than the 'C' option. If the structure of the */
/* >              input matrix is not known, and relative accuracy is */
/* >              desirable, then this option is advisable. The input matrix A */
/* >              is preprocessed with QR factorization with FULL (row and */
/* >              column) pivoting. */
/* >       = 'G'  Computation as with 'F' with an additional estimate of the */
/* >              condition number of B, where A=D*B. If A has heavily weighted */
/* >              rows, then using this condition number gives too pessimistic */
/* >              error bound. */
/* >       = 'A': Small singular values are the noise and the matrix is treated */
/* >              as numerically rank defficient. The error in the computed */
/* >              singular values is bounded by f(m,n)*epsilon*||A||. */
/* >              The computed SVD A = U * S * V^* restores A up to */
/* >              f(m,n)*epsilon*||A||. */
/* >              This gives the procedure the licence to discard (set to zero) */
/* >              all singular values below N*epsilon*||A||. */
/* >       = 'R': Similar as in 'A'. Rank revealing property of the initial */
/* >              QR factorization is used do reveal (using triangular factor) */
/* >              a gap sigma_{r+1} < epsilon * sigma_r in which case the */
/* >              numerical RANK is declared to be r. The SVD is computed with */
/* >              absolute error bounds, but more accurately than with 'A'. */
/* > \endverbatim */
/* > */
/* > \param[in] JOBU */
/* > \verbatim */
/* >          JOBU is CHARACTER*1 */
/* >         Specifies whether to compute the columns of U: */
/* >       = 'U': N columns of U are returned in the array U. */
/* >       = 'F': full set of M left sing. vectors is returned in the array U. */
/* >       = 'W': U may be used as workspace of length M*N. See the description */
/* >              of U. */
/* >       = 'N': U is not computed. */
/* > \endverbatim */
/* > */
/* > \param[in] JOBV */
/* > \verbatim */
/* >          JOBV is CHARACTER*1 */
/* >         Specifies whether to compute the matrix V: */
/* >       = 'V': N columns of V are returned in the array V; Jacobi rotations */
/* >              are not explicitly accumulated. */
/* >       = 'J': N columns of V are returned in the array V, but they are */
/* >              computed as the product of Jacobi rotations. This option is */
/* >              allowed only if JOBU .NE. 'N', i.e. in computing the full SVD. */
/* >       = 'W': V may be used as workspace of length N*N. See the description */
/* >              of V. */
/* >       = 'N': V is not computed. */
/* > \endverbatim */
/* > */
/* > \param[in] JOBR */
/* > \verbatim */
/* >          JOBR is CHARACTER*1 */
/* >         Specifies the RANGE for the singular values. Issues the licence to */
/* >         set to zero small positive singular values if they are outside */
/* >         specified range. If A .NE. 0 is scaled so that the largest singular */
/* >         value of c*A is around SQRT(BIG), BIG=DLAMCH('O'), then JOBR issues */
/* >         the licence to kill columns of A whose norm in c*A is less than */
/* >         SQRT(SFMIN) (for JOBR.EQ.'R'), or less than SMALL=SFMIN/EPSLN, */
/* >         where SFMIN=DLAMCH('S'), EPSLN=DLAMCH('E'). */
/* >       = 'N': Do not kill small columns of c*A. This option assumes that */
/* >              BLAS and QR factorizations and triangular solvers are */
/* >              implemented to work in that range. If the condition of A */
/* >              is greater than BIG, use ZGESVJ. */
/* >       = 'R': RESTRICTED range for sigma(c*A) is [SQRT(SFMIN), SQRT(BIG)] */
/* >              (roughly, as described above). This option is recommended. */
/* >                                             =========================== */
/* >         For computing the singular values in the FULL range [SFMIN,BIG] */
/* >         use ZGESVJ. */
/* > \endverbatim */
/* > */
/* > \param[in] JOBT */
/* > \verbatim */
/* >          JOBT is CHARACTER*1 */
/* >         If the matrix is square then the procedure may determine to use */
/* >         transposed A if A^* seems to be better with respect to convergence. */
/* >         If the matrix is not square, JOBT is ignored. This is subject to */
/* >         changes in the future. */
/* >         The decision is based on two values of entropy over the adjoint */
/* >         orbit of A^* * A. See the descriptions of WORK(6) and WORK(7). */
/* >       = 'T': transpose if entropy test indicates possibly faster */
/* >         convergence of Jacobi process if A^* is taken as input. If A is */
/* >         replaced with A^*, then the row pivoting is included automatically. */
/* >       = 'N': do not speculate. */
/* >         This option can be used to compute only the singular values, or the */
/* >         full SVD (U, SIGMA and V). For only one set of singular vectors */
/* >         (U or V), the caller should provide both U and V, as one of the */
/* >         matrices is used as workspace if the matrix A is transposed. */
/* >         The implementer can easily remove this constraint and make the */
/* >         code more complicated. See the descriptions of U and V. */
/* > \endverbatim */
/* > */
/* > \param[in] JOBP */
/* > \verbatim */
/* >          JOBP is CHARACTER*1 */
/* >         Issues the licence to introduce structured perturbations to drown */
/* >         denormalized numbers. This licence should be active if the */
/* >         denormals are poorly implemented, causing slow computation, */
/* >         especially in cases of fast convergence (!). For details see [1,2]. */
/* >         For the sake of simplicity, this perturbations are included only */
/* >         when the full SVD or only the singular values are requested. The */
/* >         implementer/user can easily add the perturbation for the cases of */
/* >         computing one set of singular vectors. */
/* >       = 'P': introduce perturbation */
/* >       = 'N': do not perturb */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >         The number of rows of the input matrix A.  M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >         The number of columns of the input matrix A. M >= N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is DOUBLE COMPLEX array, dimension (LDA,N) */
/* >          On entry, the M-by-N matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A.  LDA >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] SVA */
/* > \verbatim */
/* >          SVA is DOUBLE PRECISION array, dimension (N) */
/* >          On exit, */
/* >          - For WORK(1)/WORK(2) = ONE: The singular values of A. During the */
/* >            computation SVA contains Euclidean column norms of the */
/* >            iterated matrices in the array A. */
/* >          - For WORK(1) .NE. WORK(2): The singular values of A are */
/* >            (WORK(1)/WORK(2)) * SVA(1:N). This factored form is used if */
/* >            sigma_max(A) overflows or if small singular values have been */
/* >            saved from underflow by scaling the input matrix A. */
/* >          - If JOBR='R' then some of the singular values may be returned */
/* >            as exact zeros obtained by "set to zero" because they are */
/* >            below the numerical rank threshold or are denormalized numbers. */
/* > \endverbatim */
/* > */
/* > \param[out] U */
/* > \verbatim */
/* >          U is DOUBLE COMPLEX array, dimension ( LDU, N ) */
/* >          If JOBU = 'U', then U contains on exit the M-by-N matrix of */
/* >                         the left singular vectors. */
/* >          If JOBU = 'F', then U contains on exit the M-by-M matrix of */
/* >                         the left singular vectors, including an ONB */
/* >                         of the orthogonal complement of the Range(A). */
/* >          If JOBU = 'W'  .AND. (JOBV.EQ.'V' .AND. JOBT.EQ.'T' .AND. M.EQ.N), */
/* >                         then U is used as workspace if the procedure */
/* >                         replaces A with A^*. In that case, [V] is computed */
/* >                         in U as left singular vectors of A^* and then */
/* >                         copied back to the V array. This 'W' option is just */
/* >                         a reminder to the caller that in this case U is */
/* >                         reserved as workspace of length N*N. */
/* >          If JOBU = 'N'  U is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDU */
/* > \verbatim */
/* >          LDU is INTEGER */
/* >          The leading dimension of the array U,  LDU >= 1. */
/* >          IF  JOBU = 'U' or 'F' or 'W',  then LDU >= M. */
/* > \endverbatim */
/* > */
/* > \param[out] V */
/* > \verbatim */
/* >          V is DOUBLE COMPLEX array, dimension ( LDV, N ) */
/* >          If JOBV = 'V', 'J' then V contains on exit the N-by-N matrix of */
/* >                         the right singular vectors; */
/* >          If JOBV = 'W', AND (JOBU.EQ.'U' AND JOBT.EQ.'T' AND M.EQ.N), */
/* >                         then V is used as workspace if the pprocedure */
/* >                         replaces A with A^*. In that case, [U] is computed */
/* >                         in V as right singular vectors of A^* and then */
/* >                         copied back to the U array. This 'W' option is just */
/* >                         a reminder to the caller that in this case V is */
/* >                         reserved as workspace of length N*N. */
/* >          If JOBV = 'N'  V is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDV */
/* > \verbatim */
/* >          LDV is INTEGER */
/* >          The leading dimension of the array V,  LDV >= 1. */
/* >          If JOBV = 'V' or 'J' or 'W', then LDV >= N. */
/* > \endverbatim */
/* > */
/* > \param[out] CWORK */
/* > \verbatim */
/* > CWORK (workspace) */
/* >          CWORK is DOUBLE COMPLEX array, dimension at least LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          Length of CWORK to confirm proper allocation of workspace. */
/* >          LWORK depends on the job: */
/* > */
/* >          1. If only SIGMA is needed ( JOBU.EQ.'N', JOBV.EQ.'N' ) and */
/* >            1.1 .. no scaled condition estimate required (JOBE.EQ.'N'): */
/* >               LWORK >= 2*N+1. This is the minimal requirement. */
/* >               ->> For optimal performance (blocked code) the optimal value */
/* >               is LWORK >= N + (N+1)*NB. Here NB is the optimal */
/* >               block size for ZGEQP3 and ZGEQRF. */
/* >               In general, optimal LWORK is computed as */
/* >               LWORK >= max(N+LWORK(ZGEQP3),N+LWORK(ZGEQRF)). */
/* >            1.2. .. an estimate of the scaled condition number of A is */
/* >               required (JOBA='E', or 'G'). In this case, LWORK the minimal */
/* >               requirement is LWORK >= N*N + 3*N. */
/* >               ->> For optimal performance (blocked code) the optimal value */
/* >               is LWORK >= max(N+(N+1)*NB, N*N+3*N). */
/* >               In general, the optimal length LWORK is computed as */
/* >               LWORK >= max(N+LWORK(ZGEQP3),N+LWORK(ZGEQRF), */
/* >                                                     N+N*N+LWORK(CPOCON)). */
/* > */
/* >          2. If SIGMA and the right singular vectors are needed (JOBV.EQ.'V'), */
/* >             (JOBU.EQ.'N') */
/* >            -> the minimal requirement is LWORK >= 3*N. */
/* >            -> For optimal performance, LWORK >= max(N+(N+1)*NB, 3*N,2*N+N*NB), */
/* >               where NB is the optimal block size for ZGEQP3, ZGEQRF, ZGELQ, */
/* >               CUNMLQ. In general, the optimal length LWORK is computed as */
/* >               LWORK >= max(N+LWORK(ZGEQP3), N+LWORK(CPOCON), N+LWORK(ZGESVJ), */
/* >                       N+LWORK(ZGELQF), 2*N+LWORK(ZGEQRF), N+LWORK(CUNMLQ)). */
/* > */
/* >          3. If SIGMA and the left singular vectors are needed */
/* >            -> the minimal requirement is LWORK >= 3*N. */
/* >            -> For optimal performance: */
/* >               if JOBU.EQ.'U' :: LWORK >= max(3*N, N+(N+1)*NB, 2*N+N*NB), */
/* >               where NB is the optimal block size for ZGEQP3, ZGEQRF, CUNMQR. */
/* >               In general, the optimal length LWORK is computed as */
/* >               LWORK >= max(N+LWORK(ZGEQP3),N+LWORK(CPOCON), */
/* >                        2*N+LWORK(ZGEQRF), N+LWORK(CUNMQR)). */
/* > */
/* >          4. If the full SVD is needed: (JOBU.EQ.'U' or JOBU.EQ.'F') and */
/* >            4.1. if JOBV.EQ.'V' */
/* >               the minimal requirement is LWORK >= 5*N+2*N*N. */
/* >            4.2. if JOBV.EQ.'J' the minimal requirement is */
/* >               LWORK >= 4*N+N*N. */
/* >            In both cases, the allocated CWORK can accomodate blocked runs */
/* >            of ZGEQP3, ZGEQRF, ZGELQF, SUNMQR, CUNMLQ. */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is DOUBLE PRECISION array, dimension at least LRWORK. */
/* >          On exit, */
/* >          RWORK(1) = Determines the scaling factor SCALE = RWORK(2) / RWORK(1) */
/* >                    such that SCALE*SVA(1:N) are the computed singular values */
/* >                    of A. (See the description of SVA().) */
/* >          RWORK(2) = See the description of RWORK(1). */
/* >          RWORK(3) = SCONDA is an estimate for the condition number of */
/* >                    column equilibrated A. (If JOBA .EQ. 'E' or 'G') */
/* >                    SCONDA is an estimate of SQRT(||(R^* * R)^(-1)||_1). */
/* >                    It is computed using SPOCON. It holds */
/* >                    N^(-1/4) * SCONDA <= ||R^(-1)||_2 <= N^(1/4) * SCONDA */
/* >                    where R is the triangular factor from the QRF of A. */
/* >                    However, if R is truncated and the numerical rank is */
/* >                    determined to be strictly smaller than N, SCONDA is */
/* >                    returned as -1, thus indicating that the smallest */
/* >                    singular values might be lost. */
/* > */
/* >          If full SVD is needed, the following two condition numbers are */
/* >          useful for the analysis of the algorithm. They are provied for */
/* >          a developer/implementer who is familiar with the details of */
/* >          the method. */
/* > */
/* >          RWORK(4) = an estimate of the scaled condition number of the */
/* >                    triangular factor in the first QR factorization. */
/* >          RWORK(5) = an estimate of the scaled condition number of the */
/* >                    triangular factor in the second QR factorization. */
/* >          The following two parameters are computed if JOBT .EQ. 'T'. */
/* >          They are provided for a developer/implementer who is familiar */
/* >          with the details of the method. */
/* >          RWORK(6) = the entropy of A^* * A :: this is the Shannon entropy */
/* >                    of diag(A^* * A) / Trace(A^* * A) taken as point in the */
/* >                    probability simplex. */
/* >          RWORK(7) = the entropy of A * A^*. (See the description of RWORK(6).) */
/* > \endverbatim */
/* > */
/* > \param[in] LRWORK */
/* > \verbatim */
/* >          LRWORK is INTEGER */
/* >          Length of RWORK to confirm proper allocation of workspace. */
/* >          LRWORK depends on the job: */
/* > */
/* >       1. If only singular values are requested i.e. if */
/* >          LSAME(JOBU,'N') .AND. LSAME(JOBV,'N') */
/* >          then: */
/* >          1.1. If LSAME(JOBT,'T') .OR. LSAME(JOBA,'F') .OR. LSAME(JOBA,'G'), */
/* >          then LRWORK = max( 7, N + 2 * M ). */
/* >          1.2. Otherwise, LRWORK  = max( 7, 2 * N ). */
/* >       2. If singular values with the right singular vectors are requested */
/* >          i.e. if */
/* >          (LSAME(JOBV,'V').OR.LSAME(JOBV,'J')) .AND. */
/* >          .NOT.(LSAME(JOBU,'U').OR.LSAME(JOBU,'F')) */
/* >          then: */
/* >          2.1. If LSAME(JOBT,'T') .OR. LSAME(JOBA,'F') .OR. LSAME(JOBA,'G'), */
/* >          then LRWORK = max( 7, N + 2 * M ). */
/* >          2.2. Otherwise, LRWORK  = max( 7, 2 * N ). */
/* >       3. If singular values with the left singular vectors are requested, i.e. if */
/* >          (LSAME(JOBU,'U').OR.LSAME(JOBU,'F')) .AND. */
/* >          .NOT.(LSAME(JOBV,'V').OR.LSAME(JOBV,'J')) */
/* >          then: */
/* >          3.1. If LSAME(JOBT,'T') .OR. LSAME(JOBA,'F') .OR. LSAME(JOBA,'G'), */
/* >          then LRWORK = max( 7, N + 2 * M ). */
/* >          3.2. Otherwise, LRWORK  = max( 7, 2 * N ). */
/* >       4. If singular values with both the left and the right singular vectors */
/* >          are requested, i.e. if */
/* >          (LSAME(JOBU,'U').OR.LSAME(JOBU,'F')) .AND. */
/* >          (LSAME(JOBV,'V').OR.LSAME(JOBV,'J')) */
/* >          then: */
/* >          4.1. If LSAME(JOBT,'T') .OR. LSAME(JOBA,'F') .OR. LSAME(JOBA,'G'), */
/* >          then LRWORK = max( 7, N + 2 * M ). */
/* >          4.2. Otherwise, LRWORK  = max( 7, 2 * N ). */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* >          IWORK is INTEGER array, of dimension: */
/* >                If LSAME(JOBA,'F') .OR. LSAME(JOBA,'G'), then */
/* >                the dimension of IWORK is max( 3, 2 * N + M ). */
/* >                Otherwise, the dimension of IWORK is */
/* >                -> max( 3, 2*N ) for full SVD */
/* >                -> max( 3, N ) for singular values only or singular */
/* >                   values with one set of singular vectors (left or right) */
/* >          On exit, */
/* >          IWORK(1) = the numerical rank determined after the initial */
/* >                     QR factorization with pivoting. See the descriptions */
/* >                     of JOBA and JOBR. */
/* >          IWORK(2) = the number of the computed nonzero singular values */
/* >          IWORK(3) = if nonzero, a warning message: */
/* >                     If IWORK(3).EQ.1 then some of the column norms of A */
/* >                     were denormalized floats. The requested high accuracy */
/* >                     is not warranted by the data. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >           < 0  : if INFO = -i, then the i-th argument had an illegal value. */
/* >           = 0 :  successfull exit; */
/* >           > 0 :  ZGEJSV  did not converge in the maximal allowed number */
/* >                  of sweeps. The computed values may be inaccurate. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2015 */

/* > \ingroup complex16GEsing */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  ZGEJSV implements a preconditioned Jacobi SVD algorithm. It uses ZGEQP3, */
/* >  ZGEQRF, and ZGELQF as preprocessors and preconditioners. Optionally, an */
/* >  additional row pivoting can be used as a preprocessor, which in some */
/* >  cases results in much higher accuracy. An example is matrix A with the */
/* >  structure A = D1 * C * D2, where D1, D2 are arbitrarily ill-conditioned */
/* >  diagonal matrices and C is well-conditioned matrix. In that case, complete */
/* >  pivoting in the first QR factorizations provides accuracy dependent on the */
/* >  condition number of C, and independent of D1, D2. Such higher accuracy is */
/* >  not completely understood theoretically, but it works well in practice. */
/* >  Further, if A can be written as A = B*D, with well-conditioned B and some */
/* >  diagonal D, then the high accuracy is guaranteed, both theoretically and */
/* >  in software, independent of D. For more details see [1], [2]. */
/* >     The computational range for the singular values can be the full range */
/* >  ( UNDERFLOW,OVERFLOW ), provided that the machine arithmetic and the BLAS */
/* >  & LAPACK routines called by ZGEJSV are implemented to work in that range. */
/* >  If that is not the case, then the restriction for safe computation with */
/* >  the singular values in the range of normalized IEEE numbers is that the */
/* >  spectral condition number kappa(A)=sigma_max(A)/sigma_min(A) does not */
/* >  overflow. This code (ZGEJSV) is best used in this restricted range, */
/* >  meaning that singular values of magnitude below ||A||_2 / DLAMCH('O') are */
/* >  returned as zeros. See JOBR for details on this. */
/* >     Further, this implementation is somewhat slower than the one described */
/* >  in [1,2] due to replacement of some non-LAPACK components, and because */
/* >  the choice of some tuning parameters in the iterative part (ZGESVJ) is */
/* >  left to the implementer on a particular machine. */
/* >     The rank revealing QR factorization (in this code: ZGEQP3) should be */
/* >  implemented as in [3]. We have a new version of ZGEQP3 under development */
/* >  that is more robust than the current one in LAPACK, with a cleaner cut in */
/* >  rank defficient cases. It will be available in the SIGMA library [4]. */
/* >  If M is much larger than N, it is obvious that the inital QRF with */
/* >  column pivoting can be preprocessed by the QRF without pivoting. That */
/* >  well known trick is not used in ZGEJSV because in some cases heavy row */
/* >  weighting can be treated with complete pivoting. The overhead in cases */
/* >  M much larger than N is then only due to pivoting, but the benefits in */
/* >  terms of accuracy have prevailed. The implementer/user can incorporate */
/* >  this extra QRF step easily. The implementer can also improve data movement */
/* >  (matrix transpose, matrix copy, matrix transposed copy) - this */
/* >  implementation of ZGEJSV uses only the simplest, naive data movement. */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >  Zlatko Drmac (Zagreb, Croatia) and Kresimir Veselic (Hagen, Germany) */

/* > \par References: */
/*  ================ */
/* > */
/* > \verbatim */
/* > */
/* [1] Z. Drmac and K. Veselic: New fast and accurate Jacobi SVD algorithm I. */
/*     SIAM J. Matrix Anal. Appl. Vol. 35, No. 2 (2008), pp. 1322-1342. */
/*     LAPACK Working note 169. */
/* [2] Z. Drmac and K. Veselic: New fast and accurate Jacobi SVD algorithm II. */
/*     SIAM J. Matrix Anal. Appl. Vol. 35, No. 2 (2008), pp. 1343-1362. */
/*     LAPACK Working note 170. */
/* [3] Z. Drmac and Z. Bujanovic: On the failure of rank-revealing QR */
/*     factorization software - a case study. */
/*     ACM Trans. math. Softw. Vol. 35, No 2 (2008), pp. 1-28. */
/*     LAPACK Working note 176. */
/* [4] Z. Drmac: SIGMA - mathematical software library for accurate SVD, PSV, */
/*     QSVD, (H,K)-SVD computations. */
/*     Department of Mathematics, University of Zagreb, 2008. */
/* > \endverbatim */

/* >  \par Bugs, examples and comments: */
/*   ================================= */
/* > */
/* >  Please report all bugs and send interesting examples and/or comments to */
/* >  drmac@math.hr. Thank you. */
/* > */
/*  ===================================================================== */
/* Subroutine */ int zgejsv_(char *joba, char *jobu, char *jobv, char *jobr, 
	char *jobt, char *jobp, integer *m, integer *n, doublecomplex *a, 
	integer *lda, doublereal *sva, doublecomplex *u, integer *ldu, 
	doublecomplex *v, integer *ldv, doublecomplex *cwork, integer *lwork, 
	doublereal *rwork, integer *lrwork, integer *iwork, integer *info, 
	ftnlen joba_len, ftnlen jobu_len, ftnlen jobv_len, ftnlen jobr_len, 
	ftnlen jobt_len, ftnlen jobp_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, u_dim1, u_offset, v_dim1, v_offset, i__1, i__2, 
	    i__3, i__4;
    doublereal d__1, d__2, d__3;
    doublecomplex z__1;

    /* Builtin functions */
    double sqrt(doublereal), z_abs(doublecomplex *), log(doublereal);
    void d_cnjg(doublecomplex *, doublecomplex *), z_div(doublecomplex *, 
	    doublecomplex *, doublecomplex *);
    integer i_dnnt(doublereal *);

    /* Local variables */
    static integer p, q, n1, nr;
    static doublereal big, xsc, big1;
    static logical defr;
    static doublereal aapp, aaqq;
    static logical kill;
    static integer ierr;
    static doublereal temp1;
    static logical jracc;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublecomplex ctemp;
    static doublereal entra, small, sfmin;
    static logical lsvec;
    static doublereal epsln;
    static logical rsvec;
    extern /* Subroutine */ int zcopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), zswap_(integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *);
    static logical l2aber;
    extern /* Subroutine */ int ztrsm_(char *, char *, char *, char *, 
	    integer *, integer *, doublecomplex *, doublecomplex *, integer *,
	     doublecomplex *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static doublereal condr1, condr2, uscal1, uscal2;
    static logical l2kill, l2rank, l2tran, l2pert;
    extern /* Subroutine */ int zgeqp3_(integer *, integer *, doublecomplex *,
	     integer *, integer *, doublecomplex *, doublecomplex *, integer *
	    , doublereal *, integer *);
    extern doublereal dznrm2_(integer *, doublecomplex *, integer *), dlamch_(
	    char *, ftnlen);
    extern integer idamax_();
    static doublereal scalem, sconda;
    static logical goscal;
    static doublereal aatmin, aatmax;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static logical noscal;
    extern /* Subroutine */ int zdscal_(integer *, doublereal *, 
	    doublecomplex *, integer *), zlacgv_(integer *, doublecomplex *, 
	    integer *), zgelqf_(integer *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, integer *)
	    , zlascl_();
    static doublereal entrat;
    static logical almort;
    extern /* Subroutine */ int zgeqrf_(integer *, integer *, doublecomplex *,
	     integer *, doublecomplex *, doublecomplex *, integer *, integer *
	    );
    static doublereal maxprj;
    extern /* Subroutine */ int zlacpy_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen), 
	    zlaset_();
    static logical errest, transp;
    extern /* Subroutine */ int zpocon_(char *, integer *, doublecomplex *, 
	    integer *, doublereal *, doublereal *, doublecomplex *, 
	    doublereal *, integer *, ftnlen), zgesvj_(char *, char *, char *, 
	    integer *, integer *, doublecomplex *, integer *, doublereal *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *,
	     doublereal *, integer *, integer *, ftnlen, ftnlen, ftnlen), 
	    zlassq_(), zlaswp_(integer *, doublecomplex *, integer *, integer 
	    *, integer *, integer *, integer *);
    static logical rowpiv;
    extern /* Subroutine */ int zungqr_(integer *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, integer *), zunmlq_(char *, char *, integer *, integer 
	    *, integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     ftnlen, ftnlen), zunmqr_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     ftnlen, ftnlen);
    static doublereal cond_ok__;
    static integer warning, numrank;


/*  -- LAPACK computational routine (version 3.6.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2015 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  =========================================================================== */

/*     .. Local Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */

/*     .. */

/*     Test the input arguments */

#line 577 "zgejsv.f"
    /* Parameter adjustments */
#line 577 "zgejsv.f"
    --sva;
#line 577 "zgejsv.f"
    a_dim1 = *lda;
#line 577 "zgejsv.f"
    a_offset = 1 + a_dim1;
#line 577 "zgejsv.f"
    a -= a_offset;
#line 577 "zgejsv.f"
    u_dim1 = *ldu;
#line 577 "zgejsv.f"
    u_offset = 1 + u_dim1;
#line 577 "zgejsv.f"
    u -= u_offset;
#line 577 "zgejsv.f"
    v_dim1 = *ldv;
#line 577 "zgejsv.f"
    v_offset = 1 + v_dim1;
#line 577 "zgejsv.f"
    v -= v_offset;
#line 577 "zgejsv.f"
    --cwork;
#line 577 "zgejsv.f"
    --rwork;
#line 577 "zgejsv.f"
    --iwork;
#line 577 "zgejsv.f"

#line 577 "zgejsv.f"
    /* Function Body */
#line 577 "zgejsv.f"
    lsvec = lsame_(jobu, "U", (ftnlen)1, (ftnlen)1) || lsame_(jobu, "F", (
	    ftnlen)1, (ftnlen)1);
#line 578 "zgejsv.f"
    jracc = lsame_(jobv, "J", (ftnlen)1, (ftnlen)1);
#line 579 "zgejsv.f"
    rsvec = lsame_(jobv, "V", (ftnlen)1, (ftnlen)1) || jracc;
#line 580 "zgejsv.f"
    rowpiv = lsame_(joba, "F", (ftnlen)1, (ftnlen)1) || lsame_(joba, "G", (
	    ftnlen)1, (ftnlen)1);
#line 581 "zgejsv.f"
    l2rank = lsame_(joba, "R", (ftnlen)1, (ftnlen)1);
#line 582 "zgejsv.f"
    l2aber = lsame_(joba, "A", (ftnlen)1, (ftnlen)1);
#line 583 "zgejsv.f"
    errest = lsame_(joba, "E", (ftnlen)1, (ftnlen)1) || lsame_(joba, "G", (
	    ftnlen)1, (ftnlen)1);
#line 584 "zgejsv.f"
    l2tran = lsame_(jobt, "T", (ftnlen)1, (ftnlen)1);
#line 585 "zgejsv.f"
    l2kill = lsame_(jobr, "R", (ftnlen)1, (ftnlen)1);
#line 586 "zgejsv.f"
    defr = lsame_(jobr, "N", (ftnlen)1, (ftnlen)1);
#line 587 "zgejsv.f"
    l2pert = lsame_(jobp, "P", (ftnlen)1, (ftnlen)1);

#line 589 "zgejsv.f"
    if (! (rowpiv || l2rank || l2aber || errest || lsame_(joba, "C", (ftnlen)
	    1, (ftnlen)1))) {
#line 591 "zgejsv.f"
	*info = -1;
#line 592 "zgejsv.f"
    } else if (! (lsvec || lsame_(jobu, "N", (ftnlen)1, (ftnlen)1) || lsame_(
	    jobu, "W", (ftnlen)1, (ftnlen)1))) {
#line 594 "zgejsv.f"
	*info = -2;
#line 595 "zgejsv.f"
    } else if (! (rsvec || lsame_(jobv, "N", (ftnlen)1, (ftnlen)1) || lsame_(
	    jobv, "W", (ftnlen)1, (ftnlen)1)) || jracc && ! lsvec) {
#line 597 "zgejsv.f"
	*info = -3;
#line 598 "zgejsv.f"
    } else if (! (l2kill || defr)) {
#line 599 "zgejsv.f"
	*info = -4;
#line 600 "zgejsv.f"
    } else if (! (l2tran || lsame_(jobt, "N", (ftnlen)1, (ftnlen)1))) {
#line 601 "zgejsv.f"
	*info = -5;
#line 602 "zgejsv.f"
    } else if (! (l2pert || lsame_(jobp, "N", (ftnlen)1, (ftnlen)1))) {
#line 603 "zgejsv.f"
	*info = -6;
#line 604 "zgejsv.f"
    } else if (*m < 0) {
#line 605 "zgejsv.f"
	*info = -7;
#line 606 "zgejsv.f"
    } else if (*n < 0 || *n > *m) {
#line 607 "zgejsv.f"
	*info = -8;
#line 608 "zgejsv.f"
    } else if (*lda < *m) {
#line 609 "zgejsv.f"
	*info = -10;
#line 610 "zgejsv.f"
    } else if (lsvec && *ldu < *m) {
#line 611 "zgejsv.f"
	*info = -13;
#line 612 "zgejsv.f"
    } else if (rsvec && *ldv < *n) {
#line 613 "zgejsv.f"
	*info = -15;
#line 614 "zgejsv.f"
    } else if (! (lsvec || rsvec || errest) && *lwork < (*n << 1) + 1 || ! (
	    lsvec || rsvec) && errest && *lwork < *n * *n + *n * 3 || lsvec &&
	     ! rsvec && *lwork < *n * 3 || rsvec && ! lsvec && *lwork < *n * 
	    3 || lsvec && rsvec && ! jracc && *lwork < *n * 5 + (*n << 1) * *
	    n || lsvec && rsvec && jracc && *lwork < (*n << 2) + *n * *n) {
#line 627 "zgejsv.f"
	*info = -17;
#line 628 "zgejsv.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 628 "zgejsv.f"
	i__1 = *n + (*m << 1);
#line 628 "zgejsv.f"
	if (*lrwork < max(i__1,7)) {
#line 629 "zgejsv.f"
	    *info = -19;
#line 630 "zgejsv.f"
	} else {
/*        #:) */
#line 632 "zgejsv.f"
	    *info = 0;
#line 633 "zgejsv.f"
	}
#line 633 "zgejsv.f"
    }

#line 635 "zgejsv.f"
    if (*info != 0) {
/*       #:( */
#line 637 "zgejsv.f"
	i__1 = -(*info);
#line 637 "zgejsv.f"
	xerbla_("ZGEJSV", &i__1, (ftnlen)6);
#line 638 "zgejsv.f"
	return 0;
#line 639 "zgejsv.f"
    }

/*     Quick return for void matrix (Y3K safe) */
/* #:) */
#line 643 "zgejsv.f"
    if (*m == 0 || *n == 0) {
#line 643 "zgejsv.f"
	return 0;
#line 643 "zgejsv.f"
    }

/*     Determine whether the matrix U should be M x N or M x M */

#line 647 "zgejsv.f"
    if (lsvec) {
#line 648 "zgejsv.f"
	n1 = *n;
#line 649 "zgejsv.f"
	if (lsame_(jobu, "F", (ftnlen)1, (ftnlen)1)) {
#line 649 "zgejsv.f"
	    n1 = *m;
#line 649 "zgejsv.f"
	}
#line 650 "zgejsv.f"
    }

/*     Set numerical parameters */

/* !    NOTE: Make sure DLAMCH() does not fail on the target architecture. */

#line 656 "zgejsv.f"
    epsln = dlamch_("Epsilon", (ftnlen)7);
#line 657 "zgejsv.f"
    sfmin = dlamch_("SafeMinimum", (ftnlen)11);
#line 658 "zgejsv.f"
    small = sfmin / epsln;
#line 659 "zgejsv.f"
    big = dlamch_("O", (ftnlen)1);
/*     BIG   = ONE / SFMIN */

/*     Initialize SVA(1:N) = diag( ||A e_i||_2 )_1^N */

/* (!)  If necessary, scale SVA() to protect the largest norm from */
/*     overflow. It is possible that this scaling pushes the smallest */
/*     column norm left from the underflow threshold (extreme case). */

#line 668 "zgejsv.f"
    scalem = 1. / sqrt((doublereal) (*m) * (doublereal) (*n));
#line 669 "zgejsv.f"
    noscal = TRUE_;
#line 670 "zgejsv.f"
    goscal = TRUE_;
#line 671 "zgejsv.f"
    i__1 = *n;
#line 671 "zgejsv.f"
    for (p = 1; p <= i__1; ++p) {
#line 672 "zgejsv.f"
	aapp = 0.;
#line 673 "zgejsv.f"
	aaqq = 1.;
#line 674 "zgejsv.f"
	zlassq_(m, &a[p * a_dim1 + 1], &c__1, &aapp, &aaqq);
#line 675 "zgejsv.f"
	if (aapp > big) {
#line 676 "zgejsv.f"
	    *info = -9;
#line 677 "zgejsv.f"
	    i__2 = -(*info);
#line 677 "zgejsv.f"
	    xerbla_("ZGEJSV", &i__2, (ftnlen)6);
#line 678 "zgejsv.f"
	    return 0;
#line 679 "zgejsv.f"
	}
#line 680 "zgejsv.f"
	aaqq = sqrt(aaqq);
#line 681 "zgejsv.f"
	if (aapp < big / aaqq && noscal) {
#line 682 "zgejsv.f"
	    sva[p] = aapp * aaqq;
#line 683 "zgejsv.f"
	} else {
#line 684 "zgejsv.f"
	    noscal = FALSE_;
#line 685 "zgejsv.f"
	    sva[p] = aapp * (aaqq * scalem);
#line 686 "zgejsv.f"
	    if (goscal) {
#line 687 "zgejsv.f"
		goscal = FALSE_;
#line 688 "zgejsv.f"
		i__2 = p - 1;
#line 688 "zgejsv.f"
		dscal_(&i__2, &scalem, &sva[1], &c__1);
#line 689 "zgejsv.f"
	    }
#line 690 "zgejsv.f"
	}
#line 691 "zgejsv.f"
/* L1874: */
#line 691 "zgejsv.f"
    }

#line 693 "zgejsv.f"
    if (noscal) {
#line 693 "zgejsv.f"
	scalem = 1.;
#line 693 "zgejsv.f"
    }

#line 695 "zgejsv.f"
    aapp = 0.;
#line 696 "zgejsv.f"
    aaqq = big;
#line 697 "zgejsv.f"
    i__1 = *n;
#line 697 "zgejsv.f"
    for (p = 1; p <= i__1; ++p) {
/* Computing MAX */
#line 698 "zgejsv.f"
	d__1 = aapp, d__2 = sva[p];
#line 698 "zgejsv.f"
	aapp = max(d__1,d__2);
#line 699 "zgejsv.f"
	if (sva[p] != 0.) {
/* Computing MIN */
#line 699 "zgejsv.f"
	    d__1 = aaqq, d__2 = sva[p];
#line 699 "zgejsv.f"
	    aaqq = min(d__1,d__2);
#line 699 "zgejsv.f"
	}
#line 700 "zgejsv.f"
/* L4781: */
#line 700 "zgejsv.f"
    }

/*     Quick return for zero M x N matrix */
/* #:) */
#line 704 "zgejsv.f"
    if (aapp == 0.) {
#line 705 "zgejsv.f"
	if (lsvec) {
#line 705 "zgejsv.f"
	    zlaset_("G", m, &n1, &c_b1, &c_b2, &u[u_offset], ldu, (ftnlen)1);
#line 705 "zgejsv.f"
	}
#line 706 "zgejsv.f"
	if (rsvec) {
#line 706 "zgejsv.f"
	    zlaset_("G", n, n, &c_b1, &c_b2, &v[v_offset], ldv, (ftnlen)1);
#line 706 "zgejsv.f"
	}
#line 707 "zgejsv.f"
	rwork[1] = 1.;
#line 708 "zgejsv.f"
	rwork[2] = 1.;
#line 709 "zgejsv.f"
	if (errest) {
#line 709 "zgejsv.f"
	    rwork[3] = 1.;
#line 709 "zgejsv.f"
	}
#line 710 "zgejsv.f"
	if (lsvec && rsvec) {
#line 711 "zgejsv.f"
	    rwork[4] = 1.;
#line 712 "zgejsv.f"
	    rwork[5] = 1.;
#line 713 "zgejsv.f"
	}
#line 714 "zgejsv.f"
	if (l2tran) {
#line 715 "zgejsv.f"
	    rwork[6] = 0.;
#line 716 "zgejsv.f"
	    rwork[7] = 0.;
#line 717 "zgejsv.f"
	}
#line 718 "zgejsv.f"
	iwork[1] = 0;
#line 719 "zgejsv.f"
	iwork[2] = 0;
#line 720 "zgejsv.f"
	iwork[3] = 0;
#line 721 "zgejsv.f"
	return 0;
#line 722 "zgejsv.f"
    }

/*     Issue warning if denormalized column norms detected. Override the */
/*     high relative accuracy request. Issue licence to kill columns */
/*     (set them to zero) whose norm is less than sigma_max / BIG (roughly). */
/* #:( */
#line 728 "zgejsv.f"
    warning = 0;
#line 729 "zgejsv.f"
    if (aaqq <= sfmin) {
#line 730 "zgejsv.f"
	l2rank = TRUE_;
#line 731 "zgejsv.f"
	l2kill = TRUE_;
#line 732 "zgejsv.f"
	warning = 1;
#line 733 "zgejsv.f"
    }

/*     Quick return for one-column matrix */
/* #:) */
#line 737 "zgejsv.f"
    if (*n == 1) {

#line 739 "zgejsv.f"
	if (lsvec) {
#line 740 "zgejsv.f"
	    zlascl_("G", &c__0, &c__0, &sva[1], &scalem, m, &c__1, &a[a_dim1 
		    + 1], lda, &ierr, (ftnlen)1);
#line 741 "zgejsv.f"
	    zlacpy_("A", m, &c__1, &a[a_offset], lda, &u[u_offset], ldu, (
		    ftnlen)1);
/*           computing all M left singular vectors of the M x 1 matrix */
#line 743 "zgejsv.f"
	    if (n1 != *n) {
#line 744 "zgejsv.f"
		i__1 = *lwork - *n;
#line 744 "zgejsv.f"
		zgeqrf_(m, n, &u[u_offset], ldu, &cwork[1], &cwork[*n + 1], &
			i__1, &ierr);
#line 745 "zgejsv.f"
		i__1 = *lwork - *n;
#line 745 "zgejsv.f"
		zungqr_(m, &n1, &c__1, &u[u_offset], ldu, &cwork[1], &cwork[*
			n + 1], &i__1, &ierr);
#line 746 "zgejsv.f"
		zcopy_(m, &a[a_dim1 + 1], &c__1, &u[u_dim1 + 1], &c__1);
#line 747 "zgejsv.f"
	    }
#line 748 "zgejsv.f"
	}
#line 749 "zgejsv.f"
	if (rsvec) {
#line 750 "zgejsv.f"
	    i__1 = v_dim1 + 1;
#line 750 "zgejsv.f"
	    v[i__1].r = 1., v[i__1].i = 0.;
#line 751 "zgejsv.f"
	}
#line 752 "zgejsv.f"
	if (sva[1] < big * scalem) {
#line 753 "zgejsv.f"
	    sva[1] /= scalem;
#line 754 "zgejsv.f"
	    scalem = 1.;
#line 755 "zgejsv.f"
	}
#line 756 "zgejsv.f"
	rwork[1] = 1. / scalem;
#line 757 "zgejsv.f"
	rwork[2] = 1.;
#line 758 "zgejsv.f"
	if (sva[1] != 0.) {
#line 759 "zgejsv.f"
	    iwork[1] = 1;
#line 760 "zgejsv.f"
	    if (sva[1] / scalem >= sfmin) {
#line 761 "zgejsv.f"
		iwork[2] = 1;
#line 762 "zgejsv.f"
	    } else {
#line 763 "zgejsv.f"
		iwork[2] = 0;
#line 764 "zgejsv.f"
	    }
#line 765 "zgejsv.f"
	} else {
#line 766 "zgejsv.f"
	    iwork[1] = 0;
#line 767 "zgejsv.f"
	    iwork[2] = 0;
#line 768 "zgejsv.f"
	}
#line 769 "zgejsv.f"
	iwork[3] = 0;
#line 770 "zgejsv.f"
	if (errest) {
#line 770 "zgejsv.f"
	    rwork[3] = 1.;
#line 770 "zgejsv.f"
	}
#line 771 "zgejsv.f"
	if (lsvec && rsvec) {
#line 772 "zgejsv.f"
	    rwork[4] = 1.;
#line 773 "zgejsv.f"
	    rwork[5] = 1.;
#line 774 "zgejsv.f"
	}
#line 775 "zgejsv.f"
	if (l2tran) {
#line 776 "zgejsv.f"
	    rwork[6] = 0.;
#line 777 "zgejsv.f"
	    rwork[7] = 0.;
#line 778 "zgejsv.f"
	}
#line 779 "zgejsv.f"
	return 0;

#line 781 "zgejsv.f"
    }

#line 783 "zgejsv.f"
    transp = FALSE_;
#line 784 "zgejsv.f"
    l2tran = l2tran && *m == *n;

#line 786 "zgejsv.f"
    aatmax = -1.;
#line 787 "zgejsv.f"
    aatmin = big;
#line 788 "zgejsv.f"
    if (rowpiv || l2tran) {

/*     Compute the row norms, needed to determine row pivoting sequence */
/*     (in the case of heavily row weighted A, row pivoting is strongly */
/*     advised) and to collect information needed to compare the */
/*     structures of A * A^* and A^* * A (in the case L2TRAN.EQ..TRUE.). */

#line 795 "zgejsv.f"
	if (l2tran) {
#line 796 "zgejsv.f"
	    i__1 = *m;
#line 796 "zgejsv.f"
	    for (p = 1; p <= i__1; ++p) {
#line 797 "zgejsv.f"
		xsc = 0.;
#line 798 "zgejsv.f"
		temp1 = 1.;
#line 799 "zgejsv.f"
		zlassq_(n, &a[p + a_dim1], lda, &xsc, &temp1);
/*              ZLASSQ gets both the ell_2 and the ell_infinity norm */
/*              in one pass through the vector */
#line 802 "zgejsv.f"
		rwork[*m + *n + p] = xsc * scalem;
#line 803 "zgejsv.f"
		rwork[*n + p] = xsc * (scalem * sqrt(temp1));
/* Computing MAX */
#line 804 "zgejsv.f"
		d__1 = aatmax, d__2 = rwork[*n + p];
#line 804 "zgejsv.f"
		aatmax = max(d__1,d__2);
#line 805 "zgejsv.f"
		if (rwork[*n + p] != 0.) {
/* Computing MIN */
#line 805 "zgejsv.f"
		    d__1 = aatmin, d__2 = rwork[*n + p];
#line 805 "zgejsv.f"
		    aatmin = min(d__1,d__2);
#line 805 "zgejsv.f"
		}
#line 807 "zgejsv.f"
/* L1950: */
#line 807 "zgejsv.f"
	    }
#line 808 "zgejsv.f"
	} else {
#line 809 "zgejsv.f"
	    i__1 = *m;
#line 809 "zgejsv.f"
	    for (p = 1; p <= i__1; ++p) {
#line 810 "zgejsv.f"
		rwork[*m + *n + p] = scalem * z_abs(&a[p + idamax_(n, &a[p + 
			a_dim1], lda) * a_dim1]);
/* Computing MAX */
#line 811 "zgejsv.f"
		d__1 = aatmax, d__2 = rwork[*m + *n + p];
#line 811 "zgejsv.f"
		aatmax = max(d__1,d__2);
/* Computing MIN */
#line 812 "zgejsv.f"
		d__1 = aatmin, d__2 = rwork[*m + *n + p];
#line 812 "zgejsv.f"
		aatmin = min(d__1,d__2);
#line 813 "zgejsv.f"
/* L1904: */
#line 813 "zgejsv.f"
	    }
#line 814 "zgejsv.f"
	}

#line 816 "zgejsv.f"
    }

/*     For square matrix A try to determine whether A^*  would be  better */
/*     input for the preconditioned Jacobi SVD, with faster convergence. */
/*     The decision is based on an O(N) function of the vector of column */
/*     and row norms of A, based on the Shannon entropy. This should give */
/*     the right choice in most cases when the difference actually matters. */
/*     It may fail and pick the slower converging side. */

#line 825 "zgejsv.f"
    entra = 0.;
#line 826 "zgejsv.f"
    entrat = 0.;
#line 827 "zgejsv.f"
    if (l2tran) {

#line 829 "zgejsv.f"
	xsc = 0.;
#line 830 "zgejsv.f"
	temp1 = 1.;
#line 831 "zgejsv.f"
	zlassq_(n, &sva[1], &c__1, &xsc, &temp1);
#line 832 "zgejsv.f"
	temp1 = 1. / temp1;

#line 834 "zgejsv.f"
	entra = 0.;
#line 835 "zgejsv.f"
	i__1 = *n;
#line 835 "zgejsv.f"
	for (p = 1; p <= i__1; ++p) {
/* Computing 2nd power */
#line 836 "zgejsv.f"
	    d__1 = sva[p] / xsc;
#line 836 "zgejsv.f"
	    big1 = d__1 * d__1 * temp1;
#line 837 "zgejsv.f"
	    if (big1 != 0.) {
#line 837 "zgejsv.f"
		entra += big1 * log(big1);
#line 837 "zgejsv.f"
	    }
#line 838 "zgejsv.f"
/* L1113: */
#line 838 "zgejsv.f"
	}
#line 839 "zgejsv.f"
	entra = -entra / log((doublereal) (*n));

/*        Now, SVA().^2/Trace(A^* * A) is a point in the probability simplex. */
/*        It is derived from the diagonal of  A^* * A.  Do the same with the */
/*        diagonal of A * A^*, compute the entropy of the corresponding */
/*        probability distribution. Note that A * A^* and A^* * A have the */
/*        same trace. */

#line 847 "zgejsv.f"
	entrat = 0.;
#line 848 "zgejsv.f"
	i__1 = *n + *m;
#line 848 "zgejsv.f"
	for (p = *n + 1; p <= i__1; ++p) {
/* Computing 2nd power */
#line 849 "zgejsv.f"
	    d__1 = rwork[p] / xsc;
#line 849 "zgejsv.f"
	    big1 = d__1 * d__1 * temp1;
#line 850 "zgejsv.f"
	    if (big1 != 0.) {
#line 850 "zgejsv.f"
		entrat += big1 * log(big1);
#line 850 "zgejsv.f"
	    }
#line 851 "zgejsv.f"
/* L1114: */
#line 851 "zgejsv.f"
	}
#line 852 "zgejsv.f"
	entrat = -entrat / log((doublereal) (*m));

/*        Analyze the entropies and decide A or A^*. Smaller entropy */
/*        usually means better input for the algorithm. */

#line 857 "zgejsv.f"
	transp = entrat < entra;
#line 858 "zgejsv.f"
	transp = TRUE_;

/*        If A^* is better than A, take the adjoint of A. */

#line 862 "zgejsv.f"
	if (transp) {
/*           In an optimal implementation, this trivial transpose */
/*           should be replaced with faster transpose. */
#line 865 "zgejsv.f"
	    i__1 = *n - 1;
#line 865 "zgejsv.f"
	    for (p = 1; p <= i__1; ++p) {
#line 866 "zgejsv.f"
		i__2 = p + p * a_dim1;
#line 866 "zgejsv.f"
		d_cnjg(&z__1, &a[p + p * a_dim1]);
#line 866 "zgejsv.f"
		a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 867 "zgejsv.f"
		i__2 = *n;
#line 867 "zgejsv.f"
		for (q = p + 1; q <= i__2; ++q) {
#line 868 "zgejsv.f"
		    d_cnjg(&z__1, &a[q + p * a_dim1]);
#line 868 "zgejsv.f"
		    ctemp.r = z__1.r, ctemp.i = z__1.i;
#line 869 "zgejsv.f"
		    i__3 = q + p * a_dim1;
#line 869 "zgejsv.f"
		    d_cnjg(&z__1, &a[p + q * a_dim1]);
#line 869 "zgejsv.f"
		    a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 870 "zgejsv.f"
		    i__3 = p + q * a_dim1;
#line 870 "zgejsv.f"
		    a[i__3].r = ctemp.r, a[i__3].i = ctemp.i;
#line 871 "zgejsv.f"
/* L1116: */
#line 871 "zgejsv.f"
		}
#line 872 "zgejsv.f"
/* L1115: */
#line 872 "zgejsv.f"
	    }
#line 873 "zgejsv.f"
	    i__1 = *n + *n * a_dim1;
#line 873 "zgejsv.f"
	    d_cnjg(&z__1, &a[*n + *n * a_dim1]);
#line 873 "zgejsv.f"
	    a[i__1].r = z__1.r, a[i__1].i = z__1.i;
#line 874 "zgejsv.f"
	    i__1 = *n;
#line 874 "zgejsv.f"
	    for (p = 1; p <= i__1; ++p) {
#line 875 "zgejsv.f"
		rwork[*m + *n + p] = sva[p];
#line 876 "zgejsv.f"
		sva[p] = rwork[*n + p];
/*              previously computed row 2-norms are now column 2-norms */
/*              of the transposed matrix */
#line 879 "zgejsv.f"
/* L1117: */
#line 879 "zgejsv.f"
	    }
#line 880 "zgejsv.f"
	    temp1 = aapp;
#line 881 "zgejsv.f"
	    aapp = aatmax;
#line 882 "zgejsv.f"
	    aatmax = temp1;
#line 883 "zgejsv.f"
	    temp1 = aaqq;
#line 884 "zgejsv.f"
	    aaqq = aatmin;
#line 885 "zgejsv.f"
	    aatmin = temp1;
#line 886 "zgejsv.f"
	    kill = lsvec;
#line 887 "zgejsv.f"
	    lsvec = rsvec;
#line 888 "zgejsv.f"
	    rsvec = kill;
#line 889 "zgejsv.f"
	    if (lsvec) {
#line 889 "zgejsv.f"
		n1 = *n;
#line 889 "zgejsv.f"
	    }

#line 891 "zgejsv.f"
	    rowpiv = TRUE_;
#line 892 "zgejsv.f"
	}

#line 894 "zgejsv.f"
    }
/*     END IF L2TRAN */

/*     Scale the matrix so that its maximal singular value remains less */
/*     than SQRT(BIG) -- the matrix is scaled so that its maximal column */
/*     has Euclidean norm equal to SQRT(BIG/N). The only reason to keep */
/*     SQRT(BIG) instead of BIG is the fact that ZGEJSV uses LAPACK and */
/*     BLAS routines that, in some implementations, are not capable of */
/*     working in the full interval [SFMIN,BIG] and that they may provoke */
/*     overflows in the intermediate results. If the singular values spread */
/*     from SFMIN to BIG, then ZGESVJ will compute them. So, in that case, */
/*     one should use ZGESVJ instead of ZGEJSV. */

#line 907 "zgejsv.f"
    big1 = sqrt(big);
#line 908 "zgejsv.f"
    temp1 = sqrt(big / (doublereal) (*n));

#line 910 "zgejsv.f"
    zlascl_("G", &c__0, &c__0, &aapp, &temp1, n, &c__1, &sva[1], n, &ierr, (
	    ftnlen)1);
#line 911 "zgejsv.f"
    if (aaqq > aapp * sfmin) {
#line 912 "zgejsv.f"
	aaqq = aaqq / aapp * temp1;
#line 913 "zgejsv.f"
    } else {
#line 914 "zgejsv.f"
	aaqq = aaqq * temp1 / aapp;
#line 915 "zgejsv.f"
    }
#line 916 "zgejsv.f"
    temp1 *= scalem;
#line 917 "zgejsv.f"
    zlascl_("G", &c__0, &c__0, &aapp, &temp1, m, n, &a[a_offset], lda, &ierr, 
	    (ftnlen)1);

/*     To undo scaling at the end of this procedure, multiply the */
/*     computed singular values with USCAL2 / USCAL1. */

#line 922 "zgejsv.f"
    uscal1 = temp1;
#line 923 "zgejsv.f"
    uscal2 = aapp;

#line 925 "zgejsv.f"
    if (l2kill) {
/*        L2KILL enforces computation of nonzero singular values in */
/*        the restricted range of condition number of the initial A, */
/*        sigma_max(A) / sigma_min(A) approx. SQRT(BIG)/SQRT(SFMIN). */
#line 929 "zgejsv.f"
	xsc = sqrt(sfmin);
#line 930 "zgejsv.f"
    } else {
#line 931 "zgejsv.f"
	xsc = small;

/*        Now, if the condition number of A is too big, */
/*        sigma_max(A) / sigma_min(A) .GT. SQRT(BIG/N) * EPSLN / SFMIN, */
/*        as a precaution measure, the full SVD is computed using ZGESVJ */
/*        with accumulated Jacobi rotations. This provides numerically */
/*        more robust computation, at the cost of slightly increased run */
/*        time. Depending on the concrete implementation of BLAS and LAPACK */
/*        (i.e. how they behave in presence of extreme ill-conditioning) the */
/*        implementor may decide to remove this switch. */
#line 941 "zgejsv.f"
	if (aaqq < sqrt(sfmin) && lsvec && rsvec) {
#line 942 "zgejsv.f"
	    jracc = TRUE_;
#line 943 "zgejsv.f"
	}

#line 945 "zgejsv.f"
    }
#line 946 "zgejsv.f"
    if (aaqq < xsc) {
#line 947 "zgejsv.f"
	i__1 = *n;
#line 947 "zgejsv.f"
	for (p = 1; p <= i__1; ++p) {
#line 948 "zgejsv.f"
	    if (sva[p] < xsc) {
#line 949 "zgejsv.f"
		zlaset_("A", m, &c__1, &c_b1, &c_b1, &a[p * a_dim1 + 1], lda, 
			(ftnlen)1);
#line 950 "zgejsv.f"
		sva[p] = 0.;
#line 951 "zgejsv.f"
	    }
#line 952 "zgejsv.f"
/* L700: */
#line 952 "zgejsv.f"
	}
#line 953 "zgejsv.f"
    }

/*     Preconditioning using QR factorization with pivoting */

#line 957 "zgejsv.f"
    if (rowpiv) {
/*        Optional row permutation (Bjoerck row pivoting): */
/*        A result by Cox and Higham shows that the Bjoerck's */
/*        row pivoting combined with standard column pivoting */
/*        has similar effect as Powell-Reid complete pivoting. */
/*        The ell-infinity norms of A are made nonincreasing. */
#line 963 "zgejsv.f"
	i__1 = *m - 1;
#line 963 "zgejsv.f"
	for (p = 1; p <= i__1; ++p) {
#line 964 "zgejsv.f"
	    i__2 = *m - p + 1;
#line 964 "zgejsv.f"
	    q = idamax_(&i__2, &rwork[*m + *n + p], &c__1) + p - 1;
#line 965 "zgejsv.f"
	    iwork[(*n << 1) + p] = q;
#line 966 "zgejsv.f"
	    if (p != q) {
#line 967 "zgejsv.f"
		temp1 = rwork[*m + *n + p];
#line 968 "zgejsv.f"
		rwork[*m + *n + p] = rwork[*m + *n + q];
#line 969 "zgejsv.f"
		rwork[*m + *n + q] = temp1;
#line 970 "zgejsv.f"
	    }
#line 971 "zgejsv.f"
/* L1952: */
#line 971 "zgejsv.f"
	}
#line 972 "zgejsv.f"
	i__1 = *m - 1;
#line 972 "zgejsv.f"
	zlaswp_(n, &a[a_offset], lda, &c__1, &i__1, &iwork[(*n << 1) + 1], &
		c__1);
#line 973 "zgejsv.f"
    }

/*     End of the preparation phase (scaling, optional sorting and */
/*     transposing, optional flushing of small columns). */

/*     Preconditioning */

/*     If the full SVD is needed, the right singular vectors are computed */
/*     from a matrix equation, and for that we need theoretical analysis */
/*     of the Businger-Golub pivoting. So we use ZGEQP3 as the first RR QRF. */
/*     In all other cases the first RR QRF can be chosen by other criteria */
/*     (eg speed by replacing global with restricted window pivoting, such */
/*     as in xGEQPX from TOMS # 782). Good results will be obtained using */
/*     xGEQPX with properly (!) chosen numerical parameters. */
/*     Any improvement of ZGEQP3 improves overal performance of ZGEJSV. */

/*     A * P1 = Q1 * [ R1^* 0]^*: */
#line 991 "zgejsv.f"
    i__1 = *n;
#line 991 "zgejsv.f"
    for (p = 1; p <= i__1; ++p) {
/*        .. all columns are free columns */
#line 993 "zgejsv.f"
	iwork[p] = 0;
#line 994 "zgejsv.f"
/* L1963: */
#line 994 "zgejsv.f"
    }
#line 995 "zgejsv.f"
    i__1 = *lwork - *n;
#line 995 "zgejsv.f"
    zgeqp3_(m, n, &a[a_offset], lda, &iwork[1], &cwork[1], &cwork[*n + 1], &
	    i__1, &rwork[1], &ierr);

/*     The upper triangular matrix R1 from the first QRF is inspected for */
/*     rank deficiency and possibilities for deflation, or possible */
/*     ill-conditioning. Depending on the user specified flag L2RANK, */
/*     the procedure explores possibilities to reduce the numerical */
/*     rank by inspecting the computed upper triangular factor. If */
/*     L2RANK or L2ABER are up, then ZGEJSV will compute the SVD of */
/*     A + dA, where ||dA|| <= f(M,N)*EPSLN. */

#line 1006 "zgejsv.f"
    nr = 1;
#line 1007 "zgejsv.f"
    if (l2aber) {
/*        Standard absolute error bound suffices. All sigma_i with */
/*        sigma_i < N*EPSLN*||A|| are flushed to zero. This is an */
/*        agressive enforcement of lower numerical rank by introducing a */
/*        backward error of the order of N*EPSLN*||A||. */
#line 1012 "zgejsv.f"
	temp1 = sqrt((doublereal) (*n)) * epsln;
#line 1013 "zgejsv.f"
	i__1 = *n;
#line 1013 "zgejsv.f"
	for (p = 2; p <= i__1; ++p) {
#line 1014 "zgejsv.f"
	    if (z_abs(&a[p + p * a_dim1]) >= temp1 * z_abs(&a[a_dim1 + 1])) {
#line 1015 "zgejsv.f"
		++nr;
#line 1016 "zgejsv.f"
	    } else {
#line 1017 "zgejsv.f"
		goto L3002;
#line 1018 "zgejsv.f"
	    }
#line 1019 "zgejsv.f"
/* L3001: */
#line 1019 "zgejsv.f"
	}
#line 1020 "zgejsv.f"
L3002:
#line 1021 "zgejsv.f"
	;
#line 1021 "zgejsv.f"
    } else if (l2rank) {
/*        .. similarly as above, only slightly more gentle (less agressive). */
/*        Sudden drop on the diagonal of R1 is used as the criterion for */
/*        close-to-rank-defficient. */
#line 1025 "zgejsv.f"
	temp1 = sqrt(sfmin);
#line 1026 "zgejsv.f"
	i__1 = *n;
#line 1026 "zgejsv.f"
	for (p = 2; p <= i__1; ++p) {
#line 1027 "zgejsv.f"
	    if (z_abs(&a[p + p * a_dim1]) < epsln * z_abs(&a[p - 1 + (p - 1) *
		     a_dim1]) || z_abs(&a[p + p * a_dim1]) < small || l2kill 
		    && z_abs(&a[p + p * a_dim1]) < temp1) {
#line 1027 "zgejsv.f"
		goto L3402;
#line 1027 "zgejsv.f"
	    }
#line 1030 "zgejsv.f"
	    ++nr;
#line 1031 "zgejsv.f"
/* L3401: */
#line 1031 "zgejsv.f"
	}
#line 1032 "zgejsv.f"
L3402:

#line 1034 "zgejsv.f"
	;
#line 1034 "zgejsv.f"
    } else {
/*        The goal is high relative accuracy. However, if the matrix */
/*        has high scaled condition number the relative accuracy is in */
/*        general not feasible. Later on, a condition number estimator */
/*        will be deployed to estimate the scaled condition number. */
/*        Here we just remove the underflowed part of the triangular */
/*        factor. This prevents the situation in which the code is */
/*        working hard to get the accuracy not warranted by the data. */
#line 1042 "zgejsv.f"
	temp1 = sqrt(sfmin);
#line 1043 "zgejsv.f"
	i__1 = *n;
#line 1043 "zgejsv.f"
	for (p = 2; p <= i__1; ++p) {
#line 1044 "zgejsv.f"
	    if (z_abs(&a[p + p * a_dim1]) < small || l2kill && z_abs(&a[p + p 
		    * a_dim1]) < temp1) {
#line 1044 "zgejsv.f"
		goto L3302;
#line 1044 "zgejsv.f"
	    }
#line 1046 "zgejsv.f"
	    ++nr;
#line 1047 "zgejsv.f"
/* L3301: */
#line 1047 "zgejsv.f"
	}
#line 1048 "zgejsv.f"
L3302:

#line 1050 "zgejsv.f"
	;
#line 1050 "zgejsv.f"
    }

#line 1052 "zgejsv.f"
    almort = FALSE_;
#line 1053 "zgejsv.f"
    if (nr == *n) {
#line 1054 "zgejsv.f"
	maxprj = 1.;
#line 1055 "zgejsv.f"
	i__1 = *n;
#line 1055 "zgejsv.f"
	for (p = 2; p <= i__1; ++p) {
#line 1056 "zgejsv.f"
	    temp1 = z_abs(&a[p + p * a_dim1]) / sva[iwork[p]];
#line 1057 "zgejsv.f"
	    maxprj = min(maxprj,temp1);
#line 1058 "zgejsv.f"
/* L3051: */
#line 1058 "zgejsv.f"
	}
/* Computing 2nd power */
#line 1059 "zgejsv.f"
	d__1 = maxprj;
#line 1059 "zgejsv.f"
	if (d__1 * d__1 >= 1. - (doublereal) (*n) * epsln) {
#line 1059 "zgejsv.f"
	    almort = TRUE_;
#line 1059 "zgejsv.f"
	}
#line 1060 "zgejsv.f"
    }


#line 1063 "zgejsv.f"
    sconda = -1.;
#line 1064 "zgejsv.f"
    condr1 = -1.;
#line 1065 "zgejsv.f"
    condr2 = -1.;

#line 1067 "zgejsv.f"
    if (errest) {
#line 1068 "zgejsv.f"
	if (*n == nr) {
#line 1069 "zgejsv.f"
	    if (rsvec) {
/*              .. V is available as workspace */
#line 1071 "zgejsv.f"
		zlacpy_("U", n, n, &a[a_offset], lda, &v[v_offset], ldv, (
			ftnlen)1);
#line 1072 "zgejsv.f"
		i__1 = *n;
#line 1072 "zgejsv.f"
		for (p = 1; p <= i__1; ++p) {
#line 1073 "zgejsv.f"
		    temp1 = sva[iwork[p]];
#line 1074 "zgejsv.f"
		    d__1 = 1. / temp1;
#line 1074 "zgejsv.f"
		    zdscal_(&p, &d__1, &v[p * v_dim1 + 1], &c__1);
#line 1075 "zgejsv.f"
/* L3053: */
#line 1075 "zgejsv.f"
		}
#line 1076 "zgejsv.f"
		zpocon_("U", n, &v[v_offset], ldv, &c_b80, &temp1, &cwork[*n 
			+ 1], &rwork[1], &ierr, (ftnlen)1);

#line 1079 "zgejsv.f"
	    } else if (lsvec) {
/*              .. U is available as workspace */
#line 1081 "zgejsv.f"
		zlacpy_("U", n, n, &a[a_offset], lda, &u[u_offset], ldu, (
			ftnlen)1);
#line 1082 "zgejsv.f"
		i__1 = *n;
#line 1082 "zgejsv.f"
		for (p = 1; p <= i__1; ++p) {
#line 1083 "zgejsv.f"
		    temp1 = sva[iwork[p]];
#line 1084 "zgejsv.f"
		    d__1 = 1. / temp1;
#line 1084 "zgejsv.f"
		    zdscal_(&p, &d__1, &u[p * u_dim1 + 1], &c__1);
#line 1085 "zgejsv.f"
/* L3054: */
#line 1085 "zgejsv.f"
		}
#line 1086 "zgejsv.f"
		zpocon_("U", n, &u[u_offset], ldu, &c_b80, &temp1, &cwork[*n 
			+ 1], &rwork[1], &ierr, (ftnlen)1);
#line 1088 "zgejsv.f"
	    } else {
#line 1089 "zgejsv.f"
		zlacpy_("U", n, n, &a[a_offset], lda, &cwork[*n + 1], n, (
			ftnlen)1);
#line 1090 "zgejsv.f"
		i__1 = *n;
#line 1090 "zgejsv.f"
		for (p = 1; p <= i__1; ++p) {
#line 1091 "zgejsv.f"
		    temp1 = sva[iwork[p]];
#line 1092 "zgejsv.f"
		    d__1 = 1. / temp1;
#line 1092 "zgejsv.f"
		    zdscal_(&p, &d__1, &cwork[*n + (p - 1) * *n + 1], &c__1);
#line 1093 "zgejsv.f"
/* L3052: */
#line 1093 "zgejsv.f"
		}
/*           .. the columns of R are scaled to have unit Euclidean lengths. */
#line 1095 "zgejsv.f"
		zpocon_("U", n, &cwork[*n + 1], n, &c_b80, &temp1, &cwork[*n 
			+ *n * *n + 1], &rwork[1], &ierr, (ftnlen)1);

#line 1098 "zgejsv.f"
	    }
#line 1099 "zgejsv.f"
	    sconda = 1. / sqrt(temp1);
/*           SCONDA is an estimate of SQRT(||(R^* * R)^(-1)||_1). */
/*           N^(-1/4) * SCONDA <= ||R^(-1)||_2 <= N^(1/4) * SCONDA */
#line 1102 "zgejsv.f"
	} else {
#line 1103 "zgejsv.f"
	    sconda = -1.;
#line 1104 "zgejsv.f"
	}
#line 1105 "zgejsv.f"
    }

#line 1107 "zgejsv.f"
    z_div(&z__1, &a[a_dim1 + 1], &a[nr + nr * a_dim1]);
#line 1107 "zgejsv.f"
    l2pert = l2pert && z_abs(&z__1) > sqrt(big1);
/*     If there is no violent scaling, artificial perturbation is not needed. */

/*     Phase 3: */

#line 1112 "zgejsv.f"
    if (! (rsvec || lsvec)) {

/*         Singular Values only */

/*         .. transpose A(1:NR,1:N) */
/* Computing MIN */
#line 1117 "zgejsv.f"
	i__2 = *n - 1;
#line 1117 "zgejsv.f"
	i__1 = min(i__2,nr);
#line 1117 "zgejsv.f"
	for (p = 1; p <= i__1; ++p) {
#line 1118 "zgejsv.f"
	    i__2 = *n - p;
#line 1118 "zgejsv.f"
	    zcopy_(&i__2, &a[p + (p + 1) * a_dim1], lda, &a[p + 1 + p * 
		    a_dim1], &c__1);
#line 1119 "zgejsv.f"
	    i__2 = *n - p + 1;
#line 1119 "zgejsv.f"
	    zlacgv_(&i__2, &a[p + p * a_dim1], &c__1);
#line 1120 "zgejsv.f"
/* L1946: */
#line 1120 "zgejsv.f"
	}
#line 1121 "zgejsv.f"
	if (nr == *n) {
#line 1121 "zgejsv.f"
	    i__1 = *n + *n * a_dim1;
#line 1121 "zgejsv.f"
	    d_cnjg(&z__1, &a[*n + *n * a_dim1]);
#line 1121 "zgejsv.f"
	    a[i__1].r = z__1.r, a[i__1].i = z__1.i;
#line 1121 "zgejsv.f"
	}

/*        The following two DO-loops introduce small relative perturbation */
/*        into the strict upper triangle of the lower triangular matrix. */
/*        Small entries below the main diagonal are also changed. */
/*        This modification is useful if the computing environment does not */
/*        provide/allow FLUSH TO ZERO underflow, for it prevents many */
/*        annoying denormalized numbers in case of strongly scaled matrices. */
/*        The perturbation is structured so that it does not introduce any */
/*        new perturbation of the singular values, and it does not destroy */
/*        the job done by the preconditioner. */
/*        The licence for this perturbation is in the variable L2PERT, which */
/*        should be .FALSE. if FLUSH TO ZERO underflow is active. */

#line 1135 "zgejsv.f"
	if (! almort) {

#line 1137 "zgejsv.f"
	    if (l2pert) {
/*              XSC = SQRT(SMALL) */
#line 1139 "zgejsv.f"
		xsc = epsln / (doublereal) (*n);
#line 1140 "zgejsv.f"
		i__1 = nr;
#line 1140 "zgejsv.f"
		for (q = 1; q <= i__1; ++q) {
#line 1141 "zgejsv.f"
		    d__1 = xsc * z_abs(&a[q + q * a_dim1]);
#line 1141 "zgejsv.f"
		    z__1.r = d__1, z__1.i = 0.;
#line 1141 "zgejsv.f"
		    ctemp.r = z__1.r, ctemp.i = z__1.i;
#line 1142 "zgejsv.f"
		    i__2 = *n;
#line 1142 "zgejsv.f"
		    for (p = 1; p <= i__2; ++p) {
#line 1143 "zgejsv.f"
			if (p > q && z_abs(&a[p + q * a_dim1]) <= temp1 || p <
				 q) {
#line 1143 "zgejsv.f"
			    i__3 = p + q * a_dim1;
#line 1143 "zgejsv.f"
			    a[i__3].r = ctemp.r, a[i__3].i = ctemp.i;
#line 1143 "zgejsv.f"
			}
/*     $                     A(p,q) = TEMP1 * ( A(p,q) / ABS(A(p,q)) ) */
#line 1147 "zgejsv.f"
/* L4949: */
#line 1147 "zgejsv.f"
		    }
#line 1148 "zgejsv.f"
/* L4947: */
#line 1148 "zgejsv.f"
		}
#line 1149 "zgejsv.f"
	    } else {
#line 1150 "zgejsv.f"
		i__1 = nr - 1;
#line 1150 "zgejsv.f"
		i__2 = nr - 1;
#line 1150 "zgejsv.f"
		zlaset_("U", &i__1, &i__2, &c_b1, &c_b1, &a[(a_dim1 << 1) + 1]
			, lda, (ftnlen)1);
#line 1151 "zgejsv.f"
	    }

/*            .. second preconditioning using the QR factorization */

#line 1155 "zgejsv.f"
	    i__1 = *lwork - *n;
#line 1155 "zgejsv.f"
	    zgeqrf_(n, &nr, &a[a_offset], lda, &cwork[1], &cwork[*n + 1], &
		    i__1, &ierr);

/*           .. and transpose upper to lower triangular */
#line 1158 "zgejsv.f"
	    i__1 = nr - 1;
#line 1158 "zgejsv.f"
	    for (p = 1; p <= i__1; ++p) {
#line 1159 "zgejsv.f"
		i__2 = nr - p;
#line 1159 "zgejsv.f"
		zcopy_(&i__2, &a[p + (p + 1) * a_dim1], lda, &a[p + 1 + p * 
			a_dim1], &c__1);
#line 1160 "zgejsv.f"
		i__2 = nr - p + 1;
#line 1160 "zgejsv.f"
		zlacgv_(&i__2, &a[p + p * a_dim1], &c__1);
#line 1161 "zgejsv.f"
/* L1948: */
#line 1161 "zgejsv.f"
	    }

#line 1163 "zgejsv.f"
	}

/*           Row-cyclic Jacobi SVD algorithm with column pivoting */

/*           .. again some perturbation (a "background noise") is added */
/*           to drown denormals */
#line 1169 "zgejsv.f"
	if (l2pert) {
/*              XSC = SQRT(SMALL) */
#line 1171 "zgejsv.f"
	    xsc = epsln / (doublereal) (*n);
#line 1172 "zgejsv.f"
	    i__1 = nr;
#line 1172 "zgejsv.f"
	    for (q = 1; q <= i__1; ++q) {
#line 1173 "zgejsv.f"
		d__1 = xsc * z_abs(&a[q + q * a_dim1]);
#line 1173 "zgejsv.f"
		z__1.r = d__1, z__1.i = 0.;
#line 1173 "zgejsv.f"
		ctemp.r = z__1.r, ctemp.i = z__1.i;
#line 1174 "zgejsv.f"
		i__2 = nr;
#line 1174 "zgejsv.f"
		for (p = 1; p <= i__2; ++p) {
#line 1175 "zgejsv.f"
		    if (p > q && z_abs(&a[p + q * a_dim1]) <= temp1 || p < q) 
			    {
#line 1175 "zgejsv.f"
			i__3 = p + q * a_dim1;
#line 1175 "zgejsv.f"
			a[i__3].r = ctemp.r, a[i__3].i = ctemp.i;
#line 1175 "zgejsv.f"
		    }
/*     $                   A(p,q) = TEMP1 * ( A(p,q) / ABS(A(p,q)) ) */
#line 1179 "zgejsv.f"
/* L1949: */
#line 1179 "zgejsv.f"
		}
#line 1180 "zgejsv.f"
/* L1947: */
#line 1180 "zgejsv.f"
	    }
#line 1181 "zgejsv.f"
	} else {
#line 1182 "zgejsv.f"
	    i__1 = nr - 1;
#line 1182 "zgejsv.f"
	    i__2 = nr - 1;
#line 1182 "zgejsv.f"
	    zlaset_("U", &i__1, &i__2, &c_b1, &c_b1, &a[(a_dim1 << 1) + 1], 
		    lda, (ftnlen)1);
#line 1183 "zgejsv.f"
	}

/*           .. and one-sided Jacobi rotations are started on a lower */
/*           triangular matrix (plus perturbation which is ignored in */
/*           the part which destroys triangular form (confusing?!)) */

#line 1189 "zgejsv.f"
	zgesvj_("L", "NoU", "NoV", &nr, &nr, &a[a_offset], lda, &sva[1], n, &
		v[v_offset], ldv, &cwork[1], lwork, &rwork[1], lrwork, info, (
		ftnlen)1, (ftnlen)3, (ftnlen)3);

#line 1192 "zgejsv.f"
	scalem = rwork[1];
#line 1193 "zgejsv.f"
	numrank = i_dnnt(&rwork[2]);


#line 1196 "zgejsv.f"
    } else if (rsvec && ! lsvec) {

/*        -> Singular Values and Right Singular Vectors <- */

#line 1200 "zgejsv.f"
	if (almort) {

/*           .. in this case NR equals N */
#line 1203 "zgejsv.f"
	    i__1 = nr;
#line 1203 "zgejsv.f"
	    for (p = 1; p <= i__1; ++p) {
#line 1204 "zgejsv.f"
		i__2 = *n - p + 1;
#line 1204 "zgejsv.f"
		zcopy_(&i__2, &a[p + p * a_dim1], lda, &v[p + p * v_dim1], &
			c__1);
#line 1205 "zgejsv.f"
		i__2 = *n - p + 1;
#line 1205 "zgejsv.f"
		zlacgv_(&i__2, &v[p + p * v_dim1], &c__1);
#line 1206 "zgejsv.f"
/* L1998: */
#line 1206 "zgejsv.f"
	    }
#line 1207 "zgejsv.f"
	    i__1 = nr - 1;
#line 1207 "zgejsv.f"
	    i__2 = nr - 1;
#line 1207 "zgejsv.f"
	    zlaset_("Upper", &i__1, &i__2, &c_b1, &c_b1, &v[(v_dim1 << 1) + 1]
		    , ldv, (ftnlen)5);

#line 1209 "zgejsv.f"
	    zgesvj_("L", "U", "N", n, &nr, &v[v_offset], ldv, &sva[1], &nr, &
		    a[a_offset], lda, &cwork[1], lwork, &rwork[1], lrwork, 
		    info, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1211 "zgejsv.f"
	    scalem = rwork[1];
#line 1212 "zgejsv.f"
	    numrank = i_dnnt(&rwork[2]);
#line 1214 "zgejsv.f"
	} else {

/*        .. two more QR factorizations ( one QRF is not enough, two require */
/*        accumulated product of Jacobi rotations, three are perfect ) */

#line 1219 "zgejsv.f"
	    i__1 = nr - 1;
#line 1219 "zgejsv.f"
	    i__2 = nr - 1;
#line 1219 "zgejsv.f"
	    zlaset_("Lower", &i__1, &i__2, &c_b1, &c_b1, &a[a_dim1 + 2], lda, 
		    (ftnlen)5);
#line 1220 "zgejsv.f"
	    i__1 = *lwork - *n;
#line 1220 "zgejsv.f"
	    zgelqf_(&nr, n, &a[a_offset], lda, &cwork[1], &cwork[*n + 1], &
		    i__1, &ierr);
#line 1221 "zgejsv.f"
	    zlacpy_("Lower", &nr, &nr, &a[a_offset], lda, &v[v_offset], ldv, (
		    ftnlen)5);
#line 1222 "zgejsv.f"
	    i__1 = nr - 1;
#line 1222 "zgejsv.f"
	    i__2 = nr - 1;
#line 1222 "zgejsv.f"
	    zlaset_("Upper", &i__1, &i__2, &c_b1, &c_b1, &v[(v_dim1 << 1) + 1]
		    , ldv, (ftnlen)5);
#line 1223 "zgejsv.f"
	    i__1 = *lwork - (*n << 1);
#line 1223 "zgejsv.f"
	    zgeqrf_(&nr, &nr, &v[v_offset], ldv, &cwork[*n + 1], &cwork[(*n <<
		     1) + 1], &i__1, &ierr);
#line 1225 "zgejsv.f"
	    i__1 = nr;
#line 1225 "zgejsv.f"
	    for (p = 1; p <= i__1; ++p) {
#line 1226 "zgejsv.f"
		i__2 = nr - p + 1;
#line 1226 "zgejsv.f"
		zcopy_(&i__2, &v[p + p * v_dim1], ldv, &v[p + p * v_dim1], &
			c__1);
#line 1227 "zgejsv.f"
		i__2 = nr - p + 1;
#line 1227 "zgejsv.f"
		zlacgv_(&i__2, &v[p + p * v_dim1], &c__1);
#line 1228 "zgejsv.f"
/* L8998: */
#line 1228 "zgejsv.f"
	    }
#line 1229 "zgejsv.f"
	    i__1 = nr - 1;
#line 1229 "zgejsv.f"
	    i__2 = nr - 1;
#line 1229 "zgejsv.f"
	    zlaset_("Upper", &i__1, &i__2, &c_b120, &c_b120, &v[(v_dim1 << 1) 
		    + 1], ldv, (ftnlen)5);

#line 1231 "zgejsv.f"
	    i__1 = *lwork - *n;
#line 1231 "zgejsv.f"
	    zgesvj_("Lower", "U", "N", &nr, &nr, &v[v_offset], ldv, &sva[1], &
		    nr, &u[u_offset], ldu, &cwork[*n + 1], &i__1, &rwork[1], 
		    lrwork, info, (ftnlen)5, (ftnlen)1, (ftnlen)1);
#line 1233 "zgejsv.f"
	    scalem = rwork[1];
#line 1234 "zgejsv.f"
	    numrank = i_dnnt(&rwork[2]);
#line 1235 "zgejsv.f"
	    if (nr < *n) {
#line 1236 "zgejsv.f"
		i__1 = *n - nr;
#line 1236 "zgejsv.f"
		zlaset_("A", &i__1, &nr, &c_b1, &c_b1, &v[nr + 1 + v_dim1], 
			ldv, (ftnlen)1);
#line 1237 "zgejsv.f"
		i__1 = *n - nr;
#line 1237 "zgejsv.f"
		zlaset_("A", &nr, &i__1, &c_b1, &c_b1, &v[(nr + 1) * v_dim1 + 
			1], ldv, (ftnlen)1);
#line 1238 "zgejsv.f"
		i__1 = *n - nr;
#line 1238 "zgejsv.f"
		i__2 = *n - nr;
#line 1238 "zgejsv.f"
		zlaset_("A", &i__1, &i__2, &c_b1, &c_b2, &v[nr + 1 + (nr + 1) 
			* v_dim1], ldv, (ftnlen)1);
#line 1239 "zgejsv.f"
	    }

#line 1241 "zgejsv.f"
	    i__1 = *lwork - *n;
#line 1241 "zgejsv.f"
	    zunmlq_("Left", "C", n, n, &nr, &a[a_offset], lda, &cwork[1], &v[
		    v_offset], ldv, &cwork[*n + 1], &i__1, &ierr, (ftnlen)4, (
		    ftnlen)1);

#line 1244 "zgejsv.f"
	}

#line 1246 "zgejsv.f"
	i__1 = *n;
#line 1246 "zgejsv.f"
	for (p = 1; p <= i__1; ++p) {
#line 1247 "zgejsv.f"
	    zcopy_(n, &v[p + v_dim1], ldv, &a[iwork[p] + a_dim1], lda);
#line 1248 "zgejsv.f"
/* L8991: */
#line 1248 "zgejsv.f"
	}
#line 1249 "zgejsv.f"
	zlacpy_("All", n, n, &a[a_offset], lda, &v[v_offset], ldv, (ftnlen)3);

#line 1251 "zgejsv.f"
	if (transp) {
#line 1252 "zgejsv.f"
	    zlacpy_("All", n, n, &v[v_offset], ldv, &u[u_offset], ldu, (
		    ftnlen)3);
#line 1253 "zgejsv.f"
	}

#line 1255 "zgejsv.f"
    } else if (lsvec && ! rsvec) {

/*        .. Singular Values and Left Singular Vectors                 .. */

/*        .. second preconditioning step to avoid need to accumulate */
/*        Jacobi rotations in the Jacobi iterations. */
#line 1261 "zgejsv.f"
	i__1 = nr;
#line 1261 "zgejsv.f"
	for (p = 1; p <= i__1; ++p) {
#line 1262 "zgejsv.f"
	    i__2 = *n - p + 1;
#line 1262 "zgejsv.f"
	    zcopy_(&i__2, &a[p + p * a_dim1], lda, &u[p + p * u_dim1], &c__1);
#line 1263 "zgejsv.f"
	    i__2 = *n - p + 1;
#line 1263 "zgejsv.f"
	    zlacgv_(&i__2, &u[p + p * u_dim1], &c__1);
#line 1264 "zgejsv.f"
/* L1965: */
#line 1264 "zgejsv.f"
	}
#line 1265 "zgejsv.f"
	i__1 = nr - 1;
#line 1265 "zgejsv.f"
	i__2 = nr - 1;
#line 1265 "zgejsv.f"
	zlaset_("Upper", &i__1, &i__2, &c_b1, &c_b1, &u[(u_dim1 << 1) + 1], 
		ldu, (ftnlen)5);

#line 1267 "zgejsv.f"
	i__1 = *lwork - (*n << 1);
#line 1267 "zgejsv.f"
	zgeqrf_(n, &nr, &u[u_offset], ldu, &cwork[*n + 1], &cwork[(*n << 1) + 
		1], &i__1, &ierr);

#line 1270 "zgejsv.f"
	i__1 = nr - 1;
#line 1270 "zgejsv.f"
	for (p = 1; p <= i__1; ++p) {
#line 1271 "zgejsv.f"
	    i__2 = nr - p;
#line 1271 "zgejsv.f"
	    zcopy_(&i__2, &u[p + (p + 1) * u_dim1], ldu, &u[p + 1 + p * 
		    u_dim1], &c__1);
#line 1272 "zgejsv.f"
	    i__2 = *n - p + 1;
#line 1272 "zgejsv.f"
	    zlacgv_(&i__2, &u[p + p * u_dim1], &c__1);
#line 1273 "zgejsv.f"
/* L1967: */
#line 1273 "zgejsv.f"
	}
#line 1274 "zgejsv.f"
	i__1 = nr - 1;
#line 1274 "zgejsv.f"
	i__2 = nr - 1;
#line 1274 "zgejsv.f"
	zlaset_("Upper", &i__1, &i__2, &c_b1, &c_b1, &u[(u_dim1 << 1) + 1], 
		ldu, (ftnlen)5);

#line 1276 "zgejsv.f"
	i__1 = *lwork - *n;
#line 1276 "zgejsv.f"
	zgesvj_("Lower", "U", "N", &nr, &nr, &u[u_offset], ldu, &sva[1], &nr, 
		&a[a_offset], lda, &cwork[*n + 1], &i__1, &rwork[1], lrwork, 
		info, (ftnlen)5, (ftnlen)1, (ftnlen)1);
#line 1278 "zgejsv.f"
	scalem = rwork[1];
#line 1279 "zgejsv.f"
	numrank = i_dnnt(&rwork[2]);

#line 1281 "zgejsv.f"
	if (nr < *m) {
#line 1282 "zgejsv.f"
	    i__1 = *m - nr;
#line 1282 "zgejsv.f"
	    zlaset_("A", &i__1, &nr, &c_b1, &c_b1, &u[nr + 1 + u_dim1], ldu, (
		    ftnlen)1);
#line 1283 "zgejsv.f"
	    if (nr < n1) {
#line 1284 "zgejsv.f"
		i__1 = n1 - nr;
#line 1284 "zgejsv.f"
		zlaset_("A", &nr, &i__1, &c_b1, &c_b1, &u[(nr + 1) * u_dim1 + 
			1], ldu, (ftnlen)1);
#line 1285 "zgejsv.f"
		i__1 = *m - nr;
#line 1285 "zgejsv.f"
		i__2 = n1 - nr;
#line 1285 "zgejsv.f"
		zlaset_("A", &i__1, &i__2, &c_b1, &c_b2, &u[nr + 1 + (nr + 1) 
			* u_dim1], ldu, (ftnlen)1);
#line 1286 "zgejsv.f"
	    }
#line 1287 "zgejsv.f"
	}

#line 1289 "zgejsv.f"
	i__1 = *lwork - *n;
#line 1289 "zgejsv.f"
	zunmqr_("Left", "No Tr", m, &n1, n, &a[a_offset], lda, &cwork[1], &u[
		u_offset], ldu, &cwork[*n + 1], &i__1, &ierr, (ftnlen)4, (
		ftnlen)5);

#line 1292 "zgejsv.f"
	if (rowpiv) {
#line 1292 "zgejsv.f"
	    i__1 = *m - 1;
#line 1292 "zgejsv.f"
	    zlaswp_(&n1, &u[u_offset], ldu, &c__1, &i__1, &iwork[(*n << 1) + 
		    1], &c_n1);
#line 1292 "zgejsv.f"
	}

#line 1295 "zgejsv.f"
	i__1 = n1;
#line 1295 "zgejsv.f"
	for (p = 1; p <= i__1; ++p) {
#line 1296 "zgejsv.f"
	    xsc = 1. / dznrm2_(m, &u[p * u_dim1 + 1], &c__1);
#line 1297 "zgejsv.f"
	    zdscal_(m, &xsc, &u[p * u_dim1 + 1], &c__1);
#line 1298 "zgejsv.f"
/* L1974: */
#line 1298 "zgejsv.f"
	}

#line 1300 "zgejsv.f"
	if (transp) {
#line 1301 "zgejsv.f"
	    zlacpy_("All", n, n, &u[u_offset], ldu, &v[v_offset], ldv, (
		    ftnlen)3);
#line 1302 "zgejsv.f"
	}

#line 1304 "zgejsv.f"
    } else {

/*        .. Full SVD .. */

#line 1308 "zgejsv.f"
	if (! jracc) {

#line 1310 "zgejsv.f"
	    if (! almort) {

/*           Second Preconditioning Step (QRF [with pivoting]) */
/*           Note that the composition of TRANSPOSE, QRF and TRANSPOSE is */
/*           equivalent to an LQF CALL. Since in many libraries the QRF */
/*           seems to be better optimized than the LQF, we do explicit */
/*           transpose and use the QRF. This is subject to changes in an */
/*           optimized implementation of ZGEJSV. */

#line 1319 "zgejsv.f"
		i__1 = nr;
#line 1319 "zgejsv.f"
		for (p = 1; p <= i__1; ++p) {
#line 1320 "zgejsv.f"
		    i__2 = *n - p + 1;
#line 1320 "zgejsv.f"
		    zcopy_(&i__2, &a[p + p * a_dim1], lda, &v[p + p * v_dim1],
			     &c__1);
#line 1321 "zgejsv.f"
		    i__2 = *n - p + 1;
#line 1321 "zgejsv.f"
		    zlacgv_(&i__2, &v[p + p * v_dim1], &c__1);
#line 1322 "zgejsv.f"
/* L1968: */
#line 1322 "zgejsv.f"
		}

/*           .. the following two loops perturb small entries to avoid */
/*           denormals in the second QR factorization, where they are */
/*           as good as zeros. This is done to avoid painfully slow */
/*           computation with denormals. The relative size of the perturbation */
/*           is a parameter that can be changed by the implementer. */
/*           This perturbation device will be obsolete on machines with */
/*           properly implemented arithmetic. */
/*           To switch it off, set L2PERT=.FALSE. To remove it from  the */
/*           code, remove the action under L2PERT=.TRUE., leave the ELSE part. */
/*           The following two loops should be blocked and fused with the */
/*           transposed copy above. */

#line 1336 "zgejsv.f"
		if (l2pert) {
#line 1337 "zgejsv.f"
		    xsc = sqrt(small);
#line 1338 "zgejsv.f"
		    i__1 = nr;
#line 1338 "zgejsv.f"
		    for (q = 1; q <= i__1; ++q) {
#line 1339 "zgejsv.f"
			d__1 = xsc * z_abs(&v[q + q * v_dim1]);
#line 1339 "zgejsv.f"
			z__1.r = d__1, z__1.i = 0.;
#line 1339 "zgejsv.f"
			ctemp.r = z__1.r, ctemp.i = z__1.i;
#line 1340 "zgejsv.f"
			i__2 = *n;
#line 1340 "zgejsv.f"
			for (p = 1; p <= i__2; ++p) {
#line 1341 "zgejsv.f"
			    if (p > q && z_abs(&v[p + q * v_dim1]) <= temp1 ||
				     p < q) {
#line 1341 "zgejsv.f"
				i__3 = p + q * v_dim1;
#line 1341 "zgejsv.f"
				v[i__3].r = ctemp.r, v[i__3].i = ctemp.i;
#line 1341 "zgejsv.f"
			    }
/*     $                   V(p,q) = TEMP1 * ( V(p,q) / ABS(V(p,q)) ) */
#line 1345 "zgejsv.f"
			    if (p < q) {
#line 1345 "zgejsv.f"
				i__3 = p + q * v_dim1;
#line 1345 "zgejsv.f"
				i__4 = p + q * v_dim1;
#line 1345 "zgejsv.f"
				z__1.r = -v[i__4].r, z__1.i = -v[i__4].i;
#line 1345 "zgejsv.f"
				v[i__3].r = z__1.r, v[i__3].i = z__1.i;
#line 1345 "zgejsv.f"
			    }
#line 1346 "zgejsv.f"
/* L2968: */
#line 1346 "zgejsv.f"
			}
#line 1347 "zgejsv.f"
/* L2969: */
#line 1347 "zgejsv.f"
		    }
#line 1348 "zgejsv.f"
		} else {
#line 1349 "zgejsv.f"
		    i__1 = nr - 1;
#line 1349 "zgejsv.f"
		    i__2 = nr - 1;
#line 1349 "zgejsv.f"
		    zlaset_("U", &i__1, &i__2, &c_b1, &c_b1, &v[(v_dim1 << 1) 
			    + 1], ldv, (ftnlen)1);
#line 1350 "zgejsv.f"
		}

/*           Estimate the row scaled condition number of R1 */
/*           (If R1 is rectangular, N > NR, then the condition number */
/*           of the leading NR x NR submatrix is estimated.) */

#line 1356 "zgejsv.f"
		zlacpy_("L", &nr, &nr, &v[v_offset], ldv, &cwork[(*n << 1) + 
			1], &nr, (ftnlen)1);
#line 1357 "zgejsv.f"
		i__1 = nr;
#line 1357 "zgejsv.f"
		for (p = 1; p <= i__1; ++p) {
#line 1358 "zgejsv.f"
		    i__2 = nr - p + 1;
#line 1358 "zgejsv.f"
		    temp1 = dznrm2_(&i__2, &cwork[(*n << 1) + (p - 1) * nr + 
			    p], &c__1);
#line 1359 "zgejsv.f"
		    i__2 = nr - p + 1;
#line 1359 "zgejsv.f"
		    d__1 = 1. / temp1;
#line 1359 "zgejsv.f"
		    zdscal_(&i__2, &d__1, &cwork[(*n << 1) + (p - 1) * nr + p]
			    , &c__1);
#line 1360 "zgejsv.f"
/* L3950: */
#line 1360 "zgejsv.f"
		}
#line 1361 "zgejsv.f"
		zpocon_("Lower", &nr, &cwork[(*n << 1) + 1], &nr, &c_b80, &
			temp1, &cwork[(*n << 1) + nr * nr + 1], &rwork[1], &
			ierr, (ftnlen)5);
#line 1363 "zgejsv.f"
		condr1 = 1. / sqrt(temp1);
/*           .. here need a second oppinion on the condition number */
/*           .. then assume worst case scenario */
/*           R1 is OK for inverse <=> CONDR1 .LT. DFLOAT(N) */
/*           more conservative    <=> CONDR1 .LT. SQRT(DFLOAT(N)) */

#line 1369 "zgejsv.f"
		cond_ok__ = sqrt(sqrt((doublereal) nr));
/* [TP]       COND_OK is a tuning parameter. */

#line 1372 "zgejsv.f"
		if (condr1 < cond_ok__) {
/*              .. the second QRF without pivoting. Note: in an optimized */
/*              implementation, this QRF should be implemented as the QRF */
/*              of a lower triangular matrix. */
/*              R1^* = Q2 * R2 */
#line 1377 "zgejsv.f"
		    i__1 = *lwork - (*n << 1);
#line 1377 "zgejsv.f"
		    zgeqrf_(n, &nr, &v[v_offset], ldv, &cwork[*n + 1], &cwork[
			    (*n << 1) + 1], &i__1, &ierr);

#line 1380 "zgejsv.f"
		    if (l2pert) {
#line 1381 "zgejsv.f"
			xsc = sqrt(small) / epsln;
#line 1382 "zgejsv.f"
			i__1 = nr;
#line 1382 "zgejsv.f"
			for (p = 2; p <= i__1; ++p) {
#line 1383 "zgejsv.f"
			    i__2 = p - 1;
#line 1383 "zgejsv.f"
			    for (q = 1; q <= i__2; ++q) {
/* Computing MIN */
#line 1384 "zgejsv.f"
				d__2 = z_abs(&v[p + p * v_dim1]), d__3 = 
					z_abs(&v[q + q * v_dim1]);
#line 1384 "zgejsv.f"
				d__1 = xsc * min(d__2,d__3);
#line 1384 "zgejsv.f"
				z__1.r = d__1, z__1.i = 0.;
#line 1384 "zgejsv.f"
				ctemp.r = z__1.r, ctemp.i = z__1.i;
#line 1386 "zgejsv.f"
				if (z_abs(&v[q + p * v_dim1]) <= temp1) {
#line 1386 "zgejsv.f"
				    i__3 = q + p * v_dim1;
#line 1386 "zgejsv.f"
				    v[i__3].r = ctemp.r, v[i__3].i = ctemp.i;
#line 1386 "zgejsv.f"
				}
/*     $                     V(q,p) = TEMP1 * ( V(q,p) / ABS(V(q,p)) ) */
#line 1389 "zgejsv.f"
/* L3958: */
#line 1389 "zgejsv.f"
			    }
#line 1390 "zgejsv.f"
/* L3959: */
#line 1390 "zgejsv.f"
			}
#line 1391 "zgejsv.f"
		    }

#line 1393 "zgejsv.f"
		    if (nr != *n) {
#line 1393 "zgejsv.f"
			zlacpy_("A", n, &nr, &v[v_offset], ldv, &cwork[(*n << 
				1) + 1], n, (ftnlen)1);
#line 1393 "zgejsv.f"
		    }
/*              .. save ... */

/*           .. this transposed copy should be better than naive */
#line 1398 "zgejsv.f"
		    i__1 = nr - 1;
#line 1398 "zgejsv.f"
		    for (p = 1; p <= i__1; ++p) {
#line 1399 "zgejsv.f"
			i__2 = nr - p;
#line 1399 "zgejsv.f"
			zcopy_(&i__2, &v[p + (p + 1) * v_dim1], ldv, &v[p + 1 
				+ p * v_dim1], &c__1);
#line 1400 "zgejsv.f"
			i__2 = nr - p + 1;
#line 1400 "zgejsv.f"
			zlacgv_(&i__2, &v[p + p * v_dim1], &c__1);
#line 1401 "zgejsv.f"
/* L1969: */
#line 1401 "zgejsv.f"
		    }
#line 1402 "zgejsv.f"
		    i__1 = nr + nr * v_dim1;
#line 1402 "zgejsv.f"
		    d_cnjg(&z__1, &v[nr + nr * v_dim1]);
#line 1402 "zgejsv.f"
		    v[i__1].r = z__1.r, v[i__1].i = z__1.i;

#line 1404 "zgejsv.f"
		    condr2 = condr1;

#line 1406 "zgejsv.f"
		} else {

/*              .. ill-conditioned case: second QRF with pivoting */
/*              Note that windowed pivoting would be equaly good */
/*              numerically, and more run-time efficient. So, in */
/*              an optimal implementation, the next call to ZGEQP3 */
/*              should be replaced with eg. CALL ZGEQPX (ACM TOMS #782) */
/*              with properly (carefully) chosen parameters. */

/*              R1^* * P2 = Q2 * R2 */
#line 1416 "zgejsv.f"
		    i__1 = nr;
#line 1416 "zgejsv.f"
		    for (p = 1; p <= i__1; ++p) {
#line 1417 "zgejsv.f"
			iwork[*n + p] = 0;
#line 1418 "zgejsv.f"
/* L3003: */
#line 1418 "zgejsv.f"
		    }
#line 1419 "zgejsv.f"
		    i__1 = *lwork - (*n << 1);
#line 1419 "zgejsv.f"
		    zgeqp3_(n, &nr, &v[v_offset], ldv, &iwork[*n + 1], &cwork[
			    *n + 1], &cwork[(*n << 1) + 1], &i__1, &rwork[1], 
			    &ierr);
/* *               CALL ZGEQRF( N, NR, V, LDV, CWORK(N+1), CWORK(2*N+1), */
/* *     $              LWORK-2*N, IERR ) */
#line 1423 "zgejsv.f"
		    if (l2pert) {
#line 1424 "zgejsv.f"
			xsc = sqrt(small);
#line 1425 "zgejsv.f"
			i__1 = nr;
#line 1425 "zgejsv.f"
			for (p = 2; p <= i__1; ++p) {
#line 1426 "zgejsv.f"
			    i__2 = p - 1;
#line 1426 "zgejsv.f"
			    for (q = 1; q <= i__2; ++q) {
/* Computing MIN */
#line 1427 "zgejsv.f"
				d__2 = z_abs(&v[p + p * v_dim1]), d__3 = 
					z_abs(&v[q + q * v_dim1]);
#line 1427 "zgejsv.f"
				d__1 = xsc * min(d__2,d__3);
#line 1427 "zgejsv.f"
				z__1.r = d__1, z__1.i = 0.;
#line 1427 "zgejsv.f"
				ctemp.r = z__1.r, ctemp.i = z__1.i;
#line 1429 "zgejsv.f"
				if (z_abs(&v[q + p * v_dim1]) <= temp1) {
#line 1429 "zgejsv.f"
				    i__3 = q + p * v_dim1;
#line 1429 "zgejsv.f"
				    v[i__3].r = ctemp.r, v[i__3].i = ctemp.i;
#line 1429 "zgejsv.f"
				}
/*     $                     V(q,p) = TEMP1 * ( V(q,p) / ABS(V(q,p)) ) */
#line 1432 "zgejsv.f"
/* L3968: */
#line 1432 "zgejsv.f"
			    }
#line 1433 "zgejsv.f"
/* L3969: */
#line 1433 "zgejsv.f"
			}
#line 1434 "zgejsv.f"
		    }

#line 1436 "zgejsv.f"
		    zlacpy_("A", n, &nr, &v[v_offset], ldv, &cwork[(*n << 1) 
			    + 1], n, (ftnlen)1);

#line 1438 "zgejsv.f"
		    if (l2pert) {
#line 1439 "zgejsv.f"
			xsc = sqrt(small);
#line 1440 "zgejsv.f"
			i__1 = nr;
#line 1440 "zgejsv.f"
			for (p = 2; p <= i__1; ++p) {
#line 1441 "zgejsv.f"
			    i__2 = p - 1;
#line 1441 "zgejsv.f"
			    for (q = 1; q <= i__2; ++q) {
/* Computing MIN */
#line 1442 "zgejsv.f"
				d__2 = z_abs(&v[p + p * v_dim1]), d__3 = 
					z_abs(&v[q + q * v_dim1]);
#line 1442 "zgejsv.f"
				d__1 = xsc * min(d__2,d__3);
#line 1442 "zgejsv.f"
				z__1.r = d__1, z__1.i = 0.;
#line 1442 "zgejsv.f"
				ctemp.r = z__1.r, ctemp.i = z__1.i;
/*                        V(p,q) = - TEMP1*( V(q,p) / ABS(V(q,p)) ) */
#line 1445 "zgejsv.f"
				i__3 = p + q * v_dim1;
#line 1445 "zgejsv.f"
				z__1.r = -ctemp.r, z__1.i = -ctemp.i;
#line 1445 "zgejsv.f"
				v[i__3].r = z__1.r, v[i__3].i = z__1.i;
#line 1446 "zgejsv.f"
/* L8971: */
#line 1446 "zgejsv.f"
			    }
#line 1447 "zgejsv.f"
/* L8970: */
#line 1447 "zgejsv.f"
			}
#line 1448 "zgejsv.f"
		    } else {
#line 1449 "zgejsv.f"
			i__1 = nr - 1;
#line 1449 "zgejsv.f"
			i__2 = nr - 1;
#line 1449 "zgejsv.f"
			zlaset_("L", &i__1, &i__2, &c_b1, &c_b1, &v[v_dim1 + 
				2], ldv, (ftnlen)1);
#line 1450 "zgejsv.f"
		    }
/*              Now, compute R2 = L3 * Q3, the LQ factorization. */
#line 1452 "zgejsv.f"
		    i__1 = *lwork - (*n << 1) - *n * nr - nr;
#line 1452 "zgejsv.f"
		    zgelqf_(&nr, &nr, &v[v_offset], ldv, &cwork[(*n << 1) + *
			    n * nr + 1], &cwork[(*n << 1) + *n * nr + nr + 1],
			     &i__1, &ierr);
/*              .. and estimate the condition number */
#line 1455 "zgejsv.f"
		    zlacpy_("L", &nr, &nr, &v[v_offset], ldv, &cwork[(*n << 1)
			     + *n * nr + nr + 1], &nr, (ftnlen)1);
#line 1456 "zgejsv.f"
		    i__1 = nr;
#line 1456 "zgejsv.f"
		    for (p = 1; p <= i__1; ++p) {
#line 1457 "zgejsv.f"
			temp1 = dznrm2_(&p, &cwork[(*n << 1) + *n * nr + nr + 
				p], &nr);
#line 1458 "zgejsv.f"
			d__1 = 1. / temp1;
#line 1458 "zgejsv.f"
			zdscal_(&p, &d__1, &cwork[(*n << 1) + *n * nr + nr + 
				p], &nr);
#line 1459 "zgejsv.f"
/* L4950: */
#line 1459 "zgejsv.f"
		    }
#line 1460 "zgejsv.f"
		    zpocon_("L", &nr, &cwork[(*n << 1) + *n * nr + nr + 1], &
			    nr, &c_b80, &temp1, &cwork[(*n << 1) + *n * nr + 
			    nr + nr * nr + 1], &rwork[1], &ierr, (ftnlen)1);
#line 1462 "zgejsv.f"
		    condr2 = 1. / sqrt(temp1);


#line 1465 "zgejsv.f"
		    if (condr2 >= cond_ok__) {
/*                 .. save the Householder vectors used for Q3 */
/*                 (this overwrittes the copy of R2, as it will not be */
/*                 needed in this branch, but it does not overwritte the */
/*                 Huseholder vectors of Q2.). */
#line 1470 "zgejsv.f"
			zlacpy_("U", &nr, &nr, &v[v_offset], ldv, &cwork[(*n 
				<< 1) + 1], n, (ftnlen)1);
/*                 .. and the rest of the information on Q3 is in */
/*                 WORK(2*N+N*NR+1:2*N+N*NR+N) */
#line 1473 "zgejsv.f"
		    }

#line 1475 "zgejsv.f"
		}

#line 1477 "zgejsv.f"
		if (l2pert) {
#line 1478 "zgejsv.f"
		    xsc = sqrt(small);
#line 1479 "zgejsv.f"
		    i__1 = nr;
#line 1479 "zgejsv.f"
		    for (q = 2; q <= i__1; ++q) {
#line 1480 "zgejsv.f"
			i__2 = q + q * v_dim1;
#line 1480 "zgejsv.f"
			z__1.r = xsc * v[i__2].r, z__1.i = xsc * v[i__2].i;
#line 1480 "zgejsv.f"
			ctemp.r = z__1.r, ctemp.i = z__1.i;
#line 1481 "zgejsv.f"
			i__2 = q - 1;
#line 1481 "zgejsv.f"
			for (p = 1; p <= i__2; ++p) {
/*                     V(p,q) = - TEMP1*( V(p,q) / ABS(V(p,q)) ) */
#line 1483 "zgejsv.f"
			    i__3 = p + q * v_dim1;
#line 1483 "zgejsv.f"
			    z__1.r = -ctemp.r, z__1.i = -ctemp.i;
#line 1483 "zgejsv.f"
			    v[i__3].r = z__1.r, v[i__3].i = z__1.i;
#line 1484 "zgejsv.f"
/* L4969: */
#line 1484 "zgejsv.f"
			}
#line 1485 "zgejsv.f"
/* L4968: */
#line 1485 "zgejsv.f"
		    }
#line 1486 "zgejsv.f"
		} else {
#line 1487 "zgejsv.f"
		    i__1 = nr - 1;
#line 1487 "zgejsv.f"
		    i__2 = nr - 1;
#line 1487 "zgejsv.f"
		    zlaset_("U", &i__1, &i__2, &c_b1, &c_b1, &v[(v_dim1 << 1) 
			    + 1], ldv, (ftnlen)1);
#line 1488 "zgejsv.f"
		}

/*        Second preconditioning finished; continue with Jacobi SVD */
/*        The input matrix is lower trinagular. */

/*        Recover the right singular vectors as solution of a well */
/*        conditioned triangular matrix equation. */

#line 1496 "zgejsv.f"
		if (condr1 < cond_ok__) {

#line 1498 "zgejsv.f"
		    i__1 = *lwork - (*n << 1) - *n * nr - nr;
#line 1498 "zgejsv.f"
		    zgesvj_("L", "U", "N", &nr, &nr, &v[v_offset], ldv, &sva[
			    1], &nr, &u[u_offset], ldu, &cwork[(*n << 1) + *n 
			    * nr + nr + 1], &i__1, &rwork[1], lrwork, info, (
			    ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1501 "zgejsv.f"
		    scalem = rwork[1];
#line 1502 "zgejsv.f"
		    numrank = i_dnnt(&rwork[2]);
#line 1503 "zgejsv.f"
		    i__1 = nr;
#line 1503 "zgejsv.f"
		    for (p = 1; p <= i__1; ++p) {
#line 1504 "zgejsv.f"
			zcopy_(&nr, &v[p * v_dim1 + 1], &c__1, &u[p * u_dim1 
				+ 1], &c__1);
#line 1505 "zgejsv.f"
			zdscal_(&nr, &sva[p], &v[p * v_dim1 + 1], &c__1);
#line 1506 "zgejsv.f"
/* L3970: */
#line 1506 "zgejsv.f"
		    }
/*        .. pick the right matrix equation and solve it */

#line 1510 "zgejsv.f"
		    if (nr == *n) {
/* :))             .. best case, R1 is inverted. The solution of this matrix */
/*                 equation is Q2*V2 = the product of the Jacobi rotations */
/*                 used in ZGESVJ, premultiplied with the orthogonal matrix */
/*                 from the second QR factorization. */
#line 1515 "zgejsv.f"
			ztrsm_("L", "U", "N", "N", &nr, &nr, &c_b2, &a[
				a_offset], lda, &v[v_offset], ldv, (ftnlen)1, 
				(ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1516 "zgejsv.f"
		    } else {
/*                 .. R1 is well conditioned, but non-square. Adjoint of R2 */
/*                 is inverted to get the product of the Jacobi rotations */
/*                 used in ZGESVJ. The Q-factor from the second QR */
/*                 factorization is then built in explicitly. */
#line 1521 "zgejsv.f"
			ztrsm_("L", "U", "C", "N", &nr, &nr, &c_b2, &cwork[(*
				n << 1) + 1], n, &v[v_offset], ldv, (ftnlen)1,
				 (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1523 "zgejsv.f"
			if (nr < *n) {
#line 1524 "zgejsv.f"
			    i__1 = *n - nr;
#line 1524 "zgejsv.f"
			    zlaset_("A", &i__1, &nr, &c_b120, &c_b1, &v[nr + 
				    1 + v_dim1], ldv, (ftnlen)1);
#line 1525 "zgejsv.f"
			    i__1 = *n - nr;
#line 1525 "zgejsv.f"
			    zlaset_("A", &nr, &i__1, &c_b120, &c_b1, &v[(nr + 
				    1) * v_dim1 + 1], ldv, (ftnlen)1);
#line 1526 "zgejsv.f"
			    i__1 = *n - nr;
#line 1526 "zgejsv.f"
			    i__2 = *n - nr;
#line 1526 "zgejsv.f"
			    zlaset_("A", &i__1, &i__2, &c_b120, &c_b2, &v[nr 
				    + 1 + (nr + 1) * v_dim1], ldv, (ftnlen)1);
#line 1527 "zgejsv.f"
			}
#line 1528 "zgejsv.f"
			i__1 = *lwork - (*n << 1) - *n * nr - nr;
#line 1528 "zgejsv.f"
			zunmqr_("L", "N", n, n, &nr, &cwork[(*n << 1) + 1], n,
				 &cwork[*n + 1], &v[v_offset], ldv, &cwork[(*
				n << 1) + *n * nr + nr + 1], &i__1, &ierr, (
				ftnlen)1, (ftnlen)1);
#line 1530 "zgejsv.f"
		    }

#line 1532 "zgejsv.f"
		} else if (condr2 < cond_ok__) {

/*              The matrix R2 is inverted. The solution of the matrix equation */
/*              is Q3^* * V3 = the product of the Jacobi rotations (appplied to */
/*              the lower triangular L3 from the LQ factorization of */
/*              R2=L3*Q3), pre-multiplied with the transposed Q3. */
#line 1538 "zgejsv.f"
		    i__1 = *lwork - (*n << 1) - *n * nr - nr;
#line 1538 "zgejsv.f"
		    zgesvj_("L", "U", "N", &nr, &nr, &v[v_offset], ldv, &sva[
			    1], &nr, &u[u_offset], ldu, &cwork[(*n << 1) + *n 
			    * nr + nr + 1], &i__1, &rwork[1], lrwork, info, (
			    ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1541 "zgejsv.f"
		    scalem = rwork[1];
#line 1542 "zgejsv.f"
		    numrank = i_dnnt(&rwork[2]);
#line 1543 "zgejsv.f"
		    i__1 = nr;
#line 1543 "zgejsv.f"
		    for (p = 1; p <= i__1; ++p) {
#line 1544 "zgejsv.f"
			zcopy_(&nr, &v[p * v_dim1 + 1], &c__1, &u[p * u_dim1 
				+ 1], &c__1);
#line 1545 "zgejsv.f"
			zdscal_(&nr, &sva[p], &u[p * u_dim1 + 1], &c__1);
#line 1546 "zgejsv.f"
/* L3870: */
#line 1546 "zgejsv.f"
		    }
#line 1547 "zgejsv.f"
		    ztrsm_("L", "U", "N", "N", &nr, &nr, &c_b2, &cwork[(*n << 
			    1) + 1], n, &u[u_offset], ldu, (ftnlen)1, (ftnlen)
			    1, (ftnlen)1, (ftnlen)1);
/*              .. apply the permutation from the second QR factorization */
#line 1550 "zgejsv.f"
		    i__1 = nr;
#line 1550 "zgejsv.f"
		    for (q = 1; q <= i__1; ++q) {
#line 1551 "zgejsv.f"
			i__2 = nr;
#line 1551 "zgejsv.f"
			for (p = 1; p <= i__2; ++p) {
#line 1552 "zgejsv.f"
			    i__3 = (*n << 1) + *n * nr + nr + iwork[*n + p];
#line 1552 "zgejsv.f"
			    i__4 = p + q * u_dim1;
#line 1552 "zgejsv.f"
			    cwork[i__3].r = u[i__4].r, cwork[i__3].i = u[i__4]
				    .i;
#line 1553 "zgejsv.f"
/* L872: */
#line 1553 "zgejsv.f"
			}
#line 1554 "zgejsv.f"
			i__2 = nr;
#line 1554 "zgejsv.f"
			for (p = 1; p <= i__2; ++p) {
#line 1555 "zgejsv.f"
			    i__3 = p + q * u_dim1;
#line 1555 "zgejsv.f"
			    i__4 = (*n << 1) + *n * nr + nr + p;
#line 1555 "zgejsv.f"
			    u[i__3].r = cwork[i__4].r, u[i__3].i = cwork[i__4]
				    .i;
#line 1556 "zgejsv.f"
/* L874: */
#line 1556 "zgejsv.f"
			}
#line 1557 "zgejsv.f"
/* L873: */
#line 1557 "zgejsv.f"
		    }
#line 1558 "zgejsv.f"
		    if (nr < *n) {
#line 1559 "zgejsv.f"
			i__1 = *n - nr;
#line 1559 "zgejsv.f"
			zlaset_("A", &i__1, &nr, &c_b1, &c_b1, &v[nr + 1 + 
				v_dim1], ldv, (ftnlen)1);
#line 1560 "zgejsv.f"
			i__1 = *n - nr;
#line 1560 "zgejsv.f"
			zlaset_("A", &nr, &i__1, &c_b1, &c_b1, &v[(nr + 1) * 
				v_dim1 + 1], ldv, (ftnlen)1);
#line 1561 "zgejsv.f"
			i__1 = *n - nr;
#line 1561 "zgejsv.f"
			i__2 = *n - nr;
#line 1561 "zgejsv.f"
			zlaset_("A", &i__1, &i__2, &c_b1, &c_b2, &v[nr + 1 + (
				nr + 1) * v_dim1], ldv, (ftnlen)1);
#line 1562 "zgejsv.f"
		    }
#line 1563 "zgejsv.f"
		    i__1 = *lwork - (*n << 1) - *n * nr - nr;
#line 1563 "zgejsv.f"
		    zunmqr_("L", "N", n, n, &nr, &cwork[(*n << 1) + 1], n, &
			    cwork[*n + 1], &v[v_offset], ldv, &cwork[(*n << 1)
			     + *n * nr + nr + 1], &i__1, &ierr, (ftnlen)1, (
			    ftnlen)1);
#line 1565 "zgejsv.f"
		} else {
/*              Last line of defense. */
/* #:(          This is a rather pathological case: no scaled condition */
/*              improvement after two pivoted QR factorizations. Other */
/*              possibility is that the rank revealing QR factorization */
/*              or the condition estimator has failed, or the COND_OK */
/*              is set very close to ONE (which is unnecessary). Normally, */
/*              this branch should never be executed, but in rare cases of */
/*              failure of the RRQR or condition estimator, the last line of */
/*              defense ensures that ZGEJSV completes the task. */
/*              Compute the full SVD of L3 using ZGESVJ with explicit */
/*              accumulation of Jacobi rotations. */
#line 1577 "zgejsv.f"
		    i__1 = *lwork - (*n << 1) - *n * nr - nr;
#line 1577 "zgejsv.f"
		    zgesvj_("L", "U", "V", &nr, &nr, &v[v_offset], ldv, &sva[
			    1], &nr, &u[u_offset], ldu, &cwork[(*n << 1) + *n 
			    * nr + nr + 1], &i__1, &rwork[1], lrwork, info, (
			    ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1580 "zgejsv.f"
		    scalem = rwork[1];
#line 1581 "zgejsv.f"
		    numrank = i_dnnt(&rwork[2]);
#line 1582 "zgejsv.f"
		    if (nr < *n) {
#line 1583 "zgejsv.f"
			i__1 = *n - nr;
#line 1583 "zgejsv.f"
			zlaset_("A", &i__1, &nr, &c_b1, &c_b1, &v[nr + 1 + 
				v_dim1], ldv, (ftnlen)1);
#line 1584 "zgejsv.f"
			i__1 = *n - nr;
#line 1584 "zgejsv.f"
			zlaset_("A", &nr, &i__1, &c_b1, &c_b1, &v[(nr + 1) * 
				v_dim1 + 1], ldv, (ftnlen)1);
#line 1585 "zgejsv.f"
			i__1 = *n - nr;
#line 1585 "zgejsv.f"
			i__2 = *n - nr;
#line 1585 "zgejsv.f"
			zlaset_("A", &i__1, &i__2, &c_b1, &c_b2, &v[nr + 1 + (
				nr + 1) * v_dim1], ldv, (ftnlen)1);
#line 1586 "zgejsv.f"
		    }
#line 1587 "zgejsv.f"
		    i__1 = *lwork - (*n << 1) - *n * nr - nr;
#line 1587 "zgejsv.f"
		    zunmqr_("L", "N", n, n, &nr, &cwork[(*n << 1) + 1], n, &
			    cwork[*n + 1], &v[v_offset], ldv, &cwork[(*n << 1)
			     + *n * nr + nr + 1], &i__1, &ierr, (ftnlen)1, (
			    ftnlen)1);

#line 1590 "zgejsv.f"
		    i__1 = *lwork - (*n << 1) - *n * nr - nr;
#line 1590 "zgejsv.f"
		    zunmlq_("L", "C", &nr, &nr, &nr, &cwork[(*n << 1) + 1], n,
			     &cwork[(*n << 1) + *n * nr + 1], &u[u_offset], 
			    ldu, &cwork[(*n << 1) + *n * nr + nr + 1], &i__1, 
			    &ierr, (ftnlen)1, (ftnlen)1);
#line 1593 "zgejsv.f"
		    i__1 = nr;
#line 1593 "zgejsv.f"
		    for (q = 1; q <= i__1; ++q) {
#line 1594 "zgejsv.f"
			i__2 = nr;
#line 1594 "zgejsv.f"
			for (p = 1; p <= i__2; ++p) {
#line 1595 "zgejsv.f"
			    i__3 = (*n << 1) + *n * nr + nr + iwork[*n + p];
#line 1595 "zgejsv.f"
			    i__4 = p + q * u_dim1;
#line 1595 "zgejsv.f"
			    cwork[i__3].r = u[i__4].r, cwork[i__3].i = u[i__4]
				    .i;
#line 1596 "zgejsv.f"
/* L772: */
#line 1596 "zgejsv.f"
			}
#line 1597 "zgejsv.f"
			i__2 = nr;
#line 1597 "zgejsv.f"
			for (p = 1; p <= i__2; ++p) {
#line 1598 "zgejsv.f"
			    i__3 = p + q * u_dim1;
#line 1598 "zgejsv.f"
			    i__4 = (*n << 1) + *n * nr + nr + p;
#line 1598 "zgejsv.f"
			    u[i__3].r = cwork[i__4].r, u[i__3].i = cwork[i__4]
				    .i;
#line 1599 "zgejsv.f"
/* L774: */
#line 1599 "zgejsv.f"
			}
#line 1600 "zgejsv.f"
/* L773: */
#line 1600 "zgejsv.f"
		    }

#line 1602 "zgejsv.f"
		}

/*           Permute the rows of V using the (column) permutation from the */
/*           first QRF. Also, scale the columns to make them unit in */
/*           Euclidean norm. This applies to all cases. */

#line 1608 "zgejsv.f"
		temp1 = sqrt((doublereal) (*n)) * epsln;
#line 1609 "zgejsv.f"
		i__1 = *n;
#line 1609 "zgejsv.f"
		for (q = 1; q <= i__1; ++q) {
#line 1610 "zgejsv.f"
		    i__2 = *n;
#line 1610 "zgejsv.f"
		    for (p = 1; p <= i__2; ++p) {
#line 1611 "zgejsv.f"
			i__3 = (*n << 1) + *n * nr + nr + iwork[p];
#line 1611 "zgejsv.f"
			i__4 = p + q * v_dim1;
#line 1611 "zgejsv.f"
			cwork[i__3].r = v[i__4].r, cwork[i__3].i = v[i__4].i;
#line 1612 "zgejsv.f"
/* L972: */
#line 1612 "zgejsv.f"
		    }
#line 1613 "zgejsv.f"
		    i__2 = *n;
#line 1613 "zgejsv.f"
		    for (p = 1; p <= i__2; ++p) {
#line 1614 "zgejsv.f"
			i__3 = p + q * v_dim1;
#line 1614 "zgejsv.f"
			i__4 = (*n << 1) + *n * nr + nr + p;
#line 1614 "zgejsv.f"
			v[i__3].r = cwork[i__4].r, v[i__3].i = cwork[i__4].i;
#line 1615 "zgejsv.f"
/* L973: */
#line 1615 "zgejsv.f"
		    }
#line 1616 "zgejsv.f"
		    xsc = 1. / dznrm2_(n, &v[q * v_dim1 + 1], &c__1);
#line 1617 "zgejsv.f"
		    if (xsc < 1. - temp1 || xsc > temp1 + 1.) {
#line 1617 "zgejsv.f"
			zdscal_(n, &xsc, &v[q * v_dim1 + 1], &c__1);
#line 1617 "zgejsv.f"
		    }
#line 1619 "zgejsv.f"
/* L1972: */
#line 1619 "zgejsv.f"
		}
/*           At this moment, V contains the right singular vectors of A. */
/*           Next, assemble the left singular vector matrix U (M x N). */
#line 1622 "zgejsv.f"
		if (nr < *m) {
#line 1623 "zgejsv.f"
		    i__1 = *m - nr;
#line 1623 "zgejsv.f"
		    zlaset_("A", &i__1, &nr, &c_b1, &c_b1, &u[nr + 1 + u_dim1]
			    , ldu, (ftnlen)1);
#line 1624 "zgejsv.f"
		    if (nr < n1) {
#line 1625 "zgejsv.f"
			i__1 = n1 - nr;
#line 1625 "zgejsv.f"
			zlaset_("A", &nr, &i__1, &c_b1, &c_b1, &u[(nr + 1) * 
				u_dim1 + 1], ldu, (ftnlen)1);
#line 1626 "zgejsv.f"
			i__1 = *m - nr;
#line 1626 "zgejsv.f"
			i__2 = n1 - nr;
#line 1626 "zgejsv.f"
			zlaset_("A", &i__1, &i__2, &c_b1, &c_b2, &u[nr + 1 + (
				nr + 1) * u_dim1], ldu, (ftnlen)1);
#line 1628 "zgejsv.f"
		    }
#line 1629 "zgejsv.f"
		}

/*           The Q matrix from the first QRF is built into the left singular */
/*           matrix U. This applies to all cases. */

#line 1634 "zgejsv.f"
		i__1 = *lwork - *n;
#line 1634 "zgejsv.f"
		zunmqr_("Left", "No_Tr", m, &n1, n, &a[a_offset], lda, &cwork[
			1], &u[u_offset], ldu, &cwork[*n + 1], &i__1, &ierr, (
			ftnlen)4, (ftnlen)5);
/*           The columns of U are normalized. The cost is O(M*N) flops. */
#line 1638 "zgejsv.f"
		temp1 = sqrt((doublereal) (*m)) * epsln;
#line 1639 "zgejsv.f"
		i__1 = nr;
#line 1639 "zgejsv.f"
		for (p = 1; p <= i__1; ++p) {
#line 1640 "zgejsv.f"
		    xsc = 1. / dznrm2_(m, &u[p * u_dim1 + 1], &c__1);
#line 1641 "zgejsv.f"
		    if (xsc < 1. - temp1 || xsc > temp1 + 1.) {
#line 1641 "zgejsv.f"
			zdscal_(m, &xsc, &u[p * u_dim1 + 1], &c__1);
#line 1641 "zgejsv.f"
		    }
#line 1643 "zgejsv.f"
/* L1973: */
#line 1643 "zgejsv.f"
		}

/*           If the initial QRF is computed with row pivoting, the left */
/*           singular vectors must be adjusted. */

#line 1648 "zgejsv.f"
		if (rowpiv) {
#line 1648 "zgejsv.f"
		    i__1 = *m - 1;
#line 1648 "zgejsv.f"
		    zlaswp_(&n1, &u[u_offset], ldu, &c__1, &i__1, &iwork[(*n 
			    << 1) + 1], &c_n1);
#line 1648 "zgejsv.f"
		}

#line 1651 "zgejsv.f"
	    } else {

/*        .. the initial matrix A has almost orthogonal columns and */
/*        the second QRF is not needed */

#line 1656 "zgejsv.f"
		zlacpy_("Upper", n, n, &a[a_offset], lda, &cwork[*n + 1], n, (
			ftnlen)5);
#line 1657 "zgejsv.f"
		if (l2pert) {
#line 1658 "zgejsv.f"
		    xsc = sqrt(small);
#line 1659 "zgejsv.f"
		    i__1 = *n;
#line 1659 "zgejsv.f"
		    for (p = 2; p <= i__1; ++p) {
#line 1660 "zgejsv.f"
			i__2 = *n + (p - 1) * *n + p;
#line 1660 "zgejsv.f"
			z__1.r = xsc * cwork[i__2].r, z__1.i = xsc * cwork[
				i__2].i;
#line 1660 "zgejsv.f"
			ctemp.r = z__1.r, ctemp.i = z__1.i;
#line 1661 "zgejsv.f"
			i__2 = p - 1;
#line 1661 "zgejsv.f"
			for (q = 1; q <= i__2; ++q) {
/*                     CWORK(N+(q-1)*N+p)=-TEMP1 * ( CWORK(N+(p-1)*N+q) / */
/*     $                                        ABS(CWORK(N+(p-1)*N+q)) ) */
#line 1664 "zgejsv.f"
			    i__3 = *n + (q - 1) * *n + p;
#line 1664 "zgejsv.f"
			    z__1.r = -ctemp.r, z__1.i = -ctemp.i;
#line 1664 "zgejsv.f"
			    cwork[i__3].r = z__1.r, cwork[i__3].i = z__1.i;
#line 1665 "zgejsv.f"
/* L5971: */
#line 1665 "zgejsv.f"
			}
#line 1666 "zgejsv.f"
/* L5970: */
#line 1666 "zgejsv.f"
		    }
#line 1667 "zgejsv.f"
		} else {
#line 1668 "zgejsv.f"
		    i__1 = *n - 1;
#line 1668 "zgejsv.f"
		    i__2 = *n - 1;
#line 1668 "zgejsv.f"
		    zlaset_("Lower", &i__1, &i__2, &c_b1, &c_b1, &cwork[*n + 
			    2], n, (ftnlen)5);
#line 1669 "zgejsv.f"
		}

#line 1671 "zgejsv.f"
		i__1 = *lwork - *n - *n * *n;
#line 1671 "zgejsv.f"
		zgesvj_("Upper", "U", "N", n, n, &cwork[*n + 1], n, &sva[1], 
			n, &u[u_offset], ldu, &cwork[*n + *n * *n + 1], &i__1,
			 &rwork[1], lrwork, info, (ftnlen)5, (ftnlen)1, (
			ftnlen)1);

#line 1675 "zgejsv.f"
		scalem = rwork[1];
#line 1676 "zgejsv.f"
		numrank = i_dnnt(&rwork[2]);
#line 1677 "zgejsv.f"
		i__1 = *n;
#line 1677 "zgejsv.f"
		for (p = 1; p <= i__1; ++p) {
#line 1678 "zgejsv.f"
		    zcopy_(n, &cwork[*n + (p - 1) * *n + 1], &c__1, &u[p * 
			    u_dim1 + 1], &c__1);
#line 1679 "zgejsv.f"
		    zdscal_(n, &sva[p], &cwork[*n + (p - 1) * *n + 1], &c__1);
#line 1680 "zgejsv.f"
/* L6970: */
#line 1680 "zgejsv.f"
		}

#line 1682 "zgejsv.f"
		ztrsm_("Left", "Upper", "NoTrans", "No UD", n, n, &c_b2, &a[
			a_offset], lda, &cwork[*n + 1], n, (ftnlen)4, (ftnlen)
			5, (ftnlen)7, (ftnlen)5);
#line 1684 "zgejsv.f"
		i__1 = *n;
#line 1684 "zgejsv.f"
		for (p = 1; p <= i__1; ++p) {
#line 1685 "zgejsv.f"
		    zcopy_(n, &cwork[*n + p], n, &v[iwork[p] + v_dim1], ldv);
#line 1686 "zgejsv.f"
/* L6972: */
#line 1686 "zgejsv.f"
		}
#line 1687 "zgejsv.f"
		temp1 = sqrt((doublereal) (*n)) * epsln;
#line 1688 "zgejsv.f"
		i__1 = *n;
#line 1688 "zgejsv.f"
		for (p = 1; p <= i__1; ++p) {
#line 1689 "zgejsv.f"
		    xsc = 1. / dznrm2_(n, &v[p * v_dim1 + 1], &c__1);
#line 1690 "zgejsv.f"
		    if (xsc < 1. - temp1 || xsc > temp1 + 1.) {
#line 1690 "zgejsv.f"
			zdscal_(n, &xsc, &v[p * v_dim1 + 1], &c__1);
#line 1690 "zgejsv.f"
		    }
#line 1692 "zgejsv.f"
/* L6971: */
#line 1692 "zgejsv.f"
		}

/*           Assemble the left singular vector matrix U (M x N). */

#line 1696 "zgejsv.f"
		if (*n < *m) {
#line 1697 "zgejsv.f"
		    i__1 = *m - *n;
#line 1697 "zgejsv.f"
		    zlaset_("A", &i__1, n, &c_b1, &c_b1, &u[*n + 1 + u_dim1], 
			    ldu, (ftnlen)1);
#line 1698 "zgejsv.f"
		    if (*n < n1) {
#line 1699 "zgejsv.f"
			i__1 = n1 - *n;
#line 1699 "zgejsv.f"
			zlaset_("A", n, &i__1, &c_b1, &c_b1, &u[(*n + 1) * 
				u_dim1 + 1], ldu, (ftnlen)1);
#line 1700 "zgejsv.f"
			i__1 = *m - *n;
#line 1700 "zgejsv.f"
			i__2 = n1 - *n;
#line 1700 "zgejsv.f"
			zlaset_("A", &i__1, &i__2, &c_b1, &c_b2, &u[*n + 1 + (
				*n + 1) * u_dim1], ldu, (ftnlen)1);
#line 1701 "zgejsv.f"
		    }
#line 1702 "zgejsv.f"
		}
#line 1703 "zgejsv.f"
		i__1 = *lwork - *n;
#line 1703 "zgejsv.f"
		zunmqr_("Left", "No Tr", m, &n1, n, &a[a_offset], lda, &cwork[
			1], &u[u_offset], ldu, &cwork[*n + 1], &i__1, &ierr, (
			ftnlen)4, (ftnlen)5);
#line 1705 "zgejsv.f"
		temp1 = sqrt((doublereal) (*m)) * epsln;
#line 1706 "zgejsv.f"
		i__1 = n1;
#line 1706 "zgejsv.f"
		for (p = 1; p <= i__1; ++p) {
#line 1707 "zgejsv.f"
		    xsc = 1. / dznrm2_(m, &u[p * u_dim1 + 1], &c__1);
#line 1708 "zgejsv.f"
		    if (xsc < 1. - temp1 || xsc > temp1 + 1.) {
#line 1708 "zgejsv.f"
			zdscal_(m, &xsc, &u[p * u_dim1 + 1], &c__1);
#line 1708 "zgejsv.f"
		    }
#line 1710 "zgejsv.f"
/* L6973: */
#line 1710 "zgejsv.f"
		}

#line 1712 "zgejsv.f"
		if (rowpiv) {
#line 1712 "zgejsv.f"
		    i__1 = *m - 1;
#line 1712 "zgejsv.f"
		    zlaswp_(&n1, &u[u_offset], ldu, &c__1, &i__1, &iwork[(*n 
			    << 1) + 1], &c_n1);
#line 1712 "zgejsv.f"
		}

#line 1715 "zgejsv.f"
	    }

/*        end of the  >> almost orthogonal case <<  in the full SVD */

#line 1719 "zgejsv.f"
	} else {

/*        This branch deploys a preconditioned Jacobi SVD with explicitly */
/*        accumulated rotations. It is included as optional, mainly for */
/*        experimental purposes. It does perfom well, and can also be used. */
/*        In this implementation, this branch will be automatically activated */
/*        if the  condition number sigma_max(A) / sigma_min(A) is predicted */
/*        to be greater than the overflow threshold. This is because the */
/*        a posteriori computation of the singular vectors assumes robust */
/*        implementation of BLAS and some LAPACK procedures, capable of working */
/*        in presence of extreme values. Since that is not always the case, ... */

#line 1731 "zgejsv.f"
	    i__1 = nr;
#line 1731 "zgejsv.f"
	    for (p = 1; p <= i__1; ++p) {
#line 1732 "zgejsv.f"
		i__2 = *n - p + 1;
#line 1732 "zgejsv.f"
		zcopy_(&i__2, &a[p + p * a_dim1], lda, &v[p + p * v_dim1], &
			c__1);
#line 1733 "zgejsv.f"
		i__2 = *n - p + 1;
#line 1733 "zgejsv.f"
		zlacgv_(&i__2, &v[p + p * v_dim1], &c__1);
#line 1734 "zgejsv.f"
/* L7968: */
#line 1734 "zgejsv.f"
	    }

#line 1736 "zgejsv.f"
	    if (l2pert) {
#line 1737 "zgejsv.f"
		xsc = sqrt(small / epsln);
#line 1738 "zgejsv.f"
		i__1 = nr;
#line 1738 "zgejsv.f"
		for (q = 1; q <= i__1; ++q) {
#line 1739 "zgejsv.f"
		    d__1 = xsc * z_abs(&v[q + q * v_dim1]);
#line 1739 "zgejsv.f"
		    z__1.r = d__1, z__1.i = 0.;
#line 1739 "zgejsv.f"
		    ctemp.r = z__1.r, ctemp.i = z__1.i;
#line 1740 "zgejsv.f"
		    i__2 = *n;
#line 1740 "zgejsv.f"
		    for (p = 1; p <= i__2; ++p) {
#line 1741 "zgejsv.f"
			if (p > q && z_abs(&v[p + q * v_dim1]) <= temp1 || p <
				 q) {
#line 1741 "zgejsv.f"
			    i__3 = p + q * v_dim1;
#line 1741 "zgejsv.f"
			    v[i__3].r = ctemp.r, v[i__3].i = ctemp.i;
#line 1741 "zgejsv.f"
			}
/*     $                V(p,q) = TEMP1 * ( V(p,q) / ABS(V(p,q)) ) */
#line 1745 "zgejsv.f"
			if (p < q) {
#line 1745 "zgejsv.f"
			    i__3 = p + q * v_dim1;
#line 1745 "zgejsv.f"
			    i__4 = p + q * v_dim1;
#line 1745 "zgejsv.f"
			    z__1.r = -v[i__4].r, z__1.i = -v[i__4].i;
#line 1745 "zgejsv.f"
			    v[i__3].r = z__1.r, v[i__3].i = z__1.i;
#line 1745 "zgejsv.f"
			}
#line 1746 "zgejsv.f"
/* L5968: */
#line 1746 "zgejsv.f"
		    }
#line 1747 "zgejsv.f"
/* L5969: */
#line 1747 "zgejsv.f"
		}
#line 1748 "zgejsv.f"
	    } else {
#line 1749 "zgejsv.f"
		i__1 = nr - 1;
#line 1749 "zgejsv.f"
		i__2 = nr - 1;
#line 1749 "zgejsv.f"
		zlaset_("U", &i__1, &i__2, &c_b1, &c_b1, &v[(v_dim1 << 1) + 1]
			, ldv, (ftnlen)1);
#line 1750 "zgejsv.f"
	    }
#line 1752 "zgejsv.f"
	    i__1 = *lwork - (*n << 1);
#line 1752 "zgejsv.f"
	    zgeqrf_(n, &nr, &v[v_offset], ldv, &cwork[*n + 1], &cwork[(*n << 
		    1) + 1], &i__1, &ierr);
#line 1754 "zgejsv.f"
	    zlacpy_("L", n, &nr, &v[v_offset], ldv, &cwork[(*n << 1) + 1], n, 
		    (ftnlen)1);

#line 1756 "zgejsv.f"
	    i__1 = nr;
#line 1756 "zgejsv.f"
	    for (p = 1; p <= i__1; ++p) {
#line 1757 "zgejsv.f"
		i__2 = nr - p + 1;
#line 1757 "zgejsv.f"
		zcopy_(&i__2, &v[p + p * v_dim1], ldv, &u[p + p * u_dim1], &
			c__1);
#line 1758 "zgejsv.f"
		i__2 = nr - p + 1;
#line 1758 "zgejsv.f"
		zlacgv_(&i__2, &u[p + p * u_dim1], &c__1);
#line 1759 "zgejsv.f"
/* L7969: */
#line 1759 "zgejsv.f"
	    }
#line 1761 "zgejsv.f"
	    if (l2pert) {
#line 1762 "zgejsv.f"
		xsc = sqrt(small / epsln);
#line 1763 "zgejsv.f"
		i__1 = nr;
#line 1763 "zgejsv.f"
		for (q = 2; q <= i__1; ++q) {
#line 1764 "zgejsv.f"
		    i__2 = q - 1;
#line 1764 "zgejsv.f"
		    for (p = 1; p <= i__2; ++p) {
/* Computing MIN */
#line 1765 "zgejsv.f"
			d__2 = z_abs(&u[p + p * u_dim1]), d__3 = z_abs(&u[q + 
				q * u_dim1]);
#line 1765 "zgejsv.f"
			d__1 = xsc * min(d__2,d__3);
#line 1765 "zgejsv.f"
			z__1.r = d__1, z__1.i = 0.;
#line 1765 "zgejsv.f"
			ctemp.r = z__1.r, ctemp.i = z__1.i;
/*                  U(p,q) = - TEMP1 * ( U(q,p) / ABS(U(q,p)) ) */
#line 1768 "zgejsv.f"
			i__3 = p + q * u_dim1;
#line 1768 "zgejsv.f"
			z__1.r = -ctemp.r, z__1.i = -ctemp.i;
#line 1768 "zgejsv.f"
			u[i__3].r = z__1.r, u[i__3].i = z__1.i;
#line 1769 "zgejsv.f"
/* L9971: */
#line 1769 "zgejsv.f"
		    }
#line 1770 "zgejsv.f"
/* L9970: */
#line 1770 "zgejsv.f"
		}
#line 1771 "zgejsv.f"
	    } else {
#line 1772 "zgejsv.f"
		i__1 = nr - 1;
#line 1772 "zgejsv.f"
		i__2 = nr - 1;
#line 1772 "zgejsv.f"
		zlaset_("U", &i__1, &i__2, &c_b1, &c_b1, &u[(u_dim1 << 1) + 1]
			, ldu, (ftnlen)1);
#line 1773 "zgejsv.f"
	    }
#line 1775 "zgejsv.f"
	    i__1 = *lwork - (*n << 1) - *n * nr;
#line 1775 "zgejsv.f"
	    zgesvj_("L", "U", "V", &nr, &nr, &u[u_offset], ldu, &sva[1], n, &
		    v[v_offset], ldv, &cwork[(*n << 1) + *n * nr + 1], &i__1, 
		    &rwork[1], lrwork, info, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1778 "zgejsv.f"
	    scalem = rwork[1];
#line 1779 "zgejsv.f"
	    numrank = i_dnnt(&rwork[2]);
#line 1781 "zgejsv.f"
	    if (nr < *n) {
#line 1782 "zgejsv.f"
		i__1 = *n - nr;
#line 1782 "zgejsv.f"
		zlaset_("A", &i__1, &nr, &c_b120, &c_b120, &v[nr + 1 + v_dim1]
			, ldv, (ftnlen)1);
#line 1783 "zgejsv.f"
		i__1 = *n - nr;
#line 1783 "zgejsv.f"
		zlaset_("A", &nr, &i__1, &c_b120, &c_b120, &v[(nr + 1) * 
			v_dim1 + 1], ldv, (ftnlen)1);
#line 1784 "zgejsv.f"
		i__1 = *n - nr;
#line 1784 "zgejsv.f"
		i__2 = *n - nr;
#line 1784 "zgejsv.f"
		zlaset_("A", &i__1, &i__2, &c_b120, &c_b80, &v[nr + 1 + (nr + 
			1) * v_dim1], ldv, (ftnlen)1);
#line 1785 "zgejsv.f"
	    }
#line 1787 "zgejsv.f"
	    i__1 = *lwork - (*n << 1) - *n * nr - nr;
#line 1787 "zgejsv.f"
	    zunmqr_("L", "N", n, n, &nr, &cwork[(*n << 1) + 1], n, &cwork[*n 
		    + 1], &v[v_offset], ldv, &cwork[(*n << 1) + *n * nr + nr 
		    + 1], &i__1, &ierr, (ftnlen)1, (ftnlen)1);

/*           Permute the rows of V using the (column) permutation from the */
/*           first QRF. Also, scale the columns to make them unit in */
/*           Euclidean norm. This applies to all cases. */

#line 1794 "zgejsv.f"
	    temp1 = sqrt((doublereal) (*n)) * epsln;
#line 1795 "zgejsv.f"
	    i__1 = *n;
#line 1795 "zgejsv.f"
	    for (q = 1; q <= i__1; ++q) {
#line 1796 "zgejsv.f"
		i__2 = *n;
#line 1796 "zgejsv.f"
		for (p = 1; p <= i__2; ++p) {
#line 1797 "zgejsv.f"
		    i__3 = (*n << 1) + *n * nr + nr + iwork[p];
#line 1797 "zgejsv.f"
		    i__4 = p + q * v_dim1;
#line 1797 "zgejsv.f"
		    cwork[i__3].r = v[i__4].r, cwork[i__3].i = v[i__4].i;
#line 1798 "zgejsv.f"
/* L8972: */
#line 1798 "zgejsv.f"
		}
#line 1799 "zgejsv.f"
		i__2 = *n;
#line 1799 "zgejsv.f"
		for (p = 1; p <= i__2; ++p) {
#line 1800 "zgejsv.f"
		    i__3 = p + q * v_dim1;
#line 1800 "zgejsv.f"
		    i__4 = (*n << 1) + *n * nr + nr + p;
#line 1800 "zgejsv.f"
		    v[i__3].r = cwork[i__4].r, v[i__3].i = cwork[i__4].i;
#line 1801 "zgejsv.f"
/* L8973: */
#line 1801 "zgejsv.f"
		}
#line 1802 "zgejsv.f"
		xsc = 1. / dznrm2_(n, &v[q * v_dim1 + 1], &c__1);
#line 1803 "zgejsv.f"
		if (xsc < 1. - temp1 || xsc > temp1 + 1.) {
#line 1803 "zgejsv.f"
		    zdscal_(n, &xsc, &v[q * v_dim1 + 1], &c__1);
#line 1803 "zgejsv.f"
		}
#line 1805 "zgejsv.f"
/* L7972: */
#line 1805 "zgejsv.f"
	    }

/*           At this moment, V contains the right singular vectors of A. */
/*           Next, assemble the left singular vector matrix U (M x N). */

#line 1810 "zgejsv.f"
	    if (nr < *m) {
#line 1811 "zgejsv.f"
		i__1 = *m - nr;
#line 1811 "zgejsv.f"
		zlaset_("A", &i__1, &nr, &c_b1, &c_b1, &u[nr + 1 + u_dim1], 
			ldu, (ftnlen)1);
#line 1812 "zgejsv.f"
		if (nr < n1) {
#line 1813 "zgejsv.f"
		    i__1 = n1 - nr;
#line 1813 "zgejsv.f"
		    zlaset_("A", &nr, &i__1, &c_b1, &c_b1, &u[(nr + 1) * 
			    u_dim1 + 1], ldu, (ftnlen)1);
#line 1814 "zgejsv.f"
		    i__1 = *m - nr;
#line 1814 "zgejsv.f"
		    i__2 = n1 - nr;
#line 1814 "zgejsv.f"
		    zlaset_("A", &i__1, &i__2, &c_b1, &c_b2, &u[nr + 1 + (nr 
			    + 1) * u_dim1], ldu, (ftnlen)1);
#line 1815 "zgejsv.f"
		}
#line 1816 "zgejsv.f"
	    }

#line 1818 "zgejsv.f"
	    i__1 = *lwork - *n;
#line 1818 "zgejsv.f"
	    zunmqr_("Left", "No Tr", m, &n1, n, &a[a_offset], lda, &cwork[1], 
		    &u[u_offset], ldu, &cwork[*n + 1], &i__1, &ierr, (ftnlen)
		    4, (ftnlen)5);

#line 1821 "zgejsv.f"
	    if (rowpiv) {
#line 1821 "zgejsv.f"
		i__1 = *m - 1;
#line 1821 "zgejsv.f"
		zlaswp_(&n1, &u[u_offset], ldu, &c__1, &i__1, &iwork[(*n << 1)
			 + 1], &c_n1);
#line 1821 "zgejsv.f"
	    }


#line 1825 "zgejsv.f"
	}
#line 1826 "zgejsv.f"
	if (transp) {
/*           .. swap U and V because the procedure worked on A^* */
#line 1828 "zgejsv.f"
	    i__1 = *n;
#line 1828 "zgejsv.f"
	    for (p = 1; p <= i__1; ++p) {
#line 1829 "zgejsv.f"
		zswap_(n, &u[p * u_dim1 + 1], &c__1, &v[p * v_dim1 + 1], &
			c__1);
#line 1830 "zgejsv.f"
/* L6974: */
#line 1830 "zgejsv.f"
	    }
#line 1831 "zgejsv.f"
	}

#line 1833 "zgejsv.f"
    }
/*     end of the full SVD */

/*     Undo scaling, if necessary (and possible) */

#line 1838 "zgejsv.f"
    if (uscal2 <= big / sva[1] * uscal1) {
#line 1839 "zgejsv.f"
	zlascl_("G", &c__0, &c__0, &uscal1, &uscal2, &nr, &c__1, &sva[1], n, &
		ierr, (ftnlen)1);
#line 1840 "zgejsv.f"
	uscal1 = 1.;
#line 1841 "zgejsv.f"
	uscal2 = 1.;
#line 1842 "zgejsv.f"
    }

#line 1844 "zgejsv.f"
    if (nr < *n) {
#line 1845 "zgejsv.f"
	i__1 = *n;
#line 1845 "zgejsv.f"
	for (p = nr + 1; p <= i__1; ++p) {
#line 1846 "zgejsv.f"
	    sva[p] = 0.;
#line 1847 "zgejsv.f"
/* L3004: */
#line 1847 "zgejsv.f"
	}
#line 1848 "zgejsv.f"
    }

#line 1850 "zgejsv.f"
    rwork[1] = uscal2 * scalem;
#line 1851 "zgejsv.f"
    rwork[2] = uscal1;
#line 1852 "zgejsv.f"
    if (errest) {
#line 1852 "zgejsv.f"
	rwork[3] = sconda;
#line 1852 "zgejsv.f"
    }
#line 1853 "zgejsv.f"
    if (lsvec && rsvec) {
#line 1854 "zgejsv.f"
	rwork[4] = condr1;
#line 1855 "zgejsv.f"
	rwork[5] = condr2;
#line 1856 "zgejsv.f"
    }
#line 1857 "zgejsv.f"
    if (l2tran) {
#line 1858 "zgejsv.f"
	rwork[6] = entra;
#line 1859 "zgejsv.f"
	rwork[7] = entrat;
#line 1860 "zgejsv.f"
    }

#line 1862 "zgejsv.f"
    iwork[1] = nr;
#line 1863 "zgejsv.f"
    iwork[2] = numrank;
#line 1864 "zgejsv.f"
    iwork[3] = warning;

#line 1866 "zgejsv.f"
    return 0;
/*     .. */
/*     .. END OF ZGEJSV */
/*     .. */
} /* zgejsv_ */

