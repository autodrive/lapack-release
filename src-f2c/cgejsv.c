#line 1 "cgejsv.f"
/* cgejsv.f -- translated by f2c (version 20100827).
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

#line 1 "cgejsv.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};
static integer c__1 = 1;
static integer c__0 = 0;
static doublereal c_b80 = 1.;
static doublereal c_b120 = 0.;
static integer c_n1 = -1;

/* > \brief \b CGEJSV */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CGEJSV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgejsv.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgejsv.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgejsv.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*     SUBROUTINE CGEJSV( JOBA, JOBU, JOBV, JOBR, JOBT, JOBP, */
/*                         M, N, A, LDA, SVA, U, LDU, V, LDV, */
/*                         CWORK, LWORK, RWORK, LRWORK, IWORK, INFO ) */

/*     .. Scalar Arguments .. */
/*     IMPLICIT    NONE */
/*     INTEGER     INFO, LDA, LDU, LDV, LWORK, M, N */
/*     .. */
/*     .. Array Arguments .. */
/*     COMPLEX     A( LDA, * ),  U( LDU, * ), V( LDV, * ), CWORK( LWORK ) */
/*     REAL        SVA( N ), RWORK( LRWORK ) */
/*     INTEGER     IWORK( * ) */
/*     CHARACTER*1 JOBA, JOBP, JOBR, JOBT, JOBU, JOBV */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CGEJSV computes the singular value decomposition (SVD) of a real M-by-N */
/* > matrix [A], where M >= N. The SVD of [A] is written as */
/* > */
/* >              [A] = [U] * [SIGMA] * [V]^*, */
/* > */
/* > where [SIGMA] is an N-by-N (M-by-N) matrix which is zero except for its N */
/* > diagonal elements, [U] is an M-by-N (or M-by-M) orthonormal matrix, and */
/* > [V] is an N-by-N orthogonal matrix. The diagonal elements of [SIGMA] are */
/* > the singular values of [A]. The columns of [U] and [V] are the left and */
/* > the right singular vectors of [A], respectively. The matrices [U] and [V] */
/* > are computed and stored in the arrays U and V, respectively. The diagonal */
/* > of [SIGMA] is computed and stored in the array SVA. */
/* > */
/* >  Arguments: */
/* >  ========== */
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
/* >         value of c*A is around SQRT(BIG), BIG=SLAMCH('O'), then JOBR issues */
/* >         the licence to kill columns of A whose norm in c*A is less than */
/* >         SQRT(SFMIN) (for JOBR.EQ.'R'), or less than SMALL=SFMIN/EPSLN, */
/* >         where SFMIN=SLAMCH('S'), EPSLN=SLAMCH('E'). */
/* >       = 'N': Do not kill small columns of c*A. This option assumes that */
/* >              BLAS and QR factorizations and triangular solvers are */
/* >              implemented to work in that range. If the condition of A */
/* >              is greater than BIG, use CGESVJ. */
/* >       = 'R': RESTRICTED range for sigma(c*A) is [SQRT(SFMIN), SQRT(BIG)] */
/* >              (roughly, as described above). This option is recommended. */
/* >                                             =========================== */
/* >         For computing the singular values in the FULL range [SFMIN,BIG] */
/* >         use CGESVJ. */
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
/* >          A is COMPLEX array, dimension (LDA,N) */
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
/* >          SVA is REAL array, dimension (N) */
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
/* >          U is COMPLEX array, dimension ( LDU, N ) */
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
/* >          V is COMPLEX array, dimension ( LDV, N ) */
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
/* >          CWORK is COMPLEX array, dimension at least LWORK. */
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
/* >               block size for CGEQP3 and CGEQRF. */
/* >               In general, optimal LWORK is computed as */
/* >               LWORK >= max(N+LWORK(CGEQP3),N+LWORK(CGEQRF)). */
/* >            1.2. .. an estimate of the scaled condition number of A is */
/* >               required (JOBA='E', or 'G'). In this case, LWORK the minimal */
/* >               requirement is LWORK >= N*N + 3*N. */
/* >               ->> For optimal performance (blocked code) the optimal value */
/* >               is LWORK >= max(N+(N+1)*NB, N*N+3*N). */
/* >               In general, the optimal length LWORK is computed as */
/* >               LWORK >= max(N+LWORK(CGEQP3),N+LWORK(CGEQRF), */
/* >                                                     N+N*N+LWORK(CPOCON)). */
/* > */
/* >          2. If SIGMA and the right singular vectors are needed (JOBV.EQ.'V'), */
/* >             (JOBU.EQ.'N') */
/* >            -> the minimal requirement is LWORK >= 3*N. */
/* >            -> For optimal performance, LWORK >= max(N+(N+1)*NB, 3*N,2*N+N*NB), */
/* >               where NB is the optimal block size for CGEQP3, CGEQRF, CGELQ, */
/* >               CUNMLQ. In general, the optimal length LWORK is computed as */
/* >               LWORK >= max(N+LWORK(CGEQP3), N+LWORK(CPOCON), N+LWORK(CGESVJ), */
/* >                       N+LWORK(CGELQF), 2*N+LWORK(CGEQRF), N+LWORK(CUNMLQ)). */
/* > */
/* >          3. If SIGMA and the left singular vectors are needed */
/* >            -> the minimal requirement is LWORK >= 3*N. */
/* >            -> For optimal performance: */
/* >               if JOBU.EQ.'U' :: LWORK >= max(3*N, N+(N+1)*NB, 2*N+N*NB), */
/* >               where NB is the optimal block size for CGEQP3, CGEQRF, CUNMQR. */
/* >               In general, the optimal length LWORK is computed as */
/* >               LWORK >= max(N+LWORK(CGEQP3),N+LWORK(CPOCON), */
/* >                        2*N+LWORK(CGEQRF), N+LWORK(CUNMQR)). */
/* > */
/* >          4. If the full SVD is needed: (JOBU.EQ.'U' or JOBU.EQ.'F') and */
/* >            4.1. if JOBV.EQ.'V' */
/* >               the minimal requirement is LWORK >= 5*N+2*N*N. */
/* >            4.2. if JOBV.EQ.'J' the minimal requirement is */
/* >               LWORK >= 4*N+N*N. */
/* >            In both cases, the allocated CWORK can accomodate blocked runs */
/* >            of CGEQP3, CGEQRF, CGELQF, SUNMQR, CUNMLQ. */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is REAL array, dimension at least LRWORK. */
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
/* >           > 0 :  CGEJSV  did not converge in the maximal allowed number */
/* >                  of sweeps. The computed values may be inaccurate. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2015 */

/* > \ingroup complexGEsing */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  CGEJSV implements a preconditioned Jacobi SVD algorithm. It uses CGEQP3, */
/* >  CGEQRF, and CGELQF as preprocessors and preconditioners. Optionally, an */
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
/* >  & LAPACK routines called by CGEJSV are implemented to work in that range. */
/* >  If that is not the case, then the restriction for safe computation with */
/* >  the singular values in the range of normalized IEEE numbers is that the */
/* >  spectral condition number kappa(A)=sigma_max(A)/sigma_min(A) does not */
/* >  overflow. This code (CGEJSV) is best used in this restricted range, */
/* >  meaning that singular values of magnitude below ||A||_2 / SLAMCH('O') are */
/* >  returned as zeros. See JOBR for details on this. */
/* >     Further, this implementation is somewhat slower than the one described */
/* >  in [1,2] due to replacement of some non-LAPACK components, and because */
/* >  the choice of some tuning parameters in the iterative part (CGESVJ) is */
/* >  left to the implementer on a particular machine. */
/* >     The rank revealing QR factorization (in this code: CGEQP3) should be */
/* >  implemented as in [3]. We have a new version of CGEQP3 under development */
/* >  that is more robust than the current one in LAPACK, with a cleaner cut in */
/* >  rank defficient cases. It will be available in the SIGMA library [4]. */
/* >  If M is much larger than N, it is obvious that the inital QRF with */
/* >  column pivoting can be preprocessed by the QRF without pivoting. That */
/* >  well known trick is not used in CGEJSV because in some cases heavy row */
/* >  weighting can be treated with complete pivoting. The overhead in cases */
/* >  M much larger than N is then only due to pivoting, but the benefits in */
/* >  terms of accuracy have prevailed. The implementer/user can incorporate */
/* >  this extra QRF step easily. The implementer can also improve data movement */
/* >  (matrix transpose, matrix copy, matrix transposed copy) - this */
/* >  implementation of CGEJSV uses only the simplest, naive data movement. */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >  Zlatko Drmac (Zagreb, Croatia) and Kresimir Veselic (Hagen, Germany) */

/* > \par References: */
/*  ================ */
/* > */
/* > \verbatim */
/* > */
/* > [1] Z. Drmac and K. Veselic: New fast and accurate Jacobi SVD algorithm I. */
/* >     SIAM J. Matrix Anal. Appl. Vol. 35, No. 2 (2008), pp. 1322-1342. */
/* >     LAPACK Working note 169. */
/* > [2] Z. Drmac and K. Veselic: New fast and accurate Jacobi SVD algorithm II. */
/* >     SIAM J. Matrix Anal. Appl. Vol. 35, No. 2 (2008), pp. 1343-1362. */
/* >     LAPACK Working note 170. */
/* > [3] Z. Drmac and Z. Bujanovic: On the failure of rank-revealing QR */
/* >     factorization software - a case study. */
/* >     ACM Trans. math. Softw. Vol. 35, No 2 (2008), pp. 1-28. */
/* >     LAPACK Working note 176. */
/* > [4] Z. Drmac: SIGMA - mathematical software library for accurate SVD, PSV, */
/* >     QSVD, (H,K)-SVD computations. */
/* >     Department of Mathematics, University of Zagreb, 2008. */
/* > \endverbatim */

/* >  \par Bugs, examples and comments: */
/*   ================================= */
/* > */
/* >  Please report all bugs and send interesting examples and/or comments to */
/* >  drmac@math.hr. Thank you. */
/* > */
/*  ===================================================================== */
/* Subroutine */ int cgejsv_(char *joba, char *jobu, char *jobv, char *jobr, 
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
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static doublecomplex ctemp;
    static doublereal entra, small, sfmin;
    static logical lsvec;
    extern /* Subroutine */ int ccopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), cswap_(integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *);
    static doublereal epsln;
    static logical rsvec;
    extern /* Subroutine */ int ctrsm_(char *, char *, char *, char *, 
	    integer *, integer *, doublecomplex *, doublecomplex *, integer *,
	     doublecomplex *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static logical l2aber;
    extern /* Subroutine */ int cgeqp3_(integer *, integer *, doublecomplex *,
	     integer *, integer *, doublecomplex *, doublecomplex *, integer *
	    , doublereal *, integer *);
    static doublereal condr1, condr2, uscal1, uscal2;
    static logical l2kill, l2rank, l2tran;
    extern doublereal scnrm2_(integer *, doublecomplex *, integer *);
    static logical l2pert;
    extern /* Subroutine */ int clacgv_(integer *, doublecomplex *, integer *)
	    , cgelqf_(integer *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, integer *), clascl_()
	    ;
    static doublereal scalem, sconda;
    static logical goscal;
    static doublereal aatmin;
    extern doublereal slamch_(char *, ftnlen);
    static doublereal aatmax;
    extern /* Subroutine */ int cgeqrf_(integer *, integer *, doublecomplex *,
	     integer *, doublecomplex *, doublecomplex *, integer *, integer *
	    ), clacpy_(char *, integer *, integer *, doublecomplex *, integer 
	    *, doublecomplex *, integer *, ftnlen), claset_();
    static logical noscal;
    extern /* Subroutine */ int csscal_(integer *, doublereal *, 
	    doublecomplex *, integer *);
    extern integer isamax_();
    extern /* Subroutine */ int cpocon_(char *, integer *, doublecomplex *, 
	    integer *, doublereal *, doublereal *, doublecomplex *, 
	    doublereal *, integer *, ftnlen), xerbla_(char *, integer *, 
	    ftnlen), cgesvj_(char *, char *, char *, integer *, integer *, 
	    doublecomplex *, integer *, doublereal *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublereal *, integer *, integer *, ftnlen, ftnlen, ftnlen), 
	    classq_(), claswp_(integer *, doublecomplex *, integer *, integer 
	    *, integer *, integer *, integer *);
    static doublereal entrat;
    static logical almort;
    extern /* Subroutine */ int cungqr_(integer *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, integer *), cunmlq_(char *, char *, integer *, integer 
	    *, integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     ftnlen, ftnlen);
    static doublereal maxprj;
    static logical errest;
    extern /* Subroutine */ int cunmqr_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     ftnlen, ftnlen);
    static logical transp, rowpiv;
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

#line 573 "cgejsv.f"
    /* Parameter adjustments */
#line 573 "cgejsv.f"
    --sva;
#line 573 "cgejsv.f"
    a_dim1 = *lda;
#line 573 "cgejsv.f"
    a_offset = 1 + a_dim1;
#line 573 "cgejsv.f"
    a -= a_offset;
#line 573 "cgejsv.f"
    u_dim1 = *ldu;
#line 573 "cgejsv.f"
    u_offset = 1 + u_dim1;
#line 573 "cgejsv.f"
    u -= u_offset;
#line 573 "cgejsv.f"
    v_dim1 = *ldv;
#line 573 "cgejsv.f"
    v_offset = 1 + v_dim1;
#line 573 "cgejsv.f"
    v -= v_offset;
#line 573 "cgejsv.f"
    --cwork;
#line 573 "cgejsv.f"
    --rwork;
#line 573 "cgejsv.f"
    --iwork;
#line 573 "cgejsv.f"

#line 573 "cgejsv.f"
    /* Function Body */
#line 573 "cgejsv.f"
    lsvec = lsame_(jobu, "U", (ftnlen)1, (ftnlen)1) || lsame_(jobu, "F", (
	    ftnlen)1, (ftnlen)1);
#line 574 "cgejsv.f"
    jracc = lsame_(jobv, "J", (ftnlen)1, (ftnlen)1);
#line 575 "cgejsv.f"
    rsvec = lsame_(jobv, "V", (ftnlen)1, (ftnlen)1) || jracc;
#line 576 "cgejsv.f"
    rowpiv = lsame_(joba, "F", (ftnlen)1, (ftnlen)1) || lsame_(joba, "G", (
	    ftnlen)1, (ftnlen)1);
#line 577 "cgejsv.f"
    l2rank = lsame_(joba, "R", (ftnlen)1, (ftnlen)1);
#line 578 "cgejsv.f"
    l2aber = lsame_(joba, "A", (ftnlen)1, (ftnlen)1);
#line 579 "cgejsv.f"
    errest = lsame_(joba, "E", (ftnlen)1, (ftnlen)1) || lsame_(joba, "G", (
	    ftnlen)1, (ftnlen)1);
#line 580 "cgejsv.f"
    l2tran = lsame_(jobt, "T", (ftnlen)1, (ftnlen)1);
#line 581 "cgejsv.f"
    l2kill = lsame_(jobr, "R", (ftnlen)1, (ftnlen)1);
#line 582 "cgejsv.f"
    defr = lsame_(jobr, "N", (ftnlen)1, (ftnlen)1);
#line 583 "cgejsv.f"
    l2pert = lsame_(jobp, "P", (ftnlen)1, (ftnlen)1);

#line 585 "cgejsv.f"
    if (! (rowpiv || l2rank || l2aber || errest || lsame_(joba, "C", (ftnlen)
	    1, (ftnlen)1))) {
#line 587 "cgejsv.f"
	*info = -1;
#line 588 "cgejsv.f"
    } else if (! (lsvec || lsame_(jobu, "N", (ftnlen)1, (ftnlen)1) || lsame_(
	    jobu, "W", (ftnlen)1, (ftnlen)1))) {
#line 590 "cgejsv.f"
	*info = -2;
#line 591 "cgejsv.f"
    } else if (! (rsvec || lsame_(jobv, "N", (ftnlen)1, (ftnlen)1) || lsame_(
	    jobv, "W", (ftnlen)1, (ftnlen)1)) || jracc && ! lsvec) {
#line 593 "cgejsv.f"
	*info = -3;
#line 594 "cgejsv.f"
    } else if (! (l2kill || defr)) {
#line 595 "cgejsv.f"
	*info = -4;
#line 596 "cgejsv.f"
    } else if (! (l2tran || lsame_(jobt, "N", (ftnlen)1, (ftnlen)1))) {
#line 597 "cgejsv.f"
	*info = -5;
#line 598 "cgejsv.f"
    } else if (! (l2pert || lsame_(jobp, "N", (ftnlen)1, (ftnlen)1))) {
#line 599 "cgejsv.f"
	*info = -6;
#line 600 "cgejsv.f"
    } else if (*m < 0) {
#line 601 "cgejsv.f"
	*info = -7;
#line 602 "cgejsv.f"
    } else if (*n < 0 || *n > *m) {
#line 603 "cgejsv.f"
	*info = -8;
#line 604 "cgejsv.f"
    } else if (*lda < *m) {
#line 605 "cgejsv.f"
	*info = -10;
#line 606 "cgejsv.f"
    } else if (lsvec && *ldu < *m) {
#line 607 "cgejsv.f"
	*info = -13;
#line 608 "cgejsv.f"
    } else if (rsvec && *ldv < *n) {
#line 609 "cgejsv.f"
	*info = -15;
#line 610 "cgejsv.f"
    } else if (! (lsvec || rsvec || errest) && *lwork < (*n << 1) + 1 || ! (
	    lsvec || rsvec) && errest && *lwork < *n * *n + *n * 3 || lsvec &&
	     ! rsvec && *lwork < *n * 3 || rsvec && ! lsvec && *lwork < *n * 
	    3 || lsvec && rsvec && ! jracc && *lwork < *n * 5 + (*n << 1) * *
	    n || lsvec && rsvec && jracc && *lwork < (*n << 2) + *n * *n) {
#line 623 "cgejsv.f"
	*info = -17;
#line 624 "cgejsv.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 624 "cgejsv.f"
	i__1 = *n + (*m << 1);
#line 624 "cgejsv.f"
	if (*lrwork < max(i__1,7)) {
#line 625 "cgejsv.f"
	    *info = -19;
#line 626 "cgejsv.f"
	} else {
/*        #:) */
#line 628 "cgejsv.f"
	    *info = 0;
#line 629 "cgejsv.f"
	}
#line 629 "cgejsv.f"
    }

#line 631 "cgejsv.f"
    if (*info != 0) {
/*       #:( */
#line 633 "cgejsv.f"
	i__1 = -(*info);
#line 633 "cgejsv.f"
	xerbla_("CGEJSV", &i__1, (ftnlen)6);
#line 634 "cgejsv.f"
	return 0;
#line 635 "cgejsv.f"
    }

/*     Quick return for void matrix (Y3K safe) */
/* #:) */
#line 639 "cgejsv.f"
    if (*m == 0 || *n == 0) {
#line 639 "cgejsv.f"
	return 0;
#line 639 "cgejsv.f"
    }

/*     Determine whether the matrix U should be M x N or M x M */

#line 643 "cgejsv.f"
    if (lsvec) {
#line 644 "cgejsv.f"
	n1 = *n;
#line 645 "cgejsv.f"
	if (lsame_(jobu, "F", (ftnlen)1, (ftnlen)1)) {
#line 645 "cgejsv.f"
	    n1 = *m;
#line 645 "cgejsv.f"
	}
#line 646 "cgejsv.f"
    }

/*     Set numerical parameters */

/* !    NOTE: Make sure SLAMCH() does not fail on the target architecture. */

#line 652 "cgejsv.f"
    epsln = slamch_("Epsilon", (ftnlen)7);
#line 653 "cgejsv.f"
    sfmin = slamch_("SafeMinimum", (ftnlen)11);
#line 654 "cgejsv.f"
    small = sfmin / epsln;
#line 655 "cgejsv.f"
    big = slamch_("O", (ftnlen)1);
/*     BIG   = ONE / SFMIN */

/*     Initialize SVA(1:N) = diag( ||A e_i||_2 )_1^N */

/* (!)  If necessary, scale SVA() to protect the largest norm from */
/*     overflow. It is possible that this scaling pushes the smallest */
/*     column norm left from the underflow threshold (extreme case). */

#line 664 "cgejsv.f"
    scalem = 1. / sqrt((doublereal) (*m) * (doublereal) (*n));
#line 665 "cgejsv.f"
    noscal = TRUE_;
#line 666 "cgejsv.f"
    goscal = TRUE_;
#line 667 "cgejsv.f"
    i__1 = *n;
#line 667 "cgejsv.f"
    for (p = 1; p <= i__1; ++p) {
#line 668 "cgejsv.f"
	aapp = 0.;
#line 669 "cgejsv.f"
	aaqq = 1.;
#line 670 "cgejsv.f"
	classq_(m, &a[p * a_dim1 + 1], &c__1, &aapp, &aaqq);
#line 671 "cgejsv.f"
	if (aapp > big) {
#line 672 "cgejsv.f"
	    *info = -9;
#line 673 "cgejsv.f"
	    i__2 = -(*info);
#line 673 "cgejsv.f"
	    xerbla_("CGEJSV", &i__2, (ftnlen)6);
#line 674 "cgejsv.f"
	    return 0;
#line 675 "cgejsv.f"
	}
#line 676 "cgejsv.f"
	aaqq = sqrt(aaqq);
#line 677 "cgejsv.f"
	if (aapp < big / aaqq && noscal) {
#line 678 "cgejsv.f"
	    sva[p] = aapp * aaqq;
#line 679 "cgejsv.f"
	} else {
#line 680 "cgejsv.f"
	    noscal = FALSE_;
#line 681 "cgejsv.f"
	    sva[p] = aapp * (aaqq * scalem);
#line 682 "cgejsv.f"
	    if (goscal) {
#line 683 "cgejsv.f"
		goscal = FALSE_;
#line 684 "cgejsv.f"
		i__2 = p - 1;
#line 684 "cgejsv.f"
		sscal_(&i__2, &scalem, &sva[1], &c__1);
#line 685 "cgejsv.f"
	    }
#line 686 "cgejsv.f"
	}
#line 687 "cgejsv.f"
/* L1874: */
#line 687 "cgejsv.f"
    }

#line 689 "cgejsv.f"
    if (noscal) {
#line 689 "cgejsv.f"
	scalem = 1.;
#line 689 "cgejsv.f"
    }

#line 691 "cgejsv.f"
    aapp = 0.;
#line 692 "cgejsv.f"
    aaqq = big;
#line 693 "cgejsv.f"
    i__1 = *n;
#line 693 "cgejsv.f"
    for (p = 1; p <= i__1; ++p) {
/* Computing MAX */
#line 694 "cgejsv.f"
	d__1 = aapp, d__2 = sva[p];
#line 694 "cgejsv.f"
	aapp = max(d__1,d__2);
#line 695 "cgejsv.f"
	if (sva[p] != 0.) {
/* Computing MIN */
#line 695 "cgejsv.f"
	    d__1 = aaqq, d__2 = sva[p];
#line 695 "cgejsv.f"
	    aaqq = min(d__1,d__2);
#line 695 "cgejsv.f"
	}
#line 696 "cgejsv.f"
/* L4781: */
#line 696 "cgejsv.f"
    }

/*     Quick return for zero M x N matrix */
/* #:) */
#line 700 "cgejsv.f"
    if (aapp == 0.) {
#line 701 "cgejsv.f"
	if (lsvec) {
#line 701 "cgejsv.f"
	    claset_("G", m, &n1, &c_b1, &c_b2, &u[u_offset], ldu, (ftnlen)1);
#line 701 "cgejsv.f"
	}
#line 702 "cgejsv.f"
	if (rsvec) {
#line 702 "cgejsv.f"
	    claset_("G", n, n, &c_b1, &c_b2, &v[v_offset], ldv, (ftnlen)1);
#line 702 "cgejsv.f"
	}
#line 703 "cgejsv.f"
	rwork[1] = 1.;
#line 704 "cgejsv.f"
	rwork[2] = 1.;
#line 705 "cgejsv.f"
	if (errest) {
#line 705 "cgejsv.f"
	    rwork[3] = 1.;
#line 705 "cgejsv.f"
	}
#line 706 "cgejsv.f"
	if (lsvec && rsvec) {
#line 707 "cgejsv.f"
	    rwork[4] = 1.;
#line 708 "cgejsv.f"
	    rwork[5] = 1.;
#line 709 "cgejsv.f"
	}
#line 710 "cgejsv.f"
	if (l2tran) {
#line 711 "cgejsv.f"
	    rwork[6] = 0.;
#line 712 "cgejsv.f"
	    rwork[7] = 0.;
#line 713 "cgejsv.f"
	}
#line 714 "cgejsv.f"
	iwork[1] = 0;
#line 715 "cgejsv.f"
	iwork[2] = 0;
#line 716 "cgejsv.f"
	iwork[3] = 0;
#line 717 "cgejsv.f"
	return 0;
#line 718 "cgejsv.f"
    }

/*     Issue warning if denormalized column norms detected. Override the */
/*     high relative accuracy request. Issue licence to kill columns */
/*     (set them to zero) whose norm is less than sigma_max / BIG (roughly). */
/* #:( */
#line 724 "cgejsv.f"
    warning = 0;
#line 725 "cgejsv.f"
    if (aaqq <= sfmin) {
#line 726 "cgejsv.f"
	l2rank = TRUE_;
#line 727 "cgejsv.f"
	l2kill = TRUE_;
#line 728 "cgejsv.f"
	warning = 1;
#line 729 "cgejsv.f"
    }

/*     Quick return for one-column matrix */
/* #:) */
#line 733 "cgejsv.f"
    if (*n == 1) {

#line 735 "cgejsv.f"
	if (lsvec) {
#line 736 "cgejsv.f"
	    clascl_("G", &c__0, &c__0, &sva[1], &scalem, m, &c__1, &a[a_dim1 
		    + 1], lda, &ierr, (ftnlen)1);
#line 737 "cgejsv.f"
	    clacpy_("A", m, &c__1, &a[a_offset], lda, &u[u_offset], ldu, (
		    ftnlen)1);
/*           computing all M left singular vectors of the M x 1 matrix */
#line 739 "cgejsv.f"
	    if (n1 != *n) {
#line 740 "cgejsv.f"
		i__1 = *lwork - *n;
#line 740 "cgejsv.f"
		cgeqrf_(m, n, &u[u_offset], ldu, &cwork[1], &cwork[*n + 1], &
			i__1, &ierr);
#line 741 "cgejsv.f"
		i__1 = *lwork - *n;
#line 741 "cgejsv.f"
		cungqr_(m, &n1, &c__1, &u[u_offset], ldu, &cwork[1], &cwork[*
			n + 1], &i__1, &ierr);
#line 742 "cgejsv.f"
		ccopy_(m, &a[a_dim1 + 1], &c__1, &u[u_dim1 + 1], &c__1);
#line 743 "cgejsv.f"
	    }
#line 744 "cgejsv.f"
	}
#line 745 "cgejsv.f"
	if (rsvec) {
#line 746 "cgejsv.f"
	    i__1 = v_dim1 + 1;
#line 746 "cgejsv.f"
	    v[i__1].r = 1., v[i__1].i = 0.;
#line 747 "cgejsv.f"
	}
#line 748 "cgejsv.f"
	if (sva[1] < big * scalem) {
#line 749 "cgejsv.f"
	    sva[1] /= scalem;
#line 750 "cgejsv.f"
	    scalem = 1.;
#line 751 "cgejsv.f"
	}
#line 752 "cgejsv.f"
	rwork[1] = 1. / scalem;
#line 753 "cgejsv.f"
	rwork[2] = 1.;
#line 754 "cgejsv.f"
	if (sva[1] != 0.) {
#line 755 "cgejsv.f"
	    iwork[1] = 1;
#line 756 "cgejsv.f"
	    if (sva[1] / scalem >= sfmin) {
#line 757 "cgejsv.f"
		iwork[2] = 1;
#line 758 "cgejsv.f"
	    } else {
#line 759 "cgejsv.f"
		iwork[2] = 0;
#line 760 "cgejsv.f"
	    }
#line 761 "cgejsv.f"
	} else {
#line 762 "cgejsv.f"
	    iwork[1] = 0;
#line 763 "cgejsv.f"
	    iwork[2] = 0;
#line 764 "cgejsv.f"
	}
#line 765 "cgejsv.f"
	iwork[3] = 0;
#line 766 "cgejsv.f"
	if (errest) {
#line 766 "cgejsv.f"
	    rwork[3] = 1.;
#line 766 "cgejsv.f"
	}
#line 767 "cgejsv.f"
	if (lsvec && rsvec) {
#line 768 "cgejsv.f"
	    rwork[4] = 1.;
#line 769 "cgejsv.f"
	    rwork[5] = 1.;
#line 770 "cgejsv.f"
	}
#line 771 "cgejsv.f"
	if (l2tran) {
#line 772 "cgejsv.f"
	    rwork[6] = 0.;
#line 773 "cgejsv.f"
	    rwork[7] = 0.;
#line 774 "cgejsv.f"
	}
#line 775 "cgejsv.f"
	return 0;

#line 777 "cgejsv.f"
    }

#line 779 "cgejsv.f"
    transp = FALSE_;
#line 780 "cgejsv.f"
    l2tran = l2tran && *m == *n;

#line 782 "cgejsv.f"
    aatmax = -1.;
#line 783 "cgejsv.f"
    aatmin = big;
#line 784 "cgejsv.f"
    if (rowpiv || l2tran) {

/*     Compute the row norms, needed to determine row pivoting sequence */
/*     (in the case of heavily row weighted A, row pivoting is strongly */
/*     advised) and to collect information needed to compare the */
/*     structures of A * A^* and A^* * A (in the case L2TRAN.EQ..TRUE.). */

#line 791 "cgejsv.f"
	if (l2tran) {
#line 792 "cgejsv.f"
	    i__1 = *m;
#line 792 "cgejsv.f"
	    for (p = 1; p <= i__1; ++p) {
#line 793 "cgejsv.f"
		xsc = 0.;
#line 794 "cgejsv.f"
		temp1 = 1.;
#line 795 "cgejsv.f"
		classq_(n, &a[p + a_dim1], lda, &xsc, &temp1);
/*              CLASSQ gets both the ell_2 and the ell_infinity norm */
/*              in one pass through the vector */
#line 798 "cgejsv.f"
		rwork[*m + *n + p] = xsc * scalem;
#line 799 "cgejsv.f"
		rwork[*n + p] = xsc * (scalem * sqrt(temp1));
/* Computing MAX */
#line 800 "cgejsv.f"
		d__1 = aatmax, d__2 = rwork[*n + p];
#line 800 "cgejsv.f"
		aatmax = max(d__1,d__2);
#line 801 "cgejsv.f"
		if (rwork[*n + p] != 0.) {
/* Computing MIN */
#line 801 "cgejsv.f"
		    d__1 = aatmin, d__2 = rwork[*n + p];
#line 801 "cgejsv.f"
		    aatmin = min(d__1,d__2);
#line 801 "cgejsv.f"
		}
#line 803 "cgejsv.f"
/* L1950: */
#line 803 "cgejsv.f"
	    }
#line 804 "cgejsv.f"
	} else {
#line 805 "cgejsv.f"
	    i__1 = *m;
#line 805 "cgejsv.f"
	    for (p = 1; p <= i__1; ++p) {
#line 806 "cgejsv.f"
		rwork[*m + *n + p] = scalem * z_abs(&a[p + isamax_(n, &a[p + 
			a_dim1], lda) * a_dim1]);
/* Computing MAX */
#line 807 "cgejsv.f"
		d__1 = aatmax, d__2 = rwork[*m + *n + p];
#line 807 "cgejsv.f"
		aatmax = max(d__1,d__2);
/* Computing MIN */
#line 808 "cgejsv.f"
		d__1 = aatmin, d__2 = rwork[*m + *n + p];
#line 808 "cgejsv.f"
		aatmin = min(d__1,d__2);
#line 809 "cgejsv.f"
/* L1904: */
#line 809 "cgejsv.f"
	    }
#line 810 "cgejsv.f"
	}

#line 812 "cgejsv.f"
    }

/*     For square matrix A try to determine whether A^*  would be  better */
/*     input for the preconditioned Jacobi SVD, with faster convergence. */
/*     The decision is based on an O(N) function of the vector of column */
/*     and row norms of A, based on the Shannon entropy. This should give */
/*     the right choice in most cases when the difference actually matters. */
/*     It may fail and pick the slower converging side. */

#line 821 "cgejsv.f"
    entra = 0.;
#line 822 "cgejsv.f"
    entrat = 0.;
#line 823 "cgejsv.f"
    if (l2tran) {

#line 825 "cgejsv.f"
	xsc = 0.;
#line 826 "cgejsv.f"
	temp1 = 1.;
#line 827 "cgejsv.f"
	classq_(n, &sva[1], &c__1, &xsc, &temp1);
#line 828 "cgejsv.f"
	temp1 = 1. / temp1;

#line 830 "cgejsv.f"
	entra = 0.;
#line 831 "cgejsv.f"
	i__1 = *n;
#line 831 "cgejsv.f"
	for (p = 1; p <= i__1; ++p) {
/* Computing 2nd power */
#line 832 "cgejsv.f"
	    d__1 = sva[p] / xsc;
#line 832 "cgejsv.f"
	    big1 = d__1 * d__1 * temp1;
#line 833 "cgejsv.f"
	    if (big1 != 0.) {
#line 833 "cgejsv.f"
		entra += big1 * log(big1);
#line 833 "cgejsv.f"
	    }
#line 834 "cgejsv.f"
/* L1113: */
#line 834 "cgejsv.f"
	}
#line 835 "cgejsv.f"
	entra = -entra / log((doublereal) (*n));

/*        Now, SVA().^2/Trace(A^* * A) is a point in the probability simplex. */
/*        It is derived from the diagonal of  A^* * A.  Do the same with the */
/*        diagonal of A * A^*, compute the entropy of the corresponding */
/*        probability distribution. Note that A * A^* and A^* * A have the */
/*        same trace. */

#line 843 "cgejsv.f"
	entrat = 0.;
#line 844 "cgejsv.f"
	i__1 = *n + *m;
#line 844 "cgejsv.f"
	for (p = *n + 1; p <= i__1; ++p) {
/* Computing 2nd power */
#line 845 "cgejsv.f"
	    d__1 = rwork[p] / xsc;
#line 845 "cgejsv.f"
	    big1 = d__1 * d__1 * temp1;
#line 846 "cgejsv.f"
	    if (big1 != 0.) {
#line 846 "cgejsv.f"
		entrat += big1 * log(big1);
#line 846 "cgejsv.f"
	    }
#line 847 "cgejsv.f"
/* L1114: */
#line 847 "cgejsv.f"
	}
#line 848 "cgejsv.f"
	entrat = -entrat / log((doublereal) (*m));

/*        Analyze the entropies and decide A or A^*. Smaller entropy */
/*        usually means better input for the algorithm. */

#line 853 "cgejsv.f"
	transp = entrat < entra;
#line 854 "cgejsv.f"
	transp = TRUE_;

/*        If A^* is better than A, take the adjoint of A. */

#line 858 "cgejsv.f"
	if (transp) {
/*           In an optimal implementation, this trivial transpose */
/*           should be replaced with faster transpose. */
#line 861 "cgejsv.f"
	    i__1 = *n - 1;
#line 861 "cgejsv.f"
	    for (p = 1; p <= i__1; ++p) {
#line 862 "cgejsv.f"
		i__2 = p + p * a_dim1;
#line 862 "cgejsv.f"
		d_cnjg(&z__1, &a[p + p * a_dim1]);
#line 862 "cgejsv.f"
		a[i__2].r = z__1.r, a[i__2].i = z__1.i;
#line 863 "cgejsv.f"
		i__2 = *n;
#line 863 "cgejsv.f"
		for (q = p + 1; q <= i__2; ++q) {
#line 864 "cgejsv.f"
		    d_cnjg(&z__1, &a[q + p * a_dim1]);
#line 864 "cgejsv.f"
		    ctemp.r = z__1.r, ctemp.i = z__1.i;
#line 865 "cgejsv.f"
		    i__3 = q + p * a_dim1;
#line 865 "cgejsv.f"
		    d_cnjg(&z__1, &a[p + q * a_dim1]);
#line 865 "cgejsv.f"
		    a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 866 "cgejsv.f"
		    i__3 = p + q * a_dim1;
#line 866 "cgejsv.f"
		    a[i__3].r = ctemp.r, a[i__3].i = ctemp.i;
#line 867 "cgejsv.f"
/* L1116: */
#line 867 "cgejsv.f"
		}
#line 868 "cgejsv.f"
/* L1115: */
#line 868 "cgejsv.f"
	    }
#line 869 "cgejsv.f"
	    i__1 = *n + *n * a_dim1;
#line 869 "cgejsv.f"
	    d_cnjg(&z__1, &a[*n + *n * a_dim1]);
#line 869 "cgejsv.f"
	    a[i__1].r = z__1.r, a[i__1].i = z__1.i;
#line 870 "cgejsv.f"
	    i__1 = *n;
#line 870 "cgejsv.f"
	    for (p = 1; p <= i__1; ++p) {
#line 871 "cgejsv.f"
		rwork[*m + *n + p] = sva[p];
#line 872 "cgejsv.f"
		sva[p] = rwork[*n + p];
/*              previously computed row 2-norms are now column 2-norms */
/*              of the transposed matrix */
#line 875 "cgejsv.f"
/* L1117: */
#line 875 "cgejsv.f"
	    }
#line 876 "cgejsv.f"
	    temp1 = aapp;
#line 877 "cgejsv.f"
	    aapp = aatmax;
#line 878 "cgejsv.f"
	    aatmax = temp1;
#line 879 "cgejsv.f"
	    temp1 = aaqq;
#line 880 "cgejsv.f"
	    aaqq = aatmin;
#line 881 "cgejsv.f"
	    aatmin = temp1;
#line 882 "cgejsv.f"
	    kill = lsvec;
#line 883 "cgejsv.f"
	    lsvec = rsvec;
#line 884 "cgejsv.f"
	    rsvec = kill;
#line 885 "cgejsv.f"
	    if (lsvec) {
#line 885 "cgejsv.f"
		n1 = *n;
#line 885 "cgejsv.f"
	    }

#line 887 "cgejsv.f"
	    rowpiv = TRUE_;
#line 888 "cgejsv.f"
	}

#line 890 "cgejsv.f"
    }
/*     END IF L2TRAN */

/*     Scale the matrix so that its maximal singular value remains less */
/*     than SQRT(BIG) -- the matrix is scaled so that its maximal column */
/*     has Euclidean norm equal to SQRT(BIG/N). The only reason to keep */
/*     SQRT(BIG) instead of BIG is the fact that CGEJSV uses LAPACK and */
/*     BLAS routines that, in some implementations, are not capable of */
/*     working in the full interval [SFMIN,BIG] and that they may provoke */
/*     overflows in the intermediate results. If the singular values spread */
/*     from SFMIN to BIG, then CGESVJ will compute them. So, in that case, */
/*     one should use CGESVJ instead of CGEJSV. */

#line 903 "cgejsv.f"
    big1 = sqrt(big);
#line 904 "cgejsv.f"
    temp1 = sqrt(big / (doublereal) (*n));

#line 906 "cgejsv.f"
    clascl_("G", &c__0, &c__0, &aapp, &temp1, n, &c__1, &sva[1], n, &ierr, (
	    ftnlen)1);
#line 907 "cgejsv.f"
    if (aaqq > aapp * sfmin) {
#line 908 "cgejsv.f"
	aaqq = aaqq / aapp * temp1;
#line 909 "cgejsv.f"
    } else {
#line 910 "cgejsv.f"
	aaqq = aaqq * temp1 / aapp;
#line 911 "cgejsv.f"
    }
#line 912 "cgejsv.f"
    temp1 *= scalem;
#line 913 "cgejsv.f"
    clascl_("G", &c__0, &c__0, &aapp, &temp1, m, n, &a[a_offset], lda, &ierr, 
	    (ftnlen)1);

/*     To undo scaling at the end of this procedure, multiply the */
/*     computed singular values with USCAL2 / USCAL1. */

#line 918 "cgejsv.f"
    uscal1 = temp1;
#line 919 "cgejsv.f"
    uscal2 = aapp;

#line 921 "cgejsv.f"
    if (l2kill) {
/*        L2KILL enforces computation of nonzero singular values in */
/*        the restricted range of condition number of the initial A, */
/*        sigma_max(A) / sigma_min(A) approx. SQRT(BIG)/SQRT(SFMIN). */
#line 925 "cgejsv.f"
	xsc = sqrt(sfmin);
#line 926 "cgejsv.f"
    } else {
#line 927 "cgejsv.f"
	xsc = small;

/*        Now, if the condition number of A is too big, */
/*        sigma_max(A) / sigma_min(A) .GT. SQRT(BIG/N) * EPSLN / SFMIN, */
/*        as a precaution measure, the full SVD is computed using CGESVJ */
/*        with accumulated Jacobi rotations. This provides numerically */
/*        more robust computation, at the cost of slightly increased run */
/*        time. Depending on the concrete implementation of BLAS and LAPACK */
/*        (i.e. how they behave in presence of extreme ill-conditioning) the */
/*        implementor may decide to remove this switch. */
#line 937 "cgejsv.f"
	if (aaqq < sqrt(sfmin) && lsvec && rsvec) {
#line 938 "cgejsv.f"
	    jracc = TRUE_;
#line 939 "cgejsv.f"
	}

#line 941 "cgejsv.f"
    }
#line 942 "cgejsv.f"
    if (aaqq < xsc) {
#line 943 "cgejsv.f"
	i__1 = *n;
#line 943 "cgejsv.f"
	for (p = 1; p <= i__1; ++p) {
#line 944 "cgejsv.f"
	    if (sva[p] < xsc) {
#line 945 "cgejsv.f"
		claset_("A", m, &c__1, &c_b1, &c_b1, &a[p * a_dim1 + 1], lda, 
			(ftnlen)1);
#line 946 "cgejsv.f"
		sva[p] = 0.;
#line 947 "cgejsv.f"
	    }
#line 948 "cgejsv.f"
/* L700: */
#line 948 "cgejsv.f"
	}
#line 949 "cgejsv.f"
    }

/*     Preconditioning using QR factorization with pivoting */

#line 953 "cgejsv.f"
    if (rowpiv) {
/*        Optional row permutation (Bjoerck row pivoting): */
/*        A result by Cox and Higham shows that the Bjoerck's */
/*        row pivoting combined with standard column pivoting */
/*        has similar effect as Powell-Reid complete pivoting. */
/*        The ell-infinity norms of A are made nonincreasing. */
#line 959 "cgejsv.f"
	i__1 = *m - 1;
#line 959 "cgejsv.f"
	for (p = 1; p <= i__1; ++p) {
#line 960 "cgejsv.f"
	    i__2 = *m - p + 1;
#line 960 "cgejsv.f"
	    q = isamax_(&i__2, &rwork[*m + *n + p], &c__1) + p - 1;
#line 961 "cgejsv.f"
	    iwork[(*n << 1) + p] = q;
#line 962 "cgejsv.f"
	    if (p != q) {
#line 963 "cgejsv.f"
		temp1 = rwork[*m + *n + p];
#line 964 "cgejsv.f"
		rwork[*m + *n + p] = rwork[*m + *n + q];
#line 965 "cgejsv.f"
		rwork[*m + *n + q] = temp1;
#line 966 "cgejsv.f"
	    }
#line 967 "cgejsv.f"
/* L1952: */
#line 967 "cgejsv.f"
	}
#line 968 "cgejsv.f"
	i__1 = *m - 1;
#line 968 "cgejsv.f"
	claswp_(n, &a[a_offset], lda, &c__1, &i__1, &iwork[(*n << 1) + 1], &
		c__1);
#line 969 "cgejsv.f"
    }

/*     End of the preparation phase (scaling, optional sorting and */
/*     transposing, optional flushing of small columns). */

/*     Preconditioning */

/*     If the full SVD is needed, the right singular vectors are computed */
/*     from a matrix equation, and for that we need theoretical analysis */
/*     of the Businger-Golub pivoting. So we use CGEQP3 as the first RR QRF. */
/*     In all other cases the first RR QRF can be chosen by other criteria */
/*     (eg speed by replacing global with restricted window pivoting, such */
/*     as in xGEQPX from TOMS # 782). Good results will be obtained using */
/*     xGEQPX with properly (!) chosen numerical parameters. */
/*     Any improvement of CGEQP3 improves overal performance of CGEJSV. */

/*     A * P1 = Q1 * [ R1^* 0]^*: */
#line 986 "cgejsv.f"
    i__1 = *n;
#line 986 "cgejsv.f"
    for (p = 1; p <= i__1; ++p) {
/*        .. all columns are free columns */
#line 988 "cgejsv.f"
	iwork[p] = 0;
#line 989 "cgejsv.f"
/* L1963: */
#line 989 "cgejsv.f"
    }
#line 990 "cgejsv.f"
    i__1 = *lwork - *n;
#line 990 "cgejsv.f"
    cgeqp3_(m, n, &a[a_offset], lda, &iwork[1], &cwork[1], &cwork[*n + 1], &
	    i__1, &rwork[1], &ierr);

/*     The upper triangular matrix R1 from the first QRF is inspected for */
/*     rank deficiency and possibilities for deflation, or possible */
/*     ill-conditioning. Depending on the user specified flag L2RANK, */
/*     the procedure explores possibilities to reduce the numerical */
/*     rank by inspecting the computed upper triangular factor. If */
/*     L2RANK or L2ABER are up, then CGEJSV will compute the SVD of */
/*     A + dA, where ||dA|| <= f(M,N)*EPSLN. */

#line 1001 "cgejsv.f"
    nr = 1;
#line 1002 "cgejsv.f"
    if (l2aber) {
/*        Standard absolute error bound suffices. All sigma_i with */
/*        sigma_i < N*EPSLN*||A|| are flushed to zero. This is an */
/*        agressive enforcement of lower numerical rank by introducing a */
/*        backward error of the order of N*EPSLN*||A||. */
#line 1007 "cgejsv.f"
	temp1 = sqrt((doublereal) (*n)) * epsln;
#line 1008 "cgejsv.f"
	i__1 = *n;
#line 1008 "cgejsv.f"
	for (p = 2; p <= i__1; ++p) {
#line 1009 "cgejsv.f"
	    if (z_abs(&a[p + p * a_dim1]) >= temp1 * z_abs(&a[a_dim1 + 1])) {
#line 1010 "cgejsv.f"
		++nr;
#line 1011 "cgejsv.f"
	    } else {
#line 1012 "cgejsv.f"
		goto L3002;
#line 1013 "cgejsv.f"
	    }
#line 1014 "cgejsv.f"
/* L3001: */
#line 1014 "cgejsv.f"
	}
#line 1015 "cgejsv.f"
L3002:
#line 1016 "cgejsv.f"
	;
#line 1016 "cgejsv.f"
    } else if (l2rank) {
/*        .. similarly as above, only slightly more gentle (less agressive). */
/*        Sudden drop on the diagonal of R1 is used as the criterion for */
/*        close-to-rank-defficient. */
#line 1020 "cgejsv.f"
	temp1 = sqrt(sfmin);
#line 1021 "cgejsv.f"
	i__1 = *n;
#line 1021 "cgejsv.f"
	for (p = 2; p <= i__1; ++p) {
#line 1022 "cgejsv.f"
	    if (z_abs(&a[p + p * a_dim1]) < epsln * z_abs(&a[p - 1 + (p - 1) *
		     a_dim1]) || z_abs(&a[p + p * a_dim1]) < small || l2kill 
		    && z_abs(&a[p + p * a_dim1]) < temp1) {
#line 1022 "cgejsv.f"
		goto L3402;
#line 1022 "cgejsv.f"
	    }
#line 1025 "cgejsv.f"
	    ++nr;
#line 1026 "cgejsv.f"
/* L3401: */
#line 1026 "cgejsv.f"
	}
#line 1027 "cgejsv.f"
L3402:

#line 1029 "cgejsv.f"
	;
#line 1029 "cgejsv.f"
    } else {
/*        The goal is high relative accuracy. However, if the matrix */
/*        has high scaled condition number the relative accuracy is in */
/*        general not feasible. Later on, a condition number estimator */
/*        will be deployed to estimate the scaled condition number. */
/*        Here we just remove the underflowed part of the triangular */
/*        factor. This prevents the situation in which the code is */
/*        working hard to get the accuracy not warranted by the data. */
#line 1037 "cgejsv.f"
	temp1 = sqrt(sfmin);
#line 1038 "cgejsv.f"
	i__1 = *n;
#line 1038 "cgejsv.f"
	for (p = 2; p <= i__1; ++p) {
#line 1039 "cgejsv.f"
	    if (z_abs(&a[p + p * a_dim1]) < small || l2kill && z_abs(&a[p + p 
		    * a_dim1]) < temp1) {
#line 1039 "cgejsv.f"
		goto L3302;
#line 1039 "cgejsv.f"
	    }
#line 1041 "cgejsv.f"
	    ++nr;
#line 1042 "cgejsv.f"
/* L3301: */
#line 1042 "cgejsv.f"
	}
#line 1043 "cgejsv.f"
L3302:

#line 1045 "cgejsv.f"
	;
#line 1045 "cgejsv.f"
    }

#line 1047 "cgejsv.f"
    almort = FALSE_;
#line 1048 "cgejsv.f"
    if (nr == *n) {
#line 1049 "cgejsv.f"
	maxprj = 1.;
#line 1050 "cgejsv.f"
	i__1 = *n;
#line 1050 "cgejsv.f"
	for (p = 2; p <= i__1; ++p) {
#line 1051 "cgejsv.f"
	    temp1 = z_abs(&a[p + p * a_dim1]) / sva[iwork[p]];
#line 1052 "cgejsv.f"
	    maxprj = min(maxprj,temp1);
#line 1053 "cgejsv.f"
/* L3051: */
#line 1053 "cgejsv.f"
	}
/* Computing 2nd power */
#line 1054 "cgejsv.f"
	d__1 = maxprj;
#line 1054 "cgejsv.f"
	if (d__1 * d__1 >= 1. - (doublereal) (*n) * epsln) {
#line 1054 "cgejsv.f"
	    almort = TRUE_;
#line 1054 "cgejsv.f"
	}
#line 1055 "cgejsv.f"
    }


#line 1058 "cgejsv.f"
    sconda = -1.;
#line 1059 "cgejsv.f"
    condr1 = -1.;
#line 1060 "cgejsv.f"
    condr2 = -1.;

#line 1062 "cgejsv.f"
    if (errest) {
#line 1063 "cgejsv.f"
	if (*n == nr) {
#line 1064 "cgejsv.f"
	    if (rsvec) {
/*              .. V is available as workspace */
#line 1066 "cgejsv.f"
		clacpy_("U", n, n, &a[a_offset], lda, &v[v_offset], ldv, (
			ftnlen)1);
#line 1067 "cgejsv.f"
		i__1 = *n;
#line 1067 "cgejsv.f"
		for (p = 1; p <= i__1; ++p) {
#line 1068 "cgejsv.f"
		    temp1 = sva[iwork[p]];
#line 1069 "cgejsv.f"
		    d__1 = 1. / temp1;
#line 1069 "cgejsv.f"
		    csscal_(&p, &d__1, &v[p * v_dim1 + 1], &c__1);
#line 1070 "cgejsv.f"
/* L3053: */
#line 1070 "cgejsv.f"
		}
#line 1071 "cgejsv.f"
		cpocon_("U", n, &v[v_offset], ldv, &c_b80, &temp1, &cwork[*n 
			+ 1], &rwork[1], &ierr, (ftnlen)1);

#line 1074 "cgejsv.f"
	    } else if (lsvec) {
/*              .. U is available as workspace */
#line 1076 "cgejsv.f"
		clacpy_("U", n, n, &a[a_offset], lda, &u[u_offset], ldu, (
			ftnlen)1);
#line 1077 "cgejsv.f"
		i__1 = *n;
#line 1077 "cgejsv.f"
		for (p = 1; p <= i__1; ++p) {
#line 1078 "cgejsv.f"
		    temp1 = sva[iwork[p]];
#line 1079 "cgejsv.f"
		    d__1 = 1. / temp1;
#line 1079 "cgejsv.f"
		    csscal_(&p, &d__1, &u[p * u_dim1 + 1], &c__1);
#line 1080 "cgejsv.f"
/* L3054: */
#line 1080 "cgejsv.f"
		}
#line 1081 "cgejsv.f"
		cpocon_("U", n, &u[u_offset], ldu, &c_b80, &temp1, &cwork[*n 
			+ 1], &rwork[1], &ierr, (ftnlen)1);
#line 1083 "cgejsv.f"
	    } else {
#line 1084 "cgejsv.f"
		clacpy_("U", n, n, &a[a_offset], lda, &cwork[*n + 1], n, (
			ftnlen)1);
#line 1085 "cgejsv.f"
		i__1 = *n;
#line 1085 "cgejsv.f"
		for (p = 1; p <= i__1; ++p) {
#line 1086 "cgejsv.f"
		    temp1 = sva[iwork[p]];
#line 1087 "cgejsv.f"
		    d__1 = 1. / temp1;
#line 1087 "cgejsv.f"
		    csscal_(&p, &d__1, &cwork[*n + (p - 1) * *n + 1], &c__1);
#line 1088 "cgejsv.f"
/* L3052: */
#line 1088 "cgejsv.f"
		}
/*           .. the columns of R are scaled to have unit Euclidean lengths. */
#line 1090 "cgejsv.f"
		cpocon_("U", n, &cwork[*n + 1], n, &c_b80, &temp1, &cwork[*n 
			+ *n * *n + 1], &rwork[1], &ierr, (ftnlen)1);

#line 1093 "cgejsv.f"
	    }
#line 1094 "cgejsv.f"
	    sconda = 1. / sqrt(temp1);
/*           SCONDA is an estimate of SQRT(||(R^* * R)^(-1)||_1). */
/*           N^(-1/4) * SCONDA <= ||R^(-1)||_2 <= N^(1/4) * SCONDA */
#line 1097 "cgejsv.f"
	} else {
#line 1098 "cgejsv.f"
	    sconda = -1.;
#line 1099 "cgejsv.f"
	}
#line 1100 "cgejsv.f"
    }

#line 1102 "cgejsv.f"
    z_div(&z__1, &a[a_dim1 + 1], &a[nr + nr * a_dim1]);
#line 1102 "cgejsv.f"
    l2pert = l2pert && z_abs(&z__1) > sqrt(big1);
/*     If there is no violent scaling, artificial perturbation is not needed. */

/*     Phase 3: */

#line 1107 "cgejsv.f"
    if (! (rsvec || lsvec)) {

/*         Singular Values only */

/*         .. transpose A(1:NR,1:N) */
/* Computing MIN */
#line 1112 "cgejsv.f"
	i__2 = *n - 1;
#line 1112 "cgejsv.f"
	i__1 = min(i__2,nr);
#line 1112 "cgejsv.f"
	for (p = 1; p <= i__1; ++p) {
#line 1113 "cgejsv.f"
	    i__2 = *n - p;
#line 1113 "cgejsv.f"
	    ccopy_(&i__2, &a[p + (p + 1) * a_dim1], lda, &a[p + 1 + p * 
		    a_dim1], &c__1);
#line 1114 "cgejsv.f"
	    i__2 = *n - p + 1;
#line 1114 "cgejsv.f"
	    clacgv_(&i__2, &a[p + p * a_dim1], &c__1);
#line 1115 "cgejsv.f"
/* L1946: */
#line 1115 "cgejsv.f"
	}
#line 1116 "cgejsv.f"
	if (nr == *n) {
#line 1116 "cgejsv.f"
	    i__1 = *n + *n * a_dim1;
#line 1116 "cgejsv.f"
	    d_cnjg(&z__1, &a[*n + *n * a_dim1]);
#line 1116 "cgejsv.f"
	    a[i__1].r = z__1.r, a[i__1].i = z__1.i;
#line 1116 "cgejsv.f"
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

#line 1130 "cgejsv.f"
	if (! almort) {

#line 1132 "cgejsv.f"
	    if (l2pert) {
/*              XSC = SQRT(SMALL) */
#line 1134 "cgejsv.f"
		xsc = epsln / (doublereal) (*n);
#line 1135 "cgejsv.f"
		i__1 = nr;
#line 1135 "cgejsv.f"
		for (q = 1; q <= i__1; ++q) {
#line 1136 "cgejsv.f"
		    d__1 = xsc * z_abs(&a[q + q * a_dim1]);
#line 1136 "cgejsv.f"
		    z__1.r = d__1, z__1.i = 0.;
#line 1136 "cgejsv.f"
		    ctemp.r = z__1.r, ctemp.i = z__1.i;
#line 1137 "cgejsv.f"
		    i__2 = *n;
#line 1137 "cgejsv.f"
		    for (p = 1; p <= i__2; ++p) {
#line 1138 "cgejsv.f"
			if (p > q && z_abs(&a[p + q * a_dim1]) <= temp1 || p <
				 q) {
#line 1138 "cgejsv.f"
			    i__3 = p + q * a_dim1;
#line 1138 "cgejsv.f"
			    a[i__3].r = ctemp.r, a[i__3].i = ctemp.i;
#line 1138 "cgejsv.f"
			}
/*     $                     A(p,q) = TEMP1 * ( A(p,q) / ABS(A(p,q)) ) */
#line 1142 "cgejsv.f"
/* L4949: */
#line 1142 "cgejsv.f"
		    }
#line 1143 "cgejsv.f"
/* L4947: */
#line 1143 "cgejsv.f"
		}
#line 1144 "cgejsv.f"
	    } else {
#line 1145 "cgejsv.f"
		i__1 = nr - 1;
#line 1145 "cgejsv.f"
		i__2 = nr - 1;
#line 1145 "cgejsv.f"
		claset_("U", &i__1, &i__2, &c_b1, &c_b1, &a[(a_dim1 << 1) + 1]
			, lda, (ftnlen)1);
#line 1146 "cgejsv.f"
	    }

/*            .. second preconditioning using the QR factorization */

#line 1150 "cgejsv.f"
	    i__1 = *lwork - *n;
#line 1150 "cgejsv.f"
	    cgeqrf_(n, &nr, &a[a_offset], lda, &cwork[1], &cwork[*n + 1], &
		    i__1, &ierr);

/*           .. and transpose upper to lower triangular */
#line 1153 "cgejsv.f"
	    i__1 = nr - 1;
#line 1153 "cgejsv.f"
	    for (p = 1; p <= i__1; ++p) {
#line 1154 "cgejsv.f"
		i__2 = nr - p;
#line 1154 "cgejsv.f"
		ccopy_(&i__2, &a[p + (p + 1) * a_dim1], lda, &a[p + 1 + p * 
			a_dim1], &c__1);
#line 1155 "cgejsv.f"
		i__2 = nr - p + 1;
#line 1155 "cgejsv.f"
		clacgv_(&i__2, &a[p + p * a_dim1], &c__1);
#line 1156 "cgejsv.f"
/* L1948: */
#line 1156 "cgejsv.f"
	    }

#line 1158 "cgejsv.f"
	}

/*           Row-cyclic Jacobi SVD algorithm with column pivoting */

/*           .. again some perturbation (a "background noise") is added */
/*           to drown denormals */
#line 1164 "cgejsv.f"
	if (l2pert) {
/*              XSC = SQRT(SMALL) */
#line 1166 "cgejsv.f"
	    xsc = epsln / (doublereal) (*n);
#line 1167 "cgejsv.f"
	    i__1 = nr;
#line 1167 "cgejsv.f"
	    for (q = 1; q <= i__1; ++q) {
#line 1168 "cgejsv.f"
		d__1 = xsc * z_abs(&a[q + q * a_dim1]);
#line 1168 "cgejsv.f"
		z__1.r = d__1, z__1.i = 0.;
#line 1168 "cgejsv.f"
		ctemp.r = z__1.r, ctemp.i = z__1.i;
#line 1169 "cgejsv.f"
		i__2 = nr;
#line 1169 "cgejsv.f"
		for (p = 1; p <= i__2; ++p) {
#line 1170 "cgejsv.f"
		    if (p > q && z_abs(&a[p + q * a_dim1]) <= temp1 || p < q) 
			    {
#line 1170 "cgejsv.f"
			i__3 = p + q * a_dim1;
#line 1170 "cgejsv.f"
			a[i__3].r = ctemp.r, a[i__3].i = ctemp.i;
#line 1170 "cgejsv.f"
		    }
/*     $                   A(p,q) = TEMP1 * ( A(p,q) / ABS(A(p,q)) ) */
#line 1174 "cgejsv.f"
/* L1949: */
#line 1174 "cgejsv.f"
		}
#line 1175 "cgejsv.f"
/* L1947: */
#line 1175 "cgejsv.f"
	    }
#line 1176 "cgejsv.f"
	} else {
#line 1177 "cgejsv.f"
	    i__1 = nr - 1;
#line 1177 "cgejsv.f"
	    i__2 = nr - 1;
#line 1177 "cgejsv.f"
	    claset_("U", &i__1, &i__2, &c_b1, &c_b1, &a[(a_dim1 << 1) + 1], 
		    lda, (ftnlen)1);
#line 1178 "cgejsv.f"
	}

/*           .. and one-sided Jacobi rotations are started on a lower */
/*           triangular matrix (plus perturbation which is ignored in */
/*           the part which destroys triangular form (confusing?!)) */

#line 1184 "cgejsv.f"
	cgesvj_("L", "NoU", "NoV", &nr, &nr, &a[a_offset], lda, &sva[1], n, &
		v[v_offset], ldv, &cwork[1], lwork, &rwork[1], lrwork, info, (
		ftnlen)1, (ftnlen)3, (ftnlen)3);

#line 1187 "cgejsv.f"
	scalem = rwork[1];
#line 1188 "cgejsv.f"
	numrank = i_dnnt(&rwork[2]);


#line 1191 "cgejsv.f"
    } else if (rsvec && ! lsvec) {

/*        -> Singular Values and Right Singular Vectors <- */

#line 1195 "cgejsv.f"
	if (almort) {

/*           .. in this case NR equals N */
#line 1198 "cgejsv.f"
	    i__1 = nr;
#line 1198 "cgejsv.f"
	    for (p = 1; p <= i__1; ++p) {
#line 1199 "cgejsv.f"
		i__2 = *n - p + 1;
#line 1199 "cgejsv.f"
		ccopy_(&i__2, &a[p + p * a_dim1], lda, &v[p + p * v_dim1], &
			c__1);
#line 1200 "cgejsv.f"
		i__2 = *n - p + 1;
#line 1200 "cgejsv.f"
		clacgv_(&i__2, &v[p + p * v_dim1], &c__1);
#line 1201 "cgejsv.f"
/* L1998: */
#line 1201 "cgejsv.f"
	    }
#line 1202 "cgejsv.f"
	    i__1 = nr - 1;
#line 1202 "cgejsv.f"
	    i__2 = nr - 1;
#line 1202 "cgejsv.f"
	    claset_("Upper", &i__1, &i__2, &c_b1, &c_b1, &v[(v_dim1 << 1) + 1]
		    , ldv, (ftnlen)5);

#line 1204 "cgejsv.f"
	    cgesvj_("L", "U", "N", n, &nr, &v[v_offset], ldv, &sva[1], &nr, &
		    a[a_offset], lda, &cwork[1], lwork, &rwork[1], lrwork, 
		    info, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1206 "cgejsv.f"
	    scalem = rwork[1];
#line 1207 "cgejsv.f"
	    numrank = i_dnnt(&rwork[2]);
#line 1209 "cgejsv.f"
	} else {

/*        .. two more QR factorizations ( one QRF is not enough, two require */
/*        accumulated product of Jacobi rotations, three are perfect ) */

#line 1214 "cgejsv.f"
	    i__1 = nr - 1;
#line 1214 "cgejsv.f"
	    i__2 = nr - 1;
#line 1214 "cgejsv.f"
	    claset_("Lower", &i__1, &i__2, &c_b1, &c_b1, &a[a_dim1 + 2], lda, 
		    (ftnlen)5);
#line 1215 "cgejsv.f"
	    i__1 = *lwork - *n;
#line 1215 "cgejsv.f"
	    cgelqf_(&nr, n, &a[a_offset], lda, &cwork[1], &cwork[*n + 1], &
		    i__1, &ierr);
#line 1216 "cgejsv.f"
	    clacpy_("Lower", &nr, &nr, &a[a_offset], lda, &v[v_offset], ldv, (
		    ftnlen)5);
#line 1217 "cgejsv.f"
	    i__1 = nr - 1;
#line 1217 "cgejsv.f"
	    i__2 = nr - 1;
#line 1217 "cgejsv.f"
	    claset_("Upper", &i__1, &i__2, &c_b1, &c_b1, &v[(v_dim1 << 1) + 1]
		    , ldv, (ftnlen)5);
#line 1218 "cgejsv.f"
	    i__1 = *lwork - (*n << 1);
#line 1218 "cgejsv.f"
	    cgeqrf_(&nr, &nr, &v[v_offset], ldv, &cwork[*n + 1], &cwork[(*n <<
		     1) + 1], &i__1, &ierr);
#line 1220 "cgejsv.f"
	    i__1 = nr;
#line 1220 "cgejsv.f"
	    for (p = 1; p <= i__1; ++p) {
#line 1221 "cgejsv.f"
		i__2 = nr - p + 1;
#line 1221 "cgejsv.f"
		ccopy_(&i__2, &v[p + p * v_dim1], ldv, &v[p + p * v_dim1], &
			c__1);
#line 1222 "cgejsv.f"
		i__2 = nr - p + 1;
#line 1222 "cgejsv.f"
		clacgv_(&i__2, &v[p + p * v_dim1], &c__1);
#line 1223 "cgejsv.f"
/* L8998: */
#line 1223 "cgejsv.f"
	    }
#line 1224 "cgejsv.f"
	    i__1 = nr - 1;
#line 1224 "cgejsv.f"
	    i__2 = nr - 1;
#line 1224 "cgejsv.f"
	    claset_("Upper", &i__1, &i__2, &c_b120, &c_b120, &v[(v_dim1 << 1) 
		    + 1], ldv, (ftnlen)5);

#line 1226 "cgejsv.f"
	    i__1 = *lwork - *n;
#line 1226 "cgejsv.f"
	    cgesvj_("Lower", "U", "N", &nr, &nr, &v[v_offset], ldv, &sva[1], &
		    nr, &u[u_offset], ldu, &cwork[*n + 1], &i__1, &rwork[1], 
		    lrwork, info, (ftnlen)5, (ftnlen)1, (ftnlen)1);
#line 1228 "cgejsv.f"
	    scalem = rwork[1];
#line 1229 "cgejsv.f"
	    numrank = i_dnnt(&rwork[2]);
#line 1230 "cgejsv.f"
	    if (nr < *n) {
#line 1231 "cgejsv.f"
		i__1 = *n - nr;
#line 1231 "cgejsv.f"
		claset_("A", &i__1, &nr, &c_b1, &c_b1, &v[nr + 1 + v_dim1], 
			ldv, (ftnlen)1);
#line 1232 "cgejsv.f"
		i__1 = *n - nr;
#line 1232 "cgejsv.f"
		claset_("A", &nr, &i__1, &c_b1, &c_b1, &v[(nr + 1) * v_dim1 + 
			1], ldv, (ftnlen)1);
#line 1233 "cgejsv.f"
		i__1 = *n - nr;
#line 1233 "cgejsv.f"
		i__2 = *n - nr;
#line 1233 "cgejsv.f"
		claset_("A", &i__1, &i__2, &c_b1, &c_b2, &v[nr + 1 + (nr + 1) 
			* v_dim1], ldv, (ftnlen)1);
#line 1234 "cgejsv.f"
	    }

#line 1236 "cgejsv.f"
	    i__1 = *lwork - *n;
#line 1236 "cgejsv.f"
	    cunmlq_("Left", "C", n, n, &nr, &a[a_offset], lda, &cwork[1], &v[
		    v_offset], ldv, &cwork[*n + 1], &i__1, &ierr, (ftnlen)4, (
		    ftnlen)1);

#line 1239 "cgejsv.f"
	}

#line 1241 "cgejsv.f"
	i__1 = *n;
#line 1241 "cgejsv.f"
	for (p = 1; p <= i__1; ++p) {
#line 1242 "cgejsv.f"
	    ccopy_(n, &v[p + v_dim1], ldv, &a[iwork[p] + a_dim1], lda);
#line 1243 "cgejsv.f"
/* L8991: */
#line 1243 "cgejsv.f"
	}
#line 1244 "cgejsv.f"
	clacpy_("All", n, n, &a[a_offset], lda, &v[v_offset], ldv, (ftnlen)3);

#line 1246 "cgejsv.f"
	if (transp) {
#line 1247 "cgejsv.f"
	    clacpy_("All", n, n, &v[v_offset], ldv, &u[u_offset], ldu, (
		    ftnlen)3);
#line 1248 "cgejsv.f"
	}

#line 1250 "cgejsv.f"
    } else if (lsvec && ! rsvec) {

/*        .. Singular Values and Left Singular Vectors                 .. */

/*        .. second preconditioning step to avoid need to accumulate */
/*        Jacobi rotations in the Jacobi iterations. */
#line 1256 "cgejsv.f"
	i__1 = nr;
#line 1256 "cgejsv.f"
	for (p = 1; p <= i__1; ++p) {
#line 1257 "cgejsv.f"
	    i__2 = *n - p + 1;
#line 1257 "cgejsv.f"
	    ccopy_(&i__2, &a[p + p * a_dim1], lda, &u[p + p * u_dim1], &c__1);
#line 1258 "cgejsv.f"
	    i__2 = *n - p + 1;
#line 1258 "cgejsv.f"
	    clacgv_(&i__2, &u[p + p * u_dim1], &c__1);
#line 1259 "cgejsv.f"
/* L1965: */
#line 1259 "cgejsv.f"
	}
#line 1260 "cgejsv.f"
	i__1 = nr - 1;
#line 1260 "cgejsv.f"
	i__2 = nr - 1;
#line 1260 "cgejsv.f"
	claset_("Upper", &i__1, &i__2, &c_b1, &c_b1, &u[(u_dim1 << 1) + 1], 
		ldu, (ftnlen)5);

#line 1262 "cgejsv.f"
	i__1 = *lwork - (*n << 1);
#line 1262 "cgejsv.f"
	cgeqrf_(n, &nr, &u[u_offset], ldu, &cwork[*n + 1], &cwork[(*n << 1) + 
		1], &i__1, &ierr);

#line 1265 "cgejsv.f"
	i__1 = nr - 1;
#line 1265 "cgejsv.f"
	for (p = 1; p <= i__1; ++p) {
#line 1266 "cgejsv.f"
	    i__2 = nr - p;
#line 1266 "cgejsv.f"
	    ccopy_(&i__2, &u[p + (p + 1) * u_dim1], ldu, &u[p + 1 + p * 
		    u_dim1], &c__1);
#line 1267 "cgejsv.f"
	    i__2 = *n - p + 1;
#line 1267 "cgejsv.f"
	    clacgv_(&i__2, &u[p + p * u_dim1], &c__1);
#line 1268 "cgejsv.f"
/* L1967: */
#line 1268 "cgejsv.f"
	}
#line 1269 "cgejsv.f"
	i__1 = nr - 1;
#line 1269 "cgejsv.f"
	i__2 = nr - 1;
#line 1269 "cgejsv.f"
	claset_("Upper", &i__1, &i__2, &c_b1, &c_b1, &u[(u_dim1 << 1) + 1], 
		ldu, (ftnlen)5);

#line 1271 "cgejsv.f"
	i__1 = *lwork - *n;
#line 1271 "cgejsv.f"
	cgesvj_("Lower", "U", "N", &nr, &nr, &u[u_offset], ldu, &sva[1], &nr, 
		&a[a_offset], lda, &cwork[*n + 1], &i__1, &rwork[1], lrwork, 
		info, (ftnlen)5, (ftnlen)1, (ftnlen)1);
#line 1273 "cgejsv.f"
	scalem = rwork[1];
#line 1274 "cgejsv.f"
	numrank = i_dnnt(&rwork[2]);

#line 1276 "cgejsv.f"
	if (nr < *m) {
#line 1277 "cgejsv.f"
	    i__1 = *m - nr;
#line 1277 "cgejsv.f"
	    claset_("A", &i__1, &nr, &c_b1, &c_b1, &u[nr + 1 + u_dim1], ldu, (
		    ftnlen)1);
#line 1278 "cgejsv.f"
	    if (nr < n1) {
#line 1279 "cgejsv.f"
		i__1 = n1 - nr;
#line 1279 "cgejsv.f"
		claset_("A", &nr, &i__1, &c_b1, &c_b1, &u[(nr + 1) * u_dim1 + 
			1], ldu, (ftnlen)1);
#line 1280 "cgejsv.f"
		i__1 = *m - nr;
#line 1280 "cgejsv.f"
		i__2 = n1 - nr;
#line 1280 "cgejsv.f"
		claset_("A", &i__1, &i__2, &c_b1, &c_b2, &u[nr + 1 + (nr + 1) 
			* u_dim1], ldu, (ftnlen)1);
#line 1281 "cgejsv.f"
	    }
#line 1282 "cgejsv.f"
	}

#line 1284 "cgejsv.f"
	i__1 = *lwork - *n;
#line 1284 "cgejsv.f"
	cunmqr_("Left", "No Tr", m, &n1, n, &a[a_offset], lda, &cwork[1], &u[
		u_offset], ldu, &cwork[*n + 1], &i__1, &ierr, (ftnlen)4, (
		ftnlen)5);

#line 1287 "cgejsv.f"
	if (rowpiv) {
#line 1287 "cgejsv.f"
	    i__1 = *m - 1;
#line 1287 "cgejsv.f"
	    claswp_(&n1, &u[u_offset], ldu, &c__1, &i__1, &iwork[(*n << 1) + 
		    1], &c_n1);
#line 1287 "cgejsv.f"
	}

#line 1290 "cgejsv.f"
	i__1 = n1;
#line 1290 "cgejsv.f"
	for (p = 1; p <= i__1; ++p) {
#line 1291 "cgejsv.f"
	    xsc = 1. / scnrm2_(m, &u[p * u_dim1 + 1], &c__1);
#line 1292 "cgejsv.f"
	    csscal_(m, &xsc, &u[p * u_dim1 + 1], &c__1);
#line 1293 "cgejsv.f"
/* L1974: */
#line 1293 "cgejsv.f"
	}

#line 1295 "cgejsv.f"
	if (transp) {
#line 1296 "cgejsv.f"
	    clacpy_("All", n, n, &u[u_offset], ldu, &v[v_offset], ldv, (
		    ftnlen)3);
#line 1297 "cgejsv.f"
	}

#line 1299 "cgejsv.f"
    } else {

/*        .. Full SVD .. */

#line 1303 "cgejsv.f"
	if (! jracc) {

#line 1305 "cgejsv.f"
	    if (! almort) {

/*           Second Preconditioning Step (QRF [with pivoting]) */
/*           Note that the composition of TRANSPOSE, QRF and TRANSPOSE is */
/*           equivalent to an LQF CALL. Since in many libraries the QRF */
/*           seems to be better optimized than the LQF, we do explicit */
/*           transpose and use the QRF. This is subject to changes in an */
/*           optimized implementation of CGEJSV. */

#line 1314 "cgejsv.f"
		i__1 = nr;
#line 1314 "cgejsv.f"
		for (p = 1; p <= i__1; ++p) {
#line 1315 "cgejsv.f"
		    i__2 = *n - p + 1;
#line 1315 "cgejsv.f"
		    ccopy_(&i__2, &a[p + p * a_dim1], lda, &v[p + p * v_dim1],
			     &c__1);
#line 1316 "cgejsv.f"
		    i__2 = *n - p + 1;
#line 1316 "cgejsv.f"
		    clacgv_(&i__2, &v[p + p * v_dim1], &c__1);
#line 1317 "cgejsv.f"
/* L1968: */
#line 1317 "cgejsv.f"
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

#line 1331 "cgejsv.f"
		if (l2pert) {
#line 1332 "cgejsv.f"
		    xsc = sqrt(small);
#line 1333 "cgejsv.f"
		    i__1 = nr;
#line 1333 "cgejsv.f"
		    for (q = 1; q <= i__1; ++q) {
#line 1334 "cgejsv.f"
			d__1 = xsc * z_abs(&v[q + q * v_dim1]);
#line 1334 "cgejsv.f"
			z__1.r = d__1, z__1.i = 0.;
#line 1334 "cgejsv.f"
			ctemp.r = z__1.r, ctemp.i = z__1.i;
#line 1335 "cgejsv.f"
			i__2 = *n;
#line 1335 "cgejsv.f"
			for (p = 1; p <= i__2; ++p) {
#line 1336 "cgejsv.f"
			    if (p > q && z_abs(&v[p + q * v_dim1]) <= temp1 ||
				     p < q) {
#line 1336 "cgejsv.f"
				i__3 = p + q * v_dim1;
#line 1336 "cgejsv.f"
				v[i__3].r = ctemp.r, v[i__3].i = ctemp.i;
#line 1336 "cgejsv.f"
			    }
/*     $                   V(p,q) = TEMP1 * ( V(p,q) / ABS(V(p,q)) ) */
#line 1340 "cgejsv.f"
			    if (p < q) {
#line 1340 "cgejsv.f"
				i__3 = p + q * v_dim1;
#line 1340 "cgejsv.f"
				i__4 = p + q * v_dim1;
#line 1340 "cgejsv.f"
				z__1.r = -v[i__4].r, z__1.i = -v[i__4].i;
#line 1340 "cgejsv.f"
				v[i__3].r = z__1.r, v[i__3].i = z__1.i;
#line 1340 "cgejsv.f"
			    }
#line 1341 "cgejsv.f"
/* L2968: */
#line 1341 "cgejsv.f"
			}
#line 1342 "cgejsv.f"
/* L2969: */
#line 1342 "cgejsv.f"
		    }
#line 1343 "cgejsv.f"
		} else {
#line 1344 "cgejsv.f"
		    i__1 = nr - 1;
#line 1344 "cgejsv.f"
		    i__2 = nr - 1;
#line 1344 "cgejsv.f"
		    claset_("U", &i__1, &i__2, &c_b1, &c_b1, &v[(v_dim1 << 1) 
			    + 1], ldv, (ftnlen)1);
#line 1345 "cgejsv.f"
		}

/*           Estimate the row scaled condition number of R1 */
/*           (If R1 is rectangular, N > NR, then the condition number */
/*           of the leading NR x NR submatrix is estimated.) */

#line 1351 "cgejsv.f"
		clacpy_("L", &nr, &nr, &v[v_offset], ldv, &cwork[(*n << 1) + 
			1], &nr, (ftnlen)1);
#line 1352 "cgejsv.f"
		i__1 = nr;
#line 1352 "cgejsv.f"
		for (p = 1; p <= i__1; ++p) {
#line 1353 "cgejsv.f"
		    i__2 = nr - p + 1;
#line 1353 "cgejsv.f"
		    temp1 = scnrm2_(&i__2, &cwork[(*n << 1) + (p - 1) * nr + 
			    p], &c__1);
#line 1354 "cgejsv.f"
		    i__2 = nr - p + 1;
#line 1354 "cgejsv.f"
		    d__1 = 1. / temp1;
#line 1354 "cgejsv.f"
		    csscal_(&i__2, &d__1, &cwork[(*n << 1) + (p - 1) * nr + p]
			    , &c__1);
#line 1355 "cgejsv.f"
/* L3950: */
#line 1355 "cgejsv.f"
		}
#line 1356 "cgejsv.f"
		cpocon_("Lower", &nr, &cwork[(*n << 1) + 1], &nr, &c_b80, &
			temp1, &cwork[(*n << 1) + nr * nr + 1], &rwork[1], &
			ierr, (ftnlen)5);
#line 1358 "cgejsv.f"
		condr1 = 1. / sqrt(temp1);
/*           .. here need a second oppinion on the condition number */
/*           .. then assume worst case scenario */
/*           R1 is OK for inverse <=> CONDR1 .LT. FLOAT(N) */
/*           more conservative    <=> CONDR1 .LT. SQRT(FLOAT(N)) */

#line 1364 "cgejsv.f"
		cond_ok__ = sqrt(sqrt((doublereal) nr));
/* [TP]       COND_OK is a tuning parameter. */

#line 1367 "cgejsv.f"
		if (condr1 < cond_ok__) {
/*              .. the second QRF without pivoting. Note: in an optimized */
/*              implementation, this QRF should be implemented as the QRF */
/*              of a lower triangular matrix. */
/*              R1^* = Q2 * R2 */
#line 1372 "cgejsv.f"
		    i__1 = *lwork - (*n << 1);
#line 1372 "cgejsv.f"
		    cgeqrf_(n, &nr, &v[v_offset], ldv, &cwork[*n + 1], &cwork[
			    (*n << 1) + 1], &i__1, &ierr);

#line 1375 "cgejsv.f"
		    if (l2pert) {
#line 1376 "cgejsv.f"
			xsc = sqrt(small) / epsln;
#line 1377 "cgejsv.f"
			i__1 = nr;
#line 1377 "cgejsv.f"
			for (p = 2; p <= i__1; ++p) {
#line 1378 "cgejsv.f"
			    i__2 = p - 1;
#line 1378 "cgejsv.f"
			    for (q = 1; q <= i__2; ++q) {
/* Computing MIN */
#line 1379 "cgejsv.f"
				d__2 = z_abs(&v[p + p * v_dim1]), d__3 = 
					z_abs(&v[q + q * v_dim1]);
#line 1379 "cgejsv.f"
				d__1 = xsc * min(d__2,d__3);
#line 1379 "cgejsv.f"
				z__1.r = d__1, z__1.i = 0.;
#line 1379 "cgejsv.f"
				ctemp.r = z__1.r, ctemp.i = z__1.i;
#line 1381 "cgejsv.f"
				if (z_abs(&v[q + p * v_dim1]) <= temp1) {
#line 1381 "cgejsv.f"
				    i__3 = q + p * v_dim1;
#line 1381 "cgejsv.f"
				    v[i__3].r = ctemp.r, v[i__3].i = ctemp.i;
#line 1381 "cgejsv.f"
				}
/*     $                     V(q,p) = TEMP1 * ( V(q,p) / ABS(V(q,p)) ) */
#line 1384 "cgejsv.f"
/* L3958: */
#line 1384 "cgejsv.f"
			    }
#line 1385 "cgejsv.f"
/* L3959: */
#line 1385 "cgejsv.f"
			}
#line 1386 "cgejsv.f"
		    }

#line 1388 "cgejsv.f"
		    if (nr != *n) {
#line 1388 "cgejsv.f"
			clacpy_("A", n, &nr, &v[v_offset], ldv, &cwork[(*n << 
				1) + 1], n, (ftnlen)1);
#line 1388 "cgejsv.f"
		    }
/*              .. save ... */

/*           .. this transposed copy should be better than naive */
#line 1393 "cgejsv.f"
		    i__1 = nr - 1;
#line 1393 "cgejsv.f"
		    for (p = 1; p <= i__1; ++p) {
#line 1394 "cgejsv.f"
			i__2 = nr - p;
#line 1394 "cgejsv.f"
			ccopy_(&i__2, &v[p + (p + 1) * v_dim1], ldv, &v[p + 1 
				+ p * v_dim1], &c__1);
#line 1395 "cgejsv.f"
			i__2 = nr - p + 1;
#line 1395 "cgejsv.f"
			clacgv_(&i__2, &v[p + p * v_dim1], &c__1);
#line 1396 "cgejsv.f"
/* L1969: */
#line 1396 "cgejsv.f"
		    }
#line 1397 "cgejsv.f"
		    i__1 = nr + nr * v_dim1;
#line 1397 "cgejsv.f"
		    d_cnjg(&z__1, &v[nr + nr * v_dim1]);
#line 1397 "cgejsv.f"
		    v[i__1].r = z__1.r, v[i__1].i = z__1.i;

#line 1399 "cgejsv.f"
		    condr2 = condr1;

#line 1401 "cgejsv.f"
		} else {

/*              .. ill-conditioned case: second QRF with pivoting */
/*              Note that windowed pivoting would be equaly good */
/*              numerically, and more run-time efficient. So, in */
/*              an optimal implementation, the next call to CGEQP3 */
/*              should be replaced with eg. CALL CGEQPX (ACM TOMS #782) */
/*              with properly (carefully) chosen parameters. */

/*              R1^* * P2 = Q2 * R2 */
#line 1411 "cgejsv.f"
		    i__1 = nr;
#line 1411 "cgejsv.f"
		    for (p = 1; p <= i__1; ++p) {
#line 1412 "cgejsv.f"
			iwork[*n + p] = 0;
#line 1413 "cgejsv.f"
/* L3003: */
#line 1413 "cgejsv.f"
		    }
#line 1414 "cgejsv.f"
		    i__1 = *lwork - (*n << 1);
#line 1414 "cgejsv.f"
		    cgeqp3_(n, &nr, &v[v_offset], ldv, &iwork[*n + 1], &cwork[
			    *n + 1], &cwork[(*n << 1) + 1], &i__1, &rwork[1], 
			    &ierr);
/* *               CALL CGEQRF( N, NR, V, LDV, CWORK(N+1), CWORK(2*N+1), */
/* *     $              LWORK-2*N, IERR ) */
#line 1418 "cgejsv.f"
		    if (l2pert) {
#line 1419 "cgejsv.f"
			xsc = sqrt(small);
#line 1420 "cgejsv.f"
			i__1 = nr;
#line 1420 "cgejsv.f"
			for (p = 2; p <= i__1; ++p) {
#line 1421 "cgejsv.f"
			    i__2 = p - 1;
#line 1421 "cgejsv.f"
			    for (q = 1; q <= i__2; ++q) {
/* Computing MIN */
#line 1422 "cgejsv.f"
				d__2 = z_abs(&v[p + p * v_dim1]), d__3 = 
					z_abs(&v[q + q * v_dim1]);
#line 1422 "cgejsv.f"
				d__1 = xsc * min(d__2,d__3);
#line 1422 "cgejsv.f"
				z__1.r = d__1, z__1.i = 0.;
#line 1422 "cgejsv.f"
				ctemp.r = z__1.r, ctemp.i = z__1.i;
#line 1424 "cgejsv.f"
				if (z_abs(&v[q + p * v_dim1]) <= temp1) {
#line 1424 "cgejsv.f"
				    i__3 = q + p * v_dim1;
#line 1424 "cgejsv.f"
				    v[i__3].r = ctemp.r, v[i__3].i = ctemp.i;
#line 1424 "cgejsv.f"
				}
/*     $                     V(q,p) = TEMP1 * ( V(q,p) / ABS(V(q,p)) ) */
#line 1427 "cgejsv.f"
/* L3968: */
#line 1427 "cgejsv.f"
			    }
#line 1428 "cgejsv.f"
/* L3969: */
#line 1428 "cgejsv.f"
			}
#line 1429 "cgejsv.f"
		    }

#line 1431 "cgejsv.f"
		    clacpy_("A", n, &nr, &v[v_offset], ldv, &cwork[(*n << 1) 
			    + 1], n, (ftnlen)1);

#line 1433 "cgejsv.f"
		    if (l2pert) {
#line 1434 "cgejsv.f"
			xsc = sqrt(small);
#line 1435 "cgejsv.f"
			i__1 = nr;
#line 1435 "cgejsv.f"
			for (p = 2; p <= i__1; ++p) {
#line 1436 "cgejsv.f"
			    i__2 = p - 1;
#line 1436 "cgejsv.f"
			    for (q = 1; q <= i__2; ++q) {
/* Computing MIN */
#line 1437 "cgejsv.f"
				d__2 = z_abs(&v[p + p * v_dim1]), d__3 = 
					z_abs(&v[q + q * v_dim1]);
#line 1437 "cgejsv.f"
				d__1 = xsc * min(d__2,d__3);
#line 1437 "cgejsv.f"
				z__1.r = d__1, z__1.i = 0.;
#line 1437 "cgejsv.f"
				ctemp.r = z__1.r, ctemp.i = z__1.i;
/*                        V(p,q) = - TEMP1*( V(q,p) / ABS(V(q,p)) ) */
#line 1440 "cgejsv.f"
				i__3 = p + q * v_dim1;
#line 1440 "cgejsv.f"
				z__1.r = -ctemp.r, z__1.i = -ctemp.i;
#line 1440 "cgejsv.f"
				v[i__3].r = z__1.r, v[i__3].i = z__1.i;
#line 1441 "cgejsv.f"
/* L8971: */
#line 1441 "cgejsv.f"
			    }
#line 1442 "cgejsv.f"
/* L8970: */
#line 1442 "cgejsv.f"
			}
#line 1443 "cgejsv.f"
		    } else {
#line 1444 "cgejsv.f"
			i__1 = nr - 1;
#line 1444 "cgejsv.f"
			i__2 = nr - 1;
#line 1444 "cgejsv.f"
			claset_("L", &i__1, &i__2, &c_b1, &c_b1, &v[v_dim1 + 
				2], ldv, (ftnlen)1);
#line 1445 "cgejsv.f"
		    }
/*              Now, compute R2 = L3 * Q3, the LQ factorization. */
#line 1447 "cgejsv.f"
		    i__1 = *lwork - (*n << 1) - *n * nr - nr;
#line 1447 "cgejsv.f"
		    cgelqf_(&nr, &nr, &v[v_offset], ldv, &cwork[(*n << 1) + *
			    n * nr + 1], &cwork[(*n << 1) + *n * nr + nr + 1],
			     &i__1, &ierr);
/*              .. and estimate the condition number */
#line 1450 "cgejsv.f"
		    clacpy_("L", &nr, &nr, &v[v_offset], ldv, &cwork[(*n << 1)
			     + *n * nr + nr + 1], &nr, (ftnlen)1);
#line 1451 "cgejsv.f"
		    i__1 = nr;
#line 1451 "cgejsv.f"
		    for (p = 1; p <= i__1; ++p) {
#line 1452 "cgejsv.f"
			temp1 = scnrm2_(&p, &cwork[(*n << 1) + *n * nr + nr + 
				p], &nr);
#line 1453 "cgejsv.f"
			d__1 = 1. / temp1;
#line 1453 "cgejsv.f"
			csscal_(&p, &d__1, &cwork[(*n << 1) + *n * nr + nr + 
				p], &nr);
#line 1454 "cgejsv.f"
/* L4950: */
#line 1454 "cgejsv.f"
		    }
#line 1455 "cgejsv.f"
		    cpocon_("L", &nr, &cwork[(*n << 1) + *n * nr + nr + 1], &
			    nr, &c_b80, &temp1, &cwork[(*n << 1) + *n * nr + 
			    nr + nr * nr + 1], &rwork[1], &ierr, (ftnlen)1);
#line 1457 "cgejsv.f"
		    condr2 = 1. / sqrt(temp1);


#line 1460 "cgejsv.f"
		    if (condr2 >= cond_ok__) {
/*                 .. save the Householder vectors used for Q3 */
/*                 (this overwrittes the copy of R2, as it will not be */
/*                 needed in this branch, but it does not overwritte the */
/*                 Huseholder vectors of Q2.). */
#line 1465 "cgejsv.f"
			clacpy_("U", &nr, &nr, &v[v_offset], ldv, &cwork[(*n 
				<< 1) + 1], n, (ftnlen)1);
/*                 .. and the rest of the information on Q3 is in */
/*                 WORK(2*N+N*NR+1:2*N+N*NR+N) */
#line 1468 "cgejsv.f"
		    }

#line 1470 "cgejsv.f"
		}

#line 1472 "cgejsv.f"
		if (l2pert) {
#line 1473 "cgejsv.f"
		    xsc = sqrt(small);
#line 1474 "cgejsv.f"
		    i__1 = nr;
#line 1474 "cgejsv.f"
		    for (q = 2; q <= i__1; ++q) {
#line 1475 "cgejsv.f"
			i__2 = q + q * v_dim1;
#line 1475 "cgejsv.f"
			z__1.r = xsc * v[i__2].r, z__1.i = xsc * v[i__2].i;
#line 1475 "cgejsv.f"
			ctemp.r = z__1.r, ctemp.i = z__1.i;
#line 1476 "cgejsv.f"
			i__2 = q - 1;
#line 1476 "cgejsv.f"
			for (p = 1; p <= i__2; ++p) {
/*                    V(p,q) = - SIGN( TEMP1, V(q,p) ) */
/*                     V(p,q) = - TEMP1*( V(p,q) / ABS(V(p,q)) ) */
#line 1479 "cgejsv.f"
			    i__3 = p + q * v_dim1;
#line 1479 "cgejsv.f"
			    z__1.r = -ctemp.r, z__1.i = -ctemp.i;
#line 1479 "cgejsv.f"
			    v[i__3].r = z__1.r, v[i__3].i = z__1.i;
#line 1480 "cgejsv.f"
/* L4969: */
#line 1480 "cgejsv.f"
			}
#line 1481 "cgejsv.f"
/* L4968: */
#line 1481 "cgejsv.f"
		    }
#line 1482 "cgejsv.f"
		} else {
#line 1483 "cgejsv.f"
		    i__1 = nr - 1;
#line 1483 "cgejsv.f"
		    i__2 = nr - 1;
#line 1483 "cgejsv.f"
		    claset_("U", &i__1, &i__2, &c_b1, &c_b1, &v[(v_dim1 << 1) 
			    + 1], ldv, (ftnlen)1);
#line 1484 "cgejsv.f"
		}

/*        Second preconditioning finished; continue with Jacobi SVD */
/*        The input matrix is lower trinagular. */

/*        Recover the right singular vectors as solution of a well */
/*        conditioned triangular matrix equation. */

#line 1492 "cgejsv.f"
		if (condr1 < cond_ok__) {

#line 1494 "cgejsv.f"
		    i__1 = *lwork - (*n << 1) - *n * nr - nr;
#line 1494 "cgejsv.f"
		    cgesvj_("L", "U", "N", &nr, &nr, &v[v_offset], ldv, &sva[
			    1], &nr, &u[u_offset], ldu, &cwork[(*n << 1) + *n 
			    * nr + nr + 1], &i__1, &rwork[1], lrwork, info, (
			    ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1497 "cgejsv.f"
		    scalem = rwork[1];
#line 1498 "cgejsv.f"
		    numrank = i_dnnt(&rwork[2]);
#line 1499 "cgejsv.f"
		    i__1 = nr;
#line 1499 "cgejsv.f"
		    for (p = 1; p <= i__1; ++p) {
#line 1500 "cgejsv.f"
			ccopy_(&nr, &v[p * v_dim1 + 1], &c__1, &u[p * u_dim1 
				+ 1], &c__1);
#line 1501 "cgejsv.f"
			csscal_(&nr, &sva[p], &v[p * v_dim1 + 1], &c__1);
#line 1502 "cgejsv.f"
/* L3970: */
#line 1502 "cgejsv.f"
		    }
/*        .. pick the right matrix equation and solve it */

#line 1506 "cgejsv.f"
		    if (nr == *n) {
/* :))             .. best case, R1 is inverted. The solution of this matrix */
/*                 equation is Q2*V2 = the product of the Jacobi rotations */
/*                 used in CGESVJ, premultiplied with the orthogonal matrix */
/*                 from the second QR factorization. */
#line 1511 "cgejsv.f"
			ctrsm_("L", "U", "N", "N", &nr, &nr, &c_b2, &a[
				a_offset], lda, &v[v_offset], ldv, (ftnlen)1, 
				(ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1512 "cgejsv.f"
		    } else {
/*                 .. R1 is well conditioned, but non-square. Adjoint of R2 */
/*                 is inverted to get the product of the Jacobi rotations */
/*                 used in CGESVJ. The Q-factor from the second QR */
/*                 factorization is then built in explicitly. */
#line 1517 "cgejsv.f"
			ctrsm_("L", "U", "C", "N", &nr, &nr, &c_b2, &cwork[(*
				n << 1) + 1], n, &v[v_offset], ldv, (ftnlen)1,
				 (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1519 "cgejsv.f"
			if (nr < *n) {
#line 1520 "cgejsv.f"
			    i__1 = *n - nr;
#line 1520 "cgejsv.f"
			    claset_("A", &i__1, &nr, &c_b120, &c_b1, &v[nr + 
				    1 + v_dim1], ldv, (ftnlen)1);
#line 1521 "cgejsv.f"
			    i__1 = *n - nr;
#line 1521 "cgejsv.f"
			    claset_("A", &nr, &i__1, &c_b120, &c_b1, &v[(nr + 
				    1) * v_dim1 + 1], ldv, (ftnlen)1);
#line 1522 "cgejsv.f"
			    i__1 = *n - nr;
#line 1522 "cgejsv.f"
			    i__2 = *n - nr;
#line 1522 "cgejsv.f"
			    claset_("A", &i__1, &i__2, &c_b120, &c_b2, &v[nr 
				    + 1 + (nr + 1) * v_dim1], ldv, (ftnlen)1);
#line 1523 "cgejsv.f"
			}
#line 1524 "cgejsv.f"
			i__1 = *lwork - (*n << 1) - *n * nr - nr;
#line 1524 "cgejsv.f"
			cunmqr_("L", "N", n, n, &nr, &cwork[(*n << 1) + 1], n,
				 &cwork[*n + 1], &v[v_offset], ldv, &cwork[(*
				n << 1) + *n * nr + nr + 1], &i__1, &ierr, (
				ftnlen)1, (ftnlen)1);
#line 1526 "cgejsv.f"
		    }

#line 1528 "cgejsv.f"
		} else if (condr2 < cond_ok__) {

/*              The matrix R2 is inverted. The solution of the matrix equation */
/*              is Q3^* * V3 = the product of the Jacobi rotations (appplied to */
/*              the lower triangular L3 from the LQ factorization of */
/*              R2=L3*Q3), pre-multiplied with the transposed Q3. */
#line 1534 "cgejsv.f"
		    i__1 = *lwork - (*n << 1) - *n * nr - nr;
#line 1534 "cgejsv.f"
		    cgesvj_("L", "U", "N", &nr, &nr, &v[v_offset], ldv, &sva[
			    1], &nr, &u[u_offset], ldu, &cwork[(*n << 1) + *n 
			    * nr + nr + 1], &i__1, &rwork[1], lrwork, info, (
			    ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1537 "cgejsv.f"
		    scalem = rwork[1];
#line 1538 "cgejsv.f"
		    numrank = i_dnnt(&rwork[2]);
#line 1539 "cgejsv.f"
		    i__1 = nr;
#line 1539 "cgejsv.f"
		    for (p = 1; p <= i__1; ++p) {
#line 1540 "cgejsv.f"
			ccopy_(&nr, &v[p * v_dim1 + 1], &c__1, &u[p * u_dim1 
				+ 1], &c__1);
#line 1541 "cgejsv.f"
			csscal_(&nr, &sva[p], &u[p * u_dim1 + 1], &c__1);
#line 1542 "cgejsv.f"
/* L3870: */
#line 1542 "cgejsv.f"
		    }
#line 1543 "cgejsv.f"
		    ctrsm_("L", "U", "N", "N", &nr, &nr, &c_b2, &cwork[(*n << 
			    1) + 1], n, &u[u_offset], ldu, (ftnlen)1, (ftnlen)
			    1, (ftnlen)1, (ftnlen)1);
/*              .. apply the permutation from the second QR factorization */
#line 1546 "cgejsv.f"
		    i__1 = nr;
#line 1546 "cgejsv.f"
		    for (q = 1; q <= i__1; ++q) {
#line 1547 "cgejsv.f"
			i__2 = nr;
#line 1547 "cgejsv.f"
			for (p = 1; p <= i__2; ++p) {
#line 1548 "cgejsv.f"
			    i__3 = (*n << 1) + *n * nr + nr + iwork[*n + p];
#line 1548 "cgejsv.f"
			    i__4 = p + q * u_dim1;
#line 1548 "cgejsv.f"
			    cwork[i__3].r = u[i__4].r, cwork[i__3].i = u[i__4]
				    .i;
#line 1549 "cgejsv.f"
/* L872: */
#line 1549 "cgejsv.f"
			}
#line 1550 "cgejsv.f"
			i__2 = nr;
#line 1550 "cgejsv.f"
			for (p = 1; p <= i__2; ++p) {
#line 1551 "cgejsv.f"
			    i__3 = p + q * u_dim1;
#line 1551 "cgejsv.f"
			    i__4 = (*n << 1) + *n * nr + nr + p;
#line 1551 "cgejsv.f"
			    u[i__3].r = cwork[i__4].r, u[i__3].i = cwork[i__4]
				    .i;
#line 1552 "cgejsv.f"
/* L874: */
#line 1552 "cgejsv.f"
			}
#line 1553 "cgejsv.f"
/* L873: */
#line 1553 "cgejsv.f"
		    }
#line 1554 "cgejsv.f"
		    if (nr < *n) {
#line 1555 "cgejsv.f"
			i__1 = *n - nr;
#line 1555 "cgejsv.f"
			claset_("A", &i__1, &nr, &c_b1, &c_b1, &v[nr + 1 + 
				v_dim1], ldv, (ftnlen)1);
#line 1556 "cgejsv.f"
			i__1 = *n - nr;
#line 1556 "cgejsv.f"
			claset_("A", &nr, &i__1, &c_b1, &c_b1, &v[(nr + 1) * 
				v_dim1 + 1], ldv, (ftnlen)1);
#line 1557 "cgejsv.f"
			i__1 = *n - nr;
#line 1557 "cgejsv.f"
			i__2 = *n - nr;
#line 1557 "cgejsv.f"
			claset_("A", &i__1, &i__2, &c_b1, &c_b2, &v[nr + 1 + (
				nr + 1) * v_dim1], ldv, (ftnlen)1);
#line 1558 "cgejsv.f"
		    }
#line 1559 "cgejsv.f"
		    i__1 = *lwork - (*n << 1) - *n * nr - nr;
#line 1559 "cgejsv.f"
		    cunmqr_("L", "N", n, n, &nr, &cwork[(*n << 1) + 1], n, &
			    cwork[*n + 1], &v[v_offset], ldv, &cwork[(*n << 1)
			     + *n * nr + nr + 1], &i__1, &ierr, (ftnlen)1, (
			    ftnlen)1);
#line 1561 "cgejsv.f"
		} else {
/*              Last line of defense. */
/* #:(          This is a rather pathological case: no scaled condition */
/*              improvement after two pivoted QR factorizations. Other */
/*              possibility is that the rank revealing QR factorization */
/*              or the condition estimator has failed, or the COND_OK */
/*              is set very close to ONE (which is unnecessary). Normally, */
/*              this branch should never be executed, but in rare cases of */
/*              failure of the RRQR or condition estimator, the last line of */
/*              defense ensures that CGEJSV completes the task. */
/*              Compute the full SVD of L3 using CGESVJ with explicit */
/*              accumulation of Jacobi rotations. */
#line 1573 "cgejsv.f"
		    i__1 = *lwork - (*n << 1) - *n * nr - nr;
#line 1573 "cgejsv.f"
		    cgesvj_("L", "U", "V", &nr, &nr, &v[v_offset], ldv, &sva[
			    1], &nr, &u[u_offset], ldu, &cwork[(*n << 1) + *n 
			    * nr + nr + 1], &i__1, &rwork[1], lrwork, info, (
			    ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1576 "cgejsv.f"
		    scalem = rwork[1];
#line 1577 "cgejsv.f"
		    numrank = i_dnnt(&rwork[2]);
#line 1578 "cgejsv.f"
		    if (nr < *n) {
#line 1579 "cgejsv.f"
			i__1 = *n - nr;
#line 1579 "cgejsv.f"
			claset_("A", &i__1, &nr, &c_b1, &c_b1, &v[nr + 1 + 
				v_dim1], ldv, (ftnlen)1);
#line 1580 "cgejsv.f"
			i__1 = *n - nr;
#line 1580 "cgejsv.f"
			claset_("A", &nr, &i__1, &c_b1, &c_b1, &v[(nr + 1) * 
				v_dim1 + 1], ldv, (ftnlen)1);
#line 1581 "cgejsv.f"
			i__1 = *n - nr;
#line 1581 "cgejsv.f"
			i__2 = *n - nr;
#line 1581 "cgejsv.f"
			claset_("A", &i__1, &i__2, &c_b1, &c_b2, &v[nr + 1 + (
				nr + 1) * v_dim1], ldv, (ftnlen)1);
#line 1582 "cgejsv.f"
		    }
#line 1583 "cgejsv.f"
		    i__1 = *lwork - (*n << 1) - *n * nr - nr;
#line 1583 "cgejsv.f"
		    cunmqr_("L", "N", n, n, &nr, &cwork[(*n << 1) + 1], n, &
			    cwork[*n + 1], &v[v_offset], ldv, &cwork[(*n << 1)
			     + *n * nr + nr + 1], &i__1, &ierr, (ftnlen)1, (
			    ftnlen)1);

#line 1586 "cgejsv.f"
		    i__1 = *lwork - (*n << 1) - *n * nr - nr;
#line 1586 "cgejsv.f"
		    cunmlq_("L", "C", &nr, &nr, &nr, &cwork[(*n << 1) + 1], n,
			     &cwork[(*n << 1) + *n * nr + 1], &u[u_offset], 
			    ldu, &cwork[(*n << 1) + *n * nr + nr + 1], &i__1, 
			    &ierr, (ftnlen)1, (ftnlen)1);
#line 1589 "cgejsv.f"
		    i__1 = nr;
#line 1589 "cgejsv.f"
		    for (q = 1; q <= i__1; ++q) {
#line 1590 "cgejsv.f"
			i__2 = nr;
#line 1590 "cgejsv.f"
			for (p = 1; p <= i__2; ++p) {
#line 1591 "cgejsv.f"
			    i__3 = (*n << 1) + *n * nr + nr + iwork[*n + p];
#line 1591 "cgejsv.f"
			    i__4 = p + q * u_dim1;
#line 1591 "cgejsv.f"
			    cwork[i__3].r = u[i__4].r, cwork[i__3].i = u[i__4]
				    .i;
#line 1592 "cgejsv.f"
/* L772: */
#line 1592 "cgejsv.f"
			}
#line 1593 "cgejsv.f"
			i__2 = nr;
#line 1593 "cgejsv.f"
			for (p = 1; p <= i__2; ++p) {
#line 1594 "cgejsv.f"
			    i__3 = p + q * u_dim1;
#line 1594 "cgejsv.f"
			    i__4 = (*n << 1) + *n * nr + nr + p;
#line 1594 "cgejsv.f"
			    u[i__3].r = cwork[i__4].r, u[i__3].i = cwork[i__4]
				    .i;
#line 1595 "cgejsv.f"
/* L774: */
#line 1595 "cgejsv.f"
			}
#line 1596 "cgejsv.f"
/* L773: */
#line 1596 "cgejsv.f"
		    }

#line 1598 "cgejsv.f"
		}

/*           Permute the rows of V using the (column) permutation from the */
/*           first QRF. Also, scale the columns to make them unit in */
/*           Euclidean norm. This applies to all cases. */

#line 1604 "cgejsv.f"
		temp1 = sqrt((doublereal) (*n)) * epsln;
#line 1605 "cgejsv.f"
		i__1 = *n;
#line 1605 "cgejsv.f"
		for (q = 1; q <= i__1; ++q) {
#line 1606 "cgejsv.f"
		    i__2 = *n;
#line 1606 "cgejsv.f"
		    for (p = 1; p <= i__2; ++p) {
#line 1607 "cgejsv.f"
			i__3 = (*n << 1) + *n * nr + nr + iwork[p];
#line 1607 "cgejsv.f"
			i__4 = p + q * v_dim1;
#line 1607 "cgejsv.f"
			cwork[i__3].r = v[i__4].r, cwork[i__3].i = v[i__4].i;
#line 1608 "cgejsv.f"
/* L972: */
#line 1608 "cgejsv.f"
		    }
#line 1609 "cgejsv.f"
		    i__2 = *n;
#line 1609 "cgejsv.f"
		    for (p = 1; p <= i__2; ++p) {
#line 1610 "cgejsv.f"
			i__3 = p + q * v_dim1;
#line 1610 "cgejsv.f"
			i__4 = (*n << 1) + *n * nr + nr + p;
#line 1610 "cgejsv.f"
			v[i__3].r = cwork[i__4].r, v[i__3].i = cwork[i__4].i;
#line 1611 "cgejsv.f"
/* L973: */
#line 1611 "cgejsv.f"
		    }
#line 1612 "cgejsv.f"
		    xsc = 1. / scnrm2_(n, &v[q * v_dim1 + 1], &c__1);
#line 1613 "cgejsv.f"
		    if (xsc < 1. - temp1 || xsc > temp1 + 1.) {
#line 1613 "cgejsv.f"
			csscal_(n, &xsc, &v[q * v_dim1 + 1], &c__1);
#line 1613 "cgejsv.f"
		    }
#line 1615 "cgejsv.f"
/* L1972: */
#line 1615 "cgejsv.f"
		}
/*           At this moment, V contains the right singular vectors of A. */
/*           Next, assemble the left singular vector matrix U (M x N). */
#line 1618 "cgejsv.f"
		if (nr < *m) {
#line 1619 "cgejsv.f"
		    i__1 = *m - nr;
#line 1619 "cgejsv.f"
		    claset_("A", &i__1, &nr, &c_b1, &c_b1, &u[nr + 1 + u_dim1]
			    , ldu, (ftnlen)1);
#line 1620 "cgejsv.f"
		    if (nr < n1) {
#line 1621 "cgejsv.f"
			i__1 = n1 - nr;
#line 1621 "cgejsv.f"
			claset_("A", &nr, &i__1, &c_b1, &c_b1, &u[(nr + 1) * 
				u_dim1 + 1], ldu, (ftnlen)1);
#line 1622 "cgejsv.f"
			i__1 = *m - nr;
#line 1622 "cgejsv.f"
			i__2 = n1 - nr;
#line 1622 "cgejsv.f"
			claset_("A", &i__1, &i__2, &c_b1, &c_b2, &u[nr + 1 + (
				nr + 1) * u_dim1], ldu, (ftnlen)1);
#line 1624 "cgejsv.f"
		    }
#line 1625 "cgejsv.f"
		}

/*           The Q matrix from the first QRF is built into the left singular */
/*           matrix U. This applies to all cases. */

#line 1630 "cgejsv.f"
		i__1 = *lwork - *n;
#line 1630 "cgejsv.f"
		cunmqr_("Left", "No_Tr", m, &n1, n, &a[a_offset], lda, &cwork[
			1], &u[u_offset], ldu, &cwork[*n + 1], &i__1, &ierr, (
			ftnlen)4, (ftnlen)5);
/*           The columns of U are normalized. The cost is O(M*N) flops. */
#line 1634 "cgejsv.f"
		temp1 = sqrt((doublereal) (*m)) * epsln;
#line 1635 "cgejsv.f"
		i__1 = nr;
#line 1635 "cgejsv.f"
		for (p = 1; p <= i__1; ++p) {
#line 1636 "cgejsv.f"
		    xsc = 1. / scnrm2_(m, &u[p * u_dim1 + 1], &c__1);
#line 1637 "cgejsv.f"
		    if (xsc < 1. - temp1 || xsc > temp1 + 1.) {
#line 1637 "cgejsv.f"
			csscal_(m, &xsc, &u[p * u_dim1 + 1], &c__1);
#line 1637 "cgejsv.f"
		    }
#line 1639 "cgejsv.f"
/* L1973: */
#line 1639 "cgejsv.f"
		}

/*           If the initial QRF is computed with row pivoting, the left */
/*           singular vectors must be adjusted. */

#line 1644 "cgejsv.f"
		if (rowpiv) {
#line 1644 "cgejsv.f"
		    i__1 = *m - 1;
#line 1644 "cgejsv.f"
		    claswp_(&n1, &u[u_offset], ldu, &c__1, &i__1, &iwork[(*n 
			    << 1) + 1], &c_n1);
#line 1644 "cgejsv.f"
		}

#line 1647 "cgejsv.f"
	    } else {

/*        .. the initial matrix A has almost orthogonal columns and */
/*        the second QRF is not needed */

#line 1652 "cgejsv.f"
		clacpy_("Upper", n, n, &a[a_offset], lda, &cwork[*n + 1], n, (
			ftnlen)5);
#line 1653 "cgejsv.f"
		if (l2pert) {
#line 1654 "cgejsv.f"
		    xsc = sqrt(small);
#line 1655 "cgejsv.f"
		    i__1 = *n;
#line 1655 "cgejsv.f"
		    for (p = 2; p <= i__1; ++p) {
#line 1656 "cgejsv.f"
			i__2 = *n + (p - 1) * *n + p;
#line 1656 "cgejsv.f"
			z__1.r = xsc * cwork[i__2].r, z__1.i = xsc * cwork[
				i__2].i;
#line 1656 "cgejsv.f"
			ctemp.r = z__1.r, ctemp.i = z__1.i;
#line 1657 "cgejsv.f"
			i__2 = p - 1;
#line 1657 "cgejsv.f"
			for (q = 1; q <= i__2; ++q) {
/*                     CWORK(N+(q-1)*N+p)=-TEMP1 * ( CWORK(N+(p-1)*N+q) / */
/*     $                                        ABS(CWORK(N+(p-1)*N+q)) ) */
#line 1660 "cgejsv.f"
			    i__3 = *n + (q - 1) * *n + p;
#line 1660 "cgejsv.f"
			    z__1.r = -ctemp.r, z__1.i = -ctemp.i;
#line 1660 "cgejsv.f"
			    cwork[i__3].r = z__1.r, cwork[i__3].i = z__1.i;
#line 1661 "cgejsv.f"
/* L5971: */
#line 1661 "cgejsv.f"
			}
#line 1662 "cgejsv.f"
/* L5970: */
#line 1662 "cgejsv.f"
		    }
#line 1663 "cgejsv.f"
		} else {
#line 1664 "cgejsv.f"
		    i__1 = *n - 1;
#line 1664 "cgejsv.f"
		    i__2 = *n - 1;
#line 1664 "cgejsv.f"
		    claset_("Lower", &i__1, &i__2, &c_b1, &c_b1, &cwork[*n + 
			    2], n, (ftnlen)5);
#line 1665 "cgejsv.f"
		}

#line 1667 "cgejsv.f"
		i__1 = *lwork - *n - *n * *n;
#line 1667 "cgejsv.f"
		cgesvj_("Upper", "U", "N", n, n, &cwork[*n + 1], n, &sva[1], 
			n, &u[u_offset], ldu, &cwork[*n + *n * *n + 1], &i__1,
			 &rwork[1], lrwork, info, (ftnlen)5, (ftnlen)1, (
			ftnlen)1);

#line 1671 "cgejsv.f"
		scalem = rwork[1];
#line 1672 "cgejsv.f"
		numrank = i_dnnt(&rwork[2]);
#line 1673 "cgejsv.f"
		i__1 = *n;
#line 1673 "cgejsv.f"
		for (p = 1; p <= i__1; ++p) {
#line 1674 "cgejsv.f"
		    ccopy_(n, &cwork[*n + (p - 1) * *n + 1], &c__1, &u[p * 
			    u_dim1 + 1], &c__1);
#line 1675 "cgejsv.f"
		    csscal_(n, &sva[p], &cwork[*n + (p - 1) * *n + 1], &c__1);
#line 1676 "cgejsv.f"
/* L6970: */
#line 1676 "cgejsv.f"
		}

#line 1678 "cgejsv.f"
		ctrsm_("Left", "Upper", "NoTrans", "No UD", n, n, &c_b2, &a[
			a_offset], lda, &cwork[*n + 1], n, (ftnlen)4, (ftnlen)
			5, (ftnlen)7, (ftnlen)5);
#line 1680 "cgejsv.f"
		i__1 = *n;
#line 1680 "cgejsv.f"
		for (p = 1; p <= i__1; ++p) {
#line 1681 "cgejsv.f"
		    ccopy_(n, &cwork[*n + p], n, &v[iwork[p] + v_dim1], ldv);
#line 1682 "cgejsv.f"
/* L6972: */
#line 1682 "cgejsv.f"
		}
#line 1683 "cgejsv.f"
		temp1 = sqrt((doublereal) (*n)) * epsln;
#line 1684 "cgejsv.f"
		i__1 = *n;
#line 1684 "cgejsv.f"
		for (p = 1; p <= i__1; ++p) {
#line 1685 "cgejsv.f"
		    xsc = 1. / scnrm2_(n, &v[p * v_dim1 + 1], &c__1);
#line 1686 "cgejsv.f"
		    if (xsc < 1. - temp1 || xsc > temp1 + 1.) {
#line 1686 "cgejsv.f"
			csscal_(n, &xsc, &v[p * v_dim1 + 1], &c__1);
#line 1686 "cgejsv.f"
		    }
#line 1688 "cgejsv.f"
/* L6971: */
#line 1688 "cgejsv.f"
		}

/*           Assemble the left singular vector matrix U (M x N). */

#line 1692 "cgejsv.f"
		if (*n < *m) {
#line 1693 "cgejsv.f"
		    i__1 = *m - *n;
#line 1693 "cgejsv.f"
		    claset_("A", &i__1, n, &c_b1, &c_b1, &u[*n + 1 + u_dim1], 
			    ldu, (ftnlen)1);
#line 1694 "cgejsv.f"
		    if (*n < n1) {
#line 1695 "cgejsv.f"
			i__1 = n1 - *n;
#line 1695 "cgejsv.f"
			claset_("A", n, &i__1, &c_b1, &c_b1, &u[(*n + 1) * 
				u_dim1 + 1], ldu, (ftnlen)1);
#line 1696 "cgejsv.f"
			i__1 = *m - *n;
#line 1696 "cgejsv.f"
			i__2 = n1 - *n;
#line 1696 "cgejsv.f"
			claset_("A", &i__1, &i__2, &c_b1, &c_b2, &u[*n + 1 + (
				*n + 1) * u_dim1], ldu, (ftnlen)1);
#line 1697 "cgejsv.f"
		    }
#line 1698 "cgejsv.f"
		}
#line 1699 "cgejsv.f"
		i__1 = *lwork - *n;
#line 1699 "cgejsv.f"
		cunmqr_("Left", "No Tr", m, &n1, n, &a[a_offset], lda, &cwork[
			1], &u[u_offset], ldu, &cwork[*n + 1], &i__1, &ierr, (
			ftnlen)4, (ftnlen)5);
#line 1701 "cgejsv.f"
		temp1 = sqrt((doublereal) (*m)) * epsln;
#line 1702 "cgejsv.f"
		i__1 = n1;
#line 1702 "cgejsv.f"
		for (p = 1; p <= i__1; ++p) {
#line 1703 "cgejsv.f"
		    xsc = 1. / scnrm2_(m, &u[p * u_dim1 + 1], &c__1);
#line 1704 "cgejsv.f"
		    if (xsc < 1. - temp1 || xsc > temp1 + 1.) {
#line 1704 "cgejsv.f"
			csscal_(m, &xsc, &u[p * u_dim1 + 1], &c__1);
#line 1704 "cgejsv.f"
		    }
#line 1706 "cgejsv.f"
/* L6973: */
#line 1706 "cgejsv.f"
		}

#line 1708 "cgejsv.f"
		if (rowpiv) {
#line 1708 "cgejsv.f"
		    i__1 = *m - 1;
#line 1708 "cgejsv.f"
		    claswp_(&n1, &u[u_offset], ldu, &c__1, &i__1, &iwork[(*n 
			    << 1) + 1], &c_n1);
#line 1708 "cgejsv.f"
		}

#line 1711 "cgejsv.f"
	    }

/*        end of the  >> almost orthogonal case <<  in the full SVD */

#line 1715 "cgejsv.f"
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

#line 1727 "cgejsv.f"
	    i__1 = nr;
#line 1727 "cgejsv.f"
	    for (p = 1; p <= i__1; ++p) {
#line 1728 "cgejsv.f"
		i__2 = *n - p + 1;
#line 1728 "cgejsv.f"
		ccopy_(&i__2, &a[p + p * a_dim1], lda, &v[p + p * v_dim1], &
			c__1);
#line 1729 "cgejsv.f"
		i__2 = *n - p + 1;
#line 1729 "cgejsv.f"
		clacgv_(&i__2, &v[p + p * v_dim1], &c__1);
#line 1730 "cgejsv.f"
/* L7968: */
#line 1730 "cgejsv.f"
	    }

#line 1732 "cgejsv.f"
	    if (l2pert) {
#line 1733 "cgejsv.f"
		xsc = sqrt(small / epsln);
#line 1734 "cgejsv.f"
		i__1 = nr;
#line 1734 "cgejsv.f"
		for (q = 1; q <= i__1; ++q) {
#line 1735 "cgejsv.f"
		    d__1 = xsc * z_abs(&v[q + q * v_dim1]);
#line 1735 "cgejsv.f"
		    z__1.r = d__1, z__1.i = 0.;
#line 1735 "cgejsv.f"
		    ctemp.r = z__1.r, ctemp.i = z__1.i;
#line 1736 "cgejsv.f"
		    i__2 = *n;
#line 1736 "cgejsv.f"
		    for (p = 1; p <= i__2; ++p) {
#line 1737 "cgejsv.f"
			if (p > q && z_abs(&v[p + q * v_dim1]) <= temp1 || p <
				 q) {
#line 1737 "cgejsv.f"
			    i__3 = p + q * v_dim1;
#line 1737 "cgejsv.f"
			    v[i__3].r = ctemp.r, v[i__3].i = ctemp.i;
#line 1737 "cgejsv.f"
			}
/*     $                V(p,q) = TEMP1 * ( V(p,q) / ABS(V(p,q)) ) */
#line 1741 "cgejsv.f"
			if (p < q) {
#line 1741 "cgejsv.f"
			    i__3 = p + q * v_dim1;
#line 1741 "cgejsv.f"
			    i__4 = p + q * v_dim1;
#line 1741 "cgejsv.f"
			    z__1.r = -v[i__4].r, z__1.i = -v[i__4].i;
#line 1741 "cgejsv.f"
			    v[i__3].r = z__1.r, v[i__3].i = z__1.i;
#line 1741 "cgejsv.f"
			}
#line 1742 "cgejsv.f"
/* L5968: */
#line 1742 "cgejsv.f"
		    }
#line 1743 "cgejsv.f"
/* L5969: */
#line 1743 "cgejsv.f"
		}
#line 1744 "cgejsv.f"
	    } else {
#line 1745 "cgejsv.f"
		i__1 = nr - 1;
#line 1745 "cgejsv.f"
		i__2 = nr - 1;
#line 1745 "cgejsv.f"
		claset_("U", &i__1, &i__2, &c_b1, &c_b1, &v[(v_dim1 << 1) + 1]
			, ldv, (ftnlen)1);
#line 1746 "cgejsv.f"
	    }
#line 1748 "cgejsv.f"
	    i__1 = *lwork - (*n << 1);
#line 1748 "cgejsv.f"
	    cgeqrf_(n, &nr, &v[v_offset], ldv, &cwork[*n + 1], &cwork[(*n << 
		    1) + 1], &i__1, &ierr);
#line 1750 "cgejsv.f"
	    clacpy_("L", n, &nr, &v[v_offset], ldv, &cwork[(*n << 1) + 1], n, 
		    (ftnlen)1);

#line 1752 "cgejsv.f"
	    i__1 = nr;
#line 1752 "cgejsv.f"
	    for (p = 1; p <= i__1; ++p) {
#line 1753 "cgejsv.f"
		i__2 = nr - p + 1;
#line 1753 "cgejsv.f"
		ccopy_(&i__2, &v[p + p * v_dim1], ldv, &u[p + p * u_dim1], &
			c__1);
#line 1754 "cgejsv.f"
		i__2 = nr - p + 1;
#line 1754 "cgejsv.f"
		clacgv_(&i__2, &u[p + p * u_dim1], &c__1);
#line 1755 "cgejsv.f"
/* L7969: */
#line 1755 "cgejsv.f"
	    }
#line 1757 "cgejsv.f"
	    if (l2pert) {
#line 1758 "cgejsv.f"
		xsc = sqrt(small / epsln);
#line 1759 "cgejsv.f"
		i__1 = nr;
#line 1759 "cgejsv.f"
		for (q = 2; q <= i__1; ++q) {
#line 1760 "cgejsv.f"
		    i__2 = q - 1;
#line 1760 "cgejsv.f"
		    for (p = 1; p <= i__2; ++p) {
/* Computing MIN */
#line 1761 "cgejsv.f"
			d__2 = z_abs(&u[p + p * u_dim1]), d__3 = z_abs(&u[q + 
				q * u_dim1]);
#line 1761 "cgejsv.f"
			d__1 = xsc * min(d__2,d__3);
#line 1761 "cgejsv.f"
			z__1.r = d__1, z__1.i = 0.;
#line 1761 "cgejsv.f"
			ctemp.r = z__1.r, ctemp.i = z__1.i;
/*                  U(p,q) = - TEMP1 * ( U(q,p) / ABS(U(q,p)) ) */
#line 1764 "cgejsv.f"
			i__3 = p + q * u_dim1;
#line 1764 "cgejsv.f"
			z__1.r = -ctemp.r, z__1.i = -ctemp.i;
#line 1764 "cgejsv.f"
			u[i__3].r = z__1.r, u[i__3].i = z__1.i;
#line 1765 "cgejsv.f"
/* L9971: */
#line 1765 "cgejsv.f"
		    }
#line 1766 "cgejsv.f"
/* L9970: */
#line 1766 "cgejsv.f"
		}
#line 1767 "cgejsv.f"
	    } else {
#line 1768 "cgejsv.f"
		i__1 = nr - 1;
#line 1768 "cgejsv.f"
		i__2 = nr - 1;
#line 1768 "cgejsv.f"
		claset_("U", &i__1, &i__2, &c_b1, &c_b1, &u[(u_dim1 << 1) + 1]
			, ldu, (ftnlen)1);
#line 1769 "cgejsv.f"
	    }
#line 1771 "cgejsv.f"
	    i__1 = *lwork - (*n << 1) - *n * nr;
#line 1771 "cgejsv.f"
	    cgesvj_("L", "U", "V", &nr, &nr, &u[u_offset], ldu, &sva[1], n, &
		    v[v_offset], ldv, &cwork[(*n << 1) + *n * nr + 1], &i__1, 
		    &rwork[1], lrwork, info, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1774 "cgejsv.f"
	    scalem = rwork[1];
#line 1775 "cgejsv.f"
	    numrank = i_dnnt(&rwork[2]);
#line 1777 "cgejsv.f"
	    if (nr < *n) {
#line 1778 "cgejsv.f"
		i__1 = *n - nr;
#line 1778 "cgejsv.f"
		claset_("A", &i__1, &nr, &c_b120, &c_b120, &v[nr + 1 + v_dim1]
			, ldv, (ftnlen)1);
#line 1779 "cgejsv.f"
		i__1 = *n - nr;
#line 1779 "cgejsv.f"
		claset_("A", &nr, &i__1, &c_b120, &c_b120, &v[(nr + 1) * 
			v_dim1 + 1], ldv, (ftnlen)1);
#line 1780 "cgejsv.f"
		i__1 = *n - nr;
#line 1780 "cgejsv.f"
		i__2 = *n - nr;
#line 1780 "cgejsv.f"
		claset_("A", &i__1, &i__2, &c_b120, &c_b80, &v[nr + 1 + (nr + 
			1) * v_dim1], ldv, (ftnlen)1);
#line 1781 "cgejsv.f"
	    }
#line 1783 "cgejsv.f"
	    i__1 = *lwork - (*n << 1) - *n * nr - nr;
#line 1783 "cgejsv.f"
	    cunmqr_("L", "N", n, n, &nr, &cwork[(*n << 1) + 1], n, &cwork[*n 
		    + 1], &v[v_offset], ldv, &cwork[(*n << 1) + *n * nr + nr 
		    + 1], &i__1, &ierr, (ftnlen)1, (ftnlen)1);

/*           Permute the rows of V using the (column) permutation from the */
/*           first QRF. Also, scale the columns to make them unit in */
/*           Euclidean norm. This applies to all cases. */

#line 1790 "cgejsv.f"
	    temp1 = sqrt((doublereal) (*n)) * epsln;
#line 1791 "cgejsv.f"
	    i__1 = *n;
#line 1791 "cgejsv.f"
	    for (q = 1; q <= i__1; ++q) {
#line 1792 "cgejsv.f"
		i__2 = *n;
#line 1792 "cgejsv.f"
		for (p = 1; p <= i__2; ++p) {
#line 1793 "cgejsv.f"
		    i__3 = (*n << 1) + *n * nr + nr + iwork[p];
#line 1793 "cgejsv.f"
		    i__4 = p + q * v_dim1;
#line 1793 "cgejsv.f"
		    cwork[i__3].r = v[i__4].r, cwork[i__3].i = v[i__4].i;
#line 1794 "cgejsv.f"
/* L8972: */
#line 1794 "cgejsv.f"
		}
#line 1795 "cgejsv.f"
		i__2 = *n;
#line 1795 "cgejsv.f"
		for (p = 1; p <= i__2; ++p) {
#line 1796 "cgejsv.f"
		    i__3 = p + q * v_dim1;
#line 1796 "cgejsv.f"
		    i__4 = (*n << 1) + *n * nr + nr + p;
#line 1796 "cgejsv.f"
		    v[i__3].r = cwork[i__4].r, v[i__3].i = cwork[i__4].i;
#line 1797 "cgejsv.f"
/* L8973: */
#line 1797 "cgejsv.f"
		}
#line 1798 "cgejsv.f"
		xsc = 1. / scnrm2_(n, &v[q * v_dim1 + 1], &c__1);
#line 1799 "cgejsv.f"
		if (xsc < 1. - temp1 || xsc > temp1 + 1.) {
#line 1799 "cgejsv.f"
		    csscal_(n, &xsc, &v[q * v_dim1 + 1], &c__1);
#line 1799 "cgejsv.f"
		}
#line 1801 "cgejsv.f"
/* L7972: */
#line 1801 "cgejsv.f"
	    }

/*           At this moment, V contains the right singular vectors of A. */
/*           Next, assemble the left singular vector matrix U (M x N). */

#line 1806 "cgejsv.f"
	    if (nr < *m) {
#line 1807 "cgejsv.f"
		i__1 = *m - nr;
#line 1807 "cgejsv.f"
		claset_("A", &i__1, &nr, &c_b1, &c_b1, &u[nr + 1 + u_dim1], 
			ldu, (ftnlen)1);
#line 1808 "cgejsv.f"
		if (nr < n1) {
#line 1809 "cgejsv.f"
		    i__1 = n1 - nr;
#line 1809 "cgejsv.f"
		    claset_("A", &nr, &i__1, &c_b1, &c_b1, &u[(nr + 1) * 
			    u_dim1 + 1], ldu, (ftnlen)1);
#line 1810 "cgejsv.f"
		    i__1 = *m - nr;
#line 1810 "cgejsv.f"
		    i__2 = n1 - nr;
#line 1810 "cgejsv.f"
		    claset_("A", &i__1, &i__2, &c_b1, &c_b2, &u[nr + 1 + (nr 
			    + 1) * u_dim1], ldu, (ftnlen)1);
#line 1811 "cgejsv.f"
		}
#line 1812 "cgejsv.f"
	    }

#line 1814 "cgejsv.f"
	    i__1 = *lwork - *n;
#line 1814 "cgejsv.f"
	    cunmqr_("Left", "No Tr", m, &n1, n, &a[a_offset], lda, &cwork[1], 
		    &u[u_offset], ldu, &cwork[*n + 1], &i__1, &ierr, (ftnlen)
		    4, (ftnlen)5);

#line 1817 "cgejsv.f"
	    if (rowpiv) {
#line 1817 "cgejsv.f"
		i__1 = *m - 1;
#line 1817 "cgejsv.f"
		claswp_(&n1, &u[u_offset], ldu, &c__1, &i__1, &iwork[(*n << 1)
			 + 1], &c_n1);
#line 1817 "cgejsv.f"
	    }


#line 1821 "cgejsv.f"
	}
#line 1822 "cgejsv.f"
	if (transp) {
/*           .. swap U and V because the procedure worked on A^* */
#line 1824 "cgejsv.f"
	    i__1 = *n;
#line 1824 "cgejsv.f"
	    for (p = 1; p <= i__1; ++p) {
#line 1825 "cgejsv.f"
		cswap_(n, &u[p * u_dim1 + 1], &c__1, &v[p * v_dim1 + 1], &
			c__1);
#line 1826 "cgejsv.f"
/* L6974: */
#line 1826 "cgejsv.f"
	    }
#line 1827 "cgejsv.f"
	}

#line 1829 "cgejsv.f"
    }
/*     end of the full SVD */

/*     Undo scaling, if necessary (and possible) */

#line 1834 "cgejsv.f"
    if (uscal2 <= big / sva[1] * uscal1) {
#line 1835 "cgejsv.f"
	clascl_("G", &c__0, &c__0, &uscal1, &uscal2, &nr, &c__1, &sva[1], n, &
		ierr, (ftnlen)1);
#line 1836 "cgejsv.f"
	uscal1 = 1.;
#line 1837 "cgejsv.f"
	uscal2 = 1.;
#line 1838 "cgejsv.f"
    }

#line 1840 "cgejsv.f"
    if (nr < *n) {
#line 1841 "cgejsv.f"
	i__1 = *n;
#line 1841 "cgejsv.f"
	for (p = nr + 1; p <= i__1; ++p) {
#line 1842 "cgejsv.f"
	    sva[p] = 0.;
#line 1843 "cgejsv.f"
/* L3004: */
#line 1843 "cgejsv.f"
	}
#line 1844 "cgejsv.f"
    }

#line 1846 "cgejsv.f"
    rwork[1] = uscal2 * scalem;
#line 1847 "cgejsv.f"
    rwork[2] = uscal1;
#line 1848 "cgejsv.f"
    if (errest) {
#line 1848 "cgejsv.f"
	rwork[3] = sconda;
#line 1848 "cgejsv.f"
    }
#line 1849 "cgejsv.f"
    if (lsvec && rsvec) {
#line 1850 "cgejsv.f"
	rwork[4] = condr1;
#line 1851 "cgejsv.f"
	rwork[5] = condr2;
#line 1852 "cgejsv.f"
    }
#line 1853 "cgejsv.f"
    if (l2tran) {
#line 1854 "cgejsv.f"
	rwork[6] = entra;
#line 1855 "cgejsv.f"
	rwork[7] = entrat;
#line 1856 "cgejsv.f"
    }

#line 1858 "cgejsv.f"
    iwork[1] = nr;
#line 1859 "cgejsv.f"
    iwork[2] = numrank;
#line 1860 "cgejsv.f"
    iwork[3] = warning;

#line 1862 "cgejsv.f"
    return 0;
/*     .. */
/*     .. END OF CGEJSV */
/*     .. */
} /* cgejsv_ */

