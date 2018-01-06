#line 1 "sgejsv.f"
/* sgejsv.f -- translated by f2c (version 20100827).
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

#line 1 "sgejsv.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b34 = 0.;
static doublereal c_b35 = 1.;
static integer c__0 = 0;
static integer c_n1 = -1;

/* > \brief \b SGEJSV */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SGEJSV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgejsv.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgejsv.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgejsv.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SGEJSV( JOBA, JOBU, JOBV, JOBR, JOBT, JOBP, */
/*                          M, N, A, LDA, SVA, U, LDU, V, LDV, */
/*                          WORK, LWORK, IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       IMPLICIT    NONE */
/*       INTEGER     INFO, LDA, LDU, LDV, LWORK, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL        A( LDA, * ), SVA( N ), U( LDU, * ), V( LDV, * ), */
/*      $            WORK( LWORK ) */
/*       INTEGER     IWORK( * ) */
/*       CHARACTER*1 JOBA, JOBP, JOBR, JOBT, JOBU, JOBV */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SGEJSV computes the singular value decomposition (SVD) of a real M-by-N */
/* > matrix [A], where M >= N. The SVD of [A] is written as */
/* > */
/* >              [A] = [U] * [SIGMA] * [V]^t, */
/* > */
/* > where [SIGMA] is an N-by-N (M-by-N) matrix which is zero except for its N */
/* > diagonal elements, [U] is an M-by-N (or M-by-M) orthonormal matrix, and */
/* > [V] is an N-by-N orthogonal matrix. The diagonal elements of [SIGMA] are */
/* > the singular values of [A]. The columns of [U] and [V] are the left and */
/* > the right singular vectors of [A], respectively. The matrices [U] and [V] */
/* > are computed and stored in the arrays U and V, respectively. The diagonal */
/* > of [SIGMA] is computed and stored in the array SVA. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

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
/* >              The computed SVD A = U * S * V^t restores A up to */
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
/* >              is greater than BIG, use SGESVJ. */
/* >       = 'R': RESTRICTED range for sigma(c*A) is [SQRT(SFMIN), SQRT(BIG)] */
/* >              (roughly, as described above). This option is recommended. */
/* >                                             =========================== */
/* >         For computing the singular values in the FULL range [SFMIN,BIG] */
/* >         use SGESVJ. */
/* > \endverbatim */
/* > */
/* > \param[in] JOBT */
/* > \verbatim */
/* >          JOBT is CHARACTER*1 */
/* >         If the matrix is square then the procedure may determine to use */
/* >         transposed A if A^t seems to be better with respect to convergence. */
/* >         If the matrix is not square, JOBT is ignored. This is subject to */
/* >         changes in the future. */
/* >         The decision is based on two values of entropy over the adjoint */
/* >         orbit of A^t * A. See the descriptions of WORK(6) and WORK(7). */
/* >       = 'T': transpose if entropy test indicates possibly faster */
/* >         convergence of Jacobi process if A^t is taken as input. If A is */
/* >         replaced with A^t, then the row pivoting is included automatically. */
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
/* >          A is REAL array, dimension (LDA,N) */
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
/* >          U is REAL array, dimension ( LDU, N ) */
/* >          If JOBU = 'U', then U contains on exit the M-by-N matrix of */
/* >                         the left singular vectors. */
/* >          If JOBU = 'F', then U contains on exit the M-by-M matrix of */
/* >                         the left singular vectors, including an ONB */
/* >                         of the orthogonal complement of the Range(A). */
/* >          If JOBU = 'W'  .AND. (JOBV.EQ.'V' .AND. JOBT.EQ.'T' .AND. M.EQ.N), */
/* >                         then U is used as workspace if the procedure */
/* >                         replaces A with A^t. In that case, [V] is computed */
/* >                         in U as left singular vectors of A^t and then */
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
/* >          V is REAL array, dimension ( LDV, N ) */
/* >          If JOBV = 'V', 'J' then V contains on exit the N-by-N matrix of */
/* >                         the right singular vectors; */
/* >          If JOBV = 'W', AND (JOBU.EQ.'U' AND JOBT.EQ.'T' AND M.EQ.N), */
/* >                         then V is used as workspace if the pprocedure */
/* >                         replaces A with A^t. In that case, [U] is computed */
/* >                         in V as right singular vectors of A^t and then */
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
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is REAL array, dimension at least LWORK. */
/* >          On exit, */
/* >          WORK(1) = SCALE = WORK(2) / WORK(1) is the scaling factor such */
/* >                    that SCALE*SVA(1:N) are the computed singular values */
/* >                    of A. (See the description of SVA().) */
/* >          WORK(2) = See the description of WORK(1). */
/* >          WORK(3) = SCONDA is an estimate for the condition number of */
/* >                    column equilibrated A. (If JOBA .EQ. 'E' or 'G') */
/* >                    SCONDA is an estimate of SQRT(||(R^t * R)^(-1)||_1). */
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
/* >          WORK(4) = an estimate of the scaled condition number of the */
/* >                    triangular factor in the first QR factorization. */
/* >          WORK(5) = an estimate of the scaled condition number of the */
/* >                    triangular factor in the second QR factorization. */
/* >          The following two parameters are computed if JOBT .EQ. 'T'. */
/* >          They are provided for a developer/implementer who is familiar */
/* >          with the details of the method. */
/* > */
/* >          WORK(6) = the entropy of A^t*A :: this is the Shannon entropy */
/* >                    of diag(A^t*A) / Trace(A^t*A) taken as point in the */
/* >                    probability simplex. */
/* >          WORK(7) = the entropy of A*A^t. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          Length of WORK to confirm proper allocation of work space. */
/* >          LWORK depends on the job: */
/* > */
/* >          If only SIGMA is needed ( JOBU.EQ.'N', JOBV.EQ.'N' ) and */
/* >            -> .. no scaled condition estimate required (JOBE.EQ.'N'): */
/* >               LWORK >= max(2*M+N,4*N+1,7). This is the minimal requirement. */
/* >               ->> For optimal performance (blocked code) the optimal value */
/* >               is LWORK >= max(2*M+N,3*N+(N+1)*NB,7). Here NB is the optimal */
/* >               block size for DGEQP3 and DGEQRF. */
/* >               In general, optimal LWORK is computed as */
/* >               LWORK >= max(2*M+N,N+LWORK(DGEQP3),N+LWORK(DGEQRF), 7). */
/* >            -> .. an estimate of the scaled condition number of A is */
/* >               required (JOBA='E', 'G'). In this case, LWORK is the maximum */
/* >               of the above and N*N+4*N, i.e. LWORK >= max(2*M+N,N*N+4*N,7). */
/* >               ->> For optimal performance (blocked code) the optimal value */
/* >               is LWORK >= max(2*M+N,3*N+(N+1)*NB, N*N+4*N, 7). */
/* >               In general, the optimal length LWORK is computed as */
/* >               LWORK >= max(2*M+N,N+LWORK(DGEQP3),N+LWORK(DGEQRF), */
/* >                                                     N+N*N+LWORK(DPOCON),7). */
/* > */
/* >          If SIGMA and the right singular vectors are needed (JOBV.EQ.'V'), */
/* >            -> the minimal requirement is LWORK >= max(2*M+N,4*N+1,7). */
/* >            -> For optimal performance, LWORK >= max(2*M+N,3*N+(N+1)*NB,7), */
/* >               where NB is the optimal block size for DGEQP3, DGEQRF, DGELQ, */
/* >               DORMLQ. In general, the optimal length LWORK is computed as */
/* >               LWORK >= max(2*M+N,N+LWORK(DGEQP3), N+LWORK(DPOCON), */
/* >                       N+LWORK(DGELQ), 2*N+LWORK(DGEQRF), N+LWORK(DORMLQ)). */
/* > */
/* >          If SIGMA and the left singular vectors are needed */
/* >            -> the minimal requirement is LWORK >= max(2*M+N,4*N+1,7). */
/* >            -> For optimal performance: */
/* >               if JOBU.EQ.'U' :: LWORK >= max(2*M+N,3*N+(N+1)*NB,7), */
/* >               if JOBU.EQ.'F' :: LWORK >= max(2*M+N,3*N+(N+1)*NB,N+M*NB,7), */
/* >               where NB is the optimal block size for DGEQP3, DGEQRF, DORMQR. */
/* >               In general, the optimal length LWORK is computed as */
/* >               LWORK >= max(2*M+N,N+LWORK(DGEQP3),N+LWORK(DPOCON), */
/* >                        2*N+LWORK(DGEQRF), N+LWORK(DORMQR)). */
/* >               Here LWORK(DORMQR) equals N*NB (for JOBU.EQ.'U') or */
/* >               M*NB (for JOBU.EQ.'F'). */
/* > */
/* >          If the full SVD is needed: (JOBU.EQ.'U' or JOBU.EQ.'F') and */
/* >            -> if JOBV.EQ.'V' */
/* >               the minimal requirement is LWORK >= max(2*M+N,6*N+2*N*N). */
/* >            -> if JOBV.EQ.'J' the minimal requirement is */
/* >               LWORK >= max(2*M+N, 4*N+N*N,2*N+N*N+6). */
/* >            -> For optimal performance, LWORK should be additionally */
/* >               larger than N+M*NB, where NB is the optimal block size */
/* >               for DORMQR. */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* >          IWORK is INTEGER array, dimension M+3*N. */
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
/* >           > 0 :  SGEJSV  did not converge in the maximal allowed number */
/* >                  of sweeps. The computed values may be inaccurate. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup realGEsing */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  SGEJSV implements a preconditioned Jacobi SVD algorithm. It uses SGEQP3, */
/* >  SGEQRF, and SGELQF as preprocessors and preconditioners. Optionally, an */
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
/* >  & LAPACK routines called by SGEJSV are implemented to work in that range. */
/* >  If that is not the case, then the restriction for safe computation with */
/* >  the singular values in the range of normalized IEEE numbers is that the */
/* >  spectral condition number kappa(A)=sigma_max(A)/sigma_min(A) does not */
/* >  overflow. This code (SGEJSV) is best used in this restricted range, */
/* >  meaning that singular values of magnitude below ||A||_2 / SLAMCH('O') are */
/* >  returned as zeros. See JOBR for details on this. */
/* >     Further, this implementation is somewhat slower than the one described */
/* >  in [1,2] due to replacement of some non-LAPACK components, and because */
/* >  the choice of some tuning parameters in the iterative part (SGESVJ) is */
/* >  left to the implementer on a particular machine. */
/* >     The rank revealing QR factorization (in this code: SGEQP3) should be */
/* >  implemented as in [3]. We have a new version of SGEQP3 under development */
/* >  that is more robust than the current one in LAPACK, with a cleaner cut in */
/* >  rank defficient cases. It will be available in the SIGMA library [4]. */
/* >  If M is much larger than N, it is obvious that the inital QRF with */
/* >  column pivoting can be preprocessed by the QRF without pivoting. That */
/* >  well known trick is not used in SGEJSV because in some cases heavy row */
/* >  weighting can be treated with complete pivoting. The overhead in cases */
/* >  M much larger than N is then only due to pivoting, but the benefits in */
/* >  terms of accuracy have prevailed. The implementer/user can incorporate */
/* >  this extra QRF step easily. The implementer can also improve data movement */
/* >  (matrix transpose, matrix copy, matrix transposed copy) - this */
/* >  implementation of SGEJSV uses only the simplest, naive data movement. */
/* > \endverbatim */

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
/* Subroutine */ int sgejsv_(char *joba, char *jobu, char *jobv, char *jobr, 
	char *jobt, char *jobp, integer *m, integer *n, doublereal *a, 
	integer *lda, doublereal *sva, doublereal *u, integer *ldu, 
	doublereal *v, integer *ldv, doublereal *work, integer *lwork, 
	integer *iwork, integer *info, ftnlen joba_len, ftnlen jobu_len, 
	ftnlen jobv_len, ftnlen jobr_len, ftnlen jobt_len, ftnlen jobp_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, u_dim1, u_offset, v_dim1, v_offset, i__1, i__2, 
	    i__3, i__4, i__5, i__6, i__7, i__8, i__9, i__10, i__11, i__12;
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double sqrt(doublereal), log(doublereal), d_sign(doublereal *, doublereal 
	    *);
    integer i_dnnt(doublereal *);

    /* Local variables */
    static integer p, q, n1, nr;
    static doublereal big, xsc, big1;
    static logical defr;
    static doublereal aapp, aaqq;
    static logical kill;
    static integer ierr;
    static doublereal temp1;
    extern doublereal snrm2_(integer *, doublereal *, integer *);
    static logical jracc;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal small, entra, sfmin;
    static logical lsvec;
    static doublereal epsln;
    static logical rsvec;
    extern /* Subroutine */ int scopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), sswap_(integer *, doublereal *, integer 
	    *, doublereal *, integer *);
    static logical l2aber;
    extern /* Subroutine */ int strsm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static doublereal condr1, condr2, uscal1, uscal2;
    static logical l2kill, l2rank, l2tran;
    extern /* Subroutine */ int sgeqp3_(integer *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    integer *);
    static logical l2pert;
    static doublereal scalem, sconda;
    static logical goscal;
    static doublereal aatmin;
    extern doublereal slamch_(char *, ftnlen);
    static doublereal aatmax;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static logical noscal;
    extern /* Subroutine */ int sgelqf_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *);
    extern integer isamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int slascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen), sgeqrf_(integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *), slacpy_(char *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, ftnlen), slaset_(char *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, ftnlen);
    static doublereal entrat;
    static logical almort;
    static doublereal maxprj;
    extern /* Subroutine */ int spocon_(char *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    integer *, ftnlen);
    static logical errest;
    extern /* Subroutine */ int sgesvj_(char *, char *, char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen, ftnlen), slassq_(integer *, doublereal *, integer 
	    *, doublereal *, doublereal *);
    static logical transp;
    extern /* Subroutine */ int slaswp_(integer *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *), sorgqr_(integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, integer *), sormlq_(char *, char *, 
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen), sormqr_(char *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen);
    static logical rowpiv;
    static doublereal cond_ok__;
    static integer warning, numrank;


/*  -- LAPACK computational routine (version 3.4.2) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     September 2012 */

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

#line 528 "sgejsv.f"
    /* Parameter adjustments */
#line 528 "sgejsv.f"
    --sva;
#line 528 "sgejsv.f"
    a_dim1 = *lda;
#line 528 "sgejsv.f"
    a_offset = 1 + a_dim1;
#line 528 "sgejsv.f"
    a -= a_offset;
#line 528 "sgejsv.f"
    u_dim1 = *ldu;
#line 528 "sgejsv.f"
    u_offset = 1 + u_dim1;
#line 528 "sgejsv.f"
    u -= u_offset;
#line 528 "sgejsv.f"
    v_dim1 = *ldv;
#line 528 "sgejsv.f"
    v_offset = 1 + v_dim1;
#line 528 "sgejsv.f"
    v -= v_offset;
#line 528 "sgejsv.f"
    --work;
#line 528 "sgejsv.f"
    --iwork;
#line 528 "sgejsv.f"

#line 528 "sgejsv.f"
    /* Function Body */
#line 528 "sgejsv.f"
    lsvec = lsame_(jobu, "U", (ftnlen)1, (ftnlen)1) || lsame_(jobu, "F", (
	    ftnlen)1, (ftnlen)1);
#line 529 "sgejsv.f"
    jracc = lsame_(jobv, "J", (ftnlen)1, (ftnlen)1);
#line 530 "sgejsv.f"
    rsvec = lsame_(jobv, "V", (ftnlen)1, (ftnlen)1) || jracc;
#line 531 "sgejsv.f"
    rowpiv = lsame_(joba, "F", (ftnlen)1, (ftnlen)1) || lsame_(joba, "G", (
	    ftnlen)1, (ftnlen)1);
#line 532 "sgejsv.f"
    l2rank = lsame_(joba, "R", (ftnlen)1, (ftnlen)1);
#line 533 "sgejsv.f"
    l2aber = lsame_(joba, "A", (ftnlen)1, (ftnlen)1);
#line 534 "sgejsv.f"
    errest = lsame_(joba, "E", (ftnlen)1, (ftnlen)1) || lsame_(joba, "G", (
	    ftnlen)1, (ftnlen)1);
#line 535 "sgejsv.f"
    l2tran = lsame_(jobt, "T", (ftnlen)1, (ftnlen)1);
#line 536 "sgejsv.f"
    l2kill = lsame_(jobr, "R", (ftnlen)1, (ftnlen)1);
#line 537 "sgejsv.f"
    defr = lsame_(jobr, "N", (ftnlen)1, (ftnlen)1);
#line 538 "sgejsv.f"
    l2pert = lsame_(jobp, "P", (ftnlen)1, (ftnlen)1);

#line 540 "sgejsv.f"
    if (! (rowpiv || l2rank || l2aber || errest || lsame_(joba, "C", (ftnlen)
	    1, (ftnlen)1))) {
#line 542 "sgejsv.f"
	*info = -1;
#line 543 "sgejsv.f"
    } else if (! (lsvec || lsame_(jobu, "N", (ftnlen)1, (ftnlen)1) || lsame_(
	    jobu, "W", (ftnlen)1, (ftnlen)1))) {
#line 545 "sgejsv.f"
	*info = -2;
#line 546 "sgejsv.f"
    } else if (! (rsvec || lsame_(jobv, "N", (ftnlen)1, (ftnlen)1) || lsame_(
	    jobv, "W", (ftnlen)1, (ftnlen)1)) || jracc && ! lsvec) {
#line 548 "sgejsv.f"
	*info = -3;
#line 549 "sgejsv.f"
    } else if (! (l2kill || defr)) {
#line 550 "sgejsv.f"
	*info = -4;
#line 551 "sgejsv.f"
    } else if (! (l2tran || lsame_(jobt, "N", (ftnlen)1, (ftnlen)1))) {
#line 552 "sgejsv.f"
	*info = -5;
#line 553 "sgejsv.f"
    } else if (! (l2pert || lsame_(jobp, "N", (ftnlen)1, (ftnlen)1))) {
#line 554 "sgejsv.f"
	*info = -6;
#line 555 "sgejsv.f"
    } else if (*m < 0) {
#line 556 "sgejsv.f"
	*info = -7;
#line 557 "sgejsv.f"
    } else if (*n < 0 || *n > *m) {
#line 558 "sgejsv.f"
	*info = -8;
#line 559 "sgejsv.f"
    } else if (*lda < *m) {
#line 560 "sgejsv.f"
	*info = -10;
#line 561 "sgejsv.f"
    } else if (lsvec && *ldu < *m) {
#line 562 "sgejsv.f"
	*info = -13;
#line 563 "sgejsv.f"
    } else if (rsvec && *ldv < *n) {
#line 564 "sgejsv.f"
	*info = -14;
#line 565 "sgejsv.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 565 "sgejsv.f"
	i__1 = 7, i__2 = (*n << 2) + 1, i__1 = max(i__1,i__2), i__2 = (*m << 
		1) + *n;
/* Computing MAX */
#line 565 "sgejsv.f"
	i__3 = 7, i__4 = (*n << 2) + *n * *n, i__3 = max(i__3,i__4), i__4 = (*
		m << 1) + *n;
/* Computing MAX */
#line 565 "sgejsv.f"
	i__5 = 7, i__6 = (*m << 1) + *n, i__5 = max(i__5,i__6), i__6 = (*n << 
		2) + 1;
/* Computing MAX */
#line 565 "sgejsv.f"
	i__7 = 7, i__8 = (*m << 1) + *n, i__7 = max(i__7,i__8), i__8 = (*n << 
		2) + 1;
/* Computing MAX */
#line 565 "sgejsv.f"
	i__9 = (*m << 1) + *n, i__10 = *n * 6 + (*n << 1) * *n;
/* Computing MAX */
#line 565 "sgejsv.f"
	i__11 = (*m << 1) + *n, i__12 = (*n << 2) + *n * *n, i__11 = max(
		i__11,i__12), i__12 = (*n << 1) + *n * *n + 6;
#line 565 "sgejsv.f"
	if (! (lsvec || rsvec || errest) && *lwork < max(i__1,i__2) || ! (
		lsvec || rsvec) && errest && *lwork < max(i__3,i__4) || lsvec 
		&& ! rsvec && *lwork < max(i__5,i__6) || rsvec && ! lsvec && *
		lwork < max(i__7,i__8) || lsvec && rsvec && ! jracc && *lwork 
		< max(i__9,i__10) || lsvec && rsvec && jracc && *lwork < max(
		i__11,i__12)) {
#line 578 "sgejsv.f"
	    *info = -17;
#line 579 "sgejsv.f"
	} else {
/*        #:) */
#line 581 "sgejsv.f"
	    *info = 0;
#line 582 "sgejsv.f"
	}
#line 582 "sgejsv.f"
    }

#line 584 "sgejsv.f"
    if (*info != 0) {
/*       #:( */
#line 586 "sgejsv.f"
	i__1 = -(*info);
#line 586 "sgejsv.f"
	xerbla_("SGEJSV", &i__1, (ftnlen)6);
#line 587 "sgejsv.f"
	return 0;
#line 588 "sgejsv.f"
    }

/*     Quick return for void matrix (Y3K safe) */
/* #:) */
#line 592 "sgejsv.f"
    if (*m == 0 || *n == 0) {
#line 592 "sgejsv.f"
	return 0;
#line 592 "sgejsv.f"
    }

/*     Determine whether the matrix U should be M x N or M x M */

#line 596 "sgejsv.f"
    if (lsvec) {
#line 597 "sgejsv.f"
	n1 = *n;
#line 598 "sgejsv.f"
	if (lsame_(jobu, "F", (ftnlen)1, (ftnlen)1)) {
#line 598 "sgejsv.f"
	    n1 = *m;
#line 598 "sgejsv.f"
	}
#line 599 "sgejsv.f"
    }

/*     Set numerical parameters */

/* !    NOTE: Make sure SLAMCH() does not fail on the target architecture. */

#line 605 "sgejsv.f"
    epsln = slamch_("Epsilon", (ftnlen)7);
#line 606 "sgejsv.f"
    sfmin = slamch_("SafeMinimum", (ftnlen)11);
#line 607 "sgejsv.f"
    small = sfmin / epsln;
#line 608 "sgejsv.f"
    big = slamch_("O", (ftnlen)1);
/*     BIG   = ONE / SFMIN */

/*     Initialize SVA(1:N) = diag( ||A e_i||_2 )_1^N */

/* (!)  If necessary, scale SVA() to protect the largest norm from */
/*     overflow. It is possible that this scaling pushes the smallest */
/*     column norm left from the underflow threshold (extreme case). */

#line 617 "sgejsv.f"
    scalem = 1. / sqrt((doublereal) (*m) * (doublereal) (*n));
#line 618 "sgejsv.f"
    noscal = TRUE_;
#line 619 "sgejsv.f"
    goscal = TRUE_;
#line 620 "sgejsv.f"
    i__1 = *n;
#line 620 "sgejsv.f"
    for (p = 1; p <= i__1; ++p) {
#line 621 "sgejsv.f"
	aapp = 0.;
#line 622 "sgejsv.f"
	aaqq = 1.;
#line 623 "sgejsv.f"
	slassq_(m, &a[p * a_dim1 + 1], &c__1, &aapp, &aaqq);
#line 624 "sgejsv.f"
	if (aapp > big) {
#line 625 "sgejsv.f"
	    *info = -9;
#line 626 "sgejsv.f"
	    i__2 = -(*info);
#line 626 "sgejsv.f"
	    xerbla_("SGEJSV", &i__2, (ftnlen)6);
#line 627 "sgejsv.f"
	    return 0;
#line 628 "sgejsv.f"
	}
#line 629 "sgejsv.f"
	aaqq = sqrt(aaqq);
#line 630 "sgejsv.f"
	if (aapp < big / aaqq && noscal) {
#line 631 "sgejsv.f"
	    sva[p] = aapp * aaqq;
#line 632 "sgejsv.f"
	} else {
#line 633 "sgejsv.f"
	    noscal = FALSE_;
#line 634 "sgejsv.f"
	    sva[p] = aapp * (aaqq * scalem);
#line 635 "sgejsv.f"
	    if (goscal) {
#line 636 "sgejsv.f"
		goscal = FALSE_;
#line 637 "sgejsv.f"
		i__2 = p - 1;
#line 637 "sgejsv.f"
		sscal_(&i__2, &scalem, &sva[1], &c__1);
#line 638 "sgejsv.f"
	    }
#line 639 "sgejsv.f"
	}
#line 640 "sgejsv.f"
/* L1874: */
#line 640 "sgejsv.f"
    }

#line 642 "sgejsv.f"
    if (noscal) {
#line 642 "sgejsv.f"
	scalem = 1.;
#line 642 "sgejsv.f"
    }

#line 644 "sgejsv.f"
    aapp = 0.;
#line 645 "sgejsv.f"
    aaqq = big;
#line 646 "sgejsv.f"
    i__1 = *n;
#line 646 "sgejsv.f"
    for (p = 1; p <= i__1; ++p) {
/* Computing MAX */
#line 647 "sgejsv.f"
	d__1 = aapp, d__2 = sva[p];
#line 647 "sgejsv.f"
	aapp = max(d__1,d__2);
#line 648 "sgejsv.f"
	if (sva[p] != 0.) {
/* Computing MIN */
#line 648 "sgejsv.f"
	    d__1 = aaqq, d__2 = sva[p];
#line 648 "sgejsv.f"
	    aaqq = min(d__1,d__2);
#line 648 "sgejsv.f"
	}
#line 649 "sgejsv.f"
/* L4781: */
#line 649 "sgejsv.f"
    }

/*     Quick return for zero M x N matrix */
/* #:) */
#line 653 "sgejsv.f"
    if (aapp == 0.) {
#line 654 "sgejsv.f"
	if (lsvec) {
#line 654 "sgejsv.f"
	    slaset_("G", m, &n1, &c_b34, &c_b35, &u[u_offset], ldu, (ftnlen)1)
		    ;
#line 654 "sgejsv.f"
	}
#line 655 "sgejsv.f"
	if (rsvec) {
#line 655 "sgejsv.f"
	    slaset_("G", n, n, &c_b34, &c_b35, &v[v_offset], ldv, (ftnlen)1);
#line 655 "sgejsv.f"
	}
#line 656 "sgejsv.f"
	work[1] = 1.;
#line 657 "sgejsv.f"
	work[2] = 1.;
#line 658 "sgejsv.f"
	if (errest) {
#line 658 "sgejsv.f"
	    work[3] = 1.;
#line 658 "sgejsv.f"
	}
#line 659 "sgejsv.f"
	if (lsvec && rsvec) {
#line 660 "sgejsv.f"
	    work[4] = 1.;
#line 661 "sgejsv.f"
	    work[5] = 1.;
#line 662 "sgejsv.f"
	}
#line 663 "sgejsv.f"
	if (l2tran) {
#line 664 "sgejsv.f"
	    work[6] = 0.;
#line 665 "sgejsv.f"
	    work[7] = 0.;
#line 666 "sgejsv.f"
	}
#line 667 "sgejsv.f"
	iwork[1] = 0;
#line 668 "sgejsv.f"
	iwork[2] = 0;
#line 669 "sgejsv.f"
	iwork[3] = 0;
#line 670 "sgejsv.f"
	return 0;
#line 671 "sgejsv.f"
    }

/*     Issue warning if denormalized column norms detected. Override the */
/*     high relative accuracy request. Issue licence to kill columns */
/*     (set them to zero) whose norm is less than sigma_max / BIG (roughly). */
/* #:( */
#line 677 "sgejsv.f"
    warning = 0;
#line 678 "sgejsv.f"
    if (aaqq <= sfmin) {
#line 679 "sgejsv.f"
	l2rank = TRUE_;
#line 680 "sgejsv.f"
	l2kill = TRUE_;
#line 681 "sgejsv.f"
	warning = 1;
#line 682 "sgejsv.f"
    }

/*     Quick return for one-column matrix */
/* #:) */
#line 686 "sgejsv.f"
    if (*n == 1) {

#line 688 "sgejsv.f"
	if (lsvec) {
#line 689 "sgejsv.f"
	    slascl_("G", &c__0, &c__0, &sva[1], &scalem, m, &c__1, &a[a_dim1 
		    + 1], lda, &ierr, (ftnlen)1);
#line 690 "sgejsv.f"
	    slacpy_("A", m, &c__1, &a[a_offset], lda, &u[u_offset], ldu, (
		    ftnlen)1);
/*           computing all M left singular vectors of the M x 1 matrix */
#line 692 "sgejsv.f"
	    if (n1 != *n) {
#line 693 "sgejsv.f"
		i__1 = *lwork - *n;
#line 693 "sgejsv.f"
		sgeqrf_(m, n, &u[u_offset], ldu, &work[1], &work[*n + 1], &
			i__1, &ierr);
#line 694 "sgejsv.f"
		i__1 = *lwork - *n;
#line 694 "sgejsv.f"
		sorgqr_(m, &n1, &c__1, &u[u_offset], ldu, &work[1], &work[*n 
			+ 1], &i__1, &ierr);
#line 695 "sgejsv.f"
		scopy_(m, &a[a_dim1 + 1], &c__1, &u[u_dim1 + 1], &c__1);
#line 696 "sgejsv.f"
	    }
#line 697 "sgejsv.f"
	}
#line 698 "sgejsv.f"
	if (rsvec) {
#line 699 "sgejsv.f"
	    v[v_dim1 + 1] = 1.;
#line 700 "sgejsv.f"
	}
#line 701 "sgejsv.f"
	if (sva[1] < big * scalem) {
#line 702 "sgejsv.f"
	    sva[1] /= scalem;
#line 703 "sgejsv.f"
	    scalem = 1.;
#line 704 "sgejsv.f"
	}
#line 705 "sgejsv.f"
	work[1] = 1. / scalem;
#line 706 "sgejsv.f"
	work[2] = 1.;
#line 707 "sgejsv.f"
	if (sva[1] != 0.) {
#line 708 "sgejsv.f"
	    iwork[1] = 1;
#line 709 "sgejsv.f"
	    if (sva[1] / scalem >= sfmin) {
#line 710 "sgejsv.f"
		iwork[2] = 1;
#line 711 "sgejsv.f"
	    } else {
#line 712 "sgejsv.f"
		iwork[2] = 0;
#line 713 "sgejsv.f"
	    }
#line 714 "sgejsv.f"
	} else {
#line 715 "sgejsv.f"
	    iwork[1] = 0;
#line 716 "sgejsv.f"
	    iwork[2] = 0;
#line 717 "sgejsv.f"
	}
#line 718 "sgejsv.f"
	if (errest) {
#line 718 "sgejsv.f"
	    work[3] = 1.;
#line 718 "sgejsv.f"
	}
#line 719 "sgejsv.f"
	if (lsvec && rsvec) {
#line 720 "sgejsv.f"
	    work[4] = 1.;
#line 721 "sgejsv.f"
	    work[5] = 1.;
#line 722 "sgejsv.f"
	}
#line 723 "sgejsv.f"
	if (l2tran) {
#line 724 "sgejsv.f"
	    work[6] = 0.;
#line 725 "sgejsv.f"
	    work[7] = 0.;
#line 726 "sgejsv.f"
	}
#line 727 "sgejsv.f"
	return 0;

#line 729 "sgejsv.f"
    }

#line 731 "sgejsv.f"
    transp = FALSE_;
#line 732 "sgejsv.f"
    l2tran = l2tran && *m == *n;

#line 734 "sgejsv.f"
    aatmax = -1.;
#line 735 "sgejsv.f"
    aatmin = big;
#line 736 "sgejsv.f"
    if (rowpiv || l2tran) {

/*     Compute the row norms, needed to determine row pivoting sequence */
/*     (in the case of heavily row weighted A, row pivoting is strongly */
/*     advised) and to collect information needed to compare the */
/*     structures of A * A^t and A^t * A (in the case L2TRAN.EQ..TRUE.). */

#line 743 "sgejsv.f"
	if (l2tran) {
#line 744 "sgejsv.f"
	    i__1 = *m;
#line 744 "sgejsv.f"
	    for (p = 1; p <= i__1; ++p) {
#line 745 "sgejsv.f"
		xsc = 0.;
#line 746 "sgejsv.f"
		temp1 = 1.;
#line 747 "sgejsv.f"
		slassq_(n, &a[p + a_dim1], lda, &xsc, &temp1);
/*              SLASSQ gets both the ell_2 and the ell_infinity norm */
/*              in one pass through the vector */
#line 750 "sgejsv.f"
		work[*m + *n + p] = xsc * scalem;
#line 751 "sgejsv.f"
		work[*n + p] = xsc * (scalem * sqrt(temp1));
/* Computing MAX */
#line 752 "sgejsv.f"
		d__1 = aatmax, d__2 = work[*n + p];
#line 752 "sgejsv.f"
		aatmax = max(d__1,d__2);
#line 753 "sgejsv.f"
		if (work[*n + p] != 0.) {
/* Computing MIN */
#line 753 "sgejsv.f"
		    d__1 = aatmin, d__2 = work[*n + p];
#line 753 "sgejsv.f"
		    aatmin = min(d__1,d__2);
#line 753 "sgejsv.f"
		}
#line 754 "sgejsv.f"
/* L1950: */
#line 754 "sgejsv.f"
	    }
#line 755 "sgejsv.f"
	} else {
#line 756 "sgejsv.f"
	    i__1 = *m;
#line 756 "sgejsv.f"
	    for (p = 1; p <= i__1; ++p) {
#line 757 "sgejsv.f"
		work[*m + *n + p] = scalem * (d__1 = a[p + isamax_(n, &a[p + 
			a_dim1], lda) * a_dim1], abs(d__1));
/* Computing MAX */
#line 758 "sgejsv.f"
		d__1 = aatmax, d__2 = work[*m + *n + p];
#line 758 "sgejsv.f"
		aatmax = max(d__1,d__2);
/* Computing MIN */
#line 759 "sgejsv.f"
		d__1 = aatmin, d__2 = work[*m + *n + p];
#line 759 "sgejsv.f"
		aatmin = min(d__1,d__2);
#line 760 "sgejsv.f"
/* L1904: */
#line 760 "sgejsv.f"
	    }
#line 761 "sgejsv.f"
	}

#line 763 "sgejsv.f"
    }

/*     For square matrix A try to determine whether A^t  would be  better */
/*     input for the preconditioned Jacobi SVD, with faster convergence. */
/*     The decision is based on an O(N) function of the vector of column */
/*     and row norms of A, based on the Shannon entropy. This should give */
/*     the right choice in most cases when the difference actually matters. */
/*     It may fail and pick the slower converging side. */

#line 772 "sgejsv.f"
    entra = 0.;
#line 773 "sgejsv.f"
    entrat = 0.;
#line 774 "sgejsv.f"
    if (l2tran) {

#line 776 "sgejsv.f"
	xsc = 0.;
#line 777 "sgejsv.f"
	temp1 = 1.;
#line 778 "sgejsv.f"
	slassq_(n, &sva[1], &c__1, &xsc, &temp1);
#line 779 "sgejsv.f"
	temp1 = 1. / temp1;

#line 781 "sgejsv.f"
	entra = 0.;
#line 782 "sgejsv.f"
	i__1 = *n;
#line 782 "sgejsv.f"
	for (p = 1; p <= i__1; ++p) {
/* Computing 2nd power */
#line 783 "sgejsv.f"
	    d__1 = sva[p] / xsc;
#line 783 "sgejsv.f"
	    big1 = d__1 * d__1 * temp1;
#line 784 "sgejsv.f"
	    if (big1 != 0.) {
#line 784 "sgejsv.f"
		entra += big1 * log(big1);
#line 784 "sgejsv.f"
	    }
#line 785 "sgejsv.f"
/* L1113: */
#line 785 "sgejsv.f"
	}
#line 786 "sgejsv.f"
	entra = -entra / log((doublereal) (*n));

/*        Now, SVA().^2/Trace(A^t * A) is a point in the probability simplex. */
/*        It is derived from the diagonal of  A^t * A.  Do the same with the */
/*        diagonal of A * A^t, compute the entropy of the corresponding */
/*        probability distribution. Note that A * A^t and A^t * A have the */
/*        same trace. */

#line 794 "sgejsv.f"
	entrat = 0.;
#line 795 "sgejsv.f"
	i__1 = *n + *m;
#line 795 "sgejsv.f"
	for (p = *n + 1; p <= i__1; ++p) {
/* Computing 2nd power */
#line 796 "sgejsv.f"
	    d__1 = work[p] / xsc;
#line 796 "sgejsv.f"
	    big1 = d__1 * d__1 * temp1;
#line 797 "sgejsv.f"
	    if (big1 != 0.) {
#line 797 "sgejsv.f"
		entrat += big1 * log(big1);
#line 797 "sgejsv.f"
	    }
#line 798 "sgejsv.f"
/* L1114: */
#line 798 "sgejsv.f"
	}
#line 799 "sgejsv.f"
	entrat = -entrat / log((doublereal) (*m));

/*        Analyze the entropies and decide A or A^t. Smaller entropy */
/*        usually means better input for the algorithm. */

#line 804 "sgejsv.f"
	transp = entrat < entra;

/*        If A^t is better than A, transpose A. */

#line 808 "sgejsv.f"
	if (transp) {
/*           In an optimal implementation, this trivial transpose */
/*           should be replaced with faster transpose. */
#line 811 "sgejsv.f"
	    i__1 = *n - 1;
#line 811 "sgejsv.f"
	    for (p = 1; p <= i__1; ++p) {
#line 812 "sgejsv.f"
		i__2 = *n;
#line 812 "sgejsv.f"
		for (q = p + 1; q <= i__2; ++q) {
#line 813 "sgejsv.f"
		    temp1 = a[q + p * a_dim1];
#line 814 "sgejsv.f"
		    a[q + p * a_dim1] = a[p + q * a_dim1];
#line 815 "sgejsv.f"
		    a[p + q * a_dim1] = temp1;
#line 816 "sgejsv.f"
/* L1116: */
#line 816 "sgejsv.f"
		}
#line 817 "sgejsv.f"
/* L1115: */
#line 817 "sgejsv.f"
	    }
#line 818 "sgejsv.f"
	    i__1 = *n;
#line 818 "sgejsv.f"
	    for (p = 1; p <= i__1; ++p) {
#line 819 "sgejsv.f"
		work[*m + *n + p] = sva[p];
#line 820 "sgejsv.f"
		sva[p] = work[*n + p];
#line 821 "sgejsv.f"
/* L1117: */
#line 821 "sgejsv.f"
	    }
#line 822 "sgejsv.f"
	    temp1 = aapp;
#line 823 "sgejsv.f"
	    aapp = aatmax;
#line 824 "sgejsv.f"
	    aatmax = temp1;
#line 825 "sgejsv.f"
	    temp1 = aaqq;
#line 826 "sgejsv.f"
	    aaqq = aatmin;
#line 827 "sgejsv.f"
	    aatmin = temp1;
#line 828 "sgejsv.f"
	    kill = lsvec;
#line 829 "sgejsv.f"
	    lsvec = rsvec;
#line 830 "sgejsv.f"
	    rsvec = kill;
#line 831 "sgejsv.f"
	    if (lsvec) {
#line 831 "sgejsv.f"
		n1 = *n;
#line 831 "sgejsv.f"
	    }

#line 833 "sgejsv.f"
	    rowpiv = TRUE_;
#line 834 "sgejsv.f"
	}

#line 836 "sgejsv.f"
    }
/*     END IF L2TRAN */

/*     Scale the matrix so that its maximal singular value remains less */
/*     than SQRT(BIG) -- the matrix is scaled so that its maximal column */
/*     has Euclidean norm equal to SQRT(BIG/N). The only reason to keep */
/*     SQRT(BIG) instead of BIG is the fact that SGEJSV uses LAPACK and */
/*     BLAS routines that, in some implementations, are not capable of */
/*     working in the full interval [SFMIN,BIG] and that they may provoke */
/*     overflows in the intermediate results. If the singular values spread */
/*     from SFMIN to BIG, then SGESVJ will compute them. So, in that case, */
/*     one should use SGESVJ instead of SGEJSV. */

#line 849 "sgejsv.f"
    big1 = sqrt(big);
#line 850 "sgejsv.f"
    temp1 = sqrt(big / (doublereal) (*n));

#line 852 "sgejsv.f"
    slascl_("G", &c__0, &c__0, &aapp, &temp1, n, &c__1, &sva[1], n, &ierr, (
	    ftnlen)1);
#line 853 "sgejsv.f"
    if (aaqq > aapp * sfmin) {
#line 854 "sgejsv.f"
	aaqq = aaqq / aapp * temp1;
#line 855 "sgejsv.f"
    } else {
#line 856 "sgejsv.f"
	aaqq = aaqq * temp1 / aapp;
#line 857 "sgejsv.f"
    }
#line 858 "sgejsv.f"
    temp1 *= scalem;
#line 859 "sgejsv.f"
    slascl_("G", &c__0, &c__0, &aapp, &temp1, m, n, &a[a_offset], lda, &ierr, 
	    (ftnlen)1);

/*     To undo scaling at the end of this procedure, multiply the */
/*     computed singular values with USCAL2 / USCAL1. */

#line 864 "sgejsv.f"
    uscal1 = temp1;
#line 865 "sgejsv.f"
    uscal2 = aapp;

#line 867 "sgejsv.f"
    if (l2kill) {
/*        L2KILL enforces computation of nonzero singular values in */
/*        the restricted range of condition number of the initial A, */
/*        sigma_max(A) / sigma_min(A) approx. SQRT(BIG)/SQRT(SFMIN). */
#line 871 "sgejsv.f"
	xsc = sqrt(sfmin);
#line 872 "sgejsv.f"
    } else {
#line 873 "sgejsv.f"
	xsc = small;

/*        Now, if the condition number of A is too big, */
/*        sigma_max(A) / sigma_min(A) .GT. SQRT(BIG/N) * EPSLN / SFMIN, */
/*        as a precaution measure, the full SVD is computed using SGESVJ */
/*        with accumulated Jacobi rotations. This provides numerically */
/*        more robust computation, at the cost of slightly increased run */
/*        time. Depending on the concrete implementation of BLAS and LAPACK */
/*        (i.e. how they behave in presence of extreme ill-conditioning) the */
/*        implementor may decide to remove this switch. */
#line 883 "sgejsv.f"
	if (aaqq < sqrt(sfmin) && lsvec && rsvec) {
#line 884 "sgejsv.f"
	    jracc = TRUE_;
#line 885 "sgejsv.f"
	}

#line 887 "sgejsv.f"
    }
#line 888 "sgejsv.f"
    if (aaqq < xsc) {
#line 889 "sgejsv.f"
	i__1 = *n;
#line 889 "sgejsv.f"
	for (p = 1; p <= i__1; ++p) {
#line 890 "sgejsv.f"
	    if (sva[p] < xsc) {
#line 891 "sgejsv.f"
		slaset_("A", m, &c__1, &c_b34, &c_b34, &a[p * a_dim1 + 1], 
			lda, (ftnlen)1);
#line 892 "sgejsv.f"
		sva[p] = 0.;
#line 893 "sgejsv.f"
	    }
#line 894 "sgejsv.f"
/* L700: */
#line 894 "sgejsv.f"
	}
#line 895 "sgejsv.f"
    }

/*     Preconditioning using QR factorization with pivoting */

#line 899 "sgejsv.f"
    if (rowpiv) {
/*        Optional row permutation (Bjoerck row pivoting): */
/*        A result by Cox and Higham shows that the Bjoerck's */
/*        row pivoting combined with standard column pivoting */
/*        has similar effect as Powell-Reid complete pivoting. */
/*        The ell-infinity norms of A are made nonincreasing. */
#line 905 "sgejsv.f"
	i__1 = *m - 1;
#line 905 "sgejsv.f"
	for (p = 1; p <= i__1; ++p) {
#line 906 "sgejsv.f"
	    i__2 = *m - p + 1;
#line 906 "sgejsv.f"
	    q = isamax_(&i__2, &work[*m + *n + p], &c__1) + p - 1;
#line 907 "sgejsv.f"
	    iwork[(*n << 1) + p] = q;
#line 908 "sgejsv.f"
	    if (p != q) {
#line 909 "sgejsv.f"
		temp1 = work[*m + *n + p];
#line 910 "sgejsv.f"
		work[*m + *n + p] = work[*m + *n + q];
#line 911 "sgejsv.f"
		work[*m + *n + q] = temp1;
#line 912 "sgejsv.f"
	    }
#line 913 "sgejsv.f"
/* L1952: */
#line 913 "sgejsv.f"
	}
#line 914 "sgejsv.f"
	i__1 = *m - 1;
#line 914 "sgejsv.f"
	slaswp_(n, &a[a_offset], lda, &c__1, &i__1, &iwork[(*n << 1) + 1], &
		c__1);
#line 915 "sgejsv.f"
    }

/*     End of the preparation phase (scaling, optional sorting and */
/*     transposing, optional flushing of small columns). */

/*     Preconditioning */

/*     If the full SVD is needed, the right singular vectors are computed */
/*     from a matrix equation, and for that we need theoretical analysis */
/*     of the Businger-Golub pivoting. So we use SGEQP3 as the first RR QRF. */
/*     In all other cases the first RR QRF can be chosen by other criteria */
/*     (eg speed by replacing global with restricted window pivoting, such */
/*     as in SGEQPX from TOMS # 782). Good results will be obtained using */
/*     SGEQPX with properly (!) chosen numerical parameters. */
/*     Any improvement of SGEQP3 improves overal performance of SGEJSV. */

/*     A * P1 = Q1 * [ R1^t 0]^t: */
#line 932 "sgejsv.f"
    i__1 = *n;
#line 932 "sgejsv.f"
    for (p = 1; p <= i__1; ++p) {
/*        .. all columns are free columns */
#line 934 "sgejsv.f"
	iwork[p] = 0;
#line 935 "sgejsv.f"
/* L1963: */
#line 935 "sgejsv.f"
    }
#line 936 "sgejsv.f"
    i__1 = *lwork - *n;
#line 936 "sgejsv.f"
    sgeqp3_(m, n, &a[a_offset], lda, &iwork[1], &work[1], &work[*n + 1], &
	    i__1, &ierr);

/*     The upper triangular matrix R1 from the first QRF is inspected for */
/*     rank deficiency and possibilities for deflation, or possible */
/*     ill-conditioning. Depending on the user specified flag L2RANK, */
/*     the procedure explores possibilities to reduce the numerical */
/*     rank by inspecting the computed upper triangular factor. If */
/*     L2RANK or L2ABER are up, then SGEJSV will compute the SVD of */
/*     A + dA, where ||dA|| <= f(M,N)*EPSLN. */

#line 946 "sgejsv.f"
    nr = 1;
#line 947 "sgejsv.f"
    if (l2aber) {
/*        Standard absolute error bound suffices. All sigma_i with */
/*        sigma_i < N*EPSLN*||A|| are flushed to zero. This is an */
/*        agressive enforcement of lower numerical rank by introducing a */
/*        backward error of the order of N*EPSLN*||A||. */
#line 952 "sgejsv.f"
	temp1 = sqrt((doublereal) (*n)) * epsln;
#line 953 "sgejsv.f"
	i__1 = *n;
#line 953 "sgejsv.f"
	for (p = 2; p <= i__1; ++p) {
#line 954 "sgejsv.f"
	    if ((d__2 = a[p + p * a_dim1], abs(d__2)) >= temp1 * (d__1 = a[
		    a_dim1 + 1], abs(d__1))) {
#line 955 "sgejsv.f"
		++nr;
#line 956 "sgejsv.f"
	    } else {
#line 957 "sgejsv.f"
		goto L3002;
#line 958 "sgejsv.f"
	    }
#line 959 "sgejsv.f"
/* L3001: */
#line 959 "sgejsv.f"
	}
#line 960 "sgejsv.f"
L3002:
#line 961 "sgejsv.f"
	;
#line 961 "sgejsv.f"
    } else if (l2rank) {
/*        .. similarly as above, only slightly more gentle (less agressive). */
/*        Sudden drop on the diagonal of R1 is used as the criterion for */
/*        close-to-rank-defficient. */
#line 965 "sgejsv.f"
	temp1 = sqrt(sfmin);
#line 966 "sgejsv.f"
	i__1 = *n;
#line 966 "sgejsv.f"
	for (p = 2; p <= i__1; ++p) {
#line 967 "sgejsv.f"
	    if ((d__2 = a[p + p * a_dim1], abs(d__2)) < epsln * (d__1 = a[p - 
		    1 + (p - 1) * a_dim1], abs(d__1)) || (d__3 = a[p + p * 
		    a_dim1], abs(d__3)) < small || l2kill && (d__4 = a[p + p *
		     a_dim1], abs(d__4)) < temp1) {
#line 967 "sgejsv.f"
		goto L3402;
#line 967 "sgejsv.f"
	    }
#line 970 "sgejsv.f"
	    ++nr;
#line 971 "sgejsv.f"
/* L3401: */
#line 971 "sgejsv.f"
	}
#line 972 "sgejsv.f"
L3402:

#line 974 "sgejsv.f"
	;
#line 974 "sgejsv.f"
    } else {
/*        The goal is high relative accuracy. However, if the matrix */
/*        has high scaled condition number the relative accuracy is in */
/*        general not feasible. Later on, a condition number estimator */
/*        will be deployed to estimate the scaled condition number. */
/*        Here we just remove the underflowed part of the triangular */
/*        factor. This prevents the situation in which the code is */
/*        working hard to get the accuracy not warranted by the data. */
#line 982 "sgejsv.f"
	temp1 = sqrt(sfmin);
#line 983 "sgejsv.f"
	i__1 = *n;
#line 983 "sgejsv.f"
	for (p = 2; p <= i__1; ++p) {
#line 984 "sgejsv.f"
	    if ((d__1 = a[p + p * a_dim1], abs(d__1)) < small || l2kill && (
		    d__2 = a[p + p * a_dim1], abs(d__2)) < temp1) {
#line 984 "sgejsv.f"
		goto L3302;
#line 984 "sgejsv.f"
	    }
#line 986 "sgejsv.f"
	    ++nr;
#line 987 "sgejsv.f"
/* L3301: */
#line 987 "sgejsv.f"
	}
#line 988 "sgejsv.f"
L3302:

#line 990 "sgejsv.f"
	;
#line 990 "sgejsv.f"
    }

#line 992 "sgejsv.f"
    almort = FALSE_;
#line 993 "sgejsv.f"
    if (nr == *n) {
#line 994 "sgejsv.f"
	maxprj = 1.;
#line 995 "sgejsv.f"
	i__1 = *n;
#line 995 "sgejsv.f"
	for (p = 2; p <= i__1; ++p) {
#line 996 "sgejsv.f"
	    temp1 = (d__1 = a[p + p * a_dim1], abs(d__1)) / sva[iwork[p]];
#line 997 "sgejsv.f"
	    maxprj = min(maxprj,temp1);
#line 998 "sgejsv.f"
/* L3051: */
#line 998 "sgejsv.f"
	}
/* Computing 2nd power */
#line 999 "sgejsv.f"
	d__1 = maxprj;
#line 999 "sgejsv.f"
	if (d__1 * d__1 >= 1. - (doublereal) (*n) * epsln) {
#line 999 "sgejsv.f"
	    almort = TRUE_;
#line 999 "sgejsv.f"
	}
#line 1000 "sgejsv.f"
    }


#line 1003 "sgejsv.f"
    sconda = -1.;
#line 1004 "sgejsv.f"
    condr1 = -1.;
#line 1005 "sgejsv.f"
    condr2 = -1.;

#line 1007 "sgejsv.f"
    if (errest) {
#line 1008 "sgejsv.f"
	if (*n == nr) {
#line 1009 "sgejsv.f"
	    if (rsvec) {
/*              .. V is available as workspace */
#line 1011 "sgejsv.f"
		slacpy_("U", n, n, &a[a_offset], lda, &v[v_offset], ldv, (
			ftnlen)1);
#line 1012 "sgejsv.f"
		i__1 = *n;
#line 1012 "sgejsv.f"
		for (p = 1; p <= i__1; ++p) {
#line 1013 "sgejsv.f"
		    temp1 = sva[iwork[p]];
#line 1014 "sgejsv.f"
		    d__1 = 1. / temp1;
#line 1014 "sgejsv.f"
		    sscal_(&p, &d__1, &v[p * v_dim1 + 1], &c__1);
#line 1015 "sgejsv.f"
/* L3053: */
#line 1015 "sgejsv.f"
		}
#line 1016 "sgejsv.f"
		spocon_("U", n, &v[v_offset], ldv, &c_b35, &temp1, &work[*n + 
			1], &iwork[(*n << 1) + *m + 1], &ierr, (ftnlen)1);
#line 1018 "sgejsv.f"
	    } else if (lsvec) {
/*              .. U is available as workspace */
#line 1020 "sgejsv.f"
		slacpy_("U", n, n, &a[a_offset], lda, &u[u_offset], ldu, (
			ftnlen)1);
#line 1021 "sgejsv.f"
		i__1 = *n;
#line 1021 "sgejsv.f"
		for (p = 1; p <= i__1; ++p) {
#line 1022 "sgejsv.f"
		    temp1 = sva[iwork[p]];
#line 1023 "sgejsv.f"
		    d__1 = 1. / temp1;
#line 1023 "sgejsv.f"
		    sscal_(&p, &d__1, &u[p * u_dim1 + 1], &c__1);
#line 1024 "sgejsv.f"
/* L3054: */
#line 1024 "sgejsv.f"
		}
#line 1025 "sgejsv.f"
		spocon_("U", n, &u[u_offset], ldu, &c_b35, &temp1, &work[*n + 
			1], &iwork[(*n << 1) + *m + 1], &ierr, (ftnlen)1);
#line 1027 "sgejsv.f"
	    } else {
#line 1028 "sgejsv.f"
		slacpy_("U", n, n, &a[a_offset], lda, &work[*n + 1], n, (
			ftnlen)1);
#line 1029 "sgejsv.f"
		i__1 = *n;
#line 1029 "sgejsv.f"
		for (p = 1; p <= i__1; ++p) {
#line 1030 "sgejsv.f"
		    temp1 = sva[iwork[p]];
#line 1031 "sgejsv.f"
		    d__1 = 1. / temp1;
#line 1031 "sgejsv.f"
		    sscal_(&p, &d__1, &work[*n + (p - 1) * *n + 1], &c__1);
#line 1032 "sgejsv.f"
/* L3052: */
#line 1032 "sgejsv.f"
		}
/*           .. the columns of R are scaled to have unit Euclidean lengths. */
#line 1034 "sgejsv.f"
		spocon_("U", n, &work[*n + 1], n, &c_b35, &temp1, &work[*n + *
			n * *n + 1], &iwork[(*n << 1) + *m + 1], &ierr, (
			ftnlen)1);
#line 1036 "sgejsv.f"
	    }
#line 1037 "sgejsv.f"
	    sconda = 1. / sqrt(temp1);
/*           SCONDA is an estimate of SQRT(||(R^t * R)^(-1)||_1). */
/*           N^(-1/4) * SCONDA <= ||R^(-1)||_2 <= N^(1/4) * SCONDA */
#line 1040 "sgejsv.f"
	} else {
#line 1041 "sgejsv.f"
	    sconda = -1.;
#line 1042 "sgejsv.f"
	}
#line 1043 "sgejsv.f"
    }

#line 1045 "sgejsv.f"
    l2pert = l2pert && (d__1 = a[a_dim1 + 1] / a[nr + nr * a_dim1], abs(d__1))
	     > sqrt(big1);
/*     If there is no violent scaling, artificial perturbation is not needed. */

/*     Phase 3: */

#line 1050 "sgejsv.f"
    if (! (rsvec || lsvec)) {

/*         Singular Values only */

/*         .. transpose A(1:NR,1:N) */
/* Computing MIN */
#line 1055 "sgejsv.f"
	i__2 = *n - 1;
#line 1055 "sgejsv.f"
	i__1 = min(i__2,nr);
#line 1055 "sgejsv.f"
	for (p = 1; p <= i__1; ++p) {
#line 1056 "sgejsv.f"
	    i__2 = *n - p;
#line 1056 "sgejsv.f"
	    scopy_(&i__2, &a[p + (p + 1) * a_dim1], lda, &a[p + 1 + p * 
		    a_dim1], &c__1);
#line 1057 "sgejsv.f"
/* L1946: */
#line 1057 "sgejsv.f"
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

#line 1071 "sgejsv.f"
	if (! almort) {

#line 1073 "sgejsv.f"
	    if (l2pert) {
/*              XSC = SQRT(SMALL) */
#line 1075 "sgejsv.f"
		xsc = epsln / (doublereal) (*n);
#line 1076 "sgejsv.f"
		i__1 = nr;
#line 1076 "sgejsv.f"
		for (q = 1; q <= i__1; ++q) {
#line 1077 "sgejsv.f"
		    temp1 = xsc * (d__1 = a[q + q * a_dim1], abs(d__1));
#line 1078 "sgejsv.f"
		    i__2 = *n;
#line 1078 "sgejsv.f"
		    for (p = 1; p <= i__2; ++p) {
#line 1079 "sgejsv.f"
			if (p > q && (d__1 = a[p + q * a_dim1], abs(d__1)) <= 
				temp1 || p < q) {
#line 1079 "sgejsv.f"
			    a[p + q * a_dim1] = d_sign(&temp1, &a[p + q * 
				    a_dim1]);
#line 1079 "sgejsv.f"
			}
#line 1082 "sgejsv.f"
/* L4949: */
#line 1082 "sgejsv.f"
		    }
#line 1083 "sgejsv.f"
/* L4947: */
#line 1083 "sgejsv.f"
		}
#line 1084 "sgejsv.f"
	    } else {
#line 1085 "sgejsv.f"
		i__1 = nr - 1;
#line 1085 "sgejsv.f"
		i__2 = nr - 1;
#line 1085 "sgejsv.f"
		slaset_("U", &i__1, &i__2, &c_b34, &c_b34, &a[(a_dim1 << 1) + 
			1], lda, (ftnlen)1);
#line 1086 "sgejsv.f"
	    }

/*            .. second preconditioning using the QR factorization */

#line 1090 "sgejsv.f"
	    i__1 = *lwork - *n;
#line 1090 "sgejsv.f"
	    sgeqrf_(n, &nr, &a[a_offset], lda, &work[1], &work[*n + 1], &i__1,
		     &ierr);

/*           .. and transpose upper to lower triangular */
#line 1093 "sgejsv.f"
	    i__1 = nr - 1;
#line 1093 "sgejsv.f"
	    for (p = 1; p <= i__1; ++p) {
#line 1094 "sgejsv.f"
		i__2 = nr - p;
#line 1094 "sgejsv.f"
		scopy_(&i__2, &a[p + (p + 1) * a_dim1], lda, &a[p + 1 + p * 
			a_dim1], &c__1);
#line 1095 "sgejsv.f"
/* L1948: */
#line 1095 "sgejsv.f"
	    }

#line 1097 "sgejsv.f"
	}

/*           Row-cyclic Jacobi SVD algorithm with column pivoting */

/*           .. again some perturbation (a "background noise") is added */
/*           to drown denormals */
#line 1103 "sgejsv.f"
	if (l2pert) {
/*              XSC = SQRT(SMALL) */
#line 1105 "sgejsv.f"
	    xsc = epsln / (doublereal) (*n);
#line 1106 "sgejsv.f"
	    i__1 = nr;
#line 1106 "sgejsv.f"
	    for (q = 1; q <= i__1; ++q) {
#line 1107 "sgejsv.f"
		temp1 = xsc * (d__1 = a[q + q * a_dim1], abs(d__1));
#line 1108 "sgejsv.f"
		i__2 = nr;
#line 1108 "sgejsv.f"
		for (p = 1; p <= i__2; ++p) {
#line 1109 "sgejsv.f"
		    if (p > q && (d__1 = a[p + q * a_dim1], abs(d__1)) <= 
			    temp1 || p < q) {
#line 1109 "sgejsv.f"
			a[p + q * a_dim1] = d_sign(&temp1, &a[p + q * a_dim1])
				;
#line 1109 "sgejsv.f"
		    }
#line 1112 "sgejsv.f"
/* L1949: */
#line 1112 "sgejsv.f"
		}
#line 1113 "sgejsv.f"
/* L1947: */
#line 1113 "sgejsv.f"
	    }
#line 1114 "sgejsv.f"
	} else {
#line 1115 "sgejsv.f"
	    i__1 = nr - 1;
#line 1115 "sgejsv.f"
	    i__2 = nr - 1;
#line 1115 "sgejsv.f"
	    slaset_("U", &i__1, &i__2, &c_b34, &c_b34, &a[(a_dim1 << 1) + 1], 
		    lda, (ftnlen)1);
#line 1116 "sgejsv.f"
	}

/*           .. and one-sided Jacobi rotations are started on a lower */
/*           triangular matrix (plus perturbation which is ignored in */
/*           the part which destroys triangular form (confusing?!)) */

#line 1122 "sgejsv.f"
	sgesvj_("L", "NoU", "NoV", &nr, &nr, &a[a_offset], lda, &sva[1], n, &
		v[v_offset], ldv, &work[1], lwork, info, (ftnlen)1, (ftnlen)3,
		 (ftnlen)3);

#line 1125 "sgejsv.f"
	scalem = work[1];
#line 1126 "sgejsv.f"
	numrank = i_dnnt(&work[2]);


#line 1129 "sgejsv.f"
    } else if (rsvec && ! lsvec) {

/*        -> Singular Values and Right Singular Vectors <- */

#line 1133 "sgejsv.f"
	if (almort) {

/*           .. in this case NR equals N */
#line 1136 "sgejsv.f"
	    i__1 = nr;
#line 1136 "sgejsv.f"
	    for (p = 1; p <= i__1; ++p) {
#line 1137 "sgejsv.f"
		i__2 = *n - p + 1;
#line 1137 "sgejsv.f"
		scopy_(&i__2, &a[p + p * a_dim1], lda, &v[p + p * v_dim1], &
			c__1);
#line 1138 "sgejsv.f"
/* L1998: */
#line 1138 "sgejsv.f"
	    }
#line 1139 "sgejsv.f"
	    i__1 = nr - 1;
#line 1139 "sgejsv.f"
	    i__2 = nr - 1;
#line 1139 "sgejsv.f"
	    slaset_("Upper", &i__1, &i__2, &c_b34, &c_b34, &v[(v_dim1 << 1) + 
		    1], ldv, (ftnlen)5);

#line 1141 "sgejsv.f"
	    sgesvj_("L", "U", "N", n, &nr, &v[v_offset], ldv, &sva[1], &nr, &
		    a[a_offset], lda, &work[1], lwork, info, (ftnlen)1, (
		    ftnlen)1, (ftnlen)1);
#line 1143 "sgejsv.f"
	    scalem = work[1];
#line 1144 "sgejsv.f"
	    numrank = i_dnnt(&work[2]);
#line 1146 "sgejsv.f"
	} else {

/*        .. two more QR factorizations ( one QRF is not enough, two require */
/*        accumulated product of Jacobi rotations, three are perfect ) */

#line 1151 "sgejsv.f"
	    i__1 = nr - 1;
#line 1151 "sgejsv.f"
	    i__2 = nr - 1;
#line 1151 "sgejsv.f"
	    slaset_("Lower", &i__1, &i__2, &c_b34, &c_b34, &a[a_dim1 + 2], 
		    lda, (ftnlen)5);
#line 1152 "sgejsv.f"
	    i__1 = *lwork - *n;
#line 1152 "sgejsv.f"
	    sgelqf_(&nr, n, &a[a_offset], lda, &work[1], &work[*n + 1], &i__1,
		     &ierr);
#line 1153 "sgejsv.f"
	    slacpy_("Lower", &nr, &nr, &a[a_offset], lda, &v[v_offset], ldv, (
		    ftnlen)5);
#line 1154 "sgejsv.f"
	    i__1 = nr - 1;
#line 1154 "sgejsv.f"
	    i__2 = nr - 1;
#line 1154 "sgejsv.f"
	    slaset_("Upper", &i__1, &i__2, &c_b34, &c_b34, &v[(v_dim1 << 1) + 
		    1], ldv, (ftnlen)5);
#line 1155 "sgejsv.f"
	    i__1 = *lwork - (*n << 1);
#line 1155 "sgejsv.f"
	    sgeqrf_(&nr, &nr, &v[v_offset], ldv, &work[*n + 1], &work[(*n << 
		    1) + 1], &i__1, &ierr);
#line 1157 "sgejsv.f"
	    i__1 = nr;
#line 1157 "sgejsv.f"
	    for (p = 1; p <= i__1; ++p) {
#line 1158 "sgejsv.f"
		i__2 = nr - p + 1;
#line 1158 "sgejsv.f"
		scopy_(&i__2, &v[p + p * v_dim1], ldv, &v[p + p * v_dim1], &
			c__1);
#line 1159 "sgejsv.f"
/* L8998: */
#line 1159 "sgejsv.f"
	    }
#line 1160 "sgejsv.f"
	    i__1 = nr - 1;
#line 1160 "sgejsv.f"
	    i__2 = nr - 1;
#line 1160 "sgejsv.f"
	    slaset_("Upper", &i__1, &i__2, &c_b34, &c_b34, &v[(v_dim1 << 1) + 
		    1], ldv, (ftnlen)5);

#line 1162 "sgejsv.f"
	    i__1 = *lwork - *n;
#line 1162 "sgejsv.f"
	    sgesvj_("Lower", "U", "N", &nr, &nr, &v[v_offset], ldv, &sva[1], &
		    nr, &u[u_offset], ldu, &work[*n + 1], &i__1, info, (
		    ftnlen)5, (ftnlen)1, (ftnlen)1);
#line 1164 "sgejsv.f"
	    scalem = work[*n + 1];
#line 1165 "sgejsv.f"
	    numrank = i_dnnt(&work[*n + 2]);
#line 1166 "sgejsv.f"
	    if (nr < *n) {
#line 1167 "sgejsv.f"
		i__1 = *n - nr;
#line 1167 "sgejsv.f"
		slaset_("A", &i__1, &nr, &c_b34, &c_b34, &v[nr + 1 + v_dim1], 
			ldv, (ftnlen)1);
#line 1168 "sgejsv.f"
		i__1 = *n - nr;
#line 1168 "sgejsv.f"
		slaset_("A", &nr, &i__1, &c_b34, &c_b34, &v[(nr + 1) * v_dim1 
			+ 1], ldv, (ftnlen)1);
#line 1169 "sgejsv.f"
		i__1 = *n - nr;
#line 1169 "sgejsv.f"
		i__2 = *n - nr;
#line 1169 "sgejsv.f"
		slaset_("A", &i__1, &i__2, &c_b34, &c_b35, &v[nr + 1 + (nr + 
			1) * v_dim1], ldv, (ftnlen)1);
#line 1170 "sgejsv.f"
	    }

#line 1172 "sgejsv.f"
	    i__1 = *lwork - *n;
#line 1172 "sgejsv.f"
	    sormlq_("Left", "Transpose", n, n, &nr, &a[a_offset], lda, &work[
		    1], &v[v_offset], ldv, &work[*n + 1], &i__1, &ierr, (
		    ftnlen)4, (ftnlen)9);

#line 1175 "sgejsv.f"
	}

#line 1177 "sgejsv.f"
	i__1 = *n;
#line 1177 "sgejsv.f"
	for (p = 1; p <= i__1; ++p) {
#line 1178 "sgejsv.f"
	    scopy_(n, &v[p + v_dim1], ldv, &a[iwork[p] + a_dim1], lda);
#line 1179 "sgejsv.f"
/* L8991: */
#line 1179 "sgejsv.f"
	}
#line 1180 "sgejsv.f"
	slacpy_("All", n, n, &a[a_offset], lda, &v[v_offset], ldv, (ftnlen)3);

#line 1182 "sgejsv.f"
	if (transp) {
#line 1183 "sgejsv.f"
	    slacpy_("All", n, n, &v[v_offset], ldv, &u[u_offset], ldu, (
		    ftnlen)3);
#line 1184 "sgejsv.f"
	}

#line 1186 "sgejsv.f"
    } else if (lsvec && ! rsvec) {

/*        .. Singular Values and Left Singular Vectors                 .. */

/*        .. second preconditioning step to avoid need to accumulate */
/*        Jacobi rotations in the Jacobi iterations. */
#line 1192 "sgejsv.f"
	i__1 = nr;
#line 1192 "sgejsv.f"
	for (p = 1; p <= i__1; ++p) {
#line 1193 "sgejsv.f"
	    i__2 = *n - p + 1;
#line 1193 "sgejsv.f"
	    scopy_(&i__2, &a[p + p * a_dim1], lda, &u[p + p * u_dim1], &c__1);
#line 1194 "sgejsv.f"
/* L1965: */
#line 1194 "sgejsv.f"
	}
#line 1195 "sgejsv.f"
	i__1 = nr - 1;
#line 1195 "sgejsv.f"
	i__2 = nr - 1;
#line 1195 "sgejsv.f"
	slaset_("Upper", &i__1, &i__2, &c_b34, &c_b34, &u[(u_dim1 << 1) + 1], 
		ldu, (ftnlen)5);

#line 1197 "sgejsv.f"
	i__1 = *lwork - (*n << 1);
#line 1197 "sgejsv.f"
	sgeqrf_(n, &nr, &u[u_offset], ldu, &work[*n + 1], &work[(*n << 1) + 1]
		, &i__1, &ierr);

#line 1200 "sgejsv.f"
	i__1 = nr - 1;
#line 1200 "sgejsv.f"
	for (p = 1; p <= i__1; ++p) {
#line 1201 "sgejsv.f"
	    i__2 = nr - p;
#line 1201 "sgejsv.f"
	    scopy_(&i__2, &u[p + (p + 1) * u_dim1], ldu, &u[p + 1 + p * 
		    u_dim1], &c__1);
#line 1202 "sgejsv.f"
/* L1967: */
#line 1202 "sgejsv.f"
	}
#line 1203 "sgejsv.f"
	i__1 = nr - 1;
#line 1203 "sgejsv.f"
	i__2 = nr - 1;
#line 1203 "sgejsv.f"
	slaset_("Upper", &i__1, &i__2, &c_b34, &c_b34, &u[(u_dim1 << 1) + 1], 
		ldu, (ftnlen)5);

#line 1205 "sgejsv.f"
	i__1 = *lwork - *n;
#line 1205 "sgejsv.f"
	sgesvj_("Lower", "U", "N", &nr, &nr, &u[u_offset], ldu, &sva[1], &nr, 
		&a[a_offset], lda, &work[*n + 1], &i__1, info, (ftnlen)5, (
		ftnlen)1, (ftnlen)1);
#line 1207 "sgejsv.f"
	scalem = work[*n + 1];
#line 1208 "sgejsv.f"
	numrank = i_dnnt(&work[*n + 2]);

#line 1210 "sgejsv.f"
	if (nr < *m) {
#line 1211 "sgejsv.f"
	    i__1 = *m - nr;
#line 1211 "sgejsv.f"
	    slaset_("A", &i__1, &nr, &c_b34, &c_b34, &u[nr + 1 + u_dim1], ldu,
		     (ftnlen)1);
#line 1212 "sgejsv.f"
	    if (nr < n1) {
#line 1213 "sgejsv.f"
		i__1 = n1 - nr;
#line 1213 "sgejsv.f"
		slaset_("A", &nr, &i__1, &c_b34, &c_b34, &u[(nr + 1) * u_dim1 
			+ 1], ldu, (ftnlen)1);
#line 1214 "sgejsv.f"
		i__1 = *m - nr;
#line 1214 "sgejsv.f"
		i__2 = n1 - nr;
#line 1214 "sgejsv.f"
		slaset_("A", &i__1, &i__2, &c_b34, &c_b35, &u[nr + 1 + (nr + 
			1) * u_dim1], ldu, (ftnlen)1);
#line 1215 "sgejsv.f"
	    }
#line 1216 "sgejsv.f"
	}

#line 1218 "sgejsv.f"
	i__1 = *lwork - *n;
#line 1218 "sgejsv.f"
	sormqr_("Left", "No Tr", m, &n1, n, &a[a_offset], lda, &work[1], &u[
		u_offset], ldu, &work[*n + 1], &i__1, &ierr, (ftnlen)4, (
		ftnlen)5);

#line 1221 "sgejsv.f"
	if (rowpiv) {
#line 1221 "sgejsv.f"
	    i__1 = *m - 1;
#line 1221 "sgejsv.f"
	    slaswp_(&n1, &u[u_offset], ldu, &c__1, &i__1, &iwork[(*n << 1) + 
		    1], &c_n1);
#line 1221 "sgejsv.f"
	}

#line 1224 "sgejsv.f"
	i__1 = n1;
#line 1224 "sgejsv.f"
	for (p = 1; p <= i__1; ++p) {
#line 1225 "sgejsv.f"
	    xsc = 1. / snrm2_(m, &u[p * u_dim1 + 1], &c__1);
#line 1226 "sgejsv.f"
	    sscal_(m, &xsc, &u[p * u_dim1 + 1], &c__1);
#line 1227 "sgejsv.f"
/* L1974: */
#line 1227 "sgejsv.f"
	}

#line 1229 "sgejsv.f"
	if (transp) {
#line 1230 "sgejsv.f"
	    slacpy_("All", n, n, &u[u_offset], ldu, &v[v_offset], ldv, (
		    ftnlen)3);
#line 1231 "sgejsv.f"
	}

#line 1233 "sgejsv.f"
    } else {

/*        .. Full SVD .. */

#line 1237 "sgejsv.f"
	if (! jracc) {

#line 1239 "sgejsv.f"
	    if (! almort) {

/*           Second Preconditioning Step (QRF [with pivoting]) */
/*           Note that the composition of TRANSPOSE, QRF and TRANSPOSE is */
/*           equivalent to an LQF CALL. Since in many libraries the QRF */
/*           seems to be better optimized than the LQF, we do explicit */
/*           transpose and use the QRF. This is subject to changes in an */
/*           optimized implementation of SGEJSV. */

#line 1248 "sgejsv.f"
		i__1 = nr;
#line 1248 "sgejsv.f"
		for (p = 1; p <= i__1; ++p) {
#line 1249 "sgejsv.f"
		    i__2 = *n - p + 1;
#line 1249 "sgejsv.f"
		    scopy_(&i__2, &a[p + p * a_dim1], lda, &v[p + p * v_dim1],
			     &c__1);
#line 1250 "sgejsv.f"
/* L1968: */
#line 1250 "sgejsv.f"
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

#line 1264 "sgejsv.f"
		if (l2pert) {
#line 1265 "sgejsv.f"
		    xsc = sqrt(small);
#line 1266 "sgejsv.f"
		    i__1 = nr;
#line 1266 "sgejsv.f"
		    for (q = 1; q <= i__1; ++q) {
#line 1267 "sgejsv.f"
			temp1 = xsc * (d__1 = v[q + q * v_dim1], abs(d__1));
#line 1268 "sgejsv.f"
			i__2 = *n;
#line 1268 "sgejsv.f"
			for (p = 1; p <= i__2; ++p) {
#line 1269 "sgejsv.f"
			    if (p > q && (d__1 = v[p + q * v_dim1], abs(d__1))
				     <= temp1 || p < q) {
#line 1269 "sgejsv.f"
				v[p + q * v_dim1] = d_sign(&temp1, &v[p + q * 
					v_dim1]);
#line 1269 "sgejsv.f"
			    }
#line 1272 "sgejsv.f"
			    if (p < q) {
#line 1272 "sgejsv.f"
				v[p + q * v_dim1] = -v[p + q * v_dim1];
#line 1272 "sgejsv.f"
			    }
#line 1273 "sgejsv.f"
/* L2968: */
#line 1273 "sgejsv.f"
			}
#line 1274 "sgejsv.f"
/* L2969: */
#line 1274 "sgejsv.f"
		    }
#line 1275 "sgejsv.f"
		} else {
#line 1276 "sgejsv.f"
		    i__1 = nr - 1;
#line 1276 "sgejsv.f"
		    i__2 = nr - 1;
#line 1276 "sgejsv.f"
		    slaset_("U", &i__1, &i__2, &c_b34, &c_b34, &v[(v_dim1 << 
			    1) + 1], ldv, (ftnlen)1);
#line 1277 "sgejsv.f"
		}

/*           Estimate the row scaled condition number of R1 */
/*           (If R1 is rectangular, N > NR, then the condition number */
/*           of the leading NR x NR submatrix is estimated.) */

#line 1283 "sgejsv.f"
		slacpy_("L", &nr, &nr, &v[v_offset], ldv, &work[(*n << 1) + 1]
			, &nr, (ftnlen)1);
#line 1284 "sgejsv.f"
		i__1 = nr;
#line 1284 "sgejsv.f"
		for (p = 1; p <= i__1; ++p) {
#line 1285 "sgejsv.f"
		    i__2 = nr - p + 1;
#line 1285 "sgejsv.f"
		    temp1 = snrm2_(&i__2, &work[(*n << 1) + (p - 1) * nr + p],
			     &c__1);
#line 1286 "sgejsv.f"
		    i__2 = nr - p + 1;
#line 1286 "sgejsv.f"
		    d__1 = 1. / temp1;
#line 1286 "sgejsv.f"
		    sscal_(&i__2, &d__1, &work[(*n << 1) + (p - 1) * nr + p], 
			    &c__1);
#line 1287 "sgejsv.f"
/* L3950: */
#line 1287 "sgejsv.f"
		}
#line 1288 "sgejsv.f"
		spocon_("Lower", &nr, &work[(*n << 1) + 1], &nr, &c_b35, &
			temp1, &work[(*n << 1) + nr * nr + 1], &iwork[*m + (*
			n << 1) + 1], &ierr, (ftnlen)5);
#line 1290 "sgejsv.f"
		condr1 = 1. / sqrt(temp1);
/*           .. here need a second oppinion on the condition number */
/*           .. then assume worst case scenario */
/*           R1 is OK for inverse <=> CONDR1 .LT. FLOAT(N) */
/*           more conservative    <=> CONDR1 .LT. SQRT(FLOAT(N)) */

#line 1296 "sgejsv.f"
		cond_ok__ = sqrt((doublereal) nr);
/* [TP]       COND_OK is a tuning parameter. */
#line 1299 "sgejsv.f"
		if (condr1 < cond_ok__) {
/*              .. the second QRF without pivoting. Note: in an optimized */
/*              implementation, this QRF should be implemented as the QRF */
/*              of a lower triangular matrix. */
/*              R1^t = Q2 * R2 */
#line 1304 "sgejsv.f"
		    i__1 = *lwork - (*n << 1);
#line 1304 "sgejsv.f"
		    sgeqrf_(n, &nr, &v[v_offset], ldv, &work[*n + 1], &work[(*
			    n << 1) + 1], &i__1, &ierr);

#line 1307 "sgejsv.f"
		    if (l2pert) {
#line 1308 "sgejsv.f"
			xsc = sqrt(small) / epsln;
#line 1309 "sgejsv.f"
			i__1 = nr;
#line 1309 "sgejsv.f"
			for (p = 2; p <= i__1; ++p) {
#line 1310 "sgejsv.f"
			    i__2 = p - 1;
#line 1310 "sgejsv.f"
			    for (q = 1; q <= i__2; ++q) {
/* Computing MIN */
#line 1311 "sgejsv.f"
				d__3 = (d__1 = v[p + p * v_dim1], abs(d__1)), 
					d__4 = (d__2 = v[q + q * v_dim1], abs(
					d__2));
#line 1311 "sgejsv.f"
				temp1 = xsc * min(d__3,d__4);
#line 1312 "sgejsv.f"
				if ((d__1 = v[q + p * v_dim1], abs(d__1)) <= 
					temp1) {
#line 1312 "sgejsv.f"
				    v[q + p * v_dim1] = d_sign(&temp1, &v[q + 
					    p * v_dim1]);
#line 1312 "sgejsv.f"
				}
#line 1314 "sgejsv.f"
/* L3958: */
#line 1314 "sgejsv.f"
			    }
#line 1315 "sgejsv.f"
/* L3959: */
#line 1315 "sgejsv.f"
			}
#line 1316 "sgejsv.f"
		    }

#line 1318 "sgejsv.f"
		    if (nr != *n) {
#line 1318 "sgejsv.f"
			slacpy_("A", n, &nr, &v[v_offset], ldv, &work[(*n << 
				1) + 1], n, (ftnlen)1);
#line 1318 "sgejsv.f"
		    }
/*              .. save ... */

/*           .. this transposed copy should be better than naive */
#line 1323 "sgejsv.f"
		    i__1 = nr - 1;
#line 1323 "sgejsv.f"
		    for (p = 1; p <= i__1; ++p) {
#line 1324 "sgejsv.f"
			i__2 = nr - p;
#line 1324 "sgejsv.f"
			scopy_(&i__2, &v[p + (p + 1) * v_dim1], ldv, &v[p + 1 
				+ p * v_dim1], &c__1);
#line 1325 "sgejsv.f"
/* L1969: */
#line 1325 "sgejsv.f"
		    }

#line 1327 "sgejsv.f"
		    condr2 = condr1;

#line 1329 "sgejsv.f"
		} else {

/*              .. ill-conditioned case: second QRF with pivoting */
/*              Note that windowed pivoting would be equaly good */
/*              numerically, and more run-time efficient. So, in */
/*              an optimal implementation, the next call to SGEQP3 */
/*              should be replaced with eg. CALL SGEQPX (ACM TOMS #782) */
/*              with properly (carefully) chosen parameters. */

/*              R1^t * P2 = Q2 * R2 */
#line 1339 "sgejsv.f"
		    i__1 = nr;
#line 1339 "sgejsv.f"
		    for (p = 1; p <= i__1; ++p) {
#line 1340 "sgejsv.f"
			iwork[*n + p] = 0;
#line 1341 "sgejsv.f"
/* L3003: */
#line 1341 "sgejsv.f"
		    }
#line 1342 "sgejsv.f"
		    i__1 = *lwork - (*n << 1);
#line 1342 "sgejsv.f"
		    sgeqp3_(n, &nr, &v[v_offset], ldv, &iwork[*n + 1], &work[*
			    n + 1], &work[(*n << 1) + 1], &i__1, &ierr);
/* *               CALL SGEQRF( N, NR, V, LDV, WORK(N+1), WORK(2*N+1), */
/* *     $              LWORK-2*N, IERR ) */
#line 1346 "sgejsv.f"
		    if (l2pert) {
#line 1347 "sgejsv.f"
			xsc = sqrt(small);
#line 1348 "sgejsv.f"
			i__1 = nr;
#line 1348 "sgejsv.f"
			for (p = 2; p <= i__1; ++p) {
#line 1349 "sgejsv.f"
			    i__2 = p - 1;
#line 1349 "sgejsv.f"
			    for (q = 1; q <= i__2; ++q) {
/* Computing MIN */
#line 1350 "sgejsv.f"
				d__3 = (d__1 = v[p + p * v_dim1], abs(d__1)), 
					d__4 = (d__2 = v[q + q * v_dim1], abs(
					d__2));
#line 1350 "sgejsv.f"
				temp1 = xsc * min(d__3,d__4);
#line 1351 "sgejsv.f"
				if ((d__1 = v[q + p * v_dim1], abs(d__1)) <= 
					temp1) {
#line 1351 "sgejsv.f"
				    v[q + p * v_dim1] = d_sign(&temp1, &v[q + 
					    p * v_dim1]);
#line 1351 "sgejsv.f"
				}
#line 1353 "sgejsv.f"
/* L3968: */
#line 1353 "sgejsv.f"
			    }
#line 1354 "sgejsv.f"
/* L3969: */
#line 1354 "sgejsv.f"
			}
#line 1355 "sgejsv.f"
		    }

#line 1357 "sgejsv.f"
		    slacpy_("A", n, &nr, &v[v_offset], ldv, &work[(*n << 1) + 
			    1], n, (ftnlen)1);

#line 1359 "sgejsv.f"
		    if (l2pert) {
#line 1360 "sgejsv.f"
			xsc = sqrt(small);
#line 1361 "sgejsv.f"
			i__1 = nr;
#line 1361 "sgejsv.f"
			for (p = 2; p <= i__1; ++p) {
#line 1362 "sgejsv.f"
			    i__2 = p - 1;
#line 1362 "sgejsv.f"
			    for (q = 1; q <= i__2; ++q) {
/* Computing MIN */
#line 1363 "sgejsv.f"
				d__3 = (d__1 = v[p + p * v_dim1], abs(d__1)), 
					d__4 = (d__2 = v[q + q * v_dim1], abs(
					d__2));
#line 1363 "sgejsv.f"
				temp1 = xsc * min(d__3,d__4);
#line 1364 "sgejsv.f"
				v[p + q * v_dim1] = -d_sign(&temp1, &v[q + p *
					 v_dim1]);
#line 1365 "sgejsv.f"
/* L8971: */
#line 1365 "sgejsv.f"
			    }
#line 1366 "sgejsv.f"
/* L8970: */
#line 1366 "sgejsv.f"
			}
#line 1367 "sgejsv.f"
		    } else {
#line 1368 "sgejsv.f"
			i__1 = nr - 1;
#line 1368 "sgejsv.f"
			i__2 = nr - 1;
#line 1368 "sgejsv.f"
			slaset_("L", &i__1, &i__2, &c_b34, &c_b34, &v[v_dim1 
				+ 2], ldv, (ftnlen)1);
#line 1369 "sgejsv.f"
		    }
/*              Now, compute R2 = L3 * Q3, the LQ factorization. */
#line 1371 "sgejsv.f"
		    i__1 = *lwork - (*n << 1) - *n * nr - nr;
#line 1371 "sgejsv.f"
		    sgelqf_(&nr, &nr, &v[v_offset], ldv, &work[(*n << 1) + *n 
			    * nr + 1], &work[(*n << 1) + *n * nr + nr + 1], &
			    i__1, &ierr);
/*              .. and estimate the condition number */
#line 1374 "sgejsv.f"
		    slacpy_("L", &nr, &nr, &v[v_offset], ldv, &work[(*n << 1) 
			    + *n * nr + nr + 1], &nr, (ftnlen)1);
#line 1375 "sgejsv.f"
		    i__1 = nr;
#line 1375 "sgejsv.f"
		    for (p = 1; p <= i__1; ++p) {
#line 1376 "sgejsv.f"
			temp1 = snrm2_(&p, &work[(*n << 1) + *n * nr + nr + p]
				, &nr);
#line 1377 "sgejsv.f"
			d__1 = 1. / temp1;
#line 1377 "sgejsv.f"
			sscal_(&p, &d__1, &work[(*n << 1) + *n * nr + nr + p],
				 &nr);
#line 1378 "sgejsv.f"
/* L4950: */
#line 1378 "sgejsv.f"
		    }
#line 1379 "sgejsv.f"
		    spocon_("L", &nr, &work[(*n << 1) + *n * nr + nr + 1], &
			    nr, &c_b35, &temp1, &work[(*n << 1) + *n * nr + 
			    nr + nr * nr + 1], &iwork[*m + (*n << 1) + 1], &
			    ierr, (ftnlen)1);
#line 1381 "sgejsv.f"
		    condr2 = 1. / sqrt(temp1);

#line 1383 "sgejsv.f"
		    if (condr2 >= cond_ok__) {
/*                 .. save the Householder vectors used for Q3 */
/*                 (this overwrittes the copy of R2, as it will not be */
/*                 needed in this branch, but it does not overwritte the */
/*                 Huseholder vectors of Q2.). */
#line 1388 "sgejsv.f"
			slacpy_("U", &nr, &nr, &v[v_offset], ldv, &work[(*n <<
				 1) + 1], n, (ftnlen)1);
/*                 .. and the rest of the information on Q3 is in */
/*                 WORK(2*N+N*NR+1:2*N+N*NR+N) */
#line 1391 "sgejsv.f"
		    }

#line 1393 "sgejsv.f"
		}

#line 1395 "sgejsv.f"
		if (l2pert) {
#line 1396 "sgejsv.f"
		    xsc = sqrt(small);
#line 1397 "sgejsv.f"
		    i__1 = nr;
#line 1397 "sgejsv.f"
		    for (q = 2; q <= i__1; ++q) {
#line 1398 "sgejsv.f"
			temp1 = xsc * v[q + q * v_dim1];
#line 1399 "sgejsv.f"
			i__2 = q - 1;
#line 1399 "sgejsv.f"
			for (p = 1; p <= i__2; ++p) {
/*                    V(p,q) = - SIGN( TEMP1, V(q,p) ) */
#line 1401 "sgejsv.f"
			    v[p + q * v_dim1] = -d_sign(&temp1, &v[p + q * 
				    v_dim1]);
#line 1402 "sgejsv.f"
/* L4969: */
#line 1402 "sgejsv.f"
			}
#line 1403 "sgejsv.f"
/* L4968: */
#line 1403 "sgejsv.f"
		    }
#line 1404 "sgejsv.f"
		} else {
#line 1405 "sgejsv.f"
		    i__1 = nr - 1;
#line 1405 "sgejsv.f"
		    i__2 = nr - 1;
#line 1405 "sgejsv.f"
		    slaset_("U", &i__1, &i__2, &c_b34, &c_b34, &v[(v_dim1 << 
			    1) + 1], ldv, (ftnlen)1);
#line 1406 "sgejsv.f"
		}

/*        Second preconditioning finished; continue with Jacobi SVD */
/*        The input matrix is lower trinagular. */

/*        Recover the right singular vectors as solution of a well */
/*        conditioned triangular matrix equation. */

#line 1414 "sgejsv.f"
		if (condr1 < cond_ok__) {

#line 1416 "sgejsv.f"
		    i__1 = *lwork - (*n << 1) - *n * nr - nr;
#line 1416 "sgejsv.f"
		    sgesvj_("L", "U", "N", &nr, &nr, &v[v_offset], ldv, &sva[
			    1], &nr, &u[u_offset], ldu, &work[(*n << 1) + *n *
			     nr + nr + 1], &i__1, info, (ftnlen)1, (ftnlen)1, 
			    (ftnlen)1);
#line 1418 "sgejsv.f"
		    scalem = work[(*n << 1) + *n * nr + nr + 1];
#line 1419 "sgejsv.f"
		    numrank = i_dnnt(&work[(*n << 1) + *n * nr + nr + 2]);
#line 1420 "sgejsv.f"
		    i__1 = nr;
#line 1420 "sgejsv.f"
		    for (p = 1; p <= i__1; ++p) {
#line 1421 "sgejsv.f"
			scopy_(&nr, &v[p * v_dim1 + 1], &c__1, &u[p * u_dim1 
				+ 1], &c__1);
#line 1422 "sgejsv.f"
			sscal_(&nr, &sva[p], &v[p * v_dim1 + 1], &c__1);
#line 1423 "sgejsv.f"
/* L3970: */
#line 1423 "sgejsv.f"
		    }
/*        .. pick the right matrix equation and solve it */

#line 1427 "sgejsv.f"
		    if (nr == *n) {
/* :))             .. best case, R1 is inverted. The solution of this matrix */
/*                 equation is Q2*V2 = the product of the Jacobi rotations */
/*                 used in SGESVJ, premultiplied with the orthogonal matrix */
/*                 from the second QR factorization. */
#line 1432 "sgejsv.f"
			strsm_("L", "U", "N", "N", &nr, &nr, &c_b35, &a[
				a_offset], lda, &v[v_offset], ldv, (ftnlen)1, 
				(ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1433 "sgejsv.f"
		    } else {
/*                 .. R1 is well conditioned, but non-square. Transpose(R2) */
/*                 is inverted to get the product of the Jacobi rotations */
/*                 used in SGESVJ. The Q-factor from the second QR */
/*                 factorization is then built in explicitly. */
#line 1438 "sgejsv.f"
			strsm_("L", "U", "T", "N", &nr, &nr, &c_b35, &work[(*
				n << 1) + 1], n, &v[v_offset], ldv, (ftnlen)1,
				 (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1440 "sgejsv.f"
			if (nr < *n) {
#line 1441 "sgejsv.f"
			    i__1 = *n - nr;
#line 1441 "sgejsv.f"
			    slaset_("A", &i__1, &nr, &c_b34, &c_b34, &v[nr + 
				    1 + v_dim1], ldv, (ftnlen)1);
#line 1442 "sgejsv.f"
			    i__1 = *n - nr;
#line 1442 "sgejsv.f"
			    slaset_("A", &nr, &i__1, &c_b34, &c_b34, &v[(nr + 
				    1) * v_dim1 + 1], ldv, (ftnlen)1);
#line 1443 "sgejsv.f"
			    i__1 = *n - nr;
#line 1443 "sgejsv.f"
			    i__2 = *n - nr;
#line 1443 "sgejsv.f"
			    slaset_("A", &i__1, &i__2, &c_b34, &c_b35, &v[nr 
				    + 1 + (nr + 1) * v_dim1], ldv, (ftnlen)1);
#line 1444 "sgejsv.f"
			}
#line 1445 "sgejsv.f"
			i__1 = *lwork - (*n << 1) - *n * nr - nr;
#line 1445 "sgejsv.f"
			sormqr_("L", "N", n, n, &nr, &work[(*n << 1) + 1], n, 
				&work[*n + 1], &v[v_offset], ldv, &work[(*n <<
				 1) + *n * nr + nr + 1], &i__1, &ierr, (
				ftnlen)1, (ftnlen)1);
#line 1447 "sgejsv.f"
		    }

#line 1449 "sgejsv.f"
		} else if (condr2 < cond_ok__) {

/* :)           .. the input matrix A is very likely a relative of */
/*              the Kahan matrix :) */
/*              The matrix R2 is inverted. The solution of the matrix equation */
/*              is Q3^T*V3 = the product of the Jacobi rotations (appplied to */
/*              the lower triangular L3 from the LQ factorization of */
/*              R2=L3*Q3), pre-multiplied with the transposed Q3. */
#line 1457 "sgejsv.f"
		    i__1 = *lwork - (*n << 1) - *n * nr - nr;
#line 1457 "sgejsv.f"
		    sgesvj_("L", "U", "N", &nr, &nr, &v[v_offset], ldv, &sva[
			    1], &nr, &u[u_offset], ldu, &work[(*n << 1) + *n *
			     nr + nr + 1], &i__1, info, (ftnlen)1, (ftnlen)1, 
			    (ftnlen)1);
#line 1459 "sgejsv.f"
		    scalem = work[(*n << 1) + *n * nr + nr + 1];
#line 1460 "sgejsv.f"
		    numrank = i_dnnt(&work[(*n << 1) + *n * nr + nr + 2]);
#line 1461 "sgejsv.f"
		    i__1 = nr;
#line 1461 "sgejsv.f"
		    for (p = 1; p <= i__1; ++p) {
#line 1462 "sgejsv.f"
			scopy_(&nr, &v[p * v_dim1 + 1], &c__1, &u[p * u_dim1 
				+ 1], &c__1);
#line 1463 "sgejsv.f"
			sscal_(&nr, &sva[p], &u[p * u_dim1 + 1], &c__1);
#line 1464 "sgejsv.f"
/* L3870: */
#line 1464 "sgejsv.f"
		    }
#line 1465 "sgejsv.f"
		    strsm_("L", "U", "N", "N", &nr, &nr, &c_b35, &work[(*n << 
			    1) + 1], n, &u[u_offset], ldu, (ftnlen)1, (ftnlen)
			    1, (ftnlen)1, (ftnlen)1);
/*              .. apply the permutation from the second QR factorization */
#line 1467 "sgejsv.f"
		    i__1 = nr;
#line 1467 "sgejsv.f"
		    for (q = 1; q <= i__1; ++q) {
#line 1468 "sgejsv.f"
			i__2 = nr;
#line 1468 "sgejsv.f"
			for (p = 1; p <= i__2; ++p) {
#line 1469 "sgejsv.f"
			    work[(*n << 1) + *n * nr + nr + iwork[*n + p]] = 
				    u[p + q * u_dim1];
#line 1470 "sgejsv.f"
/* L872: */
#line 1470 "sgejsv.f"
			}
#line 1471 "sgejsv.f"
			i__2 = nr;
#line 1471 "sgejsv.f"
			for (p = 1; p <= i__2; ++p) {
#line 1472 "sgejsv.f"
			    u[p + q * u_dim1] = work[(*n << 1) + *n * nr + nr 
				    + p];
#line 1473 "sgejsv.f"
/* L874: */
#line 1473 "sgejsv.f"
			}
#line 1474 "sgejsv.f"
/* L873: */
#line 1474 "sgejsv.f"
		    }
#line 1475 "sgejsv.f"
		    if (nr < *n) {
#line 1476 "sgejsv.f"
			i__1 = *n - nr;
#line 1476 "sgejsv.f"
			slaset_("A", &i__1, &nr, &c_b34, &c_b34, &v[nr + 1 + 
				v_dim1], ldv, (ftnlen)1);
#line 1477 "sgejsv.f"
			i__1 = *n - nr;
#line 1477 "sgejsv.f"
			slaset_("A", &nr, &i__1, &c_b34, &c_b34, &v[(nr + 1) *
				 v_dim1 + 1], ldv, (ftnlen)1);
#line 1478 "sgejsv.f"
			i__1 = *n - nr;
#line 1478 "sgejsv.f"
			i__2 = *n - nr;
#line 1478 "sgejsv.f"
			slaset_("A", &i__1, &i__2, &c_b34, &c_b35, &v[nr + 1 
				+ (nr + 1) * v_dim1], ldv, (ftnlen)1);
#line 1479 "sgejsv.f"
		    }
#line 1480 "sgejsv.f"
		    i__1 = *lwork - (*n << 1) - *n * nr - nr;
#line 1480 "sgejsv.f"
		    sormqr_("L", "N", n, n, &nr, &work[(*n << 1) + 1], n, &
			    work[*n + 1], &v[v_offset], ldv, &work[(*n << 1) 
			    + *n * nr + nr + 1], &i__1, &ierr, (ftnlen)1, (
			    ftnlen)1);
#line 1482 "sgejsv.f"
		} else {
/*              Last line of defense. */
/* #:(          This is a rather pathological case: no scaled condition */
/*              improvement after two pivoted QR factorizations. Other */
/*              possibility is that the rank revealing QR factorization */
/*              or the condition estimator has failed, or the COND_OK */
/*              is set very close to ONE (which is unnecessary). Normally, */
/*              this branch should never be executed, but in rare cases of */
/*              failure of the RRQR or condition estimator, the last line of */
/*              defense ensures that SGEJSV completes the task. */
/*              Compute the full SVD of L3 using SGESVJ with explicit */
/*              accumulation of Jacobi rotations. */
#line 1494 "sgejsv.f"
		    i__1 = *lwork - (*n << 1) - *n * nr - nr;
#line 1494 "sgejsv.f"
		    sgesvj_("L", "U", "V", &nr, &nr, &v[v_offset], ldv, &sva[
			    1], &nr, &u[u_offset], ldu, &work[(*n << 1) + *n *
			     nr + nr + 1], &i__1, info, (ftnlen)1, (ftnlen)1, 
			    (ftnlen)1);
#line 1496 "sgejsv.f"
		    scalem = work[(*n << 1) + *n * nr + nr + 1];
#line 1497 "sgejsv.f"
		    numrank = i_dnnt(&work[(*n << 1) + *n * nr + nr + 2]);
#line 1498 "sgejsv.f"
		    if (nr < *n) {
#line 1499 "sgejsv.f"
			i__1 = *n - nr;
#line 1499 "sgejsv.f"
			slaset_("A", &i__1, &nr, &c_b34, &c_b34, &v[nr + 1 + 
				v_dim1], ldv, (ftnlen)1);
#line 1500 "sgejsv.f"
			i__1 = *n - nr;
#line 1500 "sgejsv.f"
			slaset_("A", &nr, &i__1, &c_b34, &c_b34, &v[(nr + 1) *
				 v_dim1 + 1], ldv, (ftnlen)1);
#line 1501 "sgejsv.f"
			i__1 = *n - nr;
#line 1501 "sgejsv.f"
			i__2 = *n - nr;
#line 1501 "sgejsv.f"
			slaset_("A", &i__1, &i__2, &c_b34, &c_b35, &v[nr + 1 
				+ (nr + 1) * v_dim1], ldv, (ftnlen)1);
#line 1502 "sgejsv.f"
		    }
#line 1503 "sgejsv.f"
		    i__1 = *lwork - (*n << 1) - *n * nr - nr;
#line 1503 "sgejsv.f"
		    sormqr_("L", "N", n, n, &nr, &work[(*n << 1) + 1], n, &
			    work[*n + 1], &v[v_offset], ldv, &work[(*n << 1) 
			    + *n * nr + nr + 1], &i__1, &ierr, (ftnlen)1, (
			    ftnlen)1);

#line 1506 "sgejsv.f"
		    i__1 = *lwork - (*n << 1) - *n * nr - nr;
#line 1506 "sgejsv.f"
		    sormlq_("L", "T", &nr, &nr, &nr, &work[(*n << 1) + 1], n, 
			    &work[(*n << 1) + *n * nr + 1], &u[u_offset], ldu,
			     &work[(*n << 1) + *n * nr + nr + 1], &i__1, &
			    ierr, (ftnlen)1, (ftnlen)1);
#line 1509 "sgejsv.f"
		    i__1 = nr;
#line 1509 "sgejsv.f"
		    for (q = 1; q <= i__1; ++q) {
#line 1510 "sgejsv.f"
			i__2 = nr;
#line 1510 "sgejsv.f"
			for (p = 1; p <= i__2; ++p) {
#line 1511 "sgejsv.f"
			    work[(*n << 1) + *n * nr + nr + iwork[*n + p]] = 
				    u[p + q * u_dim1];
#line 1512 "sgejsv.f"
/* L772: */
#line 1512 "sgejsv.f"
			}
#line 1513 "sgejsv.f"
			i__2 = nr;
#line 1513 "sgejsv.f"
			for (p = 1; p <= i__2; ++p) {
#line 1514 "sgejsv.f"
			    u[p + q * u_dim1] = work[(*n << 1) + *n * nr + nr 
				    + p];
#line 1515 "sgejsv.f"
/* L774: */
#line 1515 "sgejsv.f"
			}
#line 1516 "sgejsv.f"
/* L773: */
#line 1516 "sgejsv.f"
		    }

#line 1518 "sgejsv.f"
		}

/*           Permute the rows of V using the (column) permutation from the */
/*           first QRF. Also, scale the columns to make them unit in */
/*           Euclidean norm. This applies to all cases. */

#line 1524 "sgejsv.f"
		temp1 = sqrt((doublereal) (*n)) * epsln;
#line 1525 "sgejsv.f"
		i__1 = *n;
#line 1525 "sgejsv.f"
		for (q = 1; q <= i__1; ++q) {
#line 1526 "sgejsv.f"
		    i__2 = *n;
#line 1526 "sgejsv.f"
		    for (p = 1; p <= i__2; ++p) {
#line 1527 "sgejsv.f"
			work[(*n << 1) + *n * nr + nr + iwork[p]] = v[p + q * 
				v_dim1];
#line 1528 "sgejsv.f"
/* L972: */
#line 1528 "sgejsv.f"
		    }
#line 1529 "sgejsv.f"
		    i__2 = *n;
#line 1529 "sgejsv.f"
		    for (p = 1; p <= i__2; ++p) {
#line 1530 "sgejsv.f"
			v[p + q * v_dim1] = work[(*n << 1) + *n * nr + nr + p]
				;
#line 1531 "sgejsv.f"
/* L973: */
#line 1531 "sgejsv.f"
		    }
#line 1532 "sgejsv.f"
		    xsc = 1. / snrm2_(n, &v[q * v_dim1 + 1], &c__1);
#line 1533 "sgejsv.f"
		    if (xsc < 1. - temp1 || xsc > temp1 + 1.) {
#line 1533 "sgejsv.f"
			sscal_(n, &xsc, &v[q * v_dim1 + 1], &c__1);
#line 1533 "sgejsv.f"
		    }
#line 1535 "sgejsv.f"
/* L1972: */
#line 1535 "sgejsv.f"
		}
/*           At this moment, V contains the right singular vectors of A. */
/*           Next, assemble the left singular vector matrix U (M x N). */
#line 1538 "sgejsv.f"
		if (nr < *m) {
#line 1539 "sgejsv.f"
		    i__1 = *m - nr;
#line 1539 "sgejsv.f"
		    slaset_("A", &i__1, &nr, &c_b34, &c_b34, &u[nr + 1 + 
			    u_dim1], ldu, (ftnlen)1);
#line 1540 "sgejsv.f"
		    if (nr < n1) {
#line 1541 "sgejsv.f"
			i__1 = n1 - nr;
#line 1541 "sgejsv.f"
			slaset_("A", &nr, &i__1, &c_b34, &c_b34, &u[(nr + 1) *
				 u_dim1 + 1], ldu, (ftnlen)1);
#line 1542 "sgejsv.f"
			i__1 = *m - nr;
#line 1542 "sgejsv.f"
			i__2 = n1 - nr;
#line 1542 "sgejsv.f"
			slaset_("A", &i__1, &i__2, &c_b34, &c_b35, &u[nr + 1 
				+ (nr + 1) * u_dim1], ldu, (ftnlen)1);
#line 1543 "sgejsv.f"
		    }
#line 1544 "sgejsv.f"
		}

/*           The Q matrix from the first QRF is built into the left singular */
/*           matrix U. This applies to all cases. */

#line 1549 "sgejsv.f"
		i__1 = *lwork - *n;
#line 1549 "sgejsv.f"
		sormqr_("Left", "No_Tr", m, &n1, n, &a[a_offset], lda, &work[
			1], &u[u_offset], ldu, &work[*n + 1], &i__1, &ierr, (
			ftnlen)4, (ftnlen)5);
/*           The columns of U are normalized. The cost is O(M*N) flops. */
#line 1553 "sgejsv.f"
		temp1 = sqrt((doublereal) (*m)) * epsln;
#line 1554 "sgejsv.f"
		i__1 = nr;
#line 1554 "sgejsv.f"
		for (p = 1; p <= i__1; ++p) {
#line 1555 "sgejsv.f"
		    xsc = 1. / snrm2_(m, &u[p * u_dim1 + 1], &c__1);
#line 1556 "sgejsv.f"
		    if (xsc < 1. - temp1 || xsc > temp1 + 1.) {
#line 1556 "sgejsv.f"
			sscal_(m, &xsc, &u[p * u_dim1 + 1], &c__1);
#line 1556 "sgejsv.f"
		    }
#line 1558 "sgejsv.f"
/* L1973: */
#line 1558 "sgejsv.f"
		}

/*           If the initial QRF is computed with row pivoting, the left */
/*           singular vectors must be adjusted. */

#line 1563 "sgejsv.f"
		if (rowpiv) {
#line 1563 "sgejsv.f"
		    i__1 = *m - 1;
#line 1563 "sgejsv.f"
		    slaswp_(&n1, &u[u_offset], ldu, &c__1, &i__1, &iwork[(*n 
			    << 1) + 1], &c_n1);
#line 1563 "sgejsv.f"
		}

#line 1566 "sgejsv.f"
	    } else {

/*        .. the initial matrix A has almost orthogonal columns and */
/*        the second QRF is not needed */

#line 1571 "sgejsv.f"
		slacpy_("Upper", n, n, &a[a_offset], lda, &work[*n + 1], n, (
			ftnlen)5);
#line 1572 "sgejsv.f"
		if (l2pert) {
#line 1573 "sgejsv.f"
		    xsc = sqrt(small);
#line 1574 "sgejsv.f"
		    i__1 = *n;
#line 1574 "sgejsv.f"
		    for (p = 2; p <= i__1; ++p) {
#line 1575 "sgejsv.f"
			temp1 = xsc * work[*n + (p - 1) * *n + p];
#line 1576 "sgejsv.f"
			i__2 = p - 1;
#line 1576 "sgejsv.f"
			for (q = 1; q <= i__2; ++q) {
#line 1577 "sgejsv.f"
			    work[*n + (q - 1) * *n + p] = -d_sign(&temp1, &
				    work[*n + (p - 1) * *n + q]);
#line 1578 "sgejsv.f"
/* L5971: */
#line 1578 "sgejsv.f"
			}
#line 1579 "sgejsv.f"
/* L5970: */
#line 1579 "sgejsv.f"
		    }
#line 1580 "sgejsv.f"
		} else {
#line 1581 "sgejsv.f"
		    i__1 = *n - 1;
#line 1581 "sgejsv.f"
		    i__2 = *n - 1;
#line 1581 "sgejsv.f"
		    slaset_("Lower", &i__1, &i__2, &c_b34, &c_b34, &work[*n + 
			    2], n, (ftnlen)5);
#line 1582 "sgejsv.f"
		}

#line 1584 "sgejsv.f"
		i__1 = *lwork - *n - *n * *n;
#line 1584 "sgejsv.f"
		sgesvj_("Upper", "U", "N", n, n, &work[*n + 1], n, &sva[1], n,
			 &u[u_offset], ldu, &work[*n + *n * *n + 1], &i__1, 
			info, (ftnlen)5, (ftnlen)1, (ftnlen)1);

#line 1587 "sgejsv.f"
		scalem = work[*n + *n * *n + 1];
#line 1588 "sgejsv.f"
		numrank = i_dnnt(&work[*n + *n * *n + 2]);
#line 1589 "sgejsv.f"
		i__1 = *n;
#line 1589 "sgejsv.f"
		for (p = 1; p <= i__1; ++p) {
#line 1590 "sgejsv.f"
		    scopy_(n, &work[*n + (p - 1) * *n + 1], &c__1, &u[p * 
			    u_dim1 + 1], &c__1);
#line 1591 "sgejsv.f"
		    sscal_(n, &sva[p], &work[*n + (p - 1) * *n + 1], &c__1);
#line 1592 "sgejsv.f"
/* L6970: */
#line 1592 "sgejsv.f"
		}

#line 1594 "sgejsv.f"
		strsm_("Left", "Upper", "NoTrans", "No UD", n, n, &c_b35, &a[
			a_offset], lda, &work[*n + 1], n, (ftnlen)4, (ftnlen)
			5, (ftnlen)7, (ftnlen)5);
#line 1596 "sgejsv.f"
		i__1 = *n;
#line 1596 "sgejsv.f"
		for (p = 1; p <= i__1; ++p) {
#line 1597 "sgejsv.f"
		    scopy_(n, &work[*n + p], n, &v[iwork[p] + v_dim1], ldv);
#line 1598 "sgejsv.f"
/* L6972: */
#line 1598 "sgejsv.f"
		}
#line 1599 "sgejsv.f"
		temp1 = sqrt((doublereal) (*n)) * epsln;
#line 1600 "sgejsv.f"
		i__1 = *n;
#line 1600 "sgejsv.f"
		for (p = 1; p <= i__1; ++p) {
#line 1601 "sgejsv.f"
		    xsc = 1. / snrm2_(n, &v[p * v_dim1 + 1], &c__1);
#line 1602 "sgejsv.f"
		    if (xsc < 1. - temp1 || xsc > temp1 + 1.) {
#line 1602 "sgejsv.f"
			sscal_(n, &xsc, &v[p * v_dim1 + 1], &c__1);
#line 1602 "sgejsv.f"
		    }
#line 1604 "sgejsv.f"
/* L6971: */
#line 1604 "sgejsv.f"
		}

/*           Assemble the left singular vector matrix U (M x N). */

#line 1608 "sgejsv.f"
		if (*n < *m) {
#line 1609 "sgejsv.f"
		    i__1 = *m - *n;
#line 1609 "sgejsv.f"
		    slaset_("A", &i__1, n, &c_b34, &c_b34, &u[*n + 1 + u_dim1]
			    , ldu, (ftnlen)1);
#line 1610 "sgejsv.f"
		    if (*n < n1) {
#line 1611 "sgejsv.f"
			i__1 = n1 - *n;
#line 1611 "sgejsv.f"
			slaset_("A", n, &i__1, &c_b34, &c_b34, &u[(*n + 1) * 
				u_dim1 + 1], ldu, (ftnlen)1);
#line 1612 "sgejsv.f"
			i__1 = *m - *n;
#line 1612 "sgejsv.f"
			i__2 = n1 - *n;
#line 1612 "sgejsv.f"
			slaset_("A", &i__1, &i__2, &c_b34, &c_b35, &u[*n + 1 
				+ (*n + 1) * u_dim1], ldu, (ftnlen)1);
#line 1613 "sgejsv.f"
		    }
#line 1614 "sgejsv.f"
		}
#line 1615 "sgejsv.f"
		i__1 = *lwork - *n;
#line 1615 "sgejsv.f"
		sormqr_("Left", "No Tr", m, &n1, n, &a[a_offset], lda, &work[
			1], &u[u_offset], ldu, &work[*n + 1], &i__1, &ierr, (
			ftnlen)4, (ftnlen)5);
#line 1617 "sgejsv.f"
		temp1 = sqrt((doublereal) (*m)) * epsln;
#line 1618 "sgejsv.f"
		i__1 = n1;
#line 1618 "sgejsv.f"
		for (p = 1; p <= i__1; ++p) {
#line 1619 "sgejsv.f"
		    xsc = 1. / snrm2_(m, &u[p * u_dim1 + 1], &c__1);
#line 1620 "sgejsv.f"
		    if (xsc < 1. - temp1 || xsc > temp1 + 1.) {
#line 1620 "sgejsv.f"
			sscal_(m, &xsc, &u[p * u_dim1 + 1], &c__1);
#line 1620 "sgejsv.f"
		    }
#line 1622 "sgejsv.f"
/* L6973: */
#line 1622 "sgejsv.f"
		}

#line 1624 "sgejsv.f"
		if (rowpiv) {
#line 1624 "sgejsv.f"
		    i__1 = *m - 1;
#line 1624 "sgejsv.f"
		    slaswp_(&n1, &u[u_offset], ldu, &c__1, &i__1, &iwork[(*n 
			    << 1) + 1], &c_n1);
#line 1624 "sgejsv.f"
		}

#line 1627 "sgejsv.f"
	    }

/*        end of the  >> almost orthogonal case <<  in the full SVD */

#line 1631 "sgejsv.f"
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

#line 1643 "sgejsv.f"
	    i__1 = nr;
#line 1643 "sgejsv.f"
	    for (p = 1; p <= i__1; ++p) {
#line 1644 "sgejsv.f"
		i__2 = *n - p + 1;
#line 1644 "sgejsv.f"
		scopy_(&i__2, &a[p + p * a_dim1], lda, &v[p + p * v_dim1], &
			c__1);
#line 1645 "sgejsv.f"
/* L7968: */
#line 1645 "sgejsv.f"
	    }

#line 1647 "sgejsv.f"
	    if (l2pert) {
#line 1648 "sgejsv.f"
		xsc = sqrt(small / epsln);
#line 1649 "sgejsv.f"
		i__1 = nr;
#line 1649 "sgejsv.f"
		for (q = 1; q <= i__1; ++q) {
#line 1650 "sgejsv.f"
		    temp1 = xsc * (d__1 = v[q + q * v_dim1], abs(d__1));
#line 1651 "sgejsv.f"
		    i__2 = *n;
#line 1651 "sgejsv.f"
		    for (p = 1; p <= i__2; ++p) {
#line 1652 "sgejsv.f"
			if (p > q && (d__1 = v[p + q * v_dim1], abs(d__1)) <= 
				temp1 || p < q) {
#line 1652 "sgejsv.f"
			    v[p + q * v_dim1] = d_sign(&temp1, &v[p + q * 
				    v_dim1]);
#line 1652 "sgejsv.f"
			}
#line 1655 "sgejsv.f"
			if (p < q) {
#line 1655 "sgejsv.f"
			    v[p + q * v_dim1] = -v[p + q * v_dim1];
#line 1655 "sgejsv.f"
			}
#line 1656 "sgejsv.f"
/* L5968: */
#line 1656 "sgejsv.f"
		    }
#line 1657 "sgejsv.f"
/* L5969: */
#line 1657 "sgejsv.f"
		}
#line 1658 "sgejsv.f"
	    } else {
#line 1659 "sgejsv.f"
		i__1 = nr - 1;
#line 1659 "sgejsv.f"
		i__2 = nr - 1;
#line 1659 "sgejsv.f"
		slaset_("U", &i__1, &i__2, &c_b34, &c_b34, &v[(v_dim1 << 1) + 
			1], ldv, (ftnlen)1);
#line 1660 "sgejsv.f"
	    }
#line 1662 "sgejsv.f"
	    i__1 = *lwork - (*n << 1);
#line 1662 "sgejsv.f"
	    sgeqrf_(n, &nr, &v[v_offset], ldv, &work[*n + 1], &work[(*n << 1) 
		    + 1], &i__1, &ierr);
#line 1664 "sgejsv.f"
	    slacpy_("L", n, &nr, &v[v_offset], ldv, &work[(*n << 1) + 1], n, (
		    ftnlen)1);

#line 1666 "sgejsv.f"
	    i__1 = nr;
#line 1666 "sgejsv.f"
	    for (p = 1; p <= i__1; ++p) {
#line 1667 "sgejsv.f"
		i__2 = nr - p + 1;
#line 1667 "sgejsv.f"
		scopy_(&i__2, &v[p + p * v_dim1], ldv, &u[p + p * u_dim1], &
			c__1);
#line 1668 "sgejsv.f"
/* L7969: */
#line 1668 "sgejsv.f"
	    }
#line 1670 "sgejsv.f"
	    if (l2pert) {
#line 1671 "sgejsv.f"
		xsc = sqrt(small / epsln);
#line 1672 "sgejsv.f"
		i__1 = nr;
#line 1672 "sgejsv.f"
		for (q = 2; q <= i__1; ++q) {
#line 1673 "sgejsv.f"
		    i__2 = q - 1;
#line 1673 "sgejsv.f"
		    for (p = 1; p <= i__2; ++p) {
/* Computing MIN */
#line 1674 "sgejsv.f"
			d__3 = (d__1 = u[p + p * u_dim1], abs(d__1)), d__4 = (
				d__2 = u[q + q * u_dim1], abs(d__2));
#line 1674 "sgejsv.f"
			temp1 = xsc * min(d__3,d__4);
#line 1675 "sgejsv.f"
			u[p + q * u_dim1] = -d_sign(&temp1, &u[q + p * u_dim1]
				);
#line 1676 "sgejsv.f"
/* L9971: */
#line 1676 "sgejsv.f"
		    }
#line 1677 "sgejsv.f"
/* L9970: */
#line 1677 "sgejsv.f"
		}
#line 1678 "sgejsv.f"
	    } else {
#line 1679 "sgejsv.f"
		i__1 = nr - 1;
#line 1679 "sgejsv.f"
		i__2 = nr - 1;
#line 1679 "sgejsv.f"
		slaset_("U", &i__1, &i__2, &c_b34, &c_b34, &u[(u_dim1 << 1) + 
			1], ldu, (ftnlen)1);
#line 1680 "sgejsv.f"
	    }
#line 1682 "sgejsv.f"
	    i__1 = *lwork - (*n << 1) - *n * nr;
#line 1682 "sgejsv.f"
	    sgesvj_("L", "U", "V", &nr, &nr, &u[u_offset], ldu, &sva[1], n, &
		    v[v_offset], ldv, &work[(*n << 1) + *n * nr + 1], &i__1, 
		    info, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1684 "sgejsv.f"
	    scalem = work[(*n << 1) + *n * nr + 1];
#line 1685 "sgejsv.f"
	    numrank = i_dnnt(&work[(*n << 1) + *n * nr + 2]);
#line 1687 "sgejsv.f"
	    if (nr < *n) {
#line 1688 "sgejsv.f"
		i__1 = *n - nr;
#line 1688 "sgejsv.f"
		slaset_("A", &i__1, &nr, &c_b34, &c_b34, &v[nr + 1 + v_dim1], 
			ldv, (ftnlen)1);
#line 1689 "sgejsv.f"
		i__1 = *n - nr;
#line 1689 "sgejsv.f"
		slaset_("A", &nr, &i__1, &c_b34, &c_b34, &v[(nr + 1) * v_dim1 
			+ 1], ldv, (ftnlen)1);
#line 1690 "sgejsv.f"
		i__1 = *n - nr;
#line 1690 "sgejsv.f"
		i__2 = *n - nr;
#line 1690 "sgejsv.f"
		slaset_("A", &i__1, &i__2, &c_b34, &c_b35, &v[nr + 1 + (nr + 
			1) * v_dim1], ldv, (ftnlen)1);
#line 1691 "sgejsv.f"
	    }
#line 1693 "sgejsv.f"
	    i__1 = *lwork - (*n << 1) - *n * nr - nr;
#line 1693 "sgejsv.f"
	    sormqr_("L", "N", n, n, &nr, &work[(*n << 1) + 1], n, &work[*n + 
		    1], &v[v_offset], ldv, &work[(*n << 1) + *n * nr + nr + 1]
		    , &i__1, &ierr, (ftnlen)1, (ftnlen)1);

/*           Permute the rows of V using the (column) permutation from the */
/*           first QRF. Also, scale the columns to make them unit in */
/*           Euclidean norm. This applies to all cases. */

#line 1700 "sgejsv.f"
	    temp1 = sqrt((doublereal) (*n)) * epsln;
#line 1701 "sgejsv.f"
	    i__1 = *n;
#line 1701 "sgejsv.f"
	    for (q = 1; q <= i__1; ++q) {
#line 1702 "sgejsv.f"
		i__2 = *n;
#line 1702 "sgejsv.f"
		for (p = 1; p <= i__2; ++p) {
#line 1703 "sgejsv.f"
		    work[(*n << 1) + *n * nr + nr + iwork[p]] = v[p + q * 
			    v_dim1];
#line 1704 "sgejsv.f"
/* L8972: */
#line 1704 "sgejsv.f"
		}
#line 1705 "sgejsv.f"
		i__2 = *n;
#line 1705 "sgejsv.f"
		for (p = 1; p <= i__2; ++p) {
#line 1706 "sgejsv.f"
		    v[p + q * v_dim1] = work[(*n << 1) + *n * nr + nr + p];
#line 1707 "sgejsv.f"
/* L8973: */
#line 1707 "sgejsv.f"
		}
#line 1708 "sgejsv.f"
		xsc = 1. / snrm2_(n, &v[q * v_dim1 + 1], &c__1);
#line 1709 "sgejsv.f"
		if (xsc < 1. - temp1 || xsc > temp1 + 1.) {
#line 1709 "sgejsv.f"
		    sscal_(n, &xsc, &v[q * v_dim1 + 1], &c__1);
#line 1709 "sgejsv.f"
		}
#line 1711 "sgejsv.f"
/* L7972: */
#line 1711 "sgejsv.f"
	    }

/*           At this moment, V contains the right singular vectors of A. */
/*           Next, assemble the left singular vector matrix U (M x N). */

#line 1716 "sgejsv.f"
	    if (nr < *m) {
#line 1717 "sgejsv.f"
		i__1 = *m - nr;
#line 1717 "sgejsv.f"
		slaset_("A", &i__1, &nr, &c_b34, &c_b34, &u[nr + 1 + u_dim1], 
			ldu, (ftnlen)1);
#line 1718 "sgejsv.f"
		if (nr < n1) {
#line 1719 "sgejsv.f"
		    i__1 = n1 - nr;
#line 1719 "sgejsv.f"
		    slaset_("A", &nr, &i__1, &c_b34, &c_b34, &u[(nr + 1) * 
			    u_dim1 + 1], ldu, (ftnlen)1);
#line 1720 "sgejsv.f"
		    i__1 = *m - nr;
#line 1720 "sgejsv.f"
		    i__2 = n1 - nr;
#line 1720 "sgejsv.f"
		    slaset_("A", &i__1, &i__2, &c_b34, &c_b35, &u[nr + 1 + (
			    nr + 1) * u_dim1], ldu, (ftnlen)1);
#line 1721 "sgejsv.f"
		}
#line 1722 "sgejsv.f"
	    }

#line 1724 "sgejsv.f"
	    i__1 = *lwork - *n;
#line 1724 "sgejsv.f"
	    sormqr_("Left", "No Tr", m, &n1, n, &a[a_offset], lda, &work[1], &
		    u[u_offset], ldu, &work[*n + 1], &i__1, &ierr, (ftnlen)4, 
		    (ftnlen)5);

#line 1727 "sgejsv.f"
	    if (rowpiv) {
#line 1727 "sgejsv.f"
		i__1 = *m - 1;
#line 1727 "sgejsv.f"
		slaswp_(&n1, &u[u_offset], ldu, &c__1, &i__1, &iwork[(*n << 1)
			 + 1], &c_n1);
#line 1727 "sgejsv.f"
	    }


#line 1731 "sgejsv.f"
	}
#line 1732 "sgejsv.f"
	if (transp) {
/*           .. swap U and V because the procedure worked on A^t */
#line 1734 "sgejsv.f"
	    i__1 = *n;
#line 1734 "sgejsv.f"
	    for (p = 1; p <= i__1; ++p) {
#line 1735 "sgejsv.f"
		sswap_(n, &u[p * u_dim1 + 1], &c__1, &v[p * v_dim1 + 1], &
			c__1);
#line 1736 "sgejsv.f"
/* L6974: */
#line 1736 "sgejsv.f"
	    }
#line 1737 "sgejsv.f"
	}

#line 1739 "sgejsv.f"
    }
/*     end of the full SVD */

/*     Undo scaling, if necessary (and possible) */

#line 1744 "sgejsv.f"
    if (uscal2 <= big / sva[1] * uscal1) {
#line 1745 "sgejsv.f"
	slascl_("G", &c__0, &c__0, &uscal1, &uscal2, &nr, &c__1, &sva[1], n, &
		ierr, (ftnlen)1);
#line 1746 "sgejsv.f"
	uscal1 = 1.;
#line 1747 "sgejsv.f"
	uscal2 = 1.;
#line 1748 "sgejsv.f"
    }

#line 1750 "sgejsv.f"
    if (nr < *n) {
#line 1751 "sgejsv.f"
	i__1 = *n;
#line 1751 "sgejsv.f"
	for (p = nr + 1; p <= i__1; ++p) {
#line 1752 "sgejsv.f"
	    sva[p] = 0.;
#line 1753 "sgejsv.f"
/* L3004: */
#line 1753 "sgejsv.f"
	}
#line 1754 "sgejsv.f"
    }

#line 1756 "sgejsv.f"
    work[1] = uscal2 * scalem;
#line 1757 "sgejsv.f"
    work[2] = uscal1;
#line 1758 "sgejsv.f"
    if (errest) {
#line 1758 "sgejsv.f"
	work[3] = sconda;
#line 1758 "sgejsv.f"
    }
#line 1759 "sgejsv.f"
    if (lsvec && rsvec) {
#line 1760 "sgejsv.f"
	work[4] = condr1;
#line 1761 "sgejsv.f"
	work[5] = condr2;
#line 1762 "sgejsv.f"
    }
#line 1763 "sgejsv.f"
    if (l2tran) {
#line 1764 "sgejsv.f"
	work[6] = entra;
#line 1765 "sgejsv.f"
	work[7] = entrat;
#line 1766 "sgejsv.f"
    }

#line 1768 "sgejsv.f"
    iwork[1] = nr;
#line 1769 "sgejsv.f"
    iwork[2] = numrank;
#line 1770 "sgejsv.f"
    iwork[3] = warning;

#line 1772 "sgejsv.f"
    return 0;
/*     .. */
/*     .. END OF SGEJSV */
/*     .. */
} /* sgejsv_ */

