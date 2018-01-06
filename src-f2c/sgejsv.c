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
/* > SGEJSV can sometimes compute tiny singular values and their singular vectors much */
/* > more accurately than other SVD routines, see below under Further Details. */
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

/* > \date November 2015 */

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

#line 530 "sgejsv.f"
    /* Parameter adjustments */
#line 530 "sgejsv.f"
    --sva;
#line 530 "sgejsv.f"
    a_dim1 = *lda;
#line 530 "sgejsv.f"
    a_offset = 1 + a_dim1;
#line 530 "sgejsv.f"
    a -= a_offset;
#line 530 "sgejsv.f"
    u_dim1 = *ldu;
#line 530 "sgejsv.f"
    u_offset = 1 + u_dim1;
#line 530 "sgejsv.f"
    u -= u_offset;
#line 530 "sgejsv.f"
    v_dim1 = *ldv;
#line 530 "sgejsv.f"
    v_offset = 1 + v_dim1;
#line 530 "sgejsv.f"
    v -= v_offset;
#line 530 "sgejsv.f"
    --work;
#line 530 "sgejsv.f"
    --iwork;
#line 530 "sgejsv.f"

#line 530 "sgejsv.f"
    /* Function Body */
#line 530 "sgejsv.f"
    lsvec = lsame_(jobu, "U", (ftnlen)1, (ftnlen)1) || lsame_(jobu, "F", (
	    ftnlen)1, (ftnlen)1);
#line 531 "sgejsv.f"
    jracc = lsame_(jobv, "J", (ftnlen)1, (ftnlen)1);
#line 532 "sgejsv.f"
    rsvec = lsame_(jobv, "V", (ftnlen)1, (ftnlen)1) || jracc;
#line 533 "sgejsv.f"
    rowpiv = lsame_(joba, "F", (ftnlen)1, (ftnlen)1) || lsame_(joba, "G", (
	    ftnlen)1, (ftnlen)1);
#line 534 "sgejsv.f"
    l2rank = lsame_(joba, "R", (ftnlen)1, (ftnlen)1);
#line 535 "sgejsv.f"
    l2aber = lsame_(joba, "A", (ftnlen)1, (ftnlen)1);
#line 536 "sgejsv.f"
    errest = lsame_(joba, "E", (ftnlen)1, (ftnlen)1) || lsame_(joba, "G", (
	    ftnlen)1, (ftnlen)1);
#line 537 "sgejsv.f"
    l2tran = lsame_(jobt, "T", (ftnlen)1, (ftnlen)1);
#line 538 "sgejsv.f"
    l2kill = lsame_(jobr, "R", (ftnlen)1, (ftnlen)1);
#line 539 "sgejsv.f"
    defr = lsame_(jobr, "N", (ftnlen)1, (ftnlen)1);
#line 540 "sgejsv.f"
    l2pert = lsame_(jobp, "P", (ftnlen)1, (ftnlen)1);

#line 542 "sgejsv.f"
    if (! (rowpiv || l2rank || l2aber || errest || lsame_(joba, "C", (ftnlen)
	    1, (ftnlen)1))) {
#line 544 "sgejsv.f"
	*info = -1;
#line 545 "sgejsv.f"
    } else if (! (lsvec || lsame_(jobu, "N", (ftnlen)1, (ftnlen)1) || lsame_(
	    jobu, "W", (ftnlen)1, (ftnlen)1))) {
#line 547 "sgejsv.f"
	*info = -2;
#line 548 "sgejsv.f"
    } else if (! (rsvec || lsame_(jobv, "N", (ftnlen)1, (ftnlen)1) || lsame_(
	    jobv, "W", (ftnlen)1, (ftnlen)1)) || jracc && ! lsvec) {
#line 550 "sgejsv.f"
	*info = -3;
#line 551 "sgejsv.f"
    } else if (! (l2kill || defr)) {
#line 552 "sgejsv.f"
	*info = -4;
#line 553 "sgejsv.f"
    } else if (! (l2tran || lsame_(jobt, "N", (ftnlen)1, (ftnlen)1))) {
#line 554 "sgejsv.f"
	*info = -5;
#line 555 "sgejsv.f"
    } else if (! (l2pert || lsame_(jobp, "N", (ftnlen)1, (ftnlen)1))) {
#line 556 "sgejsv.f"
	*info = -6;
#line 557 "sgejsv.f"
    } else if (*m < 0) {
#line 558 "sgejsv.f"
	*info = -7;
#line 559 "sgejsv.f"
    } else if (*n < 0 || *n > *m) {
#line 560 "sgejsv.f"
	*info = -8;
#line 561 "sgejsv.f"
    } else if (*lda < *m) {
#line 562 "sgejsv.f"
	*info = -10;
#line 563 "sgejsv.f"
    } else if (lsvec && *ldu < *m) {
#line 564 "sgejsv.f"
	*info = -13;
#line 565 "sgejsv.f"
    } else if (rsvec && *ldv < *n) {
#line 566 "sgejsv.f"
	*info = -14;
#line 567 "sgejsv.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 567 "sgejsv.f"
	i__1 = 7, i__2 = (*n << 2) + 1, i__1 = max(i__1,i__2), i__2 = (*m << 
		1) + *n;
/* Computing MAX */
#line 567 "sgejsv.f"
	i__3 = 7, i__4 = (*n << 2) + *n * *n, i__3 = max(i__3,i__4), i__4 = (*
		m << 1) + *n;
/* Computing MAX */
#line 567 "sgejsv.f"
	i__5 = 7, i__6 = (*m << 1) + *n, i__5 = max(i__5,i__6), i__6 = (*n << 
		2) + 1;
/* Computing MAX */
#line 567 "sgejsv.f"
	i__7 = 7, i__8 = (*m << 1) + *n, i__7 = max(i__7,i__8), i__8 = (*n << 
		2) + 1;
/* Computing MAX */
#line 567 "sgejsv.f"
	i__9 = (*m << 1) + *n, i__10 = *n * 6 + (*n << 1) * *n;
/* Computing MAX */
#line 567 "sgejsv.f"
	i__11 = (*m << 1) + *n, i__12 = (*n << 2) + *n * *n, i__11 = max(
		i__11,i__12), i__12 = (*n << 1) + *n * *n + 6;
#line 567 "sgejsv.f"
	if (! (lsvec || rsvec || errest) && *lwork < max(i__1,i__2) || ! (
		lsvec || rsvec) && errest && *lwork < max(i__3,i__4) || lsvec 
		&& ! rsvec && *lwork < max(i__5,i__6) || rsvec && ! lsvec && *
		lwork < max(i__7,i__8) || lsvec && rsvec && ! jracc && *lwork 
		< max(i__9,i__10) || lsvec && rsvec && jracc && *lwork < max(
		i__11,i__12)) {
#line 580 "sgejsv.f"
	    *info = -17;
#line 581 "sgejsv.f"
	} else {
/*        #:) */
#line 583 "sgejsv.f"
	    *info = 0;
#line 584 "sgejsv.f"
	}
#line 584 "sgejsv.f"
    }

#line 586 "sgejsv.f"
    if (*info != 0) {
/*       #:( */
#line 588 "sgejsv.f"
	i__1 = -(*info);
#line 588 "sgejsv.f"
	xerbla_("SGEJSV", &i__1, (ftnlen)6);
#line 589 "sgejsv.f"
	return 0;
#line 590 "sgejsv.f"
    }

/*     Quick return for void matrix (Y3K safe) */
/* #:) */
#line 594 "sgejsv.f"
    if (*m == 0 || *n == 0) {
#line 594 "sgejsv.f"
	return 0;
#line 594 "sgejsv.f"
    }

/*     Determine whether the matrix U should be M x N or M x M */

#line 598 "sgejsv.f"
    if (lsvec) {
#line 599 "sgejsv.f"
	n1 = *n;
#line 600 "sgejsv.f"
	if (lsame_(jobu, "F", (ftnlen)1, (ftnlen)1)) {
#line 600 "sgejsv.f"
	    n1 = *m;
#line 600 "sgejsv.f"
	}
#line 601 "sgejsv.f"
    }

/*     Set numerical parameters */

/* !    NOTE: Make sure SLAMCH() does not fail on the target architecture. */

#line 607 "sgejsv.f"
    epsln = slamch_("Epsilon", (ftnlen)7);
#line 608 "sgejsv.f"
    sfmin = slamch_("SafeMinimum", (ftnlen)11);
#line 609 "sgejsv.f"
    small = sfmin / epsln;
#line 610 "sgejsv.f"
    big = slamch_("O", (ftnlen)1);
/*     BIG   = ONE / SFMIN */

/*     Initialize SVA(1:N) = diag( ||A e_i||_2 )_1^N */

/* (!)  If necessary, scale SVA() to protect the largest norm from */
/*     overflow. It is possible that this scaling pushes the smallest */
/*     column norm left from the underflow threshold (extreme case). */

#line 619 "sgejsv.f"
    scalem = 1. / sqrt((doublereal) (*m) * (doublereal) (*n));
#line 620 "sgejsv.f"
    noscal = TRUE_;
#line 621 "sgejsv.f"
    goscal = TRUE_;
#line 622 "sgejsv.f"
    i__1 = *n;
#line 622 "sgejsv.f"
    for (p = 1; p <= i__1; ++p) {
#line 623 "sgejsv.f"
	aapp = 0.;
#line 624 "sgejsv.f"
	aaqq = 1.;
#line 625 "sgejsv.f"
	slassq_(m, &a[p * a_dim1 + 1], &c__1, &aapp, &aaqq);
#line 626 "sgejsv.f"
	if (aapp > big) {
#line 627 "sgejsv.f"
	    *info = -9;
#line 628 "sgejsv.f"
	    i__2 = -(*info);
#line 628 "sgejsv.f"
	    xerbla_("SGEJSV", &i__2, (ftnlen)6);
#line 629 "sgejsv.f"
	    return 0;
#line 630 "sgejsv.f"
	}
#line 631 "sgejsv.f"
	aaqq = sqrt(aaqq);
#line 632 "sgejsv.f"
	if (aapp < big / aaqq && noscal) {
#line 633 "sgejsv.f"
	    sva[p] = aapp * aaqq;
#line 634 "sgejsv.f"
	} else {
#line 635 "sgejsv.f"
	    noscal = FALSE_;
#line 636 "sgejsv.f"
	    sva[p] = aapp * (aaqq * scalem);
#line 637 "sgejsv.f"
	    if (goscal) {
#line 638 "sgejsv.f"
		goscal = FALSE_;
#line 639 "sgejsv.f"
		i__2 = p - 1;
#line 639 "sgejsv.f"
		sscal_(&i__2, &scalem, &sva[1], &c__1);
#line 640 "sgejsv.f"
	    }
#line 641 "sgejsv.f"
	}
#line 642 "sgejsv.f"
/* L1874: */
#line 642 "sgejsv.f"
    }

#line 644 "sgejsv.f"
    if (noscal) {
#line 644 "sgejsv.f"
	scalem = 1.;
#line 644 "sgejsv.f"
    }

#line 646 "sgejsv.f"
    aapp = 0.;
#line 647 "sgejsv.f"
    aaqq = big;
#line 648 "sgejsv.f"
    i__1 = *n;
#line 648 "sgejsv.f"
    for (p = 1; p <= i__1; ++p) {
/* Computing MAX */
#line 649 "sgejsv.f"
	d__1 = aapp, d__2 = sva[p];
#line 649 "sgejsv.f"
	aapp = max(d__1,d__2);
#line 650 "sgejsv.f"
	if (sva[p] != 0.) {
/* Computing MIN */
#line 650 "sgejsv.f"
	    d__1 = aaqq, d__2 = sva[p];
#line 650 "sgejsv.f"
	    aaqq = min(d__1,d__2);
#line 650 "sgejsv.f"
	}
#line 651 "sgejsv.f"
/* L4781: */
#line 651 "sgejsv.f"
    }

/*     Quick return for zero M x N matrix */
/* #:) */
#line 655 "sgejsv.f"
    if (aapp == 0.) {
#line 656 "sgejsv.f"
	if (lsvec) {
#line 656 "sgejsv.f"
	    slaset_("G", m, &n1, &c_b34, &c_b35, &u[u_offset], ldu, (ftnlen)1)
		    ;
#line 656 "sgejsv.f"
	}
#line 657 "sgejsv.f"
	if (rsvec) {
#line 657 "sgejsv.f"
	    slaset_("G", n, n, &c_b34, &c_b35, &v[v_offset], ldv, (ftnlen)1);
#line 657 "sgejsv.f"
	}
#line 658 "sgejsv.f"
	work[1] = 1.;
#line 659 "sgejsv.f"
	work[2] = 1.;
#line 660 "sgejsv.f"
	if (errest) {
#line 660 "sgejsv.f"
	    work[3] = 1.;
#line 660 "sgejsv.f"
	}
#line 661 "sgejsv.f"
	if (lsvec && rsvec) {
#line 662 "sgejsv.f"
	    work[4] = 1.;
#line 663 "sgejsv.f"
	    work[5] = 1.;
#line 664 "sgejsv.f"
	}
#line 665 "sgejsv.f"
	if (l2tran) {
#line 666 "sgejsv.f"
	    work[6] = 0.;
#line 667 "sgejsv.f"
	    work[7] = 0.;
#line 668 "sgejsv.f"
	}
#line 669 "sgejsv.f"
	iwork[1] = 0;
#line 670 "sgejsv.f"
	iwork[2] = 0;
#line 671 "sgejsv.f"
	iwork[3] = 0;
#line 672 "sgejsv.f"
	return 0;
#line 673 "sgejsv.f"
    }

/*     Issue warning if denormalized column norms detected. Override the */
/*     high relative accuracy request. Issue licence to kill columns */
/*     (set them to zero) whose norm is less than sigma_max / BIG (roughly). */
/* #:( */
#line 679 "sgejsv.f"
    warning = 0;
#line 680 "sgejsv.f"
    if (aaqq <= sfmin) {
#line 681 "sgejsv.f"
	l2rank = TRUE_;
#line 682 "sgejsv.f"
	l2kill = TRUE_;
#line 683 "sgejsv.f"
	warning = 1;
#line 684 "sgejsv.f"
    }

/*     Quick return for one-column matrix */
/* #:) */
#line 688 "sgejsv.f"
    if (*n == 1) {

#line 690 "sgejsv.f"
	if (lsvec) {
#line 691 "sgejsv.f"
	    slascl_("G", &c__0, &c__0, &sva[1], &scalem, m, &c__1, &a[a_dim1 
		    + 1], lda, &ierr, (ftnlen)1);
#line 692 "sgejsv.f"
	    slacpy_("A", m, &c__1, &a[a_offset], lda, &u[u_offset], ldu, (
		    ftnlen)1);
/*           computing all M left singular vectors of the M x 1 matrix */
#line 694 "sgejsv.f"
	    if (n1 != *n) {
#line 695 "sgejsv.f"
		i__1 = *lwork - *n;
#line 695 "sgejsv.f"
		sgeqrf_(m, n, &u[u_offset], ldu, &work[1], &work[*n + 1], &
			i__1, &ierr);
#line 696 "sgejsv.f"
		i__1 = *lwork - *n;
#line 696 "sgejsv.f"
		sorgqr_(m, &n1, &c__1, &u[u_offset], ldu, &work[1], &work[*n 
			+ 1], &i__1, &ierr);
#line 697 "sgejsv.f"
		scopy_(m, &a[a_dim1 + 1], &c__1, &u[u_dim1 + 1], &c__1);
#line 698 "sgejsv.f"
	    }
#line 699 "sgejsv.f"
	}
#line 700 "sgejsv.f"
	if (rsvec) {
#line 701 "sgejsv.f"
	    v[v_dim1 + 1] = 1.;
#line 702 "sgejsv.f"
	}
#line 703 "sgejsv.f"
	if (sva[1] < big * scalem) {
#line 704 "sgejsv.f"
	    sva[1] /= scalem;
#line 705 "sgejsv.f"
	    scalem = 1.;
#line 706 "sgejsv.f"
	}
#line 707 "sgejsv.f"
	work[1] = 1. / scalem;
#line 708 "sgejsv.f"
	work[2] = 1.;
#line 709 "sgejsv.f"
	if (sva[1] != 0.) {
#line 710 "sgejsv.f"
	    iwork[1] = 1;
#line 711 "sgejsv.f"
	    if (sva[1] / scalem >= sfmin) {
#line 712 "sgejsv.f"
		iwork[2] = 1;
#line 713 "sgejsv.f"
	    } else {
#line 714 "sgejsv.f"
		iwork[2] = 0;
#line 715 "sgejsv.f"
	    }
#line 716 "sgejsv.f"
	} else {
#line 717 "sgejsv.f"
	    iwork[1] = 0;
#line 718 "sgejsv.f"
	    iwork[2] = 0;
#line 719 "sgejsv.f"
	}
#line 720 "sgejsv.f"
	if (errest) {
#line 720 "sgejsv.f"
	    work[3] = 1.;
#line 720 "sgejsv.f"
	}
#line 721 "sgejsv.f"
	if (lsvec && rsvec) {
#line 722 "sgejsv.f"
	    work[4] = 1.;
#line 723 "sgejsv.f"
	    work[5] = 1.;
#line 724 "sgejsv.f"
	}
#line 725 "sgejsv.f"
	if (l2tran) {
#line 726 "sgejsv.f"
	    work[6] = 0.;
#line 727 "sgejsv.f"
	    work[7] = 0.;
#line 728 "sgejsv.f"
	}
#line 729 "sgejsv.f"
	return 0;

#line 731 "sgejsv.f"
    }

#line 733 "sgejsv.f"
    transp = FALSE_;
#line 734 "sgejsv.f"
    l2tran = l2tran && *m == *n;

#line 736 "sgejsv.f"
    aatmax = -1.;
#line 737 "sgejsv.f"
    aatmin = big;
#line 738 "sgejsv.f"
    if (rowpiv || l2tran) {

/*     Compute the row norms, needed to determine row pivoting sequence */
/*     (in the case of heavily row weighted A, row pivoting is strongly */
/*     advised) and to collect information needed to compare the */
/*     structures of A * A^t and A^t * A (in the case L2TRAN.EQ..TRUE.). */

#line 745 "sgejsv.f"
	if (l2tran) {
#line 746 "sgejsv.f"
	    i__1 = *m;
#line 746 "sgejsv.f"
	    for (p = 1; p <= i__1; ++p) {
#line 747 "sgejsv.f"
		xsc = 0.;
#line 748 "sgejsv.f"
		temp1 = 1.;
#line 749 "sgejsv.f"
		slassq_(n, &a[p + a_dim1], lda, &xsc, &temp1);
/*              SLASSQ gets both the ell_2 and the ell_infinity norm */
/*              in one pass through the vector */
#line 752 "sgejsv.f"
		work[*m + *n + p] = xsc * scalem;
#line 753 "sgejsv.f"
		work[*n + p] = xsc * (scalem * sqrt(temp1));
/* Computing MAX */
#line 754 "sgejsv.f"
		d__1 = aatmax, d__2 = work[*n + p];
#line 754 "sgejsv.f"
		aatmax = max(d__1,d__2);
#line 755 "sgejsv.f"
		if (work[*n + p] != 0.) {
/* Computing MIN */
#line 755 "sgejsv.f"
		    d__1 = aatmin, d__2 = work[*n + p];
#line 755 "sgejsv.f"
		    aatmin = min(d__1,d__2);
#line 755 "sgejsv.f"
		}
#line 756 "sgejsv.f"
/* L1950: */
#line 756 "sgejsv.f"
	    }
#line 757 "sgejsv.f"
	} else {
#line 758 "sgejsv.f"
	    i__1 = *m;
#line 758 "sgejsv.f"
	    for (p = 1; p <= i__1; ++p) {
#line 759 "sgejsv.f"
		work[*m + *n + p] = scalem * (d__1 = a[p + isamax_(n, &a[p + 
			a_dim1], lda) * a_dim1], abs(d__1));
/* Computing MAX */
#line 760 "sgejsv.f"
		d__1 = aatmax, d__2 = work[*m + *n + p];
#line 760 "sgejsv.f"
		aatmax = max(d__1,d__2);
/* Computing MIN */
#line 761 "sgejsv.f"
		d__1 = aatmin, d__2 = work[*m + *n + p];
#line 761 "sgejsv.f"
		aatmin = min(d__1,d__2);
#line 762 "sgejsv.f"
/* L1904: */
#line 762 "sgejsv.f"
	    }
#line 763 "sgejsv.f"
	}

#line 765 "sgejsv.f"
    }

/*     For square matrix A try to determine whether A^t  would be  better */
/*     input for the preconditioned Jacobi SVD, with faster convergence. */
/*     The decision is based on an O(N) function of the vector of column */
/*     and row norms of A, based on the Shannon entropy. This should give */
/*     the right choice in most cases when the difference actually matters. */
/*     It may fail and pick the slower converging side. */

#line 774 "sgejsv.f"
    entra = 0.;
#line 775 "sgejsv.f"
    entrat = 0.;
#line 776 "sgejsv.f"
    if (l2tran) {

#line 778 "sgejsv.f"
	xsc = 0.;
#line 779 "sgejsv.f"
	temp1 = 1.;
#line 780 "sgejsv.f"
	slassq_(n, &sva[1], &c__1, &xsc, &temp1);
#line 781 "sgejsv.f"
	temp1 = 1. / temp1;

#line 783 "sgejsv.f"
	entra = 0.;
#line 784 "sgejsv.f"
	i__1 = *n;
#line 784 "sgejsv.f"
	for (p = 1; p <= i__1; ++p) {
/* Computing 2nd power */
#line 785 "sgejsv.f"
	    d__1 = sva[p] / xsc;
#line 785 "sgejsv.f"
	    big1 = d__1 * d__1 * temp1;
#line 786 "sgejsv.f"
	    if (big1 != 0.) {
#line 786 "sgejsv.f"
		entra += big1 * log(big1);
#line 786 "sgejsv.f"
	    }
#line 787 "sgejsv.f"
/* L1113: */
#line 787 "sgejsv.f"
	}
#line 788 "sgejsv.f"
	entra = -entra / log((doublereal) (*n));

/*        Now, SVA().^2/Trace(A^t * A) is a point in the probability simplex. */
/*        It is derived from the diagonal of  A^t * A.  Do the same with the */
/*        diagonal of A * A^t, compute the entropy of the corresponding */
/*        probability distribution. Note that A * A^t and A^t * A have the */
/*        same trace. */

#line 796 "sgejsv.f"
	entrat = 0.;
#line 797 "sgejsv.f"
	i__1 = *n + *m;
#line 797 "sgejsv.f"
	for (p = *n + 1; p <= i__1; ++p) {
/* Computing 2nd power */
#line 798 "sgejsv.f"
	    d__1 = work[p] / xsc;
#line 798 "sgejsv.f"
	    big1 = d__1 * d__1 * temp1;
#line 799 "sgejsv.f"
	    if (big1 != 0.) {
#line 799 "sgejsv.f"
		entrat += big1 * log(big1);
#line 799 "sgejsv.f"
	    }
#line 800 "sgejsv.f"
/* L1114: */
#line 800 "sgejsv.f"
	}
#line 801 "sgejsv.f"
	entrat = -entrat / log((doublereal) (*m));

/*        Analyze the entropies and decide A or A^t. Smaller entropy */
/*        usually means better input for the algorithm. */

#line 806 "sgejsv.f"
	transp = entrat < entra;

/*        If A^t is better than A, transpose A. */

#line 810 "sgejsv.f"
	if (transp) {
/*           In an optimal implementation, this trivial transpose */
/*           should be replaced with faster transpose. */
#line 813 "sgejsv.f"
	    i__1 = *n - 1;
#line 813 "sgejsv.f"
	    for (p = 1; p <= i__1; ++p) {
#line 814 "sgejsv.f"
		i__2 = *n;
#line 814 "sgejsv.f"
		for (q = p + 1; q <= i__2; ++q) {
#line 815 "sgejsv.f"
		    temp1 = a[q + p * a_dim1];
#line 816 "sgejsv.f"
		    a[q + p * a_dim1] = a[p + q * a_dim1];
#line 817 "sgejsv.f"
		    a[p + q * a_dim1] = temp1;
#line 818 "sgejsv.f"
/* L1116: */
#line 818 "sgejsv.f"
		}
#line 819 "sgejsv.f"
/* L1115: */
#line 819 "sgejsv.f"
	    }
#line 820 "sgejsv.f"
	    i__1 = *n;
#line 820 "sgejsv.f"
	    for (p = 1; p <= i__1; ++p) {
#line 821 "sgejsv.f"
		work[*m + *n + p] = sva[p];
#line 822 "sgejsv.f"
		sva[p] = work[*n + p];
#line 823 "sgejsv.f"
/* L1117: */
#line 823 "sgejsv.f"
	    }
#line 824 "sgejsv.f"
	    temp1 = aapp;
#line 825 "sgejsv.f"
	    aapp = aatmax;
#line 826 "sgejsv.f"
	    aatmax = temp1;
#line 827 "sgejsv.f"
	    temp1 = aaqq;
#line 828 "sgejsv.f"
	    aaqq = aatmin;
#line 829 "sgejsv.f"
	    aatmin = temp1;
#line 830 "sgejsv.f"
	    kill = lsvec;
#line 831 "sgejsv.f"
	    lsvec = rsvec;
#line 832 "sgejsv.f"
	    rsvec = kill;
#line 833 "sgejsv.f"
	    if (lsvec) {
#line 833 "sgejsv.f"
		n1 = *n;
#line 833 "sgejsv.f"
	    }

#line 835 "sgejsv.f"
	    rowpiv = TRUE_;
#line 836 "sgejsv.f"
	}

#line 838 "sgejsv.f"
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

#line 851 "sgejsv.f"
    big1 = sqrt(big);
#line 852 "sgejsv.f"
    temp1 = sqrt(big / (doublereal) (*n));

#line 854 "sgejsv.f"
    slascl_("G", &c__0, &c__0, &aapp, &temp1, n, &c__1, &sva[1], n, &ierr, (
	    ftnlen)1);
#line 855 "sgejsv.f"
    if (aaqq > aapp * sfmin) {
#line 856 "sgejsv.f"
	aaqq = aaqq / aapp * temp1;
#line 857 "sgejsv.f"
    } else {
#line 858 "sgejsv.f"
	aaqq = aaqq * temp1 / aapp;
#line 859 "sgejsv.f"
    }
#line 860 "sgejsv.f"
    temp1 *= scalem;
#line 861 "sgejsv.f"
    slascl_("G", &c__0, &c__0, &aapp, &temp1, m, n, &a[a_offset], lda, &ierr, 
	    (ftnlen)1);

/*     To undo scaling at the end of this procedure, multiply the */
/*     computed singular values with USCAL2 / USCAL1. */

#line 866 "sgejsv.f"
    uscal1 = temp1;
#line 867 "sgejsv.f"
    uscal2 = aapp;

#line 869 "sgejsv.f"
    if (l2kill) {
/*        L2KILL enforces computation of nonzero singular values in */
/*        the restricted range of condition number of the initial A, */
/*        sigma_max(A) / sigma_min(A) approx. SQRT(BIG)/SQRT(SFMIN). */
#line 873 "sgejsv.f"
	xsc = sqrt(sfmin);
#line 874 "sgejsv.f"
    } else {
#line 875 "sgejsv.f"
	xsc = small;

/*        Now, if the condition number of A is too big, */
/*        sigma_max(A) / sigma_min(A) .GT. SQRT(BIG/N) * EPSLN / SFMIN, */
/*        as a precaution measure, the full SVD is computed using SGESVJ */
/*        with accumulated Jacobi rotations. This provides numerically */
/*        more robust computation, at the cost of slightly increased run */
/*        time. Depending on the concrete implementation of BLAS and LAPACK */
/*        (i.e. how they behave in presence of extreme ill-conditioning) the */
/*        implementor may decide to remove this switch. */
#line 885 "sgejsv.f"
	if (aaqq < sqrt(sfmin) && lsvec && rsvec) {
#line 886 "sgejsv.f"
	    jracc = TRUE_;
#line 887 "sgejsv.f"
	}

#line 889 "sgejsv.f"
    }
#line 890 "sgejsv.f"
    if (aaqq < xsc) {
#line 891 "sgejsv.f"
	i__1 = *n;
#line 891 "sgejsv.f"
	for (p = 1; p <= i__1; ++p) {
#line 892 "sgejsv.f"
	    if (sva[p] < xsc) {
#line 893 "sgejsv.f"
		slaset_("A", m, &c__1, &c_b34, &c_b34, &a[p * a_dim1 + 1], 
			lda, (ftnlen)1);
#line 894 "sgejsv.f"
		sva[p] = 0.;
#line 895 "sgejsv.f"
	    }
#line 896 "sgejsv.f"
/* L700: */
#line 896 "sgejsv.f"
	}
#line 897 "sgejsv.f"
    }

/*     Preconditioning using QR factorization with pivoting */

#line 901 "sgejsv.f"
    if (rowpiv) {
/*        Optional row permutation (Bjoerck row pivoting): */
/*        A result by Cox and Higham shows that the Bjoerck's */
/*        row pivoting combined with standard column pivoting */
/*        has similar effect as Powell-Reid complete pivoting. */
/*        The ell-infinity norms of A are made nonincreasing. */
#line 907 "sgejsv.f"
	i__1 = *m - 1;
#line 907 "sgejsv.f"
	for (p = 1; p <= i__1; ++p) {
#line 908 "sgejsv.f"
	    i__2 = *m - p + 1;
#line 908 "sgejsv.f"
	    q = isamax_(&i__2, &work[*m + *n + p], &c__1) + p - 1;
#line 909 "sgejsv.f"
	    iwork[(*n << 1) + p] = q;
#line 910 "sgejsv.f"
	    if (p != q) {
#line 911 "sgejsv.f"
		temp1 = work[*m + *n + p];
#line 912 "sgejsv.f"
		work[*m + *n + p] = work[*m + *n + q];
#line 913 "sgejsv.f"
		work[*m + *n + q] = temp1;
#line 914 "sgejsv.f"
	    }
#line 915 "sgejsv.f"
/* L1952: */
#line 915 "sgejsv.f"
	}
#line 916 "sgejsv.f"
	i__1 = *m - 1;
#line 916 "sgejsv.f"
	slaswp_(n, &a[a_offset], lda, &c__1, &i__1, &iwork[(*n << 1) + 1], &
		c__1);
#line 917 "sgejsv.f"
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
#line 934 "sgejsv.f"
    i__1 = *n;
#line 934 "sgejsv.f"
    for (p = 1; p <= i__1; ++p) {
/*        .. all columns are free columns */
#line 936 "sgejsv.f"
	iwork[p] = 0;
#line 937 "sgejsv.f"
/* L1963: */
#line 937 "sgejsv.f"
    }
#line 938 "sgejsv.f"
    i__1 = *lwork - *n;
#line 938 "sgejsv.f"
    sgeqp3_(m, n, &a[a_offset], lda, &iwork[1], &work[1], &work[*n + 1], &
	    i__1, &ierr);

/*     The upper triangular matrix R1 from the first QRF is inspected for */
/*     rank deficiency and possibilities for deflation, or possible */
/*     ill-conditioning. Depending on the user specified flag L2RANK, */
/*     the procedure explores possibilities to reduce the numerical */
/*     rank by inspecting the computed upper triangular factor. If */
/*     L2RANK or L2ABER are up, then SGEJSV will compute the SVD of */
/*     A + dA, where ||dA|| <= f(M,N)*EPSLN. */

#line 948 "sgejsv.f"
    nr = 1;
#line 949 "sgejsv.f"
    if (l2aber) {
/*        Standard absolute error bound suffices. All sigma_i with */
/*        sigma_i < N*EPSLN*||A|| are flushed to zero. This is an */
/*        agressive enforcement of lower numerical rank by introducing a */
/*        backward error of the order of N*EPSLN*||A||. */
#line 954 "sgejsv.f"
	temp1 = sqrt((doublereal) (*n)) * epsln;
#line 955 "sgejsv.f"
	i__1 = *n;
#line 955 "sgejsv.f"
	for (p = 2; p <= i__1; ++p) {
#line 956 "sgejsv.f"
	    if ((d__2 = a[p + p * a_dim1], abs(d__2)) >= temp1 * (d__1 = a[
		    a_dim1 + 1], abs(d__1))) {
#line 957 "sgejsv.f"
		++nr;
#line 958 "sgejsv.f"
	    } else {
#line 959 "sgejsv.f"
		goto L3002;
#line 960 "sgejsv.f"
	    }
#line 961 "sgejsv.f"
/* L3001: */
#line 961 "sgejsv.f"
	}
#line 962 "sgejsv.f"
L3002:
#line 963 "sgejsv.f"
	;
#line 963 "sgejsv.f"
    } else if (l2rank) {
/*        .. similarly as above, only slightly more gentle (less agressive). */
/*        Sudden drop on the diagonal of R1 is used as the criterion for */
/*        close-to-rank-defficient. */
#line 967 "sgejsv.f"
	temp1 = sqrt(sfmin);
#line 968 "sgejsv.f"
	i__1 = *n;
#line 968 "sgejsv.f"
	for (p = 2; p <= i__1; ++p) {
#line 969 "sgejsv.f"
	    if ((d__2 = a[p + p * a_dim1], abs(d__2)) < epsln * (d__1 = a[p - 
		    1 + (p - 1) * a_dim1], abs(d__1)) || (d__3 = a[p + p * 
		    a_dim1], abs(d__3)) < small || l2kill && (d__4 = a[p + p *
		     a_dim1], abs(d__4)) < temp1) {
#line 969 "sgejsv.f"
		goto L3402;
#line 969 "sgejsv.f"
	    }
#line 972 "sgejsv.f"
	    ++nr;
#line 973 "sgejsv.f"
/* L3401: */
#line 973 "sgejsv.f"
	}
#line 974 "sgejsv.f"
L3402:

#line 976 "sgejsv.f"
	;
#line 976 "sgejsv.f"
    } else {
/*        The goal is high relative accuracy. However, if the matrix */
/*        has high scaled condition number the relative accuracy is in */
/*        general not feasible. Later on, a condition number estimator */
/*        will be deployed to estimate the scaled condition number. */
/*        Here we just remove the underflowed part of the triangular */
/*        factor. This prevents the situation in which the code is */
/*        working hard to get the accuracy not warranted by the data. */
#line 984 "sgejsv.f"
	temp1 = sqrt(sfmin);
#line 985 "sgejsv.f"
	i__1 = *n;
#line 985 "sgejsv.f"
	for (p = 2; p <= i__1; ++p) {
#line 986 "sgejsv.f"
	    if ((d__1 = a[p + p * a_dim1], abs(d__1)) < small || l2kill && (
		    d__2 = a[p + p * a_dim1], abs(d__2)) < temp1) {
#line 986 "sgejsv.f"
		goto L3302;
#line 986 "sgejsv.f"
	    }
#line 988 "sgejsv.f"
	    ++nr;
#line 989 "sgejsv.f"
/* L3301: */
#line 989 "sgejsv.f"
	}
#line 990 "sgejsv.f"
L3302:

#line 992 "sgejsv.f"
	;
#line 992 "sgejsv.f"
    }

#line 994 "sgejsv.f"
    almort = FALSE_;
#line 995 "sgejsv.f"
    if (nr == *n) {
#line 996 "sgejsv.f"
	maxprj = 1.;
#line 997 "sgejsv.f"
	i__1 = *n;
#line 997 "sgejsv.f"
	for (p = 2; p <= i__1; ++p) {
#line 998 "sgejsv.f"
	    temp1 = (d__1 = a[p + p * a_dim1], abs(d__1)) / sva[iwork[p]];
#line 999 "sgejsv.f"
	    maxprj = min(maxprj,temp1);
#line 1000 "sgejsv.f"
/* L3051: */
#line 1000 "sgejsv.f"
	}
/* Computing 2nd power */
#line 1001 "sgejsv.f"
	d__1 = maxprj;
#line 1001 "sgejsv.f"
	if (d__1 * d__1 >= 1. - (doublereal) (*n) * epsln) {
#line 1001 "sgejsv.f"
	    almort = TRUE_;
#line 1001 "sgejsv.f"
	}
#line 1002 "sgejsv.f"
    }


#line 1005 "sgejsv.f"
    sconda = -1.;
#line 1006 "sgejsv.f"
    condr1 = -1.;
#line 1007 "sgejsv.f"
    condr2 = -1.;

#line 1009 "sgejsv.f"
    if (errest) {
#line 1010 "sgejsv.f"
	if (*n == nr) {
#line 1011 "sgejsv.f"
	    if (rsvec) {
/*              .. V is available as workspace */
#line 1013 "sgejsv.f"
		slacpy_("U", n, n, &a[a_offset], lda, &v[v_offset], ldv, (
			ftnlen)1);
#line 1014 "sgejsv.f"
		i__1 = *n;
#line 1014 "sgejsv.f"
		for (p = 1; p <= i__1; ++p) {
#line 1015 "sgejsv.f"
		    temp1 = sva[iwork[p]];
#line 1016 "sgejsv.f"
		    d__1 = 1. / temp1;
#line 1016 "sgejsv.f"
		    sscal_(&p, &d__1, &v[p * v_dim1 + 1], &c__1);
#line 1017 "sgejsv.f"
/* L3053: */
#line 1017 "sgejsv.f"
		}
#line 1018 "sgejsv.f"
		spocon_("U", n, &v[v_offset], ldv, &c_b35, &temp1, &work[*n + 
			1], &iwork[(*n << 1) + *m + 1], &ierr, (ftnlen)1);
#line 1020 "sgejsv.f"
	    } else if (lsvec) {
/*              .. U is available as workspace */
#line 1022 "sgejsv.f"
		slacpy_("U", n, n, &a[a_offset], lda, &u[u_offset], ldu, (
			ftnlen)1);
#line 1023 "sgejsv.f"
		i__1 = *n;
#line 1023 "sgejsv.f"
		for (p = 1; p <= i__1; ++p) {
#line 1024 "sgejsv.f"
		    temp1 = sva[iwork[p]];
#line 1025 "sgejsv.f"
		    d__1 = 1. / temp1;
#line 1025 "sgejsv.f"
		    sscal_(&p, &d__1, &u[p * u_dim1 + 1], &c__1);
#line 1026 "sgejsv.f"
/* L3054: */
#line 1026 "sgejsv.f"
		}
#line 1027 "sgejsv.f"
		spocon_("U", n, &u[u_offset], ldu, &c_b35, &temp1, &work[*n + 
			1], &iwork[(*n << 1) + *m + 1], &ierr, (ftnlen)1);
#line 1029 "sgejsv.f"
	    } else {
#line 1030 "sgejsv.f"
		slacpy_("U", n, n, &a[a_offset], lda, &work[*n + 1], n, (
			ftnlen)1);
#line 1031 "sgejsv.f"
		i__1 = *n;
#line 1031 "sgejsv.f"
		for (p = 1; p <= i__1; ++p) {
#line 1032 "sgejsv.f"
		    temp1 = sva[iwork[p]];
#line 1033 "sgejsv.f"
		    d__1 = 1. / temp1;
#line 1033 "sgejsv.f"
		    sscal_(&p, &d__1, &work[*n + (p - 1) * *n + 1], &c__1);
#line 1034 "sgejsv.f"
/* L3052: */
#line 1034 "sgejsv.f"
		}
/*           .. the columns of R are scaled to have unit Euclidean lengths. */
#line 1036 "sgejsv.f"
		spocon_("U", n, &work[*n + 1], n, &c_b35, &temp1, &work[*n + *
			n * *n + 1], &iwork[(*n << 1) + *m + 1], &ierr, (
			ftnlen)1);
#line 1038 "sgejsv.f"
	    }
#line 1039 "sgejsv.f"
	    sconda = 1. / sqrt(temp1);
/*           SCONDA is an estimate of SQRT(||(R^t * R)^(-1)||_1). */
/*           N^(-1/4) * SCONDA <= ||R^(-1)||_2 <= N^(1/4) * SCONDA */
#line 1042 "sgejsv.f"
	} else {
#line 1043 "sgejsv.f"
	    sconda = -1.;
#line 1044 "sgejsv.f"
	}
#line 1045 "sgejsv.f"
    }

#line 1047 "sgejsv.f"
    l2pert = l2pert && (d__1 = a[a_dim1 + 1] / a[nr + nr * a_dim1], abs(d__1))
	     > sqrt(big1);
/*     If there is no violent scaling, artificial perturbation is not needed. */

/*     Phase 3: */

#line 1052 "sgejsv.f"
    if (! (rsvec || lsvec)) {

/*         Singular Values only */

/*         .. transpose A(1:NR,1:N) */
/* Computing MIN */
#line 1057 "sgejsv.f"
	i__2 = *n - 1;
#line 1057 "sgejsv.f"
	i__1 = min(i__2,nr);
#line 1057 "sgejsv.f"
	for (p = 1; p <= i__1; ++p) {
#line 1058 "sgejsv.f"
	    i__2 = *n - p;
#line 1058 "sgejsv.f"
	    scopy_(&i__2, &a[p + (p + 1) * a_dim1], lda, &a[p + 1 + p * 
		    a_dim1], &c__1);
#line 1059 "sgejsv.f"
/* L1946: */
#line 1059 "sgejsv.f"
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

#line 1073 "sgejsv.f"
	if (! almort) {

#line 1075 "sgejsv.f"
	    if (l2pert) {
/*              XSC = SQRT(SMALL) */
#line 1077 "sgejsv.f"
		xsc = epsln / (doublereal) (*n);
#line 1078 "sgejsv.f"
		i__1 = nr;
#line 1078 "sgejsv.f"
		for (q = 1; q <= i__1; ++q) {
#line 1079 "sgejsv.f"
		    temp1 = xsc * (d__1 = a[q + q * a_dim1], abs(d__1));
#line 1080 "sgejsv.f"
		    i__2 = *n;
#line 1080 "sgejsv.f"
		    for (p = 1; p <= i__2; ++p) {
#line 1081 "sgejsv.f"
			if (p > q && (d__1 = a[p + q * a_dim1], abs(d__1)) <= 
				temp1 || p < q) {
#line 1081 "sgejsv.f"
			    a[p + q * a_dim1] = d_sign(&temp1, &a[p + q * 
				    a_dim1]);
#line 1081 "sgejsv.f"
			}
#line 1084 "sgejsv.f"
/* L4949: */
#line 1084 "sgejsv.f"
		    }
#line 1085 "sgejsv.f"
/* L4947: */
#line 1085 "sgejsv.f"
		}
#line 1086 "sgejsv.f"
	    } else {
#line 1087 "sgejsv.f"
		i__1 = nr - 1;
#line 1087 "sgejsv.f"
		i__2 = nr - 1;
#line 1087 "sgejsv.f"
		slaset_("U", &i__1, &i__2, &c_b34, &c_b34, &a[(a_dim1 << 1) + 
			1], lda, (ftnlen)1);
#line 1088 "sgejsv.f"
	    }

/*            .. second preconditioning using the QR factorization */

#line 1092 "sgejsv.f"
	    i__1 = *lwork - *n;
#line 1092 "sgejsv.f"
	    sgeqrf_(n, &nr, &a[a_offset], lda, &work[1], &work[*n + 1], &i__1,
		     &ierr);

/*           .. and transpose upper to lower triangular */
#line 1095 "sgejsv.f"
	    i__1 = nr - 1;
#line 1095 "sgejsv.f"
	    for (p = 1; p <= i__1; ++p) {
#line 1096 "sgejsv.f"
		i__2 = nr - p;
#line 1096 "sgejsv.f"
		scopy_(&i__2, &a[p + (p + 1) * a_dim1], lda, &a[p + 1 + p * 
			a_dim1], &c__1);
#line 1097 "sgejsv.f"
/* L1948: */
#line 1097 "sgejsv.f"
	    }

#line 1099 "sgejsv.f"
	}

/*           Row-cyclic Jacobi SVD algorithm with column pivoting */

/*           .. again some perturbation (a "background noise") is added */
/*           to drown denormals */
#line 1105 "sgejsv.f"
	if (l2pert) {
/*              XSC = SQRT(SMALL) */
#line 1107 "sgejsv.f"
	    xsc = epsln / (doublereal) (*n);
#line 1108 "sgejsv.f"
	    i__1 = nr;
#line 1108 "sgejsv.f"
	    for (q = 1; q <= i__1; ++q) {
#line 1109 "sgejsv.f"
		temp1 = xsc * (d__1 = a[q + q * a_dim1], abs(d__1));
#line 1110 "sgejsv.f"
		i__2 = nr;
#line 1110 "sgejsv.f"
		for (p = 1; p <= i__2; ++p) {
#line 1111 "sgejsv.f"
		    if (p > q && (d__1 = a[p + q * a_dim1], abs(d__1)) <= 
			    temp1 || p < q) {
#line 1111 "sgejsv.f"
			a[p + q * a_dim1] = d_sign(&temp1, &a[p + q * a_dim1])
				;
#line 1111 "sgejsv.f"
		    }
#line 1114 "sgejsv.f"
/* L1949: */
#line 1114 "sgejsv.f"
		}
#line 1115 "sgejsv.f"
/* L1947: */
#line 1115 "sgejsv.f"
	    }
#line 1116 "sgejsv.f"
	} else {
#line 1117 "sgejsv.f"
	    i__1 = nr - 1;
#line 1117 "sgejsv.f"
	    i__2 = nr - 1;
#line 1117 "sgejsv.f"
	    slaset_("U", &i__1, &i__2, &c_b34, &c_b34, &a[(a_dim1 << 1) + 1], 
		    lda, (ftnlen)1);
#line 1118 "sgejsv.f"
	}

/*           .. and one-sided Jacobi rotations are started on a lower */
/*           triangular matrix (plus perturbation which is ignored in */
/*           the part which destroys triangular form (confusing?!)) */

#line 1124 "sgejsv.f"
	sgesvj_("L", "NoU", "NoV", &nr, &nr, &a[a_offset], lda, &sva[1], n, &
		v[v_offset], ldv, &work[1], lwork, info, (ftnlen)1, (ftnlen)3,
		 (ftnlen)3);

#line 1127 "sgejsv.f"
	scalem = work[1];
#line 1128 "sgejsv.f"
	numrank = i_dnnt(&work[2]);


#line 1131 "sgejsv.f"
    } else if (rsvec && ! lsvec) {

/*        -> Singular Values and Right Singular Vectors <- */

#line 1135 "sgejsv.f"
	if (almort) {

/*           .. in this case NR equals N */
#line 1138 "sgejsv.f"
	    i__1 = nr;
#line 1138 "sgejsv.f"
	    for (p = 1; p <= i__1; ++p) {
#line 1139 "sgejsv.f"
		i__2 = *n - p + 1;
#line 1139 "sgejsv.f"
		scopy_(&i__2, &a[p + p * a_dim1], lda, &v[p + p * v_dim1], &
			c__1);
#line 1140 "sgejsv.f"
/* L1998: */
#line 1140 "sgejsv.f"
	    }
#line 1141 "sgejsv.f"
	    i__1 = nr - 1;
#line 1141 "sgejsv.f"
	    i__2 = nr - 1;
#line 1141 "sgejsv.f"
	    slaset_("Upper", &i__1, &i__2, &c_b34, &c_b34, &v[(v_dim1 << 1) + 
		    1], ldv, (ftnlen)5);

#line 1143 "sgejsv.f"
	    sgesvj_("L", "U", "N", n, &nr, &v[v_offset], ldv, &sva[1], &nr, &
		    a[a_offset], lda, &work[1], lwork, info, (ftnlen)1, (
		    ftnlen)1, (ftnlen)1);
#line 1145 "sgejsv.f"
	    scalem = work[1];
#line 1146 "sgejsv.f"
	    numrank = i_dnnt(&work[2]);
#line 1148 "sgejsv.f"
	} else {

/*        .. two more QR factorizations ( one QRF is not enough, two require */
/*        accumulated product of Jacobi rotations, three are perfect ) */

#line 1153 "sgejsv.f"
	    i__1 = nr - 1;
#line 1153 "sgejsv.f"
	    i__2 = nr - 1;
#line 1153 "sgejsv.f"
	    slaset_("Lower", &i__1, &i__2, &c_b34, &c_b34, &a[a_dim1 + 2], 
		    lda, (ftnlen)5);
#line 1154 "sgejsv.f"
	    i__1 = *lwork - *n;
#line 1154 "sgejsv.f"
	    sgelqf_(&nr, n, &a[a_offset], lda, &work[1], &work[*n + 1], &i__1,
		     &ierr);
#line 1155 "sgejsv.f"
	    slacpy_("Lower", &nr, &nr, &a[a_offset], lda, &v[v_offset], ldv, (
		    ftnlen)5);
#line 1156 "sgejsv.f"
	    i__1 = nr - 1;
#line 1156 "sgejsv.f"
	    i__2 = nr - 1;
#line 1156 "sgejsv.f"
	    slaset_("Upper", &i__1, &i__2, &c_b34, &c_b34, &v[(v_dim1 << 1) + 
		    1], ldv, (ftnlen)5);
#line 1157 "sgejsv.f"
	    i__1 = *lwork - (*n << 1);
#line 1157 "sgejsv.f"
	    sgeqrf_(&nr, &nr, &v[v_offset], ldv, &work[*n + 1], &work[(*n << 
		    1) + 1], &i__1, &ierr);
#line 1159 "sgejsv.f"
	    i__1 = nr;
#line 1159 "sgejsv.f"
	    for (p = 1; p <= i__1; ++p) {
#line 1160 "sgejsv.f"
		i__2 = nr - p + 1;
#line 1160 "sgejsv.f"
		scopy_(&i__2, &v[p + p * v_dim1], ldv, &v[p + p * v_dim1], &
			c__1);
#line 1161 "sgejsv.f"
/* L8998: */
#line 1161 "sgejsv.f"
	    }
#line 1162 "sgejsv.f"
	    i__1 = nr - 1;
#line 1162 "sgejsv.f"
	    i__2 = nr - 1;
#line 1162 "sgejsv.f"
	    slaset_("Upper", &i__1, &i__2, &c_b34, &c_b34, &v[(v_dim1 << 1) + 
		    1], ldv, (ftnlen)5);

#line 1164 "sgejsv.f"
	    i__1 = *lwork - *n;
#line 1164 "sgejsv.f"
	    sgesvj_("Lower", "U", "N", &nr, &nr, &v[v_offset], ldv, &sva[1], &
		    nr, &u[u_offset], ldu, &work[*n + 1], &i__1, info, (
		    ftnlen)5, (ftnlen)1, (ftnlen)1);
#line 1166 "sgejsv.f"
	    scalem = work[*n + 1];
#line 1167 "sgejsv.f"
	    numrank = i_dnnt(&work[*n + 2]);
#line 1168 "sgejsv.f"
	    if (nr < *n) {
#line 1169 "sgejsv.f"
		i__1 = *n - nr;
#line 1169 "sgejsv.f"
		slaset_("A", &i__1, &nr, &c_b34, &c_b34, &v[nr + 1 + v_dim1], 
			ldv, (ftnlen)1);
#line 1170 "sgejsv.f"
		i__1 = *n - nr;
#line 1170 "sgejsv.f"
		slaset_("A", &nr, &i__1, &c_b34, &c_b34, &v[(nr + 1) * v_dim1 
			+ 1], ldv, (ftnlen)1);
#line 1171 "sgejsv.f"
		i__1 = *n - nr;
#line 1171 "sgejsv.f"
		i__2 = *n - nr;
#line 1171 "sgejsv.f"
		slaset_("A", &i__1, &i__2, &c_b34, &c_b35, &v[nr + 1 + (nr + 
			1) * v_dim1], ldv, (ftnlen)1);
#line 1172 "sgejsv.f"
	    }

#line 1174 "sgejsv.f"
	    i__1 = *lwork - *n;
#line 1174 "sgejsv.f"
	    sormlq_("Left", "Transpose", n, n, &nr, &a[a_offset], lda, &work[
		    1], &v[v_offset], ldv, &work[*n + 1], &i__1, &ierr, (
		    ftnlen)4, (ftnlen)9);

#line 1177 "sgejsv.f"
	}

#line 1179 "sgejsv.f"
	i__1 = *n;
#line 1179 "sgejsv.f"
	for (p = 1; p <= i__1; ++p) {
#line 1180 "sgejsv.f"
	    scopy_(n, &v[p + v_dim1], ldv, &a[iwork[p] + a_dim1], lda);
#line 1181 "sgejsv.f"
/* L8991: */
#line 1181 "sgejsv.f"
	}
#line 1182 "sgejsv.f"
	slacpy_("All", n, n, &a[a_offset], lda, &v[v_offset], ldv, (ftnlen)3);

#line 1184 "sgejsv.f"
	if (transp) {
#line 1185 "sgejsv.f"
	    slacpy_("All", n, n, &v[v_offset], ldv, &u[u_offset], ldu, (
		    ftnlen)3);
#line 1186 "sgejsv.f"
	}

#line 1188 "sgejsv.f"
    } else if (lsvec && ! rsvec) {

/*        .. Singular Values and Left Singular Vectors                 .. */

/*        .. second preconditioning step to avoid need to accumulate */
/*        Jacobi rotations in the Jacobi iterations. */
#line 1194 "sgejsv.f"
	i__1 = nr;
#line 1194 "sgejsv.f"
	for (p = 1; p <= i__1; ++p) {
#line 1195 "sgejsv.f"
	    i__2 = *n - p + 1;
#line 1195 "sgejsv.f"
	    scopy_(&i__2, &a[p + p * a_dim1], lda, &u[p + p * u_dim1], &c__1);
#line 1196 "sgejsv.f"
/* L1965: */
#line 1196 "sgejsv.f"
	}
#line 1197 "sgejsv.f"
	i__1 = nr - 1;
#line 1197 "sgejsv.f"
	i__2 = nr - 1;
#line 1197 "sgejsv.f"
	slaset_("Upper", &i__1, &i__2, &c_b34, &c_b34, &u[(u_dim1 << 1) + 1], 
		ldu, (ftnlen)5);

#line 1199 "sgejsv.f"
	i__1 = *lwork - (*n << 1);
#line 1199 "sgejsv.f"
	sgeqrf_(n, &nr, &u[u_offset], ldu, &work[*n + 1], &work[(*n << 1) + 1]
		, &i__1, &ierr);

#line 1202 "sgejsv.f"
	i__1 = nr - 1;
#line 1202 "sgejsv.f"
	for (p = 1; p <= i__1; ++p) {
#line 1203 "sgejsv.f"
	    i__2 = nr - p;
#line 1203 "sgejsv.f"
	    scopy_(&i__2, &u[p + (p + 1) * u_dim1], ldu, &u[p + 1 + p * 
		    u_dim1], &c__1);
#line 1204 "sgejsv.f"
/* L1967: */
#line 1204 "sgejsv.f"
	}
#line 1205 "sgejsv.f"
	i__1 = nr - 1;
#line 1205 "sgejsv.f"
	i__2 = nr - 1;
#line 1205 "sgejsv.f"
	slaset_("Upper", &i__1, &i__2, &c_b34, &c_b34, &u[(u_dim1 << 1) + 1], 
		ldu, (ftnlen)5);

#line 1207 "sgejsv.f"
	i__1 = *lwork - *n;
#line 1207 "sgejsv.f"
	sgesvj_("Lower", "U", "N", &nr, &nr, &u[u_offset], ldu, &sva[1], &nr, 
		&a[a_offset], lda, &work[*n + 1], &i__1, info, (ftnlen)5, (
		ftnlen)1, (ftnlen)1);
#line 1209 "sgejsv.f"
	scalem = work[*n + 1];
#line 1210 "sgejsv.f"
	numrank = i_dnnt(&work[*n + 2]);

#line 1212 "sgejsv.f"
	if (nr < *m) {
#line 1213 "sgejsv.f"
	    i__1 = *m - nr;
#line 1213 "sgejsv.f"
	    slaset_("A", &i__1, &nr, &c_b34, &c_b34, &u[nr + 1 + u_dim1], ldu,
		     (ftnlen)1);
#line 1214 "sgejsv.f"
	    if (nr < n1) {
#line 1215 "sgejsv.f"
		i__1 = n1 - nr;
#line 1215 "sgejsv.f"
		slaset_("A", &nr, &i__1, &c_b34, &c_b34, &u[(nr + 1) * u_dim1 
			+ 1], ldu, (ftnlen)1);
#line 1216 "sgejsv.f"
		i__1 = *m - nr;
#line 1216 "sgejsv.f"
		i__2 = n1 - nr;
#line 1216 "sgejsv.f"
		slaset_("A", &i__1, &i__2, &c_b34, &c_b35, &u[nr + 1 + (nr + 
			1) * u_dim1], ldu, (ftnlen)1);
#line 1217 "sgejsv.f"
	    }
#line 1218 "sgejsv.f"
	}

#line 1220 "sgejsv.f"
	i__1 = *lwork - *n;
#line 1220 "sgejsv.f"
	sormqr_("Left", "No Tr", m, &n1, n, &a[a_offset], lda, &work[1], &u[
		u_offset], ldu, &work[*n + 1], &i__1, &ierr, (ftnlen)4, (
		ftnlen)5);

#line 1223 "sgejsv.f"
	if (rowpiv) {
#line 1223 "sgejsv.f"
	    i__1 = *m - 1;
#line 1223 "sgejsv.f"
	    slaswp_(&n1, &u[u_offset], ldu, &c__1, &i__1, &iwork[(*n << 1) + 
		    1], &c_n1);
#line 1223 "sgejsv.f"
	}

#line 1226 "sgejsv.f"
	i__1 = n1;
#line 1226 "sgejsv.f"
	for (p = 1; p <= i__1; ++p) {
#line 1227 "sgejsv.f"
	    xsc = 1. / snrm2_(m, &u[p * u_dim1 + 1], &c__1);
#line 1228 "sgejsv.f"
	    sscal_(m, &xsc, &u[p * u_dim1 + 1], &c__1);
#line 1229 "sgejsv.f"
/* L1974: */
#line 1229 "sgejsv.f"
	}

#line 1231 "sgejsv.f"
	if (transp) {
#line 1232 "sgejsv.f"
	    slacpy_("All", n, n, &u[u_offset], ldu, &v[v_offset], ldv, (
		    ftnlen)3);
#line 1233 "sgejsv.f"
	}

#line 1235 "sgejsv.f"
    } else {

/*        .. Full SVD .. */

#line 1239 "sgejsv.f"
	if (! jracc) {

#line 1241 "sgejsv.f"
	    if (! almort) {

/*           Second Preconditioning Step (QRF [with pivoting]) */
/*           Note that the composition of TRANSPOSE, QRF and TRANSPOSE is */
/*           equivalent to an LQF CALL. Since in many libraries the QRF */
/*           seems to be better optimized than the LQF, we do explicit */
/*           transpose and use the QRF. This is subject to changes in an */
/*           optimized implementation of SGEJSV. */

#line 1250 "sgejsv.f"
		i__1 = nr;
#line 1250 "sgejsv.f"
		for (p = 1; p <= i__1; ++p) {
#line 1251 "sgejsv.f"
		    i__2 = *n - p + 1;
#line 1251 "sgejsv.f"
		    scopy_(&i__2, &a[p + p * a_dim1], lda, &v[p + p * v_dim1],
			     &c__1);
#line 1252 "sgejsv.f"
/* L1968: */
#line 1252 "sgejsv.f"
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

#line 1266 "sgejsv.f"
		if (l2pert) {
#line 1267 "sgejsv.f"
		    xsc = sqrt(small);
#line 1268 "sgejsv.f"
		    i__1 = nr;
#line 1268 "sgejsv.f"
		    for (q = 1; q <= i__1; ++q) {
#line 1269 "sgejsv.f"
			temp1 = xsc * (d__1 = v[q + q * v_dim1], abs(d__1));
#line 1270 "sgejsv.f"
			i__2 = *n;
#line 1270 "sgejsv.f"
			for (p = 1; p <= i__2; ++p) {
#line 1271 "sgejsv.f"
			    if (p > q && (d__1 = v[p + q * v_dim1], abs(d__1))
				     <= temp1 || p < q) {
#line 1271 "sgejsv.f"
				v[p + q * v_dim1] = d_sign(&temp1, &v[p + q * 
					v_dim1]);
#line 1271 "sgejsv.f"
			    }
#line 1274 "sgejsv.f"
			    if (p < q) {
#line 1274 "sgejsv.f"
				v[p + q * v_dim1] = -v[p + q * v_dim1];
#line 1274 "sgejsv.f"
			    }
#line 1275 "sgejsv.f"
/* L2968: */
#line 1275 "sgejsv.f"
			}
#line 1276 "sgejsv.f"
/* L2969: */
#line 1276 "sgejsv.f"
		    }
#line 1277 "sgejsv.f"
		} else {
#line 1278 "sgejsv.f"
		    i__1 = nr - 1;
#line 1278 "sgejsv.f"
		    i__2 = nr - 1;
#line 1278 "sgejsv.f"
		    slaset_("U", &i__1, &i__2, &c_b34, &c_b34, &v[(v_dim1 << 
			    1) + 1], ldv, (ftnlen)1);
#line 1279 "sgejsv.f"
		}

/*           Estimate the row scaled condition number of R1 */
/*           (If R1 is rectangular, N > NR, then the condition number */
/*           of the leading NR x NR submatrix is estimated.) */

#line 1285 "sgejsv.f"
		slacpy_("L", &nr, &nr, &v[v_offset], ldv, &work[(*n << 1) + 1]
			, &nr, (ftnlen)1);
#line 1286 "sgejsv.f"
		i__1 = nr;
#line 1286 "sgejsv.f"
		for (p = 1; p <= i__1; ++p) {
#line 1287 "sgejsv.f"
		    i__2 = nr - p + 1;
#line 1287 "sgejsv.f"
		    temp1 = snrm2_(&i__2, &work[(*n << 1) + (p - 1) * nr + p],
			     &c__1);
#line 1288 "sgejsv.f"
		    i__2 = nr - p + 1;
#line 1288 "sgejsv.f"
		    d__1 = 1. / temp1;
#line 1288 "sgejsv.f"
		    sscal_(&i__2, &d__1, &work[(*n << 1) + (p - 1) * nr + p], 
			    &c__1);
#line 1289 "sgejsv.f"
/* L3950: */
#line 1289 "sgejsv.f"
		}
#line 1290 "sgejsv.f"
		spocon_("Lower", &nr, &work[(*n << 1) + 1], &nr, &c_b35, &
			temp1, &work[(*n << 1) + nr * nr + 1], &iwork[*m + (*
			n << 1) + 1], &ierr, (ftnlen)5);
#line 1292 "sgejsv.f"
		condr1 = 1. / sqrt(temp1);
/*           .. here need a second oppinion on the condition number */
/*           .. then assume worst case scenario */
/*           R1 is OK for inverse <=> CONDR1 .LT. FLOAT(N) */
/*           more conservative    <=> CONDR1 .LT. SQRT(FLOAT(N)) */

#line 1298 "sgejsv.f"
		cond_ok__ = sqrt((doublereal) nr);
/* [TP]       COND_OK is a tuning parameter. */
#line 1301 "sgejsv.f"
		if (condr1 < cond_ok__) {
/*              .. the second QRF without pivoting. Note: in an optimized */
/*              implementation, this QRF should be implemented as the QRF */
/*              of a lower triangular matrix. */
/*              R1^t = Q2 * R2 */
#line 1306 "sgejsv.f"
		    i__1 = *lwork - (*n << 1);
#line 1306 "sgejsv.f"
		    sgeqrf_(n, &nr, &v[v_offset], ldv, &work[*n + 1], &work[(*
			    n << 1) + 1], &i__1, &ierr);

#line 1309 "sgejsv.f"
		    if (l2pert) {
#line 1310 "sgejsv.f"
			xsc = sqrt(small) / epsln;
#line 1311 "sgejsv.f"
			i__1 = nr;
#line 1311 "sgejsv.f"
			for (p = 2; p <= i__1; ++p) {
#line 1312 "sgejsv.f"
			    i__2 = p - 1;
#line 1312 "sgejsv.f"
			    for (q = 1; q <= i__2; ++q) {
/* Computing MIN */
#line 1313 "sgejsv.f"
				d__3 = (d__1 = v[p + p * v_dim1], abs(d__1)), 
					d__4 = (d__2 = v[q + q * v_dim1], abs(
					d__2));
#line 1313 "sgejsv.f"
				temp1 = xsc * min(d__3,d__4);
#line 1314 "sgejsv.f"
				if ((d__1 = v[q + p * v_dim1], abs(d__1)) <= 
					temp1) {
#line 1314 "sgejsv.f"
				    v[q + p * v_dim1] = d_sign(&temp1, &v[q + 
					    p * v_dim1]);
#line 1314 "sgejsv.f"
				}
#line 1316 "sgejsv.f"
/* L3958: */
#line 1316 "sgejsv.f"
			    }
#line 1317 "sgejsv.f"
/* L3959: */
#line 1317 "sgejsv.f"
			}
#line 1318 "sgejsv.f"
		    }

#line 1320 "sgejsv.f"
		    if (nr != *n) {
#line 1320 "sgejsv.f"
			slacpy_("A", n, &nr, &v[v_offset], ldv, &work[(*n << 
				1) + 1], n, (ftnlen)1);
#line 1320 "sgejsv.f"
		    }
/*              .. save ... */

/*           .. this transposed copy should be better than naive */
#line 1325 "sgejsv.f"
		    i__1 = nr - 1;
#line 1325 "sgejsv.f"
		    for (p = 1; p <= i__1; ++p) {
#line 1326 "sgejsv.f"
			i__2 = nr - p;
#line 1326 "sgejsv.f"
			scopy_(&i__2, &v[p + (p + 1) * v_dim1], ldv, &v[p + 1 
				+ p * v_dim1], &c__1);
#line 1327 "sgejsv.f"
/* L1969: */
#line 1327 "sgejsv.f"
		    }

#line 1329 "sgejsv.f"
		    condr2 = condr1;

#line 1331 "sgejsv.f"
		} else {

/*              .. ill-conditioned case: second QRF with pivoting */
/*              Note that windowed pivoting would be equaly good */
/*              numerically, and more run-time efficient. So, in */
/*              an optimal implementation, the next call to SGEQP3 */
/*              should be replaced with eg. CALL SGEQPX (ACM TOMS #782) */
/*              with properly (carefully) chosen parameters. */

/*              R1^t * P2 = Q2 * R2 */
#line 1341 "sgejsv.f"
		    i__1 = nr;
#line 1341 "sgejsv.f"
		    for (p = 1; p <= i__1; ++p) {
#line 1342 "sgejsv.f"
			iwork[*n + p] = 0;
#line 1343 "sgejsv.f"
/* L3003: */
#line 1343 "sgejsv.f"
		    }
#line 1344 "sgejsv.f"
		    i__1 = *lwork - (*n << 1);
#line 1344 "sgejsv.f"
		    sgeqp3_(n, &nr, &v[v_offset], ldv, &iwork[*n + 1], &work[*
			    n + 1], &work[(*n << 1) + 1], &i__1, &ierr);
/* *               CALL SGEQRF( N, NR, V, LDV, WORK(N+1), WORK(2*N+1), */
/* *     $              LWORK-2*N, IERR ) */
#line 1348 "sgejsv.f"
		    if (l2pert) {
#line 1349 "sgejsv.f"
			xsc = sqrt(small);
#line 1350 "sgejsv.f"
			i__1 = nr;
#line 1350 "sgejsv.f"
			for (p = 2; p <= i__1; ++p) {
#line 1351 "sgejsv.f"
			    i__2 = p - 1;
#line 1351 "sgejsv.f"
			    for (q = 1; q <= i__2; ++q) {
/* Computing MIN */
#line 1352 "sgejsv.f"
				d__3 = (d__1 = v[p + p * v_dim1], abs(d__1)), 
					d__4 = (d__2 = v[q + q * v_dim1], abs(
					d__2));
#line 1352 "sgejsv.f"
				temp1 = xsc * min(d__3,d__4);
#line 1353 "sgejsv.f"
				if ((d__1 = v[q + p * v_dim1], abs(d__1)) <= 
					temp1) {
#line 1353 "sgejsv.f"
				    v[q + p * v_dim1] = d_sign(&temp1, &v[q + 
					    p * v_dim1]);
#line 1353 "sgejsv.f"
				}
#line 1355 "sgejsv.f"
/* L3968: */
#line 1355 "sgejsv.f"
			    }
#line 1356 "sgejsv.f"
/* L3969: */
#line 1356 "sgejsv.f"
			}
#line 1357 "sgejsv.f"
		    }

#line 1359 "sgejsv.f"
		    slacpy_("A", n, &nr, &v[v_offset], ldv, &work[(*n << 1) + 
			    1], n, (ftnlen)1);

#line 1361 "sgejsv.f"
		    if (l2pert) {
#line 1362 "sgejsv.f"
			xsc = sqrt(small);
#line 1363 "sgejsv.f"
			i__1 = nr;
#line 1363 "sgejsv.f"
			for (p = 2; p <= i__1; ++p) {
#line 1364 "sgejsv.f"
			    i__2 = p - 1;
#line 1364 "sgejsv.f"
			    for (q = 1; q <= i__2; ++q) {
/* Computing MIN */
#line 1365 "sgejsv.f"
				d__3 = (d__1 = v[p + p * v_dim1], abs(d__1)), 
					d__4 = (d__2 = v[q + q * v_dim1], abs(
					d__2));
#line 1365 "sgejsv.f"
				temp1 = xsc * min(d__3,d__4);
#line 1366 "sgejsv.f"
				v[p + q * v_dim1] = -d_sign(&temp1, &v[q + p *
					 v_dim1]);
#line 1367 "sgejsv.f"
/* L8971: */
#line 1367 "sgejsv.f"
			    }
#line 1368 "sgejsv.f"
/* L8970: */
#line 1368 "sgejsv.f"
			}
#line 1369 "sgejsv.f"
		    } else {
#line 1370 "sgejsv.f"
			i__1 = nr - 1;
#line 1370 "sgejsv.f"
			i__2 = nr - 1;
#line 1370 "sgejsv.f"
			slaset_("L", &i__1, &i__2, &c_b34, &c_b34, &v[v_dim1 
				+ 2], ldv, (ftnlen)1);
#line 1371 "sgejsv.f"
		    }
/*              Now, compute R2 = L3 * Q3, the LQ factorization. */
#line 1373 "sgejsv.f"
		    i__1 = *lwork - (*n << 1) - *n * nr - nr;
#line 1373 "sgejsv.f"
		    sgelqf_(&nr, &nr, &v[v_offset], ldv, &work[(*n << 1) + *n 
			    * nr + 1], &work[(*n << 1) + *n * nr + nr + 1], &
			    i__1, &ierr);
/*              .. and estimate the condition number */
#line 1376 "sgejsv.f"
		    slacpy_("L", &nr, &nr, &v[v_offset], ldv, &work[(*n << 1) 
			    + *n * nr + nr + 1], &nr, (ftnlen)1);
#line 1377 "sgejsv.f"
		    i__1 = nr;
#line 1377 "sgejsv.f"
		    for (p = 1; p <= i__1; ++p) {
#line 1378 "sgejsv.f"
			temp1 = snrm2_(&p, &work[(*n << 1) + *n * nr + nr + p]
				, &nr);
#line 1379 "sgejsv.f"
			d__1 = 1. / temp1;
#line 1379 "sgejsv.f"
			sscal_(&p, &d__1, &work[(*n << 1) + *n * nr + nr + p],
				 &nr);
#line 1380 "sgejsv.f"
/* L4950: */
#line 1380 "sgejsv.f"
		    }
#line 1381 "sgejsv.f"
		    spocon_("L", &nr, &work[(*n << 1) + *n * nr + nr + 1], &
			    nr, &c_b35, &temp1, &work[(*n << 1) + *n * nr + 
			    nr + nr * nr + 1], &iwork[*m + (*n << 1) + 1], &
			    ierr, (ftnlen)1);
#line 1383 "sgejsv.f"
		    condr2 = 1. / sqrt(temp1);

#line 1385 "sgejsv.f"
		    if (condr2 >= cond_ok__) {
/*                 .. save the Householder vectors used for Q3 */
/*                 (this overwrittes the copy of R2, as it will not be */
/*                 needed in this branch, but it does not overwritte the */
/*                 Huseholder vectors of Q2.). */
#line 1390 "sgejsv.f"
			slacpy_("U", &nr, &nr, &v[v_offset], ldv, &work[(*n <<
				 1) + 1], n, (ftnlen)1);
/*                 .. and the rest of the information on Q3 is in */
/*                 WORK(2*N+N*NR+1:2*N+N*NR+N) */
#line 1393 "sgejsv.f"
		    }

#line 1395 "sgejsv.f"
		}

#line 1397 "sgejsv.f"
		if (l2pert) {
#line 1398 "sgejsv.f"
		    xsc = sqrt(small);
#line 1399 "sgejsv.f"
		    i__1 = nr;
#line 1399 "sgejsv.f"
		    for (q = 2; q <= i__1; ++q) {
#line 1400 "sgejsv.f"
			temp1 = xsc * v[q + q * v_dim1];
#line 1401 "sgejsv.f"
			i__2 = q - 1;
#line 1401 "sgejsv.f"
			for (p = 1; p <= i__2; ++p) {
/*                    V(p,q) = - SIGN( TEMP1, V(q,p) ) */
#line 1403 "sgejsv.f"
			    v[p + q * v_dim1] = -d_sign(&temp1, &v[p + q * 
				    v_dim1]);
#line 1404 "sgejsv.f"
/* L4969: */
#line 1404 "sgejsv.f"
			}
#line 1405 "sgejsv.f"
/* L4968: */
#line 1405 "sgejsv.f"
		    }
#line 1406 "sgejsv.f"
		} else {
#line 1407 "sgejsv.f"
		    i__1 = nr - 1;
#line 1407 "sgejsv.f"
		    i__2 = nr - 1;
#line 1407 "sgejsv.f"
		    slaset_("U", &i__1, &i__2, &c_b34, &c_b34, &v[(v_dim1 << 
			    1) + 1], ldv, (ftnlen)1);
#line 1408 "sgejsv.f"
		}

/*        Second preconditioning finished; continue with Jacobi SVD */
/*        The input matrix is lower trinagular. */

/*        Recover the right singular vectors as solution of a well */
/*        conditioned triangular matrix equation. */

#line 1416 "sgejsv.f"
		if (condr1 < cond_ok__) {

#line 1418 "sgejsv.f"
		    i__1 = *lwork - (*n << 1) - *n * nr - nr;
#line 1418 "sgejsv.f"
		    sgesvj_("L", "U", "N", &nr, &nr, &v[v_offset], ldv, &sva[
			    1], &nr, &u[u_offset], ldu, &work[(*n << 1) + *n *
			     nr + nr + 1], &i__1, info, (ftnlen)1, (ftnlen)1, 
			    (ftnlen)1);
#line 1420 "sgejsv.f"
		    scalem = work[(*n << 1) + *n * nr + nr + 1];
#line 1421 "sgejsv.f"
		    numrank = i_dnnt(&work[(*n << 1) + *n * nr + nr + 2]);
#line 1422 "sgejsv.f"
		    i__1 = nr;
#line 1422 "sgejsv.f"
		    for (p = 1; p <= i__1; ++p) {
#line 1423 "sgejsv.f"
			scopy_(&nr, &v[p * v_dim1 + 1], &c__1, &u[p * u_dim1 
				+ 1], &c__1);
#line 1424 "sgejsv.f"
			sscal_(&nr, &sva[p], &v[p * v_dim1 + 1], &c__1);
#line 1425 "sgejsv.f"
/* L3970: */
#line 1425 "sgejsv.f"
		    }
/*        .. pick the right matrix equation and solve it */

#line 1429 "sgejsv.f"
		    if (nr == *n) {
/* :))             .. best case, R1 is inverted. The solution of this matrix */
/*                 equation is Q2*V2 = the product of the Jacobi rotations */
/*                 used in SGESVJ, premultiplied with the orthogonal matrix */
/*                 from the second QR factorization. */
#line 1434 "sgejsv.f"
			strsm_("L", "U", "N", "N", &nr, &nr, &c_b35, &a[
				a_offset], lda, &v[v_offset], ldv, (ftnlen)1, 
				(ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1435 "sgejsv.f"
		    } else {
/*                 .. R1 is well conditioned, but non-square. Transpose(R2) */
/*                 is inverted to get the product of the Jacobi rotations */
/*                 used in SGESVJ. The Q-factor from the second QR */
/*                 factorization is then built in explicitly. */
#line 1440 "sgejsv.f"
			strsm_("L", "U", "T", "N", &nr, &nr, &c_b35, &work[(*
				n << 1) + 1], n, &v[v_offset], ldv, (ftnlen)1,
				 (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1442 "sgejsv.f"
			if (nr < *n) {
#line 1443 "sgejsv.f"
			    i__1 = *n - nr;
#line 1443 "sgejsv.f"
			    slaset_("A", &i__1, &nr, &c_b34, &c_b34, &v[nr + 
				    1 + v_dim1], ldv, (ftnlen)1);
#line 1444 "sgejsv.f"
			    i__1 = *n - nr;
#line 1444 "sgejsv.f"
			    slaset_("A", &nr, &i__1, &c_b34, &c_b34, &v[(nr + 
				    1) * v_dim1 + 1], ldv, (ftnlen)1);
#line 1445 "sgejsv.f"
			    i__1 = *n - nr;
#line 1445 "sgejsv.f"
			    i__2 = *n - nr;
#line 1445 "sgejsv.f"
			    slaset_("A", &i__1, &i__2, &c_b34, &c_b35, &v[nr 
				    + 1 + (nr + 1) * v_dim1], ldv, (ftnlen)1);
#line 1446 "sgejsv.f"
			}
#line 1447 "sgejsv.f"
			i__1 = *lwork - (*n << 1) - *n * nr - nr;
#line 1447 "sgejsv.f"
			sormqr_("L", "N", n, n, &nr, &work[(*n << 1) + 1], n, 
				&work[*n + 1], &v[v_offset], ldv, &work[(*n <<
				 1) + *n * nr + nr + 1], &i__1, &ierr, (
				ftnlen)1, (ftnlen)1);
#line 1449 "sgejsv.f"
		    }

#line 1451 "sgejsv.f"
		} else if (condr2 < cond_ok__) {

/* :)           .. the input matrix A is very likely a relative of */
/*              the Kahan matrix :) */
/*              The matrix R2 is inverted. The solution of the matrix equation */
/*              is Q3^T*V3 = the product of the Jacobi rotations (appplied to */
/*              the lower triangular L3 from the LQ factorization of */
/*              R2=L3*Q3), pre-multiplied with the transposed Q3. */
#line 1459 "sgejsv.f"
		    i__1 = *lwork - (*n << 1) - *n * nr - nr;
#line 1459 "sgejsv.f"
		    sgesvj_("L", "U", "N", &nr, &nr, &v[v_offset], ldv, &sva[
			    1], &nr, &u[u_offset], ldu, &work[(*n << 1) + *n *
			     nr + nr + 1], &i__1, info, (ftnlen)1, (ftnlen)1, 
			    (ftnlen)1);
#line 1461 "sgejsv.f"
		    scalem = work[(*n << 1) + *n * nr + nr + 1];
#line 1462 "sgejsv.f"
		    numrank = i_dnnt(&work[(*n << 1) + *n * nr + nr + 2]);
#line 1463 "sgejsv.f"
		    i__1 = nr;
#line 1463 "sgejsv.f"
		    for (p = 1; p <= i__1; ++p) {
#line 1464 "sgejsv.f"
			scopy_(&nr, &v[p * v_dim1 + 1], &c__1, &u[p * u_dim1 
				+ 1], &c__1);
#line 1465 "sgejsv.f"
			sscal_(&nr, &sva[p], &u[p * u_dim1 + 1], &c__1);
#line 1466 "sgejsv.f"
/* L3870: */
#line 1466 "sgejsv.f"
		    }
#line 1467 "sgejsv.f"
		    strsm_("L", "U", "N", "N", &nr, &nr, &c_b35, &work[(*n << 
			    1) + 1], n, &u[u_offset], ldu, (ftnlen)1, (ftnlen)
			    1, (ftnlen)1, (ftnlen)1);
/*              .. apply the permutation from the second QR factorization */
#line 1469 "sgejsv.f"
		    i__1 = nr;
#line 1469 "sgejsv.f"
		    for (q = 1; q <= i__1; ++q) {
#line 1470 "sgejsv.f"
			i__2 = nr;
#line 1470 "sgejsv.f"
			for (p = 1; p <= i__2; ++p) {
#line 1471 "sgejsv.f"
			    work[(*n << 1) + *n * nr + nr + iwork[*n + p]] = 
				    u[p + q * u_dim1];
#line 1472 "sgejsv.f"
/* L872: */
#line 1472 "sgejsv.f"
			}
#line 1473 "sgejsv.f"
			i__2 = nr;
#line 1473 "sgejsv.f"
			for (p = 1; p <= i__2; ++p) {
#line 1474 "sgejsv.f"
			    u[p + q * u_dim1] = work[(*n << 1) + *n * nr + nr 
				    + p];
#line 1475 "sgejsv.f"
/* L874: */
#line 1475 "sgejsv.f"
			}
#line 1476 "sgejsv.f"
/* L873: */
#line 1476 "sgejsv.f"
		    }
#line 1477 "sgejsv.f"
		    if (nr < *n) {
#line 1478 "sgejsv.f"
			i__1 = *n - nr;
#line 1478 "sgejsv.f"
			slaset_("A", &i__1, &nr, &c_b34, &c_b34, &v[nr + 1 + 
				v_dim1], ldv, (ftnlen)1);
#line 1479 "sgejsv.f"
			i__1 = *n - nr;
#line 1479 "sgejsv.f"
			slaset_("A", &nr, &i__1, &c_b34, &c_b34, &v[(nr + 1) *
				 v_dim1 + 1], ldv, (ftnlen)1);
#line 1480 "sgejsv.f"
			i__1 = *n - nr;
#line 1480 "sgejsv.f"
			i__2 = *n - nr;
#line 1480 "sgejsv.f"
			slaset_("A", &i__1, &i__2, &c_b34, &c_b35, &v[nr + 1 
				+ (nr + 1) * v_dim1], ldv, (ftnlen)1);
#line 1481 "sgejsv.f"
		    }
#line 1482 "sgejsv.f"
		    i__1 = *lwork - (*n << 1) - *n * nr - nr;
#line 1482 "sgejsv.f"
		    sormqr_("L", "N", n, n, &nr, &work[(*n << 1) + 1], n, &
			    work[*n + 1], &v[v_offset], ldv, &work[(*n << 1) 
			    + *n * nr + nr + 1], &i__1, &ierr, (ftnlen)1, (
			    ftnlen)1);
#line 1484 "sgejsv.f"
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
#line 1496 "sgejsv.f"
		    i__1 = *lwork - (*n << 1) - *n * nr - nr;
#line 1496 "sgejsv.f"
		    sgesvj_("L", "U", "V", &nr, &nr, &v[v_offset], ldv, &sva[
			    1], &nr, &u[u_offset], ldu, &work[(*n << 1) + *n *
			     nr + nr + 1], &i__1, info, (ftnlen)1, (ftnlen)1, 
			    (ftnlen)1);
#line 1498 "sgejsv.f"
		    scalem = work[(*n << 1) + *n * nr + nr + 1];
#line 1499 "sgejsv.f"
		    numrank = i_dnnt(&work[(*n << 1) + *n * nr + nr + 2]);
#line 1500 "sgejsv.f"
		    if (nr < *n) {
#line 1501 "sgejsv.f"
			i__1 = *n - nr;
#line 1501 "sgejsv.f"
			slaset_("A", &i__1, &nr, &c_b34, &c_b34, &v[nr + 1 + 
				v_dim1], ldv, (ftnlen)1);
#line 1502 "sgejsv.f"
			i__1 = *n - nr;
#line 1502 "sgejsv.f"
			slaset_("A", &nr, &i__1, &c_b34, &c_b34, &v[(nr + 1) *
				 v_dim1 + 1], ldv, (ftnlen)1);
#line 1503 "sgejsv.f"
			i__1 = *n - nr;
#line 1503 "sgejsv.f"
			i__2 = *n - nr;
#line 1503 "sgejsv.f"
			slaset_("A", &i__1, &i__2, &c_b34, &c_b35, &v[nr + 1 
				+ (nr + 1) * v_dim1], ldv, (ftnlen)1);
#line 1504 "sgejsv.f"
		    }
#line 1505 "sgejsv.f"
		    i__1 = *lwork - (*n << 1) - *n * nr - nr;
#line 1505 "sgejsv.f"
		    sormqr_("L", "N", n, n, &nr, &work[(*n << 1) + 1], n, &
			    work[*n + 1], &v[v_offset], ldv, &work[(*n << 1) 
			    + *n * nr + nr + 1], &i__1, &ierr, (ftnlen)1, (
			    ftnlen)1);

#line 1508 "sgejsv.f"
		    i__1 = *lwork - (*n << 1) - *n * nr - nr;
#line 1508 "sgejsv.f"
		    sormlq_("L", "T", &nr, &nr, &nr, &work[(*n << 1) + 1], n, 
			    &work[(*n << 1) + *n * nr + 1], &u[u_offset], ldu,
			     &work[(*n << 1) + *n * nr + nr + 1], &i__1, &
			    ierr, (ftnlen)1, (ftnlen)1);
#line 1511 "sgejsv.f"
		    i__1 = nr;
#line 1511 "sgejsv.f"
		    for (q = 1; q <= i__1; ++q) {
#line 1512 "sgejsv.f"
			i__2 = nr;
#line 1512 "sgejsv.f"
			for (p = 1; p <= i__2; ++p) {
#line 1513 "sgejsv.f"
			    work[(*n << 1) + *n * nr + nr + iwork[*n + p]] = 
				    u[p + q * u_dim1];
#line 1514 "sgejsv.f"
/* L772: */
#line 1514 "sgejsv.f"
			}
#line 1515 "sgejsv.f"
			i__2 = nr;
#line 1515 "sgejsv.f"
			for (p = 1; p <= i__2; ++p) {
#line 1516 "sgejsv.f"
			    u[p + q * u_dim1] = work[(*n << 1) + *n * nr + nr 
				    + p];
#line 1517 "sgejsv.f"
/* L774: */
#line 1517 "sgejsv.f"
			}
#line 1518 "sgejsv.f"
/* L773: */
#line 1518 "sgejsv.f"
		    }

#line 1520 "sgejsv.f"
		}

/*           Permute the rows of V using the (column) permutation from the */
/*           first QRF. Also, scale the columns to make them unit in */
/*           Euclidean norm. This applies to all cases. */

#line 1526 "sgejsv.f"
		temp1 = sqrt((doublereal) (*n)) * epsln;
#line 1527 "sgejsv.f"
		i__1 = *n;
#line 1527 "sgejsv.f"
		for (q = 1; q <= i__1; ++q) {
#line 1528 "sgejsv.f"
		    i__2 = *n;
#line 1528 "sgejsv.f"
		    for (p = 1; p <= i__2; ++p) {
#line 1529 "sgejsv.f"
			work[(*n << 1) + *n * nr + nr + iwork[p]] = v[p + q * 
				v_dim1];
#line 1530 "sgejsv.f"
/* L972: */
#line 1530 "sgejsv.f"
		    }
#line 1531 "sgejsv.f"
		    i__2 = *n;
#line 1531 "sgejsv.f"
		    for (p = 1; p <= i__2; ++p) {
#line 1532 "sgejsv.f"
			v[p + q * v_dim1] = work[(*n << 1) + *n * nr + nr + p]
				;
#line 1533 "sgejsv.f"
/* L973: */
#line 1533 "sgejsv.f"
		    }
#line 1534 "sgejsv.f"
		    xsc = 1. / snrm2_(n, &v[q * v_dim1 + 1], &c__1);
#line 1535 "sgejsv.f"
		    if (xsc < 1. - temp1 || xsc > temp1 + 1.) {
#line 1535 "sgejsv.f"
			sscal_(n, &xsc, &v[q * v_dim1 + 1], &c__1);
#line 1535 "sgejsv.f"
		    }
#line 1537 "sgejsv.f"
/* L1972: */
#line 1537 "sgejsv.f"
		}
/*           At this moment, V contains the right singular vectors of A. */
/*           Next, assemble the left singular vector matrix U (M x N). */
#line 1540 "sgejsv.f"
		if (nr < *m) {
#line 1541 "sgejsv.f"
		    i__1 = *m - nr;
#line 1541 "sgejsv.f"
		    slaset_("A", &i__1, &nr, &c_b34, &c_b34, &u[nr + 1 + 
			    u_dim1], ldu, (ftnlen)1);
#line 1542 "sgejsv.f"
		    if (nr < n1) {
#line 1543 "sgejsv.f"
			i__1 = n1 - nr;
#line 1543 "sgejsv.f"
			slaset_("A", &nr, &i__1, &c_b34, &c_b34, &u[(nr + 1) *
				 u_dim1 + 1], ldu, (ftnlen)1);
#line 1544 "sgejsv.f"
			i__1 = *m - nr;
#line 1544 "sgejsv.f"
			i__2 = n1 - nr;
#line 1544 "sgejsv.f"
			slaset_("A", &i__1, &i__2, &c_b34, &c_b35, &u[nr + 1 
				+ (nr + 1) * u_dim1], ldu, (ftnlen)1);
#line 1545 "sgejsv.f"
		    }
#line 1546 "sgejsv.f"
		}

/*           The Q matrix from the first QRF is built into the left singular */
/*           matrix U. This applies to all cases. */

#line 1551 "sgejsv.f"
		i__1 = *lwork - *n;
#line 1551 "sgejsv.f"
		sormqr_("Left", "No_Tr", m, &n1, n, &a[a_offset], lda, &work[
			1], &u[u_offset], ldu, &work[*n + 1], &i__1, &ierr, (
			ftnlen)4, (ftnlen)5);
/*           The columns of U are normalized. The cost is O(M*N) flops. */
#line 1555 "sgejsv.f"
		temp1 = sqrt((doublereal) (*m)) * epsln;
#line 1556 "sgejsv.f"
		i__1 = nr;
#line 1556 "sgejsv.f"
		for (p = 1; p <= i__1; ++p) {
#line 1557 "sgejsv.f"
		    xsc = 1. / snrm2_(m, &u[p * u_dim1 + 1], &c__1);
#line 1558 "sgejsv.f"
		    if (xsc < 1. - temp1 || xsc > temp1 + 1.) {
#line 1558 "sgejsv.f"
			sscal_(m, &xsc, &u[p * u_dim1 + 1], &c__1);
#line 1558 "sgejsv.f"
		    }
#line 1560 "sgejsv.f"
/* L1973: */
#line 1560 "sgejsv.f"
		}

/*           If the initial QRF is computed with row pivoting, the left */
/*           singular vectors must be adjusted. */

#line 1565 "sgejsv.f"
		if (rowpiv) {
#line 1565 "sgejsv.f"
		    i__1 = *m - 1;
#line 1565 "sgejsv.f"
		    slaswp_(&n1, &u[u_offset], ldu, &c__1, &i__1, &iwork[(*n 
			    << 1) + 1], &c_n1);
#line 1565 "sgejsv.f"
		}

#line 1568 "sgejsv.f"
	    } else {

/*        .. the initial matrix A has almost orthogonal columns and */
/*        the second QRF is not needed */

#line 1573 "sgejsv.f"
		slacpy_("Upper", n, n, &a[a_offset], lda, &work[*n + 1], n, (
			ftnlen)5);
#line 1574 "sgejsv.f"
		if (l2pert) {
#line 1575 "sgejsv.f"
		    xsc = sqrt(small);
#line 1576 "sgejsv.f"
		    i__1 = *n;
#line 1576 "sgejsv.f"
		    for (p = 2; p <= i__1; ++p) {
#line 1577 "sgejsv.f"
			temp1 = xsc * work[*n + (p - 1) * *n + p];
#line 1578 "sgejsv.f"
			i__2 = p - 1;
#line 1578 "sgejsv.f"
			for (q = 1; q <= i__2; ++q) {
#line 1579 "sgejsv.f"
			    work[*n + (q - 1) * *n + p] = -d_sign(&temp1, &
				    work[*n + (p - 1) * *n + q]);
#line 1580 "sgejsv.f"
/* L5971: */
#line 1580 "sgejsv.f"
			}
#line 1581 "sgejsv.f"
/* L5970: */
#line 1581 "sgejsv.f"
		    }
#line 1582 "sgejsv.f"
		} else {
#line 1583 "sgejsv.f"
		    i__1 = *n - 1;
#line 1583 "sgejsv.f"
		    i__2 = *n - 1;
#line 1583 "sgejsv.f"
		    slaset_("Lower", &i__1, &i__2, &c_b34, &c_b34, &work[*n + 
			    2], n, (ftnlen)5);
#line 1584 "sgejsv.f"
		}

#line 1586 "sgejsv.f"
		i__1 = *lwork - *n - *n * *n;
#line 1586 "sgejsv.f"
		sgesvj_("Upper", "U", "N", n, n, &work[*n + 1], n, &sva[1], n,
			 &u[u_offset], ldu, &work[*n + *n * *n + 1], &i__1, 
			info, (ftnlen)5, (ftnlen)1, (ftnlen)1);

#line 1589 "sgejsv.f"
		scalem = work[*n + *n * *n + 1];
#line 1590 "sgejsv.f"
		numrank = i_dnnt(&work[*n + *n * *n + 2]);
#line 1591 "sgejsv.f"
		i__1 = *n;
#line 1591 "sgejsv.f"
		for (p = 1; p <= i__1; ++p) {
#line 1592 "sgejsv.f"
		    scopy_(n, &work[*n + (p - 1) * *n + 1], &c__1, &u[p * 
			    u_dim1 + 1], &c__1);
#line 1593 "sgejsv.f"
		    sscal_(n, &sva[p], &work[*n + (p - 1) * *n + 1], &c__1);
#line 1594 "sgejsv.f"
/* L6970: */
#line 1594 "sgejsv.f"
		}

#line 1596 "sgejsv.f"
		strsm_("Left", "Upper", "NoTrans", "No UD", n, n, &c_b35, &a[
			a_offset], lda, &work[*n + 1], n, (ftnlen)4, (ftnlen)
			5, (ftnlen)7, (ftnlen)5);
#line 1598 "sgejsv.f"
		i__1 = *n;
#line 1598 "sgejsv.f"
		for (p = 1; p <= i__1; ++p) {
#line 1599 "sgejsv.f"
		    scopy_(n, &work[*n + p], n, &v[iwork[p] + v_dim1], ldv);
#line 1600 "sgejsv.f"
/* L6972: */
#line 1600 "sgejsv.f"
		}
#line 1601 "sgejsv.f"
		temp1 = sqrt((doublereal) (*n)) * epsln;
#line 1602 "sgejsv.f"
		i__1 = *n;
#line 1602 "sgejsv.f"
		for (p = 1; p <= i__1; ++p) {
#line 1603 "sgejsv.f"
		    xsc = 1. / snrm2_(n, &v[p * v_dim1 + 1], &c__1);
#line 1604 "sgejsv.f"
		    if (xsc < 1. - temp1 || xsc > temp1 + 1.) {
#line 1604 "sgejsv.f"
			sscal_(n, &xsc, &v[p * v_dim1 + 1], &c__1);
#line 1604 "sgejsv.f"
		    }
#line 1606 "sgejsv.f"
/* L6971: */
#line 1606 "sgejsv.f"
		}

/*           Assemble the left singular vector matrix U (M x N). */

#line 1610 "sgejsv.f"
		if (*n < *m) {
#line 1611 "sgejsv.f"
		    i__1 = *m - *n;
#line 1611 "sgejsv.f"
		    slaset_("A", &i__1, n, &c_b34, &c_b34, &u[*n + 1 + u_dim1]
			    , ldu, (ftnlen)1);
#line 1612 "sgejsv.f"
		    if (*n < n1) {
#line 1613 "sgejsv.f"
			i__1 = n1 - *n;
#line 1613 "sgejsv.f"
			slaset_("A", n, &i__1, &c_b34, &c_b34, &u[(*n + 1) * 
				u_dim1 + 1], ldu, (ftnlen)1);
#line 1614 "sgejsv.f"
			i__1 = *m - *n;
#line 1614 "sgejsv.f"
			i__2 = n1 - *n;
#line 1614 "sgejsv.f"
			slaset_("A", &i__1, &i__2, &c_b34, &c_b35, &u[*n + 1 
				+ (*n + 1) * u_dim1], ldu, (ftnlen)1);
#line 1615 "sgejsv.f"
		    }
#line 1616 "sgejsv.f"
		}
#line 1617 "sgejsv.f"
		i__1 = *lwork - *n;
#line 1617 "sgejsv.f"
		sormqr_("Left", "No Tr", m, &n1, n, &a[a_offset], lda, &work[
			1], &u[u_offset], ldu, &work[*n + 1], &i__1, &ierr, (
			ftnlen)4, (ftnlen)5);
#line 1619 "sgejsv.f"
		temp1 = sqrt((doublereal) (*m)) * epsln;
#line 1620 "sgejsv.f"
		i__1 = n1;
#line 1620 "sgejsv.f"
		for (p = 1; p <= i__1; ++p) {
#line 1621 "sgejsv.f"
		    xsc = 1. / snrm2_(m, &u[p * u_dim1 + 1], &c__1);
#line 1622 "sgejsv.f"
		    if (xsc < 1. - temp1 || xsc > temp1 + 1.) {
#line 1622 "sgejsv.f"
			sscal_(m, &xsc, &u[p * u_dim1 + 1], &c__1);
#line 1622 "sgejsv.f"
		    }
#line 1624 "sgejsv.f"
/* L6973: */
#line 1624 "sgejsv.f"
		}

#line 1626 "sgejsv.f"
		if (rowpiv) {
#line 1626 "sgejsv.f"
		    i__1 = *m - 1;
#line 1626 "sgejsv.f"
		    slaswp_(&n1, &u[u_offset], ldu, &c__1, &i__1, &iwork[(*n 
			    << 1) + 1], &c_n1);
#line 1626 "sgejsv.f"
		}

#line 1629 "sgejsv.f"
	    }

/*        end of the  >> almost orthogonal case <<  in the full SVD */

#line 1633 "sgejsv.f"
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

#line 1645 "sgejsv.f"
	    i__1 = nr;
#line 1645 "sgejsv.f"
	    for (p = 1; p <= i__1; ++p) {
#line 1646 "sgejsv.f"
		i__2 = *n - p + 1;
#line 1646 "sgejsv.f"
		scopy_(&i__2, &a[p + p * a_dim1], lda, &v[p + p * v_dim1], &
			c__1);
#line 1647 "sgejsv.f"
/* L7968: */
#line 1647 "sgejsv.f"
	    }

#line 1649 "sgejsv.f"
	    if (l2pert) {
#line 1650 "sgejsv.f"
		xsc = sqrt(small / epsln);
#line 1651 "sgejsv.f"
		i__1 = nr;
#line 1651 "sgejsv.f"
		for (q = 1; q <= i__1; ++q) {
#line 1652 "sgejsv.f"
		    temp1 = xsc * (d__1 = v[q + q * v_dim1], abs(d__1));
#line 1653 "sgejsv.f"
		    i__2 = *n;
#line 1653 "sgejsv.f"
		    for (p = 1; p <= i__2; ++p) {
#line 1654 "sgejsv.f"
			if (p > q && (d__1 = v[p + q * v_dim1], abs(d__1)) <= 
				temp1 || p < q) {
#line 1654 "sgejsv.f"
			    v[p + q * v_dim1] = d_sign(&temp1, &v[p + q * 
				    v_dim1]);
#line 1654 "sgejsv.f"
			}
#line 1657 "sgejsv.f"
			if (p < q) {
#line 1657 "sgejsv.f"
			    v[p + q * v_dim1] = -v[p + q * v_dim1];
#line 1657 "sgejsv.f"
			}
#line 1658 "sgejsv.f"
/* L5968: */
#line 1658 "sgejsv.f"
		    }
#line 1659 "sgejsv.f"
/* L5969: */
#line 1659 "sgejsv.f"
		}
#line 1660 "sgejsv.f"
	    } else {
#line 1661 "sgejsv.f"
		i__1 = nr - 1;
#line 1661 "sgejsv.f"
		i__2 = nr - 1;
#line 1661 "sgejsv.f"
		slaset_("U", &i__1, &i__2, &c_b34, &c_b34, &v[(v_dim1 << 1) + 
			1], ldv, (ftnlen)1);
#line 1662 "sgejsv.f"
	    }
#line 1664 "sgejsv.f"
	    i__1 = *lwork - (*n << 1);
#line 1664 "sgejsv.f"
	    sgeqrf_(n, &nr, &v[v_offset], ldv, &work[*n + 1], &work[(*n << 1) 
		    + 1], &i__1, &ierr);
#line 1666 "sgejsv.f"
	    slacpy_("L", n, &nr, &v[v_offset], ldv, &work[(*n << 1) + 1], n, (
		    ftnlen)1);

#line 1668 "sgejsv.f"
	    i__1 = nr;
#line 1668 "sgejsv.f"
	    for (p = 1; p <= i__1; ++p) {
#line 1669 "sgejsv.f"
		i__2 = nr - p + 1;
#line 1669 "sgejsv.f"
		scopy_(&i__2, &v[p + p * v_dim1], ldv, &u[p + p * u_dim1], &
			c__1);
#line 1670 "sgejsv.f"
/* L7969: */
#line 1670 "sgejsv.f"
	    }
#line 1672 "sgejsv.f"
	    if (l2pert) {
#line 1673 "sgejsv.f"
		xsc = sqrt(small / epsln);
#line 1674 "sgejsv.f"
		i__1 = nr;
#line 1674 "sgejsv.f"
		for (q = 2; q <= i__1; ++q) {
#line 1675 "sgejsv.f"
		    i__2 = q - 1;
#line 1675 "sgejsv.f"
		    for (p = 1; p <= i__2; ++p) {
/* Computing MIN */
#line 1676 "sgejsv.f"
			d__3 = (d__1 = u[p + p * u_dim1], abs(d__1)), d__4 = (
				d__2 = u[q + q * u_dim1], abs(d__2));
#line 1676 "sgejsv.f"
			temp1 = xsc * min(d__3,d__4);
#line 1677 "sgejsv.f"
			u[p + q * u_dim1] = -d_sign(&temp1, &u[q + p * u_dim1]
				);
#line 1678 "sgejsv.f"
/* L9971: */
#line 1678 "sgejsv.f"
		    }
#line 1679 "sgejsv.f"
/* L9970: */
#line 1679 "sgejsv.f"
		}
#line 1680 "sgejsv.f"
	    } else {
#line 1681 "sgejsv.f"
		i__1 = nr - 1;
#line 1681 "sgejsv.f"
		i__2 = nr - 1;
#line 1681 "sgejsv.f"
		slaset_("U", &i__1, &i__2, &c_b34, &c_b34, &u[(u_dim1 << 1) + 
			1], ldu, (ftnlen)1);
#line 1682 "sgejsv.f"
	    }
#line 1684 "sgejsv.f"
	    i__1 = *lwork - (*n << 1) - *n * nr;
#line 1684 "sgejsv.f"
	    sgesvj_("L", "U", "V", &nr, &nr, &u[u_offset], ldu, &sva[1], n, &
		    v[v_offset], ldv, &work[(*n << 1) + *n * nr + 1], &i__1, 
		    info, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1686 "sgejsv.f"
	    scalem = work[(*n << 1) + *n * nr + 1];
#line 1687 "sgejsv.f"
	    numrank = i_dnnt(&work[(*n << 1) + *n * nr + 2]);
#line 1689 "sgejsv.f"
	    if (nr < *n) {
#line 1690 "sgejsv.f"
		i__1 = *n - nr;
#line 1690 "sgejsv.f"
		slaset_("A", &i__1, &nr, &c_b34, &c_b34, &v[nr + 1 + v_dim1], 
			ldv, (ftnlen)1);
#line 1691 "sgejsv.f"
		i__1 = *n - nr;
#line 1691 "sgejsv.f"
		slaset_("A", &nr, &i__1, &c_b34, &c_b34, &v[(nr + 1) * v_dim1 
			+ 1], ldv, (ftnlen)1);
#line 1692 "sgejsv.f"
		i__1 = *n - nr;
#line 1692 "sgejsv.f"
		i__2 = *n - nr;
#line 1692 "sgejsv.f"
		slaset_("A", &i__1, &i__2, &c_b34, &c_b35, &v[nr + 1 + (nr + 
			1) * v_dim1], ldv, (ftnlen)1);
#line 1693 "sgejsv.f"
	    }
#line 1695 "sgejsv.f"
	    i__1 = *lwork - (*n << 1) - *n * nr - nr;
#line 1695 "sgejsv.f"
	    sormqr_("L", "N", n, n, &nr, &work[(*n << 1) + 1], n, &work[*n + 
		    1], &v[v_offset], ldv, &work[(*n << 1) + *n * nr + nr + 1]
		    , &i__1, &ierr, (ftnlen)1, (ftnlen)1);

/*           Permute the rows of V using the (column) permutation from the */
/*           first QRF. Also, scale the columns to make them unit in */
/*           Euclidean norm. This applies to all cases. */

#line 1702 "sgejsv.f"
	    temp1 = sqrt((doublereal) (*n)) * epsln;
#line 1703 "sgejsv.f"
	    i__1 = *n;
#line 1703 "sgejsv.f"
	    for (q = 1; q <= i__1; ++q) {
#line 1704 "sgejsv.f"
		i__2 = *n;
#line 1704 "sgejsv.f"
		for (p = 1; p <= i__2; ++p) {
#line 1705 "sgejsv.f"
		    work[(*n << 1) + *n * nr + nr + iwork[p]] = v[p + q * 
			    v_dim1];
#line 1706 "sgejsv.f"
/* L8972: */
#line 1706 "sgejsv.f"
		}
#line 1707 "sgejsv.f"
		i__2 = *n;
#line 1707 "sgejsv.f"
		for (p = 1; p <= i__2; ++p) {
#line 1708 "sgejsv.f"
		    v[p + q * v_dim1] = work[(*n << 1) + *n * nr + nr + p];
#line 1709 "sgejsv.f"
/* L8973: */
#line 1709 "sgejsv.f"
		}
#line 1710 "sgejsv.f"
		xsc = 1. / snrm2_(n, &v[q * v_dim1 + 1], &c__1);
#line 1711 "sgejsv.f"
		if (xsc < 1. - temp1 || xsc > temp1 + 1.) {
#line 1711 "sgejsv.f"
		    sscal_(n, &xsc, &v[q * v_dim1 + 1], &c__1);
#line 1711 "sgejsv.f"
		}
#line 1713 "sgejsv.f"
/* L7972: */
#line 1713 "sgejsv.f"
	    }

/*           At this moment, V contains the right singular vectors of A. */
/*           Next, assemble the left singular vector matrix U (M x N). */

#line 1718 "sgejsv.f"
	    if (nr < *m) {
#line 1719 "sgejsv.f"
		i__1 = *m - nr;
#line 1719 "sgejsv.f"
		slaset_("A", &i__1, &nr, &c_b34, &c_b34, &u[nr + 1 + u_dim1], 
			ldu, (ftnlen)1);
#line 1720 "sgejsv.f"
		if (nr < n1) {
#line 1721 "sgejsv.f"
		    i__1 = n1 - nr;
#line 1721 "sgejsv.f"
		    slaset_("A", &nr, &i__1, &c_b34, &c_b34, &u[(nr + 1) * 
			    u_dim1 + 1], ldu, (ftnlen)1);
#line 1722 "sgejsv.f"
		    i__1 = *m - nr;
#line 1722 "sgejsv.f"
		    i__2 = n1 - nr;
#line 1722 "sgejsv.f"
		    slaset_("A", &i__1, &i__2, &c_b34, &c_b35, &u[nr + 1 + (
			    nr + 1) * u_dim1], ldu, (ftnlen)1);
#line 1723 "sgejsv.f"
		}
#line 1724 "sgejsv.f"
	    }

#line 1726 "sgejsv.f"
	    i__1 = *lwork - *n;
#line 1726 "sgejsv.f"
	    sormqr_("Left", "No Tr", m, &n1, n, &a[a_offset], lda, &work[1], &
		    u[u_offset], ldu, &work[*n + 1], &i__1, &ierr, (ftnlen)4, 
		    (ftnlen)5);

#line 1729 "sgejsv.f"
	    if (rowpiv) {
#line 1729 "sgejsv.f"
		i__1 = *m - 1;
#line 1729 "sgejsv.f"
		slaswp_(&n1, &u[u_offset], ldu, &c__1, &i__1, &iwork[(*n << 1)
			 + 1], &c_n1);
#line 1729 "sgejsv.f"
	    }


#line 1733 "sgejsv.f"
	}
#line 1734 "sgejsv.f"
	if (transp) {
/*           .. swap U and V because the procedure worked on A^t */
#line 1736 "sgejsv.f"
	    i__1 = *n;
#line 1736 "sgejsv.f"
	    for (p = 1; p <= i__1; ++p) {
#line 1737 "sgejsv.f"
		sswap_(n, &u[p * u_dim1 + 1], &c__1, &v[p * v_dim1 + 1], &
			c__1);
#line 1738 "sgejsv.f"
/* L6974: */
#line 1738 "sgejsv.f"
	    }
#line 1739 "sgejsv.f"
	}

#line 1741 "sgejsv.f"
    }
/*     end of the full SVD */

/*     Undo scaling, if necessary (and possible) */

#line 1746 "sgejsv.f"
    if (uscal2 <= big / sva[1] * uscal1) {
#line 1747 "sgejsv.f"
	slascl_("G", &c__0, &c__0, &uscal1, &uscal2, &nr, &c__1, &sva[1], n, &
		ierr, (ftnlen)1);
#line 1748 "sgejsv.f"
	uscal1 = 1.;
#line 1749 "sgejsv.f"
	uscal2 = 1.;
#line 1750 "sgejsv.f"
    }

#line 1752 "sgejsv.f"
    if (nr < *n) {
#line 1753 "sgejsv.f"
	i__1 = *n;
#line 1753 "sgejsv.f"
	for (p = nr + 1; p <= i__1; ++p) {
#line 1754 "sgejsv.f"
	    sva[p] = 0.;
#line 1755 "sgejsv.f"
/* L3004: */
#line 1755 "sgejsv.f"
	}
#line 1756 "sgejsv.f"
    }

#line 1758 "sgejsv.f"
    work[1] = uscal2 * scalem;
#line 1759 "sgejsv.f"
    work[2] = uscal1;
#line 1760 "sgejsv.f"
    if (errest) {
#line 1760 "sgejsv.f"
	work[3] = sconda;
#line 1760 "sgejsv.f"
    }
#line 1761 "sgejsv.f"
    if (lsvec && rsvec) {
#line 1762 "sgejsv.f"
	work[4] = condr1;
#line 1763 "sgejsv.f"
	work[5] = condr2;
#line 1764 "sgejsv.f"
    }
#line 1765 "sgejsv.f"
    if (l2tran) {
#line 1766 "sgejsv.f"
	work[6] = entra;
#line 1767 "sgejsv.f"
	work[7] = entrat;
#line 1768 "sgejsv.f"
    }

#line 1770 "sgejsv.f"
    iwork[1] = nr;
#line 1771 "sgejsv.f"
    iwork[2] = numrank;
#line 1772 "sgejsv.f"
    iwork[3] = warning;

#line 1774 "sgejsv.f"
    return 0;
/*     .. */
/*     .. END OF SGEJSV */
/*     .. */
} /* sgejsv_ */

