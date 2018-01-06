#line 1 "dgejsv.f"
/* dgejsv.f -- translated by f2c (version 20100827).
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

#line 1 "dgejsv.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b34 = 0.;
static doublereal c_b35 = 1.;
static integer c__0 = 0;
static integer c_n1 = -1;

/* > \brief \b DGEJSV */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DGEJSV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgejsv.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgejsv.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgejsv.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DGEJSV( JOBA, JOBU, JOBV, JOBR, JOBT, JOBP, */
/*                          M, N, A, LDA, SVA, U, LDU, V, LDV, */
/*                          WORK, LWORK, IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       IMPLICIT    NONE */
/*       INTEGER     INFO, LDA, LDU, LDV, LWORK, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION A( LDA, * ), SVA( N ), U( LDU, * ), V( LDV, * ), */
/*      $            WORK( LWORK ) */
/*       INTEGER     IWORK( * ) */
/*       CHARACTER*1 JOBA, JOBP, JOBR, JOBT, JOBU, JOBV */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DGEJSV computes the singular value decomposition (SVD) of a real M-by-N */
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
/* > DGEJSV can sometimes compute tiny singular values and their singular vectors much */
/* > more accurately than other SVD routines, see below under Further Details.*> \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] JOBA */
/* > \verbatim */
/* >          JOBA is CHARACTER*1 */
/* >        Specifies the level of accuracy: */
/* >       = 'C': This option works well (high relative accuracy) if A = B * D, */
/* >             with well-conditioned B and arbitrary diagonal matrix D. */
/* >             The accuracy cannot be spoiled by COLUMN scaling. The */
/* >             accuracy of the computed output depends on the condition of */
/* >             B, and the procedure aims at the best theoretical accuracy. */
/* >             The relative error max_{i=1:N}|d sigma_i| / sigma_i is */
/* >             bounded by f(M,N)*epsilon* cond(B), independent of D. */
/* >             The input matrix is preprocessed with the QRF with column */
/* >             pivoting. This initial preprocessing and preconditioning by */
/* >             a rank revealing QR factorization is common for all values of */
/* >             JOBA. Additional actions are specified as follows: */
/* >       = 'E': Computation as with 'C' with an additional estimate of the */
/* >             condition number of B. It provides a realistic error bound. */
/* >       = 'F': If A = D1 * C * D2 with ill-conditioned diagonal scalings */
/* >             D1, D2, and well-conditioned matrix C, this option gives */
/* >             higher accuracy than the 'C' option. If the structure of the */
/* >             input matrix is not known, and relative accuracy is */
/* >             desirable, then this option is advisable. The input matrix A */
/* >             is preprocessed with QR factorization with FULL (row and */
/* >             column) pivoting. */
/* >       = 'G'  Computation as with 'F' with an additional estimate of the */
/* >             condition number of B, where A=D*B. If A has heavily weighted */
/* >             rows, then using this condition number gives too pessimistic */
/* >             error bound. */
/* >       = 'A': Small singular values are the noise and the matrix is treated */
/* >             as numerically rank defficient. The error in the computed */
/* >             singular values is bounded by f(m,n)*epsilon*||A||. */
/* >             The computed SVD A = U * S * V^t restores A up to */
/* >             f(m,n)*epsilon*||A||. */
/* >             This gives the procedure the licence to discard (set to zero) */
/* >             all singular values below N*epsilon*||A||. */
/* >       = 'R': Similar as in 'A'. Rank revealing property of the initial */
/* >             QR factorization is used do reveal (using triangular factor) */
/* >             a gap sigma_{r+1} < epsilon * sigma_r in which case the */
/* >             numerical RANK is declared to be r. The SVD is computed with */
/* >             absolute error bounds, but more accurately than with 'A'. */
/* > \endverbatim */
/* > */
/* > \param[in] JOBU */
/* > \verbatim */
/* >          JOBU is CHARACTER*1 */
/* >        Specifies whether to compute the columns of U: */
/* >       = 'U': N columns of U are returned in the array U. */
/* >       = 'F': full set of M left sing. vectors is returned in the array U. */
/* >       = 'W': U may be used as workspace of length M*N. See the description */
/* >             of U. */
/* >       = 'N': U is not computed. */
/* > \endverbatim */
/* > */
/* > \param[in] JOBV */
/* > \verbatim */
/* >          JOBV is CHARACTER*1 */
/* >        Specifies whether to compute the matrix V: */
/* >       = 'V': N columns of V are returned in the array V; Jacobi rotations */
/* >             are not explicitly accumulated. */
/* >       = 'J': N columns of V are returned in the array V, but they are */
/* >             computed as the product of Jacobi rotations. This option is */
/* >             allowed only if JOBU .NE. 'N', i.e. in computing the full SVD. */
/* >       = 'W': V may be used as workspace of length N*N. See the description */
/* >             of V. */
/* >       = 'N': V is not computed. */
/* > \endverbatim */
/* > */
/* > \param[in] JOBR */
/* > \verbatim */
/* >          JOBR is CHARACTER*1 */
/* >        Specifies the RANGE for the singular values. Issues the licence to */
/* >        set to zero small positive singular values if they are outside */
/* >        specified range. If A .NE. 0 is scaled so that the largest singular */
/* >        value of c*A is around DSQRT(BIG), BIG=SLAMCH('O'), then JOBR issues */
/* >        the licence to kill columns of A whose norm in c*A is less than */
/* >        DSQRT(SFMIN) (for JOBR.EQ.'R'), or less than SMALL=SFMIN/EPSLN, */
/* >        where SFMIN=SLAMCH('S'), EPSLN=SLAMCH('E'). */
/* >       = 'N': Do not kill small columns of c*A. This option assumes that */
/* >             BLAS and QR factorizations and triangular solvers are */
/* >             implemented to work in that range. If the condition of A */
/* >             is greater than BIG, use DGESVJ. */
/* >       = 'R': RESTRICTED range for sigma(c*A) is [DSQRT(SFMIN), DSQRT(BIG)] */
/* >             (roughly, as described above). This option is recommended. */
/* >                                            ~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/* >        For computing the singular values in the FULL range [SFMIN,BIG] */
/* >        use DGESVJ. */
/* > \endverbatim */
/* > */
/* > \param[in] JOBT */
/* > \verbatim */
/* >          JOBT is CHARACTER*1 */
/* >        If the matrix is square then the procedure may determine to use */
/* >        transposed A if A^t seems to be better with respect to convergence. */
/* >        If the matrix is not square, JOBT is ignored. This is subject to */
/* >        changes in the future. */
/* >        The decision is based on two values of entropy over the adjoint */
/* >        orbit of A^t * A. See the descriptions of WORK(6) and WORK(7). */
/* >       = 'T': transpose if entropy test indicates possibly faster */
/* >        convergence of Jacobi process if A^t is taken as input. If A is */
/* >        replaced with A^t, then the row pivoting is included automatically. */
/* >       = 'N': do not speculate. */
/* >        This option can be used to compute only the singular values, or the */
/* >        full SVD (U, SIGMA and V). For only one set of singular vectors */
/* >        (U or V), the caller should provide both U and V, as one of the */
/* >        matrices is used as workspace if the matrix A is transposed. */
/* >        The implementer can easily remove this constraint and make the */
/* >        code more complicated. See the descriptions of U and V. */
/* > \endverbatim */
/* > */
/* > \param[in] JOBP */
/* > \verbatim */
/* >          JOBP is CHARACTER*1 */
/* >        Issues the licence to introduce structured perturbations to drown */
/* >        denormalized numbers. This licence should be active if the */
/* >        denormals are poorly implemented, causing slow computation, */
/* >        especially in cases of fast convergence (!). For details see [1,2]. */
/* >        For the sake of simplicity, this perturbations are included only */
/* >        when the full SVD or only the singular values are requested. The */
/* >        implementer/user can easily add the perturbation for the cases of */
/* >        computing one set of singular vectors. */
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
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
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
/* >          U is DOUBLE PRECISION array, dimension ( LDU, N ) */
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
/* >          V is DOUBLE PRECISION array, dimension ( LDV, N ) */
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
/* >          WORK is DOUBLE PRECISION array, dimension at least LWORK. */
/* >          On exit, if N.GT.0 .AND. M.GT.0 (else not referenced), */
/* >          WORK(1) = SCALE = WORK(2) / WORK(1) is the scaling factor such */
/* >                    that SCALE*SVA(1:N) are the computed singular values */
/* >                    of A. (See the description of SVA().) */
/* >          WORK(2) = See the description of WORK(1). */
/* >          WORK(3) = SCONDA is an estimate for the condition number of */
/* >                    column equilibrated A. (If JOBA .EQ. 'E' or 'G') */
/* >                    SCONDA is an estimate of DSQRT(||(R^t * R)^(-1)||_1). */
/* >                    It is computed using DPOCON. It holds */
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
/* >           > 0 :  DGEJSV  did not converge in the maximal allowed number */
/* >                  of sweeps. The computed values may be inaccurate. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2015 */

/* > \ingroup doubleGEsing */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  DGEJSV implements a preconditioned Jacobi SVD algorithm. It uses DGEQP3, */
/* >  DGEQRF, and DGELQF as preprocessors and preconditioners. Optionally, an */
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
/* >  & LAPACK routines called by DGEJSV are implemented to work in that range. */
/* >  If that is not the case, then the restriction for safe computation with */
/* >  the singular values in the range of normalized IEEE numbers is that the */
/* >  spectral condition number kappa(A)=sigma_max(A)/sigma_min(A) does not */
/* >  overflow. This code (DGEJSV) is best used in this restricted range, */
/* >  meaning that singular values of magnitude below ||A||_2 / DLAMCH('O') are */
/* >  returned as zeros. See JOBR for details on this. */
/* >     Further, this implementation is somewhat slower than the one described */
/* >  in [1,2] due to replacement of some non-LAPACK components, and because */
/* >  the choice of some tuning parameters in the iterative part (DGESVJ) is */
/* >  left to the implementer on a particular machine. */
/* >     The rank revealing QR factorization (in this code: DGEQP3) should be */
/* >  implemented as in [3]. We have a new version of DGEQP3 under development */
/* >  that is more robust than the current one in LAPACK, with a cleaner cut in */
/* >  rank defficient cases. It will be available in the SIGMA library [4]. */
/* >  If M is much larger than N, it is obvious that the inital QRF with */
/* >  column pivoting can be preprocessed by the QRF without pivoting. That */
/* >  well known trick is not used in DGEJSV because in some cases heavy row */
/* >  weighting can be treated with complete pivoting. The overhead in cases */
/* >  M much larger than N is then only due to pivoting, but the benefits in */
/* >  terms of accuracy have prevailed. The implementer/user can incorporate */
/* >  this extra QRF step easily. The implementer can also improve data movement */
/* >  (matrix transpose, matrix copy, matrix transposed copy) - this */
/* >  implementation of DGEJSV uses only the simplest, naive data movement. */
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
/* >     ACM Trans. Math. Softw. Vol. 35, No 2 (2008), pp. 1-28. */
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
/* Subroutine */ int dgejsv_(char *joba, char *jobu, char *jobv, char *jobr, 
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
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    static doublereal temp1;
    static logical jracc;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal small, entra, sfmin;
    static logical lsvec;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dswap_(integer *, doublereal *, integer 
	    *, doublereal *, integer *);
    static doublereal epsln;
    static logical rsvec;
    extern /* Subroutine */ int dtrsm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static logical l2aber;
    extern /* Subroutine */ int dgeqp3_(integer *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    integer *);
    static doublereal condr1, condr2, uscal1, uscal2;
    static logical l2kill, l2rank, l2tran, l2pert;
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int dgelqf_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *);
    extern integer idamax_(integer *, doublereal *, integer *);
    static doublereal scalem;
    extern /* Subroutine */ int dlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen);
    static doublereal sconda;
    static logical goscal;
    static doublereal aatmin;
    extern /* Subroutine */ int dgeqrf_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *);
    static doublereal aatmax;
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    dlaset_(char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen), xerbla_(char *, integer *, 
	    ftnlen);
    static logical noscal;
    extern /* Subroutine */ int dpocon_(char *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    integer *, ftnlen), dgesvj_(char *, char *, char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen, ftnlen), dlassq_(integer *, doublereal *, integer 
	    *, doublereal *, doublereal *), dlaswp_(integer *, doublereal *, 
	    integer *, integer *, integer *, integer *, integer *);
    static doublereal entrat;
    static logical almort;
    extern /* Subroutine */ int dorgqr_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *), dormlq_(char *, char *, integer *, integer *, integer 
	    *, doublereal *, integer *, doublereal *, doublereal *, integer *,
	     doublereal *, integer *, integer *, ftnlen, ftnlen);
    static doublereal maxprj;
    static logical errest;
    extern /* Subroutine */ int dormqr_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);
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

#line 528 "dgejsv.f"
    /* Parameter adjustments */
#line 528 "dgejsv.f"
    --sva;
#line 528 "dgejsv.f"
    a_dim1 = *lda;
#line 528 "dgejsv.f"
    a_offset = 1 + a_dim1;
#line 528 "dgejsv.f"
    a -= a_offset;
#line 528 "dgejsv.f"
    u_dim1 = *ldu;
#line 528 "dgejsv.f"
    u_offset = 1 + u_dim1;
#line 528 "dgejsv.f"
    u -= u_offset;
#line 528 "dgejsv.f"
    v_dim1 = *ldv;
#line 528 "dgejsv.f"
    v_offset = 1 + v_dim1;
#line 528 "dgejsv.f"
    v -= v_offset;
#line 528 "dgejsv.f"
    --work;
#line 528 "dgejsv.f"
    --iwork;
#line 528 "dgejsv.f"

#line 528 "dgejsv.f"
    /* Function Body */
#line 528 "dgejsv.f"
    lsvec = lsame_(jobu, "U", (ftnlen)1, (ftnlen)1) || lsame_(jobu, "F", (
	    ftnlen)1, (ftnlen)1);
#line 529 "dgejsv.f"
    jracc = lsame_(jobv, "J", (ftnlen)1, (ftnlen)1);
#line 530 "dgejsv.f"
    rsvec = lsame_(jobv, "V", (ftnlen)1, (ftnlen)1) || jracc;
#line 531 "dgejsv.f"
    rowpiv = lsame_(joba, "F", (ftnlen)1, (ftnlen)1) || lsame_(joba, "G", (
	    ftnlen)1, (ftnlen)1);
#line 532 "dgejsv.f"
    l2rank = lsame_(joba, "R", (ftnlen)1, (ftnlen)1);
#line 533 "dgejsv.f"
    l2aber = lsame_(joba, "A", (ftnlen)1, (ftnlen)1);
#line 534 "dgejsv.f"
    errest = lsame_(joba, "E", (ftnlen)1, (ftnlen)1) || lsame_(joba, "G", (
	    ftnlen)1, (ftnlen)1);
#line 535 "dgejsv.f"
    l2tran = lsame_(jobt, "T", (ftnlen)1, (ftnlen)1);
#line 536 "dgejsv.f"
    l2kill = lsame_(jobr, "R", (ftnlen)1, (ftnlen)1);
#line 537 "dgejsv.f"
    defr = lsame_(jobr, "N", (ftnlen)1, (ftnlen)1);
#line 538 "dgejsv.f"
    l2pert = lsame_(jobp, "P", (ftnlen)1, (ftnlen)1);

#line 540 "dgejsv.f"
    if (! (rowpiv || l2rank || l2aber || errest || lsame_(joba, "C", (ftnlen)
	    1, (ftnlen)1))) {
#line 542 "dgejsv.f"
	*info = -1;
#line 543 "dgejsv.f"
    } else if (! (lsvec || lsame_(jobu, "N", (ftnlen)1, (ftnlen)1) || lsame_(
	    jobu, "W", (ftnlen)1, (ftnlen)1))) {
#line 545 "dgejsv.f"
	*info = -2;
#line 546 "dgejsv.f"
    } else if (! (rsvec || lsame_(jobv, "N", (ftnlen)1, (ftnlen)1) || lsame_(
	    jobv, "W", (ftnlen)1, (ftnlen)1)) || jracc && ! lsvec) {
#line 548 "dgejsv.f"
	*info = -3;
#line 549 "dgejsv.f"
    } else if (! (l2kill || defr)) {
#line 550 "dgejsv.f"
	*info = -4;
#line 551 "dgejsv.f"
    } else if (! (l2tran || lsame_(jobt, "N", (ftnlen)1, (ftnlen)1))) {
#line 552 "dgejsv.f"
	*info = -5;
#line 553 "dgejsv.f"
    } else if (! (l2pert || lsame_(jobp, "N", (ftnlen)1, (ftnlen)1))) {
#line 554 "dgejsv.f"
	*info = -6;
#line 555 "dgejsv.f"
    } else if (*m < 0) {
#line 556 "dgejsv.f"
	*info = -7;
#line 557 "dgejsv.f"
    } else if (*n < 0 || *n > *m) {
#line 558 "dgejsv.f"
	*info = -8;
#line 559 "dgejsv.f"
    } else if (*lda < *m) {
#line 560 "dgejsv.f"
	*info = -10;
#line 561 "dgejsv.f"
    } else if (lsvec && *ldu < *m) {
#line 562 "dgejsv.f"
	*info = -13;
#line 563 "dgejsv.f"
    } else if (rsvec && *ldv < *n) {
#line 564 "dgejsv.f"
	*info = -14;
#line 565 "dgejsv.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 565 "dgejsv.f"
	i__1 = 7, i__2 = (*n << 2) + 1, i__1 = max(i__1,i__2), i__2 = (*m << 
		1) + *n;
/* Computing MAX */
#line 565 "dgejsv.f"
	i__3 = 7, i__4 = (*n << 2) + *n * *n, i__3 = max(i__3,i__4), i__4 = (*
		m << 1) + *n;
/* Computing MAX */
#line 565 "dgejsv.f"
	i__5 = 7, i__6 = (*m << 1) + *n, i__5 = max(i__5,i__6), i__6 = (*n << 
		2) + 1;
/* Computing MAX */
#line 565 "dgejsv.f"
	i__7 = 7, i__8 = (*m << 1) + *n, i__7 = max(i__7,i__8), i__8 = (*n << 
		2) + 1;
/* Computing MAX */
#line 565 "dgejsv.f"
	i__9 = (*m << 1) + *n, i__10 = *n * 6 + (*n << 1) * *n;
/* Computing MAX */
#line 565 "dgejsv.f"
	i__11 = (*m << 1) + *n, i__12 = (*n << 2) + *n * *n, i__11 = max(
		i__11,i__12), i__12 = (*n << 1) + *n * *n + 6;
#line 565 "dgejsv.f"
	if (! (lsvec || rsvec || errest) && *lwork < max(i__1,i__2) || ! (
		lsvec || rsvec) && errest && *lwork < max(i__3,i__4) || lsvec 
		&& ! rsvec && *lwork < max(i__5,i__6) || rsvec && ! lsvec && *
		lwork < max(i__7,i__8) || lsvec && rsvec && ! jracc && *lwork 
		< max(i__9,i__10) || lsvec && rsvec && jracc && *lwork < max(
		i__11,i__12)) {
#line 578 "dgejsv.f"
	    *info = -17;
#line 579 "dgejsv.f"
	} else {
/*        #:) */
#line 581 "dgejsv.f"
	    *info = 0;
#line 582 "dgejsv.f"
	}
#line 582 "dgejsv.f"
    }

#line 584 "dgejsv.f"
    if (*info != 0) {
/*       #:( */
#line 586 "dgejsv.f"
	i__1 = -(*info);
#line 586 "dgejsv.f"
	xerbla_("DGEJSV", &i__1, (ftnlen)6);
#line 587 "dgejsv.f"
	return 0;
#line 588 "dgejsv.f"
    }

/*     Quick return for void matrix (Y3K safe) */
/* #:) */
#line 592 "dgejsv.f"
    if (*m == 0 || *n == 0) {
#line 592 "dgejsv.f"
	return 0;
#line 592 "dgejsv.f"
    }

/*     Determine whether the matrix U should be M x N or M x M */

#line 596 "dgejsv.f"
    if (lsvec) {
#line 597 "dgejsv.f"
	n1 = *n;
#line 598 "dgejsv.f"
	if (lsame_(jobu, "F", (ftnlen)1, (ftnlen)1)) {
#line 598 "dgejsv.f"
	    n1 = *m;
#line 598 "dgejsv.f"
	}
#line 599 "dgejsv.f"
    }

/*     Set numerical parameters */

/* !    NOTE: Make sure DLAMCH() does not fail on the target architecture. */

#line 605 "dgejsv.f"
    epsln = dlamch_("Epsilon", (ftnlen)7);
#line 606 "dgejsv.f"
    sfmin = dlamch_("SafeMinimum", (ftnlen)11);
#line 607 "dgejsv.f"
    small = sfmin / epsln;
#line 608 "dgejsv.f"
    big = dlamch_("O", (ftnlen)1);
/*     BIG   = ONE / SFMIN */

/*     Initialize SVA(1:N) = diag( ||A e_i||_2 )_1^N */

/* (!)  If necessary, scale SVA() to protect the largest norm from */
/*     overflow. It is possible that this scaling pushes the smallest */
/*     column norm left from the underflow threshold (extreme case). */

#line 617 "dgejsv.f"
    scalem = 1. / sqrt((doublereal) (*m) * (doublereal) (*n));
#line 618 "dgejsv.f"
    noscal = TRUE_;
#line 619 "dgejsv.f"
    goscal = TRUE_;
#line 620 "dgejsv.f"
    i__1 = *n;
#line 620 "dgejsv.f"
    for (p = 1; p <= i__1; ++p) {
#line 621 "dgejsv.f"
	aapp = 0.;
#line 622 "dgejsv.f"
	aaqq = 1.;
#line 623 "dgejsv.f"
	dlassq_(m, &a[p * a_dim1 + 1], &c__1, &aapp, &aaqq);
#line 624 "dgejsv.f"
	if (aapp > big) {
#line 625 "dgejsv.f"
	    *info = -9;
#line 626 "dgejsv.f"
	    i__2 = -(*info);
#line 626 "dgejsv.f"
	    xerbla_("DGEJSV", &i__2, (ftnlen)6);
#line 627 "dgejsv.f"
	    return 0;
#line 628 "dgejsv.f"
	}
#line 629 "dgejsv.f"
	aaqq = sqrt(aaqq);
#line 630 "dgejsv.f"
	if (aapp < big / aaqq && noscal) {
#line 631 "dgejsv.f"
	    sva[p] = aapp * aaqq;
#line 632 "dgejsv.f"
	} else {
#line 633 "dgejsv.f"
	    noscal = FALSE_;
#line 634 "dgejsv.f"
	    sva[p] = aapp * (aaqq * scalem);
#line 635 "dgejsv.f"
	    if (goscal) {
#line 636 "dgejsv.f"
		goscal = FALSE_;
#line 637 "dgejsv.f"
		i__2 = p - 1;
#line 637 "dgejsv.f"
		dscal_(&i__2, &scalem, &sva[1], &c__1);
#line 638 "dgejsv.f"
	    }
#line 639 "dgejsv.f"
	}
#line 640 "dgejsv.f"
/* L1874: */
#line 640 "dgejsv.f"
    }

#line 642 "dgejsv.f"
    if (noscal) {
#line 642 "dgejsv.f"
	scalem = 1.;
#line 642 "dgejsv.f"
    }

#line 644 "dgejsv.f"
    aapp = 0.;
#line 645 "dgejsv.f"
    aaqq = big;
#line 646 "dgejsv.f"
    i__1 = *n;
#line 646 "dgejsv.f"
    for (p = 1; p <= i__1; ++p) {
/* Computing MAX */
#line 647 "dgejsv.f"
	d__1 = aapp, d__2 = sva[p];
#line 647 "dgejsv.f"
	aapp = max(d__1,d__2);
#line 648 "dgejsv.f"
	if (sva[p] != 0.) {
/* Computing MIN */
#line 648 "dgejsv.f"
	    d__1 = aaqq, d__2 = sva[p];
#line 648 "dgejsv.f"
	    aaqq = min(d__1,d__2);
#line 648 "dgejsv.f"
	}
#line 649 "dgejsv.f"
/* L4781: */
#line 649 "dgejsv.f"
    }

/*     Quick return for zero M x N matrix */
/* #:) */
#line 653 "dgejsv.f"
    if (aapp == 0.) {
#line 654 "dgejsv.f"
	if (lsvec) {
#line 654 "dgejsv.f"
	    dlaset_("G", m, &n1, &c_b34, &c_b35, &u[u_offset], ldu, (ftnlen)1)
		    ;
#line 654 "dgejsv.f"
	}
#line 655 "dgejsv.f"
	if (rsvec) {
#line 655 "dgejsv.f"
	    dlaset_("G", n, n, &c_b34, &c_b35, &v[v_offset], ldv, (ftnlen)1);
#line 655 "dgejsv.f"
	}
#line 656 "dgejsv.f"
	work[1] = 1.;
#line 657 "dgejsv.f"
	work[2] = 1.;
#line 658 "dgejsv.f"
	if (errest) {
#line 658 "dgejsv.f"
	    work[3] = 1.;
#line 658 "dgejsv.f"
	}
#line 659 "dgejsv.f"
	if (lsvec && rsvec) {
#line 660 "dgejsv.f"
	    work[4] = 1.;
#line 661 "dgejsv.f"
	    work[5] = 1.;
#line 662 "dgejsv.f"
	}
#line 663 "dgejsv.f"
	if (l2tran) {
#line 664 "dgejsv.f"
	    work[6] = 0.;
#line 665 "dgejsv.f"
	    work[7] = 0.;
#line 666 "dgejsv.f"
	}
#line 667 "dgejsv.f"
	iwork[1] = 0;
#line 668 "dgejsv.f"
	iwork[2] = 0;
#line 669 "dgejsv.f"
	iwork[3] = 0;
#line 670 "dgejsv.f"
	return 0;
#line 671 "dgejsv.f"
    }

/*     Issue warning if denormalized column norms detected. Override the */
/*     high relative accuracy request. Issue licence to kill columns */
/*     (set them to zero) whose norm is less than sigma_max / BIG (roughly). */
/* #:( */
#line 677 "dgejsv.f"
    warning = 0;
#line 678 "dgejsv.f"
    if (aaqq <= sfmin) {
#line 679 "dgejsv.f"
	l2rank = TRUE_;
#line 680 "dgejsv.f"
	l2kill = TRUE_;
#line 681 "dgejsv.f"
	warning = 1;
#line 682 "dgejsv.f"
    }

/*     Quick return for one-column matrix */
/* #:) */
#line 686 "dgejsv.f"
    if (*n == 1) {

#line 688 "dgejsv.f"
	if (lsvec) {
#line 689 "dgejsv.f"
	    dlascl_("G", &c__0, &c__0, &sva[1], &scalem, m, &c__1, &a[a_dim1 
		    + 1], lda, &ierr, (ftnlen)1);
#line 690 "dgejsv.f"
	    dlacpy_("A", m, &c__1, &a[a_offset], lda, &u[u_offset], ldu, (
		    ftnlen)1);
/*           computing all M left singular vectors of the M x 1 matrix */
#line 692 "dgejsv.f"
	    if (n1 != *n) {
#line 693 "dgejsv.f"
		i__1 = *lwork - *n;
#line 693 "dgejsv.f"
		dgeqrf_(m, n, &u[u_offset], ldu, &work[1], &work[*n + 1], &
			i__1, &ierr);
#line 694 "dgejsv.f"
		i__1 = *lwork - *n;
#line 694 "dgejsv.f"
		dorgqr_(m, &n1, &c__1, &u[u_offset], ldu, &work[1], &work[*n 
			+ 1], &i__1, &ierr);
#line 695 "dgejsv.f"
		dcopy_(m, &a[a_dim1 + 1], &c__1, &u[u_dim1 + 1], &c__1);
#line 696 "dgejsv.f"
	    }
#line 697 "dgejsv.f"
	}
#line 698 "dgejsv.f"
	if (rsvec) {
#line 699 "dgejsv.f"
	    v[v_dim1 + 1] = 1.;
#line 700 "dgejsv.f"
	}
#line 701 "dgejsv.f"
	if (sva[1] < big * scalem) {
#line 702 "dgejsv.f"
	    sva[1] /= scalem;
#line 703 "dgejsv.f"
	    scalem = 1.;
#line 704 "dgejsv.f"
	}
#line 705 "dgejsv.f"
	work[1] = 1. / scalem;
#line 706 "dgejsv.f"
	work[2] = 1.;
#line 707 "dgejsv.f"
	if (sva[1] != 0.) {
#line 708 "dgejsv.f"
	    iwork[1] = 1;
#line 709 "dgejsv.f"
	    if (sva[1] / scalem >= sfmin) {
#line 710 "dgejsv.f"
		iwork[2] = 1;
#line 711 "dgejsv.f"
	    } else {
#line 712 "dgejsv.f"
		iwork[2] = 0;
#line 713 "dgejsv.f"
	    }
#line 714 "dgejsv.f"
	} else {
#line 715 "dgejsv.f"
	    iwork[1] = 0;
#line 716 "dgejsv.f"
	    iwork[2] = 0;
#line 717 "dgejsv.f"
	}
#line 718 "dgejsv.f"
	if (errest) {
#line 718 "dgejsv.f"
	    work[3] = 1.;
#line 718 "dgejsv.f"
	}
#line 719 "dgejsv.f"
	if (lsvec && rsvec) {
#line 720 "dgejsv.f"
	    work[4] = 1.;
#line 721 "dgejsv.f"
	    work[5] = 1.;
#line 722 "dgejsv.f"
	}
#line 723 "dgejsv.f"
	if (l2tran) {
#line 724 "dgejsv.f"
	    work[6] = 0.;
#line 725 "dgejsv.f"
	    work[7] = 0.;
#line 726 "dgejsv.f"
	}
#line 727 "dgejsv.f"
	return 0;

#line 729 "dgejsv.f"
    }

#line 731 "dgejsv.f"
    transp = FALSE_;
#line 732 "dgejsv.f"
    l2tran = l2tran && *m == *n;

#line 734 "dgejsv.f"
    aatmax = -1.;
#line 735 "dgejsv.f"
    aatmin = big;
#line 736 "dgejsv.f"
    if (rowpiv || l2tran) {

/*     Compute the row norms, needed to determine row pivoting sequence */
/*     (in the case of heavily row weighted A, row pivoting is strongly */
/*     advised) and to collect information needed to compare the */
/*     structures of A * A^t and A^t * A (in the case L2TRAN.EQ..TRUE.). */

#line 743 "dgejsv.f"
	if (l2tran) {
#line 744 "dgejsv.f"
	    i__1 = *m;
#line 744 "dgejsv.f"
	    for (p = 1; p <= i__1; ++p) {
#line 745 "dgejsv.f"
		xsc = 0.;
#line 746 "dgejsv.f"
		temp1 = 1.;
#line 747 "dgejsv.f"
		dlassq_(n, &a[p + a_dim1], lda, &xsc, &temp1);
/*              DLASSQ gets both the ell_2 and the ell_infinity norm */
/*              in one pass through the vector */
#line 750 "dgejsv.f"
		work[*m + *n + p] = xsc * scalem;
#line 751 "dgejsv.f"
		work[*n + p] = xsc * (scalem * sqrt(temp1));
/* Computing MAX */
#line 752 "dgejsv.f"
		d__1 = aatmax, d__2 = work[*n + p];
#line 752 "dgejsv.f"
		aatmax = max(d__1,d__2);
#line 753 "dgejsv.f"
		if (work[*n + p] != 0.) {
/* Computing MIN */
#line 753 "dgejsv.f"
		    d__1 = aatmin, d__2 = work[*n + p];
#line 753 "dgejsv.f"
		    aatmin = min(d__1,d__2);
#line 753 "dgejsv.f"
		}
#line 754 "dgejsv.f"
/* L1950: */
#line 754 "dgejsv.f"
	    }
#line 755 "dgejsv.f"
	} else {
#line 756 "dgejsv.f"
	    i__1 = *m;
#line 756 "dgejsv.f"
	    for (p = 1; p <= i__1; ++p) {
#line 757 "dgejsv.f"
		work[*m + *n + p] = scalem * (d__1 = a[p + idamax_(n, &a[p + 
			a_dim1], lda) * a_dim1], abs(d__1));
/* Computing MAX */
#line 758 "dgejsv.f"
		d__1 = aatmax, d__2 = work[*m + *n + p];
#line 758 "dgejsv.f"
		aatmax = max(d__1,d__2);
/* Computing MIN */
#line 759 "dgejsv.f"
		d__1 = aatmin, d__2 = work[*m + *n + p];
#line 759 "dgejsv.f"
		aatmin = min(d__1,d__2);
#line 760 "dgejsv.f"
/* L1904: */
#line 760 "dgejsv.f"
	    }
#line 761 "dgejsv.f"
	}

#line 763 "dgejsv.f"
    }

/*     For square matrix A try to determine whether A^t  would be  better */
/*     input for the preconditioned Jacobi SVD, with faster convergence. */
/*     The decision is based on an O(N) function of the vector of column */
/*     and row norms of A, based on the Shannon entropy. This should give */
/*     the right choice in most cases when the difference actually matters. */
/*     It may fail and pick the slower converging side. */

#line 772 "dgejsv.f"
    entra = 0.;
#line 773 "dgejsv.f"
    entrat = 0.;
#line 774 "dgejsv.f"
    if (l2tran) {

#line 776 "dgejsv.f"
	xsc = 0.;
#line 777 "dgejsv.f"
	temp1 = 1.;
#line 778 "dgejsv.f"
	dlassq_(n, &sva[1], &c__1, &xsc, &temp1);
#line 779 "dgejsv.f"
	temp1 = 1. / temp1;

#line 781 "dgejsv.f"
	entra = 0.;
#line 782 "dgejsv.f"
	i__1 = *n;
#line 782 "dgejsv.f"
	for (p = 1; p <= i__1; ++p) {
/* Computing 2nd power */
#line 783 "dgejsv.f"
	    d__1 = sva[p] / xsc;
#line 783 "dgejsv.f"
	    big1 = d__1 * d__1 * temp1;
#line 784 "dgejsv.f"
	    if (big1 != 0.) {
#line 784 "dgejsv.f"
		entra += big1 * log(big1);
#line 784 "dgejsv.f"
	    }
#line 785 "dgejsv.f"
/* L1113: */
#line 785 "dgejsv.f"
	}
#line 786 "dgejsv.f"
	entra = -entra / log((doublereal) (*n));

/*        Now, SVA().^2/Trace(A^t * A) is a point in the probability simplex. */
/*        It is derived from the diagonal of  A^t * A.  Do the same with the */
/*        diagonal of A * A^t, compute the entropy of the corresponding */
/*        probability distribution. Note that A * A^t and A^t * A have the */
/*        same trace. */

#line 794 "dgejsv.f"
	entrat = 0.;
#line 795 "dgejsv.f"
	i__1 = *n + *m;
#line 795 "dgejsv.f"
	for (p = *n + 1; p <= i__1; ++p) {
/* Computing 2nd power */
#line 796 "dgejsv.f"
	    d__1 = work[p] / xsc;
#line 796 "dgejsv.f"
	    big1 = d__1 * d__1 * temp1;
#line 797 "dgejsv.f"
	    if (big1 != 0.) {
#line 797 "dgejsv.f"
		entrat += big1 * log(big1);
#line 797 "dgejsv.f"
	    }
#line 798 "dgejsv.f"
/* L1114: */
#line 798 "dgejsv.f"
	}
#line 799 "dgejsv.f"
	entrat = -entrat / log((doublereal) (*m));

/*        Analyze the entropies and decide A or A^t. Smaller entropy */
/*        usually means better input for the algorithm. */

#line 804 "dgejsv.f"
	transp = entrat < entra;

/*        If A^t is better than A, transpose A. */

#line 808 "dgejsv.f"
	if (transp) {
/*           In an optimal implementation, this trivial transpose */
/*           should be replaced with faster transpose. */
#line 811 "dgejsv.f"
	    i__1 = *n - 1;
#line 811 "dgejsv.f"
	    for (p = 1; p <= i__1; ++p) {
#line 812 "dgejsv.f"
		i__2 = *n;
#line 812 "dgejsv.f"
		for (q = p + 1; q <= i__2; ++q) {
#line 813 "dgejsv.f"
		    temp1 = a[q + p * a_dim1];
#line 814 "dgejsv.f"
		    a[q + p * a_dim1] = a[p + q * a_dim1];
#line 815 "dgejsv.f"
		    a[p + q * a_dim1] = temp1;
#line 816 "dgejsv.f"
/* L1116: */
#line 816 "dgejsv.f"
		}
#line 817 "dgejsv.f"
/* L1115: */
#line 817 "dgejsv.f"
	    }
#line 818 "dgejsv.f"
	    i__1 = *n;
#line 818 "dgejsv.f"
	    for (p = 1; p <= i__1; ++p) {
#line 819 "dgejsv.f"
		work[*m + *n + p] = sva[p];
#line 820 "dgejsv.f"
		sva[p] = work[*n + p];
#line 821 "dgejsv.f"
/* L1117: */
#line 821 "dgejsv.f"
	    }
#line 822 "dgejsv.f"
	    temp1 = aapp;
#line 823 "dgejsv.f"
	    aapp = aatmax;
#line 824 "dgejsv.f"
	    aatmax = temp1;
#line 825 "dgejsv.f"
	    temp1 = aaqq;
#line 826 "dgejsv.f"
	    aaqq = aatmin;
#line 827 "dgejsv.f"
	    aatmin = temp1;
#line 828 "dgejsv.f"
	    kill = lsvec;
#line 829 "dgejsv.f"
	    lsvec = rsvec;
#line 830 "dgejsv.f"
	    rsvec = kill;
#line 831 "dgejsv.f"
	    if (lsvec) {
#line 831 "dgejsv.f"
		n1 = *n;
#line 831 "dgejsv.f"
	    }

#line 833 "dgejsv.f"
	    rowpiv = TRUE_;
#line 834 "dgejsv.f"
	}

#line 836 "dgejsv.f"
    }
/*     END IF L2TRAN */

/*     Scale the matrix so that its maximal singular value remains less */
/*     than DSQRT(BIG) -- the matrix is scaled so that its maximal column */
/*     has Euclidean norm equal to DSQRT(BIG/N). The only reason to keep */
/*     DSQRT(BIG) instead of BIG is the fact that DGEJSV uses LAPACK and */
/*     BLAS routines that, in some implementations, are not capable of */
/*     working in the full interval [SFMIN,BIG] and that they may provoke */
/*     overflows in the intermediate results. If the singular values spread */
/*     from SFMIN to BIG, then DGESVJ will compute them. So, in that case, */
/*     one should use DGESVJ instead of DGEJSV. */

#line 849 "dgejsv.f"
    big1 = sqrt(big);
#line 850 "dgejsv.f"
    temp1 = sqrt(big / (doublereal) (*n));

#line 852 "dgejsv.f"
    dlascl_("G", &c__0, &c__0, &aapp, &temp1, n, &c__1, &sva[1], n, &ierr, (
	    ftnlen)1);
#line 853 "dgejsv.f"
    if (aaqq > aapp * sfmin) {
#line 854 "dgejsv.f"
	aaqq = aaqq / aapp * temp1;
#line 855 "dgejsv.f"
    } else {
#line 856 "dgejsv.f"
	aaqq = aaqq * temp1 / aapp;
#line 857 "dgejsv.f"
    }
#line 858 "dgejsv.f"
    temp1 *= scalem;
#line 859 "dgejsv.f"
    dlascl_("G", &c__0, &c__0, &aapp, &temp1, m, n, &a[a_offset], lda, &ierr, 
	    (ftnlen)1);

/*     To undo scaling at the end of this procedure, multiply the */
/*     computed singular values with USCAL2 / USCAL1. */

#line 864 "dgejsv.f"
    uscal1 = temp1;
#line 865 "dgejsv.f"
    uscal2 = aapp;

#line 867 "dgejsv.f"
    if (l2kill) {
/*        L2KILL enforces computation of nonzero singular values in */
/*        the restricted range of condition number of the initial A, */
/*        sigma_max(A) / sigma_min(A) approx. DSQRT(BIG)/DSQRT(SFMIN). */
#line 871 "dgejsv.f"
	xsc = sqrt(sfmin);
#line 872 "dgejsv.f"
    } else {
#line 873 "dgejsv.f"
	xsc = small;

/*        Now, if the condition number of A is too big, */
/*        sigma_max(A) / sigma_min(A) .GT. DSQRT(BIG/N) * EPSLN / SFMIN, */
/*        as a precaution measure, the full SVD is computed using DGESVJ */
/*        with accumulated Jacobi rotations. This provides numerically */
/*        more robust computation, at the cost of slightly increased run */
/*        time. Depending on the concrete implementation of BLAS and LAPACK */
/*        (i.e. how they behave in presence of extreme ill-conditioning) the */
/*        implementor may decide to remove this switch. */
#line 883 "dgejsv.f"
	if (aaqq < sqrt(sfmin) && lsvec && rsvec) {
#line 884 "dgejsv.f"
	    jracc = TRUE_;
#line 885 "dgejsv.f"
	}

#line 887 "dgejsv.f"
    }
#line 888 "dgejsv.f"
    if (aaqq < xsc) {
#line 889 "dgejsv.f"
	i__1 = *n;
#line 889 "dgejsv.f"
	for (p = 1; p <= i__1; ++p) {
#line 890 "dgejsv.f"
	    if (sva[p] < xsc) {
#line 891 "dgejsv.f"
		dlaset_("A", m, &c__1, &c_b34, &c_b34, &a[p * a_dim1 + 1], 
			lda, (ftnlen)1);
#line 892 "dgejsv.f"
		sva[p] = 0.;
#line 893 "dgejsv.f"
	    }
#line 894 "dgejsv.f"
/* L700: */
#line 894 "dgejsv.f"
	}
#line 895 "dgejsv.f"
    }

/*     Preconditioning using QR factorization with pivoting */

#line 899 "dgejsv.f"
    if (rowpiv) {
/*        Optional row permutation (Bjoerck row pivoting): */
/*        A result by Cox and Higham shows that the Bjoerck's */
/*        row pivoting combined with standard column pivoting */
/*        has similar effect as Powell-Reid complete pivoting. */
/*        The ell-infinity norms of A are made nonincreasing. */
#line 905 "dgejsv.f"
	i__1 = *m - 1;
#line 905 "dgejsv.f"
	for (p = 1; p <= i__1; ++p) {
#line 906 "dgejsv.f"
	    i__2 = *m - p + 1;
#line 906 "dgejsv.f"
	    q = idamax_(&i__2, &work[*m + *n + p], &c__1) + p - 1;
#line 907 "dgejsv.f"
	    iwork[(*n << 1) + p] = q;
#line 908 "dgejsv.f"
	    if (p != q) {
#line 909 "dgejsv.f"
		temp1 = work[*m + *n + p];
#line 910 "dgejsv.f"
		work[*m + *n + p] = work[*m + *n + q];
#line 911 "dgejsv.f"
		work[*m + *n + q] = temp1;
#line 912 "dgejsv.f"
	    }
#line 913 "dgejsv.f"
/* L1952: */
#line 913 "dgejsv.f"
	}
#line 914 "dgejsv.f"
	i__1 = *m - 1;
#line 914 "dgejsv.f"
	dlaswp_(n, &a[a_offset], lda, &c__1, &i__1, &iwork[(*n << 1) + 1], &
		c__1);
#line 915 "dgejsv.f"
    }

/*     End of the preparation phase (scaling, optional sorting and */
/*     transposing, optional flushing of small columns). */

/*     Preconditioning */

/*     If the full SVD is needed, the right singular vectors are computed */
/*     from a matrix equation, and for that we need theoretical analysis */
/*     of the Businger-Golub pivoting. So we use DGEQP3 as the first RR QRF. */
/*     In all other cases the first RR QRF can be chosen by other criteria */
/*     (eg speed by replacing global with restricted window pivoting, such */
/*     as in SGEQPX from TOMS # 782). Good results will be obtained using */
/*     SGEQPX with properly (!) chosen numerical parameters. */
/*     Any improvement of DGEQP3 improves overal performance of DGEJSV. */

/*     A * P1 = Q1 * [ R1^t 0]^t: */
#line 932 "dgejsv.f"
    i__1 = *n;
#line 932 "dgejsv.f"
    for (p = 1; p <= i__1; ++p) {
/*        .. all columns are free columns */
#line 934 "dgejsv.f"
	iwork[p] = 0;
#line 935 "dgejsv.f"
/* L1963: */
#line 935 "dgejsv.f"
    }
#line 936 "dgejsv.f"
    i__1 = *lwork - *n;
#line 936 "dgejsv.f"
    dgeqp3_(m, n, &a[a_offset], lda, &iwork[1], &work[1], &work[*n + 1], &
	    i__1, &ierr);

/*     The upper triangular matrix R1 from the first QRF is inspected for */
/*     rank deficiency and possibilities for deflation, or possible */
/*     ill-conditioning. Depending on the user specified flag L2RANK, */
/*     the procedure explores possibilities to reduce the numerical */
/*     rank by inspecting the computed upper triangular factor. If */
/*     L2RANK or L2ABER are up, then DGEJSV will compute the SVD of */
/*     A + dA, where ||dA|| <= f(M,N)*EPSLN. */

#line 946 "dgejsv.f"
    nr = 1;
#line 947 "dgejsv.f"
    if (l2aber) {
/*        Standard absolute error bound suffices. All sigma_i with */
/*        sigma_i < N*EPSLN*||A|| are flushed to zero. This is an */
/*        agressive enforcement of lower numerical rank by introducing a */
/*        backward error of the order of N*EPSLN*||A||. */
#line 952 "dgejsv.f"
	temp1 = sqrt((doublereal) (*n)) * epsln;
#line 953 "dgejsv.f"
	i__1 = *n;
#line 953 "dgejsv.f"
	for (p = 2; p <= i__1; ++p) {
#line 954 "dgejsv.f"
	    if ((d__2 = a[p + p * a_dim1], abs(d__2)) >= temp1 * (d__1 = a[
		    a_dim1 + 1], abs(d__1))) {
#line 955 "dgejsv.f"
		++nr;
#line 956 "dgejsv.f"
	    } else {
#line 957 "dgejsv.f"
		goto L3002;
#line 958 "dgejsv.f"
	    }
#line 959 "dgejsv.f"
/* L3001: */
#line 959 "dgejsv.f"
	}
#line 960 "dgejsv.f"
L3002:
#line 961 "dgejsv.f"
	;
#line 961 "dgejsv.f"
    } else if (l2rank) {
/*        .. similarly as above, only slightly more gentle (less agressive). */
/*        Sudden drop on the diagonal of R1 is used as the criterion for */
/*        close-to-rank-defficient. */
#line 965 "dgejsv.f"
	temp1 = sqrt(sfmin);
#line 966 "dgejsv.f"
	i__1 = *n;
#line 966 "dgejsv.f"
	for (p = 2; p <= i__1; ++p) {
#line 967 "dgejsv.f"
	    if ((d__2 = a[p + p * a_dim1], abs(d__2)) < epsln * (d__1 = a[p - 
		    1 + (p - 1) * a_dim1], abs(d__1)) || (d__3 = a[p + p * 
		    a_dim1], abs(d__3)) < small || l2kill && (d__4 = a[p + p *
		     a_dim1], abs(d__4)) < temp1) {
#line 967 "dgejsv.f"
		goto L3402;
#line 967 "dgejsv.f"
	    }
#line 970 "dgejsv.f"
	    ++nr;
#line 971 "dgejsv.f"
/* L3401: */
#line 971 "dgejsv.f"
	}
#line 972 "dgejsv.f"
L3402:

#line 974 "dgejsv.f"
	;
#line 974 "dgejsv.f"
    } else {
/*        The goal is high relative accuracy. However, if the matrix */
/*        has high scaled condition number the relative accuracy is in */
/*        general not feasible. Later on, a condition number estimator */
/*        will be deployed to estimate the scaled condition number. */
/*        Here we just remove the underflowed part of the triangular */
/*        factor. This prevents the situation in which the code is */
/*        working hard to get the accuracy not warranted by the data. */
#line 982 "dgejsv.f"
	temp1 = sqrt(sfmin);
#line 983 "dgejsv.f"
	i__1 = *n;
#line 983 "dgejsv.f"
	for (p = 2; p <= i__1; ++p) {
#line 984 "dgejsv.f"
	    if ((d__1 = a[p + p * a_dim1], abs(d__1)) < small || l2kill && (
		    d__2 = a[p + p * a_dim1], abs(d__2)) < temp1) {
#line 984 "dgejsv.f"
		goto L3302;
#line 984 "dgejsv.f"
	    }
#line 986 "dgejsv.f"
	    ++nr;
#line 987 "dgejsv.f"
/* L3301: */
#line 987 "dgejsv.f"
	}
#line 988 "dgejsv.f"
L3302:

#line 990 "dgejsv.f"
	;
#line 990 "dgejsv.f"
    }

#line 992 "dgejsv.f"
    almort = FALSE_;
#line 993 "dgejsv.f"
    if (nr == *n) {
#line 994 "dgejsv.f"
	maxprj = 1.;
#line 995 "dgejsv.f"
	i__1 = *n;
#line 995 "dgejsv.f"
	for (p = 2; p <= i__1; ++p) {
#line 996 "dgejsv.f"
	    temp1 = (d__1 = a[p + p * a_dim1], abs(d__1)) / sva[iwork[p]];
#line 997 "dgejsv.f"
	    maxprj = min(maxprj,temp1);
#line 998 "dgejsv.f"
/* L3051: */
#line 998 "dgejsv.f"
	}
/* Computing 2nd power */
#line 999 "dgejsv.f"
	d__1 = maxprj;
#line 999 "dgejsv.f"
	if (d__1 * d__1 >= 1. - (doublereal) (*n) * epsln) {
#line 999 "dgejsv.f"
	    almort = TRUE_;
#line 999 "dgejsv.f"
	}
#line 1000 "dgejsv.f"
    }


#line 1003 "dgejsv.f"
    sconda = -1.;
#line 1004 "dgejsv.f"
    condr1 = -1.;
#line 1005 "dgejsv.f"
    condr2 = -1.;

#line 1007 "dgejsv.f"
    if (errest) {
#line 1008 "dgejsv.f"
	if (*n == nr) {
#line 1009 "dgejsv.f"
	    if (rsvec) {
/*              .. V is available as workspace */
#line 1011 "dgejsv.f"
		dlacpy_("U", n, n, &a[a_offset], lda, &v[v_offset], ldv, (
			ftnlen)1);
#line 1012 "dgejsv.f"
		i__1 = *n;
#line 1012 "dgejsv.f"
		for (p = 1; p <= i__1; ++p) {
#line 1013 "dgejsv.f"
		    temp1 = sva[iwork[p]];
#line 1014 "dgejsv.f"
		    d__1 = 1. / temp1;
#line 1014 "dgejsv.f"
		    dscal_(&p, &d__1, &v[p * v_dim1 + 1], &c__1);
#line 1015 "dgejsv.f"
/* L3053: */
#line 1015 "dgejsv.f"
		}
#line 1016 "dgejsv.f"
		dpocon_("U", n, &v[v_offset], ldv, &c_b35, &temp1, &work[*n + 
			1], &iwork[(*n << 1) + *m + 1], &ierr, (ftnlen)1);
#line 1018 "dgejsv.f"
	    } else if (lsvec) {
/*              .. U is available as workspace */
#line 1020 "dgejsv.f"
		dlacpy_("U", n, n, &a[a_offset], lda, &u[u_offset], ldu, (
			ftnlen)1);
#line 1021 "dgejsv.f"
		i__1 = *n;
#line 1021 "dgejsv.f"
		for (p = 1; p <= i__1; ++p) {
#line 1022 "dgejsv.f"
		    temp1 = sva[iwork[p]];
#line 1023 "dgejsv.f"
		    d__1 = 1. / temp1;
#line 1023 "dgejsv.f"
		    dscal_(&p, &d__1, &u[p * u_dim1 + 1], &c__1);
#line 1024 "dgejsv.f"
/* L3054: */
#line 1024 "dgejsv.f"
		}
#line 1025 "dgejsv.f"
		dpocon_("U", n, &u[u_offset], ldu, &c_b35, &temp1, &work[*n + 
			1], &iwork[(*n << 1) + *m + 1], &ierr, (ftnlen)1);
#line 1027 "dgejsv.f"
	    } else {
#line 1028 "dgejsv.f"
		dlacpy_("U", n, n, &a[a_offset], lda, &work[*n + 1], n, (
			ftnlen)1);
#line 1029 "dgejsv.f"
		i__1 = *n;
#line 1029 "dgejsv.f"
		for (p = 1; p <= i__1; ++p) {
#line 1030 "dgejsv.f"
		    temp1 = sva[iwork[p]];
#line 1031 "dgejsv.f"
		    d__1 = 1. / temp1;
#line 1031 "dgejsv.f"
		    dscal_(&p, &d__1, &work[*n + (p - 1) * *n + 1], &c__1);
#line 1032 "dgejsv.f"
/* L3052: */
#line 1032 "dgejsv.f"
		}
/*           .. the columns of R are scaled to have unit Euclidean lengths. */
#line 1034 "dgejsv.f"
		dpocon_("U", n, &work[*n + 1], n, &c_b35, &temp1, &work[*n + *
			n * *n + 1], &iwork[(*n << 1) + *m + 1], &ierr, (
			ftnlen)1);
#line 1036 "dgejsv.f"
	    }
#line 1037 "dgejsv.f"
	    sconda = 1. / sqrt(temp1);
/*           SCONDA is an estimate of DSQRT(||(R^t * R)^(-1)||_1). */
/*           N^(-1/4) * SCONDA <= ||R^(-1)||_2 <= N^(1/4) * SCONDA */
#line 1040 "dgejsv.f"
	} else {
#line 1041 "dgejsv.f"
	    sconda = -1.;
#line 1042 "dgejsv.f"
	}
#line 1043 "dgejsv.f"
    }

#line 1045 "dgejsv.f"
    l2pert = l2pert && (d__1 = a[a_dim1 + 1] / a[nr + nr * a_dim1], abs(d__1))
	     > sqrt(big1);
/*     If there is no violent scaling, artificial perturbation is not needed. */

/*     Phase 3: */

#line 1050 "dgejsv.f"
    if (! (rsvec || lsvec)) {

/*         Singular Values only */

/*         .. transpose A(1:NR,1:N) */
/* Computing MIN */
#line 1055 "dgejsv.f"
	i__2 = *n - 1;
#line 1055 "dgejsv.f"
	i__1 = min(i__2,nr);
#line 1055 "dgejsv.f"
	for (p = 1; p <= i__1; ++p) {
#line 1056 "dgejsv.f"
	    i__2 = *n - p;
#line 1056 "dgejsv.f"
	    dcopy_(&i__2, &a[p + (p + 1) * a_dim1], lda, &a[p + 1 + p * 
		    a_dim1], &c__1);
#line 1057 "dgejsv.f"
/* L1946: */
#line 1057 "dgejsv.f"
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

#line 1071 "dgejsv.f"
	if (! almort) {

#line 1073 "dgejsv.f"
	    if (l2pert) {
/*              XSC = DSQRT(SMALL) */
#line 1075 "dgejsv.f"
		xsc = epsln / (doublereal) (*n);
#line 1076 "dgejsv.f"
		i__1 = nr;
#line 1076 "dgejsv.f"
		for (q = 1; q <= i__1; ++q) {
#line 1077 "dgejsv.f"
		    temp1 = xsc * (d__1 = a[q + q * a_dim1], abs(d__1));
#line 1078 "dgejsv.f"
		    i__2 = *n;
#line 1078 "dgejsv.f"
		    for (p = 1; p <= i__2; ++p) {
#line 1079 "dgejsv.f"
			if (p > q && (d__1 = a[p + q * a_dim1], abs(d__1)) <= 
				temp1 || p < q) {
#line 1079 "dgejsv.f"
			    a[p + q * a_dim1] = d_sign(&temp1, &a[p + q * 
				    a_dim1]);
#line 1079 "dgejsv.f"
			}
#line 1082 "dgejsv.f"
/* L4949: */
#line 1082 "dgejsv.f"
		    }
#line 1083 "dgejsv.f"
/* L4947: */
#line 1083 "dgejsv.f"
		}
#line 1084 "dgejsv.f"
	    } else {
#line 1085 "dgejsv.f"
		i__1 = nr - 1;
#line 1085 "dgejsv.f"
		i__2 = nr - 1;
#line 1085 "dgejsv.f"
		dlaset_("U", &i__1, &i__2, &c_b34, &c_b34, &a[(a_dim1 << 1) + 
			1], lda, (ftnlen)1);
#line 1086 "dgejsv.f"
	    }

/*            .. second preconditioning using the QR factorization */

#line 1090 "dgejsv.f"
	    i__1 = *lwork - *n;
#line 1090 "dgejsv.f"
	    dgeqrf_(n, &nr, &a[a_offset], lda, &work[1], &work[*n + 1], &i__1,
		     &ierr);

/*           .. and transpose upper to lower triangular */
#line 1093 "dgejsv.f"
	    i__1 = nr - 1;
#line 1093 "dgejsv.f"
	    for (p = 1; p <= i__1; ++p) {
#line 1094 "dgejsv.f"
		i__2 = nr - p;
#line 1094 "dgejsv.f"
		dcopy_(&i__2, &a[p + (p + 1) * a_dim1], lda, &a[p + 1 + p * 
			a_dim1], &c__1);
#line 1095 "dgejsv.f"
/* L1948: */
#line 1095 "dgejsv.f"
	    }

#line 1097 "dgejsv.f"
	}

/*           Row-cyclic Jacobi SVD algorithm with column pivoting */

/*           .. again some perturbation (a "background noise") is added */
/*           to drown denormals */
#line 1103 "dgejsv.f"
	if (l2pert) {
/*              XSC = DSQRT(SMALL) */
#line 1105 "dgejsv.f"
	    xsc = epsln / (doublereal) (*n);
#line 1106 "dgejsv.f"
	    i__1 = nr;
#line 1106 "dgejsv.f"
	    for (q = 1; q <= i__1; ++q) {
#line 1107 "dgejsv.f"
		temp1 = xsc * (d__1 = a[q + q * a_dim1], abs(d__1));
#line 1108 "dgejsv.f"
		i__2 = nr;
#line 1108 "dgejsv.f"
		for (p = 1; p <= i__2; ++p) {
#line 1109 "dgejsv.f"
		    if (p > q && (d__1 = a[p + q * a_dim1], abs(d__1)) <= 
			    temp1 || p < q) {
#line 1109 "dgejsv.f"
			a[p + q * a_dim1] = d_sign(&temp1, &a[p + q * a_dim1])
				;
#line 1109 "dgejsv.f"
		    }
#line 1112 "dgejsv.f"
/* L1949: */
#line 1112 "dgejsv.f"
		}
#line 1113 "dgejsv.f"
/* L1947: */
#line 1113 "dgejsv.f"
	    }
#line 1114 "dgejsv.f"
	} else {
#line 1115 "dgejsv.f"
	    i__1 = nr - 1;
#line 1115 "dgejsv.f"
	    i__2 = nr - 1;
#line 1115 "dgejsv.f"
	    dlaset_("U", &i__1, &i__2, &c_b34, &c_b34, &a[(a_dim1 << 1) + 1], 
		    lda, (ftnlen)1);
#line 1116 "dgejsv.f"
	}

/*           .. and one-sided Jacobi rotations are started on a lower */
/*           triangular matrix (plus perturbation which is ignored in */
/*           the part which destroys triangular form (confusing?!)) */

#line 1122 "dgejsv.f"
	dgesvj_("L", "NoU", "NoV", &nr, &nr, &a[a_offset], lda, &sva[1], n, &
		v[v_offset], ldv, &work[1], lwork, info, (ftnlen)1, (ftnlen)3,
		 (ftnlen)3);

#line 1125 "dgejsv.f"
	scalem = work[1];
#line 1126 "dgejsv.f"
	numrank = i_dnnt(&work[2]);


#line 1129 "dgejsv.f"
    } else if (rsvec && ! lsvec) {

/*        -> Singular Values and Right Singular Vectors <- */

#line 1133 "dgejsv.f"
	if (almort) {

/*           .. in this case NR equals N */
#line 1136 "dgejsv.f"
	    i__1 = nr;
#line 1136 "dgejsv.f"
	    for (p = 1; p <= i__1; ++p) {
#line 1137 "dgejsv.f"
		i__2 = *n - p + 1;
#line 1137 "dgejsv.f"
		dcopy_(&i__2, &a[p + p * a_dim1], lda, &v[p + p * v_dim1], &
			c__1);
#line 1138 "dgejsv.f"
/* L1998: */
#line 1138 "dgejsv.f"
	    }
#line 1139 "dgejsv.f"
	    i__1 = nr - 1;
#line 1139 "dgejsv.f"
	    i__2 = nr - 1;
#line 1139 "dgejsv.f"
	    dlaset_("Upper", &i__1, &i__2, &c_b34, &c_b34, &v[(v_dim1 << 1) + 
		    1], ldv, (ftnlen)5);

#line 1141 "dgejsv.f"
	    dgesvj_("L", "U", "N", n, &nr, &v[v_offset], ldv, &sva[1], &nr, &
		    a[a_offset], lda, &work[1], lwork, info, (ftnlen)1, (
		    ftnlen)1, (ftnlen)1);
#line 1143 "dgejsv.f"
	    scalem = work[1];
#line 1144 "dgejsv.f"
	    numrank = i_dnnt(&work[2]);
#line 1146 "dgejsv.f"
	} else {

/*        .. two more QR factorizations ( one QRF is not enough, two require */
/*        accumulated product of Jacobi rotations, three are perfect ) */

#line 1151 "dgejsv.f"
	    i__1 = nr - 1;
#line 1151 "dgejsv.f"
	    i__2 = nr - 1;
#line 1151 "dgejsv.f"
	    dlaset_("Lower", &i__1, &i__2, &c_b34, &c_b34, &a[a_dim1 + 2], 
		    lda, (ftnlen)5);
#line 1152 "dgejsv.f"
	    i__1 = *lwork - *n;
#line 1152 "dgejsv.f"
	    dgelqf_(&nr, n, &a[a_offset], lda, &work[1], &work[*n + 1], &i__1,
		     &ierr);
#line 1153 "dgejsv.f"
	    dlacpy_("Lower", &nr, &nr, &a[a_offset], lda, &v[v_offset], ldv, (
		    ftnlen)5);
#line 1154 "dgejsv.f"
	    i__1 = nr - 1;
#line 1154 "dgejsv.f"
	    i__2 = nr - 1;
#line 1154 "dgejsv.f"
	    dlaset_("Upper", &i__1, &i__2, &c_b34, &c_b34, &v[(v_dim1 << 1) + 
		    1], ldv, (ftnlen)5);
#line 1155 "dgejsv.f"
	    i__1 = *lwork - (*n << 1);
#line 1155 "dgejsv.f"
	    dgeqrf_(&nr, &nr, &v[v_offset], ldv, &work[*n + 1], &work[(*n << 
		    1) + 1], &i__1, &ierr);
#line 1157 "dgejsv.f"
	    i__1 = nr;
#line 1157 "dgejsv.f"
	    for (p = 1; p <= i__1; ++p) {
#line 1158 "dgejsv.f"
		i__2 = nr - p + 1;
#line 1158 "dgejsv.f"
		dcopy_(&i__2, &v[p + p * v_dim1], ldv, &v[p + p * v_dim1], &
			c__1);
#line 1159 "dgejsv.f"
/* L8998: */
#line 1159 "dgejsv.f"
	    }
#line 1160 "dgejsv.f"
	    i__1 = nr - 1;
#line 1160 "dgejsv.f"
	    i__2 = nr - 1;
#line 1160 "dgejsv.f"
	    dlaset_("Upper", &i__1, &i__2, &c_b34, &c_b34, &v[(v_dim1 << 1) + 
		    1], ldv, (ftnlen)5);

#line 1162 "dgejsv.f"
	    dgesvj_("Lower", "U", "N", &nr, &nr, &v[v_offset], ldv, &sva[1], &
		    nr, &u[u_offset], ldu, &work[*n + 1], lwork, info, (
		    ftnlen)5, (ftnlen)1, (ftnlen)1);
#line 1164 "dgejsv.f"
	    scalem = work[*n + 1];
#line 1165 "dgejsv.f"
	    numrank = i_dnnt(&work[*n + 2]);
#line 1166 "dgejsv.f"
	    if (nr < *n) {
#line 1167 "dgejsv.f"
		i__1 = *n - nr;
#line 1167 "dgejsv.f"
		dlaset_("A", &i__1, &nr, &c_b34, &c_b34, &v[nr + 1 + v_dim1], 
			ldv, (ftnlen)1);
#line 1168 "dgejsv.f"
		i__1 = *n - nr;
#line 1168 "dgejsv.f"
		dlaset_("A", &nr, &i__1, &c_b34, &c_b34, &v[(nr + 1) * v_dim1 
			+ 1], ldv, (ftnlen)1);
#line 1169 "dgejsv.f"
		i__1 = *n - nr;
#line 1169 "dgejsv.f"
		i__2 = *n - nr;
#line 1169 "dgejsv.f"
		dlaset_("A", &i__1, &i__2, &c_b34, &c_b35, &v[nr + 1 + (nr + 
			1) * v_dim1], ldv, (ftnlen)1);
#line 1170 "dgejsv.f"
	    }

#line 1172 "dgejsv.f"
	    i__1 = *lwork - *n;
#line 1172 "dgejsv.f"
	    dormlq_("Left", "Transpose", n, n, &nr, &a[a_offset], lda, &work[
		    1], &v[v_offset], ldv, &work[*n + 1], &i__1, &ierr, (
		    ftnlen)4, (ftnlen)9);

#line 1175 "dgejsv.f"
	}

#line 1177 "dgejsv.f"
	i__1 = *n;
#line 1177 "dgejsv.f"
	for (p = 1; p <= i__1; ++p) {
#line 1178 "dgejsv.f"
	    dcopy_(n, &v[p + v_dim1], ldv, &a[iwork[p] + a_dim1], lda);
#line 1179 "dgejsv.f"
/* L8991: */
#line 1179 "dgejsv.f"
	}
#line 1180 "dgejsv.f"
	dlacpy_("All", n, n, &a[a_offset], lda, &v[v_offset], ldv, (ftnlen)3);

#line 1182 "dgejsv.f"
	if (transp) {
#line 1183 "dgejsv.f"
	    dlacpy_("All", n, n, &v[v_offset], ldv, &u[u_offset], ldu, (
		    ftnlen)3);
#line 1184 "dgejsv.f"
	}

#line 1186 "dgejsv.f"
    } else if (lsvec && ! rsvec) {

/*        .. Singular Values and Left Singular Vectors                 .. */

/*        .. second preconditioning step to avoid need to accumulate */
/*        Jacobi rotations in the Jacobi iterations. */
#line 1192 "dgejsv.f"
	i__1 = nr;
#line 1192 "dgejsv.f"
	for (p = 1; p <= i__1; ++p) {
#line 1193 "dgejsv.f"
	    i__2 = *n - p + 1;
#line 1193 "dgejsv.f"
	    dcopy_(&i__2, &a[p + p * a_dim1], lda, &u[p + p * u_dim1], &c__1);
#line 1194 "dgejsv.f"
/* L1965: */
#line 1194 "dgejsv.f"
	}
#line 1195 "dgejsv.f"
	i__1 = nr - 1;
#line 1195 "dgejsv.f"
	i__2 = nr - 1;
#line 1195 "dgejsv.f"
	dlaset_("Upper", &i__1, &i__2, &c_b34, &c_b34, &u[(u_dim1 << 1) + 1], 
		ldu, (ftnlen)5);

#line 1197 "dgejsv.f"
	i__1 = *lwork - (*n << 1);
#line 1197 "dgejsv.f"
	dgeqrf_(n, &nr, &u[u_offset], ldu, &work[*n + 1], &work[(*n << 1) + 1]
		, &i__1, &ierr);

#line 1200 "dgejsv.f"
	i__1 = nr - 1;
#line 1200 "dgejsv.f"
	for (p = 1; p <= i__1; ++p) {
#line 1201 "dgejsv.f"
	    i__2 = nr - p;
#line 1201 "dgejsv.f"
	    dcopy_(&i__2, &u[p + (p + 1) * u_dim1], ldu, &u[p + 1 + p * 
		    u_dim1], &c__1);
#line 1202 "dgejsv.f"
/* L1967: */
#line 1202 "dgejsv.f"
	}
#line 1203 "dgejsv.f"
	i__1 = nr - 1;
#line 1203 "dgejsv.f"
	i__2 = nr - 1;
#line 1203 "dgejsv.f"
	dlaset_("Upper", &i__1, &i__2, &c_b34, &c_b34, &u[(u_dim1 << 1) + 1], 
		ldu, (ftnlen)5);

#line 1205 "dgejsv.f"
	i__1 = *lwork - *n;
#line 1205 "dgejsv.f"
	dgesvj_("Lower", "U", "N", &nr, &nr, &u[u_offset], ldu, &sva[1], &nr, 
		&a[a_offset], lda, &work[*n + 1], &i__1, info, (ftnlen)5, (
		ftnlen)1, (ftnlen)1);
#line 1207 "dgejsv.f"
	scalem = work[*n + 1];
#line 1208 "dgejsv.f"
	numrank = i_dnnt(&work[*n + 2]);

#line 1210 "dgejsv.f"
	if (nr < *m) {
#line 1211 "dgejsv.f"
	    i__1 = *m - nr;
#line 1211 "dgejsv.f"
	    dlaset_("A", &i__1, &nr, &c_b34, &c_b34, &u[nr + 1 + u_dim1], ldu,
		     (ftnlen)1);
#line 1212 "dgejsv.f"
	    if (nr < n1) {
#line 1213 "dgejsv.f"
		i__1 = n1 - nr;
#line 1213 "dgejsv.f"
		dlaset_("A", &nr, &i__1, &c_b34, &c_b34, &u[(nr + 1) * u_dim1 
			+ 1], ldu, (ftnlen)1);
#line 1214 "dgejsv.f"
		i__1 = *m - nr;
#line 1214 "dgejsv.f"
		i__2 = n1 - nr;
#line 1214 "dgejsv.f"
		dlaset_("A", &i__1, &i__2, &c_b34, &c_b35, &u[nr + 1 + (nr + 
			1) * u_dim1], ldu, (ftnlen)1);
#line 1215 "dgejsv.f"
	    }
#line 1216 "dgejsv.f"
	}

#line 1218 "dgejsv.f"
	i__1 = *lwork - *n;
#line 1218 "dgejsv.f"
	dormqr_("Left", "No Tr", m, &n1, n, &a[a_offset], lda, &work[1], &u[
		u_offset], ldu, &work[*n + 1], &i__1, &ierr, (ftnlen)4, (
		ftnlen)5);

#line 1221 "dgejsv.f"
	if (rowpiv) {
#line 1221 "dgejsv.f"
	    i__1 = *m - 1;
#line 1221 "dgejsv.f"
	    dlaswp_(&n1, &u[u_offset], ldu, &c__1, &i__1, &iwork[(*n << 1) + 
		    1], &c_n1);
#line 1221 "dgejsv.f"
	}

#line 1224 "dgejsv.f"
	i__1 = n1;
#line 1224 "dgejsv.f"
	for (p = 1; p <= i__1; ++p) {
#line 1225 "dgejsv.f"
	    xsc = 1. / dnrm2_(m, &u[p * u_dim1 + 1], &c__1);
#line 1226 "dgejsv.f"
	    dscal_(m, &xsc, &u[p * u_dim1 + 1], &c__1);
#line 1227 "dgejsv.f"
/* L1974: */
#line 1227 "dgejsv.f"
	}

#line 1229 "dgejsv.f"
	if (transp) {
#line 1230 "dgejsv.f"
	    dlacpy_("All", n, n, &u[u_offset], ldu, &v[v_offset], ldv, (
		    ftnlen)3);
#line 1231 "dgejsv.f"
	}

#line 1233 "dgejsv.f"
    } else {

/*        .. Full SVD .. */

#line 1237 "dgejsv.f"
	if (! jracc) {

#line 1239 "dgejsv.f"
	    if (! almort) {

/*           Second Preconditioning Step (QRF [with pivoting]) */
/*           Note that the composition of TRANSPOSE, QRF and TRANSPOSE is */
/*           equivalent to an LQF CALL. Since in many libraries the QRF */
/*           seems to be better optimized than the LQF, we do explicit */
/*           transpose and use the QRF. This is subject to changes in an */
/*           optimized implementation of DGEJSV. */

#line 1248 "dgejsv.f"
		i__1 = nr;
#line 1248 "dgejsv.f"
		for (p = 1; p <= i__1; ++p) {
#line 1249 "dgejsv.f"
		    i__2 = *n - p + 1;
#line 1249 "dgejsv.f"
		    dcopy_(&i__2, &a[p + p * a_dim1], lda, &v[p + p * v_dim1],
			     &c__1);
#line 1250 "dgejsv.f"
/* L1968: */
#line 1250 "dgejsv.f"
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

#line 1264 "dgejsv.f"
		if (l2pert) {
#line 1265 "dgejsv.f"
		    xsc = sqrt(small);
#line 1266 "dgejsv.f"
		    i__1 = nr;
#line 1266 "dgejsv.f"
		    for (q = 1; q <= i__1; ++q) {
#line 1267 "dgejsv.f"
			temp1 = xsc * (d__1 = v[q + q * v_dim1], abs(d__1));
#line 1268 "dgejsv.f"
			i__2 = *n;
#line 1268 "dgejsv.f"
			for (p = 1; p <= i__2; ++p) {
#line 1269 "dgejsv.f"
			    if (p > q && (d__1 = v[p + q * v_dim1], abs(d__1))
				     <= temp1 || p < q) {
#line 1269 "dgejsv.f"
				v[p + q * v_dim1] = d_sign(&temp1, &v[p + q * 
					v_dim1]);
#line 1269 "dgejsv.f"
			    }
#line 1272 "dgejsv.f"
			    if (p < q) {
#line 1272 "dgejsv.f"
				v[p + q * v_dim1] = -v[p + q * v_dim1];
#line 1272 "dgejsv.f"
			    }
#line 1273 "dgejsv.f"
/* L2968: */
#line 1273 "dgejsv.f"
			}
#line 1274 "dgejsv.f"
/* L2969: */
#line 1274 "dgejsv.f"
		    }
#line 1275 "dgejsv.f"
		} else {
#line 1276 "dgejsv.f"
		    i__1 = nr - 1;
#line 1276 "dgejsv.f"
		    i__2 = nr - 1;
#line 1276 "dgejsv.f"
		    dlaset_("U", &i__1, &i__2, &c_b34, &c_b34, &v[(v_dim1 << 
			    1) + 1], ldv, (ftnlen)1);
#line 1277 "dgejsv.f"
		}

/*           Estimate the row scaled condition number of R1 */
/*           (If R1 is rectangular, N > NR, then the condition number */
/*           of the leading NR x NR submatrix is estimated.) */

#line 1283 "dgejsv.f"
		dlacpy_("L", &nr, &nr, &v[v_offset], ldv, &work[(*n << 1) + 1]
			, &nr, (ftnlen)1);
#line 1284 "dgejsv.f"
		i__1 = nr;
#line 1284 "dgejsv.f"
		for (p = 1; p <= i__1; ++p) {
#line 1285 "dgejsv.f"
		    i__2 = nr - p + 1;
#line 1285 "dgejsv.f"
		    temp1 = dnrm2_(&i__2, &work[(*n << 1) + (p - 1) * nr + p],
			     &c__1);
#line 1286 "dgejsv.f"
		    i__2 = nr - p + 1;
#line 1286 "dgejsv.f"
		    d__1 = 1. / temp1;
#line 1286 "dgejsv.f"
		    dscal_(&i__2, &d__1, &work[(*n << 1) + (p - 1) * nr + p], 
			    &c__1);
#line 1287 "dgejsv.f"
/* L3950: */
#line 1287 "dgejsv.f"
		}
#line 1288 "dgejsv.f"
		dpocon_("Lower", &nr, &work[(*n << 1) + 1], &nr, &c_b35, &
			temp1, &work[(*n << 1) + nr * nr + 1], &iwork[*m + (*
			n << 1) + 1], &ierr, (ftnlen)5);
#line 1290 "dgejsv.f"
		condr1 = 1. / sqrt(temp1);
/*           .. here need a second oppinion on the condition number */
/*           .. then assume worst case scenario */
/*           R1 is OK for inverse <=> CONDR1 .LT. DBLE(N) */
/*           more conservative    <=> CONDR1 .LT. DSQRT(DBLE(N)) */

#line 1296 "dgejsv.f"
		cond_ok__ = sqrt((doublereal) nr);
/* [TP]       COND_OK is a tuning parameter. */
#line 1299 "dgejsv.f"
		if (condr1 < cond_ok__) {
/*              .. the second QRF without pivoting. Note: in an optimized */
/*              implementation, this QRF should be implemented as the QRF */
/*              of a lower triangular matrix. */
/*              R1^t = Q2 * R2 */
#line 1304 "dgejsv.f"
		    i__1 = *lwork - (*n << 1);
#line 1304 "dgejsv.f"
		    dgeqrf_(n, &nr, &v[v_offset], ldv, &work[*n + 1], &work[(*
			    n << 1) + 1], &i__1, &ierr);

#line 1307 "dgejsv.f"
		    if (l2pert) {
#line 1308 "dgejsv.f"
			xsc = sqrt(small) / epsln;
#line 1309 "dgejsv.f"
			i__1 = nr;
#line 1309 "dgejsv.f"
			for (p = 2; p <= i__1; ++p) {
#line 1310 "dgejsv.f"
			    i__2 = p - 1;
#line 1310 "dgejsv.f"
			    for (q = 1; q <= i__2; ++q) {
/* Computing MIN */
#line 1311 "dgejsv.f"
				d__3 = (d__1 = v[p + p * v_dim1], abs(d__1)), 
					d__4 = (d__2 = v[q + q * v_dim1], abs(
					d__2));
#line 1311 "dgejsv.f"
				temp1 = xsc * min(d__3,d__4);
#line 1312 "dgejsv.f"
				if ((d__1 = v[q + p * v_dim1], abs(d__1)) <= 
					temp1) {
#line 1312 "dgejsv.f"
				    v[q + p * v_dim1] = d_sign(&temp1, &v[q + 
					    p * v_dim1]);
#line 1312 "dgejsv.f"
				}
#line 1314 "dgejsv.f"
/* L3958: */
#line 1314 "dgejsv.f"
			    }
#line 1315 "dgejsv.f"
/* L3959: */
#line 1315 "dgejsv.f"
			}
#line 1316 "dgejsv.f"
		    }

#line 1318 "dgejsv.f"
		    if (nr != *n) {
#line 1318 "dgejsv.f"
			dlacpy_("A", n, &nr, &v[v_offset], ldv, &work[(*n << 
				1) + 1], n, (ftnlen)1);
#line 1318 "dgejsv.f"
		    }
/*              .. save ... */

/*           .. this transposed copy should be better than naive */
#line 1323 "dgejsv.f"
		    i__1 = nr - 1;
#line 1323 "dgejsv.f"
		    for (p = 1; p <= i__1; ++p) {
#line 1324 "dgejsv.f"
			i__2 = nr - p;
#line 1324 "dgejsv.f"
			dcopy_(&i__2, &v[p + (p + 1) * v_dim1], ldv, &v[p + 1 
				+ p * v_dim1], &c__1);
#line 1325 "dgejsv.f"
/* L1969: */
#line 1325 "dgejsv.f"
		    }

#line 1327 "dgejsv.f"
		    condr2 = condr1;

#line 1329 "dgejsv.f"
		} else {

/*              .. ill-conditioned case: second QRF with pivoting */
/*              Note that windowed pivoting would be equaly good */
/*              numerically, and more run-time efficient. So, in */
/*              an optimal implementation, the next call to DGEQP3 */
/*              should be replaced with eg. CALL SGEQPX (ACM TOMS #782) */
/*              with properly (carefully) chosen parameters. */

/*              R1^t * P2 = Q2 * R2 */
#line 1339 "dgejsv.f"
		    i__1 = nr;
#line 1339 "dgejsv.f"
		    for (p = 1; p <= i__1; ++p) {
#line 1340 "dgejsv.f"
			iwork[*n + p] = 0;
#line 1341 "dgejsv.f"
/* L3003: */
#line 1341 "dgejsv.f"
		    }
#line 1342 "dgejsv.f"
		    i__1 = *lwork - (*n << 1);
#line 1342 "dgejsv.f"
		    dgeqp3_(n, &nr, &v[v_offset], ldv, &iwork[*n + 1], &work[*
			    n + 1], &work[(*n << 1) + 1], &i__1, &ierr);
/* *               CALL DGEQRF( N, NR, V, LDV, WORK(N+1), WORK(2*N+1), */
/* *     $              LWORK-2*N, IERR ) */
#line 1346 "dgejsv.f"
		    if (l2pert) {
#line 1347 "dgejsv.f"
			xsc = sqrt(small);
#line 1348 "dgejsv.f"
			i__1 = nr;
#line 1348 "dgejsv.f"
			for (p = 2; p <= i__1; ++p) {
#line 1349 "dgejsv.f"
			    i__2 = p - 1;
#line 1349 "dgejsv.f"
			    for (q = 1; q <= i__2; ++q) {
/* Computing MIN */
#line 1350 "dgejsv.f"
				d__3 = (d__1 = v[p + p * v_dim1], abs(d__1)), 
					d__4 = (d__2 = v[q + q * v_dim1], abs(
					d__2));
#line 1350 "dgejsv.f"
				temp1 = xsc * min(d__3,d__4);
#line 1351 "dgejsv.f"
				if ((d__1 = v[q + p * v_dim1], abs(d__1)) <= 
					temp1) {
#line 1351 "dgejsv.f"
				    v[q + p * v_dim1] = d_sign(&temp1, &v[q + 
					    p * v_dim1]);
#line 1351 "dgejsv.f"
				}
#line 1353 "dgejsv.f"
/* L3968: */
#line 1353 "dgejsv.f"
			    }
#line 1354 "dgejsv.f"
/* L3969: */
#line 1354 "dgejsv.f"
			}
#line 1355 "dgejsv.f"
		    }

#line 1357 "dgejsv.f"
		    dlacpy_("A", n, &nr, &v[v_offset], ldv, &work[(*n << 1) + 
			    1], n, (ftnlen)1);

#line 1359 "dgejsv.f"
		    if (l2pert) {
#line 1360 "dgejsv.f"
			xsc = sqrt(small);
#line 1361 "dgejsv.f"
			i__1 = nr;
#line 1361 "dgejsv.f"
			for (p = 2; p <= i__1; ++p) {
#line 1362 "dgejsv.f"
			    i__2 = p - 1;
#line 1362 "dgejsv.f"
			    for (q = 1; q <= i__2; ++q) {
/* Computing MIN */
#line 1363 "dgejsv.f"
				d__3 = (d__1 = v[p + p * v_dim1], abs(d__1)), 
					d__4 = (d__2 = v[q + q * v_dim1], abs(
					d__2));
#line 1363 "dgejsv.f"
				temp1 = xsc * min(d__3,d__4);
#line 1364 "dgejsv.f"
				v[p + q * v_dim1] = -d_sign(&temp1, &v[q + p *
					 v_dim1]);
#line 1365 "dgejsv.f"
/* L8971: */
#line 1365 "dgejsv.f"
			    }
#line 1366 "dgejsv.f"
/* L8970: */
#line 1366 "dgejsv.f"
			}
#line 1367 "dgejsv.f"
		    } else {
#line 1368 "dgejsv.f"
			i__1 = nr - 1;
#line 1368 "dgejsv.f"
			i__2 = nr - 1;
#line 1368 "dgejsv.f"
			dlaset_("L", &i__1, &i__2, &c_b34, &c_b34, &v[v_dim1 
				+ 2], ldv, (ftnlen)1);
#line 1369 "dgejsv.f"
		    }
/*              Now, compute R2 = L3 * Q3, the LQ factorization. */
#line 1371 "dgejsv.f"
		    i__1 = *lwork - (*n << 1) - *n * nr - nr;
#line 1371 "dgejsv.f"
		    dgelqf_(&nr, &nr, &v[v_offset], ldv, &work[(*n << 1) + *n 
			    * nr + 1], &work[(*n << 1) + *n * nr + nr + 1], &
			    i__1, &ierr);
/*              .. and estimate the condition number */
#line 1374 "dgejsv.f"
		    dlacpy_("L", &nr, &nr, &v[v_offset], ldv, &work[(*n << 1) 
			    + *n * nr + nr + 1], &nr, (ftnlen)1);
#line 1375 "dgejsv.f"
		    i__1 = nr;
#line 1375 "dgejsv.f"
		    for (p = 1; p <= i__1; ++p) {
#line 1376 "dgejsv.f"
			temp1 = dnrm2_(&p, &work[(*n << 1) + *n * nr + nr + p]
				, &nr);
#line 1377 "dgejsv.f"
			d__1 = 1. / temp1;
#line 1377 "dgejsv.f"
			dscal_(&p, &d__1, &work[(*n << 1) + *n * nr + nr + p],
				 &nr);
#line 1378 "dgejsv.f"
/* L4950: */
#line 1378 "dgejsv.f"
		    }
#line 1379 "dgejsv.f"
		    dpocon_("L", &nr, &work[(*n << 1) + *n * nr + nr + 1], &
			    nr, &c_b35, &temp1, &work[(*n << 1) + *n * nr + 
			    nr + nr * nr + 1], &iwork[*m + (*n << 1) + 1], &
			    ierr, (ftnlen)1);
#line 1381 "dgejsv.f"
		    condr2 = 1. / sqrt(temp1);

#line 1383 "dgejsv.f"
		    if (condr2 >= cond_ok__) {
/*                 .. save the Householder vectors used for Q3 */
/*                 (this overwrittes the copy of R2, as it will not be */
/*                 needed in this branch, but it does not overwritte the */
/*                 Huseholder vectors of Q2.). */
#line 1388 "dgejsv.f"
			dlacpy_("U", &nr, &nr, &v[v_offset], ldv, &work[(*n <<
				 1) + 1], n, (ftnlen)1);
/*                 .. and the rest of the information on Q3 is in */
/*                 WORK(2*N+N*NR+1:2*N+N*NR+N) */
#line 1391 "dgejsv.f"
		    }

#line 1393 "dgejsv.f"
		}

#line 1395 "dgejsv.f"
		if (l2pert) {
#line 1396 "dgejsv.f"
		    xsc = sqrt(small);
#line 1397 "dgejsv.f"
		    i__1 = nr;
#line 1397 "dgejsv.f"
		    for (q = 2; q <= i__1; ++q) {
#line 1398 "dgejsv.f"
			temp1 = xsc * v[q + q * v_dim1];
#line 1399 "dgejsv.f"
			i__2 = q - 1;
#line 1399 "dgejsv.f"
			for (p = 1; p <= i__2; ++p) {
/*                    V(p,q) = - DSIGN( TEMP1, V(q,p) ) */
#line 1401 "dgejsv.f"
			    v[p + q * v_dim1] = -d_sign(&temp1, &v[p + q * 
				    v_dim1]);
#line 1402 "dgejsv.f"
/* L4969: */
#line 1402 "dgejsv.f"
			}
#line 1403 "dgejsv.f"
/* L4968: */
#line 1403 "dgejsv.f"
		    }
#line 1404 "dgejsv.f"
		} else {
#line 1405 "dgejsv.f"
		    i__1 = nr - 1;
#line 1405 "dgejsv.f"
		    i__2 = nr - 1;
#line 1405 "dgejsv.f"
		    dlaset_("U", &i__1, &i__2, &c_b34, &c_b34, &v[(v_dim1 << 
			    1) + 1], ldv, (ftnlen)1);
#line 1406 "dgejsv.f"
		}

/*        Second preconditioning finished; continue with Jacobi SVD */
/*        The input matrix is lower trinagular. */

/*        Recover the right singular vectors as solution of a well */
/*        conditioned triangular matrix equation. */

#line 1414 "dgejsv.f"
		if (condr1 < cond_ok__) {

#line 1416 "dgejsv.f"
		    i__1 = *lwork - (*n << 1) - *n * nr - nr;
#line 1416 "dgejsv.f"
		    dgesvj_("L", "U", "N", &nr, &nr, &v[v_offset], ldv, &sva[
			    1], &nr, &u[u_offset], ldu, &work[(*n << 1) + *n *
			     nr + nr + 1], &i__1, info, (ftnlen)1, (ftnlen)1, 
			    (ftnlen)1);
#line 1418 "dgejsv.f"
		    scalem = work[(*n << 1) + *n * nr + nr + 1];
#line 1419 "dgejsv.f"
		    numrank = i_dnnt(&work[(*n << 1) + *n * nr + nr + 2]);
#line 1420 "dgejsv.f"
		    i__1 = nr;
#line 1420 "dgejsv.f"
		    for (p = 1; p <= i__1; ++p) {
#line 1421 "dgejsv.f"
			dcopy_(&nr, &v[p * v_dim1 + 1], &c__1, &u[p * u_dim1 
				+ 1], &c__1);
#line 1422 "dgejsv.f"
			dscal_(&nr, &sva[p], &v[p * v_dim1 + 1], &c__1);
#line 1423 "dgejsv.f"
/* L3970: */
#line 1423 "dgejsv.f"
		    }
/*        .. pick the right matrix equation and solve it */

#line 1427 "dgejsv.f"
		    if (nr == *n) {
/* :))             .. best case, R1 is inverted. The solution of this matrix */
/*                 equation is Q2*V2 = the product of the Jacobi rotations */
/*                 used in DGESVJ, premultiplied with the orthogonal matrix */
/*                 from the second QR factorization. */
#line 1432 "dgejsv.f"
			dtrsm_("L", "U", "N", "N", &nr, &nr, &c_b35, &a[
				a_offset], lda, &v[v_offset], ldv, (ftnlen)1, 
				(ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1433 "dgejsv.f"
		    } else {
/*                 .. R1 is well conditioned, but non-square. Transpose(R2) */
/*                 is inverted to get the product of the Jacobi rotations */
/*                 used in DGESVJ. The Q-factor from the second QR */
/*                 factorization is then built in explicitly. */
#line 1438 "dgejsv.f"
			dtrsm_("L", "U", "T", "N", &nr, &nr, &c_b35, &work[(*
				n << 1) + 1], n, &v[v_offset], ldv, (ftnlen)1,
				 (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1440 "dgejsv.f"
			if (nr < *n) {
#line 1441 "dgejsv.f"
			    i__1 = *n - nr;
#line 1441 "dgejsv.f"
			    dlaset_("A", &i__1, &nr, &c_b34, &c_b34, &v[nr + 
				    1 + v_dim1], ldv, (ftnlen)1);
#line 1442 "dgejsv.f"
			    i__1 = *n - nr;
#line 1442 "dgejsv.f"
			    dlaset_("A", &nr, &i__1, &c_b34, &c_b34, &v[(nr + 
				    1) * v_dim1 + 1], ldv, (ftnlen)1);
#line 1443 "dgejsv.f"
			    i__1 = *n - nr;
#line 1443 "dgejsv.f"
			    i__2 = *n - nr;
#line 1443 "dgejsv.f"
			    dlaset_("A", &i__1, &i__2, &c_b34, &c_b35, &v[nr 
				    + 1 + (nr + 1) * v_dim1], ldv, (ftnlen)1);
#line 1444 "dgejsv.f"
			}
#line 1445 "dgejsv.f"
			i__1 = *lwork - (*n << 1) - *n * nr - nr;
#line 1445 "dgejsv.f"
			dormqr_("L", "N", n, n, &nr, &work[(*n << 1) + 1], n, 
				&work[*n + 1], &v[v_offset], ldv, &work[(*n <<
				 1) + *n * nr + nr + 1], &i__1, &ierr, (
				ftnlen)1, (ftnlen)1);
#line 1447 "dgejsv.f"
		    }

#line 1449 "dgejsv.f"
		} else if (condr2 < cond_ok__) {

/* :)           .. the input matrix A is very likely a relative of */
/*              the Kahan matrix :) */
/*              The matrix R2 is inverted. The solution of the matrix equation */
/*              is Q3^T*V3 = the product of the Jacobi rotations (appplied to */
/*              the lower triangular L3 from the LQ factorization of */
/*              R2=L3*Q3), pre-multiplied with the transposed Q3. */
#line 1457 "dgejsv.f"
		    i__1 = *lwork - (*n << 1) - *n * nr - nr;
#line 1457 "dgejsv.f"
		    dgesvj_("L", "U", "N", &nr, &nr, &v[v_offset], ldv, &sva[
			    1], &nr, &u[u_offset], ldu, &work[(*n << 1) + *n *
			     nr + nr + 1], &i__1, info, (ftnlen)1, (ftnlen)1, 
			    (ftnlen)1);
#line 1459 "dgejsv.f"
		    scalem = work[(*n << 1) + *n * nr + nr + 1];
#line 1460 "dgejsv.f"
		    numrank = i_dnnt(&work[(*n << 1) + *n * nr + nr + 2]);
#line 1461 "dgejsv.f"
		    i__1 = nr;
#line 1461 "dgejsv.f"
		    for (p = 1; p <= i__1; ++p) {
#line 1462 "dgejsv.f"
			dcopy_(&nr, &v[p * v_dim1 + 1], &c__1, &u[p * u_dim1 
				+ 1], &c__1);
#line 1463 "dgejsv.f"
			dscal_(&nr, &sva[p], &u[p * u_dim1 + 1], &c__1);
#line 1464 "dgejsv.f"
/* L3870: */
#line 1464 "dgejsv.f"
		    }
#line 1465 "dgejsv.f"
		    dtrsm_("L", "U", "N", "N", &nr, &nr, &c_b35, &work[(*n << 
			    1) + 1], n, &u[u_offset], ldu, (ftnlen)1, (ftnlen)
			    1, (ftnlen)1, (ftnlen)1);
/*              .. apply the permutation from the second QR factorization */
#line 1467 "dgejsv.f"
		    i__1 = nr;
#line 1467 "dgejsv.f"
		    for (q = 1; q <= i__1; ++q) {
#line 1468 "dgejsv.f"
			i__2 = nr;
#line 1468 "dgejsv.f"
			for (p = 1; p <= i__2; ++p) {
#line 1469 "dgejsv.f"
			    work[(*n << 1) + *n * nr + nr + iwork[*n + p]] = 
				    u[p + q * u_dim1];
#line 1470 "dgejsv.f"
/* L872: */
#line 1470 "dgejsv.f"
			}
#line 1471 "dgejsv.f"
			i__2 = nr;
#line 1471 "dgejsv.f"
			for (p = 1; p <= i__2; ++p) {
#line 1472 "dgejsv.f"
			    u[p + q * u_dim1] = work[(*n << 1) + *n * nr + nr 
				    + p];
#line 1473 "dgejsv.f"
/* L874: */
#line 1473 "dgejsv.f"
			}
#line 1474 "dgejsv.f"
/* L873: */
#line 1474 "dgejsv.f"
		    }
#line 1475 "dgejsv.f"
		    if (nr < *n) {
#line 1476 "dgejsv.f"
			i__1 = *n - nr;
#line 1476 "dgejsv.f"
			dlaset_("A", &i__1, &nr, &c_b34, &c_b34, &v[nr + 1 + 
				v_dim1], ldv, (ftnlen)1);
#line 1477 "dgejsv.f"
			i__1 = *n - nr;
#line 1477 "dgejsv.f"
			dlaset_("A", &nr, &i__1, &c_b34, &c_b34, &v[(nr + 1) *
				 v_dim1 + 1], ldv, (ftnlen)1);
#line 1478 "dgejsv.f"
			i__1 = *n - nr;
#line 1478 "dgejsv.f"
			i__2 = *n - nr;
#line 1478 "dgejsv.f"
			dlaset_("A", &i__1, &i__2, &c_b34, &c_b35, &v[nr + 1 
				+ (nr + 1) * v_dim1], ldv, (ftnlen)1);
#line 1479 "dgejsv.f"
		    }
#line 1480 "dgejsv.f"
		    i__1 = *lwork - (*n << 1) - *n * nr - nr;
#line 1480 "dgejsv.f"
		    dormqr_("L", "N", n, n, &nr, &work[(*n << 1) + 1], n, &
			    work[*n + 1], &v[v_offset], ldv, &work[(*n << 1) 
			    + *n * nr + nr + 1], &i__1, &ierr, (ftnlen)1, (
			    ftnlen)1);
#line 1482 "dgejsv.f"
		} else {
/*              Last line of defense. */
/* #:(          This is a rather pathological case: no scaled condition */
/*              improvement after two pivoted QR factorizations. Other */
/*              possibility is that the rank revealing QR factorization */
/*              or the condition estimator has failed, or the COND_OK */
/*              is set very close to ONE (which is unnecessary). Normally, */
/*              this branch should never be executed, but in rare cases of */
/*              failure of the RRQR or condition estimator, the last line of */
/*              defense ensures that DGEJSV completes the task. */
/*              Compute the full SVD of L3 using DGESVJ with explicit */
/*              accumulation of Jacobi rotations. */
#line 1494 "dgejsv.f"
		    i__1 = *lwork - (*n << 1) - *n * nr - nr;
#line 1494 "dgejsv.f"
		    dgesvj_("L", "U", "V", &nr, &nr, &v[v_offset], ldv, &sva[
			    1], &nr, &u[u_offset], ldu, &work[(*n << 1) + *n *
			     nr + nr + 1], &i__1, info, (ftnlen)1, (ftnlen)1, 
			    (ftnlen)1);
#line 1496 "dgejsv.f"
		    scalem = work[(*n << 1) + *n * nr + nr + 1];
#line 1497 "dgejsv.f"
		    numrank = i_dnnt(&work[(*n << 1) + *n * nr + nr + 2]);
#line 1498 "dgejsv.f"
		    if (nr < *n) {
#line 1499 "dgejsv.f"
			i__1 = *n - nr;
#line 1499 "dgejsv.f"
			dlaset_("A", &i__1, &nr, &c_b34, &c_b34, &v[nr + 1 + 
				v_dim1], ldv, (ftnlen)1);
#line 1500 "dgejsv.f"
			i__1 = *n - nr;
#line 1500 "dgejsv.f"
			dlaset_("A", &nr, &i__1, &c_b34, &c_b34, &v[(nr + 1) *
				 v_dim1 + 1], ldv, (ftnlen)1);
#line 1501 "dgejsv.f"
			i__1 = *n - nr;
#line 1501 "dgejsv.f"
			i__2 = *n - nr;
#line 1501 "dgejsv.f"
			dlaset_("A", &i__1, &i__2, &c_b34, &c_b35, &v[nr + 1 
				+ (nr + 1) * v_dim1], ldv, (ftnlen)1);
#line 1502 "dgejsv.f"
		    }
#line 1503 "dgejsv.f"
		    i__1 = *lwork - (*n << 1) - *n * nr - nr;
#line 1503 "dgejsv.f"
		    dormqr_("L", "N", n, n, &nr, &work[(*n << 1) + 1], n, &
			    work[*n + 1], &v[v_offset], ldv, &work[(*n << 1) 
			    + *n * nr + nr + 1], &i__1, &ierr, (ftnlen)1, (
			    ftnlen)1);

#line 1506 "dgejsv.f"
		    i__1 = *lwork - (*n << 1) - *n * nr - nr;
#line 1506 "dgejsv.f"
		    dormlq_("L", "T", &nr, &nr, &nr, &work[(*n << 1) + 1], n, 
			    &work[(*n << 1) + *n * nr + 1], &u[u_offset], ldu,
			     &work[(*n << 1) + *n * nr + nr + 1], &i__1, &
			    ierr, (ftnlen)1, (ftnlen)1);
#line 1509 "dgejsv.f"
		    i__1 = nr;
#line 1509 "dgejsv.f"
		    for (q = 1; q <= i__1; ++q) {
#line 1510 "dgejsv.f"
			i__2 = nr;
#line 1510 "dgejsv.f"
			for (p = 1; p <= i__2; ++p) {
#line 1511 "dgejsv.f"
			    work[(*n << 1) + *n * nr + nr + iwork[*n + p]] = 
				    u[p + q * u_dim1];
#line 1512 "dgejsv.f"
/* L772: */
#line 1512 "dgejsv.f"
			}
#line 1513 "dgejsv.f"
			i__2 = nr;
#line 1513 "dgejsv.f"
			for (p = 1; p <= i__2; ++p) {
#line 1514 "dgejsv.f"
			    u[p + q * u_dim1] = work[(*n << 1) + *n * nr + nr 
				    + p];
#line 1515 "dgejsv.f"
/* L774: */
#line 1515 "dgejsv.f"
			}
#line 1516 "dgejsv.f"
/* L773: */
#line 1516 "dgejsv.f"
		    }

#line 1518 "dgejsv.f"
		}

/*           Permute the rows of V using the (column) permutation from the */
/*           first QRF. Also, scale the columns to make them unit in */
/*           Euclidean norm. This applies to all cases. */

#line 1524 "dgejsv.f"
		temp1 = sqrt((doublereal) (*n)) * epsln;
#line 1525 "dgejsv.f"
		i__1 = *n;
#line 1525 "dgejsv.f"
		for (q = 1; q <= i__1; ++q) {
#line 1526 "dgejsv.f"
		    i__2 = *n;
#line 1526 "dgejsv.f"
		    for (p = 1; p <= i__2; ++p) {
#line 1527 "dgejsv.f"
			work[(*n << 1) + *n * nr + nr + iwork[p]] = v[p + q * 
				v_dim1];
#line 1528 "dgejsv.f"
/* L972: */
#line 1528 "dgejsv.f"
		    }
#line 1529 "dgejsv.f"
		    i__2 = *n;
#line 1529 "dgejsv.f"
		    for (p = 1; p <= i__2; ++p) {
#line 1530 "dgejsv.f"
			v[p + q * v_dim1] = work[(*n << 1) + *n * nr + nr + p]
				;
#line 1531 "dgejsv.f"
/* L973: */
#line 1531 "dgejsv.f"
		    }
#line 1532 "dgejsv.f"
		    xsc = 1. / dnrm2_(n, &v[q * v_dim1 + 1], &c__1);
#line 1533 "dgejsv.f"
		    if (xsc < 1. - temp1 || xsc > temp1 + 1.) {
#line 1533 "dgejsv.f"
			dscal_(n, &xsc, &v[q * v_dim1 + 1], &c__1);
#line 1533 "dgejsv.f"
		    }
#line 1535 "dgejsv.f"
/* L1972: */
#line 1535 "dgejsv.f"
		}
/*           At this moment, V contains the right singular vectors of A. */
/*           Next, assemble the left singular vector matrix U (M x N). */
#line 1538 "dgejsv.f"
		if (nr < *m) {
#line 1539 "dgejsv.f"
		    i__1 = *m - nr;
#line 1539 "dgejsv.f"
		    dlaset_("A", &i__1, &nr, &c_b34, &c_b34, &u[nr + 1 + 
			    u_dim1], ldu, (ftnlen)1);
#line 1540 "dgejsv.f"
		    if (nr < n1) {
#line 1541 "dgejsv.f"
			i__1 = n1 - nr;
#line 1541 "dgejsv.f"
			dlaset_("A", &nr, &i__1, &c_b34, &c_b34, &u[(nr + 1) *
				 u_dim1 + 1], ldu, (ftnlen)1);
#line 1542 "dgejsv.f"
			i__1 = *m - nr;
#line 1542 "dgejsv.f"
			i__2 = n1 - nr;
#line 1542 "dgejsv.f"
			dlaset_("A", &i__1, &i__2, &c_b34, &c_b35, &u[nr + 1 
				+ (nr + 1) * u_dim1], ldu, (ftnlen)1);
#line 1543 "dgejsv.f"
		    }
#line 1544 "dgejsv.f"
		}

/*           The Q matrix from the first QRF is built into the left singular */
/*           matrix U. This applies to all cases. */

#line 1549 "dgejsv.f"
		i__1 = *lwork - *n;
#line 1549 "dgejsv.f"
		dormqr_("Left", "No_Tr", m, &n1, n, &a[a_offset], lda, &work[
			1], &u[u_offset], ldu, &work[*n + 1], &i__1, &ierr, (
			ftnlen)4, (ftnlen)5);
/*           The columns of U are normalized. The cost is O(M*N) flops. */
#line 1553 "dgejsv.f"
		temp1 = sqrt((doublereal) (*m)) * epsln;
#line 1554 "dgejsv.f"
		i__1 = nr;
#line 1554 "dgejsv.f"
		for (p = 1; p <= i__1; ++p) {
#line 1555 "dgejsv.f"
		    xsc = 1. / dnrm2_(m, &u[p * u_dim1 + 1], &c__1);
#line 1556 "dgejsv.f"
		    if (xsc < 1. - temp1 || xsc > temp1 + 1.) {
#line 1556 "dgejsv.f"
			dscal_(m, &xsc, &u[p * u_dim1 + 1], &c__1);
#line 1556 "dgejsv.f"
		    }
#line 1558 "dgejsv.f"
/* L1973: */
#line 1558 "dgejsv.f"
		}

/*           If the initial QRF is computed with row pivoting, the left */
/*           singular vectors must be adjusted. */

#line 1563 "dgejsv.f"
		if (rowpiv) {
#line 1563 "dgejsv.f"
		    i__1 = *m - 1;
#line 1563 "dgejsv.f"
		    dlaswp_(&n1, &u[u_offset], ldu, &c__1, &i__1, &iwork[(*n 
			    << 1) + 1], &c_n1);
#line 1563 "dgejsv.f"
		}

#line 1566 "dgejsv.f"
	    } else {

/*        .. the initial matrix A has almost orthogonal columns and */
/*        the second QRF is not needed */

#line 1571 "dgejsv.f"
		dlacpy_("Upper", n, n, &a[a_offset], lda, &work[*n + 1], n, (
			ftnlen)5);
#line 1572 "dgejsv.f"
		if (l2pert) {
#line 1573 "dgejsv.f"
		    xsc = sqrt(small);
#line 1574 "dgejsv.f"
		    i__1 = *n;
#line 1574 "dgejsv.f"
		    for (p = 2; p <= i__1; ++p) {
#line 1575 "dgejsv.f"
			temp1 = xsc * work[*n + (p - 1) * *n + p];
#line 1576 "dgejsv.f"
			i__2 = p - 1;
#line 1576 "dgejsv.f"
			for (q = 1; q <= i__2; ++q) {
#line 1577 "dgejsv.f"
			    work[*n + (q - 1) * *n + p] = -d_sign(&temp1, &
				    work[*n + (p - 1) * *n + q]);
#line 1578 "dgejsv.f"
/* L5971: */
#line 1578 "dgejsv.f"
			}
#line 1579 "dgejsv.f"
/* L5970: */
#line 1579 "dgejsv.f"
		    }
#line 1580 "dgejsv.f"
		} else {
#line 1581 "dgejsv.f"
		    i__1 = *n - 1;
#line 1581 "dgejsv.f"
		    i__2 = *n - 1;
#line 1581 "dgejsv.f"
		    dlaset_("Lower", &i__1, &i__2, &c_b34, &c_b34, &work[*n + 
			    2], n, (ftnlen)5);
#line 1582 "dgejsv.f"
		}

#line 1584 "dgejsv.f"
		i__1 = *lwork - *n - *n * *n;
#line 1584 "dgejsv.f"
		dgesvj_("Upper", "U", "N", n, n, &work[*n + 1], n, &sva[1], n,
			 &u[u_offset], ldu, &work[*n + *n * *n + 1], &i__1, 
			info, (ftnlen)5, (ftnlen)1, (ftnlen)1);

#line 1587 "dgejsv.f"
		scalem = work[*n + *n * *n + 1];
#line 1588 "dgejsv.f"
		numrank = i_dnnt(&work[*n + *n * *n + 2]);
#line 1589 "dgejsv.f"
		i__1 = *n;
#line 1589 "dgejsv.f"
		for (p = 1; p <= i__1; ++p) {
#line 1590 "dgejsv.f"
		    dcopy_(n, &work[*n + (p - 1) * *n + 1], &c__1, &u[p * 
			    u_dim1 + 1], &c__1);
#line 1591 "dgejsv.f"
		    dscal_(n, &sva[p], &work[*n + (p - 1) * *n + 1], &c__1);
#line 1592 "dgejsv.f"
/* L6970: */
#line 1592 "dgejsv.f"
		}

#line 1594 "dgejsv.f"
		dtrsm_("Left", "Upper", "NoTrans", "No UD", n, n, &c_b35, &a[
			a_offset], lda, &work[*n + 1], n, (ftnlen)4, (ftnlen)
			5, (ftnlen)7, (ftnlen)5);
#line 1596 "dgejsv.f"
		i__1 = *n;
#line 1596 "dgejsv.f"
		for (p = 1; p <= i__1; ++p) {
#line 1597 "dgejsv.f"
		    dcopy_(n, &work[*n + p], n, &v[iwork[p] + v_dim1], ldv);
#line 1598 "dgejsv.f"
/* L6972: */
#line 1598 "dgejsv.f"
		}
#line 1599 "dgejsv.f"
		temp1 = sqrt((doublereal) (*n)) * epsln;
#line 1600 "dgejsv.f"
		i__1 = *n;
#line 1600 "dgejsv.f"
		for (p = 1; p <= i__1; ++p) {
#line 1601 "dgejsv.f"
		    xsc = 1. / dnrm2_(n, &v[p * v_dim1 + 1], &c__1);
#line 1602 "dgejsv.f"
		    if (xsc < 1. - temp1 || xsc > temp1 + 1.) {
#line 1602 "dgejsv.f"
			dscal_(n, &xsc, &v[p * v_dim1 + 1], &c__1);
#line 1602 "dgejsv.f"
		    }
#line 1604 "dgejsv.f"
/* L6971: */
#line 1604 "dgejsv.f"
		}

/*           Assemble the left singular vector matrix U (M x N). */

#line 1608 "dgejsv.f"
		if (*n < *m) {
#line 1609 "dgejsv.f"
		    i__1 = *m - *n;
#line 1609 "dgejsv.f"
		    dlaset_("A", &i__1, n, &c_b34, &c_b34, &u[*n + 1 + u_dim1]
			    , ldu, (ftnlen)1);
#line 1610 "dgejsv.f"
		    if (*n < n1) {
#line 1611 "dgejsv.f"
			i__1 = n1 - *n;
#line 1611 "dgejsv.f"
			dlaset_("A", n, &i__1, &c_b34, &c_b34, &u[(*n + 1) * 
				u_dim1 + 1], ldu, (ftnlen)1);
#line 1612 "dgejsv.f"
			i__1 = *m - *n;
#line 1612 "dgejsv.f"
			i__2 = n1 - *n;
#line 1612 "dgejsv.f"
			dlaset_("A", &i__1, &i__2, &c_b34, &c_b35, &u[*n + 1 
				+ (*n + 1) * u_dim1], ldu, (ftnlen)1);
#line 1613 "dgejsv.f"
		    }
#line 1614 "dgejsv.f"
		}
#line 1615 "dgejsv.f"
		i__1 = *lwork - *n;
#line 1615 "dgejsv.f"
		dormqr_("Left", "No Tr", m, &n1, n, &a[a_offset], lda, &work[
			1], &u[u_offset], ldu, &work[*n + 1], &i__1, &ierr, (
			ftnlen)4, (ftnlen)5);
#line 1617 "dgejsv.f"
		temp1 = sqrt((doublereal) (*m)) * epsln;
#line 1618 "dgejsv.f"
		i__1 = n1;
#line 1618 "dgejsv.f"
		for (p = 1; p <= i__1; ++p) {
#line 1619 "dgejsv.f"
		    xsc = 1. / dnrm2_(m, &u[p * u_dim1 + 1], &c__1);
#line 1620 "dgejsv.f"
		    if (xsc < 1. - temp1 || xsc > temp1 + 1.) {
#line 1620 "dgejsv.f"
			dscal_(m, &xsc, &u[p * u_dim1 + 1], &c__1);
#line 1620 "dgejsv.f"
		    }
#line 1622 "dgejsv.f"
/* L6973: */
#line 1622 "dgejsv.f"
		}

#line 1624 "dgejsv.f"
		if (rowpiv) {
#line 1624 "dgejsv.f"
		    i__1 = *m - 1;
#line 1624 "dgejsv.f"
		    dlaswp_(&n1, &u[u_offset], ldu, &c__1, &i__1, &iwork[(*n 
			    << 1) + 1], &c_n1);
#line 1624 "dgejsv.f"
		}

#line 1627 "dgejsv.f"
	    }

/*        end of the  >> almost orthogonal case <<  in the full SVD */

#line 1631 "dgejsv.f"
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

#line 1643 "dgejsv.f"
	    i__1 = nr;
#line 1643 "dgejsv.f"
	    for (p = 1; p <= i__1; ++p) {
#line 1644 "dgejsv.f"
		i__2 = *n - p + 1;
#line 1644 "dgejsv.f"
		dcopy_(&i__2, &a[p + p * a_dim1], lda, &v[p + p * v_dim1], &
			c__1);
#line 1645 "dgejsv.f"
/* L7968: */
#line 1645 "dgejsv.f"
	    }

#line 1647 "dgejsv.f"
	    if (l2pert) {
#line 1648 "dgejsv.f"
		xsc = sqrt(small / epsln);
#line 1649 "dgejsv.f"
		i__1 = nr;
#line 1649 "dgejsv.f"
		for (q = 1; q <= i__1; ++q) {
#line 1650 "dgejsv.f"
		    temp1 = xsc * (d__1 = v[q + q * v_dim1], abs(d__1));
#line 1651 "dgejsv.f"
		    i__2 = *n;
#line 1651 "dgejsv.f"
		    for (p = 1; p <= i__2; ++p) {
#line 1652 "dgejsv.f"
			if (p > q && (d__1 = v[p + q * v_dim1], abs(d__1)) <= 
				temp1 || p < q) {
#line 1652 "dgejsv.f"
			    v[p + q * v_dim1] = d_sign(&temp1, &v[p + q * 
				    v_dim1]);
#line 1652 "dgejsv.f"
			}
#line 1655 "dgejsv.f"
			if (p < q) {
#line 1655 "dgejsv.f"
			    v[p + q * v_dim1] = -v[p + q * v_dim1];
#line 1655 "dgejsv.f"
			}
#line 1656 "dgejsv.f"
/* L5968: */
#line 1656 "dgejsv.f"
		    }
#line 1657 "dgejsv.f"
/* L5969: */
#line 1657 "dgejsv.f"
		}
#line 1658 "dgejsv.f"
	    } else {
#line 1659 "dgejsv.f"
		i__1 = nr - 1;
#line 1659 "dgejsv.f"
		i__2 = nr - 1;
#line 1659 "dgejsv.f"
		dlaset_("U", &i__1, &i__2, &c_b34, &c_b34, &v[(v_dim1 << 1) + 
			1], ldv, (ftnlen)1);
#line 1660 "dgejsv.f"
	    }
#line 1662 "dgejsv.f"
	    i__1 = *lwork - (*n << 1);
#line 1662 "dgejsv.f"
	    dgeqrf_(n, &nr, &v[v_offset], ldv, &work[*n + 1], &work[(*n << 1) 
		    + 1], &i__1, &ierr);
#line 1664 "dgejsv.f"
	    dlacpy_("L", n, &nr, &v[v_offset], ldv, &work[(*n << 1) + 1], n, (
		    ftnlen)1);

#line 1666 "dgejsv.f"
	    i__1 = nr;
#line 1666 "dgejsv.f"
	    for (p = 1; p <= i__1; ++p) {
#line 1667 "dgejsv.f"
		i__2 = nr - p + 1;
#line 1667 "dgejsv.f"
		dcopy_(&i__2, &v[p + p * v_dim1], ldv, &u[p + p * u_dim1], &
			c__1);
#line 1668 "dgejsv.f"
/* L7969: */
#line 1668 "dgejsv.f"
	    }
#line 1670 "dgejsv.f"
	    if (l2pert) {
#line 1671 "dgejsv.f"
		xsc = sqrt(small / epsln);
#line 1672 "dgejsv.f"
		i__1 = nr;
#line 1672 "dgejsv.f"
		for (q = 2; q <= i__1; ++q) {
#line 1673 "dgejsv.f"
		    i__2 = q - 1;
#line 1673 "dgejsv.f"
		    for (p = 1; p <= i__2; ++p) {
/* Computing MIN */
#line 1674 "dgejsv.f"
			d__3 = (d__1 = u[p + p * u_dim1], abs(d__1)), d__4 = (
				d__2 = u[q + q * u_dim1], abs(d__2));
#line 1674 "dgejsv.f"
			temp1 = xsc * min(d__3,d__4);
#line 1675 "dgejsv.f"
			u[p + q * u_dim1] = -d_sign(&temp1, &u[q + p * u_dim1]
				);
#line 1676 "dgejsv.f"
/* L9971: */
#line 1676 "dgejsv.f"
		    }
#line 1677 "dgejsv.f"
/* L9970: */
#line 1677 "dgejsv.f"
		}
#line 1678 "dgejsv.f"
	    } else {
#line 1679 "dgejsv.f"
		i__1 = nr - 1;
#line 1679 "dgejsv.f"
		i__2 = nr - 1;
#line 1679 "dgejsv.f"
		dlaset_("U", &i__1, &i__2, &c_b34, &c_b34, &u[(u_dim1 << 1) + 
			1], ldu, (ftnlen)1);
#line 1680 "dgejsv.f"
	    }
#line 1682 "dgejsv.f"
	    i__1 = *lwork - (*n << 1) - *n * nr;
#line 1682 "dgejsv.f"
	    dgesvj_("G", "U", "V", &nr, &nr, &u[u_offset], ldu, &sva[1], n, &
		    v[v_offset], ldv, &work[(*n << 1) + *n * nr + 1], &i__1, 
		    info, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1684 "dgejsv.f"
	    scalem = work[(*n << 1) + *n * nr + 1];
#line 1685 "dgejsv.f"
	    numrank = i_dnnt(&work[(*n << 1) + *n * nr + 2]);
#line 1687 "dgejsv.f"
	    if (nr < *n) {
#line 1688 "dgejsv.f"
		i__1 = *n - nr;
#line 1688 "dgejsv.f"
		dlaset_("A", &i__1, &nr, &c_b34, &c_b34, &v[nr + 1 + v_dim1], 
			ldv, (ftnlen)1);
#line 1689 "dgejsv.f"
		i__1 = *n - nr;
#line 1689 "dgejsv.f"
		dlaset_("A", &nr, &i__1, &c_b34, &c_b34, &v[(nr + 1) * v_dim1 
			+ 1], ldv, (ftnlen)1);
#line 1690 "dgejsv.f"
		i__1 = *n - nr;
#line 1690 "dgejsv.f"
		i__2 = *n - nr;
#line 1690 "dgejsv.f"
		dlaset_("A", &i__1, &i__2, &c_b34, &c_b35, &v[nr + 1 + (nr + 
			1) * v_dim1], ldv, (ftnlen)1);
#line 1691 "dgejsv.f"
	    }
#line 1693 "dgejsv.f"
	    i__1 = *lwork - (*n << 1) - *n * nr - nr;
#line 1693 "dgejsv.f"
	    dormqr_("L", "N", n, n, &nr, &work[(*n << 1) + 1], n, &work[*n + 
		    1], &v[v_offset], ldv, &work[(*n << 1) + *n * nr + nr + 1]
		    , &i__1, &ierr, (ftnlen)1, (ftnlen)1);

/*           Permute the rows of V using the (column) permutation from the */
/*           first QRF. Also, scale the columns to make them unit in */
/*           Euclidean norm. This applies to all cases. */

#line 1700 "dgejsv.f"
	    temp1 = sqrt((doublereal) (*n)) * epsln;
#line 1701 "dgejsv.f"
	    i__1 = *n;
#line 1701 "dgejsv.f"
	    for (q = 1; q <= i__1; ++q) {
#line 1702 "dgejsv.f"
		i__2 = *n;
#line 1702 "dgejsv.f"
		for (p = 1; p <= i__2; ++p) {
#line 1703 "dgejsv.f"
		    work[(*n << 1) + *n * nr + nr + iwork[p]] = v[p + q * 
			    v_dim1];
#line 1704 "dgejsv.f"
/* L8972: */
#line 1704 "dgejsv.f"
		}
#line 1705 "dgejsv.f"
		i__2 = *n;
#line 1705 "dgejsv.f"
		for (p = 1; p <= i__2; ++p) {
#line 1706 "dgejsv.f"
		    v[p + q * v_dim1] = work[(*n << 1) + *n * nr + nr + p];
#line 1707 "dgejsv.f"
/* L8973: */
#line 1707 "dgejsv.f"
		}
#line 1708 "dgejsv.f"
		xsc = 1. / dnrm2_(n, &v[q * v_dim1 + 1], &c__1);
#line 1709 "dgejsv.f"
		if (xsc < 1. - temp1 || xsc > temp1 + 1.) {
#line 1709 "dgejsv.f"
		    dscal_(n, &xsc, &v[q * v_dim1 + 1], &c__1);
#line 1709 "dgejsv.f"
		}
#line 1711 "dgejsv.f"
/* L7972: */
#line 1711 "dgejsv.f"
	    }

/*           At this moment, V contains the right singular vectors of A. */
/*           Next, assemble the left singular vector matrix U (M x N). */

#line 1716 "dgejsv.f"
	    if (nr < *m) {
#line 1717 "dgejsv.f"
		i__1 = *m - nr;
#line 1717 "dgejsv.f"
		dlaset_("A", &i__1, &nr, &c_b34, &c_b34, &u[nr + 1 + u_dim1], 
			ldu, (ftnlen)1);
#line 1718 "dgejsv.f"
		if (nr < n1) {
#line 1719 "dgejsv.f"
		    i__1 = n1 - nr;
#line 1719 "dgejsv.f"
		    dlaset_("A", &nr, &i__1, &c_b34, &c_b34, &u[(nr + 1) * 
			    u_dim1 + 1], ldu, (ftnlen)1);
#line 1720 "dgejsv.f"
		    i__1 = *m - nr;
#line 1720 "dgejsv.f"
		    i__2 = n1 - nr;
#line 1720 "dgejsv.f"
		    dlaset_("A", &i__1, &i__2, &c_b34, &c_b35, &u[nr + 1 + (
			    nr + 1) * u_dim1], ldu, (ftnlen)1);
#line 1721 "dgejsv.f"
		}
#line 1722 "dgejsv.f"
	    }

#line 1724 "dgejsv.f"
	    i__1 = *lwork - *n;
#line 1724 "dgejsv.f"
	    dormqr_("Left", "No Tr", m, &n1, n, &a[a_offset], lda, &work[1], &
		    u[u_offset], ldu, &work[*n + 1], &i__1, &ierr, (ftnlen)4, 
		    (ftnlen)5);

#line 1727 "dgejsv.f"
	    if (rowpiv) {
#line 1727 "dgejsv.f"
		i__1 = *m - 1;
#line 1727 "dgejsv.f"
		dlaswp_(&n1, &u[u_offset], ldu, &c__1, &i__1, &iwork[(*n << 1)
			 + 1], &c_n1);
#line 1727 "dgejsv.f"
	    }


#line 1731 "dgejsv.f"
	}
#line 1732 "dgejsv.f"
	if (transp) {
/*           .. swap U and V because the procedure worked on A^t */
#line 1734 "dgejsv.f"
	    i__1 = *n;
#line 1734 "dgejsv.f"
	    for (p = 1; p <= i__1; ++p) {
#line 1735 "dgejsv.f"
		dswap_(n, &u[p * u_dim1 + 1], &c__1, &v[p * v_dim1 + 1], &
			c__1);
#line 1736 "dgejsv.f"
/* L6974: */
#line 1736 "dgejsv.f"
	    }
#line 1737 "dgejsv.f"
	}

#line 1739 "dgejsv.f"
    }
/*     end of the full SVD */

/*     Undo scaling, if necessary (and possible) */

#line 1744 "dgejsv.f"
    if (uscal2 <= big / sva[1] * uscal1) {
#line 1745 "dgejsv.f"
	dlascl_("G", &c__0, &c__0, &uscal1, &uscal2, &nr, &c__1, &sva[1], n, &
		ierr, (ftnlen)1);
#line 1746 "dgejsv.f"
	uscal1 = 1.;
#line 1747 "dgejsv.f"
	uscal2 = 1.;
#line 1748 "dgejsv.f"
    }

#line 1750 "dgejsv.f"
    if (nr < *n) {
#line 1751 "dgejsv.f"
	i__1 = *n;
#line 1751 "dgejsv.f"
	for (p = nr + 1; p <= i__1; ++p) {
#line 1752 "dgejsv.f"
	    sva[p] = 0.;
#line 1753 "dgejsv.f"
/* L3004: */
#line 1753 "dgejsv.f"
	}
#line 1754 "dgejsv.f"
    }

#line 1756 "dgejsv.f"
    work[1] = uscal2 * scalem;
#line 1757 "dgejsv.f"
    work[2] = uscal1;
#line 1758 "dgejsv.f"
    if (errest) {
#line 1758 "dgejsv.f"
	work[3] = sconda;
#line 1758 "dgejsv.f"
    }
#line 1759 "dgejsv.f"
    if (lsvec && rsvec) {
#line 1760 "dgejsv.f"
	work[4] = condr1;
#line 1761 "dgejsv.f"
	work[5] = condr2;
#line 1762 "dgejsv.f"
    }
#line 1763 "dgejsv.f"
    if (l2tran) {
#line 1764 "dgejsv.f"
	work[6] = entra;
#line 1765 "dgejsv.f"
	work[7] = entrat;
#line 1766 "dgejsv.f"
    }

#line 1768 "dgejsv.f"
    iwork[1] = nr;
#line 1769 "dgejsv.f"
    iwork[2] = numrank;
#line 1770 "dgejsv.f"
    iwork[3] = warning;

#line 1772 "dgejsv.f"
    return 0;
/*     .. */
/*     .. END OF DGEJSV */
/*     .. */
} /* dgejsv_ */

