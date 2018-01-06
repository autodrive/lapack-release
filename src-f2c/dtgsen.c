#line 1 "dtgsen.f"
/* dtgsen.f -- translated by f2c (version 20100827).
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

#line 1 "dtgsen.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c__2 = 2;
static doublereal c_b28 = 1.;

/* > \brief \b DTGSEN */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DTGSEN + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dtgsen.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dtgsen.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dtgsen.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DTGSEN( IJOB, WANTQ, WANTZ, SELECT, N, A, LDA, B, LDB, */
/*                          ALPHAR, ALPHAI, BETA, Q, LDQ, Z, LDZ, M, PL, */
/*                          PR, DIF, WORK, LWORK, IWORK, LIWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       LOGICAL            WANTQ, WANTZ */
/*       INTEGER            IJOB, INFO, LDA, LDB, LDQ, LDZ, LIWORK, LWORK, */
/*      $                   M, N */
/*       DOUBLE PRECISION   PL, PR */
/*       .. */
/*       .. Array Arguments .. */
/*       LOGICAL            SELECT( * ) */
/*       INTEGER            IWORK( * ) */
/*       DOUBLE PRECISION   A( LDA, * ), ALPHAI( * ), ALPHAR( * ), */
/*      $                   B( LDB, * ), BETA( * ), DIF( * ), Q( LDQ, * ), */
/*      $                   WORK( * ), Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DTGSEN reorders the generalized real Schur decomposition of a real */
/* > matrix pair (A, B) (in terms of an orthonormal equivalence trans- */
/* > formation Q**T * (A, B) * Z), so that a selected cluster of eigenvalues */
/* > appears in the leading diagonal blocks of the upper quasi-triangular */
/* > matrix A and the upper triangular B. The leading columns of Q and */
/* > Z form orthonormal bases of the corresponding left and right eigen- */
/* > spaces (deflating subspaces). (A, B) must be in generalized real */
/* > Schur canonical form (as returned by DGGES), i.e. A is block upper */
/* > triangular with 1-by-1 and 2-by-2 diagonal blocks. B is upper */
/* > triangular. */
/* > */
/* > DTGSEN also computes the generalized eigenvalues */
/* > */
/* >             w(j) = (ALPHAR(j) + i*ALPHAI(j))/BETA(j) */
/* > */
/* > of the reordered matrix pair (A, B). */
/* > */
/* > Optionally, DTGSEN computes the estimates of reciprocal condition */
/* > numbers for eigenvalues and eigenspaces. These are Difu[(A11,B11), */
/* > (A22,B22)] and Difl[(A11,B11), (A22,B22)], i.e. the separation(s) */
/* > between the matrix pairs (A11, B11) and (A22,B22) that correspond to */
/* > the selected cluster and the eigenvalues outside the cluster, resp., */
/* > and norms of "projections" onto left and right eigenspaces w.r.t. */
/* > the selected cluster in the (1,1)-block. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] IJOB */
/* > \verbatim */
/* >          IJOB is INTEGER */
/* >          Specifies whether condition numbers are required for the */
/* >          cluster of eigenvalues (PL and PR) or the deflating subspaces */
/* >          (Difu and Difl): */
/* >           =0: Only reorder w.r.t. SELECT. No extras. */
/* >           =1: Reciprocal of norms of "projections" onto left and right */
/* >               eigenspaces w.r.t. the selected cluster (PL and PR). */
/* >           =2: Upper bounds on Difu and Difl. F-norm-based estimate */
/* >               (DIF(1:2)). */
/* >           =3: Estimate of Difu and Difl. 1-norm-based estimate */
/* >               (DIF(1:2)). */
/* >               About 5 times as expensive as IJOB = 2. */
/* >           =4: Compute PL, PR and DIF (i.e. 0, 1 and 2 above): Economic */
/* >               version to get it all. */
/* >           =5: Compute PL, PR and DIF (i.e. 0, 1 and 3 above) */
/* > \endverbatim */
/* > */
/* > \param[in] WANTQ */
/* > \verbatim */
/* >          WANTQ is LOGICAL */
/* >          .TRUE. : update the left transformation matrix Q; */
/* >          .FALSE.: do not update Q. */
/* > \endverbatim */
/* > */
/* > \param[in] WANTZ */
/* > \verbatim */
/* >          WANTZ is LOGICAL */
/* >          .TRUE. : update the right transformation matrix Z; */
/* >          .FALSE.: do not update Z. */
/* > \endverbatim */
/* > */
/* > \param[in] SELECT */
/* > \verbatim */
/* >          SELECT is LOGICAL array, dimension (N) */
/* >          SELECT specifies the eigenvalues in the selected cluster. */
/* >          To select a real eigenvalue w(j), SELECT(j) must be set to */
/* >          .TRUE.. To select a complex conjugate pair of eigenvalues */
/* >          w(j) and w(j+1), corresponding to a 2-by-2 diagonal block, */
/* >          either SELECT(j) or SELECT(j+1) or both must be set to */
/* >          .TRUE.; a complex conjugate pair of eigenvalues must be */
/* >          either both included in the cluster or both excluded. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrices A and B. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension(LDA,N) */
/* >          On entry, the upper quasi-triangular matrix A, with (A, B) in */
/* >          generalized real Schur canonical form. */
/* >          On exit, A is overwritten by the reordered matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A. LDA >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* >          B is DOUBLE PRECISION array, dimension(LDB,N) */
/* >          On entry, the upper triangular matrix B, with (A, B) in */
/* >          generalized real Schur canonical form. */
/* >          On exit, B is overwritten by the reordered matrix B. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >          The leading dimension of the array B. LDB >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] ALPHAR */
/* > \verbatim */
/* >          ALPHAR is DOUBLE PRECISION array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] ALPHAI */
/* > \verbatim */
/* >          ALPHAI is DOUBLE PRECISION array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] BETA */
/* > \verbatim */
/* >          BETA is DOUBLE PRECISION array, dimension (N) */
/* > */
/* >          On exit, (ALPHAR(j) + ALPHAI(j)*i)/BETA(j), j=1,...,N, will */
/* >          be the generalized eigenvalues.  ALPHAR(j) + ALPHAI(j)*i */
/* >          and BETA(j),j=1,...,N  are the diagonals of the complex Schur */
/* >          form (S,T) that would result if the 2-by-2 diagonal blocks of */
/* >          the real generalized Schur form of (A,B) were further reduced */
/* >          to triangular form using complex unitary transformations. */
/* >          If ALPHAI(j) is zero, then the j-th eigenvalue is real; if */
/* >          positive, then the j-th and (j+1)-st eigenvalues are a */
/* >          complex conjugate pair, with ALPHAI(j+1) negative. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Q */
/* > \verbatim */
/* >          Q is DOUBLE PRECISION array, dimension (LDQ,N) */
/* >          On entry, if WANTQ = .TRUE., Q is an N-by-N matrix. */
/* >          On exit, Q has been postmultiplied by the left orthogonal */
/* >          transformation matrix which reorder (A, B); The leading M */
/* >          columns of Q form orthonormal bases for the specified pair of */
/* >          left eigenspaces (deflating subspaces). */
/* >          If WANTQ = .FALSE., Q is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDQ */
/* > \verbatim */
/* >          LDQ is INTEGER */
/* >          The leading dimension of the array Q.  LDQ >= 1; */
/* >          and if WANTQ = .TRUE., LDQ >= N. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Z */
/* > \verbatim */
/* >          Z is DOUBLE PRECISION array, dimension (LDZ,N) */
/* >          On entry, if WANTZ = .TRUE., Z is an N-by-N matrix. */
/* >          On exit, Z has been postmultiplied by the left orthogonal */
/* >          transformation matrix which reorder (A, B); The leading M */
/* >          columns of Z form orthonormal bases for the specified pair of */
/* >          left eigenspaces (deflating subspaces). */
/* >          If WANTZ = .FALSE., Z is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDZ */
/* > \verbatim */
/* >          LDZ is INTEGER */
/* >          The leading dimension of the array Z. LDZ >= 1; */
/* >          If WANTZ = .TRUE., LDZ >= N. */
/* > \endverbatim */
/* > */
/* > \param[out] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The dimension of the specified pair of left and right eigen- */
/* >          spaces (deflating subspaces). 0 <= M <= N. */
/* > \endverbatim */
/* > */
/* > \param[out] PL */
/* > \verbatim */
/* >          PL is DOUBLE PRECISION */
/* > \endverbatim */
/* > \param[out] PR */
/* > \verbatim */
/* >          PR is DOUBLE PRECISION */
/* > */
/* >          If IJOB = 1, 4 or 5, PL, PR are lower bounds on the */
/* >          reciprocal of the norm of "projections" onto left and right */
/* >          eigenspaces with respect to the selected cluster. */
/* >          0 < PL, PR <= 1. */
/* >          If M = 0 or M = N, PL = PR  = 1. */
/* >          If IJOB = 0, 2 or 3, PL and PR are not referenced. */
/* > \endverbatim */
/* > */
/* > \param[out] DIF */
/* > \verbatim */
/* >          DIF is DOUBLE PRECISION array, dimension (2). */
/* >          If IJOB >= 2, DIF(1:2) store the estimates of Difu and Difl. */
/* >          If IJOB = 2 or 4, DIF(1:2) are F-norm-based upper bounds on */
/* >          Difu and Difl. If IJOB = 3 or 5, DIF(1:2) are 1-norm-based */
/* >          estimates of Difu and Difl. */
/* >          If M = 0 or N, DIF(1:2) = F-norm([A, B]). */
/* >          If IJOB = 0 or 1, DIF is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, */
/* >          dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The dimension of the array WORK. LWORK >=  4*N+16. */
/* >          If IJOB = 1, 2 or 4, LWORK >= MAX(4*N+16, 2*M*(N-M)). */
/* >          If IJOB = 3 or 5, LWORK >= MAX(4*N+16, 4*M*(N-M)). */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the optimal size of the WORK array, returns */
/* >          this value as the first entry of the WORK array, and no error */
/* >          message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* >          IWORK is INTEGER array, dimension (MAX(1,LIWORK)) */
/* >          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LIWORK */
/* > \verbatim */
/* >          LIWORK is INTEGER */
/* >          The dimension of the array IWORK. LIWORK >= 1. */
/* >          If IJOB = 1, 2 or 4, LIWORK >=  N+6. */
/* >          If IJOB = 3 or 5, LIWORK >= MAX(2*M*(N-M), N+6). */
/* > */
/* >          If LIWORK = -1, then a workspace query is assumed; the */
/* >          routine only calculates the optimal size of the IWORK array, */
/* >          returns this value as the first entry of the IWORK array, and */
/* >          no error message related to LIWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >            =0: Successful exit. */
/* >            <0: If INFO = -i, the i-th argument had an illegal value. */
/* >            =1: Reordering of (A, B) failed because the transformed */
/* >                matrix pair (A, B) would be too far from generalized */
/* >                Schur form; the problem is very ill-conditioned. */
/* >                (A, B) may have been partially reordered. */
/* >                If requested, 0 is returned in DIF(*), PL and PR. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup doubleOTHERcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  DTGSEN first collects the selected eigenvalues by computing */
/* >  orthogonal U and W that move them to the top left corner of (A, B). */
/* >  In other words, the selected eigenvalues are the eigenvalues of */
/* >  (A11, B11) in: */
/* > */
/* >              U**T*(A, B)*W = (A11 A12) (B11 B12) n1 */
/* >                              ( 0  A22),( 0  B22) n2 */
/* >                                n1  n2    n1  n2 */
/* > */
/* >  where N = n1+n2 and U**T means the transpose of U. The first n1 columns */
/* >  of U and W span the specified pair of left and right eigenspaces */
/* >  (deflating subspaces) of (A, B). */
/* > */
/* >  If (A, B) has been obtained from the generalized real Schur */
/* >  decomposition of a matrix pair (C, D) = Q*(A, B)*Z**T, then the */
/* >  reordered generalized real Schur form of (C, D) is given by */
/* > */
/* >           (C, D) = (Q*U)*(U**T*(A, B)*W)*(Z*W)**T, */
/* > */
/* >  and the first n1 columns of Q*U and Z*W span the corresponding */
/* >  deflating subspaces of (C, D) (Q and Z store Q*U and Z*W, resp.). */
/* > */
/* >  Note that if the selected eigenvalue is sufficiently ill-conditioned, */
/* >  then its value may differ significantly from its value before */
/* >  reordering. */
/* > */
/* >  The reciprocal condition numbers of the left and right eigenspaces */
/* >  spanned by the first n1 columns of U and W (or Q*U and Z*W) may */
/* >  be returned in DIF(1:2), corresponding to Difu and Difl, resp. */
/* > */
/* >  The Difu and Difl are defined as: */
/* > */
/* >       Difu[(A11, B11), (A22, B22)] = sigma-min( Zu ) */
/* >  and */
/* >       Difl[(A11, B11), (A22, B22)] = Difu[(A22, B22), (A11, B11)], */
/* > */
/* >  where sigma-min(Zu) is the smallest singular value of the */
/* >  (2*n1*n2)-by-(2*n1*n2) matrix */
/* > */
/* >       Zu = [ kron(In2, A11)  -kron(A22**T, In1) ] */
/* >            [ kron(In2, B11)  -kron(B22**T, In1) ]. */
/* > */
/* >  Here, Inx is the identity matrix of size nx and A22**T is the */
/* >  transpose of A22. kron(X, Y) is the Kronecker product between */
/* >  the matrices X and Y. */
/* > */
/* >  When DIF(2) is small, small changes in (A, B) can cause large changes */
/* >  in the deflating subspace. An approximate (asymptotic) bound on the */
/* >  maximum angular error in the computed deflating subspaces is */
/* > */
/* >       EPS * norm((A, B)) / DIF(2), */
/* > */
/* >  where EPS is the machine precision. */
/* > */
/* >  The reciprocal norm of the projectors on the left and right */
/* >  eigenspaces associated with (A11, B11) may be returned in PL and PR. */
/* >  They are computed as follows. First we compute L and R so that */
/* >  P*(A, B)*Q is block diagonal, where */
/* > */
/* >       P = ( I -L ) n1           Q = ( I R ) n1 */
/* >           ( 0  I ) n2    and        ( 0 I ) n2 */
/* >             n1 n2                    n1 n2 */
/* > */
/* >  and (L, R) is the solution to the generalized Sylvester equation */
/* > */
/* >       A11*R - L*A22 = -A12 */
/* >       B11*R - L*B22 = -B12 */
/* > */
/* >  Then PL = (F-norm(L)**2+1)**(-1/2) and PR = (F-norm(R)**2+1)**(-1/2). */
/* >  An approximate (asymptotic) bound on the average absolute error of */
/* >  the selected eigenvalues is */
/* > */
/* >       EPS * norm((A, B)) / PL. */
/* > */
/* >  There are also global error bounds which valid for perturbations up */
/* >  to a certain restriction:  A lower bound (x) on the smallest */
/* >  F-norm(E,F) for which an eigenvalue of (A11, B11) may move and */
/* >  coalesce with an eigenvalue of (A22, B22) under perturbation (E,F), */
/* >  (i.e. (A + E, B + F), is */
/* > */
/* >   x = min(Difu,Difl)/((1/(PL*PL)+1/(PR*PR))**(1/2)+2*max(1/PL,1/PR)). */
/* > */
/* >  An approximate bound on x can be computed from DIF(1:2), PL and PR. */
/* > */
/* >  If y = ( F-norm(E,F) / x) <= 1, the angles between the perturbed */
/* >  (L', R') and unperturbed (L, R) left and right deflating subspaces */
/* >  associated with the selected cluster in the (1,1)-blocks can be */
/* >  bounded as */
/* > */
/* >   max-angle(L, L') <= arctan( y * PL / (1 - y * (1 - PL * PL)**(1/2)) */
/* >   max-angle(R, R') <= arctan( y * PR / (1 - y * (1 - PR * PR)**(1/2)) */
/* > */
/* >  See LAPACK User's Guide section 4.11 or the following references */
/* >  for more information. */
/* > */
/* >  Note that if the default method for computing the Frobenius-norm- */
/* >  based estimate DIF is not wanted (see DLATDF), then the parameter */
/* >  IDIFJB (see below) should be changed from 3 to 4 (routine DLATDF */
/* >  (IJOB = 2 will be used)). See DTGSYL for more details. */
/* > \endverbatim */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Bo Kagstrom and Peter Poromaa, Department of Computing Science, */
/* >     Umea University, S-901 87 Umea, Sweden. */

/* > \par References: */
/*  ================ */
/* > */
/* > \verbatim */
/* > */
/* >  [1] B. Kagstrom; A Direct Method for Reordering Eigenvalues in the */
/* >      Generalized Real Schur Form of a Regular Matrix Pair (A, B), in */
/* >      M.S. Moonen et al (eds), Linear Algebra for Large Scale and */
/* >      Real-Time Applications, Kluwer Academic Publ. 1993, pp 195-218. */
/* > */
/* >  [2] B. Kagstrom and P. Poromaa; Computing Eigenspaces with Specified */
/* >      Eigenvalues of a Regular Matrix Pair (A, B) and Condition */
/* >      Estimation: Theory, Algorithms and Software, */
/* >      Report UMINF - 94.04, Department of Computing Science, Umea */
/* >      University, S-901 87 Umea, Sweden, 1994. Also as LAPACK Working */
/* >      Note 87. To appear in Numerical Algorithms, 1996. */
/* > */
/* >  [3] B. Kagstrom and P. Poromaa, LAPACK-Style Algorithms and Software */
/* >      for Solving the Generalized Sylvester Equation and Estimating the */
/* >      Separation between Regular Matrix Pairs, Report UMINF - 93.23, */
/* >      Department of Computing Science, Umea University, S-901 87 Umea, */
/* >      Sweden, December 1993, Revised April 1994, Also as LAPACK Working */
/* >      Note 75. To appear in ACM Trans. on Math. Software, Vol 22, No 1, */
/* >      1996. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dtgsen_(integer *ijob, logical *wantq, logical *wantz, 
	logical *select, integer *n, doublereal *a, integer *lda, doublereal *
	b, integer *ldb, doublereal *alphar, doublereal *alphai, doublereal *
	beta, doublereal *q, integer *ldq, doublereal *z__, integer *ldz, 
	integer *m, doublereal *pl, doublereal *pr, doublereal *dif, 
	doublereal *work, integer *lwork, integer *iwork, integer *liwork, 
	integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, q_dim1, q_offset, z_dim1, 
	    z_offset, i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal), d_sign(doublereal *, doublereal *);

    /* Local variables */
    static integer i__, k, n1, n2, kk, ks, mn2, ijb;
    static doublereal eps;
    static integer kase;
    static logical pair;
    static integer ierr;
    static doublereal dsum;
    static logical swap;
    extern /* Subroutine */ int dlag2_(doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *);
    static integer isave[3];
    static logical wantd;
    static integer lwmin;
    static logical wantp;
    extern /* Subroutine */ int dlacn2_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, integer *);
    static logical wantd1, wantd2;
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal dscale, rdscal;
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen), dtgexc_(logical *, logical *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    integer *, doublereal *, integer *, integer *), dlassq_(integer *,
	     doublereal *, integer *, doublereal *, doublereal *);
    static integer liwmin;
    extern /* Subroutine */ int dtgsyl_(char *, integer *, integer *, integer 
	    *, doublereal *, integer *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     integer *, integer *, integer *, ftnlen);
    static doublereal smlnum;
    static logical lquery;


/*  -- LAPACK computational routine (version 3.4.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Decode and test the input parameters */

#line 507 "dtgsen.f"
    /* Parameter adjustments */
#line 507 "dtgsen.f"
    --select;
#line 507 "dtgsen.f"
    a_dim1 = *lda;
#line 507 "dtgsen.f"
    a_offset = 1 + a_dim1;
#line 507 "dtgsen.f"
    a -= a_offset;
#line 507 "dtgsen.f"
    b_dim1 = *ldb;
#line 507 "dtgsen.f"
    b_offset = 1 + b_dim1;
#line 507 "dtgsen.f"
    b -= b_offset;
#line 507 "dtgsen.f"
    --alphar;
#line 507 "dtgsen.f"
    --alphai;
#line 507 "dtgsen.f"
    --beta;
#line 507 "dtgsen.f"
    q_dim1 = *ldq;
#line 507 "dtgsen.f"
    q_offset = 1 + q_dim1;
#line 507 "dtgsen.f"
    q -= q_offset;
#line 507 "dtgsen.f"
    z_dim1 = *ldz;
#line 507 "dtgsen.f"
    z_offset = 1 + z_dim1;
#line 507 "dtgsen.f"
    z__ -= z_offset;
#line 507 "dtgsen.f"
    --dif;
#line 507 "dtgsen.f"
    --work;
#line 507 "dtgsen.f"
    --iwork;
#line 507 "dtgsen.f"

#line 507 "dtgsen.f"
    /* Function Body */
#line 507 "dtgsen.f"
    *info = 0;
#line 508 "dtgsen.f"
    lquery = *lwork == -1 || *liwork == -1;

#line 510 "dtgsen.f"
    if (*ijob < 0 || *ijob > 5) {
#line 511 "dtgsen.f"
	*info = -1;
#line 512 "dtgsen.f"
    } else if (*n < 0) {
#line 513 "dtgsen.f"
	*info = -5;
#line 514 "dtgsen.f"
    } else if (*lda < max(1,*n)) {
#line 515 "dtgsen.f"
	*info = -7;
#line 516 "dtgsen.f"
    } else if (*ldb < max(1,*n)) {
#line 517 "dtgsen.f"
	*info = -9;
#line 518 "dtgsen.f"
    } else if (*ldq < 1 || *wantq && *ldq < *n) {
#line 519 "dtgsen.f"
	*info = -14;
#line 520 "dtgsen.f"
    } else if (*ldz < 1 || *wantz && *ldz < *n) {
#line 521 "dtgsen.f"
	*info = -16;
#line 522 "dtgsen.f"
    }

#line 524 "dtgsen.f"
    if (*info != 0) {
#line 525 "dtgsen.f"
	i__1 = -(*info);
#line 525 "dtgsen.f"
	xerbla_("DTGSEN", &i__1, (ftnlen)6);
#line 526 "dtgsen.f"
	return 0;
#line 527 "dtgsen.f"
    }

/*     Get machine constants */

#line 531 "dtgsen.f"
    eps = dlamch_("P", (ftnlen)1);
#line 532 "dtgsen.f"
    smlnum = dlamch_("S", (ftnlen)1) / eps;
#line 533 "dtgsen.f"
    ierr = 0;

#line 535 "dtgsen.f"
    wantp = *ijob == 1 || *ijob >= 4;
#line 536 "dtgsen.f"
    wantd1 = *ijob == 2 || *ijob == 4;
#line 537 "dtgsen.f"
    wantd2 = *ijob == 3 || *ijob == 5;
#line 538 "dtgsen.f"
    wantd = wantd1 || wantd2;

/*     Set M to the dimension of the specified pair of deflating */
/*     subspaces. */

#line 543 "dtgsen.f"
    *m = 0;
#line 544 "dtgsen.f"
    pair = FALSE_;
#line 545 "dtgsen.f"
    i__1 = *n;
#line 545 "dtgsen.f"
    for (k = 1; k <= i__1; ++k) {
#line 546 "dtgsen.f"
	if (pair) {
#line 547 "dtgsen.f"
	    pair = FALSE_;
#line 548 "dtgsen.f"
	} else {
#line 549 "dtgsen.f"
	    if (k < *n) {
#line 550 "dtgsen.f"
		if (a[k + 1 + k * a_dim1] == 0.) {
#line 551 "dtgsen.f"
		    if (select[k]) {
#line 551 "dtgsen.f"
			++(*m);
#line 551 "dtgsen.f"
		    }
#line 553 "dtgsen.f"
		} else {
#line 554 "dtgsen.f"
		    pair = TRUE_;
#line 555 "dtgsen.f"
		    if (select[k] || select[k + 1]) {
#line 555 "dtgsen.f"
			*m += 2;
#line 555 "dtgsen.f"
		    }
#line 557 "dtgsen.f"
		}
#line 558 "dtgsen.f"
	    } else {
#line 559 "dtgsen.f"
		if (select[*n]) {
#line 559 "dtgsen.f"
		    ++(*m);
#line 559 "dtgsen.f"
		}
#line 561 "dtgsen.f"
	    }
#line 562 "dtgsen.f"
	}
#line 563 "dtgsen.f"
/* L10: */
#line 563 "dtgsen.f"
    }

#line 565 "dtgsen.f"
    if (*ijob == 1 || *ijob == 2 || *ijob == 4) {
/* Computing MAX */
#line 566 "dtgsen.f"
	i__1 = 1, i__2 = (*n << 2) + 16, i__1 = max(i__1,i__2), i__2 = (*m << 
		1) * (*n - *m);
#line 566 "dtgsen.f"
	lwmin = max(i__1,i__2);
/* Computing MAX */
#line 567 "dtgsen.f"
	i__1 = 1, i__2 = *n + 6;
#line 567 "dtgsen.f"
	liwmin = max(i__1,i__2);
#line 568 "dtgsen.f"
    } else if (*ijob == 3 || *ijob == 5) {
/* Computing MAX */
#line 569 "dtgsen.f"
	i__1 = 1, i__2 = (*n << 2) + 16, i__1 = max(i__1,i__2), i__2 = (*m << 
		2) * (*n - *m);
#line 569 "dtgsen.f"
	lwmin = max(i__1,i__2);
/* Computing MAX */
#line 570 "dtgsen.f"
	i__1 = 1, i__2 = (*m << 1) * (*n - *m), i__1 = max(i__1,i__2), i__2 = 
		*n + 6;
#line 570 "dtgsen.f"
	liwmin = max(i__1,i__2);
#line 571 "dtgsen.f"
    } else {
/* Computing MAX */
#line 572 "dtgsen.f"
	i__1 = 1, i__2 = (*n << 2) + 16;
#line 572 "dtgsen.f"
	lwmin = max(i__1,i__2);
#line 573 "dtgsen.f"
	liwmin = 1;
#line 574 "dtgsen.f"
    }

#line 576 "dtgsen.f"
    work[1] = (doublereal) lwmin;
#line 577 "dtgsen.f"
    iwork[1] = liwmin;

#line 579 "dtgsen.f"
    if (*lwork < lwmin && ! lquery) {
#line 580 "dtgsen.f"
	*info = -22;
#line 581 "dtgsen.f"
    } else if (*liwork < liwmin && ! lquery) {
#line 582 "dtgsen.f"
	*info = -24;
#line 583 "dtgsen.f"
    }

#line 585 "dtgsen.f"
    if (*info != 0) {
#line 586 "dtgsen.f"
	i__1 = -(*info);
#line 586 "dtgsen.f"
	xerbla_("DTGSEN", &i__1, (ftnlen)6);
#line 587 "dtgsen.f"
	return 0;
#line 588 "dtgsen.f"
    } else if (lquery) {
#line 589 "dtgsen.f"
	return 0;
#line 590 "dtgsen.f"
    }

/*     Quick return if possible. */

#line 594 "dtgsen.f"
    if (*m == *n || *m == 0) {
#line 595 "dtgsen.f"
	if (wantp) {
#line 596 "dtgsen.f"
	    *pl = 1.;
#line 597 "dtgsen.f"
	    *pr = 1.;
#line 598 "dtgsen.f"
	}
#line 599 "dtgsen.f"
	if (wantd) {
#line 600 "dtgsen.f"
	    dscale = 0.;
#line 601 "dtgsen.f"
	    dsum = 1.;
#line 602 "dtgsen.f"
	    i__1 = *n;
#line 602 "dtgsen.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 603 "dtgsen.f"
		dlassq_(n, &a[i__ * a_dim1 + 1], &c__1, &dscale, &dsum);
#line 604 "dtgsen.f"
		dlassq_(n, &b[i__ * b_dim1 + 1], &c__1, &dscale, &dsum);
#line 605 "dtgsen.f"
/* L20: */
#line 605 "dtgsen.f"
	    }
#line 606 "dtgsen.f"
	    dif[1] = dscale * sqrt(dsum);
#line 607 "dtgsen.f"
	    dif[2] = dif[1];
#line 608 "dtgsen.f"
	}
#line 609 "dtgsen.f"
	goto L60;
#line 610 "dtgsen.f"
    }

/*     Collect the selected blocks at the top-left corner of (A, B). */

#line 614 "dtgsen.f"
    ks = 0;
#line 615 "dtgsen.f"
    pair = FALSE_;
#line 616 "dtgsen.f"
    i__1 = *n;
#line 616 "dtgsen.f"
    for (k = 1; k <= i__1; ++k) {
#line 617 "dtgsen.f"
	if (pair) {
#line 618 "dtgsen.f"
	    pair = FALSE_;
#line 619 "dtgsen.f"
	} else {

#line 621 "dtgsen.f"
	    swap = select[k];
#line 622 "dtgsen.f"
	    if (k < *n) {
#line 623 "dtgsen.f"
		if (a[k + 1 + k * a_dim1] != 0.) {
#line 624 "dtgsen.f"
		    pair = TRUE_;
#line 625 "dtgsen.f"
		    swap = swap || select[k + 1];
#line 626 "dtgsen.f"
		}
#line 627 "dtgsen.f"
	    }

#line 629 "dtgsen.f"
	    if (swap) {
#line 630 "dtgsen.f"
		++ks;

/*              Swap the K-th block to position KS. */
/*              Perform the reordering of diagonal blocks in (A, B) */
/*              by orthogonal transformation matrices and update */
/*              Q and Z accordingly (if requested): */

#line 637 "dtgsen.f"
		kk = k;
#line 638 "dtgsen.f"
		if (k != ks) {
#line 638 "dtgsen.f"
		    dtgexc_(wantq, wantz, n, &a[a_offset], lda, &b[b_offset], 
			    ldb, &q[q_offset], ldq, &z__[z_offset], ldz, &kk, 
			    &ks, &work[1], lwork, &ierr);
#line 638 "dtgsen.f"
		}

#line 642 "dtgsen.f"
		if (ierr > 0) {

/*                 Swap is rejected: exit. */

#line 646 "dtgsen.f"
		    *info = 1;
#line 647 "dtgsen.f"
		    if (wantp) {
#line 648 "dtgsen.f"
			*pl = 0.;
#line 649 "dtgsen.f"
			*pr = 0.;
#line 650 "dtgsen.f"
		    }
#line 651 "dtgsen.f"
		    if (wantd) {
#line 652 "dtgsen.f"
			dif[1] = 0.;
#line 653 "dtgsen.f"
			dif[2] = 0.;
#line 654 "dtgsen.f"
		    }
#line 655 "dtgsen.f"
		    goto L60;
#line 656 "dtgsen.f"
		}

#line 658 "dtgsen.f"
		if (pair) {
#line 658 "dtgsen.f"
		    ++ks;
#line 658 "dtgsen.f"
		}
#line 660 "dtgsen.f"
	    }
#line 661 "dtgsen.f"
	}
#line 662 "dtgsen.f"
/* L30: */
#line 662 "dtgsen.f"
    }
#line 663 "dtgsen.f"
    if (wantp) {

/*        Solve generalized Sylvester equation for R and L */
/*        and compute PL and PR. */

#line 668 "dtgsen.f"
	n1 = *m;
#line 669 "dtgsen.f"
	n2 = *n - *m;
#line 670 "dtgsen.f"
	i__ = n1 + 1;
#line 671 "dtgsen.f"
	ijb = 0;
#line 672 "dtgsen.f"
	dlacpy_("Full", &n1, &n2, &a[i__ * a_dim1 + 1], lda, &work[1], &n1, (
		ftnlen)4);
#line 673 "dtgsen.f"
	dlacpy_("Full", &n1, &n2, &b[i__ * b_dim1 + 1], ldb, &work[n1 * n2 + 
		1], &n1, (ftnlen)4);
#line 675 "dtgsen.f"
	i__1 = *lwork - (n1 << 1) * n2;
#line 675 "dtgsen.f"
	dtgsyl_("N", &ijb, &n1, &n2, &a[a_offset], lda, &a[i__ + i__ * a_dim1]
		, lda, &work[1], &n1, &b[b_offset], ldb, &b[i__ + i__ * 
		b_dim1], ldb, &work[n1 * n2 + 1], &n1, &dscale, &dif[1], &
		work[(n1 * n2 << 1) + 1], &i__1, &iwork[1], &ierr, (ftnlen)1);

/*        Estimate the reciprocal of norms of "projections" onto left */
/*        and right eigenspaces. */

#line 683 "dtgsen.f"
	rdscal = 0.;
#line 684 "dtgsen.f"
	dsum = 1.;
#line 685 "dtgsen.f"
	i__1 = n1 * n2;
#line 685 "dtgsen.f"
	dlassq_(&i__1, &work[1], &c__1, &rdscal, &dsum);
#line 686 "dtgsen.f"
	*pl = rdscal * sqrt(dsum);
#line 687 "dtgsen.f"
	if (*pl == 0.) {
#line 688 "dtgsen.f"
	    *pl = 1.;
#line 689 "dtgsen.f"
	} else {
#line 690 "dtgsen.f"
	    *pl = dscale / (sqrt(dscale * dscale / *pl + *pl) * sqrt(*pl));
#line 691 "dtgsen.f"
	}
#line 692 "dtgsen.f"
	rdscal = 0.;
#line 693 "dtgsen.f"
	dsum = 1.;
#line 694 "dtgsen.f"
	i__1 = n1 * n2;
#line 694 "dtgsen.f"
	dlassq_(&i__1, &work[n1 * n2 + 1], &c__1, &rdscal, &dsum);
#line 695 "dtgsen.f"
	*pr = rdscal * sqrt(dsum);
#line 696 "dtgsen.f"
	if (*pr == 0.) {
#line 697 "dtgsen.f"
	    *pr = 1.;
#line 698 "dtgsen.f"
	} else {
#line 699 "dtgsen.f"
	    *pr = dscale / (sqrt(dscale * dscale / *pr + *pr) * sqrt(*pr));
#line 700 "dtgsen.f"
	}
#line 701 "dtgsen.f"
    }

#line 703 "dtgsen.f"
    if (wantd) {

/*        Compute estimates of Difu and Difl. */

#line 707 "dtgsen.f"
	if (wantd1) {
#line 708 "dtgsen.f"
	    n1 = *m;
#line 709 "dtgsen.f"
	    n2 = *n - *m;
#line 710 "dtgsen.f"
	    i__ = n1 + 1;
#line 711 "dtgsen.f"
	    ijb = 3;

/*           Frobenius norm-based Difu-estimate. */

#line 715 "dtgsen.f"
	    i__1 = *lwork - (n1 << 1) * n2;
#line 715 "dtgsen.f"
	    dtgsyl_("N", &ijb, &n1, &n2, &a[a_offset], lda, &a[i__ + i__ * 
		    a_dim1], lda, &work[1], &n1, &b[b_offset], ldb, &b[i__ + 
		    i__ * b_dim1], ldb, &work[n1 * n2 + 1], &n1, &dscale, &
		    dif[1], &work[(n1 << 1) * n2 + 1], &i__1, &iwork[1], &
		    ierr, (ftnlen)1);

/*           Frobenius norm-based Difl-estimate. */

#line 722 "dtgsen.f"
	    i__1 = *lwork - (n1 << 1) * n2;
#line 722 "dtgsen.f"
	    dtgsyl_("N", &ijb, &n2, &n1, &a[i__ + i__ * a_dim1], lda, &a[
		    a_offset], lda, &work[1], &n2, &b[i__ + i__ * b_dim1], 
		    ldb, &b[b_offset], ldb, &work[n1 * n2 + 1], &n2, &dscale, 
		    &dif[2], &work[(n1 << 1) * n2 + 1], &i__1, &iwork[1], &
		    ierr, (ftnlen)1);
#line 726 "dtgsen.f"
	} else {


/*           Compute 1-norm-based estimates of Difu and Difl using */
/*           reversed communication with DLACN2. In each step a */
/*           generalized Sylvester equation or a transposed variant */
/*           is solved. */

#line 734 "dtgsen.f"
	    kase = 0;
#line 735 "dtgsen.f"
	    n1 = *m;
#line 736 "dtgsen.f"
	    n2 = *n - *m;
#line 737 "dtgsen.f"
	    i__ = n1 + 1;
#line 738 "dtgsen.f"
	    ijb = 0;
#line 739 "dtgsen.f"
	    mn2 = (n1 << 1) * n2;

/*           1-norm-based estimate of Difu. */

#line 743 "dtgsen.f"
L40:
#line 744 "dtgsen.f"
	    dlacn2_(&mn2, &work[mn2 + 1], &work[1], &iwork[1], &dif[1], &kase,
		     isave);
#line 746 "dtgsen.f"
	    if (kase != 0) {
#line 747 "dtgsen.f"
		if (kase == 1) {

/*                 Solve generalized Sylvester equation. */

#line 751 "dtgsen.f"
		    i__1 = *lwork - (n1 << 1) * n2;
#line 751 "dtgsen.f"
		    dtgsyl_("N", &ijb, &n1, &n2, &a[a_offset], lda, &a[i__ + 
			    i__ * a_dim1], lda, &work[1], &n1, &b[b_offset], 
			    ldb, &b[i__ + i__ * b_dim1], ldb, &work[n1 * n2 + 
			    1], &n1, &dscale, &dif[1], &work[(n1 << 1) * n2 + 
			    1], &i__1, &iwork[1], &ierr, (ftnlen)1);
#line 756 "dtgsen.f"
		} else {

/*                 Solve the transposed variant. */

#line 760 "dtgsen.f"
		    i__1 = *lwork - (n1 << 1) * n2;
#line 760 "dtgsen.f"
		    dtgsyl_("T", &ijb, &n1, &n2, &a[a_offset], lda, &a[i__ + 
			    i__ * a_dim1], lda, &work[1], &n1, &b[b_offset], 
			    ldb, &b[i__ + i__ * b_dim1], ldb, &work[n1 * n2 + 
			    1], &n1, &dscale, &dif[1], &work[(n1 << 1) * n2 + 
			    1], &i__1, &iwork[1], &ierr, (ftnlen)1);
#line 765 "dtgsen.f"
		}
#line 766 "dtgsen.f"
		goto L40;
#line 767 "dtgsen.f"
	    }
#line 768 "dtgsen.f"
	    dif[1] = dscale / dif[1];

/*           1-norm-based estimate of Difl. */

#line 772 "dtgsen.f"
L50:
#line 773 "dtgsen.f"
	    dlacn2_(&mn2, &work[mn2 + 1], &work[1], &iwork[1], &dif[2], &kase,
		     isave);
#line 775 "dtgsen.f"
	    if (kase != 0) {
#line 776 "dtgsen.f"
		if (kase == 1) {

/*                 Solve generalized Sylvester equation. */

#line 780 "dtgsen.f"
		    i__1 = *lwork - (n1 << 1) * n2;
#line 780 "dtgsen.f"
		    dtgsyl_("N", &ijb, &n2, &n1, &a[i__ + i__ * a_dim1], lda, 
			    &a[a_offset], lda, &work[1], &n2, &b[i__ + i__ * 
			    b_dim1], ldb, &b[b_offset], ldb, &work[n1 * n2 + 
			    1], &n2, &dscale, &dif[2], &work[(n1 << 1) * n2 + 
			    1], &i__1, &iwork[1], &ierr, (ftnlen)1);
#line 785 "dtgsen.f"
		} else {

/*                 Solve the transposed variant. */

#line 789 "dtgsen.f"
		    i__1 = *lwork - (n1 << 1) * n2;
#line 789 "dtgsen.f"
		    dtgsyl_("T", &ijb, &n2, &n1, &a[i__ + i__ * a_dim1], lda, 
			    &a[a_offset], lda, &work[1], &n2, &b[i__ + i__ * 
			    b_dim1], ldb, &b[b_offset], ldb, &work[n1 * n2 + 
			    1], &n2, &dscale, &dif[2], &work[(n1 << 1) * n2 + 
			    1], &i__1, &iwork[1], &ierr, (ftnlen)1);
#line 794 "dtgsen.f"
		}
#line 795 "dtgsen.f"
		goto L50;
#line 796 "dtgsen.f"
	    }
#line 797 "dtgsen.f"
	    dif[2] = dscale / dif[2];

#line 799 "dtgsen.f"
	}
#line 800 "dtgsen.f"
    }

#line 802 "dtgsen.f"
L60:

/*     Compute generalized eigenvalues of reordered pair (A, B) and */
/*     normalize the generalized Schur form. */

#line 807 "dtgsen.f"
    pair = FALSE_;
#line 808 "dtgsen.f"
    i__1 = *n;
#line 808 "dtgsen.f"
    for (k = 1; k <= i__1; ++k) {
#line 809 "dtgsen.f"
	if (pair) {
#line 810 "dtgsen.f"
	    pair = FALSE_;
#line 811 "dtgsen.f"
	} else {

#line 813 "dtgsen.f"
	    if (k < *n) {
#line 814 "dtgsen.f"
		if (a[k + 1 + k * a_dim1] != 0.) {
#line 815 "dtgsen.f"
		    pair = TRUE_;
#line 816 "dtgsen.f"
		}
#line 817 "dtgsen.f"
	    }

#line 819 "dtgsen.f"
	    if (pair) {

/*             Compute the eigenvalue(s) at position K. */

#line 823 "dtgsen.f"
		work[1] = a[k + k * a_dim1];
#line 824 "dtgsen.f"
		work[2] = a[k + 1 + k * a_dim1];
#line 825 "dtgsen.f"
		work[3] = a[k + (k + 1) * a_dim1];
#line 826 "dtgsen.f"
		work[4] = a[k + 1 + (k + 1) * a_dim1];
#line 827 "dtgsen.f"
		work[5] = b[k + k * b_dim1];
#line 828 "dtgsen.f"
		work[6] = b[k + 1 + k * b_dim1];
#line 829 "dtgsen.f"
		work[7] = b[k + (k + 1) * b_dim1];
#line 830 "dtgsen.f"
		work[8] = b[k + 1 + (k + 1) * b_dim1];
#line 831 "dtgsen.f"
		d__1 = smlnum * eps;
#line 831 "dtgsen.f"
		dlag2_(&work[1], &c__2, &work[5], &c__2, &d__1, &beta[k], &
			beta[k + 1], &alphar[k], &alphar[k + 1], &alphai[k]);
#line 834 "dtgsen.f"
		alphai[k + 1] = -alphai[k];

#line 836 "dtgsen.f"
	    } else {

#line 838 "dtgsen.f"
		if (d_sign(&c_b28, &b[k + k * b_dim1]) < 0.) {

/*                 If B(K,K) is negative, make it positive */

#line 842 "dtgsen.f"
		    i__2 = *n;
#line 842 "dtgsen.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 843 "dtgsen.f"
			a[k + i__ * a_dim1] = -a[k + i__ * a_dim1];
#line 844 "dtgsen.f"
			b[k + i__ * b_dim1] = -b[k + i__ * b_dim1];
#line 845 "dtgsen.f"
			if (*wantq) {
#line 845 "dtgsen.f"
			    q[i__ + k * q_dim1] = -q[i__ + k * q_dim1];
#line 845 "dtgsen.f"
			}
#line 846 "dtgsen.f"
/* L70: */
#line 846 "dtgsen.f"
		    }
#line 847 "dtgsen.f"
		}

#line 849 "dtgsen.f"
		alphar[k] = a[k + k * a_dim1];
#line 850 "dtgsen.f"
		alphai[k] = 0.;
#line 851 "dtgsen.f"
		beta[k] = b[k + k * b_dim1];

#line 853 "dtgsen.f"
	    }
#line 854 "dtgsen.f"
	}
#line 855 "dtgsen.f"
/* L80: */
#line 855 "dtgsen.f"
    }

#line 857 "dtgsen.f"
    work[1] = (doublereal) lwmin;
#line 858 "dtgsen.f"
    iwork[1] = liwmin;

#line 860 "dtgsen.f"
    return 0;

/*     End of DTGSEN */

} /* dtgsen_ */

