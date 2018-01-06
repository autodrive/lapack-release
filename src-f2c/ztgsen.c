#line 1 "ztgsen.f"
/* ztgsen.f -- translated by f2c (version 20100827).
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

#line 1 "ztgsen.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b ZTGSEN */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZTGSEN + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztgsen.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztgsen.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztgsen.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZTGSEN( IJOB, WANTQ, WANTZ, SELECT, N, A, LDA, B, LDB, */
/*                          ALPHA, BETA, Q, LDQ, Z, LDZ, M, PL, PR, DIF, */
/*                          WORK, LWORK, IWORK, LIWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       LOGICAL            WANTQ, WANTZ */
/*       INTEGER            IJOB, INFO, LDA, LDB, LDQ, LDZ, LIWORK, LWORK, */
/*      $                   M, N */
/*       DOUBLE PRECISION   PL, PR */
/*       .. */
/*       .. Array Arguments .. */
/*       LOGICAL            SELECT( * ) */
/*       INTEGER            IWORK( * ) */
/*       DOUBLE PRECISION   DIF( * ) */
/*       COMPLEX*16         A( LDA, * ), ALPHA( * ), B( LDB, * ), */
/*      $                   BETA( * ), Q( LDQ, * ), WORK( * ), Z( LDZ, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZTGSEN reorders the generalized Schur decomposition of a complex */
/* > matrix pair (A, B) (in terms of an unitary equivalence trans- */
/* > formation Q**H * (A, B) * Z), so that a selected cluster of eigenvalues */
/* > appears in the leading diagonal blocks of the pair (A,B). The leading */
/* > columns of Q and Z form unitary bases of the corresponding left and */
/* > right eigenspaces (deflating subspaces). (A, B) must be in */
/* > generalized Schur canonical form, that is, A and B are both upper */
/* > triangular. */
/* > */
/* > ZTGSEN also computes the generalized eigenvalues */
/* > */
/* >          w(j)= ALPHA(j) / BETA(j) */
/* > */
/* > of the reordered matrix pair (A, B). */
/* > */
/* > Optionally, the routine computes estimates of reciprocal condition */
/* > numbers for eigenvalues and eigenspaces. These are Difu[(A11,B11), */
/* > (A22,B22)] and Difl[(A11,B11), (A22,B22)], i.e. the separation(s) */
/* > between the matrix pairs (A11, B11) and (A22,B22) that correspond to */
/* > the selected cluster and the eigenvalues outside the cluster, resp., */
/* > and norms of "projections" onto left and right eigenspaces w.r.t. */
/* > the selected cluster in the (1,1)-block. */
/* > */
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
/* >          SELECT specifies the eigenvalues in the selected cluster. To */
/* >          select an eigenvalue w(j), SELECT(j) must be set to */
/* >          .TRUE.. */
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
/* >          A is COMPLEX*16 array, dimension(LDA,N) */
/* >          On entry, the upper triangular matrix A, in generalized */
/* >          Schur canonical form. */
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
/* >          B is COMPLEX*16 array, dimension(LDB,N) */
/* >          On entry, the upper triangular matrix B, in generalized */
/* >          Schur canonical form. */
/* >          On exit, B is overwritten by the reordered matrix B. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >          The leading dimension of the array B. LDB >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] ALPHA */
/* > \verbatim */
/* >          ALPHA is COMPLEX*16 array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] BETA */
/* > \verbatim */
/* >          BETA is COMPLEX*16 array, dimension (N) */
/* > */
/* >          The diagonal elements of A and B, respectively, */
/* >          when the pair (A,B) has been reduced to generalized Schur */
/* >          form.  ALPHA(i)/BETA(i) i=1,...,N are the generalized */
/* >          eigenvalues. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Q */
/* > \verbatim */
/* >          Q is COMPLEX*16 array, dimension (LDQ,N) */
/* >          On entry, if WANTQ = .TRUE., Q is an N-by-N matrix. */
/* >          On exit, Q has been postmultiplied by the left unitary */
/* >          transformation matrix which reorder (A, B); The leading M */
/* >          columns of Q form orthonormal bases for the specified pair of */
/* >          left eigenspaces (deflating subspaces). */
/* >          If WANTQ = .FALSE., Q is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDQ */
/* > \verbatim */
/* >          LDQ is INTEGER */
/* >          The leading dimension of the array Q. LDQ >= 1. */
/* >          If WANTQ = .TRUE., LDQ >= N. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Z */
/* > \verbatim */
/* >          Z is COMPLEX*16 array, dimension (LDZ,N) */
/* >          On entry, if WANTZ = .TRUE., Z is an N-by-N matrix. */
/* >          On exit, Z has been postmultiplied by the left unitary */
/* >          transformation matrix which reorder (A, B); The leading M */
/* >          columns of Z form orthonormal bases for the specified pair of */
/* >          left eigenspaces (deflating subspaces). */
/* >          If WANTZ = .FALSE., Z is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDZ */
/* > \verbatim */
/* >          LDZ is INTEGER */
/* >          The leading dimension of the array Z. LDZ >= 1. */
/* >          If WANTZ = .TRUE., LDZ >= N. */
/* > \endverbatim */
/* > */
/* > \param[out] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The dimension of the specified pair of left and right */
/* >          eigenspaces, (deflating subspaces) 0 <= M <= N. */
/* > \endverbatim */
/* > */
/* > \param[out] PL */
/* > \verbatim */
/* >          PL is DOUBLE PRECISION */
/* > \endverbatim */
/* > */
/* > \param[out] PR */
/* > \verbatim */
/* >          PR is DOUBLE PRECISION */
/* > */
/* >          If IJOB = 1, 4 or 5, PL, PR are lower bounds on the */
/* >          reciprocal  of the norm of "projections" onto left and right */
/* >          eigenspace with respect to the selected cluster. */
/* >          0 < PL, PR <= 1. */
/* >          If M = 0 or M = N, PL = PR  = 1. */
/* >          If IJOB = 0, 2 or 3 PL, PR are not referenced. */
/* > \endverbatim */
/* > */
/* > \param[out] DIF */
/* > \verbatim */
/* >          DIF is DOUBLE PRECISION array, dimension (2). */
/* >          If IJOB >= 2, DIF(1:2) store the estimates of Difu and Difl. */
/* >          If IJOB = 2 or 4, DIF(1:2) are F-norm-based upper bounds on */
/* >          Difu and Difl. If IJOB = 3 or 5, DIF(1:2) are 1-norm-based */
/* >          estimates of Difu and Difl, computed using reversed */
/* >          communication with ZLACN2. */
/* >          If M = 0 or N, DIF(1:2) = F-norm([A, B]). */
/* >          If IJOB = 0 or 1, DIF is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The dimension of the array WORK. LWORK >=  1 */
/* >          If IJOB = 1, 2 or 4, LWORK >=  2*M*(N-M) */
/* >          If IJOB = 3 or 5, LWORK >=  4*M*(N-M) */
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
/* >          If IJOB = 1, 2 or 4, LIWORK >=  N+2; */
/* >          If IJOB = 3 or 5, LIWORK >= MAX(N+2, 2*M*(N-M)); */
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

/* > \date June 2016 */

/* > \ingroup complex16OTHERcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  ZTGSEN first collects the selected eigenvalues by computing unitary */
/* >  U and W that move them to the top left corner of (A, B). In other */
/* >  words, the selected eigenvalues are the eigenvalues of (A11, B11) in */
/* > */
/* >              U**H*(A, B)*W = (A11 A12) (B11 B12) n1 */
/* >                              ( 0  A22),( 0  B22) n2 */
/* >                                n1  n2    n1  n2 */
/* > */
/* >  where N = n1+n2 and U**H means the conjugate transpose of U. The first */
/* >  n1 columns of U and W span the specified pair of left and right */
/* >  eigenspaces (deflating subspaces) of (A, B). */
/* > */
/* >  If (A, B) has been obtained from the generalized real Schur */
/* >  decomposition of a matrix pair (C, D) = Q*(A, B)*Z**H, then the */
/* >  reordered generalized Schur form of (C, D) is given by */
/* > */
/* >           (C, D) = (Q*U)*(U**H *(A, B)*W)*(Z*W)**H, */
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
/* >       Zu = [ kron(In2, A11)  -kron(A22**H, In1) ] */
/* >            [ kron(In2, B11)  -kron(B22**H, In1) ]. */
/* > */
/* >  Here, Inx is the identity matrix of size nx and A22**H is the */
/* >  conjugate transpose of A22. kron(X, Y) is the Kronecker product between */
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
/* >  based estimate DIF is not wanted (see ZLATDF), then the parameter */
/* >  IDIFJB (see below) should be changed from 3 to 4 (routine ZLATDF */
/* >  (IJOB = 2 will be used)). See ZTGSYL for more details. */
/* > \endverbatim */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Bo Kagstrom and Peter Poromaa, Department of Computing Science, */
/* >     Umea University, S-901 87 Umea, Sweden. */

/* > \par References: */
/*  ================ */
/* > */
/* >  [1] B. Kagstrom; A Direct Method for Reordering Eigenvalues in the */
/* >      Generalized Real Schur Form of a Regular Matrix Pair (A, B), in */
/* >      M.S. Moonen et al (eds), Linear Algebra for Large Scale and */
/* >      Real-Time Applications, Kluwer Academic Publ. 1993, pp 195-218. */
/* > \n */
/* >  [2] B. Kagstrom and P. Poromaa; Computing Eigenspaces with Specified */
/* >      Eigenvalues of a Regular Matrix Pair (A, B) and Condition */
/* >      Estimation: Theory, Algorithms and Software, Report */
/* >      UMINF - 94.04, Department of Computing Science, Umea University, */
/* >      S-901 87 Umea, Sweden, 1994. Also as LAPACK Working Note 87. */
/* >      To appear in Numerical Algorithms, 1996. */
/* > \n */
/* >  [3] B. Kagstrom and P. Poromaa, LAPACK-Style Algorithms and Software */
/* >      for Solving the Generalized Sylvester Equation and Estimating the */
/* >      Separation between Regular Matrix Pairs, Report UMINF - 93.23, */
/* >      Department of Computing Science, Umea University, S-901 87 Umea, */
/* >      Sweden, December 1993, Revised April 1994, Also as LAPACK working */
/* >      Note 75. To appear in ACM Trans. on Math. Software, Vol 22, No 1, */
/* >      1996. */
/* > */
/*  ===================================================================== */
/* Subroutine */ int ztgsen_(integer *ijob, logical *wantq, logical *wantz, 
	logical *select, integer *n, doublecomplex *a, integer *lda, 
	doublecomplex *b, integer *ldb, doublecomplex *alpha, doublecomplex *
	beta, doublecomplex *q, integer *ldq, doublecomplex *z__, integer *
	ldz, integer *m, doublereal *pl, doublereal *pr, doublereal *dif, 
	doublecomplex *work, integer *lwork, integer *iwork, integer *liwork, 
	integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, q_dim1, q_offset, z_dim1, 
	    z_offset, i__1, i__2, i__3;
    doublecomplex z__1, z__2;

    /* Builtin functions */
    double sqrt(doublereal), z_abs(doublecomplex *);
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, k, n1, n2, ks, mn2, ijb, kase, ierr;
    static doublereal dsum;
    static logical swap;
    static doublecomplex temp1, temp2;
    static integer isave[3];
    extern /* Subroutine */ int zscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *);
    static logical wantd;
    static integer lwmin;
    static logical wantp;
    extern /* Subroutine */ int zlacn2_(integer *, doublecomplex *, 
	    doublecomplex *, doublereal *, integer *, integer *);
    static logical wantd1, wantd2;
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal dscale, rdscal, safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static integer liwmin;
    extern /* Subroutine */ int zlacpy_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen), 
	    ztgexc_(logical *, logical *, integer *, doublecomplex *, integer 
	    *, doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, integer *, integer *, integer *), 
	    zlassq_(integer *, doublecomplex *, integer *, doublereal *, 
	    doublereal *);
    static logical lquery;
    extern /* Subroutine */ int ztgsyl_(char *, integer *, integer *, integer 
	    *, doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublereal *, doublereal *, doublecomplex *, integer *, integer *,
	     integer *, ftnlen);


/*  -- LAPACK computational routine (version 3.7.1) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     June 2016 */

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
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Decode and test the input parameters */

#line 488 "ztgsen.f"
    /* Parameter adjustments */
#line 488 "ztgsen.f"
    --select;
#line 488 "ztgsen.f"
    a_dim1 = *lda;
#line 488 "ztgsen.f"
    a_offset = 1 + a_dim1;
#line 488 "ztgsen.f"
    a -= a_offset;
#line 488 "ztgsen.f"
    b_dim1 = *ldb;
#line 488 "ztgsen.f"
    b_offset = 1 + b_dim1;
#line 488 "ztgsen.f"
    b -= b_offset;
#line 488 "ztgsen.f"
    --alpha;
#line 488 "ztgsen.f"
    --beta;
#line 488 "ztgsen.f"
    q_dim1 = *ldq;
#line 488 "ztgsen.f"
    q_offset = 1 + q_dim1;
#line 488 "ztgsen.f"
    q -= q_offset;
#line 488 "ztgsen.f"
    z_dim1 = *ldz;
#line 488 "ztgsen.f"
    z_offset = 1 + z_dim1;
#line 488 "ztgsen.f"
    z__ -= z_offset;
#line 488 "ztgsen.f"
    --dif;
#line 488 "ztgsen.f"
    --work;
#line 488 "ztgsen.f"
    --iwork;
#line 488 "ztgsen.f"

#line 488 "ztgsen.f"
    /* Function Body */
#line 488 "ztgsen.f"
    *info = 0;
#line 489 "ztgsen.f"
    lquery = *lwork == -1 || *liwork == -1;

#line 491 "ztgsen.f"
    if (*ijob < 0 || *ijob > 5) {
#line 492 "ztgsen.f"
	*info = -1;
#line 493 "ztgsen.f"
    } else if (*n < 0) {
#line 494 "ztgsen.f"
	*info = -5;
#line 495 "ztgsen.f"
    } else if (*lda < max(1,*n)) {
#line 496 "ztgsen.f"
	*info = -7;
#line 497 "ztgsen.f"
    } else if (*ldb < max(1,*n)) {
#line 498 "ztgsen.f"
	*info = -9;
#line 499 "ztgsen.f"
    } else if (*ldq < 1 || *wantq && *ldq < *n) {
#line 500 "ztgsen.f"
	*info = -13;
#line 501 "ztgsen.f"
    } else if (*ldz < 1 || *wantz && *ldz < *n) {
#line 502 "ztgsen.f"
	*info = -15;
#line 503 "ztgsen.f"
    }

#line 505 "ztgsen.f"
    if (*info != 0) {
#line 506 "ztgsen.f"
	i__1 = -(*info);
#line 506 "ztgsen.f"
	xerbla_("ZTGSEN", &i__1, (ftnlen)6);
#line 507 "ztgsen.f"
	return 0;
#line 508 "ztgsen.f"
    }

#line 510 "ztgsen.f"
    ierr = 0;

#line 512 "ztgsen.f"
    wantp = *ijob == 1 || *ijob >= 4;
#line 513 "ztgsen.f"
    wantd1 = *ijob == 2 || *ijob == 4;
#line 514 "ztgsen.f"
    wantd2 = *ijob == 3 || *ijob == 5;
#line 515 "ztgsen.f"
    wantd = wantd1 || wantd2;

/*     Set M to the dimension of the specified pair of deflating */
/*     subspaces. */

#line 520 "ztgsen.f"
    *m = 0;
#line 521 "ztgsen.f"
    if (! lquery || *ijob != 0) {
#line 522 "ztgsen.f"
	i__1 = *n;
#line 522 "ztgsen.f"
	for (k = 1; k <= i__1; ++k) {
#line 523 "ztgsen.f"
	    i__2 = k;
#line 523 "ztgsen.f"
	    i__3 = k + k * a_dim1;
#line 523 "ztgsen.f"
	    alpha[i__2].r = a[i__3].r, alpha[i__2].i = a[i__3].i;
#line 524 "ztgsen.f"
	    i__2 = k;
#line 524 "ztgsen.f"
	    i__3 = k + k * b_dim1;
#line 524 "ztgsen.f"
	    beta[i__2].r = b[i__3].r, beta[i__2].i = b[i__3].i;
#line 525 "ztgsen.f"
	    if (k < *n) {
#line 526 "ztgsen.f"
		if (select[k]) {
#line 526 "ztgsen.f"
		    ++(*m);
#line 526 "ztgsen.f"
		}
#line 528 "ztgsen.f"
	    } else {
#line 529 "ztgsen.f"
		if (select[*n]) {
#line 529 "ztgsen.f"
		    ++(*m);
#line 529 "ztgsen.f"
		}
#line 531 "ztgsen.f"
	    }
#line 532 "ztgsen.f"
/* L10: */
#line 532 "ztgsen.f"
	}
#line 533 "ztgsen.f"
    }

#line 535 "ztgsen.f"
    if (*ijob == 1 || *ijob == 2 || *ijob == 4) {
/* Computing MAX */
#line 536 "ztgsen.f"
	i__1 = 1, i__2 = (*m << 1) * (*n - *m);
#line 536 "ztgsen.f"
	lwmin = max(i__1,i__2);
/* Computing MAX */
#line 537 "ztgsen.f"
	i__1 = 1, i__2 = *n + 2;
#line 537 "ztgsen.f"
	liwmin = max(i__1,i__2);
#line 538 "ztgsen.f"
    } else if (*ijob == 3 || *ijob == 5) {
/* Computing MAX */
#line 539 "ztgsen.f"
	i__1 = 1, i__2 = (*m << 2) * (*n - *m);
#line 539 "ztgsen.f"
	lwmin = max(i__1,i__2);
/* Computing MAX */
#line 540 "ztgsen.f"
	i__1 = 1, i__2 = (*m << 1) * (*n - *m), i__1 = max(i__1,i__2), i__2 = 
		*n + 2;
#line 540 "ztgsen.f"
	liwmin = max(i__1,i__2);
#line 541 "ztgsen.f"
    } else {
#line 542 "ztgsen.f"
	lwmin = 1;
#line 543 "ztgsen.f"
	liwmin = 1;
#line 544 "ztgsen.f"
    }

#line 546 "ztgsen.f"
    work[1].r = (doublereal) lwmin, work[1].i = 0.;
#line 547 "ztgsen.f"
    iwork[1] = liwmin;

#line 549 "ztgsen.f"
    if (*lwork < lwmin && ! lquery) {
#line 550 "ztgsen.f"
	*info = -21;
#line 551 "ztgsen.f"
    } else if (*liwork < liwmin && ! lquery) {
#line 552 "ztgsen.f"
	*info = -23;
#line 553 "ztgsen.f"
    }

#line 555 "ztgsen.f"
    if (*info != 0) {
#line 556 "ztgsen.f"
	i__1 = -(*info);
#line 556 "ztgsen.f"
	xerbla_("ZTGSEN", &i__1, (ftnlen)6);
#line 557 "ztgsen.f"
	return 0;
#line 558 "ztgsen.f"
    } else if (lquery) {
#line 559 "ztgsen.f"
	return 0;
#line 560 "ztgsen.f"
    }

/*     Quick return if possible. */

#line 564 "ztgsen.f"
    if (*m == *n || *m == 0) {
#line 565 "ztgsen.f"
	if (wantp) {
#line 566 "ztgsen.f"
	    *pl = 1.;
#line 567 "ztgsen.f"
	    *pr = 1.;
#line 568 "ztgsen.f"
	}
#line 569 "ztgsen.f"
	if (wantd) {
#line 570 "ztgsen.f"
	    dscale = 0.;
#line 571 "ztgsen.f"
	    dsum = 1.;
#line 572 "ztgsen.f"
	    i__1 = *n;
#line 572 "ztgsen.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 573 "ztgsen.f"
		zlassq_(n, &a[i__ * a_dim1 + 1], &c__1, &dscale, &dsum);
#line 574 "ztgsen.f"
		zlassq_(n, &b[i__ * b_dim1 + 1], &c__1, &dscale, &dsum);
#line 575 "ztgsen.f"
/* L20: */
#line 575 "ztgsen.f"
	    }
#line 576 "ztgsen.f"
	    dif[1] = dscale * sqrt(dsum);
#line 577 "ztgsen.f"
	    dif[2] = dif[1];
#line 578 "ztgsen.f"
	}
#line 579 "ztgsen.f"
	goto L70;
#line 580 "ztgsen.f"
    }

/*     Get machine constant */

#line 584 "ztgsen.f"
    safmin = dlamch_("S", (ftnlen)1);

/*     Collect the selected blocks at the top-left corner of (A, B). */

#line 588 "ztgsen.f"
    ks = 0;
#line 589 "ztgsen.f"
    i__1 = *n;
#line 589 "ztgsen.f"
    for (k = 1; k <= i__1; ++k) {
#line 590 "ztgsen.f"
	swap = select[k];
#line 591 "ztgsen.f"
	if (swap) {
#line 592 "ztgsen.f"
	    ++ks;

/*           Swap the K-th block to position KS. Compute unitary Q */
/*           and Z that will swap adjacent diagonal blocks in (A, B). */

#line 597 "ztgsen.f"
	    if (k != ks) {
#line 597 "ztgsen.f"
		ztgexc_(wantq, wantz, n, &a[a_offset], lda, &b[b_offset], ldb,
			 &q[q_offset], ldq, &z__[z_offset], ldz, &k, &ks, &
			ierr);
#line 597 "ztgsen.f"
	    }

#line 601 "ztgsen.f"
	    if (ierr > 0) {

/*              Swap is rejected: exit. */

#line 605 "ztgsen.f"
		*info = 1;
#line 606 "ztgsen.f"
		if (wantp) {
#line 607 "ztgsen.f"
		    *pl = 0.;
#line 608 "ztgsen.f"
		    *pr = 0.;
#line 609 "ztgsen.f"
		}
#line 610 "ztgsen.f"
		if (wantd) {
#line 611 "ztgsen.f"
		    dif[1] = 0.;
#line 612 "ztgsen.f"
		    dif[2] = 0.;
#line 613 "ztgsen.f"
		}
#line 614 "ztgsen.f"
		goto L70;
#line 615 "ztgsen.f"
	    }
#line 616 "ztgsen.f"
	}
#line 617 "ztgsen.f"
/* L30: */
#line 617 "ztgsen.f"
    }
#line 618 "ztgsen.f"
    if (wantp) {

/*        Solve generalized Sylvester equation for R and L: */
/*                   A11 * R - L * A22 = A12 */
/*                   B11 * R - L * B22 = B12 */

#line 624 "ztgsen.f"
	n1 = *m;
#line 625 "ztgsen.f"
	n2 = *n - *m;
#line 626 "ztgsen.f"
	i__ = n1 + 1;
#line 627 "ztgsen.f"
	zlacpy_("Full", &n1, &n2, &a[i__ * a_dim1 + 1], lda, &work[1], &n1, (
		ftnlen)4);
#line 628 "ztgsen.f"
	zlacpy_("Full", &n1, &n2, &b[i__ * b_dim1 + 1], ldb, &work[n1 * n2 + 
		1], &n1, (ftnlen)4);
#line 630 "ztgsen.f"
	ijb = 0;
#line 631 "ztgsen.f"
	i__1 = *lwork - (n1 << 1) * n2;
#line 631 "ztgsen.f"
	ztgsyl_("N", &ijb, &n1, &n2, &a[a_offset], lda, &a[i__ + i__ * a_dim1]
		, lda, &work[1], &n1, &b[b_offset], ldb, &b[i__ + i__ * 
		b_dim1], ldb, &work[n1 * n2 + 1], &n1, &dscale, &dif[1], &
		work[(n1 * n2 << 1) + 1], &i__1, &iwork[1], &ierr, (ftnlen)1);

/*        Estimate the reciprocal of norms of "projections" onto */
/*        left and right eigenspaces */

#line 639 "ztgsen.f"
	rdscal = 0.;
#line 640 "ztgsen.f"
	dsum = 1.;
#line 641 "ztgsen.f"
	i__1 = n1 * n2;
#line 641 "ztgsen.f"
	zlassq_(&i__1, &work[1], &c__1, &rdscal, &dsum);
#line 642 "ztgsen.f"
	*pl = rdscal * sqrt(dsum);
#line 643 "ztgsen.f"
	if (*pl == 0.) {
#line 644 "ztgsen.f"
	    *pl = 1.;
#line 645 "ztgsen.f"
	} else {
#line 646 "ztgsen.f"
	    *pl = dscale / (sqrt(dscale * dscale / *pl + *pl) * sqrt(*pl));
#line 647 "ztgsen.f"
	}
#line 648 "ztgsen.f"
	rdscal = 0.;
#line 649 "ztgsen.f"
	dsum = 1.;
#line 650 "ztgsen.f"
	i__1 = n1 * n2;
#line 650 "ztgsen.f"
	zlassq_(&i__1, &work[n1 * n2 + 1], &c__1, &rdscal, &dsum);
#line 651 "ztgsen.f"
	*pr = rdscal * sqrt(dsum);
#line 652 "ztgsen.f"
	if (*pr == 0.) {
#line 653 "ztgsen.f"
	    *pr = 1.;
#line 654 "ztgsen.f"
	} else {
#line 655 "ztgsen.f"
	    *pr = dscale / (sqrt(dscale * dscale / *pr + *pr) * sqrt(*pr));
#line 656 "ztgsen.f"
	}
#line 657 "ztgsen.f"
    }
#line 658 "ztgsen.f"
    if (wantd) {

/*        Compute estimates Difu and Difl. */

#line 662 "ztgsen.f"
	if (wantd1) {
#line 663 "ztgsen.f"
	    n1 = *m;
#line 664 "ztgsen.f"
	    n2 = *n - *m;
#line 665 "ztgsen.f"
	    i__ = n1 + 1;
#line 666 "ztgsen.f"
	    ijb = 3;

/*           Frobenius norm-based Difu estimate. */

#line 670 "ztgsen.f"
	    i__1 = *lwork - (n1 << 1) * n2;
#line 670 "ztgsen.f"
	    ztgsyl_("N", &ijb, &n1, &n2, &a[a_offset], lda, &a[i__ + i__ * 
		    a_dim1], lda, &work[1], &n1, &b[b_offset], ldb, &b[i__ + 
		    i__ * b_dim1], ldb, &work[n1 * n2 + 1], &n1, &dscale, &
		    dif[1], &work[(n1 * n2 << 1) + 1], &i__1, &iwork[1], &
		    ierr, (ftnlen)1);

/*           Frobenius norm-based Difl estimate. */

#line 677 "ztgsen.f"
	    i__1 = *lwork - (n1 << 1) * n2;
#line 677 "ztgsen.f"
	    ztgsyl_("N", &ijb, &n2, &n1, &a[i__ + i__ * a_dim1], lda, &a[
		    a_offset], lda, &work[1], &n2, &b[i__ + i__ * b_dim1], 
		    ldb, &b[b_offset], ldb, &work[n1 * n2 + 1], &n2, &dscale, 
		    &dif[2], &work[(n1 * n2 << 1) + 1], &i__1, &iwork[1], &
		    ierr, (ftnlen)1);
#line 681 "ztgsen.f"
	} else {

/*           Compute 1-norm-based estimates of Difu and Difl using */
/*           reversed communication with ZLACN2. In each step a */
/*           generalized Sylvester equation or a transposed variant */
/*           is solved. */

#line 688 "ztgsen.f"
	    kase = 0;
#line 689 "ztgsen.f"
	    n1 = *m;
#line 690 "ztgsen.f"
	    n2 = *n - *m;
#line 691 "ztgsen.f"
	    i__ = n1 + 1;
#line 692 "ztgsen.f"
	    ijb = 0;
#line 693 "ztgsen.f"
	    mn2 = (n1 << 1) * n2;

/*           1-norm-based estimate of Difu. */

#line 697 "ztgsen.f"
L40:
#line 698 "ztgsen.f"
	    zlacn2_(&mn2, &work[mn2 + 1], &work[1], &dif[1], &kase, isave);
#line 700 "ztgsen.f"
	    if (kase != 0) {
#line 701 "ztgsen.f"
		if (kase == 1) {

/*                 Solve generalized Sylvester equation */

#line 705 "ztgsen.f"
		    i__1 = *lwork - (n1 << 1) * n2;
#line 705 "ztgsen.f"
		    ztgsyl_("N", &ijb, &n1, &n2, &a[a_offset], lda, &a[i__ + 
			    i__ * a_dim1], lda, &work[1], &n1, &b[b_offset], 
			    ldb, &b[i__ + i__ * b_dim1], ldb, &work[n1 * n2 + 
			    1], &n1, &dscale, &dif[1], &work[(n1 * n2 << 1) + 
			    1], &i__1, &iwork[1], &ierr, (ftnlen)1);
#line 710 "ztgsen.f"
		} else {

/*                 Solve the transposed variant. */

#line 714 "ztgsen.f"
		    i__1 = *lwork - (n1 << 1) * n2;
#line 714 "ztgsen.f"
		    ztgsyl_("C", &ijb, &n1, &n2, &a[a_offset], lda, &a[i__ + 
			    i__ * a_dim1], lda, &work[1], &n1, &b[b_offset], 
			    ldb, &b[i__ + i__ * b_dim1], ldb, &work[n1 * n2 + 
			    1], &n1, &dscale, &dif[1], &work[(n1 * n2 << 1) + 
			    1], &i__1, &iwork[1], &ierr, (ftnlen)1);
#line 719 "ztgsen.f"
		}
#line 720 "ztgsen.f"
		goto L40;
#line 721 "ztgsen.f"
	    }
#line 722 "ztgsen.f"
	    dif[1] = dscale / dif[1];

/*           1-norm-based estimate of Difl. */

#line 726 "ztgsen.f"
L50:
#line 727 "ztgsen.f"
	    zlacn2_(&mn2, &work[mn2 + 1], &work[1], &dif[2], &kase, isave);
#line 729 "ztgsen.f"
	    if (kase != 0) {
#line 730 "ztgsen.f"
		if (kase == 1) {

/*                 Solve generalized Sylvester equation */

#line 734 "ztgsen.f"
		    i__1 = *lwork - (n1 << 1) * n2;
#line 734 "ztgsen.f"
		    ztgsyl_("N", &ijb, &n2, &n1, &a[i__ + i__ * a_dim1], lda, 
			    &a[a_offset], lda, &work[1], &n2, &b[i__ + i__ * 
			    b_dim1], ldb, &b[b_offset], ldb, &work[n1 * n2 + 
			    1], &n2, &dscale, &dif[2], &work[(n1 * n2 << 1) + 
			    1], &i__1, &iwork[1], &ierr, (ftnlen)1);
#line 739 "ztgsen.f"
		} else {

/*                 Solve the transposed variant. */

#line 743 "ztgsen.f"
		    i__1 = *lwork - (n1 << 1) * n2;
#line 743 "ztgsen.f"
		    ztgsyl_("C", &ijb, &n2, &n1, &a[i__ + i__ * a_dim1], lda, 
			    &a[a_offset], lda, &work[1], &n2, &b[b_offset], 
			    ldb, &b[i__ + i__ * b_dim1], ldb, &work[n1 * n2 + 
			    1], &n2, &dscale, &dif[2], &work[(n1 * n2 << 1) + 
			    1], &i__1, &iwork[1], &ierr, (ftnlen)1);
#line 748 "ztgsen.f"
		}
#line 749 "ztgsen.f"
		goto L50;
#line 750 "ztgsen.f"
	    }
#line 751 "ztgsen.f"
	    dif[2] = dscale / dif[2];
#line 752 "ztgsen.f"
	}
#line 753 "ztgsen.f"
    }

/*     If B(K,K) is complex, make it real and positive (normalization */
/*     of the generalized Schur form) and Store the generalized */
/*     eigenvalues of reordered pair (A, B) */

#line 759 "ztgsen.f"
    i__1 = *n;
#line 759 "ztgsen.f"
    for (k = 1; k <= i__1; ++k) {
#line 760 "ztgsen.f"
	dscale = z_abs(&b[k + k * b_dim1]);
#line 761 "ztgsen.f"
	if (dscale > safmin) {
#line 762 "ztgsen.f"
	    i__2 = k + k * b_dim1;
#line 762 "ztgsen.f"
	    z__2.r = b[i__2].r / dscale, z__2.i = b[i__2].i / dscale;
#line 762 "ztgsen.f"
	    d_cnjg(&z__1, &z__2);
#line 762 "ztgsen.f"
	    temp1.r = z__1.r, temp1.i = z__1.i;
#line 763 "ztgsen.f"
	    i__2 = k + k * b_dim1;
#line 763 "ztgsen.f"
	    z__1.r = b[i__2].r / dscale, z__1.i = b[i__2].i / dscale;
#line 763 "ztgsen.f"
	    temp2.r = z__1.r, temp2.i = z__1.i;
#line 764 "ztgsen.f"
	    i__2 = k + k * b_dim1;
#line 764 "ztgsen.f"
	    b[i__2].r = dscale, b[i__2].i = 0.;
#line 765 "ztgsen.f"
	    i__2 = *n - k;
#line 765 "ztgsen.f"
	    zscal_(&i__2, &temp1, &b[k + (k + 1) * b_dim1], ldb);
#line 766 "ztgsen.f"
	    i__2 = *n - k + 1;
#line 766 "ztgsen.f"
	    zscal_(&i__2, &temp1, &a[k + k * a_dim1], lda);
#line 767 "ztgsen.f"
	    if (*wantq) {
#line 767 "ztgsen.f"
		zscal_(n, &temp2, &q[k * q_dim1 + 1], &c__1);
#line 767 "ztgsen.f"
	    }
#line 769 "ztgsen.f"
	} else {
#line 770 "ztgsen.f"
	    i__2 = k + k * b_dim1;
#line 770 "ztgsen.f"
	    b[i__2].r = 0., b[i__2].i = 0.;
#line 771 "ztgsen.f"
	}

#line 773 "ztgsen.f"
	i__2 = k;
#line 773 "ztgsen.f"
	i__3 = k + k * a_dim1;
#line 773 "ztgsen.f"
	alpha[i__2].r = a[i__3].r, alpha[i__2].i = a[i__3].i;
#line 774 "ztgsen.f"
	i__2 = k;
#line 774 "ztgsen.f"
	i__3 = k + k * b_dim1;
#line 774 "ztgsen.f"
	beta[i__2].r = b[i__3].r, beta[i__2].i = b[i__3].i;

#line 776 "ztgsen.f"
/* L60: */
#line 776 "ztgsen.f"
    }

#line 778 "ztgsen.f"
L70:

#line 780 "ztgsen.f"
    work[1].r = (doublereal) lwmin, work[1].i = 0.;
#line 781 "ztgsen.f"
    iwork[1] = liwmin;

#line 783 "ztgsen.f"
    return 0;

/*     End of ZTGSEN */

} /* ztgsen_ */

