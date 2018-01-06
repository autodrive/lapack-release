#line 1 "ctgsja.f"
/* ctgsja.f -- translated by f2c (version 20100827).
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

#line 1 "ctgsja.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};
static integer c__1 = 1;
static doublereal c_b39 = -1.;
static doublereal c_b42 = 1.;

/* > \brief \b CTGSJA */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CTGSJA + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ctgsja.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ctgsja.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ctgsja.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CTGSJA( JOBU, JOBV, JOBQ, M, P, N, K, L, A, LDA, B, */
/*                          LDB, TOLA, TOLB, ALPHA, BETA, U, LDU, V, LDV, */
/*                          Q, LDQ, WORK, NCYCLE, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBQ, JOBU, JOBV */
/*       INTEGER            INFO, K, L, LDA, LDB, LDQ, LDU, LDV, M, N, */
/*      $                   NCYCLE, P */
/*       REAL               TOLA, TOLB */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               ALPHA( * ), BETA( * ) */
/*       COMPLEX            A( LDA, * ), B( LDB, * ), Q( LDQ, * ), */
/*      $                   U( LDU, * ), V( LDV, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CTGSJA computes the generalized singular value decomposition (GSVD) */
/* > of two complex upper triangular (or trapezoidal) matrices A and B. */
/* > */
/* > On entry, it is assumed that matrices A and B have the following */
/* > forms, which may be obtained by the preprocessing subroutine CGGSVP */
/* > from a general M-by-N matrix A and P-by-N matrix B: */
/* > */
/* >              N-K-L  K    L */
/* >    A =    K ( 0    A12  A13 ) if M-K-L >= 0; */
/* >           L ( 0     0   A23 ) */
/* >       M-K-L ( 0     0    0  ) */
/* > */
/* >            N-K-L  K    L */
/* >    A =  K ( 0    A12  A13 ) if M-K-L < 0; */
/* >       M-K ( 0     0   A23 ) */
/* > */
/* >            N-K-L  K    L */
/* >    B =  L ( 0     0   B13 ) */
/* >       P-L ( 0     0    0  ) */
/* > */
/* > where the K-by-K matrix A12 and L-by-L matrix B13 are nonsingular */
/* > upper triangular; A23 is L-by-L upper triangular if M-K-L >= 0, */
/* > otherwise A23 is (M-K)-by-L upper trapezoidal. */
/* > */
/* > On exit, */
/* > */
/* >        U**H *A*Q = D1*( 0 R ),    V**H *B*Q = D2*( 0 R ), */
/* > */
/* > where U, V and Q are unitary matrices. */
/* > R is a nonsingular upper triangular matrix, and D1 */
/* > and D2 are ``diagonal'' matrices, which are of the following */
/* > structures: */
/* > */
/* > If M-K-L >= 0, */
/* > */
/* >                     K  L */
/* >        D1 =     K ( I  0 ) */
/* >                 L ( 0  C ) */
/* >             M-K-L ( 0  0 ) */
/* > */
/* >                    K  L */
/* >        D2 = L   ( 0  S ) */
/* >             P-L ( 0  0 ) */
/* > */
/* >                N-K-L  K    L */
/* >   ( 0 R ) = K (  0   R11  R12 ) K */
/* >             L (  0    0   R22 ) L */
/* > */
/* > where */
/* > */
/* >   C = diag( ALPHA(K+1), ... , ALPHA(K+L) ), */
/* >   S = diag( BETA(K+1),  ... , BETA(K+L) ), */
/* >   C**2 + S**2 = I. */
/* > */
/* >   R is stored in A(1:K+L,N-K-L+1:N) on exit. */
/* > */
/* > If M-K-L < 0, */
/* > */
/* >                K M-K K+L-M */
/* >     D1 =   K ( I  0    0   ) */
/* >          M-K ( 0  C    0   ) */
/* > */
/* >                  K M-K K+L-M */
/* >     D2 =   M-K ( 0  S    0   ) */
/* >          K+L-M ( 0  0    I   ) */
/* >            P-L ( 0  0    0   ) */
/* > */
/* >                N-K-L  K   M-K  K+L-M */
/* > ( 0 R ) =    K ( 0    R11  R12  R13  ) */
/* >           M-K ( 0     0   R22  R23  ) */
/* >         K+L-M ( 0     0    0   R33  ) */
/* > */
/* > where */
/* > C = diag( ALPHA(K+1), ... , ALPHA(M) ), */
/* > S = diag( BETA(K+1),  ... , BETA(M) ), */
/* > C**2 + S**2 = I. */
/* > */
/* > R = ( R11 R12 R13 ) is stored in A(1:M, N-K-L+1:N) and R33 is stored */
/* >     (  0  R22 R23 ) */
/* > in B(M-K+1:L,N+M-K-L+1:N) on exit. */
/* > */
/* > The computation of the unitary transformation matrices U, V or Q */
/* > is optional.  These matrices may either be formed explicitly, or they */
/* > may be postmultiplied into input matrices U1, V1, or Q1. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] JOBU */
/* > \verbatim */
/* >          JOBU is CHARACTER*1 */
/* >          = 'U':  U must contain a unitary matrix U1 on entry, and */
/* >                  the product U1*U is returned; */
/* >          = 'I':  U is initialized to the unit matrix, and the */
/* >                  unitary matrix U is returned; */
/* >          = 'N':  U is not computed. */
/* > \endverbatim */
/* > */
/* > \param[in] JOBV */
/* > \verbatim */
/* >          JOBV is CHARACTER*1 */
/* >          = 'V':  V must contain a unitary matrix V1 on entry, and */
/* >                  the product V1*V is returned; */
/* >          = 'I':  V is initialized to the unit matrix, and the */
/* >                  unitary matrix V is returned; */
/* >          = 'N':  V is not computed. */
/* > \endverbatim */
/* > */
/* > \param[in] JOBQ */
/* > \verbatim */
/* >          JOBQ is CHARACTER*1 */
/* >          = 'Q':  Q must contain a unitary matrix Q1 on entry, and */
/* >                  the product Q1*Q is returned; */
/* >          = 'I':  Q is initialized to the unit matrix, and the */
/* >                  unitary matrix Q is returned; */
/* >          = 'N':  Q is not computed. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the matrix A.  M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] P */
/* > \verbatim */
/* >          P is INTEGER */
/* >          The number of rows of the matrix B.  P >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrices A and B.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] K */
/* > \verbatim */
/* >          K is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[in] L */
/* > \verbatim */
/* >          L is INTEGER */
/* > */
/* >          K and L specify the subblocks in the input matrices A and B: */
/* >          A23 = A(K+1:MIN(K+L,M),N-L+1:N) and B13 = B(1:L,,N-L+1:N) */
/* >          of A and B, whose GSVD is going to be computed by CTGSJA. */
/* >          See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX array, dimension (LDA,N) */
/* >          On entry, the M-by-N matrix A. */
/* >          On exit, A(N-K+1:N,1:MIN(K+L,M) ) contains the triangular */
/* >          matrix R or part of R.  See Purpose for details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A. LDA >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* >          B is COMPLEX array, dimension (LDB,N) */
/* >          On entry, the P-by-N matrix B. */
/* >          On exit, if necessary, B(M-K+1:L,N+M-K-L+1:N) contains */
/* >          a part of R.  See Purpose for details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >          The leading dimension of the array B. LDB >= max(1,P). */
/* > \endverbatim */
/* > */
/* > \param[in] TOLA */
/* > \verbatim */
/* >          TOLA is REAL */
/* > \endverbatim */
/* > */
/* > \param[in] TOLB */
/* > \verbatim */
/* >          TOLB is REAL */
/* > */
/* >          TOLA and TOLB are the convergence criteria for the Jacobi- */
/* >          Kogbetliantz iteration procedure. Generally, they are the */
/* >          same as used in the preprocessing step, say */
/* >              TOLA = MAX(M,N)*norm(A)*MACHEPS, */
/* >              TOLB = MAX(P,N)*norm(B)*MACHEPS. */
/* > \endverbatim */
/* > */
/* > \param[out] ALPHA */
/* > \verbatim */
/* >          ALPHA is REAL array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] BETA */
/* > \verbatim */
/* >          BETA is REAL array, dimension (N) */
/* > */
/* >          On exit, ALPHA and BETA contain the generalized singular */
/* >          value pairs of A and B; */
/* >            ALPHA(1:K) = 1, */
/* >            BETA(1:K)  = 0, */
/* >          and if M-K-L >= 0, */
/* >            ALPHA(K+1:K+L) = diag(C), */
/* >            BETA(K+1:K+L)  = diag(S), */
/* >          or if M-K-L < 0, */
/* >            ALPHA(K+1:M)= C, ALPHA(M+1:K+L)= 0 */
/* >            BETA(K+1:M) = S, BETA(M+1:K+L) = 1. */
/* >          Furthermore, if K+L < N, */
/* >            ALPHA(K+L+1:N) = 0 */
/* >            BETA(K+L+1:N)  = 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] U */
/* > \verbatim */
/* >          U is COMPLEX array, dimension (LDU,M) */
/* >          On entry, if JOBU = 'U', U must contain a matrix U1 (usually */
/* >          the unitary matrix returned by CGGSVP). */
/* >          On exit, */
/* >          if JOBU = 'I', U contains the unitary matrix U; */
/* >          if JOBU = 'U', U contains the product U1*U. */
/* >          If JOBU = 'N', U is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDU */
/* > \verbatim */
/* >          LDU is INTEGER */
/* >          The leading dimension of the array U. LDU >= max(1,M) if */
/* >          JOBU = 'U'; LDU >= 1 otherwise. */
/* > \endverbatim */
/* > */
/* > \param[in,out] V */
/* > \verbatim */
/* >          V is COMPLEX array, dimension (LDV,P) */
/* >          On entry, if JOBV = 'V', V must contain a matrix V1 (usually */
/* >          the unitary matrix returned by CGGSVP). */
/* >          On exit, */
/* >          if JOBV = 'I', V contains the unitary matrix V; */
/* >          if JOBV = 'V', V contains the product V1*V. */
/* >          If JOBV = 'N', V is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDV */
/* > \verbatim */
/* >          LDV is INTEGER */
/* >          The leading dimension of the array V. LDV >= max(1,P) if */
/* >          JOBV = 'V'; LDV >= 1 otherwise. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Q */
/* > \verbatim */
/* >          Q is COMPLEX array, dimension (LDQ,N) */
/* >          On entry, if JOBQ = 'Q', Q must contain a matrix Q1 (usually */
/* >          the unitary matrix returned by CGGSVP). */
/* >          On exit, */
/* >          if JOBQ = 'I', Q contains the unitary matrix Q; */
/* >          if JOBQ = 'Q', Q contains the product Q1*Q. */
/* >          If JOBQ = 'N', Q is not referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDQ */
/* > \verbatim */
/* >          LDQ is INTEGER */
/* >          The leading dimension of the array Q. LDQ >= max(1,N) if */
/* >          JOBQ = 'Q'; LDQ >= 1 otherwise. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX array, dimension (2*N) */
/* > \endverbatim */
/* > */
/* > \param[out] NCYCLE */
/* > \verbatim */
/* >          NCYCLE is INTEGER */
/* >          The number of cycles required for convergence. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value. */
/* >          = 1:  the procedure does not converge after MAXIT cycles. */
/* > \endverbatim */

/* > \par Internal Parameters: */
/*  ========================= */
/* > */
/* > \verbatim */
/* >  MAXIT   INTEGER */
/* >          MAXIT specifies the total loops that the iterative procedure */
/* >          may take. If after MAXIT cycles, the routine fails to */
/* >          converge, we return INFO = 1. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup complexOTHERcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  CTGSJA essentially uses a variant of Kogbetliantz algorithm to reduce */
/* >  min(L,M-K)-by-L triangular (or trapezoidal) matrix A23 and L-by-L */
/* >  matrix B13 to the form: */
/* > */
/* >           U1**H *A13*Q1 = C1*R1; V1**H *B13*Q1 = S1*R1, */
/* > */
/* >  where U1, V1 and Q1 are unitary matrix. */
/* >  C1 and S1 are diagonal matrices satisfying */
/* > */
/* >                C1**2 + S1**2 = I, */
/* > */
/* >  and R1 is an L-by-L nonsingular upper triangular matrix. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int ctgsja_(char *jobu, char *jobv, char *jobq, integer *m, 
	integer *p, integer *n, integer *k, integer *l, doublecomplex *a, 
	integer *lda, doublecomplex *b, integer *ldb, doublereal *tola, 
	doublereal *tolb, doublereal *alpha, doublereal *beta, doublecomplex *
	u, integer *ldu, doublecomplex *v, integer *ldv, doublecomplex *q, 
	integer *ldq, doublecomplex *work, integer *ncycle, integer *info, 
	ftnlen jobu_len, ftnlen jobv_len, ftnlen jobq_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, q_dim1, q_offset, u_dim1, 
	    u_offset, v_dim1, v_offset, i__1, i__2, i__3, i__4;
    doublereal d__1;
    doublecomplex z__1;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, j;
    static doublereal a1, b1, a3, b3;
    static doublecomplex a2, b2;
    static doublereal csq, csu, csv;
    static doublecomplex snq;
    static doublereal rwk;
    static doublecomplex snu, snv;
    extern /* Subroutine */ int crot_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublereal *, doublecomplex *);
    static doublereal gamma;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int ccopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    static logical initq, initu, initv, wantq, upper;
    static doublereal error, ssmin;
    static logical wantu, wantv;
    extern /* Subroutine */ int clags2_(logical *, doublereal *, 
	    doublecomplex *, doublereal *, doublereal *, doublecomplex *, 
	    doublereal *, doublereal *, doublecomplex *, doublereal *, 
	    doublecomplex *, doublereal *, doublecomplex *), clapll_(integer *
	    , doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublereal *), csscal_(integer *, doublereal *, doublecomplex *, 
	    integer *);
    static integer kcycle;
    extern /* Subroutine */ int claset_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, doublecomplex *, integer *, 
	    ftnlen), xerbla_(char *, integer *, ftnlen), slartg_(doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *);


/*  -- LAPACK computational routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

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

/*     Decode and test the input parameters */

#line 433 "ctgsja.f"
    /* Parameter adjustments */
#line 433 "ctgsja.f"
    a_dim1 = *lda;
#line 433 "ctgsja.f"
    a_offset = 1 + a_dim1;
#line 433 "ctgsja.f"
    a -= a_offset;
#line 433 "ctgsja.f"
    b_dim1 = *ldb;
#line 433 "ctgsja.f"
    b_offset = 1 + b_dim1;
#line 433 "ctgsja.f"
    b -= b_offset;
#line 433 "ctgsja.f"
    --alpha;
#line 433 "ctgsja.f"
    --beta;
#line 433 "ctgsja.f"
    u_dim1 = *ldu;
#line 433 "ctgsja.f"
    u_offset = 1 + u_dim1;
#line 433 "ctgsja.f"
    u -= u_offset;
#line 433 "ctgsja.f"
    v_dim1 = *ldv;
#line 433 "ctgsja.f"
    v_offset = 1 + v_dim1;
#line 433 "ctgsja.f"
    v -= v_offset;
#line 433 "ctgsja.f"
    q_dim1 = *ldq;
#line 433 "ctgsja.f"
    q_offset = 1 + q_dim1;
#line 433 "ctgsja.f"
    q -= q_offset;
#line 433 "ctgsja.f"
    --work;
#line 433 "ctgsja.f"

#line 433 "ctgsja.f"
    /* Function Body */
#line 433 "ctgsja.f"
    initu = lsame_(jobu, "I", (ftnlen)1, (ftnlen)1);
#line 434 "ctgsja.f"
    wantu = initu || lsame_(jobu, "U", (ftnlen)1, (ftnlen)1);

#line 436 "ctgsja.f"
    initv = lsame_(jobv, "I", (ftnlen)1, (ftnlen)1);
#line 437 "ctgsja.f"
    wantv = initv || lsame_(jobv, "V", (ftnlen)1, (ftnlen)1);

#line 439 "ctgsja.f"
    initq = lsame_(jobq, "I", (ftnlen)1, (ftnlen)1);
#line 440 "ctgsja.f"
    wantq = initq || lsame_(jobq, "Q", (ftnlen)1, (ftnlen)1);

#line 442 "ctgsja.f"
    *info = 0;
#line 443 "ctgsja.f"
    if (! (initu || wantu || lsame_(jobu, "N", (ftnlen)1, (ftnlen)1))) {
#line 444 "ctgsja.f"
	*info = -1;
#line 445 "ctgsja.f"
    } else if (! (initv || wantv || lsame_(jobv, "N", (ftnlen)1, (ftnlen)1))) 
	    {
#line 446 "ctgsja.f"
	*info = -2;
#line 447 "ctgsja.f"
    } else if (! (initq || wantq || lsame_(jobq, "N", (ftnlen)1, (ftnlen)1))) 
	    {
#line 448 "ctgsja.f"
	*info = -3;
#line 449 "ctgsja.f"
    } else if (*m < 0) {
#line 450 "ctgsja.f"
	*info = -4;
#line 451 "ctgsja.f"
    } else if (*p < 0) {
#line 452 "ctgsja.f"
	*info = -5;
#line 453 "ctgsja.f"
    } else if (*n < 0) {
#line 454 "ctgsja.f"
	*info = -6;
#line 455 "ctgsja.f"
    } else if (*lda < max(1,*m)) {
#line 456 "ctgsja.f"
	*info = -10;
#line 457 "ctgsja.f"
    } else if (*ldb < max(1,*p)) {
#line 458 "ctgsja.f"
	*info = -12;
#line 459 "ctgsja.f"
    } else if (*ldu < 1 || wantu && *ldu < *m) {
#line 460 "ctgsja.f"
	*info = -18;
#line 461 "ctgsja.f"
    } else if (*ldv < 1 || wantv && *ldv < *p) {
#line 462 "ctgsja.f"
	*info = -20;
#line 463 "ctgsja.f"
    } else if (*ldq < 1 || wantq && *ldq < *n) {
#line 464 "ctgsja.f"
	*info = -22;
#line 465 "ctgsja.f"
    }
#line 466 "ctgsja.f"
    if (*info != 0) {
#line 467 "ctgsja.f"
	i__1 = -(*info);
#line 467 "ctgsja.f"
	xerbla_("CTGSJA", &i__1, (ftnlen)6);
#line 468 "ctgsja.f"
	return 0;
#line 469 "ctgsja.f"
    }

/*     Initialize U, V and Q, if necessary */

#line 473 "ctgsja.f"
    if (initu) {
#line 473 "ctgsja.f"
	claset_("Full", m, m, &c_b1, &c_b2, &u[u_offset], ldu, (ftnlen)4);
#line 473 "ctgsja.f"
    }
#line 475 "ctgsja.f"
    if (initv) {
#line 475 "ctgsja.f"
	claset_("Full", p, p, &c_b1, &c_b2, &v[v_offset], ldv, (ftnlen)4);
#line 475 "ctgsja.f"
    }
#line 477 "ctgsja.f"
    if (initq) {
#line 477 "ctgsja.f"
	claset_("Full", n, n, &c_b1, &c_b2, &q[q_offset], ldq, (ftnlen)4);
#line 477 "ctgsja.f"
    }

/*     Loop until convergence */

#line 482 "ctgsja.f"
    upper = FALSE_;
#line 483 "ctgsja.f"
    for (kcycle = 1; kcycle <= 40; ++kcycle) {

#line 485 "ctgsja.f"
	upper = ! upper;

#line 487 "ctgsja.f"
	i__1 = *l - 1;
#line 487 "ctgsja.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 488 "ctgsja.f"
	    i__2 = *l;
#line 488 "ctgsja.f"
	    for (j = i__ + 1; j <= i__2; ++j) {

#line 490 "ctgsja.f"
		a1 = 0.;
#line 491 "ctgsja.f"
		a2.r = 0., a2.i = 0.;
#line 492 "ctgsja.f"
		a3 = 0.;
#line 493 "ctgsja.f"
		if (*k + i__ <= *m) {
#line 493 "ctgsja.f"
		    i__3 = *k + i__ + (*n - *l + i__) * a_dim1;
#line 493 "ctgsja.f"
		    a1 = a[i__3].r;
#line 493 "ctgsja.f"
		}
#line 495 "ctgsja.f"
		if (*k + j <= *m) {
#line 495 "ctgsja.f"
		    i__3 = *k + j + (*n - *l + j) * a_dim1;
#line 495 "ctgsja.f"
		    a3 = a[i__3].r;
#line 495 "ctgsja.f"
		}

#line 498 "ctgsja.f"
		i__3 = i__ + (*n - *l + i__) * b_dim1;
#line 498 "ctgsja.f"
		b1 = b[i__3].r;
#line 499 "ctgsja.f"
		i__3 = j + (*n - *l + j) * b_dim1;
#line 499 "ctgsja.f"
		b3 = b[i__3].r;

#line 501 "ctgsja.f"
		if (upper) {
#line 502 "ctgsja.f"
		    if (*k + i__ <= *m) {
#line 502 "ctgsja.f"
			i__3 = *k + i__ + (*n - *l + j) * a_dim1;
#line 502 "ctgsja.f"
			a2.r = a[i__3].r, a2.i = a[i__3].i;
#line 502 "ctgsja.f"
		    }
#line 504 "ctgsja.f"
		    i__3 = i__ + (*n - *l + j) * b_dim1;
#line 504 "ctgsja.f"
		    b2.r = b[i__3].r, b2.i = b[i__3].i;
#line 505 "ctgsja.f"
		} else {
#line 506 "ctgsja.f"
		    if (*k + j <= *m) {
#line 506 "ctgsja.f"
			i__3 = *k + j + (*n - *l + i__) * a_dim1;
#line 506 "ctgsja.f"
			a2.r = a[i__3].r, a2.i = a[i__3].i;
#line 506 "ctgsja.f"
		    }
#line 508 "ctgsja.f"
		    i__3 = j + (*n - *l + i__) * b_dim1;
#line 508 "ctgsja.f"
		    b2.r = b[i__3].r, b2.i = b[i__3].i;
#line 509 "ctgsja.f"
		}

#line 511 "ctgsja.f"
		clags2_(&upper, &a1, &a2, &a3, &b1, &b2, &b3, &csu, &snu, &
			csv, &snv, &csq, &snq);

/*              Update (K+I)-th and (K+J)-th rows of matrix A: U**H *A */

#line 516 "ctgsja.f"
		if (*k + j <= *m) {
#line 516 "ctgsja.f"
		    d_cnjg(&z__1, &snu);
#line 516 "ctgsja.f"
		    crot_(l, &a[*k + j + (*n - *l + 1) * a_dim1], lda, &a[*k 
			    + i__ + (*n - *l + 1) * a_dim1], lda, &csu, &z__1)
			    ;
#line 516 "ctgsja.f"
		}

/*              Update I-th and J-th rows of matrix B: V**H *B */

#line 522 "ctgsja.f"
		d_cnjg(&z__1, &snv);
#line 522 "ctgsja.f"
		crot_(l, &b[j + (*n - *l + 1) * b_dim1], ldb, &b[i__ + (*n - *
			l + 1) * b_dim1], ldb, &csv, &z__1);

/*              Update (N-L+I)-th and (N-L+J)-th columns of matrices */
/*              A and B: A*Q and B*Q */

/* Computing MIN */
#line 528 "ctgsja.f"
		i__4 = *k + *l;
#line 528 "ctgsja.f"
		i__3 = min(i__4,*m);
#line 528 "ctgsja.f"
		crot_(&i__3, &a[(*n - *l + j) * a_dim1 + 1], &c__1, &a[(*n - *
			l + i__) * a_dim1 + 1], &c__1, &csq, &snq);

#line 531 "ctgsja.f"
		crot_(l, &b[(*n - *l + j) * b_dim1 + 1], &c__1, &b[(*n - *l + 
			i__) * b_dim1 + 1], &c__1, &csq, &snq);

#line 534 "ctgsja.f"
		if (upper) {
#line 535 "ctgsja.f"
		    if (*k + i__ <= *m) {
#line 535 "ctgsja.f"
			i__3 = *k + i__ + (*n - *l + j) * a_dim1;
#line 535 "ctgsja.f"
			a[i__3].r = 0., a[i__3].i = 0.;
#line 535 "ctgsja.f"
		    }
#line 537 "ctgsja.f"
		    i__3 = i__ + (*n - *l + j) * b_dim1;
#line 537 "ctgsja.f"
		    b[i__3].r = 0., b[i__3].i = 0.;
#line 538 "ctgsja.f"
		} else {
#line 539 "ctgsja.f"
		    if (*k + j <= *m) {
#line 539 "ctgsja.f"
			i__3 = *k + j + (*n - *l + i__) * a_dim1;
#line 539 "ctgsja.f"
			a[i__3].r = 0., a[i__3].i = 0.;
#line 539 "ctgsja.f"
		    }
#line 541 "ctgsja.f"
		    i__3 = j + (*n - *l + i__) * b_dim1;
#line 541 "ctgsja.f"
		    b[i__3].r = 0., b[i__3].i = 0.;
#line 542 "ctgsja.f"
		}

/*              Ensure that the diagonal elements of A and B are real. */

#line 546 "ctgsja.f"
		if (*k + i__ <= *m) {
#line 546 "ctgsja.f"
		    i__3 = *k + i__ + (*n - *l + i__) * a_dim1;
#line 546 "ctgsja.f"
		    i__4 = *k + i__ + (*n - *l + i__) * a_dim1;
#line 546 "ctgsja.f"
		    d__1 = a[i__4].r;
#line 546 "ctgsja.f"
		    a[i__3].r = d__1, a[i__3].i = 0.;
#line 546 "ctgsja.f"
		}
#line 548 "ctgsja.f"
		if (*k + j <= *m) {
#line 548 "ctgsja.f"
		    i__3 = *k + j + (*n - *l + j) * a_dim1;
#line 548 "ctgsja.f"
		    i__4 = *k + j + (*n - *l + j) * a_dim1;
#line 548 "ctgsja.f"
		    d__1 = a[i__4].r;
#line 548 "ctgsja.f"
		    a[i__3].r = d__1, a[i__3].i = 0.;
#line 548 "ctgsja.f"
		}
#line 550 "ctgsja.f"
		i__3 = i__ + (*n - *l + i__) * b_dim1;
#line 550 "ctgsja.f"
		i__4 = i__ + (*n - *l + i__) * b_dim1;
#line 550 "ctgsja.f"
		d__1 = b[i__4].r;
#line 550 "ctgsja.f"
		b[i__3].r = d__1, b[i__3].i = 0.;
#line 551 "ctgsja.f"
		i__3 = j + (*n - *l + j) * b_dim1;
#line 551 "ctgsja.f"
		i__4 = j + (*n - *l + j) * b_dim1;
#line 551 "ctgsja.f"
		d__1 = b[i__4].r;
#line 551 "ctgsja.f"
		b[i__3].r = d__1, b[i__3].i = 0.;

/*              Update unitary matrices U, V, Q, if desired. */

#line 555 "ctgsja.f"
		if (wantu && *k + j <= *m) {
#line 555 "ctgsja.f"
		    crot_(m, &u[(*k + j) * u_dim1 + 1], &c__1, &u[(*k + i__) *
			     u_dim1 + 1], &c__1, &csu, &snu);
#line 555 "ctgsja.f"
		}

#line 559 "ctgsja.f"
		if (wantv) {
#line 559 "ctgsja.f"
		    crot_(p, &v[j * v_dim1 + 1], &c__1, &v[i__ * v_dim1 + 1], 
			    &c__1, &csv, &snv);
#line 559 "ctgsja.f"
		}

#line 562 "ctgsja.f"
		if (wantq) {
#line 562 "ctgsja.f"
		    crot_(n, &q[(*n - *l + j) * q_dim1 + 1], &c__1, &q[(*n - *
			    l + i__) * q_dim1 + 1], &c__1, &csq, &snq);
#line 562 "ctgsja.f"
		}

#line 566 "ctgsja.f"
/* L10: */
#line 566 "ctgsja.f"
	    }
#line 567 "ctgsja.f"
/* L20: */
#line 567 "ctgsja.f"
	}

#line 569 "ctgsja.f"
	if (! upper) {

/*           The matrices A13 and B13 were lower triangular at the start */
/*           of the cycle, and are now upper triangular. */

/*           Convergence test: test the parallelism of the corresponding */
/*           rows of A and B. */

#line 577 "ctgsja.f"
	    error = 0.;
/* Computing MIN */
#line 578 "ctgsja.f"
	    i__2 = *l, i__3 = *m - *k;
#line 578 "ctgsja.f"
	    i__1 = min(i__2,i__3);
#line 578 "ctgsja.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 579 "ctgsja.f"
		i__2 = *l - i__ + 1;
#line 579 "ctgsja.f"
		ccopy_(&i__2, &a[*k + i__ + (*n - *l + i__) * a_dim1], lda, &
			work[1], &c__1);
#line 580 "ctgsja.f"
		i__2 = *l - i__ + 1;
#line 580 "ctgsja.f"
		ccopy_(&i__2, &b[i__ + (*n - *l + i__) * b_dim1], ldb, &work[*
			l + 1], &c__1);
#line 581 "ctgsja.f"
		i__2 = *l - i__ + 1;
#line 581 "ctgsja.f"
		clapll_(&i__2, &work[1], &c__1, &work[*l + 1], &c__1, &ssmin);
#line 582 "ctgsja.f"
		error = max(error,ssmin);
#line 583 "ctgsja.f"
/* L30: */
#line 583 "ctgsja.f"
	    }

#line 585 "ctgsja.f"
	    if (abs(error) <= min(*tola,*tolb)) {
#line 585 "ctgsja.f"
		goto L50;
#line 585 "ctgsja.f"
	    }
#line 587 "ctgsja.f"
	}

/*        End of cycle loop */

#line 591 "ctgsja.f"
/* L40: */
#line 591 "ctgsja.f"
    }

/*     The algorithm has not converged after MAXIT cycles. */

#line 595 "ctgsja.f"
    *info = 1;
#line 596 "ctgsja.f"
    goto L100;

#line 598 "ctgsja.f"
L50:

/*     If ERROR <= MIN(TOLA,TOLB), then the algorithm has converged. */
/*     Compute the generalized singular value pairs (ALPHA, BETA), and */
/*     set the triangular matrix R to array A. */

#line 604 "ctgsja.f"
    i__1 = *k;
#line 604 "ctgsja.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 605 "ctgsja.f"
	alpha[i__] = 1.;
#line 606 "ctgsja.f"
	beta[i__] = 0.;
#line 607 "ctgsja.f"
/* L60: */
#line 607 "ctgsja.f"
    }

/* Computing MIN */
#line 609 "ctgsja.f"
    i__2 = *l, i__3 = *m - *k;
#line 609 "ctgsja.f"
    i__1 = min(i__2,i__3);
#line 609 "ctgsja.f"
    for (i__ = 1; i__ <= i__1; ++i__) {

#line 611 "ctgsja.f"
	i__2 = *k + i__ + (*n - *l + i__) * a_dim1;
#line 611 "ctgsja.f"
	a1 = a[i__2].r;
#line 612 "ctgsja.f"
	i__2 = i__ + (*n - *l + i__) * b_dim1;
#line 612 "ctgsja.f"
	b1 = b[i__2].r;

#line 614 "ctgsja.f"
	if (a1 != 0.) {
#line 615 "ctgsja.f"
	    gamma = b1 / a1;

#line 617 "ctgsja.f"
	    if (gamma < 0.) {
#line 618 "ctgsja.f"
		i__2 = *l - i__ + 1;
#line 618 "ctgsja.f"
		csscal_(&i__2, &c_b39, &b[i__ + (*n - *l + i__) * b_dim1], 
			ldb);
#line 619 "ctgsja.f"
		if (wantv) {
#line 619 "ctgsja.f"
		    csscal_(p, &c_b39, &v[i__ * v_dim1 + 1], &c__1);
#line 619 "ctgsja.f"
		}
#line 621 "ctgsja.f"
	    }

#line 623 "ctgsja.f"
	    d__1 = abs(gamma);
#line 623 "ctgsja.f"
	    slartg_(&d__1, &c_b42, &beta[*k + i__], &alpha[*k + i__], &rwk);

#line 626 "ctgsja.f"
	    if (alpha[*k + i__] >= beta[*k + i__]) {
#line 627 "ctgsja.f"
		i__2 = *l - i__ + 1;
#line 627 "ctgsja.f"
		d__1 = 1. / alpha[*k + i__];
#line 627 "ctgsja.f"
		csscal_(&i__2, &d__1, &a[*k + i__ + (*n - *l + i__) * a_dim1],
			 lda);
#line 629 "ctgsja.f"
	    } else {
#line 630 "ctgsja.f"
		i__2 = *l - i__ + 1;
#line 630 "ctgsja.f"
		d__1 = 1. / beta[*k + i__];
#line 630 "ctgsja.f"
		csscal_(&i__2, &d__1, &b[i__ + (*n - *l + i__) * b_dim1], ldb)
			;
#line 632 "ctgsja.f"
		i__2 = *l - i__ + 1;
#line 632 "ctgsja.f"
		ccopy_(&i__2, &b[i__ + (*n - *l + i__) * b_dim1], ldb, &a[*k 
			+ i__ + (*n - *l + i__) * a_dim1], lda);
#line 634 "ctgsja.f"
	    }

#line 636 "ctgsja.f"
	} else {
#line 637 "ctgsja.f"
	    alpha[*k + i__] = 0.;
#line 638 "ctgsja.f"
	    beta[*k + i__] = 1.;
#line 639 "ctgsja.f"
	    i__2 = *l - i__ + 1;
#line 639 "ctgsja.f"
	    ccopy_(&i__2, &b[i__ + (*n - *l + i__) * b_dim1], ldb, &a[*k + 
		    i__ + (*n - *l + i__) * a_dim1], lda);
#line 641 "ctgsja.f"
	}
#line 642 "ctgsja.f"
/* L70: */
#line 642 "ctgsja.f"
    }

/*     Post-assignment */

#line 646 "ctgsja.f"
    i__1 = *k + *l;
#line 646 "ctgsja.f"
    for (i__ = *m + 1; i__ <= i__1; ++i__) {
#line 647 "ctgsja.f"
	alpha[i__] = 0.;
#line 648 "ctgsja.f"
	beta[i__] = 1.;
#line 649 "ctgsja.f"
/* L80: */
#line 649 "ctgsja.f"
    }

#line 651 "ctgsja.f"
    if (*k + *l < *n) {
#line 652 "ctgsja.f"
	i__1 = *n;
#line 652 "ctgsja.f"
	for (i__ = *k + *l + 1; i__ <= i__1; ++i__) {
#line 653 "ctgsja.f"
	    alpha[i__] = 0.;
#line 654 "ctgsja.f"
	    beta[i__] = 0.;
#line 655 "ctgsja.f"
/* L90: */
#line 655 "ctgsja.f"
	}
#line 656 "ctgsja.f"
    }

#line 658 "ctgsja.f"
L100:
#line 659 "ctgsja.f"
    *ncycle = kcycle;

#line 661 "ctgsja.f"
    return 0;

/*     End of CTGSJA */

} /* ctgsja_ */

