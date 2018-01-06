#line 1 "dtgsja.f"
/* dtgsja.f -- translated by f2c (version 20100827).
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

#line 1 "dtgsja.f"
/* Table of constant values */

static doublereal c_b13 = 0.;
static doublereal c_b14 = 1.;
static integer c__1 = 1;
static doublereal c_b43 = -1.;

/* > \brief \b DTGSJA */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DTGSJA + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dtgsja.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dtgsja.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dtgsja.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DTGSJA( JOBU, JOBV, JOBQ, M, P, N, K, L, A, LDA, B, */
/*                          LDB, TOLA, TOLB, ALPHA, BETA, U, LDU, V, LDV, */
/*                          Q, LDQ, WORK, NCYCLE, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBQ, JOBU, JOBV */
/*       INTEGER            INFO, K, L, LDA, LDB, LDQ, LDU, LDV, M, N, */
/*      $                   NCYCLE, P */
/*       DOUBLE PRECISION   TOLA, TOLB */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   A( LDA, * ), ALPHA( * ), B( LDB, * ), */
/*      $                   BETA( * ), Q( LDQ, * ), U( LDU, * ), */
/*      $                   V( LDV, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DTGSJA computes the generalized singular value decomposition (GSVD) */
/* > of two real upper triangular (or trapezoidal) matrices A and B. */
/* > */
/* > On entry, it is assumed that matrices A and B have the following */
/* > forms, which may be obtained by the preprocessing subroutine DGGSVP */
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
/* >        U**T *A*Q = D1*( 0 R ),    V**T *B*Q = D2*( 0 R ), */
/* > */
/* > where U, V and Q are orthogonal matrices. */
/* > R is a nonsingular upper triangular matrix, and D1 and D2 are */
/* > ``diagonal'' matrices, which are of the following structures: */
/* > */
/* > If M-K-L >= 0, */
/* > */
/* >                     K  L */
/* >        D1 =     K ( I  0 ) */
/* >                 L ( 0  C ) */
/* >             M-K-L ( 0  0 ) */
/* > */
/* >                   K  L */
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
/* > The computation of the orthogonal transformation matrices U, V or Q */
/* > is optional.  These matrices may either be formed explicitly, or they */
/* > may be postmultiplied into input matrices U1, V1, or Q1. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] JOBU */
/* > \verbatim */
/* >          JOBU is CHARACTER*1 */
/* >          = 'U':  U must contain an orthogonal matrix U1 on entry, and */
/* >                  the product U1*U is returned; */
/* >          = 'I':  U is initialized to the unit matrix, and the */
/* >                  orthogonal matrix U is returned; */
/* >          = 'N':  U is not computed. */
/* > \endverbatim */
/* > */
/* > \param[in] JOBV */
/* > \verbatim */
/* >          JOBV is CHARACTER*1 */
/* >          = 'V':  V must contain an orthogonal matrix V1 on entry, and */
/* >                  the product V1*V is returned; */
/* >          = 'I':  V is initialized to the unit matrix, and the */
/* >                  orthogonal matrix V is returned; */
/* >          = 'N':  V is not computed. */
/* > \endverbatim */
/* > */
/* > \param[in] JOBQ */
/* > \verbatim */
/* >          JOBQ is CHARACTER*1 */
/* >          = 'Q':  Q must contain an orthogonal matrix Q1 on entry, and */
/* >                  the product Q1*Q is returned; */
/* >          = 'I':  Q is initialized to the unit matrix, and the */
/* >                  orthogonal matrix Q is returned; */
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
/* >          A23 = A(K+1:MIN(K+L,M),N-L+1:N) and B13 = B(1:L,N-L+1:N) */
/* >          of A and B, whose GSVD is going to be computed by DTGSJA. */
/* >          See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
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
/* >          B is DOUBLE PRECISION array, dimension (LDB,N) */
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
/* >          TOLA is DOUBLE PRECISION */
/* > \endverbatim */
/* > */
/* > \param[in] TOLB */
/* > \verbatim */
/* >          TOLB is DOUBLE PRECISION */
/* > */
/* >          TOLA and TOLB are the convergence criteria for the Jacobi- */
/* >          Kogbetliantz iteration procedure. Generally, they are the */
/* >          same as used in the preprocessing step, say */
/* >              TOLA = max(M,N)*norm(A)*MAZHEPS, */
/* >              TOLB = max(P,N)*norm(B)*MAZHEPS. */
/* > \endverbatim */
/* > */
/* > \param[out] ALPHA */
/* > \verbatim */
/* >          ALPHA is DOUBLE PRECISION array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] BETA */
/* > \verbatim */
/* >          BETA is DOUBLE PRECISION array, dimension (N) */
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
/* >            ALPHA(K+L+1:N) = 0 and */
/* >            BETA(K+L+1:N)  = 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] U */
/* > \verbatim */
/* >          U is DOUBLE PRECISION array, dimension (LDU,M) */
/* >          On entry, if JOBU = 'U', U must contain a matrix U1 (usually */
/* >          the orthogonal matrix returned by DGGSVP). */
/* >          On exit, */
/* >          if JOBU = 'I', U contains the orthogonal matrix U; */
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
/* >          V is DOUBLE PRECISION array, dimension (LDV,P) */
/* >          On entry, if JOBV = 'V', V must contain a matrix V1 (usually */
/* >          the orthogonal matrix returned by DGGSVP). */
/* >          On exit, */
/* >          if JOBV = 'I', V contains the orthogonal matrix V; */
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
/* >          Q is DOUBLE PRECISION array, dimension (LDQ,N) */
/* >          On entry, if JOBQ = 'Q', Q must contain a matrix Q1 (usually */
/* >          the orthogonal matrix returned by DGGSVP). */
/* >          On exit, */
/* >          if JOBQ = 'I', Q contains the orthogonal matrix Q; */
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
/* >          WORK is DOUBLE PRECISION array, dimension (2*N) */
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
/* > */
/* > \verbatim */
/* >  Internal Parameters */
/* >  =================== */
/* > */
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

/* > \date November 2011 */

/* > \ingroup doubleOTHERcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  DTGSJA essentially uses a variant of Kogbetliantz algorithm to reduce */
/* >  min(L,M-K)-by-L triangular (or trapezoidal) matrix A23 and L-by-L */
/* >  matrix B13 to the form: */
/* > */
/* >           U1**T *A13*Q1 = C1*R1; V1**T *B13*Q1 = S1*R1, */
/* > */
/* >  where U1, V1 and Q1 are orthogonal matrix, and Z**T is the transpose */
/* >  of Z.  C1 and S1 are diagonal matrices satisfying */
/* > */
/* >                C1**2 + S1**2 = I, */
/* > */
/* >  and R1 is an L-by-L nonsingular upper triangular matrix. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dtgsja_(char *jobu, char *jobv, char *jobq, integer *m, 
	integer *p, integer *n, integer *k, integer *l, doublereal *a, 
	integer *lda, doublereal *b, integer *ldb, doublereal *tola, 
	doublereal *tolb, doublereal *alpha, doublereal *beta, doublereal *u, 
	integer *ldu, doublereal *v, integer *ldv, doublereal *q, integer *
	ldq, doublereal *work, integer *ncycle, integer *info, ftnlen 
	jobu_len, ftnlen jobv_len, ftnlen jobq_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, q_dim1, q_offset, u_dim1, 
	    u_offset, v_dim1, v_offset, i__1, i__2, i__3, i__4;
    doublereal d__1;

    /* Local variables */
    static integer i__, j;
    static doublereal a1, a2, a3, b1, b2, b3, csq, csu, csv, snq, rwk, snu, 
	    snv;
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *);
    static doublereal gamma;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static logical initq, initu, initv, wantq, upper;
    static doublereal error, ssmin;
    static logical wantu, wantv;
    extern /* Subroutine */ int dlags2_(logical *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *), dlapll_(integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *);
    static integer kcycle;
    extern /* Subroutine */ int dlartg_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *), dlaset_(char *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, ftnlen), xerbla_(char *, integer *, ftnlen);


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
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Decode and test the input parameters */

#line 428 "dtgsja.f"
    /* Parameter adjustments */
#line 428 "dtgsja.f"
    a_dim1 = *lda;
#line 428 "dtgsja.f"
    a_offset = 1 + a_dim1;
#line 428 "dtgsja.f"
    a -= a_offset;
#line 428 "dtgsja.f"
    b_dim1 = *ldb;
#line 428 "dtgsja.f"
    b_offset = 1 + b_dim1;
#line 428 "dtgsja.f"
    b -= b_offset;
#line 428 "dtgsja.f"
    --alpha;
#line 428 "dtgsja.f"
    --beta;
#line 428 "dtgsja.f"
    u_dim1 = *ldu;
#line 428 "dtgsja.f"
    u_offset = 1 + u_dim1;
#line 428 "dtgsja.f"
    u -= u_offset;
#line 428 "dtgsja.f"
    v_dim1 = *ldv;
#line 428 "dtgsja.f"
    v_offset = 1 + v_dim1;
#line 428 "dtgsja.f"
    v -= v_offset;
#line 428 "dtgsja.f"
    q_dim1 = *ldq;
#line 428 "dtgsja.f"
    q_offset = 1 + q_dim1;
#line 428 "dtgsja.f"
    q -= q_offset;
#line 428 "dtgsja.f"
    --work;
#line 428 "dtgsja.f"

#line 428 "dtgsja.f"
    /* Function Body */
#line 428 "dtgsja.f"
    initu = lsame_(jobu, "I", (ftnlen)1, (ftnlen)1);
#line 429 "dtgsja.f"
    wantu = initu || lsame_(jobu, "U", (ftnlen)1, (ftnlen)1);

#line 431 "dtgsja.f"
    initv = lsame_(jobv, "I", (ftnlen)1, (ftnlen)1);
#line 432 "dtgsja.f"
    wantv = initv || lsame_(jobv, "V", (ftnlen)1, (ftnlen)1);

#line 434 "dtgsja.f"
    initq = lsame_(jobq, "I", (ftnlen)1, (ftnlen)1);
#line 435 "dtgsja.f"
    wantq = initq || lsame_(jobq, "Q", (ftnlen)1, (ftnlen)1);

#line 437 "dtgsja.f"
    *info = 0;
#line 438 "dtgsja.f"
    if (! (initu || wantu || lsame_(jobu, "N", (ftnlen)1, (ftnlen)1))) {
#line 439 "dtgsja.f"
	*info = -1;
#line 440 "dtgsja.f"
    } else if (! (initv || wantv || lsame_(jobv, "N", (ftnlen)1, (ftnlen)1))) 
	    {
#line 441 "dtgsja.f"
	*info = -2;
#line 442 "dtgsja.f"
    } else if (! (initq || wantq || lsame_(jobq, "N", (ftnlen)1, (ftnlen)1))) 
	    {
#line 443 "dtgsja.f"
	*info = -3;
#line 444 "dtgsja.f"
    } else if (*m < 0) {
#line 445 "dtgsja.f"
	*info = -4;
#line 446 "dtgsja.f"
    } else if (*p < 0) {
#line 447 "dtgsja.f"
	*info = -5;
#line 448 "dtgsja.f"
    } else if (*n < 0) {
#line 449 "dtgsja.f"
	*info = -6;
#line 450 "dtgsja.f"
    } else if (*lda < max(1,*m)) {
#line 451 "dtgsja.f"
	*info = -10;
#line 452 "dtgsja.f"
    } else if (*ldb < max(1,*p)) {
#line 453 "dtgsja.f"
	*info = -12;
#line 454 "dtgsja.f"
    } else if (*ldu < 1 || wantu && *ldu < *m) {
#line 455 "dtgsja.f"
	*info = -18;
#line 456 "dtgsja.f"
    } else if (*ldv < 1 || wantv && *ldv < *p) {
#line 457 "dtgsja.f"
	*info = -20;
#line 458 "dtgsja.f"
    } else if (*ldq < 1 || wantq && *ldq < *n) {
#line 459 "dtgsja.f"
	*info = -22;
#line 460 "dtgsja.f"
    }
#line 461 "dtgsja.f"
    if (*info != 0) {
#line 462 "dtgsja.f"
	i__1 = -(*info);
#line 462 "dtgsja.f"
	xerbla_("DTGSJA", &i__1, (ftnlen)6);
#line 463 "dtgsja.f"
	return 0;
#line 464 "dtgsja.f"
    }

/*     Initialize U, V and Q, if necessary */

#line 468 "dtgsja.f"
    if (initu) {
#line 468 "dtgsja.f"
	dlaset_("Full", m, m, &c_b13, &c_b14, &u[u_offset], ldu, (ftnlen)4);
#line 468 "dtgsja.f"
    }
#line 470 "dtgsja.f"
    if (initv) {
#line 470 "dtgsja.f"
	dlaset_("Full", p, p, &c_b13, &c_b14, &v[v_offset], ldv, (ftnlen)4);
#line 470 "dtgsja.f"
    }
#line 472 "dtgsja.f"
    if (initq) {
#line 472 "dtgsja.f"
	dlaset_("Full", n, n, &c_b13, &c_b14, &q[q_offset], ldq, (ftnlen)4);
#line 472 "dtgsja.f"
    }

/*     Loop until convergence */

#line 477 "dtgsja.f"
    upper = FALSE_;
#line 478 "dtgsja.f"
    for (kcycle = 1; kcycle <= 40; ++kcycle) {

#line 480 "dtgsja.f"
	upper = ! upper;

#line 482 "dtgsja.f"
	i__1 = *l - 1;
#line 482 "dtgsja.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 483 "dtgsja.f"
	    i__2 = *l;
#line 483 "dtgsja.f"
	    for (j = i__ + 1; j <= i__2; ++j) {

#line 485 "dtgsja.f"
		a1 = 0.;
#line 486 "dtgsja.f"
		a2 = 0.;
#line 487 "dtgsja.f"
		a3 = 0.;
#line 488 "dtgsja.f"
		if (*k + i__ <= *m) {
#line 488 "dtgsja.f"
		    a1 = a[*k + i__ + (*n - *l + i__) * a_dim1];
#line 488 "dtgsja.f"
		}
#line 490 "dtgsja.f"
		if (*k + j <= *m) {
#line 490 "dtgsja.f"
		    a3 = a[*k + j + (*n - *l + j) * a_dim1];
#line 490 "dtgsja.f"
		}

#line 493 "dtgsja.f"
		b1 = b[i__ + (*n - *l + i__) * b_dim1];
#line 494 "dtgsja.f"
		b3 = b[j + (*n - *l + j) * b_dim1];

#line 496 "dtgsja.f"
		if (upper) {
#line 497 "dtgsja.f"
		    if (*k + i__ <= *m) {
#line 497 "dtgsja.f"
			a2 = a[*k + i__ + (*n - *l + j) * a_dim1];
#line 497 "dtgsja.f"
		    }
#line 499 "dtgsja.f"
		    b2 = b[i__ + (*n - *l + j) * b_dim1];
#line 500 "dtgsja.f"
		} else {
#line 501 "dtgsja.f"
		    if (*k + j <= *m) {
#line 501 "dtgsja.f"
			a2 = a[*k + j + (*n - *l + i__) * a_dim1];
#line 501 "dtgsja.f"
		    }
#line 503 "dtgsja.f"
		    b2 = b[j + (*n - *l + i__) * b_dim1];
#line 504 "dtgsja.f"
		}

#line 506 "dtgsja.f"
		dlags2_(&upper, &a1, &a2, &a3, &b1, &b2, &b3, &csu, &snu, &
			csv, &snv, &csq, &snq);

/*              Update (K+I)-th and (K+J)-th rows of matrix A: U**T *A */

#line 511 "dtgsja.f"
		if (*k + j <= *m) {
#line 511 "dtgsja.f"
		    drot_(l, &a[*k + j + (*n - *l + 1) * a_dim1], lda, &a[*k 
			    + i__ + (*n - *l + 1) * a_dim1], lda, &csu, &snu);
#line 511 "dtgsja.f"
		}

/*              Update I-th and J-th rows of matrix B: V**T *B */

#line 517 "dtgsja.f"
		drot_(l, &b[j + (*n - *l + 1) * b_dim1], ldb, &b[i__ + (*n - *
			l + 1) * b_dim1], ldb, &csv, &snv);

/*              Update (N-L+I)-th and (N-L+J)-th columns of matrices */
/*              A and B: A*Q and B*Q */

/* Computing MIN */
#line 523 "dtgsja.f"
		i__4 = *k + *l;
#line 523 "dtgsja.f"
		i__3 = min(i__4,*m);
#line 523 "dtgsja.f"
		drot_(&i__3, &a[(*n - *l + j) * a_dim1 + 1], &c__1, &a[(*n - *
			l + i__) * a_dim1 + 1], &c__1, &csq, &snq);

#line 526 "dtgsja.f"
		drot_(l, &b[(*n - *l + j) * b_dim1 + 1], &c__1, &b[(*n - *l + 
			i__) * b_dim1 + 1], &c__1, &csq, &snq);

#line 529 "dtgsja.f"
		if (upper) {
#line 530 "dtgsja.f"
		    if (*k + i__ <= *m) {
#line 530 "dtgsja.f"
			a[*k + i__ + (*n - *l + j) * a_dim1] = 0.;
#line 530 "dtgsja.f"
		    }
#line 532 "dtgsja.f"
		    b[i__ + (*n - *l + j) * b_dim1] = 0.;
#line 533 "dtgsja.f"
		} else {
#line 534 "dtgsja.f"
		    if (*k + j <= *m) {
#line 534 "dtgsja.f"
			a[*k + j + (*n - *l + i__) * a_dim1] = 0.;
#line 534 "dtgsja.f"
		    }
#line 536 "dtgsja.f"
		    b[j + (*n - *l + i__) * b_dim1] = 0.;
#line 537 "dtgsja.f"
		}

/*              Update orthogonal matrices U, V, Q, if desired. */

#line 541 "dtgsja.f"
		if (wantu && *k + j <= *m) {
#line 541 "dtgsja.f"
		    drot_(m, &u[(*k + j) * u_dim1 + 1], &c__1, &u[(*k + i__) *
			     u_dim1 + 1], &c__1, &csu, &snu);
#line 541 "dtgsja.f"
		}

#line 545 "dtgsja.f"
		if (wantv) {
#line 545 "dtgsja.f"
		    drot_(p, &v[j * v_dim1 + 1], &c__1, &v[i__ * v_dim1 + 1], 
			    &c__1, &csv, &snv);
#line 545 "dtgsja.f"
		}

#line 548 "dtgsja.f"
		if (wantq) {
#line 548 "dtgsja.f"
		    drot_(n, &q[(*n - *l + j) * q_dim1 + 1], &c__1, &q[(*n - *
			    l + i__) * q_dim1 + 1], &c__1, &csq, &snq);
#line 548 "dtgsja.f"
		}

#line 552 "dtgsja.f"
/* L10: */
#line 552 "dtgsja.f"
	    }
#line 553 "dtgsja.f"
/* L20: */
#line 553 "dtgsja.f"
	}

#line 555 "dtgsja.f"
	if (! upper) {

/*           The matrices A13 and B13 were lower triangular at the start */
/*           of the cycle, and are now upper triangular. */

/*           Convergence test: test the parallelism of the corresponding */
/*           rows of A and B. */

#line 563 "dtgsja.f"
	    error = 0.;
/* Computing MIN */
#line 564 "dtgsja.f"
	    i__2 = *l, i__3 = *m - *k;
#line 564 "dtgsja.f"
	    i__1 = min(i__2,i__3);
#line 564 "dtgsja.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 565 "dtgsja.f"
		i__2 = *l - i__ + 1;
#line 565 "dtgsja.f"
		dcopy_(&i__2, &a[*k + i__ + (*n - *l + i__) * a_dim1], lda, &
			work[1], &c__1);
#line 566 "dtgsja.f"
		i__2 = *l - i__ + 1;
#line 566 "dtgsja.f"
		dcopy_(&i__2, &b[i__ + (*n - *l + i__) * b_dim1], ldb, &work[*
			l + 1], &c__1);
#line 567 "dtgsja.f"
		i__2 = *l - i__ + 1;
#line 567 "dtgsja.f"
		dlapll_(&i__2, &work[1], &c__1, &work[*l + 1], &c__1, &ssmin);
#line 568 "dtgsja.f"
		error = max(error,ssmin);
#line 569 "dtgsja.f"
/* L30: */
#line 569 "dtgsja.f"
	    }

#line 571 "dtgsja.f"
	    if (abs(error) <= min(*tola,*tolb)) {
#line 571 "dtgsja.f"
		goto L50;
#line 571 "dtgsja.f"
	    }
#line 573 "dtgsja.f"
	}

/*        End of cycle loop */

#line 577 "dtgsja.f"
/* L40: */
#line 577 "dtgsja.f"
    }

/*     The algorithm has not converged after MAXIT cycles. */

#line 581 "dtgsja.f"
    *info = 1;
#line 582 "dtgsja.f"
    goto L100;

#line 584 "dtgsja.f"
L50:

/*     If ERROR <= MIN(TOLA,TOLB), then the algorithm has converged. */
/*     Compute the generalized singular value pairs (ALPHA, BETA), and */
/*     set the triangular matrix R to array A. */

#line 590 "dtgsja.f"
    i__1 = *k;
#line 590 "dtgsja.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 591 "dtgsja.f"
	alpha[i__] = 1.;
#line 592 "dtgsja.f"
	beta[i__] = 0.;
#line 593 "dtgsja.f"
/* L60: */
#line 593 "dtgsja.f"
    }

/* Computing MIN */
#line 595 "dtgsja.f"
    i__2 = *l, i__3 = *m - *k;
#line 595 "dtgsja.f"
    i__1 = min(i__2,i__3);
#line 595 "dtgsja.f"
    for (i__ = 1; i__ <= i__1; ++i__) {

#line 597 "dtgsja.f"
	a1 = a[*k + i__ + (*n - *l + i__) * a_dim1];
#line 598 "dtgsja.f"
	b1 = b[i__ + (*n - *l + i__) * b_dim1];

#line 600 "dtgsja.f"
	if (a1 != 0.) {
#line 601 "dtgsja.f"
	    gamma = b1 / a1;

/*           change sign if necessary */

#line 605 "dtgsja.f"
	    if (gamma < 0.) {
#line 606 "dtgsja.f"
		i__2 = *l - i__ + 1;
#line 606 "dtgsja.f"
		dscal_(&i__2, &c_b43, &b[i__ + (*n - *l + i__) * b_dim1], ldb)
			;
#line 607 "dtgsja.f"
		if (wantv) {
#line 607 "dtgsja.f"
		    dscal_(p, &c_b43, &v[i__ * v_dim1 + 1], &c__1);
#line 607 "dtgsja.f"
		}
#line 609 "dtgsja.f"
	    }

#line 611 "dtgsja.f"
	    d__1 = abs(gamma);
#line 611 "dtgsja.f"
	    dlartg_(&d__1, &c_b14, &beta[*k + i__], &alpha[*k + i__], &rwk);

#line 614 "dtgsja.f"
	    if (alpha[*k + i__] >= beta[*k + i__]) {
#line 615 "dtgsja.f"
		i__2 = *l - i__ + 1;
#line 615 "dtgsja.f"
		d__1 = 1. / alpha[*k + i__];
#line 615 "dtgsja.f"
		dscal_(&i__2, &d__1, &a[*k + i__ + (*n - *l + i__) * a_dim1], 
			lda);
#line 617 "dtgsja.f"
	    } else {
#line 618 "dtgsja.f"
		i__2 = *l - i__ + 1;
#line 618 "dtgsja.f"
		d__1 = 1. / beta[*k + i__];
#line 618 "dtgsja.f"
		dscal_(&i__2, &d__1, &b[i__ + (*n - *l + i__) * b_dim1], ldb);
#line 620 "dtgsja.f"
		i__2 = *l - i__ + 1;
#line 620 "dtgsja.f"
		dcopy_(&i__2, &b[i__ + (*n - *l + i__) * b_dim1], ldb, &a[*k 
			+ i__ + (*n - *l + i__) * a_dim1], lda);
#line 622 "dtgsja.f"
	    }

#line 624 "dtgsja.f"
	} else {

#line 626 "dtgsja.f"
	    alpha[*k + i__] = 0.;
#line 627 "dtgsja.f"
	    beta[*k + i__] = 1.;
#line 628 "dtgsja.f"
	    i__2 = *l - i__ + 1;
#line 628 "dtgsja.f"
	    dcopy_(&i__2, &b[i__ + (*n - *l + i__) * b_dim1], ldb, &a[*k + 
		    i__ + (*n - *l + i__) * a_dim1], lda);

#line 631 "dtgsja.f"
	}

#line 633 "dtgsja.f"
/* L70: */
#line 633 "dtgsja.f"
    }

/*     Post-assignment */

#line 637 "dtgsja.f"
    i__1 = *k + *l;
#line 637 "dtgsja.f"
    for (i__ = *m + 1; i__ <= i__1; ++i__) {
#line 638 "dtgsja.f"
	alpha[i__] = 0.;
#line 639 "dtgsja.f"
	beta[i__] = 1.;
#line 640 "dtgsja.f"
/* L80: */
#line 640 "dtgsja.f"
    }

#line 642 "dtgsja.f"
    if (*k + *l < *n) {
#line 643 "dtgsja.f"
	i__1 = *n;
#line 643 "dtgsja.f"
	for (i__ = *k + *l + 1; i__ <= i__1; ++i__) {
#line 644 "dtgsja.f"
	    alpha[i__] = 0.;
#line 645 "dtgsja.f"
	    beta[i__] = 0.;
#line 646 "dtgsja.f"
/* L90: */
#line 646 "dtgsja.f"
	}
#line 647 "dtgsja.f"
    }

#line 649 "dtgsja.f"
L100:
#line 650 "dtgsja.f"
    *ncycle = kcycle;
#line 651 "dtgsja.f"
    return 0;

/*     End of DTGSJA */

} /* dtgsja_ */

