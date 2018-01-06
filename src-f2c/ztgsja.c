#line 1 "ztgsja.f"
/* ztgsja.f -- translated by f2c (version 20100827).
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

#line 1 "ztgsja.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};
static integer c__1 = 1;
static doublereal c_b39 = -1.;
static doublereal c_b42 = 1.;

/* > \brief \b ZTGSJA */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZTGSJA + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztgsja.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztgsja.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztgsja.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZTGSJA( JOBU, JOBV, JOBQ, M, P, N, K, L, A, LDA, B, */
/*                          LDB, TOLA, TOLB, ALPHA, BETA, U, LDU, V, LDV, */
/*                          Q, LDQ, WORK, NCYCLE, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBQ, JOBU, JOBV */
/*       INTEGER            INFO, K, L, LDA, LDB, LDQ, LDU, LDV, M, N, */
/*      $                   NCYCLE, P */
/*       DOUBLE PRECISION   TOLA, TOLB */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   ALPHA( * ), BETA( * ) */
/*       COMPLEX*16         A( LDA, * ), B( LDB, * ), Q( LDQ, * ), */
/*      $                   U( LDU, * ), V( LDV, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZTGSJA computes the generalized singular value decomposition (GSVD) */
/* > of two complex upper triangular (or trapezoidal) matrices A and B. */
/* > */
/* > On entry, it is assumed that matrices A and B have the following */
/* > forms, which may be obtained by the preprocessing subroutine ZGGSVP */
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
/* >          of A and B, whose GSVD is going to be computed by ZTGSJA. */
/* >          See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
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
/* >          B is COMPLEX*16 array, dimension (LDB,N) */
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
/* >              TOLA = MAX(M,N)*norm(A)*MAZHEPS, */
/* >              TOLB = MAX(P,N)*norm(B)*MAZHEPS. */
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
/* >          U is COMPLEX*16 array, dimension (LDU,M) */
/* >          On entry, if JOBU = 'U', U must contain a matrix U1 (usually */
/* >          the unitary matrix returned by ZGGSVP). */
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
/* >          V is COMPLEX*16 array, dimension (LDV,P) */
/* >          On entry, if JOBV = 'V', V must contain a matrix V1 (usually */
/* >          the unitary matrix returned by ZGGSVP). */
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
/* >          Q is COMPLEX*16 array, dimension (LDQ,N) */
/* >          On entry, if JOBQ = 'Q', Q must contain a matrix Q1 (usually */
/* >          the unitary matrix returned by ZGGSVP). */
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
/* >          WORK is COMPLEX*16 array, dimension (2*N) */
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

/* > \date November 2011 */

/* > \ingroup complex16OTHERcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  ZTGSJA essentially uses a variant of Kogbetliantz algorithm to reduce */
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
/* Subroutine */ int ztgsja_(char *jobu, char *jobv, char *jobq, integer *m, 
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
    extern /* Subroutine */ int zrot_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublereal *, doublecomplex *);
    static doublereal gamma;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical initq, initu, initv, wantq, upper;
    static doublereal error, ssmin;
    static logical wantu, wantv;
    extern /* Subroutine */ int zcopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), zlags2_(logical *, doublereal *, 
	    doublecomplex *, doublereal *, doublereal *, doublecomplex *, 
	    doublereal *, doublereal *, doublecomplex *, doublereal *, 
	    doublecomplex *, doublereal *, doublecomplex *);
    static integer kcycle;
    extern /* Subroutine */ int dlartg_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *), xerbla_(char *, 
	    integer *, ftnlen), zdscal_(integer *, doublereal *, 
	    doublecomplex *, integer *), zlapll_(integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, doublereal *), zlaset_(
	    char *, integer *, integer *, doublecomplex *, doublecomplex *, 
	    doublecomplex *, integer *, ftnlen);


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

#line 433 "ztgsja.f"
    /* Parameter adjustments */
#line 433 "ztgsja.f"
    a_dim1 = *lda;
#line 433 "ztgsja.f"
    a_offset = 1 + a_dim1;
#line 433 "ztgsja.f"
    a -= a_offset;
#line 433 "ztgsja.f"
    b_dim1 = *ldb;
#line 433 "ztgsja.f"
    b_offset = 1 + b_dim1;
#line 433 "ztgsja.f"
    b -= b_offset;
#line 433 "ztgsja.f"
    --alpha;
#line 433 "ztgsja.f"
    --beta;
#line 433 "ztgsja.f"
    u_dim1 = *ldu;
#line 433 "ztgsja.f"
    u_offset = 1 + u_dim1;
#line 433 "ztgsja.f"
    u -= u_offset;
#line 433 "ztgsja.f"
    v_dim1 = *ldv;
#line 433 "ztgsja.f"
    v_offset = 1 + v_dim1;
#line 433 "ztgsja.f"
    v -= v_offset;
#line 433 "ztgsja.f"
    q_dim1 = *ldq;
#line 433 "ztgsja.f"
    q_offset = 1 + q_dim1;
#line 433 "ztgsja.f"
    q -= q_offset;
#line 433 "ztgsja.f"
    --work;
#line 433 "ztgsja.f"

#line 433 "ztgsja.f"
    /* Function Body */
#line 433 "ztgsja.f"
    initu = lsame_(jobu, "I", (ftnlen)1, (ftnlen)1);
#line 434 "ztgsja.f"
    wantu = initu || lsame_(jobu, "U", (ftnlen)1, (ftnlen)1);

#line 436 "ztgsja.f"
    initv = lsame_(jobv, "I", (ftnlen)1, (ftnlen)1);
#line 437 "ztgsja.f"
    wantv = initv || lsame_(jobv, "V", (ftnlen)1, (ftnlen)1);

#line 439 "ztgsja.f"
    initq = lsame_(jobq, "I", (ftnlen)1, (ftnlen)1);
#line 440 "ztgsja.f"
    wantq = initq || lsame_(jobq, "Q", (ftnlen)1, (ftnlen)1);

#line 442 "ztgsja.f"
    *info = 0;
#line 443 "ztgsja.f"
    if (! (initu || wantu || lsame_(jobu, "N", (ftnlen)1, (ftnlen)1))) {
#line 444 "ztgsja.f"
	*info = -1;
#line 445 "ztgsja.f"
    } else if (! (initv || wantv || lsame_(jobv, "N", (ftnlen)1, (ftnlen)1))) 
	    {
#line 446 "ztgsja.f"
	*info = -2;
#line 447 "ztgsja.f"
    } else if (! (initq || wantq || lsame_(jobq, "N", (ftnlen)1, (ftnlen)1))) 
	    {
#line 448 "ztgsja.f"
	*info = -3;
#line 449 "ztgsja.f"
    } else if (*m < 0) {
#line 450 "ztgsja.f"
	*info = -4;
#line 451 "ztgsja.f"
    } else if (*p < 0) {
#line 452 "ztgsja.f"
	*info = -5;
#line 453 "ztgsja.f"
    } else if (*n < 0) {
#line 454 "ztgsja.f"
	*info = -6;
#line 455 "ztgsja.f"
    } else if (*lda < max(1,*m)) {
#line 456 "ztgsja.f"
	*info = -10;
#line 457 "ztgsja.f"
    } else if (*ldb < max(1,*p)) {
#line 458 "ztgsja.f"
	*info = -12;
#line 459 "ztgsja.f"
    } else if (*ldu < 1 || wantu && *ldu < *m) {
#line 460 "ztgsja.f"
	*info = -18;
#line 461 "ztgsja.f"
    } else if (*ldv < 1 || wantv && *ldv < *p) {
#line 462 "ztgsja.f"
	*info = -20;
#line 463 "ztgsja.f"
    } else if (*ldq < 1 || wantq && *ldq < *n) {
#line 464 "ztgsja.f"
	*info = -22;
#line 465 "ztgsja.f"
    }
#line 466 "ztgsja.f"
    if (*info != 0) {
#line 467 "ztgsja.f"
	i__1 = -(*info);
#line 467 "ztgsja.f"
	xerbla_("ZTGSJA", &i__1, (ftnlen)6);
#line 468 "ztgsja.f"
	return 0;
#line 469 "ztgsja.f"
    }

/*     Initialize U, V and Q, if necessary */

#line 473 "ztgsja.f"
    if (initu) {
#line 473 "ztgsja.f"
	zlaset_("Full", m, m, &c_b1, &c_b2, &u[u_offset], ldu, (ftnlen)4);
#line 473 "ztgsja.f"
    }
#line 475 "ztgsja.f"
    if (initv) {
#line 475 "ztgsja.f"
	zlaset_("Full", p, p, &c_b1, &c_b2, &v[v_offset], ldv, (ftnlen)4);
#line 475 "ztgsja.f"
    }
#line 477 "ztgsja.f"
    if (initq) {
#line 477 "ztgsja.f"
	zlaset_("Full", n, n, &c_b1, &c_b2, &q[q_offset], ldq, (ftnlen)4);
#line 477 "ztgsja.f"
    }

/*     Loop until convergence */

#line 482 "ztgsja.f"
    upper = FALSE_;
#line 483 "ztgsja.f"
    for (kcycle = 1; kcycle <= 40; ++kcycle) {

#line 485 "ztgsja.f"
	upper = ! upper;

#line 487 "ztgsja.f"
	i__1 = *l - 1;
#line 487 "ztgsja.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 488 "ztgsja.f"
	    i__2 = *l;
#line 488 "ztgsja.f"
	    for (j = i__ + 1; j <= i__2; ++j) {

#line 490 "ztgsja.f"
		a1 = 0.;
#line 491 "ztgsja.f"
		a2.r = 0., a2.i = 0.;
#line 492 "ztgsja.f"
		a3 = 0.;
#line 493 "ztgsja.f"
		if (*k + i__ <= *m) {
#line 493 "ztgsja.f"
		    i__3 = *k + i__ + (*n - *l + i__) * a_dim1;
#line 493 "ztgsja.f"
		    a1 = a[i__3].r;
#line 493 "ztgsja.f"
		}
#line 495 "ztgsja.f"
		if (*k + j <= *m) {
#line 495 "ztgsja.f"
		    i__3 = *k + j + (*n - *l + j) * a_dim1;
#line 495 "ztgsja.f"
		    a3 = a[i__3].r;
#line 495 "ztgsja.f"
		}

#line 498 "ztgsja.f"
		i__3 = i__ + (*n - *l + i__) * b_dim1;
#line 498 "ztgsja.f"
		b1 = b[i__3].r;
#line 499 "ztgsja.f"
		i__3 = j + (*n - *l + j) * b_dim1;
#line 499 "ztgsja.f"
		b3 = b[i__3].r;

#line 501 "ztgsja.f"
		if (upper) {
#line 502 "ztgsja.f"
		    if (*k + i__ <= *m) {
#line 502 "ztgsja.f"
			i__3 = *k + i__ + (*n - *l + j) * a_dim1;
#line 502 "ztgsja.f"
			a2.r = a[i__3].r, a2.i = a[i__3].i;
#line 502 "ztgsja.f"
		    }
#line 504 "ztgsja.f"
		    i__3 = i__ + (*n - *l + j) * b_dim1;
#line 504 "ztgsja.f"
		    b2.r = b[i__3].r, b2.i = b[i__3].i;
#line 505 "ztgsja.f"
		} else {
#line 506 "ztgsja.f"
		    if (*k + j <= *m) {
#line 506 "ztgsja.f"
			i__3 = *k + j + (*n - *l + i__) * a_dim1;
#line 506 "ztgsja.f"
			a2.r = a[i__3].r, a2.i = a[i__3].i;
#line 506 "ztgsja.f"
		    }
#line 508 "ztgsja.f"
		    i__3 = j + (*n - *l + i__) * b_dim1;
#line 508 "ztgsja.f"
		    b2.r = b[i__3].r, b2.i = b[i__3].i;
#line 509 "ztgsja.f"
		}

#line 511 "ztgsja.f"
		zlags2_(&upper, &a1, &a2, &a3, &b1, &b2, &b3, &csu, &snu, &
			csv, &snv, &csq, &snq);

/*              Update (K+I)-th and (K+J)-th rows of matrix A: U**H *A */

#line 516 "ztgsja.f"
		if (*k + j <= *m) {
#line 516 "ztgsja.f"
		    d_cnjg(&z__1, &snu);
#line 516 "ztgsja.f"
		    zrot_(l, &a[*k + j + (*n - *l + 1) * a_dim1], lda, &a[*k 
			    + i__ + (*n - *l + 1) * a_dim1], lda, &csu, &z__1)
			    ;
#line 516 "ztgsja.f"
		}

/*              Update I-th and J-th rows of matrix B: V**H *B */

#line 522 "ztgsja.f"
		d_cnjg(&z__1, &snv);
#line 522 "ztgsja.f"
		zrot_(l, &b[j + (*n - *l + 1) * b_dim1], ldb, &b[i__ + (*n - *
			l + 1) * b_dim1], ldb, &csv, &z__1);

/*              Update (N-L+I)-th and (N-L+J)-th columns of matrices */
/*              A and B: A*Q and B*Q */

/* Computing MIN */
#line 528 "ztgsja.f"
		i__4 = *k + *l;
#line 528 "ztgsja.f"
		i__3 = min(i__4,*m);
#line 528 "ztgsja.f"
		zrot_(&i__3, &a[(*n - *l + j) * a_dim1 + 1], &c__1, &a[(*n - *
			l + i__) * a_dim1 + 1], &c__1, &csq, &snq);

#line 531 "ztgsja.f"
		zrot_(l, &b[(*n - *l + j) * b_dim1 + 1], &c__1, &b[(*n - *l + 
			i__) * b_dim1 + 1], &c__1, &csq, &snq);

#line 534 "ztgsja.f"
		if (upper) {
#line 535 "ztgsja.f"
		    if (*k + i__ <= *m) {
#line 535 "ztgsja.f"
			i__3 = *k + i__ + (*n - *l + j) * a_dim1;
#line 535 "ztgsja.f"
			a[i__3].r = 0., a[i__3].i = 0.;
#line 535 "ztgsja.f"
		    }
#line 537 "ztgsja.f"
		    i__3 = i__ + (*n - *l + j) * b_dim1;
#line 537 "ztgsja.f"
		    b[i__3].r = 0., b[i__3].i = 0.;
#line 538 "ztgsja.f"
		} else {
#line 539 "ztgsja.f"
		    if (*k + j <= *m) {
#line 539 "ztgsja.f"
			i__3 = *k + j + (*n - *l + i__) * a_dim1;
#line 539 "ztgsja.f"
			a[i__3].r = 0., a[i__3].i = 0.;
#line 539 "ztgsja.f"
		    }
#line 541 "ztgsja.f"
		    i__3 = j + (*n - *l + i__) * b_dim1;
#line 541 "ztgsja.f"
		    b[i__3].r = 0., b[i__3].i = 0.;
#line 542 "ztgsja.f"
		}

/*              Ensure that the diagonal elements of A and B are real. */

#line 546 "ztgsja.f"
		if (*k + i__ <= *m) {
#line 546 "ztgsja.f"
		    i__3 = *k + i__ + (*n - *l + i__) * a_dim1;
#line 546 "ztgsja.f"
		    i__4 = *k + i__ + (*n - *l + i__) * a_dim1;
#line 546 "ztgsja.f"
		    d__1 = a[i__4].r;
#line 546 "ztgsja.f"
		    a[i__3].r = d__1, a[i__3].i = 0.;
#line 546 "ztgsja.f"
		}
#line 548 "ztgsja.f"
		if (*k + j <= *m) {
#line 548 "ztgsja.f"
		    i__3 = *k + j + (*n - *l + j) * a_dim1;
#line 548 "ztgsja.f"
		    i__4 = *k + j + (*n - *l + j) * a_dim1;
#line 548 "ztgsja.f"
		    d__1 = a[i__4].r;
#line 548 "ztgsja.f"
		    a[i__3].r = d__1, a[i__3].i = 0.;
#line 548 "ztgsja.f"
		}
#line 550 "ztgsja.f"
		i__3 = i__ + (*n - *l + i__) * b_dim1;
#line 550 "ztgsja.f"
		i__4 = i__ + (*n - *l + i__) * b_dim1;
#line 550 "ztgsja.f"
		d__1 = b[i__4].r;
#line 550 "ztgsja.f"
		b[i__3].r = d__1, b[i__3].i = 0.;
#line 551 "ztgsja.f"
		i__3 = j + (*n - *l + j) * b_dim1;
#line 551 "ztgsja.f"
		i__4 = j + (*n - *l + j) * b_dim1;
#line 551 "ztgsja.f"
		d__1 = b[i__4].r;
#line 551 "ztgsja.f"
		b[i__3].r = d__1, b[i__3].i = 0.;

/*              Update unitary matrices U, V, Q, if desired. */

#line 555 "ztgsja.f"
		if (wantu && *k + j <= *m) {
#line 555 "ztgsja.f"
		    zrot_(m, &u[(*k + j) * u_dim1 + 1], &c__1, &u[(*k + i__) *
			     u_dim1 + 1], &c__1, &csu, &snu);
#line 555 "ztgsja.f"
		}

#line 559 "ztgsja.f"
		if (wantv) {
#line 559 "ztgsja.f"
		    zrot_(p, &v[j * v_dim1 + 1], &c__1, &v[i__ * v_dim1 + 1], 
			    &c__1, &csv, &snv);
#line 559 "ztgsja.f"
		}

#line 562 "ztgsja.f"
		if (wantq) {
#line 562 "ztgsja.f"
		    zrot_(n, &q[(*n - *l + j) * q_dim1 + 1], &c__1, &q[(*n - *
			    l + i__) * q_dim1 + 1], &c__1, &csq, &snq);
#line 562 "ztgsja.f"
		}

#line 566 "ztgsja.f"
/* L10: */
#line 566 "ztgsja.f"
	    }
#line 567 "ztgsja.f"
/* L20: */
#line 567 "ztgsja.f"
	}

#line 569 "ztgsja.f"
	if (! upper) {

/*           The matrices A13 and B13 were lower triangular at the start */
/*           of the cycle, and are now upper triangular. */

/*           Convergence test: test the parallelism of the corresponding */
/*           rows of A and B. */

#line 577 "ztgsja.f"
	    error = 0.;
/* Computing MIN */
#line 578 "ztgsja.f"
	    i__2 = *l, i__3 = *m - *k;
#line 578 "ztgsja.f"
	    i__1 = min(i__2,i__3);
#line 578 "ztgsja.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 579 "ztgsja.f"
		i__2 = *l - i__ + 1;
#line 579 "ztgsja.f"
		zcopy_(&i__2, &a[*k + i__ + (*n - *l + i__) * a_dim1], lda, &
			work[1], &c__1);
#line 580 "ztgsja.f"
		i__2 = *l - i__ + 1;
#line 580 "ztgsja.f"
		zcopy_(&i__2, &b[i__ + (*n - *l + i__) * b_dim1], ldb, &work[*
			l + 1], &c__1);
#line 581 "ztgsja.f"
		i__2 = *l - i__ + 1;
#line 581 "ztgsja.f"
		zlapll_(&i__2, &work[1], &c__1, &work[*l + 1], &c__1, &ssmin);
#line 582 "ztgsja.f"
		error = max(error,ssmin);
#line 583 "ztgsja.f"
/* L30: */
#line 583 "ztgsja.f"
	    }

#line 585 "ztgsja.f"
	    if (abs(error) <= min(*tola,*tolb)) {
#line 585 "ztgsja.f"
		goto L50;
#line 585 "ztgsja.f"
	    }
#line 587 "ztgsja.f"
	}

/*        End of cycle loop */

#line 591 "ztgsja.f"
/* L40: */
#line 591 "ztgsja.f"
    }

/*     The algorithm has not converged after MAXIT cycles. */

#line 595 "ztgsja.f"
    *info = 1;
#line 596 "ztgsja.f"
    goto L100;

#line 598 "ztgsja.f"
L50:

/*     If ERROR <= MIN(TOLA,TOLB), then the algorithm has converged. */
/*     Compute the generalized singular value pairs (ALPHA, BETA), and */
/*     set the triangular matrix R to array A. */

#line 604 "ztgsja.f"
    i__1 = *k;
#line 604 "ztgsja.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 605 "ztgsja.f"
	alpha[i__] = 1.;
#line 606 "ztgsja.f"
	beta[i__] = 0.;
#line 607 "ztgsja.f"
/* L60: */
#line 607 "ztgsja.f"
    }

/* Computing MIN */
#line 609 "ztgsja.f"
    i__2 = *l, i__3 = *m - *k;
#line 609 "ztgsja.f"
    i__1 = min(i__2,i__3);
#line 609 "ztgsja.f"
    for (i__ = 1; i__ <= i__1; ++i__) {

#line 611 "ztgsja.f"
	i__2 = *k + i__ + (*n - *l + i__) * a_dim1;
#line 611 "ztgsja.f"
	a1 = a[i__2].r;
#line 612 "ztgsja.f"
	i__2 = i__ + (*n - *l + i__) * b_dim1;
#line 612 "ztgsja.f"
	b1 = b[i__2].r;

#line 614 "ztgsja.f"
	if (a1 != 0.) {
#line 615 "ztgsja.f"
	    gamma = b1 / a1;

#line 617 "ztgsja.f"
	    if (gamma < 0.) {
#line 618 "ztgsja.f"
		i__2 = *l - i__ + 1;
#line 618 "ztgsja.f"
		zdscal_(&i__2, &c_b39, &b[i__ + (*n - *l + i__) * b_dim1], 
			ldb);
#line 619 "ztgsja.f"
		if (wantv) {
#line 619 "ztgsja.f"
		    zdscal_(p, &c_b39, &v[i__ * v_dim1 + 1], &c__1);
#line 619 "ztgsja.f"
		}
#line 621 "ztgsja.f"
	    }

#line 623 "ztgsja.f"
	    d__1 = abs(gamma);
#line 623 "ztgsja.f"
	    dlartg_(&d__1, &c_b42, &beta[*k + i__], &alpha[*k + i__], &rwk);

#line 626 "ztgsja.f"
	    if (alpha[*k + i__] >= beta[*k + i__]) {
#line 627 "ztgsja.f"
		i__2 = *l - i__ + 1;
#line 627 "ztgsja.f"
		d__1 = 1. / alpha[*k + i__];
#line 627 "ztgsja.f"
		zdscal_(&i__2, &d__1, &a[*k + i__ + (*n - *l + i__) * a_dim1],
			 lda);
#line 629 "ztgsja.f"
	    } else {
#line 630 "ztgsja.f"
		i__2 = *l - i__ + 1;
#line 630 "ztgsja.f"
		d__1 = 1. / beta[*k + i__];
#line 630 "ztgsja.f"
		zdscal_(&i__2, &d__1, &b[i__ + (*n - *l + i__) * b_dim1], ldb)
			;
#line 632 "ztgsja.f"
		i__2 = *l - i__ + 1;
#line 632 "ztgsja.f"
		zcopy_(&i__2, &b[i__ + (*n - *l + i__) * b_dim1], ldb, &a[*k 
			+ i__ + (*n - *l + i__) * a_dim1], lda);
#line 634 "ztgsja.f"
	    }

#line 636 "ztgsja.f"
	} else {

#line 638 "ztgsja.f"
	    alpha[*k + i__] = 0.;
#line 639 "ztgsja.f"
	    beta[*k + i__] = 1.;
#line 640 "ztgsja.f"
	    i__2 = *l - i__ + 1;
#line 640 "ztgsja.f"
	    zcopy_(&i__2, &b[i__ + (*n - *l + i__) * b_dim1], ldb, &a[*k + 
		    i__ + (*n - *l + i__) * a_dim1], lda);
#line 642 "ztgsja.f"
	}
#line 643 "ztgsja.f"
/* L70: */
#line 643 "ztgsja.f"
    }

/*     Post-assignment */

#line 647 "ztgsja.f"
    i__1 = *k + *l;
#line 647 "ztgsja.f"
    for (i__ = *m + 1; i__ <= i__1; ++i__) {
#line 648 "ztgsja.f"
	alpha[i__] = 0.;
#line 649 "ztgsja.f"
	beta[i__] = 1.;
#line 650 "ztgsja.f"
/* L80: */
#line 650 "ztgsja.f"
    }

#line 652 "ztgsja.f"
    if (*k + *l < *n) {
#line 653 "ztgsja.f"
	i__1 = *n;
#line 653 "ztgsja.f"
	for (i__ = *k + *l + 1; i__ <= i__1; ++i__) {
#line 654 "ztgsja.f"
	    alpha[i__] = 0.;
#line 655 "ztgsja.f"
	    beta[i__] = 0.;
#line 656 "ztgsja.f"
/* L90: */
#line 656 "ztgsja.f"
	}
#line 657 "ztgsja.f"
    }

#line 659 "ztgsja.f"
L100:
#line 660 "ztgsja.f"
    *ncycle = kcycle;

#line 662 "ztgsja.f"
    return 0;

/*     End of ZTGSJA */

} /* ztgsja_ */

