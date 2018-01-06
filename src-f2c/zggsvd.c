#line 1 "zggsvd.f"
/* zggsvd.f -- translated by f2c (version 20100827).
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

#line 1 "zggsvd.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief <b> ZGGSVD computes the singular value decomposition (SVD) for OTHER matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZGGSVD + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zggsvd.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zggsvd.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zggsvd.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZGGSVD( JOBU, JOBV, JOBQ, M, N, P, K, L, A, LDA, B, */
/*                          LDB, ALPHA, BETA, U, LDU, V, LDV, Q, LDQ, WORK, */
/*                          RWORK, IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBQ, JOBU, JOBV */
/*       INTEGER            INFO, K, L, LDA, LDB, LDQ, LDU, LDV, M, N, P */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       DOUBLE PRECISION   ALPHA( * ), BETA( * ), RWORK( * ) */
/*       COMPLEX*16         A( LDA, * ), B( LDB, * ), Q( LDQ, * ), */
/*      $                   U( LDU, * ), V( LDV, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZGGSVD computes the generalized singular value decomposition (GSVD) */
/* > of an M-by-N complex matrix A and P-by-N complex matrix B: */
/* > */
/* >       U**H*A*Q = D1*( 0 R ),    V**H*B*Q = D2*( 0 R ) */
/* > */
/* > where U, V and Q are unitary matrices. */
/* > Let K+L = the effective numerical rank of the */
/* > matrix (A**H,B**H)**H, then R is a (K+L)-by-(K+L) nonsingular upper */
/* > triangular matrix, D1 and D2 are M-by-(K+L) and P-by-(K+L) "diagonal" */
/* > matrices and of the following structures, respectively: */
/* > */
/* > If M-K-L >= 0, */
/* > */
/* >                     K  L */
/* >        D1 =     K ( I  0 ) */
/* >                 L ( 0  C ) */
/* >             M-K-L ( 0  0 ) */
/* > */
/* >                   K  L */
/* >        D2 =   L ( 0  S ) */
/* >             P-L ( 0  0 ) */
/* > */
/* >                 N-K-L  K    L */
/* >   ( 0 R ) = K (  0   R11  R12 ) */
/* >             L (  0    0   R22 ) */
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
/* >                   K M-K K+L-M */
/* >        D1 =   K ( I  0    0   ) */
/* >             M-K ( 0  C    0   ) */
/* > */
/* >                     K M-K K+L-M */
/* >        D2 =   M-K ( 0  S    0  ) */
/* >             K+L-M ( 0  0    I  ) */
/* >               P-L ( 0  0    0  ) */
/* > */
/* >                    N-K-L  K   M-K  K+L-M */
/* >   ( 0 R ) =     K ( 0    R11  R12  R13  ) */
/* >               M-K ( 0     0   R22  R23  ) */
/* >             K+L-M ( 0     0    0   R33  ) */
/* > */
/* > where */
/* > */
/* >   C = diag( ALPHA(K+1), ... , ALPHA(M) ), */
/* >   S = diag( BETA(K+1),  ... , BETA(M) ), */
/* >   C**2 + S**2 = I. */
/* > */
/* >   (R11 R12 R13 ) is stored in A(1:M, N-K-L+1:N), and R33 is stored */
/* >   ( 0  R22 R23 ) */
/* >   in B(M-K+1:L,N+M-K-L+1:N) on exit. */
/* > */
/* > The routine computes C, S, R, and optionally the unitary */
/* > transformation matrices U, V and Q. */
/* > */
/* > In particular, if B is an N-by-N nonsingular matrix, then the GSVD of */
/* > A and B implicitly gives the SVD of A*inv(B): */
/* >                      A*inv(B) = U*(D1*inv(D2))*V**H. */
/* > If ( A**H,B**H)**H has orthnormal columns, then the GSVD of A and B is also */
/* > equal to the CS decomposition of A and B. Furthermore, the GSVD can */
/* > be used to derive the solution of the eigenvalue problem: */
/* >                      A**H*A x = lambda* B**H*B x. */
/* > In some literature, the GSVD of A and B is presented in the form */
/* >                  U**H*A*X = ( 0 D1 ),   V**H*B*X = ( 0 D2 ) */
/* > where U and V are orthogonal and X is nonsingular, and D1 and D2 are */
/* > ``diagonal''.  The former GSVD form can be converted to the latter */
/* > form by taking the nonsingular matrix X as */
/* > */
/* >                       X = Q*(  I   0    ) */
/* >                             (  0 inv(R) ) */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] JOBU */
/* > \verbatim */
/* >          JOBU is CHARACTER*1 */
/* >          = 'U':  Unitary matrix U is computed; */
/* >          = 'N':  U is not computed. */
/* > \endverbatim */
/* > */
/* > \param[in] JOBV */
/* > \verbatim */
/* >          JOBV is CHARACTER*1 */
/* >          = 'V':  Unitary matrix V is computed; */
/* >          = 'N':  V is not computed. */
/* > \endverbatim */
/* > */
/* > \param[in] JOBQ */
/* > \verbatim */
/* >          JOBQ is CHARACTER*1 */
/* >          = 'Q':  Unitary matrix Q is computed; */
/* >          = 'N':  Q is not computed. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the matrix A.  M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrices A and B.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] P */
/* > \verbatim */
/* >          P is INTEGER */
/* >          The number of rows of the matrix B.  P >= 0. */
/* > \endverbatim */
/* > */
/* > \param[out] K */
/* > \verbatim */
/* >          K is INTEGER */
/* > \endverbatim */
/* > */
/* > \param[out] L */
/* > \verbatim */
/* >          L is INTEGER */
/* > */
/* >          On exit, K and L specify the dimension of the subblocks */
/* >          described in Purpose. */
/* >          K + L = effective numerical rank of (A**H,B**H)**H. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
/* >          On entry, the M-by-N matrix A. */
/* >          On exit, A contains the triangular matrix R, or part of R. */
/* >          See Purpose for details. */
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
/* >          On exit, B contains part of the triangular matrix R if */
/* >          M-K-L < 0.  See Purpose for details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >          The leading dimension of the array B. LDB >= max(1,P). */
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
/* >            ALPHA(K+1:K+L) = C, */
/* >            BETA(K+1:K+L)  = S, */
/* >          or if M-K-L < 0, */
/* >            ALPHA(K+1:M)=C, ALPHA(M+1:K+L)=0 */
/* >            BETA(K+1:M) =S, BETA(M+1:K+L) =1 */
/* >          and */
/* >            ALPHA(K+L+1:N) = 0 */
/* >            BETA(K+L+1:N)  = 0 */
/* > \endverbatim */
/* > */
/* > \param[out] U */
/* > \verbatim */
/* >          U is COMPLEX*16 array, dimension (LDU,M) */
/* >          If JOBU = 'U', U contains the M-by-M unitary matrix U. */
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
/* > \param[out] V */
/* > \verbatim */
/* >          V is COMPLEX*16 array, dimension (LDV,P) */
/* >          If JOBV = 'V', V contains the P-by-P unitary matrix V. */
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
/* > \param[out] Q */
/* > \verbatim */
/* >          Q is COMPLEX*16 array, dimension (LDQ,N) */
/* >          If JOBQ = 'Q', Q contains the N-by-N unitary matrix Q. */
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
/* >          WORK is COMPLEX*16 array, dimension (max(3*N,M,P)+N) */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is DOUBLE PRECISION array, dimension (2*N) */
/* > \endverbatim */
/* > */
/* > \param[out] IWORK */
/* > \verbatim */
/* >          IWORK is INTEGER array, dimension (N) */
/* >          On exit, IWORK stores the sorting information. More */
/* >          precisely, the following loop will sort ALPHA */
/* >             for I = K+1, min(M,K+L) */
/* >                 swap ALPHA(I) and ALPHA(IWORK(I)) */
/* >             endfor */
/* >          such that ALPHA(1) >= ALPHA(2) >= ... >= ALPHA(N). */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit. */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value. */
/* >          > 0:  if INFO = 1, the Jacobi-type procedure failed to */
/* >                converge.  For further details, see subroutine ZTGSJA. */
/* > \endverbatim */

/* > \par Internal Parameters: */
/*  ========================= */
/* > */
/* > \verbatim */
/* >  TOLA    DOUBLE PRECISION */
/* >  TOLB    DOUBLE PRECISION */
/* >          TOLA and TOLB are the thresholds to determine the effective */
/* >          rank of (A**H,B**H)**H. Generally, they are set to */
/* >                   TOLA = MAX(M,N)*norm(A)*MAZHEPS, */
/* >                   TOLB = MAX(P,N)*norm(B)*MAZHEPS. */
/* >          The size of TOLA and TOLB may affect the size of backward */
/* >          errors of the decomposition. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup complex16OTHERsing */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Ming Gu and Huan Ren, Computer Science Division, University of */
/* >     California at Berkeley, USA */
/* > */
/*  ===================================================================== */
/* Subroutine */ int zggsvd_(char *jobu, char *jobv, char *jobq, integer *m, 
	integer *n, integer *p, integer *k, integer *l, doublecomplex *a, 
	integer *lda, doublecomplex *b, integer *ldb, doublereal *alpha, 
	doublereal *beta, doublecomplex *u, integer *ldu, doublecomplex *v, 
	integer *ldv, doublecomplex *q, integer *ldq, doublecomplex *work, 
	doublereal *rwork, integer *iwork, integer *info, ftnlen jobu_len, 
	ftnlen jobv_len, ftnlen jobq_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, q_dim1, q_offset, u_dim1, 
	    u_offset, v_dim1, v_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j;
    static doublereal ulp;
    static integer ibnd;
    static doublereal tola;
    static integer isub;
    static doublereal tolb, unfl, temp, smax;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal anorm, bnorm;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static logical wantq, wantu, wantv;
    extern doublereal dlamch_(char *, ftnlen);
    static integer ncycle;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern doublereal zlange_(char *, integer *, integer *, doublecomplex *, 
	    integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int ztgsja_(char *, char *, char *, integer *, 
	    integer *, integer *, integer *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, integer *, ftnlen, ftnlen, ftnlen), 
	    zggsvp_(char *, char *, char *, integer *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublecomplex *,
	     integer *, doublecomplex *, integer *, doublecomplex *, integer *
	    , integer *, doublereal *, doublecomplex *, doublecomplex *, 
	    integer *, ftnlen, ftnlen, ftnlen);


/*  -- LAPACK driver routine (version 3.4.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

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

#line 376 "zggsvd.f"
    /* Parameter adjustments */
#line 376 "zggsvd.f"
    a_dim1 = *lda;
#line 376 "zggsvd.f"
    a_offset = 1 + a_dim1;
#line 376 "zggsvd.f"
    a -= a_offset;
#line 376 "zggsvd.f"
    b_dim1 = *ldb;
#line 376 "zggsvd.f"
    b_offset = 1 + b_dim1;
#line 376 "zggsvd.f"
    b -= b_offset;
#line 376 "zggsvd.f"
    --alpha;
#line 376 "zggsvd.f"
    --beta;
#line 376 "zggsvd.f"
    u_dim1 = *ldu;
#line 376 "zggsvd.f"
    u_offset = 1 + u_dim1;
#line 376 "zggsvd.f"
    u -= u_offset;
#line 376 "zggsvd.f"
    v_dim1 = *ldv;
#line 376 "zggsvd.f"
    v_offset = 1 + v_dim1;
#line 376 "zggsvd.f"
    v -= v_offset;
#line 376 "zggsvd.f"
    q_dim1 = *ldq;
#line 376 "zggsvd.f"
    q_offset = 1 + q_dim1;
#line 376 "zggsvd.f"
    q -= q_offset;
#line 376 "zggsvd.f"
    --work;
#line 376 "zggsvd.f"
    --rwork;
#line 376 "zggsvd.f"
    --iwork;
#line 376 "zggsvd.f"

#line 376 "zggsvd.f"
    /* Function Body */
#line 376 "zggsvd.f"
    wantu = lsame_(jobu, "U", (ftnlen)1, (ftnlen)1);
#line 377 "zggsvd.f"
    wantv = lsame_(jobv, "V", (ftnlen)1, (ftnlen)1);
#line 378 "zggsvd.f"
    wantq = lsame_(jobq, "Q", (ftnlen)1, (ftnlen)1);

#line 380 "zggsvd.f"
    *info = 0;
#line 381 "zggsvd.f"
    if (! (wantu || lsame_(jobu, "N", (ftnlen)1, (ftnlen)1))) {
#line 382 "zggsvd.f"
	*info = -1;
#line 383 "zggsvd.f"
    } else if (! (wantv || lsame_(jobv, "N", (ftnlen)1, (ftnlen)1))) {
#line 384 "zggsvd.f"
	*info = -2;
#line 385 "zggsvd.f"
    } else if (! (wantq || lsame_(jobq, "N", (ftnlen)1, (ftnlen)1))) {
#line 386 "zggsvd.f"
	*info = -3;
#line 387 "zggsvd.f"
    } else if (*m < 0) {
#line 388 "zggsvd.f"
	*info = -4;
#line 389 "zggsvd.f"
    } else if (*n < 0) {
#line 390 "zggsvd.f"
	*info = -5;
#line 391 "zggsvd.f"
    } else if (*p < 0) {
#line 392 "zggsvd.f"
	*info = -6;
#line 393 "zggsvd.f"
    } else if (*lda < max(1,*m)) {
#line 394 "zggsvd.f"
	*info = -10;
#line 395 "zggsvd.f"
    } else if (*ldb < max(1,*p)) {
#line 396 "zggsvd.f"
	*info = -12;
#line 397 "zggsvd.f"
    } else if (*ldu < 1 || wantu && *ldu < *m) {
#line 398 "zggsvd.f"
	*info = -16;
#line 399 "zggsvd.f"
    } else if (*ldv < 1 || wantv && *ldv < *p) {
#line 400 "zggsvd.f"
	*info = -18;
#line 401 "zggsvd.f"
    } else if (*ldq < 1 || wantq && *ldq < *n) {
#line 402 "zggsvd.f"
	*info = -20;
#line 403 "zggsvd.f"
    }
#line 404 "zggsvd.f"
    if (*info != 0) {
#line 405 "zggsvd.f"
	i__1 = -(*info);
#line 405 "zggsvd.f"
	xerbla_("ZGGSVD", &i__1, (ftnlen)6);
#line 406 "zggsvd.f"
	return 0;
#line 407 "zggsvd.f"
    }

/*     Compute the Frobenius norm of matrices A and B */

#line 411 "zggsvd.f"
    anorm = zlange_("1", m, n, &a[a_offset], lda, &rwork[1], (ftnlen)1);
#line 412 "zggsvd.f"
    bnorm = zlange_("1", p, n, &b[b_offset], ldb, &rwork[1], (ftnlen)1);

/*     Get machine precision and set up threshold for determining */
/*     the effective numerical rank of the matrices A and B. */

#line 417 "zggsvd.f"
    ulp = dlamch_("Precision", (ftnlen)9);
#line 418 "zggsvd.f"
    unfl = dlamch_("Safe Minimum", (ftnlen)12);
#line 419 "zggsvd.f"
    tola = max(*m,*n) * max(anorm,unfl) * ulp;
#line 420 "zggsvd.f"
    tolb = max(*p,*n) * max(bnorm,unfl) * ulp;

#line 422 "zggsvd.f"
    zggsvp_(jobu, jobv, jobq, m, p, n, &a[a_offset], lda, &b[b_offset], ldb, &
	    tola, &tolb, k, l, &u[u_offset], ldu, &v[v_offset], ldv, &q[
	    q_offset], ldq, &iwork[1], &rwork[1], &work[1], &work[*n + 1], 
	    info, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*     Compute the GSVD of two upper "triangular" matrices */

#line 428 "zggsvd.f"
    ztgsja_(jobu, jobv, jobq, m, p, n, k, l, &a[a_offset], lda, &b[b_offset], 
	    ldb, &tola, &tolb, &alpha[1], &beta[1], &u[u_offset], ldu, &v[
	    v_offset], ldv, &q[q_offset], ldq, &work[1], &ncycle, info, (
	    ftnlen)1, (ftnlen)1, (ftnlen)1);

/*     Sort the singular values and store the pivot indices in IWORK */
/*     Copy ALPHA to RWORK, then sort ALPHA in RWORK */

#line 435 "zggsvd.f"
    dcopy_(n, &alpha[1], &c__1, &rwork[1], &c__1);
/* Computing MIN */
#line 436 "zggsvd.f"
    i__1 = *l, i__2 = *m - *k;
#line 436 "zggsvd.f"
    ibnd = min(i__1,i__2);
#line 437 "zggsvd.f"
    i__1 = ibnd;
#line 437 "zggsvd.f"
    for (i__ = 1; i__ <= i__1; ++i__) {

/*        Scan for largest ALPHA(K+I) */

#line 441 "zggsvd.f"
	isub = i__;
#line 442 "zggsvd.f"
	smax = rwork[*k + i__];
#line 443 "zggsvd.f"
	i__2 = ibnd;
#line 443 "zggsvd.f"
	for (j = i__ + 1; j <= i__2; ++j) {
#line 444 "zggsvd.f"
	    temp = rwork[*k + j];
#line 445 "zggsvd.f"
	    if (temp > smax) {
#line 446 "zggsvd.f"
		isub = j;
#line 447 "zggsvd.f"
		smax = temp;
#line 448 "zggsvd.f"
	    }
#line 449 "zggsvd.f"
/* L10: */
#line 449 "zggsvd.f"
	}
#line 450 "zggsvd.f"
	if (isub != i__) {
#line 451 "zggsvd.f"
	    rwork[*k + isub] = rwork[*k + i__];
#line 452 "zggsvd.f"
	    rwork[*k + i__] = smax;
#line 453 "zggsvd.f"
	    iwork[*k + i__] = *k + isub;
#line 454 "zggsvd.f"
	} else {
#line 455 "zggsvd.f"
	    iwork[*k + i__] = *k + i__;
#line 456 "zggsvd.f"
	}
#line 457 "zggsvd.f"
/* L20: */
#line 457 "zggsvd.f"
    }

#line 459 "zggsvd.f"
    return 0;

/*     End of ZGGSVD */

} /* zggsvd_ */

