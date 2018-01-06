#line 1 "cggsvd3.f"
/* cggsvd3.f -- translated by f2c (version 20100827).
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

#line 1 "cggsvd3.f"
/* Table of constant values */

static integer c_n1 = -1;
static integer c__1 = 1;

/* > \brief <b> CGGSVD3 computes the singular value decomposition (SVD) for OTHER matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CGGSVD3 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cggsvd3
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cggsvd3
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cggsvd3
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CGGSVD3( JOBU, JOBV, JOBQ, M, N, P, K, L, A, LDA, B, */
/*                           LDB, ALPHA, BETA, U, LDU, V, LDV, Q, LDQ, WORK, */
/*                           LWORK, RWORK, IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBQ, JOBU, JOBV */
/*       INTEGER            INFO, K, L, LDA, LDB, LDQ, LDU, LDV, M, N, P, LWORK */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       REAL               ALPHA( * ), BETA( * ), RWORK( * ) */
/*       COMPLEX            A( LDA, * ), B( LDB, * ), Q( LDQ, * ), */
/*      $                   U( LDU, * ), V( LDV, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CGGSVD3 computes the generalized singular value decomposition (GSVD) */
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
/* > If ( A**H,B**H)**H has orthonormal columns, then the GSVD of A and B is also */
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
/* >          A is COMPLEX array, dimension (LDA,N) */
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
/* >          B is COMPLEX array, dimension (LDB,N) */
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
/* >          U is COMPLEX array, dimension (LDU,M) */
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
/* >          V is COMPLEX array, dimension (LDV,P) */
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
/* >          Q is COMPLEX array, dimension (LDQ,N) */
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
/* >          WORK is COMPLEX array, dimension (MAX(1,LWORK)) */
/* >          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */
/* > \endverbatim */
/* > */
/* > \param[in] LWORK */
/* > \verbatim */
/* >          LWORK is INTEGER */
/* >          The dimension of the array WORK. */
/* > */
/* >          If LWORK = -1, then a workspace query is assumed; the routine */
/* >          only calculates the optimal size of the WORK array, returns */
/* >          this value as the first entry of the WORK array, and no error */
/* >          message related to LWORK is issued by XERBLA. */
/* > \endverbatim */
/* > */
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is REAL array, dimension (2*N) */
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
/* >                converge.  For further details, see subroutine CTGSJA. */
/* > \endverbatim */

/* > \par Internal Parameters: */
/*  ========================= */
/* > */
/* > \verbatim */
/* >  TOLA    REAL */
/* >  TOLB    REAL */
/* >          TOLA and TOLB are the thresholds to determine the effective */
/* >          rank of (A**H,B**H)**H. Generally, they are set to */
/* >                   TOLA = MAX(M,N)*norm(A)*MACHEPS, */
/* >                   TOLB = MAX(P,N)*norm(B)*MACHEPS. */
/* >          The size of TOLA and TOLB may affect the size of backward */
/* >          errors of the decomposition. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date August 2015 */

/* > \ingroup complexOTHERsing */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Ming Gu and Huan Ren, Computer Science Division, University of */
/* >     California at Berkeley, USA */
/* > */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* >  CGGSVD3 replaces the deprecated subroutine CGGSVD. */
/* > */
/*  ===================================================================== */
/* Subroutine */ int cggsvd3_(char *jobu, char *jobv, char *jobq, integer *m, 
	integer *n, integer *p, integer *k, integer *l, doublecomplex *a, 
	integer *lda, doublecomplex *b, integer *ldb, doublereal *alpha, 
	doublereal *beta, doublecomplex *u, integer *ldu, doublecomplex *v, 
	integer *ldv, doublecomplex *q, integer *ldq, doublecomplex *work, 
	integer *lwork, doublereal *rwork, integer *iwork, integer *info, 
	ftnlen jobu_len, ftnlen jobv_len, ftnlen jobq_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, q_dim1, q_offset, u_dim1, 
	    u_offset, v_dim1, v_offset, i__1, i__2;
    doublecomplex z__1;

    /* Local variables */
    static integer i__, j;
    static doublereal ulp;
    static integer ibnd;
    static doublereal tola;
    static integer isub;
    static doublereal tolb, unfl, temp, smax;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal anorm, bnorm;
    static logical wantq;
    extern /* Subroutine */ int scopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static logical wantu, wantv;
    extern doublereal clange_(char *, integer *, integer *, doublecomplex *, 
	    integer *, doublereal *, ftnlen), slamch_(char *, ftnlen);
    extern /* Subroutine */ int ctgsja_(char *, char *, char *, integer *, 
	    integer *, integer *, integer *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, integer *, ftnlen, ftnlen, ftnlen);
    static integer ncycle;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static integer lwkopt;
    static logical lquery;
    extern /* Subroutine */ int cggsvp3_(char *, char *, char *, integer *, 
	    integer *, integer *, doublecomplex *, integer *, doublecomplex *,
	     integer *, doublereal *, doublereal *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, integer *, doublereal *, 
	    doublecomplex *, doublecomplex *, integer *, integer *, ftnlen, 
	    ftnlen, ftnlen);


/*  -- LAPACK driver routine (version 3.6.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     August 2015 */

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

#line 396 "cggsvd3.f"
    /* Parameter adjustments */
#line 396 "cggsvd3.f"
    a_dim1 = *lda;
#line 396 "cggsvd3.f"
    a_offset = 1 + a_dim1;
#line 396 "cggsvd3.f"
    a -= a_offset;
#line 396 "cggsvd3.f"
    b_dim1 = *ldb;
#line 396 "cggsvd3.f"
    b_offset = 1 + b_dim1;
#line 396 "cggsvd3.f"
    b -= b_offset;
#line 396 "cggsvd3.f"
    --alpha;
#line 396 "cggsvd3.f"
    --beta;
#line 396 "cggsvd3.f"
    u_dim1 = *ldu;
#line 396 "cggsvd3.f"
    u_offset = 1 + u_dim1;
#line 396 "cggsvd3.f"
    u -= u_offset;
#line 396 "cggsvd3.f"
    v_dim1 = *ldv;
#line 396 "cggsvd3.f"
    v_offset = 1 + v_dim1;
#line 396 "cggsvd3.f"
    v -= v_offset;
#line 396 "cggsvd3.f"
    q_dim1 = *ldq;
#line 396 "cggsvd3.f"
    q_offset = 1 + q_dim1;
#line 396 "cggsvd3.f"
    q -= q_offset;
#line 396 "cggsvd3.f"
    --work;
#line 396 "cggsvd3.f"
    --rwork;
#line 396 "cggsvd3.f"
    --iwork;
#line 396 "cggsvd3.f"

#line 396 "cggsvd3.f"
    /* Function Body */
#line 396 "cggsvd3.f"
    wantu = lsame_(jobu, "U", (ftnlen)1, (ftnlen)1);
#line 397 "cggsvd3.f"
    wantv = lsame_(jobv, "V", (ftnlen)1, (ftnlen)1);
#line 398 "cggsvd3.f"
    wantq = lsame_(jobq, "Q", (ftnlen)1, (ftnlen)1);
#line 399 "cggsvd3.f"
    lquery = *lwork == -1;
#line 400 "cggsvd3.f"
    lwkopt = 1;

/*     Test the input arguments */

#line 404 "cggsvd3.f"
    *info = 0;
#line 405 "cggsvd3.f"
    if (! (wantu || lsame_(jobu, "N", (ftnlen)1, (ftnlen)1))) {
#line 406 "cggsvd3.f"
	*info = -1;
#line 407 "cggsvd3.f"
    } else if (! (wantv || lsame_(jobv, "N", (ftnlen)1, (ftnlen)1))) {
#line 408 "cggsvd3.f"
	*info = -2;
#line 409 "cggsvd3.f"
    } else if (! (wantq || lsame_(jobq, "N", (ftnlen)1, (ftnlen)1))) {
#line 410 "cggsvd3.f"
	*info = -3;
#line 411 "cggsvd3.f"
    } else if (*m < 0) {
#line 412 "cggsvd3.f"
	*info = -4;
#line 413 "cggsvd3.f"
    } else if (*n < 0) {
#line 414 "cggsvd3.f"
	*info = -5;
#line 415 "cggsvd3.f"
    } else if (*p < 0) {
#line 416 "cggsvd3.f"
	*info = -6;
#line 417 "cggsvd3.f"
    } else if (*lda < max(1,*m)) {
#line 418 "cggsvd3.f"
	*info = -10;
#line 419 "cggsvd3.f"
    } else if (*ldb < max(1,*p)) {
#line 420 "cggsvd3.f"
	*info = -12;
#line 421 "cggsvd3.f"
    } else if (*ldu < 1 || wantu && *ldu < *m) {
#line 422 "cggsvd3.f"
	*info = -16;
#line 423 "cggsvd3.f"
    } else if (*ldv < 1 || wantv && *ldv < *p) {
#line 424 "cggsvd3.f"
	*info = -18;
#line 425 "cggsvd3.f"
    } else if (*ldq < 1 || wantq && *ldq < *n) {
#line 426 "cggsvd3.f"
	*info = -20;
#line 427 "cggsvd3.f"
    } else if (*lwork < 1 && ! lquery) {
#line 428 "cggsvd3.f"
	*info = -24;
#line 429 "cggsvd3.f"
    }

/*     Compute workspace */

#line 433 "cggsvd3.f"
    if (*info == 0) {
#line 434 "cggsvd3.f"
	cggsvp3_(jobu, jobv, jobq, m, p, n, &a[a_offset], lda, &b[b_offset], 
		ldb, &tola, &tolb, k, l, &u[u_offset], ldu, &v[v_offset], ldv,
		 &q[q_offset], ldq, &iwork[1], &rwork[1], &work[1], &work[1], 
		&c_n1, info, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 437 "cggsvd3.f"
	lwkopt = *n + (integer) work[1].r;
/* Computing MAX */
#line 438 "cggsvd3.f"
	i__1 = *n << 1;
#line 438 "cggsvd3.f"
	lwkopt = max(i__1,lwkopt);
#line 439 "cggsvd3.f"
	lwkopt = max(1,lwkopt);
#line 440 "cggsvd3.f"
	z__1.r = (doublereal) lwkopt, z__1.i = 0.;
#line 440 "cggsvd3.f"
	work[1].r = z__1.r, work[1].i = z__1.i;
#line 441 "cggsvd3.f"
    }

#line 443 "cggsvd3.f"
    if (*info != 0) {
#line 444 "cggsvd3.f"
	i__1 = -(*info);
#line 444 "cggsvd3.f"
	xerbla_("CGGSVD3", &i__1, (ftnlen)7);
#line 445 "cggsvd3.f"
	return 0;
#line 446 "cggsvd3.f"
    }
#line 447 "cggsvd3.f"
    if (lquery) {
#line 448 "cggsvd3.f"
	return 0;
#line 449 "cggsvd3.f"
    }

/*     Compute the Frobenius norm of matrices A and B */

#line 453 "cggsvd3.f"
    anorm = clange_("1", m, n, &a[a_offset], lda, &rwork[1], (ftnlen)1);
#line 454 "cggsvd3.f"
    bnorm = clange_("1", p, n, &b[b_offset], ldb, &rwork[1], (ftnlen)1);

/*     Get machine precision and set up threshold for determining */
/*     the effective numerical rank of the matrices A and B. */

#line 459 "cggsvd3.f"
    ulp = slamch_("Precision", (ftnlen)9);
#line 460 "cggsvd3.f"
    unfl = slamch_("Safe Minimum", (ftnlen)12);
#line 461 "cggsvd3.f"
    tola = max(*m,*n) * max(anorm,unfl) * ulp;
#line 462 "cggsvd3.f"
    tolb = max(*p,*n) * max(bnorm,unfl) * ulp;

#line 464 "cggsvd3.f"
    i__1 = *lwork - *n;
#line 464 "cggsvd3.f"
    cggsvp3_(jobu, jobv, jobq, m, p, n, &a[a_offset], lda, &b[b_offset], ldb, 
	    &tola, &tolb, k, l, &u[u_offset], ldu, &v[v_offset], ldv, &q[
	    q_offset], ldq, &iwork[1], &rwork[1], &work[1], &work[*n + 1], &
	    i__1, info, (ftnlen)1, (ftnlen)1, (ftnlen)1);

/*     Compute the GSVD of two upper "triangular" matrices */

#line 470 "cggsvd3.f"
    ctgsja_(jobu, jobv, jobq, m, p, n, k, l, &a[a_offset], lda, &b[b_offset], 
	    ldb, &tola, &tolb, &alpha[1], &beta[1], &u[u_offset], ldu, &v[
	    v_offset], ldv, &q[q_offset], ldq, &work[1], &ncycle, info, (
	    ftnlen)1, (ftnlen)1, (ftnlen)1);

/*     Sort the singular values and store the pivot indices in IWORK */
/*     Copy ALPHA to RWORK, then sort ALPHA in RWORK */

#line 477 "cggsvd3.f"
    scopy_(n, &alpha[1], &c__1, &rwork[1], &c__1);
/* Computing MIN */
#line 478 "cggsvd3.f"
    i__1 = *l, i__2 = *m - *k;
#line 478 "cggsvd3.f"
    ibnd = min(i__1,i__2);
#line 479 "cggsvd3.f"
    i__1 = ibnd;
#line 479 "cggsvd3.f"
    for (i__ = 1; i__ <= i__1; ++i__) {

/*        Scan for largest ALPHA(K+I) */

#line 483 "cggsvd3.f"
	isub = i__;
#line 484 "cggsvd3.f"
	smax = rwork[*k + i__];
#line 485 "cggsvd3.f"
	i__2 = ibnd;
#line 485 "cggsvd3.f"
	for (j = i__ + 1; j <= i__2; ++j) {
#line 486 "cggsvd3.f"
	    temp = rwork[*k + j];
#line 487 "cggsvd3.f"
	    if (temp > smax) {
#line 488 "cggsvd3.f"
		isub = j;
#line 489 "cggsvd3.f"
		smax = temp;
#line 490 "cggsvd3.f"
	    }
#line 491 "cggsvd3.f"
/* L10: */
#line 491 "cggsvd3.f"
	}
#line 492 "cggsvd3.f"
	if (isub != i__) {
#line 493 "cggsvd3.f"
	    rwork[*k + isub] = rwork[*k + i__];
#line 494 "cggsvd3.f"
	    rwork[*k + i__] = smax;
#line 495 "cggsvd3.f"
	    iwork[*k + i__] = *k + isub;
#line 496 "cggsvd3.f"
	} else {
#line 497 "cggsvd3.f"
	    iwork[*k + i__] = *k + i__;
#line 498 "cggsvd3.f"
	}
#line 499 "cggsvd3.f"
/* L20: */
#line 499 "cggsvd3.f"
    }

#line 501 "cggsvd3.f"
    z__1.r = (doublereal) lwkopt, z__1.i = 0.;
#line 501 "cggsvd3.f"
    work[1].r = z__1.r, work[1].i = z__1.i;
#line 502 "cggsvd3.f"
    return 0;

/*     End of CGGSVD3 */

} /* cggsvd3_ */

