#line 1 "dggsvd3.f"
/* dggsvd3.f -- translated by f2c (version 20100827).
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

#line 1 "dggsvd3.f"
/* Table of constant values */

static integer c_n1 = -1;
static integer c__1 = 1;

/* > \brief <b> DGGSVD3 computes the singular value decomposition (SVD) for OTHER matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DGGSVD3 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dggsvd3
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dggsvd3
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dggsvd3
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DGGSVD3( JOBU, JOBV, JOBQ, M, N, P, K, L, A, LDA, B, */
/*                           LDB, ALPHA, BETA, U, LDU, V, LDV, Q, LDQ, WORK, */
/*                           LWORK, IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBQ, JOBU, JOBV */
/*       INTEGER            INFO, K, L, LDA, LDB, LDQ, LDU, LDV, M, N, P, LWORK */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       DOUBLE PRECISION   A( LDA, * ), ALPHA( * ), B( LDB, * ), */
/*      $                   BETA( * ), Q( LDQ, * ), U( LDU, * ), */
/*      $                   V( LDV, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DGGSVD3 computes the generalized singular value decomposition (GSVD) */
/* > of an M-by-N real matrix A and P-by-N real matrix B: */
/* > */
/* >       U**T*A*Q = D1*( 0 R ),    V**T*B*Q = D2*( 0 R ) */
/* > */
/* > where U, V and Q are orthogonal matrices. */
/* > Let K+L = the effective numerical rank of the matrix (A**T,B**T)**T, */
/* > then R is a K+L-by-K+L nonsingular upper triangular matrix, D1 and */
/* > D2 are M-by-(K+L) and P-by-(K+L) "diagonal" matrices and of the */
/* > following structures, respectively: */
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
/* > The routine computes C, S, R, and optionally the orthogonal */
/* > transformation matrices U, V and Q. */
/* > */
/* > In particular, if B is an N-by-N nonsingular matrix, then the GSVD of */
/* > A and B implicitly gives the SVD of A*inv(B): */
/* >                      A*inv(B) = U*(D1*inv(D2))*V**T. */
/* > If ( A**T,B**T)**T  has orthonormal columns, then the GSVD of A and B is */
/* > also equal to the CS decomposition of A and B. Furthermore, the GSVD */
/* > can be used to derive the solution of the eigenvalue problem: */
/* >                      A**T*A x = lambda* B**T*B x. */
/* > In some literature, the GSVD of A and B is presented in the form */
/* >                  U**T*A*X = ( 0 D1 ),   V**T*B*X = ( 0 D2 ) */
/* > where U and V are orthogonal and X is nonsingular, D1 and D2 are */
/* > ``diagonal''.  The former GSVD form can be converted to the latter */
/* > form by taking the nonsingular matrix X as */
/* > */
/* >                      X = Q*( I   0    ) */
/* >                            ( 0 inv(R) ). */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] JOBU */
/* > \verbatim */
/* >          JOBU is CHARACTER*1 */
/* >          = 'U':  Orthogonal matrix U is computed; */
/* >          = 'N':  U is not computed. */
/* > \endverbatim */
/* > */
/* > \param[in] JOBV */
/* > \verbatim */
/* >          JOBV is CHARACTER*1 */
/* >          = 'V':  Orthogonal matrix V is computed; */
/* >          = 'N':  V is not computed. */
/* > \endverbatim */
/* > */
/* > \param[in] JOBQ */
/* > \verbatim */
/* >          JOBQ is CHARACTER*1 */
/* >          = 'Q':  Orthogonal matrix Q is computed; */
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
/* >          K + L = effective numerical rank of (A**T,B**T)**T. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
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
/* >          B is DOUBLE PRECISION array, dimension (LDB,N) */
/* >          On entry, the P-by-N matrix B. */
/* >          On exit, B contains the triangular matrix R if M-K-L < 0. */
/* >          See Purpose for details. */
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
/* >          U is DOUBLE PRECISION array, dimension (LDU,M) */
/* >          If JOBU = 'U', U contains the M-by-M orthogonal matrix U. */
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
/* >          V is DOUBLE PRECISION array, dimension (LDV,P) */
/* >          If JOBV = 'V', V contains the P-by-P orthogonal matrix V. */
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
/* >          Q is DOUBLE PRECISION array, dimension (LDQ,N) */
/* >          If JOBQ = 'Q', Q contains the N-by-N orthogonal matrix Q. */
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
/* >          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)) */
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
/* >                converge.  For further details, see subroutine DTGSJA. */
/* > \endverbatim */

/* > \par Internal Parameters: */
/*  ========================= */
/* > */
/* > \verbatim */
/* >  TOLA    DOUBLE PRECISION */
/* >  TOLB    DOUBLE PRECISION */
/* >          TOLA and TOLB are the thresholds to determine the effective */
/* >          rank of (A**T,B**T)**T. Generally, they are set to */
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

/* > \ingroup doubleOTHERsing */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Ming Gu and Huan Ren, Computer Science Division, University of */
/* >     California at Berkeley, USA */
/* > */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* >  DGGSVD3 replaces the deprecated subroutine DGGSVD. */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dggsvd3_(char *jobu, char *jobv, char *jobq, integer *m, 
	integer *n, integer *p, integer *k, integer *l, doublereal *a, 
	integer *lda, doublereal *b, integer *ldb, doublereal *alpha, 
	doublereal *beta, doublereal *u, integer *ldu, doublereal *v, integer 
	*ldv, doublereal *q, integer *ldq, doublereal *work, integer *lwork, 
	integer *iwork, integer *info, ftnlen jobu_len, ftnlen jobv_len, 
	ftnlen jobq_len)
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
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int dtgsja_(char *, char *, char *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, integer 
	    *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);
    static integer ncycle;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static integer lwkopt;
    static logical lquery;
    extern /* Subroutine */ int dggsvp3_(char *, char *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);


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

#line 391 "dggsvd3.f"
    /* Parameter adjustments */
#line 391 "dggsvd3.f"
    a_dim1 = *lda;
#line 391 "dggsvd3.f"
    a_offset = 1 + a_dim1;
#line 391 "dggsvd3.f"
    a -= a_offset;
#line 391 "dggsvd3.f"
    b_dim1 = *ldb;
#line 391 "dggsvd3.f"
    b_offset = 1 + b_dim1;
#line 391 "dggsvd3.f"
    b -= b_offset;
#line 391 "dggsvd3.f"
    --alpha;
#line 391 "dggsvd3.f"
    --beta;
#line 391 "dggsvd3.f"
    u_dim1 = *ldu;
#line 391 "dggsvd3.f"
    u_offset = 1 + u_dim1;
#line 391 "dggsvd3.f"
    u -= u_offset;
#line 391 "dggsvd3.f"
    v_dim1 = *ldv;
#line 391 "dggsvd3.f"
    v_offset = 1 + v_dim1;
#line 391 "dggsvd3.f"
    v -= v_offset;
#line 391 "dggsvd3.f"
    q_dim1 = *ldq;
#line 391 "dggsvd3.f"
    q_offset = 1 + q_dim1;
#line 391 "dggsvd3.f"
    q -= q_offset;
#line 391 "dggsvd3.f"
    --work;
#line 391 "dggsvd3.f"
    --iwork;
#line 391 "dggsvd3.f"

#line 391 "dggsvd3.f"
    /* Function Body */
#line 391 "dggsvd3.f"
    wantu = lsame_(jobu, "U", (ftnlen)1, (ftnlen)1);
#line 392 "dggsvd3.f"
    wantv = lsame_(jobv, "V", (ftnlen)1, (ftnlen)1);
#line 393 "dggsvd3.f"
    wantq = lsame_(jobq, "Q", (ftnlen)1, (ftnlen)1);
#line 394 "dggsvd3.f"
    lquery = *lwork == -1;
#line 395 "dggsvd3.f"
    lwkopt = 1;

/*     Test the input arguments */

#line 399 "dggsvd3.f"
    *info = 0;
#line 400 "dggsvd3.f"
    if (! (wantu || lsame_(jobu, "N", (ftnlen)1, (ftnlen)1))) {
#line 401 "dggsvd3.f"
	*info = -1;
#line 402 "dggsvd3.f"
    } else if (! (wantv || lsame_(jobv, "N", (ftnlen)1, (ftnlen)1))) {
#line 403 "dggsvd3.f"
	*info = -2;
#line 404 "dggsvd3.f"
    } else if (! (wantq || lsame_(jobq, "N", (ftnlen)1, (ftnlen)1))) {
#line 405 "dggsvd3.f"
	*info = -3;
#line 406 "dggsvd3.f"
    } else if (*m < 0) {
#line 407 "dggsvd3.f"
	*info = -4;
#line 408 "dggsvd3.f"
    } else if (*n < 0) {
#line 409 "dggsvd3.f"
	*info = -5;
#line 410 "dggsvd3.f"
    } else if (*p < 0) {
#line 411 "dggsvd3.f"
	*info = -6;
#line 412 "dggsvd3.f"
    } else if (*lda < max(1,*m)) {
#line 413 "dggsvd3.f"
	*info = -10;
#line 414 "dggsvd3.f"
    } else if (*ldb < max(1,*p)) {
#line 415 "dggsvd3.f"
	*info = -12;
#line 416 "dggsvd3.f"
    } else if (*ldu < 1 || wantu && *ldu < *m) {
#line 417 "dggsvd3.f"
	*info = -16;
#line 418 "dggsvd3.f"
    } else if (*ldv < 1 || wantv && *ldv < *p) {
#line 419 "dggsvd3.f"
	*info = -18;
#line 420 "dggsvd3.f"
    } else if (*ldq < 1 || wantq && *ldq < *n) {
#line 421 "dggsvd3.f"
	*info = -20;
#line 422 "dggsvd3.f"
    } else if (*lwork < 1 && ! lquery) {
#line 423 "dggsvd3.f"
	*info = -24;
#line 424 "dggsvd3.f"
    }

/*     Compute workspace */

#line 428 "dggsvd3.f"
    if (*info == 0) {
#line 429 "dggsvd3.f"
	dggsvp3_(jobu, jobv, jobq, m, p, n, &a[a_offset], lda, &b[b_offset], 
		ldb, &tola, &tolb, k, l, &u[u_offset], ldu, &v[v_offset], ldv,
		 &q[q_offset], ldq, &iwork[1], &work[1], &work[1], &c_n1, 
		info, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 432 "dggsvd3.f"
	lwkopt = *n + (integer) work[1];
/* Computing MAX */
#line 433 "dggsvd3.f"
	i__1 = *n << 1;
#line 433 "dggsvd3.f"
	lwkopt = max(i__1,lwkopt);
#line 434 "dggsvd3.f"
	lwkopt = max(1,lwkopt);
#line 435 "dggsvd3.f"
	work[1] = (doublereal) lwkopt;
#line 436 "dggsvd3.f"
    }

#line 438 "dggsvd3.f"
    if (*info != 0) {
#line 439 "dggsvd3.f"
	i__1 = -(*info);
#line 439 "dggsvd3.f"
	xerbla_("DGGSVD3", &i__1, (ftnlen)7);
#line 440 "dggsvd3.f"
	return 0;
#line 441 "dggsvd3.f"
    }
#line 442 "dggsvd3.f"
    if (lquery) {
#line 443 "dggsvd3.f"
	return 0;
#line 444 "dggsvd3.f"
    }

/*     Compute the Frobenius norm of matrices A and B */

#line 448 "dggsvd3.f"
    anorm = dlange_("1", m, n, &a[a_offset], lda, &work[1], (ftnlen)1);
#line 449 "dggsvd3.f"
    bnorm = dlange_("1", p, n, &b[b_offset], ldb, &work[1], (ftnlen)1);

/*     Get machine precision and set up threshold for determining */
/*     the effective numerical rank of the matrices A and B. */

#line 454 "dggsvd3.f"
    ulp = dlamch_("Precision", (ftnlen)9);
#line 455 "dggsvd3.f"
    unfl = dlamch_("Safe Minimum", (ftnlen)12);
#line 456 "dggsvd3.f"
    tola = max(*m,*n) * max(anorm,unfl) * ulp;
#line 457 "dggsvd3.f"
    tolb = max(*p,*n) * max(bnorm,unfl) * ulp;

/*     Preprocessing */

#line 461 "dggsvd3.f"
    i__1 = *lwork - *n;
#line 461 "dggsvd3.f"
    dggsvp3_(jobu, jobv, jobq, m, p, n, &a[a_offset], lda, &b[b_offset], ldb, 
	    &tola, &tolb, k, l, &u[u_offset], ldu, &v[v_offset], ldv, &q[
	    q_offset], ldq, &iwork[1], &work[1], &work[*n + 1], &i__1, info, (
	    ftnlen)1, (ftnlen)1, (ftnlen)1);

/*     Compute the GSVD of two upper "triangular" matrices */

#line 467 "dggsvd3.f"
    dtgsja_(jobu, jobv, jobq, m, p, n, k, l, &a[a_offset], lda, &b[b_offset], 
	    ldb, &tola, &tolb, &alpha[1], &beta[1], &u[u_offset], ldu, &v[
	    v_offset], ldv, &q[q_offset], ldq, &work[1], &ncycle, info, (
	    ftnlen)1, (ftnlen)1, (ftnlen)1);

/*     Sort the singular values and store the pivot indices in IWORK */
/*     Copy ALPHA to WORK, then sort ALPHA in WORK */

#line 474 "dggsvd3.f"
    dcopy_(n, &alpha[1], &c__1, &work[1], &c__1);
/* Computing MIN */
#line 475 "dggsvd3.f"
    i__1 = *l, i__2 = *m - *k;
#line 475 "dggsvd3.f"
    ibnd = min(i__1,i__2);
#line 476 "dggsvd3.f"
    i__1 = ibnd;
#line 476 "dggsvd3.f"
    for (i__ = 1; i__ <= i__1; ++i__) {

/*        Scan for largest ALPHA(K+I) */

#line 480 "dggsvd3.f"
	isub = i__;
#line 481 "dggsvd3.f"
	smax = work[*k + i__];
#line 482 "dggsvd3.f"
	i__2 = ibnd;
#line 482 "dggsvd3.f"
	for (j = i__ + 1; j <= i__2; ++j) {
#line 483 "dggsvd3.f"
	    temp = work[*k + j];
#line 484 "dggsvd3.f"
	    if (temp > smax) {
#line 485 "dggsvd3.f"
		isub = j;
#line 486 "dggsvd3.f"
		smax = temp;
#line 487 "dggsvd3.f"
	    }
#line 488 "dggsvd3.f"
/* L10: */
#line 488 "dggsvd3.f"
	}
#line 489 "dggsvd3.f"
	if (isub != i__) {
#line 490 "dggsvd3.f"
	    work[*k + isub] = work[*k + i__];
#line 491 "dggsvd3.f"
	    work[*k + i__] = smax;
#line 492 "dggsvd3.f"
	    iwork[*k + i__] = *k + isub;
#line 493 "dggsvd3.f"
	} else {
#line 494 "dggsvd3.f"
	    iwork[*k + i__] = *k + i__;
#line 495 "dggsvd3.f"
	}
#line 496 "dggsvd3.f"
/* L20: */
#line 496 "dggsvd3.f"
    }

#line 498 "dggsvd3.f"
    work[1] = (doublereal) lwkopt;
#line 499 "dggsvd3.f"
    return 0;

/*     End of DGGSVD3 */

} /* dggsvd3_ */

