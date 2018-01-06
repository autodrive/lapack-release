#line 1 "dggsvd.f"
/* dggsvd.f -- translated by f2c (version 20100827).
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

#line 1 "dggsvd.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief <b> DGGSVD computes the singular value decomposition (SVD) for OTHER matrices</b> */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DGGSVD + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dggsvd.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dggsvd.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dggsvd.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DGGSVD( JOBU, JOBV, JOBQ, M, N, P, K, L, A, LDA, B, */
/*                          LDB, ALPHA, BETA, U, LDU, V, LDV, Q, LDQ, WORK, */
/*                          IWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBQ, JOBU, JOBV */
/*       INTEGER            INFO, K, L, LDA, LDB, LDQ, LDU, LDV, M, N, P */
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
/* > DGGSVD computes the generalized singular value decomposition (GSVD) */
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
/* >          WORK is DOUBLE PRECISION array, */
/* >                      dimension (max(3*N,M,P)+N) */
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
/* >          = 0:  successful exit */
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
/* >          rank of (A',B')**T. Generally, they are set to */
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

/* > \ingroup doubleOTHERsing */

/* > \par Contributors: */
/*  ================== */
/* > */
/* >     Ming Gu and Huan Ren, Computer Science Division, University of */
/* >     California at Berkeley, USA */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dggsvd_(char *jobu, char *jobv, char *jobq, integer *m, 
	integer *n, integer *p, integer *k, integer *l, doublereal *a, 
	integer *lda, doublereal *b, integer *ldb, doublereal *alpha, 
	doublereal *beta, doublereal *u, integer *ldu, doublereal *v, integer 
	*ldv, doublereal *q, integer *ldq, doublereal *work, integer *iwork, 
	integer *info, ftnlen jobu_len, ftnlen jobv_len, ftnlen jobq_len)
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
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), dggsvp_(
	    char *, char *, char *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen, ftnlen, ftnlen);


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

/*     Test the input parameters */

#line 373 "dggsvd.f"
    /* Parameter adjustments */
#line 373 "dggsvd.f"
    a_dim1 = *lda;
#line 373 "dggsvd.f"
    a_offset = 1 + a_dim1;
#line 373 "dggsvd.f"
    a -= a_offset;
#line 373 "dggsvd.f"
    b_dim1 = *ldb;
#line 373 "dggsvd.f"
    b_offset = 1 + b_dim1;
#line 373 "dggsvd.f"
    b -= b_offset;
#line 373 "dggsvd.f"
    --alpha;
#line 373 "dggsvd.f"
    --beta;
#line 373 "dggsvd.f"
    u_dim1 = *ldu;
#line 373 "dggsvd.f"
    u_offset = 1 + u_dim1;
#line 373 "dggsvd.f"
    u -= u_offset;
#line 373 "dggsvd.f"
    v_dim1 = *ldv;
#line 373 "dggsvd.f"
    v_offset = 1 + v_dim1;
#line 373 "dggsvd.f"
    v -= v_offset;
#line 373 "dggsvd.f"
    q_dim1 = *ldq;
#line 373 "dggsvd.f"
    q_offset = 1 + q_dim1;
#line 373 "dggsvd.f"
    q -= q_offset;
#line 373 "dggsvd.f"
    --work;
#line 373 "dggsvd.f"
    --iwork;
#line 373 "dggsvd.f"

#line 373 "dggsvd.f"
    /* Function Body */
#line 373 "dggsvd.f"
    wantu = lsame_(jobu, "U", (ftnlen)1, (ftnlen)1);
#line 374 "dggsvd.f"
    wantv = lsame_(jobv, "V", (ftnlen)1, (ftnlen)1);
#line 375 "dggsvd.f"
    wantq = lsame_(jobq, "Q", (ftnlen)1, (ftnlen)1);

#line 377 "dggsvd.f"
    *info = 0;
#line 378 "dggsvd.f"
    if (! (wantu || lsame_(jobu, "N", (ftnlen)1, (ftnlen)1))) {
#line 379 "dggsvd.f"
	*info = -1;
#line 380 "dggsvd.f"
    } else if (! (wantv || lsame_(jobv, "N", (ftnlen)1, (ftnlen)1))) {
#line 381 "dggsvd.f"
	*info = -2;
#line 382 "dggsvd.f"
    } else if (! (wantq || lsame_(jobq, "N", (ftnlen)1, (ftnlen)1))) {
#line 383 "dggsvd.f"
	*info = -3;
#line 384 "dggsvd.f"
    } else if (*m < 0) {
#line 385 "dggsvd.f"
	*info = -4;
#line 386 "dggsvd.f"
    } else if (*n < 0) {
#line 387 "dggsvd.f"
	*info = -5;
#line 388 "dggsvd.f"
    } else if (*p < 0) {
#line 389 "dggsvd.f"
	*info = -6;
#line 390 "dggsvd.f"
    } else if (*lda < max(1,*m)) {
#line 391 "dggsvd.f"
	*info = -10;
#line 392 "dggsvd.f"
    } else if (*ldb < max(1,*p)) {
#line 393 "dggsvd.f"
	*info = -12;
#line 394 "dggsvd.f"
    } else if (*ldu < 1 || wantu && *ldu < *m) {
#line 395 "dggsvd.f"
	*info = -16;
#line 396 "dggsvd.f"
    } else if (*ldv < 1 || wantv && *ldv < *p) {
#line 397 "dggsvd.f"
	*info = -18;
#line 398 "dggsvd.f"
    } else if (*ldq < 1 || wantq && *ldq < *n) {
#line 399 "dggsvd.f"
	*info = -20;
#line 400 "dggsvd.f"
    }
#line 401 "dggsvd.f"
    if (*info != 0) {
#line 402 "dggsvd.f"
	i__1 = -(*info);
#line 402 "dggsvd.f"
	xerbla_("DGGSVD", &i__1, (ftnlen)6);
#line 403 "dggsvd.f"
	return 0;
#line 404 "dggsvd.f"
    }

/*     Compute the Frobenius norm of matrices A and B */

#line 408 "dggsvd.f"
    anorm = dlange_("1", m, n, &a[a_offset], lda, &work[1], (ftnlen)1);
#line 409 "dggsvd.f"
    bnorm = dlange_("1", p, n, &b[b_offset], ldb, &work[1], (ftnlen)1);

/*     Get machine precision and set up threshold for determining */
/*     the effective numerical rank of the matrices A and B. */

#line 414 "dggsvd.f"
    ulp = dlamch_("Precision", (ftnlen)9);
#line 415 "dggsvd.f"
    unfl = dlamch_("Safe Minimum", (ftnlen)12);
#line 416 "dggsvd.f"
    tola = max(*m,*n) * max(anorm,unfl) * ulp;
#line 417 "dggsvd.f"
    tolb = max(*p,*n) * max(bnorm,unfl) * ulp;

/*     Preprocessing */

#line 421 "dggsvd.f"
    dggsvp_(jobu, jobv, jobq, m, p, n, &a[a_offset], lda, &b[b_offset], ldb, &
	    tola, &tolb, k, l, &u[u_offset], ldu, &v[v_offset], ldv, &q[
	    q_offset], ldq, &iwork[1], &work[1], &work[*n + 1], info, (ftnlen)
	    1, (ftnlen)1, (ftnlen)1);

/*     Compute the GSVD of two upper "triangular" matrices */

#line 427 "dggsvd.f"
    dtgsja_(jobu, jobv, jobq, m, p, n, k, l, &a[a_offset], lda, &b[b_offset], 
	    ldb, &tola, &tolb, &alpha[1], &beta[1], &u[u_offset], ldu, &v[
	    v_offset], ldv, &q[q_offset], ldq, &work[1], &ncycle, info, (
	    ftnlen)1, (ftnlen)1, (ftnlen)1);

/*     Sort the singular values and store the pivot indices in IWORK */
/*     Copy ALPHA to WORK, then sort ALPHA in WORK */

#line 434 "dggsvd.f"
    dcopy_(n, &alpha[1], &c__1, &work[1], &c__1);
/* Computing MIN */
#line 435 "dggsvd.f"
    i__1 = *l, i__2 = *m - *k;
#line 435 "dggsvd.f"
    ibnd = min(i__1,i__2);
#line 436 "dggsvd.f"
    i__1 = ibnd;
#line 436 "dggsvd.f"
    for (i__ = 1; i__ <= i__1; ++i__) {

/*        Scan for largest ALPHA(K+I) */

#line 440 "dggsvd.f"
	isub = i__;
#line 441 "dggsvd.f"
	smax = work[*k + i__];
#line 442 "dggsvd.f"
	i__2 = ibnd;
#line 442 "dggsvd.f"
	for (j = i__ + 1; j <= i__2; ++j) {
#line 443 "dggsvd.f"
	    temp = work[*k + j];
#line 444 "dggsvd.f"
	    if (temp > smax) {
#line 445 "dggsvd.f"
		isub = j;
#line 446 "dggsvd.f"
		smax = temp;
#line 447 "dggsvd.f"
	    }
#line 448 "dggsvd.f"
/* L10: */
#line 448 "dggsvd.f"
	}
#line 449 "dggsvd.f"
	if (isub != i__) {
#line 450 "dggsvd.f"
	    work[*k + isub] = work[*k + i__];
#line 451 "dggsvd.f"
	    work[*k + i__] = smax;
#line 452 "dggsvd.f"
	    iwork[*k + i__] = *k + isub;
#line 453 "dggsvd.f"
	} else {
#line 454 "dggsvd.f"
	    iwork[*k + i__] = *k + i__;
#line 455 "dggsvd.f"
	}
#line 456 "dggsvd.f"
/* L20: */
#line 456 "dggsvd.f"
    }

#line 458 "dggsvd.f"
    return 0;

/*     End of DGGSVD */

} /* dggsvd_ */

