#line 1 "sggsvp3.f"
/* sggsvp3.f -- translated by f2c (version 20100827).
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

#line 1 "sggsvp3.f"
/* Table of constant values */

static integer c_n1 = -1;
static doublereal c_b14 = 0.;
static doublereal c_b24 = 1.;

/* > \brief \b SGGSVP3 */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SGGSVP3 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sggsvp3
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sggsvp3
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sggsvp3
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SGGSVP3( JOBU, JOBV, JOBQ, M, P, N, A, LDA, B, LDB, */
/*                           TOLA, TOLB, K, L, U, LDU, V, LDV, Q, LDQ, */
/*                           IWORK, TAU, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBQ, JOBU, JOBV */
/*       INTEGER            INFO, K, L, LDA, LDB, LDQ, LDU, LDV, M, N, P, LWORK */
/*       REAL               TOLA, TOLB */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       REAL               A( LDA, * ), B( LDB, * ), Q( LDQ, * ), */
/*      $                   TAU( * ), U( LDU, * ), V( LDV, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SGGSVP3 computes orthogonal matrices U, V and Q such that */
/* > */
/* >                    N-K-L  K    L */
/* >  U**T*A*Q =     K ( 0    A12  A13 )  if M-K-L >= 0; */
/* >                 L ( 0     0   A23 ) */
/* >             M-K-L ( 0     0    0  ) */
/* > */
/* >                  N-K-L  K    L */
/* >         =     K ( 0    A12  A13 )  if M-K-L < 0; */
/* >             M-K ( 0     0   A23 ) */
/* > */
/* >                  N-K-L  K    L */
/* >  V**T*B*Q =   L ( 0     0   B13 ) */
/* >             P-L ( 0     0    0  ) */
/* > */
/* > where the K-by-K matrix A12 and L-by-L matrix B13 are nonsingular */
/* > upper triangular; A23 is L-by-L upper triangular if M-K-L >= 0, */
/* > otherwise A23 is (M-K)-by-L upper trapezoidal.  K+L = the effective */
/* > numerical rank of the (M+P)-by-N matrix (A**T,B**T)**T. */
/* > */
/* > This decomposition is the preprocessing step for computing the */
/* > Generalized Singular Value Decomposition (GSVD), see subroutine */
/* > SGGSVD3. */
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
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is REAL array, dimension (LDA,N) */
/* >          On entry, the M-by-N matrix A. */
/* >          On exit, A contains the triangular (or trapezoidal) matrix */
/* >          described in the Purpose section. */
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
/* >          B is REAL array, dimension (LDB,N) */
/* >          On entry, the P-by-N matrix B. */
/* >          On exit, B contains the triangular matrix described in */
/* >          the Purpose section. */
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
/* >          TOLA and TOLB are the thresholds to determine the effective */
/* >          numerical rank of matrix B and a subblock of A. Generally, */
/* >          they are set to */
/* >             TOLA = MAX(M,N)*norm(A)*MACHEPS, */
/* >             TOLB = MAX(P,N)*norm(B)*MACHEPS. */
/* >          The size of TOLA and TOLB may affect the size of backward */
/* >          errors of the decomposition. */
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
/* >          described in Purpose section. */
/* >          K + L = effective numerical rank of (A**T,B**T)**T. */
/* > \endverbatim */
/* > */
/* > \param[out] U */
/* > \verbatim */
/* >          U is REAL array, dimension (LDU,M) */
/* >          If JOBU = 'U', U contains the orthogonal matrix U. */
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
/* >          V is REAL array, dimension (LDV,P) */
/* >          If JOBV = 'V', V contains the orthogonal matrix V. */
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
/* >          Q is REAL array, dimension (LDQ,N) */
/* >          If JOBQ = 'Q', Q contains the orthogonal matrix Q. */
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
/* > \param[out] IWORK */
/* > \verbatim */
/* >          IWORK is INTEGER array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] TAU */
/* > \verbatim */
/* >          TAU is REAL array, dimension (N) */
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
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date August 2015 */

/* > \ingroup realOTHERcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  The subroutine uses LAPACK subroutine SGEQP3 for the QR factorization */
/* >  with column pivoting to detect the effective numerical rank of the */
/* >  a matrix. It may be replaced by a better rank determination strategy. */
/* > */
/* >  SGGSVP3 replaces the deprecated subroutine SGGSVP. */
/* > */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int sggsvp3_(char *jobu, char *jobv, char *jobq, integer *m, 
	integer *p, integer *n, doublereal *a, integer *lda, doublereal *b, 
	integer *ldb, doublereal *tola, doublereal *tolb, integer *k, integer 
	*l, doublereal *u, integer *ldu, doublereal *v, integer *ldv, 
	doublereal *q, integer *ldq, integer *iwork, doublereal *tau, 
	doublereal *work, integer *lwork, integer *info, ftnlen jobu_len, 
	ftnlen jobv_len, ftnlen jobq_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, q_dim1, q_offset, u_dim1, 
	    u_offset, v_dim1, v_offset, i__1, i__2, i__3;
    doublereal d__1;

    /* Local variables */
    static integer i__, j;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical wantq, wantu, wantv;
    extern /* Subroutine */ int sgeqp3_(integer *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    integer *), sgeqr2_(integer *, integer *, doublereal *, integer *,
	     doublereal *, doublereal *, integer *), sgerq2_(integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *), sorg2r_(integer *, integer *, integer *, doublereal *,
	     integer *, doublereal *, doublereal *, integer *), sorm2r_(char *
	    , char *, integer *, integer *, integer *, doublereal *, integer *
	    , doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen, ftnlen), sormr2_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, ftnlen, ftnlen), xerbla_(char 
	    *, integer *, ftnlen), slacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    slaset_(char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen), slapmt_(logical *, integer *, 
	    integer *, doublereal *, integer *, integer *);
    static logical forwrd;
    static integer lwkopt;
    static logical lquery;


/*  -- LAPACK computational routine (version 3.6.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     August 2015 */


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

/*     Test the input parameters */

#line 319 "sggsvp3.f"
    /* Parameter adjustments */
#line 319 "sggsvp3.f"
    a_dim1 = *lda;
#line 319 "sggsvp3.f"
    a_offset = 1 + a_dim1;
#line 319 "sggsvp3.f"
    a -= a_offset;
#line 319 "sggsvp3.f"
    b_dim1 = *ldb;
#line 319 "sggsvp3.f"
    b_offset = 1 + b_dim1;
#line 319 "sggsvp3.f"
    b -= b_offset;
#line 319 "sggsvp3.f"
    u_dim1 = *ldu;
#line 319 "sggsvp3.f"
    u_offset = 1 + u_dim1;
#line 319 "sggsvp3.f"
    u -= u_offset;
#line 319 "sggsvp3.f"
    v_dim1 = *ldv;
#line 319 "sggsvp3.f"
    v_offset = 1 + v_dim1;
#line 319 "sggsvp3.f"
    v -= v_offset;
#line 319 "sggsvp3.f"
    q_dim1 = *ldq;
#line 319 "sggsvp3.f"
    q_offset = 1 + q_dim1;
#line 319 "sggsvp3.f"
    q -= q_offset;
#line 319 "sggsvp3.f"
    --iwork;
#line 319 "sggsvp3.f"
    --tau;
#line 319 "sggsvp3.f"
    --work;
#line 319 "sggsvp3.f"

#line 319 "sggsvp3.f"
    /* Function Body */
#line 319 "sggsvp3.f"
    wantu = lsame_(jobu, "U", (ftnlen)1, (ftnlen)1);
#line 320 "sggsvp3.f"
    wantv = lsame_(jobv, "V", (ftnlen)1, (ftnlen)1);
#line 321 "sggsvp3.f"
    wantq = lsame_(jobq, "Q", (ftnlen)1, (ftnlen)1);
#line 322 "sggsvp3.f"
    forwrd = TRUE_;
#line 323 "sggsvp3.f"
    lquery = *lwork == -1;
#line 324 "sggsvp3.f"
    lwkopt = 1;

/*     Test the input arguments */

#line 328 "sggsvp3.f"
    *info = 0;
#line 329 "sggsvp3.f"
    if (! (wantu || lsame_(jobu, "N", (ftnlen)1, (ftnlen)1))) {
#line 330 "sggsvp3.f"
	*info = -1;
#line 331 "sggsvp3.f"
    } else if (! (wantv || lsame_(jobv, "N", (ftnlen)1, (ftnlen)1))) {
#line 332 "sggsvp3.f"
	*info = -2;
#line 333 "sggsvp3.f"
    } else if (! (wantq || lsame_(jobq, "N", (ftnlen)1, (ftnlen)1))) {
#line 334 "sggsvp3.f"
	*info = -3;
#line 335 "sggsvp3.f"
    } else if (*m < 0) {
#line 336 "sggsvp3.f"
	*info = -4;
#line 337 "sggsvp3.f"
    } else if (*p < 0) {
#line 338 "sggsvp3.f"
	*info = -5;
#line 339 "sggsvp3.f"
    } else if (*n < 0) {
#line 340 "sggsvp3.f"
	*info = -6;
#line 341 "sggsvp3.f"
    } else if (*lda < max(1,*m)) {
#line 342 "sggsvp3.f"
	*info = -8;
#line 343 "sggsvp3.f"
    } else if (*ldb < max(1,*p)) {
#line 344 "sggsvp3.f"
	*info = -10;
#line 345 "sggsvp3.f"
    } else if (*ldu < 1 || wantu && *ldu < *m) {
#line 346 "sggsvp3.f"
	*info = -16;
#line 347 "sggsvp3.f"
    } else if (*ldv < 1 || wantv && *ldv < *p) {
#line 348 "sggsvp3.f"
	*info = -18;
#line 349 "sggsvp3.f"
    } else if (*ldq < 1 || wantq && *ldq < *n) {
#line 350 "sggsvp3.f"
	*info = -20;
#line 351 "sggsvp3.f"
    } else if (*lwork < 1 && ! lquery) {
#line 352 "sggsvp3.f"
	*info = -24;
#line 353 "sggsvp3.f"
    }

/*     Compute workspace */

#line 357 "sggsvp3.f"
    if (*info == 0) {
#line 358 "sggsvp3.f"
	sgeqp3_(p, n, &b[b_offset], ldb, &iwork[1], &tau[1], &work[1], &c_n1, 
		info);
#line 359 "sggsvp3.f"
	lwkopt = (integer) work[1];
#line 360 "sggsvp3.f"
	if (wantv) {
#line 361 "sggsvp3.f"
	    lwkopt = max(lwkopt,*p);
#line 362 "sggsvp3.f"
	}
/* Computing MAX */
#line 363 "sggsvp3.f"
	i__1 = lwkopt, i__2 = min(*n,*p);
#line 363 "sggsvp3.f"
	lwkopt = max(i__1,i__2);
#line 364 "sggsvp3.f"
	lwkopt = max(lwkopt,*m);
#line 365 "sggsvp3.f"
	if (wantq) {
#line 366 "sggsvp3.f"
	    lwkopt = max(lwkopt,*n);
#line 367 "sggsvp3.f"
	}
#line 368 "sggsvp3.f"
	sgeqp3_(m, n, &a[a_offset], lda, &iwork[1], &tau[1], &work[1], &c_n1, 
		info);
/* Computing MAX */
#line 369 "sggsvp3.f"
	i__1 = lwkopt, i__2 = (integer) work[1];
#line 369 "sggsvp3.f"
	lwkopt = max(i__1,i__2);
#line 370 "sggsvp3.f"
	lwkopt = max(1,lwkopt);
#line 371 "sggsvp3.f"
	work[1] = (doublereal) lwkopt;
#line 372 "sggsvp3.f"
    }

#line 374 "sggsvp3.f"
    if (*info != 0) {
#line 375 "sggsvp3.f"
	i__1 = -(*info);
#line 375 "sggsvp3.f"
	xerbla_("SGGSVP3", &i__1, (ftnlen)7);
#line 376 "sggsvp3.f"
	return 0;
#line 377 "sggsvp3.f"
    }
#line 378 "sggsvp3.f"
    if (lquery) {
#line 379 "sggsvp3.f"
	return 0;
#line 380 "sggsvp3.f"
    }

/*     QR with column pivoting of B: B*P = V*( S11 S12 ) */
/*                                           (  0   0  ) */

#line 385 "sggsvp3.f"
    i__1 = *n;
#line 385 "sggsvp3.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 386 "sggsvp3.f"
	iwork[i__] = 0;
#line 387 "sggsvp3.f"
/* L10: */
#line 387 "sggsvp3.f"
    }
#line 388 "sggsvp3.f"
    sgeqp3_(p, n, &b[b_offset], ldb, &iwork[1], &tau[1], &work[1], lwork, 
	    info);

/*     Update A := A*P */

#line 392 "sggsvp3.f"
    slapmt_(&forwrd, m, n, &a[a_offset], lda, &iwork[1]);

/*     Determine the effective rank of matrix B. */

#line 396 "sggsvp3.f"
    *l = 0;
#line 397 "sggsvp3.f"
    i__1 = min(*p,*n);
#line 397 "sggsvp3.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 398 "sggsvp3.f"
	if ((d__1 = b[i__ + i__ * b_dim1], abs(d__1)) > *tolb) {
#line 398 "sggsvp3.f"
	    ++(*l);
#line 398 "sggsvp3.f"
	}
#line 400 "sggsvp3.f"
/* L20: */
#line 400 "sggsvp3.f"
    }

#line 402 "sggsvp3.f"
    if (wantv) {

/*        Copy the details of V, and form V. */

#line 406 "sggsvp3.f"
	slaset_("Full", p, p, &c_b14, &c_b14, &v[v_offset], ldv, (ftnlen)4);
#line 407 "sggsvp3.f"
	if (*p > 1) {
#line 407 "sggsvp3.f"
	    i__1 = *p - 1;
#line 407 "sggsvp3.f"
	    slacpy_("Lower", &i__1, n, &b[b_dim1 + 2], ldb, &v[v_dim1 + 2], 
		    ldv, (ftnlen)5);
#line 407 "sggsvp3.f"
	}
#line 410 "sggsvp3.f"
	i__1 = min(*p,*n);
#line 410 "sggsvp3.f"
	sorg2r_(p, p, &i__1, &v[v_offset], ldv, &tau[1], &work[1], info);
#line 411 "sggsvp3.f"
    }

/*     Clean up B */

#line 415 "sggsvp3.f"
    i__1 = *l - 1;
#line 415 "sggsvp3.f"
    for (j = 1; j <= i__1; ++j) {
#line 416 "sggsvp3.f"
	i__2 = *l;
#line 416 "sggsvp3.f"
	for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 417 "sggsvp3.f"
	    b[i__ + j * b_dim1] = 0.;
#line 418 "sggsvp3.f"
/* L30: */
#line 418 "sggsvp3.f"
	}
#line 419 "sggsvp3.f"
/* L40: */
#line 419 "sggsvp3.f"
    }
#line 420 "sggsvp3.f"
    if (*p > *l) {
#line 420 "sggsvp3.f"
	i__1 = *p - *l;
#line 420 "sggsvp3.f"
	slaset_("Full", &i__1, n, &c_b14, &c_b14, &b[*l + 1 + b_dim1], ldb, (
		ftnlen)4);
#line 420 "sggsvp3.f"
    }

#line 423 "sggsvp3.f"
    if (wantq) {

/*        Set Q = I and Update Q := Q*P */

#line 427 "sggsvp3.f"
	slaset_("Full", n, n, &c_b14, &c_b24, &q[q_offset], ldq, (ftnlen)4);
#line 428 "sggsvp3.f"
	slapmt_(&forwrd, n, n, &q[q_offset], ldq, &iwork[1]);
#line 429 "sggsvp3.f"
    }

#line 431 "sggsvp3.f"
    if (*p >= *l && *n != *l) {

/*        RQ factorization of (S11 S12): ( S11 S12 ) = ( 0 S12 )*Z */

#line 435 "sggsvp3.f"
	sgerq2_(l, n, &b[b_offset], ldb, &tau[1], &work[1], info);

/*        Update A := A*Z**T */

#line 439 "sggsvp3.f"
	sormr2_("Right", "Transpose", m, n, l, &b[b_offset], ldb, &tau[1], &a[
		a_offset], lda, &work[1], info, (ftnlen)5, (ftnlen)9);

#line 442 "sggsvp3.f"
	if (wantq) {

/*           Update Q := Q*Z**T */

#line 446 "sggsvp3.f"
	    sormr2_("Right", "Transpose", n, n, l, &b[b_offset], ldb, &tau[1],
		     &q[q_offset], ldq, &work[1], info, (ftnlen)5, (ftnlen)9);
#line 448 "sggsvp3.f"
	}

/*        Clean up B */

#line 452 "sggsvp3.f"
	i__1 = *n - *l;
#line 452 "sggsvp3.f"
	slaset_("Full", l, &i__1, &c_b14, &c_b14, &b[b_offset], ldb, (ftnlen)
		4);
#line 453 "sggsvp3.f"
	i__1 = *n;
#line 453 "sggsvp3.f"
	for (j = *n - *l + 1; j <= i__1; ++j) {
#line 454 "sggsvp3.f"
	    i__2 = *l;
#line 454 "sggsvp3.f"
	    for (i__ = j - *n + *l + 1; i__ <= i__2; ++i__) {
#line 455 "sggsvp3.f"
		b[i__ + j * b_dim1] = 0.;
#line 456 "sggsvp3.f"
/* L50: */
#line 456 "sggsvp3.f"
	    }
#line 457 "sggsvp3.f"
/* L60: */
#line 457 "sggsvp3.f"
	}

#line 459 "sggsvp3.f"
    }

/*     Let              N-L     L */
/*                A = ( A11    A12 ) M, */

/*     then the following does the complete QR decomposition of A11: */

/*              A11 = U*(  0  T12 )*P1**T */
/*                      (  0   0  ) */

#line 469 "sggsvp3.f"
    i__1 = *n - *l;
#line 469 "sggsvp3.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 470 "sggsvp3.f"
	iwork[i__] = 0;
#line 471 "sggsvp3.f"
/* L70: */
#line 471 "sggsvp3.f"
    }
#line 472 "sggsvp3.f"
    i__1 = *n - *l;
#line 472 "sggsvp3.f"
    sgeqp3_(m, &i__1, &a[a_offset], lda, &iwork[1], &tau[1], &work[1], lwork, 
	    info);

/*     Determine the effective rank of A11 */

#line 476 "sggsvp3.f"
    *k = 0;
/* Computing MIN */
#line 477 "sggsvp3.f"
    i__2 = *m, i__3 = *n - *l;
#line 477 "sggsvp3.f"
    i__1 = min(i__2,i__3);
#line 477 "sggsvp3.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 478 "sggsvp3.f"
	if ((d__1 = a[i__ + i__ * a_dim1], abs(d__1)) > *tola) {
#line 478 "sggsvp3.f"
	    ++(*k);
#line 478 "sggsvp3.f"
	}
#line 480 "sggsvp3.f"
/* L80: */
#line 480 "sggsvp3.f"
    }

/*     Update A12 := U**T*A12, where A12 = A( 1:M, N-L+1:N ) */

/* Computing MIN */
#line 484 "sggsvp3.f"
    i__2 = *m, i__3 = *n - *l;
#line 484 "sggsvp3.f"
    i__1 = min(i__2,i__3);
#line 484 "sggsvp3.f"
    sorm2r_("Left", "Transpose", m, l, &i__1, &a[a_offset], lda, &tau[1], &a[(
	    *n - *l + 1) * a_dim1 + 1], lda, &work[1], info, (ftnlen)4, (
	    ftnlen)9);

#line 487 "sggsvp3.f"
    if (wantu) {

/*        Copy the details of U, and form U */

#line 491 "sggsvp3.f"
	slaset_("Full", m, m, &c_b14, &c_b14, &u[u_offset], ldu, (ftnlen)4);
#line 492 "sggsvp3.f"
	if (*m > 1) {
#line 492 "sggsvp3.f"
	    i__1 = *m - 1;
#line 492 "sggsvp3.f"
	    i__2 = *n - *l;
#line 492 "sggsvp3.f"
	    slacpy_("Lower", &i__1, &i__2, &a[a_dim1 + 2], lda, &u[u_dim1 + 2]
		    , ldu, (ftnlen)5);
#line 492 "sggsvp3.f"
	}
/* Computing MIN */
#line 495 "sggsvp3.f"
	i__2 = *m, i__3 = *n - *l;
#line 495 "sggsvp3.f"
	i__1 = min(i__2,i__3);
#line 495 "sggsvp3.f"
	sorg2r_(m, m, &i__1, &u[u_offset], ldu, &tau[1], &work[1], info);
#line 496 "sggsvp3.f"
    }

#line 498 "sggsvp3.f"
    if (wantq) {

/*        Update Q( 1:N, 1:N-L )  = Q( 1:N, 1:N-L )*P1 */

#line 502 "sggsvp3.f"
	i__1 = *n - *l;
#line 502 "sggsvp3.f"
	slapmt_(&forwrd, n, &i__1, &q[q_offset], ldq, &iwork[1]);
#line 503 "sggsvp3.f"
    }

/*     Clean up A: set the strictly lower triangular part of */
/*     A(1:K, 1:K) = 0, and A( K+1:M, 1:N-L ) = 0. */

#line 508 "sggsvp3.f"
    i__1 = *k - 1;
#line 508 "sggsvp3.f"
    for (j = 1; j <= i__1; ++j) {
#line 509 "sggsvp3.f"
	i__2 = *k;
#line 509 "sggsvp3.f"
	for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 510 "sggsvp3.f"
	    a[i__ + j * a_dim1] = 0.;
#line 511 "sggsvp3.f"
/* L90: */
#line 511 "sggsvp3.f"
	}
#line 512 "sggsvp3.f"
/* L100: */
#line 512 "sggsvp3.f"
    }
#line 513 "sggsvp3.f"
    if (*m > *k) {
#line 513 "sggsvp3.f"
	i__1 = *m - *k;
#line 513 "sggsvp3.f"
	i__2 = *n - *l;
#line 513 "sggsvp3.f"
	slaset_("Full", &i__1, &i__2, &c_b14, &c_b14, &a[*k + 1 + a_dim1], 
		lda, (ftnlen)4);
#line 513 "sggsvp3.f"
    }

#line 516 "sggsvp3.f"
    if (*n - *l > *k) {

/*        RQ factorization of ( T11 T12 ) = ( 0 T12 )*Z1 */

#line 520 "sggsvp3.f"
	i__1 = *n - *l;
#line 520 "sggsvp3.f"
	sgerq2_(k, &i__1, &a[a_offset], lda, &tau[1], &work[1], info);

#line 522 "sggsvp3.f"
	if (wantq) {

/*           Update Q( 1:N,1:N-L ) = Q( 1:N,1:N-L )*Z1**T */

#line 526 "sggsvp3.f"
	    i__1 = *n - *l;
#line 526 "sggsvp3.f"
	    sormr2_("Right", "Transpose", n, &i__1, k, &a[a_offset], lda, &
		    tau[1], &q[q_offset], ldq, &work[1], info, (ftnlen)5, (
		    ftnlen)9);
#line 528 "sggsvp3.f"
	}

/*        Clean up A */

#line 532 "sggsvp3.f"
	i__1 = *n - *l - *k;
#line 532 "sggsvp3.f"
	slaset_("Full", k, &i__1, &c_b14, &c_b14, &a[a_offset], lda, (ftnlen)
		4);
#line 533 "sggsvp3.f"
	i__1 = *n - *l;
#line 533 "sggsvp3.f"
	for (j = *n - *l - *k + 1; j <= i__1; ++j) {
#line 534 "sggsvp3.f"
	    i__2 = *k;
#line 534 "sggsvp3.f"
	    for (i__ = j - *n + *l + *k + 1; i__ <= i__2; ++i__) {
#line 535 "sggsvp3.f"
		a[i__ + j * a_dim1] = 0.;
#line 536 "sggsvp3.f"
/* L110: */
#line 536 "sggsvp3.f"
	    }
#line 537 "sggsvp3.f"
/* L120: */
#line 537 "sggsvp3.f"
	}

#line 539 "sggsvp3.f"
    }

#line 541 "sggsvp3.f"
    if (*m > *k) {

/*        QR factorization of A( K+1:M,N-L+1:N ) */

#line 545 "sggsvp3.f"
	i__1 = *m - *k;
#line 545 "sggsvp3.f"
	sgeqr2_(&i__1, l, &a[*k + 1 + (*n - *l + 1) * a_dim1], lda, &tau[1], &
		work[1], info);

#line 547 "sggsvp3.f"
	if (wantu) {

/*           Update U(:,K+1:M) := U(:,K+1:M)*U1 */

#line 551 "sggsvp3.f"
	    i__1 = *m - *k;
/* Computing MIN */
#line 551 "sggsvp3.f"
	    i__3 = *m - *k;
#line 551 "sggsvp3.f"
	    i__2 = min(i__3,*l);
#line 551 "sggsvp3.f"
	    sorm2r_("Right", "No transpose", m, &i__1, &i__2, &a[*k + 1 + (*n 
		    - *l + 1) * a_dim1], lda, &tau[1], &u[(*k + 1) * u_dim1 + 
		    1], ldu, &work[1], info, (ftnlen)5, (ftnlen)12);
#line 554 "sggsvp3.f"
	}

/*        Clean up */

#line 558 "sggsvp3.f"
	i__1 = *n;
#line 558 "sggsvp3.f"
	for (j = *n - *l + 1; j <= i__1; ++j) {
#line 559 "sggsvp3.f"
	    i__2 = *m;
#line 559 "sggsvp3.f"
	    for (i__ = j - *n + *k + *l + 1; i__ <= i__2; ++i__) {
#line 560 "sggsvp3.f"
		a[i__ + j * a_dim1] = 0.;
#line 561 "sggsvp3.f"
/* L130: */
#line 561 "sggsvp3.f"
	    }
#line 562 "sggsvp3.f"
/* L140: */
#line 562 "sggsvp3.f"
	}

#line 564 "sggsvp3.f"
    }

#line 566 "sggsvp3.f"
    work[1] = (doublereal) lwkopt;
#line 567 "sggsvp3.f"
    return 0;

/*     End of SGGSVP3 */

} /* sggsvp3_ */

