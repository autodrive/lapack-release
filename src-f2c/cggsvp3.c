#line 1 "cggsvp3.f"
/* cggsvp3.f -- translated by f2c (version 20100827).
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

#line 1 "cggsvp3.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};
static integer c_n1 = -1;

/* > \brief \b CGGSVP3 */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CGGSVP3 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cggsvp3
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cggsvp3
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cggsvp3
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CGGSVP3( JOBU, JOBV, JOBQ, M, P, N, A, LDA, B, LDB, */
/*                           TOLA, TOLB, K, L, U, LDU, V, LDV, Q, LDQ, */
/*                           IWORK, RWORK, TAU, WORK, LWORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBQ, JOBU, JOBV */
/*       INTEGER            INFO, K, L, LDA, LDB, LDQ, LDU, LDV, M, N, P, LWORK */
/*       REAL               TOLA, TOLB */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       REAL               RWORK( * ) */
/*       COMPLEX            A( LDA, * ), B( LDB, * ), Q( LDQ, * ), */
/*      $                   TAU( * ), U( LDU, * ), V( LDV, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CGGSVP3 computes unitary matrices U, V and Q such that */
/* > */
/* >                    N-K-L  K    L */
/* >  U**H*A*Q =     K ( 0    A12  A13 )  if M-K-L >= 0; */
/* >                 L ( 0     0   A23 ) */
/* >             M-K-L ( 0     0    0  ) */
/* > */
/* >                  N-K-L  K    L */
/* >         =     K ( 0    A12  A13 )  if M-K-L < 0; */
/* >             M-K ( 0     0   A23 ) */
/* > */
/* >                  N-K-L  K    L */
/* >  V**H*B*Q =   L ( 0     0   B13 ) */
/* >             P-L ( 0     0    0  ) */
/* > */
/* > where the K-by-K matrix A12 and L-by-L matrix B13 are nonsingular */
/* > upper triangular; A23 is L-by-L upper triangular if M-K-L >= 0, */
/* > otherwise A23 is (M-K)-by-L upper trapezoidal.  K+L = the effective */
/* > numerical rank of the (M+P)-by-N matrix (A**H,B**H)**H. */
/* > */
/* > This decomposition is the preprocessing step for computing the */
/* > Generalized Singular Value Decomposition (GSVD), see subroutine */
/* > CGGSVD3. */
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
/* >          A is COMPLEX array, dimension (LDA,N) */
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
/* >          B is COMPLEX array, dimension (LDB,N) */
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
/* >          K + L = effective numerical rank of (A**H,B**H)**H. */
/* > \endverbatim */
/* > */
/* > \param[out] U */
/* > \verbatim */
/* >          U is COMPLEX array, dimension (LDU,M) */
/* >          If JOBU = 'U', U contains the unitary matrix U. */
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
/* >          If JOBV = 'V', V contains the unitary matrix V. */
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
/* >          If JOBQ = 'Q', Q contains the unitary matrix Q. */
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
/* > \param[out] RWORK */
/* > \verbatim */
/* >          RWORK is REAL array, dimension (2*N) */
/* > \endverbatim */
/* > */
/* > \param[out] TAU */
/* > \verbatim */
/* >          TAU is COMPLEX array, dimension (N) */
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

/* > \ingroup complexOTHERcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  The subroutine uses LAPACK subroutine CGEQP3 for the QR factorization */
/* >  with column pivoting to detect the effective numerical rank of the */
/* >  a matrix. It may be replaced by a better rank determination strategy. */
/* > */
/* >  CGGSVP3 replaces the deprecated subroutine CGGSVP. */
/* > */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int cggsvp3_(char *jobu, char *jobv, char *jobq, integer *m, 
	integer *p, integer *n, doublecomplex *a, integer *lda, doublecomplex 
	*b, integer *ldb, doublereal *tola, doublereal *tolb, integer *k, 
	integer *l, doublecomplex *u, integer *ldu, doublecomplex *v, integer 
	*ldv, doublecomplex *q, integer *ldq, integer *iwork, doublereal *
	rwork, doublecomplex *tau, doublecomplex *work, integer *lwork, 
	integer *info, ftnlen jobu_len, ftnlen jobv_len, ftnlen jobq_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, q_dim1, q_offset, u_dim1, 
	    u_offset, v_dim1, v_offset, i__1, i__2, i__3;
    doublecomplex z__1;

    /* Builtin functions */
    double z_abs(doublecomplex *);

    /* Local variables */
    static integer i__, j;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical wantq, wantu, wantv;
    extern /* Subroutine */ int cgeqp3_(integer *, integer *, doublecomplex *,
	     integer *, integer *, doublecomplex *, doublecomplex *, integer *
	    , doublereal *, integer *), cgeqr2_(integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *), cgerq2_(integer *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *), cung2r_(
	    integer *, integer *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, doublecomplex *, integer *), cunm2r_(char *, 
	    char *, integer *, integer *, integer *, doublecomplex *, integer 
	    *, doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, ftnlen, ftnlen), cunmr2_(char *, char *, integer *, 
	    integer *, integer *, doublecomplex *, integer *, doublecomplex *,
	     doublecomplex *, integer *, doublecomplex *, integer *, ftnlen, 
	    ftnlen), clacpy_(char *, integer *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, ftnlen), claset_(char *, 
	    integer *, integer *, doublecomplex *, doublecomplex *, 
	    doublecomplex *, integer *, ftnlen), xerbla_(char *, integer *, 
	    ftnlen), clapmt_(logical *, integer *, integer *, doublecomplex *,
	     integer *, integer *);
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

#line 328 "cggsvp3.f"
    /* Parameter adjustments */
#line 328 "cggsvp3.f"
    a_dim1 = *lda;
#line 328 "cggsvp3.f"
    a_offset = 1 + a_dim1;
#line 328 "cggsvp3.f"
    a -= a_offset;
#line 328 "cggsvp3.f"
    b_dim1 = *ldb;
#line 328 "cggsvp3.f"
    b_offset = 1 + b_dim1;
#line 328 "cggsvp3.f"
    b -= b_offset;
#line 328 "cggsvp3.f"
    u_dim1 = *ldu;
#line 328 "cggsvp3.f"
    u_offset = 1 + u_dim1;
#line 328 "cggsvp3.f"
    u -= u_offset;
#line 328 "cggsvp3.f"
    v_dim1 = *ldv;
#line 328 "cggsvp3.f"
    v_offset = 1 + v_dim1;
#line 328 "cggsvp3.f"
    v -= v_offset;
#line 328 "cggsvp3.f"
    q_dim1 = *ldq;
#line 328 "cggsvp3.f"
    q_offset = 1 + q_dim1;
#line 328 "cggsvp3.f"
    q -= q_offset;
#line 328 "cggsvp3.f"
    --iwork;
#line 328 "cggsvp3.f"
    --rwork;
#line 328 "cggsvp3.f"
    --tau;
#line 328 "cggsvp3.f"
    --work;
#line 328 "cggsvp3.f"

#line 328 "cggsvp3.f"
    /* Function Body */
#line 328 "cggsvp3.f"
    wantu = lsame_(jobu, "U", (ftnlen)1, (ftnlen)1);
#line 329 "cggsvp3.f"
    wantv = lsame_(jobv, "V", (ftnlen)1, (ftnlen)1);
#line 330 "cggsvp3.f"
    wantq = lsame_(jobq, "Q", (ftnlen)1, (ftnlen)1);
#line 331 "cggsvp3.f"
    forwrd = TRUE_;
#line 332 "cggsvp3.f"
    lquery = *lwork == -1;
#line 333 "cggsvp3.f"
    lwkopt = 1;

/*     Test the input arguments */

#line 337 "cggsvp3.f"
    *info = 0;
#line 338 "cggsvp3.f"
    if (! (wantu || lsame_(jobu, "N", (ftnlen)1, (ftnlen)1))) {
#line 339 "cggsvp3.f"
	*info = -1;
#line 340 "cggsvp3.f"
    } else if (! (wantv || lsame_(jobv, "N", (ftnlen)1, (ftnlen)1))) {
#line 341 "cggsvp3.f"
	*info = -2;
#line 342 "cggsvp3.f"
    } else if (! (wantq || lsame_(jobq, "N", (ftnlen)1, (ftnlen)1))) {
#line 343 "cggsvp3.f"
	*info = -3;
#line 344 "cggsvp3.f"
    } else if (*m < 0) {
#line 345 "cggsvp3.f"
	*info = -4;
#line 346 "cggsvp3.f"
    } else if (*p < 0) {
#line 347 "cggsvp3.f"
	*info = -5;
#line 348 "cggsvp3.f"
    } else if (*n < 0) {
#line 349 "cggsvp3.f"
	*info = -6;
#line 350 "cggsvp3.f"
    } else if (*lda < max(1,*m)) {
#line 351 "cggsvp3.f"
	*info = -8;
#line 352 "cggsvp3.f"
    } else if (*ldb < max(1,*p)) {
#line 353 "cggsvp3.f"
	*info = -10;
#line 354 "cggsvp3.f"
    } else if (*ldu < 1 || wantu && *ldu < *m) {
#line 355 "cggsvp3.f"
	*info = -16;
#line 356 "cggsvp3.f"
    } else if (*ldv < 1 || wantv && *ldv < *p) {
#line 357 "cggsvp3.f"
	*info = -18;
#line 358 "cggsvp3.f"
    } else if (*ldq < 1 || wantq && *ldq < *n) {
#line 359 "cggsvp3.f"
	*info = -20;
#line 360 "cggsvp3.f"
    } else if (*lwork < 1 && ! lquery) {
#line 361 "cggsvp3.f"
	*info = -24;
#line 362 "cggsvp3.f"
    }

/*     Compute workspace */

#line 366 "cggsvp3.f"
    if (*info == 0) {
#line 367 "cggsvp3.f"
	cgeqp3_(p, n, &b[b_offset], ldb, &iwork[1], &tau[1], &work[1], &c_n1, 
		&rwork[1], info);
#line 368 "cggsvp3.f"
	lwkopt = (integer) work[1].r;
#line 369 "cggsvp3.f"
	if (wantv) {
#line 370 "cggsvp3.f"
	    lwkopt = max(lwkopt,*p);
#line 371 "cggsvp3.f"
	}
/* Computing MAX */
#line 372 "cggsvp3.f"
	i__1 = lwkopt, i__2 = min(*n,*p);
#line 372 "cggsvp3.f"
	lwkopt = max(i__1,i__2);
#line 373 "cggsvp3.f"
	lwkopt = max(lwkopt,*m);
#line 374 "cggsvp3.f"
	if (wantq) {
#line 375 "cggsvp3.f"
	    lwkopt = max(lwkopt,*n);
#line 376 "cggsvp3.f"
	}
#line 377 "cggsvp3.f"
	cgeqp3_(m, n, &a[a_offset], lda, &iwork[1], &tau[1], &work[1], &c_n1, 
		&rwork[1], info);
/* Computing MAX */
#line 378 "cggsvp3.f"
	i__1 = lwkopt, i__2 = (integer) work[1].r;
#line 378 "cggsvp3.f"
	lwkopt = max(i__1,i__2);
#line 379 "cggsvp3.f"
	lwkopt = max(1,lwkopt);
#line 380 "cggsvp3.f"
	z__1.r = (doublereal) lwkopt, z__1.i = 0.;
#line 380 "cggsvp3.f"
	work[1].r = z__1.r, work[1].i = z__1.i;
#line 381 "cggsvp3.f"
    }

#line 383 "cggsvp3.f"
    if (*info != 0) {
#line 384 "cggsvp3.f"
	i__1 = -(*info);
#line 384 "cggsvp3.f"
	xerbla_("CGGSVP3", &i__1, (ftnlen)7);
#line 385 "cggsvp3.f"
	return 0;
#line 386 "cggsvp3.f"
    }
#line 387 "cggsvp3.f"
    if (lquery) {
#line 388 "cggsvp3.f"
	return 0;
#line 389 "cggsvp3.f"
    }

/*     QR with column pivoting of B: B*P = V*( S11 S12 ) */
/*                                           (  0   0  ) */

#line 394 "cggsvp3.f"
    i__1 = *n;
#line 394 "cggsvp3.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 395 "cggsvp3.f"
	iwork[i__] = 0;
#line 396 "cggsvp3.f"
/* L10: */
#line 396 "cggsvp3.f"
    }
#line 397 "cggsvp3.f"
    cgeqp3_(p, n, &b[b_offset], ldb, &iwork[1], &tau[1], &work[1], lwork, &
	    rwork[1], info);

/*     Update A := A*P */

#line 401 "cggsvp3.f"
    clapmt_(&forwrd, m, n, &a[a_offset], lda, &iwork[1]);

/*     Determine the effective rank of matrix B. */

#line 405 "cggsvp3.f"
    *l = 0;
#line 406 "cggsvp3.f"
    i__1 = min(*p,*n);
#line 406 "cggsvp3.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 407 "cggsvp3.f"
	if (z_abs(&b[i__ + i__ * b_dim1]) > *tolb) {
#line 407 "cggsvp3.f"
	    ++(*l);
#line 407 "cggsvp3.f"
	}
#line 409 "cggsvp3.f"
/* L20: */
#line 409 "cggsvp3.f"
    }

#line 411 "cggsvp3.f"
    if (wantv) {

/*        Copy the details of V, and form V. */

#line 415 "cggsvp3.f"
	claset_("Full", p, p, &c_b1, &c_b1, &v[v_offset], ldv, (ftnlen)4);
#line 416 "cggsvp3.f"
	if (*p > 1) {
#line 416 "cggsvp3.f"
	    i__1 = *p - 1;
#line 416 "cggsvp3.f"
	    clacpy_("Lower", &i__1, n, &b[b_dim1 + 2], ldb, &v[v_dim1 + 2], 
		    ldv, (ftnlen)5);
#line 416 "cggsvp3.f"
	}
#line 419 "cggsvp3.f"
	i__1 = min(*p,*n);
#line 419 "cggsvp3.f"
	cung2r_(p, p, &i__1, &v[v_offset], ldv, &tau[1], &work[1], info);
#line 420 "cggsvp3.f"
    }

/*     Clean up B */

#line 424 "cggsvp3.f"
    i__1 = *l - 1;
#line 424 "cggsvp3.f"
    for (j = 1; j <= i__1; ++j) {
#line 425 "cggsvp3.f"
	i__2 = *l;
#line 425 "cggsvp3.f"
	for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 426 "cggsvp3.f"
	    i__3 = i__ + j * b_dim1;
#line 426 "cggsvp3.f"
	    b[i__3].r = 0., b[i__3].i = 0.;
#line 427 "cggsvp3.f"
/* L30: */
#line 427 "cggsvp3.f"
	}
#line 428 "cggsvp3.f"
/* L40: */
#line 428 "cggsvp3.f"
    }
#line 429 "cggsvp3.f"
    if (*p > *l) {
#line 429 "cggsvp3.f"
	i__1 = *p - *l;
#line 429 "cggsvp3.f"
	claset_("Full", &i__1, n, &c_b1, &c_b1, &b[*l + 1 + b_dim1], ldb, (
		ftnlen)4);
#line 429 "cggsvp3.f"
    }

#line 432 "cggsvp3.f"
    if (wantq) {

/*        Set Q = I and Update Q := Q*P */

#line 436 "cggsvp3.f"
	claset_("Full", n, n, &c_b1, &c_b2, &q[q_offset], ldq, (ftnlen)4);
#line 437 "cggsvp3.f"
	clapmt_(&forwrd, n, n, &q[q_offset], ldq, &iwork[1]);
#line 438 "cggsvp3.f"
    }

#line 440 "cggsvp3.f"
    if (*p >= *l && *n != *l) {

/*        RQ factorization of ( S11 S12 ) = ( 0 S12 )*Z */

#line 444 "cggsvp3.f"
	cgerq2_(l, n, &b[b_offset], ldb, &tau[1], &work[1], info);

/*        Update A := A*Z**H */

#line 448 "cggsvp3.f"
	cunmr2_("Right", "Conjugate transpose", m, n, l, &b[b_offset], ldb, &
		tau[1], &a[a_offset], lda, &work[1], info, (ftnlen)5, (ftnlen)
		19);
#line 450 "cggsvp3.f"
	if (wantq) {

/*           Update Q := Q*Z**H */

#line 454 "cggsvp3.f"
	    cunmr2_("Right", "Conjugate transpose", n, n, l, &b[b_offset], 
		    ldb, &tau[1], &q[q_offset], ldq, &work[1], info, (ftnlen)
		    5, (ftnlen)19);
#line 456 "cggsvp3.f"
	}

/*        Clean up B */

#line 460 "cggsvp3.f"
	i__1 = *n - *l;
#line 460 "cggsvp3.f"
	claset_("Full", l, &i__1, &c_b1, &c_b1, &b[b_offset], ldb, (ftnlen)4);
#line 461 "cggsvp3.f"
	i__1 = *n;
#line 461 "cggsvp3.f"
	for (j = *n - *l + 1; j <= i__1; ++j) {
#line 462 "cggsvp3.f"
	    i__2 = *l;
#line 462 "cggsvp3.f"
	    for (i__ = j - *n + *l + 1; i__ <= i__2; ++i__) {
#line 463 "cggsvp3.f"
		i__3 = i__ + j * b_dim1;
#line 463 "cggsvp3.f"
		b[i__3].r = 0., b[i__3].i = 0.;
#line 464 "cggsvp3.f"
/* L50: */
#line 464 "cggsvp3.f"
	    }
#line 465 "cggsvp3.f"
/* L60: */
#line 465 "cggsvp3.f"
	}

#line 467 "cggsvp3.f"
    }

/*     Let              N-L     L */
/*                A = ( A11    A12 ) M, */

/*     then the following does the complete QR decomposition of A11: */

/*              A11 = U*(  0  T12 )*P1**H */
/*                      (  0   0  ) */

#line 477 "cggsvp3.f"
    i__1 = *n - *l;
#line 477 "cggsvp3.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 478 "cggsvp3.f"
	iwork[i__] = 0;
#line 479 "cggsvp3.f"
/* L70: */
#line 479 "cggsvp3.f"
    }
#line 480 "cggsvp3.f"
    i__1 = *n - *l;
#line 480 "cggsvp3.f"
    cgeqp3_(m, &i__1, &a[a_offset], lda, &iwork[1], &tau[1], &work[1], lwork, 
	    &rwork[1], info);

/*     Determine the effective rank of A11 */

#line 485 "cggsvp3.f"
    *k = 0;
/* Computing MIN */
#line 486 "cggsvp3.f"
    i__2 = *m, i__3 = *n - *l;
#line 486 "cggsvp3.f"
    i__1 = min(i__2,i__3);
#line 486 "cggsvp3.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 487 "cggsvp3.f"
	if (z_abs(&a[i__ + i__ * a_dim1]) > *tola) {
#line 487 "cggsvp3.f"
	    ++(*k);
#line 487 "cggsvp3.f"
	}
#line 489 "cggsvp3.f"
/* L80: */
#line 489 "cggsvp3.f"
    }

/*     Update A12 := U**H*A12, where A12 = A( 1:M, N-L+1:N ) */

/* Computing MIN */
#line 493 "cggsvp3.f"
    i__2 = *m, i__3 = *n - *l;
#line 493 "cggsvp3.f"
    i__1 = min(i__2,i__3);
#line 493 "cggsvp3.f"
    cunm2r_("Left", "Conjugate transpose", m, l, &i__1, &a[a_offset], lda, &
	    tau[1], &a[(*n - *l + 1) * a_dim1 + 1], lda, &work[1], info, (
	    ftnlen)4, (ftnlen)19);

#line 496 "cggsvp3.f"
    if (wantu) {

/*        Copy the details of U, and form U */

#line 500 "cggsvp3.f"
	claset_("Full", m, m, &c_b1, &c_b1, &u[u_offset], ldu, (ftnlen)4);
#line 501 "cggsvp3.f"
	if (*m > 1) {
#line 501 "cggsvp3.f"
	    i__1 = *m - 1;
#line 501 "cggsvp3.f"
	    i__2 = *n - *l;
#line 501 "cggsvp3.f"
	    clacpy_("Lower", &i__1, &i__2, &a[a_dim1 + 2], lda, &u[u_dim1 + 2]
		    , ldu, (ftnlen)5);
#line 501 "cggsvp3.f"
	}
/* Computing MIN */
#line 504 "cggsvp3.f"
	i__2 = *m, i__3 = *n - *l;
#line 504 "cggsvp3.f"
	i__1 = min(i__2,i__3);
#line 504 "cggsvp3.f"
	cung2r_(m, m, &i__1, &u[u_offset], ldu, &tau[1], &work[1], info);
#line 505 "cggsvp3.f"
    }

#line 507 "cggsvp3.f"
    if (wantq) {

/*        Update Q( 1:N, 1:N-L )  = Q( 1:N, 1:N-L )*P1 */

#line 511 "cggsvp3.f"
	i__1 = *n - *l;
#line 511 "cggsvp3.f"
	clapmt_(&forwrd, n, &i__1, &q[q_offset], ldq, &iwork[1]);
#line 512 "cggsvp3.f"
    }

/*     Clean up A: set the strictly lower triangular part of */
/*     A(1:K, 1:K) = 0, and A( K+1:M, 1:N-L ) = 0. */

#line 517 "cggsvp3.f"
    i__1 = *k - 1;
#line 517 "cggsvp3.f"
    for (j = 1; j <= i__1; ++j) {
#line 518 "cggsvp3.f"
	i__2 = *k;
#line 518 "cggsvp3.f"
	for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 519 "cggsvp3.f"
	    i__3 = i__ + j * a_dim1;
#line 519 "cggsvp3.f"
	    a[i__3].r = 0., a[i__3].i = 0.;
#line 520 "cggsvp3.f"
/* L90: */
#line 520 "cggsvp3.f"
	}
#line 521 "cggsvp3.f"
/* L100: */
#line 521 "cggsvp3.f"
    }
#line 522 "cggsvp3.f"
    if (*m > *k) {
#line 522 "cggsvp3.f"
	i__1 = *m - *k;
#line 522 "cggsvp3.f"
	i__2 = *n - *l;
#line 522 "cggsvp3.f"
	claset_("Full", &i__1, &i__2, &c_b1, &c_b1, &a[*k + 1 + a_dim1], lda, 
		(ftnlen)4);
#line 522 "cggsvp3.f"
    }

#line 525 "cggsvp3.f"
    if (*n - *l > *k) {

/*        RQ factorization of ( T11 T12 ) = ( 0 T12 )*Z1 */

#line 529 "cggsvp3.f"
	i__1 = *n - *l;
#line 529 "cggsvp3.f"
	cgerq2_(k, &i__1, &a[a_offset], lda, &tau[1], &work[1], info);

#line 531 "cggsvp3.f"
	if (wantq) {

/*           Update Q( 1:N,1:N-L ) = Q( 1:N,1:N-L )*Z1**H */

#line 535 "cggsvp3.f"
	    i__1 = *n - *l;
#line 535 "cggsvp3.f"
	    cunmr2_("Right", "Conjugate transpose", n, &i__1, k, &a[a_offset],
		     lda, &tau[1], &q[q_offset], ldq, &work[1], info, (ftnlen)
		    5, (ftnlen)19);
#line 537 "cggsvp3.f"
	}

/*        Clean up A */

#line 541 "cggsvp3.f"
	i__1 = *n - *l - *k;
#line 541 "cggsvp3.f"
	claset_("Full", k, &i__1, &c_b1, &c_b1, &a[a_offset], lda, (ftnlen)4);
#line 542 "cggsvp3.f"
	i__1 = *n - *l;
#line 542 "cggsvp3.f"
	for (j = *n - *l - *k + 1; j <= i__1; ++j) {
#line 543 "cggsvp3.f"
	    i__2 = *k;
#line 543 "cggsvp3.f"
	    for (i__ = j - *n + *l + *k + 1; i__ <= i__2; ++i__) {
#line 544 "cggsvp3.f"
		i__3 = i__ + j * a_dim1;
#line 544 "cggsvp3.f"
		a[i__3].r = 0., a[i__3].i = 0.;
#line 545 "cggsvp3.f"
/* L110: */
#line 545 "cggsvp3.f"
	    }
#line 546 "cggsvp3.f"
/* L120: */
#line 546 "cggsvp3.f"
	}

#line 548 "cggsvp3.f"
    }

#line 550 "cggsvp3.f"
    if (*m > *k) {

/*        QR factorization of A( K+1:M,N-L+1:N ) */

#line 554 "cggsvp3.f"
	i__1 = *m - *k;
#line 554 "cggsvp3.f"
	cgeqr2_(&i__1, l, &a[*k + 1 + (*n - *l + 1) * a_dim1], lda, &tau[1], &
		work[1], info);

#line 556 "cggsvp3.f"
	if (wantu) {

/*           Update U(:,K+1:M) := U(:,K+1:M)*U1 */

#line 560 "cggsvp3.f"
	    i__1 = *m - *k;
/* Computing MIN */
#line 560 "cggsvp3.f"
	    i__3 = *m - *k;
#line 560 "cggsvp3.f"
	    i__2 = min(i__3,*l);
#line 560 "cggsvp3.f"
	    cunm2r_("Right", "No transpose", m, &i__1, &i__2, &a[*k + 1 + (*n 
		    - *l + 1) * a_dim1], lda, &tau[1], &u[(*k + 1) * u_dim1 + 
		    1], ldu, &work[1], info, (ftnlen)5, (ftnlen)12);
#line 563 "cggsvp3.f"
	}

/*        Clean up */

#line 567 "cggsvp3.f"
	i__1 = *n;
#line 567 "cggsvp3.f"
	for (j = *n - *l + 1; j <= i__1; ++j) {
#line 568 "cggsvp3.f"
	    i__2 = *m;
#line 568 "cggsvp3.f"
	    for (i__ = j - *n + *k + *l + 1; i__ <= i__2; ++i__) {
#line 569 "cggsvp3.f"
		i__3 = i__ + j * a_dim1;
#line 569 "cggsvp3.f"
		a[i__3].r = 0., a[i__3].i = 0.;
#line 570 "cggsvp3.f"
/* L130: */
#line 570 "cggsvp3.f"
	    }
#line 571 "cggsvp3.f"
/* L140: */
#line 571 "cggsvp3.f"
	}

#line 573 "cggsvp3.f"
    }

#line 575 "cggsvp3.f"
    z__1.r = (doublereal) lwkopt, z__1.i = 0.;
#line 575 "cggsvp3.f"
    work[1].r = z__1.r, work[1].i = z__1.i;
#line 576 "cggsvp3.f"
    return 0;

/*     End of CGGSVP3 */

} /* cggsvp3_ */

