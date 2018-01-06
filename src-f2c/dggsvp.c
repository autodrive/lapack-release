#line 1 "dggsvp.f"
/* dggsvp.f -- translated by f2c (version 20100827).
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

#line 1 "dggsvp.f"
/* Table of constant values */

static doublereal c_b12 = 0.;
static doublereal c_b22 = 1.;

/* > \brief \b DGGSVP */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DGGSVP + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dggsvp.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dggsvp.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dggsvp.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DGGSVP( JOBU, JOBV, JOBQ, M, P, N, A, LDA, B, LDB, */
/*                          TOLA, TOLB, K, L, U, LDU, V, LDV, Q, LDQ, */
/*                          IWORK, TAU, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBQ, JOBU, JOBV */
/*       INTEGER            INFO, K, L, LDA, LDB, LDQ, LDU, LDV, M, N, P */
/*       DOUBLE PRECISION   TOLA, TOLB */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), Q( LDQ, * ), */
/*      $                   TAU( * ), U( LDU, * ), V( LDV, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DGGSVP computes orthogonal matrices U, V and Q such that */
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
/* > DGGSVD. */
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
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
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
/* >          B is DOUBLE PRECISION array, dimension (LDB,N) */
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
/* >          TOLA is DOUBLE PRECISION */
/* > \endverbatim */
/* > */
/* > \param[in] TOLB */
/* > \verbatim */
/* >          TOLB is DOUBLE PRECISION */
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
/* >          U is DOUBLE PRECISION array, dimension (LDU,M) */
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
/* >          V is DOUBLE PRECISION array, dimension (LDV,P) */
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
/* >          Q is DOUBLE PRECISION array, dimension (LDQ,N) */
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
/* >          TAU is DOUBLE PRECISION array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension (max(3*N,M,P)) */
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

/* > \date November 2011 */

/* > \ingroup doubleOTHERcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* >  The subroutine uses LAPACK subroutine DGEQPF for the QR factorization */
/* >  with column pivoting to detect the effective numerical rank of the */
/* >  a matrix. It may be replaced by a better rank determination strategy. */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dggsvp_(char *jobu, char *jobv, char *jobq, integer *m, 
	integer *p, integer *n, doublereal *a, integer *lda, doublereal *b, 
	integer *ldb, doublereal *tola, doublereal *tolb, integer *k, integer 
	*l, doublereal *u, integer *ldu, doublereal *v, integer *ldv, 
	doublereal *q, integer *ldq, integer *iwork, doublereal *tau, 
	doublereal *work, integer *info, ftnlen jobu_len, ftnlen jobv_len, 
	ftnlen jobq_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, q_dim1, q_offset, u_dim1, 
	    u_offset, v_dim1, v_offset, i__1, i__2, i__3;
    doublereal d__1;

    /* Local variables */
    static integer i__, j;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical wantq, wantu, wantv;
    extern /* Subroutine */ int dgeqr2_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *), dgerq2_(
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *), dorg2r_(integer *, integer *, integer *,
	     doublereal *, integer *, doublereal *, doublereal *, integer *), 
	    dorm2r_(char *, char *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen), dormr2_(char *, char *, 
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen, ftnlen), dgeqpf_(integer *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, doublereal *, integer *), 
	    dlacpy_(char *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen), dlaset_(char *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    ftnlen), xerbla_(char *, integer *, ftnlen), dlapmt_(logical *, 
	    integer *, integer *, doublereal *, integer *, integer *);
    static logical forwrd;


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

/*     Test the input parameters */

#line 298 "dggsvp.f"
    /* Parameter adjustments */
#line 298 "dggsvp.f"
    a_dim1 = *lda;
#line 298 "dggsvp.f"
    a_offset = 1 + a_dim1;
#line 298 "dggsvp.f"
    a -= a_offset;
#line 298 "dggsvp.f"
    b_dim1 = *ldb;
#line 298 "dggsvp.f"
    b_offset = 1 + b_dim1;
#line 298 "dggsvp.f"
    b -= b_offset;
#line 298 "dggsvp.f"
    u_dim1 = *ldu;
#line 298 "dggsvp.f"
    u_offset = 1 + u_dim1;
#line 298 "dggsvp.f"
    u -= u_offset;
#line 298 "dggsvp.f"
    v_dim1 = *ldv;
#line 298 "dggsvp.f"
    v_offset = 1 + v_dim1;
#line 298 "dggsvp.f"
    v -= v_offset;
#line 298 "dggsvp.f"
    q_dim1 = *ldq;
#line 298 "dggsvp.f"
    q_offset = 1 + q_dim1;
#line 298 "dggsvp.f"
    q -= q_offset;
#line 298 "dggsvp.f"
    --iwork;
#line 298 "dggsvp.f"
    --tau;
#line 298 "dggsvp.f"
    --work;
#line 298 "dggsvp.f"

#line 298 "dggsvp.f"
    /* Function Body */
#line 298 "dggsvp.f"
    wantu = lsame_(jobu, "U", (ftnlen)1, (ftnlen)1);
#line 299 "dggsvp.f"
    wantv = lsame_(jobv, "V", (ftnlen)1, (ftnlen)1);
#line 300 "dggsvp.f"
    wantq = lsame_(jobq, "Q", (ftnlen)1, (ftnlen)1);
#line 301 "dggsvp.f"
    forwrd = TRUE_;

#line 303 "dggsvp.f"
    *info = 0;
#line 304 "dggsvp.f"
    if (! (wantu || lsame_(jobu, "N", (ftnlen)1, (ftnlen)1))) {
#line 305 "dggsvp.f"
	*info = -1;
#line 306 "dggsvp.f"
    } else if (! (wantv || lsame_(jobv, "N", (ftnlen)1, (ftnlen)1))) {
#line 307 "dggsvp.f"
	*info = -2;
#line 308 "dggsvp.f"
    } else if (! (wantq || lsame_(jobq, "N", (ftnlen)1, (ftnlen)1))) {
#line 309 "dggsvp.f"
	*info = -3;
#line 310 "dggsvp.f"
    } else if (*m < 0) {
#line 311 "dggsvp.f"
	*info = -4;
#line 312 "dggsvp.f"
    } else if (*p < 0) {
#line 313 "dggsvp.f"
	*info = -5;
#line 314 "dggsvp.f"
    } else if (*n < 0) {
#line 315 "dggsvp.f"
	*info = -6;
#line 316 "dggsvp.f"
    } else if (*lda < max(1,*m)) {
#line 317 "dggsvp.f"
	*info = -8;
#line 318 "dggsvp.f"
    } else if (*ldb < max(1,*p)) {
#line 319 "dggsvp.f"
	*info = -10;
#line 320 "dggsvp.f"
    } else if (*ldu < 1 || wantu && *ldu < *m) {
#line 321 "dggsvp.f"
	*info = -16;
#line 322 "dggsvp.f"
    } else if (*ldv < 1 || wantv && *ldv < *p) {
#line 323 "dggsvp.f"
	*info = -18;
#line 324 "dggsvp.f"
    } else if (*ldq < 1 || wantq && *ldq < *n) {
#line 325 "dggsvp.f"
	*info = -20;
#line 326 "dggsvp.f"
    }
#line 327 "dggsvp.f"
    if (*info != 0) {
#line 328 "dggsvp.f"
	i__1 = -(*info);
#line 328 "dggsvp.f"
	xerbla_("DGGSVP", &i__1, (ftnlen)6);
#line 329 "dggsvp.f"
	return 0;
#line 330 "dggsvp.f"
    }

/*     QR with column pivoting of B: B*P = V*( S11 S12 ) */
/*                                           (  0   0  ) */

#line 335 "dggsvp.f"
    i__1 = *n;
#line 335 "dggsvp.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 336 "dggsvp.f"
	iwork[i__] = 0;
#line 337 "dggsvp.f"
/* L10: */
#line 337 "dggsvp.f"
    }
#line 338 "dggsvp.f"
    dgeqpf_(p, n, &b[b_offset], ldb, &iwork[1], &tau[1], &work[1], info);

/*     Update A := A*P */

#line 342 "dggsvp.f"
    dlapmt_(&forwrd, m, n, &a[a_offset], lda, &iwork[1]);

/*     Determine the effective rank of matrix B. */

#line 346 "dggsvp.f"
    *l = 0;
#line 347 "dggsvp.f"
    i__1 = min(*p,*n);
#line 347 "dggsvp.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 348 "dggsvp.f"
	if ((d__1 = b[i__ + i__ * b_dim1], abs(d__1)) > *tolb) {
#line 348 "dggsvp.f"
	    ++(*l);
#line 348 "dggsvp.f"
	}
#line 350 "dggsvp.f"
/* L20: */
#line 350 "dggsvp.f"
    }

#line 352 "dggsvp.f"
    if (wantv) {

/*        Copy the details of V, and form V. */

#line 356 "dggsvp.f"
	dlaset_("Full", p, p, &c_b12, &c_b12, &v[v_offset], ldv, (ftnlen)4);
#line 357 "dggsvp.f"
	if (*p > 1) {
#line 357 "dggsvp.f"
	    i__1 = *p - 1;
#line 357 "dggsvp.f"
	    dlacpy_("Lower", &i__1, n, &b[b_dim1 + 2], ldb, &v[v_dim1 + 2], 
		    ldv, (ftnlen)5);
#line 357 "dggsvp.f"
	}
#line 360 "dggsvp.f"
	i__1 = min(*p,*n);
#line 360 "dggsvp.f"
	dorg2r_(p, p, &i__1, &v[v_offset], ldv, &tau[1], &work[1], info);
#line 361 "dggsvp.f"
    }

/*     Clean up B */

#line 365 "dggsvp.f"
    i__1 = *l - 1;
#line 365 "dggsvp.f"
    for (j = 1; j <= i__1; ++j) {
#line 366 "dggsvp.f"
	i__2 = *l;
#line 366 "dggsvp.f"
	for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 367 "dggsvp.f"
	    b[i__ + j * b_dim1] = 0.;
#line 368 "dggsvp.f"
/* L30: */
#line 368 "dggsvp.f"
	}
#line 369 "dggsvp.f"
/* L40: */
#line 369 "dggsvp.f"
    }
#line 370 "dggsvp.f"
    if (*p > *l) {
#line 370 "dggsvp.f"
	i__1 = *p - *l;
#line 370 "dggsvp.f"
	dlaset_("Full", &i__1, n, &c_b12, &c_b12, &b[*l + 1 + b_dim1], ldb, (
		ftnlen)4);
#line 370 "dggsvp.f"
    }

#line 373 "dggsvp.f"
    if (wantq) {

/*        Set Q = I and Update Q := Q*P */

#line 377 "dggsvp.f"
	dlaset_("Full", n, n, &c_b12, &c_b22, &q[q_offset], ldq, (ftnlen)4);
#line 378 "dggsvp.f"
	dlapmt_(&forwrd, n, n, &q[q_offset], ldq, &iwork[1]);
#line 379 "dggsvp.f"
    }

#line 381 "dggsvp.f"
    if (*p >= *l && *n != *l) {

/*        RQ factorization of (S11 S12): ( S11 S12 ) = ( 0 S12 )*Z */

#line 385 "dggsvp.f"
	dgerq2_(l, n, &b[b_offset], ldb, &tau[1], &work[1], info);

/*        Update A := A*Z**T */

#line 389 "dggsvp.f"
	dormr2_("Right", "Transpose", m, n, l, &b[b_offset], ldb, &tau[1], &a[
		a_offset], lda, &work[1], info, (ftnlen)5, (ftnlen)9);

#line 392 "dggsvp.f"
	if (wantq) {

/*           Update Q := Q*Z**T */

#line 396 "dggsvp.f"
	    dormr2_("Right", "Transpose", n, n, l, &b[b_offset], ldb, &tau[1],
		     &q[q_offset], ldq, &work[1], info, (ftnlen)5, (ftnlen)9);
#line 398 "dggsvp.f"
	}

/*        Clean up B */

#line 402 "dggsvp.f"
	i__1 = *n - *l;
#line 402 "dggsvp.f"
	dlaset_("Full", l, &i__1, &c_b12, &c_b12, &b[b_offset], ldb, (ftnlen)
		4);
#line 403 "dggsvp.f"
	i__1 = *n;
#line 403 "dggsvp.f"
	for (j = *n - *l + 1; j <= i__1; ++j) {
#line 404 "dggsvp.f"
	    i__2 = *l;
#line 404 "dggsvp.f"
	    for (i__ = j - *n + *l + 1; i__ <= i__2; ++i__) {
#line 405 "dggsvp.f"
		b[i__ + j * b_dim1] = 0.;
#line 406 "dggsvp.f"
/* L50: */
#line 406 "dggsvp.f"
	    }
#line 407 "dggsvp.f"
/* L60: */
#line 407 "dggsvp.f"
	}

#line 409 "dggsvp.f"
    }

/*     Let              N-L     L */
/*                A = ( A11    A12 ) M, */

/*     then the following does the complete QR decomposition of A11: */

/*              A11 = U*(  0  T12 )*P1**T */
/*                      (  0   0  ) */

#line 419 "dggsvp.f"
    i__1 = *n - *l;
#line 419 "dggsvp.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 420 "dggsvp.f"
	iwork[i__] = 0;
#line 421 "dggsvp.f"
/* L70: */
#line 421 "dggsvp.f"
    }
#line 422 "dggsvp.f"
    i__1 = *n - *l;
#line 422 "dggsvp.f"
    dgeqpf_(m, &i__1, &a[a_offset], lda, &iwork[1], &tau[1], &work[1], info);

/*     Determine the effective rank of A11 */

#line 426 "dggsvp.f"
    *k = 0;
/* Computing MIN */
#line 427 "dggsvp.f"
    i__2 = *m, i__3 = *n - *l;
#line 427 "dggsvp.f"
    i__1 = min(i__2,i__3);
#line 427 "dggsvp.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 428 "dggsvp.f"
	if ((d__1 = a[i__ + i__ * a_dim1], abs(d__1)) > *tola) {
#line 428 "dggsvp.f"
	    ++(*k);
#line 428 "dggsvp.f"
	}
#line 430 "dggsvp.f"
/* L80: */
#line 430 "dggsvp.f"
    }

/*     Update A12 := U**T*A12, where A12 = A( 1:M, N-L+1:N ) */

/* Computing MIN */
#line 434 "dggsvp.f"
    i__2 = *m, i__3 = *n - *l;
#line 434 "dggsvp.f"
    i__1 = min(i__2,i__3);
#line 434 "dggsvp.f"
    dorm2r_("Left", "Transpose", m, l, &i__1, &a[a_offset], lda, &tau[1], &a[(
	    *n - *l + 1) * a_dim1 + 1], lda, &work[1], info, (ftnlen)4, (
	    ftnlen)9);

#line 437 "dggsvp.f"
    if (wantu) {

/*        Copy the details of U, and form U */

#line 441 "dggsvp.f"
	dlaset_("Full", m, m, &c_b12, &c_b12, &u[u_offset], ldu, (ftnlen)4);
#line 442 "dggsvp.f"
	if (*m > 1) {
#line 442 "dggsvp.f"
	    i__1 = *m - 1;
#line 442 "dggsvp.f"
	    i__2 = *n - *l;
#line 442 "dggsvp.f"
	    dlacpy_("Lower", &i__1, &i__2, &a[a_dim1 + 2], lda, &u[u_dim1 + 2]
		    , ldu, (ftnlen)5);
#line 442 "dggsvp.f"
	}
/* Computing MIN */
#line 445 "dggsvp.f"
	i__2 = *m, i__3 = *n - *l;
#line 445 "dggsvp.f"
	i__1 = min(i__2,i__3);
#line 445 "dggsvp.f"
	dorg2r_(m, m, &i__1, &u[u_offset], ldu, &tau[1], &work[1], info);
#line 446 "dggsvp.f"
    }

#line 448 "dggsvp.f"
    if (wantq) {

/*        Update Q( 1:N, 1:N-L )  = Q( 1:N, 1:N-L )*P1 */

#line 452 "dggsvp.f"
	i__1 = *n - *l;
#line 452 "dggsvp.f"
	dlapmt_(&forwrd, n, &i__1, &q[q_offset], ldq, &iwork[1]);
#line 453 "dggsvp.f"
    }

/*     Clean up A: set the strictly lower triangular part of */
/*     A(1:K, 1:K) = 0, and A( K+1:M, 1:N-L ) = 0. */

#line 458 "dggsvp.f"
    i__1 = *k - 1;
#line 458 "dggsvp.f"
    for (j = 1; j <= i__1; ++j) {
#line 459 "dggsvp.f"
	i__2 = *k;
#line 459 "dggsvp.f"
	for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 460 "dggsvp.f"
	    a[i__ + j * a_dim1] = 0.;
#line 461 "dggsvp.f"
/* L90: */
#line 461 "dggsvp.f"
	}
#line 462 "dggsvp.f"
/* L100: */
#line 462 "dggsvp.f"
    }
#line 463 "dggsvp.f"
    if (*m > *k) {
#line 463 "dggsvp.f"
	i__1 = *m - *k;
#line 463 "dggsvp.f"
	i__2 = *n - *l;
#line 463 "dggsvp.f"
	dlaset_("Full", &i__1, &i__2, &c_b12, &c_b12, &a[*k + 1 + a_dim1], 
		lda, (ftnlen)4);
#line 463 "dggsvp.f"
    }

#line 466 "dggsvp.f"
    if (*n - *l > *k) {

/*        RQ factorization of ( T11 T12 ) = ( 0 T12 )*Z1 */

#line 470 "dggsvp.f"
	i__1 = *n - *l;
#line 470 "dggsvp.f"
	dgerq2_(k, &i__1, &a[a_offset], lda, &tau[1], &work[1], info);

#line 472 "dggsvp.f"
	if (wantq) {

/*           Update Q( 1:N,1:N-L ) = Q( 1:N,1:N-L )*Z1**T */

#line 476 "dggsvp.f"
	    i__1 = *n - *l;
#line 476 "dggsvp.f"
	    dormr2_("Right", "Transpose", n, &i__1, k, &a[a_offset], lda, &
		    tau[1], &q[q_offset], ldq, &work[1], info, (ftnlen)5, (
		    ftnlen)9);
#line 478 "dggsvp.f"
	}

/*        Clean up A */

#line 482 "dggsvp.f"
	i__1 = *n - *l - *k;
#line 482 "dggsvp.f"
	dlaset_("Full", k, &i__1, &c_b12, &c_b12, &a[a_offset], lda, (ftnlen)
		4);
#line 483 "dggsvp.f"
	i__1 = *n - *l;
#line 483 "dggsvp.f"
	for (j = *n - *l - *k + 1; j <= i__1; ++j) {
#line 484 "dggsvp.f"
	    i__2 = *k;
#line 484 "dggsvp.f"
	    for (i__ = j - *n + *l + *k + 1; i__ <= i__2; ++i__) {
#line 485 "dggsvp.f"
		a[i__ + j * a_dim1] = 0.;
#line 486 "dggsvp.f"
/* L110: */
#line 486 "dggsvp.f"
	    }
#line 487 "dggsvp.f"
/* L120: */
#line 487 "dggsvp.f"
	}

#line 489 "dggsvp.f"
    }

#line 491 "dggsvp.f"
    if (*m > *k) {

/*        QR factorization of A( K+1:M,N-L+1:N ) */

#line 495 "dggsvp.f"
	i__1 = *m - *k;
#line 495 "dggsvp.f"
	dgeqr2_(&i__1, l, &a[*k + 1 + (*n - *l + 1) * a_dim1], lda, &tau[1], &
		work[1], info);

#line 497 "dggsvp.f"
	if (wantu) {

/*           Update U(:,K+1:M) := U(:,K+1:M)*U1 */

#line 501 "dggsvp.f"
	    i__1 = *m - *k;
/* Computing MIN */
#line 501 "dggsvp.f"
	    i__3 = *m - *k;
#line 501 "dggsvp.f"
	    i__2 = min(i__3,*l);
#line 501 "dggsvp.f"
	    dorm2r_("Right", "No transpose", m, &i__1, &i__2, &a[*k + 1 + (*n 
		    - *l + 1) * a_dim1], lda, &tau[1], &u[(*k + 1) * u_dim1 + 
		    1], ldu, &work[1], info, (ftnlen)5, (ftnlen)12);
#line 504 "dggsvp.f"
	}

/*        Clean up */

#line 508 "dggsvp.f"
	i__1 = *n;
#line 508 "dggsvp.f"
	for (j = *n - *l + 1; j <= i__1; ++j) {
#line 509 "dggsvp.f"
	    i__2 = *m;
#line 509 "dggsvp.f"
	    for (i__ = j - *n + *k + *l + 1; i__ <= i__2; ++i__) {
#line 510 "dggsvp.f"
		a[i__ + j * a_dim1] = 0.;
#line 511 "dggsvp.f"
/* L130: */
#line 511 "dggsvp.f"
	    }
#line 512 "dggsvp.f"
/* L140: */
#line 512 "dggsvp.f"
	}

#line 514 "dggsvp.f"
    }

#line 516 "dggsvp.f"
    return 0;

/*     End of DGGSVP */

} /* dggsvp_ */

