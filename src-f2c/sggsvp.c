#line 1 "sggsvp.f"
/* sggsvp.f -- translated by f2c (version 20100827).
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

#line 1 "sggsvp.f"
/* Table of constant values */

static doublereal c_b12 = 0.;
static doublereal c_b22 = 1.;

/* > \brief \b SGGSVP */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SGGSVP + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sggsvp.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sggsvp.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sggsvp.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SGGSVP( JOBU, JOBV, JOBQ, M, P, N, A, LDA, B, LDB, */
/*                          TOLA, TOLB, K, L, U, LDU, V, LDV, Q, LDQ, */
/*                          IWORK, TAU, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBQ, JOBU, JOBV */
/*       INTEGER            INFO, K, L, LDA, LDB, LDQ, LDU, LDV, M, N, P */
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
/* > SGGSVP computes orthogonal matrices U, V and Q such that */
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
/* > SGGSVD. */
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
/* >          WORK is REAL array, dimension (max(3*N,M,P)) */
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

/* > \ingroup realOTHERcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* >  The subroutine uses LAPACK subroutine SGEQPF for the QR factorization */
/* >  with column pivoting to detect the effective numerical rank of the */
/* >  a matrix. It may be replaced by a better rank determination strategy. */
/* > */
/*  ===================================================================== */
/* Subroutine */ int sggsvp_(char *jobu, char *jobv, char *jobq, integer *m, 
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
    extern /* Subroutine */ int sgeqr2_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *), sgerq2_(
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *), sorg2r_(integer *, integer *, integer *,
	     doublereal *, integer *, doublereal *, doublereal *, integer *), 
	    sorm2r_(char *, char *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen), sormr2_(char *, char *, 
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen, ftnlen), xerbla_(char *, integer *, ftnlen), sgeqpf_(
	    integer *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, integer *), slacpy_(char *, integer *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen), slaset_(char *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, ftnlen), slapmt_(logical *,
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

#line 298 "sggsvp.f"
    /* Parameter adjustments */
#line 298 "sggsvp.f"
    a_dim1 = *lda;
#line 298 "sggsvp.f"
    a_offset = 1 + a_dim1;
#line 298 "sggsvp.f"
    a -= a_offset;
#line 298 "sggsvp.f"
    b_dim1 = *ldb;
#line 298 "sggsvp.f"
    b_offset = 1 + b_dim1;
#line 298 "sggsvp.f"
    b -= b_offset;
#line 298 "sggsvp.f"
    u_dim1 = *ldu;
#line 298 "sggsvp.f"
    u_offset = 1 + u_dim1;
#line 298 "sggsvp.f"
    u -= u_offset;
#line 298 "sggsvp.f"
    v_dim1 = *ldv;
#line 298 "sggsvp.f"
    v_offset = 1 + v_dim1;
#line 298 "sggsvp.f"
    v -= v_offset;
#line 298 "sggsvp.f"
    q_dim1 = *ldq;
#line 298 "sggsvp.f"
    q_offset = 1 + q_dim1;
#line 298 "sggsvp.f"
    q -= q_offset;
#line 298 "sggsvp.f"
    --iwork;
#line 298 "sggsvp.f"
    --tau;
#line 298 "sggsvp.f"
    --work;
#line 298 "sggsvp.f"

#line 298 "sggsvp.f"
    /* Function Body */
#line 298 "sggsvp.f"
    wantu = lsame_(jobu, "U", (ftnlen)1, (ftnlen)1);
#line 299 "sggsvp.f"
    wantv = lsame_(jobv, "V", (ftnlen)1, (ftnlen)1);
#line 300 "sggsvp.f"
    wantq = lsame_(jobq, "Q", (ftnlen)1, (ftnlen)1);
#line 301 "sggsvp.f"
    forwrd = TRUE_;

#line 303 "sggsvp.f"
    *info = 0;
#line 304 "sggsvp.f"
    if (! (wantu || lsame_(jobu, "N", (ftnlen)1, (ftnlen)1))) {
#line 305 "sggsvp.f"
	*info = -1;
#line 306 "sggsvp.f"
    } else if (! (wantv || lsame_(jobv, "N", (ftnlen)1, (ftnlen)1))) {
#line 307 "sggsvp.f"
	*info = -2;
#line 308 "sggsvp.f"
    } else if (! (wantq || lsame_(jobq, "N", (ftnlen)1, (ftnlen)1))) {
#line 309 "sggsvp.f"
	*info = -3;
#line 310 "sggsvp.f"
    } else if (*m < 0) {
#line 311 "sggsvp.f"
	*info = -4;
#line 312 "sggsvp.f"
    } else if (*p < 0) {
#line 313 "sggsvp.f"
	*info = -5;
#line 314 "sggsvp.f"
    } else if (*n < 0) {
#line 315 "sggsvp.f"
	*info = -6;
#line 316 "sggsvp.f"
    } else if (*lda < max(1,*m)) {
#line 317 "sggsvp.f"
	*info = -8;
#line 318 "sggsvp.f"
    } else if (*ldb < max(1,*p)) {
#line 319 "sggsvp.f"
	*info = -10;
#line 320 "sggsvp.f"
    } else if (*ldu < 1 || wantu && *ldu < *m) {
#line 321 "sggsvp.f"
	*info = -16;
#line 322 "sggsvp.f"
    } else if (*ldv < 1 || wantv && *ldv < *p) {
#line 323 "sggsvp.f"
	*info = -18;
#line 324 "sggsvp.f"
    } else if (*ldq < 1 || wantq && *ldq < *n) {
#line 325 "sggsvp.f"
	*info = -20;
#line 326 "sggsvp.f"
    }
#line 327 "sggsvp.f"
    if (*info != 0) {
#line 328 "sggsvp.f"
	i__1 = -(*info);
#line 328 "sggsvp.f"
	xerbla_("SGGSVP", &i__1, (ftnlen)6);
#line 329 "sggsvp.f"
	return 0;
#line 330 "sggsvp.f"
    }

/*     QR with column pivoting of B: B*P = V*( S11 S12 ) */
/*                                           (  0   0  ) */

#line 335 "sggsvp.f"
    i__1 = *n;
#line 335 "sggsvp.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 336 "sggsvp.f"
	iwork[i__] = 0;
#line 337 "sggsvp.f"
/* L10: */
#line 337 "sggsvp.f"
    }
#line 338 "sggsvp.f"
    sgeqpf_(p, n, &b[b_offset], ldb, &iwork[1], &tau[1], &work[1], info);

/*     Update A := A*P */

#line 342 "sggsvp.f"
    slapmt_(&forwrd, m, n, &a[a_offset], lda, &iwork[1]);

/*     Determine the effective rank of matrix B. */

#line 346 "sggsvp.f"
    *l = 0;
#line 347 "sggsvp.f"
    i__1 = min(*p,*n);
#line 347 "sggsvp.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 348 "sggsvp.f"
	if ((d__1 = b[i__ + i__ * b_dim1], abs(d__1)) > *tolb) {
#line 348 "sggsvp.f"
	    ++(*l);
#line 348 "sggsvp.f"
	}
#line 350 "sggsvp.f"
/* L20: */
#line 350 "sggsvp.f"
    }

#line 352 "sggsvp.f"
    if (wantv) {

/*        Copy the details of V, and form V. */

#line 356 "sggsvp.f"
	slaset_("Full", p, p, &c_b12, &c_b12, &v[v_offset], ldv, (ftnlen)4);
#line 357 "sggsvp.f"
	if (*p > 1) {
#line 357 "sggsvp.f"
	    i__1 = *p - 1;
#line 357 "sggsvp.f"
	    slacpy_("Lower", &i__1, n, &b[b_dim1 + 2], ldb, &v[v_dim1 + 2], 
		    ldv, (ftnlen)5);
#line 357 "sggsvp.f"
	}
#line 360 "sggsvp.f"
	i__1 = min(*p,*n);
#line 360 "sggsvp.f"
	sorg2r_(p, p, &i__1, &v[v_offset], ldv, &tau[1], &work[1], info);
#line 361 "sggsvp.f"
    }

/*     Clean up B */

#line 365 "sggsvp.f"
    i__1 = *l - 1;
#line 365 "sggsvp.f"
    for (j = 1; j <= i__1; ++j) {
#line 366 "sggsvp.f"
	i__2 = *l;
#line 366 "sggsvp.f"
	for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 367 "sggsvp.f"
	    b[i__ + j * b_dim1] = 0.;
#line 368 "sggsvp.f"
/* L30: */
#line 368 "sggsvp.f"
	}
#line 369 "sggsvp.f"
/* L40: */
#line 369 "sggsvp.f"
    }
#line 370 "sggsvp.f"
    if (*p > *l) {
#line 370 "sggsvp.f"
	i__1 = *p - *l;
#line 370 "sggsvp.f"
	slaset_("Full", &i__1, n, &c_b12, &c_b12, &b[*l + 1 + b_dim1], ldb, (
		ftnlen)4);
#line 370 "sggsvp.f"
    }

#line 373 "sggsvp.f"
    if (wantq) {

/*        Set Q = I and Update Q := Q*P */

#line 377 "sggsvp.f"
	slaset_("Full", n, n, &c_b12, &c_b22, &q[q_offset], ldq, (ftnlen)4);
#line 378 "sggsvp.f"
	slapmt_(&forwrd, n, n, &q[q_offset], ldq, &iwork[1]);
#line 379 "sggsvp.f"
    }

#line 381 "sggsvp.f"
    if (*p >= *l && *n != *l) {

/*        RQ factorization of (S11 S12): ( S11 S12 ) = ( 0 S12 )*Z */

#line 385 "sggsvp.f"
	sgerq2_(l, n, &b[b_offset], ldb, &tau[1], &work[1], info);

/*        Update A := A*Z**T */

#line 389 "sggsvp.f"
	sormr2_("Right", "Transpose", m, n, l, &b[b_offset], ldb, &tau[1], &a[
		a_offset], lda, &work[1], info, (ftnlen)5, (ftnlen)9);

#line 392 "sggsvp.f"
	if (wantq) {

/*           Update Q := Q*Z**T */

#line 396 "sggsvp.f"
	    sormr2_("Right", "Transpose", n, n, l, &b[b_offset], ldb, &tau[1],
		     &q[q_offset], ldq, &work[1], info, (ftnlen)5, (ftnlen)9);
#line 398 "sggsvp.f"
	}

/*        Clean up B */

#line 402 "sggsvp.f"
	i__1 = *n - *l;
#line 402 "sggsvp.f"
	slaset_("Full", l, &i__1, &c_b12, &c_b12, &b[b_offset], ldb, (ftnlen)
		4);
#line 403 "sggsvp.f"
	i__1 = *n;
#line 403 "sggsvp.f"
	for (j = *n - *l + 1; j <= i__1; ++j) {
#line 404 "sggsvp.f"
	    i__2 = *l;
#line 404 "sggsvp.f"
	    for (i__ = j - *n + *l + 1; i__ <= i__2; ++i__) {
#line 405 "sggsvp.f"
		b[i__ + j * b_dim1] = 0.;
#line 406 "sggsvp.f"
/* L50: */
#line 406 "sggsvp.f"
	    }
#line 407 "sggsvp.f"
/* L60: */
#line 407 "sggsvp.f"
	}

#line 409 "sggsvp.f"
    }

/*     Let              N-L     L */
/*                A = ( A11    A12 ) M, */

/*     then the following does the complete QR decomposition of A11: */

/*              A11 = U*(  0  T12 )*P1**T */
/*                      (  0   0  ) */

#line 419 "sggsvp.f"
    i__1 = *n - *l;
#line 419 "sggsvp.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 420 "sggsvp.f"
	iwork[i__] = 0;
#line 421 "sggsvp.f"
/* L70: */
#line 421 "sggsvp.f"
    }
#line 422 "sggsvp.f"
    i__1 = *n - *l;
#line 422 "sggsvp.f"
    sgeqpf_(m, &i__1, &a[a_offset], lda, &iwork[1], &tau[1], &work[1], info);

/*     Determine the effective rank of A11 */

#line 426 "sggsvp.f"
    *k = 0;
/* Computing MIN */
#line 427 "sggsvp.f"
    i__2 = *m, i__3 = *n - *l;
#line 427 "sggsvp.f"
    i__1 = min(i__2,i__3);
#line 427 "sggsvp.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 428 "sggsvp.f"
	if ((d__1 = a[i__ + i__ * a_dim1], abs(d__1)) > *tola) {
#line 428 "sggsvp.f"
	    ++(*k);
#line 428 "sggsvp.f"
	}
#line 430 "sggsvp.f"
/* L80: */
#line 430 "sggsvp.f"
    }

/*     Update A12 := U**T*A12, where A12 = A( 1:M, N-L+1:N ) */

/* Computing MIN */
#line 434 "sggsvp.f"
    i__2 = *m, i__3 = *n - *l;
#line 434 "sggsvp.f"
    i__1 = min(i__2,i__3);
#line 434 "sggsvp.f"
    sorm2r_("Left", "Transpose", m, l, &i__1, &a[a_offset], lda, &tau[1], &a[(
	    *n - *l + 1) * a_dim1 + 1], lda, &work[1], info, (ftnlen)4, (
	    ftnlen)9);

#line 437 "sggsvp.f"
    if (wantu) {

/*        Copy the details of U, and form U */

#line 441 "sggsvp.f"
	slaset_("Full", m, m, &c_b12, &c_b12, &u[u_offset], ldu, (ftnlen)4);
#line 442 "sggsvp.f"
	if (*m > 1) {
#line 442 "sggsvp.f"
	    i__1 = *m - 1;
#line 442 "sggsvp.f"
	    i__2 = *n - *l;
#line 442 "sggsvp.f"
	    slacpy_("Lower", &i__1, &i__2, &a[a_dim1 + 2], lda, &u[u_dim1 + 2]
		    , ldu, (ftnlen)5);
#line 442 "sggsvp.f"
	}
/* Computing MIN */
#line 445 "sggsvp.f"
	i__2 = *m, i__3 = *n - *l;
#line 445 "sggsvp.f"
	i__1 = min(i__2,i__3);
#line 445 "sggsvp.f"
	sorg2r_(m, m, &i__1, &u[u_offset], ldu, &tau[1], &work[1], info);
#line 446 "sggsvp.f"
    }

#line 448 "sggsvp.f"
    if (wantq) {

/*        Update Q( 1:N, 1:N-L )  = Q( 1:N, 1:N-L )*P1 */

#line 452 "sggsvp.f"
	i__1 = *n - *l;
#line 452 "sggsvp.f"
	slapmt_(&forwrd, n, &i__1, &q[q_offset], ldq, &iwork[1]);
#line 453 "sggsvp.f"
    }

/*     Clean up A: set the strictly lower triangular part of */
/*     A(1:K, 1:K) = 0, and A( K+1:M, 1:N-L ) = 0. */

#line 458 "sggsvp.f"
    i__1 = *k - 1;
#line 458 "sggsvp.f"
    for (j = 1; j <= i__1; ++j) {
#line 459 "sggsvp.f"
	i__2 = *k;
#line 459 "sggsvp.f"
	for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 460 "sggsvp.f"
	    a[i__ + j * a_dim1] = 0.;
#line 461 "sggsvp.f"
/* L90: */
#line 461 "sggsvp.f"
	}
#line 462 "sggsvp.f"
/* L100: */
#line 462 "sggsvp.f"
    }
#line 463 "sggsvp.f"
    if (*m > *k) {
#line 463 "sggsvp.f"
	i__1 = *m - *k;
#line 463 "sggsvp.f"
	i__2 = *n - *l;
#line 463 "sggsvp.f"
	slaset_("Full", &i__1, &i__2, &c_b12, &c_b12, &a[*k + 1 + a_dim1], 
		lda, (ftnlen)4);
#line 463 "sggsvp.f"
    }

#line 466 "sggsvp.f"
    if (*n - *l > *k) {

/*        RQ factorization of ( T11 T12 ) = ( 0 T12 )*Z1 */

#line 470 "sggsvp.f"
	i__1 = *n - *l;
#line 470 "sggsvp.f"
	sgerq2_(k, &i__1, &a[a_offset], lda, &tau[1], &work[1], info);

#line 472 "sggsvp.f"
	if (wantq) {

/*           Update Q( 1:N,1:N-L ) = Q( 1:N,1:N-L )*Z1**T */

#line 476 "sggsvp.f"
	    i__1 = *n - *l;
#line 476 "sggsvp.f"
	    sormr2_("Right", "Transpose", n, &i__1, k, &a[a_offset], lda, &
		    tau[1], &q[q_offset], ldq, &work[1], info, (ftnlen)5, (
		    ftnlen)9);
#line 478 "sggsvp.f"
	}

/*        Clean up A */

#line 482 "sggsvp.f"
	i__1 = *n - *l - *k;
#line 482 "sggsvp.f"
	slaset_("Full", k, &i__1, &c_b12, &c_b12, &a[a_offset], lda, (ftnlen)
		4);
#line 483 "sggsvp.f"
	i__1 = *n - *l;
#line 483 "sggsvp.f"
	for (j = *n - *l - *k + 1; j <= i__1; ++j) {
#line 484 "sggsvp.f"
	    i__2 = *k;
#line 484 "sggsvp.f"
	    for (i__ = j - *n + *l + *k + 1; i__ <= i__2; ++i__) {
#line 485 "sggsvp.f"
		a[i__ + j * a_dim1] = 0.;
#line 486 "sggsvp.f"
/* L110: */
#line 486 "sggsvp.f"
	    }
#line 487 "sggsvp.f"
/* L120: */
#line 487 "sggsvp.f"
	}

#line 489 "sggsvp.f"
    }

#line 491 "sggsvp.f"
    if (*m > *k) {

/*        QR factorization of A( K+1:M,N-L+1:N ) */

#line 495 "sggsvp.f"
	i__1 = *m - *k;
#line 495 "sggsvp.f"
	sgeqr2_(&i__1, l, &a[*k + 1 + (*n - *l + 1) * a_dim1], lda, &tau[1], &
		work[1], info);

#line 497 "sggsvp.f"
	if (wantu) {

/*           Update U(:,K+1:M) := U(:,K+1:M)*U1 */

#line 501 "sggsvp.f"
	    i__1 = *m - *k;
/* Computing MIN */
#line 501 "sggsvp.f"
	    i__3 = *m - *k;
#line 501 "sggsvp.f"
	    i__2 = min(i__3,*l);
#line 501 "sggsvp.f"
	    sorm2r_("Right", "No transpose", m, &i__1, &i__2, &a[*k + 1 + (*n 
		    - *l + 1) * a_dim1], lda, &tau[1], &u[(*k + 1) * u_dim1 + 
		    1], ldu, &work[1], info, (ftnlen)5, (ftnlen)12);
#line 504 "sggsvp.f"
	}

/*        Clean up */

#line 508 "sggsvp.f"
	i__1 = *n;
#line 508 "sggsvp.f"
	for (j = *n - *l + 1; j <= i__1; ++j) {
#line 509 "sggsvp.f"
	    i__2 = *m;
#line 509 "sggsvp.f"
	    for (i__ = j - *n + *k + *l + 1; i__ <= i__2; ++i__) {
#line 510 "sggsvp.f"
		a[i__ + j * a_dim1] = 0.;
#line 511 "sggsvp.f"
/* L130: */
#line 511 "sggsvp.f"
	    }
#line 512 "sggsvp.f"
/* L140: */
#line 512 "sggsvp.f"
	}

#line 514 "sggsvp.f"
    }

#line 516 "sggsvp.f"
    return 0;

/*     End of SGGSVP */

} /* sggsvp_ */

