#line 1 "cggsvp.f"
/* cggsvp.f -- translated by f2c (version 20100827).
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

#line 1 "cggsvp.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};

/* > \brief \b CGGSVP */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CGGSVP + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cggsvp.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cggsvp.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cggsvp.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CGGSVP( JOBU, JOBV, JOBQ, M, P, N, A, LDA, B, LDB, */
/*                          TOLA, TOLB, K, L, U, LDU, V, LDV, Q, LDQ, */
/*                          IWORK, RWORK, TAU, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBQ, JOBU, JOBV */
/*       INTEGER            INFO, K, L, LDA, LDB, LDQ, LDU, LDV, M, N, P */
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
/* > CGGSVP computes unitary matrices U, V and Q such that */
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
/* > CGGSVD. */
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
/* >          WORK is COMPLEX array, dimension (max(3*N,M,P)) */
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

/* > \ingroup complexOTHERcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* >  The subroutine uses LAPACK subroutine CGEQPF for the QR factorization */
/* >  with column pivoting to detect the effective numerical rank of the */
/* >  a matrix. It may be replaced by a better rank determination strategy. */
/* > */
/*  ===================================================================== */
/* Subroutine */ int cggsvp_(char *jobu, char *jobv, char *jobq, integer *m, 
	integer *p, integer *n, doublecomplex *a, integer *lda, doublecomplex 
	*b, integer *ldb, doublereal *tola, doublereal *tolb, integer *k, 
	integer *l, doublecomplex *u, integer *ldu, doublecomplex *v, integer 
	*ldv, doublecomplex *q, integer *ldq, integer *iwork, doublereal *
	rwork, doublecomplex *tau, doublecomplex *work, integer *info, ftnlen 
	jobu_len, ftnlen jobv_len, ftnlen jobq_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, q_dim1, q_offset, u_dim1, 
	    u_offset, v_dim1, v_offset, i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Builtin functions */
    double d_imag(doublecomplex *);

    /* Local variables */
    static integer i__, j;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical wantq, wantu, wantv;
    extern /* Subroutine */ int cgeqr2_(integer *, integer *, doublecomplex *,
	     integer *, doublecomplex *, doublecomplex *, integer *), cgerq2_(
	    integer *, integer *, doublecomplex *, integer *, doublecomplex *,
	     doublecomplex *, integer *), cung2r_(integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *), cunm2r_(char *, char *, integer *, 
	    integer *, integer *, doublecomplex *, integer *, doublecomplex *,
	     doublecomplex *, integer *, doublecomplex *, integer *, ftnlen, 
	    ftnlen), cunmr2_(char *, char *, integer *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, ftnlen, ftnlen), cgeqpf_(
	    integer *, integer *, doublecomplex *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, doublereal *, integer *), 
	    clacpy_(char *, integer *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, ftnlen), claset_(char *, integer *, 
	    integer *, doublecomplex *, doublecomplex *, doublecomplex *, 
	    integer *, ftnlen), xerbla_(char *, integer *, ftnlen), clapmt_(
	    logical *, integer *, integer *, doublecomplex *, integer *, 
	    integer *);
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
/*     .. Statement Functions .. */
/*     .. */
/*     .. Statement Function definitions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters */

#line 313 "cggsvp.f"
    /* Parameter adjustments */
#line 313 "cggsvp.f"
    a_dim1 = *lda;
#line 313 "cggsvp.f"
    a_offset = 1 + a_dim1;
#line 313 "cggsvp.f"
    a -= a_offset;
#line 313 "cggsvp.f"
    b_dim1 = *ldb;
#line 313 "cggsvp.f"
    b_offset = 1 + b_dim1;
#line 313 "cggsvp.f"
    b -= b_offset;
#line 313 "cggsvp.f"
    u_dim1 = *ldu;
#line 313 "cggsvp.f"
    u_offset = 1 + u_dim1;
#line 313 "cggsvp.f"
    u -= u_offset;
#line 313 "cggsvp.f"
    v_dim1 = *ldv;
#line 313 "cggsvp.f"
    v_offset = 1 + v_dim1;
#line 313 "cggsvp.f"
    v -= v_offset;
#line 313 "cggsvp.f"
    q_dim1 = *ldq;
#line 313 "cggsvp.f"
    q_offset = 1 + q_dim1;
#line 313 "cggsvp.f"
    q -= q_offset;
#line 313 "cggsvp.f"
    --iwork;
#line 313 "cggsvp.f"
    --rwork;
#line 313 "cggsvp.f"
    --tau;
#line 313 "cggsvp.f"
    --work;
#line 313 "cggsvp.f"

#line 313 "cggsvp.f"
    /* Function Body */
#line 313 "cggsvp.f"
    wantu = lsame_(jobu, "U", (ftnlen)1, (ftnlen)1);
#line 314 "cggsvp.f"
    wantv = lsame_(jobv, "V", (ftnlen)1, (ftnlen)1);
#line 315 "cggsvp.f"
    wantq = lsame_(jobq, "Q", (ftnlen)1, (ftnlen)1);
#line 316 "cggsvp.f"
    forwrd = TRUE_;

#line 318 "cggsvp.f"
    *info = 0;
#line 319 "cggsvp.f"
    if (! (wantu || lsame_(jobu, "N", (ftnlen)1, (ftnlen)1))) {
#line 320 "cggsvp.f"
	*info = -1;
#line 321 "cggsvp.f"
    } else if (! (wantv || lsame_(jobv, "N", (ftnlen)1, (ftnlen)1))) {
#line 322 "cggsvp.f"
	*info = -2;
#line 323 "cggsvp.f"
    } else if (! (wantq || lsame_(jobq, "N", (ftnlen)1, (ftnlen)1))) {
#line 324 "cggsvp.f"
	*info = -3;
#line 325 "cggsvp.f"
    } else if (*m < 0) {
#line 326 "cggsvp.f"
	*info = -4;
#line 327 "cggsvp.f"
    } else if (*p < 0) {
#line 328 "cggsvp.f"
	*info = -5;
#line 329 "cggsvp.f"
    } else if (*n < 0) {
#line 330 "cggsvp.f"
	*info = -6;
#line 331 "cggsvp.f"
    } else if (*lda < max(1,*m)) {
#line 332 "cggsvp.f"
	*info = -8;
#line 333 "cggsvp.f"
    } else if (*ldb < max(1,*p)) {
#line 334 "cggsvp.f"
	*info = -10;
#line 335 "cggsvp.f"
    } else if (*ldu < 1 || wantu && *ldu < *m) {
#line 336 "cggsvp.f"
	*info = -16;
#line 337 "cggsvp.f"
    } else if (*ldv < 1 || wantv && *ldv < *p) {
#line 338 "cggsvp.f"
	*info = -18;
#line 339 "cggsvp.f"
    } else if (*ldq < 1 || wantq && *ldq < *n) {
#line 340 "cggsvp.f"
	*info = -20;
#line 341 "cggsvp.f"
    }
#line 342 "cggsvp.f"
    if (*info != 0) {
#line 343 "cggsvp.f"
	i__1 = -(*info);
#line 343 "cggsvp.f"
	xerbla_("CGGSVP", &i__1, (ftnlen)6);
#line 344 "cggsvp.f"
	return 0;
#line 345 "cggsvp.f"
    }

/*     QR with column pivoting of B: B*P = V*( S11 S12 ) */
/*                                           (  0   0  ) */

#line 350 "cggsvp.f"
    i__1 = *n;
#line 350 "cggsvp.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 351 "cggsvp.f"
	iwork[i__] = 0;
#line 352 "cggsvp.f"
/* L10: */
#line 352 "cggsvp.f"
    }
#line 353 "cggsvp.f"
    cgeqpf_(p, n, &b[b_offset], ldb, &iwork[1], &tau[1], &work[1], &rwork[1], 
	    info);

/*     Update A := A*P */

#line 357 "cggsvp.f"
    clapmt_(&forwrd, m, n, &a[a_offset], lda, &iwork[1]);

/*     Determine the effective rank of matrix B. */

#line 361 "cggsvp.f"
    *l = 0;
#line 362 "cggsvp.f"
    i__1 = min(*p,*n);
#line 362 "cggsvp.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 363 "cggsvp.f"
	i__2 = i__ + i__ * b_dim1;
#line 363 "cggsvp.f"
	if ((d__1 = b[i__2].r, abs(d__1)) + (d__2 = d_imag(&b[i__ + i__ * 
		b_dim1]), abs(d__2)) > *tolb) {
#line 363 "cggsvp.f"
	    ++(*l);
#line 363 "cggsvp.f"
	}
#line 365 "cggsvp.f"
/* L20: */
#line 365 "cggsvp.f"
    }

#line 367 "cggsvp.f"
    if (wantv) {

/*        Copy the details of V, and form V. */

#line 371 "cggsvp.f"
	claset_("Full", p, p, &c_b1, &c_b1, &v[v_offset], ldv, (ftnlen)4);
#line 372 "cggsvp.f"
	if (*p > 1) {
#line 372 "cggsvp.f"
	    i__1 = *p - 1;
#line 372 "cggsvp.f"
	    clacpy_("Lower", &i__1, n, &b[b_dim1 + 2], ldb, &v[v_dim1 + 2], 
		    ldv, (ftnlen)5);
#line 372 "cggsvp.f"
	}
#line 375 "cggsvp.f"
	i__1 = min(*p,*n);
#line 375 "cggsvp.f"
	cung2r_(p, p, &i__1, &v[v_offset], ldv, &tau[1], &work[1], info);
#line 376 "cggsvp.f"
    }

/*     Clean up B */

#line 380 "cggsvp.f"
    i__1 = *l - 1;
#line 380 "cggsvp.f"
    for (j = 1; j <= i__1; ++j) {
#line 381 "cggsvp.f"
	i__2 = *l;
#line 381 "cggsvp.f"
	for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 382 "cggsvp.f"
	    i__3 = i__ + j * b_dim1;
#line 382 "cggsvp.f"
	    b[i__3].r = 0., b[i__3].i = 0.;
#line 383 "cggsvp.f"
/* L30: */
#line 383 "cggsvp.f"
	}
#line 384 "cggsvp.f"
/* L40: */
#line 384 "cggsvp.f"
    }
#line 385 "cggsvp.f"
    if (*p > *l) {
#line 385 "cggsvp.f"
	i__1 = *p - *l;
#line 385 "cggsvp.f"
	claset_("Full", &i__1, n, &c_b1, &c_b1, &b[*l + 1 + b_dim1], ldb, (
		ftnlen)4);
#line 385 "cggsvp.f"
    }

#line 388 "cggsvp.f"
    if (wantq) {

/*        Set Q = I and Update Q := Q*P */

#line 392 "cggsvp.f"
	claset_("Full", n, n, &c_b1, &c_b2, &q[q_offset], ldq, (ftnlen)4);
#line 393 "cggsvp.f"
	clapmt_(&forwrd, n, n, &q[q_offset], ldq, &iwork[1]);
#line 394 "cggsvp.f"
    }

#line 396 "cggsvp.f"
    if (*p >= *l && *n != *l) {

/*        RQ factorization of ( S11 S12 ) = ( 0 S12 )*Z */

#line 400 "cggsvp.f"
	cgerq2_(l, n, &b[b_offset], ldb, &tau[1], &work[1], info);

/*        Update A := A*Z**H */

#line 404 "cggsvp.f"
	cunmr2_("Right", "Conjugate transpose", m, n, l, &b[b_offset], ldb, &
		tau[1], &a[a_offset], lda, &work[1], info, (ftnlen)5, (ftnlen)
		19);
#line 406 "cggsvp.f"
	if (wantq) {

/*           Update Q := Q*Z**H */

#line 410 "cggsvp.f"
	    cunmr2_("Right", "Conjugate transpose", n, n, l, &b[b_offset], 
		    ldb, &tau[1], &q[q_offset], ldq, &work[1], info, (ftnlen)
		    5, (ftnlen)19);
#line 412 "cggsvp.f"
	}

/*        Clean up B */

#line 416 "cggsvp.f"
	i__1 = *n - *l;
#line 416 "cggsvp.f"
	claset_("Full", l, &i__1, &c_b1, &c_b1, &b[b_offset], ldb, (ftnlen)4);
#line 417 "cggsvp.f"
	i__1 = *n;
#line 417 "cggsvp.f"
	for (j = *n - *l + 1; j <= i__1; ++j) {
#line 418 "cggsvp.f"
	    i__2 = *l;
#line 418 "cggsvp.f"
	    for (i__ = j - *n + *l + 1; i__ <= i__2; ++i__) {
#line 419 "cggsvp.f"
		i__3 = i__ + j * b_dim1;
#line 419 "cggsvp.f"
		b[i__3].r = 0., b[i__3].i = 0.;
#line 420 "cggsvp.f"
/* L50: */
#line 420 "cggsvp.f"
	    }
#line 421 "cggsvp.f"
/* L60: */
#line 421 "cggsvp.f"
	}

#line 423 "cggsvp.f"
    }

/*     Let              N-L     L */
/*                A = ( A11    A12 ) M, */

/*     then the following does the complete QR decomposition of A11: */

/*              A11 = U*(  0  T12 )*P1**H */
/*                      (  0   0  ) */

#line 433 "cggsvp.f"
    i__1 = *n - *l;
#line 433 "cggsvp.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 434 "cggsvp.f"
	iwork[i__] = 0;
#line 435 "cggsvp.f"
/* L70: */
#line 435 "cggsvp.f"
    }
#line 436 "cggsvp.f"
    i__1 = *n - *l;
#line 436 "cggsvp.f"
    cgeqpf_(m, &i__1, &a[a_offset], lda, &iwork[1], &tau[1], &work[1], &rwork[
	    1], info);

/*     Determine the effective rank of A11 */

#line 440 "cggsvp.f"
    *k = 0;
/* Computing MIN */
#line 441 "cggsvp.f"
    i__2 = *m, i__3 = *n - *l;
#line 441 "cggsvp.f"
    i__1 = min(i__2,i__3);
#line 441 "cggsvp.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 442 "cggsvp.f"
	i__2 = i__ + i__ * a_dim1;
#line 442 "cggsvp.f"
	if ((d__1 = a[i__2].r, abs(d__1)) + (d__2 = d_imag(&a[i__ + i__ * 
		a_dim1]), abs(d__2)) > *tola) {
#line 442 "cggsvp.f"
	    ++(*k);
#line 442 "cggsvp.f"
	}
#line 444 "cggsvp.f"
/* L80: */
#line 444 "cggsvp.f"
    }

/*     Update A12 := U**H*A12, where A12 = A( 1:M, N-L+1:N ) */

/* Computing MIN */
#line 448 "cggsvp.f"
    i__2 = *m, i__3 = *n - *l;
#line 448 "cggsvp.f"
    i__1 = min(i__2,i__3);
#line 448 "cggsvp.f"
    cunm2r_("Left", "Conjugate transpose", m, l, &i__1, &a[a_offset], lda, &
	    tau[1], &a[(*n - *l + 1) * a_dim1 + 1], lda, &work[1], info, (
	    ftnlen)4, (ftnlen)19);

#line 451 "cggsvp.f"
    if (wantu) {

/*        Copy the details of U, and form U */

#line 455 "cggsvp.f"
	claset_("Full", m, m, &c_b1, &c_b1, &u[u_offset], ldu, (ftnlen)4);
#line 456 "cggsvp.f"
	if (*m > 1) {
#line 456 "cggsvp.f"
	    i__1 = *m - 1;
#line 456 "cggsvp.f"
	    i__2 = *n - *l;
#line 456 "cggsvp.f"
	    clacpy_("Lower", &i__1, &i__2, &a[a_dim1 + 2], lda, &u[u_dim1 + 2]
		    , ldu, (ftnlen)5);
#line 456 "cggsvp.f"
	}
/* Computing MIN */
#line 459 "cggsvp.f"
	i__2 = *m, i__3 = *n - *l;
#line 459 "cggsvp.f"
	i__1 = min(i__2,i__3);
#line 459 "cggsvp.f"
	cung2r_(m, m, &i__1, &u[u_offset], ldu, &tau[1], &work[1], info);
#line 460 "cggsvp.f"
    }

#line 462 "cggsvp.f"
    if (wantq) {

/*        Update Q( 1:N, 1:N-L )  = Q( 1:N, 1:N-L )*P1 */

#line 466 "cggsvp.f"
	i__1 = *n - *l;
#line 466 "cggsvp.f"
	clapmt_(&forwrd, n, &i__1, &q[q_offset], ldq, &iwork[1]);
#line 467 "cggsvp.f"
    }

/*     Clean up A: set the strictly lower triangular part of */
/*     A(1:K, 1:K) = 0, and A( K+1:M, 1:N-L ) = 0. */

#line 472 "cggsvp.f"
    i__1 = *k - 1;
#line 472 "cggsvp.f"
    for (j = 1; j <= i__1; ++j) {
#line 473 "cggsvp.f"
	i__2 = *k;
#line 473 "cggsvp.f"
	for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 474 "cggsvp.f"
	    i__3 = i__ + j * a_dim1;
#line 474 "cggsvp.f"
	    a[i__3].r = 0., a[i__3].i = 0.;
#line 475 "cggsvp.f"
/* L90: */
#line 475 "cggsvp.f"
	}
#line 476 "cggsvp.f"
/* L100: */
#line 476 "cggsvp.f"
    }
#line 477 "cggsvp.f"
    if (*m > *k) {
#line 477 "cggsvp.f"
	i__1 = *m - *k;
#line 477 "cggsvp.f"
	i__2 = *n - *l;
#line 477 "cggsvp.f"
	claset_("Full", &i__1, &i__2, &c_b1, &c_b1, &a[*k + 1 + a_dim1], lda, 
		(ftnlen)4);
#line 477 "cggsvp.f"
    }

#line 480 "cggsvp.f"
    if (*n - *l > *k) {

/*        RQ factorization of ( T11 T12 ) = ( 0 T12 )*Z1 */

#line 484 "cggsvp.f"
	i__1 = *n - *l;
#line 484 "cggsvp.f"
	cgerq2_(k, &i__1, &a[a_offset], lda, &tau[1], &work[1], info);

#line 486 "cggsvp.f"
	if (wantq) {

/*           Update Q( 1:N,1:N-L ) = Q( 1:N,1:N-L )*Z1**H */

#line 490 "cggsvp.f"
	    i__1 = *n - *l;
#line 490 "cggsvp.f"
	    cunmr2_("Right", "Conjugate transpose", n, &i__1, k, &a[a_offset],
		     lda, &tau[1], &q[q_offset], ldq, &work[1], info, (ftnlen)
		    5, (ftnlen)19);
#line 492 "cggsvp.f"
	}

/*        Clean up A */

#line 496 "cggsvp.f"
	i__1 = *n - *l - *k;
#line 496 "cggsvp.f"
	claset_("Full", k, &i__1, &c_b1, &c_b1, &a[a_offset], lda, (ftnlen)4);
#line 497 "cggsvp.f"
	i__1 = *n - *l;
#line 497 "cggsvp.f"
	for (j = *n - *l - *k + 1; j <= i__1; ++j) {
#line 498 "cggsvp.f"
	    i__2 = *k;
#line 498 "cggsvp.f"
	    for (i__ = j - *n + *l + *k + 1; i__ <= i__2; ++i__) {
#line 499 "cggsvp.f"
		i__3 = i__ + j * a_dim1;
#line 499 "cggsvp.f"
		a[i__3].r = 0., a[i__3].i = 0.;
#line 500 "cggsvp.f"
/* L110: */
#line 500 "cggsvp.f"
	    }
#line 501 "cggsvp.f"
/* L120: */
#line 501 "cggsvp.f"
	}

#line 503 "cggsvp.f"
    }

#line 505 "cggsvp.f"
    if (*m > *k) {

/*        QR factorization of A( K+1:M,N-L+1:N ) */

#line 509 "cggsvp.f"
	i__1 = *m - *k;
#line 509 "cggsvp.f"
	cgeqr2_(&i__1, l, &a[*k + 1 + (*n - *l + 1) * a_dim1], lda, &tau[1], &
		work[1], info);

#line 511 "cggsvp.f"
	if (wantu) {

/*           Update U(:,K+1:M) := U(:,K+1:M)*U1 */

#line 515 "cggsvp.f"
	    i__1 = *m - *k;
/* Computing MIN */
#line 515 "cggsvp.f"
	    i__3 = *m - *k;
#line 515 "cggsvp.f"
	    i__2 = min(i__3,*l);
#line 515 "cggsvp.f"
	    cunm2r_("Right", "No transpose", m, &i__1, &i__2, &a[*k + 1 + (*n 
		    - *l + 1) * a_dim1], lda, &tau[1], &u[(*k + 1) * u_dim1 + 
		    1], ldu, &work[1], info, (ftnlen)5, (ftnlen)12);
#line 518 "cggsvp.f"
	}

/*        Clean up */

#line 522 "cggsvp.f"
	i__1 = *n;
#line 522 "cggsvp.f"
	for (j = *n - *l + 1; j <= i__1; ++j) {
#line 523 "cggsvp.f"
	    i__2 = *m;
#line 523 "cggsvp.f"
	    for (i__ = j - *n + *k + *l + 1; i__ <= i__2; ++i__) {
#line 524 "cggsvp.f"
		i__3 = i__ + j * a_dim1;
#line 524 "cggsvp.f"
		a[i__3].r = 0., a[i__3].i = 0.;
#line 525 "cggsvp.f"
/* L130: */
#line 525 "cggsvp.f"
	    }
#line 526 "cggsvp.f"
/* L140: */
#line 526 "cggsvp.f"
	}

#line 528 "cggsvp.f"
    }

#line 530 "cggsvp.f"
    return 0;

/*     End of CGGSVP */

} /* cggsvp_ */

