#line 1 "zggsvp.f"
/* zggsvp.f -- translated by f2c (version 20100827).
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

#line 1 "zggsvp.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};

/* > \brief \b ZGGSVP */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZGGSVP + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zggsvp.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zggsvp.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zggsvp.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZGGSVP( JOBU, JOBV, JOBQ, M, P, N, A, LDA, B, LDB, */
/*                          TOLA, TOLB, K, L, U, LDU, V, LDV, Q, LDQ, */
/*                          IWORK, RWORK, TAU, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOBQ, JOBU, JOBV */
/*       INTEGER            INFO, K, L, LDA, LDB, LDQ, LDU, LDV, M, N, P */
/*       DOUBLE PRECISION   TOLA, TOLB */
/*       .. */
/*       .. Array Arguments .. */
/*       INTEGER            IWORK( * ) */
/*       DOUBLE PRECISION   RWORK( * ) */
/*       COMPLEX*16         A( LDA, * ), B( LDB, * ), Q( LDQ, * ), */
/*      $                   TAU( * ), U( LDU, * ), V( LDV, * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZGGSVP computes unitary matrices U, V and Q such that */
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
/* > ZGGSVD. */
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
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
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
/* >          B is COMPLEX*16 array, dimension (LDB,N) */
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
/* >             TOLA = MAX(M,N)*norm(A)*MAZHEPS, */
/* >             TOLB = MAX(P,N)*norm(B)*MAZHEPS. */
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
/* >          U is COMPLEX*16 array, dimension (LDU,M) */
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
/* >          V is COMPLEX*16 array, dimension (LDV,P) */
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
/* >          Q is COMPLEX*16 array, dimension (LDQ,N) */
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
/* >          RWORK is DOUBLE PRECISION array, dimension (2*N) */
/* > \endverbatim */
/* > */
/* > \param[out] TAU */
/* > \verbatim */
/* >          TAU is COMPLEX*16 array, dimension (N) */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX*16 array, dimension (max(3*N,M,P)) */
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

/* > \ingroup complex16OTHERcomputational */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  The subroutine uses LAPACK subroutine ZGEQPF for the QR factorization */
/* >  with column pivoting to detect the effective numerical rank of the */
/* >  a matrix. It may be replaced by a better rank determination strategy. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int zggsvp_(char *jobu, char *jobv, char *jobq, integer *m, 
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
    extern /* Subroutine */ int zgeqr2_(integer *, integer *, doublecomplex *,
	     integer *, doublecomplex *, doublecomplex *, integer *), zgerq2_(
	    integer *, integer *, doublecomplex *, integer *, doublecomplex *,
	     doublecomplex *, integer *), zung2r_(integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *), zunm2r_(char *, char *, integer *, 
	    integer *, integer *, doublecomplex *, integer *, doublecomplex *,
	     doublecomplex *, integer *, doublecomplex *, integer *, ftnlen, 
	    ftnlen), zunmr2_(char *, char *, integer *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, ftnlen, ftnlen), xerbla_(
	    char *, integer *, ftnlen), zgeqpf_(integer *, integer *, 
	    doublecomplex *, integer *, integer *, doublecomplex *, 
	    doublecomplex *, doublereal *, integer *), zlacpy_(char *, 
	    integer *, integer *, doublecomplex *, integer *, doublecomplex *,
	     integer *, ftnlen);
    static logical forwrd;
    extern /* Subroutine */ int zlaset_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, doublecomplex *, integer *, 
	    ftnlen), zlapmt_(logical *, integer *, integer *, doublecomplex *,
	     integer *, integer *);


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

#line 316 "zggsvp.f"
    /* Parameter adjustments */
#line 316 "zggsvp.f"
    a_dim1 = *lda;
#line 316 "zggsvp.f"
    a_offset = 1 + a_dim1;
#line 316 "zggsvp.f"
    a -= a_offset;
#line 316 "zggsvp.f"
    b_dim1 = *ldb;
#line 316 "zggsvp.f"
    b_offset = 1 + b_dim1;
#line 316 "zggsvp.f"
    b -= b_offset;
#line 316 "zggsvp.f"
    u_dim1 = *ldu;
#line 316 "zggsvp.f"
    u_offset = 1 + u_dim1;
#line 316 "zggsvp.f"
    u -= u_offset;
#line 316 "zggsvp.f"
    v_dim1 = *ldv;
#line 316 "zggsvp.f"
    v_offset = 1 + v_dim1;
#line 316 "zggsvp.f"
    v -= v_offset;
#line 316 "zggsvp.f"
    q_dim1 = *ldq;
#line 316 "zggsvp.f"
    q_offset = 1 + q_dim1;
#line 316 "zggsvp.f"
    q -= q_offset;
#line 316 "zggsvp.f"
    --iwork;
#line 316 "zggsvp.f"
    --rwork;
#line 316 "zggsvp.f"
    --tau;
#line 316 "zggsvp.f"
    --work;
#line 316 "zggsvp.f"

#line 316 "zggsvp.f"
    /* Function Body */
#line 316 "zggsvp.f"
    wantu = lsame_(jobu, "U", (ftnlen)1, (ftnlen)1);
#line 317 "zggsvp.f"
    wantv = lsame_(jobv, "V", (ftnlen)1, (ftnlen)1);
#line 318 "zggsvp.f"
    wantq = lsame_(jobq, "Q", (ftnlen)1, (ftnlen)1);
#line 319 "zggsvp.f"
    forwrd = TRUE_;

#line 321 "zggsvp.f"
    *info = 0;
#line 322 "zggsvp.f"
    if (! (wantu || lsame_(jobu, "N", (ftnlen)1, (ftnlen)1))) {
#line 323 "zggsvp.f"
	*info = -1;
#line 324 "zggsvp.f"
    } else if (! (wantv || lsame_(jobv, "N", (ftnlen)1, (ftnlen)1))) {
#line 325 "zggsvp.f"
	*info = -2;
#line 326 "zggsvp.f"
    } else if (! (wantq || lsame_(jobq, "N", (ftnlen)1, (ftnlen)1))) {
#line 327 "zggsvp.f"
	*info = -3;
#line 328 "zggsvp.f"
    } else if (*m < 0) {
#line 329 "zggsvp.f"
	*info = -4;
#line 330 "zggsvp.f"
    } else if (*p < 0) {
#line 331 "zggsvp.f"
	*info = -5;
#line 332 "zggsvp.f"
    } else if (*n < 0) {
#line 333 "zggsvp.f"
	*info = -6;
#line 334 "zggsvp.f"
    } else if (*lda < max(1,*m)) {
#line 335 "zggsvp.f"
	*info = -8;
#line 336 "zggsvp.f"
    } else if (*ldb < max(1,*p)) {
#line 337 "zggsvp.f"
	*info = -10;
#line 338 "zggsvp.f"
    } else if (*ldu < 1 || wantu && *ldu < *m) {
#line 339 "zggsvp.f"
	*info = -16;
#line 340 "zggsvp.f"
    } else if (*ldv < 1 || wantv && *ldv < *p) {
#line 341 "zggsvp.f"
	*info = -18;
#line 342 "zggsvp.f"
    } else if (*ldq < 1 || wantq && *ldq < *n) {
#line 343 "zggsvp.f"
	*info = -20;
#line 344 "zggsvp.f"
    }
#line 345 "zggsvp.f"
    if (*info != 0) {
#line 346 "zggsvp.f"
	i__1 = -(*info);
#line 346 "zggsvp.f"
	xerbla_("ZGGSVP", &i__1, (ftnlen)6);
#line 347 "zggsvp.f"
	return 0;
#line 348 "zggsvp.f"
    }

/*     QR with column pivoting of B: B*P = V*( S11 S12 ) */
/*                                           (  0   0  ) */

#line 353 "zggsvp.f"
    i__1 = *n;
#line 353 "zggsvp.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 354 "zggsvp.f"
	iwork[i__] = 0;
#line 355 "zggsvp.f"
/* L10: */
#line 355 "zggsvp.f"
    }
#line 356 "zggsvp.f"
    zgeqpf_(p, n, &b[b_offset], ldb, &iwork[1], &tau[1], &work[1], &rwork[1], 
	    info);

/*     Update A := A*P */

#line 360 "zggsvp.f"
    zlapmt_(&forwrd, m, n, &a[a_offset], lda, &iwork[1]);

/*     Determine the effective rank of matrix B. */

#line 364 "zggsvp.f"
    *l = 0;
#line 365 "zggsvp.f"
    i__1 = min(*p,*n);
#line 365 "zggsvp.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 366 "zggsvp.f"
	i__2 = i__ + i__ * b_dim1;
#line 366 "zggsvp.f"
	if ((d__1 = b[i__2].r, abs(d__1)) + (d__2 = d_imag(&b[i__ + i__ * 
		b_dim1]), abs(d__2)) > *tolb) {
#line 366 "zggsvp.f"
	    ++(*l);
#line 366 "zggsvp.f"
	}
#line 368 "zggsvp.f"
/* L20: */
#line 368 "zggsvp.f"
    }

#line 370 "zggsvp.f"
    if (wantv) {

/*        Copy the details of V, and form V. */

#line 374 "zggsvp.f"
	zlaset_("Full", p, p, &c_b1, &c_b1, &v[v_offset], ldv, (ftnlen)4);
#line 375 "zggsvp.f"
	if (*p > 1) {
#line 375 "zggsvp.f"
	    i__1 = *p - 1;
#line 375 "zggsvp.f"
	    zlacpy_("Lower", &i__1, n, &b[b_dim1 + 2], ldb, &v[v_dim1 + 2], 
		    ldv, (ftnlen)5);
#line 375 "zggsvp.f"
	}
#line 378 "zggsvp.f"
	i__1 = min(*p,*n);
#line 378 "zggsvp.f"
	zung2r_(p, p, &i__1, &v[v_offset], ldv, &tau[1], &work[1], info);
#line 379 "zggsvp.f"
    }

/*     Clean up B */

#line 383 "zggsvp.f"
    i__1 = *l - 1;
#line 383 "zggsvp.f"
    for (j = 1; j <= i__1; ++j) {
#line 384 "zggsvp.f"
	i__2 = *l;
#line 384 "zggsvp.f"
	for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 385 "zggsvp.f"
	    i__3 = i__ + j * b_dim1;
#line 385 "zggsvp.f"
	    b[i__3].r = 0., b[i__3].i = 0.;
#line 386 "zggsvp.f"
/* L30: */
#line 386 "zggsvp.f"
	}
#line 387 "zggsvp.f"
/* L40: */
#line 387 "zggsvp.f"
    }
#line 388 "zggsvp.f"
    if (*p > *l) {
#line 388 "zggsvp.f"
	i__1 = *p - *l;
#line 388 "zggsvp.f"
	zlaset_("Full", &i__1, n, &c_b1, &c_b1, &b[*l + 1 + b_dim1], ldb, (
		ftnlen)4);
#line 388 "zggsvp.f"
    }

#line 391 "zggsvp.f"
    if (wantq) {

/*        Set Q = I and Update Q := Q*P */

#line 395 "zggsvp.f"
	zlaset_("Full", n, n, &c_b1, &c_b2, &q[q_offset], ldq, (ftnlen)4);
#line 396 "zggsvp.f"
	zlapmt_(&forwrd, n, n, &q[q_offset], ldq, &iwork[1]);
#line 397 "zggsvp.f"
    }

#line 399 "zggsvp.f"
    if (*p >= *l && *n != *l) {

/*        RQ factorization of ( S11 S12 ) = ( 0 S12 )*Z */

#line 403 "zggsvp.f"
	zgerq2_(l, n, &b[b_offset], ldb, &tau[1], &work[1], info);

/*        Update A := A*Z**H */

#line 407 "zggsvp.f"
	zunmr2_("Right", "Conjugate transpose", m, n, l, &b[b_offset], ldb, &
		tau[1], &a[a_offset], lda, &work[1], info, (ftnlen)5, (ftnlen)
		19);
#line 409 "zggsvp.f"
	if (wantq) {

/*           Update Q := Q*Z**H */

#line 413 "zggsvp.f"
	    zunmr2_("Right", "Conjugate transpose", n, n, l, &b[b_offset], 
		    ldb, &tau[1], &q[q_offset], ldq, &work[1], info, (ftnlen)
		    5, (ftnlen)19);
#line 415 "zggsvp.f"
	}

/*        Clean up B */

#line 419 "zggsvp.f"
	i__1 = *n - *l;
#line 419 "zggsvp.f"
	zlaset_("Full", l, &i__1, &c_b1, &c_b1, &b[b_offset], ldb, (ftnlen)4);
#line 420 "zggsvp.f"
	i__1 = *n;
#line 420 "zggsvp.f"
	for (j = *n - *l + 1; j <= i__1; ++j) {
#line 421 "zggsvp.f"
	    i__2 = *l;
#line 421 "zggsvp.f"
	    for (i__ = j - *n + *l + 1; i__ <= i__2; ++i__) {
#line 422 "zggsvp.f"
		i__3 = i__ + j * b_dim1;
#line 422 "zggsvp.f"
		b[i__3].r = 0., b[i__3].i = 0.;
#line 423 "zggsvp.f"
/* L50: */
#line 423 "zggsvp.f"
	    }
#line 424 "zggsvp.f"
/* L60: */
#line 424 "zggsvp.f"
	}

#line 426 "zggsvp.f"
    }

/*     Let              N-L     L */
/*                A = ( A11    A12 ) M, */

/*     then the following does the complete QR decomposition of A11: */

/*              A11 = U*(  0  T12 )*P1**H */
/*                      (  0   0  ) */

#line 436 "zggsvp.f"
    i__1 = *n - *l;
#line 436 "zggsvp.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 437 "zggsvp.f"
	iwork[i__] = 0;
#line 438 "zggsvp.f"
/* L70: */
#line 438 "zggsvp.f"
    }
#line 439 "zggsvp.f"
    i__1 = *n - *l;
#line 439 "zggsvp.f"
    zgeqpf_(m, &i__1, &a[a_offset], lda, &iwork[1], &tau[1], &work[1], &rwork[
	    1], info);

/*     Determine the effective rank of A11 */

#line 443 "zggsvp.f"
    *k = 0;
/* Computing MIN */
#line 444 "zggsvp.f"
    i__2 = *m, i__3 = *n - *l;
#line 444 "zggsvp.f"
    i__1 = min(i__2,i__3);
#line 444 "zggsvp.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 445 "zggsvp.f"
	i__2 = i__ + i__ * a_dim1;
#line 445 "zggsvp.f"
	if ((d__1 = a[i__2].r, abs(d__1)) + (d__2 = d_imag(&a[i__ + i__ * 
		a_dim1]), abs(d__2)) > *tola) {
#line 445 "zggsvp.f"
	    ++(*k);
#line 445 "zggsvp.f"
	}
#line 447 "zggsvp.f"
/* L80: */
#line 447 "zggsvp.f"
    }

/*     Update A12 := U**H*A12, where A12 = A( 1:M, N-L+1:N ) */

/* Computing MIN */
#line 451 "zggsvp.f"
    i__2 = *m, i__3 = *n - *l;
#line 451 "zggsvp.f"
    i__1 = min(i__2,i__3);
#line 451 "zggsvp.f"
    zunm2r_("Left", "Conjugate transpose", m, l, &i__1, &a[a_offset], lda, &
	    tau[1], &a[(*n - *l + 1) * a_dim1 + 1], lda, &work[1], info, (
	    ftnlen)4, (ftnlen)19);

#line 454 "zggsvp.f"
    if (wantu) {

/*        Copy the details of U, and form U */

#line 458 "zggsvp.f"
	zlaset_("Full", m, m, &c_b1, &c_b1, &u[u_offset], ldu, (ftnlen)4);
#line 459 "zggsvp.f"
	if (*m > 1) {
#line 459 "zggsvp.f"
	    i__1 = *m - 1;
#line 459 "zggsvp.f"
	    i__2 = *n - *l;
#line 459 "zggsvp.f"
	    zlacpy_("Lower", &i__1, &i__2, &a[a_dim1 + 2], lda, &u[u_dim1 + 2]
		    , ldu, (ftnlen)5);
#line 459 "zggsvp.f"
	}
/* Computing MIN */
#line 462 "zggsvp.f"
	i__2 = *m, i__3 = *n - *l;
#line 462 "zggsvp.f"
	i__1 = min(i__2,i__3);
#line 462 "zggsvp.f"
	zung2r_(m, m, &i__1, &u[u_offset], ldu, &tau[1], &work[1], info);
#line 463 "zggsvp.f"
    }

#line 465 "zggsvp.f"
    if (wantq) {

/*        Update Q( 1:N, 1:N-L )  = Q( 1:N, 1:N-L )*P1 */

#line 469 "zggsvp.f"
	i__1 = *n - *l;
#line 469 "zggsvp.f"
	zlapmt_(&forwrd, n, &i__1, &q[q_offset], ldq, &iwork[1]);
#line 470 "zggsvp.f"
    }

/*     Clean up A: set the strictly lower triangular part of */
/*     A(1:K, 1:K) = 0, and A( K+1:M, 1:N-L ) = 0. */

#line 475 "zggsvp.f"
    i__1 = *k - 1;
#line 475 "zggsvp.f"
    for (j = 1; j <= i__1; ++j) {
#line 476 "zggsvp.f"
	i__2 = *k;
#line 476 "zggsvp.f"
	for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 477 "zggsvp.f"
	    i__3 = i__ + j * a_dim1;
#line 477 "zggsvp.f"
	    a[i__3].r = 0., a[i__3].i = 0.;
#line 478 "zggsvp.f"
/* L90: */
#line 478 "zggsvp.f"
	}
#line 479 "zggsvp.f"
/* L100: */
#line 479 "zggsvp.f"
    }
#line 480 "zggsvp.f"
    if (*m > *k) {
#line 480 "zggsvp.f"
	i__1 = *m - *k;
#line 480 "zggsvp.f"
	i__2 = *n - *l;
#line 480 "zggsvp.f"
	zlaset_("Full", &i__1, &i__2, &c_b1, &c_b1, &a[*k + 1 + a_dim1], lda, 
		(ftnlen)4);
#line 480 "zggsvp.f"
    }

#line 483 "zggsvp.f"
    if (*n - *l > *k) {

/*        RQ factorization of ( T11 T12 ) = ( 0 T12 )*Z1 */

#line 487 "zggsvp.f"
	i__1 = *n - *l;
#line 487 "zggsvp.f"
	zgerq2_(k, &i__1, &a[a_offset], lda, &tau[1], &work[1], info);

#line 489 "zggsvp.f"
	if (wantq) {

/*           Update Q( 1:N,1:N-L ) = Q( 1:N,1:N-L )*Z1**H */

#line 493 "zggsvp.f"
	    i__1 = *n - *l;
#line 493 "zggsvp.f"
	    zunmr2_("Right", "Conjugate transpose", n, &i__1, k, &a[a_offset],
		     lda, &tau[1], &q[q_offset], ldq, &work[1], info, (ftnlen)
		    5, (ftnlen)19);
#line 495 "zggsvp.f"
	}

/*        Clean up A */

#line 499 "zggsvp.f"
	i__1 = *n - *l - *k;
#line 499 "zggsvp.f"
	zlaset_("Full", k, &i__1, &c_b1, &c_b1, &a[a_offset], lda, (ftnlen)4);
#line 500 "zggsvp.f"
	i__1 = *n - *l;
#line 500 "zggsvp.f"
	for (j = *n - *l - *k + 1; j <= i__1; ++j) {
#line 501 "zggsvp.f"
	    i__2 = *k;
#line 501 "zggsvp.f"
	    for (i__ = j - *n + *l + *k + 1; i__ <= i__2; ++i__) {
#line 502 "zggsvp.f"
		i__3 = i__ + j * a_dim1;
#line 502 "zggsvp.f"
		a[i__3].r = 0., a[i__3].i = 0.;
#line 503 "zggsvp.f"
/* L110: */
#line 503 "zggsvp.f"
	    }
#line 504 "zggsvp.f"
/* L120: */
#line 504 "zggsvp.f"
	}

#line 506 "zggsvp.f"
    }

#line 508 "zggsvp.f"
    if (*m > *k) {

/*        QR factorization of A( K+1:M,N-L+1:N ) */

#line 512 "zggsvp.f"
	i__1 = *m - *k;
#line 512 "zggsvp.f"
	zgeqr2_(&i__1, l, &a[*k + 1 + (*n - *l + 1) * a_dim1], lda, &tau[1], &
		work[1], info);

#line 514 "zggsvp.f"
	if (wantu) {

/*           Update U(:,K+1:M) := U(:,K+1:M)*U1 */

#line 518 "zggsvp.f"
	    i__1 = *m - *k;
/* Computing MIN */
#line 518 "zggsvp.f"
	    i__3 = *m - *k;
#line 518 "zggsvp.f"
	    i__2 = min(i__3,*l);
#line 518 "zggsvp.f"
	    zunm2r_("Right", "No transpose", m, &i__1, &i__2, &a[*k + 1 + (*n 
		    - *l + 1) * a_dim1], lda, &tau[1], &u[(*k + 1) * u_dim1 + 
		    1], ldu, &work[1], info, (ftnlen)5, (ftnlen)12);
#line 521 "zggsvp.f"
	}

/*        Clean up */

#line 525 "zggsvp.f"
	i__1 = *n;
#line 525 "zggsvp.f"
	for (j = *n - *l + 1; j <= i__1; ++j) {
#line 526 "zggsvp.f"
	    i__2 = *m;
#line 526 "zggsvp.f"
	    for (i__ = j - *n + *k + *l + 1; i__ <= i__2; ++i__) {
#line 527 "zggsvp.f"
		i__3 = i__ + j * a_dim1;
#line 527 "zggsvp.f"
		a[i__3].r = 0., a[i__3].i = 0.;
#line 528 "zggsvp.f"
/* L130: */
#line 528 "zggsvp.f"
	    }
#line 529 "zggsvp.f"
/* L140: */
#line 529 "zggsvp.f"
	}

#line 531 "zggsvp.f"
    }

#line 533 "zggsvp.f"
    return 0;

/*     End of ZGGSVP */

} /* zggsvp_ */

