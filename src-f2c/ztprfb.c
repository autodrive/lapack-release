#line 1 "ztprfb.f"
/* ztprfb.f -- translated by f2c (version 20100827).
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

#line 1 "ztprfb.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static doublecomplex c_b2 = {0.,0.};

/* > \brief \b ZTPRFB applies a real or complex "triangular-pentagonal" blocked reflector to a real or complex
 matrix, which is composed of two blocks. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZTPRFB + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztprfb.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztprfb.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztprfb.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZTPRFB( SIDE, TRANS, DIRECT, STOREV, M, N, K, L, */
/*                          V, LDV, T, LDT, A, LDA, B, LDB, WORK, LDWORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER DIRECT, SIDE, STOREV, TRANS */
/*       INTEGER   K, L, LDA, LDB, LDT, LDV, LDWORK, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16   A( LDA, * ), B( LDB, * ), T( LDT, * ), */
/*      $          V( LDV, * ), WORK( LDWORK, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZTPRFB applies a complex "triangular-pentagonal" block reflector H or its */
/* > conjugate transpose H**H to a complex matrix C, which is composed of two */
/* > blocks A and B, either from the left or right. */
/* > */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] SIDE */
/* > \verbatim */
/* >          SIDE is CHARACTER*1 */
/* >          = 'L': apply H or H**H from the Left */
/* >          = 'R': apply H or H**H from the Right */
/* > \endverbatim */
/* > */
/* > \param[in] TRANS */
/* > \verbatim */
/* >          TRANS is CHARACTER*1 */
/* >          = 'N': apply H (No transpose) */
/* >          = 'C': apply H**H (Conjugate transpose) */
/* > \endverbatim */
/* > */
/* > \param[in] DIRECT */
/* > \verbatim */
/* >          DIRECT is CHARACTER*1 */
/* >          Indicates how H is formed from a product of elementary */
/* >          reflectors */
/* >          = 'F': H = H(1) H(2) . . . H(k) (Forward) */
/* >          = 'B': H = H(k) . . . H(2) H(1) (Backward) */
/* > \endverbatim */
/* > */
/* > \param[in] STOREV */
/* > \verbatim */
/* >          STOREV is CHARACTER*1 */
/* >          Indicates how the vectors which define the elementary */
/* >          reflectors are stored: */
/* >          = 'C': Columns */
/* >          = 'R': Rows */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the matrix B. */
/* >          M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix B. */
/* >          N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] K */
/* > \verbatim */
/* >          K is INTEGER */
/* >          The order of the matrix T, i.e. the number of elementary */
/* >          reflectors whose product defines the block reflector. */
/* >          K >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] L */
/* > \verbatim */
/* >          L is INTEGER */
/* >          The order of the trapezoidal part of V. */
/* >          K >= L >= 0.  See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] V */
/* > \verbatim */
/* >          V is COMPLEX*16 array, dimension */
/* >                                (LDV,K) if STOREV = 'C' */
/* >                                (LDV,M) if STOREV = 'R' and SIDE = 'L' */
/* >                                (LDV,N) if STOREV = 'R' and SIDE = 'R' */
/* >          The pentagonal matrix V, which contains the elementary reflectors */
/* >          H(1), H(2), ..., H(K).  See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDV */
/* > \verbatim */
/* >          LDV is INTEGER */
/* >          The leading dimension of the array V. */
/* >          If STOREV = 'C' and SIDE = 'L', LDV >= max(1,M); */
/* >          if STOREV = 'C' and SIDE = 'R', LDV >= max(1,N); */
/* >          if STOREV = 'R', LDV >= K. */
/* > \endverbatim */
/* > */
/* > \param[in] T */
/* > \verbatim */
/* >          T is COMPLEX*16 array, dimension (LDT,K) */
/* >          The triangular K-by-K matrix T in the representation of the */
/* >          block reflector. */
/* > \endverbatim */
/* > */
/* > \param[in] LDT */
/* > \verbatim */
/* >          LDT is INTEGER */
/* >          The leading dimension of the array T. */
/* >          LDT >= K. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension */
/* >          (LDA,N) if SIDE = 'L' or (LDA,K) if SIDE = 'R' */
/* >          On entry, the K-by-N or M-by-K matrix A. */
/* >          On exit, A is overwritten by the corresponding block of */
/* >          H*C or H**H*C or C*H or C*H**H.  See Futher Details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The leading dimension of the array A. */
/* >          If SIDE = 'L', LDC >= max(1,K); */
/* >          If SIDE = 'R', LDC >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* >          B is COMPLEX*16 array, dimension (LDB,N) */
/* >          On entry, the M-by-N matrix B. */
/* >          On exit, B is overwritten by the corresponding block of */
/* >          H*C or H**H*C or C*H or C*H**H.  See Further Details. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >          The leading dimension of the array B. */
/* >          LDB >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX*16 array, dimension */
/* >          (LDWORK,N) if SIDE = 'L', */
/* >          (LDWORK,K) if SIDE = 'R'. */
/* > \endverbatim */
/* > */
/* > \param[in] LDWORK */
/* > \verbatim */
/* >          LDWORK is INTEGER */
/* >          The leading dimension of the array WORK. */
/* >          If SIDE = 'L', LDWORK >= K; */
/* >          if SIDE = 'R', LDWORK >= M. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup complex16OTHERauxiliary */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  The matrix C is a composite matrix formed from blocks A and B. */
/* >  The block B is of size M-by-N; if SIDE = 'R', A is of size M-by-K, */
/* >  and if SIDE = 'L', A is of size K-by-N. */
/* > */
/* >  If SIDE = 'R' and DIRECT = 'F', C = [A B]. */
/* > */
/* >  If SIDE = 'L' and DIRECT = 'F', C = [A] */
/* >                                      [B]. */
/* > */
/* >  If SIDE = 'R' and DIRECT = 'B', C = [B A]. */
/* > */
/* >  If SIDE = 'L' and DIRECT = 'B', C = [B] */
/* >                                      [A]. */
/* > */
/* >  The pentagonal matrix V is composed of a rectangular block V1 and a */
/* >  trapezoidal block V2.  The size of the trapezoidal block is determined by */
/* >  the parameter L, where 0<=L<=K.  If L=K, the V2 block of V is triangular; */
/* >  if L=0, there is no trapezoidal block, thus V = V1 is rectangular. */
/* > */
/* >  If DIRECT = 'F' and STOREV = 'C':  V = [V1] */
/* >                                         [V2] */
/* >     - V2 is upper trapezoidal (first L rows of K-by-K upper triangular) */
/* > */
/* >  If DIRECT = 'F' and STOREV = 'R':  V = [V1 V2] */
/* > */
/* >     - V2 is lower trapezoidal (first L columns of K-by-K lower triangular) */
/* > */
/* >  If DIRECT = 'B' and STOREV = 'C':  V = [V2] */
/* >                                         [V1] */
/* >     - V2 is lower trapezoidal (last L rows of K-by-K lower triangular) */
/* > */
/* >  If DIRECT = 'B' and STOREV = 'R':  V = [V2 V1] */
/* > */
/* >     - V2 is upper trapezoidal (last L columns of K-by-K upper triangular) */
/* > */
/* >  If STOREV = 'C' and SIDE = 'L', V is M-by-K with V2 L-by-K. */
/* > */
/* >  If STOREV = 'C' and SIDE = 'R', V is N-by-K with V2 L-by-K. */
/* > */
/* >  If STOREV = 'R' and SIDE = 'L', V is K-by-M with V2 K-by-L. */
/* > */
/* >  If STOREV = 'R' and SIDE = 'R', V is K-by-N with V2 K-by-L. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int ztprfb_(char *side, char *trans, char *direct, char *
	storev, integer *m, integer *n, integer *k, integer *l, doublecomplex 
	*v, integer *ldv, doublecomplex *t, integer *ldt, doublecomplex *a, 
	integer *lda, doublecomplex *b, integer *ldb, doublecomplex *work, 
	integer *ldwork, ftnlen side_len, ftnlen trans_len, ftnlen direct_len,
	 ftnlen storev_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, t_dim1, t_offset, v_dim1, 
	    v_offset, work_dim1, work_offset, i__1, i__2, i__3, i__4, i__5;
    doublecomplex z__1;

    /* Local variables */
    static logical backward;
    static integer i__, j, kp, mp, np;
    static logical row, left;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical right;
    extern /* Subroutine */ int zgemm_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, ftnlen, ftnlen), ztrmm_(char *, char *, char *, char *,
	     integer *, integer *, doublecomplex *, doublecomplex *, integer *
	    , doublecomplex *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static logical column, forward;


/*  -- LAPACK auxiliary routine (version 3.4.2) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     September 2012 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ========================================================================== */

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

/*     Quick return if possible */

#line 292 "ztprfb.f"
    /* Parameter adjustments */
#line 292 "ztprfb.f"
    v_dim1 = *ldv;
#line 292 "ztprfb.f"
    v_offset = 1 + v_dim1;
#line 292 "ztprfb.f"
    v -= v_offset;
#line 292 "ztprfb.f"
    t_dim1 = *ldt;
#line 292 "ztprfb.f"
    t_offset = 1 + t_dim1;
#line 292 "ztprfb.f"
    t -= t_offset;
#line 292 "ztprfb.f"
    a_dim1 = *lda;
#line 292 "ztprfb.f"
    a_offset = 1 + a_dim1;
#line 292 "ztprfb.f"
    a -= a_offset;
#line 292 "ztprfb.f"
    b_dim1 = *ldb;
#line 292 "ztprfb.f"
    b_offset = 1 + b_dim1;
#line 292 "ztprfb.f"
    b -= b_offset;
#line 292 "ztprfb.f"
    work_dim1 = *ldwork;
#line 292 "ztprfb.f"
    work_offset = 1 + work_dim1;
#line 292 "ztprfb.f"
    work -= work_offset;
#line 292 "ztprfb.f"

#line 292 "ztprfb.f"
    /* Function Body */
#line 292 "ztprfb.f"
    if (*m <= 0 || *n <= 0 || *k <= 0 || *l < 0) {
#line 292 "ztprfb.f"
	return 0;
#line 292 "ztprfb.f"
    }

#line 294 "ztprfb.f"
    if (lsame_(storev, "C", (ftnlen)1, (ftnlen)1)) {
#line 295 "ztprfb.f"
	column = TRUE_;
#line 296 "ztprfb.f"
	row = FALSE_;
#line 297 "ztprfb.f"
    } else if (lsame_(storev, "R", (ftnlen)1, (ftnlen)1)) {
#line 298 "ztprfb.f"
	column = FALSE_;
#line 299 "ztprfb.f"
	row = TRUE_;
#line 300 "ztprfb.f"
    } else {
#line 301 "ztprfb.f"
	column = FALSE_;
#line 302 "ztprfb.f"
	row = FALSE_;
#line 303 "ztprfb.f"
    }

#line 305 "ztprfb.f"
    if (lsame_(side, "L", (ftnlen)1, (ftnlen)1)) {
#line 306 "ztprfb.f"
	left = TRUE_;
#line 307 "ztprfb.f"
	right = FALSE_;
#line 308 "ztprfb.f"
    } else if (lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {
#line 309 "ztprfb.f"
	left = FALSE_;
#line 310 "ztprfb.f"
	right = TRUE_;
#line 311 "ztprfb.f"
    } else {
#line 312 "ztprfb.f"
	left = FALSE_;
#line 313 "ztprfb.f"
	right = FALSE_;
#line 314 "ztprfb.f"
    }

#line 316 "ztprfb.f"
    if (lsame_(direct, "F", (ftnlen)1, (ftnlen)1)) {
#line 317 "ztprfb.f"
	forward = TRUE_;
#line 318 "ztprfb.f"
	backward = FALSE_;
#line 319 "ztprfb.f"
    } else if (lsame_(direct, "B", (ftnlen)1, (ftnlen)1)) {
#line 320 "ztprfb.f"
	forward = FALSE_;
#line 321 "ztprfb.f"
	backward = TRUE_;
#line 322 "ztprfb.f"
    } else {
#line 323 "ztprfb.f"
	forward = FALSE_;
#line 324 "ztprfb.f"
	backward = FALSE_;
#line 325 "ztprfb.f"
    }

/* --------------------------------------------------------------------------- */

#line 329 "ztprfb.f"
    if (column && forward && left) {

/* --------------------------------------------------------------------------- */

/*        Let  W =  [ I ]    (K-by-K) */
/*                  [ V ]    (M-by-K) */

/*        Form  H C  or  H**H C  where  C = [ A ]  (K-by-N) */
/*                                          [ B ]  (M-by-N) */

/*        H = I - W T W**H          or  H**H = I - W T**H W**H */

/*        A = A -   T (A + V**H B)  or  A = A -   T**H (A + V**H B) */
/*        B = B - V T (A + V**H B)  or  B = B - V T**H (A + V**H B) */

/* --------------------------------------------------------------------------- */

/* Computing MIN */
#line 346 "ztprfb.f"
	i__1 = *m - *l + 1;
#line 346 "ztprfb.f"
	mp = min(i__1,*m);
/* Computing MIN */
#line 347 "ztprfb.f"
	i__1 = *l + 1;
#line 347 "ztprfb.f"
	kp = min(i__1,*k);

#line 349 "ztprfb.f"
	i__1 = *n;
#line 349 "ztprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 350 "ztprfb.f"
	    i__2 = *l;
#line 350 "ztprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 351 "ztprfb.f"
		i__3 = i__ + j * work_dim1;
#line 351 "ztprfb.f"
		i__4 = *m - *l + i__ + j * b_dim1;
#line 351 "ztprfb.f"
		work[i__3].r = b[i__4].r, work[i__3].i = b[i__4].i;
#line 352 "ztprfb.f"
	    }
#line 353 "ztprfb.f"
	}
#line 354 "ztprfb.f"
	ztrmm_("L", "U", "C", "N", l, n, &c_b1, &v[mp + v_dim1], ldv, &work[
		work_offset], ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		ftnlen)1);
#line 356 "ztprfb.f"
	i__1 = *m - *l;
#line 356 "ztprfb.f"
	zgemm_("C", "N", l, n, &i__1, &c_b1, &v[v_offset], ldv, &b[b_offset], 
		ldb, &c_b1, &work[work_offset], ldwork, (ftnlen)1, (ftnlen)1);
#line 358 "ztprfb.f"
	i__1 = *k - *l;
#line 358 "ztprfb.f"
	zgemm_("C", "N", &i__1, n, m, &c_b1, &v[kp * v_dim1 + 1], ldv, &b[
		b_offset], ldb, &c_b2, &work[kp + work_dim1], ldwork, (ftnlen)
		1, (ftnlen)1);

#line 361 "ztprfb.f"
	i__1 = *n;
#line 361 "ztprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 362 "ztprfb.f"
	    i__2 = *k;
#line 362 "ztprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 363 "ztprfb.f"
		i__3 = i__ + j * work_dim1;
#line 363 "ztprfb.f"
		i__4 = i__ + j * work_dim1;
#line 363 "ztprfb.f"
		i__5 = i__ + j * a_dim1;
#line 363 "ztprfb.f"
		z__1.r = work[i__4].r + a[i__5].r, z__1.i = work[i__4].i + a[
			i__5].i;
#line 363 "ztprfb.f"
		work[i__3].r = z__1.r, work[i__3].i = z__1.i;
#line 364 "ztprfb.f"
	    }
#line 365 "ztprfb.f"
	}

#line 367 "ztprfb.f"
	ztrmm_("L", "U", trans, "N", k, n, &c_b1, &t[t_offset], ldt, &work[
		work_offset], ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		ftnlen)1);

#line 370 "ztprfb.f"
	i__1 = *n;
#line 370 "ztprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 371 "ztprfb.f"
	    i__2 = *k;
#line 371 "ztprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 372 "ztprfb.f"
		i__3 = i__ + j * a_dim1;
#line 372 "ztprfb.f"
		i__4 = i__ + j * a_dim1;
#line 372 "ztprfb.f"
		i__5 = i__ + j * work_dim1;
#line 372 "ztprfb.f"
		z__1.r = a[i__4].r - work[i__5].r, z__1.i = a[i__4].i - work[
			i__5].i;
#line 372 "ztprfb.f"
		a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 373 "ztprfb.f"
	    }
#line 374 "ztprfb.f"
	}

#line 376 "ztprfb.f"
	i__1 = *m - *l;
#line 376 "ztprfb.f"
	z__1.r = -1., z__1.i = -0.;
#line 376 "ztprfb.f"
	zgemm_("N", "N", &i__1, n, k, &z__1, &v[v_offset], ldv, &work[
		work_offset], ldwork, &c_b1, &b[b_offset], ldb, (ftnlen)1, (
		ftnlen)1);
#line 378 "ztprfb.f"
	i__1 = *k - *l;
#line 378 "ztprfb.f"
	z__1.r = -1., z__1.i = -0.;
#line 378 "ztprfb.f"
	zgemm_("N", "N", l, n, &i__1, &z__1, &v[mp + kp * v_dim1], ldv, &work[
		kp + work_dim1], ldwork, &c_b1, &b[mp + b_dim1], ldb, (ftnlen)
		1, (ftnlen)1);
#line 380 "ztprfb.f"
	ztrmm_("L", "U", "N", "N", l, n, &c_b1, &v[mp + v_dim1], ldv, &work[
		work_offset], ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		ftnlen)1);
#line 382 "ztprfb.f"
	i__1 = *n;
#line 382 "ztprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 383 "ztprfb.f"
	    i__2 = *l;
#line 383 "ztprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 384 "ztprfb.f"
		i__3 = *m - *l + i__ + j * b_dim1;
#line 384 "ztprfb.f"
		i__4 = *m - *l + i__ + j * b_dim1;
#line 384 "ztprfb.f"
		i__5 = i__ + j * work_dim1;
#line 384 "ztprfb.f"
		z__1.r = b[i__4].r - work[i__5].r, z__1.i = b[i__4].i - work[
			i__5].i;
#line 384 "ztprfb.f"
		b[i__3].r = z__1.r, b[i__3].i = z__1.i;
#line 385 "ztprfb.f"
	    }
#line 386 "ztprfb.f"
	}

/* --------------------------------------------------------------------------- */

#line 390 "ztprfb.f"
    } else if (column && forward && right) {

/* --------------------------------------------------------------------------- */

/*        Let  W =  [ I ]    (K-by-K) */
/*                  [ V ]    (N-by-K) */

/*        Form  C H or  C H**H  where  C = [ A B ] (A is M-by-K, B is M-by-N) */

/*        H = I - W T W**H          or  H**H = I - W T**H W**H */

/*        A = A - (A + B V) T      or  A = A - (A + B V) T**H */
/*        B = B - (A + B V) T V**H  or  B = B - (A + B V) T**H V**H */

/* --------------------------------------------------------------------------- */

/* Computing MIN */
#line 406 "ztprfb.f"
	i__1 = *n - *l + 1;
#line 406 "ztprfb.f"
	np = min(i__1,*n);
/* Computing MIN */
#line 407 "ztprfb.f"
	i__1 = *l + 1;
#line 407 "ztprfb.f"
	kp = min(i__1,*k);

#line 409 "ztprfb.f"
	i__1 = *l;
#line 409 "ztprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 410 "ztprfb.f"
	    i__2 = *m;
#line 410 "ztprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 411 "ztprfb.f"
		i__3 = i__ + j * work_dim1;
#line 411 "ztprfb.f"
		i__4 = i__ + (*n - *l + j) * b_dim1;
#line 411 "ztprfb.f"
		work[i__3].r = b[i__4].r, work[i__3].i = b[i__4].i;
#line 412 "ztprfb.f"
	    }
#line 413 "ztprfb.f"
	}
#line 414 "ztprfb.f"
	ztrmm_("R", "U", "N", "N", m, l, &c_b1, &v[np + v_dim1], ldv, &work[
		work_offset], ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		ftnlen)1);
#line 416 "ztprfb.f"
	i__1 = *n - *l;
#line 416 "ztprfb.f"
	zgemm_("N", "N", m, l, &i__1, &c_b1, &b[b_offset], ldb, &v[v_offset], 
		ldv, &c_b1, &work[work_offset], ldwork, (ftnlen)1, (ftnlen)1);
#line 418 "ztprfb.f"
	i__1 = *k - *l;
#line 418 "ztprfb.f"
	zgemm_("N", "N", m, &i__1, n, &c_b1, &b[b_offset], ldb, &v[kp * 
		v_dim1 + 1], ldv, &c_b2, &work[kp * work_dim1 + 1], ldwork, (
		ftnlen)1, (ftnlen)1);

#line 421 "ztprfb.f"
	i__1 = *k;
#line 421 "ztprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 422 "ztprfb.f"
	    i__2 = *m;
#line 422 "ztprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 423 "ztprfb.f"
		i__3 = i__ + j * work_dim1;
#line 423 "ztprfb.f"
		i__4 = i__ + j * work_dim1;
#line 423 "ztprfb.f"
		i__5 = i__ + j * a_dim1;
#line 423 "ztprfb.f"
		z__1.r = work[i__4].r + a[i__5].r, z__1.i = work[i__4].i + a[
			i__5].i;
#line 423 "ztprfb.f"
		work[i__3].r = z__1.r, work[i__3].i = z__1.i;
#line 424 "ztprfb.f"
	    }
#line 425 "ztprfb.f"
	}

#line 427 "ztprfb.f"
	ztrmm_("R", "U", trans, "N", m, k, &c_b1, &t[t_offset], ldt, &work[
		work_offset], ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		ftnlen)1);

#line 430 "ztprfb.f"
	i__1 = *k;
#line 430 "ztprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 431 "ztprfb.f"
	    i__2 = *m;
#line 431 "ztprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 432 "ztprfb.f"
		i__3 = i__ + j * a_dim1;
#line 432 "ztprfb.f"
		i__4 = i__ + j * a_dim1;
#line 432 "ztprfb.f"
		i__5 = i__ + j * work_dim1;
#line 432 "ztprfb.f"
		z__1.r = a[i__4].r - work[i__5].r, z__1.i = a[i__4].i - work[
			i__5].i;
#line 432 "ztprfb.f"
		a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 433 "ztprfb.f"
	    }
#line 434 "ztprfb.f"
	}

#line 436 "ztprfb.f"
	i__1 = *n - *l;
#line 436 "ztprfb.f"
	z__1.r = -1., z__1.i = -0.;
#line 436 "ztprfb.f"
	zgemm_("N", "C", m, &i__1, k, &z__1, &work[work_offset], ldwork, &v[
		v_offset], ldv, &c_b1, &b[b_offset], ldb, (ftnlen)1, (ftnlen)
		1);
#line 438 "ztprfb.f"
	i__1 = *k - *l;
#line 438 "ztprfb.f"
	z__1.r = -1., z__1.i = -0.;
#line 438 "ztprfb.f"
	zgemm_("N", "C", m, l, &i__1, &z__1, &work[kp * work_dim1 + 1], 
		ldwork, &v[np + kp * v_dim1], ldv, &c_b1, &b[np * b_dim1 + 1],
		 ldb, (ftnlen)1, (ftnlen)1);
#line 440 "ztprfb.f"
	ztrmm_("R", "U", "C", "N", m, l, &c_b1, &v[np + v_dim1], ldv, &work[
		work_offset], ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		ftnlen)1);
#line 442 "ztprfb.f"
	i__1 = *l;
#line 442 "ztprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 443 "ztprfb.f"
	    i__2 = *m;
#line 443 "ztprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 444 "ztprfb.f"
		i__3 = i__ + (*n - *l + j) * b_dim1;
#line 444 "ztprfb.f"
		i__4 = i__ + (*n - *l + j) * b_dim1;
#line 444 "ztprfb.f"
		i__5 = i__ + j * work_dim1;
#line 444 "ztprfb.f"
		z__1.r = b[i__4].r - work[i__5].r, z__1.i = b[i__4].i - work[
			i__5].i;
#line 444 "ztprfb.f"
		b[i__3].r = z__1.r, b[i__3].i = z__1.i;
#line 445 "ztprfb.f"
	    }
#line 446 "ztprfb.f"
	}

/* --------------------------------------------------------------------------- */

#line 450 "ztprfb.f"
    } else if (column && backward && left) {

/* --------------------------------------------------------------------------- */

/*        Let  W =  [ V ]    (M-by-K) */
/*                  [ I ]    (K-by-K) */

/*        Form  H C  or  H**H C  where  C = [ B ]  (M-by-N) */
/*                                          [ A ]  (K-by-N) */

/*        H = I - W T W**H          or  H**H = I - W T**H W**H */

/*        A = A -   T (A + V**H B)  or  A = A -   T**H (A + V**H B) */
/*        B = B - V T (A + V**H B)  or  B = B - V T**H (A + V**H B) */

/* --------------------------------------------------------------------------- */

/* Computing MIN */
#line 467 "ztprfb.f"
	i__1 = *l + 1;
#line 467 "ztprfb.f"
	mp = min(i__1,*m);
/* Computing MIN */
#line 468 "ztprfb.f"
	i__1 = *k - *l + 1;
#line 468 "ztprfb.f"
	kp = min(i__1,*k);

#line 470 "ztprfb.f"
	i__1 = *n;
#line 470 "ztprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 471 "ztprfb.f"
	    i__2 = *l;
#line 471 "ztprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 472 "ztprfb.f"
		i__3 = *k - *l + i__ + j * work_dim1;
#line 472 "ztprfb.f"
		i__4 = i__ + j * b_dim1;
#line 472 "ztprfb.f"
		work[i__3].r = b[i__4].r, work[i__3].i = b[i__4].i;
#line 473 "ztprfb.f"
	    }
#line 474 "ztprfb.f"
	}

#line 476 "ztprfb.f"
	ztrmm_("L", "L", "C", "N", l, n, &c_b1, &v[kp * v_dim1 + 1], ldv, &
		work[kp + work_dim1], ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1,
		 (ftnlen)1);
#line 478 "ztprfb.f"
	i__1 = *m - *l;
#line 478 "ztprfb.f"
	zgemm_("C", "N", l, n, &i__1, &c_b1, &v[mp + kp * v_dim1], ldv, &b[mp 
		+ b_dim1], ldb, &c_b1, &work[kp + work_dim1], ldwork, (ftnlen)
		1, (ftnlen)1);
#line 480 "ztprfb.f"
	i__1 = *k - *l;
#line 480 "ztprfb.f"
	zgemm_("C", "N", &i__1, n, m, &c_b1, &v[v_offset], ldv, &b[b_offset], 
		ldb, &c_b2, &work[work_offset], ldwork, (ftnlen)1, (ftnlen)1);

#line 483 "ztprfb.f"
	i__1 = *n;
#line 483 "ztprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 484 "ztprfb.f"
	    i__2 = *k;
#line 484 "ztprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 485 "ztprfb.f"
		i__3 = i__ + j * work_dim1;
#line 485 "ztprfb.f"
		i__4 = i__ + j * work_dim1;
#line 485 "ztprfb.f"
		i__5 = i__ + j * a_dim1;
#line 485 "ztprfb.f"
		z__1.r = work[i__4].r + a[i__5].r, z__1.i = work[i__4].i + a[
			i__5].i;
#line 485 "ztprfb.f"
		work[i__3].r = z__1.r, work[i__3].i = z__1.i;
#line 486 "ztprfb.f"
	    }
#line 487 "ztprfb.f"
	}

#line 489 "ztprfb.f"
	ztrmm_("L", "L", trans, "N", k, n, &c_b1, &t[t_offset], ldt, &work[
		work_offset], ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		ftnlen)1);

#line 492 "ztprfb.f"
	i__1 = *n;
#line 492 "ztprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 493 "ztprfb.f"
	    i__2 = *k;
#line 493 "ztprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 494 "ztprfb.f"
		i__3 = i__ + j * a_dim1;
#line 494 "ztprfb.f"
		i__4 = i__ + j * a_dim1;
#line 494 "ztprfb.f"
		i__5 = i__ + j * work_dim1;
#line 494 "ztprfb.f"
		z__1.r = a[i__4].r - work[i__5].r, z__1.i = a[i__4].i - work[
			i__5].i;
#line 494 "ztprfb.f"
		a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 495 "ztprfb.f"
	    }
#line 496 "ztprfb.f"
	}

#line 498 "ztprfb.f"
	i__1 = *m - *l;
#line 498 "ztprfb.f"
	z__1.r = -1., z__1.i = -0.;
#line 498 "ztprfb.f"
	zgemm_("N", "N", &i__1, n, k, &z__1, &v[mp + v_dim1], ldv, &work[
		work_offset], ldwork, &c_b1, &b[mp + b_dim1], ldb, (ftnlen)1, 
		(ftnlen)1);
#line 500 "ztprfb.f"
	i__1 = *k - *l;
#line 500 "ztprfb.f"
	z__1.r = -1., z__1.i = -0.;
#line 500 "ztprfb.f"
	zgemm_("N", "N", l, n, &i__1, &z__1, &v[v_offset], ldv, &work[
		work_offset], ldwork, &c_b1, &b[b_offset], ldb, (ftnlen)1, (
		ftnlen)1);
#line 502 "ztprfb.f"
	ztrmm_("L", "L", "N", "N", l, n, &c_b1, &v[kp * v_dim1 + 1], ldv, &
		work[kp + work_dim1], ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1,
		 (ftnlen)1);
#line 504 "ztprfb.f"
	i__1 = *n;
#line 504 "ztprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 505 "ztprfb.f"
	    i__2 = *l;
#line 505 "ztprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 506 "ztprfb.f"
		i__3 = i__ + j * b_dim1;
#line 506 "ztprfb.f"
		i__4 = i__ + j * b_dim1;
#line 506 "ztprfb.f"
		i__5 = *k - *l + i__ + j * work_dim1;
#line 506 "ztprfb.f"
		z__1.r = b[i__4].r - work[i__5].r, z__1.i = b[i__4].i - work[
			i__5].i;
#line 506 "ztprfb.f"
		b[i__3].r = z__1.r, b[i__3].i = z__1.i;
#line 507 "ztprfb.f"
	    }
#line 508 "ztprfb.f"
	}

/* --------------------------------------------------------------------------- */

#line 512 "ztprfb.f"
    } else if (column && backward && right) {

/* --------------------------------------------------------------------------- */

/*        Let  W =  [ V ]    (N-by-K) */
/*                  [ I ]    (K-by-K) */

/*        Form  C H  or  C H**H  where  C = [ B A ] (B is M-by-N, A is M-by-K) */

/*        H = I - W T W**H          or  H**H = I - W T**H W**H */

/*        A = A - (A + B V) T      or  A = A - (A + B V) T**H */
/*        B = B - (A + B V) T V**H  or  B = B - (A + B V) T**H V**H */

/* --------------------------------------------------------------------------- */

/* Computing MIN */
#line 528 "ztprfb.f"
	i__1 = *l + 1;
#line 528 "ztprfb.f"
	np = min(i__1,*n);
/* Computing MIN */
#line 529 "ztprfb.f"
	i__1 = *k - *l + 1;
#line 529 "ztprfb.f"
	kp = min(i__1,*k);

#line 531 "ztprfb.f"
	i__1 = *l;
#line 531 "ztprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 532 "ztprfb.f"
	    i__2 = *m;
#line 532 "ztprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 533 "ztprfb.f"
		i__3 = i__ + (*k - *l + j) * work_dim1;
#line 533 "ztprfb.f"
		i__4 = i__ + j * b_dim1;
#line 533 "ztprfb.f"
		work[i__3].r = b[i__4].r, work[i__3].i = b[i__4].i;
#line 534 "ztprfb.f"
	    }
#line 535 "ztprfb.f"
	}
#line 536 "ztprfb.f"
	ztrmm_("R", "L", "N", "N", m, l, &c_b1, &v[kp * v_dim1 + 1], ldv, &
		work[kp * work_dim1 + 1], ldwork, (ftnlen)1, (ftnlen)1, (
		ftnlen)1, (ftnlen)1);
#line 538 "ztprfb.f"
	i__1 = *n - *l;
#line 538 "ztprfb.f"
	zgemm_("N", "N", m, l, &i__1, &c_b1, &b[np * b_dim1 + 1], ldb, &v[np 
		+ kp * v_dim1], ldv, &c_b1, &work[kp * work_dim1 + 1], ldwork,
		 (ftnlen)1, (ftnlen)1);
#line 540 "ztprfb.f"
	i__1 = *k - *l;
#line 540 "ztprfb.f"
	zgemm_("N", "N", m, &i__1, n, &c_b1, &b[b_offset], ldb, &v[v_offset], 
		ldv, &c_b2, &work[work_offset], ldwork, (ftnlen)1, (ftnlen)1);

#line 543 "ztprfb.f"
	i__1 = *k;
#line 543 "ztprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 544 "ztprfb.f"
	    i__2 = *m;
#line 544 "ztprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 545 "ztprfb.f"
		i__3 = i__ + j * work_dim1;
#line 545 "ztprfb.f"
		i__4 = i__ + j * work_dim1;
#line 545 "ztprfb.f"
		i__5 = i__ + j * a_dim1;
#line 545 "ztprfb.f"
		z__1.r = work[i__4].r + a[i__5].r, z__1.i = work[i__4].i + a[
			i__5].i;
#line 545 "ztprfb.f"
		work[i__3].r = z__1.r, work[i__3].i = z__1.i;
#line 546 "ztprfb.f"
	    }
#line 547 "ztprfb.f"
	}

#line 549 "ztprfb.f"
	ztrmm_("R", "L", trans, "N", m, k, &c_b1, &t[t_offset], ldt, &work[
		work_offset], ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		ftnlen)1);

#line 552 "ztprfb.f"
	i__1 = *k;
#line 552 "ztprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 553 "ztprfb.f"
	    i__2 = *m;
#line 553 "ztprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 554 "ztprfb.f"
		i__3 = i__ + j * a_dim1;
#line 554 "ztprfb.f"
		i__4 = i__ + j * a_dim1;
#line 554 "ztprfb.f"
		i__5 = i__ + j * work_dim1;
#line 554 "ztprfb.f"
		z__1.r = a[i__4].r - work[i__5].r, z__1.i = a[i__4].i - work[
			i__5].i;
#line 554 "ztprfb.f"
		a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 555 "ztprfb.f"
	    }
#line 556 "ztprfb.f"
	}

#line 558 "ztprfb.f"
	i__1 = *n - *l;
#line 558 "ztprfb.f"
	z__1.r = -1., z__1.i = -0.;
#line 558 "ztprfb.f"
	zgemm_("N", "C", m, &i__1, k, &z__1, &work[work_offset], ldwork, &v[
		np + v_dim1], ldv, &c_b1, &b[np * b_dim1 + 1], ldb, (ftnlen)1,
		 (ftnlen)1);
#line 560 "ztprfb.f"
	i__1 = *k - *l;
#line 560 "ztprfb.f"
	z__1.r = -1., z__1.i = -0.;
#line 560 "ztprfb.f"
	zgemm_("N", "C", m, l, &i__1, &z__1, &work[work_offset], ldwork, &v[
		v_offset], ldv, &c_b1, &b[b_offset], ldb, (ftnlen)1, (ftnlen)
		1);
#line 562 "ztprfb.f"
	ztrmm_("R", "L", "C", "N", m, l, &c_b1, &v[kp * v_dim1 + 1], ldv, &
		work[kp * work_dim1 + 1], ldwork, (ftnlen)1, (ftnlen)1, (
		ftnlen)1, (ftnlen)1);
#line 564 "ztprfb.f"
	i__1 = *l;
#line 564 "ztprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 565 "ztprfb.f"
	    i__2 = *m;
#line 565 "ztprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 566 "ztprfb.f"
		i__3 = i__ + j * b_dim1;
#line 566 "ztprfb.f"
		i__4 = i__ + j * b_dim1;
#line 566 "ztprfb.f"
		i__5 = i__ + (*k - *l + j) * work_dim1;
#line 566 "ztprfb.f"
		z__1.r = b[i__4].r - work[i__5].r, z__1.i = b[i__4].i - work[
			i__5].i;
#line 566 "ztprfb.f"
		b[i__3].r = z__1.r, b[i__3].i = z__1.i;
#line 567 "ztprfb.f"
	    }
#line 568 "ztprfb.f"
	}

/* --------------------------------------------------------------------------- */

#line 572 "ztprfb.f"
    } else if (row && forward && left) {

/* --------------------------------------------------------------------------- */

/*        Let  W =  [ I V ] ( I is K-by-K, V is K-by-M ) */

/*        Form  H C  or  H**H C  where  C = [ A ]  (K-by-N) */
/*                                          [ B ]  (M-by-N) */

/*        H = I - W**H T W          or  H**H = I - W**H T**H W */

/*        A = A -     T (A + V B)  or  A = A -     T**H (A + V B) */
/*        B = B - V**H T (A + V B)  or  B = B - V**H T**H (A + V B) */

/* --------------------------------------------------------------------------- */

/* Computing MIN */
#line 588 "ztprfb.f"
	i__1 = *m - *l + 1;
#line 588 "ztprfb.f"
	mp = min(i__1,*m);
/* Computing MIN */
#line 589 "ztprfb.f"
	i__1 = *l + 1;
#line 589 "ztprfb.f"
	kp = min(i__1,*k);

#line 591 "ztprfb.f"
	i__1 = *n;
#line 591 "ztprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 592 "ztprfb.f"
	    i__2 = *l;
#line 592 "ztprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 593 "ztprfb.f"
		i__3 = i__ + j * work_dim1;
#line 593 "ztprfb.f"
		i__4 = *m - *l + i__ + j * b_dim1;
#line 593 "ztprfb.f"
		work[i__3].r = b[i__4].r, work[i__3].i = b[i__4].i;
#line 594 "ztprfb.f"
	    }
#line 595 "ztprfb.f"
	}
#line 596 "ztprfb.f"
	ztrmm_("L", "L", "N", "N", l, n, &c_b1, &v[mp * v_dim1 + 1], ldv, &
		work[work_offset], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		ftnlen)1);
#line 598 "ztprfb.f"
	i__1 = *m - *l;
#line 598 "ztprfb.f"
	zgemm_("N", "N", l, n, &i__1, &c_b1, &v[v_offset], ldv, &b[b_offset], 
		ldb, &c_b1, &work[work_offset], ldwork, (ftnlen)1, (ftnlen)1);
#line 600 "ztprfb.f"
	i__1 = *k - *l;
#line 600 "ztprfb.f"
	zgemm_("N", "N", &i__1, n, m, &c_b1, &v[kp + v_dim1], ldv, &b[
		b_offset], ldb, &c_b2, &work[kp + work_dim1], ldwork, (ftnlen)
		1, (ftnlen)1);

#line 603 "ztprfb.f"
	i__1 = *n;
#line 603 "ztprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 604 "ztprfb.f"
	    i__2 = *k;
#line 604 "ztprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 605 "ztprfb.f"
		i__3 = i__ + j * work_dim1;
#line 605 "ztprfb.f"
		i__4 = i__ + j * work_dim1;
#line 605 "ztprfb.f"
		i__5 = i__ + j * a_dim1;
#line 605 "ztprfb.f"
		z__1.r = work[i__4].r + a[i__5].r, z__1.i = work[i__4].i + a[
			i__5].i;
#line 605 "ztprfb.f"
		work[i__3].r = z__1.r, work[i__3].i = z__1.i;
#line 606 "ztprfb.f"
	    }
#line 607 "ztprfb.f"
	}

#line 609 "ztprfb.f"
	ztrmm_("L", "U", trans, "N", k, n, &c_b1, &t[t_offset], ldt, &work[
		work_offset], ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		ftnlen)1);

#line 612 "ztprfb.f"
	i__1 = *n;
#line 612 "ztprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 613 "ztprfb.f"
	    i__2 = *k;
#line 613 "ztprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 614 "ztprfb.f"
		i__3 = i__ + j * a_dim1;
#line 614 "ztprfb.f"
		i__4 = i__ + j * a_dim1;
#line 614 "ztprfb.f"
		i__5 = i__ + j * work_dim1;
#line 614 "ztprfb.f"
		z__1.r = a[i__4].r - work[i__5].r, z__1.i = a[i__4].i - work[
			i__5].i;
#line 614 "ztprfb.f"
		a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 615 "ztprfb.f"
	    }
#line 616 "ztprfb.f"
	}

#line 618 "ztprfb.f"
	i__1 = *m - *l;
#line 618 "ztprfb.f"
	z__1.r = -1., z__1.i = -0.;
#line 618 "ztprfb.f"
	zgemm_("C", "N", &i__1, n, k, &z__1, &v[v_offset], ldv, &work[
		work_offset], ldwork, &c_b1, &b[b_offset], ldb, (ftnlen)1, (
		ftnlen)1);
#line 620 "ztprfb.f"
	i__1 = *k - *l;
#line 620 "ztprfb.f"
	z__1.r = -1., z__1.i = -0.;
#line 620 "ztprfb.f"
	zgemm_("C", "N", l, n, &i__1, &z__1, &v[kp + mp * v_dim1], ldv, &work[
		kp + work_dim1], ldwork, &c_b1, &b[mp + b_dim1], ldb, (ftnlen)
		1, (ftnlen)1);
#line 622 "ztprfb.f"
	ztrmm_("L", "L", "C", "N", l, n, &c_b1, &v[mp * v_dim1 + 1], ldv, &
		work[work_offset], ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		ftnlen)1);
#line 624 "ztprfb.f"
	i__1 = *n;
#line 624 "ztprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 625 "ztprfb.f"
	    i__2 = *l;
#line 625 "ztprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 626 "ztprfb.f"
		i__3 = *m - *l + i__ + j * b_dim1;
#line 626 "ztprfb.f"
		i__4 = *m - *l + i__ + j * b_dim1;
#line 626 "ztprfb.f"
		i__5 = i__ + j * work_dim1;
#line 626 "ztprfb.f"
		z__1.r = b[i__4].r - work[i__5].r, z__1.i = b[i__4].i - work[
			i__5].i;
#line 626 "ztprfb.f"
		b[i__3].r = z__1.r, b[i__3].i = z__1.i;
#line 627 "ztprfb.f"
	    }
#line 628 "ztprfb.f"
	}

/* --------------------------------------------------------------------------- */

#line 632 "ztprfb.f"
    } else if (row && forward && right) {

/* --------------------------------------------------------------------------- */

/*        Let  W =  [ I V ] ( I is K-by-K, V is K-by-N ) */

/*        Form  C H  or  C H**H  where  C = [ A B ] (A is M-by-K, B is M-by-N) */

/*        H = I - W**H T W            or  H**H = I - W**H T**H W */

/*        A = A - (A + B V**H) T      or  A = A - (A + B V**H) T**H */
/*        B = B - (A + B V**H) T V    or  B = B - (A + B V**H) T**H V */

/* --------------------------------------------------------------------------- */

/* Computing MIN */
#line 647 "ztprfb.f"
	i__1 = *n - *l + 1;
#line 647 "ztprfb.f"
	np = min(i__1,*n);
/* Computing MIN */
#line 648 "ztprfb.f"
	i__1 = *l + 1;
#line 648 "ztprfb.f"
	kp = min(i__1,*k);

#line 650 "ztprfb.f"
	i__1 = *l;
#line 650 "ztprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 651 "ztprfb.f"
	    i__2 = *m;
#line 651 "ztprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 652 "ztprfb.f"
		i__3 = i__ + j * work_dim1;
#line 652 "ztprfb.f"
		i__4 = i__ + (*n - *l + j) * b_dim1;
#line 652 "ztprfb.f"
		work[i__3].r = b[i__4].r, work[i__3].i = b[i__4].i;
#line 653 "ztprfb.f"
	    }
#line 654 "ztprfb.f"
	}
#line 655 "ztprfb.f"
	ztrmm_("R", "L", "C", "N", m, l, &c_b1, &v[np * v_dim1 + 1], ldv, &
		work[work_offset], ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		ftnlen)1);
#line 657 "ztprfb.f"
	i__1 = *n - *l;
#line 657 "ztprfb.f"
	zgemm_("N", "C", m, l, &i__1, &c_b1, &b[b_offset], ldb, &v[v_offset], 
		ldv, &c_b1, &work[work_offset], ldwork, (ftnlen)1, (ftnlen)1);
#line 659 "ztprfb.f"
	i__1 = *k - *l;
#line 659 "ztprfb.f"
	zgemm_("N", "C", m, &i__1, n, &c_b1, &b[b_offset], ldb, &v[kp + 
		v_dim1], ldv, &c_b2, &work[kp * work_dim1 + 1], ldwork, (
		ftnlen)1, (ftnlen)1);

#line 662 "ztprfb.f"
	i__1 = *k;
#line 662 "ztprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 663 "ztprfb.f"
	    i__2 = *m;
#line 663 "ztprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 664 "ztprfb.f"
		i__3 = i__ + j * work_dim1;
#line 664 "ztprfb.f"
		i__4 = i__ + j * work_dim1;
#line 664 "ztprfb.f"
		i__5 = i__ + j * a_dim1;
#line 664 "ztprfb.f"
		z__1.r = work[i__4].r + a[i__5].r, z__1.i = work[i__4].i + a[
			i__5].i;
#line 664 "ztprfb.f"
		work[i__3].r = z__1.r, work[i__3].i = z__1.i;
#line 665 "ztprfb.f"
	    }
#line 666 "ztprfb.f"
	}

#line 668 "ztprfb.f"
	ztrmm_("R", "U", trans, "N", m, k, &c_b1, &t[t_offset], ldt, &work[
		work_offset], ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		ftnlen)1);

#line 671 "ztprfb.f"
	i__1 = *k;
#line 671 "ztprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 672 "ztprfb.f"
	    i__2 = *m;
#line 672 "ztprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 673 "ztprfb.f"
		i__3 = i__ + j * a_dim1;
#line 673 "ztprfb.f"
		i__4 = i__ + j * a_dim1;
#line 673 "ztprfb.f"
		i__5 = i__ + j * work_dim1;
#line 673 "ztprfb.f"
		z__1.r = a[i__4].r - work[i__5].r, z__1.i = a[i__4].i - work[
			i__5].i;
#line 673 "ztprfb.f"
		a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 674 "ztprfb.f"
	    }
#line 675 "ztprfb.f"
	}

#line 677 "ztprfb.f"
	i__1 = *n - *l;
#line 677 "ztprfb.f"
	z__1.r = -1., z__1.i = -0.;
#line 677 "ztprfb.f"
	zgemm_("N", "N", m, &i__1, k, &z__1, &work[work_offset], ldwork, &v[
		v_offset], ldv, &c_b1, &b[b_offset], ldb, (ftnlen)1, (ftnlen)
		1);
#line 679 "ztprfb.f"
	i__1 = *k - *l;
#line 679 "ztprfb.f"
	z__1.r = -1., z__1.i = -0.;
#line 679 "ztprfb.f"
	zgemm_("N", "N", m, l, &i__1, &z__1, &work[kp * work_dim1 + 1], 
		ldwork, &v[kp + np * v_dim1], ldv, &c_b1, &b[np * b_dim1 + 1],
		 ldb, (ftnlen)1, (ftnlen)1);
#line 681 "ztprfb.f"
	ztrmm_("R", "L", "N", "N", m, l, &c_b1, &v[np * v_dim1 + 1], ldv, &
		work[work_offset], ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		ftnlen)1);
#line 683 "ztprfb.f"
	i__1 = *l;
#line 683 "ztprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 684 "ztprfb.f"
	    i__2 = *m;
#line 684 "ztprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 685 "ztprfb.f"
		i__3 = i__ + (*n - *l + j) * b_dim1;
#line 685 "ztprfb.f"
		i__4 = i__ + (*n - *l + j) * b_dim1;
#line 685 "ztprfb.f"
		i__5 = i__ + j * work_dim1;
#line 685 "ztprfb.f"
		z__1.r = b[i__4].r - work[i__5].r, z__1.i = b[i__4].i - work[
			i__5].i;
#line 685 "ztprfb.f"
		b[i__3].r = z__1.r, b[i__3].i = z__1.i;
#line 686 "ztprfb.f"
	    }
#line 687 "ztprfb.f"
	}

/* --------------------------------------------------------------------------- */

#line 691 "ztprfb.f"
    } else if (row && backward && left) {

/* --------------------------------------------------------------------------- */

/*        Let  W =  [ V I ] ( I is K-by-K, V is K-by-M ) */

/*        Form  H C  or  H**H C  where  C = [ B ]  (M-by-N) */
/*                                          [ A ]  (K-by-N) */

/*        H = I - W**H T W          or  H**H = I - W**H T**H W */

/*        A = A -     T (A + V B)  or  A = A -     T**H (A + V B) */
/*        B = B - V**H T (A + V B)  or  B = B - V**H T**H (A + V B) */

/* --------------------------------------------------------------------------- */

/* Computing MIN */
#line 707 "ztprfb.f"
	i__1 = *l + 1;
#line 707 "ztprfb.f"
	mp = min(i__1,*m);
/* Computing MIN */
#line 708 "ztprfb.f"
	i__1 = *k - *l + 1;
#line 708 "ztprfb.f"
	kp = min(i__1,*k);

#line 710 "ztprfb.f"
	i__1 = *n;
#line 710 "ztprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 711 "ztprfb.f"
	    i__2 = *l;
#line 711 "ztprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 712 "ztprfb.f"
		i__3 = *k - *l + i__ + j * work_dim1;
#line 712 "ztprfb.f"
		i__4 = i__ + j * b_dim1;
#line 712 "ztprfb.f"
		work[i__3].r = b[i__4].r, work[i__3].i = b[i__4].i;
#line 713 "ztprfb.f"
	    }
#line 714 "ztprfb.f"
	}
#line 715 "ztprfb.f"
	ztrmm_("L", "U", "N", "N", l, n, &c_b1, &v[kp + v_dim1], ldv, &work[
		kp + work_dim1], ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		ftnlen)1);
#line 717 "ztprfb.f"
	i__1 = *m - *l;
#line 717 "ztprfb.f"
	zgemm_("N", "N", l, n, &i__1, &c_b1, &v[kp + mp * v_dim1], ldv, &b[mp 
		+ b_dim1], ldb, &c_b1, &work[kp + work_dim1], ldwork, (ftnlen)
		1, (ftnlen)1);
#line 719 "ztprfb.f"
	i__1 = *k - *l;
#line 719 "ztprfb.f"
	zgemm_("N", "N", &i__1, n, m, &c_b1, &v[v_offset], ldv, &b[b_offset], 
		ldb, &c_b2, &work[work_offset], ldwork, (ftnlen)1, (ftnlen)1);

#line 722 "ztprfb.f"
	i__1 = *n;
#line 722 "ztprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 723 "ztprfb.f"
	    i__2 = *k;
#line 723 "ztprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 724 "ztprfb.f"
		i__3 = i__ + j * work_dim1;
#line 724 "ztprfb.f"
		i__4 = i__ + j * work_dim1;
#line 724 "ztprfb.f"
		i__5 = i__ + j * a_dim1;
#line 724 "ztprfb.f"
		z__1.r = work[i__4].r + a[i__5].r, z__1.i = work[i__4].i + a[
			i__5].i;
#line 724 "ztprfb.f"
		work[i__3].r = z__1.r, work[i__3].i = z__1.i;
#line 725 "ztprfb.f"
	    }
#line 726 "ztprfb.f"
	}

#line 728 "ztprfb.f"
	ztrmm_("L", "L ", trans, "N", k, n, &c_b1, &t[t_offset], ldt, &work[
		work_offset], ldwork, (ftnlen)1, (ftnlen)2, (ftnlen)1, (
		ftnlen)1);

#line 731 "ztprfb.f"
	i__1 = *n;
#line 731 "ztprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 732 "ztprfb.f"
	    i__2 = *k;
#line 732 "ztprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 733 "ztprfb.f"
		i__3 = i__ + j * a_dim1;
#line 733 "ztprfb.f"
		i__4 = i__ + j * a_dim1;
#line 733 "ztprfb.f"
		i__5 = i__ + j * work_dim1;
#line 733 "ztprfb.f"
		z__1.r = a[i__4].r - work[i__5].r, z__1.i = a[i__4].i - work[
			i__5].i;
#line 733 "ztprfb.f"
		a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 734 "ztprfb.f"
	    }
#line 735 "ztprfb.f"
	}

#line 737 "ztprfb.f"
	i__1 = *m - *l;
#line 737 "ztprfb.f"
	z__1.r = -1., z__1.i = -0.;
#line 737 "ztprfb.f"
	zgemm_("C", "N", &i__1, n, k, &z__1, &v[mp * v_dim1 + 1], ldv, &work[
		work_offset], ldwork, &c_b1, &b[mp + b_dim1], ldb, (ftnlen)1, 
		(ftnlen)1);
#line 739 "ztprfb.f"
	i__1 = *k - *l;
#line 739 "ztprfb.f"
	z__1.r = -1., z__1.i = -0.;
#line 739 "ztprfb.f"
	zgemm_("C", "N", l, n, &i__1, &z__1, &v[v_offset], ldv, &work[
		work_offset], ldwork, &c_b1, &b[b_offset], ldb, (ftnlen)1, (
		ftnlen)1);
#line 741 "ztprfb.f"
	ztrmm_("L", "U", "C", "N", l, n, &c_b1, &v[kp + v_dim1], ldv, &work[
		kp + work_dim1], ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		ftnlen)1);
#line 743 "ztprfb.f"
	i__1 = *n;
#line 743 "ztprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 744 "ztprfb.f"
	    i__2 = *l;
#line 744 "ztprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 745 "ztprfb.f"
		i__3 = i__ + j * b_dim1;
#line 745 "ztprfb.f"
		i__4 = i__ + j * b_dim1;
#line 745 "ztprfb.f"
		i__5 = *k - *l + i__ + j * work_dim1;
#line 745 "ztprfb.f"
		z__1.r = b[i__4].r - work[i__5].r, z__1.i = b[i__4].i - work[
			i__5].i;
#line 745 "ztprfb.f"
		b[i__3].r = z__1.r, b[i__3].i = z__1.i;
#line 746 "ztprfb.f"
	    }
#line 747 "ztprfb.f"
	}

/* --------------------------------------------------------------------------- */

#line 751 "ztprfb.f"
    } else if (row && backward && right) {

/* --------------------------------------------------------------------------- */

/*        Let  W =  [ V I ] ( I is K-by-K, V is K-by-N ) */

/*        Form  C H  or  C H**H  where  C = [ B A ] (A is M-by-K, B is M-by-N) */

/*        H = I - W**H T W            or  H**H = I - W**H T**H W */

/*        A = A - (A + B V**H) T      or  A = A - (A + B V**H) T**H */
/*        B = B - (A + B V**H) T V    or  B = B - (A + B V**H) T**H V */

/* --------------------------------------------------------------------------- */

/* Computing MIN */
#line 766 "ztprfb.f"
	i__1 = *l + 1;
#line 766 "ztprfb.f"
	np = min(i__1,*n);
/* Computing MIN */
#line 767 "ztprfb.f"
	i__1 = *k - *l + 1;
#line 767 "ztprfb.f"
	kp = min(i__1,*k);

#line 769 "ztprfb.f"
	i__1 = *l;
#line 769 "ztprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 770 "ztprfb.f"
	    i__2 = *m;
#line 770 "ztprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 771 "ztprfb.f"
		i__3 = i__ + (*k - *l + j) * work_dim1;
#line 771 "ztprfb.f"
		i__4 = i__ + j * b_dim1;
#line 771 "ztprfb.f"
		work[i__3].r = b[i__4].r, work[i__3].i = b[i__4].i;
#line 772 "ztprfb.f"
	    }
#line 773 "ztprfb.f"
	}
#line 774 "ztprfb.f"
	ztrmm_("R", "U", "C", "N", m, l, &c_b1, &v[kp + v_dim1], ldv, &work[
		kp * work_dim1 + 1], ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, 
		(ftnlen)1);
#line 776 "ztprfb.f"
	i__1 = *n - *l;
#line 776 "ztprfb.f"
	zgemm_("N", "C", m, l, &i__1, &c_b1, &b[np * b_dim1 + 1], ldb, &v[kp 
		+ np * v_dim1], ldv, &c_b1, &work[kp * work_dim1 + 1], ldwork,
		 (ftnlen)1, (ftnlen)1);
#line 778 "ztprfb.f"
	i__1 = *k - *l;
#line 778 "ztprfb.f"
	zgemm_("N", "C", m, &i__1, n, &c_b1, &b[b_offset], ldb, &v[v_offset], 
		ldv, &c_b2, &work[work_offset], ldwork, (ftnlen)1, (ftnlen)1);

#line 781 "ztprfb.f"
	i__1 = *k;
#line 781 "ztprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 782 "ztprfb.f"
	    i__2 = *m;
#line 782 "ztprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 783 "ztprfb.f"
		i__3 = i__ + j * work_dim1;
#line 783 "ztprfb.f"
		i__4 = i__ + j * work_dim1;
#line 783 "ztprfb.f"
		i__5 = i__ + j * a_dim1;
#line 783 "ztprfb.f"
		z__1.r = work[i__4].r + a[i__5].r, z__1.i = work[i__4].i + a[
			i__5].i;
#line 783 "ztprfb.f"
		work[i__3].r = z__1.r, work[i__3].i = z__1.i;
#line 784 "ztprfb.f"
	    }
#line 785 "ztprfb.f"
	}

#line 787 "ztprfb.f"
	ztrmm_("R", "L", trans, "N", m, k, &c_b1, &t[t_offset], ldt, &work[
		work_offset], ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		ftnlen)1);

#line 790 "ztprfb.f"
	i__1 = *k;
#line 790 "ztprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 791 "ztprfb.f"
	    i__2 = *m;
#line 791 "ztprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 792 "ztprfb.f"
		i__3 = i__ + j * a_dim1;
#line 792 "ztprfb.f"
		i__4 = i__ + j * a_dim1;
#line 792 "ztprfb.f"
		i__5 = i__ + j * work_dim1;
#line 792 "ztprfb.f"
		z__1.r = a[i__4].r - work[i__5].r, z__1.i = a[i__4].i - work[
			i__5].i;
#line 792 "ztprfb.f"
		a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 793 "ztprfb.f"
	    }
#line 794 "ztprfb.f"
	}

#line 796 "ztprfb.f"
	i__1 = *n - *l;
#line 796 "ztprfb.f"
	z__1.r = -1., z__1.i = -0.;
#line 796 "ztprfb.f"
	zgemm_("N", "N", m, &i__1, k, &z__1, &work[work_offset], ldwork, &v[
		np * v_dim1 + 1], ldv, &c_b1, &b[np * b_dim1 + 1], ldb, (
		ftnlen)1, (ftnlen)1);
#line 798 "ztprfb.f"
	i__1 = *k - *l;
#line 798 "ztprfb.f"
	z__1.r = -1., z__1.i = -0.;
#line 798 "ztprfb.f"
	zgemm_("N", "N", m, l, &i__1, &z__1, &work[work_offset], ldwork, &v[
		v_offset], ldv, &c_b1, &b[b_offset], ldb, (ftnlen)1, (ftnlen)
		1);
#line 800 "ztprfb.f"
	ztrmm_("R", "U", "N", "N", m, l, &c_b1, &v[kp + v_dim1], ldv, &work[
		kp * work_dim1 + 1], ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, 
		(ftnlen)1);
#line 802 "ztprfb.f"
	i__1 = *l;
#line 802 "ztprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 803 "ztprfb.f"
	    i__2 = *m;
#line 803 "ztprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 804 "ztprfb.f"
		i__3 = i__ + j * b_dim1;
#line 804 "ztprfb.f"
		i__4 = i__ + j * b_dim1;
#line 804 "ztprfb.f"
		i__5 = i__ + (*k - *l + j) * work_dim1;
#line 804 "ztprfb.f"
		z__1.r = b[i__4].r - work[i__5].r, z__1.i = b[i__4].i - work[
			i__5].i;
#line 804 "ztprfb.f"
		b[i__3].r = z__1.r, b[i__3].i = z__1.i;
#line 805 "ztprfb.f"
	    }
#line 806 "ztprfb.f"
	}

#line 808 "ztprfb.f"
    }

#line 810 "ztprfb.f"
    return 0;

/*     End of ZTPRFB */

} /* ztprfb_ */

