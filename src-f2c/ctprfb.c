#line 1 "ctprfb.f"
/* ctprfb.f -- translated by f2c (version 20100827).
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

#line 1 "ctprfb.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static doublecomplex c_b2 = {0.,0.};

/* > \brief \b CTPRFB applies a real or complex "triangular-pentagonal" blocked reflector to a real or complex
 matrix, which is composed of two blocks. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CTPRFB + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ctprfb.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ctprfb.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ctprfb.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CTPRFB( SIDE, TRANS, DIRECT, STOREV, M, N, K, L, */
/*                          V, LDV, T, LDT, A, LDA, B, LDB, WORK, LDWORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER DIRECT, SIDE, STOREV, TRANS */
/*       INTEGER   K, L, LDA, LDB, LDT, LDV, LDWORK, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX   A( LDA, * ), B( LDB, * ), T( LDT, * ), */
/*      $          V( LDV, * ), WORK( LDWORK, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CTPRFB applies a complex "triangular-pentagonal" block reflector H or its */
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
/* >          V is COMPLEX array, dimension */
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
/* >          T is COMPLEX array, dimension (LDT,K) */
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
/* >          A is COMPLEX array, dimension */
/* >          (LDA,N) if SIDE = 'L' or (LDA,K) if SIDE = 'R' */
/* >          On entry, the K-by-N or M-by-K matrix A. */
/* >          On exit, A is overwritten by the corresponding block of */
/* >          H*C or H**H*C or C*H or C*H**H.  See Further Details. */
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
/* >          B is COMPLEX array, dimension (LDB,N) */
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
/* >          WORK is COMPLEX array, dimension */
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

/* > \date December 2016 */

/* > \ingroup complexOTHERauxiliary */

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
/* Subroutine */ int ctprfb_(char *side, char *trans, char *direct, char *
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
    extern /* Subroutine */ int cgemm_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical right;
    extern /* Subroutine */ int ctrmm_(char *, char *, char *, char *, 
	    integer *, integer *, doublecomplex *, doublecomplex *, integer *,
	     doublecomplex *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static logical column, forward;


/*  -- LAPACK auxiliary routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

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

#line 292 "ctprfb.f"
    /* Parameter adjustments */
#line 292 "ctprfb.f"
    v_dim1 = *ldv;
#line 292 "ctprfb.f"
    v_offset = 1 + v_dim1;
#line 292 "ctprfb.f"
    v -= v_offset;
#line 292 "ctprfb.f"
    t_dim1 = *ldt;
#line 292 "ctprfb.f"
    t_offset = 1 + t_dim1;
#line 292 "ctprfb.f"
    t -= t_offset;
#line 292 "ctprfb.f"
    a_dim1 = *lda;
#line 292 "ctprfb.f"
    a_offset = 1 + a_dim1;
#line 292 "ctprfb.f"
    a -= a_offset;
#line 292 "ctprfb.f"
    b_dim1 = *ldb;
#line 292 "ctprfb.f"
    b_offset = 1 + b_dim1;
#line 292 "ctprfb.f"
    b -= b_offset;
#line 292 "ctprfb.f"
    work_dim1 = *ldwork;
#line 292 "ctprfb.f"
    work_offset = 1 + work_dim1;
#line 292 "ctprfb.f"
    work -= work_offset;
#line 292 "ctprfb.f"

#line 292 "ctprfb.f"
    /* Function Body */
#line 292 "ctprfb.f"
    if (*m <= 0 || *n <= 0 || *k <= 0 || *l < 0) {
#line 292 "ctprfb.f"
	return 0;
#line 292 "ctprfb.f"
    }

#line 294 "ctprfb.f"
    if (lsame_(storev, "C", (ftnlen)1, (ftnlen)1)) {
#line 295 "ctprfb.f"
	column = TRUE_;
#line 296 "ctprfb.f"
	row = FALSE_;
#line 297 "ctprfb.f"
    } else if (lsame_(storev, "R", (ftnlen)1, (ftnlen)1)) {
#line 298 "ctprfb.f"
	column = FALSE_;
#line 299 "ctprfb.f"
	row = TRUE_;
#line 300 "ctprfb.f"
    } else {
#line 301 "ctprfb.f"
	column = FALSE_;
#line 302 "ctprfb.f"
	row = FALSE_;
#line 303 "ctprfb.f"
    }

#line 305 "ctprfb.f"
    if (lsame_(side, "L", (ftnlen)1, (ftnlen)1)) {
#line 306 "ctprfb.f"
	left = TRUE_;
#line 307 "ctprfb.f"
	right = FALSE_;
#line 308 "ctprfb.f"
    } else if (lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {
#line 309 "ctprfb.f"
	left = FALSE_;
#line 310 "ctprfb.f"
	right = TRUE_;
#line 311 "ctprfb.f"
    } else {
#line 312 "ctprfb.f"
	left = FALSE_;
#line 313 "ctprfb.f"
	right = FALSE_;
#line 314 "ctprfb.f"
    }

#line 316 "ctprfb.f"
    if (lsame_(direct, "F", (ftnlen)1, (ftnlen)1)) {
#line 317 "ctprfb.f"
	forward = TRUE_;
#line 318 "ctprfb.f"
	backward = FALSE_;
#line 319 "ctprfb.f"
    } else if (lsame_(direct, "B", (ftnlen)1, (ftnlen)1)) {
#line 320 "ctprfb.f"
	forward = FALSE_;
#line 321 "ctprfb.f"
	backward = TRUE_;
#line 322 "ctprfb.f"
    } else {
#line 323 "ctprfb.f"
	forward = FALSE_;
#line 324 "ctprfb.f"
	backward = FALSE_;
#line 325 "ctprfb.f"
    }

/* --------------------------------------------------------------------------- */

#line 329 "ctprfb.f"
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
#line 346 "ctprfb.f"
	i__1 = *m - *l + 1;
#line 346 "ctprfb.f"
	mp = min(i__1,*m);
/* Computing MIN */
#line 347 "ctprfb.f"
	i__1 = *l + 1;
#line 347 "ctprfb.f"
	kp = min(i__1,*k);

#line 349 "ctprfb.f"
	i__1 = *n;
#line 349 "ctprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 350 "ctprfb.f"
	    i__2 = *l;
#line 350 "ctprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 351 "ctprfb.f"
		i__3 = i__ + j * work_dim1;
#line 351 "ctprfb.f"
		i__4 = *m - *l + i__ + j * b_dim1;
#line 351 "ctprfb.f"
		work[i__3].r = b[i__4].r, work[i__3].i = b[i__4].i;
#line 352 "ctprfb.f"
	    }
#line 353 "ctprfb.f"
	}
#line 354 "ctprfb.f"
	ctrmm_("L", "U", "C", "N", l, n, &c_b1, &v[mp + v_dim1], ldv, &work[
		work_offset], ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		ftnlen)1);
#line 356 "ctprfb.f"
	i__1 = *m - *l;
#line 356 "ctprfb.f"
	cgemm_("C", "N", l, n, &i__1, &c_b1, &v[v_offset], ldv, &b[b_offset], 
		ldb, &c_b1, &work[work_offset], ldwork, (ftnlen)1, (ftnlen)1);
#line 358 "ctprfb.f"
	i__1 = *k - *l;
#line 358 "ctprfb.f"
	cgemm_("C", "N", &i__1, n, m, &c_b1, &v[kp * v_dim1 + 1], ldv, &b[
		b_offset], ldb, &c_b2, &work[kp + work_dim1], ldwork, (ftnlen)
		1, (ftnlen)1);

#line 361 "ctprfb.f"
	i__1 = *n;
#line 361 "ctprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 362 "ctprfb.f"
	    i__2 = *k;
#line 362 "ctprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 363 "ctprfb.f"
		i__3 = i__ + j * work_dim1;
#line 363 "ctprfb.f"
		i__4 = i__ + j * work_dim1;
#line 363 "ctprfb.f"
		i__5 = i__ + j * a_dim1;
#line 363 "ctprfb.f"
		z__1.r = work[i__4].r + a[i__5].r, z__1.i = work[i__4].i + a[
			i__5].i;
#line 363 "ctprfb.f"
		work[i__3].r = z__1.r, work[i__3].i = z__1.i;
#line 364 "ctprfb.f"
	    }
#line 365 "ctprfb.f"
	}

#line 367 "ctprfb.f"
	ctrmm_("L", "U", trans, "N", k, n, &c_b1, &t[t_offset], ldt, &work[
		work_offset], ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		ftnlen)1);

#line 370 "ctprfb.f"
	i__1 = *n;
#line 370 "ctprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 371 "ctprfb.f"
	    i__2 = *k;
#line 371 "ctprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 372 "ctprfb.f"
		i__3 = i__ + j * a_dim1;
#line 372 "ctprfb.f"
		i__4 = i__ + j * a_dim1;
#line 372 "ctprfb.f"
		i__5 = i__ + j * work_dim1;
#line 372 "ctprfb.f"
		z__1.r = a[i__4].r - work[i__5].r, z__1.i = a[i__4].i - work[
			i__5].i;
#line 372 "ctprfb.f"
		a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 373 "ctprfb.f"
	    }
#line 374 "ctprfb.f"
	}

#line 376 "ctprfb.f"
	i__1 = *m - *l;
#line 376 "ctprfb.f"
	z__1.r = -1., z__1.i = -0.;
#line 376 "ctprfb.f"
	cgemm_("N", "N", &i__1, n, k, &z__1, &v[v_offset], ldv, &work[
		work_offset], ldwork, &c_b1, &b[b_offset], ldb, (ftnlen)1, (
		ftnlen)1);
#line 378 "ctprfb.f"
	i__1 = *k - *l;
#line 378 "ctprfb.f"
	z__1.r = -1., z__1.i = -0.;
#line 378 "ctprfb.f"
	cgemm_("N", "N", l, n, &i__1, &z__1, &v[mp + kp * v_dim1], ldv, &work[
		kp + work_dim1], ldwork, &c_b1, &b[mp + b_dim1], ldb, (ftnlen)
		1, (ftnlen)1);
#line 380 "ctprfb.f"
	ctrmm_("L", "U", "N", "N", l, n, &c_b1, &v[mp + v_dim1], ldv, &work[
		work_offset], ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		ftnlen)1);
#line 382 "ctprfb.f"
	i__1 = *n;
#line 382 "ctprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 383 "ctprfb.f"
	    i__2 = *l;
#line 383 "ctprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 384 "ctprfb.f"
		i__3 = *m - *l + i__ + j * b_dim1;
#line 384 "ctprfb.f"
		i__4 = *m - *l + i__ + j * b_dim1;
#line 384 "ctprfb.f"
		i__5 = i__ + j * work_dim1;
#line 384 "ctprfb.f"
		z__1.r = b[i__4].r - work[i__5].r, z__1.i = b[i__4].i - work[
			i__5].i;
#line 384 "ctprfb.f"
		b[i__3].r = z__1.r, b[i__3].i = z__1.i;
#line 385 "ctprfb.f"
	    }
#line 386 "ctprfb.f"
	}

/* --------------------------------------------------------------------------- */

#line 390 "ctprfb.f"
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
#line 406 "ctprfb.f"
	i__1 = *n - *l + 1;
#line 406 "ctprfb.f"
	np = min(i__1,*n);
/* Computing MIN */
#line 407 "ctprfb.f"
	i__1 = *l + 1;
#line 407 "ctprfb.f"
	kp = min(i__1,*k);

#line 409 "ctprfb.f"
	i__1 = *l;
#line 409 "ctprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 410 "ctprfb.f"
	    i__2 = *m;
#line 410 "ctprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 411 "ctprfb.f"
		i__3 = i__ + j * work_dim1;
#line 411 "ctprfb.f"
		i__4 = i__ + (*n - *l + j) * b_dim1;
#line 411 "ctprfb.f"
		work[i__3].r = b[i__4].r, work[i__3].i = b[i__4].i;
#line 412 "ctprfb.f"
	    }
#line 413 "ctprfb.f"
	}
#line 414 "ctprfb.f"
	ctrmm_("R", "U", "N", "N", m, l, &c_b1, &v[np + v_dim1], ldv, &work[
		work_offset], ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		ftnlen)1);
#line 416 "ctprfb.f"
	i__1 = *n - *l;
#line 416 "ctprfb.f"
	cgemm_("N", "N", m, l, &i__1, &c_b1, &b[b_offset], ldb, &v[v_offset], 
		ldv, &c_b1, &work[work_offset], ldwork, (ftnlen)1, (ftnlen)1);
#line 418 "ctprfb.f"
	i__1 = *k - *l;
#line 418 "ctprfb.f"
	cgemm_("N", "N", m, &i__1, n, &c_b1, &b[b_offset], ldb, &v[kp * 
		v_dim1 + 1], ldv, &c_b2, &work[kp * work_dim1 + 1], ldwork, (
		ftnlen)1, (ftnlen)1);

#line 421 "ctprfb.f"
	i__1 = *k;
#line 421 "ctprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 422 "ctprfb.f"
	    i__2 = *m;
#line 422 "ctprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 423 "ctprfb.f"
		i__3 = i__ + j * work_dim1;
#line 423 "ctprfb.f"
		i__4 = i__ + j * work_dim1;
#line 423 "ctprfb.f"
		i__5 = i__ + j * a_dim1;
#line 423 "ctprfb.f"
		z__1.r = work[i__4].r + a[i__5].r, z__1.i = work[i__4].i + a[
			i__5].i;
#line 423 "ctprfb.f"
		work[i__3].r = z__1.r, work[i__3].i = z__1.i;
#line 424 "ctprfb.f"
	    }
#line 425 "ctprfb.f"
	}

#line 427 "ctprfb.f"
	ctrmm_("R", "U", trans, "N", m, k, &c_b1, &t[t_offset], ldt, &work[
		work_offset], ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		ftnlen)1);

#line 430 "ctprfb.f"
	i__1 = *k;
#line 430 "ctprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 431 "ctprfb.f"
	    i__2 = *m;
#line 431 "ctprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 432 "ctprfb.f"
		i__3 = i__ + j * a_dim1;
#line 432 "ctprfb.f"
		i__4 = i__ + j * a_dim1;
#line 432 "ctprfb.f"
		i__5 = i__ + j * work_dim1;
#line 432 "ctprfb.f"
		z__1.r = a[i__4].r - work[i__5].r, z__1.i = a[i__4].i - work[
			i__5].i;
#line 432 "ctprfb.f"
		a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 433 "ctprfb.f"
	    }
#line 434 "ctprfb.f"
	}

#line 436 "ctprfb.f"
	i__1 = *n - *l;
#line 436 "ctprfb.f"
	z__1.r = -1., z__1.i = -0.;
#line 436 "ctprfb.f"
	cgemm_("N", "C", m, &i__1, k, &z__1, &work[work_offset], ldwork, &v[
		v_offset], ldv, &c_b1, &b[b_offset], ldb, (ftnlen)1, (ftnlen)
		1);
#line 438 "ctprfb.f"
	i__1 = *k - *l;
#line 438 "ctprfb.f"
	z__1.r = -1., z__1.i = -0.;
#line 438 "ctprfb.f"
	cgemm_("N", "C", m, l, &i__1, &z__1, &work[kp * work_dim1 + 1], 
		ldwork, &v[np + kp * v_dim1], ldv, &c_b1, &b[np * b_dim1 + 1],
		 ldb, (ftnlen)1, (ftnlen)1);
#line 440 "ctprfb.f"
	ctrmm_("R", "U", "C", "N", m, l, &c_b1, &v[np + v_dim1], ldv, &work[
		work_offset], ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		ftnlen)1);
#line 442 "ctprfb.f"
	i__1 = *l;
#line 442 "ctprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 443 "ctprfb.f"
	    i__2 = *m;
#line 443 "ctprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 444 "ctprfb.f"
		i__3 = i__ + (*n - *l + j) * b_dim1;
#line 444 "ctprfb.f"
		i__4 = i__ + (*n - *l + j) * b_dim1;
#line 444 "ctprfb.f"
		i__5 = i__ + j * work_dim1;
#line 444 "ctprfb.f"
		z__1.r = b[i__4].r - work[i__5].r, z__1.i = b[i__4].i - work[
			i__5].i;
#line 444 "ctprfb.f"
		b[i__3].r = z__1.r, b[i__3].i = z__1.i;
#line 445 "ctprfb.f"
	    }
#line 446 "ctprfb.f"
	}

/* --------------------------------------------------------------------------- */

#line 450 "ctprfb.f"
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
#line 467 "ctprfb.f"
	i__1 = *l + 1;
#line 467 "ctprfb.f"
	mp = min(i__1,*m);
/* Computing MIN */
#line 468 "ctprfb.f"
	i__1 = *k - *l + 1;
#line 468 "ctprfb.f"
	kp = min(i__1,*k);

#line 470 "ctprfb.f"
	i__1 = *n;
#line 470 "ctprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 471 "ctprfb.f"
	    i__2 = *l;
#line 471 "ctprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 472 "ctprfb.f"
		i__3 = *k - *l + i__ + j * work_dim1;
#line 472 "ctprfb.f"
		i__4 = i__ + j * b_dim1;
#line 472 "ctprfb.f"
		work[i__3].r = b[i__4].r, work[i__3].i = b[i__4].i;
#line 473 "ctprfb.f"
	    }
#line 474 "ctprfb.f"
	}

#line 476 "ctprfb.f"
	ctrmm_("L", "L", "C", "N", l, n, &c_b1, &v[kp * v_dim1 + 1], ldv, &
		work[kp + work_dim1], ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1,
		 (ftnlen)1);
#line 478 "ctprfb.f"
	i__1 = *m - *l;
#line 478 "ctprfb.f"
	cgemm_("C", "N", l, n, &i__1, &c_b1, &v[mp + kp * v_dim1], ldv, &b[mp 
		+ b_dim1], ldb, &c_b1, &work[kp + work_dim1], ldwork, (ftnlen)
		1, (ftnlen)1);
#line 480 "ctprfb.f"
	i__1 = *k - *l;
#line 480 "ctprfb.f"
	cgemm_("C", "N", &i__1, n, m, &c_b1, &v[v_offset], ldv, &b[b_offset], 
		ldb, &c_b2, &work[work_offset], ldwork, (ftnlen)1, (ftnlen)1);

#line 483 "ctprfb.f"
	i__1 = *n;
#line 483 "ctprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 484 "ctprfb.f"
	    i__2 = *k;
#line 484 "ctprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 485 "ctprfb.f"
		i__3 = i__ + j * work_dim1;
#line 485 "ctprfb.f"
		i__4 = i__ + j * work_dim1;
#line 485 "ctprfb.f"
		i__5 = i__ + j * a_dim1;
#line 485 "ctprfb.f"
		z__1.r = work[i__4].r + a[i__5].r, z__1.i = work[i__4].i + a[
			i__5].i;
#line 485 "ctprfb.f"
		work[i__3].r = z__1.r, work[i__3].i = z__1.i;
#line 486 "ctprfb.f"
	    }
#line 487 "ctprfb.f"
	}

#line 489 "ctprfb.f"
	ctrmm_("L", "L", trans, "N", k, n, &c_b1, &t[t_offset], ldt, &work[
		work_offset], ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		ftnlen)1);

#line 492 "ctprfb.f"
	i__1 = *n;
#line 492 "ctprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 493 "ctprfb.f"
	    i__2 = *k;
#line 493 "ctprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 494 "ctprfb.f"
		i__3 = i__ + j * a_dim1;
#line 494 "ctprfb.f"
		i__4 = i__ + j * a_dim1;
#line 494 "ctprfb.f"
		i__5 = i__ + j * work_dim1;
#line 494 "ctprfb.f"
		z__1.r = a[i__4].r - work[i__5].r, z__1.i = a[i__4].i - work[
			i__5].i;
#line 494 "ctprfb.f"
		a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 495 "ctprfb.f"
	    }
#line 496 "ctprfb.f"
	}

#line 498 "ctprfb.f"
	i__1 = *m - *l;
#line 498 "ctprfb.f"
	z__1.r = -1., z__1.i = -0.;
#line 498 "ctprfb.f"
	cgemm_("N", "N", &i__1, n, k, &z__1, &v[mp + v_dim1], ldv, &work[
		work_offset], ldwork, &c_b1, &b[mp + b_dim1], ldb, (ftnlen)1, 
		(ftnlen)1);
#line 500 "ctprfb.f"
	i__1 = *k - *l;
#line 500 "ctprfb.f"
	z__1.r = -1., z__1.i = -0.;
#line 500 "ctprfb.f"
	cgemm_("N", "N", l, n, &i__1, &z__1, &v[v_offset], ldv, &work[
		work_offset], ldwork, &c_b1, &b[b_offset], ldb, (ftnlen)1, (
		ftnlen)1);
#line 502 "ctprfb.f"
	ctrmm_("L", "L", "N", "N", l, n, &c_b1, &v[kp * v_dim1 + 1], ldv, &
		work[kp + work_dim1], ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1,
		 (ftnlen)1);
#line 504 "ctprfb.f"
	i__1 = *n;
#line 504 "ctprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 505 "ctprfb.f"
	    i__2 = *l;
#line 505 "ctprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 506 "ctprfb.f"
		i__3 = i__ + j * b_dim1;
#line 506 "ctprfb.f"
		i__4 = i__ + j * b_dim1;
#line 506 "ctprfb.f"
		i__5 = *k - *l + i__ + j * work_dim1;
#line 506 "ctprfb.f"
		z__1.r = b[i__4].r - work[i__5].r, z__1.i = b[i__4].i - work[
			i__5].i;
#line 506 "ctprfb.f"
		b[i__3].r = z__1.r, b[i__3].i = z__1.i;
#line 507 "ctprfb.f"
	    }
#line 508 "ctprfb.f"
	}

/* --------------------------------------------------------------------------- */

#line 512 "ctprfb.f"
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
#line 528 "ctprfb.f"
	i__1 = *l + 1;
#line 528 "ctprfb.f"
	np = min(i__1,*n);
/* Computing MIN */
#line 529 "ctprfb.f"
	i__1 = *k - *l + 1;
#line 529 "ctprfb.f"
	kp = min(i__1,*k);

#line 531 "ctprfb.f"
	i__1 = *l;
#line 531 "ctprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 532 "ctprfb.f"
	    i__2 = *m;
#line 532 "ctprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 533 "ctprfb.f"
		i__3 = i__ + (*k - *l + j) * work_dim1;
#line 533 "ctprfb.f"
		i__4 = i__ + j * b_dim1;
#line 533 "ctprfb.f"
		work[i__3].r = b[i__4].r, work[i__3].i = b[i__4].i;
#line 534 "ctprfb.f"
	    }
#line 535 "ctprfb.f"
	}
#line 536 "ctprfb.f"
	ctrmm_("R", "L", "N", "N", m, l, &c_b1, &v[kp * v_dim1 + 1], ldv, &
		work[kp * work_dim1 + 1], ldwork, (ftnlen)1, (ftnlen)1, (
		ftnlen)1, (ftnlen)1);
#line 538 "ctprfb.f"
	i__1 = *n - *l;
#line 538 "ctprfb.f"
	cgemm_("N", "N", m, l, &i__1, &c_b1, &b[np * b_dim1 + 1], ldb, &v[np 
		+ kp * v_dim1], ldv, &c_b1, &work[kp * work_dim1 + 1], ldwork,
		 (ftnlen)1, (ftnlen)1);
#line 540 "ctprfb.f"
	i__1 = *k - *l;
#line 540 "ctprfb.f"
	cgemm_("N", "N", m, &i__1, n, &c_b1, &b[b_offset], ldb, &v[v_offset], 
		ldv, &c_b2, &work[work_offset], ldwork, (ftnlen)1, (ftnlen)1);

#line 543 "ctprfb.f"
	i__1 = *k;
#line 543 "ctprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 544 "ctprfb.f"
	    i__2 = *m;
#line 544 "ctprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 545 "ctprfb.f"
		i__3 = i__ + j * work_dim1;
#line 545 "ctprfb.f"
		i__4 = i__ + j * work_dim1;
#line 545 "ctprfb.f"
		i__5 = i__ + j * a_dim1;
#line 545 "ctprfb.f"
		z__1.r = work[i__4].r + a[i__5].r, z__1.i = work[i__4].i + a[
			i__5].i;
#line 545 "ctprfb.f"
		work[i__3].r = z__1.r, work[i__3].i = z__1.i;
#line 546 "ctprfb.f"
	    }
#line 547 "ctprfb.f"
	}

#line 549 "ctprfb.f"
	ctrmm_("R", "L", trans, "N", m, k, &c_b1, &t[t_offset], ldt, &work[
		work_offset], ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		ftnlen)1);

#line 552 "ctprfb.f"
	i__1 = *k;
#line 552 "ctprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 553 "ctprfb.f"
	    i__2 = *m;
#line 553 "ctprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 554 "ctprfb.f"
		i__3 = i__ + j * a_dim1;
#line 554 "ctprfb.f"
		i__4 = i__ + j * a_dim1;
#line 554 "ctprfb.f"
		i__5 = i__ + j * work_dim1;
#line 554 "ctprfb.f"
		z__1.r = a[i__4].r - work[i__5].r, z__1.i = a[i__4].i - work[
			i__5].i;
#line 554 "ctprfb.f"
		a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 555 "ctprfb.f"
	    }
#line 556 "ctprfb.f"
	}

#line 558 "ctprfb.f"
	i__1 = *n - *l;
#line 558 "ctprfb.f"
	z__1.r = -1., z__1.i = -0.;
#line 558 "ctprfb.f"
	cgemm_("N", "C", m, &i__1, k, &z__1, &work[work_offset], ldwork, &v[
		np + v_dim1], ldv, &c_b1, &b[np * b_dim1 + 1], ldb, (ftnlen)1,
		 (ftnlen)1);
#line 560 "ctprfb.f"
	i__1 = *k - *l;
#line 560 "ctprfb.f"
	z__1.r = -1., z__1.i = -0.;
#line 560 "ctprfb.f"
	cgemm_("N", "C", m, l, &i__1, &z__1, &work[work_offset], ldwork, &v[
		v_offset], ldv, &c_b1, &b[b_offset], ldb, (ftnlen)1, (ftnlen)
		1);
#line 562 "ctprfb.f"
	ctrmm_("R", "L", "C", "N", m, l, &c_b1, &v[kp * v_dim1 + 1], ldv, &
		work[kp * work_dim1 + 1], ldwork, (ftnlen)1, (ftnlen)1, (
		ftnlen)1, (ftnlen)1);
#line 564 "ctprfb.f"
	i__1 = *l;
#line 564 "ctprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 565 "ctprfb.f"
	    i__2 = *m;
#line 565 "ctprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 566 "ctprfb.f"
		i__3 = i__ + j * b_dim1;
#line 566 "ctprfb.f"
		i__4 = i__ + j * b_dim1;
#line 566 "ctprfb.f"
		i__5 = i__ + (*k - *l + j) * work_dim1;
#line 566 "ctprfb.f"
		z__1.r = b[i__4].r - work[i__5].r, z__1.i = b[i__4].i - work[
			i__5].i;
#line 566 "ctprfb.f"
		b[i__3].r = z__1.r, b[i__3].i = z__1.i;
#line 567 "ctprfb.f"
	    }
#line 568 "ctprfb.f"
	}

/* --------------------------------------------------------------------------- */

#line 572 "ctprfb.f"
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
#line 588 "ctprfb.f"
	i__1 = *m - *l + 1;
#line 588 "ctprfb.f"
	mp = min(i__1,*m);
/* Computing MIN */
#line 589 "ctprfb.f"
	i__1 = *l + 1;
#line 589 "ctprfb.f"
	kp = min(i__1,*k);

#line 591 "ctprfb.f"
	i__1 = *n;
#line 591 "ctprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 592 "ctprfb.f"
	    i__2 = *l;
#line 592 "ctprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 593 "ctprfb.f"
		i__3 = i__ + j * work_dim1;
#line 593 "ctprfb.f"
		i__4 = *m - *l + i__ + j * b_dim1;
#line 593 "ctprfb.f"
		work[i__3].r = b[i__4].r, work[i__3].i = b[i__4].i;
#line 594 "ctprfb.f"
	    }
#line 595 "ctprfb.f"
	}
#line 596 "ctprfb.f"
	ctrmm_("L", "L", "N", "N", l, n, &c_b1, &v[mp * v_dim1 + 1], ldv, &
		work[work_offset], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		ftnlen)1);
#line 598 "ctprfb.f"
	i__1 = *m - *l;
#line 598 "ctprfb.f"
	cgemm_("N", "N", l, n, &i__1, &c_b1, &v[v_offset], ldv, &b[b_offset], 
		ldb, &c_b1, &work[work_offset], ldwork, (ftnlen)1, (ftnlen)1);
#line 600 "ctprfb.f"
	i__1 = *k - *l;
#line 600 "ctprfb.f"
	cgemm_("N", "N", &i__1, n, m, &c_b1, &v[kp + v_dim1], ldv, &b[
		b_offset], ldb, &c_b2, &work[kp + work_dim1], ldwork, (ftnlen)
		1, (ftnlen)1);

#line 603 "ctprfb.f"
	i__1 = *n;
#line 603 "ctprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 604 "ctprfb.f"
	    i__2 = *k;
#line 604 "ctprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 605 "ctprfb.f"
		i__3 = i__ + j * work_dim1;
#line 605 "ctprfb.f"
		i__4 = i__ + j * work_dim1;
#line 605 "ctprfb.f"
		i__5 = i__ + j * a_dim1;
#line 605 "ctprfb.f"
		z__1.r = work[i__4].r + a[i__5].r, z__1.i = work[i__4].i + a[
			i__5].i;
#line 605 "ctprfb.f"
		work[i__3].r = z__1.r, work[i__3].i = z__1.i;
#line 606 "ctprfb.f"
	    }
#line 607 "ctprfb.f"
	}

#line 609 "ctprfb.f"
	ctrmm_("L", "U", trans, "N", k, n, &c_b1, &t[t_offset], ldt, &work[
		work_offset], ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		ftnlen)1);

#line 612 "ctprfb.f"
	i__1 = *n;
#line 612 "ctprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 613 "ctprfb.f"
	    i__2 = *k;
#line 613 "ctprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 614 "ctprfb.f"
		i__3 = i__ + j * a_dim1;
#line 614 "ctprfb.f"
		i__4 = i__ + j * a_dim1;
#line 614 "ctprfb.f"
		i__5 = i__ + j * work_dim1;
#line 614 "ctprfb.f"
		z__1.r = a[i__4].r - work[i__5].r, z__1.i = a[i__4].i - work[
			i__5].i;
#line 614 "ctprfb.f"
		a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 615 "ctprfb.f"
	    }
#line 616 "ctprfb.f"
	}

#line 618 "ctprfb.f"
	i__1 = *m - *l;
#line 618 "ctprfb.f"
	z__1.r = -1., z__1.i = -0.;
#line 618 "ctprfb.f"
	cgemm_("C", "N", &i__1, n, k, &z__1, &v[v_offset], ldv, &work[
		work_offset], ldwork, &c_b1, &b[b_offset], ldb, (ftnlen)1, (
		ftnlen)1);
#line 620 "ctprfb.f"
	i__1 = *k - *l;
#line 620 "ctprfb.f"
	z__1.r = -1., z__1.i = -0.;
#line 620 "ctprfb.f"
	cgemm_("C", "N", l, n, &i__1, &z__1, &v[kp + mp * v_dim1], ldv, &work[
		kp + work_dim1], ldwork, &c_b1, &b[mp + b_dim1], ldb, (ftnlen)
		1, (ftnlen)1);
#line 622 "ctprfb.f"
	ctrmm_("L", "L", "C", "N", l, n, &c_b1, &v[mp * v_dim1 + 1], ldv, &
		work[work_offset], ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		ftnlen)1);
#line 624 "ctprfb.f"
	i__1 = *n;
#line 624 "ctprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 625 "ctprfb.f"
	    i__2 = *l;
#line 625 "ctprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 626 "ctprfb.f"
		i__3 = *m - *l + i__ + j * b_dim1;
#line 626 "ctprfb.f"
		i__4 = *m - *l + i__ + j * b_dim1;
#line 626 "ctprfb.f"
		i__5 = i__ + j * work_dim1;
#line 626 "ctprfb.f"
		z__1.r = b[i__4].r - work[i__5].r, z__1.i = b[i__4].i - work[
			i__5].i;
#line 626 "ctprfb.f"
		b[i__3].r = z__1.r, b[i__3].i = z__1.i;
#line 627 "ctprfb.f"
	    }
#line 628 "ctprfb.f"
	}

/* --------------------------------------------------------------------------- */

#line 632 "ctprfb.f"
    } else if (row && forward && right) {

/* --------------------------------------------------------------------------- */

/*        Let  W =  [ I V ] ( I is K-by-K, V is K-by-N ) */

/*        Form  C H  or  C H**H  where  C = [ A B ] (A is M-by-K, B is M-by-N) */

/*        H = I - W**H T W            or  H**H = I - W**H T**H W */

/*        A = A - (A + B V**H) T      or  A = A - (A + B V**H) T**H */
/*        B = B - (A + B V**H) T V    or  B = B - (A + B V**H) T**H V */

/* --------------------------------------------------------------------------- */

/* Computing MIN */
#line 647 "ctprfb.f"
	i__1 = *n - *l + 1;
#line 647 "ctprfb.f"
	np = min(i__1,*n);
/* Computing MIN */
#line 648 "ctprfb.f"
	i__1 = *l + 1;
#line 648 "ctprfb.f"
	kp = min(i__1,*k);

#line 650 "ctprfb.f"
	i__1 = *l;
#line 650 "ctprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 651 "ctprfb.f"
	    i__2 = *m;
#line 651 "ctprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 652 "ctprfb.f"
		i__3 = i__ + j * work_dim1;
#line 652 "ctprfb.f"
		i__4 = i__ + (*n - *l + j) * b_dim1;
#line 652 "ctprfb.f"
		work[i__3].r = b[i__4].r, work[i__3].i = b[i__4].i;
#line 653 "ctprfb.f"
	    }
#line 654 "ctprfb.f"
	}
#line 655 "ctprfb.f"
	ctrmm_("R", "L", "C", "N", m, l, &c_b1, &v[np * v_dim1 + 1], ldv, &
		work[work_offset], ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		ftnlen)1);
#line 657 "ctprfb.f"
	i__1 = *n - *l;
#line 657 "ctprfb.f"
	cgemm_("N", "C", m, l, &i__1, &c_b1, &b[b_offset], ldb, &v[v_offset], 
		ldv, &c_b1, &work[work_offset], ldwork, (ftnlen)1, (ftnlen)1);
#line 659 "ctprfb.f"
	i__1 = *k - *l;
#line 659 "ctprfb.f"
	cgemm_("N", "C", m, &i__1, n, &c_b1, &b[b_offset], ldb, &v[kp + 
		v_dim1], ldv, &c_b2, &work[kp * work_dim1 + 1], ldwork, (
		ftnlen)1, (ftnlen)1);

#line 662 "ctprfb.f"
	i__1 = *k;
#line 662 "ctprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 663 "ctprfb.f"
	    i__2 = *m;
#line 663 "ctprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 664 "ctprfb.f"
		i__3 = i__ + j * work_dim1;
#line 664 "ctprfb.f"
		i__4 = i__ + j * work_dim1;
#line 664 "ctprfb.f"
		i__5 = i__ + j * a_dim1;
#line 664 "ctprfb.f"
		z__1.r = work[i__4].r + a[i__5].r, z__1.i = work[i__4].i + a[
			i__5].i;
#line 664 "ctprfb.f"
		work[i__3].r = z__1.r, work[i__3].i = z__1.i;
#line 665 "ctprfb.f"
	    }
#line 666 "ctprfb.f"
	}

#line 668 "ctprfb.f"
	ctrmm_("R", "U", trans, "N", m, k, &c_b1, &t[t_offset], ldt, &work[
		work_offset], ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		ftnlen)1);

#line 671 "ctprfb.f"
	i__1 = *k;
#line 671 "ctprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 672 "ctprfb.f"
	    i__2 = *m;
#line 672 "ctprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 673 "ctprfb.f"
		i__3 = i__ + j * a_dim1;
#line 673 "ctprfb.f"
		i__4 = i__ + j * a_dim1;
#line 673 "ctprfb.f"
		i__5 = i__ + j * work_dim1;
#line 673 "ctprfb.f"
		z__1.r = a[i__4].r - work[i__5].r, z__1.i = a[i__4].i - work[
			i__5].i;
#line 673 "ctprfb.f"
		a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 674 "ctprfb.f"
	    }
#line 675 "ctprfb.f"
	}

#line 677 "ctprfb.f"
	i__1 = *n - *l;
#line 677 "ctprfb.f"
	z__1.r = -1., z__1.i = -0.;
#line 677 "ctprfb.f"
	cgemm_("N", "N", m, &i__1, k, &z__1, &work[work_offset], ldwork, &v[
		v_offset], ldv, &c_b1, &b[b_offset], ldb, (ftnlen)1, (ftnlen)
		1);
#line 679 "ctprfb.f"
	i__1 = *k - *l;
#line 679 "ctprfb.f"
	z__1.r = -1., z__1.i = -0.;
#line 679 "ctprfb.f"
	cgemm_("N", "N", m, l, &i__1, &z__1, &work[kp * work_dim1 + 1], 
		ldwork, &v[kp + np * v_dim1], ldv, &c_b1, &b[np * b_dim1 + 1],
		 ldb, (ftnlen)1, (ftnlen)1);
#line 681 "ctprfb.f"
	ctrmm_("R", "L", "N", "N", m, l, &c_b1, &v[np * v_dim1 + 1], ldv, &
		work[work_offset], ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		ftnlen)1);
#line 683 "ctprfb.f"
	i__1 = *l;
#line 683 "ctprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 684 "ctprfb.f"
	    i__2 = *m;
#line 684 "ctprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 685 "ctprfb.f"
		i__3 = i__ + (*n - *l + j) * b_dim1;
#line 685 "ctprfb.f"
		i__4 = i__ + (*n - *l + j) * b_dim1;
#line 685 "ctprfb.f"
		i__5 = i__ + j * work_dim1;
#line 685 "ctprfb.f"
		z__1.r = b[i__4].r - work[i__5].r, z__1.i = b[i__4].i - work[
			i__5].i;
#line 685 "ctprfb.f"
		b[i__3].r = z__1.r, b[i__3].i = z__1.i;
#line 686 "ctprfb.f"
	    }
#line 687 "ctprfb.f"
	}

/* --------------------------------------------------------------------------- */

#line 691 "ctprfb.f"
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
#line 707 "ctprfb.f"
	i__1 = *l + 1;
#line 707 "ctprfb.f"
	mp = min(i__1,*m);
/* Computing MIN */
#line 708 "ctprfb.f"
	i__1 = *k - *l + 1;
#line 708 "ctprfb.f"
	kp = min(i__1,*k);

#line 710 "ctprfb.f"
	i__1 = *n;
#line 710 "ctprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 711 "ctprfb.f"
	    i__2 = *l;
#line 711 "ctprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 712 "ctprfb.f"
		i__3 = *k - *l + i__ + j * work_dim1;
#line 712 "ctprfb.f"
		i__4 = i__ + j * b_dim1;
#line 712 "ctprfb.f"
		work[i__3].r = b[i__4].r, work[i__3].i = b[i__4].i;
#line 713 "ctprfb.f"
	    }
#line 714 "ctprfb.f"
	}
#line 715 "ctprfb.f"
	ctrmm_("L", "U", "N", "N", l, n, &c_b1, &v[kp + v_dim1], ldv, &work[
		kp + work_dim1], ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		ftnlen)1);
#line 717 "ctprfb.f"
	i__1 = *m - *l;
#line 717 "ctprfb.f"
	cgemm_("N", "N", l, n, &i__1, &c_b1, &v[kp + mp * v_dim1], ldv, &b[mp 
		+ b_dim1], ldb, &c_b1, &work[kp + work_dim1], ldwork, (ftnlen)
		1, (ftnlen)1);
#line 719 "ctprfb.f"
	i__1 = *k - *l;
#line 719 "ctprfb.f"
	cgemm_("N", "N", &i__1, n, m, &c_b1, &v[v_offset], ldv, &b[b_offset], 
		ldb, &c_b2, &work[work_offset], ldwork, (ftnlen)1, (ftnlen)1);

#line 722 "ctprfb.f"
	i__1 = *n;
#line 722 "ctprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 723 "ctprfb.f"
	    i__2 = *k;
#line 723 "ctprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 724 "ctprfb.f"
		i__3 = i__ + j * work_dim1;
#line 724 "ctprfb.f"
		i__4 = i__ + j * work_dim1;
#line 724 "ctprfb.f"
		i__5 = i__ + j * a_dim1;
#line 724 "ctprfb.f"
		z__1.r = work[i__4].r + a[i__5].r, z__1.i = work[i__4].i + a[
			i__5].i;
#line 724 "ctprfb.f"
		work[i__3].r = z__1.r, work[i__3].i = z__1.i;
#line 725 "ctprfb.f"
	    }
#line 726 "ctprfb.f"
	}

#line 728 "ctprfb.f"
	ctrmm_("L", "L ", trans, "N", k, n, &c_b1, &t[t_offset], ldt, &work[
		work_offset], ldwork, (ftnlen)1, (ftnlen)2, (ftnlen)1, (
		ftnlen)1);

#line 731 "ctprfb.f"
	i__1 = *n;
#line 731 "ctprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 732 "ctprfb.f"
	    i__2 = *k;
#line 732 "ctprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 733 "ctprfb.f"
		i__3 = i__ + j * a_dim1;
#line 733 "ctprfb.f"
		i__4 = i__ + j * a_dim1;
#line 733 "ctprfb.f"
		i__5 = i__ + j * work_dim1;
#line 733 "ctprfb.f"
		z__1.r = a[i__4].r - work[i__5].r, z__1.i = a[i__4].i - work[
			i__5].i;
#line 733 "ctprfb.f"
		a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 734 "ctprfb.f"
	    }
#line 735 "ctprfb.f"
	}

#line 737 "ctprfb.f"
	i__1 = *m - *l;
#line 737 "ctprfb.f"
	z__1.r = -1., z__1.i = -0.;
#line 737 "ctprfb.f"
	cgemm_("C", "N", &i__1, n, k, &z__1, &v[mp * v_dim1 + 1], ldv, &work[
		work_offset], ldwork, &c_b1, &b[mp + b_dim1], ldb, (ftnlen)1, 
		(ftnlen)1);
#line 739 "ctprfb.f"
	i__1 = *k - *l;
#line 739 "ctprfb.f"
	z__1.r = -1., z__1.i = -0.;
#line 739 "ctprfb.f"
	cgemm_("C", "N", l, n, &i__1, &z__1, &v[v_offset], ldv, &work[
		work_offset], ldwork, &c_b1, &b[b_offset], ldb, (ftnlen)1, (
		ftnlen)1);
#line 741 "ctprfb.f"
	ctrmm_("L", "U", "C", "N", l, n, &c_b1, &v[kp + v_dim1], ldv, &work[
		kp + work_dim1], ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		ftnlen)1);
#line 743 "ctprfb.f"
	i__1 = *n;
#line 743 "ctprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 744 "ctprfb.f"
	    i__2 = *l;
#line 744 "ctprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 745 "ctprfb.f"
		i__3 = i__ + j * b_dim1;
#line 745 "ctprfb.f"
		i__4 = i__ + j * b_dim1;
#line 745 "ctprfb.f"
		i__5 = *k - *l + i__ + j * work_dim1;
#line 745 "ctprfb.f"
		z__1.r = b[i__4].r - work[i__5].r, z__1.i = b[i__4].i - work[
			i__5].i;
#line 745 "ctprfb.f"
		b[i__3].r = z__1.r, b[i__3].i = z__1.i;
#line 746 "ctprfb.f"
	    }
#line 747 "ctprfb.f"
	}

/* --------------------------------------------------------------------------- */

#line 751 "ctprfb.f"
    } else if (row && backward && right) {

/* --------------------------------------------------------------------------- */

/*        Let  W =  [ V I ] ( I is K-by-K, V is K-by-N ) */

/*        Form  C H  or  C H**H  where  C = [ B A ] (A is M-by-K, B is M-by-N) */

/*        H = I - W**H T W            or  H**H = I - W**H T**H W */

/*        A = A - (A + B V**H) T      or  A = A - (A + B V**H) T**H */
/*        B = B - (A + B V**H) T V    or  B = B - (A + B V**H) T**H V */

/* --------------------------------------------------------------------------- */

/* Computing MIN */
#line 766 "ctprfb.f"
	i__1 = *l + 1;
#line 766 "ctprfb.f"
	np = min(i__1,*n);
/* Computing MIN */
#line 767 "ctprfb.f"
	i__1 = *k - *l + 1;
#line 767 "ctprfb.f"
	kp = min(i__1,*k);

#line 769 "ctprfb.f"
	i__1 = *l;
#line 769 "ctprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 770 "ctprfb.f"
	    i__2 = *m;
#line 770 "ctprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 771 "ctprfb.f"
		i__3 = i__ + (*k - *l + j) * work_dim1;
#line 771 "ctprfb.f"
		i__4 = i__ + j * b_dim1;
#line 771 "ctprfb.f"
		work[i__3].r = b[i__4].r, work[i__3].i = b[i__4].i;
#line 772 "ctprfb.f"
	    }
#line 773 "ctprfb.f"
	}
#line 774 "ctprfb.f"
	ctrmm_("R", "U", "C", "N", m, l, &c_b1, &v[kp + v_dim1], ldv, &work[
		kp * work_dim1 + 1], ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, 
		(ftnlen)1);
#line 776 "ctprfb.f"
	i__1 = *n - *l;
#line 776 "ctprfb.f"
	cgemm_("N", "C", m, l, &i__1, &c_b1, &b[np * b_dim1 + 1], ldb, &v[kp 
		+ np * v_dim1], ldv, &c_b1, &work[kp * work_dim1 + 1], ldwork,
		 (ftnlen)1, (ftnlen)1);
#line 778 "ctprfb.f"
	i__1 = *k - *l;
#line 778 "ctprfb.f"
	cgemm_("N", "C", m, &i__1, n, &c_b1, &b[b_offset], ldb, &v[v_offset], 
		ldv, &c_b2, &work[work_offset], ldwork, (ftnlen)1, (ftnlen)1);

#line 781 "ctprfb.f"
	i__1 = *k;
#line 781 "ctprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 782 "ctprfb.f"
	    i__2 = *m;
#line 782 "ctprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 783 "ctprfb.f"
		i__3 = i__ + j * work_dim1;
#line 783 "ctprfb.f"
		i__4 = i__ + j * work_dim1;
#line 783 "ctprfb.f"
		i__5 = i__ + j * a_dim1;
#line 783 "ctprfb.f"
		z__1.r = work[i__4].r + a[i__5].r, z__1.i = work[i__4].i + a[
			i__5].i;
#line 783 "ctprfb.f"
		work[i__3].r = z__1.r, work[i__3].i = z__1.i;
#line 784 "ctprfb.f"
	    }
#line 785 "ctprfb.f"
	}

#line 787 "ctprfb.f"
	ctrmm_("R", "L", trans, "N", m, k, &c_b1, &t[t_offset], ldt, &work[
		work_offset], ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		ftnlen)1);

#line 790 "ctprfb.f"
	i__1 = *k;
#line 790 "ctprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 791 "ctprfb.f"
	    i__2 = *m;
#line 791 "ctprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 792 "ctprfb.f"
		i__3 = i__ + j * a_dim1;
#line 792 "ctprfb.f"
		i__4 = i__ + j * a_dim1;
#line 792 "ctprfb.f"
		i__5 = i__ + j * work_dim1;
#line 792 "ctprfb.f"
		z__1.r = a[i__4].r - work[i__5].r, z__1.i = a[i__4].i - work[
			i__5].i;
#line 792 "ctprfb.f"
		a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 793 "ctprfb.f"
	    }
#line 794 "ctprfb.f"
	}

#line 796 "ctprfb.f"
	i__1 = *n - *l;
#line 796 "ctprfb.f"
	z__1.r = -1., z__1.i = -0.;
#line 796 "ctprfb.f"
	cgemm_("N", "N", m, &i__1, k, &z__1, &work[work_offset], ldwork, &v[
		np * v_dim1 + 1], ldv, &c_b1, &b[np * b_dim1 + 1], ldb, (
		ftnlen)1, (ftnlen)1);
#line 798 "ctprfb.f"
	i__1 = *k - *l;
#line 798 "ctprfb.f"
	z__1.r = -1., z__1.i = -0.;
#line 798 "ctprfb.f"
	cgemm_("N", "N", m, l, &i__1, &z__1, &work[work_offset], ldwork, &v[
		v_offset], ldv, &c_b1, &b[b_offset], ldb, (ftnlen)1, (ftnlen)
		1);
#line 800 "ctprfb.f"
	ctrmm_("R", "U", "N", "N", m, l, &c_b1, &v[kp + v_dim1], ldv, &work[
		kp * work_dim1 + 1], ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, 
		(ftnlen)1);
#line 802 "ctprfb.f"
	i__1 = *l;
#line 802 "ctprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 803 "ctprfb.f"
	    i__2 = *m;
#line 803 "ctprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 804 "ctprfb.f"
		i__3 = i__ + j * b_dim1;
#line 804 "ctprfb.f"
		i__4 = i__ + j * b_dim1;
#line 804 "ctprfb.f"
		i__5 = i__ + (*k - *l + j) * work_dim1;
#line 804 "ctprfb.f"
		z__1.r = b[i__4].r - work[i__5].r, z__1.i = b[i__4].i - work[
			i__5].i;
#line 804 "ctprfb.f"
		b[i__3].r = z__1.r, b[i__3].i = z__1.i;
#line 805 "ctprfb.f"
	    }
#line 806 "ctprfb.f"
	}

#line 808 "ctprfb.f"
    }

#line 810 "ctprfb.f"
    return 0;

/*     End of CTPRFB */

} /* ctprfb_ */

