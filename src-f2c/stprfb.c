#line 1 "stprfb.f"
/* stprfb.f -- translated by f2c (version 20100827).
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

#line 1 "stprfb.f"
/* Table of constant values */

static doublereal c_b12 = 1.;
static doublereal c_b20 = 0.;
static doublereal c_b27 = -1.;

/* > \brief \b STPRFB applies a real or complex "triangular-pentagonal" blocked reflector to a real or complex
 matrix, which is composed of two blocks. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download STPRFB + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/stprfb.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/stprfb.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/stprfb.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE STPRFB( SIDE, TRANS, DIRECT, STOREV, M, N, K, L, */
/*                          V, LDV, T, LDT, A, LDA, B, LDB, WORK, LDWORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER DIRECT, SIDE, STOREV, TRANS */
/*       INTEGER   K, L, LDA, LDB, LDT, LDV, LDWORK, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL   A( LDA, * ), B( LDB, * ), T( LDT, * ), */
/*      $          V( LDV, * ), WORK( LDWORK, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > STPRFB applies a real "triangular-pentagonal" block reflector H or its */
/* > conjugate transpose H^H to a real matrix C, which is composed of two */
/* > blocks A and B, either from the left or right. */
/* > */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] SIDE */
/* > \verbatim */
/* >          SIDE is CHARACTER*1 */
/* >          = 'L': apply H or H^H from the Left */
/* >          = 'R': apply H or H^H from the Right */
/* > \endverbatim */
/* > */
/* > \param[in] TRANS */
/* > \verbatim */
/* >          TRANS is CHARACTER*1 */
/* >          = 'N': apply H (No transpose) */
/* >          = 'C': apply H^H (Conjugate transpose) */
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
/* >          V is REAL array, dimension */
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
/* >          T is REAL array, dimension (LDT,K) */
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
/* >          A is REAL array, dimension */
/* >          (LDA,N) if SIDE = 'L' or (LDA,K) if SIDE = 'R' */
/* >          On entry, the K-by-N or M-by-K matrix A. */
/* >          On exit, A is overwritten by the corresponding block of */
/* >          H*C or H^H*C or C*H or C*H^H.  See Futher Details. */
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
/* >          B is REAL array, dimension (LDB,N) */
/* >          On entry, the M-by-N matrix B. */
/* >          On exit, B is overwritten by the corresponding block of */
/* >          H*C or H^H*C or C*H or C*H^H.  See Further Details. */
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
/* >          WORK is REAL array, dimension */
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

/* > \ingroup realOTHERauxiliary */

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
/* Subroutine */ int stprfb_(char *side, char *trans, char *direct, char *
	storev, integer *m, integer *n, integer *k, integer *l, doublereal *v,
	 integer *ldv, doublereal *t, integer *ldt, doublereal *a, integer *
	lda, doublereal *b, integer *ldb, doublereal *work, integer *ldwork, 
	ftnlen side_len, ftnlen trans_len, ftnlen direct_len, ftnlen 
	storev_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, t_dim1, t_offset, v_dim1, 
	    v_offset, work_dim1, work_offset, i__1, i__2;

    /* Local variables */
    static logical backward;
    static integer i__, j, kp, mp, np;
    static logical row, left;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int sgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    static logical right;
    extern /* Subroutine */ int strmm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
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
/*     .. Executable Statements .. */

/*     Quick return if possible */

#line 289 "stprfb.f"
    /* Parameter adjustments */
#line 289 "stprfb.f"
    v_dim1 = *ldv;
#line 289 "stprfb.f"
    v_offset = 1 + v_dim1;
#line 289 "stprfb.f"
    v -= v_offset;
#line 289 "stprfb.f"
    t_dim1 = *ldt;
#line 289 "stprfb.f"
    t_offset = 1 + t_dim1;
#line 289 "stprfb.f"
    t -= t_offset;
#line 289 "stprfb.f"
    a_dim1 = *lda;
#line 289 "stprfb.f"
    a_offset = 1 + a_dim1;
#line 289 "stprfb.f"
    a -= a_offset;
#line 289 "stprfb.f"
    b_dim1 = *ldb;
#line 289 "stprfb.f"
    b_offset = 1 + b_dim1;
#line 289 "stprfb.f"
    b -= b_offset;
#line 289 "stprfb.f"
    work_dim1 = *ldwork;
#line 289 "stprfb.f"
    work_offset = 1 + work_dim1;
#line 289 "stprfb.f"
    work -= work_offset;
#line 289 "stprfb.f"

#line 289 "stprfb.f"
    /* Function Body */
#line 289 "stprfb.f"
    if (*m <= 0 || *n <= 0 || *k <= 0 || *l < 0) {
#line 289 "stprfb.f"
	return 0;
#line 289 "stprfb.f"
    }

#line 291 "stprfb.f"
    if (lsame_(storev, "C", (ftnlen)1, (ftnlen)1)) {
#line 292 "stprfb.f"
	column = TRUE_;
#line 293 "stprfb.f"
	row = FALSE_;
#line 294 "stprfb.f"
    } else if (lsame_(storev, "R", (ftnlen)1, (ftnlen)1)) {
#line 295 "stprfb.f"
	column = FALSE_;
#line 296 "stprfb.f"
	row = TRUE_;
#line 297 "stprfb.f"
    } else {
#line 298 "stprfb.f"
	column = FALSE_;
#line 299 "stprfb.f"
	row = FALSE_;
#line 300 "stprfb.f"
    }

#line 302 "stprfb.f"
    if (lsame_(side, "L", (ftnlen)1, (ftnlen)1)) {
#line 303 "stprfb.f"
	left = TRUE_;
#line 304 "stprfb.f"
	right = FALSE_;
#line 305 "stprfb.f"
    } else if (lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {
#line 306 "stprfb.f"
	left = FALSE_;
#line 307 "stprfb.f"
	right = TRUE_;
#line 308 "stprfb.f"
    } else {
#line 309 "stprfb.f"
	left = FALSE_;
#line 310 "stprfb.f"
	right = FALSE_;
#line 311 "stprfb.f"
    }

#line 313 "stprfb.f"
    if (lsame_(direct, "F", (ftnlen)1, (ftnlen)1)) {
#line 314 "stprfb.f"
	forward = TRUE_;
#line 315 "stprfb.f"
	backward = FALSE_;
#line 316 "stprfb.f"
    } else if (lsame_(direct, "B", (ftnlen)1, (ftnlen)1)) {
#line 317 "stprfb.f"
	forward = FALSE_;
#line 318 "stprfb.f"
	backward = TRUE_;
#line 319 "stprfb.f"
    } else {
#line 320 "stprfb.f"
	forward = FALSE_;
#line 321 "stprfb.f"
	backward = FALSE_;
#line 322 "stprfb.f"
    }

/* --------------------------------------------------------------------------- */

#line 326 "stprfb.f"
    if (column && forward && left) {

/* --------------------------------------------------------------------------- */

/*        Let  W =  [ I ]    (K-by-K) */
/*                  [ V ]    (M-by-K) */

/*        Form  H C  or  H^H C  where  C = [ A ]  (K-by-N) */
/*                                         [ B ]  (M-by-N) */

/*        H = I - W T W^H          or  H^H = I - W T^H W^H */

/*        A = A -   T (A + V^H B)  or  A = A -   T^H (A + V^H B) */
/*        B = B - V T (A + V^H B)  or  B = B - V T^H (A + V^H B) */

/* --------------------------------------------------------------------------- */

/* Computing MIN */
#line 343 "stprfb.f"
	i__1 = *m - *l + 1;
#line 343 "stprfb.f"
	mp = min(i__1,*m);
/* Computing MIN */
#line 344 "stprfb.f"
	i__1 = *l + 1;
#line 344 "stprfb.f"
	kp = min(i__1,*k);

#line 346 "stprfb.f"
	i__1 = *n;
#line 346 "stprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 347 "stprfb.f"
	    i__2 = *l;
#line 347 "stprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 348 "stprfb.f"
		work[i__ + j * work_dim1] = b[*m - *l + i__ + j * b_dim1];
#line 349 "stprfb.f"
	    }
#line 350 "stprfb.f"
	}
#line 351 "stprfb.f"
	strmm_("L", "U", "T", "N", l, n, &c_b12, &v[mp + v_dim1], ldv, &work[
		work_offset], ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		ftnlen)1);
#line 353 "stprfb.f"
	i__1 = *m - *l;
#line 353 "stprfb.f"
	sgemm_("T", "N", l, n, &i__1, &c_b12, &v[v_offset], ldv, &b[b_offset],
		 ldb, &c_b12, &work[work_offset], ldwork, (ftnlen)1, (ftnlen)
		1);
#line 355 "stprfb.f"
	i__1 = *k - *l;
#line 355 "stprfb.f"
	sgemm_("T", "N", &i__1, n, m, &c_b12, &v[kp * v_dim1 + 1], ldv, &b[
		b_offset], ldb, &c_b20, &work[kp + work_dim1], ldwork, (
		ftnlen)1, (ftnlen)1);

#line 358 "stprfb.f"
	i__1 = *n;
#line 358 "stprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 359 "stprfb.f"
	    i__2 = *k;
#line 359 "stprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 360 "stprfb.f"
		work[i__ + j * work_dim1] += a[i__ + j * a_dim1];
#line 361 "stprfb.f"
	    }
#line 362 "stprfb.f"
	}

#line 364 "stprfb.f"
	strmm_("L", "U", trans, "N", k, n, &c_b12, &t[t_offset], ldt, &work[
		work_offset], ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		ftnlen)1);

#line 367 "stprfb.f"
	i__1 = *n;
#line 367 "stprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 368 "stprfb.f"
	    i__2 = *k;
#line 368 "stprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 369 "stprfb.f"
		a[i__ + j * a_dim1] -= work[i__ + j * work_dim1];
#line 370 "stprfb.f"
	    }
#line 371 "stprfb.f"
	}

#line 373 "stprfb.f"
	i__1 = *m - *l;
#line 373 "stprfb.f"
	sgemm_("N", "N", &i__1, n, k, &c_b27, &v[v_offset], ldv, &work[
		work_offset], ldwork, &c_b12, &b[b_offset], ldb, (ftnlen)1, (
		ftnlen)1);
#line 375 "stprfb.f"
	i__1 = *k - *l;
#line 375 "stprfb.f"
	sgemm_("N", "N", l, n, &i__1, &c_b27, &v[mp + kp * v_dim1], ldv, &
		work[kp + work_dim1], ldwork, &c_b12, &b[mp + b_dim1], ldb, (
		ftnlen)1, (ftnlen)1);
#line 377 "stprfb.f"
	strmm_("L", "U", "N", "N", l, n, &c_b12, &v[mp + v_dim1], ldv, &work[
		work_offset], ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		ftnlen)1);
#line 379 "stprfb.f"
	i__1 = *n;
#line 379 "stprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 380 "stprfb.f"
	    i__2 = *l;
#line 380 "stprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 381 "stprfb.f"
		b[*m - *l + i__ + j * b_dim1] -= work[i__ + j * work_dim1];
#line 382 "stprfb.f"
	    }
#line 383 "stprfb.f"
	}

/* --------------------------------------------------------------------------- */

#line 387 "stprfb.f"
    } else if (column && forward && right) {

/* --------------------------------------------------------------------------- */

/*        Let  W =  [ I ]    (K-by-K) */
/*                  [ V ]    (N-by-K) */

/*        Form  C H or  C H^H  where  C = [ A B ] (A is M-by-K, B is M-by-N) */

/*        H = I - W T W^H          or  H^H = I - W T^H W^H */

/*        A = A - (A + B V) T      or  A = A - (A + B V) T^H */
/*        B = B - (A + B V) T V^H  or  B = B - (A + B V) T^H V^H */

/* --------------------------------------------------------------------------- */

/* Computing MIN */
#line 403 "stprfb.f"
	i__1 = *n - *l + 1;
#line 403 "stprfb.f"
	np = min(i__1,*n);
/* Computing MIN */
#line 404 "stprfb.f"
	i__1 = *l + 1;
#line 404 "stprfb.f"
	kp = min(i__1,*k);

#line 406 "stprfb.f"
	i__1 = *l;
#line 406 "stprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 407 "stprfb.f"
	    i__2 = *m;
#line 407 "stprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 408 "stprfb.f"
		work[i__ + j * work_dim1] = b[i__ + (*n - *l + j) * b_dim1];
#line 409 "stprfb.f"
	    }
#line 410 "stprfb.f"
	}
#line 411 "stprfb.f"
	strmm_("R", "U", "N", "N", m, l, &c_b12, &v[np + v_dim1], ldv, &work[
		work_offset], ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		ftnlen)1);
#line 413 "stprfb.f"
	i__1 = *n - *l;
#line 413 "stprfb.f"
	sgemm_("N", "N", m, l, &i__1, &c_b12, &b[b_offset], ldb, &v[v_offset],
		 ldv, &c_b12, &work[work_offset], ldwork, (ftnlen)1, (ftnlen)
		1);
#line 415 "stprfb.f"
	i__1 = *k - *l;
#line 415 "stprfb.f"
	sgemm_("N", "N", m, &i__1, n, &c_b12, &b[b_offset], ldb, &v[kp * 
		v_dim1 + 1], ldv, &c_b20, &work[kp * work_dim1 + 1], ldwork, (
		ftnlen)1, (ftnlen)1);

#line 418 "stprfb.f"
	i__1 = *k;
#line 418 "stprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 419 "stprfb.f"
	    i__2 = *m;
#line 419 "stprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 420 "stprfb.f"
		work[i__ + j * work_dim1] += a[i__ + j * a_dim1];
#line 421 "stprfb.f"
	    }
#line 422 "stprfb.f"
	}

#line 424 "stprfb.f"
	strmm_("R", "U", trans, "N", m, k, &c_b12, &t[t_offset], ldt, &work[
		work_offset], ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		ftnlen)1);

#line 427 "stprfb.f"
	i__1 = *k;
#line 427 "stprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 428 "stprfb.f"
	    i__2 = *m;
#line 428 "stprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 429 "stprfb.f"
		a[i__ + j * a_dim1] -= work[i__ + j * work_dim1];
#line 430 "stprfb.f"
	    }
#line 431 "stprfb.f"
	}

#line 433 "stprfb.f"
	i__1 = *n - *l;
#line 433 "stprfb.f"
	sgemm_("N", "T", m, &i__1, k, &c_b27, &work[work_offset], ldwork, &v[
		v_offset], ldv, &c_b12, &b[b_offset], ldb, (ftnlen)1, (ftnlen)
		1);
#line 435 "stprfb.f"
	i__1 = *k - *l;
#line 435 "stprfb.f"
	sgemm_("N", "T", m, l, &i__1, &c_b27, &work[kp * work_dim1 + 1], 
		ldwork, &v[np + kp * v_dim1], ldv, &c_b12, &b[np * b_dim1 + 1]
		, ldb, (ftnlen)1, (ftnlen)1);
#line 437 "stprfb.f"
	strmm_("R", "U", "T", "N", m, l, &c_b12, &v[np + v_dim1], ldv, &work[
		work_offset], ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		ftnlen)1);
#line 439 "stprfb.f"
	i__1 = *l;
#line 439 "stprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 440 "stprfb.f"
	    i__2 = *m;
#line 440 "stprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 441 "stprfb.f"
		b[i__ + (*n - *l + j) * b_dim1] -= work[i__ + j * work_dim1];
#line 442 "stprfb.f"
	    }
#line 443 "stprfb.f"
	}

/* --------------------------------------------------------------------------- */

#line 447 "stprfb.f"
    } else if (column && backward && left) {

/* --------------------------------------------------------------------------- */

/*        Let  W =  [ V ]    (M-by-K) */
/*                  [ I ]    (K-by-K) */

/*        Form  H C  or  H^H C  where  C = [ B ]  (M-by-N) */
/*                                         [ A ]  (K-by-N) */

/*        H = I - W T W^H          or  H^H = I - W T^H W^H */

/*        A = A -   T (A + V^H B)  or  A = A -   T^H (A + V^H B) */
/*        B = B - V T (A + V^H B)  or  B = B - V T^H (A + V^H B) */

/* --------------------------------------------------------------------------- */

/* Computing MIN */
#line 464 "stprfb.f"
	i__1 = *l + 1;
#line 464 "stprfb.f"
	mp = min(i__1,*m);
/* Computing MIN */
#line 465 "stprfb.f"
	i__1 = *k - *l + 1;
#line 465 "stprfb.f"
	kp = min(i__1,*k);

#line 467 "stprfb.f"
	i__1 = *n;
#line 467 "stprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 468 "stprfb.f"
	    i__2 = *l;
#line 468 "stprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 469 "stprfb.f"
		work[*k - *l + i__ + j * work_dim1] = b[i__ + j * b_dim1];
#line 470 "stprfb.f"
	    }
#line 471 "stprfb.f"
	}

#line 473 "stprfb.f"
	strmm_("L", "L", "T", "N", l, n, &c_b12, &v[kp * v_dim1 + 1], ldv, &
		work[kp + work_dim1], ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1,
		 (ftnlen)1);
#line 475 "stprfb.f"
	i__1 = *m - *l;
#line 475 "stprfb.f"
	sgemm_("T", "N", l, n, &i__1, &c_b12, &v[mp + kp * v_dim1], ldv, &b[
		mp + b_dim1], ldb, &c_b12, &work[kp + work_dim1], ldwork, (
		ftnlen)1, (ftnlen)1);
#line 477 "stprfb.f"
	i__1 = *k - *l;
#line 477 "stprfb.f"
	sgemm_("T", "N", &i__1, n, m, &c_b12, &v[v_offset], ldv, &b[b_offset],
		 ldb, &c_b20, &work[work_offset], ldwork, (ftnlen)1, (ftnlen)
		1);

#line 480 "stprfb.f"
	i__1 = *n;
#line 480 "stprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 481 "stprfb.f"
	    i__2 = *k;
#line 481 "stprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 482 "stprfb.f"
		work[i__ + j * work_dim1] += a[i__ + j * a_dim1];
#line 483 "stprfb.f"
	    }
#line 484 "stprfb.f"
	}

#line 486 "stprfb.f"
	strmm_("L", "L", trans, "N", k, n, &c_b12, &t[t_offset], ldt, &work[
		work_offset], ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		ftnlen)1);

#line 489 "stprfb.f"
	i__1 = *n;
#line 489 "stprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 490 "stprfb.f"
	    i__2 = *k;
#line 490 "stprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 491 "stprfb.f"
		a[i__ + j * a_dim1] -= work[i__ + j * work_dim1];
#line 492 "stprfb.f"
	    }
#line 493 "stprfb.f"
	}

#line 495 "stprfb.f"
	i__1 = *m - *l;
#line 495 "stprfb.f"
	sgemm_("N", "N", &i__1, n, k, &c_b27, &v[mp + v_dim1], ldv, &work[
		work_offset], ldwork, &c_b12, &b[mp + b_dim1], ldb, (ftnlen)1,
		 (ftnlen)1);
#line 497 "stprfb.f"
	i__1 = *k - *l;
#line 497 "stprfb.f"
	sgemm_("N", "N", l, n, &i__1, &c_b27, &v[v_offset], ldv, &work[
		work_offset], ldwork, &c_b12, &b[b_offset], ldb, (ftnlen)1, (
		ftnlen)1);
#line 499 "stprfb.f"
	strmm_("L", "L", "N", "N", l, n, &c_b12, &v[kp * v_dim1 + 1], ldv, &
		work[kp + work_dim1], ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1,
		 (ftnlen)1);
#line 501 "stprfb.f"
	i__1 = *n;
#line 501 "stprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 502 "stprfb.f"
	    i__2 = *l;
#line 502 "stprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 503 "stprfb.f"
		b[i__ + j * b_dim1] -= work[*k - *l + i__ + j * work_dim1];
#line 504 "stprfb.f"
	    }
#line 505 "stprfb.f"
	}

/* --------------------------------------------------------------------------- */

#line 509 "stprfb.f"
    } else if (column && backward && right) {

/* --------------------------------------------------------------------------- */

/*        Let  W =  [ V ]    (N-by-K) */
/*                  [ I ]    (K-by-K) */

/*        Form  C H  or  C H^H  where  C = [ B A ] (B is M-by-N, A is M-by-K) */

/*        H = I - W T W^H          or  H^H = I - W T^H W^H */

/*        A = A - (A + B V) T      or  A = A - (A + B V) T^H */
/*        B = B - (A + B V) T V^H  or  B = B - (A + B V) T^H V^H */

/* --------------------------------------------------------------------------- */

/* Computing MIN */
#line 525 "stprfb.f"
	i__1 = *l + 1;
#line 525 "stprfb.f"
	np = min(i__1,*n);
/* Computing MIN */
#line 526 "stprfb.f"
	i__1 = *k - *l + 1;
#line 526 "stprfb.f"
	kp = min(i__1,*k);

#line 528 "stprfb.f"
	i__1 = *l;
#line 528 "stprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 529 "stprfb.f"
	    i__2 = *m;
#line 529 "stprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 530 "stprfb.f"
		work[i__ + (*k - *l + j) * work_dim1] = b[i__ + j * b_dim1];
#line 531 "stprfb.f"
	    }
#line 532 "stprfb.f"
	}
#line 533 "stprfb.f"
	strmm_("R", "L", "N", "N", m, l, &c_b12, &v[kp * v_dim1 + 1], ldv, &
		work[kp * work_dim1 + 1], ldwork, (ftnlen)1, (ftnlen)1, (
		ftnlen)1, (ftnlen)1);
#line 535 "stprfb.f"
	i__1 = *n - *l;
#line 535 "stprfb.f"
	sgemm_("N", "N", m, l, &i__1, &c_b12, &b[np * b_dim1 + 1], ldb, &v[np 
		+ kp * v_dim1], ldv, &c_b12, &work[kp * work_dim1 + 1], 
		ldwork, (ftnlen)1, (ftnlen)1);
#line 537 "stprfb.f"
	i__1 = *k - *l;
#line 537 "stprfb.f"
	sgemm_("N", "N", m, &i__1, n, &c_b12, &b[b_offset], ldb, &v[v_offset],
		 ldv, &c_b20, &work[work_offset], ldwork, (ftnlen)1, (ftnlen)
		1);

#line 540 "stprfb.f"
	i__1 = *k;
#line 540 "stprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 541 "stprfb.f"
	    i__2 = *m;
#line 541 "stprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 542 "stprfb.f"
		work[i__ + j * work_dim1] += a[i__ + j * a_dim1];
#line 543 "stprfb.f"
	    }
#line 544 "stprfb.f"
	}

#line 546 "stprfb.f"
	strmm_("R", "L", trans, "N", m, k, &c_b12, &t[t_offset], ldt, &work[
		work_offset], ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		ftnlen)1);

#line 549 "stprfb.f"
	i__1 = *k;
#line 549 "stprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 550 "stprfb.f"
	    i__2 = *m;
#line 550 "stprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 551 "stprfb.f"
		a[i__ + j * a_dim1] -= work[i__ + j * work_dim1];
#line 552 "stprfb.f"
	    }
#line 553 "stprfb.f"
	}

#line 555 "stprfb.f"
	i__1 = *n - *l;
#line 555 "stprfb.f"
	sgemm_("N", "T", m, &i__1, k, &c_b27, &work[work_offset], ldwork, &v[
		np + v_dim1], ldv, &c_b12, &b[np * b_dim1 + 1], ldb, (ftnlen)
		1, (ftnlen)1);
#line 557 "stprfb.f"
	i__1 = *k - *l;
#line 557 "stprfb.f"
	sgemm_("N", "T", m, l, &i__1, &c_b27, &work[work_offset], ldwork, &v[
		v_offset], ldv, &c_b12, &b[b_offset], ldb, (ftnlen)1, (ftnlen)
		1);
#line 559 "stprfb.f"
	strmm_("R", "L", "T", "N", m, l, &c_b12, &v[kp * v_dim1 + 1], ldv, &
		work[kp * work_dim1 + 1], ldwork, (ftnlen)1, (ftnlen)1, (
		ftnlen)1, (ftnlen)1);
#line 561 "stprfb.f"
	i__1 = *l;
#line 561 "stprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 562 "stprfb.f"
	    i__2 = *m;
#line 562 "stprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 563 "stprfb.f"
		b[i__ + j * b_dim1] -= work[i__ + (*k - *l + j) * work_dim1];
#line 564 "stprfb.f"
	    }
#line 565 "stprfb.f"
	}

/* --------------------------------------------------------------------------- */

#line 569 "stprfb.f"
    } else if (row && forward && left) {

/* --------------------------------------------------------------------------- */

/*        Let  W =  [ I V ] ( I is K-by-K, V is K-by-M ) */

/*        Form  H C  or  H^H C  where  C = [ A ]  (K-by-N) */
/*                                         [ B ]  (M-by-N) */

/*        H = I - W^H T W          or  H^H = I - W^H T^H W */

/*        A = A -     T (A + V B)  or  A = A -     T^H (A + V B) */
/*        B = B - V^H T (A + V B)  or  B = B - V^H T^H (A + V B) */

/* --------------------------------------------------------------------------- */

/* Computing MIN */
#line 585 "stprfb.f"
	i__1 = *m - *l + 1;
#line 585 "stprfb.f"
	mp = min(i__1,*m);
/* Computing MIN */
#line 586 "stprfb.f"
	i__1 = *l + 1;
#line 586 "stprfb.f"
	kp = min(i__1,*k);

#line 588 "stprfb.f"
	i__1 = *n;
#line 588 "stprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 589 "stprfb.f"
	    i__2 = *l;
#line 589 "stprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 590 "stprfb.f"
		work[i__ + j * work_dim1] = b[*m - *l + i__ + j * b_dim1];
#line 591 "stprfb.f"
	    }
#line 592 "stprfb.f"
	}
#line 593 "stprfb.f"
	strmm_("L", "L", "N", "N", l, n, &c_b12, &v[mp * v_dim1 + 1], ldv, &
		work[work_offset], ldb, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		ftnlen)1);
#line 595 "stprfb.f"
	i__1 = *m - *l;
#line 595 "stprfb.f"
	sgemm_("N", "N", l, n, &i__1, &c_b12, &v[v_offset], ldv, &b[b_offset],
		 ldb, &c_b12, &work[work_offset], ldwork, (ftnlen)1, (ftnlen)
		1);
#line 597 "stprfb.f"
	i__1 = *k - *l;
#line 597 "stprfb.f"
	sgemm_("N", "N", &i__1, n, m, &c_b12, &v[kp + v_dim1], ldv, &b[
		b_offset], ldb, &c_b20, &work[kp + work_dim1], ldwork, (
		ftnlen)1, (ftnlen)1);

#line 600 "stprfb.f"
	i__1 = *n;
#line 600 "stprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 601 "stprfb.f"
	    i__2 = *k;
#line 601 "stprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 602 "stprfb.f"
		work[i__ + j * work_dim1] += a[i__ + j * a_dim1];
#line 603 "stprfb.f"
	    }
#line 604 "stprfb.f"
	}

#line 606 "stprfb.f"
	strmm_("L", "U", trans, "N", k, n, &c_b12, &t[t_offset], ldt, &work[
		work_offset], ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		ftnlen)1);

#line 609 "stprfb.f"
	i__1 = *n;
#line 609 "stprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 610 "stprfb.f"
	    i__2 = *k;
#line 610 "stprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 611 "stprfb.f"
		a[i__ + j * a_dim1] -= work[i__ + j * work_dim1];
#line 612 "stprfb.f"
	    }
#line 613 "stprfb.f"
	}

#line 615 "stprfb.f"
	i__1 = *m - *l;
#line 615 "stprfb.f"
	sgemm_("T", "N", &i__1, n, k, &c_b27, &v[v_offset], ldv, &work[
		work_offset], ldwork, &c_b12, &b[b_offset], ldb, (ftnlen)1, (
		ftnlen)1);
#line 617 "stprfb.f"
	i__1 = *k - *l;
#line 617 "stprfb.f"
	sgemm_("T", "N", l, n, &i__1, &c_b27, &v[kp + mp * v_dim1], ldv, &
		work[kp + work_dim1], ldwork, &c_b12, &b[mp + b_dim1], ldb, (
		ftnlen)1, (ftnlen)1);
#line 619 "stprfb.f"
	strmm_("L", "L", "T", "N", l, n, &c_b12, &v[mp * v_dim1 + 1], ldv, &
		work[work_offset], ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		ftnlen)1);
#line 621 "stprfb.f"
	i__1 = *n;
#line 621 "stprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 622 "stprfb.f"
	    i__2 = *l;
#line 622 "stprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 623 "stprfb.f"
		b[*m - *l + i__ + j * b_dim1] -= work[i__ + j * work_dim1];
#line 624 "stprfb.f"
	    }
#line 625 "stprfb.f"
	}

/* --------------------------------------------------------------------------- */

#line 629 "stprfb.f"
    } else if (row && forward && right) {

/* --------------------------------------------------------------------------- */

/*        Let  W =  [ I V ] ( I is K-by-K, V is K-by-N ) */

/*        Form  C H  or  C H^H  where  C = [ A B ] (A is M-by-K, B is M-by-N) */

/*        H = I - W^H T W            or  H^H = I - W^H T^H W */

/*        A = A - (A + B V^H) T      or  A = A - (A + B V^H) T^H */
/*        B = B - (A + B V^H) T V    or  B = B - (A + B V^H) T^H V */

/* --------------------------------------------------------------------------- */

/* Computing MIN */
#line 644 "stprfb.f"
	i__1 = *n - *l + 1;
#line 644 "stprfb.f"
	np = min(i__1,*n);
/* Computing MIN */
#line 645 "stprfb.f"
	i__1 = *l + 1;
#line 645 "stprfb.f"
	kp = min(i__1,*k);

#line 647 "stprfb.f"
	i__1 = *l;
#line 647 "stprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 648 "stprfb.f"
	    i__2 = *m;
#line 648 "stprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 649 "stprfb.f"
		work[i__ + j * work_dim1] = b[i__ + (*n - *l + j) * b_dim1];
#line 650 "stprfb.f"
	    }
#line 651 "stprfb.f"
	}
#line 652 "stprfb.f"
	strmm_("R", "L", "T", "N", m, l, &c_b12, &v[np * v_dim1 + 1], ldv, &
		work[work_offset], ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		ftnlen)1);
#line 654 "stprfb.f"
	i__1 = *n - *l;
#line 654 "stprfb.f"
	sgemm_("N", "T", m, l, &i__1, &c_b12, &b[b_offset], ldb, &v[v_offset],
		 ldv, &c_b12, &work[work_offset], ldwork, (ftnlen)1, (ftnlen)
		1);
#line 656 "stprfb.f"
	i__1 = *k - *l;
#line 656 "stprfb.f"
	sgemm_("N", "T", m, &i__1, n, &c_b12, &b[b_offset], ldb, &v[kp + 
		v_dim1], ldv, &c_b20, &work[kp * work_dim1 + 1], ldwork, (
		ftnlen)1, (ftnlen)1);

#line 659 "stprfb.f"
	i__1 = *k;
#line 659 "stprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 660 "stprfb.f"
	    i__2 = *m;
#line 660 "stprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 661 "stprfb.f"
		work[i__ + j * work_dim1] += a[i__ + j * a_dim1];
#line 662 "stprfb.f"
	    }
#line 663 "stprfb.f"
	}

#line 665 "stprfb.f"
	strmm_("R", "U", trans, "N", m, k, &c_b12, &t[t_offset], ldt, &work[
		work_offset], ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		ftnlen)1);

#line 668 "stprfb.f"
	i__1 = *k;
#line 668 "stprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 669 "stprfb.f"
	    i__2 = *m;
#line 669 "stprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 670 "stprfb.f"
		a[i__ + j * a_dim1] -= work[i__ + j * work_dim1];
#line 671 "stprfb.f"
	    }
#line 672 "stprfb.f"
	}

#line 674 "stprfb.f"
	i__1 = *n - *l;
#line 674 "stprfb.f"
	sgemm_("N", "N", m, &i__1, k, &c_b27, &work[work_offset], ldwork, &v[
		v_offset], ldv, &c_b12, &b[b_offset], ldb, (ftnlen)1, (ftnlen)
		1);
#line 676 "stprfb.f"
	i__1 = *k - *l;
#line 676 "stprfb.f"
	sgemm_("N", "N", m, l, &i__1, &c_b27, &work[kp * work_dim1 + 1], 
		ldwork, &v[kp + np * v_dim1], ldv, &c_b12, &b[np * b_dim1 + 1]
		, ldb, (ftnlen)1, (ftnlen)1);
#line 678 "stprfb.f"
	strmm_("R", "L", "N", "N", m, l, &c_b12, &v[np * v_dim1 + 1], ldv, &
		work[work_offset], ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		ftnlen)1);
#line 680 "stprfb.f"
	i__1 = *l;
#line 680 "stprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 681 "stprfb.f"
	    i__2 = *m;
#line 681 "stprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 682 "stprfb.f"
		b[i__ + (*n - *l + j) * b_dim1] -= work[i__ + j * work_dim1];
#line 683 "stprfb.f"
	    }
#line 684 "stprfb.f"
	}

/* --------------------------------------------------------------------------- */

#line 688 "stprfb.f"
    } else if (row && backward && left) {

/* --------------------------------------------------------------------------- */

/*        Let  W =  [ V I ] ( I is K-by-K, V is K-by-M ) */

/*        Form  H C  or  H^H C  where  C = [ B ]  (M-by-N) */
/*                                         [ A ]  (K-by-N) */

/*        H = I - W^H T W          or  H^H = I - W^H T^H W */

/*        A = A -     T (A + V B)  or  A = A -     T^H (A + V B) */
/*        B = B - V^H T (A + V B)  or  B = B - V^H T^H (A + V B) */

/* --------------------------------------------------------------------------- */

/* Computing MIN */
#line 704 "stprfb.f"
	i__1 = *l + 1;
#line 704 "stprfb.f"
	mp = min(i__1,*m);
/* Computing MIN */
#line 705 "stprfb.f"
	i__1 = *k - *l + 1;
#line 705 "stprfb.f"
	kp = min(i__1,*k);

#line 707 "stprfb.f"
	i__1 = *n;
#line 707 "stprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 708 "stprfb.f"
	    i__2 = *l;
#line 708 "stprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 709 "stprfb.f"
		work[*k - *l + i__ + j * work_dim1] = b[i__ + j * b_dim1];
#line 710 "stprfb.f"
	    }
#line 711 "stprfb.f"
	}
#line 712 "stprfb.f"
	strmm_("L", "U", "N", "N", l, n, &c_b12, &v[kp + v_dim1], ldv, &work[
		kp + work_dim1], ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		ftnlen)1);
#line 714 "stprfb.f"
	i__1 = *m - *l;
#line 714 "stprfb.f"
	sgemm_("N", "N", l, n, &i__1, &c_b12, &v[kp + mp * v_dim1], ldv, &b[
		mp + b_dim1], ldb, &c_b12, &work[kp + work_dim1], ldwork, (
		ftnlen)1, (ftnlen)1);
#line 716 "stprfb.f"
	i__1 = *k - *l;
#line 716 "stprfb.f"
	sgemm_("N", "N", &i__1, n, m, &c_b12, &v[v_offset], ldv, &b[b_offset],
		 ldb, &c_b20, &work[work_offset], ldwork, (ftnlen)1, (ftnlen)
		1);

#line 719 "stprfb.f"
	i__1 = *n;
#line 719 "stprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 720 "stprfb.f"
	    i__2 = *k;
#line 720 "stprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 721 "stprfb.f"
		work[i__ + j * work_dim1] += a[i__ + j * a_dim1];
#line 722 "stprfb.f"
	    }
#line 723 "stprfb.f"
	}

#line 725 "stprfb.f"
	strmm_("L", "L ", trans, "N", k, n, &c_b12, &t[t_offset], ldt, &work[
		work_offset], ldwork, (ftnlen)1, (ftnlen)2, (ftnlen)1, (
		ftnlen)1);

#line 728 "stprfb.f"
	i__1 = *n;
#line 728 "stprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 729 "stprfb.f"
	    i__2 = *k;
#line 729 "stprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 730 "stprfb.f"
		a[i__ + j * a_dim1] -= work[i__ + j * work_dim1];
#line 731 "stprfb.f"
	    }
#line 732 "stprfb.f"
	}

#line 734 "stprfb.f"
	i__1 = *m - *l;
#line 734 "stprfb.f"
	sgemm_("T", "N", &i__1, n, k, &c_b27, &v[mp * v_dim1 + 1], ldv, &work[
		work_offset], ldwork, &c_b12, &b[mp + b_dim1], ldb, (ftnlen)1,
		 (ftnlen)1);
#line 736 "stprfb.f"
	i__1 = *k - *l;
#line 736 "stprfb.f"
	sgemm_("T", "N", l, n, &i__1, &c_b27, &v[v_offset], ldv, &work[
		work_offset], ldwork, &c_b12, &b[b_offset], ldb, (ftnlen)1, (
		ftnlen)1);
#line 738 "stprfb.f"
	strmm_("L", "U", "T", "N", l, n, &c_b12, &v[kp + v_dim1], ldv, &work[
		kp + work_dim1], ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		ftnlen)1);
#line 740 "stprfb.f"
	i__1 = *n;
#line 740 "stprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 741 "stprfb.f"
	    i__2 = *l;
#line 741 "stprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 742 "stprfb.f"
		b[i__ + j * b_dim1] -= work[*k - *l + i__ + j * work_dim1];
#line 743 "stprfb.f"
	    }
#line 744 "stprfb.f"
	}

/* --------------------------------------------------------------------------- */

#line 748 "stprfb.f"
    } else if (row && backward && right) {

/* --------------------------------------------------------------------------- */

/*        Let  W =  [ V I ] ( I is K-by-K, V is K-by-N ) */

/*        Form  C H  or  C H^H  where  C = [ B A ] (A is M-by-K, B is M-by-N) */

/*        H = I - W^H T W            or  H^H = I - W^H T^H W */

/*        A = A - (A + B V^H) T      or  A = A - (A + B V^H) T^H */
/*        B = B - (A + B V^H) T V    or  B = B - (A + B V^H) T^H V */

/* --------------------------------------------------------------------------- */

/* Computing MIN */
#line 763 "stprfb.f"
	i__1 = *l + 1;
#line 763 "stprfb.f"
	np = min(i__1,*n);
/* Computing MIN */
#line 764 "stprfb.f"
	i__1 = *k - *l + 1;
#line 764 "stprfb.f"
	kp = min(i__1,*k);

#line 766 "stprfb.f"
	i__1 = *l;
#line 766 "stprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 767 "stprfb.f"
	    i__2 = *m;
#line 767 "stprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 768 "stprfb.f"
		work[i__ + (*k - *l + j) * work_dim1] = b[i__ + j * b_dim1];
#line 769 "stprfb.f"
	    }
#line 770 "stprfb.f"
	}
#line 771 "stprfb.f"
	strmm_("R", "U", "T", "N", m, l, &c_b12, &v[kp + v_dim1], ldv, &work[
		kp * work_dim1 + 1], ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, 
		(ftnlen)1);
#line 773 "stprfb.f"
	i__1 = *n - *l;
#line 773 "stprfb.f"
	sgemm_("N", "T", m, l, &i__1, &c_b12, &b[np * b_dim1 + 1], ldb, &v[kp 
		+ np * v_dim1], ldv, &c_b12, &work[kp * work_dim1 + 1], 
		ldwork, (ftnlen)1, (ftnlen)1);
#line 775 "stprfb.f"
	i__1 = *k - *l;
#line 775 "stprfb.f"
	sgemm_("N", "T", m, &i__1, n, &c_b12, &b[b_offset], ldb, &v[v_offset],
		 ldv, &c_b20, &work[work_offset], ldwork, (ftnlen)1, (ftnlen)
		1);

#line 778 "stprfb.f"
	i__1 = *k;
#line 778 "stprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 779 "stprfb.f"
	    i__2 = *m;
#line 779 "stprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 780 "stprfb.f"
		work[i__ + j * work_dim1] += a[i__ + j * a_dim1];
#line 781 "stprfb.f"
	    }
#line 782 "stprfb.f"
	}

#line 784 "stprfb.f"
	strmm_("R", "L", trans, "N", m, k, &c_b12, &t[t_offset], ldt, &work[
		work_offset], ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, (
		ftnlen)1);

#line 787 "stprfb.f"
	i__1 = *k;
#line 787 "stprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 788 "stprfb.f"
	    i__2 = *m;
#line 788 "stprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 789 "stprfb.f"
		a[i__ + j * a_dim1] -= work[i__ + j * work_dim1];
#line 790 "stprfb.f"
	    }
#line 791 "stprfb.f"
	}

#line 793 "stprfb.f"
	i__1 = *n - *l;
#line 793 "stprfb.f"
	sgemm_("N", "N", m, &i__1, k, &c_b27, &work[work_offset], ldwork, &v[
		np * v_dim1 + 1], ldv, &c_b12, &b[np * b_dim1 + 1], ldb, (
		ftnlen)1, (ftnlen)1);
#line 795 "stprfb.f"
	i__1 = *k - *l;
#line 795 "stprfb.f"
	sgemm_("N", "N", m, l, &i__1, &c_b27, &work[work_offset], ldwork, &v[
		v_offset], ldv, &c_b12, &b[b_offset], ldb, (ftnlen)1, (ftnlen)
		1);
#line 797 "stprfb.f"
	strmm_("R", "U", "N", "N", m, l, &c_b12, &v[kp + v_dim1], ldv, &work[
		kp * work_dim1 + 1], ldwork, (ftnlen)1, (ftnlen)1, (ftnlen)1, 
		(ftnlen)1);
#line 799 "stprfb.f"
	i__1 = *l;
#line 799 "stprfb.f"
	for (j = 1; j <= i__1; ++j) {
#line 800 "stprfb.f"
	    i__2 = *m;
#line 800 "stprfb.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 801 "stprfb.f"
		b[i__ + j * b_dim1] -= work[i__ + (*k - *l + j) * work_dim1];
#line 802 "stprfb.f"
	    }
#line 803 "stprfb.f"
	}

#line 805 "stprfb.f"
    }

#line 807 "stprfb.f"
    return 0;

/*     End of STPRFB */

} /* stprfb_ */

