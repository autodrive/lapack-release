#line 1 "dlarfb.f"
/* dlarfb.f -- translated by f2c (version 20100827).
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

#line 1 "dlarfb.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b14 = 1.;
static doublereal c_b25 = -1.;

/* > \brief \b DLARFB applies a block reflector or its transpose to a general rectangular matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLARFB + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlarfb.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlarfb.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlarfb.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLARFB( SIDE, TRANS, DIRECT, STOREV, M, N, K, V, LDV, */
/*                          T, LDT, C, LDC, WORK, LDWORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          DIRECT, SIDE, STOREV, TRANS */
/*       INTEGER            K, LDC, LDT, LDV, LDWORK, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   C( LDC, * ), T( LDT, * ), V( LDV, * ), */
/*      $                   WORK( LDWORK, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLARFB applies a real block reflector H or its transpose H**T to a */
/* > real m by n matrix C, from either the left or the right. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] SIDE */
/* > \verbatim */
/* >          SIDE is CHARACTER*1 */
/* >          = 'L': apply H or H**T from the Left */
/* >          = 'R': apply H or H**T from the Right */
/* > \endverbatim */
/* > */
/* > \param[in] TRANS */
/* > \verbatim */
/* >          TRANS is CHARACTER*1 */
/* >          = 'N': apply H (No transpose) */
/* >          = 'T': apply H**T (Transpose) */
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
/* >          = 'C': Columnwise */
/* >          = 'R': Rowwise */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the matrix C. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix C. */
/* > \endverbatim */
/* > */
/* > \param[in] K */
/* > \verbatim */
/* >          K is INTEGER */
/* >          The order of the matrix T (= the number of elementary */
/* >          reflectors whose product defines the block reflector). */
/* > \endverbatim */
/* > */
/* > \param[in] V */
/* > \verbatim */
/* >          V is DOUBLE PRECISION array, dimension */
/* >                                (LDV,K) if STOREV = 'C' */
/* >                                (LDV,M) if STOREV = 'R' and SIDE = 'L' */
/* >                                (LDV,N) if STOREV = 'R' and SIDE = 'R' */
/* >          The matrix V. See Further Details. */
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
/* >          T is DOUBLE PRECISION array, dimension (LDT,K) */
/* >          The triangular k by k matrix T in the representation of the */
/* >          block reflector. */
/* > \endverbatim */
/* > */
/* > \param[in] LDT */
/* > \verbatim */
/* >          LDT is INTEGER */
/* >          The leading dimension of the array T. LDT >= K. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* >          C is DOUBLE PRECISION array, dimension (LDC,N) */
/* >          On entry, the m by n matrix C. */
/* >          On exit, C is overwritten by H*C or H**T*C or C*H or C*H**T. */
/* > \endverbatim */
/* > */
/* > \param[in] LDC */
/* > \verbatim */
/* >          LDC is INTEGER */
/* >          The leading dimension of the array C. LDC >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension (LDWORK,K) */
/* > \endverbatim */
/* > */
/* > \param[in] LDWORK */
/* > \verbatim */
/* >          LDWORK is INTEGER */
/* >          The leading dimension of the array WORK. */
/* >          If SIDE = 'L', LDWORK >= max(1,N); */
/* >          if SIDE = 'R', LDWORK >= max(1,M). */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date June 2013 */

/* > \ingroup doubleOTHERauxiliary */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  The shape of the matrix V and the storage of the vectors which define */
/* >  the H(i) is best illustrated by the following example with n = 5 and */
/* >  k = 3. The elements equal to 1 are not stored; the corresponding */
/* >  array elements are modified but restored on exit. The rest of the */
/* >  array is not used. */
/* > */
/* >  DIRECT = 'F' and STOREV = 'C':         DIRECT = 'F' and STOREV = 'R': */
/* > */
/* >               V = (  1       )                 V = (  1 v1 v1 v1 v1 ) */
/* >                   ( v1  1    )                     (     1 v2 v2 v2 ) */
/* >                   ( v1 v2  1 )                     (        1 v3 v3 ) */
/* >                   ( v1 v2 v3 ) */
/* >                   ( v1 v2 v3 ) */
/* > */
/* >  DIRECT = 'B' and STOREV = 'C':         DIRECT = 'B' and STOREV = 'R': */
/* > */
/* >               V = ( v1 v2 v3 )                 V = ( v1 v1  1       ) */
/* >                   ( v1 v2 v3 )                     ( v2 v2 v2  1    ) */
/* >                   (  1 v2 v3 )                     ( v3 v3 v3 v3  1 ) */
/* >                   (     1 v3 ) */
/* >                   (        1 ) */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dlarfb_(char *side, char *trans, char *direct, char *
	storev, integer *m, integer *n, integer *k, doublereal *v, integer *
	ldv, doublereal *t, integer *ldt, doublereal *c__, integer *ldc, 
	doublereal *work, integer *ldwork, ftnlen side_len, ftnlen trans_len, 
	ftnlen direct_len, ftnlen storev_len)
{
    /* System generated locals */
    integer c_dim1, c_offset, t_dim1, t_offset, v_dim1, v_offset, work_dim1, 
	    work_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dtrmm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static char transt[1];


/*  -- LAPACK auxiliary routine (version 3.5.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     June 2013 */

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
/*     .. Executable Statements .. */

/*     Quick return if possible */

#line 233 "dlarfb.f"
    /* Parameter adjustments */
#line 233 "dlarfb.f"
    v_dim1 = *ldv;
#line 233 "dlarfb.f"
    v_offset = 1 + v_dim1;
#line 233 "dlarfb.f"
    v -= v_offset;
#line 233 "dlarfb.f"
    t_dim1 = *ldt;
#line 233 "dlarfb.f"
    t_offset = 1 + t_dim1;
#line 233 "dlarfb.f"
    t -= t_offset;
#line 233 "dlarfb.f"
    c_dim1 = *ldc;
#line 233 "dlarfb.f"
    c_offset = 1 + c_dim1;
#line 233 "dlarfb.f"
    c__ -= c_offset;
#line 233 "dlarfb.f"
    work_dim1 = *ldwork;
#line 233 "dlarfb.f"
    work_offset = 1 + work_dim1;
#line 233 "dlarfb.f"
    work -= work_offset;
#line 233 "dlarfb.f"

#line 233 "dlarfb.f"
    /* Function Body */
#line 233 "dlarfb.f"
    if (*m <= 0 || *n <= 0) {
#line 233 "dlarfb.f"
	return 0;
#line 233 "dlarfb.f"
    }

#line 236 "dlarfb.f"
    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {
#line 237 "dlarfb.f"
	*(unsigned char *)transt = 'T';
#line 238 "dlarfb.f"
    } else {
#line 239 "dlarfb.f"
	*(unsigned char *)transt = 'N';
#line 240 "dlarfb.f"
    }

#line 242 "dlarfb.f"
    if (lsame_(storev, "C", (ftnlen)1, (ftnlen)1)) {

#line 244 "dlarfb.f"
	if (lsame_(direct, "F", (ftnlen)1, (ftnlen)1)) {

/*           Let  V =  ( V1 )    (first K rows) */
/*                     ( V2 ) */
/*           where  V1  is unit lower triangular. */

#line 250 "dlarfb.f"
	    if (lsame_(side, "L", (ftnlen)1, (ftnlen)1)) {

/*              Form  H * C  or  H**T * C  where  C = ( C1 ) */
/*                                                    ( C2 ) */

/*              W := C**T * V  =  (C1**T * V1 + C2**T * V2)  (stored in WORK) */

/*              W := C1**T */

#line 259 "dlarfb.f"
		i__1 = *k;
#line 259 "dlarfb.f"
		for (j = 1; j <= i__1; ++j) {
#line 260 "dlarfb.f"
		    dcopy_(n, &c__[j + c_dim1], ldc, &work[j * work_dim1 + 1],
			     &c__1);
#line 261 "dlarfb.f"
/* L10: */
#line 261 "dlarfb.f"
		}

/*              W := W * V1 */

#line 265 "dlarfb.f"
		dtrmm_("Right", "Lower", "No transpose", "Unit", n, k, &c_b14,
			 &v[v_offset], ldv, &work[work_offset], ldwork, (
			ftnlen)5, (ftnlen)5, (ftnlen)12, (ftnlen)4);
#line 267 "dlarfb.f"
		if (*m > *k) {

/*                 W := W + C2**T * V2 */

#line 271 "dlarfb.f"
		    i__1 = *m - *k;
#line 271 "dlarfb.f"
		    dgemm_("Transpose", "No transpose", n, k, &i__1, &c_b14, &
			    c__[*k + 1 + c_dim1], ldc, &v[*k + 1 + v_dim1], 
			    ldv, &c_b14, &work[work_offset], ldwork, (ftnlen)
			    9, (ftnlen)12);
#line 274 "dlarfb.f"
		}

/*              W := W * T**T  or  W * T */

#line 278 "dlarfb.f"
		dtrmm_("Right", "Upper", transt, "Non-unit", n, k, &c_b14, &t[
			t_offset], ldt, &work[work_offset], ldwork, (ftnlen)5,
			 (ftnlen)5, (ftnlen)1, (ftnlen)8);

/*              C := C - V * W**T */

#line 283 "dlarfb.f"
		if (*m > *k) {

/*                 C2 := C2 - V2 * W**T */

#line 287 "dlarfb.f"
		    i__1 = *m - *k;
#line 287 "dlarfb.f"
		    dgemm_("No transpose", "Transpose", &i__1, n, k, &c_b25, &
			    v[*k + 1 + v_dim1], ldv, &work[work_offset], 
			    ldwork, &c_b14, &c__[*k + 1 + c_dim1], ldc, (
			    ftnlen)12, (ftnlen)9);
#line 290 "dlarfb.f"
		}

/*              W := W * V1**T */

#line 294 "dlarfb.f"
		dtrmm_("Right", "Lower", "Transpose", "Unit", n, k, &c_b14, &
			v[v_offset], ldv, &work[work_offset], ldwork, (ftnlen)
			5, (ftnlen)5, (ftnlen)9, (ftnlen)4);

/*              C1 := C1 - W**T */

#line 299 "dlarfb.f"
		i__1 = *k;
#line 299 "dlarfb.f"
		for (j = 1; j <= i__1; ++j) {
#line 300 "dlarfb.f"
		    i__2 = *n;
#line 300 "dlarfb.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 301 "dlarfb.f"
			c__[j + i__ * c_dim1] -= work[i__ + j * work_dim1];
#line 302 "dlarfb.f"
/* L20: */
#line 302 "dlarfb.f"
		    }
#line 303 "dlarfb.f"
/* L30: */
#line 303 "dlarfb.f"
		}

#line 305 "dlarfb.f"
	    } else if (lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {

/*              Form  C * H  or  C * H**T  where  C = ( C1  C2 ) */

/*              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK) */

/*              W := C1 */

#line 313 "dlarfb.f"
		i__1 = *k;
#line 313 "dlarfb.f"
		for (j = 1; j <= i__1; ++j) {
#line 314 "dlarfb.f"
		    dcopy_(m, &c__[j * c_dim1 + 1], &c__1, &work[j * 
			    work_dim1 + 1], &c__1);
#line 315 "dlarfb.f"
/* L40: */
#line 315 "dlarfb.f"
		}

/*              W := W * V1 */

#line 319 "dlarfb.f"
		dtrmm_("Right", "Lower", "No transpose", "Unit", m, k, &c_b14,
			 &v[v_offset], ldv, &work[work_offset], ldwork, (
			ftnlen)5, (ftnlen)5, (ftnlen)12, (ftnlen)4);
#line 321 "dlarfb.f"
		if (*n > *k) {

/*                 W := W + C2 * V2 */

#line 325 "dlarfb.f"
		    i__1 = *n - *k;
#line 325 "dlarfb.f"
		    dgemm_("No transpose", "No transpose", m, k, &i__1, &
			    c_b14, &c__[(*k + 1) * c_dim1 + 1], ldc, &v[*k + 
			    1 + v_dim1], ldv, &c_b14, &work[work_offset], 
			    ldwork, (ftnlen)12, (ftnlen)12);
#line 328 "dlarfb.f"
		}

/*              W := W * T  or  W * T**T */

#line 332 "dlarfb.f"
		dtrmm_("Right", "Upper", trans, "Non-unit", m, k, &c_b14, &t[
			t_offset], ldt, &work[work_offset], ldwork, (ftnlen)5,
			 (ftnlen)5, (ftnlen)1, (ftnlen)8);

/*              C := C - W * V**T */

#line 337 "dlarfb.f"
		if (*n > *k) {

/*                 C2 := C2 - W * V2**T */

#line 341 "dlarfb.f"
		    i__1 = *n - *k;
#line 341 "dlarfb.f"
		    dgemm_("No transpose", "Transpose", m, &i__1, k, &c_b25, &
			    work[work_offset], ldwork, &v[*k + 1 + v_dim1], 
			    ldv, &c_b14, &c__[(*k + 1) * c_dim1 + 1], ldc, (
			    ftnlen)12, (ftnlen)9);
#line 344 "dlarfb.f"
		}

/*              W := W * V1**T */

#line 348 "dlarfb.f"
		dtrmm_("Right", "Lower", "Transpose", "Unit", m, k, &c_b14, &
			v[v_offset], ldv, &work[work_offset], ldwork, (ftnlen)
			5, (ftnlen)5, (ftnlen)9, (ftnlen)4);

/*              C1 := C1 - W */

#line 353 "dlarfb.f"
		i__1 = *k;
#line 353 "dlarfb.f"
		for (j = 1; j <= i__1; ++j) {
#line 354 "dlarfb.f"
		    i__2 = *m;
#line 354 "dlarfb.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 355 "dlarfb.f"
			c__[i__ + j * c_dim1] -= work[i__ + j * work_dim1];
#line 356 "dlarfb.f"
/* L50: */
#line 356 "dlarfb.f"
		    }
#line 357 "dlarfb.f"
/* L60: */
#line 357 "dlarfb.f"
		}
#line 358 "dlarfb.f"
	    }

#line 360 "dlarfb.f"
	} else {

/*           Let  V =  ( V1 ) */
/*                     ( V2 )    (last K rows) */
/*           where  V2  is unit upper triangular. */

#line 366 "dlarfb.f"
	    if (lsame_(side, "L", (ftnlen)1, (ftnlen)1)) {

/*              Form  H * C  or  H**T * C  where  C = ( C1 ) */
/*                                                    ( C2 ) */

/*              W := C**T * V  =  (C1**T * V1 + C2**T * V2)  (stored in WORK) */

/*              W := C2**T */

#line 375 "dlarfb.f"
		i__1 = *k;
#line 375 "dlarfb.f"
		for (j = 1; j <= i__1; ++j) {
#line 376 "dlarfb.f"
		    dcopy_(n, &c__[*m - *k + j + c_dim1], ldc, &work[j * 
			    work_dim1 + 1], &c__1);
#line 377 "dlarfb.f"
/* L70: */
#line 377 "dlarfb.f"
		}

/*              W := W * V2 */

#line 381 "dlarfb.f"
		dtrmm_("Right", "Upper", "No transpose", "Unit", n, k, &c_b14,
			 &v[*m - *k + 1 + v_dim1], ldv, &work[work_offset], 
			ldwork, (ftnlen)5, (ftnlen)5, (ftnlen)12, (ftnlen)4);
#line 383 "dlarfb.f"
		if (*m > *k) {

/*                 W := W + C1**T * V1 */

#line 387 "dlarfb.f"
		    i__1 = *m - *k;
#line 387 "dlarfb.f"
		    dgemm_("Transpose", "No transpose", n, k, &i__1, &c_b14, &
			    c__[c_offset], ldc, &v[v_offset], ldv, &c_b14, &
			    work[work_offset], ldwork, (ftnlen)9, (ftnlen)12);
#line 389 "dlarfb.f"
		}

/*              W := W * T**T  or  W * T */

#line 393 "dlarfb.f"
		dtrmm_("Right", "Lower", transt, "Non-unit", n, k, &c_b14, &t[
			t_offset], ldt, &work[work_offset], ldwork, (ftnlen)5,
			 (ftnlen)5, (ftnlen)1, (ftnlen)8);

/*              C := C - V * W**T */

#line 398 "dlarfb.f"
		if (*m > *k) {

/*                 C1 := C1 - V1 * W**T */

#line 402 "dlarfb.f"
		    i__1 = *m - *k;
#line 402 "dlarfb.f"
		    dgemm_("No transpose", "Transpose", &i__1, n, k, &c_b25, &
			    v[v_offset], ldv, &work[work_offset], ldwork, &
			    c_b14, &c__[c_offset], ldc, (ftnlen)12, (ftnlen)9)
			    ;
#line 404 "dlarfb.f"
		}

/*              W := W * V2**T */

#line 408 "dlarfb.f"
		dtrmm_("Right", "Upper", "Transpose", "Unit", n, k, &c_b14, &
			v[*m - *k + 1 + v_dim1], ldv, &work[work_offset], 
			ldwork, (ftnlen)5, (ftnlen)5, (ftnlen)9, (ftnlen)4);

/*              C2 := C2 - W**T */

#line 413 "dlarfb.f"
		i__1 = *k;
#line 413 "dlarfb.f"
		for (j = 1; j <= i__1; ++j) {
#line 414 "dlarfb.f"
		    i__2 = *n;
#line 414 "dlarfb.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 415 "dlarfb.f"
			c__[*m - *k + j + i__ * c_dim1] -= work[i__ + j * 
				work_dim1];
#line 416 "dlarfb.f"
/* L80: */
#line 416 "dlarfb.f"
		    }
#line 417 "dlarfb.f"
/* L90: */
#line 417 "dlarfb.f"
		}

#line 419 "dlarfb.f"
	    } else if (lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {

/*              Form  C * H  or  C * H**T  where  C = ( C1  C2 ) */

/*              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK) */

/*              W := C2 */

#line 427 "dlarfb.f"
		i__1 = *k;
#line 427 "dlarfb.f"
		for (j = 1; j <= i__1; ++j) {
#line 428 "dlarfb.f"
		    dcopy_(m, &c__[(*n - *k + j) * c_dim1 + 1], &c__1, &work[
			    j * work_dim1 + 1], &c__1);
#line 429 "dlarfb.f"
/* L100: */
#line 429 "dlarfb.f"
		}

/*              W := W * V2 */

#line 433 "dlarfb.f"
		dtrmm_("Right", "Upper", "No transpose", "Unit", m, k, &c_b14,
			 &v[*n - *k + 1 + v_dim1], ldv, &work[work_offset], 
			ldwork, (ftnlen)5, (ftnlen)5, (ftnlen)12, (ftnlen)4);
#line 435 "dlarfb.f"
		if (*n > *k) {

/*                 W := W + C1 * V1 */

#line 439 "dlarfb.f"
		    i__1 = *n - *k;
#line 439 "dlarfb.f"
		    dgemm_("No transpose", "No transpose", m, k, &i__1, &
			    c_b14, &c__[c_offset], ldc, &v[v_offset], ldv, &
			    c_b14, &work[work_offset], ldwork, (ftnlen)12, (
			    ftnlen)12);
#line 441 "dlarfb.f"
		}

/*              W := W * T  or  W * T**T */

#line 445 "dlarfb.f"
		dtrmm_("Right", "Lower", trans, "Non-unit", m, k, &c_b14, &t[
			t_offset], ldt, &work[work_offset], ldwork, (ftnlen)5,
			 (ftnlen)5, (ftnlen)1, (ftnlen)8);

/*              C := C - W * V**T */

#line 450 "dlarfb.f"
		if (*n > *k) {

/*                 C1 := C1 - W * V1**T */

#line 454 "dlarfb.f"
		    i__1 = *n - *k;
#line 454 "dlarfb.f"
		    dgemm_("No transpose", "Transpose", m, &i__1, k, &c_b25, &
			    work[work_offset], ldwork, &v[v_offset], ldv, &
			    c_b14, &c__[c_offset], ldc, (ftnlen)12, (ftnlen)9)
			    ;
#line 456 "dlarfb.f"
		}

/*              W := W * V2**T */

#line 460 "dlarfb.f"
		dtrmm_("Right", "Upper", "Transpose", "Unit", m, k, &c_b14, &
			v[*n - *k + 1 + v_dim1], ldv, &work[work_offset], 
			ldwork, (ftnlen)5, (ftnlen)5, (ftnlen)9, (ftnlen)4);

/*              C2 := C2 - W */

#line 465 "dlarfb.f"
		i__1 = *k;
#line 465 "dlarfb.f"
		for (j = 1; j <= i__1; ++j) {
#line 466 "dlarfb.f"
		    i__2 = *m;
#line 466 "dlarfb.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 467 "dlarfb.f"
			c__[i__ + (*n - *k + j) * c_dim1] -= work[i__ + j * 
				work_dim1];
#line 468 "dlarfb.f"
/* L110: */
#line 468 "dlarfb.f"
		    }
#line 469 "dlarfb.f"
/* L120: */
#line 469 "dlarfb.f"
		}
#line 470 "dlarfb.f"
	    }
#line 471 "dlarfb.f"
	}

#line 473 "dlarfb.f"
    } else if (lsame_(storev, "R", (ftnlen)1, (ftnlen)1)) {

#line 475 "dlarfb.f"
	if (lsame_(direct, "F", (ftnlen)1, (ftnlen)1)) {

/*           Let  V =  ( V1  V2 )    (V1: first K columns) */
/*           where  V1  is unit upper triangular. */

#line 480 "dlarfb.f"
	    if (lsame_(side, "L", (ftnlen)1, (ftnlen)1)) {

/*              Form  H * C  or  H**T * C  where  C = ( C1 ) */
/*                                                    ( C2 ) */

/*              W := C**T * V**T  =  (C1**T * V1**T + C2**T * V2**T) (stored in WORK) */

/*              W := C1**T */

#line 489 "dlarfb.f"
		i__1 = *k;
#line 489 "dlarfb.f"
		for (j = 1; j <= i__1; ++j) {
#line 490 "dlarfb.f"
		    dcopy_(n, &c__[j + c_dim1], ldc, &work[j * work_dim1 + 1],
			     &c__1);
#line 491 "dlarfb.f"
/* L130: */
#line 491 "dlarfb.f"
		}

/*              W := W * V1**T */

#line 495 "dlarfb.f"
		dtrmm_("Right", "Upper", "Transpose", "Unit", n, k, &c_b14, &
			v[v_offset], ldv, &work[work_offset], ldwork, (ftnlen)
			5, (ftnlen)5, (ftnlen)9, (ftnlen)4);
#line 497 "dlarfb.f"
		if (*m > *k) {

/*                 W := W + C2**T * V2**T */

#line 501 "dlarfb.f"
		    i__1 = *m - *k;
#line 501 "dlarfb.f"
		    dgemm_("Transpose", "Transpose", n, k, &i__1, &c_b14, &
			    c__[*k + 1 + c_dim1], ldc, &v[(*k + 1) * v_dim1 + 
			    1], ldv, &c_b14, &work[work_offset], ldwork, (
			    ftnlen)9, (ftnlen)9);
#line 504 "dlarfb.f"
		}

/*              W := W * T**T  or  W * T */

#line 508 "dlarfb.f"
		dtrmm_("Right", "Upper", transt, "Non-unit", n, k, &c_b14, &t[
			t_offset], ldt, &work[work_offset], ldwork, (ftnlen)5,
			 (ftnlen)5, (ftnlen)1, (ftnlen)8);

/*              C := C - V**T * W**T */

#line 513 "dlarfb.f"
		if (*m > *k) {

/*                 C2 := C2 - V2**T * W**T */

#line 517 "dlarfb.f"
		    i__1 = *m - *k;
#line 517 "dlarfb.f"
		    dgemm_("Transpose", "Transpose", &i__1, n, k, &c_b25, &v[(
			    *k + 1) * v_dim1 + 1], ldv, &work[work_offset], 
			    ldwork, &c_b14, &c__[*k + 1 + c_dim1], ldc, (
			    ftnlen)9, (ftnlen)9);
#line 520 "dlarfb.f"
		}

/*              W := W * V1 */

#line 524 "dlarfb.f"
		dtrmm_("Right", "Upper", "No transpose", "Unit", n, k, &c_b14,
			 &v[v_offset], ldv, &work[work_offset], ldwork, (
			ftnlen)5, (ftnlen)5, (ftnlen)12, (ftnlen)4);

/*              C1 := C1 - W**T */

#line 529 "dlarfb.f"
		i__1 = *k;
#line 529 "dlarfb.f"
		for (j = 1; j <= i__1; ++j) {
#line 530 "dlarfb.f"
		    i__2 = *n;
#line 530 "dlarfb.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 531 "dlarfb.f"
			c__[j + i__ * c_dim1] -= work[i__ + j * work_dim1];
#line 532 "dlarfb.f"
/* L140: */
#line 532 "dlarfb.f"
		    }
#line 533 "dlarfb.f"
/* L150: */
#line 533 "dlarfb.f"
		}

#line 535 "dlarfb.f"
	    } else if (lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {

/*              Form  C * H  or  C * H**T  where  C = ( C1  C2 ) */

/*              W := C * V**T  =  (C1*V1**T + C2*V2**T)  (stored in WORK) */

/*              W := C1 */

#line 543 "dlarfb.f"
		i__1 = *k;
#line 543 "dlarfb.f"
		for (j = 1; j <= i__1; ++j) {
#line 544 "dlarfb.f"
		    dcopy_(m, &c__[j * c_dim1 + 1], &c__1, &work[j * 
			    work_dim1 + 1], &c__1);
#line 545 "dlarfb.f"
/* L160: */
#line 545 "dlarfb.f"
		}

/*              W := W * V1**T */

#line 549 "dlarfb.f"
		dtrmm_("Right", "Upper", "Transpose", "Unit", m, k, &c_b14, &
			v[v_offset], ldv, &work[work_offset], ldwork, (ftnlen)
			5, (ftnlen)5, (ftnlen)9, (ftnlen)4);
#line 551 "dlarfb.f"
		if (*n > *k) {

/*                 W := W + C2 * V2**T */

#line 555 "dlarfb.f"
		    i__1 = *n - *k;
#line 555 "dlarfb.f"
		    dgemm_("No transpose", "Transpose", m, k, &i__1, &c_b14, &
			    c__[(*k + 1) * c_dim1 + 1], ldc, &v[(*k + 1) * 
			    v_dim1 + 1], ldv, &c_b14, &work[work_offset], 
			    ldwork, (ftnlen)12, (ftnlen)9);
#line 558 "dlarfb.f"
		}

/*              W := W * T  or  W * T**T */

#line 562 "dlarfb.f"
		dtrmm_("Right", "Upper", trans, "Non-unit", m, k, &c_b14, &t[
			t_offset], ldt, &work[work_offset], ldwork, (ftnlen)5,
			 (ftnlen)5, (ftnlen)1, (ftnlen)8);

/*              C := C - W * V */

#line 567 "dlarfb.f"
		if (*n > *k) {

/*                 C2 := C2 - W * V2 */

#line 571 "dlarfb.f"
		    i__1 = *n - *k;
#line 571 "dlarfb.f"
		    dgemm_("No transpose", "No transpose", m, &i__1, k, &
			    c_b25, &work[work_offset], ldwork, &v[(*k + 1) * 
			    v_dim1 + 1], ldv, &c_b14, &c__[(*k + 1) * c_dim1 
			    + 1], ldc, (ftnlen)12, (ftnlen)12);
#line 574 "dlarfb.f"
		}

/*              W := W * V1 */

#line 578 "dlarfb.f"
		dtrmm_("Right", "Upper", "No transpose", "Unit", m, k, &c_b14,
			 &v[v_offset], ldv, &work[work_offset], ldwork, (
			ftnlen)5, (ftnlen)5, (ftnlen)12, (ftnlen)4);

/*              C1 := C1 - W */

#line 583 "dlarfb.f"
		i__1 = *k;
#line 583 "dlarfb.f"
		for (j = 1; j <= i__1; ++j) {
#line 584 "dlarfb.f"
		    i__2 = *m;
#line 584 "dlarfb.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 585 "dlarfb.f"
			c__[i__ + j * c_dim1] -= work[i__ + j * work_dim1];
#line 586 "dlarfb.f"
/* L170: */
#line 586 "dlarfb.f"
		    }
#line 587 "dlarfb.f"
/* L180: */
#line 587 "dlarfb.f"
		}

#line 589 "dlarfb.f"
	    }

#line 591 "dlarfb.f"
	} else {

/*           Let  V =  ( V1  V2 )    (V2: last K columns) */
/*           where  V2  is unit lower triangular. */

#line 596 "dlarfb.f"
	    if (lsame_(side, "L", (ftnlen)1, (ftnlen)1)) {

/*              Form  H * C  or  H**T * C  where  C = ( C1 ) */
/*                                                    ( C2 ) */

/*              W := C**T * V**T  =  (C1**T * V1**T + C2**T * V2**T) (stored in WORK) */

/*              W := C2**T */

#line 605 "dlarfb.f"
		i__1 = *k;
#line 605 "dlarfb.f"
		for (j = 1; j <= i__1; ++j) {
#line 606 "dlarfb.f"
		    dcopy_(n, &c__[*m - *k + j + c_dim1], ldc, &work[j * 
			    work_dim1 + 1], &c__1);
#line 607 "dlarfb.f"
/* L190: */
#line 607 "dlarfb.f"
		}

/*              W := W * V2**T */

#line 611 "dlarfb.f"
		dtrmm_("Right", "Lower", "Transpose", "Unit", n, k, &c_b14, &
			v[(*m - *k + 1) * v_dim1 + 1], ldv, &work[work_offset]
			, ldwork, (ftnlen)5, (ftnlen)5, (ftnlen)9, (ftnlen)4);
#line 613 "dlarfb.f"
		if (*m > *k) {

/*                 W := W + C1**T * V1**T */

#line 617 "dlarfb.f"
		    i__1 = *m - *k;
#line 617 "dlarfb.f"
		    dgemm_("Transpose", "Transpose", n, k, &i__1, &c_b14, &
			    c__[c_offset], ldc, &v[v_offset], ldv, &c_b14, &
			    work[work_offset], ldwork, (ftnlen)9, (ftnlen)9);
#line 619 "dlarfb.f"
		}

/*              W := W * T**T  or  W * T */

#line 623 "dlarfb.f"
		dtrmm_("Right", "Lower", transt, "Non-unit", n, k, &c_b14, &t[
			t_offset], ldt, &work[work_offset], ldwork, (ftnlen)5,
			 (ftnlen)5, (ftnlen)1, (ftnlen)8);

/*              C := C - V**T * W**T */

#line 628 "dlarfb.f"
		if (*m > *k) {

/*                 C1 := C1 - V1**T * W**T */

#line 632 "dlarfb.f"
		    i__1 = *m - *k;
#line 632 "dlarfb.f"
		    dgemm_("Transpose", "Transpose", &i__1, n, k, &c_b25, &v[
			    v_offset], ldv, &work[work_offset], ldwork, &
			    c_b14, &c__[c_offset], ldc, (ftnlen)9, (ftnlen)9);
#line 634 "dlarfb.f"
		}

/*              W := W * V2 */

#line 638 "dlarfb.f"
		dtrmm_("Right", "Lower", "No transpose", "Unit", n, k, &c_b14,
			 &v[(*m - *k + 1) * v_dim1 + 1], ldv, &work[
			work_offset], ldwork, (ftnlen)5, (ftnlen)5, (ftnlen)
			12, (ftnlen)4);

/*              C2 := C2 - W**T */

#line 643 "dlarfb.f"
		i__1 = *k;
#line 643 "dlarfb.f"
		for (j = 1; j <= i__1; ++j) {
#line 644 "dlarfb.f"
		    i__2 = *n;
#line 644 "dlarfb.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 645 "dlarfb.f"
			c__[*m - *k + j + i__ * c_dim1] -= work[i__ + j * 
				work_dim1];
#line 646 "dlarfb.f"
/* L200: */
#line 646 "dlarfb.f"
		    }
#line 647 "dlarfb.f"
/* L210: */
#line 647 "dlarfb.f"
		}

#line 649 "dlarfb.f"
	    } else if (lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {

/*              Form  C * H  or  C * H'  where  C = ( C1  C2 ) */

/*              W := C * V**T  =  (C1*V1**T + C2*V2**T)  (stored in WORK) */

/*              W := C2 */

#line 657 "dlarfb.f"
		i__1 = *k;
#line 657 "dlarfb.f"
		for (j = 1; j <= i__1; ++j) {
#line 658 "dlarfb.f"
		    dcopy_(m, &c__[(*n - *k + j) * c_dim1 + 1], &c__1, &work[
			    j * work_dim1 + 1], &c__1);
#line 659 "dlarfb.f"
/* L220: */
#line 659 "dlarfb.f"
		}

/*              W := W * V2**T */

#line 663 "dlarfb.f"
		dtrmm_("Right", "Lower", "Transpose", "Unit", m, k, &c_b14, &
			v[(*n - *k + 1) * v_dim1 + 1], ldv, &work[work_offset]
			, ldwork, (ftnlen)5, (ftnlen)5, (ftnlen)9, (ftnlen)4);
#line 665 "dlarfb.f"
		if (*n > *k) {

/*                 W := W + C1 * V1**T */

#line 669 "dlarfb.f"
		    i__1 = *n - *k;
#line 669 "dlarfb.f"
		    dgemm_("No transpose", "Transpose", m, k, &i__1, &c_b14, &
			    c__[c_offset], ldc, &v[v_offset], ldv, &c_b14, &
			    work[work_offset], ldwork, (ftnlen)12, (ftnlen)9);
#line 671 "dlarfb.f"
		}

/*              W := W * T  or  W * T**T */

#line 675 "dlarfb.f"
		dtrmm_("Right", "Lower", trans, "Non-unit", m, k, &c_b14, &t[
			t_offset], ldt, &work[work_offset], ldwork, (ftnlen)5,
			 (ftnlen)5, (ftnlen)1, (ftnlen)8);

/*              C := C - W * V */

#line 680 "dlarfb.f"
		if (*n > *k) {

/*                 C1 := C1 - W * V1 */

#line 684 "dlarfb.f"
		    i__1 = *n - *k;
#line 684 "dlarfb.f"
		    dgemm_("No transpose", "No transpose", m, &i__1, k, &
			    c_b25, &work[work_offset], ldwork, &v[v_offset], 
			    ldv, &c_b14, &c__[c_offset], ldc, (ftnlen)12, (
			    ftnlen)12);
#line 686 "dlarfb.f"
		}

/*              W := W * V2 */

#line 690 "dlarfb.f"
		dtrmm_("Right", "Lower", "No transpose", "Unit", m, k, &c_b14,
			 &v[(*n - *k + 1) * v_dim1 + 1], ldv, &work[
			work_offset], ldwork, (ftnlen)5, (ftnlen)5, (ftnlen)
			12, (ftnlen)4);

/*              C1 := C1 - W */

#line 695 "dlarfb.f"
		i__1 = *k;
#line 695 "dlarfb.f"
		for (j = 1; j <= i__1; ++j) {
#line 696 "dlarfb.f"
		    i__2 = *m;
#line 696 "dlarfb.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 697 "dlarfb.f"
			c__[i__ + (*n - *k + j) * c_dim1] -= work[i__ + j * 
				work_dim1];
#line 698 "dlarfb.f"
/* L230: */
#line 698 "dlarfb.f"
		    }
#line 699 "dlarfb.f"
/* L240: */
#line 699 "dlarfb.f"
		}

#line 701 "dlarfb.f"
	    }

#line 703 "dlarfb.f"
	}
#line 704 "dlarfb.f"
    }

#line 706 "dlarfb.f"
    return 0;

/*     End of DLARFB */

} /* dlarfb_ */

