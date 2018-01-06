#line 1 "slarfb.f"
/* slarfb.f -- translated by f2c (version 20100827).
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

#line 1 "slarfb.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b14 = 1.;
static doublereal c_b25 = -1.;

/* > \brief \b SLARFB applies a block reflector or its transpose to a general rectangular matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLARFB + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slarfb.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slarfb.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slarfb.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLARFB( SIDE, TRANS, DIRECT, STOREV, M, N, K, V, LDV, */
/*                          T, LDT, C, LDC, WORK, LDWORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          DIRECT, SIDE, STOREV, TRANS */
/*       INTEGER            K, LDC, LDT, LDV, LDWORK, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               C( LDC, * ), T( LDT, * ), V( LDV, * ), */
/*      $                   WORK( LDWORK, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLARFB applies a real block reflector H or its transpose H**T to a */
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
/* >          V is REAL array, dimension */
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
/* >          T is REAL array, dimension (LDT,K) */
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
/* >          C is REAL array, dimension (LDC,N) */
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
/* >          WORK is REAL array, dimension (LDWORK,K) */
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

/* > \ingroup realOTHERauxiliary */

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
/* Subroutine */ int slarfb_(char *side, char *trans, char *direct, char *
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
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int sgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen),
	     scopy_(integer *, doublereal *, integer *, doublereal *, integer 
	    *), strmm_(char *, char *, char *, char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen, ftnlen, ftnlen, ftnlen);
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

#line 233 "slarfb.f"
    /* Parameter adjustments */
#line 233 "slarfb.f"
    v_dim1 = *ldv;
#line 233 "slarfb.f"
    v_offset = 1 + v_dim1;
#line 233 "slarfb.f"
    v -= v_offset;
#line 233 "slarfb.f"
    t_dim1 = *ldt;
#line 233 "slarfb.f"
    t_offset = 1 + t_dim1;
#line 233 "slarfb.f"
    t -= t_offset;
#line 233 "slarfb.f"
    c_dim1 = *ldc;
#line 233 "slarfb.f"
    c_offset = 1 + c_dim1;
#line 233 "slarfb.f"
    c__ -= c_offset;
#line 233 "slarfb.f"
    work_dim1 = *ldwork;
#line 233 "slarfb.f"
    work_offset = 1 + work_dim1;
#line 233 "slarfb.f"
    work -= work_offset;
#line 233 "slarfb.f"

#line 233 "slarfb.f"
    /* Function Body */
#line 233 "slarfb.f"
    if (*m <= 0 || *n <= 0) {
#line 233 "slarfb.f"
	return 0;
#line 233 "slarfb.f"
    }

#line 236 "slarfb.f"
    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {
#line 237 "slarfb.f"
	*(unsigned char *)transt = 'T';
#line 238 "slarfb.f"
    } else {
#line 239 "slarfb.f"
	*(unsigned char *)transt = 'N';
#line 240 "slarfb.f"
    }

#line 242 "slarfb.f"
    if (lsame_(storev, "C", (ftnlen)1, (ftnlen)1)) {

#line 244 "slarfb.f"
	if (lsame_(direct, "F", (ftnlen)1, (ftnlen)1)) {

/*           Let  V =  ( V1 )    (first K rows) */
/*                     ( V2 ) */
/*           where  V1  is unit lower triangular. */

#line 250 "slarfb.f"
	    if (lsame_(side, "L", (ftnlen)1, (ftnlen)1)) {

/*              Form  H * C  or  H**T * C  where  C = ( C1 ) */
/*                                                    ( C2 ) */

/*              W := C**T * V  =  (C1**T * V1 + C2**T * V2)  (stored in WORK) */

/*              W := C1**T */

#line 259 "slarfb.f"
		i__1 = *k;
#line 259 "slarfb.f"
		for (j = 1; j <= i__1; ++j) {
#line 260 "slarfb.f"
		    scopy_(n, &c__[j + c_dim1], ldc, &work[j * work_dim1 + 1],
			     &c__1);
#line 261 "slarfb.f"
/* L10: */
#line 261 "slarfb.f"
		}

/*              W := W * V1 */

#line 265 "slarfb.f"
		strmm_("Right", "Lower", "No transpose", "Unit", n, k, &c_b14,
			 &v[v_offset], ldv, &work[work_offset], ldwork, (
			ftnlen)5, (ftnlen)5, (ftnlen)12, (ftnlen)4);
#line 267 "slarfb.f"
		if (*m > *k) {

/*                 W := W + C2**T * V2 */

#line 271 "slarfb.f"
		    i__1 = *m - *k;
#line 271 "slarfb.f"
		    sgemm_("Transpose", "No transpose", n, k, &i__1, &c_b14, &
			    c__[*k + 1 + c_dim1], ldc, &v[*k + 1 + v_dim1], 
			    ldv, &c_b14, &work[work_offset], ldwork, (ftnlen)
			    9, (ftnlen)12);
#line 274 "slarfb.f"
		}

/*              W := W * T**T  or  W * T */

#line 278 "slarfb.f"
		strmm_("Right", "Upper", transt, "Non-unit", n, k, &c_b14, &t[
			t_offset], ldt, &work[work_offset], ldwork, (ftnlen)5,
			 (ftnlen)5, (ftnlen)1, (ftnlen)8);

/*              C := C - V * W**T */

#line 283 "slarfb.f"
		if (*m > *k) {

/*                 C2 := C2 - V2 * W**T */

#line 287 "slarfb.f"
		    i__1 = *m - *k;
#line 287 "slarfb.f"
		    sgemm_("No transpose", "Transpose", &i__1, n, k, &c_b25, &
			    v[*k + 1 + v_dim1], ldv, &work[work_offset], 
			    ldwork, &c_b14, &c__[*k + 1 + c_dim1], ldc, (
			    ftnlen)12, (ftnlen)9);
#line 290 "slarfb.f"
		}

/*              W := W * V1**T */

#line 294 "slarfb.f"
		strmm_("Right", "Lower", "Transpose", "Unit", n, k, &c_b14, &
			v[v_offset], ldv, &work[work_offset], ldwork, (ftnlen)
			5, (ftnlen)5, (ftnlen)9, (ftnlen)4);

/*              C1 := C1 - W**T */

#line 299 "slarfb.f"
		i__1 = *k;
#line 299 "slarfb.f"
		for (j = 1; j <= i__1; ++j) {
#line 300 "slarfb.f"
		    i__2 = *n;
#line 300 "slarfb.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 301 "slarfb.f"
			c__[j + i__ * c_dim1] -= work[i__ + j * work_dim1];
#line 302 "slarfb.f"
/* L20: */
#line 302 "slarfb.f"
		    }
#line 303 "slarfb.f"
/* L30: */
#line 303 "slarfb.f"
		}

#line 305 "slarfb.f"
	    } else if (lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {

/*              Form  C * H  or  C * H**T  where  C = ( C1  C2 ) */

/*              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK) */

/*              W := C1 */

#line 313 "slarfb.f"
		i__1 = *k;
#line 313 "slarfb.f"
		for (j = 1; j <= i__1; ++j) {
#line 314 "slarfb.f"
		    scopy_(m, &c__[j * c_dim1 + 1], &c__1, &work[j * 
			    work_dim1 + 1], &c__1);
#line 315 "slarfb.f"
/* L40: */
#line 315 "slarfb.f"
		}

/*              W := W * V1 */

#line 319 "slarfb.f"
		strmm_("Right", "Lower", "No transpose", "Unit", m, k, &c_b14,
			 &v[v_offset], ldv, &work[work_offset], ldwork, (
			ftnlen)5, (ftnlen)5, (ftnlen)12, (ftnlen)4);
#line 321 "slarfb.f"
		if (*n > *k) {

/*                 W := W + C2 * V2 */

#line 325 "slarfb.f"
		    i__1 = *n - *k;
#line 325 "slarfb.f"
		    sgemm_("No transpose", "No transpose", m, k, &i__1, &
			    c_b14, &c__[(*k + 1) * c_dim1 + 1], ldc, &v[*k + 
			    1 + v_dim1], ldv, &c_b14, &work[work_offset], 
			    ldwork, (ftnlen)12, (ftnlen)12);
#line 328 "slarfb.f"
		}

/*              W := W * T  or  W * T**T */

#line 332 "slarfb.f"
		strmm_("Right", "Upper", trans, "Non-unit", m, k, &c_b14, &t[
			t_offset], ldt, &work[work_offset], ldwork, (ftnlen)5,
			 (ftnlen)5, (ftnlen)1, (ftnlen)8);

/*              C := C - W * V**T */

#line 337 "slarfb.f"
		if (*n > *k) {

/*                 C2 := C2 - W * V2**T */

#line 341 "slarfb.f"
		    i__1 = *n - *k;
#line 341 "slarfb.f"
		    sgemm_("No transpose", "Transpose", m, &i__1, k, &c_b25, &
			    work[work_offset], ldwork, &v[*k + 1 + v_dim1], 
			    ldv, &c_b14, &c__[(*k + 1) * c_dim1 + 1], ldc, (
			    ftnlen)12, (ftnlen)9);
#line 344 "slarfb.f"
		}

/*              W := W * V1**T */

#line 348 "slarfb.f"
		strmm_("Right", "Lower", "Transpose", "Unit", m, k, &c_b14, &
			v[v_offset], ldv, &work[work_offset], ldwork, (ftnlen)
			5, (ftnlen)5, (ftnlen)9, (ftnlen)4);

/*              C1 := C1 - W */

#line 353 "slarfb.f"
		i__1 = *k;
#line 353 "slarfb.f"
		for (j = 1; j <= i__1; ++j) {
#line 354 "slarfb.f"
		    i__2 = *m;
#line 354 "slarfb.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 355 "slarfb.f"
			c__[i__ + j * c_dim1] -= work[i__ + j * work_dim1];
#line 356 "slarfb.f"
/* L50: */
#line 356 "slarfb.f"
		    }
#line 357 "slarfb.f"
/* L60: */
#line 357 "slarfb.f"
		}
#line 358 "slarfb.f"
	    }

#line 360 "slarfb.f"
	} else {

/*           Let  V =  ( V1 ) */
/*                     ( V2 )    (last K rows) */
/*           where  V2  is unit upper triangular. */

#line 366 "slarfb.f"
	    if (lsame_(side, "L", (ftnlen)1, (ftnlen)1)) {

/*              Form  H * C  or  H**T * C  where  C = ( C1 ) */
/*                                                    ( C2 ) */

/*              W := C**T * V  =  (C1**T * V1 + C2**T * V2)  (stored in WORK) */

/*              W := C2**T */

#line 375 "slarfb.f"
		i__1 = *k;
#line 375 "slarfb.f"
		for (j = 1; j <= i__1; ++j) {
#line 376 "slarfb.f"
		    scopy_(n, &c__[*m - *k + j + c_dim1], ldc, &work[j * 
			    work_dim1 + 1], &c__1);
#line 377 "slarfb.f"
/* L70: */
#line 377 "slarfb.f"
		}

/*              W := W * V2 */

#line 381 "slarfb.f"
		strmm_("Right", "Upper", "No transpose", "Unit", n, k, &c_b14,
			 &v[*m - *k + 1 + v_dim1], ldv, &work[work_offset], 
			ldwork, (ftnlen)5, (ftnlen)5, (ftnlen)12, (ftnlen)4);
#line 383 "slarfb.f"
		if (*m > *k) {

/*                 W := W + C1**T * V1 */

#line 387 "slarfb.f"
		    i__1 = *m - *k;
#line 387 "slarfb.f"
		    sgemm_("Transpose", "No transpose", n, k, &i__1, &c_b14, &
			    c__[c_offset], ldc, &v[v_offset], ldv, &c_b14, &
			    work[work_offset], ldwork, (ftnlen)9, (ftnlen)12);
#line 389 "slarfb.f"
		}

/*              W := W * T**T  or  W * T */

#line 393 "slarfb.f"
		strmm_("Right", "Lower", transt, "Non-unit", n, k, &c_b14, &t[
			t_offset], ldt, &work[work_offset], ldwork, (ftnlen)5,
			 (ftnlen)5, (ftnlen)1, (ftnlen)8);

/*              C := C - V * W**T */

#line 398 "slarfb.f"
		if (*m > *k) {

/*                 C1 := C1 - V1 * W**T */

#line 402 "slarfb.f"
		    i__1 = *m - *k;
#line 402 "slarfb.f"
		    sgemm_("No transpose", "Transpose", &i__1, n, k, &c_b25, &
			    v[v_offset], ldv, &work[work_offset], ldwork, &
			    c_b14, &c__[c_offset], ldc, (ftnlen)12, (ftnlen)9)
			    ;
#line 404 "slarfb.f"
		}

/*              W := W * V2**T */

#line 408 "slarfb.f"
		strmm_("Right", "Upper", "Transpose", "Unit", n, k, &c_b14, &
			v[*m - *k + 1 + v_dim1], ldv, &work[work_offset], 
			ldwork, (ftnlen)5, (ftnlen)5, (ftnlen)9, (ftnlen)4);

/*              C2 := C2 - W**T */

#line 413 "slarfb.f"
		i__1 = *k;
#line 413 "slarfb.f"
		for (j = 1; j <= i__1; ++j) {
#line 414 "slarfb.f"
		    i__2 = *n;
#line 414 "slarfb.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 415 "slarfb.f"
			c__[*m - *k + j + i__ * c_dim1] -= work[i__ + j * 
				work_dim1];
#line 416 "slarfb.f"
/* L80: */
#line 416 "slarfb.f"
		    }
#line 417 "slarfb.f"
/* L90: */
#line 417 "slarfb.f"
		}

#line 419 "slarfb.f"
	    } else if (lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {

/*              Form  C * H  or  C * H'  where  C = ( C1  C2 ) */

/*              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK) */

/*              W := C2 */

#line 427 "slarfb.f"
		i__1 = *k;
#line 427 "slarfb.f"
		for (j = 1; j <= i__1; ++j) {
#line 428 "slarfb.f"
		    scopy_(m, &c__[(*n - *k + j) * c_dim1 + 1], &c__1, &work[
			    j * work_dim1 + 1], &c__1);
#line 429 "slarfb.f"
/* L100: */
#line 429 "slarfb.f"
		}

/*              W := W * V2 */

#line 433 "slarfb.f"
		strmm_("Right", "Upper", "No transpose", "Unit", m, k, &c_b14,
			 &v[*n - *k + 1 + v_dim1], ldv, &work[work_offset], 
			ldwork, (ftnlen)5, (ftnlen)5, (ftnlen)12, (ftnlen)4);
#line 435 "slarfb.f"
		if (*n > *k) {

/*                 W := W + C1 * V1 */

#line 439 "slarfb.f"
		    i__1 = *n - *k;
#line 439 "slarfb.f"
		    sgemm_("No transpose", "No transpose", m, k, &i__1, &
			    c_b14, &c__[c_offset], ldc, &v[v_offset], ldv, &
			    c_b14, &work[work_offset], ldwork, (ftnlen)12, (
			    ftnlen)12);
#line 441 "slarfb.f"
		}

/*              W := W * T  or  W * T**T */

#line 445 "slarfb.f"
		strmm_("Right", "Lower", trans, "Non-unit", m, k, &c_b14, &t[
			t_offset], ldt, &work[work_offset], ldwork, (ftnlen)5,
			 (ftnlen)5, (ftnlen)1, (ftnlen)8);

/*              C := C - W * V**T */

#line 450 "slarfb.f"
		if (*n > *k) {

/*                 C1 := C1 - W * V1**T */

#line 454 "slarfb.f"
		    i__1 = *n - *k;
#line 454 "slarfb.f"
		    sgemm_("No transpose", "Transpose", m, &i__1, k, &c_b25, &
			    work[work_offset], ldwork, &v[v_offset], ldv, &
			    c_b14, &c__[c_offset], ldc, (ftnlen)12, (ftnlen)9)
			    ;
#line 456 "slarfb.f"
		}

/*              W := W * V2**T */

#line 460 "slarfb.f"
		strmm_("Right", "Upper", "Transpose", "Unit", m, k, &c_b14, &
			v[*n - *k + 1 + v_dim1], ldv, &work[work_offset], 
			ldwork, (ftnlen)5, (ftnlen)5, (ftnlen)9, (ftnlen)4);

/*              C2 := C2 - W */

#line 465 "slarfb.f"
		i__1 = *k;
#line 465 "slarfb.f"
		for (j = 1; j <= i__1; ++j) {
#line 466 "slarfb.f"
		    i__2 = *m;
#line 466 "slarfb.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 467 "slarfb.f"
			c__[i__ + (*n - *k + j) * c_dim1] -= work[i__ + j * 
				work_dim1];
#line 468 "slarfb.f"
/* L110: */
#line 468 "slarfb.f"
		    }
#line 469 "slarfb.f"
/* L120: */
#line 469 "slarfb.f"
		}
#line 470 "slarfb.f"
	    }
#line 471 "slarfb.f"
	}

#line 473 "slarfb.f"
    } else if (lsame_(storev, "R", (ftnlen)1, (ftnlen)1)) {

#line 475 "slarfb.f"
	if (lsame_(direct, "F", (ftnlen)1, (ftnlen)1)) {

/*           Let  V =  ( V1  V2 )    (V1: first K columns) */
/*           where  V1  is unit upper triangular. */

#line 480 "slarfb.f"
	    if (lsame_(side, "L", (ftnlen)1, (ftnlen)1)) {

/*              Form  H * C  or  H**T * C  where  C = ( C1 ) */
/*                                                    ( C2 ) */

/*              W := C**T * V**T  =  (C1**T * V1**T + C2**T * V2**T) (stored in WORK) */

/*              W := C1**T */

#line 489 "slarfb.f"
		i__1 = *k;
#line 489 "slarfb.f"
		for (j = 1; j <= i__1; ++j) {
#line 490 "slarfb.f"
		    scopy_(n, &c__[j + c_dim1], ldc, &work[j * work_dim1 + 1],
			     &c__1);
#line 491 "slarfb.f"
/* L130: */
#line 491 "slarfb.f"
		}

/*              W := W * V1**T */

#line 495 "slarfb.f"
		strmm_("Right", "Upper", "Transpose", "Unit", n, k, &c_b14, &
			v[v_offset], ldv, &work[work_offset], ldwork, (ftnlen)
			5, (ftnlen)5, (ftnlen)9, (ftnlen)4);
#line 497 "slarfb.f"
		if (*m > *k) {

/*                 W := W + C2**T * V2**T */

#line 501 "slarfb.f"
		    i__1 = *m - *k;
#line 501 "slarfb.f"
		    sgemm_("Transpose", "Transpose", n, k, &i__1, &c_b14, &
			    c__[*k + 1 + c_dim1], ldc, &v[(*k + 1) * v_dim1 + 
			    1], ldv, &c_b14, &work[work_offset], ldwork, (
			    ftnlen)9, (ftnlen)9);
#line 504 "slarfb.f"
		}

/*              W := W * T**T  or  W * T */

#line 508 "slarfb.f"
		strmm_("Right", "Upper", transt, "Non-unit", n, k, &c_b14, &t[
			t_offset], ldt, &work[work_offset], ldwork, (ftnlen)5,
			 (ftnlen)5, (ftnlen)1, (ftnlen)8);

/*              C := C - V**T * W**T */

#line 513 "slarfb.f"
		if (*m > *k) {

/*                 C2 := C2 - V2**T * W**T */

#line 517 "slarfb.f"
		    i__1 = *m - *k;
#line 517 "slarfb.f"
		    sgemm_("Transpose", "Transpose", &i__1, n, k, &c_b25, &v[(
			    *k + 1) * v_dim1 + 1], ldv, &work[work_offset], 
			    ldwork, &c_b14, &c__[*k + 1 + c_dim1], ldc, (
			    ftnlen)9, (ftnlen)9);
#line 520 "slarfb.f"
		}

/*              W := W * V1 */

#line 524 "slarfb.f"
		strmm_("Right", "Upper", "No transpose", "Unit", n, k, &c_b14,
			 &v[v_offset], ldv, &work[work_offset], ldwork, (
			ftnlen)5, (ftnlen)5, (ftnlen)12, (ftnlen)4);

/*              C1 := C1 - W**T */

#line 529 "slarfb.f"
		i__1 = *k;
#line 529 "slarfb.f"
		for (j = 1; j <= i__1; ++j) {
#line 530 "slarfb.f"
		    i__2 = *n;
#line 530 "slarfb.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 531 "slarfb.f"
			c__[j + i__ * c_dim1] -= work[i__ + j * work_dim1];
#line 532 "slarfb.f"
/* L140: */
#line 532 "slarfb.f"
		    }
#line 533 "slarfb.f"
/* L150: */
#line 533 "slarfb.f"
		}

#line 535 "slarfb.f"
	    } else if (lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {

/*              Form  C * H  or  C * H**T  where  C = ( C1  C2 ) */

/*              W := C * V**T  =  (C1*V1**T + C2*V2**T)  (stored in WORK) */

/*              W := C1 */

#line 543 "slarfb.f"
		i__1 = *k;
#line 543 "slarfb.f"
		for (j = 1; j <= i__1; ++j) {
#line 544 "slarfb.f"
		    scopy_(m, &c__[j * c_dim1 + 1], &c__1, &work[j * 
			    work_dim1 + 1], &c__1);
#line 545 "slarfb.f"
/* L160: */
#line 545 "slarfb.f"
		}

/*              W := W * V1**T */

#line 549 "slarfb.f"
		strmm_("Right", "Upper", "Transpose", "Unit", m, k, &c_b14, &
			v[v_offset], ldv, &work[work_offset], ldwork, (ftnlen)
			5, (ftnlen)5, (ftnlen)9, (ftnlen)4);
#line 551 "slarfb.f"
		if (*n > *k) {

/*                 W := W + C2 * V2**T */

#line 555 "slarfb.f"
		    i__1 = *n - *k;
#line 555 "slarfb.f"
		    sgemm_("No transpose", "Transpose", m, k, &i__1, &c_b14, &
			    c__[(*k + 1) * c_dim1 + 1], ldc, &v[(*k + 1) * 
			    v_dim1 + 1], ldv, &c_b14, &work[work_offset], 
			    ldwork, (ftnlen)12, (ftnlen)9);
#line 558 "slarfb.f"
		}

/*              W := W * T  or  W * T**T */

#line 562 "slarfb.f"
		strmm_("Right", "Upper", trans, "Non-unit", m, k, &c_b14, &t[
			t_offset], ldt, &work[work_offset], ldwork, (ftnlen)5,
			 (ftnlen)5, (ftnlen)1, (ftnlen)8);

/*              C := C - W * V */

#line 567 "slarfb.f"
		if (*n > *k) {

/*                 C2 := C2 - W * V2 */

#line 571 "slarfb.f"
		    i__1 = *n - *k;
#line 571 "slarfb.f"
		    sgemm_("No transpose", "No transpose", m, &i__1, k, &
			    c_b25, &work[work_offset], ldwork, &v[(*k + 1) * 
			    v_dim1 + 1], ldv, &c_b14, &c__[(*k + 1) * c_dim1 
			    + 1], ldc, (ftnlen)12, (ftnlen)12);
#line 574 "slarfb.f"
		}

/*              W := W * V1 */

#line 578 "slarfb.f"
		strmm_("Right", "Upper", "No transpose", "Unit", m, k, &c_b14,
			 &v[v_offset], ldv, &work[work_offset], ldwork, (
			ftnlen)5, (ftnlen)5, (ftnlen)12, (ftnlen)4);

/*              C1 := C1 - W */

#line 583 "slarfb.f"
		i__1 = *k;
#line 583 "slarfb.f"
		for (j = 1; j <= i__1; ++j) {
#line 584 "slarfb.f"
		    i__2 = *m;
#line 584 "slarfb.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 585 "slarfb.f"
			c__[i__ + j * c_dim1] -= work[i__ + j * work_dim1];
#line 586 "slarfb.f"
/* L170: */
#line 586 "slarfb.f"
		    }
#line 587 "slarfb.f"
/* L180: */
#line 587 "slarfb.f"
		}

#line 589 "slarfb.f"
	    }

#line 591 "slarfb.f"
	} else {

/*           Let  V =  ( V1  V2 )    (V2: last K columns) */
/*           where  V2  is unit lower triangular. */

#line 596 "slarfb.f"
	    if (lsame_(side, "L", (ftnlen)1, (ftnlen)1)) {

/*              Form  H * C  or  H**T * C  where  C = ( C1 ) */
/*                                                    ( C2 ) */

/*              W := C**T * V**T  =  (C1**T * V1**T + C2**T * V2**T) (stored in WORK) */

/*              W := C2**T */

#line 605 "slarfb.f"
		i__1 = *k;
#line 605 "slarfb.f"
		for (j = 1; j <= i__1; ++j) {
#line 606 "slarfb.f"
		    scopy_(n, &c__[*m - *k + j + c_dim1], ldc, &work[j * 
			    work_dim1 + 1], &c__1);
#line 607 "slarfb.f"
/* L190: */
#line 607 "slarfb.f"
		}

/*              W := W * V2**T */

#line 611 "slarfb.f"
		strmm_("Right", "Lower", "Transpose", "Unit", n, k, &c_b14, &
			v[(*m - *k + 1) * v_dim1 + 1], ldv, &work[work_offset]
			, ldwork, (ftnlen)5, (ftnlen)5, (ftnlen)9, (ftnlen)4);
#line 613 "slarfb.f"
		if (*m > *k) {

/*                 W := W + C1**T * V1**T */

#line 617 "slarfb.f"
		    i__1 = *m - *k;
#line 617 "slarfb.f"
		    sgemm_("Transpose", "Transpose", n, k, &i__1, &c_b14, &
			    c__[c_offset], ldc, &v[v_offset], ldv, &c_b14, &
			    work[work_offset], ldwork, (ftnlen)9, (ftnlen)9);
#line 619 "slarfb.f"
		}

/*              W := W * T**T  or  W * T */

#line 623 "slarfb.f"
		strmm_("Right", "Lower", transt, "Non-unit", n, k, &c_b14, &t[
			t_offset], ldt, &work[work_offset], ldwork, (ftnlen)5,
			 (ftnlen)5, (ftnlen)1, (ftnlen)8);

/*              C := C - V**T * W**T */

#line 628 "slarfb.f"
		if (*m > *k) {

/*                 C1 := C1 - V1**T * W**T */

#line 632 "slarfb.f"
		    i__1 = *m - *k;
#line 632 "slarfb.f"
		    sgemm_("Transpose", "Transpose", &i__1, n, k, &c_b25, &v[
			    v_offset], ldv, &work[work_offset], ldwork, &
			    c_b14, &c__[c_offset], ldc, (ftnlen)9, (ftnlen)9);
#line 634 "slarfb.f"
		}

/*              W := W * V2 */

#line 638 "slarfb.f"
		strmm_("Right", "Lower", "No transpose", "Unit", n, k, &c_b14,
			 &v[(*m - *k + 1) * v_dim1 + 1], ldv, &work[
			work_offset], ldwork, (ftnlen)5, (ftnlen)5, (ftnlen)
			12, (ftnlen)4);

/*              C2 := C2 - W**T */

#line 643 "slarfb.f"
		i__1 = *k;
#line 643 "slarfb.f"
		for (j = 1; j <= i__1; ++j) {
#line 644 "slarfb.f"
		    i__2 = *n;
#line 644 "slarfb.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 645 "slarfb.f"
			c__[*m - *k + j + i__ * c_dim1] -= work[i__ + j * 
				work_dim1];
#line 646 "slarfb.f"
/* L200: */
#line 646 "slarfb.f"
		    }
#line 647 "slarfb.f"
/* L210: */
#line 647 "slarfb.f"
		}

#line 649 "slarfb.f"
	    } else if (lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {

/*              Form  C * H  or  C * H**T  where  C = ( C1  C2 ) */

/*              W := C * V**T  =  (C1*V1**T + C2*V2**T)  (stored in WORK) */

/*              W := C2 */

#line 657 "slarfb.f"
		i__1 = *k;
#line 657 "slarfb.f"
		for (j = 1; j <= i__1; ++j) {
#line 658 "slarfb.f"
		    scopy_(m, &c__[(*n - *k + j) * c_dim1 + 1], &c__1, &work[
			    j * work_dim1 + 1], &c__1);
#line 659 "slarfb.f"
/* L220: */
#line 659 "slarfb.f"
		}

/*              W := W * V2**T */

#line 663 "slarfb.f"
		strmm_("Right", "Lower", "Transpose", "Unit", m, k, &c_b14, &
			v[(*n - *k + 1) * v_dim1 + 1], ldv, &work[work_offset]
			, ldwork, (ftnlen)5, (ftnlen)5, (ftnlen)9, (ftnlen)4);
#line 665 "slarfb.f"
		if (*n > *k) {

/*                 W := W + C1 * V1**T */

#line 669 "slarfb.f"
		    i__1 = *n - *k;
#line 669 "slarfb.f"
		    sgemm_("No transpose", "Transpose", m, k, &i__1, &c_b14, &
			    c__[c_offset], ldc, &v[v_offset], ldv, &c_b14, &
			    work[work_offset], ldwork, (ftnlen)12, (ftnlen)9);
#line 671 "slarfb.f"
		}

/*              W := W * T  or  W * T**T */

#line 675 "slarfb.f"
		strmm_("Right", "Lower", trans, "Non-unit", m, k, &c_b14, &t[
			t_offset], ldt, &work[work_offset], ldwork, (ftnlen)5,
			 (ftnlen)5, (ftnlen)1, (ftnlen)8);

/*              C := C - W * V */

#line 680 "slarfb.f"
		if (*n > *k) {

/*                 C1 := C1 - W * V1 */

#line 684 "slarfb.f"
		    i__1 = *n - *k;
#line 684 "slarfb.f"
		    sgemm_("No transpose", "No transpose", m, &i__1, k, &
			    c_b25, &work[work_offset], ldwork, &v[v_offset], 
			    ldv, &c_b14, &c__[c_offset], ldc, (ftnlen)12, (
			    ftnlen)12);
#line 686 "slarfb.f"
		}

/*              W := W * V2 */

#line 690 "slarfb.f"
		strmm_("Right", "Lower", "No transpose", "Unit", m, k, &c_b14,
			 &v[(*n - *k + 1) * v_dim1 + 1], ldv, &work[
			work_offset], ldwork, (ftnlen)5, (ftnlen)5, (ftnlen)
			12, (ftnlen)4);

/*              C1 := C1 - W */

#line 695 "slarfb.f"
		i__1 = *k;
#line 695 "slarfb.f"
		for (j = 1; j <= i__1; ++j) {
#line 696 "slarfb.f"
		    i__2 = *m;
#line 696 "slarfb.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 697 "slarfb.f"
			c__[i__ + (*n - *k + j) * c_dim1] -= work[i__ + j * 
				work_dim1];
#line 698 "slarfb.f"
/* L230: */
#line 698 "slarfb.f"
		    }
#line 699 "slarfb.f"
/* L240: */
#line 699 "slarfb.f"
		}

#line 701 "slarfb.f"
	    }

#line 703 "slarfb.f"
	}
#line 704 "slarfb.f"
    }

#line 706 "slarfb.f"
    return 0;

/*     End of SLARFB */

} /* slarfb_ */

