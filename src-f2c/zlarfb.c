#line 1 "zlarfb.f"
/* zlarfb.f -- translated by f2c (version 20100827).
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

#line 1 "zlarfb.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static integer c__1 = 1;

/* > \brief \b ZLARFB applies a block reflector or its conjugate-transpose to a general rectangular matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZLARFB + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlarfb.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlarfb.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlarfb.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZLARFB( SIDE, TRANS, DIRECT, STOREV, M, N, K, V, LDV, */
/*                          T, LDT, C, LDC, WORK, LDWORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          DIRECT, SIDE, STOREV, TRANS */
/*       INTEGER            K, LDC, LDT, LDV, LDWORK, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         C( LDC, * ), T( LDT, * ), V( LDV, * ), */
/*      $                   WORK( LDWORK, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLARFB applies a complex block reflector H or its transpose H**H to a */
/* > complex M-by-N matrix C, from either the left or the right. */
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
/* >          V is COMPLEX*16 array, dimension */
/* >                                (LDV,K) if STOREV = 'C' */
/* >                                (LDV,M) if STOREV = 'R' and SIDE = 'L' */
/* >                                (LDV,N) if STOREV = 'R' and SIDE = 'R' */
/* >          See Further Details. */
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
/* >          The leading dimension of the array T. LDT >= K. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* >          C is COMPLEX*16 array, dimension (LDC,N) */
/* >          On entry, the M-by-N matrix C. */
/* >          On exit, C is overwritten by H*C or H**H*C or C*H or C*H**H. */
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
/* >          WORK is COMPLEX*16 array, dimension (LDWORK,K) */
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

/* > \ingroup complex16OTHERauxiliary */

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
/* Subroutine */ int zlarfb_(char *side, char *trans, char *direct, char *
	storev, integer *m, integer *n, integer *k, doublecomplex *v, integer 
	*ldv, doublecomplex *t, integer *ldt, doublecomplex *c__, integer *
	ldc, doublecomplex *work, integer *ldwork, ftnlen side_len, ftnlen 
	trans_len, ftnlen direct_len, ftnlen storev_len)
{
    /* System generated locals */
    integer c_dim1, c_offset, t_dim1, t_offset, v_dim1, v_offset, work_dim1, 
	    work_offset, i__1, i__2, i__3, i__4, i__5;
    doublecomplex z__1, z__2;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, j;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int zgemm_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, ftnlen, ftnlen), zcopy_(integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *), ztrmm_(char *, char *, 
	    char *, char *, integer *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen, 
	    ftnlen, ftnlen, ftnlen), zlacgv_(integer *, doublecomplex *, 
	    integer *);
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
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Quick return if possible */

#line 236 "zlarfb.f"
    /* Parameter adjustments */
#line 236 "zlarfb.f"
    v_dim1 = *ldv;
#line 236 "zlarfb.f"
    v_offset = 1 + v_dim1;
#line 236 "zlarfb.f"
    v -= v_offset;
#line 236 "zlarfb.f"
    t_dim1 = *ldt;
#line 236 "zlarfb.f"
    t_offset = 1 + t_dim1;
#line 236 "zlarfb.f"
    t -= t_offset;
#line 236 "zlarfb.f"
    c_dim1 = *ldc;
#line 236 "zlarfb.f"
    c_offset = 1 + c_dim1;
#line 236 "zlarfb.f"
    c__ -= c_offset;
#line 236 "zlarfb.f"
    work_dim1 = *ldwork;
#line 236 "zlarfb.f"
    work_offset = 1 + work_dim1;
#line 236 "zlarfb.f"
    work -= work_offset;
#line 236 "zlarfb.f"

#line 236 "zlarfb.f"
    /* Function Body */
#line 236 "zlarfb.f"
    if (*m <= 0 || *n <= 0) {
#line 236 "zlarfb.f"
	return 0;
#line 236 "zlarfb.f"
    }

#line 239 "zlarfb.f"
    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {
#line 240 "zlarfb.f"
	*(unsigned char *)transt = 'C';
#line 241 "zlarfb.f"
    } else {
#line 242 "zlarfb.f"
	*(unsigned char *)transt = 'N';
#line 243 "zlarfb.f"
    }

#line 245 "zlarfb.f"
    if (lsame_(storev, "C", (ftnlen)1, (ftnlen)1)) {

#line 247 "zlarfb.f"
	if (lsame_(direct, "F", (ftnlen)1, (ftnlen)1)) {

/*           Let  V =  ( V1 )    (first K rows) */
/*                     ( V2 ) */
/*           where  V1  is unit lower triangular. */

#line 253 "zlarfb.f"
	    if (lsame_(side, "L", (ftnlen)1, (ftnlen)1)) {

/*              Form  H * C  or  H**H * C  where  C = ( C1 ) */
/*                                                    ( C2 ) */

/*              W := C**H * V  =  (C1**H * V1 + C2**H * V2)  (stored in WORK) */

/*              W := C1**H */

#line 262 "zlarfb.f"
		i__1 = *k;
#line 262 "zlarfb.f"
		for (j = 1; j <= i__1; ++j) {
#line 263 "zlarfb.f"
		    zcopy_(n, &c__[j + c_dim1], ldc, &work[j * work_dim1 + 1],
			     &c__1);
#line 264 "zlarfb.f"
		    zlacgv_(n, &work[j * work_dim1 + 1], &c__1);
#line 265 "zlarfb.f"
/* L10: */
#line 265 "zlarfb.f"
		}

/*              W := W * V1 */

#line 269 "zlarfb.f"
		ztrmm_("Right", "Lower", "No transpose", "Unit", n, k, &c_b1, 
			&v[v_offset], ldv, &work[work_offset], ldwork, (
			ftnlen)5, (ftnlen)5, (ftnlen)12, (ftnlen)4);
#line 271 "zlarfb.f"
		if (*m > *k) {

/*                 W := W + C2**H * V2 */

#line 275 "zlarfb.f"
		    i__1 = *m - *k;
#line 275 "zlarfb.f"
		    zgemm_("Conjugate transpose", "No transpose", n, k, &i__1,
			     &c_b1, &c__[*k + 1 + c_dim1], ldc, &v[*k + 1 + 
			    v_dim1], ldv, &c_b1, &work[work_offset], ldwork, (
			    ftnlen)19, (ftnlen)12);
#line 278 "zlarfb.f"
		}

/*              W := W * T**H  or  W * T */

#line 282 "zlarfb.f"
		ztrmm_("Right", "Upper", transt, "Non-unit", n, k, &c_b1, &t[
			t_offset], ldt, &work[work_offset], ldwork, (ftnlen)5,
			 (ftnlen)5, (ftnlen)1, (ftnlen)8);

/*              C := C - V * W**H */

#line 287 "zlarfb.f"
		if (*m > *k) {

/*                 C2 := C2 - V2 * W**H */

#line 291 "zlarfb.f"
		    i__1 = *m - *k;
#line 291 "zlarfb.f"
		    z__1.r = -1., z__1.i = -0.;
#line 291 "zlarfb.f"
		    zgemm_("No transpose", "Conjugate transpose", &i__1, n, k,
			     &z__1, &v[*k + 1 + v_dim1], ldv, &work[
			    work_offset], ldwork, &c_b1, &c__[*k + 1 + c_dim1]
			    , ldc, (ftnlen)12, (ftnlen)19);
#line 294 "zlarfb.f"
		}

/*              W := W * V1**H */

#line 298 "zlarfb.f"
		ztrmm_("Right", "Lower", "Conjugate transpose", "Unit", n, k, 
			&c_b1, &v[v_offset], ldv, &work[work_offset], ldwork, 
			(ftnlen)5, (ftnlen)5, (ftnlen)19, (ftnlen)4);

/*              C1 := C1 - W**H */

#line 303 "zlarfb.f"
		i__1 = *k;
#line 303 "zlarfb.f"
		for (j = 1; j <= i__1; ++j) {
#line 304 "zlarfb.f"
		    i__2 = *n;
#line 304 "zlarfb.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 305 "zlarfb.f"
			i__3 = j + i__ * c_dim1;
#line 305 "zlarfb.f"
			i__4 = j + i__ * c_dim1;
#line 305 "zlarfb.f"
			d_cnjg(&z__2, &work[i__ + j * work_dim1]);
#line 305 "zlarfb.f"
			z__1.r = c__[i__4].r - z__2.r, z__1.i = c__[i__4].i - 
				z__2.i;
#line 305 "zlarfb.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 306 "zlarfb.f"
/* L20: */
#line 306 "zlarfb.f"
		    }
#line 307 "zlarfb.f"
/* L30: */
#line 307 "zlarfb.f"
		}

#line 309 "zlarfb.f"
	    } else if (lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {

/*              Form  C * H  or  C * H**H  where  C = ( C1  C2 ) */

/*              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK) */

/*              W := C1 */

#line 317 "zlarfb.f"
		i__1 = *k;
#line 317 "zlarfb.f"
		for (j = 1; j <= i__1; ++j) {
#line 318 "zlarfb.f"
		    zcopy_(m, &c__[j * c_dim1 + 1], &c__1, &work[j * 
			    work_dim1 + 1], &c__1);
#line 319 "zlarfb.f"
/* L40: */
#line 319 "zlarfb.f"
		}

/*              W := W * V1 */

#line 323 "zlarfb.f"
		ztrmm_("Right", "Lower", "No transpose", "Unit", m, k, &c_b1, 
			&v[v_offset], ldv, &work[work_offset], ldwork, (
			ftnlen)5, (ftnlen)5, (ftnlen)12, (ftnlen)4);
#line 325 "zlarfb.f"
		if (*n > *k) {

/*                 W := W + C2 * V2 */

#line 329 "zlarfb.f"
		    i__1 = *n - *k;
#line 329 "zlarfb.f"
		    zgemm_("No transpose", "No transpose", m, k, &i__1, &c_b1,
			     &c__[(*k + 1) * c_dim1 + 1], ldc, &v[*k + 1 + 
			    v_dim1], ldv, &c_b1, &work[work_offset], ldwork, (
			    ftnlen)12, (ftnlen)12);
#line 332 "zlarfb.f"
		}

/*              W := W * T  or  W * T**H */

#line 336 "zlarfb.f"
		ztrmm_("Right", "Upper", trans, "Non-unit", m, k, &c_b1, &t[
			t_offset], ldt, &work[work_offset], ldwork, (ftnlen)5,
			 (ftnlen)5, (ftnlen)1, (ftnlen)8);

/*              C := C - W * V**H */

#line 341 "zlarfb.f"
		if (*n > *k) {

/*                 C2 := C2 - W * V2**H */

#line 345 "zlarfb.f"
		    i__1 = *n - *k;
#line 345 "zlarfb.f"
		    z__1.r = -1., z__1.i = -0.;
#line 345 "zlarfb.f"
		    zgemm_("No transpose", "Conjugate transpose", m, &i__1, k,
			     &z__1, &work[work_offset], ldwork, &v[*k + 1 + 
			    v_dim1], ldv, &c_b1, &c__[(*k + 1) * c_dim1 + 1], 
			    ldc, (ftnlen)12, (ftnlen)19);
#line 348 "zlarfb.f"
		}

/*              W := W * V1**H */

#line 352 "zlarfb.f"
		ztrmm_("Right", "Lower", "Conjugate transpose", "Unit", m, k, 
			&c_b1, &v[v_offset], ldv, &work[work_offset], ldwork, 
			(ftnlen)5, (ftnlen)5, (ftnlen)19, (ftnlen)4);

/*              C1 := C1 - W */

#line 357 "zlarfb.f"
		i__1 = *k;
#line 357 "zlarfb.f"
		for (j = 1; j <= i__1; ++j) {
#line 358 "zlarfb.f"
		    i__2 = *m;
#line 358 "zlarfb.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 359 "zlarfb.f"
			i__3 = i__ + j * c_dim1;
#line 359 "zlarfb.f"
			i__4 = i__ + j * c_dim1;
#line 359 "zlarfb.f"
			i__5 = i__ + j * work_dim1;
#line 359 "zlarfb.f"
			z__1.r = c__[i__4].r - work[i__5].r, z__1.i = c__[
				i__4].i - work[i__5].i;
#line 359 "zlarfb.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 360 "zlarfb.f"
/* L50: */
#line 360 "zlarfb.f"
		    }
#line 361 "zlarfb.f"
/* L60: */
#line 361 "zlarfb.f"
		}
#line 362 "zlarfb.f"
	    }

#line 364 "zlarfb.f"
	} else {

/*           Let  V =  ( V1 ) */
/*                     ( V2 )    (last K rows) */
/*           where  V2  is unit upper triangular. */

#line 370 "zlarfb.f"
	    if (lsame_(side, "L", (ftnlen)1, (ftnlen)1)) {

/*              Form  H * C  or  H**H * C  where  C = ( C1 ) */
/*                                                    ( C2 ) */

/*              W := C**H * V  =  (C1**H * V1 + C2**H * V2)  (stored in WORK) */

/*              W := C2**H */

#line 379 "zlarfb.f"
		i__1 = *k;
#line 379 "zlarfb.f"
		for (j = 1; j <= i__1; ++j) {
#line 380 "zlarfb.f"
		    zcopy_(n, &c__[*m - *k + j + c_dim1], ldc, &work[j * 
			    work_dim1 + 1], &c__1);
#line 381 "zlarfb.f"
		    zlacgv_(n, &work[j * work_dim1 + 1], &c__1);
#line 382 "zlarfb.f"
/* L70: */
#line 382 "zlarfb.f"
		}

/*              W := W * V2 */

#line 386 "zlarfb.f"
		ztrmm_("Right", "Upper", "No transpose", "Unit", n, k, &c_b1, 
			&v[*m - *k + 1 + v_dim1], ldv, &work[work_offset], 
			ldwork, (ftnlen)5, (ftnlen)5, (ftnlen)12, (ftnlen)4);
#line 388 "zlarfb.f"
		if (*m > *k) {

/*                 W := W + C1**H * V1 */

#line 392 "zlarfb.f"
		    i__1 = *m - *k;
#line 392 "zlarfb.f"
		    zgemm_("Conjugate transpose", "No transpose", n, k, &i__1,
			     &c_b1, &c__[c_offset], ldc, &v[v_offset], ldv, &
			    c_b1, &work[work_offset], ldwork, (ftnlen)19, (
			    ftnlen)12);
#line 395 "zlarfb.f"
		}

/*              W := W * T**H  or  W * T */

#line 399 "zlarfb.f"
		ztrmm_("Right", "Lower", transt, "Non-unit", n, k, &c_b1, &t[
			t_offset], ldt, &work[work_offset], ldwork, (ftnlen)5,
			 (ftnlen)5, (ftnlen)1, (ftnlen)8);

/*              C := C - V * W**H */

#line 404 "zlarfb.f"
		if (*m > *k) {

/*                 C1 := C1 - V1 * W**H */

#line 408 "zlarfb.f"
		    i__1 = *m - *k;
#line 408 "zlarfb.f"
		    z__1.r = -1., z__1.i = -0.;
#line 408 "zlarfb.f"
		    zgemm_("No transpose", "Conjugate transpose", &i__1, n, k,
			     &z__1, &v[v_offset], ldv, &work[work_offset], 
			    ldwork, &c_b1, &c__[c_offset], ldc, (ftnlen)12, (
			    ftnlen)19);
#line 411 "zlarfb.f"
		}

/*              W := W * V2**H */

#line 415 "zlarfb.f"
		ztrmm_("Right", "Upper", "Conjugate transpose", "Unit", n, k, 
			&c_b1, &v[*m - *k + 1 + v_dim1], ldv, &work[
			work_offset], ldwork, (ftnlen)5, (ftnlen)5, (ftnlen)
			19, (ftnlen)4);

/*              C2 := C2 - W**H */

#line 421 "zlarfb.f"
		i__1 = *k;
#line 421 "zlarfb.f"
		for (j = 1; j <= i__1; ++j) {
#line 422 "zlarfb.f"
		    i__2 = *n;
#line 422 "zlarfb.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 423 "zlarfb.f"
			i__3 = *m - *k + j + i__ * c_dim1;
#line 423 "zlarfb.f"
			i__4 = *m - *k + j + i__ * c_dim1;
#line 423 "zlarfb.f"
			d_cnjg(&z__2, &work[i__ + j * work_dim1]);
#line 423 "zlarfb.f"
			z__1.r = c__[i__4].r - z__2.r, z__1.i = c__[i__4].i - 
				z__2.i;
#line 423 "zlarfb.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 425 "zlarfb.f"
/* L80: */
#line 425 "zlarfb.f"
		    }
#line 426 "zlarfb.f"
/* L90: */
#line 426 "zlarfb.f"
		}

#line 428 "zlarfb.f"
	    } else if (lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {

/*              Form  C * H  or  C * H**H  where  C = ( C1  C2 ) */

/*              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK) */

/*              W := C2 */

#line 436 "zlarfb.f"
		i__1 = *k;
#line 436 "zlarfb.f"
		for (j = 1; j <= i__1; ++j) {
#line 437 "zlarfb.f"
		    zcopy_(m, &c__[(*n - *k + j) * c_dim1 + 1], &c__1, &work[
			    j * work_dim1 + 1], &c__1);
#line 438 "zlarfb.f"
/* L100: */
#line 438 "zlarfb.f"
		}

/*              W := W * V2 */

#line 442 "zlarfb.f"
		ztrmm_("Right", "Upper", "No transpose", "Unit", m, k, &c_b1, 
			&v[*n - *k + 1 + v_dim1], ldv, &work[work_offset], 
			ldwork, (ftnlen)5, (ftnlen)5, (ftnlen)12, (ftnlen)4);
#line 444 "zlarfb.f"
		if (*n > *k) {

/*                 W := W + C1 * V1 */

#line 448 "zlarfb.f"
		    i__1 = *n - *k;
#line 448 "zlarfb.f"
		    zgemm_("No transpose", "No transpose", m, k, &i__1, &c_b1,
			     &c__[c_offset], ldc, &v[v_offset], ldv, &c_b1, &
			    work[work_offset], ldwork, (ftnlen)12, (ftnlen)12)
			    ;
#line 450 "zlarfb.f"
		}

/*              W := W * T  or  W * T**H */

#line 454 "zlarfb.f"
		ztrmm_("Right", "Lower", trans, "Non-unit", m, k, &c_b1, &t[
			t_offset], ldt, &work[work_offset], ldwork, (ftnlen)5,
			 (ftnlen)5, (ftnlen)1, (ftnlen)8);

/*              C := C - W * V**H */

#line 459 "zlarfb.f"
		if (*n > *k) {

/*                 C1 := C1 - W * V1**H */

#line 463 "zlarfb.f"
		    i__1 = *n - *k;
#line 463 "zlarfb.f"
		    z__1.r = -1., z__1.i = -0.;
#line 463 "zlarfb.f"
		    zgemm_("No transpose", "Conjugate transpose", m, &i__1, k,
			     &z__1, &work[work_offset], ldwork, &v[v_offset], 
			    ldv, &c_b1, &c__[c_offset], ldc, (ftnlen)12, (
			    ftnlen)19);
#line 466 "zlarfb.f"
		}

/*              W := W * V2**H */

#line 470 "zlarfb.f"
		ztrmm_("Right", "Upper", "Conjugate transpose", "Unit", m, k, 
			&c_b1, &v[*n - *k + 1 + v_dim1], ldv, &work[
			work_offset], ldwork, (ftnlen)5, (ftnlen)5, (ftnlen)
			19, (ftnlen)4);

/*              C2 := C2 - W */

#line 476 "zlarfb.f"
		i__1 = *k;
#line 476 "zlarfb.f"
		for (j = 1; j <= i__1; ++j) {
#line 477 "zlarfb.f"
		    i__2 = *m;
#line 477 "zlarfb.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 478 "zlarfb.f"
			i__3 = i__ + (*n - *k + j) * c_dim1;
#line 478 "zlarfb.f"
			i__4 = i__ + (*n - *k + j) * c_dim1;
#line 478 "zlarfb.f"
			i__5 = i__ + j * work_dim1;
#line 478 "zlarfb.f"
			z__1.r = c__[i__4].r - work[i__5].r, z__1.i = c__[
				i__4].i - work[i__5].i;
#line 478 "zlarfb.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 479 "zlarfb.f"
/* L110: */
#line 479 "zlarfb.f"
		    }
#line 480 "zlarfb.f"
/* L120: */
#line 480 "zlarfb.f"
		}
#line 481 "zlarfb.f"
	    }
#line 482 "zlarfb.f"
	}

#line 484 "zlarfb.f"
    } else if (lsame_(storev, "R", (ftnlen)1, (ftnlen)1)) {

#line 486 "zlarfb.f"
	if (lsame_(direct, "F", (ftnlen)1, (ftnlen)1)) {

/*           Let  V =  ( V1  V2 )    (V1: first K columns) */
/*           where  V1  is unit upper triangular. */

#line 491 "zlarfb.f"
	    if (lsame_(side, "L", (ftnlen)1, (ftnlen)1)) {

/*              Form  H * C  or  H**H * C  where  C = ( C1 ) */
/*                                                    ( C2 ) */

/*              W := C**H * V**H  =  (C1**H * V1**H + C2**H * V2**H) (stored in WORK) */

/*              W := C1**H */

#line 500 "zlarfb.f"
		i__1 = *k;
#line 500 "zlarfb.f"
		for (j = 1; j <= i__1; ++j) {
#line 501 "zlarfb.f"
		    zcopy_(n, &c__[j + c_dim1], ldc, &work[j * work_dim1 + 1],
			     &c__1);
#line 502 "zlarfb.f"
		    zlacgv_(n, &work[j * work_dim1 + 1], &c__1);
#line 503 "zlarfb.f"
/* L130: */
#line 503 "zlarfb.f"
		}

/*              W := W * V1**H */

#line 507 "zlarfb.f"
		ztrmm_("Right", "Upper", "Conjugate transpose", "Unit", n, k, 
			&c_b1, &v[v_offset], ldv, &work[work_offset], ldwork, 
			(ftnlen)5, (ftnlen)5, (ftnlen)19, (ftnlen)4);
#line 509 "zlarfb.f"
		if (*m > *k) {

/*                 W := W + C2**H * V2**H */

#line 513 "zlarfb.f"
		    i__1 = *m - *k;
#line 513 "zlarfb.f"
		    zgemm_("Conjugate transpose", "Conjugate transpose", n, k,
			     &i__1, &c_b1, &c__[*k + 1 + c_dim1], ldc, &v[(*k 
			    + 1) * v_dim1 + 1], ldv, &c_b1, &work[work_offset]
			    , ldwork, (ftnlen)19, (ftnlen)19);
#line 517 "zlarfb.f"
		}

/*              W := W * T**H  or  W * T */

#line 521 "zlarfb.f"
		ztrmm_("Right", "Upper", transt, "Non-unit", n, k, &c_b1, &t[
			t_offset], ldt, &work[work_offset], ldwork, (ftnlen)5,
			 (ftnlen)5, (ftnlen)1, (ftnlen)8);

/*              C := C - V**H * W**H */

#line 526 "zlarfb.f"
		if (*m > *k) {

/*                 C2 := C2 - V2**H * W**H */

#line 530 "zlarfb.f"
		    i__1 = *m - *k;
#line 530 "zlarfb.f"
		    z__1.r = -1., z__1.i = -0.;
#line 530 "zlarfb.f"
		    zgemm_("Conjugate transpose", "Conjugate transpose", &
			    i__1, n, k, &z__1, &v[(*k + 1) * v_dim1 + 1], ldv,
			     &work[work_offset], ldwork, &c_b1, &c__[*k + 1 + 
			    c_dim1], ldc, (ftnlen)19, (ftnlen)19);
#line 534 "zlarfb.f"
		}

/*              W := W * V1 */

#line 538 "zlarfb.f"
		ztrmm_("Right", "Upper", "No transpose", "Unit", n, k, &c_b1, 
			&v[v_offset], ldv, &work[work_offset], ldwork, (
			ftnlen)5, (ftnlen)5, (ftnlen)12, (ftnlen)4);

/*              C1 := C1 - W**H */

#line 543 "zlarfb.f"
		i__1 = *k;
#line 543 "zlarfb.f"
		for (j = 1; j <= i__1; ++j) {
#line 544 "zlarfb.f"
		    i__2 = *n;
#line 544 "zlarfb.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 545 "zlarfb.f"
			i__3 = j + i__ * c_dim1;
#line 545 "zlarfb.f"
			i__4 = j + i__ * c_dim1;
#line 545 "zlarfb.f"
			d_cnjg(&z__2, &work[i__ + j * work_dim1]);
#line 545 "zlarfb.f"
			z__1.r = c__[i__4].r - z__2.r, z__1.i = c__[i__4].i - 
				z__2.i;
#line 545 "zlarfb.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 546 "zlarfb.f"
/* L140: */
#line 546 "zlarfb.f"
		    }
#line 547 "zlarfb.f"
/* L150: */
#line 547 "zlarfb.f"
		}

#line 549 "zlarfb.f"
	    } else if (lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {

/*              Form  C * H  or  C * H**H  where  C = ( C1  C2 ) */

/*              W := C * V**H  =  (C1*V1**H + C2*V2**H)  (stored in WORK) */

/*              W := C1 */

#line 557 "zlarfb.f"
		i__1 = *k;
#line 557 "zlarfb.f"
		for (j = 1; j <= i__1; ++j) {
#line 558 "zlarfb.f"
		    zcopy_(m, &c__[j * c_dim1 + 1], &c__1, &work[j * 
			    work_dim1 + 1], &c__1);
#line 559 "zlarfb.f"
/* L160: */
#line 559 "zlarfb.f"
		}

/*              W := W * V1**H */

#line 563 "zlarfb.f"
		ztrmm_("Right", "Upper", "Conjugate transpose", "Unit", m, k, 
			&c_b1, &v[v_offset], ldv, &work[work_offset], ldwork, 
			(ftnlen)5, (ftnlen)5, (ftnlen)19, (ftnlen)4);
#line 565 "zlarfb.f"
		if (*n > *k) {

/*                 W := W + C2 * V2**H */

#line 569 "zlarfb.f"
		    i__1 = *n - *k;
#line 569 "zlarfb.f"
		    zgemm_("No transpose", "Conjugate transpose", m, k, &i__1,
			     &c_b1, &c__[(*k + 1) * c_dim1 + 1], ldc, &v[(*k 
			    + 1) * v_dim1 + 1], ldv, &c_b1, &work[work_offset]
			    , ldwork, (ftnlen)12, (ftnlen)19);
#line 572 "zlarfb.f"
		}

/*              W := W * T  or  W * T**H */

#line 576 "zlarfb.f"
		ztrmm_("Right", "Upper", trans, "Non-unit", m, k, &c_b1, &t[
			t_offset], ldt, &work[work_offset], ldwork, (ftnlen)5,
			 (ftnlen)5, (ftnlen)1, (ftnlen)8);

/*              C := C - W * V */

#line 581 "zlarfb.f"
		if (*n > *k) {

/*                 C2 := C2 - W * V2 */

#line 585 "zlarfb.f"
		    i__1 = *n - *k;
#line 585 "zlarfb.f"
		    z__1.r = -1., z__1.i = -0.;
#line 585 "zlarfb.f"
		    zgemm_("No transpose", "No transpose", m, &i__1, k, &z__1,
			     &work[work_offset], ldwork, &v[(*k + 1) * v_dim1 
			    + 1], ldv, &c_b1, &c__[(*k + 1) * c_dim1 + 1], 
			    ldc, (ftnlen)12, (ftnlen)12);
#line 588 "zlarfb.f"
		}

/*              W := W * V1 */

#line 592 "zlarfb.f"
		ztrmm_("Right", "Upper", "No transpose", "Unit", m, k, &c_b1, 
			&v[v_offset], ldv, &work[work_offset], ldwork, (
			ftnlen)5, (ftnlen)5, (ftnlen)12, (ftnlen)4);

/*              C1 := C1 - W */

#line 597 "zlarfb.f"
		i__1 = *k;
#line 597 "zlarfb.f"
		for (j = 1; j <= i__1; ++j) {
#line 598 "zlarfb.f"
		    i__2 = *m;
#line 598 "zlarfb.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 599 "zlarfb.f"
			i__3 = i__ + j * c_dim1;
#line 599 "zlarfb.f"
			i__4 = i__ + j * c_dim1;
#line 599 "zlarfb.f"
			i__5 = i__ + j * work_dim1;
#line 599 "zlarfb.f"
			z__1.r = c__[i__4].r - work[i__5].r, z__1.i = c__[
				i__4].i - work[i__5].i;
#line 599 "zlarfb.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 600 "zlarfb.f"
/* L170: */
#line 600 "zlarfb.f"
		    }
#line 601 "zlarfb.f"
/* L180: */
#line 601 "zlarfb.f"
		}

#line 603 "zlarfb.f"
	    }

#line 605 "zlarfb.f"
	} else {

/*           Let  V =  ( V1  V2 )    (V2: last K columns) */
/*           where  V2  is unit lower triangular. */

#line 610 "zlarfb.f"
	    if (lsame_(side, "L", (ftnlen)1, (ftnlen)1)) {

/*              Form  H * C  or  H**H * C  where  C = ( C1 ) */
/*                                                    ( C2 ) */

/*              W := C**H * V**H  =  (C1**H * V1**H + C2**H * V2**H) (stored in WORK) */

/*              W := C2**H */

#line 619 "zlarfb.f"
		i__1 = *k;
#line 619 "zlarfb.f"
		for (j = 1; j <= i__1; ++j) {
#line 620 "zlarfb.f"
		    zcopy_(n, &c__[*m - *k + j + c_dim1], ldc, &work[j * 
			    work_dim1 + 1], &c__1);
#line 621 "zlarfb.f"
		    zlacgv_(n, &work[j * work_dim1 + 1], &c__1);
#line 622 "zlarfb.f"
/* L190: */
#line 622 "zlarfb.f"
		}

/*              W := W * V2**H */

#line 626 "zlarfb.f"
		ztrmm_("Right", "Lower", "Conjugate transpose", "Unit", n, k, 
			&c_b1, &v[(*m - *k + 1) * v_dim1 + 1], ldv, &work[
			work_offset], ldwork, (ftnlen)5, (ftnlen)5, (ftnlen)
			19, (ftnlen)4);
#line 629 "zlarfb.f"
		if (*m > *k) {

/*                 W := W + C1**H * V1**H */

#line 633 "zlarfb.f"
		    i__1 = *m - *k;
#line 633 "zlarfb.f"
		    zgemm_("Conjugate transpose", "Conjugate transpose", n, k,
			     &i__1, &c_b1, &c__[c_offset], ldc, &v[v_offset], 
			    ldv, &c_b1, &work[work_offset], ldwork, (ftnlen)
			    19, (ftnlen)19);
#line 636 "zlarfb.f"
		}

/*              W := W * T**H  or  W * T */

#line 640 "zlarfb.f"
		ztrmm_("Right", "Lower", transt, "Non-unit", n, k, &c_b1, &t[
			t_offset], ldt, &work[work_offset], ldwork, (ftnlen)5,
			 (ftnlen)5, (ftnlen)1, (ftnlen)8);

/*              C := C - V**H * W**H */

#line 645 "zlarfb.f"
		if (*m > *k) {

/*                 C1 := C1 - V1**H * W**H */

#line 649 "zlarfb.f"
		    i__1 = *m - *k;
#line 649 "zlarfb.f"
		    z__1.r = -1., z__1.i = -0.;
#line 649 "zlarfb.f"
		    zgemm_("Conjugate transpose", "Conjugate transpose", &
			    i__1, n, k, &z__1, &v[v_offset], ldv, &work[
			    work_offset], ldwork, &c_b1, &c__[c_offset], ldc, 
			    (ftnlen)19, (ftnlen)19);
#line 652 "zlarfb.f"
		}

/*              W := W * V2 */

#line 656 "zlarfb.f"
		ztrmm_("Right", "Lower", "No transpose", "Unit", n, k, &c_b1, 
			&v[(*m - *k + 1) * v_dim1 + 1], ldv, &work[
			work_offset], ldwork, (ftnlen)5, (ftnlen)5, (ftnlen)
			12, (ftnlen)4);

/*              C2 := C2 - W**H */

#line 661 "zlarfb.f"
		i__1 = *k;
#line 661 "zlarfb.f"
		for (j = 1; j <= i__1; ++j) {
#line 662 "zlarfb.f"
		    i__2 = *n;
#line 662 "zlarfb.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 663 "zlarfb.f"
			i__3 = *m - *k + j + i__ * c_dim1;
#line 663 "zlarfb.f"
			i__4 = *m - *k + j + i__ * c_dim1;
#line 663 "zlarfb.f"
			d_cnjg(&z__2, &work[i__ + j * work_dim1]);
#line 663 "zlarfb.f"
			z__1.r = c__[i__4].r - z__2.r, z__1.i = c__[i__4].i - 
				z__2.i;
#line 663 "zlarfb.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 665 "zlarfb.f"
/* L200: */
#line 665 "zlarfb.f"
		    }
#line 666 "zlarfb.f"
/* L210: */
#line 666 "zlarfb.f"
		}

#line 668 "zlarfb.f"
	    } else if (lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {

/*              Form  C * H  or  C * H**H  where  C = ( C1  C2 ) */

/*              W := C * V**H  =  (C1*V1**H + C2*V2**H)  (stored in WORK) */

/*              W := C2 */

#line 676 "zlarfb.f"
		i__1 = *k;
#line 676 "zlarfb.f"
		for (j = 1; j <= i__1; ++j) {
#line 677 "zlarfb.f"
		    zcopy_(m, &c__[(*n - *k + j) * c_dim1 + 1], &c__1, &work[
			    j * work_dim1 + 1], &c__1);
#line 678 "zlarfb.f"
/* L220: */
#line 678 "zlarfb.f"
		}

/*              W := W * V2**H */

#line 682 "zlarfb.f"
		ztrmm_("Right", "Lower", "Conjugate transpose", "Unit", m, k, 
			&c_b1, &v[(*n - *k + 1) * v_dim1 + 1], ldv, &work[
			work_offset], ldwork, (ftnlen)5, (ftnlen)5, (ftnlen)
			19, (ftnlen)4);
#line 685 "zlarfb.f"
		if (*n > *k) {

/*                 W := W + C1 * V1**H */

#line 689 "zlarfb.f"
		    i__1 = *n - *k;
#line 689 "zlarfb.f"
		    zgemm_("No transpose", "Conjugate transpose", m, k, &i__1,
			     &c_b1, &c__[c_offset], ldc, &v[v_offset], ldv, &
			    c_b1, &work[work_offset], ldwork, (ftnlen)12, (
			    ftnlen)19);
#line 692 "zlarfb.f"
		}

/*              W := W * T  or  W * T**H */

#line 696 "zlarfb.f"
		ztrmm_("Right", "Lower", trans, "Non-unit", m, k, &c_b1, &t[
			t_offset], ldt, &work[work_offset], ldwork, (ftnlen)5,
			 (ftnlen)5, (ftnlen)1, (ftnlen)8);

/*              C := C - W * V */

#line 701 "zlarfb.f"
		if (*n > *k) {

/*                 C1 := C1 - W * V1 */

#line 705 "zlarfb.f"
		    i__1 = *n - *k;
#line 705 "zlarfb.f"
		    z__1.r = -1., z__1.i = -0.;
#line 705 "zlarfb.f"
		    zgemm_("No transpose", "No transpose", m, &i__1, k, &z__1,
			     &work[work_offset], ldwork, &v[v_offset], ldv, &
			    c_b1, &c__[c_offset], ldc, (ftnlen)12, (ftnlen)12)
			    ;
#line 707 "zlarfb.f"
		}

/*              W := W * V2 */

#line 711 "zlarfb.f"
		ztrmm_("Right", "Lower", "No transpose", "Unit", m, k, &c_b1, 
			&v[(*n - *k + 1) * v_dim1 + 1], ldv, &work[
			work_offset], ldwork, (ftnlen)5, (ftnlen)5, (ftnlen)
			12, (ftnlen)4);

/*              C1 := C1 - W */

#line 716 "zlarfb.f"
		i__1 = *k;
#line 716 "zlarfb.f"
		for (j = 1; j <= i__1; ++j) {
#line 717 "zlarfb.f"
		    i__2 = *m;
#line 717 "zlarfb.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 718 "zlarfb.f"
			i__3 = i__ + (*n - *k + j) * c_dim1;
#line 718 "zlarfb.f"
			i__4 = i__ + (*n - *k + j) * c_dim1;
#line 718 "zlarfb.f"
			i__5 = i__ + j * work_dim1;
#line 718 "zlarfb.f"
			z__1.r = c__[i__4].r - work[i__5].r, z__1.i = c__[
				i__4].i - work[i__5].i;
#line 718 "zlarfb.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 719 "zlarfb.f"
/* L230: */
#line 719 "zlarfb.f"
		    }
#line 720 "zlarfb.f"
/* L240: */
#line 720 "zlarfb.f"
		}

#line 722 "zlarfb.f"
	    }

#line 724 "zlarfb.f"
	}
#line 725 "zlarfb.f"
    }

#line 727 "zlarfb.f"
    return 0;

/*     End of ZLARFB */

} /* zlarfb_ */

