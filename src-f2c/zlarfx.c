#line 1 "zlarfx.f"
/* zlarfx.f -- translated by f2c (version 20100827).
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

#line 1 "zlarfx.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b ZLARFX applies an elementary reflector to a general rectangular matrix, with loop unrolling whe
n the reflector has order ≤ 10. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZLARFX + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlarfx.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlarfx.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlarfx.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZLARFX( SIDE, M, N, V, TAU, C, LDC, WORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          SIDE */
/*       INTEGER            LDC, M, N */
/*       COMPLEX*16         TAU */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         C( LDC, * ), V( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLARFX applies a complex elementary reflector H to a complex m by n */
/* > matrix C, from either the left or the right. H is represented in the */
/* > form */
/* > */
/* >       H = I - tau * v * v**H */
/* > */
/* > where tau is a complex scalar and v is a complex vector. */
/* > */
/* > If tau = 0, then H is taken to be the unit matrix */
/* > */
/* > This version uses inline code if H has order < 11. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] SIDE */
/* > \verbatim */
/* >          SIDE is CHARACTER*1 */
/* >          = 'L': form  H * C */
/* >          = 'R': form  C * H */
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
/* > \param[in] V */
/* > \verbatim */
/* >          V is COMPLEX*16 array, dimension (M) if SIDE = 'L' */
/* >                                        or (N) if SIDE = 'R' */
/* >          The vector v in the representation of H. */
/* > \endverbatim */
/* > */
/* > \param[in] TAU */
/* > \verbatim */
/* >          TAU is COMPLEX*16 */
/* >          The value tau in the representation of H. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* >          C is COMPLEX*16 array, dimension (LDC,N) */
/* >          On entry, the m by n matrix C. */
/* >          On exit, C is overwritten by the matrix H * C if SIDE = 'L', */
/* >          or C * H if SIDE = 'R'. */
/* > \endverbatim */
/* > */
/* > \param[in] LDC */
/* > \verbatim */
/* >          LDC is INTEGER */
/* >          The leading dimension of the array C. LDA >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX*16 array, dimension (N) if SIDE = 'L' */
/* >                                            or (M) if SIDE = 'R' */
/* >          WORK is not referenced if H has order < 11. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup complex16OTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int zlarfx_(char *side, integer *m, integer *n, 
	doublecomplex *v, doublecomplex *tau, doublecomplex *c__, integer *
	ldc, doublecomplex *work, ftnlen side_len)
{
    /* System generated locals */
    integer c_dim1, c_offset, i__1, i__2, i__3, i__4, i__5, i__6, i__7, i__8, 
	    i__9, i__10, i__11;
    doublecomplex z__1, z__2, z__3, z__4, z__5, z__6, z__7, z__8, z__9, z__10,
	     z__11, z__12, z__13, z__14, z__15, z__16, z__17, z__18, z__19;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer j;
    static doublecomplex t1, t2, t3, t4, t5, t6, t7, t8, t9, v1, v2, v3, v4, 
	    v5, v6, v7, v8, v9, t10, v10, sum;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int zlarf_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, doublecomplex *, ftnlen);


/*  -- LAPACK auxiliary routine (version 3.4.2) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     September 2012 */

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

#line 160 "zlarfx.f"
    /* Parameter adjustments */
#line 160 "zlarfx.f"
    --v;
#line 160 "zlarfx.f"
    c_dim1 = *ldc;
#line 160 "zlarfx.f"
    c_offset = 1 + c_dim1;
#line 160 "zlarfx.f"
    c__ -= c_offset;
#line 160 "zlarfx.f"
    --work;
#line 160 "zlarfx.f"

#line 160 "zlarfx.f"
    /* Function Body */
#line 160 "zlarfx.f"
    if (tau->r == 0. && tau->i == 0.) {
#line 160 "zlarfx.f"
	return 0;
#line 160 "zlarfx.f"
    }
#line 162 "zlarfx.f"
    if (lsame_(side, "L", (ftnlen)1, (ftnlen)1)) {

/*        Form  H * C, where H has order m. */

#line 166 "zlarfx.f"
	switch (*m) {
#line 166 "zlarfx.f"
	    case 1:  goto L10;
#line 166 "zlarfx.f"
	    case 2:  goto L30;
#line 166 "zlarfx.f"
	    case 3:  goto L50;
#line 166 "zlarfx.f"
	    case 4:  goto L70;
#line 166 "zlarfx.f"
	    case 5:  goto L90;
#line 166 "zlarfx.f"
	    case 6:  goto L110;
#line 166 "zlarfx.f"
	    case 7:  goto L130;
#line 166 "zlarfx.f"
	    case 8:  goto L150;
#line 166 "zlarfx.f"
	    case 9:  goto L170;
#line 166 "zlarfx.f"
	    case 10:  goto L190;
#line 166 "zlarfx.f"
	}

/*        Code for general M */

#line 171 "zlarfx.f"
	zlarf_(side, m, n, &v[1], &c__1, tau, &c__[c_offset], ldc, &work[1], (
		ftnlen)1);
#line 172 "zlarfx.f"
	goto L410;
#line 173 "zlarfx.f"
L10:

/*        Special code for 1 x 1 Householder */

#line 177 "zlarfx.f"
	z__3.r = tau->r * v[1].r - tau->i * v[1].i, z__3.i = tau->r * v[1].i 
		+ tau->i * v[1].r;
#line 177 "zlarfx.f"
	d_cnjg(&z__4, &v[1]);
#line 177 "zlarfx.f"
	z__2.r = z__3.r * z__4.r - z__3.i * z__4.i, z__2.i = z__3.r * z__4.i 
		+ z__3.i * z__4.r;
#line 177 "zlarfx.f"
	z__1.r = 1. - z__2.r, z__1.i = 0. - z__2.i;
#line 177 "zlarfx.f"
	t1.r = z__1.r, t1.i = z__1.i;
#line 178 "zlarfx.f"
	i__1 = *n;
#line 178 "zlarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 179 "zlarfx.f"
	    i__2 = j * c_dim1 + 1;
#line 179 "zlarfx.f"
	    i__3 = j * c_dim1 + 1;
#line 179 "zlarfx.f"
	    z__1.r = t1.r * c__[i__3].r - t1.i * c__[i__3].i, z__1.i = t1.r * 
		    c__[i__3].i + t1.i * c__[i__3].r;
#line 179 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 180 "zlarfx.f"
/* L20: */
#line 180 "zlarfx.f"
	}
#line 181 "zlarfx.f"
	goto L410;
#line 182 "zlarfx.f"
L30:

/*        Special code for 2 x 2 Householder */

#line 186 "zlarfx.f"
	d_cnjg(&z__1, &v[1]);
#line 186 "zlarfx.f"
	v1.r = z__1.r, v1.i = z__1.i;
#line 187 "zlarfx.f"
	d_cnjg(&z__2, &v1);
#line 187 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 187 "zlarfx.f"
	t1.r = z__1.r, t1.i = z__1.i;
#line 188 "zlarfx.f"
	d_cnjg(&z__1, &v[2]);
#line 188 "zlarfx.f"
	v2.r = z__1.r, v2.i = z__1.i;
#line 189 "zlarfx.f"
	d_cnjg(&z__2, &v2);
#line 189 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 189 "zlarfx.f"
	t2.r = z__1.r, t2.i = z__1.i;
#line 190 "zlarfx.f"
	i__1 = *n;
#line 190 "zlarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 191 "zlarfx.f"
	    i__2 = j * c_dim1 + 1;
#line 191 "zlarfx.f"
	    z__2.r = v1.r * c__[i__2].r - v1.i * c__[i__2].i, z__2.i = v1.r * 
		    c__[i__2].i + v1.i * c__[i__2].r;
#line 191 "zlarfx.f"
	    i__3 = j * c_dim1 + 2;
#line 191 "zlarfx.f"
	    z__3.r = v2.r * c__[i__3].r - v2.i * c__[i__3].i, z__3.i = v2.r * 
		    c__[i__3].i + v2.i * c__[i__3].r;
#line 191 "zlarfx.f"
	    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 191 "zlarfx.f"
	    sum.r = z__1.r, sum.i = z__1.i;
#line 192 "zlarfx.f"
	    i__2 = j * c_dim1 + 1;
#line 192 "zlarfx.f"
	    i__3 = j * c_dim1 + 1;
#line 192 "zlarfx.f"
	    z__2.r = sum.r * t1.r - sum.i * t1.i, z__2.i = sum.r * t1.i + 
		    sum.i * t1.r;
#line 192 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 192 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 193 "zlarfx.f"
	    i__2 = j * c_dim1 + 2;
#line 193 "zlarfx.f"
	    i__3 = j * c_dim1 + 2;
#line 193 "zlarfx.f"
	    z__2.r = sum.r * t2.r - sum.i * t2.i, z__2.i = sum.r * t2.i + 
		    sum.i * t2.r;
#line 193 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 193 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 194 "zlarfx.f"
/* L40: */
#line 194 "zlarfx.f"
	}
#line 195 "zlarfx.f"
	goto L410;
#line 196 "zlarfx.f"
L50:

/*        Special code for 3 x 3 Householder */

#line 200 "zlarfx.f"
	d_cnjg(&z__1, &v[1]);
#line 200 "zlarfx.f"
	v1.r = z__1.r, v1.i = z__1.i;
#line 201 "zlarfx.f"
	d_cnjg(&z__2, &v1);
#line 201 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 201 "zlarfx.f"
	t1.r = z__1.r, t1.i = z__1.i;
#line 202 "zlarfx.f"
	d_cnjg(&z__1, &v[2]);
#line 202 "zlarfx.f"
	v2.r = z__1.r, v2.i = z__1.i;
#line 203 "zlarfx.f"
	d_cnjg(&z__2, &v2);
#line 203 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 203 "zlarfx.f"
	t2.r = z__1.r, t2.i = z__1.i;
#line 204 "zlarfx.f"
	d_cnjg(&z__1, &v[3]);
#line 204 "zlarfx.f"
	v3.r = z__1.r, v3.i = z__1.i;
#line 205 "zlarfx.f"
	d_cnjg(&z__2, &v3);
#line 205 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 205 "zlarfx.f"
	t3.r = z__1.r, t3.i = z__1.i;
#line 206 "zlarfx.f"
	i__1 = *n;
#line 206 "zlarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 207 "zlarfx.f"
	    i__2 = j * c_dim1 + 1;
#line 207 "zlarfx.f"
	    z__3.r = v1.r * c__[i__2].r - v1.i * c__[i__2].i, z__3.i = v1.r * 
		    c__[i__2].i + v1.i * c__[i__2].r;
#line 207 "zlarfx.f"
	    i__3 = j * c_dim1 + 2;
#line 207 "zlarfx.f"
	    z__4.r = v2.r * c__[i__3].r - v2.i * c__[i__3].i, z__4.i = v2.r * 
		    c__[i__3].i + v2.i * c__[i__3].r;
#line 207 "zlarfx.f"
	    z__2.r = z__3.r + z__4.r, z__2.i = z__3.i + z__4.i;
#line 207 "zlarfx.f"
	    i__4 = j * c_dim1 + 3;
#line 207 "zlarfx.f"
	    z__5.r = v3.r * c__[i__4].r - v3.i * c__[i__4].i, z__5.i = v3.r * 
		    c__[i__4].i + v3.i * c__[i__4].r;
#line 207 "zlarfx.f"
	    z__1.r = z__2.r + z__5.r, z__1.i = z__2.i + z__5.i;
#line 207 "zlarfx.f"
	    sum.r = z__1.r, sum.i = z__1.i;
#line 208 "zlarfx.f"
	    i__2 = j * c_dim1 + 1;
#line 208 "zlarfx.f"
	    i__3 = j * c_dim1 + 1;
#line 208 "zlarfx.f"
	    z__2.r = sum.r * t1.r - sum.i * t1.i, z__2.i = sum.r * t1.i + 
		    sum.i * t1.r;
#line 208 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 208 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 209 "zlarfx.f"
	    i__2 = j * c_dim1 + 2;
#line 209 "zlarfx.f"
	    i__3 = j * c_dim1 + 2;
#line 209 "zlarfx.f"
	    z__2.r = sum.r * t2.r - sum.i * t2.i, z__2.i = sum.r * t2.i + 
		    sum.i * t2.r;
#line 209 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 209 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 210 "zlarfx.f"
	    i__2 = j * c_dim1 + 3;
#line 210 "zlarfx.f"
	    i__3 = j * c_dim1 + 3;
#line 210 "zlarfx.f"
	    z__2.r = sum.r * t3.r - sum.i * t3.i, z__2.i = sum.r * t3.i + 
		    sum.i * t3.r;
#line 210 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 210 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 211 "zlarfx.f"
/* L60: */
#line 211 "zlarfx.f"
	}
#line 212 "zlarfx.f"
	goto L410;
#line 213 "zlarfx.f"
L70:

/*        Special code for 4 x 4 Householder */

#line 217 "zlarfx.f"
	d_cnjg(&z__1, &v[1]);
#line 217 "zlarfx.f"
	v1.r = z__1.r, v1.i = z__1.i;
#line 218 "zlarfx.f"
	d_cnjg(&z__2, &v1);
#line 218 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 218 "zlarfx.f"
	t1.r = z__1.r, t1.i = z__1.i;
#line 219 "zlarfx.f"
	d_cnjg(&z__1, &v[2]);
#line 219 "zlarfx.f"
	v2.r = z__1.r, v2.i = z__1.i;
#line 220 "zlarfx.f"
	d_cnjg(&z__2, &v2);
#line 220 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 220 "zlarfx.f"
	t2.r = z__1.r, t2.i = z__1.i;
#line 221 "zlarfx.f"
	d_cnjg(&z__1, &v[3]);
#line 221 "zlarfx.f"
	v3.r = z__1.r, v3.i = z__1.i;
#line 222 "zlarfx.f"
	d_cnjg(&z__2, &v3);
#line 222 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 222 "zlarfx.f"
	t3.r = z__1.r, t3.i = z__1.i;
#line 223 "zlarfx.f"
	d_cnjg(&z__1, &v[4]);
#line 223 "zlarfx.f"
	v4.r = z__1.r, v4.i = z__1.i;
#line 224 "zlarfx.f"
	d_cnjg(&z__2, &v4);
#line 224 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 224 "zlarfx.f"
	t4.r = z__1.r, t4.i = z__1.i;
#line 225 "zlarfx.f"
	i__1 = *n;
#line 225 "zlarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 226 "zlarfx.f"
	    i__2 = j * c_dim1 + 1;
#line 226 "zlarfx.f"
	    z__4.r = v1.r * c__[i__2].r - v1.i * c__[i__2].i, z__4.i = v1.r * 
		    c__[i__2].i + v1.i * c__[i__2].r;
#line 226 "zlarfx.f"
	    i__3 = j * c_dim1 + 2;
#line 226 "zlarfx.f"
	    z__5.r = v2.r * c__[i__3].r - v2.i * c__[i__3].i, z__5.i = v2.r * 
		    c__[i__3].i + v2.i * c__[i__3].r;
#line 226 "zlarfx.f"
	    z__3.r = z__4.r + z__5.r, z__3.i = z__4.i + z__5.i;
#line 226 "zlarfx.f"
	    i__4 = j * c_dim1 + 3;
#line 226 "zlarfx.f"
	    z__6.r = v3.r * c__[i__4].r - v3.i * c__[i__4].i, z__6.i = v3.r * 
		    c__[i__4].i + v3.i * c__[i__4].r;
#line 226 "zlarfx.f"
	    z__2.r = z__3.r + z__6.r, z__2.i = z__3.i + z__6.i;
#line 226 "zlarfx.f"
	    i__5 = j * c_dim1 + 4;
#line 226 "zlarfx.f"
	    z__7.r = v4.r * c__[i__5].r - v4.i * c__[i__5].i, z__7.i = v4.r * 
		    c__[i__5].i + v4.i * c__[i__5].r;
#line 226 "zlarfx.f"
	    z__1.r = z__2.r + z__7.r, z__1.i = z__2.i + z__7.i;
#line 226 "zlarfx.f"
	    sum.r = z__1.r, sum.i = z__1.i;
#line 228 "zlarfx.f"
	    i__2 = j * c_dim1 + 1;
#line 228 "zlarfx.f"
	    i__3 = j * c_dim1 + 1;
#line 228 "zlarfx.f"
	    z__2.r = sum.r * t1.r - sum.i * t1.i, z__2.i = sum.r * t1.i + 
		    sum.i * t1.r;
#line 228 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 228 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 229 "zlarfx.f"
	    i__2 = j * c_dim1 + 2;
#line 229 "zlarfx.f"
	    i__3 = j * c_dim1 + 2;
#line 229 "zlarfx.f"
	    z__2.r = sum.r * t2.r - sum.i * t2.i, z__2.i = sum.r * t2.i + 
		    sum.i * t2.r;
#line 229 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 229 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 230 "zlarfx.f"
	    i__2 = j * c_dim1 + 3;
#line 230 "zlarfx.f"
	    i__3 = j * c_dim1 + 3;
#line 230 "zlarfx.f"
	    z__2.r = sum.r * t3.r - sum.i * t3.i, z__2.i = sum.r * t3.i + 
		    sum.i * t3.r;
#line 230 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 230 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 231 "zlarfx.f"
	    i__2 = j * c_dim1 + 4;
#line 231 "zlarfx.f"
	    i__3 = j * c_dim1 + 4;
#line 231 "zlarfx.f"
	    z__2.r = sum.r * t4.r - sum.i * t4.i, z__2.i = sum.r * t4.i + 
		    sum.i * t4.r;
#line 231 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 231 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 232 "zlarfx.f"
/* L80: */
#line 232 "zlarfx.f"
	}
#line 233 "zlarfx.f"
	goto L410;
#line 234 "zlarfx.f"
L90:

/*        Special code for 5 x 5 Householder */

#line 238 "zlarfx.f"
	d_cnjg(&z__1, &v[1]);
#line 238 "zlarfx.f"
	v1.r = z__1.r, v1.i = z__1.i;
#line 239 "zlarfx.f"
	d_cnjg(&z__2, &v1);
#line 239 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 239 "zlarfx.f"
	t1.r = z__1.r, t1.i = z__1.i;
#line 240 "zlarfx.f"
	d_cnjg(&z__1, &v[2]);
#line 240 "zlarfx.f"
	v2.r = z__1.r, v2.i = z__1.i;
#line 241 "zlarfx.f"
	d_cnjg(&z__2, &v2);
#line 241 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 241 "zlarfx.f"
	t2.r = z__1.r, t2.i = z__1.i;
#line 242 "zlarfx.f"
	d_cnjg(&z__1, &v[3]);
#line 242 "zlarfx.f"
	v3.r = z__1.r, v3.i = z__1.i;
#line 243 "zlarfx.f"
	d_cnjg(&z__2, &v3);
#line 243 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 243 "zlarfx.f"
	t3.r = z__1.r, t3.i = z__1.i;
#line 244 "zlarfx.f"
	d_cnjg(&z__1, &v[4]);
#line 244 "zlarfx.f"
	v4.r = z__1.r, v4.i = z__1.i;
#line 245 "zlarfx.f"
	d_cnjg(&z__2, &v4);
#line 245 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 245 "zlarfx.f"
	t4.r = z__1.r, t4.i = z__1.i;
#line 246 "zlarfx.f"
	d_cnjg(&z__1, &v[5]);
#line 246 "zlarfx.f"
	v5.r = z__1.r, v5.i = z__1.i;
#line 247 "zlarfx.f"
	d_cnjg(&z__2, &v5);
#line 247 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 247 "zlarfx.f"
	t5.r = z__1.r, t5.i = z__1.i;
#line 248 "zlarfx.f"
	i__1 = *n;
#line 248 "zlarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 249 "zlarfx.f"
	    i__2 = j * c_dim1 + 1;
#line 249 "zlarfx.f"
	    z__5.r = v1.r * c__[i__2].r - v1.i * c__[i__2].i, z__5.i = v1.r * 
		    c__[i__2].i + v1.i * c__[i__2].r;
#line 249 "zlarfx.f"
	    i__3 = j * c_dim1 + 2;
#line 249 "zlarfx.f"
	    z__6.r = v2.r * c__[i__3].r - v2.i * c__[i__3].i, z__6.i = v2.r * 
		    c__[i__3].i + v2.i * c__[i__3].r;
#line 249 "zlarfx.f"
	    z__4.r = z__5.r + z__6.r, z__4.i = z__5.i + z__6.i;
#line 249 "zlarfx.f"
	    i__4 = j * c_dim1 + 3;
#line 249 "zlarfx.f"
	    z__7.r = v3.r * c__[i__4].r - v3.i * c__[i__4].i, z__7.i = v3.r * 
		    c__[i__4].i + v3.i * c__[i__4].r;
#line 249 "zlarfx.f"
	    z__3.r = z__4.r + z__7.r, z__3.i = z__4.i + z__7.i;
#line 249 "zlarfx.f"
	    i__5 = j * c_dim1 + 4;
#line 249 "zlarfx.f"
	    z__8.r = v4.r * c__[i__5].r - v4.i * c__[i__5].i, z__8.i = v4.r * 
		    c__[i__5].i + v4.i * c__[i__5].r;
#line 249 "zlarfx.f"
	    z__2.r = z__3.r + z__8.r, z__2.i = z__3.i + z__8.i;
#line 249 "zlarfx.f"
	    i__6 = j * c_dim1 + 5;
#line 249 "zlarfx.f"
	    z__9.r = v5.r * c__[i__6].r - v5.i * c__[i__6].i, z__9.i = v5.r * 
		    c__[i__6].i + v5.i * c__[i__6].r;
#line 249 "zlarfx.f"
	    z__1.r = z__2.r + z__9.r, z__1.i = z__2.i + z__9.i;
#line 249 "zlarfx.f"
	    sum.r = z__1.r, sum.i = z__1.i;
#line 251 "zlarfx.f"
	    i__2 = j * c_dim1 + 1;
#line 251 "zlarfx.f"
	    i__3 = j * c_dim1 + 1;
#line 251 "zlarfx.f"
	    z__2.r = sum.r * t1.r - sum.i * t1.i, z__2.i = sum.r * t1.i + 
		    sum.i * t1.r;
#line 251 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 251 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 252 "zlarfx.f"
	    i__2 = j * c_dim1 + 2;
#line 252 "zlarfx.f"
	    i__3 = j * c_dim1 + 2;
#line 252 "zlarfx.f"
	    z__2.r = sum.r * t2.r - sum.i * t2.i, z__2.i = sum.r * t2.i + 
		    sum.i * t2.r;
#line 252 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 252 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 253 "zlarfx.f"
	    i__2 = j * c_dim1 + 3;
#line 253 "zlarfx.f"
	    i__3 = j * c_dim1 + 3;
#line 253 "zlarfx.f"
	    z__2.r = sum.r * t3.r - sum.i * t3.i, z__2.i = sum.r * t3.i + 
		    sum.i * t3.r;
#line 253 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 253 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 254 "zlarfx.f"
	    i__2 = j * c_dim1 + 4;
#line 254 "zlarfx.f"
	    i__3 = j * c_dim1 + 4;
#line 254 "zlarfx.f"
	    z__2.r = sum.r * t4.r - sum.i * t4.i, z__2.i = sum.r * t4.i + 
		    sum.i * t4.r;
#line 254 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 254 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 255 "zlarfx.f"
	    i__2 = j * c_dim1 + 5;
#line 255 "zlarfx.f"
	    i__3 = j * c_dim1 + 5;
#line 255 "zlarfx.f"
	    z__2.r = sum.r * t5.r - sum.i * t5.i, z__2.i = sum.r * t5.i + 
		    sum.i * t5.r;
#line 255 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 255 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 256 "zlarfx.f"
/* L100: */
#line 256 "zlarfx.f"
	}
#line 257 "zlarfx.f"
	goto L410;
#line 258 "zlarfx.f"
L110:

/*        Special code for 6 x 6 Householder */

#line 262 "zlarfx.f"
	d_cnjg(&z__1, &v[1]);
#line 262 "zlarfx.f"
	v1.r = z__1.r, v1.i = z__1.i;
#line 263 "zlarfx.f"
	d_cnjg(&z__2, &v1);
#line 263 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 263 "zlarfx.f"
	t1.r = z__1.r, t1.i = z__1.i;
#line 264 "zlarfx.f"
	d_cnjg(&z__1, &v[2]);
#line 264 "zlarfx.f"
	v2.r = z__1.r, v2.i = z__1.i;
#line 265 "zlarfx.f"
	d_cnjg(&z__2, &v2);
#line 265 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 265 "zlarfx.f"
	t2.r = z__1.r, t2.i = z__1.i;
#line 266 "zlarfx.f"
	d_cnjg(&z__1, &v[3]);
#line 266 "zlarfx.f"
	v3.r = z__1.r, v3.i = z__1.i;
#line 267 "zlarfx.f"
	d_cnjg(&z__2, &v3);
#line 267 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 267 "zlarfx.f"
	t3.r = z__1.r, t3.i = z__1.i;
#line 268 "zlarfx.f"
	d_cnjg(&z__1, &v[4]);
#line 268 "zlarfx.f"
	v4.r = z__1.r, v4.i = z__1.i;
#line 269 "zlarfx.f"
	d_cnjg(&z__2, &v4);
#line 269 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 269 "zlarfx.f"
	t4.r = z__1.r, t4.i = z__1.i;
#line 270 "zlarfx.f"
	d_cnjg(&z__1, &v[5]);
#line 270 "zlarfx.f"
	v5.r = z__1.r, v5.i = z__1.i;
#line 271 "zlarfx.f"
	d_cnjg(&z__2, &v5);
#line 271 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 271 "zlarfx.f"
	t5.r = z__1.r, t5.i = z__1.i;
#line 272 "zlarfx.f"
	d_cnjg(&z__1, &v[6]);
#line 272 "zlarfx.f"
	v6.r = z__1.r, v6.i = z__1.i;
#line 273 "zlarfx.f"
	d_cnjg(&z__2, &v6);
#line 273 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 273 "zlarfx.f"
	t6.r = z__1.r, t6.i = z__1.i;
#line 274 "zlarfx.f"
	i__1 = *n;
#line 274 "zlarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 275 "zlarfx.f"
	    i__2 = j * c_dim1 + 1;
#line 275 "zlarfx.f"
	    z__6.r = v1.r * c__[i__2].r - v1.i * c__[i__2].i, z__6.i = v1.r * 
		    c__[i__2].i + v1.i * c__[i__2].r;
#line 275 "zlarfx.f"
	    i__3 = j * c_dim1 + 2;
#line 275 "zlarfx.f"
	    z__7.r = v2.r * c__[i__3].r - v2.i * c__[i__3].i, z__7.i = v2.r * 
		    c__[i__3].i + v2.i * c__[i__3].r;
#line 275 "zlarfx.f"
	    z__5.r = z__6.r + z__7.r, z__5.i = z__6.i + z__7.i;
#line 275 "zlarfx.f"
	    i__4 = j * c_dim1 + 3;
#line 275 "zlarfx.f"
	    z__8.r = v3.r * c__[i__4].r - v3.i * c__[i__4].i, z__8.i = v3.r * 
		    c__[i__4].i + v3.i * c__[i__4].r;
#line 275 "zlarfx.f"
	    z__4.r = z__5.r + z__8.r, z__4.i = z__5.i + z__8.i;
#line 275 "zlarfx.f"
	    i__5 = j * c_dim1 + 4;
#line 275 "zlarfx.f"
	    z__9.r = v4.r * c__[i__5].r - v4.i * c__[i__5].i, z__9.i = v4.r * 
		    c__[i__5].i + v4.i * c__[i__5].r;
#line 275 "zlarfx.f"
	    z__3.r = z__4.r + z__9.r, z__3.i = z__4.i + z__9.i;
#line 275 "zlarfx.f"
	    i__6 = j * c_dim1 + 5;
#line 275 "zlarfx.f"
	    z__10.r = v5.r * c__[i__6].r - v5.i * c__[i__6].i, z__10.i = v5.r 
		    * c__[i__6].i + v5.i * c__[i__6].r;
#line 275 "zlarfx.f"
	    z__2.r = z__3.r + z__10.r, z__2.i = z__3.i + z__10.i;
#line 275 "zlarfx.f"
	    i__7 = j * c_dim1 + 6;
#line 275 "zlarfx.f"
	    z__11.r = v6.r * c__[i__7].r - v6.i * c__[i__7].i, z__11.i = v6.r 
		    * c__[i__7].i + v6.i * c__[i__7].r;
#line 275 "zlarfx.f"
	    z__1.r = z__2.r + z__11.r, z__1.i = z__2.i + z__11.i;
#line 275 "zlarfx.f"
	    sum.r = z__1.r, sum.i = z__1.i;
#line 277 "zlarfx.f"
	    i__2 = j * c_dim1 + 1;
#line 277 "zlarfx.f"
	    i__3 = j * c_dim1 + 1;
#line 277 "zlarfx.f"
	    z__2.r = sum.r * t1.r - sum.i * t1.i, z__2.i = sum.r * t1.i + 
		    sum.i * t1.r;
#line 277 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 277 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 278 "zlarfx.f"
	    i__2 = j * c_dim1 + 2;
#line 278 "zlarfx.f"
	    i__3 = j * c_dim1 + 2;
#line 278 "zlarfx.f"
	    z__2.r = sum.r * t2.r - sum.i * t2.i, z__2.i = sum.r * t2.i + 
		    sum.i * t2.r;
#line 278 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 278 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 279 "zlarfx.f"
	    i__2 = j * c_dim1 + 3;
#line 279 "zlarfx.f"
	    i__3 = j * c_dim1 + 3;
#line 279 "zlarfx.f"
	    z__2.r = sum.r * t3.r - sum.i * t3.i, z__2.i = sum.r * t3.i + 
		    sum.i * t3.r;
#line 279 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 279 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 280 "zlarfx.f"
	    i__2 = j * c_dim1 + 4;
#line 280 "zlarfx.f"
	    i__3 = j * c_dim1 + 4;
#line 280 "zlarfx.f"
	    z__2.r = sum.r * t4.r - sum.i * t4.i, z__2.i = sum.r * t4.i + 
		    sum.i * t4.r;
#line 280 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 280 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 281 "zlarfx.f"
	    i__2 = j * c_dim1 + 5;
#line 281 "zlarfx.f"
	    i__3 = j * c_dim1 + 5;
#line 281 "zlarfx.f"
	    z__2.r = sum.r * t5.r - sum.i * t5.i, z__2.i = sum.r * t5.i + 
		    sum.i * t5.r;
#line 281 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 281 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 282 "zlarfx.f"
	    i__2 = j * c_dim1 + 6;
#line 282 "zlarfx.f"
	    i__3 = j * c_dim1 + 6;
#line 282 "zlarfx.f"
	    z__2.r = sum.r * t6.r - sum.i * t6.i, z__2.i = sum.r * t6.i + 
		    sum.i * t6.r;
#line 282 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 282 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 283 "zlarfx.f"
/* L120: */
#line 283 "zlarfx.f"
	}
#line 284 "zlarfx.f"
	goto L410;
#line 285 "zlarfx.f"
L130:

/*        Special code for 7 x 7 Householder */

#line 289 "zlarfx.f"
	d_cnjg(&z__1, &v[1]);
#line 289 "zlarfx.f"
	v1.r = z__1.r, v1.i = z__1.i;
#line 290 "zlarfx.f"
	d_cnjg(&z__2, &v1);
#line 290 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 290 "zlarfx.f"
	t1.r = z__1.r, t1.i = z__1.i;
#line 291 "zlarfx.f"
	d_cnjg(&z__1, &v[2]);
#line 291 "zlarfx.f"
	v2.r = z__1.r, v2.i = z__1.i;
#line 292 "zlarfx.f"
	d_cnjg(&z__2, &v2);
#line 292 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 292 "zlarfx.f"
	t2.r = z__1.r, t2.i = z__1.i;
#line 293 "zlarfx.f"
	d_cnjg(&z__1, &v[3]);
#line 293 "zlarfx.f"
	v3.r = z__1.r, v3.i = z__1.i;
#line 294 "zlarfx.f"
	d_cnjg(&z__2, &v3);
#line 294 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 294 "zlarfx.f"
	t3.r = z__1.r, t3.i = z__1.i;
#line 295 "zlarfx.f"
	d_cnjg(&z__1, &v[4]);
#line 295 "zlarfx.f"
	v4.r = z__1.r, v4.i = z__1.i;
#line 296 "zlarfx.f"
	d_cnjg(&z__2, &v4);
#line 296 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 296 "zlarfx.f"
	t4.r = z__1.r, t4.i = z__1.i;
#line 297 "zlarfx.f"
	d_cnjg(&z__1, &v[5]);
#line 297 "zlarfx.f"
	v5.r = z__1.r, v5.i = z__1.i;
#line 298 "zlarfx.f"
	d_cnjg(&z__2, &v5);
#line 298 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 298 "zlarfx.f"
	t5.r = z__1.r, t5.i = z__1.i;
#line 299 "zlarfx.f"
	d_cnjg(&z__1, &v[6]);
#line 299 "zlarfx.f"
	v6.r = z__1.r, v6.i = z__1.i;
#line 300 "zlarfx.f"
	d_cnjg(&z__2, &v6);
#line 300 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 300 "zlarfx.f"
	t6.r = z__1.r, t6.i = z__1.i;
#line 301 "zlarfx.f"
	d_cnjg(&z__1, &v[7]);
#line 301 "zlarfx.f"
	v7.r = z__1.r, v7.i = z__1.i;
#line 302 "zlarfx.f"
	d_cnjg(&z__2, &v7);
#line 302 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 302 "zlarfx.f"
	t7.r = z__1.r, t7.i = z__1.i;
#line 303 "zlarfx.f"
	i__1 = *n;
#line 303 "zlarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 304 "zlarfx.f"
	    i__2 = j * c_dim1 + 1;
#line 304 "zlarfx.f"
	    z__7.r = v1.r * c__[i__2].r - v1.i * c__[i__2].i, z__7.i = v1.r * 
		    c__[i__2].i + v1.i * c__[i__2].r;
#line 304 "zlarfx.f"
	    i__3 = j * c_dim1 + 2;
#line 304 "zlarfx.f"
	    z__8.r = v2.r * c__[i__3].r - v2.i * c__[i__3].i, z__8.i = v2.r * 
		    c__[i__3].i + v2.i * c__[i__3].r;
#line 304 "zlarfx.f"
	    z__6.r = z__7.r + z__8.r, z__6.i = z__7.i + z__8.i;
#line 304 "zlarfx.f"
	    i__4 = j * c_dim1 + 3;
#line 304 "zlarfx.f"
	    z__9.r = v3.r * c__[i__4].r - v3.i * c__[i__4].i, z__9.i = v3.r * 
		    c__[i__4].i + v3.i * c__[i__4].r;
#line 304 "zlarfx.f"
	    z__5.r = z__6.r + z__9.r, z__5.i = z__6.i + z__9.i;
#line 304 "zlarfx.f"
	    i__5 = j * c_dim1 + 4;
#line 304 "zlarfx.f"
	    z__10.r = v4.r * c__[i__5].r - v4.i * c__[i__5].i, z__10.i = v4.r 
		    * c__[i__5].i + v4.i * c__[i__5].r;
#line 304 "zlarfx.f"
	    z__4.r = z__5.r + z__10.r, z__4.i = z__5.i + z__10.i;
#line 304 "zlarfx.f"
	    i__6 = j * c_dim1 + 5;
#line 304 "zlarfx.f"
	    z__11.r = v5.r * c__[i__6].r - v5.i * c__[i__6].i, z__11.i = v5.r 
		    * c__[i__6].i + v5.i * c__[i__6].r;
#line 304 "zlarfx.f"
	    z__3.r = z__4.r + z__11.r, z__3.i = z__4.i + z__11.i;
#line 304 "zlarfx.f"
	    i__7 = j * c_dim1 + 6;
#line 304 "zlarfx.f"
	    z__12.r = v6.r * c__[i__7].r - v6.i * c__[i__7].i, z__12.i = v6.r 
		    * c__[i__7].i + v6.i * c__[i__7].r;
#line 304 "zlarfx.f"
	    z__2.r = z__3.r + z__12.r, z__2.i = z__3.i + z__12.i;
#line 304 "zlarfx.f"
	    i__8 = j * c_dim1 + 7;
#line 304 "zlarfx.f"
	    z__13.r = v7.r * c__[i__8].r - v7.i * c__[i__8].i, z__13.i = v7.r 
		    * c__[i__8].i + v7.i * c__[i__8].r;
#line 304 "zlarfx.f"
	    z__1.r = z__2.r + z__13.r, z__1.i = z__2.i + z__13.i;
#line 304 "zlarfx.f"
	    sum.r = z__1.r, sum.i = z__1.i;
#line 307 "zlarfx.f"
	    i__2 = j * c_dim1 + 1;
#line 307 "zlarfx.f"
	    i__3 = j * c_dim1 + 1;
#line 307 "zlarfx.f"
	    z__2.r = sum.r * t1.r - sum.i * t1.i, z__2.i = sum.r * t1.i + 
		    sum.i * t1.r;
#line 307 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 307 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 308 "zlarfx.f"
	    i__2 = j * c_dim1 + 2;
#line 308 "zlarfx.f"
	    i__3 = j * c_dim1 + 2;
#line 308 "zlarfx.f"
	    z__2.r = sum.r * t2.r - sum.i * t2.i, z__2.i = sum.r * t2.i + 
		    sum.i * t2.r;
#line 308 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 308 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 309 "zlarfx.f"
	    i__2 = j * c_dim1 + 3;
#line 309 "zlarfx.f"
	    i__3 = j * c_dim1 + 3;
#line 309 "zlarfx.f"
	    z__2.r = sum.r * t3.r - sum.i * t3.i, z__2.i = sum.r * t3.i + 
		    sum.i * t3.r;
#line 309 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 309 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 310 "zlarfx.f"
	    i__2 = j * c_dim1 + 4;
#line 310 "zlarfx.f"
	    i__3 = j * c_dim1 + 4;
#line 310 "zlarfx.f"
	    z__2.r = sum.r * t4.r - sum.i * t4.i, z__2.i = sum.r * t4.i + 
		    sum.i * t4.r;
#line 310 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 310 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 311 "zlarfx.f"
	    i__2 = j * c_dim1 + 5;
#line 311 "zlarfx.f"
	    i__3 = j * c_dim1 + 5;
#line 311 "zlarfx.f"
	    z__2.r = sum.r * t5.r - sum.i * t5.i, z__2.i = sum.r * t5.i + 
		    sum.i * t5.r;
#line 311 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 311 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 312 "zlarfx.f"
	    i__2 = j * c_dim1 + 6;
#line 312 "zlarfx.f"
	    i__3 = j * c_dim1 + 6;
#line 312 "zlarfx.f"
	    z__2.r = sum.r * t6.r - sum.i * t6.i, z__2.i = sum.r * t6.i + 
		    sum.i * t6.r;
#line 312 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 312 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 313 "zlarfx.f"
	    i__2 = j * c_dim1 + 7;
#line 313 "zlarfx.f"
	    i__3 = j * c_dim1 + 7;
#line 313 "zlarfx.f"
	    z__2.r = sum.r * t7.r - sum.i * t7.i, z__2.i = sum.r * t7.i + 
		    sum.i * t7.r;
#line 313 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 313 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 314 "zlarfx.f"
/* L140: */
#line 314 "zlarfx.f"
	}
#line 315 "zlarfx.f"
	goto L410;
#line 316 "zlarfx.f"
L150:

/*        Special code for 8 x 8 Householder */

#line 320 "zlarfx.f"
	d_cnjg(&z__1, &v[1]);
#line 320 "zlarfx.f"
	v1.r = z__1.r, v1.i = z__1.i;
#line 321 "zlarfx.f"
	d_cnjg(&z__2, &v1);
#line 321 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 321 "zlarfx.f"
	t1.r = z__1.r, t1.i = z__1.i;
#line 322 "zlarfx.f"
	d_cnjg(&z__1, &v[2]);
#line 322 "zlarfx.f"
	v2.r = z__1.r, v2.i = z__1.i;
#line 323 "zlarfx.f"
	d_cnjg(&z__2, &v2);
#line 323 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 323 "zlarfx.f"
	t2.r = z__1.r, t2.i = z__1.i;
#line 324 "zlarfx.f"
	d_cnjg(&z__1, &v[3]);
#line 324 "zlarfx.f"
	v3.r = z__1.r, v3.i = z__1.i;
#line 325 "zlarfx.f"
	d_cnjg(&z__2, &v3);
#line 325 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 325 "zlarfx.f"
	t3.r = z__1.r, t3.i = z__1.i;
#line 326 "zlarfx.f"
	d_cnjg(&z__1, &v[4]);
#line 326 "zlarfx.f"
	v4.r = z__1.r, v4.i = z__1.i;
#line 327 "zlarfx.f"
	d_cnjg(&z__2, &v4);
#line 327 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 327 "zlarfx.f"
	t4.r = z__1.r, t4.i = z__1.i;
#line 328 "zlarfx.f"
	d_cnjg(&z__1, &v[5]);
#line 328 "zlarfx.f"
	v5.r = z__1.r, v5.i = z__1.i;
#line 329 "zlarfx.f"
	d_cnjg(&z__2, &v5);
#line 329 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 329 "zlarfx.f"
	t5.r = z__1.r, t5.i = z__1.i;
#line 330 "zlarfx.f"
	d_cnjg(&z__1, &v[6]);
#line 330 "zlarfx.f"
	v6.r = z__1.r, v6.i = z__1.i;
#line 331 "zlarfx.f"
	d_cnjg(&z__2, &v6);
#line 331 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 331 "zlarfx.f"
	t6.r = z__1.r, t6.i = z__1.i;
#line 332 "zlarfx.f"
	d_cnjg(&z__1, &v[7]);
#line 332 "zlarfx.f"
	v7.r = z__1.r, v7.i = z__1.i;
#line 333 "zlarfx.f"
	d_cnjg(&z__2, &v7);
#line 333 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 333 "zlarfx.f"
	t7.r = z__1.r, t7.i = z__1.i;
#line 334 "zlarfx.f"
	d_cnjg(&z__1, &v[8]);
#line 334 "zlarfx.f"
	v8.r = z__1.r, v8.i = z__1.i;
#line 335 "zlarfx.f"
	d_cnjg(&z__2, &v8);
#line 335 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 335 "zlarfx.f"
	t8.r = z__1.r, t8.i = z__1.i;
#line 336 "zlarfx.f"
	i__1 = *n;
#line 336 "zlarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 337 "zlarfx.f"
	    i__2 = j * c_dim1 + 1;
#line 337 "zlarfx.f"
	    z__8.r = v1.r * c__[i__2].r - v1.i * c__[i__2].i, z__8.i = v1.r * 
		    c__[i__2].i + v1.i * c__[i__2].r;
#line 337 "zlarfx.f"
	    i__3 = j * c_dim1 + 2;
#line 337 "zlarfx.f"
	    z__9.r = v2.r * c__[i__3].r - v2.i * c__[i__3].i, z__9.i = v2.r * 
		    c__[i__3].i + v2.i * c__[i__3].r;
#line 337 "zlarfx.f"
	    z__7.r = z__8.r + z__9.r, z__7.i = z__8.i + z__9.i;
#line 337 "zlarfx.f"
	    i__4 = j * c_dim1 + 3;
#line 337 "zlarfx.f"
	    z__10.r = v3.r * c__[i__4].r - v3.i * c__[i__4].i, z__10.i = v3.r 
		    * c__[i__4].i + v3.i * c__[i__4].r;
#line 337 "zlarfx.f"
	    z__6.r = z__7.r + z__10.r, z__6.i = z__7.i + z__10.i;
#line 337 "zlarfx.f"
	    i__5 = j * c_dim1 + 4;
#line 337 "zlarfx.f"
	    z__11.r = v4.r * c__[i__5].r - v4.i * c__[i__5].i, z__11.i = v4.r 
		    * c__[i__5].i + v4.i * c__[i__5].r;
#line 337 "zlarfx.f"
	    z__5.r = z__6.r + z__11.r, z__5.i = z__6.i + z__11.i;
#line 337 "zlarfx.f"
	    i__6 = j * c_dim1 + 5;
#line 337 "zlarfx.f"
	    z__12.r = v5.r * c__[i__6].r - v5.i * c__[i__6].i, z__12.i = v5.r 
		    * c__[i__6].i + v5.i * c__[i__6].r;
#line 337 "zlarfx.f"
	    z__4.r = z__5.r + z__12.r, z__4.i = z__5.i + z__12.i;
#line 337 "zlarfx.f"
	    i__7 = j * c_dim1 + 6;
#line 337 "zlarfx.f"
	    z__13.r = v6.r * c__[i__7].r - v6.i * c__[i__7].i, z__13.i = v6.r 
		    * c__[i__7].i + v6.i * c__[i__7].r;
#line 337 "zlarfx.f"
	    z__3.r = z__4.r + z__13.r, z__3.i = z__4.i + z__13.i;
#line 337 "zlarfx.f"
	    i__8 = j * c_dim1 + 7;
#line 337 "zlarfx.f"
	    z__14.r = v7.r * c__[i__8].r - v7.i * c__[i__8].i, z__14.i = v7.r 
		    * c__[i__8].i + v7.i * c__[i__8].r;
#line 337 "zlarfx.f"
	    z__2.r = z__3.r + z__14.r, z__2.i = z__3.i + z__14.i;
#line 337 "zlarfx.f"
	    i__9 = j * c_dim1 + 8;
#line 337 "zlarfx.f"
	    z__15.r = v8.r * c__[i__9].r - v8.i * c__[i__9].i, z__15.i = v8.r 
		    * c__[i__9].i + v8.i * c__[i__9].r;
#line 337 "zlarfx.f"
	    z__1.r = z__2.r + z__15.r, z__1.i = z__2.i + z__15.i;
#line 337 "zlarfx.f"
	    sum.r = z__1.r, sum.i = z__1.i;
#line 340 "zlarfx.f"
	    i__2 = j * c_dim1 + 1;
#line 340 "zlarfx.f"
	    i__3 = j * c_dim1 + 1;
#line 340 "zlarfx.f"
	    z__2.r = sum.r * t1.r - sum.i * t1.i, z__2.i = sum.r * t1.i + 
		    sum.i * t1.r;
#line 340 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 340 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 341 "zlarfx.f"
	    i__2 = j * c_dim1 + 2;
#line 341 "zlarfx.f"
	    i__3 = j * c_dim1 + 2;
#line 341 "zlarfx.f"
	    z__2.r = sum.r * t2.r - sum.i * t2.i, z__2.i = sum.r * t2.i + 
		    sum.i * t2.r;
#line 341 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 341 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 342 "zlarfx.f"
	    i__2 = j * c_dim1 + 3;
#line 342 "zlarfx.f"
	    i__3 = j * c_dim1 + 3;
#line 342 "zlarfx.f"
	    z__2.r = sum.r * t3.r - sum.i * t3.i, z__2.i = sum.r * t3.i + 
		    sum.i * t3.r;
#line 342 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 342 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 343 "zlarfx.f"
	    i__2 = j * c_dim1 + 4;
#line 343 "zlarfx.f"
	    i__3 = j * c_dim1 + 4;
#line 343 "zlarfx.f"
	    z__2.r = sum.r * t4.r - sum.i * t4.i, z__2.i = sum.r * t4.i + 
		    sum.i * t4.r;
#line 343 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 343 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 344 "zlarfx.f"
	    i__2 = j * c_dim1 + 5;
#line 344 "zlarfx.f"
	    i__3 = j * c_dim1 + 5;
#line 344 "zlarfx.f"
	    z__2.r = sum.r * t5.r - sum.i * t5.i, z__2.i = sum.r * t5.i + 
		    sum.i * t5.r;
#line 344 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 344 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 345 "zlarfx.f"
	    i__2 = j * c_dim1 + 6;
#line 345 "zlarfx.f"
	    i__3 = j * c_dim1 + 6;
#line 345 "zlarfx.f"
	    z__2.r = sum.r * t6.r - sum.i * t6.i, z__2.i = sum.r * t6.i + 
		    sum.i * t6.r;
#line 345 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 345 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 346 "zlarfx.f"
	    i__2 = j * c_dim1 + 7;
#line 346 "zlarfx.f"
	    i__3 = j * c_dim1 + 7;
#line 346 "zlarfx.f"
	    z__2.r = sum.r * t7.r - sum.i * t7.i, z__2.i = sum.r * t7.i + 
		    sum.i * t7.r;
#line 346 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 346 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 347 "zlarfx.f"
	    i__2 = j * c_dim1 + 8;
#line 347 "zlarfx.f"
	    i__3 = j * c_dim1 + 8;
#line 347 "zlarfx.f"
	    z__2.r = sum.r * t8.r - sum.i * t8.i, z__2.i = sum.r * t8.i + 
		    sum.i * t8.r;
#line 347 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 347 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 348 "zlarfx.f"
/* L160: */
#line 348 "zlarfx.f"
	}
#line 349 "zlarfx.f"
	goto L410;
#line 350 "zlarfx.f"
L170:

/*        Special code for 9 x 9 Householder */

#line 354 "zlarfx.f"
	d_cnjg(&z__1, &v[1]);
#line 354 "zlarfx.f"
	v1.r = z__1.r, v1.i = z__1.i;
#line 355 "zlarfx.f"
	d_cnjg(&z__2, &v1);
#line 355 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 355 "zlarfx.f"
	t1.r = z__1.r, t1.i = z__1.i;
#line 356 "zlarfx.f"
	d_cnjg(&z__1, &v[2]);
#line 356 "zlarfx.f"
	v2.r = z__1.r, v2.i = z__1.i;
#line 357 "zlarfx.f"
	d_cnjg(&z__2, &v2);
#line 357 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 357 "zlarfx.f"
	t2.r = z__1.r, t2.i = z__1.i;
#line 358 "zlarfx.f"
	d_cnjg(&z__1, &v[3]);
#line 358 "zlarfx.f"
	v3.r = z__1.r, v3.i = z__1.i;
#line 359 "zlarfx.f"
	d_cnjg(&z__2, &v3);
#line 359 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 359 "zlarfx.f"
	t3.r = z__1.r, t3.i = z__1.i;
#line 360 "zlarfx.f"
	d_cnjg(&z__1, &v[4]);
#line 360 "zlarfx.f"
	v4.r = z__1.r, v4.i = z__1.i;
#line 361 "zlarfx.f"
	d_cnjg(&z__2, &v4);
#line 361 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 361 "zlarfx.f"
	t4.r = z__1.r, t4.i = z__1.i;
#line 362 "zlarfx.f"
	d_cnjg(&z__1, &v[5]);
#line 362 "zlarfx.f"
	v5.r = z__1.r, v5.i = z__1.i;
#line 363 "zlarfx.f"
	d_cnjg(&z__2, &v5);
#line 363 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 363 "zlarfx.f"
	t5.r = z__1.r, t5.i = z__1.i;
#line 364 "zlarfx.f"
	d_cnjg(&z__1, &v[6]);
#line 364 "zlarfx.f"
	v6.r = z__1.r, v6.i = z__1.i;
#line 365 "zlarfx.f"
	d_cnjg(&z__2, &v6);
#line 365 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 365 "zlarfx.f"
	t6.r = z__1.r, t6.i = z__1.i;
#line 366 "zlarfx.f"
	d_cnjg(&z__1, &v[7]);
#line 366 "zlarfx.f"
	v7.r = z__1.r, v7.i = z__1.i;
#line 367 "zlarfx.f"
	d_cnjg(&z__2, &v7);
#line 367 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 367 "zlarfx.f"
	t7.r = z__1.r, t7.i = z__1.i;
#line 368 "zlarfx.f"
	d_cnjg(&z__1, &v[8]);
#line 368 "zlarfx.f"
	v8.r = z__1.r, v8.i = z__1.i;
#line 369 "zlarfx.f"
	d_cnjg(&z__2, &v8);
#line 369 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 369 "zlarfx.f"
	t8.r = z__1.r, t8.i = z__1.i;
#line 370 "zlarfx.f"
	d_cnjg(&z__1, &v[9]);
#line 370 "zlarfx.f"
	v9.r = z__1.r, v9.i = z__1.i;
#line 371 "zlarfx.f"
	d_cnjg(&z__2, &v9);
#line 371 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 371 "zlarfx.f"
	t9.r = z__1.r, t9.i = z__1.i;
#line 372 "zlarfx.f"
	i__1 = *n;
#line 372 "zlarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 373 "zlarfx.f"
	    i__2 = j * c_dim1 + 1;
#line 373 "zlarfx.f"
	    z__9.r = v1.r * c__[i__2].r - v1.i * c__[i__2].i, z__9.i = v1.r * 
		    c__[i__2].i + v1.i * c__[i__2].r;
#line 373 "zlarfx.f"
	    i__3 = j * c_dim1 + 2;
#line 373 "zlarfx.f"
	    z__10.r = v2.r * c__[i__3].r - v2.i * c__[i__3].i, z__10.i = v2.r 
		    * c__[i__3].i + v2.i * c__[i__3].r;
#line 373 "zlarfx.f"
	    z__8.r = z__9.r + z__10.r, z__8.i = z__9.i + z__10.i;
#line 373 "zlarfx.f"
	    i__4 = j * c_dim1 + 3;
#line 373 "zlarfx.f"
	    z__11.r = v3.r * c__[i__4].r - v3.i * c__[i__4].i, z__11.i = v3.r 
		    * c__[i__4].i + v3.i * c__[i__4].r;
#line 373 "zlarfx.f"
	    z__7.r = z__8.r + z__11.r, z__7.i = z__8.i + z__11.i;
#line 373 "zlarfx.f"
	    i__5 = j * c_dim1 + 4;
#line 373 "zlarfx.f"
	    z__12.r = v4.r * c__[i__5].r - v4.i * c__[i__5].i, z__12.i = v4.r 
		    * c__[i__5].i + v4.i * c__[i__5].r;
#line 373 "zlarfx.f"
	    z__6.r = z__7.r + z__12.r, z__6.i = z__7.i + z__12.i;
#line 373 "zlarfx.f"
	    i__6 = j * c_dim1 + 5;
#line 373 "zlarfx.f"
	    z__13.r = v5.r * c__[i__6].r - v5.i * c__[i__6].i, z__13.i = v5.r 
		    * c__[i__6].i + v5.i * c__[i__6].r;
#line 373 "zlarfx.f"
	    z__5.r = z__6.r + z__13.r, z__5.i = z__6.i + z__13.i;
#line 373 "zlarfx.f"
	    i__7 = j * c_dim1 + 6;
#line 373 "zlarfx.f"
	    z__14.r = v6.r * c__[i__7].r - v6.i * c__[i__7].i, z__14.i = v6.r 
		    * c__[i__7].i + v6.i * c__[i__7].r;
#line 373 "zlarfx.f"
	    z__4.r = z__5.r + z__14.r, z__4.i = z__5.i + z__14.i;
#line 373 "zlarfx.f"
	    i__8 = j * c_dim1 + 7;
#line 373 "zlarfx.f"
	    z__15.r = v7.r * c__[i__8].r - v7.i * c__[i__8].i, z__15.i = v7.r 
		    * c__[i__8].i + v7.i * c__[i__8].r;
#line 373 "zlarfx.f"
	    z__3.r = z__4.r + z__15.r, z__3.i = z__4.i + z__15.i;
#line 373 "zlarfx.f"
	    i__9 = j * c_dim1 + 8;
#line 373 "zlarfx.f"
	    z__16.r = v8.r * c__[i__9].r - v8.i * c__[i__9].i, z__16.i = v8.r 
		    * c__[i__9].i + v8.i * c__[i__9].r;
#line 373 "zlarfx.f"
	    z__2.r = z__3.r + z__16.r, z__2.i = z__3.i + z__16.i;
#line 373 "zlarfx.f"
	    i__10 = j * c_dim1 + 9;
#line 373 "zlarfx.f"
	    z__17.r = v9.r * c__[i__10].r - v9.i * c__[i__10].i, z__17.i = 
		    v9.r * c__[i__10].i + v9.i * c__[i__10].r;
#line 373 "zlarfx.f"
	    z__1.r = z__2.r + z__17.r, z__1.i = z__2.i + z__17.i;
#line 373 "zlarfx.f"
	    sum.r = z__1.r, sum.i = z__1.i;
#line 376 "zlarfx.f"
	    i__2 = j * c_dim1 + 1;
#line 376 "zlarfx.f"
	    i__3 = j * c_dim1 + 1;
#line 376 "zlarfx.f"
	    z__2.r = sum.r * t1.r - sum.i * t1.i, z__2.i = sum.r * t1.i + 
		    sum.i * t1.r;
#line 376 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 376 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 377 "zlarfx.f"
	    i__2 = j * c_dim1 + 2;
#line 377 "zlarfx.f"
	    i__3 = j * c_dim1 + 2;
#line 377 "zlarfx.f"
	    z__2.r = sum.r * t2.r - sum.i * t2.i, z__2.i = sum.r * t2.i + 
		    sum.i * t2.r;
#line 377 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 377 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 378 "zlarfx.f"
	    i__2 = j * c_dim1 + 3;
#line 378 "zlarfx.f"
	    i__3 = j * c_dim1 + 3;
#line 378 "zlarfx.f"
	    z__2.r = sum.r * t3.r - sum.i * t3.i, z__2.i = sum.r * t3.i + 
		    sum.i * t3.r;
#line 378 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 378 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 379 "zlarfx.f"
	    i__2 = j * c_dim1 + 4;
#line 379 "zlarfx.f"
	    i__3 = j * c_dim1 + 4;
#line 379 "zlarfx.f"
	    z__2.r = sum.r * t4.r - sum.i * t4.i, z__2.i = sum.r * t4.i + 
		    sum.i * t4.r;
#line 379 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 379 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 380 "zlarfx.f"
	    i__2 = j * c_dim1 + 5;
#line 380 "zlarfx.f"
	    i__3 = j * c_dim1 + 5;
#line 380 "zlarfx.f"
	    z__2.r = sum.r * t5.r - sum.i * t5.i, z__2.i = sum.r * t5.i + 
		    sum.i * t5.r;
#line 380 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 380 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 381 "zlarfx.f"
	    i__2 = j * c_dim1 + 6;
#line 381 "zlarfx.f"
	    i__3 = j * c_dim1 + 6;
#line 381 "zlarfx.f"
	    z__2.r = sum.r * t6.r - sum.i * t6.i, z__2.i = sum.r * t6.i + 
		    sum.i * t6.r;
#line 381 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 381 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 382 "zlarfx.f"
	    i__2 = j * c_dim1 + 7;
#line 382 "zlarfx.f"
	    i__3 = j * c_dim1 + 7;
#line 382 "zlarfx.f"
	    z__2.r = sum.r * t7.r - sum.i * t7.i, z__2.i = sum.r * t7.i + 
		    sum.i * t7.r;
#line 382 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 382 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 383 "zlarfx.f"
	    i__2 = j * c_dim1 + 8;
#line 383 "zlarfx.f"
	    i__3 = j * c_dim1 + 8;
#line 383 "zlarfx.f"
	    z__2.r = sum.r * t8.r - sum.i * t8.i, z__2.i = sum.r * t8.i + 
		    sum.i * t8.r;
#line 383 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 383 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 384 "zlarfx.f"
	    i__2 = j * c_dim1 + 9;
#line 384 "zlarfx.f"
	    i__3 = j * c_dim1 + 9;
#line 384 "zlarfx.f"
	    z__2.r = sum.r * t9.r - sum.i * t9.i, z__2.i = sum.r * t9.i + 
		    sum.i * t9.r;
#line 384 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 384 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 385 "zlarfx.f"
/* L180: */
#line 385 "zlarfx.f"
	}
#line 386 "zlarfx.f"
	goto L410;
#line 387 "zlarfx.f"
L190:

/*        Special code for 10 x 10 Householder */

#line 391 "zlarfx.f"
	d_cnjg(&z__1, &v[1]);
#line 391 "zlarfx.f"
	v1.r = z__1.r, v1.i = z__1.i;
#line 392 "zlarfx.f"
	d_cnjg(&z__2, &v1);
#line 392 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 392 "zlarfx.f"
	t1.r = z__1.r, t1.i = z__1.i;
#line 393 "zlarfx.f"
	d_cnjg(&z__1, &v[2]);
#line 393 "zlarfx.f"
	v2.r = z__1.r, v2.i = z__1.i;
#line 394 "zlarfx.f"
	d_cnjg(&z__2, &v2);
#line 394 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 394 "zlarfx.f"
	t2.r = z__1.r, t2.i = z__1.i;
#line 395 "zlarfx.f"
	d_cnjg(&z__1, &v[3]);
#line 395 "zlarfx.f"
	v3.r = z__1.r, v3.i = z__1.i;
#line 396 "zlarfx.f"
	d_cnjg(&z__2, &v3);
#line 396 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 396 "zlarfx.f"
	t3.r = z__1.r, t3.i = z__1.i;
#line 397 "zlarfx.f"
	d_cnjg(&z__1, &v[4]);
#line 397 "zlarfx.f"
	v4.r = z__1.r, v4.i = z__1.i;
#line 398 "zlarfx.f"
	d_cnjg(&z__2, &v4);
#line 398 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 398 "zlarfx.f"
	t4.r = z__1.r, t4.i = z__1.i;
#line 399 "zlarfx.f"
	d_cnjg(&z__1, &v[5]);
#line 399 "zlarfx.f"
	v5.r = z__1.r, v5.i = z__1.i;
#line 400 "zlarfx.f"
	d_cnjg(&z__2, &v5);
#line 400 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 400 "zlarfx.f"
	t5.r = z__1.r, t5.i = z__1.i;
#line 401 "zlarfx.f"
	d_cnjg(&z__1, &v[6]);
#line 401 "zlarfx.f"
	v6.r = z__1.r, v6.i = z__1.i;
#line 402 "zlarfx.f"
	d_cnjg(&z__2, &v6);
#line 402 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 402 "zlarfx.f"
	t6.r = z__1.r, t6.i = z__1.i;
#line 403 "zlarfx.f"
	d_cnjg(&z__1, &v[7]);
#line 403 "zlarfx.f"
	v7.r = z__1.r, v7.i = z__1.i;
#line 404 "zlarfx.f"
	d_cnjg(&z__2, &v7);
#line 404 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 404 "zlarfx.f"
	t7.r = z__1.r, t7.i = z__1.i;
#line 405 "zlarfx.f"
	d_cnjg(&z__1, &v[8]);
#line 405 "zlarfx.f"
	v8.r = z__1.r, v8.i = z__1.i;
#line 406 "zlarfx.f"
	d_cnjg(&z__2, &v8);
#line 406 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 406 "zlarfx.f"
	t8.r = z__1.r, t8.i = z__1.i;
#line 407 "zlarfx.f"
	d_cnjg(&z__1, &v[9]);
#line 407 "zlarfx.f"
	v9.r = z__1.r, v9.i = z__1.i;
#line 408 "zlarfx.f"
	d_cnjg(&z__2, &v9);
#line 408 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 408 "zlarfx.f"
	t9.r = z__1.r, t9.i = z__1.i;
#line 409 "zlarfx.f"
	d_cnjg(&z__1, &v[10]);
#line 409 "zlarfx.f"
	v10.r = z__1.r, v10.i = z__1.i;
#line 410 "zlarfx.f"
	d_cnjg(&z__2, &v10);
#line 410 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 410 "zlarfx.f"
	t10.r = z__1.r, t10.i = z__1.i;
#line 411 "zlarfx.f"
	i__1 = *n;
#line 411 "zlarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 412 "zlarfx.f"
	    i__2 = j * c_dim1 + 1;
#line 412 "zlarfx.f"
	    z__10.r = v1.r * c__[i__2].r - v1.i * c__[i__2].i, z__10.i = v1.r 
		    * c__[i__2].i + v1.i * c__[i__2].r;
#line 412 "zlarfx.f"
	    i__3 = j * c_dim1 + 2;
#line 412 "zlarfx.f"
	    z__11.r = v2.r * c__[i__3].r - v2.i * c__[i__3].i, z__11.i = v2.r 
		    * c__[i__3].i + v2.i * c__[i__3].r;
#line 412 "zlarfx.f"
	    z__9.r = z__10.r + z__11.r, z__9.i = z__10.i + z__11.i;
#line 412 "zlarfx.f"
	    i__4 = j * c_dim1 + 3;
#line 412 "zlarfx.f"
	    z__12.r = v3.r * c__[i__4].r - v3.i * c__[i__4].i, z__12.i = v3.r 
		    * c__[i__4].i + v3.i * c__[i__4].r;
#line 412 "zlarfx.f"
	    z__8.r = z__9.r + z__12.r, z__8.i = z__9.i + z__12.i;
#line 412 "zlarfx.f"
	    i__5 = j * c_dim1 + 4;
#line 412 "zlarfx.f"
	    z__13.r = v4.r * c__[i__5].r - v4.i * c__[i__5].i, z__13.i = v4.r 
		    * c__[i__5].i + v4.i * c__[i__5].r;
#line 412 "zlarfx.f"
	    z__7.r = z__8.r + z__13.r, z__7.i = z__8.i + z__13.i;
#line 412 "zlarfx.f"
	    i__6 = j * c_dim1 + 5;
#line 412 "zlarfx.f"
	    z__14.r = v5.r * c__[i__6].r - v5.i * c__[i__6].i, z__14.i = v5.r 
		    * c__[i__6].i + v5.i * c__[i__6].r;
#line 412 "zlarfx.f"
	    z__6.r = z__7.r + z__14.r, z__6.i = z__7.i + z__14.i;
#line 412 "zlarfx.f"
	    i__7 = j * c_dim1 + 6;
#line 412 "zlarfx.f"
	    z__15.r = v6.r * c__[i__7].r - v6.i * c__[i__7].i, z__15.i = v6.r 
		    * c__[i__7].i + v6.i * c__[i__7].r;
#line 412 "zlarfx.f"
	    z__5.r = z__6.r + z__15.r, z__5.i = z__6.i + z__15.i;
#line 412 "zlarfx.f"
	    i__8 = j * c_dim1 + 7;
#line 412 "zlarfx.f"
	    z__16.r = v7.r * c__[i__8].r - v7.i * c__[i__8].i, z__16.i = v7.r 
		    * c__[i__8].i + v7.i * c__[i__8].r;
#line 412 "zlarfx.f"
	    z__4.r = z__5.r + z__16.r, z__4.i = z__5.i + z__16.i;
#line 412 "zlarfx.f"
	    i__9 = j * c_dim1 + 8;
#line 412 "zlarfx.f"
	    z__17.r = v8.r * c__[i__9].r - v8.i * c__[i__9].i, z__17.i = v8.r 
		    * c__[i__9].i + v8.i * c__[i__9].r;
#line 412 "zlarfx.f"
	    z__3.r = z__4.r + z__17.r, z__3.i = z__4.i + z__17.i;
#line 412 "zlarfx.f"
	    i__10 = j * c_dim1 + 9;
#line 412 "zlarfx.f"
	    z__18.r = v9.r * c__[i__10].r - v9.i * c__[i__10].i, z__18.i = 
		    v9.r * c__[i__10].i + v9.i * c__[i__10].r;
#line 412 "zlarfx.f"
	    z__2.r = z__3.r + z__18.r, z__2.i = z__3.i + z__18.i;
#line 412 "zlarfx.f"
	    i__11 = j * c_dim1 + 10;
#line 412 "zlarfx.f"
	    z__19.r = v10.r * c__[i__11].r - v10.i * c__[i__11].i, z__19.i = 
		    v10.r * c__[i__11].i + v10.i * c__[i__11].r;
#line 412 "zlarfx.f"
	    z__1.r = z__2.r + z__19.r, z__1.i = z__2.i + z__19.i;
#line 412 "zlarfx.f"
	    sum.r = z__1.r, sum.i = z__1.i;
#line 416 "zlarfx.f"
	    i__2 = j * c_dim1 + 1;
#line 416 "zlarfx.f"
	    i__3 = j * c_dim1 + 1;
#line 416 "zlarfx.f"
	    z__2.r = sum.r * t1.r - sum.i * t1.i, z__2.i = sum.r * t1.i + 
		    sum.i * t1.r;
#line 416 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 416 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 417 "zlarfx.f"
	    i__2 = j * c_dim1 + 2;
#line 417 "zlarfx.f"
	    i__3 = j * c_dim1 + 2;
#line 417 "zlarfx.f"
	    z__2.r = sum.r * t2.r - sum.i * t2.i, z__2.i = sum.r * t2.i + 
		    sum.i * t2.r;
#line 417 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 417 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 418 "zlarfx.f"
	    i__2 = j * c_dim1 + 3;
#line 418 "zlarfx.f"
	    i__3 = j * c_dim1 + 3;
#line 418 "zlarfx.f"
	    z__2.r = sum.r * t3.r - sum.i * t3.i, z__2.i = sum.r * t3.i + 
		    sum.i * t3.r;
#line 418 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 418 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 419 "zlarfx.f"
	    i__2 = j * c_dim1 + 4;
#line 419 "zlarfx.f"
	    i__3 = j * c_dim1 + 4;
#line 419 "zlarfx.f"
	    z__2.r = sum.r * t4.r - sum.i * t4.i, z__2.i = sum.r * t4.i + 
		    sum.i * t4.r;
#line 419 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 419 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 420 "zlarfx.f"
	    i__2 = j * c_dim1 + 5;
#line 420 "zlarfx.f"
	    i__3 = j * c_dim1 + 5;
#line 420 "zlarfx.f"
	    z__2.r = sum.r * t5.r - sum.i * t5.i, z__2.i = sum.r * t5.i + 
		    sum.i * t5.r;
#line 420 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 420 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 421 "zlarfx.f"
	    i__2 = j * c_dim1 + 6;
#line 421 "zlarfx.f"
	    i__3 = j * c_dim1 + 6;
#line 421 "zlarfx.f"
	    z__2.r = sum.r * t6.r - sum.i * t6.i, z__2.i = sum.r * t6.i + 
		    sum.i * t6.r;
#line 421 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 421 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 422 "zlarfx.f"
	    i__2 = j * c_dim1 + 7;
#line 422 "zlarfx.f"
	    i__3 = j * c_dim1 + 7;
#line 422 "zlarfx.f"
	    z__2.r = sum.r * t7.r - sum.i * t7.i, z__2.i = sum.r * t7.i + 
		    sum.i * t7.r;
#line 422 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 422 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 423 "zlarfx.f"
	    i__2 = j * c_dim1 + 8;
#line 423 "zlarfx.f"
	    i__3 = j * c_dim1 + 8;
#line 423 "zlarfx.f"
	    z__2.r = sum.r * t8.r - sum.i * t8.i, z__2.i = sum.r * t8.i + 
		    sum.i * t8.r;
#line 423 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 423 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 424 "zlarfx.f"
	    i__2 = j * c_dim1 + 9;
#line 424 "zlarfx.f"
	    i__3 = j * c_dim1 + 9;
#line 424 "zlarfx.f"
	    z__2.r = sum.r * t9.r - sum.i * t9.i, z__2.i = sum.r * t9.i + 
		    sum.i * t9.r;
#line 424 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 424 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 425 "zlarfx.f"
	    i__2 = j * c_dim1 + 10;
#line 425 "zlarfx.f"
	    i__3 = j * c_dim1 + 10;
#line 425 "zlarfx.f"
	    z__2.r = sum.r * t10.r - sum.i * t10.i, z__2.i = sum.r * t10.i + 
		    sum.i * t10.r;
#line 425 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 425 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 426 "zlarfx.f"
/* L200: */
#line 426 "zlarfx.f"
	}
#line 427 "zlarfx.f"
	goto L410;
#line 428 "zlarfx.f"
    } else {

/*        Form  C * H, where H has order n. */

#line 432 "zlarfx.f"
	switch (*n) {
#line 432 "zlarfx.f"
	    case 1:  goto L210;
#line 432 "zlarfx.f"
	    case 2:  goto L230;
#line 432 "zlarfx.f"
	    case 3:  goto L250;
#line 432 "zlarfx.f"
	    case 4:  goto L270;
#line 432 "zlarfx.f"
	    case 5:  goto L290;
#line 432 "zlarfx.f"
	    case 6:  goto L310;
#line 432 "zlarfx.f"
	    case 7:  goto L330;
#line 432 "zlarfx.f"
	    case 8:  goto L350;
#line 432 "zlarfx.f"
	    case 9:  goto L370;
#line 432 "zlarfx.f"
	    case 10:  goto L390;
#line 432 "zlarfx.f"
	}

/*        Code for general N */

#line 437 "zlarfx.f"
	zlarf_(side, m, n, &v[1], &c__1, tau, &c__[c_offset], ldc, &work[1], (
		ftnlen)1);
#line 438 "zlarfx.f"
	goto L410;
#line 439 "zlarfx.f"
L210:

/*        Special code for 1 x 1 Householder */

#line 443 "zlarfx.f"
	z__3.r = tau->r * v[1].r - tau->i * v[1].i, z__3.i = tau->r * v[1].i 
		+ tau->i * v[1].r;
#line 443 "zlarfx.f"
	d_cnjg(&z__4, &v[1]);
#line 443 "zlarfx.f"
	z__2.r = z__3.r * z__4.r - z__3.i * z__4.i, z__2.i = z__3.r * z__4.i 
		+ z__3.i * z__4.r;
#line 443 "zlarfx.f"
	z__1.r = 1. - z__2.r, z__1.i = 0. - z__2.i;
#line 443 "zlarfx.f"
	t1.r = z__1.r, t1.i = z__1.i;
#line 444 "zlarfx.f"
	i__1 = *m;
#line 444 "zlarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 445 "zlarfx.f"
	    i__2 = j + c_dim1;
#line 445 "zlarfx.f"
	    i__3 = j + c_dim1;
#line 445 "zlarfx.f"
	    z__1.r = t1.r * c__[i__3].r - t1.i * c__[i__3].i, z__1.i = t1.r * 
		    c__[i__3].i + t1.i * c__[i__3].r;
#line 445 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 446 "zlarfx.f"
/* L220: */
#line 446 "zlarfx.f"
	}
#line 447 "zlarfx.f"
	goto L410;
#line 448 "zlarfx.f"
L230:

/*        Special code for 2 x 2 Householder */

#line 452 "zlarfx.f"
	v1.r = v[1].r, v1.i = v[1].i;
#line 453 "zlarfx.f"
	d_cnjg(&z__2, &v1);
#line 453 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 453 "zlarfx.f"
	t1.r = z__1.r, t1.i = z__1.i;
#line 454 "zlarfx.f"
	v2.r = v[2].r, v2.i = v[2].i;
#line 455 "zlarfx.f"
	d_cnjg(&z__2, &v2);
#line 455 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 455 "zlarfx.f"
	t2.r = z__1.r, t2.i = z__1.i;
#line 456 "zlarfx.f"
	i__1 = *m;
#line 456 "zlarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 457 "zlarfx.f"
	    i__2 = j + c_dim1;
#line 457 "zlarfx.f"
	    z__2.r = v1.r * c__[i__2].r - v1.i * c__[i__2].i, z__2.i = v1.r * 
		    c__[i__2].i + v1.i * c__[i__2].r;
#line 457 "zlarfx.f"
	    i__3 = j + (c_dim1 << 1);
#line 457 "zlarfx.f"
	    z__3.r = v2.r * c__[i__3].r - v2.i * c__[i__3].i, z__3.i = v2.r * 
		    c__[i__3].i + v2.i * c__[i__3].r;
#line 457 "zlarfx.f"
	    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 457 "zlarfx.f"
	    sum.r = z__1.r, sum.i = z__1.i;
#line 458 "zlarfx.f"
	    i__2 = j + c_dim1;
#line 458 "zlarfx.f"
	    i__3 = j + c_dim1;
#line 458 "zlarfx.f"
	    z__2.r = sum.r * t1.r - sum.i * t1.i, z__2.i = sum.r * t1.i + 
		    sum.i * t1.r;
#line 458 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 458 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 459 "zlarfx.f"
	    i__2 = j + (c_dim1 << 1);
#line 459 "zlarfx.f"
	    i__3 = j + (c_dim1 << 1);
#line 459 "zlarfx.f"
	    z__2.r = sum.r * t2.r - sum.i * t2.i, z__2.i = sum.r * t2.i + 
		    sum.i * t2.r;
#line 459 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 459 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 460 "zlarfx.f"
/* L240: */
#line 460 "zlarfx.f"
	}
#line 461 "zlarfx.f"
	goto L410;
#line 462 "zlarfx.f"
L250:

/*        Special code for 3 x 3 Householder */

#line 466 "zlarfx.f"
	v1.r = v[1].r, v1.i = v[1].i;
#line 467 "zlarfx.f"
	d_cnjg(&z__2, &v1);
#line 467 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 467 "zlarfx.f"
	t1.r = z__1.r, t1.i = z__1.i;
#line 468 "zlarfx.f"
	v2.r = v[2].r, v2.i = v[2].i;
#line 469 "zlarfx.f"
	d_cnjg(&z__2, &v2);
#line 469 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 469 "zlarfx.f"
	t2.r = z__1.r, t2.i = z__1.i;
#line 470 "zlarfx.f"
	v3.r = v[3].r, v3.i = v[3].i;
#line 471 "zlarfx.f"
	d_cnjg(&z__2, &v3);
#line 471 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 471 "zlarfx.f"
	t3.r = z__1.r, t3.i = z__1.i;
#line 472 "zlarfx.f"
	i__1 = *m;
#line 472 "zlarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 473 "zlarfx.f"
	    i__2 = j + c_dim1;
#line 473 "zlarfx.f"
	    z__3.r = v1.r * c__[i__2].r - v1.i * c__[i__2].i, z__3.i = v1.r * 
		    c__[i__2].i + v1.i * c__[i__2].r;
#line 473 "zlarfx.f"
	    i__3 = j + (c_dim1 << 1);
#line 473 "zlarfx.f"
	    z__4.r = v2.r * c__[i__3].r - v2.i * c__[i__3].i, z__4.i = v2.r * 
		    c__[i__3].i + v2.i * c__[i__3].r;
#line 473 "zlarfx.f"
	    z__2.r = z__3.r + z__4.r, z__2.i = z__3.i + z__4.i;
#line 473 "zlarfx.f"
	    i__4 = j + c_dim1 * 3;
#line 473 "zlarfx.f"
	    z__5.r = v3.r * c__[i__4].r - v3.i * c__[i__4].i, z__5.i = v3.r * 
		    c__[i__4].i + v3.i * c__[i__4].r;
#line 473 "zlarfx.f"
	    z__1.r = z__2.r + z__5.r, z__1.i = z__2.i + z__5.i;
#line 473 "zlarfx.f"
	    sum.r = z__1.r, sum.i = z__1.i;
#line 474 "zlarfx.f"
	    i__2 = j + c_dim1;
#line 474 "zlarfx.f"
	    i__3 = j + c_dim1;
#line 474 "zlarfx.f"
	    z__2.r = sum.r * t1.r - sum.i * t1.i, z__2.i = sum.r * t1.i + 
		    sum.i * t1.r;
#line 474 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 474 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 475 "zlarfx.f"
	    i__2 = j + (c_dim1 << 1);
#line 475 "zlarfx.f"
	    i__3 = j + (c_dim1 << 1);
#line 475 "zlarfx.f"
	    z__2.r = sum.r * t2.r - sum.i * t2.i, z__2.i = sum.r * t2.i + 
		    sum.i * t2.r;
#line 475 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 475 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 476 "zlarfx.f"
	    i__2 = j + c_dim1 * 3;
#line 476 "zlarfx.f"
	    i__3 = j + c_dim1 * 3;
#line 476 "zlarfx.f"
	    z__2.r = sum.r * t3.r - sum.i * t3.i, z__2.i = sum.r * t3.i + 
		    sum.i * t3.r;
#line 476 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 476 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 477 "zlarfx.f"
/* L260: */
#line 477 "zlarfx.f"
	}
#line 478 "zlarfx.f"
	goto L410;
#line 479 "zlarfx.f"
L270:

/*        Special code for 4 x 4 Householder */

#line 483 "zlarfx.f"
	v1.r = v[1].r, v1.i = v[1].i;
#line 484 "zlarfx.f"
	d_cnjg(&z__2, &v1);
#line 484 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 484 "zlarfx.f"
	t1.r = z__1.r, t1.i = z__1.i;
#line 485 "zlarfx.f"
	v2.r = v[2].r, v2.i = v[2].i;
#line 486 "zlarfx.f"
	d_cnjg(&z__2, &v2);
#line 486 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 486 "zlarfx.f"
	t2.r = z__1.r, t2.i = z__1.i;
#line 487 "zlarfx.f"
	v3.r = v[3].r, v3.i = v[3].i;
#line 488 "zlarfx.f"
	d_cnjg(&z__2, &v3);
#line 488 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 488 "zlarfx.f"
	t3.r = z__1.r, t3.i = z__1.i;
#line 489 "zlarfx.f"
	v4.r = v[4].r, v4.i = v[4].i;
#line 490 "zlarfx.f"
	d_cnjg(&z__2, &v4);
#line 490 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 490 "zlarfx.f"
	t4.r = z__1.r, t4.i = z__1.i;
#line 491 "zlarfx.f"
	i__1 = *m;
#line 491 "zlarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 492 "zlarfx.f"
	    i__2 = j + c_dim1;
#line 492 "zlarfx.f"
	    z__4.r = v1.r * c__[i__2].r - v1.i * c__[i__2].i, z__4.i = v1.r * 
		    c__[i__2].i + v1.i * c__[i__2].r;
#line 492 "zlarfx.f"
	    i__3 = j + (c_dim1 << 1);
#line 492 "zlarfx.f"
	    z__5.r = v2.r * c__[i__3].r - v2.i * c__[i__3].i, z__5.i = v2.r * 
		    c__[i__3].i + v2.i * c__[i__3].r;
#line 492 "zlarfx.f"
	    z__3.r = z__4.r + z__5.r, z__3.i = z__4.i + z__5.i;
#line 492 "zlarfx.f"
	    i__4 = j + c_dim1 * 3;
#line 492 "zlarfx.f"
	    z__6.r = v3.r * c__[i__4].r - v3.i * c__[i__4].i, z__6.i = v3.r * 
		    c__[i__4].i + v3.i * c__[i__4].r;
#line 492 "zlarfx.f"
	    z__2.r = z__3.r + z__6.r, z__2.i = z__3.i + z__6.i;
#line 492 "zlarfx.f"
	    i__5 = j + (c_dim1 << 2);
#line 492 "zlarfx.f"
	    z__7.r = v4.r * c__[i__5].r - v4.i * c__[i__5].i, z__7.i = v4.r * 
		    c__[i__5].i + v4.i * c__[i__5].r;
#line 492 "zlarfx.f"
	    z__1.r = z__2.r + z__7.r, z__1.i = z__2.i + z__7.i;
#line 492 "zlarfx.f"
	    sum.r = z__1.r, sum.i = z__1.i;
#line 494 "zlarfx.f"
	    i__2 = j + c_dim1;
#line 494 "zlarfx.f"
	    i__3 = j + c_dim1;
#line 494 "zlarfx.f"
	    z__2.r = sum.r * t1.r - sum.i * t1.i, z__2.i = sum.r * t1.i + 
		    sum.i * t1.r;
#line 494 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 494 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 495 "zlarfx.f"
	    i__2 = j + (c_dim1 << 1);
#line 495 "zlarfx.f"
	    i__3 = j + (c_dim1 << 1);
#line 495 "zlarfx.f"
	    z__2.r = sum.r * t2.r - sum.i * t2.i, z__2.i = sum.r * t2.i + 
		    sum.i * t2.r;
#line 495 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 495 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 496 "zlarfx.f"
	    i__2 = j + c_dim1 * 3;
#line 496 "zlarfx.f"
	    i__3 = j + c_dim1 * 3;
#line 496 "zlarfx.f"
	    z__2.r = sum.r * t3.r - sum.i * t3.i, z__2.i = sum.r * t3.i + 
		    sum.i * t3.r;
#line 496 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 496 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 497 "zlarfx.f"
	    i__2 = j + (c_dim1 << 2);
#line 497 "zlarfx.f"
	    i__3 = j + (c_dim1 << 2);
#line 497 "zlarfx.f"
	    z__2.r = sum.r * t4.r - sum.i * t4.i, z__2.i = sum.r * t4.i + 
		    sum.i * t4.r;
#line 497 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 497 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 498 "zlarfx.f"
/* L280: */
#line 498 "zlarfx.f"
	}
#line 499 "zlarfx.f"
	goto L410;
#line 500 "zlarfx.f"
L290:

/*        Special code for 5 x 5 Householder */

#line 504 "zlarfx.f"
	v1.r = v[1].r, v1.i = v[1].i;
#line 505 "zlarfx.f"
	d_cnjg(&z__2, &v1);
#line 505 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 505 "zlarfx.f"
	t1.r = z__1.r, t1.i = z__1.i;
#line 506 "zlarfx.f"
	v2.r = v[2].r, v2.i = v[2].i;
#line 507 "zlarfx.f"
	d_cnjg(&z__2, &v2);
#line 507 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 507 "zlarfx.f"
	t2.r = z__1.r, t2.i = z__1.i;
#line 508 "zlarfx.f"
	v3.r = v[3].r, v3.i = v[3].i;
#line 509 "zlarfx.f"
	d_cnjg(&z__2, &v3);
#line 509 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 509 "zlarfx.f"
	t3.r = z__1.r, t3.i = z__1.i;
#line 510 "zlarfx.f"
	v4.r = v[4].r, v4.i = v[4].i;
#line 511 "zlarfx.f"
	d_cnjg(&z__2, &v4);
#line 511 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 511 "zlarfx.f"
	t4.r = z__1.r, t4.i = z__1.i;
#line 512 "zlarfx.f"
	v5.r = v[5].r, v5.i = v[5].i;
#line 513 "zlarfx.f"
	d_cnjg(&z__2, &v5);
#line 513 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 513 "zlarfx.f"
	t5.r = z__1.r, t5.i = z__1.i;
#line 514 "zlarfx.f"
	i__1 = *m;
#line 514 "zlarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 515 "zlarfx.f"
	    i__2 = j + c_dim1;
#line 515 "zlarfx.f"
	    z__5.r = v1.r * c__[i__2].r - v1.i * c__[i__2].i, z__5.i = v1.r * 
		    c__[i__2].i + v1.i * c__[i__2].r;
#line 515 "zlarfx.f"
	    i__3 = j + (c_dim1 << 1);
#line 515 "zlarfx.f"
	    z__6.r = v2.r * c__[i__3].r - v2.i * c__[i__3].i, z__6.i = v2.r * 
		    c__[i__3].i + v2.i * c__[i__3].r;
#line 515 "zlarfx.f"
	    z__4.r = z__5.r + z__6.r, z__4.i = z__5.i + z__6.i;
#line 515 "zlarfx.f"
	    i__4 = j + c_dim1 * 3;
#line 515 "zlarfx.f"
	    z__7.r = v3.r * c__[i__4].r - v3.i * c__[i__4].i, z__7.i = v3.r * 
		    c__[i__4].i + v3.i * c__[i__4].r;
#line 515 "zlarfx.f"
	    z__3.r = z__4.r + z__7.r, z__3.i = z__4.i + z__7.i;
#line 515 "zlarfx.f"
	    i__5 = j + (c_dim1 << 2);
#line 515 "zlarfx.f"
	    z__8.r = v4.r * c__[i__5].r - v4.i * c__[i__5].i, z__8.i = v4.r * 
		    c__[i__5].i + v4.i * c__[i__5].r;
#line 515 "zlarfx.f"
	    z__2.r = z__3.r + z__8.r, z__2.i = z__3.i + z__8.i;
#line 515 "zlarfx.f"
	    i__6 = j + c_dim1 * 5;
#line 515 "zlarfx.f"
	    z__9.r = v5.r * c__[i__6].r - v5.i * c__[i__6].i, z__9.i = v5.r * 
		    c__[i__6].i + v5.i * c__[i__6].r;
#line 515 "zlarfx.f"
	    z__1.r = z__2.r + z__9.r, z__1.i = z__2.i + z__9.i;
#line 515 "zlarfx.f"
	    sum.r = z__1.r, sum.i = z__1.i;
#line 517 "zlarfx.f"
	    i__2 = j + c_dim1;
#line 517 "zlarfx.f"
	    i__3 = j + c_dim1;
#line 517 "zlarfx.f"
	    z__2.r = sum.r * t1.r - sum.i * t1.i, z__2.i = sum.r * t1.i + 
		    sum.i * t1.r;
#line 517 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 517 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 518 "zlarfx.f"
	    i__2 = j + (c_dim1 << 1);
#line 518 "zlarfx.f"
	    i__3 = j + (c_dim1 << 1);
#line 518 "zlarfx.f"
	    z__2.r = sum.r * t2.r - sum.i * t2.i, z__2.i = sum.r * t2.i + 
		    sum.i * t2.r;
#line 518 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 518 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 519 "zlarfx.f"
	    i__2 = j + c_dim1 * 3;
#line 519 "zlarfx.f"
	    i__3 = j + c_dim1 * 3;
#line 519 "zlarfx.f"
	    z__2.r = sum.r * t3.r - sum.i * t3.i, z__2.i = sum.r * t3.i + 
		    sum.i * t3.r;
#line 519 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 519 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 520 "zlarfx.f"
	    i__2 = j + (c_dim1 << 2);
#line 520 "zlarfx.f"
	    i__3 = j + (c_dim1 << 2);
#line 520 "zlarfx.f"
	    z__2.r = sum.r * t4.r - sum.i * t4.i, z__2.i = sum.r * t4.i + 
		    sum.i * t4.r;
#line 520 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 520 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 521 "zlarfx.f"
	    i__2 = j + c_dim1 * 5;
#line 521 "zlarfx.f"
	    i__3 = j + c_dim1 * 5;
#line 521 "zlarfx.f"
	    z__2.r = sum.r * t5.r - sum.i * t5.i, z__2.i = sum.r * t5.i + 
		    sum.i * t5.r;
#line 521 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 521 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 522 "zlarfx.f"
/* L300: */
#line 522 "zlarfx.f"
	}
#line 523 "zlarfx.f"
	goto L410;
#line 524 "zlarfx.f"
L310:

/*        Special code for 6 x 6 Householder */

#line 528 "zlarfx.f"
	v1.r = v[1].r, v1.i = v[1].i;
#line 529 "zlarfx.f"
	d_cnjg(&z__2, &v1);
#line 529 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 529 "zlarfx.f"
	t1.r = z__1.r, t1.i = z__1.i;
#line 530 "zlarfx.f"
	v2.r = v[2].r, v2.i = v[2].i;
#line 531 "zlarfx.f"
	d_cnjg(&z__2, &v2);
#line 531 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 531 "zlarfx.f"
	t2.r = z__1.r, t2.i = z__1.i;
#line 532 "zlarfx.f"
	v3.r = v[3].r, v3.i = v[3].i;
#line 533 "zlarfx.f"
	d_cnjg(&z__2, &v3);
#line 533 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 533 "zlarfx.f"
	t3.r = z__1.r, t3.i = z__1.i;
#line 534 "zlarfx.f"
	v4.r = v[4].r, v4.i = v[4].i;
#line 535 "zlarfx.f"
	d_cnjg(&z__2, &v4);
#line 535 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 535 "zlarfx.f"
	t4.r = z__1.r, t4.i = z__1.i;
#line 536 "zlarfx.f"
	v5.r = v[5].r, v5.i = v[5].i;
#line 537 "zlarfx.f"
	d_cnjg(&z__2, &v5);
#line 537 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 537 "zlarfx.f"
	t5.r = z__1.r, t5.i = z__1.i;
#line 538 "zlarfx.f"
	v6.r = v[6].r, v6.i = v[6].i;
#line 539 "zlarfx.f"
	d_cnjg(&z__2, &v6);
#line 539 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 539 "zlarfx.f"
	t6.r = z__1.r, t6.i = z__1.i;
#line 540 "zlarfx.f"
	i__1 = *m;
#line 540 "zlarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 541 "zlarfx.f"
	    i__2 = j + c_dim1;
#line 541 "zlarfx.f"
	    z__6.r = v1.r * c__[i__2].r - v1.i * c__[i__2].i, z__6.i = v1.r * 
		    c__[i__2].i + v1.i * c__[i__2].r;
#line 541 "zlarfx.f"
	    i__3 = j + (c_dim1 << 1);
#line 541 "zlarfx.f"
	    z__7.r = v2.r * c__[i__3].r - v2.i * c__[i__3].i, z__7.i = v2.r * 
		    c__[i__3].i + v2.i * c__[i__3].r;
#line 541 "zlarfx.f"
	    z__5.r = z__6.r + z__7.r, z__5.i = z__6.i + z__7.i;
#line 541 "zlarfx.f"
	    i__4 = j + c_dim1 * 3;
#line 541 "zlarfx.f"
	    z__8.r = v3.r * c__[i__4].r - v3.i * c__[i__4].i, z__8.i = v3.r * 
		    c__[i__4].i + v3.i * c__[i__4].r;
#line 541 "zlarfx.f"
	    z__4.r = z__5.r + z__8.r, z__4.i = z__5.i + z__8.i;
#line 541 "zlarfx.f"
	    i__5 = j + (c_dim1 << 2);
#line 541 "zlarfx.f"
	    z__9.r = v4.r * c__[i__5].r - v4.i * c__[i__5].i, z__9.i = v4.r * 
		    c__[i__5].i + v4.i * c__[i__5].r;
#line 541 "zlarfx.f"
	    z__3.r = z__4.r + z__9.r, z__3.i = z__4.i + z__9.i;
#line 541 "zlarfx.f"
	    i__6 = j + c_dim1 * 5;
#line 541 "zlarfx.f"
	    z__10.r = v5.r * c__[i__6].r - v5.i * c__[i__6].i, z__10.i = v5.r 
		    * c__[i__6].i + v5.i * c__[i__6].r;
#line 541 "zlarfx.f"
	    z__2.r = z__3.r + z__10.r, z__2.i = z__3.i + z__10.i;
#line 541 "zlarfx.f"
	    i__7 = j + c_dim1 * 6;
#line 541 "zlarfx.f"
	    z__11.r = v6.r * c__[i__7].r - v6.i * c__[i__7].i, z__11.i = v6.r 
		    * c__[i__7].i + v6.i * c__[i__7].r;
#line 541 "zlarfx.f"
	    z__1.r = z__2.r + z__11.r, z__1.i = z__2.i + z__11.i;
#line 541 "zlarfx.f"
	    sum.r = z__1.r, sum.i = z__1.i;
#line 543 "zlarfx.f"
	    i__2 = j + c_dim1;
#line 543 "zlarfx.f"
	    i__3 = j + c_dim1;
#line 543 "zlarfx.f"
	    z__2.r = sum.r * t1.r - sum.i * t1.i, z__2.i = sum.r * t1.i + 
		    sum.i * t1.r;
#line 543 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 543 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 544 "zlarfx.f"
	    i__2 = j + (c_dim1 << 1);
#line 544 "zlarfx.f"
	    i__3 = j + (c_dim1 << 1);
#line 544 "zlarfx.f"
	    z__2.r = sum.r * t2.r - sum.i * t2.i, z__2.i = sum.r * t2.i + 
		    sum.i * t2.r;
#line 544 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 544 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 545 "zlarfx.f"
	    i__2 = j + c_dim1 * 3;
#line 545 "zlarfx.f"
	    i__3 = j + c_dim1 * 3;
#line 545 "zlarfx.f"
	    z__2.r = sum.r * t3.r - sum.i * t3.i, z__2.i = sum.r * t3.i + 
		    sum.i * t3.r;
#line 545 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 545 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 546 "zlarfx.f"
	    i__2 = j + (c_dim1 << 2);
#line 546 "zlarfx.f"
	    i__3 = j + (c_dim1 << 2);
#line 546 "zlarfx.f"
	    z__2.r = sum.r * t4.r - sum.i * t4.i, z__2.i = sum.r * t4.i + 
		    sum.i * t4.r;
#line 546 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 546 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 547 "zlarfx.f"
	    i__2 = j + c_dim1 * 5;
#line 547 "zlarfx.f"
	    i__3 = j + c_dim1 * 5;
#line 547 "zlarfx.f"
	    z__2.r = sum.r * t5.r - sum.i * t5.i, z__2.i = sum.r * t5.i + 
		    sum.i * t5.r;
#line 547 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 547 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 548 "zlarfx.f"
	    i__2 = j + c_dim1 * 6;
#line 548 "zlarfx.f"
	    i__3 = j + c_dim1 * 6;
#line 548 "zlarfx.f"
	    z__2.r = sum.r * t6.r - sum.i * t6.i, z__2.i = sum.r * t6.i + 
		    sum.i * t6.r;
#line 548 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 548 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 549 "zlarfx.f"
/* L320: */
#line 549 "zlarfx.f"
	}
#line 550 "zlarfx.f"
	goto L410;
#line 551 "zlarfx.f"
L330:

/*        Special code for 7 x 7 Householder */

#line 555 "zlarfx.f"
	v1.r = v[1].r, v1.i = v[1].i;
#line 556 "zlarfx.f"
	d_cnjg(&z__2, &v1);
#line 556 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 556 "zlarfx.f"
	t1.r = z__1.r, t1.i = z__1.i;
#line 557 "zlarfx.f"
	v2.r = v[2].r, v2.i = v[2].i;
#line 558 "zlarfx.f"
	d_cnjg(&z__2, &v2);
#line 558 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 558 "zlarfx.f"
	t2.r = z__1.r, t2.i = z__1.i;
#line 559 "zlarfx.f"
	v3.r = v[3].r, v3.i = v[3].i;
#line 560 "zlarfx.f"
	d_cnjg(&z__2, &v3);
#line 560 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 560 "zlarfx.f"
	t3.r = z__1.r, t3.i = z__1.i;
#line 561 "zlarfx.f"
	v4.r = v[4].r, v4.i = v[4].i;
#line 562 "zlarfx.f"
	d_cnjg(&z__2, &v4);
#line 562 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 562 "zlarfx.f"
	t4.r = z__1.r, t4.i = z__1.i;
#line 563 "zlarfx.f"
	v5.r = v[5].r, v5.i = v[5].i;
#line 564 "zlarfx.f"
	d_cnjg(&z__2, &v5);
#line 564 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 564 "zlarfx.f"
	t5.r = z__1.r, t5.i = z__1.i;
#line 565 "zlarfx.f"
	v6.r = v[6].r, v6.i = v[6].i;
#line 566 "zlarfx.f"
	d_cnjg(&z__2, &v6);
#line 566 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 566 "zlarfx.f"
	t6.r = z__1.r, t6.i = z__1.i;
#line 567 "zlarfx.f"
	v7.r = v[7].r, v7.i = v[7].i;
#line 568 "zlarfx.f"
	d_cnjg(&z__2, &v7);
#line 568 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 568 "zlarfx.f"
	t7.r = z__1.r, t7.i = z__1.i;
#line 569 "zlarfx.f"
	i__1 = *m;
#line 569 "zlarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 570 "zlarfx.f"
	    i__2 = j + c_dim1;
#line 570 "zlarfx.f"
	    z__7.r = v1.r * c__[i__2].r - v1.i * c__[i__2].i, z__7.i = v1.r * 
		    c__[i__2].i + v1.i * c__[i__2].r;
#line 570 "zlarfx.f"
	    i__3 = j + (c_dim1 << 1);
#line 570 "zlarfx.f"
	    z__8.r = v2.r * c__[i__3].r - v2.i * c__[i__3].i, z__8.i = v2.r * 
		    c__[i__3].i + v2.i * c__[i__3].r;
#line 570 "zlarfx.f"
	    z__6.r = z__7.r + z__8.r, z__6.i = z__7.i + z__8.i;
#line 570 "zlarfx.f"
	    i__4 = j + c_dim1 * 3;
#line 570 "zlarfx.f"
	    z__9.r = v3.r * c__[i__4].r - v3.i * c__[i__4].i, z__9.i = v3.r * 
		    c__[i__4].i + v3.i * c__[i__4].r;
#line 570 "zlarfx.f"
	    z__5.r = z__6.r + z__9.r, z__5.i = z__6.i + z__9.i;
#line 570 "zlarfx.f"
	    i__5 = j + (c_dim1 << 2);
#line 570 "zlarfx.f"
	    z__10.r = v4.r * c__[i__5].r - v4.i * c__[i__5].i, z__10.i = v4.r 
		    * c__[i__5].i + v4.i * c__[i__5].r;
#line 570 "zlarfx.f"
	    z__4.r = z__5.r + z__10.r, z__4.i = z__5.i + z__10.i;
#line 570 "zlarfx.f"
	    i__6 = j + c_dim1 * 5;
#line 570 "zlarfx.f"
	    z__11.r = v5.r * c__[i__6].r - v5.i * c__[i__6].i, z__11.i = v5.r 
		    * c__[i__6].i + v5.i * c__[i__6].r;
#line 570 "zlarfx.f"
	    z__3.r = z__4.r + z__11.r, z__3.i = z__4.i + z__11.i;
#line 570 "zlarfx.f"
	    i__7 = j + c_dim1 * 6;
#line 570 "zlarfx.f"
	    z__12.r = v6.r * c__[i__7].r - v6.i * c__[i__7].i, z__12.i = v6.r 
		    * c__[i__7].i + v6.i * c__[i__7].r;
#line 570 "zlarfx.f"
	    z__2.r = z__3.r + z__12.r, z__2.i = z__3.i + z__12.i;
#line 570 "zlarfx.f"
	    i__8 = j + c_dim1 * 7;
#line 570 "zlarfx.f"
	    z__13.r = v7.r * c__[i__8].r - v7.i * c__[i__8].i, z__13.i = v7.r 
		    * c__[i__8].i + v7.i * c__[i__8].r;
#line 570 "zlarfx.f"
	    z__1.r = z__2.r + z__13.r, z__1.i = z__2.i + z__13.i;
#line 570 "zlarfx.f"
	    sum.r = z__1.r, sum.i = z__1.i;
#line 573 "zlarfx.f"
	    i__2 = j + c_dim1;
#line 573 "zlarfx.f"
	    i__3 = j + c_dim1;
#line 573 "zlarfx.f"
	    z__2.r = sum.r * t1.r - sum.i * t1.i, z__2.i = sum.r * t1.i + 
		    sum.i * t1.r;
#line 573 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 573 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 574 "zlarfx.f"
	    i__2 = j + (c_dim1 << 1);
#line 574 "zlarfx.f"
	    i__3 = j + (c_dim1 << 1);
#line 574 "zlarfx.f"
	    z__2.r = sum.r * t2.r - sum.i * t2.i, z__2.i = sum.r * t2.i + 
		    sum.i * t2.r;
#line 574 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 574 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 575 "zlarfx.f"
	    i__2 = j + c_dim1 * 3;
#line 575 "zlarfx.f"
	    i__3 = j + c_dim1 * 3;
#line 575 "zlarfx.f"
	    z__2.r = sum.r * t3.r - sum.i * t3.i, z__2.i = sum.r * t3.i + 
		    sum.i * t3.r;
#line 575 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 575 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 576 "zlarfx.f"
	    i__2 = j + (c_dim1 << 2);
#line 576 "zlarfx.f"
	    i__3 = j + (c_dim1 << 2);
#line 576 "zlarfx.f"
	    z__2.r = sum.r * t4.r - sum.i * t4.i, z__2.i = sum.r * t4.i + 
		    sum.i * t4.r;
#line 576 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 576 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 577 "zlarfx.f"
	    i__2 = j + c_dim1 * 5;
#line 577 "zlarfx.f"
	    i__3 = j + c_dim1 * 5;
#line 577 "zlarfx.f"
	    z__2.r = sum.r * t5.r - sum.i * t5.i, z__2.i = sum.r * t5.i + 
		    sum.i * t5.r;
#line 577 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 577 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 578 "zlarfx.f"
	    i__2 = j + c_dim1 * 6;
#line 578 "zlarfx.f"
	    i__3 = j + c_dim1 * 6;
#line 578 "zlarfx.f"
	    z__2.r = sum.r * t6.r - sum.i * t6.i, z__2.i = sum.r * t6.i + 
		    sum.i * t6.r;
#line 578 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 578 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 579 "zlarfx.f"
	    i__2 = j + c_dim1 * 7;
#line 579 "zlarfx.f"
	    i__3 = j + c_dim1 * 7;
#line 579 "zlarfx.f"
	    z__2.r = sum.r * t7.r - sum.i * t7.i, z__2.i = sum.r * t7.i + 
		    sum.i * t7.r;
#line 579 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 579 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 580 "zlarfx.f"
/* L340: */
#line 580 "zlarfx.f"
	}
#line 581 "zlarfx.f"
	goto L410;
#line 582 "zlarfx.f"
L350:

/*        Special code for 8 x 8 Householder */

#line 586 "zlarfx.f"
	v1.r = v[1].r, v1.i = v[1].i;
#line 587 "zlarfx.f"
	d_cnjg(&z__2, &v1);
#line 587 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 587 "zlarfx.f"
	t1.r = z__1.r, t1.i = z__1.i;
#line 588 "zlarfx.f"
	v2.r = v[2].r, v2.i = v[2].i;
#line 589 "zlarfx.f"
	d_cnjg(&z__2, &v2);
#line 589 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 589 "zlarfx.f"
	t2.r = z__1.r, t2.i = z__1.i;
#line 590 "zlarfx.f"
	v3.r = v[3].r, v3.i = v[3].i;
#line 591 "zlarfx.f"
	d_cnjg(&z__2, &v3);
#line 591 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 591 "zlarfx.f"
	t3.r = z__1.r, t3.i = z__1.i;
#line 592 "zlarfx.f"
	v4.r = v[4].r, v4.i = v[4].i;
#line 593 "zlarfx.f"
	d_cnjg(&z__2, &v4);
#line 593 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 593 "zlarfx.f"
	t4.r = z__1.r, t4.i = z__1.i;
#line 594 "zlarfx.f"
	v5.r = v[5].r, v5.i = v[5].i;
#line 595 "zlarfx.f"
	d_cnjg(&z__2, &v5);
#line 595 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 595 "zlarfx.f"
	t5.r = z__1.r, t5.i = z__1.i;
#line 596 "zlarfx.f"
	v6.r = v[6].r, v6.i = v[6].i;
#line 597 "zlarfx.f"
	d_cnjg(&z__2, &v6);
#line 597 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 597 "zlarfx.f"
	t6.r = z__1.r, t6.i = z__1.i;
#line 598 "zlarfx.f"
	v7.r = v[7].r, v7.i = v[7].i;
#line 599 "zlarfx.f"
	d_cnjg(&z__2, &v7);
#line 599 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 599 "zlarfx.f"
	t7.r = z__1.r, t7.i = z__1.i;
#line 600 "zlarfx.f"
	v8.r = v[8].r, v8.i = v[8].i;
#line 601 "zlarfx.f"
	d_cnjg(&z__2, &v8);
#line 601 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 601 "zlarfx.f"
	t8.r = z__1.r, t8.i = z__1.i;
#line 602 "zlarfx.f"
	i__1 = *m;
#line 602 "zlarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 603 "zlarfx.f"
	    i__2 = j + c_dim1;
#line 603 "zlarfx.f"
	    z__8.r = v1.r * c__[i__2].r - v1.i * c__[i__2].i, z__8.i = v1.r * 
		    c__[i__2].i + v1.i * c__[i__2].r;
#line 603 "zlarfx.f"
	    i__3 = j + (c_dim1 << 1);
#line 603 "zlarfx.f"
	    z__9.r = v2.r * c__[i__3].r - v2.i * c__[i__3].i, z__9.i = v2.r * 
		    c__[i__3].i + v2.i * c__[i__3].r;
#line 603 "zlarfx.f"
	    z__7.r = z__8.r + z__9.r, z__7.i = z__8.i + z__9.i;
#line 603 "zlarfx.f"
	    i__4 = j + c_dim1 * 3;
#line 603 "zlarfx.f"
	    z__10.r = v3.r * c__[i__4].r - v3.i * c__[i__4].i, z__10.i = v3.r 
		    * c__[i__4].i + v3.i * c__[i__4].r;
#line 603 "zlarfx.f"
	    z__6.r = z__7.r + z__10.r, z__6.i = z__7.i + z__10.i;
#line 603 "zlarfx.f"
	    i__5 = j + (c_dim1 << 2);
#line 603 "zlarfx.f"
	    z__11.r = v4.r * c__[i__5].r - v4.i * c__[i__5].i, z__11.i = v4.r 
		    * c__[i__5].i + v4.i * c__[i__5].r;
#line 603 "zlarfx.f"
	    z__5.r = z__6.r + z__11.r, z__5.i = z__6.i + z__11.i;
#line 603 "zlarfx.f"
	    i__6 = j + c_dim1 * 5;
#line 603 "zlarfx.f"
	    z__12.r = v5.r * c__[i__6].r - v5.i * c__[i__6].i, z__12.i = v5.r 
		    * c__[i__6].i + v5.i * c__[i__6].r;
#line 603 "zlarfx.f"
	    z__4.r = z__5.r + z__12.r, z__4.i = z__5.i + z__12.i;
#line 603 "zlarfx.f"
	    i__7 = j + c_dim1 * 6;
#line 603 "zlarfx.f"
	    z__13.r = v6.r * c__[i__7].r - v6.i * c__[i__7].i, z__13.i = v6.r 
		    * c__[i__7].i + v6.i * c__[i__7].r;
#line 603 "zlarfx.f"
	    z__3.r = z__4.r + z__13.r, z__3.i = z__4.i + z__13.i;
#line 603 "zlarfx.f"
	    i__8 = j + c_dim1 * 7;
#line 603 "zlarfx.f"
	    z__14.r = v7.r * c__[i__8].r - v7.i * c__[i__8].i, z__14.i = v7.r 
		    * c__[i__8].i + v7.i * c__[i__8].r;
#line 603 "zlarfx.f"
	    z__2.r = z__3.r + z__14.r, z__2.i = z__3.i + z__14.i;
#line 603 "zlarfx.f"
	    i__9 = j + (c_dim1 << 3);
#line 603 "zlarfx.f"
	    z__15.r = v8.r * c__[i__9].r - v8.i * c__[i__9].i, z__15.i = v8.r 
		    * c__[i__9].i + v8.i * c__[i__9].r;
#line 603 "zlarfx.f"
	    z__1.r = z__2.r + z__15.r, z__1.i = z__2.i + z__15.i;
#line 603 "zlarfx.f"
	    sum.r = z__1.r, sum.i = z__1.i;
#line 606 "zlarfx.f"
	    i__2 = j + c_dim1;
#line 606 "zlarfx.f"
	    i__3 = j + c_dim1;
#line 606 "zlarfx.f"
	    z__2.r = sum.r * t1.r - sum.i * t1.i, z__2.i = sum.r * t1.i + 
		    sum.i * t1.r;
#line 606 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 606 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 607 "zlarfx.f"
	    i__2 = j + (c_dim1 << 1);
#line 607 "zlarfx.f"
	    i__3 = j + (c_dim1 << 1);
#line 607 "zlarfx.f"
	    z__2.r = sum.r * t2.r - sum.i * t2.i, z__2.i = sum.r * t2.i + 
		    sum.i * t2.r;
#line 607 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 607 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 608 "zlarfx.f"
	    i__2 = j + c_dim1 * 3;
#line 608 "zlarfx.f"
	    i__3 = j + c_dim1 * 3;
#line 608 "zlarfx.f"
	    z__2.r = sum.r * t3.r - sum.i * t3.i, z__2.i = sum.r * t3.i + 
		    sum.i * t3.r;
#line 608 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 608 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 609 "zlarfx.f"
	    i__2 = j + (c_dim1 << 2);
#line 609 "zlarfx.f"
	    i__3 = j + (c_dim1 << 2);
#line 609 "zlarfx.f"
	    z__2.r = sum.r * t4.r - sum.i * t4.i, z__2.i = sum.r * t4.i + 
		    sum.i * t4.r;
#line 609 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 609 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 610 "zlarfx.f"
	    i__2 = j + c_dim1 * 5;
#line 610 "zlarfx.f"
	    i__3 = j + c_dim1 * 5;
#line 610 "zlarfx.f"
	    z__2.r = sum.r * t5.r - sum.i * t5.i, z__2.i = sum.r * t5.i + 
		    sum.i * t5.r;
#line 610 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 610 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 611 "zlarfx.f"
	    i__2 = j + c_dim1 * 6;
#line 611 "zlarfx.f"
	    i__3 = j + c_dim1 * 6;
#line 611 "zlarfx.f"
	    z__2.r = sum.r * t6.r - sum.i * t6.i, z__2.i = sum.r * t6.i + 
		    sum.i * t6.r;
#line 611 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 611 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 612 "zlarfx.f"
	    i__2 = j + c_dim1 * 7;
#line 612 "zlarfx.f"
	    i__3 = j + c_dim1 * 7;
#line 612 "zlarfx.f"
	    z__2.r = sum.r * t7.r - sum.i * t7.i, z__2.i = sum.r * t7.i + 
		    sum.i * t7.r;
#line 612 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 612 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 613 "zlarfx.f"
	    i__2 = j + (c_dim1 << 3);
#line 613 "zlarfx.f"
	    i__3 = j + (c_dim1 << 3);
#line 613 "zlarfx.f"
	    z__2.r = sum.r * t8.r - sum.i * t8.i, z__2.i = sum.r * t8.i + 
		    sum.i * t8.r;
#line 613 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 613 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 614 "zlarfx.f"
/* L360: */
#line 614 "zlarfx.f"
	}
#line 615 "zlarfx.f"
	goto L410;
#line 616 "zlarfx.f"
L370:

/*        Special code for 9 x 9 Householder */

#line 620 "zlarfx.f"
	v1.r = v[1].r, v1.i = v[1].i;
#line 621 "zlarfx.f"
	d_cnjg(&z__2, &v1);
#line 621 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 621 "zlarfx.f"
	t1.r = z__1.r, t1.i = z__1.i;
#line 622 "zlarfx.f"
	v2.r = v[2].r, v2.i = v[2].i;
#line 623 "zlarfx.f"
	d_cnjg(&z__2, &v2);
#line 623 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 623 "zlarfx.f"
	t2.r = z__1.r, t2.i = z__1.i;
#line 624 "zlarfx.f"
	v3.r = v[3].r, v3.i = v[3].i;
#line 625 "zlarfx.f"
	d_cnjg(&z__2, &v3);
#line 625 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 625 "zlarfx.f"
	t3.r = z__1.r, t3.i = z__1.i;
#line 626 "zlarfx.f"
	v4.r = v[4].r, v4.i = v[4].i;
#line 627 "zlarfx.f"
	d_cnjg(&z__2, &v4);
#line 627 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 627 "zlarfx.f"
	t4.r = z__1.r, t4.i = z__1.i;
#line 628 "zlarfx.f"
	v5.r = v[5].r, v5.i = v[5].i;
#line 629 "zlarfx.f"
	d_cnjg(&z__2, &v5);
#line 629 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 629 "zlarfx.f"
	t5.r = z__1.r, t5.i = z__1.i;
#line 630 "zlarfx.f"
	v6.r = v[6].r, v6.i = v[6].i;
#line 631 "zlarfx.f"
	d_cnjg(&z__2, &v6);
#line 631 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 631 "zlarfx.f"
	t6.r = z__1.r, t6.i = z__1.i;
#line 632 "zlarfx.f"
	v7.r = v[7].r, v7.i = v[7].i;
#line 633 "zlarfx.f"
	d_cnjg(&z__2, &v7);
#line 633 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 633 "zlarfx.f"
	t7.r = z__1.r, t7.i = z__1.i;
#line 634 "zlarfx.f"
	v8.r = v[8].r, v8.i = v[8].i;
#line 635 "zlarfx.f"
	d_cnjg(&z__2, &v8);
#line 635 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 635 "zlarfx.f"
	t8.r = z__1.r, t8.i = z__1.i;
#line 636 "zlarfx.f"
	v9.r = v[9].r, v9.i = v[9].i;
#line 637 "zlarfx.f"
	d_cnjg(&z__2, &v9);
#line 637 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 637 "zlarfx.f"
	t9.r = z__1.r, t9.i = z__1.i;
#line 638 "zlarfx.f"
	i__1 = *m;
#line 638 "zlarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 639 "zlarfx.f"
	    i__2 = j + c_dim1;
#line 639 "zlarfx.f"
	    z__9.r = v1.r * c__[i__2].r - v1.i * c__[i__2].i, z__9.i = v1.r * 
		    c__[i__2].i + v1.i * c__[i__2].r;
#line 639 "zlarfx.f"
	    i__3 = j + (c_dim1 << 1);
#line 639 "zlarfx.f"
	    z__10.r = v2.r * c__[i__3].r - v2.i * c__[i__3].i, z__10.i = v2.r 
		    * c__[i__3].i + v2.i * c__[i__3].r;
#line 639 "zlarfx.f"
	    z__8.r = z__9.r + z__10.r, z__8.i = z__9.i + z__10.i;
#line 639 "zlarfx.f"
	    i__4 = j + c_dim1 * 3;
#line 639 "zlarfx.f"
	    z__11.r = v3.r * c__[i__4].r - v3.i * c__[i__4].i, z__11.i = v3.r 
		    * c__[i__4].i + v3.i * c__[i__4].r;
#line 639 "zlarfx.f"
	    z__7.r = z__8.r + z__11.r, z__7.i = z__8.i + z__11.i;
#line 639 "zlarfx.f"
	    i__5 = j + (c_dim1 << 2);
#line 639 "zlarfx.f"
	    z__12.r = v4.r * c__[i__5].r - v4.i * c__[i__5].i, z__12.i = v4.r 
		    * c__[i__5].i + v4.i * c__[i__5].r;
#line 639 "zlarfx.f"
	    z__6.r = z__7.r + z__12.r, z__6.i = z__7.i + z__12.i;
#line 639 "zlarfx.f"
	    i__6 = j + c_dim1 * 5;
#line 639 "zlarfx.f"
	    z__13.r = v5.r * c__[i__6].r - v5.i * c__[i__6].i, z__13.i = v5.r 
		    * c__[i__6].i + v5.i * c__[i__6].r;
#line 639 "zlarfx.f"
	    z__5.r = z__6.r + z__13.r, z__5.i = z__6.i + z__13.i;
#line 639 "zlarfx.f"
	    i__7 = j + c_dim1 * 6;
#line 639 "zlarfx.f"
	    z__14.r = v6.r * c__[i__7].r - v6.i * c__[i__7].i, z__14.i = v6.r 
		    * c__[i__7].i + v6.i * c__[i__7].r;
#line 639 "zlarfx.f"
	    z__4.r = z__5.r + z__14.r, z__4.i = z__5.i + z__14.i;
#line 639 "zlarfx.f"
	    i__8 = j + c_dim1 * 7;
#line 639 "zlarfx.f"
	    z__15.r = v7.r * c__[i__8].r - v7.i * c__[i__8].i, z__15.i = v7.r 
		    * c__[i__8].i + v7.i * c__[i__8].r;
#line 639 "zlarfx.f"
	    z__3.r = z__4.r + z__15.r, z__3.i = z__4.i + z__15.i;
#line 639 "zlarfx.f"
	    i__9 = j + (c_dim1 << 3);
#line 639 "zlarfx.f"
	    z__16.r = v8.r * c__[i__9].r - v8.i * c__[i__9].i, z__16.i = v8.r 
		    * c__[i__9].i + v8.i * c__[i__9].r;
#line 639 "zlarfx.f"
	    z__2.r = z__3.r + z__16.r, z__2.i = z__3.i + z__16.i;
#line 639 "zlarfx.f"
	    i__10 = j + c_dim1 * 9;
#line 639 "zlarfx.f"
	    z__17.r = v9.r * c__[i__10].r - v9.i * c__[i__10].i, z__17.i = 
		    v9.r * c__[i__10].i + v9.i * c__[i__10].r;
#line 639 "zlarfx.f"
	    z__1.r = z__2.r + z__17.r, z__1.i = z__2.i + z__17.i;
#line 639 "zlarfx.f"
	    sum.r = z__1.r, sum.i = z__1.i;
#line 642 "zlarfx.f"
	    i__2 = j + c_dim1;
#line 642 "zlarfx.f"
	    i__3 = j + c_dim1;
#line 642 "zlarfx.f"
	    z__2.r = sum.r * t1.r - sum.i * t1.i, z__2.i = sum.r * t1.i + 
		    sum.i * t1.r;
#line 642 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 642 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 643 "zlarfx.f"
	    i__2 = j + (c_dim1 << 1);
#line 643 "zlarfx.f"
	    i__3 = j + (c_dim1 << 1);
#line 643 "zlarfx.f"
	    z__2.r = sum.r * t2.r - sum.i * t2.i, z__2.i = sum.r * t2.i + 
		    sum.i * t2.r;
#line 643 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 643 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 644 "zlarfx.f"
	    i__2 = j + c_dim1 * 3;
#line 644 "zlarfx.f"
	    i__3 = j + c_dim1 * 3;
#line 644 "zlarfx.f"
	    z__2.r = sum.r * t3.r - sum.i * t3.i, z__2.i = sum.r * t3.i + 
		    sum.i * t3.r;
#line 644 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 644 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 645 "zlarfx.f"
	    i__2 = j + (c_dim1 << 2);
#line 645 "zlarfx.f"
	    i__3 = j + (c_dim1 << 2);
#line 645 "zlarfx.f"
	    z__2.r = sum.r * t4.r - sum.i * t4.i, z__2.i = sum.r * t4.i + 
		    sum.i * t4.r;
#line 645 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 645 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 646 "zlarfx.f"
	    i__2 = j + c_dim1 * 5;
#line 646 "zlarfx.f"
	    i__3 = j + c_dim1 * 5;
#line 646 "zlarfx.f"
	    z__2.r = sum.r * t5.r - sum.i * t5.i, z__2.i = sum.r * t5.i + 
		    sum.i * t5.r;
#line 646 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 646 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 647 "zlarfx.f"
	    i__2 = j + c_dim1 * 6;
#line 647 "zlarfx.f"
	    i__3 = j + c_dim1 * 6;
#line 647 "zlarfx.f"
	    z__2.r = sum.r * t6.r - sum.i * t6.i, z__2.i = sum.r * t6.i + 
		    sum.i * t6.r;
#line 647 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 647 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 648 "zlarfx.f"
	    i__2 = j + c_dim1 * 7;
#line 648 "zlarfx.f"
	    i__3 = j + c_dim1 * 7;
#line 648 "zlarfx.f"
	    z__2.r = sum.r * t7.r - sum.i * t7.i, z__2.i = sum.r * t7.i + 
		    sum.i * t7.r;
#line 648 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 648 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 649 "zlarfx.f"
	    i__2 = j + (c_dim1 << 3);
#line 649 "zlarfx.f"
	    i__3 = j + (c_dim1 << 3);
#line 649 "zlarfx.f"
	    z__2.r = sum.r * t8.r - sum.i * t8.i, z__2.i = sum.r * t8.i + 
		    sum.i * t8.r;
#line 649 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 649 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 650 "zlarfx.f"
	    i__2 = j + c_dim1 * 9;
#line 650 "zlarfx.f"
	    i__3 = j + c_dim1 * 9;
#line 650 "zlarfx.f"
	    z__2.r = sum.r * t9.r - sum.i * t9.i, z__2.i = sum.r * t9.i + 
		    sum.i * t9.r;
#line 650 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 650 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 651 "zlarfx.f"
/* L380: */
#line 651 "zlarfx.f"
	}
#line 652 "zlarfx.f"
	goto L410;
#line 653 "zlarfx.f"
L390:

/*        Special code for 10 x 10 Householder */

#line 657 "zlarfx.f"
	v1.r = v[1].r, v1.i = v[1].i;
#line 658 "zlarfx.f"
	d_cnjg(&z__2, &v1);
#line 658 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 658 "zlarfx.f"
	t1.r = z__1.r, t1.i = z__1.i;
#line 659 "zlarfx.f"
	v2.r = v[2].r, v2.i = v[2].i;
#line 660 "zlarfx.f"
	d_cnjg(&z__2, &v2);
#line 660 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 660 "zlarfx.f"
	t2.r = z__1.r, t2.i = z__1.i;
#line 661 "zlarfx.f"
	v3.r = v[3].r, v3.i = v[3].i;
#line 662 "zlarfx.f"
	d_cnjg(&z__2, &v3);
#line 662 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 662 "zlarfx.f"
	t3.r = z__1.r, t3.i = z__1.i;
#line 663 "zlarfx.f"
	v4.r = v[4].r, v4.i = v[4].i;
#line 664 "zlarfx.f"
	d_cnjg(&z__2, &v4);
#line 664 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 664 "zlarfx.f"
	t4.r = z__1.r, t4.i = z__1.i;
#line 665 "zlarfx.f"
	v5.r = v[5].r, v5.i = v[5].i;
#line 666 "zlarfx.f"
	d_cnjg(&z__2, &v5);
#line 666 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 666 "zlarfx.f"
	t5.r = z__1.r, t5.i = z__1.i;
#line 667 "zlarfx.f"
	v6.r = v[6].r, v6.i = v[6].i;
#line 668 "zlarfx.f"
	d_cnjg(&z__2, &v6);
#line 668 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 668 "zlarfx.f"
	t6.r = z__1.r, t6.i = z__1.i;
#line 669 "zlarfx.f"
	v7.r = v[7].r, v7.i = v[7].i;
#line 670 "zlarfx.f"
	d_cnjg(&z__2, &v7);
#line 670 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 670 "zlarfx.f"
	t7.r = z__1.r, t7.i = z__1.i;
#line 671 "zlarfx.f"
	v8.r = v[8].r, v8.i = v[8].i;
#line 672 "zlarfx.f"
	d_cnjg(&z__2, &v8);
#line 672 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 672 "zlarfx.f"
	t8.r = z__1.r, t8.i = z__1.i;
#line 673 "zlarfx.f"
	v9.r = v[9].r, v9.i = v[9].i;
#line 674 "zlarfx.f"
	d_cnjg(&z__2, &v9);
#line 674 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 674 "zlarfx.f"
	t9.r = z__1.r, t9.i = z__1.i;
#line 675 "zlarfx.f"
	v10.r = v[10].r, v10.i = v[10].i;
#line 676 "zlarfx.f"
	d_cnjg(&z__2, &v10);
#line 676 "zlarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 676 "zlarfx.f"
	t10.r = z__1.r, t10.i = z__1.i;
#line 677 "zlarfx.f"
	i__1 = *m;
#line 677 "zlarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 678 "zlarfx.f"
	    i__2 = j + c_dim1;
#line 678 "zlarfx.f"
	    z__10.r = v1.r * c__[i__2].r - v1.i * c__[i__2].i, z__10.i = v1.r 
		    * c__[i__2].i + v1.i * c__[i__2].r;
#line 678 "zlarfx.f"
	    i__3 = j + (c_dim1 << 1);
#line 678 "zlarfx.f"
	    z__11.r = v2.r * c__[i__3].r - v2.i * c__[i__3].i, z__11.i = v2.r 
		    * c__[i__3].i + v2.i * c__[i__3].r;
#line 678 "zlarfx.f"
	    z__9.r = z__10.r + z__11.r, z__9.i = z__10.i + z__11.i;
#line 678 "zlarfx.f"
	    i__4 = j + c_dim1 * 3;
#line 678 "zlarfx.f"
	    z__12.r = v3.r * c__[i__4].r - v3.i * c__[i__4].i, z__12.i = v3.r 
		    * c__[i__4].i + v3.i * c__[i__4].r;
#line 678 "zlarfx.f"
	    z__8.r = z__9.r + z__12.r, z__8.i = z__9.i + z__12.i;
#line 678 "zlarfx.f"
	    i__5 = j + (c_dim1 << 2);
#line 678 "zlarfx.f"
	    z__13.r = v4.r * c__[i__5].r - v4.i * c__[i__5].i, z__13.i = v4.r 
		    * c__[i__5].i + v4.i * c__[i__5].r;
#line 678 "zlarfx.f"
	    z__7.r = z__8.r + z__13.r, z__7.i = z__8.i + z__13.i;
#line 678 "zlarfx.f"
	    i__6 = j + c_dim1 * 5;
#line 678 "zlarfx.f"
	    z__14.r = v5.r * c__[i__6].r - v5.i * c__[i__6].i, z__14.i = v5.r 
		    * c__[i__6].i + v5.i * c__[i__6].r;
#line 678 "zlarfx.f"
	    z__6.r = z__7.r + z__14.r, z__6.i = z__7.i + z__14.i;
#line 678 "zlarfx.f"
	    i__7 = j + c_dim1 * 6;
#line 678 "zlarfx.f"
	    z__15.r = v6.r * c__[i__7].r - v6.i * c__[i__7].i, z__15.i = v6.r 
		    * c__[i__7].i + v6.i * c__[i__7].r;
#line 678 "zlarfx.f"
	    z__5.r = z__6.r + z__15.r, z__5.i = z__6.i + z__15.i;
#line 678 "zlarfx.f"
	    i__8 = j + c_dim1 * 7;
#line 678 "zlarfx.f"
	    z__16.r = v7.r * c__[i__8].r - v7.i * c__[i__8].i, z__16.i = v7.r 
		    * c__[i__8].i + v7.i * c__[i__8].r;
#line 678 "zlarfx.f"
	    z__4.r = z__5.r + z__16.r, z__4.i = z__5.i + z__16.i;
#line 678 "zlarfx.f"
	    i__9 = j + (c_dim1 << 3);
#line 678 "zlarfx.f"
	    z__17.r = v8.r * c__[i__9].r - v8.i * c__[i__9].i, z__17.i = v8.r 
		    * c__[i__9].i + v8.i * c__[i__9].r;
#line 678 "zlarfx.f"
	    z__3.r = z__4.r + z__17.r, z__3.i = z__4.i + z__17.i;
#line 678 "zlarfx.f"
	    i__10 = j + c_dim1 * 9;
#line 678 "zlarfx.f"
	    z__18.r = v9.r * c__[i__10].r - v9.i * c__[i__10].i, z__18.i = 
		    v9.r * c__[i__10].i + v9.i * c__[i__10].r;
#line 678 "zlarfx.f"
	    z__2.r = z__3.r + z__18.r, z__2.i = z__3.i + z__18.i;
#line 678 "zlarfx.f"
	    i__11 = j + c_dim1 * 10;
#line 678 "zlarfx.f"
	    z__19.r = v10.r * c__[i__11].r - v10.i * c__[i__11].i, z__19.i = 
		    v10.r * c__[i__11].i + v10.i * c__[i__11].r;
#line 678 "zlarfx.f"
	    z__1.r = z__2.r + z__19.r, z__1.i = z__2.i + z__19.i;
#line 678 "zlarfx.f"
	    sum.r = z__1.r, sum.i = z__1.i;
#line 682 "zlarfx.f"
	    i__2 = j + c_dim1;
#line 682 "zlarfx.f"
	    i__3 = j + c_dim1;
#line 682 "zlarfx.f"
	    z__2.r = sum.r * t1.r - sum.i * t1.i, z__2.i = sum.r * t1.i + 
		    sum.i * t1.r;
#line 682 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 682 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 683 "zlarfx.f"
	    i__2 = j + (c_dim1 << 1);
#line 683 "zlarfx.f"
	    i__3 = j + (c_dim1 << 1);
#line 683 "zlarfx.f"
	    z__2.r = sum.r * t2.r - sum.i * t2.i, z__2.i = sum.r * t2.i + 
		    sum.i * t2.r;
#line 683 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 683 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 684 "zlarfx.f"
	    i__2 = j + c_dim1 * 3;
#line 684 "zlarfx.f"
	    i__3 = j + c_dim1 * 3;
#line 684 "zlarfx.f"
	    z__2.r = sum.r * t3.r - sum.i * t3.i, z__2.i = sum.r * t3.i + 
		    sum.i * t3.r;
#line 684 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 684 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 685 "zlarfx.f"
	    i__2 = j + (c_dim1 << 2);
#line 685 "zlarfx.f"
	    i__3 = j + (c_dim1 << 2);
#line 685 "zlarfx.f"
	    z__2.r = sum.r * t4.r - sum.i * t4.i, z__2.i = sum.r * t4.i + 
		    sum.i * t4.r;
#line 685 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 685 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 686 "zlarfx.f"
	    i__2 = j + c_dim1 * 5;
#line 686 "zlarfx.f"
	    i__3 = j + c_dim1 * 5;
#line 686 "zlarfx.f"
	    z__2.r = sum.r * t5.r - sum.i * t5.i, z__2.i = sum.r * t5.i + 
		    sum.i * t5.r;
#line 686 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 686 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 687 "zlarfx.f"
	    i__2 = j + c_dim1 * 6;
#line 687 "zlarfx.f"
	    i__3 = j + c_dim1 * 6;
#line 687 "zlarfx.f"
	    z__2.r = sum.r * t6.r - sum.i * t6.i, z__2.i = sum.r * t6.i + 
		    sum.i * t6.r;
#line 687 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 687 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 688 "zlarfx.f"
	    i__2 = j + c_dim1 * 7;
#line 688 "zlarfx.f"
	    i__3 = j + c_dim1 * 7;
#line 688 "zlarfx.f"
	    z__2.r = sum.r * t7.r - sum.i * t7.i, z__2.i = sum.r * t7.i + 
		    sum.i * t7.r;
#line 688 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 688 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 689 "zlarfx.f"
	    i__2 = j + (c_dim1 << 3);
#line 689 "zlarfx.f"
	    i__3 = j + (c_dim1 << 3);
#line 689 "zlarfx.f"
	    z__2.r = sum.r * t8.r - sum.i * t8.i, z__2.i = sum.r * t8.i + 
		    sum.i * t8.r;
#line 689 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 689 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 690 "zlarfx.f"
	    i__2 = j + c_dim1 * 9;
#line 690 "zlarfx.f"
	    i__3 = j + c_dim1 * 9;
#line 690 "zlarfx.f"
	    z__2.r = sum.r * t9.r - sum.i * t9.i, z__2.i = sum.r * t9.i + 
		    sum.i * t9.r;
#line 690 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 690 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 691 "zlarfx.f"
	    i__2 = j + c_dim1 * 10;
#line 691 "zlarfx.f"
	    i__3 = j + c_dim1 * 10;
#line 691 "zlarfx.f"
	    z__2.r = sum.r * t10.r - sum.i * t10.i, z__2.i = sum.r * t10.i + 
		    sum.i * t10.r;
#line 691 "zlarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 691 "zlarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 692 "zlarfx.f"
/* L400: */
#line 692 "zlarfx.f"
	}
#line 693 "zlarfx.f"
	goto L410;
#line 694 "zlarfx.f"
    }
#line 695 "zlarfx.f"
L410:
#line 696 "zlarfx.f"
    return 0;

/*     End of ZLARFX */

} /* zlarfx_ */

