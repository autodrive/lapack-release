#line 1 "clarfx.f"
/* clarfx.f -- translated by f2c (version 20100827).
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

#line 1 "clarfx.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b CLARFX applies an elementary reflector to a general rectangular matrix, with loop unrolling whe
n the reflector has order ≤ 10. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLARFX + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clarfx.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clarfx.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clarfx.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CLARFX( SIDE, M, N, V, TAU, C, LDC, WORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          SIDE */
/*       INTEGER            LDC, M, N */
/*       COMPLEX            TAU */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX            C( LDC, * ), V( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLARFX applies a complex elementary reflector H to a complex m by n */
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
/* >          V is COMPLEX array, dimension (M) if SIDE = 'L' */
/* >                                        or (N) if SIDE = 'R' */
/* >          The vector v in the representation of H. */
/* > \endverbatim */
/* > */
/* > \param[in] TAU */
/* > \verbatim */
/* >          TAU is COMPLEX */
/* >          The value tau in the representation of H. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* >          C is COMPLEX array, dimension (LDC,N) */
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
/* >          WORK is COMPLEX array, dimension (N) if SIDE = 'L' */
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

/* > \ingroup complexOTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int clarfx_(char *side, integer *m, integer *n, 
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
    extern /* Subroutine */ int clarf_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, doublecomplex *, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);


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

#line 160 "clarfx.f"
    /* Parameter adjustments */
#line 160 "clarfx.f"
    --v;
#line 160 "clarfx.f"
    c_dim1 = *ldc;
#line 160 "clarfx.f"
    c_offset = 1 + c_dim1;
#line 160 "clarfx.f"
    c__ -= c_offset;
#line 160 "clarfx.f"
    --work;
#line 160 "clarfx.f"

#line 160 "clarfx.f"
    /* Function Body */
#line 160 "clarfx.f"
    if (tau->r == 0. && tau->i == 0.) {
#line 160 "clarfx.f"
	return 0;
#line 160 "clarfx.f"
    }
#line 162 "clarfx.f"
    if (lsame_(side, "L", (ftnlen)1, (ftnlen)1)) {

/*        Form  H * C, where H has order m. */

#line 166 "clarfx.f"
	switch (*m) {
#line 166 "clarfx.f"
	    case 1:  goto L10;
#line 166 "clarfx.f"
	    case 2:  goto L30;
#line 166 "clarfx.f"
	    case 3:  goto L50;
#line 166 "clarfx.f"
	    case 4:  goto L70;
#line 166 "clarfx.f"
	    case 5:  goto L90;
#line 166 "clarfx.f"
	    case 6:  goto L110;
#line 166 "clarfx.f"
	    case 7:  goto L130;
#line 166 "clarfx.f"
	    case 8:  goto L150;
#line 166 "clarfx.f"
	    case 9:  goto L170;
#line 166 "clarfx.f"
	    case 10:  goto L190;
#line 166 "clarfx.f"
	}

/*        Code for general M */

#line 171 "clarfx.f"
	clarf_(side, m, n, &v[1], &c__1, tau, &c__[c_offset], ldc, &work[1], (
		ftnlen)1);
#line 172 "clarfx.f"
	goto L410;
#line 173 "clarfx.f"
L10:

/*        Special code for 1 x 1 Householder */

#line 177 "clarfx.f"
	z__3.r = tau->r * v[1].r - tau->i * v[1].i, z__3.i = tau->r * v[1].i 
		+ tau->i * v[1].r;
#line 177 "clarfx.f"
	d_cnjg(&z__4, &v[1]);
#line 177 "clarfx.f"
	z__2.r = z__3.r * z__4.r - z__3.i * z__4.i, z__2.i = z__3.r * z__4.i 
		+ z__3.i * z__4.r;
#line 177 "clarfx.f"
	z__1.r = 1. - z__2.r, z__1.i = 0. - z__2.i;
#line 177 "clarfx.f"
	t1.r = z__1.r, t1.i = z__1.i;
#line 178 "clarfx.f"
	i__1 = *n;
#line 178 "clarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 179 "clarfx.f"
	    i__2 = j * c_dim1 + 1;
#line 179 "clarfx.f"
	    i__3 = j * c_dim1 + 1;
#line 179 "clarfx.f"
	    z__1.r = t1.r * c__[i__3].r - t1.i * c__[i__3].i, z__1.i = t1.r * 
		    c__[i__3].i + t1.i * c__[i__3].r;
#line 179 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 180 "clarfx.f"
/* L20: */
#line 180 "clarfx.f"
	}
#line 181 "clarfx.f"
	goto L410;
#line 182 "clarfx.f"
L30:

/*        Special code for 2 x 2 Householder */

#line 186 "clarfx.f"
	d_cnjg(&z__1, &v[1]);
#line 186 "clarfx.f"
	v1.r = z__1.r, v1.i = z__1.i;
#line 187 "clarfx.f"
	d_cnjg(&z__2, &v1);
#line 187 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 187 "clarfx.f"
	t1.r = z__1.r, t1.i = z__1.i;
#line 188 "clarfx.f"
	d_cnjg(&z__1, &v[2]);
#line 188 "clarfx.f"
	v2.r = z__1.r, v2.i = z__1.i;
#line 189 "clarfx.f"
	d_cnjg(&z__2, &v2);
#line 189 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 189 "clarfx.f"
	t2.r = z__1.r, t2.i = z__1.i;
#line 190 "clarfx.f"
	i__1 = *n;
#line 190 "clarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 191 "clarfx.f"
	    i__2 = j * c_dim1 + 1;
#line 191 "clarfx.f"
	    z__2.r = v1.r * c__[i__2].r - v1.i * c__[i__2].i, z__2.i = v1.r * 
		    c__[i__2].i + v1.i * c__[i__2].r;
#line 191 "clarfx.f"
	    i__3 = j * c_dim1 + 2;
#line 191 "clarfx.f"
	    z__3.r = v2.r * c__[i__3].r - v2.i * c__[i__3].i, z__3.i = v2.r * 
		    c__[i__3].i + v2.i * c__[i__3].r;
#line 191 "clarfx.f"
	    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 191 "clarfx.f"
	    sum.r = z__1.r, sum.i = z__1.i;
#line 192 "clarfx.f"
	    i__2 = j * c_dim1 + 1;
#line 192 "clarfx.f"
	    i__3 = j * c_dim1 + 1;
#line 192 "clarfx.f"
	    z__2.r = sum.r * t1.r - sum.i * t1.i, z__2.i = sum.r * t1.i + 
		    sum.i * t1.r;
#line 192 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 192 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 193 "clarfx.f"
	    i__2 = j * c_dim1 + 2;
#line 193 "clarfx.f"
	    i__3 = j * c_dim1 + 2;
#line 193 "clarfx.f"
	    z__2.r = sum.r * t2.r - sum.i * t2.i, z__2.i = sum.r * t2.i + 
		    sum.i * t2.r;
#line 193 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 193 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 194 "clarfx.f"
/* L40: */
#line 194 "clarfx.f"
	}
#line 195 "clarfx.f"
	goto L410;
#line 196 "clarfx.f"
L50:

/*        Special code for 3 x 3 Householder */

#line 200 "clarfx.f"
	d_cnjg(&z__1, &v[1]);
#line 200 "clarfx.f"
	v1.r = z__1.r, v1.i = z__1.i;
#line 201 "clarfx.f"
	d_cnjg(&z__2, &v1);
#line 201 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 201 "clarfx.f"
	t1.r = z__1.r, t1.i = z__1.i;
#line 202 "clarfx.f"
	d_cnjg(&z__1, &v[2]);
#line 202 "clarfx.f"
	v2.r = z__1.r, v2.i = z__1.i;
#line 203 "clarfx.f"
	d_cnjg(&z__2, &v2);
#line 203 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 203 "clarfx.f"
	t2.r = z__1.r, t2.i = z__1.i;
#line 204 "clarfx.f"
	d_cnjg(&z__1, &v[3]);
#line 204 "clarfx.f"
	v3.r = z__1.r, v3.i = z__1.i;
#line 205 "clarfx.f"
	d_cnjg(&z__2, &v3);
#line 205 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 205 "clarfx.f"
	t3.r = z__1.r, t3.i = z__1.i;
#line 206 "clarfx.f"
	i__1 = *n;
#line 206 "clarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 207 "clarfx.f"
	    i__2 = j * c_dim1 + 1;
#line 207 "clarfx.f"
	    z__3.r = v1.r * c__[i__2].r - v1.i * c__[i__2].i, z__3.i = v1.r * 
		    c__[i__2].i + v1.i * c__[i__2].r;
#line 207 "clarfx.f"
	    i__3 = j * c_dim1 + 2;
#line 207 "clarfx.f"
	    z__4.r = v2.r * c__[i__3].r - v2.i * c__[i__3].i, z__4.i = v2.r * 
		    c__[i__3].i + v2.i * c__[i__3].r;
#line 207 "clarfx.f"
	    z__2.r = z__3.r + z__4.r, z__2.i = z__3.i + z__4.i;
#line 207 "clarfx.f"
	    i__4 = j * c_dim1 + 3;
#line 207 "clarfx.f"
	    z__5.r = v3.r * c__[i__4].r - v3.i * c__[i__4].i, z__5.i = v3.r * 
		    c__[i__4].i + v3.i * c__[i__4].r;
#line 207 "clarfx.f"
	    z__1.r = z__2.r + z__5.r, z__1.i = z__2.i + z__5.i;
#line 207 "clarfx.f"
	    sum.r = z__1.r, sum.i = z__1.i;
#line 208 "clarfx.f"
	    i__2 = j * c_dim1 + 1;
#line 208 "clarfx.f"
	    i__3 = j * c_dim1 + 1;
#line 208 "clarfx.f"
	    z__2.r = sum.r * t1.r - sum.i * t1.i, z__2.i = sum.r * t1.i + 
		    sum.i * t1.r;
#line 208 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 208 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 209 "clarfx.f"
	    i__2 = j * c_dim1 + 2;
#line 209 "clarfx.f"
	    i__3 = j * c_dim1 + 2;
#line 209 "clarfx.f"
	    z__2.r = sum.r * t2.r - sum.i * t2.i, z__2.i = sum.r * t2.i + 
		    sum.i * t2.r;
#line 209 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 209 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 210 "clarfx.f"
	    i__2 = j * c_dim1 + 3;
#line 210 "clarfx.f"
	    i__3 = j * c_dim1 + 3;
#line 210 "clarfx.f"
	    z__2.r = sum.r * t3.r - sum.i * t3.i, z__2.i = sum.r * t3.i + 
		    sum.i * t3.r;
#line 210 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 210 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 211 "clarfx.f"
/* L60: */
#line 211 "clarfx.f"
	}
#line 212 "clarfx.f"
	goto L410;
#line 213 "clarfx.f"
L70:

/*        Special code for 4 x 4 Householder */

#line 217 "clarfx.f"
	d_cnjg(&z__1, &v[1]);
#line 217 "clarfx.f"
	v1.r = z__1.r, v1.i = z__1.i;
#line 218 "clarfx.f"
	d_cnjg(&z__2, &v1);
#line 218 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 218 "clarfx.f"
	t1.r = z__1.r, t1.i = z__1.i;
#line 219 "clarfx.f"
	d_cnjg(&z__1, &v[2]);
#line 219 "clarfx.f"
	v2.r = z__1.r, v2.i = z__1.i;
#line 220 "clarfx.f"
	d_cnjg(&z__2, &v2);
#line 220 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 220 "clarfx.f"
	t2.r = z__1.r, t2.i = z__1.i;
#line 221 "clarfx.f"
	d_cnjg(&z__1, &v[3]);
#line 221 "clarfx.f"
	v3.r = z__1.r, v3.i = z__1.i;
#line 222 "clarfx.f"
	d_cnjg(&z__2, &v3);
#line 222 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 222 "clarfx.f"
	t3.r = z__1.r, t3.i = z__1.i;
#line 223 "clarfx.f"
	d_cnjg(&z__1, &v[4]);
#line 223 "clarfx.f"
	v4.r = z__1.r, v4.i = z__1.i;
#line 224 "clarfx.f"
	d_cnjg(&z__2, &v4);
#line 224 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 224 "clarfx.f"
	t4.r = z__1.r, t4.i = z__1.i;
#line 225 "clarfx.f"
	i__1 = *n;
#line 225 "clarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 226 "clarfx.f"
	    i__2 = j * c_dim1 + 1;
#line 226 "clarfx.f"
	    z__4.r = v1.r * c__[i__2].r - v1.i * c__[i__2].i, z__4.i = v1.r * 
		    c__[i__2].i + v1.i * c__[i__2].r;
#line 226 "clarfx.f"
	    i__3 = j * c_dim1 + 2;
#line 226 "clarfx.f"
	    z__5.r = v2.r * c__[i__3].r - v2.i * c__[i__3].i, z__5.i = v2.r * 
		    c__[i__3].i + v2.i * c__[i__3].r;
#line 226 "clarfx.f"
	    z__3.r = z__4.r + z__5.r, z__3.i = z__4.i + z__5.i;
#line 226 "clarfx.f"
	    i__4 = j * c_dim1 + 3;
#line 226 "clarfx.f"
	    z__6.r = v3.r * c__[i__4].r - v3.i * c__[i__4].i, z__6.i = v3.r * 
		    c__[i__4].i + v3.i * c__[i__4].r;
#line 226 "clarfx.f"
	    z__2.r = z__3.r + z__6.r, z__2.i = z__3.i + z__6.i;
#line 226 "clarfx.f"
	    i__5 = j * c_dim1 + 4;
#line 226 "clarfx.f"
	    z__7.r = v4.r * c__[i__5].r - v4.i * c__[i__5].i, z__7.i = v4.r * 
		    c__[i__5].i + v4.i * c__[i__5].r;
#line 226 "clarfx.f"
	    z__1.r = z__2.r + z__7.r, z__1.i = z__2.i + z__7.i;
#line 226 "clarfx.f"
	    sum.r = z__1.r, sum.i = z__1.i;
#line 228 "clarfx.f"
	    i__2 = j * c_dim1 + 1;
#line 228 "clarfx.f"
	    i__3 = j * c_dim1 + 1;
#line 228 "clarfx.f"
	    z__2.r = sum.r * t1.r - sum.i * t1.i, z__2.i = sum.r * t1.i + 
		    sum.i * t1.r;
#line 228 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 228 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 229 "clarfx.f"
	    i__2 = j * c_dim1 + 2;
#line 229 "clarfx.f"
	    i__3 = j * c_dim1 + 2;
#line 229 "clarfx.f"
	    z__2.r = sum.r * t2.r - sum.i * t2.i, z__2.i = sum.r * t2.i + 
		    sum.i * t2.r;
#line 229 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 229 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 230 "clarfx.f"
	    i__2 = j * c_dim1 + 3;
#line 230 "clarfx.f"
	    i__3 = j * c_dim1 + 3;
#line 230 "clarfx.f"
	    z__2.r = sum.r * t3.r - sum.i * t3.i, z__2.i = sum.r * t3.i + 
		    sum.i * t3.r;
#line 230 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 230 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 231 "clarfx.f"
	    i__2 = j * c_dim1 + 4;
#line 231 "clarfx.f"
	    i__3 = j * c_dim1 + 4;
#line 231 "clarfx.f"
	    z__2.r = sum.r * t4.r - sum.i * t4.i, z__2.i = sum.r * t4.i + 
		    sum.i * t4.r;
#line 231 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 231 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 232 "clarfx.f"
/* L80: */
#line 232 "clarfx.f"
	}
#line 233 "clarfx.f"
	goto L410;
#line 234 "clarfx.f"
L90:

/*        Special code for 5 x 5 Householder */

#line 238 "clarfx.f"
	d_cnjg(&z__1, &v[1]);
#line 238 "clarfx.f"
	v1.r = z__1.r, v1.i = z__1.i;
#line 239 "clarfx.f"
	d_cnjg(&z__2, &v1);
#line 239 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 239 "clarfx.f"
	t1.r = z__1.r, t1.i = z__1.i;
#line 240 "clarfx.f"
	d_cnjg(&z__1, &v[2]);
#line 240 "clarfx.f"
	v2.r = z__1.r, v2.i = z__1.i;
#line 241 "clarfx.f"
	d_cnjg(&z__2, &v2);
#line 241 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 241 "clarfx.f"
	t2.r = z__1.r, t2.i = z__1.i;
#line 242 "clarfx.f"
	d_cnjg(&z__1, &v[3]);
#line 242 "clarfx.f"
	v3.r = z__1.r, v3.i = z__1.i;
#line 243 "clarfx.f"
	d_cnjg(&z__2, &v3);
#line 243 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 243 "clarfx.f"
	t3.r = z__1.r, t3.i = z__1.i;
#line 244 "clarfx.f"
	d_cnjg(&z__1, &v[4]);
#line 244 "clarfx.f"
	v4.r = z__1.r, v4.i = z__1.i;
#line 245 "clarfx.f"
	d_cnjg(&z__2, &v4);
#line 245 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 245 "clarfx.f"
	t4.r = z__1.r, t4.i = z__1.i;
#line 246 "clarfx.f"
	d_cnjg(&z__1, &v[5]);
#line 246 "clarfx.f"
	v5.r = z__1.r, v5.i = z__1.i;
#line 247 "clarfx.f"
	d_cnjg(&z__2, &v5);
#line 247 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 247 "clarfx.f"
	t5.r = z__1.r, t5.i = z__1.i;
#line 248 "clarfx.f"
	i__1 = *n;
#line 248 "clarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 249 "clarfx.f"
	    i__2 = j * c_dim1 + 1;
#line 249 "clarfx.f"
	    z__5.r = v1.r * c__[i__2].r - v1.i * c__[i__2].i, z__5.i = v1.r * 
		    c__[i__2].i + v1.i * c__[i__2].r;
#line 249 "clarfx.f"
	    i__3 = j * c_dim1 + 2;
#line 249 "clarfx.f"
	    z__6.r = v2.r * c__[i__3].r - v2.i * c__[i__3].i, z__6.i = v2.r * 
		    c__[i__3].i + v2.i * c__[i__3].r;
#line 249 "clarfx.f"
	    z__4.r = z__5.r + z__6.r, z__4.i = z__5.i + z__6.i;
#line 249 "clarfx.f"
	    i__4 = j * c_dim1 + 3;
#line 249 "clarfx.f"
	    z__7.r = v3.r * c__[i__4].r - v3.i * c__[i__4].i, z__7.i = v3.r * 
		    c__[i__4].i + v3.i * c__[i__4].r;
#line 249 "clarfx.f"
	    z__3.r = z__4.r + z__7.r, z__3.i = z__4.i + z__7.i;
#line 249 "clarfx.f"
	    i__5 = j * c_dim1 + 4;
#line 249 "clarfx.f"
	    z__8.r = v4.r * c__[i__5].r - v4.i * c__[i__5].i, z__8.i = v4.r * 
		    c__[i__5].i + v4.i * c__[i__5].r;
#line 249 "clarfx.f"
	    z__2.r = z__3.r + z__8.r, z__2.i = z__3.i + z__8.i;
#line 249 "clarfx.f"
	    i__6 = j * c_dim1 + 5;
#line 249 "clarfx.f"
	    z__9.r = v5.r * c__[i__6].r - v5.i * c__[i__6].i, z__9.i = v5.r * 
		    c__[i__6].i + v5.i * c__[i__6].r;
#line 249 "clarfx.f"
	    z__1.r = z__2.r + z__9.r, z__1.i = z__2.i + z__9.i;
#line 249 "clarfx.f"
	    sum.r = z__1.r, sum.i = z__1.i;
#line 251 "clarfx.f"
	    i__2 = j * c_dim1 + 1;
#line 251 "clarfx.f"
	    i__3 = j * c_dim1 + 1;
#line 251 "clarfx.f"
	    z__2.r = sum.r * t1.r - sum.i * t1.i, z__2.i = sum.r * t1.i + 
		    sum.i * t1.r;
#line 251 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 251 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 252 "clarfx.f"
	    i__2 = j * c_dim1 + 2;
#line 252 "clarfx.f"
	    i__3 = j * c_dim1 + 2;
#line 252 "clarfx.f"
	    z__2.r = sum.r * t2.r - sum.i * t2.i, z__2.i = sum.r * t2.i + 
		    sum.i * t2.r;
#line 252 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 252 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 253 "clarfx.f"
	    i__2 = j * c_dim1 + 3;
#line 253 "clarfx.f"
	    i__3 = j * c_dim1 + 3;
#line 253 "clarfx.f"
	    z__2.r = sum.r * t3.r - sum.i * t3.i, z__2.i = sum.r * t3.i + 
		    sum.i * t3.r;
#line 253 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 253 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 254 "clarfx.f"
	    i__2 = j * c_dim1 + 4;
#line 254 "clarfx.f"
	    i__3 = j * c_dim1 + 4;
#line 254 "clarfx.f"
	    z__2.r = sum.r * t4.r - sum.i * t4.i, z__2.i = sum.r * t4.i + 
		    sum.i * t4.r;
#line 254 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 254 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 255 "clarfx.f"
	    i__2 = j * c_dim1 + 5;
#line 255 "clarfx.f"
	    i__3 = j * c_dim1 + 5;
#line 255 "clarfx.f"
	    z__2.r = sum.r * t5.r - sum.i * t5.i, z__2.i = sum.r * t5.i + 
		    sum.i * t5.r;
#line 255 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 255 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 256 "clarfx.f"
/* L100: */
#line 256 "clarfx.f"
	}
#line 257 "clarfx.f"
	goto L410;
#line 258 "clarfx.f"
L110:

/*        Special code for 6 x 6 Householder */

#line 262 "clarfx.f"
	d_cnjg(&z__1, &v[1]);
#line 262 "clarfx.f"
	v1.r = z__1.r, v1.i = z__1.i;
#line 263 "clarfx.f"
	d_cnjg(&z__2, &v1);
#line 263 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 263 "clarfx.f"
	t1.r = z__1.r, t1.i = z__1.i;
#line 264 "clarfx.f"
	d_cnjg(&z__1, &v[2]);
#line 264 "clarfx.f"
	v2.r = z__1.r, v2.i = z__1.i;
#line 265 "clarfx.f"
	d_cnjg(&z__2, &v2);
#line 265 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 265 "clarfx.f"
	t2.r = z__1.r, t2.i = z__1.i;
#line 266 "clarfx.f"
	d_cnjg(&z__1, &v[3]);
#line 266 "clarfx.f"
	v3.r = z__1.r, v3.i = z__1.i;
#line 267 "clarfx.f"
	d_cnjg(&z__2, &v3);
#line 267 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 267 "clarfx.f"
	t3.r = z__1.r, t3.i = z__1.i;
#line 268 "clarfx.f"
	d_cnjg(&z__1, &v[4]);
#line 268 "clarfx.f"
	v4.r = z__1.r, v4.i = z__1.i;
#line 269 "clarfx.f"
	d_cnjg(&z__2, &v4);
#line 269 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 269 "clarfx.f"
	t4.r = z__1.r, t4.i = z__1.i;
#line 270 "clarfx.f"
	d_cnjg(&z__1, &v[5]);
#line 270 "clarfx.f"
	v5.r = z__1.r, v5.i = z__1.i;
#line 271 "clarfx.f"
	d_cnjg(&z__2, &v5);
#line 271 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 271 "clarfx.f"
	t5.r = z__1.r, t5.i = z__1.i;
#line 272 "clarfx.f"
	d_cnjg(&z__1, &v[6]);
#line 272 "clarfx.f"
	v6.r = z__1.r, v6.i = z__1.i;
#line 273 "clarfx.f"
	d_cnjg(&z__2, &v6);
#line 273 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 273 "clarfx.f"
	t6.r = z__1.r, t6.i = z__1.i;
#line 274 "clarfx.f"
	i__1 = *n;
#line 274 "clarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 275 "clarfx.f"
	    i__2 = j * c_dim1 + 1;
#line 275 "clarfx.f"
	    z__6.r = v1.r * c__[i__2].r - v1.i * c__[i__2].i, z__6.i = v1.r * 
		    c__[i__2].i + v1.i * c__[i__2].r;
#line 275 "clarfx.f"
	    i__3 = j * c_dim1 + 2;
#line 275 "clarfx.f"
	    z__7.r = v2.r * c__[i__3].r - v2.i * c__[i__3].i, z__7.i = v2.r * 
		    c__[i__3].i + v2.i * c__[i__3].r;
#line 275 "clarfx.f"
	    z__5.r = z__6.r + z__7.r, z__5.i = z__6.i + z__7.i;
#line 275 "clarfx.f"
	    i__4 = j * c_dim1 + 3;
#line 275 "clarfx.f"
	    z__8.r = v3.r * c__[i__4].r - v3.i * c__[i__4].i, z__8.i = v3.r * 
		    c__[i__4].i + v3.i * c__[i__4].r;
#line 275 "clarfx.f"
	    z__4.r = z__5.r + z__8.r, z__4.i = z__5.i + z__8.i;
#line 275 "clarfx.f"
	    i__5 = j * c_dim1 + 4;
#line 275 "clarfx.f"
	    z__9.r = v4.r * c__[i__5].r - v4.i * c__[i__5].i, z__9.i = v4.r * 
		    c__[i__5].i + v4.i * c__[i__5].r;
#line 275 "clarfx.f"
	    z__3.r = z__4.r + z__9.r, z__3.i = z__4.i + z__9.i;
#line 275 "clarfx.f"
	    i__6 = j * c_dim1 + 5;
#line 275 "clarfx.f"
	    z__10.r = v5.r * c__[i__6].r - v5.i * c__[i__6].i, z__10.i = v5.r 
		    * c__[i__6].i + v5.i * c__[i__6].r;
#line 275 "clarfx.f"
	    z__2.r = z__3.r + z__10.r, z__2.i = z__3.i + z__10.i;
#line 275 "clarfx.f"
	    i__7 = j * c_dim1 + 6;
#line 275 "clarfx.f"
	    z__11.r = v6.r * c__[i__7].r - v6.i * c__[i__7].i, z__11.i = v6.r 
		    * c__[i__7].i + v6.i * c__[i__7].r;
#line 275 "clarfx.f"
	    z__1.r = z__2.r + z__11.r, z__1.i = z__2.i + z__11.i;
#line 275 "clarfx.f"
	    sum.r = z__1.r, sum.i = z__1.i;
#line 277 "clarfx.f"
	    i__2 = j * c_dim1 + 1;
#line 277 "clarfx.f"
	    i__3 = j * c_dim1 + 1;
#line 277 "clarfx.f"
	    z__2.r = sum.r * t1.r - sum.i * t1.i, z__2.i = sum.r * t1.i + 
		    sum.i * t1.r;
#line 277 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 277 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 278 "clarfx.f"
	    i__2 = j * c_dim1 + 2;
#line 278 "clarfx.f"
	    i__3 = j * c_dim1 + 2;
#line 278 "clarfx.f"
	    z__2.r = sum.r * t2.r - sum.i * t2.i, z__2.i = sum.r * t2.i + 
		    sum.i * t2.r;
#line 278 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 278 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 279 "clarfx.f"
	    i__2 = j * c_dim1 + 3;
#line 279 "clarfx.f"
	    i__3 = j * c_dim1 + 3;
#line 279 "clarfx.f"
	    z__2.r = sum.r * t3.r - sum.i * t3.i, z__2.i = sum.r * t3.i + 
		    sum.i * t3.r;
#line 279 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 279 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 280 "clarfx.f"
	    i__2 = j * c_dim1 + 4;
#line 280 "clarfx.f"
	    i__3 = j * c_dim1 + 4;
#line 280 "clarfx.f"
	    z__2.r = sum.r * t4.r - sum.i * t4.i, z__2.i = sum.r * t4.i + 
		    sum.i * t4.r;
#line 280 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 280 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 281 "clarfx.f"
	    i__2 = j * c_dim1 + 5;
#line 281 "clarfx.f"
	    i__3 = j * c_dim1 + 5;
#line 281 "clarfx.f"
	    z__2.r = sum.r * t5.r - sum.i * t5.i, z__2.i = sum.r * t5.i + 
		    sum.i * t5.r;
#line 281 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 281 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 282 "clarfx.f"
	    i__2 = j * c_dim1 + 6;
#line 282 "clarfx.f"
	    i__3 = j * c_dim1 + 6;
#line 282 "clarfx.f"
	    z__2.r = sum.r * t6.r - sum.i * t6.i, z__2.i = sum.r * t6.i + 
		    sum.i * t6.r;
#line 282 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 282 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 283 "clarfx.f"
/* L120: */
#line 283 "clarfx.f"
	}
#line 284 "clarfx.f"
	goto L410;
#line 285 "clarfx.f"
L130:

/*        Special code for 7 x 7 Householder */

#line 289 "clarfx.f"
	d_cnjg(&z__1, &v[1]);
#line 289 "clarfx.f"
	v1.r = z__1.r, v1.i = z__1.i;
#line 290 "clarfx.f"
	d_cnjg(&z__2, &v1);
#line 290 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 290 "clarfx.f"
	t1.r = z__1.r, t1.i = z__1.i;
#line 291 "clarfx.f"
	d_cnjg(&z__1, &v[2]);
#line 291 "clarfx.f"
	v2.r = z__1.r, v2.i = z__1.i;
#line 292 "clarfx.f"
	d_cnjg(&z__2, &v2);
#line 292 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 292 "clarfx.f"
	t2.r = z__1.r, t2.i = z__1.i;
#line 293 "clarfx.f"
	d_cnjg(&z__1, &v[3]);
#line 293 "clarfx.f"
	v3.r = z__1.r, v3.i = z__1.i;
#line 294 "clarfx.f"
	d_cnjg(&z__2, &v3);
#line 294 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 294 "clarfx.f"
	t3.r = z__1.r, t3.i = z__1.i;
#line 295 "clarfx.f"
	d_cnjg(&z__1, &v[4]);
#line 295 "clarfx.f"
	v4.r = z__1.r, v4.i = z__1.i;
#line 296 "clarfx.f"
	d_cnjg(&z__2, &v4);
#line 296 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 296 "clarfx.f"
	t4.r = z__1.r, t4.i = z__1.i;
#line 297 "clarfx.f"
	d_cnjg(&z__1, &v[5]);
#line 297 "clarfx.f"
	v5.r = z__1.r, v5.i = z__1.i;
#line 298 "clarfx.f"
	d_cnjg(&z__2, &v5);
#line 298 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 298 "clarfx.f"
	t5.r = z__1.r, t5.i = z__1.i;
#line 299 "clarfx.f"
	d_cnjg(&z__1, &v[6]);
#line 299 "clarfx.f"
	v6.r = z__1.r, v6.i = z__1.i;
#line 300 "clarfx.f"
	d_cnjg(&z__2, &v6);
#line 300 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 300 "clarfx.f"
	t6.r = z__1.r, t6.i = z__1.i;
#line 301 "clarfx.f"
	d_cnjg(&z__1, &v[7]);
#line 301 "clarfx.f"
	v7.r = z__1.r, v7.i = z__1.i;
#line 302 "clarfx.f"
	d_cnjg(&z__2, &v7);
#line 302 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 302 "clarfx.f"
	t7.r = z__1.r, t7.i = z__1.i;
#line 303 "clarfx.f"
	i__1 = *n;
#line 303 "clarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 304 "clarfx.f"
	    i__2 = j * c_dim1 + 1;
#line 304 "clarfx.f"
	    z__7.r = v1.r * c__[i__2].r - v1.i * c__[i__2].i, z__7.i = v1.r * 
		    c__[i__2].i + v1.i * c__[i__2].r;
#line 304 "clarfx.f"
	    i__3 = j * c_dim1 + 2;
#line 304 "clarfx.f"
	    z__8.r = v2.r * c__[i__3].r - v2.i * c__[i__3].i, z__8.i = v2.r * 
		    c__[i__3].i + v2.i * c__[i__3].r;
#line 304 "clarfx.f"
	    z__6.r = z__7.r + z__8.r, z__6.i = z__7.i + z__8.i;
#line 304 "clarfx.f"
	    i__4 = j * c_dim1 + 3;
#line 304 "clarfx.f"
	    z__9.r = v3.r * c__[i__4].r - v3.i * c__[i__4].i, z__9.i = v3.r * 
		    c__[i__4].i + v3.i * c__[i__4].r;
#line 304 "clarfx.f"
	    z__5.r = z__6.r + z__9.r, z__5.i = z__6.i + z__9.i;
#line 304 "clarfx.f"
	    i__5 = j * c_dim1 + 4;
#line 304 "clarfx.f"
	    z__10.r = v4.r * c__[i__5].r - v4.i * c__[i__5].i, z__10.i = v4.r 
		    * c__[i__5].i + v4.i * c__[i__5].r;
#line 304 "clarfx.f"
	    z__4.r = z__5.r + z__10.r, z__4.i = z__5.i + z__10.i;
#line 304 "clarfx.f"
	    i__6 = j * c_dim1 + 5;
#line 304 "clarfx.f"
	    z__11.r = v5.r * c__[i__6].r - v5.i * c__[i__6].i, z__11.i = v5.r 
		    * c__[i__6].i + v5.i * c__[i__6].r;
#line 304 "clarfx.f"
	    z__3.r = z__4.r + z__11.r, z__3.i = z__4.i + z__11.i;
#line 304 "clarfx.f"
	    i__7 = j * c_dim1 + 6;
#line 304 "clarfx.f"
	    z__12.r = v6.r * c__[i__7].r - v6.i * c__[i__7].i, z__12.i = v6.r 
		    * c__[i__7].i + v6.i * c__[i__7].r;
#line 304 "clarfx.f"
	    z__2.r = z__3.r + z__12.r, z__2.i = z__3.i + z__12.i;
#line 304 "clarfx.f"
	    i__8 = j * c_dim1 + 7;
#line 304 "clarfx.f"
	    z__13.r = v7.r * c__[i__8].r - v7.i * c__[i__8].i, z__13.i = v7.r 
		    * c__[i__8].i + v7.i * c__[i__8].r;
#line 304 "clarfx.f"
	    z__1.r = z__2.r + z__13.r, z__1.i = z__2.i + z__13.i;
#line 304 "clarfx.f"
	    sum.r = z__1.r, sum.i = z__1.i;
#line 307 "clarfx.f"
	    i__2 = j * c_dim1 + 1;
#line 307 "clarfx.f"
	    i__3 = j * c_dim1 + 1;
#line 307 "clarfx.f"
	    z__2.r = sum.r * t1.r - sum.i * t1.i, z__2.i = sum.r * t1.i + 
		    sum.i * t1.r;
#line 307 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 307 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 308 "clarfx.f"
	    i__2 = j * c_dim1 + 2;
#line 308 "clarfx.f"
	    i__3 = j * c_dim1 + 2;
#line 308 "clarfx.f"
	    z__2.r = sum.r * t2.r - sum.i * t2.i, z__2.i = sum.r * t2.i + 
		    sum.i * t2.r;
#line 308 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 308 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 309 "clarfx.f"
	    i__2 = j * c_dim1 + 3;
#line 309 "clarfx.f"
	    i__3 = j * c_dim1 + 3;
#line 309 "clarfx.f"
	    z__2.r = sum.r * t3.r - sum.i * t3.i, z__2.i = sum.r * t3.i + 
		    sum.i * t3.r;
#line 309 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 309 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 310 "clarfx.f"
	    i__2 = j * c_dim1 + 4;
#line 310 "clarfx.f"
	    i__3 = j * c_dim1 + 4;
#line 310 "clarfx.f"
	    z__2.r = sum.r * t4.r - sum.i * t4.i, z__2.i = sum.r * t4.i + 
		    sum.i * t4.r;
#line 310 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 310 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 311 "clarfx.f"
	    i__2 = j * c_dim1 + 5;
#line 311 "clarfx.f"
	    i__3 = j * c_dim1 + 5;
#line 311 "clarfx.f"
	    z__2.r = sum.r * t5.r - sum.i * t5.i, z__2.i = sum.r * t5.i + 
		    sum.i * t5.r;
#line 311 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 311 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 312 "clarfx.f"
	    i__2 = j * c_dim1 + 6;
#line 312 "clarfx.f"
	    i__3 = j * c_dim1 + 6;
#line 312 "clarfx.f"
	    z__2.r = sum.r * t6.r - sum.i * t6.i, z__2.i = sum.r * t6.i + 
		    sum.i * t6.r;
#line 312 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 312 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 313 "clarfx.f"
	    i__2 = j * c_dim1 + 7;
#line 313 "clarfx.f"
	    i__3 = j * c_dim1 + 7;
#line 313 "clarfx.f"
	    z__2.r = sum.r * t7.r - sum.i * t7.i, z__2.i = sum.r * t7.i + 
		    sum.i * t7.r;
#line 313 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 313 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 314 "clarfx.f"
/* L140: */
#line 314 "clarfx.f"
	}
#line 315 "clarfx.f"
	goto L410;
#line 316 "clarfx.f"
L150:

/*        Special code for 8 x 8 Householder */

#line 320 "clarfx.f"
	d_cnjg(&z__1, &v[1]);
#line 320 "clarfx.f"
	v1.r = z__1.r, v1.i = z__1.i;
#line 321 "clarfx.f"
	d_cnjg(&z__2, &v1);
#line 321 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 321 "clarfx.f"
	t1.r = z__1.r, t1.i = z__1.i;
#line 322 "clarfx.f"
	d_cnjg(&z__1, &v[2]);
#line 322 "clarfx.f"
	v2.r = z__1.r, v2.i = z__1.i;
#line 323 "clarfx.f"
	d_cnjg(&z__2, &v2);
#line 323 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 323 "clarfx.f"
	t2.r = z__1.r, t2.i = z__1.i;
#line 324 "clarfx.f"
	d_cnjg(&z__1, &v[3]);
#line 324 "clarfx.f"
	v3.r = z__1.r, v3.i = z__1.i;
#line 325 "clarfx.f"
	d_cnjg(&z__2, &v3);
#line 325 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 325 "clarfx.f"
	t3.r = z__1.r, t3.i = z__1.i;
#line 326 "clarfx.f"
	d_cnjg(&z__1, &v[4]);
#line 326 "clarfx.f"
	v4.r = z__1.r, v4.i = z__1.i;
#line 327 "clarfx.f"
	d_cnjg(&z__2, &v4);
#line 327 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 327 "clarfx.f"
	t4.r = z__1.r, t4.i = z__1.i;
#line 328 "clarfx.f"
	d_cnjg(&z__1, &v[5]);
#line 328 "clarfx.f"
	v5.r = z__1.r, v5.i = z__1.i;
#line 329 "clarfx.f"
	d_cnjg(&z__2, &v5);
#line 329 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 329 "clarfx.f"
	t5.r = z__1.r, t5.i = z__1.i;
#line 330 "clarfx.f"
	d_cnjg(&z__1, &v[6]);
#line 330 "clarfx.f"
	v6.r = z__1.r, v6.i = z__1.i;
#line 331 "clarfx.f"
	d_cnjg(&z__2, &v6);
#line 331 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 331 "clarfx.f"
	t6.r = z__1.r, t6.i = z__1.i;
#line 332 "clarfx.f"
	d_cnjg(&z__1, &v[7]);
#line 332 "clarfx.f"
	v7.r = z__1.r, v7.i = z__1.i;
#line 333 "clarfx.f"
	d_cnjg(&z__2, &v7);
#line 333 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 333 "clarfx.f"
	t7.r = z__1.r, t7.i = z__1.i;
#line 334 "clarfx.f"
	d_cnjg(&z__1, &v[8]);
#line 334 "clarfx.f"
	v8.r = z__1.r, v8.i = z__1.i;
#line 335 "clarfx.f"
	d_cnjg(&z__2, &v8);
#line 335 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 335 "clarfx.f"
	t8.r = z__1.r, t8.i = z__1.i;
#line 336 "clarfx.f"
	i__1 = *n;
#line 336 "clarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 337 "clarfx.f"
	    i__2 = j * c_dim1 + 1;
#line 337 "clarfx.f"
	    z__8.r = v1.r * c__[i__2].r - v1.i * c__[i__2].i, z__8.i = v1.r * 
		    c__[i__2].i + v1.i * c__[i__2].r;
#line 337 "clarfx.f"
	    i__3 = j * c_dim1 + 2;
#line 337 "clarfx.f"
	    z__9.r = v2.r * c__[i__3].r - v2.i * c__[i__3].i, z__9.i = v2.r * 
		    c__[i__3].i + v2.i * c__[i__3].r;
#line 337 "clarfx.f"
	    z__7.r = z__8.r + z__9.r, z__7.i = z__8.i + z__9.i;
#line 337 "clarfx.f"
	    i__4 = j * c_dim1 + 3;
#line 337 "clarfx.f"
	    z__10.r = v3.r * c__[i__4].r - v3.i * c__[i__4].i, z__10.i = v3.r 
		    * c__[i__4].i + v3.i * c__[i__4].r;
#line 337 "clarfx.f"
	    z__6.r = z__7.r + z__10.r, z__6.i = z__7.i + z__10.i;
#line 337 "clarfx.f"
	    i__5 = j * c_dim1 + 4;
#line 337 "clarfx.f"
	    z__11.r = v4.r * c__[i__5].r - v4.i * c__[i__5].i, z__11.i = v4.r 
		    * c__[i__5].i + v4.i * c__[i__5].r;
#line 337 "clarfx.f"
	    z__5.r = z__6.r + z__11.r, z__5.i = z__6.i + z__11.i;
#line 337 "clarfx.f"
	    i__6 = j * c_dim1 + 5;
#line 337 "clarfx.f"
	    z__12.r = v5.r * c__[i__6].r - v5.i * c__[i__6].i, z__12.i = v5.r 
		    * c__[i__6].i + v5.i * c__[i__6].r;
#line 337 "clarfx.f"
	    z__4.r = z__5.r + z__12.r, z__4.i = z__5.i + z__12.i;
#line 337 "clarfx.f"
	    i__7 = j * c_dim1 + 6;
#line 337 "clarfx.f"
	    z__13.r = v6.r * c__[i__7].r - v6.i * c__[i__7].i, z__13.i = v6.r 
		    * c__[i__7].i + v6.i * c__[i__7].r;
#line 337 "clarfx.f"
	    z__3.r = z__4.r + z__13.r, z__3.i = z__4.i + z__13.i;
#line 337 "clarfx.f"
	    i__8 = j * c_dim1 + 7;
#line 337 "clarfx.f"
	    z__14.r = v7.r * c__[i__8].r - v7.i * c__[i__8].i, z__14.i = v7.r 
		    * c__[i__8].i + v7.i * c__[i__8].r;
#line 337 "clarfx.f"
	    z__2.r = z__3.r + z__14.r, z__2.i = z__3.i + z__14.i;
#line 337 "clarfx.f"
	    i__9 = j * c_dim1 + 8;
#line 337 "clarfx.f"
	    z__15.r = v8.r * c__[i__9].r - v8.i * c__[i__9].i, z__15.i = v8.r 
		    * c__[i__9].i + v8.i * c__[i__9].r;
#line 337 "clarfx.f"
	    z__1.r = z__2.r + z__15.r, z__1.i = z__2.i + z__15.i;
#line 337 "clarfx.f"
	    sum.r = z__1.r, sum.i = z__1.i;
#line 340 "clarfx.f"
	    i__2 = j * c_dim1 + 1;
#line 340 "clarfx.f"
	    i__3 = j * c_dim1 + 1;
#line 340 "clarfx.f"
	    z__2.r = sum.r * t1.r - sum.i * t1.i, z__2.i = sum.r * t1.i + 
		    sum.i * t1.r;
#line 340 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 340 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 341 "clarfx.f"
	    i__2 = j * c_dim1 + 2;
#line 341 "clarfx.f"
	    i__3 = j * c_dim1 + 2;
#line 341 "clarfx.f"
	    z__2.r = sum.r * t2.r - sum.i * t2.i, z__2.i = sum.r * t2.i + 
		    sum.i * t2.r;
#line 341 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 341 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 342 "clarfx.f"
	    i__2 = j * c_dim1 + 3;
#line 342 "clarfx.f"
	    i__3 = j * c_dim1 + 3;
#line 342 "clarfx.f"
	    z__2.r = sum.r * t3.r - sum.i * t3.i, z__2.i = sum.r * t3.i + 
		    sum.i * t3.r;
#line 342 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 342 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 343 "clarfx.f"
	    i__2 = j * c_dim1 + 4;
#line 343 "clarfx.f"
	    i__3 = j * c_dim1 + 4;
#line 343 "clarfx.f"
	    z__2.r = sum.r * t4.r - sum.i * t4.i, z__2.i = sum.r * t4.i + 
		    sum.i * t4.r;
#line 343 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 343 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 344 "clarfx.f"
	    i__2 = j * c_dim1 + 5;
#line 344 "clarfx.f"
	    i__3 = j * c_dim1 + 5;
#line 344 "clarfx.f"
	    z__2.r = sum.r * t5.r - sum.i * t5.i, z__2.i = sum.r * t5.i + 
		    sum.i * t5.r;
#line 344 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 344 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 345 "clarfx.f"
	    i__2 = j * c_dim1 + 6;
#line 345 "clarfx.f"
	    i__3 = j * c_dim1 + 6;
#line 345 "clarfx.f"
	    z__2.r = sum.r * t6.r - sum.i * t6.i, z__2.i = sum.r * t6.i + 
		    sum.i * t6.r;
#line 345 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 345 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 346 "clarfx.f"
	    i__2 = j * c_dim1 + 7;
#line 346 "clarfx.f"
	    i__3 = j * c_dim1 + 7;
#line 346 "clarfx.f"
	    z__2.r = sum.r * t7.r - sum.i * t7.i, z__2.i = sum.r * t7.i + 
		    sum.i * t7.r;
#line 346 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 346 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 347 "clarfx.f"
	    i__2 = j * c_dim1 + 8;
#line 347 "clarfx.f"
	    i__3 = j * c_dim1 + 8;
#line 347 "clarfx.f"
	    z__2.r = sum.r * t8.r - sum.i * t8.i, z__2.i = sum.r * t8.i + 
		    sum.i * t8.r;
#line 347 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 347 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 348 "clarfx.f"
/* L160: */
#line 348 "clarfx.f"
	}
#line 349 "clarfx.f"
	goto L410;
#line 350 "clarfx.f"
L170:

/*        Special code for 9 x 9 Householder */

#line 354 "clarfx.f"
	d_cnjg(&z__1, &v[1]);
#line 354 "clarfx.f"
	v1.r = z__1.r, v1.i = z__1.i;
#line 355 "clarfx.f"
	d_cnjg(&z__2, &v1);
#line 355 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 355 "clarfx.f"
	t1.r = z__1.r, t1.i = z__1.i;
#line 356 "clarfx.f"
	d_cnjg(&z__1, &v[2]);
#line 356 "clarfx.f"
	v2.r = z__1.r, v2.i = z__1.i;
#line 357 "clarfx.f"
	d_cnjg(&z__2, &v2);
#line 357 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 357 "clarfx.f"
	t2.r = z__1.r, t2.i = z__1.i;
#line 358 "clarfx.f"
	d_cnjg(&z__1, &v[3]);
#line 358 "clarfx.f"
	v3.r = z__1.r, v3.i = z__1.i;
#line 359 "clarfx.f"
	d_cnjg(&z__2, &v3);
#line 359 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 359 "clarfx.f"
	t3.r = z__1.r, t3.i = z__1.i;
#line 360 "clarfx.f"
	d_cnjg(&z__1, &v[4]);
#line 360 "clarfx.f"
	v4.r = z__1.r, v4.i = z__1.i;
#line 361 "clarfx.f"
	d_cnjg(&z__2, &v4);
#line 361 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 361 "clarfx.f"
	t4.r = z__1.r, t4.i = z__1.i;
#line 362 "clarfx.f"
	d_cnjg(&z__1, &v[5]);
#line 362 "clarfx.f"
	v5.r = z__1.r, v5.i = z__1.i;
#line 363 "clarfx.f"
	d_cnjg(&z__2, &v5);
#line 363 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 363 "clarfx.f"
	t5.r = z__1.r, t5.i = z__1.i;
#line 364 "clarfx.f"
	d_cnjg(&z__1, &v[6]);
#line 364 "clarfx.f"
	v6.r = z__1.r, v6.i = z__1.i;
#line 365 "clarfx.f"
	d_cnjg(&z__2, &v6);
#line 365 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 365 "clarfx.f"
	t6.r = z__1.r, t6.i = z__1.i;
#line 366 "clarfx.f"
	d_cnjg(&z__1, &v[7]);
#line 366 "clarfx.f"
	v7.r = z__1.r, v7.i = z__1.i;
#line 367 "clarfx.f"
	d_cnjg(&z__2, &v7);
#line 367 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 367 "clarfx.f"
	t7.r = z__1.r, t7.i = z__1.i;
#line 368 "clarfx.f"
	d_cnjg(&z__1, &v[8]);
#line 368 "clarfx.f"
	v8.r = z__1.r, v8.i = z__1.i;
#line 369 "clarfx.f"
	d_cnjg(&z__2, &v8);
#line 369 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 369 "clarfx.f"
	t8.r = z__1.r, t8.i = z__1.i;
#line 370 "clarfx.f"
	d_cnjg(&z__1, &v[9]);
#line 370 "clarfx.f"
	v9.r = z__1.r, v9.i = z__1.i;
#line 371 "clarfx.f"
	d_cnjg(&z__2, &v9);
#line 371 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 371 "clarfx.f"
	t9.r = z__1.r, t9.i = z__1.i;
#line 372 "clarfx.f"
	i__1 = *n;
#line 372 "clarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 373 "clarfx.f"
	    i__2 = j * c_dim1 + 1;
#line 373 "clarfx.f"
	    z__9.r = v1.r * c__[i__2].r - v1.i * c__[i__2].i, z__9.i = v1.r * 
		    c__[i__2].i + v1.i * c__[i__2].r;
#line 373 "clarfx.f"
	    i__3 = j * c_dim1 + 2;
#line 373 "clarfx.f"
	    z__10.r = v2.r * c__[i__3].r - v2.i * c__[i__3].i, z__10.i = v2.r 
		    * c__[i__3].i + v2.i * c__[i__3].r;
#line 373 "clarfx.f"
	    z__8.r = z__9.r + z__10.r, z__8.i = z__9.i + z__10.i;
#line 373 "clarfx.f"
	    i__4 = j * c_dim1 + 3;
#line 373 "clarfx.f"
	    z__11.r = v3.r * c__[i__4].r - v3.i * c__[i__4].i, z__11.i = v3.r 
		    * c__[i__4].i + v3.i * c__[i__4].r;
#line 373 "clarfx.f"
	    z__7.r = z__8.r + z__11.r, z__7.i = z__8.i + z__11.i;
#line 373 "clarfx.f"
	    i__5 = j * c_dim1 + 4;
#line 373 "clarfx.f"
	    z__12.r = v4.r * c__[i__5].r - v4.i * c__[i__5].i, z__12.i = v4.r 
		    * c__[i__5].i + v4.i * c__[i__5].r;
#line 373 "clarfx.f"
	    z__6.r = z__7.r + z__12.r, z__6.i = z__7.i + z__12.i;
#line 373 "clarfx.f"
	    i__6 = j * c_dim1 + 5;
#line 373 "clarfx.f"
	    z__13.r = v5.r * c__[i__6].r - v5.i * c__[i__6].i, z__13.i = v5.r 
		    * c__[i__6].i + v5.i * c__[i__6].r;
#line 373 "clarfx.f"
	    z__5.r = z__6.r + z__13.r, z__5.i = z__6.i + z__13.i;
#line 373 "clarfx.f"
	    i__7 = j * c_dim1 + 6;
#line 373 "clarfx.f"
	    z__14.r = v6.r * c__[i__7].r - v6.i * c__[i__7].i, z__14.i = v6.r 
		    * c__[i__7].i + v6.i * c__[i__7].r;
#line 373 "clarfx.f"
	    z__4.r = z__5.r + z__14.r, z__4.i = z__5.i + z__14.i;
#line 373 "clarfx.f"
	    i__8 = j * c_dim1 + 7;
#line 373 "clarfx.f"
	    z__15.r = v7.r * c__[i__8].r - v7.i * c__[i__8].i, z__15.i = v7.r 
		    * c__[i__8].i + v7.i * c__[i__8].r;
#line 373 "clarfx.f"
	    z__3.r = z__4.r + z__15.r, z__3.i = z__4.i + z__15.i;
#line 373 "clarfx.f"
	    i__9 = j * c_dim1 + 8;
#line 373 "clarfx.f"
	    z__16.r = v8.r * c__[i__9].r - v8.i * c__[i__9].i, z__16.i = v8.r 
		    * c__[i__9].i + v8.i * c__[i__9].r;
#line 373 "clarfx.f"
	    z__2.r = z__3.r + z__16.r, z__2.i = z__3.i + z__16.i;
#line 373 "clarfx.f"
	    i__10 = j * c_dim1 + 9;
#line 373 "clarfx.f"
	    z__17.r = v9.r * c__[i__10].r - v9.i * c__[i__10].i, z__17.i = 
		    v9.r * c__[i__10].i + v9.i * c__[i__10].r;
#line 373 "clarfx.f"
	    z__1.r = z__2.r + z__17.r, z__1.i = z__2.i + z__17.i;
#line 373 "clarfx.f"
	    sum.r = z__1.r, sum.i = z__1.i;
#line 376 "clarfx.f"
	    i__2 = j * c_dim1 + 1;
#line 376 "clarfx.f"
	    i__3 = j * c_dim1 + 1;
#line 376 "clarfx.f"
	    z__2.r = sum.r * t1.r - sum.i * t1.i, z__2.i = sum.r * t1.i + 
		    sum.i * t1.r;
#line 376 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 376 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 377 "clarfx.f"
	    i__2 = j * c_dim1 + 2;
#line 377 "clarfx.f"
	    i__3 = j * c_dim1 + 2;
#line 377 "clarfx.f"
	    z__2.r = sum.r * t2.r - sum.i * t2.i, z__2.i = sum.r * t2.i + 
		    sum.i * t2.r;
#line 377 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 377 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 378 "clarfx.f"
	    i__2 = j * c_dim1 + 3;
#line 378 "clarfx.f"
	    i__3 = j * c_dim1 + 3;
#line 378 "clarfx.f"
	    z__2.r = sum.r * t3.r - sum.i * t3.i, z__2.i = sum.r * t3.i + 
		    sum.i * t3.r;
#line 378 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 378 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 379 "clarfx.f"
	    i__2 = j * c_dim1 + 4;
#line 379 "clarfx.f"
	    i__3 = j * c_dim1 + 4;
#line 379 "clarfx.f"
	    z__2.r = sum.r * t4.r - sum.i * t4.i, z__2.i = sum.r * t4.i + 
		    sum.i * t4.r;
#line 379 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 379 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 380 "clarfx.f"
	    i__2 = j * c_dim1 + 5;
#line 380 "clarfx.f"
	    i__3 = j * c_dim1 + 5;
#line 380 "clarfx.f"
	    z__2.r = sum.r * t5.r - sum.i * t5.i, z__2.i = sum.r * t5.i + 
		    sum.i * t5.r;
#line 380 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 380 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 381 "clarfx.f"
	    i__2 = j * c_dim1 + 6;
#line 381 "clarfx.f"
	    i__3 = j * c_dim1 + 6;
#line 381 "clarfx.f"
	    z__2.r = sum.r * t6.r - sum.i * t6.i, z__2.i = sum.r * t6.i + 
		    sum.i * t6.r;
#line 381 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 381 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 382 "clarfx.f"
	    i__2 = j * c_dim1 + 7;
#line 382 "clarfx.f"
	    i__3 = j * c_dim1 + 7;
#line 382 "clarfx.f"
	    z__2.r = sum.r * t7.r - sum.i * t7.i, z__2.i = sum.r * t7.i + 
		    sum.i * t7.r;
#line 382 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 382 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 383 "clarfx.f"
	    i__2 = j * c_dim1 + 8;
#line 383 "clarfx.f"
	    i__3 = j * c_dim1 + 8;
#line 383 "clarfx.f"
	    z__2.r = sum.r * t8.r - sum.i * t8.i, z__2.i = sum.r * t8.i + 
		    sum.i * t8.r;
#line 383 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 383 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 384 "clarfx.f"
	    i__2 = j * c_dim1 + 9;
#line 384 "clarfx.f"
	    i__3 = j * c_dim1 + 9;
#line 384 "clarfx.f"
	    z__2.r = sum.r * t9.r - sum.i * t9.i, z__2.i = sum.r * t9.i + 
		    sum.i * t9.r;
#line 384 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 384 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 385 "clarfx.f"
/* L180: */
#line 385 "clarfx.f"
	}
#line 386 "clarfx.f"
	goto L410;
#line 387 "clarfx.f"
L190:

/*        Special code for 10 x 10 Householder */

#line 391 "clarfx.f"
	d_cnjg(&z__1, &v[1]);
#line 391 "clarfx.f"
	v1.r = z__1.r, v1.i = z__1.i;
#line 392 "clarfx.f"
	d_cnjg(&z__2, &v1);
#line 392 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 392 "clarfx.f"
	t1.r = z__1.r, t1.i = z__1.i;
#line 393 "clarfx.f"
	d_cnjg(&z__1, &v[2]);
#line 393 "clarfx.f"
	v2.r = z__1.r, v2.i = z__1.i;
#line 394 "clarfx.f"
	d_cnjg(&z__2, &v2);
#line 394 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 394 "clarfx.f"
	t2.r = z__1.r, t2.i = z__1.i;
#line 395 "clarfx.f"
	d_cnjg(&z__1, &v[3]);
#line 395 "clarfx.f"
	v3.r = z__1.r, v3.i = z__1.i;
#line 396 "clarfx.f"
	d_cnjg(&z__2, &v3);
#line 396 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 396 "clarfx.f"
	t3.r = z__1.r, t3.i = z__1.i;
#line 397 "clarfx.f"
	d_cnjg(&z__1, &v[4]);
#line 397 "clarfx.f"
	v4.r = z__1.r, v4.i = z__1.i;
#line 398 "clarfx.f"
	d_cnjg(&z__2, &v4);
#line 398 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 398 "clarfx.f"
	t4.r = z__1.r, t4.i = z__1.i;
#line 399 "clarfx.f"
	d_cnjg(&z__1, &v[5]);
#line 399 "clarfx.f"
	v5.r = z__1.r, v5.i = z__1.i;
#line 400 "clarfx.f"
	d_cnjg(&z__2, &v5);
#line 400 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 400 "clarfx.f"
	t5.r = z__1.r, t5.i = z__1.i;
#line 401 "clarfx.f"
	d_cnjg(&z__1, &v[6]);
#line 401 "clarfx.f"
	v6.r = z__1.r, v6.i = z__1.i;
#line 402 "clarfx.f"
	d_cnjg(&z__2, &v6);
#line 402 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 402 "clarfx.f"
	t6.r = z__1.r, t6.i = z__1.i;
#line 403 "clarfx.f"
	d_cnjg(&z__1, &v[7]);
#line 403 "clarfx.f"
	v7.r = z__1.r, v7.i = z__1.i;
#line 404 "clarfx.f"
	d_cnjg(&z__2, &v7);
#line 404 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 404 "clarfx.f"
	t7.r = z__1.r, t7.i = z__1.i;
#line 405 "clarfx.f"
	d_cnjg(&z__1, &v[8]);
#line 405 "clarfx.f"
	v8.r = z__1.r, v8.i = z__1.i;
#line 406 "clarfx.f"
	d_cnjg(&z__2, &v8);
#line 406 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 406 "clarfx.f"
	t8.r = z__1.r, t8.i = z__1.i;
#line 407 "clarfx.f"
	d_cnjg(&z__1, &v[9]);
#line 407 "clarfx.f"
	v9.r = z__1.r, v9.i = z__1.i;
#line 408 "clarfx.f"
	d_cnjg(&z__2, &v9);
#line 408 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 408 "clarfx.f"
	t9.r = z__1.r, t9.i = z__1.i;
#line 409 "clarfx.f"
	d_cnjg(&z__1, &v[10]);
#line 409 "clarfx.f"
	v10.r = z__1.r, v10.i = z__1.i;
#line 410 "clarfx.f"
	d_cnjg(&z__2, &v10);
#line 410 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 410 "clarfx.f"
	t10.r = z__1.r, t10.i = z__1.i;
#line 411 "clarfx.f"
	i__1 = *n;
#line 411 "clarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 412 "clarfx.f"
	    i__2 = j * c_dim1 + 1;
#line 412 "clarfx.f"
	    z__10.r = v1.r * c__[i__2].r - v1.i * c__[i__2].i, z__10.i = v1.r 
		    * c__[i__2].i + v1.i * c__[i__2].r;
#line 412 "clarfx.f"
	    i__3 = j * c_dim1 + 2;
#line 412 "clarfx.f"
	    z__11.r = v2.r * c__[i__3].r - v2.i * c__[i__3].i, z__11.i = v2.r 
		    * c__[i__3].i + v2.i * c__[i__3].r;
#line 412 "clarfx.f"
	    z__9.r = z__10.r + z__11.r, z__9.i = z__10.i + z__11.i;
#line 412 "clarfx.f"
	    i__4 = j * c_dim1 + 3;
#line 412 "clarfx.f"
	    z__12.r = v3.r * c__[i__4].r - v3.i * c__[i__4].i, z__12.i = v3.r 
		    * c__[i__4].i + v3.i * c__[i__4].r;
#line 412 "clarfx.f"
	    z__8.r = z__9.r + z__12.r, z__8.i = z__9.i + z__12.i;
#line 412 "clarfx.f"
	    i__5 = j * c_dim1 + 4;
#line 412 "clarfx.f"
	    z__13.r = v4.r * c__[i__5].r - v4.i * c__[i__5].i, z__13.i = v4.r 
		    * c__[i__5].i + v4.i * c__[i__5].r;
#line 412 "clarfx.f"
	    z__7.r = z__8.r + z__13.r, z__7.i = z__8.i + z__13.i;
#line 412 "clarfx.f"
	    i__6 = j * c_dim1 + 5;
#line 412 "clarfx.f"
	    z__14.r = v5.r * c__[i__6].r - v5.i * c__[i__6].i, z__14.i = v5.r 
		    * c__[i__6].i + v5.i * c__[i__6].r;
#line 412 "clarfx.f"
	    z__6.r = z__7.r + z__14.r, z__6.i = z__7.i + z__14.i;
#line 412 "clarfx.f"
	    i__7 = j * c_dim1 + 6;
#line 412 "clarfx.f"
	    z__15.r = v6.r * c__[i__7].r - v6.i * c__[i__7].i, z__15.i = v6.r 
		    * c__[i__7].i + v6.i * c__[i__7].r;
#line 412 "clarfx.f"
	    z__5.r = z__6.r + z__15.r, z__5.i = z__6.i + z__15.i;
#line 412 "clarfx.f"
	    i__8 = j * c_dim1 + 7;
#line 412 "clarfx.f"
	    z__16.r = v7.r * c__[i__8].r - v7.i * c__[i__8].i, z__16.i = v7.r 
		    * c__[i__8].i + v7.i * c__[i__8].r;
#line 412 "clarfx.f"
	    z__4.r = z__5.r + z__16.r, z__4.i = z__5.i + z__16.i;
#line 412 "clarfx.f"
	    i__9 = j * c_dim1 + 8;
#line 412 "clarfx.f"
	    z__17.r = v8.r * c__[i__9].r - v8.i * c__[i__9].i, z__17.i = v8.r 
		    * c__[i__9].i + v8.i * c__[i__9].r;
#line 412 "clarfx.f"
	    z__3.r = z__4.r + z__17.r, z__3.i = z__4.i + z__17.i;
#line 412 "clarfx.f"
	    i__10 = j * c_dim1 + 9;
#line 412 "clarfx.f"
	    z__18.r = v9.r * c__[i__10].r - v9.i * c__[i__10].i, z__18.i = 
		    v9.r * c__[i__10].i + v9.i * c__[i__10].r;
#line 412 "clarfx.f"
	    z__2.r = z__3.r + z__18.r, z__2.i = z__3.i + z__18.i;
#line 412 "clarfx.f"
	    i__11 = j * c_dim1 + 10;
#line 412 "clarfx.f"
	    z__19.r = v10.r * c__[i__11].r - v10.i * c__[i__11].i, z__19.i = 
		    v10.r * c__[i__11].i + v10.i * c__[i__11].r;
#line 412 "clarfx.f"
	    z__1.r = z__2.r + z__19.r, z__1.i = z__2.i + z__19.i;
#line 412 "clarfx.f"
	    sum.r = z__1.r, sum.i = z__1.i;
#line 416 "clarfx.f"
	    i__2 = j * c_dim1 + 1;
#line 416 "clarfx.f"
	    i__3 = j * c_dim1 + 1;
#line 416 "clarfx.f"
	    z__2.r = sum.r * t1.r - sum.i * t1.i, z__2.i = sum.r * t1.i + 
		    sum.i * t1.r;
#line 416 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 416 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 417 "clarfx.f"
	    i__2 = j * c_dim1 + 2;
#line 417 "clarfx.f"
	    i__3 = j * c_dim1 + 2;
#line 417 "clarfx.f"
	    z__2.r = sum.r * t2.r - sum.i * t2.i, z__2.i = sum.r * t2.i + 
		    sum.i * t2.r;
#line 417 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 417 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 418 "clarfx.f"
	    i__2 = j * c_dim1 + 3;
#line 418 "clarfx.f"
	    i__3 = j * c_dim1 + 3;
#line 418 "clarfx.f"
	    z__2.r = sum.r * t3.r - sum.i * t3.i, z__2.i = sum.r * t3.i + 
		    sum.i * t3.r;
#line 418 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 418 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 419 "clarfx.f"
	    i__2 = j * c_dim1 + 4;
#line 419 "clarfx.f"
	    i__3 = j * c_dim1 + 4;
#line 419 "clarfx.f"
	    z__2.r = sum.r * t4.r - sum.i * t4.i, z__2.i = sum.r * t4.i + 
		    sum.i * t4.r;
#line 419 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 419 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 420 "clarfx.f"
	    i__2 = j * c_dim1 + 5;
#line 420 "clarfx.f"
	    i__3 = j * c_dim1 + 5;
#line 420 "clarfx.f"
	    z__2.r = sum.r * t5.r - sum.i * t5.i, z__2.i = sum.r * t5.i + 
		    sum.i * t5.r;
#line 420 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 420 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 421 "clarfx.f"
	    i__2 = j * c_dim1 + 6;
#line 421 "clarfx.f"
	    i__3 = j * c_dim1 + 6;
#line 421 "clarfx.f"
	    z__2.r = sum.r * t6.r - sum.i * t6.i, z__2.i = sum.r * t6.i + 
		    sum.i * t6.r;
#line 421 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 421 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 422 "clarfx.f"
	    i__2 = j * c_dim1 + 7;
#line 422 "clarfx.f"
	    i__3 = j * c_dim1 + 7;
#line 422 "clarfx.f"
	    z__2.r = sum.r * t7.r - sum.i * t7.i, z__2.i = sum.r * t7.i + 
		    sum.i * t7.r;
#line 422 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 422 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 423 "clarfx.f"
	    i__2 = j * c_dim1 + 8;
#line 423 "clarfx.f"
	    i__3 = j * c_dim1 + 8;
#line 423 "clarfx.f"
	    z__2.r = sum.r * t8.r - sum.i * t8.i, z__2.i = sum.r * t8.i + 
		    sum.i * t8.r;
#line 423 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 423 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 424 "clarfx.f"
	    i__2 = j * c_dim1 + 9;
#line 424 "clarfx.f"
	    i__3 = j * c_dim1 + 9;
#line 424 "clarfx.f"
	    z__2.r = sum.r * t9.r - sum.i * t9.i, z__2.i = sum.r * t9.i + 
		    sum.i * t9.r;
#line 424 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 424 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 425 "clarfx.f"
	    i__2 = j * c_dim1 + 10;
#line 425 "clarfx.f"
	    i__3 = j * c_dim1 + 10;
#line 425 "clarfx.f"
	    z__2.r = sum.r * t10.r - sum.i * t10.i, z__2.i = sum.r * t10.i + 
		    sum.i * t10.r;
#line 425 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 425 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 426 "clarfx.f"
/* L200: */
#line 426 "clarfx.f"
	}
#line 427 "clarfx.f"
	goto L410;
#line 428 "clarfx.f"
    } else {

/*        Form  C * H, where H has order n. */

#line 432 "clarfx.f"
	switch (*n) {
#line 432 "clarfx.f"
	    case 1:  goto L210;
#line 432 "clarfx.f"
	    case 2:  goto L230;
#line 432 "clarfx.f"
	    case 3:  goto L250;
#line 432 "clarfx.f"
	    case 4:  goto L270;
#line 432 "clarfx.f"
	    case 5:  goto L290;
#line 432 "clarfx.f"
	    case 6:  goto L310;
#line 432 "clarfx.f"
	    case 7:  goto L330;
#line 432 "clarfx.f"
	    case 8:  goto L350;
#line 432 "clarfx.f"
	    case 9:  goto L370;
#line 432 "clarfx.f"
	    case 10:  goto L390;
#line 432 "clarfx.f"
	}

/*        Code for general N */

#line 437 "clarfx.f"
	clarf_(side, m, n, &v[1], &c__1, tau, &c__[c_offset], ldc, &work[1], (
		ftnlen)1);
#line 438 "clarfx.f"
	goto L410;
#line 439 "clarfx.f"
L210:

/*        Special code for 1 x 1 Householder */

#line 443 "clarfx.f"
	z__3.r = tau->r * v[1].r - tau->i * v[1].i, z__3.i = tau->r * v[1].i 
		+ tau->i * v[1].r;
#line 443 "clarfx.f"
	d_cnjg(&z__4, &v[1]);
#line 443 "clarfx.f"
	z__2.r = z__3.r * z__4.r - z__3.i * z__4.i, z__2.i = z__3.r * z__4.i 
		+ z__3.i * z__4.r;
#line 443 "clarfx.f"
	z__1.r = 1. - z__2.r, z__1.i = 0. - z__2.i;
#line 443 "clarfx.f"
	t1.r = z__1.r, t1.i = z__1.i;
#line 444 "clarfx.f"
	i__1 = *m;
#line 444 "clarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 445 "clarfx.f"
	    i__2 = j + c_dim1;
#line 445 "clarfx.f"
	    i__3 = j + c_dim1;
#line 445 "clarfx.f"
	    z__1.r = t1.r * c__[i__3].r - t1.i * c__[i__3].i, z__1.i = t1.r * 
		    c__[i__3].i + t1.i * c__[i__3].r;
#line 445 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 446 "clarfx.f"
/* L220: */
#line 446 "clarfx.f"
	}
#line 447 "clarfx.f"
	goto L410;
#line 448 "clarfx.f"
L230:

/*        Special code for 2 x 2 Householder */

#line 452 "clarfx.f"
	v1.r = v[1].r, v1.i = v[1].i;
#line 453 "clarfx.f"
	d_cnjg(&z__2, &v1);
#line 453 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 453 "clarfx.f"
	t1.r = z__1.r, t1.i = z__1.i;
#line 454 "clarfx.f"
	v2.r = v[2].r, v2.i = v[2].i;
#line 455 "clarfx.f"
	d_cnjg(&z__2, &v2);
#line 455 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 455 "clarfx.f"
	t2.r = z__1.r, t2.i = z__1.i;
#line 456 "clarfx.f"
	i__1 = *m;
#line 456 "clarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 457 "clarfx.f"
	    i__2 = j + c_dim1;
#line 457 "clarfx.f"
	    z__2.r = v1.r * c__[i__2].r - v1.i * c__[i__2].i, z__2.i = v1.r * 
		    c__[i__2].i + v1.i * c__[i__2].r;
#line 457 "clarfx.f"
	    i__3 = j + (c_dim1 << 1);
#line 457 "clarfx.f"
	    z__3.r = v2.r * c__[i__3].r - v2.i * c__[i__3].i, z__3.i = v2.r * 
		    c__[i__3].i + v2.i * c__[i__3].r;
#line 457 "clarfx.f"
	    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 457 "clarfx.f"
	    sum.r = z__1.r, sum.i = z__1.i;
#line 458 "clarfx.f"
	    i__2 = j + c_dim1;
#line 458 "clarfx.f"
	    i__3 = j + c_dim1;
#line 458 "clarfx.f"
	    z__2.r = sum.r * t1.r - sum.i * t1.i, z__2.i = sum.r * t1.i + 
		    sum.i * t1.r;
#line 458 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 458 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 459 "clarfx.f"
	    i__2 = j + (c_dim1 << 1);
#line 459 "clarfx.f"
	    i__3 = j + (c_dim1 << 1);
#line 459 "clarfx.f"
	    z__2.r = sum.r * t2.r - sum.i * t2.i, z__2.i = sum.r * t2.i + 
		    sum.i * t2.r;
#line 459 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 459 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 460 "clarfx.f"
/* L240: */
#line 460 "clarfx.f"
	}
#line 461 "clarfx.f"
	goto L410;
#line 462 "clarfx.f"
L250:

/*        Special code for 3 x 3 Householder */

#line 466 "clarfx.f"
	v1.r = v[1].r, v1.i = v[1].i;
#line 467 "clarfx.f"
	d_cnjg(&z__2, &v1);
#line 467 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 467 "clarfx.f"
	t1.r = z__1.r, t1.i = z__1.i;
#line 468 "clarfx.f"
	v2.r = v[2].r, v2.i = v[2].i;
#line 469 "clarfx.f"
	d_cnjg(&z__2, &v2);
#line 469 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 469 "clarfx.f"
	t2.r = z__1.r, t2.i = z__1.i;
#line 470 "clarfx.f"
	v3.r = v[3].r, v3.i = v[3].i;
#line 471 "clarfx.f"
	d_cnjg(&z__2, &v3);
#line 471 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 471 "clarfx.f"
	t3.r = z__1.r, t3.i = z__1.i;
#line 472 "clarfx.f"
	i__1 = *m;
#line 472 "clarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 473 "clarfx.f"
	    i__2 = j + c_dim1;
#line 473 "clarfx.f"
	    z__3.r = v1.r * c__[i__2].r - v1.i * c__[i__2].i, z__3.i = v1.r * 
		    c__[i__2].i + v1.i * c__[i__2].r;
#line 473 "clarfx.f"
	    i__3 = j + (c_dim1 << 1);
#line 473 "clarfx.f"
	    z__4.r = v2.r * c__[i__3].r - v2.i * c__[i__3].i, z__4.i = v2.r * 
		    c__[i__3].i + v2.i * c__[i__3].r;
#line 473 "clarfx.f"
	    z__2.r = z__3.r + z__4.r, z__2.i = z__3.i + z__4.i;
#line 473 "clarfx.f"
	    i__4 = j + c_dim1 * 3;
#line 473 "clarfx.f"
	    z__5.r = v3.r * c__[i__4].r - v3.i * c__[i__4].i, z__5.i = v3.r * 
		    c__[i__4].i + v3.i * c__[i__4].r;
#line 473 "clarfx.f"
	    z__1.r = z__2.r + z__5.r, z__1.i = z__2.i + z__5.i;
#line 473 "clarfx.f"
	    sum.r = z__1.r, sum.i = z__1.i;
#line 474 "clarfx.f"
	    i__2 = j + c_dim1;
#line 474 "clarfx.f"
	    i__3 = j + c_dim1;
#line 474 "clarfx.f"
	    z__2.r = sum.r * t1.r - sum.i * t1.i, z__2.i = sum.r * t1.i + 
		    sum.i * t1.r;
#line 474 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 474 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 475 "clarfx.f"
	    i__2 = j + (c_dim1 << 1);
#line 475 "clarfx.f"
	    i__3 = j + (c_dim1 << 1);
#line 475 "clarfx.f"
	    z__2.r = sum.r * t2.r - sum.i * t2.i, z__2.i = sum.r * t2.i + 
		    sum.i * t2.r;
#line 475 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 475 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 476 "clarfx.f"
	    i__2 = j + c_dim1 * 3;
#line 476 "clarfx.f"
	    i__3 = j + c_dim1 * 3;
#line 476 "clarfx.f"
	    z__2.r = sum.r * t3.r - sum.i * t3.i, z__2.i = sum.r * t3.i + 
		    sum.i * t3.r;
#line 476 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 476 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 477 "clarfx.f"
/* L260: */
#line 477 "clarfx.f"
	}
#line 478 "clarfx.f"
	goto L410;
#line 479 "clarfx.f"
L270:

/*        Special code for 4 x 4 Householder */

#line 483 "clarfx.f"
	v1.r = v[1].r, v1.i = v[1].i;
#line 484 "clarfx.f"
	d_cnjg(&z__2, &v1);
#line 484 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 484 "clarfx.f"
	t1.r = z__1.r, t1.i = z__1.i;
#line 485 "clarfx.f"
	v2.r = v[2].r, v2.i = v[2].i;
#line 486 "clarfx.f"
	d_cnjg(&z__2, &v2);
#line 486 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 486 "clarfx.f"
	t2.r = z__1.r, t2.i = z__1.i;
#line 487 "clarfx.f"
	v3.r = v[3].r, v3.i = v[3].i;
#line 488 "clarfx.f"
	d_cnjg(&z__2, &v3);
#line 488 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 488 "clarfx.f"
	t3.r = z__1.r, t3.i = z__1.i;
#line 489 "clarfx.f"
	v4.r = v[4].r, v4.i = v[4].i;
#line 490 "clarfx.f"
	d_cnjg(&z__2, &v4);
#line 490 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 490 "clarfx.f"
	t4.r = z__1.r, t4.i = z__1.i;
#line 491 "clarfx.f"
	i__1 = *m;
#line 491 "clarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 492 "clarfx.f"
	    i__2 = j + c_dim1;
#line 492 "clarfx.f"
	    z__4.r = v1.r * c__[i__2].r - v1.i * c__[i__2].i, z__4.i = v1.r * 
		    c__[i__2].i + v1.i * c__[i__2].r;
#line 492 "clarfx.f"
	    i__3 = j + (c_dim1 << 1);
#line 492 "clarfx.f"
	    z__5.r = v2.r * c__[i__3].r - v2.i * c__[i__3].i, z__5.i = v2.r * 
		    c__[i__3].i + v2.i * c__[i__3].r;
#line 492 "clarfx.f"
	    z__3.r = z__4.r + z__5.r, z__3.i = z__4.i + z__5.i;
#line 492 "clarfx.f"
	    i__4 = j + c_dim1 * 3;
#line 492 "clarfx.f"
	    z__6.r = v3.r * c__[i__4].r - v3.i * c__[i__4].i, z__6.i = v3.r * 
		    c__[i__4].i + v3.i * c__[i__4].r;
#line 492 "clarfx.f"
	    z__2.r = z__3.r + z__6.r, z__2.i = z__3.i + z__6.i;
#line 492 "clarfx.f"
	    i__5 = j + (c_dim1 << 2);
#line 492 "clarfx.f"
	    z__7.r = v4.r * c__[i__5].r - v4.i * c__[i__5].i, z__7.i = v4.r * 
		    c__[i__5].i + v4.i * c__[i__5].r;
#line 492 "clarfx.f"
	    z__1.r = z__2.r + z__7.r, z__1.i = z__2.i + z__7.i;
#line 492 "clarfx.f"
	    sum.r = z__1.r, sum.i = z__1.i;
#line 494 "clarfx.f"
	    i__2 = j + c_dim1;
#line 494 "clarfx.f"
	    i__3 = j + c_dim1;
#line 494 "clarfx.f"
	    z__2.r = sum.r * t1.r - sum.i * t1.i, z__2.i = sum.r * t1.i + 
		    sum.i * t1.r;
#line 494 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 494 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 495 "clarfx.f"
	    i__2 = j + (c_dim1 << 1);
#line 495 "clarfx.f"
	    i__3 = j + (c_dim1 << 1);
#line 495 "clarfx.f"
	    z__2.r = sum.r * t2.r - sum.i * t2.i, z__2.i = sum.r * t2.i + 
		    sum.i * t2.r;
#line 495 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 495 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 496 "clarfx.f"
	    i__2 = j + c_dim1 * 3;
#line 496 "clarfx.f"
	    i__3 = j + c_dim1 * 3;
#line 496 "clarfx.f"
	    z__2.r = sum.r * t3.r - sum.i * t3.i, z__2.i = sum.r * t3.i + 
		    sum.i * t3.r;
#line 496 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 496 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 497 "clarfx.f"
	    i__2 = j + (c_dim1 << 2);
#line 497 "clarfx.f"
	    i__3 = j + (c_dim1 << 2);
#line 497 "clarfx.f"
	    z__2.r = sum.r * t4.r - sum.i * t4.i, z__2.i = sum.r * t4.i + 
		    sum.i * t4.r;
#line 497 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 497 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 498 "clarfx.f"
/* L280: */
#line 498 "clarfx.f"
	}
#line 499 "clarfx.f"
	goto L410;
#line 500 "clarfx.f"
L290:

/*        Special code for 5 x 5 Householder */

#line 504 "clarfx.f"
	v1.r = v[1].r, v1.i = v[1].i;
#line 505 "clarfx.f"
	d_cnjg(&z__2, &v1);
#line 505 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 505 "clarfx.f"
	t1.r = z__1.r, t1.i = z__1.i;
#line 506 "clarfx.f"
	v2.r = v[2].r, v2.i = v[2].i;
#line 507 "clarfx.f"
	d_cnjg(&z__2, &v2);
#line 507 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 507 "clarfx.f"
	t2.r = z__1.r, t2.i = z__1.i;
#line 508 "clarfx.f"
	v3.r = v[3].r, v3.i = v[3].i;
#line 509 "clarfx.f"
	d_cnjg(&z__2, &v3);
#line 509 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 509 "clarfx.f"
	t3.r = z__1.r, t3.i = z__1.i;
#line 510 "clarfx.f"
	v4.r = v[4].r, v4.i = v[4].i;
#line 511 "clarfx.f"
	d_cnjg(&z__2, &v4);
#line 511 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 511 "clarfx.f"
	t4.r = z__1.r, t4.i = z__1.i;
#line 512 "clarfx.f"
	v5.r = v[5].r, v5.i = v[5].i;
#line 513 "clarfx.f"
	d_cnjg(&z__2, &v5);
#line 513 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 513 "clarfx.f"
	t5.r = z__1.r, t5.i = z__1.i;
#line 514 "clarfx.f"
	i__1 = *m;
#line 514 "clarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 515 "clarfx.f"
	    i__2 = j + c_dim1;
#line 515 "clarfx.f"
	    z__5.r = v1.r * c__[i__2].r - v1.i * c__[i__2].i, z__5.i = v1.r * 
		    c__[i__2].i + v1.i * c__[i__2].r;
#line 515 "clarfx.f"
	    i__3 = j + (c_dim1 << 1);
#line 515 "clarfx.f"
	    z__6.r = v2.r * c__[i__3].r - v2.i * c__[i__3].i, z__6.i = v2.r * 
		    c__[i__3].i + v2.i * c__[i__3].r;
#line 515 "clarfx.f"
	    z__4.r = z__5.r + z__6.r, z__4.i = z__5.i + z__6.i;
#line 515 "clarfx.f"
	    i__4 = j + c_dim1 * 3;
#line 515 "clarfx.f"
	    z__7.r = v3.r * c__[i__4].r - v3.i * c__[i__4].i, z__7.i = v3.r * 
		    c__[i__4].i + v3.i * c__[i__4].r;
#line 515 "clarfx.f"
	    z__3.r = z__4.r + z__7.r, z__3.i = z__4.i + z__7.i;
#line 515 "clarfx.f"
	    i__5 = j + (c_dim1 << 2);
#line 515 "clarfx.f"
	    z__8.r = v4.r * c__[i__5].r - v4.i * c__[i__5].i, z__8.i = v4.r * 
		    c__[i__5].i + v4.i * c__[i__5].r;
#line 515 "clarfx.f"
	    z__2.r = z__3.r + z__8.r, z__2.i = z__3.i + z__8.i;
#line 515 "clarfx.f"
	    i__6 = j + c_dim1 * 5;
#line 515 "clarfx.f"
	    z__9.r = v5.r * c__[i__6].r - v5.i * c__[i__6].i, z__9.i = v5.r * 
		    c__[i__6].i + v5.i * c__[i__6].r;
#line 515 "clarfx.f"
	    z__1.r = z__2.r + z__9.r, z__1.i = z__2.i + z__9.i;
#line 515 "clarfx.f"
	    sum.r = z__1.r, sum.i = z__1.i;
#line 517 "clarfx.f"
	    i__2 = j + c_dim1;
#line 517 "clarfx.f"
	    i__3 = j + c_dim1;
#line 517 "clarfx.f"
	    z__2.r = sum.r * t1.r - sum.i * t1.i, z__2.i = sum.r * t1.i + 
		    sum.i * t1.r;
#line 517 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 517 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 518 "clarfx.f"
	    i__2 = j + (c_dim1 << 1);
#line 518 "clarfx.f"
	    i__3 = j + (c_dim1 << 1);
#line 518 "clarfx.f"
	    z__2.r = sum.r * t2.r - sum.i * t2.i, z__2.i = sum.r * t2.i + 
		    sum.i * t2.r;
#line 518 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 518 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 519 "clarfx.f"
	    i__2 = j + c_dim1 * 3;
#line 519 "clarfx.f"
	    i__3 = j + c_dim1 * 3;
#line 519 "clarfx.f"
	    z__2.r = sum.r * t3.r - sum.i * t3.i, z__2.i = sum.r * t3.i + 
		    sum.i * t3.r;
#line 519 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 519 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 520 "clarfx.f"
	    i__2 = j + (c_dim1 << 2);
#line 520 "clarfx.f"
	    i__3 = j + (c_dim1 << 2);
#line 520 "clarfx.f"
	    z__2.r = sum.r * t4.r - sum.i * t4.i, z__2.i = sum.r * t4.i + 
		    sum.i * t4.r;
#line 520 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 520 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 521 "clarfx.f"
	    i__2 = j + c_dim1 * 5;
#line 521 "clarfx.f"
	    i__3 = j + c_dim1 * 5;
#line 521 "clarfx.f"
	    z__2.r = sum.r * t5.r - sum.i * t5.i, z__2.i = sum.r * t5.i + 
		    sum.i * t5.r;
#line 521 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 521 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 522 "clarfx.f"
/* L300: */
#line 522 "clarfx.f"
	}
#line 523 "clarfx.f"
	goto L410;
#line 524 "clarfx.f"
L310:

/*        Special code for 6 x 6 Householder */

#line 528 "clarfx.f"
	v1.r = v[1].r, v1.i = v[1].i;
#line 529 "clarfx.f"
	d_cnjg(&z__2, &v1);
#line 529 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 529 "clarfx.f"
	t1.r = z__1.r, t1.i = z__1.i;
#line 530 "clarfx.f"
	v2.r = v[2].r, v2.i = v[2].i;
#line 531 "clarfx.f"
	d_cnjg(&z__2, &v2);
#line 531 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 531 "clarfx.f"
	t2.r = z__1.r, t2.i = z__1.i;
#line 532 "clarfx.f"
	v3.r = v[3].r, v3.i = v[3].i;
#line 533 "clarfx.f"
	d_cnjg(&z__2, &v3);
#line 533 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 533 "clarfx.f"
	t3.r = z__1.r, t3.i = z__1.i;
#line 534 "clarfx.f"
	v4.r = v[4].r, v4.i = v[4].i;
#line 535 "clarfx.f"
	d_cnjg(&z__2, &v4);
#line 535 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 535 "clarfx.f"
	t4.r = z__1.r, t4.i = z__1.i;
#line 536 "clarfx.f"
	v5.r = v[5].r, v5.i = v[5].i;
#line 537 "clarfx.f"
	d_cnjg(&z__2, &v5);
#line 537 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 537 "clarfx.f"
	t5.r = z__1.r, t5.i = z__1.i;
#line 538 "clarfx.f"
	v6.r = v[6].r, v6.i = v[6].i;
#line 539 "clarfx.f"
	d_cnjg(&z__2, &v6);
#line 539 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 539 "clarfx.f"
	t6.r = z__1.r, t6.i = z__1.i;
#line 540 "clarfx.f"
	i__1 = *m;
#line 540 "clarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 541 "clarfx.f"
	    i__2 = j + c_dim1;
#line 541 "clarfx.f"
	    z__6.r = v1.r * c__[i__2].r - v1.i * c__[i__2].i, z__6.i = v1.r * 
		    c__[i__2].i + v1.i * c__[i__2].r;
#line 541 "clarfx.f"
	    i__3 = j + (c_dim1 << 1);
#line 541 "clarfx.f"
	    z__7.r = v2.r * c__[i__3].r - v2.i * c__[i__3].i, z__7.i = v2.r * 
		    c__[i__3].i + v2.i * c__[i__3].r;
#line 541 "clarfx.f"
	    z__5.r = z__6.r + z__7.r, z__5.i = z__6.i + z__7.i;
#line 541 "clarfx.f"
	    i__4 = j + c_dim1 * 3;
#line 541 "clarfx.f"
	    z__8.r = v3.r * c__[i__4].r - v3.i * c__[i__4].i, z__8.i = v3.r * 
		    c__[i__4].i + v3.i * c__[i__4].r;
#line 541 "clarfx.f"
	    z__4.r = z__5.r + z__8.r, z__4.i = z__5.i + z__8.i;
#line 541 "clarfx.f"
	    i__5 = j + (c_dim1 << 2);
#line 541 "clarfx.f"
	    z__9.r = v4.r * c__[i__5].r - v4.i * c__[i__5].i, z__9.i = v4.r * 
		    c__[i__5].i + v4.i * c__[i__5].r;
#line 541 "clarfx.f"
	    z__3.r = z__4.r + z__9.r, z__3.i = z__4.i + z__9.i;
#line 541 "clarfx.f"
	    i__6 = j + c_dim1 * 5;
#line 541 "clarfx.f"
	    z__10.r = v5.r * c__[i__6].r - v5.i * c__[i__6].i, z__10.i = v5.r 
		    * c__[i__6].i + v5.i * c__[i__6].r;
#line 541 "clarfx.f"
	    z__2.r = z__3.r + z__10.r, z__2.i = z__3.i + z__10.i;
#line 541 "clarfx.f"
	    i__7 = j + c_dim1 * 6;
#line 541 "clarfx.f"
	    z__11.r = v6.r * c__[i__7].r - v6.i * c__[i__7].i, z__11.i = v6.r 
		    * c__[i__7].i + v6.i * c__[i__7].r;
#line 541 "clarfx.f"
	    z__1.r = z__2.r + z__11.r, z__1.i = z__2.i + z__11.i;
#line 541 "clarfx.f"
	    sum.r = z__1.r, sum.i = z__1.i;
#line 543 "clarfx.f"
	    i__2 = j + c_dim1;
#line 543 "clarfx.f"
	    i__3 = j + c_dim1;
#line 543 "clarfx.f"
	    z__2.r = sum.r * t1.r - sum.i * t1.i, z__2.i = sum.r * t1.i + 
		    sum.i * t1.r;
#line 543 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 543 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 544 "clarfx.f"
	    i__2 = j + (c_dim1 << 1);
#line 544 "clarfx.f"
	    i__3 = j + (c_dim1 << 1);
#line 544 "clarfx.f"
	    z__2.r = sum.r * t2.r - sum.i * t2.i, z__2.i = sum.r * t2.i + 
		    sum.i * t2.r;
#line 544 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 544 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 545 "clarfx.f"
	    i__2 = j + c_dim1 * 3;
#line 545 "clarfx.f"
	    i__3 = j + c_dim1 * 3;
#line 545 "clarfx.f"
	    z__2.r = sum.r * t3.r - sum.i * t3.i, z__2.i = sum.r * t3.i + 
		    sum.i * t3.r;
#line 545 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 545 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 546 "clarfx.f"
	    i__2 = j + (c_dim1 << 2);
#line 546 "clarfx.f"
	    i__3 = j + (c_dim1 << 2);
#line 546 "clarfx.f"
	    z__2.r = sum.r * t4.r - sum.i * t4.i, z__2.i = sum.r * t4.i + 
		    sum.i * t4.r;
#line 546 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 546 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 547 "clarfx.f"
	    i__2 = j + c_dim1 * 5;
#line 547 "clarfx.f"
	    i__3 = j + c_dim1 * 5;
#line 547 "clarfx.f"
	    z__2.r = sum.r * t5.r - sum.i * t5.i, z__2.i = sum.r * t5.i + 
		    sum.i * t5.r;
#line 547 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 547 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 548 "clarfx.f"
	    i__2 = j + c_dim1 * 6;
#line 548 "clarfx.f"
	    i__3 = j + c_dim1 * 6;
#line 548 "clarfx.f"
	    z__2.r = sum.r * t6.r - sum.i * t6.i, z__2.i = sum.r * t6.i + 
		    sum.i * t6.r;
#line 548 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 548 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 549 "clarfx.f"
/* L320: */
#line 549 "clarfx.f"
	}
#line 550 "clarfx.f"
	goto L410;
#line 551 "clarfx.f"
L330:

/*        Special code for 7 x 7 Householder */

#line 555 "clarfx.f"
	v1.r = v[1].r, v1.i = v[1].i;
#line 556 "clarfx.f"
	d_cnjg(&z__2, &v1);
#line 556 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 556 "clarfx.f"
	t1.r = z__1.r, t1.i = z__1.i;
#line 557 "clarfx.f"
	v2.r = v[2].r, v2.i = v[2].i;
#line 558 "clarfx.f"
	d_cnjg(&z__2, &v2);
#line 558 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 558 "clarfx.f"
	t2.r = z__1.r, t2.i = z__1.i;
#line 559 "clarfx.f"
	v3.r = v[3].r, v3.i = v[3].i;
#line 560 "clarfx.f"
	d_cnjg(&z__2, &v3);
#line 560 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 560 "clarfx.f"
	t3.r = z__1.r, t3.i = z__1.i;
#line 561 "clarfx.f"
	v4.r = v[4].r, v4.i = v[4].i;
#line 562 "clarfx.f"
	d_cnjg(&z__2, &v4);
#line 562 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 562 "clarfx.f"
	t4.r = z__1.r, t4.i = z__1.i;
#line 563 "clarfx.f"
	v5.r = v[5].r, v5.i = v[5].i;
#line 564 "clarfx.f"
	d_cnjg(&z__2, &v5);
#line 564 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 564 "clarfx.f"
	t5.r = z__1.r, t5.i = z__1.i;
#line 565 "clarfx.f"
	v6.r = v[6].r, v6.i = v[6].i;
#line 566 "clarfx.f"
	d_cnjg(&z__2, &v6);
#line 566 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 566 "clarfx.f"
	t6.r = z__1.r, t6.i = z__1.i;
#line 567 "clarfx.f"
	v7.r = v[7].r, v7.i = v[7].i;
#line 568 "clarfx.f"
	d_cnjg(&z__2, &v7);
#line 568 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 568 "clarfx.f"
	t7.r = z__1.r, t7.i = z__1.i;
#line 569 "clarfx.f"
	i__1 = *m;
#line 569 "clarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 570 "clarfx.f"
	    i__2 = j + c_dim1;
#line 570 "clarfx.f"
	    z__7.r = v1.r * c__[i__2].r - v1.i * c__[i__2].i, z__7.i = v1.r * 
		    c__[i__2].i + v1.i * c__[i__2].r;
#line 570 "clarfx.f"
	    i__3 = j + (c_dim1 << 1);
#line 570 "clarfx.f"
	    z__8.r = v2.r * c__[i__3].r - v2.i * c__[i__3].i, z__8.i = v2.r * 
		    c__[i__3].i + v2.i * c__[i__3].r;
#line 570 "clarfx.f"
	    z__6.r = z__7.r + z__8.r, z__6.i = z__7.i + z__8.i;
#line 570 "clarfx.f"
	    i__4 = j + c_dim1 * 3;
#line 570 "clarfx.f"
	    z__9.r = v3.r * c__[i__4].r - v3.i * c__[i__4].i, z__9.i = v3.r * 
		    c__[i__4].i + v3.i * c__[i__4].r;
#line 570 "clarfx.f"
	    z__5.r = z__6.r + z__9.r, z__5.i = z__6.i + z__9.i;
#line 570 "clarfx.f"
	    i__5 = j + (c_dim1 << 2);
#line 570 "clarfx.f"
	    z__10.r = v4.r * c__[i__5].r - v4.i * c__[i__5].i, z__10.i = v4.r 
		    * c__[i__5].i + v4.i * c__[i__5].r;
#line 570 "clarfx.f"
	    z__4.r = z__5.r + z__10.r, z__4.i = z__5.i + z__10.i;
#line 570 "clarfx.f"
	    i__6 = j + c_dim1 * 5;
#line 570 "clarfx.f"
	    z__11.r = v5.r * c__[i__6].r - v5.i * c__[i__6].i, z__11.i = v5.r 
		    * c__[i__6].i + v5.i * c__[i__6].r;
#line 570 "clarfx.f"
	    z__3.r = z__4.r + z__11.r, z__3.i = z__4.i + z__11.i;
#line 570 "clarfx.f"
	    i__7 = j + c_dim1 * 6;
#line 570 "clarfx.f"
	    z__12.r = v6.r * c__[i__7].r - v6.i * c__[i__7].i, z__12.i = v6.r 
		    * c__[i__7].i + v6.i * c__[i__7].r;
#line 570 "clarfx.f"
	    z__2.r = z__3.r + z__12.r, z__2.i = z__3.i + z__12.i;
#line 570 "clarfx.f"
	    i__8 = j + c_dim1 * 7;
#line 570 "clarfx.f"
	    z__13.r = v7.r * c__[i__8].r - v7.i * c__[i__8].i, z__13.i = v7.r 
		    * c__[i__8].i + v7.i * c__[i__8].r;
#line 570 "clarfx.f"
	    z__1.r = z__2.r + z__13.r, z__1.i = z__2.i + z__13.i;
#line 570 "clarfx.f"
	    sum.r = z__1.r, sum.i = z__1.i;
#line 573 "clarfx.f"
	    i__2 = j + c_dim1;
#line 573 "clarfx.f"
	    i__3 = j + c_dim1;
#line 573 "clarfx.f"
	    z__2.r = sum.r * t1.r - sum.i * t1.i, z__2.i = sum.r * t1.i + 
		    sum.i * t1.r;
#line 573 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 573 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 574 "clarfx.f"
	    i__2 = j + (c_dim1 << 1);
#line 574 "clarfx.f"
	    i__3 = j + (c_dim1 << 1);
#line 574 "clarfx.f"
	    z__2.r = sum.r * t2.r - sum.i * t2.i, z__2.i = sum.r * t2.i + 
		    sum.i * t2.r;
#line 574 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 574 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 575 "clarfx.f"
	    i__2 = j + c_dim1 * 3;
#line 575 "clarfx.f"
	    i__3 = j + c_dim1 * 3;
#line 575 "clarfx.f"
	    z__2.r = sum.r * t3.r - sum.i * t3.i, z__2.i = sum.r * t3.i + 
		    sum.i * t3.r;
#line 575 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 575 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 576 "clarfx.f"
	    i__2 = j + (c_dim1 << 2);
#line 576 "clarfx.f"
	    i__3 = j + (c_dim1 << 2);
#line 576 "clarfx.f"
	    z__2.r = sum.r * t4.r - sum.i * t4.i, z__2.i = sum.r * t4.i + 
		    sum.i * t4.r;
#line 576 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 576 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 577 "clarfx.f"
	    i__2 = j + c_dim1 * 5;
#line 577 "clarfx.f"
	    i__3 = j + c_dim1 * 5;
#line 577 "clarfx.f"
	    z__2.r = sum.r * t5.r - sum.i * t5.i, z__2.i = sum.r * t5.i + 
		    sum.i * t5.r;
#line 577 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 577 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 578 "clarfx.f"
	    i__2 = j + c_dim1 * 6;
#line 578 "clarfx.f"
	    i__3 = j + c_dim1 * 6;
#line 578 "clarfx.f"
	    z__2.r = sum.r * t6.r - sum.i * t6.i, z__2.i = sum.r * t6.i + 
		    sum.i * t6.r;
#line 578 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 578 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 579 "clarfx.f"
	    i__2 = j + c_dim1 * 7;
#line 579 "clarfx.f"
	    i__3 = j + c_dim1 * 7;
#line 579 "clarfx.f"
	    z__2.r = sum.r * t7.r - sum.i * t7.i, z__2.i = sum.r * t7.i + 
		    sum.i * t7.r;
#line 579 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 579 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 580 "clarfx.f"
/* L340: */
#line 580 "clarfx.f"
	}
#line 581 "clarfx.f"
	goto L410;
#line 582 "clarfx.f"
L350:

/*        Special code for 8 x 8 Householder */

#line 586 "clarfx.f"
	v1.r = v[1].r, v1.i = v[1].i;
#line 587 "clarfx.f"
	d_cnjg(&z__2, &v1);
#line 587 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 587 "clarfx.f"
	t1.r = z__1.r, t1.i = z__1.i;
#line 588 "clarfx.f"
	v2.r = v[2].r, v2.i = v[2].i;
#line 589 "clarfx.f"
	d_cnjg(&z__2, &v2);
#line 589 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 589 "clarfx.f"
	t2.r = z__1.r, t2.i = z__1.i;
#line 590 "clarfx.f"
	v3.r = v[3].r, v3.i = v[3].i;
#line 591 "clarfx.f"
	d_cnjg(&z__2, &v3);
#line 591 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 591 "clarfx.f"
	t3.r = z__1.r, t3.i = z__1.i;
#line 592 "clarfx.f"
	v4.r = v[4].r, v4.i = v[4].i;
#line 593 "clarfx.f"
	d_cnjg(&z__2, &v4);
#line 593 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 593 "clarfx.f"
	t4.r = z__1.r, t4.i = z__1.i;
#line 594 "clarfx.f"
	v5.r = v[5].r, v5.i = v[5].i;
#line 595 "clarfx.f"
	d_cnjg(&z__2, &v5);
#line 595 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 595 "clarfx.f"
	t5.r = z__1.r, t5.i = z__1.i;
#line 596 "clarfx.f"
	v6.r = v[6].r, v6.i = v[6].i;
#line 597 "clarfx.f"
	d_cnjg(&z__2, &v6);
#line 597 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 597 "clarfx.f"
	t6.r = z__1.r, t6.i = z__1.i;
#line 598 "clarfx.f"
	v7.r = v[7].r, v7.i = v[7].i;
#line 599 "clarfx.f"
	d_cnjg(&z__2, &v7);
#line 599 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 599 "clarfx.f"
	t7.r = z__1.r, t7.i = z__1.i;
#line 600 "clarfx.f"
	v8.r = v[8].r, v8.i = v[8].i;
#line 601 "clarfx.f"
	d_cnjg(&z__2, &v8);
#line 601 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 601 "clarfx.f"
	t8.r = z__1.r, t8.i = z__1.i;
#line 602 "clarfx.f"
	i__1 = *m;
#line 602 "clarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 603 "clarfx.f"
	    i__2 = j + c_dim1;
#line 603 "clarfx.f"
	    z__8.r = v1.r * c__[i__2].r - v1.i * c__[i__2].i, z__8.i = v1.r * 
		    c__[i__2].i + v1.i * c__[i__2].r;
#line 603 "clarfx.f"
	    i__3 = j + (c_dim1 << 1);
#line 603 "clarfx.f"
	    z__9.r = v2.r * c__[i__3].r - v2.i * c__[i__3].i, z__9.i = v2.r * 
		    c__[i__3].i + v2.i * c__[i__3].r;
#line 603 "clarfx.f"
	    z__7.r = z__8.r + z__9.r, z__7.i = z__8.i + z__9.i;
#line 603 "clarfx.f"
	    i__4 = j + c_dim1 * 3;
#line 603 "clarfx.f"
	    z__10.r = v3.r * c__[i__4].r - v3.i * c__[i__4].i, z__10.i = v3.r 
		    * c__[i__4].i + v3.i * c__[i__4].r;
#line 603 "clarfx.f"
	    z__6.r = z__7.r + z__10.r, z__6.i = z__7.i + z__10.i;
#line 603 "clarfx.f"
	    i__5 = j + (c_dim1 << 2);
#line 603 "clarfx.f"
	    z__11.r = v4.r * c__[i__5].r - v4.i * c__[i__5].i, z__11.i = v4.r 
		    * c__[i__5].i + v4.i * c__[i__5].r;
#line 603 "clarfx.f"
	    z__5.r = z__6.r + z__11.r, z__5.i = z__6.i + z__11.i;
#line 603 "clarfx.f"
	    i__6 = j + c_dim1 * 5;
#line 603 "clarfx.f"
	    z__12.r = v5.r * c__[i__6].r - v5.i * c__[i__6].i, z__12.i = v5.r 
		    * c__[i__6].i + v5.i * c__[i__6].r;
#line 603 "clarfx.f"
	    z__4.r = z__5.r + z__12.r, z__4.i = z__5.i + z__12.i;
#line 603 "clarfx.f"
	    i__7 = j + c_dim1 * 6;
#line 603 "clarfx.f"
	    z__13.r = v6.r * c__[i__7].r - v6.i * c__[i__7].i, z__13.i = v6.r 
		    * c__[i__7].i + v6.i * c__[i__7].r;
#line 603 "clarfx.f"
	    z__3.r = z__4.r + z__13.r, z__3.i = z__4.i + z__13.i;
#line 603 "clarfx.f"
	    i__8 = j + c_dim1 * 7;
#line 603 "clarfx.f"
	    z__14.r = v7.r * c__[i__8].r - v7.i * c__[i__8].i, z__14.i = v7.r 
		    * c__[i__8].i + v7.i * c__[i__8].r;
#line 603 "clarfx.f"
	    z__2.r = z__3.r + z__14.r, z__2.i = z__3.i + z__14.i;
#line 603 "clarfx.f"
	    i__9 = j + (c_dim1 << 3);
#line 603 "clarfx.f"
	    z__15.r = v8.r * c__[i__9].r - v8.i * c__[i__9].i, z__15.i = v8.r 
		    * c__[i__9].i + v8.i * c__[i__9].r;
#line 603 "clarfx.f"
	    z__1.r = z__2.r + z__15.r, z__1.i = z__2.i + z__15.i;
#line 603 "clarfx.f"
	    sum.r = z__1.r, sum.i = z__1.i;
#line 606 "clarfx.f"
	    i__2 = j + c_dim1;
#line 606 "clarfx.f"
	    i__3 = j + c_dim1;
#line 606 "clarfx.f"
	    z__2.r = sum.r * t1.r - sum.i * t1.i, z__2.i = sum.r * t1.i + 
		    sum.i * t1.r;
#line 606 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 606 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 607 "clarfx.f"
	    i__2 = j + (c_dim1 << 1);
#line 607 "clarfx.f"
	    i__3 = j + (c_dim1 << 1);
#line 607 "clarfx.f"
	    z__2.r = sum.r * t2.r - sum.i * t2.i, z__2.i = sum.r * t2.i + 
		    sum.i * t2.r;
#line 607 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 607 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 608 "clarfx.f"
	    i__2 = j + c_dim1 * 3;
#line 608 "clarfx.f"
	    i__3 = j + c_dim1 * 3;
#line 608 "clarfx.f"
	    z__2.r = sum.r * t3.r - sum.i * t3.i, z__2.i = sum.r * t3.i + 
		    sum.i * t3.r;
#line 608 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 608 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 609 "clarfx.f"
	    i__2 = j + (c_dim1 << 2);
#line 609 "clarfx.f"
	    i__3 = j + (c_dim1 << 2);
#line 609 "clarfx.f"
	    z__2.r = sum.r * t4.r - sum.i * t4.i, z__2.i = sum.r * t4.i + 
		    sum.i * t4.r;
#line 609 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 609 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 610 "clarfx.f"
	    i__2 = j + c_dim1 * 5;
#line 610 "clarfx.f"
	    i__3 = j + c_dim1 * 5;
#line 610 "clarfx.f"
	    z__2.r = sum.r * t5.r - sum.i * t5.i, z__2.i = sum.r * t5.i + 
		    sum.i * t5.r;
#line 610 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 610 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 611 "clarfx.f"
	    i__2 = j + c_dim1 * 6;
#line 611 "clarfx.f"
	    i__3 = j + c_dim1 * 6;
#line 611 "clarfx.f"
	    z__2.r = sum.r * t6.r - sum.i * t6.i, z__2.i = sum.r * t6.i + 
		    sum.i * t6.r;
#line 611 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 611 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 612 "clarfx.f"
	    i__2 = j + c_dim1 * 7;
#line 612 "clarfx.f"
	    i__3 = j + c_dim1 * 7;
#line 612 "clarfx.f"
	    z__2.r = sum.r * t7.r - sum.i * t7.i, z__2.i = sum.r * t7.i + 
		    sum.i * t7.r;
#line 612 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 612 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 613 "clarfx.f"
	    i__2 = j + (c_dim1 << 3);
#line 613 "clarfx.f"
	    i__3 = j + (c_dim1 << 3);
#line 613 "clarfx.f"
	    z__2.r = sum.r * t8.r - sum.i * t8.i, z__2.i = sum.r * t8.i + 
		    sum.i * t8.r;
#line 613 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 613 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 614 "clarfx.f"
/* L360: */
#line 614 "clarfx.f"
	}
#line 615 "clarfx.f"
	goto L410;
#line 616 "clarfx.f"
L370:

/*        Special code for 9 x 9 Householder */

#line 620 "clarfx.f"
	v1.r = v[1].r, v1.i = v[1].i;
#line 621 "clarfx.f"
	d_cnjg(&z__2, &v1);
#line 621 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 621 "clarfx.f"
	t1.r = z__1.r, t1.i = z__1.i;
#line 622 "clarfx.f"
	v2.r = v[2].r, v2.i = v[2].i;
#line 623 "clarfx.f"
	d_cnjg(&z__2, &v2);
#line 623 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 623 "clarfx.f"
	t2.r = z__1.r, t2.i = z__1.i;
#line 624 "clarfx.f"
	v3.r = v[3].r, v3.i = v[3].i;
#line 625 "clarfx.f"
	d_cnjg(&z__2, &v3);
#line 625 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 625 "clarfx.f"
	t3.r = z__1.r, t3.i = z__1.i;
#line 626 "clarfx.f"
	v4.r = v[4].r, v4.i = v[4].i;
#line 627 "clarfx.f"
	d_cnjg(&z__2, &v4);
#line 627 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 627 "clarfx.f"
	t4.r = z__1.r, t4.i = z__1.i;
#line 628 "clarfx.f"
	v5.r = v[5].r, v5.i = v[5].i;
#line 629 "clarfx.f"
	d_cnjg(&z__2, &v5);
#line 629 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 629 "clarfx.f"
	t5.r = z__1.r, t5.i = z__1.i;
#line 630 "clarfx.f"
	v6.r = v[6].r, v6.i = v[6].i;
#line 631 "clarfx.f"
	d_cnjg(&z__2, &v6);
#line 631 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 631 "clarfx.f"
	t6.r = z__1.r, t6.i = z__1.i;
#line 632 "clarfx.f"
	v7.r = v[7].r, v7.i = v[7].i;
#line 633 "clarfx.f"
	d_cnjg(&z__2, &v7);
#line 633 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 633 "clarfx.f"
	t7.r = z__1.r, t7.i = z__1.i;
#line 634 "clarfx.f"
	v8.r = v[8].r, v8.i = v[8].i;
#line 635 "clarfx.f"
	d_cnjg(&z__2, &v8);
#line 635 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 635 "clarfx.f"
	t8.r = z__1.r, t8.i = z__1.i;
#line 636 "clarfx.f"
	v9.r = v[9].r, v9.i = v[9].i;
#line 637 "clarfx.f"
	d_cnjg(&z__2, &v9);
#line 637 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 637 "clarfx.f"
	t9.r = z__1.r, t9.i = z__1.i;
#line 638 "clarfx.f"
	i__1 = *m;
#line 638 "clarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 639 "clarfx.f"
	    i__2 = j + c_dim1;
#line 639 "clarfx.f"
	    z__9.r = v1.r * c__[i__2].r - v1.i * c__[i__2].i, z__9.i = v1.r * 
		    c__[i__2].i + v1.i * c__[i__2].r;
#line 639 "clarfx.f"
	    i__3 = j + (c_dim1 << 1);
#line 639 "clarfx.f"
	    z__10.r = v2.r * c__[i__3].r - v2.i * c__[i__3].i, z__10.i = v2.r 
		    * c__[i__3].i + v2.i * c__[i__3].r;
#line 639 "clarfx.f"
	    z__8.r = z__9.r + z__10.r, z__8.i = z__9.i + z__10.i;
#line 639 "clarfx.f"
	    i__4 = j + c_dim1 * 3;
#line 639 "clarfx.f"
	    z__11.r = v3.r * c__[i__4].r - v3.i * c__[i__4].i, z__11.i = v3.r 
		    * c__[i__4].i + v3.i * c__[i__4].r;
#line 639 "clarfx.f"
	    z__7.r = z__8.r + z__11.r, z__7.i = z__8.i + z__11.i;
#line 639 "clarfx.f"
	    i__5 = j + (c_dim1 << 2);
#line 639 "clarfx.f"
	    z__12.r = v4.r * c__[i__5].r - v4.i * c__[i__5].i, z__12.i = v4.r 
		    * c__[i__5].i + v4.i * c__[i__5].r;
#line 639 "clarfx.f"
	    z__6.r = z__7.r + z__12.r, z__6.i = z__7.i + z__12.i;
#line 639 "clarfx.f"
	    i__6 = j + c_dim1 * 5;
#line 639 "clarfx.f"
	    z__13.r = v5.r * c__[i__6].r - v5.i * c__[i__6].i, z__13.i = v5.r 
		    * c__[i__6].i + v5.i * c__[i__6].r;
#line 639 "clarfx.f"
	    z__5.r = z__6.r + z__13.r, z__5.i = z__6.i + z__13.i;
#line 639 "clarfx.f"
	    i__7 = j + c_dim1 * 6;
#line 639 "clarfx.f"
	    z__14.r = v6.r * c__[i__7].r - v6.i * c__[i__7].i, z__14.i = v6.r 
		    * c__[i__7].i + v6.i * c__[i__7].r;
#line 639 "clarfx.f"
	    z__4.r = z__5.r + z__14.r, z__4.i = z__5.i + z__14.i;
#line 639 "clarfx.f"
	    i__8 = j + c_dim1 * 7;
#line 639 "clarfx.f"
	    z__15.r = v7.r * c__[i__8].r - v7.i * c__[i__8].i, z__15.i = v7.r 
		    * c__[i__8].i + v7.i * c__[i__8].r;
#line 639 "clarfx.f"
	    z__3.r = z__4.r + z__15.r, z__3.i = z__4.i + z__15.i;
#line 639 "clarfx.f"
	    i__9 = j + (c_dim1 << 3);
#line 639 "clarfx.f"
	    z__16.r = v8.r * c__[i__9].r - v8.i * c__[i__9].i, z__16.i = v8.r 
		    * c__[i__9].i + v8.i * c__[i__9].r;
#line 639 "clarfx.f"
	    z__2.r = z__3.r + z__16.r, z__2.i = z__3.i + z__16.i;
#line 639 "clarfx.f"
	    i__10 = j + c_dim1 * 9;
#line 639 "clarfx.f"
	    z__17.r = v9.r * c__[i__10].r - v9.i * c__[i__10].i, z__17.i = 
		    v9.r * c__[i__10].i + v9.i * c__[i__10].r;
#line 639 "clarfx.f"
	    z__1.r = z__2.r + z__17.r, z__1.i = z__2.i + z__17.i;
#line 639 "clarfx.f"
	    sum.r = z__1.r, sum.i = z__1.i;
#line 642 "clarfx.f"
	    i__2 = j + c_dim1;
#line 642 "clarfx.f"
	    i__3 = j + c_dim1;
#line 642 "clarfx.f"
	    z__2.r = sum.r * t1.r - sum.i * t1.i, z__2.i = sum.r * t1.i + 
		    sum.i * t1.r;
#line 642 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 642 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 643 "clarfx.f"
	    i__2 = j + (c_dim1 << 1);
#line 643 "clarfx.f"
	    i__3 = j + (c_dim1 << 1);
#line 643 "clarfx.f"
	    z__2.r = sum.r * t2.r - sum.i * t2.i, z__2.i = sum.r * t2.i + 
		    sum.i * t2.r;
#line 643 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 643 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 644 "clarfx.f"
	    i__2 = j + c_dim1 * 3;
#line 644 "clarfx.f"
	    i__3 = j + c_dim1 * 3;
#line 644 "clarfx.f"
	    z__2.r = sum.r * t3.r - sum.i * t3.i, z__2.i = sum.r * t3.i + 
		    sum.i * t3.r;
#line 644 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 644 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 645 "clarfx.f"
	    i__2 = j + (c_dim1 << 2);
#line 645 "clarfx.f"
	    i__3 = j + (c_dim1 << 2);
#line 645 "clarfx.f"
	    z__2.r = sum.r * t4.r - sum.i * t4.i, z__2.i = sum.r * t4.i + 
		    sum.i * t4.r;
#line 645 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 645 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 646 "clarfx.f"
	    i__2 = j + c_dim1 * 5;
#line 646 "clarfx.f"
	    i__3 = j + c_dim1 * 5;
#line 646 "clarfx.f"
	    z__2.r = sum.r * t5.r - sum.i * t5.i, z__2.i = sum.r * t5.i + 
		    sum.i * t5.r;
#line 646 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 646 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 647 "clarfx.f"
	    i__2 = j + c_dim1 * 6;
#line 647 "clarfx.f"
	    i__3 = j + c_dim1 * 6;
#line 647 "clarfx.f"
	    z__2.r = sum.r * t6.r - sum.i * t6.i, z__2.i = sum.r * t6.i + 
		    sum.i * t6.r;
#line 647 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 647 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 648 "clarfx.f"
	    i__2 = j + c_dim1 * 7;
#line 648 "clarfx.f"
	    i__3 = j + c_dim1 * 7;
#line 648 "clarfx.f"
	    z__2.r = sum.r * t7.r - sum.i * t7.i, z__2.i = sum.r * t7.i + 
		    sum.i * t7.r;
#line 648 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 648 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 649 "clarfx.f"
	    i__2 = j + (c_dim1 << 3);
#line 649 "clarfx.f"
	    i__3 = j + (c_dim1 << 3);
#line 649 "clarfx.f"
	    z__2.r = sum.r * t8.r - sum.i * t8.i, z__2.i = sum.r * t8.i + 
		    sum.i * t8.r;
#line 649 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 649 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 650 "clarfx.f"
	    i__2 = j + c_dim1 * 9;
#line 650 "clarfx.f"
	    i__3 = j + c_dim1 * 9;
#line 650 "clarfx.f"
	    z__2.r = sum.r * t9.r - sum.i * t9.i, z__2.i = sum.r * t9.i + 
		    sum.i * t9.r;
#line 650 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 650 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 651 "clarfx.f"
/* L380: */
#line 651 "clarfx.f"
	}
#line 652 "clarfx.f"
	goto L410;
#line 653 "clarfx.f"
L390:

/*        Special code for 10 x 10 Householder */

#line 657 "clarfx.f"
	v1.r = v[1].r, v1.i = v[1].i;
#line 658 "clarfx.f"
	d_cnjg(&z__2, &v1);
#line 658 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 658 "clarfx.f"
	t1.r = z__1.r, t1.i = z__1.i;
#line 659 "clarfx.f"
	v2.r = v[2].r, v2.i = v[2].i;
#line 660 "clarfx.f"
	d_cnjg(&z__2, &v2);
#line 660 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 660 "clarfx.f"
	t2.r = z__1.r, t2.i = z__1.i;
#line 661 "clarfx.f"
	v3.r = v[3].r, v3.i = v[3].i;
#line 662 "clarfx.f"
	d_cnjg(&z__2, &v3);
#line 662 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 662 "clarfx.f"
	t3.r = z__1.r, t3.i = z__1.i;
#line 663 "clarfx.f"
	v4.r = v[4].r, v4.i = v[4].i;
#line 664 "clarfx.f"
	d_cnjg(&z__2, &v4);
#line 664 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 664 "clarfx.f"
	t4.r = z__1.r, t4.i = z__1.i;
#line 665 "clarfx.f"
	v5.r = v[5].r, v5.i = v[5].i;
#line 666 "clarfx.f"
	d_cnjg(&z__2, &v5);
#line 666 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 666 "clarfx.f"
	t5.r = z__1.r, t5.i = z__1.i;
#line 667 "clarfx.f"
	v6.r = v[6].r, v6.i = v[6].i;
#line 668 "clarfx.f"
	d_cnjg(&z__2, &v6);
#line 668 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 668 "clarfx.f"
	t6.r = z__1.r, t6.i = z__1.i;
#line 669 "clarfx.f"
	v7.r = v[7].r, v7.i = v[7].i;
#line 670 "clarfx.f"
	d_cnjg(&z__2, &v7);
#line 670 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 670 "clarfx.f"
	t7.r = z__1.r, t7.i = z__1.i;
#line 671 "clarfx.f"
	v8.r = v[8].r, v8.i = v[8].i;
#line 672 "clarfx.f"
	d_cnjg(&z__2, &v8);
#line 672 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 672 "clarfx.f"
	t8.r = z__1.r, t8.i = z__1.i;
#line 673 "clarfx.f"
	v9.r = v[9].r, v9.i = v[9].i;
#line 674 "clarfx.f"
	d_cnjg(&z__2, &v9);
#line 674 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 674 "clarfx.f"
	t9.r = z__1.r, t9.i = z__1.i;
#line 675 "clarfx.f"
	v10.r = v[10].r, v10.i = v[10].i;
#line 676 "clarfx.f"
	d_cnjg(&z__2, &v10);
#line 676 "clarfx.f"
	z__1.r = tau->r * z__2.r - tau->i * z__2.i, z__1.i = tau->r * z__2.i 
		+ tau->i * z__2.r;
#line 676 "clarfx.f"
	t10.r = z__1.r, t10.i = z__1.i;
#line 677 "clarfx.f"
	i__1 = *m;
#line 677 "clarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 678 "clarfx.f"
	    i__2 = j + c_dim1;
#line 678 "clarfx.f"
	    z__10.r = v1.r * c__[i__2].r - v1.i * c__[i__2].i, z__10.i = v1.r 
		    * c__[i__2].i + v1.i * c__[i__2].r;
#line 678 "clarfx.f"
	    i__3 = j + (c_dim1 << 1);
#line 678 "clarfx.f"
	    z__11.r = v2.r * c__[i__3].r - v2.i * c__[i__3].i, z__11.i = v2.r 
		    * c__[i__3].i + v2.i * c__[i__3].r;
#line 678 "clarfx.f"
	    z__9.r = z__10.r + z__11.r, z__9.i = z__10.i + z__11.i;
#line 678 "clarfx.f"
	    i__4 = j + c_dim1 * 3;
#line 678 "clarfx.f"
	    z__12.r = v3.r * c__[i__4].r - v3.i * c__[i__4].i, z__12.i = v3.r 
		    * c__[i__4].i + v3.i * c__[i__4].r;
#line 678 "clarfx.f"
	    z__8.r = z__9.r + z__12.r, z__8.i = z__9.i + z__12.i;
#line 678 "clarfx.f"
	    i__5 = j + (c_dim1 << 2);
#line 678 "clarfx.f"
	    z__13.r = v4.r * c__[i__5].r - v4.i * c__[i__5].i, z__13.i = v4.r 
		    * c__[i__5].i + v4.i * c__[i__5].r;
#line 678 "clarfx.f"
	    z__7.r = z__8.r + z__13.r, z__7.i = z__8.i + z__13.i;
#line 678 "clarfx.f"
	    i__6 = j + c_dim1 * 5;
#line 678 "clarfx.f"
	    z__14.r = v5.r * c__[i__6].r - v5.i * c__[i__6].i, z__14.i = v5.r 
		    * c__[i__6].i + v5.i * c__[i__6].r;
#line 678 "clarfx.f"
	    z__6.r = z__7.r + z__14.r, z__6.i = z__7.i + z__14.i;
#line 678 "clarfx.f"
	    i__7 = j + c_dim1 * 6;
#line 678 "clarfx.f"
	    z__15.r = v6.r * c__[i__7].r - v6.i * c__[i__7].i, z__15.i = v6.r 
		    * c__[i__7].i + v6.i * c__[i__7].r;
#line 678 "clarfx.f"
	    z__5.r = z__6.r + z__15.r, z__5.i = z__6.i + z__15.i;
#line 678 "clarfx.f"
	    i__8 = j + c_dim1 * 7;
#line 678 "clarfx.f"
	    z__16.r = v7.r * c__[i__8].r - v7.i * c__[i__8].i, z__16.i = v7.r 
		    * c__[i__8].i + v7.i * c__[i__8].r;
#line 678 "clarfx.f"
	    z__4.r = z__5.r + z__16.r, z__4.i = z__5.i + z__16.i;
#line 678 "clarfx.f"
	    i__9 = j + (c_dim1 << 3);
#line 678 "clarfx.f"
	    z__17.r = v8.r * c__[i__9].r - v8.i * c__[i__9].i, z__17.i = v8.r 
		    * c__[i__9].i + v8.i * c__[i__9].r;
#line 678 "clarfx.f"
	    z__3.r = z__4.r + z__17.r, z__3.i = z__4.i + z__17.i;
#line 678 "clarfx.f"
	    i__10 = j + c_dim1 * 9;
#line 678 "clarfx.f"
	    z__18.r = v9.r * c__[i__10].r - v9.i * c__[i__10].i, z__18.i = 
		    v9.r * c__[i__10].i + v9.i * c__[i__10].r;
#line 678 "clarfx.f"
	    z__2.r = z__3.r + z__18.r, z__2.i = z__3.i + z__18.i;
#line 678 "clarfx.f"
	    i__11 = j + c_dim1 * 10;
#line 678 "clarfx.f"
	    z__19.r = v10.r * c__[i__11].r - v10.i * c__[i__11].i, z__19.i = 
		    v10.r * c__[i__11].i + v10.i * c__[i__11].r;
#line 678 "clarfx.f"
	    z__1.r = z__2.r + z__19.r, z__1.i = z__2.i + z__19.i;
#line 678 "clarfx.f"
	    sum.r = z__1.r, sum.i = z__1.i;
#line 682 "clarfx.f"
	    i__2 = j + c_dim1;
#line 682 "clarfx.f"
	    i__3 = j + c_dim1;
#line 682 "clarfx.f"
	    z__2.r = sum.r * t1.r - sum.i * t1.i, z__2.i = sum.r * t1.i + 
		    sum.i * t1.r;
#line 682 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 682 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 683 "clarfx.f"
	    i__2 = j + (c_dim1 << 1);
#line 683 "clarfx.f"
	    i__3 = j + (c_dim1 << 1);
#line 683 "clarfx.f"
	    z__2.r = sum.r * t2.r - sum.i * t2.i, z__2.i = sum.r * t2.i + 
		    sum.i * t2.r;
#line 683 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 683 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 684 "clarfx.f"
	    i__2 = j + c_dim1 * 3;
#line 684 "clarfx.f"
	    i__3 = j + c_dim1 * 3;
#line 684 "clarfx.f"
	    z__2.r = sum.r * t3.r - sum.i * t3.i, z__2.i = sum.r * t3.i + 
		    sum.i * t3.r;
#line 684 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 684 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 685 "clarfx.f"
	    i__2 = j + (c_dim1 << 2);
#line 685 "clarfx.f"
	    i__3 = j + (c_dim1 << 2);
#line 685 "clarfx.f"
	    z__2.r = sum.r * t4.r - sum.i * t4.i, z__2.i = sum.r * t4.i + 
		    sum.i * t4.r;
#line 685 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 685 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 686 "clarfx.f"
	    i__2 = j + c_dim1 * 5;
#line 686 "clarfx.f"
	    i__3 = j + c_dim1 * 5;
#line 686 "clarfx.f"
	    z__2.r = sum.r * t5.r - sum.i * t5.i, z__2.i = sum.r * t5.i + 
		    sum.i * t5.r;
#line 686 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 686 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 687 "clarfx.f"
	    i__2 = j + c_dim1 * 6;
#line 687 "clarfx.f"
	    i__3 = j + c_dim1 * 6;
#line 687 "clarfx.f"
	    z__2.r = sum.r * t6.r - sum.i * t6.i, z__2.i = sum.r * t6.i + 
		    sum.i * t6.r;
#line 687 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 687 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 688 "clarfx.f"
	    i__2 = j + c_dim1 * 7;
#line 688 "clarfx.f"
	    i__3 = j + c_dim1 * 7;
#line 688 "clarfx.f"
	    z__2.r = sum.r * t7.r - sum.i * t7.i, z__2.i = sum.r * t7.i + 
		    sum.i * t7.r;
#line 688 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 688 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 689 "clarfx.f"
	    i__2 = j + (c_dim1 << 3);
#line 689 "clarfx.f"
	    i__3 = j + (c_dim1 << 3);
#line 689 "clarfx.f"
	    z__2.r = sum.r * t8.r - sum.i * t8.i, z__2.i = sum.r * t8.i + 
		    sum.i * t8.r;
#line 689 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 689 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 690 "clarfx.f"
	    i__2 = j + c_dim1 * 9;
#line 690 "clarfx.f"
	    i__3 = j + c_dim1 * 9;
#line 690 "clarfx.f"
	    z__2.r = sum.r * t9.r - sum.i * t9.i, z__2.i = sum.r * t9.i + 
		    sum.i * t9.r;
#line 690 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 690 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 691 "clarfx.f"
	    i__2 = j + c_dim1 * 10;
#line 691 "clarfx.f"
	    i__3 = j + c_dim1 * 10;
#line 691 "clarfx.f"
	    z__2.r = sum.r * t10.r - sum.i * t10.i, z__2.i = sum.r * t10.i + 
		    sum.i * t10.r;
#line 691 "clarfx.f"
	    z__1.r = c__[i__3].r - z__2.r, z__1.i = c__[i__3].i - z__2.i;
#line 691 "clarfx.f"
	    c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
#line 692 "clarfx.f"
/* L400: */
#line 692 "clarfx.f"
	}
#line 693 "clarfx.f"
	goto L410;
#line 694 "clarfx.f"
    }
#line 695 "clarfx.f"
L410:
#line 695 "clarfx.f"
    return 0;

/*     End of CLARFX */

} /* clarfx_ */

