#line 1 "dlarfx.f"
/* dlarfx.f -- translated by f2c (version 20100827).
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

#line 1 "dlarfx.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b DLARFX applies an elementary reflector to a general rectangular matrix, with loop unrolling whe
n the reflector has order ≤ 10. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLARFX + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlarfx.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlarfx.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlarfx.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLARFX( SIDE, M, N, V, TAU, C, LDC, WORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          SIDE */
/*       INTEGER            LDC, M, N */
/*       DOUBLE PRECISION   TAU */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   C( LDC, * ), V( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLARFX applies a real elementary reflector H to a real m by n */
/* > matrix C, from either the left or the right. H is represented in the */
/* > form */
/* > */
/* >       H = I - tau * v * v**T */
/* > */
/* > where tau is a real scalar and v is a real vector. */
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
/* >          V is DOUBLE PRECISION array, dimension (M) if SIDE = 'L' */
/* >                                     or (N) if SIDE = 'R' */
/* >          The vector v in the representation of H. */
/* > \endverbatim */
/* > */
/* > \param[in] TAU */
/* > \verbatim */
/* >          TAU is DOUBLE PRECISION */
/* >          The value tau in the representation of H. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* >          C is DOUBLE PRECISION array, dimension (LDC,N) */
/* >          On entry, the m by n matrix C. */
/* >          On exit, C is overwritten by the matrix H * C if SIDE = 'L', */
/* >          or C * H if SIDE = 'R'. */
/* > \endverbatim */
/* > */
/* > \param[in] LDC */
/* > \verbatim */
/* >          LDC is INTEGER */
/* >          The leading dimension of the array C. LDA >= (1,M). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension */
/* >                      (N) if SIDE = 'L' */
/* >                      or (M) if SIDE = 'R' */
/* >          WORK is not referenced if H has order < 11. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup doubleOTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int dlarfx_(char *side, integer *m, integer *n, doublereal *
	v, doublereal *tau, doublereal *c__, integer *ldc, doublereal *work, 
	ftnlen side_len)
{
    /* System generated locals */
    integer c_dim1, c_offset, i__1;

    /* Local variables */
    static integer j;
    static doublereal t1, t2, t3, t4, t5, t6, t7, t8, t9, v1, v2, v3, v4, v5, 
	    v6, v7, v8, v9, t10, v10, sum;
    extern /* Subroutine */ int dlarf_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, ftnlen);
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
/*     .. Executable Statements .. */

#line 157 "dlarfx.f"
    /* Parameter adjustments */
#line 157 "dlarfx.f"
    --v;
#line 157 "dlarfx.f"
    c_dim1 = *ldc;
#line 157 "dlarfx.f"
    c_offset = 1 + c_dim1;
#line 157 "dlarfx.f"
    c__ -= c_offset;
#line 157 "dlarfx.f"
    --work;
#line 157 "dlarfx.f"

#line 157 "dlarfx.f"
    /* Function Body */
#line 157 "dlarfx.f"
    if (*tau == 0.) {
#line 157 "dlarfx.f"
	return 0;
#line 157 "dlarfx.f"
    }
#line 159 "dlarfx.f"
    if (lsame_(side, "L", (ftnlen)1, (ftnlen)1)) {

/*        Form  H * C, where H has order m. */

#line 163 "dlarfx.f"
	switch (*m) {
#line 163 "dlarfx.f"
	    case 1:  goto L10;
#line 163 "dlarfx.f"
	    case 2:  goto L30;
#line 163 "dlarfx.f"
	    case 3:  goto L50;
#line 163 "dlarfx.f"
	    case 4:  goto L70;
#line 163 "dlarfx.f"
	    case 5:  goto L90;
#line 163 "dlarfx.f"
	    case 6:  goto L110;
#line 163 "dlarfx.f"
	    case 7:  goto L130;
#line 163 "dlarfx.f"
	    case 8:  goto L150;
#line 163 "dlarfx.f"
	    case 9:  goto L170;
#line 163 "dlarfx.f"
	    case 10:  goto L190;
#line 163 "dlarfx.f"
	}

/*        Code for general M */

#line 168 "dlarfx.f"
	dlarf_(side, m, n, &v[1], &c__1, tau, &c__[c_offset], ldc, &work[1], (
		ftnlen)1);
#line 169 "dlarfx.f"
	goto L410;
#line 170 "dlarfx.f"
L10:

/*        Special code for 1 x 1 Householder */

#line 174 "dlarfx.f"
	t1 = 1. - *tau * v[1] * v[1];
#line 175 "dlarfx.f"
	i__1 = *n;
#line 175 "dlarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 176 "dlarfx.f"
	    c__[j * c_dim1 + 1] = t1 * c__[j * c_dim1 + 1];
#line 177 "dlarfx.f"
/* L20: */
#line 177 "dlarfx.f"
	}
#line 178 "dlarfx.f"
	goto L410;
#line 179 "dlarfx.f"
L30:

/*        Special code for 2 x 2 Householder */

#line 183 "dlarfx.f"
	v1 = v[1];
#line 184 "dlarfx.f"
	t1 = *tau * v1;
#line 185 "dlarfx.f"
	v2 = v[2];
#line 186 "dlarfx.f"
	t2 = *tau * v2;
#line 187 "dlarfx.f"
	i__1 = *n;
#line 187 "dlarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 188 "dlarfx.f"
	    sum = v1 * c__[j * c_dim1 + 1] + v2 * c__[j * c_dim1 + 2];
#line 189 "dlarfx.f"
	    c__[j * c_dim1 + 1] -= sum * t1;
#line 190 "dlarfx.f"
	    c__[j * c_dim1 + 2] -= sum * t2;
#line 191 "dlarfx.f"
/* L40: */
#line 191 "dlarfx.f"
	}
#line 192 "dlarfx.f"
	goto L410;
#line 193 "dlarfx.f"
L50:

/*        Special code for 3 x 3 Householder */

#line 197 "dlarfx.f"
	v1 = v[1];
#line 198 "dlarfx.f"
	t1 = *tau * v1;
#line 199 "dlarfx.f"
	v2 = v[2];
#line 200 "dlarfx.f"
	t2 = *tau * v2;
#line 201 "dlarfx.f"
	v3 = v[3];
#line 202 "dlarfx.f"
	t3 = *tau * v3;
#line 203 "dlarfx.f"
	i__1 = *n;
#line 203 "dlarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 204 "dlarfx.f"
	    sum = v1 * c__[j * c_dim1 + 1] + v2 * c__[j * c_dim1 + 2] + v3 * 
		    c__[j * c_dim1 + 3];
#line 205 "dlarfx.f"
	    c__[j * c_dim1 + 1] -= sum * t1;
#line 206 "dlarfx.f"
	    c__[j * c_dim1 + 2] -= sum * t2;
#line 207 "dlarfx.f"
	    c__[j * c_dim1 + 3] -= sum * t3;
#line 208 "dlarfx.f"
/* L60: */
#line 208 "dlarfx.f"
	}
#line 209 "dlarfx.f"
	goto L410;
#line 210 "dlarfx.f"
L70:

/*        Special code for 4 x 4 Householder */

#line 214 "dlarfx.f"
	v1 = v[1];
#line 215 "dlarfx.f"
	t1 = *tau * v1;
#line 216 "dlarfx.f"
	v2 = v[2];
#line 217 "dlarfx.f"
	t2 = *tau * v2;
#line 218 "dlarfx.f"
	v3 = v[3];
#line 219 "dlarfx.f"
	t3 = *tau * v3;
#line 220 "dlarfx.f"
	v4 = v[4];
#line 221 "dlarfx.f"
	t4 = *tau * v4;
#line 222 "dlarfx.f"
	i__1 = *n;
#line 222 "dlarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 223 "dlarfx.f"
	    sum = v1 * c__[j * c_dim1 + 1] + v2 * c__[j * c_dim1 + 2] + v3 * 
		    c__[j * c_dim1 + 3] + v4 * c__[j * c_dim1 + 4];
#line 225 "dlarfx.f"
	    c__[j * c_dim1 + 1] -= sum * t1;
#line 226 "dlarfx.f"
	    c__[j * c_dim1 + 2] -= sum * t2;
#line 227 "dlarfx.f"
	    c__[j * c_dim1 + 3] -= sum * t3;
#line 228 "dlarfx.f"
	    c__[j * c_dim1 + 4] -= sum * t4;
#line 229 "dlarfx.f"
/* L80: */
#line 229 "dlarfx.f"
	}
#line 230 "dlarfx.f"
	goto L410;
#line 231 "dlarfx.f"
L90:

/*        Special code for 5 x 5 Householder */

#line 235 "dlarfx.f"
	v1 = v[1];
#line 236 "dlarfx.f"
	t1 = *tau * v1;
#line 237 "dlarfx.f"
	v2 = v[2];
#line 238 "dlarfx.f"
	t2 = *tau * v2;
#line 239 "dlarfx.f"
	v3 = v[3];
#line 240 "dlarfx.f"
	t3 = *tau * v3;
#line 241 "dlarfx.f"
	v4 = v[4];
#line 242 "dlarfx.f"
	t4 = *tau * v4;
#line 243 "dlarfx.f"
	v5 = v[5];
#line 244 "dlarfx.f"
	t5 = *tau * v5;
#line 245 "dlarfx.f"
	i__1 = *n;
#line 245 "dlarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 246 "dlarfx.f"
	    sum = v1 * c__[j * c_dim1 + 1] + v2 * c__[j * c_dim1 + 2] + v3 * 
		    c__[j * c_dim1 + 3] + v4 * c__[j * c_dim1 + 4] + v5 * c__[
		    j * c_dim1 + 5];
#line 248 "dlarfx.f"
	    c__[j * c_dim1 + 1] -= sum * t1;
#line 249 "dlarfx.f"
	    c__[j * c_dim1 + 2] -= sum * t2;
#line 250 "dlarfx.f"
	    c__[j * c_dim1 + 3] -= sum * t3;
#line 251 "dlarfx.f"
	    c__[j * c_dim1 + 4] -= sum * t4;
#line 252 "dlarfx.f"
	    c__[j * c_dim1 + 5] -= sum * t5;
#line 253 "dlarfx.f"
/* L100: */
#line 253 "dlarfx.f"
	}
#line 254 "dlarfx.f"
	goto L410;
#line 255 "dlarfx.f"
L110:

/*        Special code for 6 x 6 Householder */

#line 259 "dlarfx.f"
	v1 = v[1];
#line 260 "dlarfx.f"
	t1 = *tau * v1;
#line 261 "dlarfx.f"
	v2 = v[2];
#line 262 "dlarfx.f"
	t2 = *tau * v2;
#line 263 "dlarfx.f"
	v3 = v[3];
#line 264 "dlarfx.f"
	t3 = *tau * v3;
#line 265 "dlarfx.f"
	v4 = v[4];
#line 266 "dlarfx.f"
	t4 = *tau * v4;
#line 267 "dlarfx.f"
	v5 = v[5];
#line 268 "dlarfx.f"
	t5 = *tau * v5;
#line 269 "dlarfx.f"
	v6 = v[6];
#line 270 "dlarfx.f"
	t6 = *tau * v6;
#line 271 "dlarfx.f"
	i__1 = *n;
#line 271 "dlarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 272 "dlarfx.f"
	    sum = v1 * c__[j * c_dim1 + 1] + v2 * c__[j * c_dim1 + 2] + v3 * 
		    c__[j * c_dim1 + 3] + v4 * c__[j * c_dim1 + 4] + v5 * c__[
		    j * c_dim1 + 5] + v6 * c__[j * c_dim1 + 6];
#line 274 "dlarfx.f"
	    c__[j * c_dim1 + 1] -= sum * t1;
#line 275 "dlarfx.f"
	    c__[j * c_dim1 + 2] -= sum * t2;
#line 276 "dlarfx.f"
	    c__[j * c_dim1 + 3] -= sum * t3;
#line 277 "dlarfx.f"
	    c__[j * c_dim1 + 4] -= sum * t4;
#line 278 "dlarfx.f"
	    c__[j * c_dim1 + 5] -= sum * t5;
#line 279 "dlarfx.f"
	    c__[j * c_dim1 + 6] -= sum * t6;
#line 280 "dlarfx.f"
/* L120: */
#line 280 "dlarfx.f"
	}
#line 281 "dlarfx.f"
	goto L410;
#line 282 "dlarfx.f"
L130:

/*        Special code for 7 x 7 Householder */

#line 286 "dlarfx.f"
	v1 = v[1];
#line 287 "dlarfx.f"
	t1 = *tau * v1;
#line 288 "dlarfx.f"
	v2 = v[2];
#line 289 "dlarfx.f"
	t2 = *tau * v2;
#line 290 "dlarfx.f"
	v3 = v[3];
#line 291 "dlarfx.f"
	t3 = *tau * v3;
#line 292 "dlarfx.f"
	v4 = v[4];
#line 293 "dlarfx.f"
	t4 = *tau * v4;
#line 294 "dlarfx.f"
	v5 = v[5];
#line 295 "dlarfx.f"
	t5 = *tau * v5;
#line 296 "dlarfx.f"
	v6 = v[6];
#line 297 "dlarfx.f"
	t6 = *tau * v6;
#line 298 "dlarfx.f"
	v7 = v[7];
#line 299 "dlarfx.f"
	t7 = *tau * v7;
#line 300 "dlarfx.f"
	i__1 = *n;
#line 300 "dlarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 301 "dlarfx.f"
	    sum = v1 * c__[j * c_dim1 + 1] + v2 * c__[j * c_dim1 + 2] + v3 * 
		    c__[j * c_dim1 + 3] + v4 * c__[j * c_dim1 + 4] + v5 * c__[
		    j * c_dim1 + 5] + v6 * c__[j * c_dim1 + 6] + v7 * c__[j * 
		    c_dim1 + 7];
#line 304 "dlarfx.f"
	    c__[j * c_dim1 + 1] -= sum * t1;
#line 305 "dlarfx.f"
	    c__[j * c_dim1 + 2] -= sum * t2;
#line 306 "dlarfx.f"
	    c__[j * c_dim1 + 3] -= sum * t3;
#line 307 "dlarfx.f"
	    c__[j * c_dim1 + 4] -= sum * t4;
#line 308 "dlarfx.f"
	    c__[j * c_dim1 + 5] -= sum * t5;
#line 309 "dlarfx.f"
	    c__[j * c_dim1 + 6] -= sum * t6;
#line 310 "dlarfx.f"
	    c__[j * c_dim1 + 7] -= sum * t7;
#line 311 "dlarfx.f"
/* L140: */
#line 311 "dlarfx.f"
	}
#line 312 "dlarfx.f"
	goto L410;
#line 313 "dlarfx.f"
L150:

/*        Special code for 8 x 8 Householder */

#line 317 "dlarfx.f"
	v1 = v[1];
#line 318 "dlarfx.f"
	t1 = *tau * v1;
#line 319 "dlarfx.f"
	v2 = v[2];
#line 320 "dlarfx.f"
	t2 = *tau * v2;
#line 321 "dlarfx.f"
	v3 = v[3];
#line 322 "dlarfx.f"
	t3 = *tau * v3;
#line 323 "dlarfx.f"
	v4 = v[4];
#line 324 "dlarfx.f"
	t4 = *tau * v4;
#line 325 "dlarfx.f"
	v5 = v[5];
#line 326 "dlarfx.f"
	t5 = *tau * v5;
#line 327 "dlarfx.f"
	v6 = v[6];
#line 328 "dlarfx.f"
	t6 = *tau * v6;
#line 329 "dlarfx.f"
	v7 = v[7];
#line 330 "dlarfx.f"
	t7 = *tau * v7;
#line 331 "dlarfx.f"
	v8 = v[8];
#line 332 "dlarfx.f"
	t8 = *tau * v8;
#line 333 "dlarfx.f"
	i__1 = *n;
#line 333 "dlarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 334 "dlarfx.f"
	    sum = v1 * c__[j * c_dim1 + 1] + v2 * c__[j * c_dim1 + 2] + v3 * 
		    c__[j * c_dim1 + 3] + v4 * c__[j * c_dim1 + 4] + v5 * c__[
		    j * c_dim1 + 5] + v6 * c__[j * c_dim1 + 6] + v7 * c__[j * 
		    c_dim1 + 7] + v8 * c__[j * c_dim1 + 8];
#line 337 "dlarfx.f"
	    c__[j * c_dim1 + 1] -= sum * t1;
#line 338 "dlarfx.f"
	    c__[j * c_dim1 + 2] -= sum * t2;
#line 339 "dlarfx.f"
	    c__[j * c_dim1 + 3] -= sum * t3;
#line 340 "dlarfx.f"
	    c__[j * c_dim1 + 4] -= sum * t4;
#line 341 "dlarfx.f"
	    c__[j * c_dim1 + 5] -= sum * t5;
#line 342 "dlarfx.f"
	    c__[j * c_dim1 + 6] -= sum * t6;
#line 343 "dlarfx.f"
	    c__[j * c_dim1 + 7] -= sum * t7;
#line 344 "dlarfx.f"
	    c__[j * c_dim1 + 8] -= sum * t8;
#line 345 "dlarfx.f"
/* L160: */
#line 345 "dlarfx.f"
	}
#line 346 "dlarfx.f"
	goto L410;
#line 347 "dlarfx.f"
L170:

/*        Special code for 9 x 9 Householder */

#line 351 "dlarfx.f"
	v1 = v[1];
#line 352 "dlarfx.f"
	t1 = *tau * v1;
#line 353 "dlarfx.f"
	v2 = v[2];
#line 354 "dlarfx.f"
	t2 = *tau * v2;
#line 355 "dlarfx.f"
	v3 = v[3];
#line 356 "dlarfx.f"
	t3 = *tau * v3;
#line 357 "dlarfx.f"
	v4 = v[4];
#line 358 "dlarfx.f"
	t4 = *tau * v4;
#line 359 "dlarfx.f"
	v5 = v[5];
#line 360 "dlarfx.f"
	t5 = *tau * v5;
#line 361 "dlarfx.f"
	v6 = v[6];
#line 362 "dlarfx.f"
	t6 = *tau * v6;
#line 363 "dlarfx.f"
	v7 = v[7];
#line 364 "dlarfx.f"
	t7 = *tau * v7;
#line 365 "dlarfx.f"
	v8 = v[8];
#line 366 "dlarfx.f"
	t8 = *tau * v8;
#line 367 "dlarfx.f"
	v9 = v[9];
#line 368 "dlarfx.f"
	t9 = *tau * v9;
#line 369 "dlarfx.f"
	i__1 = *n;
#line 369 "dlarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 370 "dlarfx.f"
	    sum = v1 * c__[j * c_dim1 + 1] + v2 * c__[j * c_dim1 + 2] + v3 * 
		    c__[j * c_dim1 + 3] + v4 * c__[j * c_dim1 + 4] + v5 * c__[
		    j * c_dim1 + 5] + v6 * c__[j * c_dim1 + 6] + v7 * c__[j * 
		    c_dim1 + 7] + v8 * c__[j * c_dim1 + 8] + v9 * c__[j * 
		    c_dim1 + 9];
#line 373 "dlarfx.f"
	    c__[j * c_dim1 + 1] -= sum * t1;
#line 374 "dlarfx.f"
	    c__[j * c_dim1 + 2] -= sum * t2;
#line 375 "dlarfx.f"
	    c__[j * c_dim1 + 3] -= sum * t3;
#line 376 "dlarfx.f"
	    c__[j * c_dim1 + 4] -= sum * t4;
#line 377 "dlarfx.f"
	    c__[j * c_dim1 + 5] -= sum * t5;
#line 378 "dlarfx.f"
	    c__[j * c_dim1 + 6] -= sum * t6;
#line 379 "dlarfx.f"
	    c__[j * c_dim1 + 7] -= sum * t7;
#line 380 "dlarfx.f"
	    c__[j * c_dim1 + 8] -= sum * t8;
#line 381 "dlarfx.f"
	    c__[j * c_dim1 + 9] -= sum * t9;
#line 382 "dlarfx.f"
/* L180: */
#line 382 "dlarfx.f"
	}
#line 383 "dlarfx.f"
	goto L410;
#line 384 "dlarfx.f"
L190:

/*        Special code for 10 x 10 Householder */

#line 388 "dlarfx.f"
	v1 = v[1];
#line 389 "dlarfx.f"
	t1 = *tau * v1;
#line 390 "dlarfx.f"
	v2 = v[2];
#line 391 "dlarfx.f"
	t2 = *tau * v2;
#line 392 "dlarfx.f"
	v3 = v[3];
#line 393 "dlarfx.f"
	t3 = *tau * v3;
#line 394 "dlarfx.f"
	v4 = v[4];
#line 395 "dlarfx.f"
	t4 = *tau * v4;
#line 396 "dlarfx.f"
	v5 = v[5];
#line 397 "dlarfx.f"
	t5 = *tau * v5;
#line 398 "dlarfx.f"
	v6 = v[6];
#line 399 "dlarfx.f"
	t6 = *tau * v6;
#line 400 "dlarfx.f"
	v7 = v[7];
#line 401 "dlarfx.f"
	t7 = *tau * v7;
#line 402 "dlarfx.f"
	v8 = v[8];
#line 403 "dlarfx.f"
	t8 = *tau * v8;
#line 404 "dlarfx.f"
	v9 = v[9];
#line 405 "dlarfx.f"
	t9 = *tau * v9;
#line 406 "dlarfx.f"
	v10 = v[10];
#line 407 "dlarfx.f"
	t10 = *tau * v10;
#line 408 "dlarfx.f"
	i__1 = *n;
#line 408 "dlarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 409 "dlarfx.f"
	    sum = v1 * c__[j * c_dim1 + 1] + v2 * c__[j * c_dim1 + 2] + v3 * 
		    c__[j * c_dim1 + 3] + v4 * c__[j * c_dim1 + 4] + v5 * c__[
		    j * c_dim1 + 5] + v6 * c__[j * c_dim1 + 6] + v7 * c__[j * 
		    c_dim1 + 7] + v8 * c__[j * c_dim1 + 8] + v9 * c__[j * 
		    c_dim1 + 9] + v10 * c__[j * c_dim1 + 10];
#line 413 "dlarfx.f"
	    c__[j * c_dim1 + 1] -= sum * t1;
#line 414 "dlarfx.f"
	    c__[j * c_dim1 + 2] -= sum * t2;
#line 415 "dlarfx.f"
	    c__[j * c_dim1 + 3] -= sum * t3;
#line 416 "dlarfx.f"
	    c__[j * c_dim1 + 4] -= sum * t4;
#line 417 "dlarfx.f"
	    c__[j * c_dim1 + 5] -= sum * t5;
#line 418 "dlarfx.f"
	    c__[j * c_dim1 + 6] -= sum * t6;
#line 419 "dlarfx.f"
	    c__[j * c_dim1 + 7] -= sum * t7;
#line 420 "dlarfx.f"
	    c__[j * c_dim1 + 8] -= sum * t8;
#line 421 "dlarfx.f"
	    c__[j * c_dim1 + 9] -= sum * t9;
#line 422 "dlarfx.f"
	    c__[j * c_dim1 + 10] -= sum * t10;
#line 423 "dlarfx.f"
/* L200: */
#line 423 "dlarfx.f"
	}
#line 424 "dlarfx.f"
	goto L410;
#line 425 "dlarfx.f"
    } else {

/*        Form  C * H, where H has order n. */

#line 429 "dlarfx.f"
	switch (*n) {
#line 429 "dlarfx.f"
	    case 1:  goto L210;
#line 429 "dlarfx.f"
	    case 2:  goto L230;
#line 429 "dlarfx.f"
	    case 3:  goto L250;
#line 429 "dlarfx.f"
	    case 4:  goto L270;
#line 429 "dlarfx.f"
	    case 5:  goto L290;
#line 429 "dlarfx.f"
	    case 6:  goto L310;
#line 429 "dlarfx.f"
	    case 7:  goto L330;
#line 429 "dlarfx.f"
	    case 8:  goto L350;
#line 429 "dlarfx.f"
	    case 9:  goto L370;
#line 429 "dlarfx.f"
	    case 10:  goto L390;
#line 429 "dlarfx.f"
	}

/*        Code for general N */

#line 434 "dlarfx.f"
	dlarf_(side, m, n, &v[1], &c__1, tau, &c__[c_offset], ldc, &work[1], (
		ftnlen)1);
#line 435 "dlarfx.f"
	goto L410;
#line 436 "dlarfx.f"
L210:

/*        Special code for 1 x 1 Householder */

#line 440 "dlarfx.f"
	t1 = 1. - *tau * v[1] * v[1];
#line 441 "dlarfx.f"
	i__1 = *m;
#line 441 "dlarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 442 "dlarfx.f"
	    c__[j + c_dim1] = t1 * c__[j + c_dim1];
#line 443 "dlarfx.f"
/* L220: */
#line 443 "dlarfx.f"
	}
#line 444 "dlarfx.f"
	goto L410;
#line 445 "dlarfx.f"
L230:

/*        Special code for 2 x 2 Householder */

#line 449 "dlarfx.f"
	v1 = v[1];
#line 450 "dlarfx.f"
	t1 = *tau * v1;
#line 451 "dlarfx.f"
	v2 = v[2];
#line 452 "dlarfx.f"
	t2 = *tau * v2;
#line 453 "dlarfx.f"
	i__1 = *m;
#line 453 "dlarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 454 "dlarfx.f"
	    sum = v1 * c__[j + c_dim1] + v2 * c__[j + (c_dim1 << 1)];
#line 455 "dlarfx.f"
	    c__[j + c_dim1] -= sum * t1;
#line 456 "dlarfx.f"
	    c__[j + (c_dim1 << 1)] -= sum * t2;
#line 457 "dlarfx.f"
/* L240: */
#line 457 "dlarfx.f"
	}
#line 458 "dlarfx.f"
	goto L410;
#line 459 "dlarfx.f"
L250:

/*        Special code for 3 x 3 Householder */

#line 463 "dlarfx.f"
	v1 = v[1];
#line 464 "dlarfx.f"
	t1 = *tau * v1;
#line 465 "dlarfx.f"
	v2 = v[2];
#line 466 "dlarfx.f"
	t2 = *tau * v2;
#line 467 "dlarfx.f"
	v3 = v[3];
#line 468 "dlarfx.f"
	t3 = *tau * v3;
#line 469 "dlarfx.f"
	i__1 = *m;
#line 469 "dlarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 470 "dlarfx.f"
	    sum = v1 * c__[j + c_dim1] + v2 * c__[j + (c_dim1 << 1)] + v3 * 
		    c__[j + c_dim1 * 3];
#line 471 "dlarfx.f"
	    c__[j + c_dim1] -= sum * t1;
#line 472 "dlarfx.f"
	    c__[j + (c_dim1 << 1)] -= sum * t2;
#line 473 "dlarfx.f"
	    c__[j + c_dim1 * 3] -= sum * t3;
#line 474 "dlarfx.f"
/* L260: */
#line 474 "dlarfx.f"
	}
#line 475 "dlarfx.f"
	goto L410;
#line 476 "dlarfx.f"
L270:

/*        Special code for 4 x 4 Householder */

#line 480 "dlarfx.f"
	v1 = v[1];
#line 481 "dlarfx.f"
	t1 = *tau * v1;
#line 482 "dlarfx.f"
	v2 = v[2];
#line 483 "dlarfx.f"
	t2 = *tau * v2;
#line 484 "dlarfx.f"
	v3 = v[3];
#line 485 "dlarfx.f"
	t3 = *tau * v3;
#line 486 "dlarfx.f"
	v4 = v[4];
#line 487 "dlarfx.f"
	t4 = *tau * v4;
#line 488 "dlarfx.f"
	i__1 = *m;
#line 488 "dlarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 489 "dlarfx.f"
	    sum = v1 * c__[j + c_dim1] + v2 * c__[j + (c_dim1 << 1)] + v3 * 
		    c__[j + c_dim1 * 3] + v4 * c__[j + (c_dim1 << 2)];
#line 491 "dlarfx.f"
	    c__[j + c_dim1] -= sum * t1;
#line 492 "dlarfx.f"
	    c__[j + (c_dim1 << 1)] -= sum * t2;
#line 493 "dlarfx.f"
	    c__[j + c_dim1 * 3] -= sum * t3;
#line 494 "dlarfx.f"
	    c__[j + (c_dim1 << 2)] -= sum * t4;
#line 495 "dlarfx.f"
/* L280: */
#line 495 "dlarfx.f"
	}
#line 496 "dlarfx.f"
	goto L410;
#line 497 "dlarfx.f"
L290:

/*        Special code for 5 x 5 Householder */

#line 501 "dlarfx.f"
	v1 = v[1];
#line 502 "dlarfx.f"
	t1 = *tau * v1;
#line 503 "dlarfx.f"
	v2 = v[2];
#line 504 "dlarfx.f"
	t2 = *tau * v2;
#line 505 "dlarfx.f"
	v3 = v[3];
#line 506 "dlarfx.f"
	t3 = *tau * v3;
#line 507 "dlarfx.f"
	v4 = v[4];
#line 508 "dlarfx.f"
	t4 = *tau * v4;
#line 509 "dlarfx.f"
	v5 = v[5];
#line 510 "dlarfx.f"
	t5 = *tau * v5;
#line 511 "dlarfx.f"
	i__1 = *m;
#line 511 "dlarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 512 "dlarfx.f"
	    sum = v1 * c__[j + c_dim1] + v2 * c__[j + (c_dim1 << 1)] + v3 * 
		    c__[j + c_dim1 * 3] + v4 * c__[j + (c_dim1 << 2)] + v5 * 
		    c__[j + c_dim1 * 5];
#line 514 "dlarfx.f"
	    c__[j + c_dim1] -= sum * t1;
#line 515 "dlarfx.f"
	    c__[j + (c_dim1 << 1)] -= sum * t2;
#line 516 "dlarfx.f"
	    c__[j + c_dim1 * 3] -= sum * t3;
#line 517 "dlarfx.f"
	    c__[j + (c_dim1 << 2)] -= sum * t4;
#line 518 "dlarfx.f"
	    c__[j + c_dim1 * 5] -= sum * t5;
#line 519 "dlarfx.f"
/* L300: */
#line 519 "dlarfx.f"
	}
#line 520 "dlarfx.f"
	goto L410;
#line 521 "dlarfx.f"
L310:

/*        Special code for 6 x 6 Householder */

#line 525 "dlarfx.f"
	v1 = v[1];
#line 526 "dlarfx.f"
	t1 = *tau * v1;
#line 527 "dlarfx.f"
	v2 = v[2];
#line 528 "dlarfx.f"
	t2 = *tau * v2;
#line 529 "dlarfx.f"
	v3 = v[3];
#line 530 "dlarfx.f"
	t3 = *tau * v3;
#line 531 "dlarfx.f"
	v4 = v[4];
#line 532 "dlarfx.f"
	t4 = *tau * v4;
#line 533 "dlarfx.f"
	v5 = v[5];
#line 534 "dlarfx.f"
	t5 = *tau * v5;
#line 535 "dlarfx.f"
	v6 = v[6];
#line 536 "dlarfx.f"
	t6 = *tau * v6;
#line 537 "dlarfx.f"
	i__1 = *m;
#line 537 "dlarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 538 "dlarfx.f"
	    sum = v1 * c__[j + c_dim1] + v2 * c__[j + (c_dim1 << 1)] + v3 * 
		    c__[j + c_dim1 * 3] + v4 * c__[j + (c_dim1 << 2)] + v5 * 
		    c__[j + c_dim1 * 5] + v6 * c__[j + c_dim1 * 6];
#line 540 "dlarfx.f"
	    c__[j + c_dim1] -= sum * t1;
#line 541 "dlarfx.f"
	    c__[j + (c_dim1 << 1)] -= sum * t2;
#line 542 "dlarfx.f"
	    c__[j + c_dim1 * 3] -= sum * t3;
#line 543 "dlarfx.f"
	    c__[j + (c_dim1 << 2)] -= sum * t4;
#line 544 "dlarfx.f"
	    c__[j + c_dim1 * 5] -= sum * t5;
#line 545 "dlarfx.f"
	    c__[j + c_dim1 * 6] -= sum * t6;
#line 546 "dlarfx.f"
/* L320: */
#line 546 "dlarfx.f"
	}
#line 547 "dlarfx.f"
	goto L410;
#line 548 "dlarfx.f"
L330:

/*        Special code for 7 x 7 Householder */

#line 552 "dlarfx.f"
	v1 = v[1];
#line 553 "dlarfx.f"
	t1 = *tau * v1;
#line 554 "dlarfx.f"
	v2 = v[2];
#line 555 "dlarfx.f"
	t2 = *tau * v2;
#line 556 "dlarfx.f"
	v3 = v[3];
#line 557 "dlarfx.f"
	t3 = *tau * v3;
#line 558 "dlarfx.f"
	v4 = v[4];
#line 559 "dlarfx.f"
	t4 = *tau * v4;
#line 560 "dlarfx.f"
	v5 = v[5];
#line 561 "dlarfx.f"
	t5 = *tau * v5;
#line 562 "dlarfx.f"
	v6 = v[6];
#line 563 "dlarfx.f"
	t6 = *tau * v6;
#line 564 "dlarfx.f"
	v7 = v[7];
#line 565 "dlarfx.f"
	t7 = *tau * v7;
#line 566 "dlarfx.f"
	i__1 = *m;
#line 566 "dlarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 567 "dlarfx.f"
	    sum = v1 * c__[j + c_dim1] + v2 * c__[j + (c_dim1 << 1)] + v3 * 
		    c__[j + c_dim1 * 3] + v4 * c__[j + (c_dim1 << 2)] + v5 * 
		    c__[j + c_dim1 * 5] + v6 * c__[j + c_dim1 * 6] + v7 * c__[
		    j + c_dim1 * 7];
#line 570 "dlarfx.f"
	    c__[j + c_dim1] -= sum * t1;
#line 571 "dlarfx.f"
	    c__[j + (c_dim1 << 1)] -= sum * t2;
#line 572 "dlarfx.f"
	    c__[j + c_dim1 * 3] -= sum * t3;
#line 573 "dlarfx.f"
	    c__[j + (c_dim1 << 2)] -= sum * t4;
#line 574 "dlarfx.f"
	    c__[j + c_dim1 * 5] -= sum * t5;
#line 575 "dlarfx.f"
	    c__[j + c_dim1 * 6] -= sum * t6;
#line 576 "dlarfx.f"
	    c__[j + c_dim1 * 7] -= sum * t7;
#line 577 "dlarfx.f"
/* L340: */
#line 577 "dlarfx.f"
	}
#line 578 "dlarfx.f"
	goto L410;
#line 579 "dlarfx.f"
L350:

/*        Special code for 8 x 8 Householder */

#line 583 "dlarfx.f"
	v1 = v[1];
#line 584 "dlarfx.f"
	t1 = *tau * v1;
#line 585 "dlarfx.f"
	v2 = v[2];
#line 586 "dlarfx.f"
	t2 = *tau * v2;
#line 587 "dlarfx.f"
	v3 = v[3];
#line 588 "dlarfx.f"
	t3 = *tau * v3;
#line 589 "dlarfx.f"
	v4 = v[4];
#line 590 "dlarfx.f"
	t4 = *tau * v4;
#line 591 "dlarfx.f"
	v5 = v[5];
#line 592 "dlarfx.f"
	t5 = *tau * v5;
#line 593 "dlarfx.f"
	v6 = v[6];
#line 594 "dlarfx.f"
	t6 = *tau * v6;
#line 595 "dlarfx.f"
	v7 = v[7];
#line 596 "dlarfx.f"
	t7 = *tau * v7;
#line 597 "dlarfx.f"
	v8 = v[8];
#line 598 "dlarfx.f"
	t8 = *tau * v8;
#line 599 "dlarfx.f"
	i__1 = *m;
#line 599 "dlarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 600 "dlarfx.f"
	    sum = v1 * c__[j + c_dim1] + v2 * c__[j + (c_dim1 << 1)] + v3 * 
		    c__[j + c_dim1 * 3] + v4 * c__[j + (c_dim1 << 2)] + v5 * 
		    c__[j + c_dim1 * 5] + v6 * c__[j + c_dim1 * 6] + v7 * c__[
		    j + c_dim1 * 7] + v8 * c__[j + (c_dim1 << 3)];
#line 603 "dlarfx.f"
	    c__[j + c_dim1] -= sum * t1;
#line 604 "dlarfx.f"
	    c__[j + (c_dim1 << 1)] -= sum * t2;
#line 605 "dlarfx.f"
	    c__[j + c_dim1 * 3] -= sum * t3;
#line 606 "dlarfx.f"
	    c__[j + (c_dim1 << 2)] -= sum * t4;
#line 607 "dlarfx.f"
	    c__[j + c_dim1 * 5] -= sum * t5;
#line 608 "dlarfx.f"
	    c__[j + c_dim1 * 6] -= sum * t6;
#line 609 "dlarfx.f"
	    c__[j + c_dim1 * 7] -= sum * t7;
#line 610 "dlarfx.f"
	    c__[j + (c_dim1 << 3)] -= sum * t8;
#line 611 "dlarfx.f"
/* L360: */
#line 611 "dlarfx.f"
	}
#line 612 "dlarfx.f"
	goto L410;
#line 613 "dlarfx.f"
L370:

/*        Special code for 9 x 9 Householder */

#line 617 "dlarfx.f"
	v1 = v[1];
#line 618 "dlarfx.f"
	t1 = *tau * v1;
#line 619 "dlarfx.f"
	v2 = v[2];
#line 620 "dlarfx.f"
	t2 = *tau * v2;
#line 621 "dlarfx.f"
	v3 = v[3];
#line 622 "dlarfx.f"
	t3 = *tau * v3;
#line 623 "dlarfx.f"
	v4 = v[4];
#line 624 "dlarfx.f"
	t4 = *tau * v4;
#line 625 "dlarfx.f"
	v5 = v[5];
#line 626 "dlarfx.f"
	t5 = *tau * v5;
#line 627 "dlarfx.f"
	v6 = v[6];
#line 628 "dlarfx.f"
	t6 = *tau * v6;
#line 629 "dlarfx.f"
	v7 = v[7];
#line 630 "dlarfx.f"
	t7 = *tau * v7;
#line 631 "dlarfx.f"
	v8 = v[8];
#line 632 "dlarfx.f"
	t8 = *tau * v8;
#line 633 "dlarfx.f"
	v9 = v[9];
#line 634 "dlarfx.f"
	t9 = *tau * v9;
#line 635 "dlarfx.f"
	i__1 = *m;
#line 635 "dlarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 636 "dlarfx.f"
	    sum = v1 * c__[j + c_dim1] + v2 * c__[j + (c_dim1 << 1)] + v3 * 
		    c__[j + c_dim1 * 3] + v4 * c__[j + (c_dim1 << 2)] + v5 * 
		    c__[j + c_dim1 * 5] + v6 * c__[j + c_dim1 * 6] + v7 * c__[
		    j + c_dim1 * 7] + v8 * c__[j + (c_dim1 << 3)] + v9 * c__[
		    j + c_dim1 * 9];
#line 639 "dlarfx.f"
	    c__[j + c_dim1] -= sum * t1;
#line 640 "dlarfx.f"
	    c__[j + (c_dim1 << 1)] -= sum * t2;
#line 641 "dlarfx.f"
	    c__[j + c_dim1 * 3] -= sum * t3;
#line 642 "dlarfx.f"
	    c__[j + (c_dim1 << 2)] -= sum * t4;
#line 643 "dlarfx.f"
	    c__[j + c_dim1 * 5] -= sum * t5;
#line 644 "dlarfx.f"
	    c__[j + c_dim1 * 6] -= sum * t6;
#line 645 "dlarfx.f"
	    c__[j + c_dim1 * 7] -= sum * t7;
#line 646 "dlarfx.f"
	    c__[j + (c_dim1 << 3)] -= sum * t8;
#line 647 "dlarfx.f"
	    c__[j + c_dim1 * 9] -= sum * t9;
#line 648 "dlarfx.f"
/* L380: */
#line 648 "dlarfx.f"
	}
#line 649 "dlarfx.f"
	goto L410;
#line 650 "dlarfx.f"
L390:

/*        Special code for 10 x 10 Householder */

#line 654 "dlarfx.f"
	v1 = v[1];
#line 655 "dlarfx.f"
	t1 = *tau * v1;
#line 656 "dlarfx.f"
	v2 = v[2];
#line 657 "dlarfx.f"
	t2 = *tau * v2;
#line 658 "dlarfx.f"
	v3 = v[3];
#line 659 "dlarfx.f"
	t3 = *tau * v3;
#line 660 "dlarfx.f"
	v4 = v[4];
#line 661 "dlarfx.f"
	t4 = *tau * v4;
#line 662 "dlarfx.f"
	v5 = v[5];
#line 663 "dlarfx.f"
	t5 = *tau * v5;
#line 664 "dlarfx.f"
	v6 = v[6];
#line 665 "dlarfx.f"
	t6 = *tau * v6;
#line 666 "dlarfx.f"
	v7 = v[7];
#line 667 "dlarfx.f"
	t7 = *tau * v7;
#line 668 "dlarfx.f"
	v8 = v[8];
#line 669 "dlarfx.f"
	t8 = *tau * v8;
#line 670 "dlarfx.f"
	v9 = v[9];
#line 671 "dlarfx.f"
	t9 = *tau * v9;
#line 672 "dlarfx.f"
	v10 = v[10];
#line 673 "dlarfx.f"
	t10 = *tau * v10;
#line 674 "dlarfx.f"
	i__1 = *m;
#line 674 "dlarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 675 "dlarfx.f"
	    sum = v1 * c__[j + c_dim1] + v2 * c__[j + (c_dim1 << 1)] + v3 * 
		    c__[j + c_dim1 * 3] + v4 * c__[j + (c_dim1 << 2)] + v5 * 
		    c__[j + c_dim1 * 5] + v6 * c__[j + c_dim1 * 6] + v7 * c__[
		    j + c_dim1 * 7] + v8 * c__[j + (c_dim1 << 3)] + v9 * c__[
		    j + c_dim1 * 9] + v10 * c__[j + c_dim1 * 10];
#line 679 "dlarfx.f"
	    c__[j + c_dim1] -= sum * t1;
#line 680 "dlarfx.f"
	    c__[j + (c_dim1 << 1)] -= sum * t2;
#line 681 "dlarfx.f"
	    c__[j + c_dim1 * 3] -= sum * t3;
#line 682 "dlarfx.f"
	    c__[j + (c_dim1 << 2)] -= sum * t4;
#line 683 "dlarfx.f"
	    c__[j + c_dim1 * 5] -= sum * t5;
#line 684 "dlarfx.f"
	    c__[j + c_dim1 * 6] -= sum * t6;
#line 685 "dlarfx.f"
	    c__[j + c_dim1 * 7] -= sum * t7;
#line 686 "dlarfx.f"
	    c__[j + (c_dim1 << 3)] -= sum * t8;
#line 687 "dlarfx.f"
	    c__[j + c_dim1 * 9] -= sum * t9;
#line 688 "dlarfx.f"
	    c__[j + c_dim1 * 10] -= sum * t10;
#line 689 "dlarfx.f"
/* L400: */
#line 689 "dlarfx.f"
	}
#line 690 "dlarfx.f"
	goto L410;
#line 691 "dlarfx.f"
    }
#line 692 "dlarfx.f"
L410:
#line 693 "dlarfx.f"
    return 0;

/*     End of DLARFX */

} /* dlarfx_ */

