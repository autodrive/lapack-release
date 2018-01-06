#line 1 "slarfx.f"
/* slarfx.f -- translated by f2c (version 20100827).
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

#line 1 "slarfx.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b SLARFX applies an elementary reflector to a general rectangular matrix, with loop unrolling whe
n the reflector has order ≤ 10. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLARFX + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slarfx.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slarfx.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slarfx.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLARFX( SIDE, M, N, V, TAU, C, LDC, WORK ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          SIDE */
/*       INTEGER            LDC, M, N */
/*       REAL               TAU */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               C( LDC, * ), V( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLARFX applies a real elementary reflector H to a real m by n */
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
/* >          V is REAL array, dimension (M) if SIDE = 'L' */
/* >                                     or (N) if SIDE = 'R' */
/* >          The vector v in the representation of H. */
/* > \endverbatim */
/* > */
/* > \param[in] TAU */
/* > \verbatim */
/* >          TAU is REAL */
/* >          The value tau in the representation of H. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* >          C is REAL array, dimension (LDC,N) */
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
/* >          WORK is REAL array, dimension */
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

/* > \date December 2016 */

/* > \ingroup realOTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int slarfx_(char *side, integer *m, integer *n, doublereal *
	v, doublereal *tau, doublereal *c__, integer *ldc, doublereal *work, 
	ftnlen side_len)
{
    /* System generated locals */
    integer c_dim1, c_offset, i__1;

    /* Local variables */
    static integer j;
    static doublereal t1, t2, t3, t4, t5, t6, t7, t8, t9, v1, v2, v3, v4, v5, 
	    v6, v7, v8, v9, t10, v10, sum;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int slarf_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, ftnlen);


/*  -- LAPACK auxiliary routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

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

#line 157 "slarfx.f"
    /* Parameter adjustments */
#line 157 "slarfx.f"
    --v;
#line 157 "slarfx.f"
    c_dim1 = *ldc;
#line 157 "slarfx.f"
    c_offset = 1 + c_dim1;
#line 157 "slarfx.f"
    c__ -= c_offset;
#line 157 "slarfx.f"
    --work;
#line 157 "slarfx.f"

#line 157 "slarfx.f"
    /* Function Body */
#line 157 "slarfx.f"
    if (*tau == 0.) {
#line 157 "slarfx.f"
	return 0;
#line 157 "slarfx.f"
    }
#line 159 "slarfx.f"
    if (lsame_(side, "L", (ftnlen)1, (ftnlen)1)) {

/*        Form  H * C, where H has order m. */

#line 163 "slarfx.f"
	switch (*m) {
#line 163 "slarfx.f"
	    case 1:  goto L10;
#line 163 "slarfx.f"
	    case 2:  goto L30;
#line 163 "slarfx.f"
	    case 3:  goto L50;
#line 163 "slarfx.f"
	    case 4:  goto L70;
#line 163 "slarfx.f"
	    case 5:  goto L90;
#line 163 "slarfx.f"
	    case 6:  goto L110;
#line 163 "slarfx.f"
	    case 7:  goto L130;
#line 163 "slarfx.f"
	    case 8:  goto L150;
#line 163 "slarfx.f"
	    case 9:  goto L170;
#line 163 "slarfx.f"
	    case 10:  goto L190;
#line 163 "slarfx.f"
	}

/*        Code for general M */

#line 168 "slarfx.f"
	slarf_(side, m, n, &v[1], &c__1, tau, &c__[c_offset], ldc, &work[1], (
		ftnlen)1);
#line 169 "slarfx.f"
	goto L410;
#line 170 "slarfx.f"
L10:

/*        Special code for 1 x 1 Householder */

#line 174 "slarfx.f"
	t1 = 1. - *tau * v[1] * v[1];
#line 175 "slarfx.f"
	i__1 = *n;
#line 175 "slarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 176 "slarfx.f"
	    c__[j * c_dim1 + 1] = t1 * c__[j * c_dim1 + 1];
#line 177 "slarfx.f"
/* L20: */
#line 177 "slarfx.f"
	}
#line 178 "slarfx.f"
	goto L410;
#line 179 "slarfx.f"
L30:

/*        Special code for 2 x 2 Householder */

#line 183 "slarfx.f"
	v1 = v[1];
#line 184 "slarfx.f"
	t1 = *tau * v1;
#line 185 "slarfx.f"
	v2 = v[2];
#line 186 "slarfx.f"
	t2 = *tau * v2;
#line 187 "slarfx.f"
	i__1 = *n;
#line 187 "slarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 188 "slarfx.f"
	    sum = v1 * c__[j * c_dim1 + 1] + v2 * c__[j * c_dim1 + 2];
#line 189 "slarfx.f"
	    c__[j * c_dim1 + 1] -= sum * t1;
#line 190 "slarfx.f"
	    c__[j * c_dim1 + 2] -= sum * t2;
#line 191 "slarfx.f"
/* L40: */
#line 191 "slarfx.f"
	}
#line 192 "slarfx.f"
	goto L410;
#line 193 "slarfx.f"
L50:

/*        Special code for 3 x 3 Householder */

#line 197 "slarfx.f"
	v1 = v[1];
#line 198 "slarfx.f"
	t1 = *tau * v1;
#line 199 "slarfx.f"
	v2 = v[2];
#line 200 "slarfx.f"
	t2 = *tau * v2;
#line 201 "slarfx.f"
	v3 = v[3];
#line 202 "slarfx.f"
	t3 = *tau * v3;
#line 203 "slarfx.f"
	i__1 = *n;
#line 203 "slarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 204 "slarfx.f"
	    sum = v1 * c__[j * c_dim1 + 1] + v2 * c__[j * c_dim1 + 2] + v3 * 
		    c__[j * c_dim1 + 3];
#line 205 "slarfx.f"
	    c__[j * c_dim1 + 1] -= sum * t1;
#line 206 "slarfx.f"
	    c__[j * c_dim1 + 2] -= sum * t2;
#line 207 "slarfx.f"
	    c__[j * c_dim1 + 3] -= sum * t3;
#line 208 "slarfx.f"
/* L60: */
#line 208 "slarfx.f"
	}
#line 209 "slarfx.f"
	goto L410;
#line 210 "slarfx.f"
L70:

/*        Special code for 4 x 4 Householder */

#line 214 "slarfx.f"
	v1 = v[1];
#line 215 "slarfx.f"
	t1 = *tau * v1;
#line 216 "slarfx.f"
	v2 = v[2];
#line 217 "slarfx.f"
	t2 = *tau * v2;
#line 218 "slarfx.f"
	v3 = v[3];
#line 219 "slarfx.f"
	t3 = *tau * v3;
#line 220 "slarfx.f"
	v4 = v[4];
#line 221 "slarfx.f"
	t4 = *tau * v4;
#line 222 "slarfx.f"
	i__1 = *n;
#line 222 "slarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 223 "slarfx.f"
	    sum = v1 * c__[j * c_dim1 + 1] + v2 * c__[j * c_dim1 + 2] + v3 * 
		    c__[j * c_dim1 + 3] + v4 * c__[j * c_dim1 + 4];
#line 225 "slarfx.f"
	    c__[j * c_dim1 + 1] -= sum * t1;
#line 226 "slarfx.f"
	    c__[j * c_dim1 + 2] -= sum * t2;
#line 227 "slarfx.f"
	    c__[j * c_dim1 + 3] -= sum * t3;
#line 228 "slarfx.f"
	    c__[j * c_dim1 + 4] -= sum * t4;
#line 229 "slarfx.f"
/* L80: */
#line 229 "slarfx.f"
	}
#line 230 "slarfx.f"
	goto L410;
#line 231 "slarfx.f"
L90:

/*        Special code for 5 x 5 Householder */

#line 235 "slarfx.f"
	v1 = v[1];
#line 236 "slarfx.f"
	t1 = *tau * v1;
#line 237 "slarfx.f"
	v2 = v[2];
#line 238 "slarfx.f"
	t2 = *tau * v2;
#line 239 "slarfx.f"
	v3 = v[3];
#line 240 "slarfx.f"
	t3 = *tau * v3;
#line 241 "slarfx.f"
	v4 = v[4];
#line 242 "slarfx.f"
	t4 = *tau * v4;
#line 243 "slarfx.f"
	v5 = v[5];
#line 244 "slarfx.f"
	t5 = *tau * v5;
#line 245 "slarfx.f"
	i__1 = *n;
#line 245 "slarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 246 "slarfx.f"
	    sum = v1 * c__[j * c_dim1 + 1] + v2 * c__[j * c_dim1 + 2] + v3 * 
		    c__[j * c_dim1 + 3] + v4 * c__[j * c_dim1 + 4] + v5 * c__[
		    j * c_dim1 + 5];
#line 248 "slarfx.f"
	    c__[j * c_dim1 + 1] -= sum * t1;
#line 249 "slarfx.f"
	    c__[j * c_dim1 + 2] -= sum * t2;
#line 250 "slarfx.f"
	    c__[j * c_dim1 + 3] -= sum * t3;
#line 251 "slarfx.f"
	    c__[j * c_dim1 + 4] -= sum * t4;
#line 252 "slarfx.f"
	    c__[j * c_dim1 + 5] -= sum * t5;
#line 253 "slarfx.f"
/* L100: */
#line 253 "slarfx.f"
	}
#line 254 "slarfx.f"
	goto L410;
#line 255 "slarfx.f"
L110:

/*        Special code for 6 x 6 Householder */

#line 259 "slarfx.f"
	v1 = v[1];
#line 260 "slarfx.f"
	t1 = *tau * v1;
#line 261 "slarfx.f"
	v2 = v[2];
#line 262 "slarfx.f"
	t2 = *tau * v2;
#line 263 "slarfx.f"
	v3 = v[3];
#line 264 "slarfx.f"
	t3 = *tau * v3;
#line 265 "slarfx.f"
	v4 = v[4];
#line 266 "slarfx.f"
	t4 = *tau * v4;
#line 267 "slarfx.f"
	v5 = v[5];
#line 268 "slarfx.f"
	t5 = *tau * v5;
#line 269 "slarfx.f"
	v6 = v[6];
#line 270 "slarfx.f"
	t6 = *tau * v6;
#line 271 "slarfx.f"
	i__1 = *n;
#line 271 "slarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 272 "slarfx.f"
	    sum = v1 * c__[j * c_dim1 + 1] + v2 * c__[j * c_dim1 + 2] + v3 * 
		    c__[j * c_dim1 + 3] + v4 * c__[j * c_dim1 + 4] + v5 * c__[
		    j * c_dim1 + 5] + v6 * c__[j * c_dim1 + 6];
#line 274 "slarfx.f"
	    c__[j * c_dim1 + 1] -= sum * t1;
#line 275 "slarfx.f"
	    c__[j * c_dim1 + 2] -= sum * t2;
#line 276 "slarfx.f"
	    c__[j * c_dim1 + 3] -= sum * t3;
#line 277 "slarfx.f"
	    c__[j * c_dim1 + 4] -= sum * t4;
#line 278 "slarfx.f"
	    c__[j * c_dim1 + 5] -= sum * t5;
#line 279 "slarfx.f"
	    c__[j * c_dim1 + 6] -= sum * t6;
#line 280 "slarfx.f"
/* L120: */
#line 280 "slarfx.f"
	}
#line 281 "slarfx.f"
	goto L410;
#line 282 "slarfx.f"
L130:

/*        Special code for 7 x 7 Householder */

#line 286 "slarfx.f"
	v1 = v[1];
#line 287 "slarfx.f"
	t1 = *tau * v1;
#line 288 "slarfx.f"
	v2 = v[2];
#line 289 "slarfx.f"
	t2 = *tau * v2;
#line 290 "slarfx.f"
	v3 = v[3];
#line 291 "slarfx.f"
	t3 = *tau * v3;
#line 292 "slarfx.f"
	v4 = v[4];
#line 293 "slarfx.f"
	t4 = *tau * v4;
#line 294 "slarfx.f"
	v5 = v[5];
#line 295 "slarfx.f"
	t5 = *tau * v5;
#line 296 "slarfx.f"
	v6 = v[6];
#line 297 "slarfx.f"
	t6 = *tau * v6;
#line 298 "slarfx.f"
	v7 = v[7];
#line 299 "slarfx.f"
	t7 = *tau * v7;
#line 300 "slarfx.f"
	i__1 = *n;
#line 300 "slarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 301 "slarfx.f"
	    sum = v1 * c__[j * c_dim1 + 1] + v2 * c__[j * c_dim1 + 2] + v3 * 
		    c__[j * c_dim1 + 3] + v4 * c__[j * c_dim1 + 4] + v5 * c__[
		    j * c_dim1 + 5] + v6 * c__[j * c_dim1 + 6] + v7 * c__[j * 
		    c_dim1 + 7];
#line 304 "slarfx.f"
	    c__[j * c_dim1 + 1] -= sum * t1;
#line 305 "slarfx.f"
	    c__[j * c_dim1 + 2] -= sum * t2;
#line 306 "slarfx.f"
	    c__[j * c_dim1 + 3] -= sum * t3;
#line 307 "slarfx.f"
	    c__[j * c_dim1 + 4] -= sum * t4;
#line 308 "slarfx.f"
	    c__[j * c_dim1 + 5] -= sum * t5;
#line 309 "slarfx.f"
	    c__[j * c_dim1 + 6] -= sum * t6;
#line 310 "slarfx.f"
	    c__[j * c_dim1 + 7] -= sum * t7;
#line 311 "slarfx.f"
/* L140: */
#line 311 "slarfx.f"
	}
#line 312 "slarfx.f"
	goto L410;
#line 313 "slarfx.f"
L150:

/*        Special code for 8 x 8 Householder */

#line 317 "slarfx.f"
	v1 = v[1];
#line 318 "slarfx.f"
	t1 = *tau * v1;
#line 319 "slarfx.f"
	v2 = v[2];
#line 320 "slarfx.f"
	t2 = *tau * v2;
#line 321 "slarfx.f"
	v3 = v[3];
#line 322 "slarfx.f"
	t3 = *tau * v3;
#line 323 "slarfx.f"
	v4 = v[4];
#line 324 "slarfx.f"
	t4 = *tau * v4;
#line 325 "slarfx.f"
	v5 = v[5];
#line 326 "slarfx.f"
	t5 = *tau * v5;
#line 327 "slarfx.f"
	v6 = v[6];
#line 328 "slarfx.f"
	t6 = *tau * v6;
#line 329 "slarfx.f"
	v7 = v[7];
#line 330 "slarfx.f"
	t7 = *tau * v7;
#line 331 "slarfx.f"
	v8 = v[8];
#line 332 "slarfx.f"
	t8 = *tau * v8;
#line 333 "slarfx.f"
	i__1 = *n;
#line 333 "slarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 334 "slarfx.f"
	    sum = v1 * c__[j * c_dim1 + 1] + v2 * c__[j * c_dim1 + 2] + v3 * 
		    c__[j * c_dim1 + 3] + v4 * c__[j * c_dim1 + 4] + v5 * c__[
		    j * c_dim1 + 5] + v6 * c__[j * c_dim1 + 6] + v7 * c__[j * 
		    c_dim1 + 7] + v8 * c__[j * c_dim1 + 8];
#line 337 "slarfx.f"
	    c__[j * c_dim1 + 1] -= sum * t1;
#line 338 "slarfx.f"
	    c__[j * c_dim1 + 2] -= sum * t2;
#line 339 "slarfx.f"
	    c__[j * c_dim1 + 3] -= sum * t3;
#line 340 "slarfx.f"
	    c__[j * c_dim1 + 4] -= sum * t4;
#line 341 "slarfx.f"
	    c__[j * c_dim1 + 5] -= sum * t5;
#line 342 "slarfx.f"
	    c__[j * c_dim1 + 6] -= sum * t6;
#line 343 "slarfx.f"
	    c__[j * c_dim1 + 7] -= sum * t7;
#line 344 "slarfx.f"
	    c__[j * c_dim1 + 8] -= sum * t8;
#line 345 "slarfx.f"
/* L160: */
#line 345 "slarfx.f"
	}
#line 346 "slarfx.f"
	goto L410;
#line 347 "slarfx.f"
L170:

/*        Special code for 9 x 9 Householder */

#line 351 "slarfx.f"
	v1 = v[1];
#line 352 "slarfx.f"
	t1 = *tau * v1;
#line 353 "slarfx.f"
	v2 = v[2];
#line 354 "slarfx.f"
	t2 = *tau * v2;
#line 355 "slarfx.f"
	v3 = v[3];
#line 356 "slarfx.f"
	t3 = *tau * v3;
#line 357 "slarfx.f"
	v4 = v[4];
#line 358 "slarfx.f"
	t4 = *tau * v4;
#line 359 "slarfx.f"
	v5 = v[5];
#line 360 "slarfx.f"
	t5 = *tau * v5;
#line 361 "slarfx.f"
	v6 = v[6];
#line 362 "slarfx.f"
	t6 = *tau * v6;
#line 363 "slarfx.f"
	v7 = v[7];
#line 364 "slarfx.f"
	t7 = *tau * v7;
#line 365 "slarfx.f"
	v8 = v[8];
#line 366 "slarfx.f"
	t8 = *tau * v8;
#line 367 "slarfx.f"
	v9 = v[9];
#line 368 "slarfx.f"
	t9 = *tau * v9;
#line 369 "slarfx.f"
	i__1 = *n;
#line 369 "slarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 370 "slarfx.f"
	    sum = v1 * c__[j * c_dim1 + 1] + v2 * c__[j * c_dim1 + 2] + v3 * 
		    c__[j * c_dim1 + 3] + v4 * c__[j * c_dim1 + 4] + v5 * c__[
		    j * c_dim1 + 5] + v6 * c__[j * c_dim1 + 6] + v7 * c__[j * 
		    c_dim1 + 7] + v8 * c__[j * c_dim1 + 8] + v9 * c__[j * 
		    c_dim1 + 9];
#line 373 "slarfx.f"
	    c__[j * c_dim1 + 1] -= sum * t1;
#line 374 "slarfx.f"
	    c__[j * c_dim1 + 2] -= sum * t2;
#line 375 "slarfx.f"
	    c__[j * c_dim1 + 3] -= sum * t3;
#line 376 "slarfx.f"
	    c__[j * c_dim1 + 4] -= sum * t4;
#line 377 "slarfx.f"
	    c__[j * c_dim1 + 5] -= sum * t5;
#line 378 "slarfx.f"
	    c__[j * c_dim1 + 6] -= sum * t6;
#line 379 "slarfx.f"
	    c__[j * c_dim1 + 7] -= sum * t7;
#line 380 "slarfx.f"
	    c__[j * c_dim1 + 8] -= sum * t8;
#line 381 "slarfx.f"
	    c__[j * c_dim1 + 9] -= sum * t9;
#line 382 "slarfx.f"
/* L180: */
#line 382 "slarfx.f"
	}
#line 383 "slarfx.f"
	goto L410;
#line 384 "slarfx.f"
L190:

/*        Special code for 10 x 10 Householder */

#line 388 "slarfx.f"
	v1 = v[1];
#line 389 "slarfx.f"
	t1 = *tau * v1;
#line 390 "slarfx.f"
	v2 = v[2];
#line 391 "slarfx.f"
	t2 = *tau * v2;
#line 392 "slarfx.f"
	v3 = v[3];
#line 393 "slarfx.f"
	t3 = *tau * v3;
#line 394 "slarfx.f"
	v4 = v[4];
#line 395 "slarfx.f"
	t4 = *tau * v4;
#line 396 "slarfx.f"
	v5 = v[5];
#line 397 "slarfx.f"
	t5 = *tau * v5;
#line 398 "slarfx.f"
	v6 = v[6];
#line 399 "slarfx.f"
	t6 = *tau * v6;
#line 400 "slarfx.f"
	v7 = v[7];
#line 401 "slarfx.f"
	t7 = *tau * v7;
#line 402 "slarfx.f"
	v8 = v[8];
#line 403 "slarfx.f"
	t8 = *tau * v8;
#line 404 "slarfx.f"
	v9 = v[9];
#line 405 "slarfx.f"
	t9 = *tau * v9;
#line 406 "slarfx.f"
	v10 = v[10];
#line 407 "slarfx.f"
	t10 = *tau * v10;
#line 408 "slarfx.f"
	i__1 = *n;
#line 408 "slarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 409 "slarfx.f"
	    sum = v1 * c__[j * c_dim1 + 1] + v2 * c__[j * c_dim1 + 2] + v3 * 
		    c__[j * c_dim1 + 3] + v4 * c__[j * c_dim1 + 4] + v5 * c__[
		    j * c_dim1 + 5] + v6 * c__[j * c_dim1 + 6] + v7 * c__[j * 
		    c_dim1 + 7] + v8 * c__[j * c_dim1 + 8] + v9 * c__[j * 
		    c_dim1 + 9] + v10 * c__[j * c_dim1 + 10];
#line 413 "slarfx.f"
	    c__[j * c_dim1 + 1] -= sum * t1;
#line 414 "slarfx.f"
	    c__[j * c_dim1 + 2] -= sum * t2;
#line 415 "slarfx.f"
	    c__[j * c_dim1 + 3] -= sum * t3;
#line 416 "slarfx.f"
	    c__[j * c_dim1 + 4] -= sum * t4;
#line 417 "slarfx.f"
	    c__[j * c_dim1 + 5] -= sum * t5;
#line 418 "slarfx.f"
	    c__[j * c_dim1 + 6] -= sum * t6;
#line 419 "slarfx.f"
	    c__[j * c_dim1 + 7] -= sum * t7;
#line 420 "slarfx.f"
	    c__[j * c_dim1 + 8] -= sum * t8;
#line 421 "slarfx.f"
	    c__[j * c_dim1 + 9] -= sum * t9;
#line 422 "slarfx.f"
	    c__[j * c_dim1 + 10] -= sum * t10;
#line 423 "slarfx.f"
/* L200: */
#line 423 "slarfx.f"
	}
#line 424 "slarfx.f"
	goto L410;
#line 425 "slarfx.f"
    } else {

/*        Form  C * H, where H has order n. */

#line 429 "slarfx.f"
	switch (*n) {
#line 429 "slarfx.f"
	    case 1:  goto L210;
#line 429 "slarfx.f"
	    case 2:  goto L230;
#line 429 "slarfx.f"
	    case 3:  goto L250;
#line 429 "slarfx.f"
	    case 4:  goto L270;
#line 429 "slarfx.f"
	    case 5:  goto L290;
#line 429 "slarfx.f"
	    case 6:  goto L310;
#line 429 "slarfx.f"
	    case 7:  goto L330;
#line 429 "slarfx.f"
	    case 8:  goto L350;
#line 429 "slarfx.f"
	    case 9:  goto L370;
#line 429 "slarfx.f"
	    case 10:  goto L390;
#line 429 "slarfx.f"
	}

/*        Code for general N */

#line 434 "slarfx.f"
	slarf_(side, m, n, &v[1], &c__1, tau, &c__[c_offset], ldc, &work[1], (
		ftnlen)1);
#line 435 "slarfx.f"
	goto L410;
#line 436 "slarfx.f"
L210:

/*        Special code for 1 x 1 Householder */

#line 440 "slarfx.f"
	t1 = 1. - *tau * v[1] * v[1];
#line 441 "slarfx.f"
	i__1 = *m;
#line 441 "slarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 442 "slarfx.f"
	    c__[j + c_dim1] = t1 * c__[j + c_dim1];
#line 443 "slarfx.f"
/* L220: */
#line 443 "slarfx.f"
	}
#line 444 "slarfx.f"
	goto L410;
#line 445 "slarfx.f"
L230:

/*        Special code for 2 x 2 Householder */

#line 449 "slarfx.f"
	v1 = v[1];
#line 450 "slarfx.f"
	t1 = *tau * v1;
#line 451 "slarfx.f"
	v2 = v[2];
#line 452 "slarfx.f"
	t2 = *tau * v2;
#line 453 "slarfx.f"
	i__1 = *m;
#line 453 "slarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 454 "slarfx.f"
	    sum = v1 * c__[j + c_dim1] + v2 * c__[j + (c_dim1 << 1)];
#line 455 "slarfx.f"
	    c__[j + c_dim1] -= sum * t1;
#line 456 "slarfx.f"
	    c__[j + (c_dim1 << 1)] -= sum * t2;
#line 457 "slarfx.f"
/* L240: */
#line 457 "slarfx.f"
	}
#line 458 "slarfx.f"
	goto L410;
#line 459 "slarfx.f"
L250:

/*        Special code for 3 x 3 Householder */

#line 463 "slarfx.f"
	v1 = v[1];
#line 464 "slarfx.f"
	t1 = *tau * v1;
#line 465 "slarfx.f"
	v2 = v[2];
#line 466 "slarfx.f"
	t2 = *tau * v2;
#line 467 "slarfx.f"
	v3 = v[3];
#line 468 "slarfx.f"
	t3 = *tau * v3;
#line 469 "slarfx.f"
	i__1 = *m;
#line 469 "slarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 470 "slarfx.f"
	    sum = v1 * c__[j + c_dim1] + v2 * c__[j + (c_dim1 << 1)] + v3 * 
		    c__[j + c_dim1 * 3];
#line 471 "slarfx.f"
	    c__[j + c_dim1] -= sum * t1;
#line 472 "slarfx.f"
	    c__[j + (c_dim1 << 1)] -= sum * t2;
#line 473 "slarfx.f"
	    c__[j + c_dim1 * 3] -= sum * t3;
#line 474 "slarfx.f"
/* L260: */
#line 474 "slarfx.f"
	}
#line 475 "slarfx.f"
	goto L410;
#line 476 "slarfx.f"
L270:

/*        Special code for 4 x 4 Householder */

#line 480 "slarfx.f"
	v1 = v[1];
#line 481 "slarfx.f"
	t1 = *tau * v1;
#line 482 "slarfx.f"
	v2 = v[2];
#line 483 "slarfx.f"
	t2 = *tau * v2;
#line 484 "slarfx.f"
	v3 = v[3];
#line 485 "slarfx.f"
	t3 = *tau * v3;
#line 486 "slarfx.f"
	v4 = v[4];
#line 487 "slarfx.f"
	t4 = *tau * v4;
#line 488 "slarfx.f"
	i__1 = *m;
#line 488 "slarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 489 "slarfx.f"
	    sum = v1 * c__[j + c_dim1] + v2 * c__[j + (c_dim1 << 1)] + v3 * 
		    c__[j + c_dim1 * 3] + v4 * c__[j + (c_dim1 << 2)];
#line 491 "slarfx.f"
	    c__[j + c_dim1] -= sum * t1;
#line 492 "slarfx.f"
	    c__[j + (c_dim1 << 1)] -= sum * t2;
#line 493 "slarfx.f"
	    c__[j + c_dim1 * 3] -= sum * t3;
#line 494 "slarfx.f"
	    c__[j + (c_dim1 << 2)] -= sum * t4;
#line 495 "slarfx.f"
/* L280: */
#line 495 "slarfx.f"
	}
#line 496 "slarfx.f"
	goto L410;
#line 497 "slarfx.f"
L290:

/*        Special code for 5 x 5 Householder */

#line 501 "slarfx.f"
	v1 = v[1];
#line 502 "slarfx.f"
	t1 = *tau * v1;
#line 503 "slarfx.f"
	v2 = v[2];
#line 504 "slarfx.f"
	t2 = *tau * v2;
#line 505 "slarfx.f"
	v3 = v[3];
#line 506 "slarfx.f"
	t3 = *tau * v3;
#line 507 "slarfx.f"
	v4 = v[4];
#line 508 "slarfx.f"
	t4 = *tau * v4;
#line 509 "slarfx.f"
	v5 = v[5];
#line 510 "slarfx.f"
	t5 = *tau * v5;
#line 511 "slarfx.f"
	i__1 = *m;
#line 511 "slarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 512 "slarfx.f"
	    sum = v1 * c__[j + c_dim1] + v2 * c__[j + (c_dim1 << 1)] + v3 * 
		    c__[j + c_dim1 * 3] + v4 * c__[j + (c_dim1 << 2)] + v5 * 
		    c__[j + c_dim1 * 5];
#line 514 "slarfx.f"
	    c__[j + c_dim1] -= sum * t1;
#line 515 "slarfx.f"
	    c__[j + (c_dim1 << 1)] -= sum * t2;
#line 516 "slarfx.f"
	    c__[j + c_dim1 * 3] -= sum * t3;
#line 517 "slarfx.f"
	    c__[j + (c_dim1 << 2)] -= sum * t4;
#line 518 "slarfx.f"
	    c__[j + c_dim1 * 5] -= sum * t5;
#line 519 "slarfx.f"
/* L300: */
#line 519 "slarfx.f"
	}
#line 520 "slarfx.f"
	goto L410;
#line 521 "slarfx.f"
L310:

/*        Special code for 6 x 6 Householder */

#line 525 "slarfx.f"
	v1 = v[1];
#line 526 "slarfx.f"
	t1 = *tau * v1;
#line 527 "slarfx.f"
	v2 = v[2];
#line 528 "slarfx.f"
	t2 = *tau * v2;
#line 529 "slarfx.f"
	v3 = v[3];
#line 530 "slarfx.f"
	t3 = *tau * v3;
#line 531 "slarfx.f"
	v4 = v[4];
#line 532 "slarfx.f"
	t4 = *tau * v4;
#line 533 "slarfx.f"
	v5 = v[5];
#line 534 "slarfx.f"
	t5 = *tau * v5;
#line 535 "slarfx.f"
	v6 = v[6];
#line 536 "slarfx.f"
	t6 = *tau * v6;
#line 537 "slarfx.f"
	i__1 = *m;
#line 537 "slarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 538 "slarfx.f"
	    sum = v1 * c__[j + c_dim1] + v2 * c__[j + (c_dim1 << 1)] + v3 * 
		    c__[j + c_dim1 * 3] + v4 * c__[j + (c_dim1 << 2)] + v5 * 
		    c__[j + c_dim1 * 5] + v6 * c__[j + c_dim1 * 6];
#line 540 "slarfx.f"
	    c__[j + c_dim1] -= sum * t1;
#line 541 "slarfx.f"
	    c__[j + (c_dim1 << 1)] -= sum * t2;
#line 542 "slarfx.f"
	    c__[j + c_dim1 * 3] -= sum * t3;
#line 543 "slarfx.f"
	    c__[j + (c_dim1 << 2)] -= sum * t4;
#line 544 "slarfx.f"
	    c__[j + c_dim1 * 5] -= sum * t5;
#line 545 "slarfx.f"
	    c__[j + c_dim1 * 6] -= sum * t6;
#line 546 "slarfx.f"
/* L320: */
#line 546 "slarfx.f"
	}
#line 547 "slarfx.f"
	goto L410;
#line 548 "slarfx.f"
L330:

/*        Special code for 7 x 7 Householder */

#line 552 "slarfx.f"
	v1 = v[1];
#line 553 "slarfx.f"
	t1 = *tau * v1;
#line 554 "slarfx.f"
	v2 = v[2];
#line 555 "slarfx.f"
	t2 = *tau * v2;
#line 556 "slarfx.f"
	v3 = v[3];
#line 557 "slarfx.f"
	t3 = *tau * v3;
#line 558 "slarfx.f"
	v4 = v[4];
#line 559 "slarfx.f"
	t4 = *tau * v4;
#line 560 "slarfx.f"
	v5 = v[5];
#line 561 "slarfx.f"
	t5 = *tau * v5;
#line 562 "slarfx.f"
	v6 = v[6];
#line 563 "slarfx.f"
	t6 = *tau * v6;
#line 564 "slarfx.f"
	v7 = v[7];
#line 565 "slarfx.f"
	t7 = *tau * v7;
#line 566 "slarfx.f"
	i__1 = *m;
#line 566 "slarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 567 "slarfx.f"
	    sum = v1 * c__[j + c_dim1] + v2 * c__[j + (c_dim1 << 1)] + v3 * 
		    c__[j + c_dim1 * 3] + v4 * c__[j + (c_dim1 << 2)] + v5 * 
		    c__[j + c_dim1 * 5] + v6 * c__[j + c_dim1 * 6] + v7 * c__[
		    j + c_dim1 * 7];
#line 570 "slarfx.f"
	    c__[j + c_dim1] -= sum * t1;
#line 571 "slarfx.f"
	    c__[j + (c_dim1 << 1)] -= sum * t2;
#line 572 "slarfx.f"
	    c__[j + c_dim1 * 3] -= sum * t3;
#line 573 "slarfx.f"
	    c__[j + (c_dim1 << 2)] -= sum * t4;
#line 574 "slarfx.f"
	    c__[j + c_dim1 * 5] -= sum * t5;
#line 575 "slarfx.f"
	    c__[j + c_dim1 * 6] -= sum * t6;
#line 576 "slarfx.f"
	    c__[j + c_dim1 * 7] -= sum * t7;
#line 577 "slarfx.f"
/* L340: */
#line 577 "slarfx.f"
	}
#line 578 "slarfx.f"
	goto L410;
#line 579 "slarfx.f"
L350:

/*        Special code for 8 x 8 Householder */

#line 583 "slarfx.f"
	v1 = v[1];
#line 584 "slarfx.f"
	t1 = *tau * v1;
#line 585 "slarfx.f"
	v2 = v[2];
#line 586 "slarfx.f"
	t2 = *tau * v2;
#line 587 "slarfx.f"
	v3 = v[3];
#line 588 "slarfx.f"
	t3 = *tau * v3;
#line 589 "slarfx.f"
	v4 = v[4];
#line 590 "slarfx.f"
	t4 = *tau * v4;
#line 591 "slarfx.f"
	v5 = v[5];
#line 592 "slarfx.f"
	t5 = *tau * v5;
#line 593 "slarfx.f"
	v6 = v[6];
#line 594 "slarfx.f"
	t6 = *tau * v6;
#line 595 "slarfx.f"
	v7 = v[7];
#line 596 "slarfx.f"
	t7 = *tau * v7;
#line 597 "slarfx.f"
	v8 = v[8];
#line 598 "slarfx.f"
	t8 = *tau * v8;
#line 599 "slarfx.f"
	i__1 = *m;
#line 599 "slarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 600 "slarfx.f"
	    sum = v1 * c__[j + c_dim1] + v2 * c__[j + (c_dim1 << 1)] + v3 * 
		    c__[j + c_dim1 * 3] + v4 * c__[j + (c_dim1 << 2)] + v5 * 
		    c__[j + c_dim1 * 5] + v6 * c__[j + c_dim1 * 6] + v7 * c__[
		    j + c_dim1 * 7] + v8 * c__[j + (c_dim1 << 3)];
#line 603 "slarfx.f"
	    c__[j + c_dim1] -= sum * t1;
#line 604 "slarfx.f"
	    c__[j + (c_dim1 << 1)] -= sum * t2;
#line 605 "slarfx.f"
	    c__[j + c_dim1 * 3] -= sum * t3;
#line 606 "slarfx.f"
	    c__[j + (c_dim1 << 2)] -= sum * t4;
#line 607 "slarfx.f"
	    c__[j + c_dim1 * 5] -= sum * t5;
#line 608 "slarfx.f"
	    c__[j + c_dim1 * 6] -= sum * t6;
#line 609 "slarfx.f"
	    c__[j + c_dim1 * 7] -= sum * t7;
#line 610 "slarfx.f"
	    c__[j + (c_dim1 << 3)] -= sum * t8;
#line 611 "slarfx.f"
/* L360: */
#line 611 "slarfx.f"
	}
#line 612 "slarfx.f"
	goto L410;
#line 613 "slarfx.f"
L370:

/*        Special code for 9 x 9 Householder */

#line 617 "slarfx.f"
	v1 = v[1];
#line 618 "slarfx.f"
	t1 = *tau * v1;
#line 619 "slarfx.f"
	v2 = v[2];
#line 620 "slarfx.f"
	t2 = *tau * v2;
#line 621 "slarfx.f"
	v3 = v[3];
#line 622 "slarfx.f"
	t3 = *tau * v3;
#line 623 "slarfx.f"
	v4 = v[4];
#line 624 "slarfx.f"
	t4 = *tau * v4;
#line 625 "slarfx.f"
	v5 = v[5];
#line 626 "slarfx.f"
	t5 = *tau * v5;
#line 627 "slarfx.f"
	v6 = v[6];
#line 628 "slarfx.f"
	t6 = *tau * v6;
#line 629 "slarfx.f"
	v7 = v[7];
#line 630 "slarfx.f"
	t7 = *tau * v7;
#line 631 "slarfx.f"
	v8 = v[8];
#line 632 "slarfx.f"
	t8 = *tau * v8;
#line 633 "slarfx.f"
	v9 = v[9];
#line 634 "slarfx.f"
	t9 = *tau * v9;
#line 635 "slarfx.f"
	i__1 = *m;
#line 635 "slarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 636 "slarfx.f"
	    sum = v1 * c__[j + c_dim1] + v2 * c__[j + (c_dim1 << 1)] + v3 * 
		    c__[j + c_dim1 * 3] + v4 * c__[j + (c_dim1 << 2)] + v5 * 
		    c__[j + c_dim1 * 5] + v6 * c__[j + c_dim1 * 6] + v7 * c__[
		    j + c_dim1 * 7] + v8 * c__[j + (c_dim1 << 3)] + v9 * c__[
		    j + c_dim1 * 9];
#line 639 "slarfx.f"
	    c__[j + c_dim1] -= sum * t1;
#line 640 "slarfx.f"
	    c__[j + (c_dim1 << 1)] -= sum * t2;
#line 641 "slarfx.f"
	    c__[j + c_dim1 * 3] -= sum * t3;
#line 642 "slarfx.f"
	    c__[j + (c_dim1 << 2)] -= sum * t4;
#line 643 "slarfx.f"
	    c__[j + c_dim1 * 5] -= sum * t5;
#line 644 "slarfx.f"
	    c__[j + c_dim1 * 6] -= sum * t6;
#line 645 "slarfx.f"
	    c__[j + c_dim1 * 7] -= sum * t7;
#line 646 "slarfx.f"
	    c__[j + (c_dim1 << 3)] -= sum * t8;
#line 647 "slarfx.f"
	    c__[j + c_dim1 * 9] -= sum * t9;
#line 648 "slarfx.f"
/* L380: */
#line 648 "slarfx.f"
	}
#line 649 "slarfx.f"
	goto L410;
#line 650 "slarfx.f"
L390:

/*        Special code for 10 x 10 Householder */

#line 654 "slarfx.f"
	v1 = v[1];
#line 655 "slarfx.f"
	t1 = *tau * v1;
#line 656 "slarfx.f"
	v2 = v[2];
#line 657 "slarfx.f"
	t2 = *tau * v2;
#line 658 "slarfx.f"
	v3 = v[3];
#line 659 "slarfx.f"
	t3 = *tau * v3;
#line 660 "slarfx.f"
	v4 = v[4];
#line 661 "slarfx.f"
	t4 = *tau * v4;
#line 662 "slarfx.f"
	v5 = v[5];
#line 663 "slarfx.f"
	t5 = *tau * v5;
#line 664 "slarfx.f"
	v6 = v[6];
#line 665 "slarfx.f"
	t6 = *tau * v6;
#line 666 "slarfx.f"
	v7 = v[7];
#line 667 "slarfx.f"
	t7 = *tau * v7;
#line 668 "slarfx.f"
	v8 = v[8];
#line 669 "slarfx.f"
	t8 = *tau * v8;
#line 670 "slarfx.f"
	v9 = v[9];
#line 671 "slarfx.f"
	t9 = *tau * v9;
#line 672 "slarfx.f"
	v10 = v[10];
#line 673 "slarfx.f"
	t10 = *tau * v10;
#line 674 "slarfx.f"
	i__1 = *m;
#line 674 "slarfx.f"
	for (j = 1; j <= i__1; ++j) {
#line 675 "slarfx.f"
	    sum = v1 * c__[j + c_dim1] + v2 * c__[j + (c_dim1 << 1)] + v3 * 
		    c__[j + c_dim1 * 3] + v4 * c__[j + (c_dim1 << 2)] + v5 * 
		    c__[j + c_dim1 * 5] + v6 * c__[j + c_dim1 * 6] + v7 * c__[
		    j + c_dim1 * 7] + v8 * c__[j + (c_dim1 << 3)] + v9 * c__[
		    j + c_dim1 * 9] + v10 * c__[j + c_dim1 * 10];
#line 679 "slarfx.f"
	    c__[j + c_dim1] -= sum * t1;
#line 680 "slarfx.f"
	    c__[j + (c_dim1 << 1)] -= sum * t2;
#line 681 "slarfx.f"
	    c__[j + c_dim1 * 3] -= sum * t3;
#line 682 "slarfx.f"
	    c__[j + (c_dim1 << 2)] -= sum * t4;
#line 683 "slarfx.f"
	    c__[j + c_dim1 * 5] -= sum * t5;
#line 684 "slarfx.f"
	    c__[j + c_dim1 * 6] -= sum * t6;
#line 685 "slarfx.f"
	    c__[j + c_dim1 * 7] -= sum * t7;
#line 686 "slarfx.f"
	    c__[j + (c_dim1 << 3)] -= sum * t8;
#line 687 "slarfx.f"
	    c__[j + c_dim1 * 9] -= sum * t9;
#line 688 "slarfx.f"
	    c__[j + c_dim1 * 10] -= sum * t10;
#line 689 "slarfx.f"
/* L400: */
#line 689 "slarfx.f"
	}
#line 690 "slarfx.f"
	goto L410;
#line 691 "slarfx.f"
    }
#line 692 "slarfx.f"
L410:
#line 692 "slarfx.f"
    return 0;

/*     End of SLARFX */

} /* slarfx_ */

