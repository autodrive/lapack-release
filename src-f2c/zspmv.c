#line 1 "zspmv.f"
/* zspmv.f -- translated by f2c (version 20100827).
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

#line 1 "zspmv.f"
/* > \brief \b ZSPMV computes a matrix-vector product for complex vectors using a complex symmetric packed mat
rix */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZSPMV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zspmv.f
"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zspmv.f
"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zspmv.f
"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZSPMV( UPLO, N, ALPHA, AP, X, INCX, BETA, Y, INCY ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INCX, INCY, N */
/*       COMPLEX*16         ALPHA, BETA */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         AP( * ), X( * ), Y( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZSPMV  performs the matrix-vector operation */
/* > */
/* >    y := alpha*A*x + beta*y, */
/* > */
/* > where alpha and beta are scalars, x and y are n element vectors and */
/* > A is an n by n symmetric matrix, supplied in packed form. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >           On entry, UPLO specifies whether the upper or lower */
/* >           triangular part of the matrix A is supplied in the packed */
/* >           array AP as follows: */
/* > */
/* >              UPLO = 'U' or 'u'   The upper triangular part of A is */
/* >                                  supplied in AP. */
/* > */
/* >              UPLO = 'L' or 'l'   The lower triangular part of A is */
/* >                                  supplied in AP. */
/* > */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >           On entry, N specifies the order of the matrix A. */
/* >           N must be at least zero. */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] ALPHA */
/* > \verbatim */
/* >          ALPHA is COMPLEX*16 */
/* >           On entry, ALPHA specifies the scalar alpha. */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] AP */
/* > \verbatim */
/* >          AP is COMPLEX*16 array, dimension at least */
/* >           ( ( N*( N + 1 ) )/2 ). */
/* >           Before entry, with UPLO = 'U' or 'u', the array AP must */
/* >           contain the upper triangular part of the symmetric matrix */
/* >           packed sequentially, column by column, so that AP( 1 ) */
/* >           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 ) */
/* >           and a( 2, 2 ) respectively, and so on. */
/* >           Before entry, with UPLO = 'L' or 'l', the array AP must */
/* >           contain the lower triangular part of the symmetric matrix */
/* >           packed sequentially, column by column, so that AP( 1 ) */
/* >           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 ) */
/* >           and a( 3, 1 ) respectively, and so on. */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] X */
/* > \verbatim */
/* >          X is COMPLEX*16 array, dimension at least */
/* >           ( 1 + ( N - 1 )*abs( INCX ) ). */
/* >           Before entry, the incremented array X must contain the N- */
/* >           element vector x. */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* >          INCX is INTEGER */
/* >           On entry, INCX specifies the increment for the elements of */
/* >           X. INCX must not be zero. */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] BETA */
/* > \verbatim */
/* >          BETA is COMPLEX*16 */
/* >           On entry, BETA specifies the scalar beta. When BETA is */
/* >           supplied as zero then Y need not be set on input. */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Y */
/* > \verbatim */
/* >          Y is COMPLEX*16 array, dimension at least */
/* >           ( 1 + ( N - 1 )*abs( INCY ) ). */
/* >           Before entry, the incremented array Y must contain the n */
/* >           element vector y. On exit, Y is overwritten by the updated */
/* >           vector y. */
/* > \endverbatim */
/* > */
/* > \param[in] INCY */
/* > \verbatim */
/* >          INCY is INTEGER */
/* >           On entry, INCY specifies the increment for the elements of */
/* >           Y. INCY must not be zero. */
/* >           Unchanged on exit. */
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
/* Subroutine */ int zspmv_(char *uplo, integer *n, doublecomplex *alpha, 
	doublecomplex *ap, doublecomplex *x, integer *incx, doublecomplex *
	beta, doublecomplex *y, integer *incy, ftnlen uplo_len)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5;
    doublecomplex z__1, z__2, z__3, z__4;

    /* Local variables */
    static integer i__, j, k, kk, ix, iy, jx, jy, kx, ky, info;
    static doublecomplex temp1, temp2;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);


/*  -- LAPACK auxiliary routine (version 3.4.2) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     September 2012 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/* ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 191 "zspmv.f"
    /* Parameter adjustments */
#line 191 "zspmv.f"
    --y;
#line 191 "zspmv.f"
    --x;
#line 191 "zspmv.f"
    --ap;
#line 191 "zspmv.f"

#line 191 "zspmv.f"
    /* Function Body */
#line 191 "zspmv.f"
    info = 0;
#line 192 "zspmv.f"
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
#line 193 "zspmv.f"
	info = 1;
#line 194 "zspmv.f"
    } else if (*n < 0) {
#line 195 "zspmv.f"
	info = 2;
#line 196 "zspmv.f"
    } else if (*incx == 0) {
#line 197 "zspmv.f"
	info = 6;
#line 198 "zspmv.f"
    } else if (*incy == 0) {
#line 199 "zspmv.f"
	info = 9;
#line 200 "zspmv.f"
    }
#line 201 "zspmv.f"
    if (info != 0) {
#line 202 "zspmv.f"
	xerbla_("ZSPMV ", &info, (ftnlen)6);
#line 203 "zspmv.f"
	return 0;
#line 204 "zspmv.f"
    }

/*     Quick return if possible. */

#line 208 "zspmv.f"
    if (*n == 0 || alpha->r == 0. && alpha->i == 0. && (beta->r == 1. && 
	    beta->i == 0.)) {
#line 208 "zspmv.f"
	return 0;
#line 208 "zspmv.f"
    }

/*     Set up the start points in  X  and  Y. */

#line 213 "zspmv.f"
    if (*incx > 0) {
#line 214 "zspmv.f"
	kx = 1;
#line 215 "zspmv.f"
    } else {
#line 216 "zspmv.f"
	kx = 1 - (*n - 1) * *incx;
#line 217 "zspmv.f"
    }
#line 218 "zspmv.f"
    if (*incy > 0) {
#line 219 "zspmv.f"
	ky = 1;
#line 220 "zspmv.f"
    } else {
#line 221 "zspmv.f"
	ky = 1 - (*n - 1) * *incy;
#line 222 "zspmv.f"
    }

/*     Start the operations. In this version the elements of the array AP */
/*     are accessed sequentially with one pass through AP. */

/*     First form  y := beta*y. */

#line 229 "zspmv.f"
    if (beta->r != 1. || beta->i != 0.) {
#line 230 "zspmv.f"
	if (*incy == 1) {
#line 231 "zspmv.f"
	    if (beta->r == 0. && beta->i == 0.) {
#line 232 "zspmv.f"
		i__1 = *n;
#line 232 "zspmv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 233 "zspmv.f"
		    i__2 = i__;
#line 233 "zspmv.f"
		    y[i__2].r = 0., y[i__2].i = 0.;
#line 234 "zspmv.f"
/* L10: */
#line 234 "zspmv.f"
		}
#line 235 "zspmv.f"
	    } else {
#line 236 "zspmv.f"
		i__1 = *n;
#line 236 "zspmv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 237 "zspmv.f"
		    i__2 = i__;
#line 237 "zspmv.f"
		    i__3 = i__;
#line 237 "zspmv.f"
		    z__1.r = beta->r * y[i__3].r - beta->i * y[i__3].i, 
			    z__1.i = beta->r * y[i__3].i + beta->i * y[i__3]
			    .r;
#line 237 "zspmv.f"
		    y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 238 "zspmv.f"
/* L20: */
#line 238 "zspmv.f"
		}
#line 239 "zspmv.f"
	    }
#line 240 "zspmv.f"
	} else {
#line 241 "zspmv.f"
	    iy = ky;
#line 242 "zspmv.f"
	    if (beta->r == 0. && beta->i == 0.) {
#line 243 "zspmv.f"
		i__1 = *n;
#line 243 "zspmv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 244 "zspmv.f"
		    i__2 = iy;
#line 244 "zspmv.f"
		    y[i__2].r = 0., y[i__2].i = 0.;
#line 245 "zspmv.f"
		    iy += *incy;
#line 246 "zspmv.f"
/* L30: */
#line 246 "zspmv.f"
		}
#line 247 "zspmv.f"
	    } else {
#line 248 "zspmv.f"
		i__1 = *n;
#line 248 "zspmv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 249 "zspmv.f"
		    i__2 = iy;
#line 249 "zspmv.f"
		    i__3 = iy;
#line 249 "zspmv.f"
		    z__1.r = beta->r * y[i__3].r - beta->i * y[i__3].i, 
			    z__1.i = beta->r * y[i__3].i + beta->i * y[i__3]
			    .r;
#line 249 "zspmv.f"
		    y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 250 "zspmv.f"
		    iy += *incy;
#line 251 "zspmv.f"
/* L40: */
#line 251 "zspmv.f"
		}
#line 252 "zspmv.f"
	    }
#line 253 "zspmv.f"
	}
#line 254 "zspmv.f"
    }
#line 255 "zspmv.f"
    if (alpha->r == 0. && alpha->i == 0.) {
#line 255 "zspmv.f"
	return 0;
#line 255 "zspmv.f"
    }
#line 257 "zspmv.f"
    kk = 1;
#line 258 "zspmv.f"
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {

/*        Form  y  when AP contains the upper triangle. */

#line 262 "zspmv.f"
	if (*incx == 1 && *incy == 1) {
#line 263 "zspmv.f"
	    i__1 = *n;
#line 263 "zspmv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 264 "zspmv.f"
		i__2 = j;
#line 264 "zspmv.f"
		z__1.r = alpha->r * x[i__2].r - alpha->i * x[i__2].i, z__1.i =
			 alpha->r * x[i__2].i + alpha->i * x[i__2].r;
#line 264 "zspmv.f"
		temp1.r = z__1.r, temp1.i = z__1.i;
#line 265 "zspmv.f"
		temp2.r = 0., temp2.i = 0.;
#line 266 "zspmv.f"
		k = kk;
#line 267 "zspmv.f"
		i__2 = j - 1;
#line 267 "zspmv.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 268 "zspmv.f"
		    i__3 = i__;
#line 268 "zspmv.f"
		    i__4 = i__;
#line 268 "zspmv.f"
		    i__5 = k;
#line 268 "zspmv.f"
		    z__2.r = temp1.r * ap[i__5].r - temp1.i * ap[i__5].i, 
			    z__2.i = temp1.r * ap[i__5].i + temp1.i * ap[i__5]
			    .r;
#line 268 "zspmv.f"
		    z__1.r = y[i__4].r + z__2.r, z__1.i = y[i__4].i + z__2.i;
#line 268 "zspmv.f"
		    y[i__3].r = z__1.r, y[i__3].i = z__1.i;
#line 269 "zspmv.f"
		    i__3 = k;
#line 269 "zspmv.f"
		    i__4 = i__;
#line 269 "zspmv.f"
		    z__2.r = ap[i__3].r * x[i__4].r - ap[i__3].i * x[i__4].i, 
			    z__2.i = ap[i__3].r * x[i__4].i + ap[i__3].i * x[
			    i__4].r;
#line 269 "zspmv.f"
		    z__1.r = temp2.r + z__2.r, z__1.i = temp2.i + z__2.i;
#line 269 "zspmv.f"
		    temp2.r = z__1.r, temp2.i = z__1.i;
#line 270 "zspmv.f"
		    ++k;
#line 271 "zspmv.f"
/* L50: */
#line 271 "zspmv.f"
		}
#line 272 "zspmv.f"
		i__2 = j;
#line 272 "zspmv.f"
		i__3 = j;
#line 272 "zspmv.f"
		i__4 = kk + j - 1;
#line 272 "zspmv.f"
		z__3.r = temp1.r * ap[i__4].r - temp1.i * ap[i__4].i, z__3.i =
			 temp1.r * ap[i__4].i + temp1.i * ap[i__4].r;
#line 272 "zspmv.f"
		z__2.r = y[i__3].r + z__3.r, z__2.i = y[i__3].i + z__3.i;
#line 272 "zspmv.f"
		z__4.r = alpha->r * temp2.r - alpha->i * temp2.i, z__4.i = 
			alpha->r * temp2.i + alpha->i * temp2.r;
#line 272 "zspmv.f"
		z__1.r = z__2.r + z__4.r, z__1.i = z__2.i + z__4.i;
#line 272 "zspmv.f"
		y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 273 "zspmv.f"
		kk += j;
#line 274 "zspmv.f"
/* L60: */
#line 274 "zspmv.f"
	    }
#line 275 "zspmv.f"
	} else {
#line 276 "zspmv.f"
	    jx = kx;
#line 277 "zspmv.f"
	    jy = ky;
#line 278 "zspmv.f"
	    i__1 = *n;
#line 278 "zspmv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 279 "zspmv.f"
		i__2 = jx;
#line 279 "zspmv.f"
		z__1.r = alpha->r * x[i__2].r - alpha->i * x[i__2].i, z__1.i =
			 alpha->r * x[i__2].i + alpha->i * x[i__2].r;
#line 279 "zspmv.f"
		temp1.r = z__1.r, temp1.i = z__1.i;
#line 280 "zspmv.f"
		temp2.r = 0., temp2.i = 0.;
#line 281 "zspmv.f"
		ix = kx;
#line 282 "zspmv.f"
		iy = ky;
#line 283 "zspmv.f"
		i__2 = kk + j - 2;
#line 283 "zspmv.f"
		for (k = kk; k <= i__2; ++k) {
#line 284 "zspmv.f"
		    i__3 = iy;
#line 284 "zspmv.f"
		    i__4 = iy;
#line 284 "zspmv.f"
		    i__5 = k;
#line 284 "zspmv.f"
		    z__2.r = temp1.r * ap[i__5].r - temp1.i * ap[i__5].i, 
			    z__2.i = temp1.r * ap[i__5].i + temp1.i * ap[i__5]
			    .r;
#line 284 "zspmv.f"
		    z__1.r = y[i__4].r + z__2.r, z__1.i = y[i__4].i + z__2.i;
#line 284 "zspmv.f"
		    y[i__3].r = z__1.r, y[i__3].i = z__1.i;
#line 285 "zspmv.f"
		    i__3 = k;
#line 285 "zspmv.f"
		    i__4 = ix;
#line 285 "zspmv.f"
		    z__2.r = ap[i__3].r * x[i__4].r - ap[i__3].i * x[i__4].i, 
			    z__2.i = ap[i__3].r * x[i__4].i + ap[i__3].i * x[
			    i__4].r;
#line 285 "zspmv.f"
		    z__1.r = temp2.r + z__2.r, z__1.i = temp2.i + z__2.i;
#line 285 "zspmv.f"
		    temp2.r = z__1.r, temp2.i = z__1.i;
#line 286 "zspmv.f"
		    ix += *incx;
#line 287 "zspmv.f"
		    iy += *incy;
#line 288 "zspmv.f"
/* L70: */
#line 288 "zspmv.f"
		}
#line 289 "zspmv.f"
		i__2 = jy;
#line 289 "zspmv.f"
		i__3 = jy;
#line 289 "zspmv.f"
		i__4 = kk + j - 1;
#line 289 "zspmv.f"
		z__3.r = temp1.r * ap[i__4].r - temp1.i * ap[i__4].i, z__3.i =
			 temp1.r * ap[i__4].i + temp1.i * ap[i__4].r;
#line 289 "zspmv.f"
		z__2.r = y[i__3].r + z__3.r, z__2.i = y[i__3].i + z__3.i;
#line 289 "zspmv.f"
		z__4.r = alpha->r * temp2.r - alpha->i * temp2.i, z__4.i = 
			alpha->r * temp2.i + alpha->i * temp2.r;
#line 289 "zspmv.f"
		z__1.r = z__2.r + z__4.r, z__1.i = z__2.i + z__4.i;
#line 289 "zspmv.f"
		y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 290 "zspmv.f"
		jx += *incx;
#line 291 "zspmv.f"
		jy += *incy;
#line 292 "zspmv.f"
		kk += j;
#line 293 "zspmv.f"
/* L80: */
#line 293 "zspmv.f"
	    }
#line 294 "zspmv.f"
	}
#line 295 "zspmv.f"
    } else {

/*        Form  y  when AP contains the lower triangle. */

#line 299 "zspmv.f"
	if (*incx == 1 && *incy == 1) {
#line 300 "zspmv.f"
	    i__1 = *n;
#line 300 "zspmv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 301 "zspmv.f"
		i__2 = j;
#line 301 "zspmv.f"
		z__1.r = alpha->r * x[i__2].r - alpha->i * x[i__2].i, z__1.i =
			 alpha->r * x[i__2].i + alpha->i * x[i__2].r;
#line 301 "zspmv.f"
		temp1.r = z__1.r, temp1.i = z__1.i;
#line 302 "zspmv.f"
		temp2.r = 0., temp2.i = 0.;
#line 303 "zspmv.f"
		i__2 = j;
#line 303 "zspmv.f"
		i__3 = j;
#line 303 "zspmv.f"
		i__4 = kk;
#line 303 "zspmv.f"
		z__2.r = temp1.r * ap[i__4].r - temp1.i * ap[i__4].i, z__2.i =
			 temp1.r * ap[i__4].i + temp1.i * ap[i__4].r;
#line 303 "zspmv.f"
		z__1.r = y[i__3].r + z__2.r, z__1.i = y[i__3].i + z__2.i;
#line 303 "zspmv.f"
		y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 304 "zspmv.f"
		k = kk + 1;
#line 305 "zspmv.f"
		i__2 = *n;
#line 305 "zspmv.f"
		for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 306 "zspmv.f"
		    i__3 = i__;
#line 306 "zspmv.f"
		    i__4 = i__;
#line 306 "zspmv.f"
		    i__5 = k;
#line 306 "zspmv.f"
		    z__2.r = temp1.r * ap[i__5].r - temp1.i * ap[i__5].i, 
			    z__2.i = temp1.r * ap[i__5].i + temp1.i * ap[i__5]
			    .r;
#line 306 "zspmv.f"
		    z__1.r = y[i__4].r + z__2.r, z__1.i = y[i__4].i + z__2.i;
#line 306 "zspmv.f"
		    y[i__3].r = z__1.r, y[i__3].i = z__1.i;
#line 307 "zspmv.f"
		    i__3 = k;
#line 307 "zspmv.f"
		    i__4 = i__;
#line 307 "zspmv.f"
		    z__2.r = ap[i__3].r * x[i__4].r - ap[i__3].i * x[i__4].i, 
			    z__2.i = ap[i__3].r * x[i__4].i + ap[i__3].i * x[
			    i__4].r;
#line 307 "zspmv.f"
		    z__1.r = temp2.r + z__2.r, z__1.i = temp2.i + z__2.i;
#line 307 "zspmv.f"
		    temp2.r = z__1.r, temp2.i = z__1.i;
#line 308 "zspmv.f"
		    ++k;
#line 309 "zspmv.f"
/* L90: */
#line 309 "zspmv.f"
		}
#line 310 "zspmv.f"
		i__2 = j;
#line 310 "zspmv.f"
		i__3 = j;
#line 310 "zspmv.f"
		z__2.r = alpha->r * temp2.r - alpha->i * temp2.i, z__2.i = 
			alpha->r * temp2.i + alpha->i * temp2.r;
#line 310 "zspmv.f"
		z__1.r = y[i__3].r + z__2.r, z__1.i = y[i__3].i + z__2.i;
#line 310 "zspmv.f"
		y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 311 "zspmv.f"
		kk += *n - j + 1;
#line 312 "zspmv.f"
/* L100: */
#line 312 "zspmv.f"
	    }
#line 313 "zspmv.f"
	} else {
#line 314 "zspmv.f"
	    jx = kx;
#line 315 "zspmv.f"
	    jy = ky;
#line 316 "zspmv.f"
	    i__1 = *n;
#line 316 "zspmv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 317 "zspmv.f"
		i__2 = jx;
#line 317 "zspmv.f"
		z__1.r = alpha->r * x[i__2].r - alpha->i * x[i__2].i, z__1.i =
			 alpha->r * x[i__2].i + alpha->i * x[i__2].r;
#line 317 "zspmv.f"
		temp1.r = z__1.r, temp1.i = z__1.i;
#line 318 "zspmv.f"
		temp2.r = 0., temp2.i = 0.;
#line 319 "zspmv.f"
		i__2 = jy;
#line 319 "zspmv.f"
		i__3 = jy;
#line 319 "zspmv.f"
		i__4 = kk;
#line 319 "zspmv.f"
		z__2.r = temp1.r * ap[i__4].r - temp1.i * ap[i__4].i, z__2.i =
			 temp1.r * ap[i__4].i + temp1.i * ap[i__4].r;
#line 319 "zspmv.f"
		z__1.r = y[i__3].r + z__2.r, z__1.i = y[i__3].i + z__2.i;
#line 319 "zspmv.f"
		y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 320 "zspmv.f"
		ix = jx;
#line 321 "zspmv.f"
		iy = jy;
#line 322 "zspmv.f"
		i__2 = kk + *n - j;
#line 322 "zspmv.f"
		for (k = kk + 1; k <= i__2; ++k) {
#line 323 "zspmv.f"
		    ix += *incx;
#line 324 "zspmv.f"
		    iy += *incy;
#line 325 "zspmv.f"
		    i__3 = iy;
#line 325 "zspmv.f"
		    i__4 = iy;
#line 325 "zspmv.f"
		    i__5 = k;
#line 325 "zspmv.f"
		    z__2.r = temp1.r * ap[i__5].r - temp1.i * ap[i__5].i, 
			    z__2.i = temp1.r * ap[i__5].i + temp1.i * ap[i__5]
			    .r;
#line 325 "zspmv.f"
		    z__1.r = y[i__4].r + z__2.r, z__1.i = y[i__4].i + z__2.i;
#line 325 "zspmv.f"
		    y[i__3].r = z__1.r, y[i__3].i = z__1.i;
#line 326 "zspmv.f"
		    i__3 = k;
#line 326 "zspmv.f"
		    i__4 = ix;
#line 326 "zspmv.f"
		    z__2.r = ap[i__3].r * x[i__4].r - ap[i__3].i * x[i__4].i, 
			    z__2.i = ap[i__3].r * x[i__4].i + ap[i__3].i * x[
			    i__4].r;
#line 326 "zspmv.f"
		    z__1.r = temp2.r + z__2.r, z__1.i = temp2.i + z__2.i;
#line 326 "zspmv.f"
		    temp2.r = z__1.r, temp2.i = z__1.i;
#line 327 "zspmv.f"
/* L110: */
#line 327 "zspmv.f"
		}
#line 328 "zspmv.f"
		i__2 = jy;
#line 328 "zspmv.f"
		i__3 = jy;
#line 328 "zspmv.f"
		z__2.r = alpha->r * temp2.r - alpha->i * temp2.i, z__2.i = 
			alpha->r * temp2.i + alpha->i * temp2.r;
#line 328 "zspmv.f"
		z__1.r = y[i__3].r + z__2.r, z__1.i = y[i__3].i + z__2.i;
#line 328 "zspmv.f"
		y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 329 "zspmv.f"
		jx += *incx;
#line 330 "zspmv.f"
		jy += *incy;
#line 331 "zspmv.f"
		kk += *n - j + 1;
#line 332 "zspmv.f"
/* L120: */
#line 332 "zspmv.f"
	    }
#line 333 "zspmv.f"
	}
#line 334 "zspmv.f"
    }

#line 336 "zspmv.f"
    return 0;

/*     End of ZSPMV */

} /* zspmv_ */

