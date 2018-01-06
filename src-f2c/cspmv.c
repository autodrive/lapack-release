#line 1 "cspmv.f"
/* cspmv.f -- translated by f2c (version 20100827).
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

#line 1 "cspmv.f"
/* > \brief \b CSPMV computes a matrix-vector product for complex vectors using a complex symmetric packed mat
rix */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CSPMV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cspmv.f
"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cspmv.f
"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cspmv.f
"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CSPMV( UPLO, N, ALPHA, AP, X, INCX, BETA, Y, INCY ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INCX, INCY, N */
/*       COMPLEX            ALPHA, BETA */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX            AP( * ), X( * ), Y( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CSPMV  performs the matrix-vector operation */
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
/* >          ALPHA is COMPLEX */
/* >           On entry, ALPHA specifies the scalar alpha. */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] AP */
/* > \verbatim */
/* >          AP is COMPLEX array, dimension at least */
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
/* >          X is COMPLEX array, dimension at least */
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
/* >          BETA is COMPLEX */
/* >           On entry, BETA specifies the scalar beta. When BETA is */
/* >           supplied as zero then Y need not be set on input. */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Y */
/* > \verbatim */
/* >          Y is COMPLEX array, dimension at least */
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

/* > \ingroup complexOTHERauxiliary */

/*  ===================================================================== */
/* Subroutine */ int cspmv_(char *uplo, integer *n, doublecomplex *alpha, 
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

#line 191 "cspmv.f"
    /* Parameter adjustments */
#line 191 "cspmv.f"
    --y;
#line 191 "cspmv.f"
    --x;
#line 191 "cspmv.f"
    --ap;
#line 191 "cspmv.f"

#line 191 "cspmv.f"
    /* Function Body */
#line 191 "cspmv.f"
    info = 0;
#line 192 "cspmv.f"
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
#line 193 "cspmv.f"
	info = 1;
#line 194 "cspmv.f"
    } else if (*n < 0) {
#line 195 "cspmv.f"
	info = 2;
#line 196 "cspmv.f"
    } else if (*incx == 0) {
#line 197 "cspmv.f"
	info = 6;
#line 198 "cspmv.f"
    } else if (*incy == 0) {
#line 199 "cspmv.f"
	info = 9;
#line 200 "cspmv.f"
    }
#line 201 "cspmv.f"
    if (info != 0) {
#line 202 "cspmv.f"
	xerbla_("CSPMV ", &info, (ftnlen)6);
#line 203 "cspmv.f"
	return 0;
#line 204 "cspmv.f"
    }

/*     Quick return if possible. */

#line 208 "cspmv.f"
    if (*n == 0 || alpha->r == 0. && alpha->i == 0. && (beta->r == 1. && 
	    beta->i == 0.)) {
#line 208 "cspmv.f"
	return 0;
#line 208 "cspmv.f"
    }

/*     Set up the start points in  X  and  Y. */

#line 213 "cspmv.f"
    if (*incx > 0) {
#line 214 "cspmv.f"
	kx = 1;
#line 215 "cspmv.f"
    } else {
#line 216 "cspmv.f"
	kx = 1 - (*n - 1) * *incx;
#line 217 "cspmv.f"
    }
#line 218 "cspmv.f"
    if (*incy > 0) {
#line 219 "cspmv.f"
	ky = 1;
#line 220 "cspmv.f"
    } else {
#line 221 "cspmv.f"
	ky = 1 - (*n - 1) * *incy;
#line 222 "cspmv.f"
    }

/*     Start the operations. In this version the elements of the array AP */
/*     are accessed sequentially with one pass through AP. */

/*     First form  y := beta*y. */

#line 229 "cspmv.f"
    if (beta->r != 1. || beta->i != 0.) {
#line 230 "cspmv.f"
	if (*incy == 1) {
#line 231 "cspmv.f"
	    if (beta->r == 0. && beta->i == 0.) {
#line 232 "cspmv.f"
		i__1 = *n;
#line 232 "cspmv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 233 "cspmv.f"
		    i__2 = i__;
#line 233 "cspmv.f"
		    y[i__2].r = 0., y[i__2].i = 0.;
#line 234 "cspmv.f"
/* L10: */
#line 234 "cspmv.f"
		}
#line 235 "cspmv.f"
	    } else {
#line 236 "cspmv.f"
		i__1 = *n;
#line 236 "cspmv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 237 "cspmv.f"
		    i__2 = i__;
#line 237 "cspmv.f"
		    i__3 = i__;
#line 237 "cspmv.f"
		    z__1.r = beta->r * y[i__3].r - beta->i * y[i__3].i, 
			    z__1.i = beta->r * y[i__3].i + beta->i * y[i__3]
			    .r;
#line 237 "cspmv.f"
		    y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 238 "cspmv.f"
/* L20: */
#line 238 "cspmv.f"
		}
#line 239 "cspmv.f"
	    }
#line 240 "cspmv.f"
	} else {
#line 241 "cspmv.f"
	    iy = ky;
#line 242 "cspmv.f"
	    if (beta->r == 0. && beta->i == 0.) {
#line 243 "cspmv.f"
		i__1 = *n;
#line 243 "cspmv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 244 "cspmv.f"
		    i__2 = iy;
#line 244 "cspmv.f"
		    y[i__2].r = 0., y[i__2].i = 0.;
#line 245 "cspmv.f"
		    iy += *incy;
#line 246 "cspmv.f"
/* L30: */
#line 246 "cspmv.f"
		}
#line 247 "cspmv.f"
	    } else {
#line 248 "cspmv.f"
		i__1 = *n;
#line 248 "cspmv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 249 "cspmv.f"
		    i__2 = iy;
#line 249 "cspmv.f"
		    i__3 = iy;
#line 249 "cspmv.f"
		    z__1.r = beta->r * y[i__3].r - beta->i * y[i__3].i, 
			    z__1.i = beta->r * y[i__3].i + beta->i * y[i__3]
			    .r;
#line 249 "cspmv.f"
		    y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 250 "cspmv.f"
		    iy += *incy;
#line 251 "cspmv.f"
/* L40: */
#line 251 "cspmv.f"
		}
#line 252 "cspmv.f"
	    }
#line 253 "cspmv.f"
	}
#line 254 "cspmv.f"
    }
#line 255 "cspmv.f"
    if (alpha->r == 0. && alpha->i == 0.) {
#line 255 "cspmv.f"
	return 0;
#line 255 "cspmv.f"
    }
#line 257 "cspmv.f"
    kk = 1;
#line 258 "cspmv.f"
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {

/*        Form  y  when AP contains the upper triangle. */

#line 262 "cspmv.f"
	if (*incx == 1 && *incy == 1) {
#line 263 "cspmv.f"
	    i__1 = *n;
#line 263 "cspmv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 264 "cspmv.f"
		i__2 = j;
#line 264 "cspmv.f"
		z__1.r = alpha->r * x[i__2].r - alpha->i * x[i__2].i, z__1.i =
			 alpha->r * x[i__2].i + alpha->i * x[i__2].r;
#line 264 "cspmv.f"
		temp1.r = z__1.r, temp1.i = z__1.i;
#line 265 "cspmv.f"
		temp2.r = 0., temp2.i = 0.;
#line 266 "cspmv.f"
		k = kk;
#line 267 "cspmv.f"
		i__2 = j - 1;
#line 267 "cspmv.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 268 "cspmv.f"
		    i__3 = i__;
#line 268 "cspmv.f"
		    i__4 = i__;
#line 268 "cspmv.f"
		    i__5 = k;
#line 268 "cspmv.f"
		    z__2.r = temp1.r * ap[i__5].r - temp1.i * ap[i__5].i, 
			    z__2.i = temp1.r * ap[i__5].i + temp1.i * ap[i__5]
			    .r;
#line 268 "cspmv.f"
		    z__1.r = y[i__4].r + z__2.r, z__1.i = y[i__4].i + z__2.i;
#line 268 "cspmv.f"
		    y[i__3].r = z__1.r, y[i__3].i = z__1.i;
#line 269 "cspmv.f"
		    i__3 = k;
#line 269 "cspmv.f"
		    i__4 = i__;
#line 269 "cspmv.f"
		    z__2.r = ap[i__3].r * x[i__4].r - ap[i__3].i * x[i__4].i, 
			    z__2.i = ap[i__3].r * x[i__4].i + ap[i__3].i * x[
			    i__4].r;
#line 269 "cspmv.f"
		    z__1.r = temp2.r + z__2.r, z__1.i = temp2.i + z__2.i;
#line 269 "cspmv.f"
		    temp2.r = z__1.r, temp2.i = z__1.i;
#line 270 "cspmv.f"
		    ++k;
#line 271 "cspmv.f"
/* L50: */
#line 271 "cspmv.f"
		}
#line 272 "cspmv.f"
		i__2 = j;
#line 272 "cspmv.f"
		i__3 = j;
#line 272 "cspmv.f"
		i__4 = kk + j - 1;
#line 272 "cspmv.f"
		z__3.r = temp1.r * ap[i__4].r - temp1.i * ap[i__4].i, z__3.i =
			 temp1.r * ap[i__4].i + temp1.i * ap[i__4].r;
#line 272 "cspmv.f"
		z__2.r = y[i__3].r + z__3.r, z__2.i = y[i__3].i + z__3.i;
#line 272 "cspmv.f"
		z__4.r = alpha->r * temp2.r - alpha->i * temp2.i, z__4.i = 
			alpha->r * temp2.i + alpha->i * temp2.r;
#line 272 "cspmv.f"
		z__1.r = z__2.r + z__4.r, z__1.i = z__2.i + z__4.i;
#line 272 "cspmv.f"
		y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 273 "cspmv.f"
		kk += j;
#line 274 "cspmv.f"
/* L60: */
#line 274 "cspmv.f"
	    }
#line 275 "cspmv.f"
	} else {
#line 276 "cspmv.f"
	    jx = kx;
#line 277 "cspmv.f"
	    jy = ky;
#line 278 "cspmv.f"
	    i__1 = *n;
#line 278 "cspmv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 279 "cspmv.f"
		i__2 = jx;
#line 279 "cspmv.f"
		z__1.r = alpha->r * x[i__2].r - alpha->i * x[i__2].i, z__1.i =
			 alpha->r * x[i__2].i + alpha->i * x[i__2].r;
#line 279 "cspmv.f"
		temp1.r = z__1.r, temp1.i = z__1.i;
#line 280 "cspmv.f"
		temp2.r = 0., temp2.i = 0.;
#line 281 "cspmv.f"
		ix = kx;
#line 282 "cspmv.f"
		iy = ky;
#line 283 "cspmv.f"
		i__2 = kk + j - 2;
#line 283 "cspmv.f"
		for (k = kk; k <= i__2; ++k) {
#line 284 "cspmv.f"
		    i__3 = iy;
#line 284 "cspmv.f"
		    i__4 = iy;
#line 284 "cspmv.f"
		    i__5 = k;
#line 284 "cspmv.f"
		    z__2.r = temp1.r * ap[i__5].r - temp1.i * ap[i__5].i, 
			    z__2.i = temp1.r * ap[i__5].i + temp1.i * ap[i__5]
			    .r;
#line 284 "cspmv.f"
		    z__1.r = y[i__4].r + z__2.r, z__1.i = y[i__4].i + z__2.i;
#line 284 "cspmv.f"
		    y[i__3].r = z__1.r, y[i__3].i = z__1.i;
#line 285 "cspmv.f"
		    i__3 = k;
#line 285 "cspmv.f"
		    i__4 = ix;
#line 285 "cspmv.f"
		    z__2.r = ap[i__3].r * x[i__4].r - ap[i__3].i * x[i__4].i, 
			    z__2.i = ap[i__3].r * x[i__4].i + ap[i__3].i * x[
			    i__4].r;
#line 285 "cspmv.f"
		    z__1.r = temp2.r + z__2.r, z__1.i = temp2.i + z__2.i;
#line 285 "cspmv.f"
		    temp2.r = z__1.r, temp2.i = z__1.i;
#line 286 "cspmv.f"
		    ix += *incx;
#line 287 "cspmv.f"
		    iy += *incy;
#line 288 "cspmv.f"
/* L70: */
#line 288 "cspmv.f"
		}
#line 289 "cspmv.f"
		i__2 = jy;
#line 289 "cspmv.f"
		i__3 = jy;
#line 289 "cspmv.f"
		i__4 = kk + j - 1;
#line 289 "cspmv.f"
		z__3.r = temp1.r * ap[i__4].r - temp1.i * ap[i__4].i, z__3.i =
			 temp1.r * ap[i__4].i + temp1.i * ap[i__4].r;
#line 289 "cspmv.f"
		z__2.r = y[i__3].r + z__3.r, z__2.i = y[i__3].i + z__3.i;
#line 289 "cspmv.f"
		z__4.r = alpha->r * temp2.r - alpha->i * temp2.i, z__4.i = 
			alpha->r * temp2.i + alpha->i * temp2.r;
#line 289 "cspmv.f"
		z__1.r = z__2.r + z__4.r, z__1.i = z__2.i + z__4.i;
#line 289 "cspmv.f"
		y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 290 "cspmv.f"
		jx += *incx;
#line 291 "cspmv.f"
		jy += *incy;
#line 292 "cspmv.f"
		kk += j;
#line 293 "cspmv.f"
/* L80: */
#line 293 "cspmv.f"
	    }
#line 294 "cspmv.f"
	}
#line 295 "cspmv.f"
    } else {

/*        Form  y  when AP contains the lower triangle. */

#line 299 "cspmv.f"
	if (*incx == 1 && *incy == 1) {
#line 300 "cspmv.f"
	    i__1 = *n;
#line 300 "cspmv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 301 "cspmv.f"
		i__2 = j;
#line 301 "cspmv.f"
		z__1.r = alpha->r * x[i__2].r - alpha->i * x[i__2].i, z__1.i =
			 alpha->r * x[i__2].i + alpha->i * x[i__2].r;
#line 301 "cspmv.f"
		temp1.r = z__1.r, temp1.i = z__1.i;
#line 302 "cspmv.f"
		temp2.r = 0., temp2.i = 0.;
#line 303 "cspmv.f"
		i__2 = j;
#line 303 "cspmv.f"
		i__3 = j;
#line 303 "cspmv.f"
		i__4 = kk;
#line 303 "cspmv.f"
		z__2.r = temp1.r * ap[i__4].r - temp1.i * ap[i__4].i, z__2.i =
			 temp1.r * ap[i__4].i + temp1.i * ap[i__4].r;
#line 303 "cspmv.f"
		z__1.r = y[i__3].r + z__2.r, z__1.i = y[i__3].i + z__2.i;
#line 303 "cspmv.f"
		y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 304 "cspmv.f"
		k = kk + 1;
#line 305 "cspmv.f"
		i__2 = *n;
#line 305 "cspmv.f"
		for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 306 "cspmv.f"
		    i__3 = i__;
#line 306 "cspmv.f"
		    i__4 = i__;
#line 306 "cspmv.f"
		    i__5 = k;
#line 306 "cspmv.f"
		    z__2.r = temp1.r * ap[i__5].r - temp1.i * ap[i__5].i, 
			    z__2.i = temp1.r * ap[i__5].i + temp1.i * ap[i__5]
			    .r;
#line 306 "cspmv.f"
		    z__1.r = y[i__4].r + z__2.r, z__1.i = y[i__4].i + z__2.i;
#line 306 "cspmv.f"
		    y[i__3].r = z__1.r, y[i__3].i = z__1.i;
#line 307 "cspmv.f"
		    i__3 = k;
#line 307 "cspmv.f"
		    i__4 = i__;
#line 307 "cspmv.f"
		    z__2.r = ap[i__3].r * x[i__4].r - ap[i__3].i * x[i__4].i, 
			    z__2.i = ap[i__3].r * x[i__4].i + ap[i__3].i * x[
			    i__4].r;
#line 307 "cspmv.f"
		    z__1.r = temp2.r + z__2.r, z__1.i = temp2.i + z__2.i;
#line 307 "cspmv.f"
		    temp2.r = z__1.r, temp2.i = z__1.i;
#line 308 "cspmv.f"
		    ++k;
#line 309 "cspmv.f"
/* L90: */
#line 309 "cspmv.f"
		}
#line 310 "cspmv.f"
		i__2 = j;
#line 310 "cspmv.f"
		i__3 = j;
#line 310 "cspmv.f"
		z__2.r = alpha->r * temp2.r - alpha->i * temp2.i, z__2.i = 
			alpha->r * temp2.i + alpha->i * temp2.r;
#line 310 "cspmv.f"
		z__1.r = y[i__3].r + z__2.r, z__1.i = y[i__3].i + z__2.i;
#line 310 "cspmv.f"
		y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 311 "cspmv.f"
		kk += *n - j + 1;
#line 312 "cspmv.f"
/* L100: */
#line 312 "cspmv.f"
	    }
#line 313 "cspmv.f"
	} else {
#line 314 "cspmv.f"
	    jx = kx;
#line 315 "cspmv.f"
	    jy = ky;
#line 316 "cspmv.f"
	    i__1 = *n;
#line 316 "cspmv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 317 "cspmv.f"
		i__2 = jx;
#line 317 "cspmv.f"
		z__1.r = alpha->r * x[i__2].r - alpha->i * x[i__2].i, z__1.i =
			 alpha->r * x[i__2].i + alpha->i * x[i__2].r;
#line 317 "cspmv.f"
		temp1.r = z__1.r, temp1.i = z__1.i;
#line 318 "cspmv.f"
		temp2.r = 0., temp2.i = 0.;
#line 319 "cspmv.f"
		i__2 = jy;
#line 319 "cspmv.f"
		i__3 = jy;
#line 319 "cspmv.f"
		i__4 = kk;
#line 319 "cspmv.f"
		z__2.r = temp1.r * ap[i__4].r - temp1.i * ap[i__4].i, z__2.i =
			 temp1.r * ap[i__4].i + temp1.i * ap[i__4].r;
#line 319 "cspmv.f"
		z__1.r = y[i__3].r + z__2.r, z__1.i = y[i__3].i + z__2.i;
#line 319 "cspmv.f"
		y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 320 "cspmv.f"
		ix = jx;
#line 321 "cspmv.f"
		iy = jy;
#line 322 "cspmv.f"
		i__2 = kk + *n - j;
#line 322 "cspmv.f"
		for (k = kk + 1; k <= i__2; ++k) {
#line 323 "cspmv.f"
		    ix += *incx;
#line 324 "cspmv.f"
		    iy += *incy;
#line 325 "cspmv.f"
		    i__3 = iy;
#line 325 "cspmv.f"
		    i__4 = iy;
#line 325 "cspmv.f"
		    i__5 = k;
#line 325 "cspmv.f"
		    z__2.r = temp1.r * ap[i__5].r - temp1.i * ap[i__5].i, 
			    z__2.i = temp1.r * ap[i__5].i + temp1.i * ap[i__5]
			    .r;
#line 325 "cspmv.f"
		    z__1.r = y[i__4].r + z__2.r, z__1.i = y[i__4].i + z__2.i;
#line 325 "cspmv.f"
		    y[i__3].r = z__1.r, y[i__3].i = z__1.i;
#line 326 "cspmv.f"
		    i__3 = k;
#line 326 "cspmv.f"
		    i__4 = ix;
#line 326 "cspmv.f"
		    z__2.r = ap[i__3].r * x[i__4].r - ap[i__3].i * x[i__4].i, 
			    z__2.i = ap[i__3].r * x[i__4].i + ap[i__3].i * x[
			    i__4].r;
#line 326 "cspmv.f"
		    z__1.r = temp2.r + z__2.r, z__1.i = temp2.i + z__2.i;
#line 326 "cspmv.f"
		    temp2.r = z__1.r, temp2.i = z__1.i;
#line 327 "cspmv.f"
/* L110: */
#line 327 "cspmv.f"
		}
#line 328 "cspmv.f"
		i__2 = jy;
#line 328 "cspmv.f"
		i__3 = jy;
#line 328 "cspmv.f"
		z__2.r = alpha->r * temp2.r - alpha->i * temp2.i, z__2.i = 
			alpha->r * temp2.i + alpha->i * temp2.r;
#line 328 "cspmv.f"
		z__1.r = y[i__3].r + z__2.r, z__1.i = y[i__3].i + z__2.i;
#line 328 "cspmv.f"
		y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 329 "cspmv.f"
		jx += *incx;
#line 330 "cspmv.f"
		jy += *incy;
#line 331 "cspmv.f"
		kk += *n - j + 1;
#line 332 "cspmv.f"
/* L120: */
#line 332 "cspmv.f"
	    }
#line 333 "cspmv.f"
	}
#line 334 "cspmv.f"
    }

#line 336 "cspmv.f"
    return 0;

/*     End of CSPMV */

} /* cspmv_ */

