#line 1 "zsymv.f"
/* zsymv.f -- translated by f2c (version 20100827).
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

#line 1 "zsymv.f"
/* > \brief \b ZSYMV computes a matrix-vector product for a complex symmetric matrix. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZSYMV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zsymv.f
"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zsymv.f
"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zsymv.f
"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZSYMV( UPLO, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INCX, INCY, LDA, N */
/*       COMPLEX*16         ALPHA, BETA */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         A( LDA, * ), X( * ), Y( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZSYMV  performs the matrix-vector  operation */
/* > */
/* >    y := alpha*A*x + beta*y, */
/* > */
/* > where alpha and beta are scalars, x and y are n element vectors and */
/* > A is an n by n symmetric matrix. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >           On entry, UPLO specifies whether the upper or lower */
/* >           triangular part of the array A is to be referenced as */
/* >           follows: */
/* > */
/* >              UPLO = 'U' or 'u'   Only the upper triangular part of A */
/* >                                  is to be referenced. */
/* > */
/* >              UPLO = 'L' or 'l'   Only the lower triangular part of A */
/* >                                  is to be referenced. */
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
/* > \param[in] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, dimension ( LDA, N ) */
/* >           Before entry, with  UPLO = 'U' or 'u', the leading n by n */
/* >           upper triangular part of the array A must contain the upper */
/* >           triangular part of the symmetric matrix and the strictly */
/* >           lower triangular part of A is not referenced. */
/* >           Before entry, with UPLO = 'L' or 'l', the leading n by n */
/* >           lower triangular part of the array A must contain the lower */
/* >           triangular part of the symmetric matrix and the strictly */
/* >           upper triangular part of A is not referenced. */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >           On entry, LDA specifies the first dimension of A as declared */
/* >           in the calling (sub) program. LDA must be at least */
/* >           max( 1, N ). */
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

/* > \ingroup complex16SYauxiliary */

/*  ===================================================================== */
/* Subroutine */ int zsymv_(char *uplo, integer *n, doublecomplex *alpha, 
	doublecomplex *a, integer *lda, doublecomplex *x, integer *incx, 
	doublecomplex *beta, doublecomplex *y, integer *incy, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;
    doublecomplex z__1, z__2, z__3, z__4;

    /* Local variables */
    static integer i__, j, ix, iy, jx, jy, kx, ky, info;
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
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 200 "zsymv.f"
    /* Parameter adjustments */
#line 200 "zsymv.f"
    a_dim1 = *lda;
#line 200 "zsymv.f"
    a_offset = 1 + a_dim1;
#line 200 "zsymv.f"
    a -= a_offset;
#line 200 "zsymv.f"
    --x;
#line 200 "zsymv.f"
    --y;
#line 200 "zsymv.f"

#line 200 "zsymv.f"
    /* Function Body */
#line 200 "zsymv.f"
    info = 0;
#line 201 "zsymv.f"
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
#line 202 "zsymv.f"
	info = 1;
#line 203 "zsymv.f"
    } else if (*n < 0) {
#line 204 "zsymv.f"
	info = 2;
#line 205 "zsymv.f"
    } else if (*lda < max(1,*n)) {
#line 206 "zsymv.f"
	info = 5;
#line 207 "zsymv.f"
    } else if (*incx == 0) {
#line 208 "zsymv.f"
	info = 7;
#line 209 "zsymv.f"
    } else if (*incy == 0) {
#line 210 "zsymv.f"
	info = 10;
#line 211 "zsymv.f"
    }
#line 212 "zsymv.f"
    if (info != 0) {
#line 213 "zsymv.f"
	xerbla_("ZSYMV ", &info, (ftnlen)6);
#line 214 "zsymv.f"
	return 0;
#line 215 "zsymv.f"
    }

/*     Quick return if possible. */

#line 219 "zsymv.f"
    if (*n == 0 || alpha->r == 0. && alpha->i == 0. && (beta->r == 1. && 
	    beta->i == 0.)) {
#line 219 "zsymv.f"
	return 0;
#line 219 "zsymv.f"
    }

/*     Set up the start points in  X  and  Y. */

#line 224 "zsymv.f"
    if (*incx > 0) {
#line 225 "zsymv.f"
	kx = 1;
#line 226 "zsymv.f"
    } else {
#line 227 "zsymv.f"
	kx = 1 - (*n - 1) * *incx;
#line 228 "zsymv.f"
    }
#line 229 "zsymv.f"
    if (*incy > 0) {
#line 230 "zsymv.f"
	ky = 1;
#line 231 "zsymv.f"
    } else {
#line 232 "zsymv.f"
	ky = 1 - (*n - 1) * *incy;
#line 233 "zsymv.f"
    }

/*     Start the operations. In this version the elements of A are */
/*     accessed sequentially with one pass through the triangular part */
/*     of A. */

/*     First form  y := beta*y. */

#line 241 "zsymv.f"
    if (beta->r != 1. || beta->i != 0.) {
#line 242 "zsymv.f"
	if (*incy == 1) {
#line 243 "zsymv.f"
	    if (beta->r == 0. && beta->i == 0.) {
#line 244 "zsymv.f"
		i__1 = *n;
#line 244 "zsymv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 245 "zsymv.f"
		    i__2 = i__;
#line 245 "zsymv.f"
		    y[i__2].r = 0., y[i__2].i = 0.;
#line 246 "zsymv.f"
/* L10: */
#line 246 "zsymv.f"
		}
#line 247 "zsymv.f"
	    } else {
#line 248 "zsymv.f"
		i__1 = *n;
#line 248 "zsymv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 249 "zsymv.f"
		    i__2 = i__;
#line 249 "zsymv.f"
		    i__3 = i__;
#line 249 "zsymv.f"
		    z__1.r = beta->r * y[i__3].r - beta->i * y[i__3].i, 
			    z__1.i = beta->r * y[i__3].i + beta->i * y[i__3]
			    .r;
#line 249 "zsymv.f"
		    y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 250 "zsymv.f"
/* L20: */
#line 250 "zsymv.f"
		}
#line 251 "zsymv.f"
	    }
#line 252 "zsymv.f"
	} else {
#line 253 "zsymv.f"
	    iy = ky;
#line 254 "zsymv.f"
	    if (beta->r == 0. && beta->i == 0.) {
#line 255 "zsymv.f"
		i__1 = *n;
#line 255 "zsymv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 256 "zsymv.f"
		    i__2 = iy;
#line 256 "zsymv.f"
		    y[i__2].r = 0., y[i__2].i = 0.;
#line 257 "zsymv.f"
		    iy += *incy;
#line 258 "zsymv.f"
/* L30: */
#line 258 "zsymv.f"
		}
#line 259 "zsymv.f"
	    } else {
#line 260 "zsymv.f"
		i__1 = *n;
#line 260 "zsymv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 261 "zsymv.f"
		    i__2 = iy;
#line 261 "zsymv.f"
		    i__3 = iy;
#line 261 "zsymv.f"
		    z__1.r = beta->r * y[i__3].r - beta->i * y[i__3].i, 
			    z__1.i = beta->r * y[i__3].i + beta->i * y[i__3]
			    .r;
#line 261 "zsymv.f"
		    y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 262 "zsymv.f"
		    iy += *incy;
#line 263 "zsymv.f"
/* L40: */
#line 263 "zsymv.f"
		}
#line 264 "zsymv.f"
	    }
#line 265 "zsymv.f"
	}
#line 266 "zsymv.f"
    }
#line 267 "zsymv.f"
    if (alpha->r == 0. && alpha->i == 0.) {
#line 267 "zsymv.f"
	return 0;
#line 267 "zsymv.f"
    }
#line 269 "zsymv.f"
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {

/*        Form  y  when A is stored in upper triangle. */

#line 273 "zsymv.f"
	if (*incx == 1 && *incy == 1) {
#line 274 "zsymv.f"
	    i__1 = *n;
#line 274 "zsymv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 275 "zsymv.f"
		i__2 = j;
#line 275 "zsymv.f"
		z__1.r = alpha->r * x[i__2].r - alpha->i * x[i__2].i, z__1.i =
			 alpha->r * x[i__2].i + alpha->i * x[i__2].r;
#line 275 "zsymv.f"
		temp1.r = z__1.r, temp1.i = z__1.i;
#line 276 "zsymv.f"
		temp2.r = 0., temp2.i = 0.;
#line 277 "zsymv.f"
		i__2 = j - 1;
#line 277 "zsymv.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 278 "zsymv.f"
		    i__3 = i__;
#line 278 "zsymv.f"
		    i__4 = i__;
#line 278 "zsymv.f"
		    i__5 = i__ + j * a_dim1;
#line 278 "zsymv.f"
		    z__2.r = temp1.r * a[i__5].r - temp1.i * a[i__5].i, 
			    z__2.i = temp1.r * a[i__5].i + temp1.i * a[i__5]
			    .r;
#line 278 "zsymv.f"
		    z__1.r = y[i__4].r + z__2.r, z__1.i = y[i__4].i + z__2.i;
#line 278 "zsymv.f"
		    y[i__3].r = z__1.r, y[i__3].i = z__1.i;
#line 279 "zsymv.f"
		    i__3 = i__ + j * a_dim1;
#line 279 "zsymv.f"
		    i__4 = i__;
#line 279 "zsymv.f"
		    z__2.r = a[i__3].r * x[i__4].r - a[i__3].i * x[i__4].i, 
			    z__2.i = a[i__3].r * x[i__4].i + a[i__3].i * x[
			    i__4].r;
#line 279 "zsymv.f"
		    z__1.r = temp2.r + z__2.r, z__1.i = temp2.i + z__2.i;
#line 279 "zsymv.f"
		    temp2.r = z__1.r, temp2.i = z__1.i;
#line 280 "zsymv.f"
/* L50: */
#line 280 "zsymv.f"
		}
#line 281 "zsymv.f"
		i__2 = j;
#line 281 "zsymv.f"
		i__3 = j;
#line 281 "zsymv.f"
		i__4 = j + j * a_dim1;
#line 281 "zsymv.f"
		z__3.r = temp1.r * a[i__4].r - temp1.i * a[i__4].i, z__3.i = 
			temp1.r * a[i__4].i + temp1.i * a[i__4].r;
#line 281 "zsymv.f"
		z__2.r = y[i__3].r + z__3.r, z__2.i = y[i__3].i + z__3.i;
#line 281 "zsymv.f"
		z__4.r = alpha->r * temp2.r - alpha->i * temp2.i, z__4.i = 
			alpha->r * temp2.i + alpha->i * temp2.r;
#line 281 "zsymv.f"
		z__1.r = z__2.r + z__4.r, z__1.i = z__2.i + z__4.i;
#line 281 "zsymv.f"
		y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 282 "zsymv.f"
/* L60: */
#line 282 "zsymv.f"
	    }
#line 283 "zsymv.f"
	} else {
#line 284 "zsymv.f"
	    jx = kx;
#line 285 "zsymv.f"
	    jy = ky;
#line 286 "zsymv.f"
	    i__1 = *n;
#line 286 "zsymv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 287 "zsymv.f"
		i__2 = jx;
#line 287 "zsymv.f"
		z__1.r = alpha->r * x[i__2].r - alpha->i * x[i__2].i, z__1.i =
			 alpha->r * x[i__2].i + alpha->i * x[i__2].r;
#line 287 "zsymv.f"
		temp1.r = z__1.r, temp1.i = z__1.i;
#line 288 "zsymv.f"
		temp2.r = 0., temp2.i = 0.;
#line 289 "zsymv.f"
		ix = kx;
#line 290 "zsymv.f"
		iy = ky;
#line 291 "zsymv.f"
		i__2 = j - 1;
#line 291 "zsymv.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 292 "zsymv.f"
		    i__3 = iy;
#line 292 "zsymv.f"
		    i__4 = iy;
#line 292 "zsymv.f"
		    i__5 = i__ + j * a_dim1;
#line 292 "zsymv.f"
		    z__2.r = temp1.r * a[i__5].r - temp1.i * a[i__5].i, 
			    z__2.i = temp1.r * a[i__5].i + temp1.i * a[i__5]
			    .r;
#line 292 "zsymv.f"
		    z__1.r = y[i__4].r + z__2.r, z__1.i = y[i__4].i + z__2.i;
#line 292 "zsymv.f"
		    y[i__3].r = z__1.r, y[i__3].i = z__1.i;
#line 293 "zsymv.f"
		    i__3 = i__ + j * a_dim1;
#line 293 "zsymv.f"
		    i__4 = ix;
#line 293 "zsymv.f"
		    z__2.r = a[i__3].r * x[i__4].r - a[i__3].i * x[i__4].i, 
			    z__2.i = a[i__3].r * x[i__4].i + a[i__3].i * x[
			    i__4].r;
#line 293 "zsymv.f"
		    z__1.r = temp2.r + z__2.r, z__1.i = temp2.i + z__2.i;
#line 293 "zsymv.f"
		    temp2.r = z__1.r, temp2.i = z__1.i;
#line 294 "zsymv.f"
		    ix += *incx;
#line 295 "zsymv.f"
		    iy += *incy;
#line 296 "zsymv.f"
/* L70: */
#line 296 "zsymv.f"
		}
#line 297 "zsymv.f"
		i__2 = jy;
#line 297 "zsymv.f"
		i__3 = jy;
#line 297 "zsymv.f"
		i__4 = j + j * a_dim1;
#line 297 "zsymv.f"
		z__3.r = temp1.r * a[i__4].r - temp1.i * a[i__4].i, z__3.i = 
			temp1.r * a[i__4].i + temp1.i * a[i__4].r;
#line 297 "zsymv.f"
		z__2.r = y[i__3].r + z__3.r, z__2.i = y[i__3].i + z__3.i;
#line 297 "zsymv.f"
		z__4.r = alpha->r * temp2.r - alpha->i * temp2.i, z__4.i = 
			alpha->r * temp2.i + alpha->i * temp2.r;
#line 297 "zsymv.f"
		z__1.r = z__2.r + z__4.r, z__1.i = z__2.i + z__4.i;
#line 297 "zsymv.f"
		y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 298 "zsymv.f"
		jx += *incx;
#line 299 "zsymv.f"
		jy += *incy;
#line 300 "zsymv.f"
/* L80: */
#line 300 "zsymv.f"
	    }
#line 301 "zsymv.f"
	}
#line 302 "zsymv.f"
    } else {

/*        Form  y  when A is stored in lower triangle. */

#line 306 "zsymv.f"
	if (*incx == 1 && *incy == 1) {
#line 307 "zsymv.f"
	    i__1 = *n;
#line 307 "zsymv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 308 "zsymv.f"
		i__2 = j;
#line 308 "zsymv.f"
		z__1.r = alpha->r * x[i__2].r - alpha->i * x[i__2].i, z__1.i =
			 alpha->r * x[i__2].i + alpha->i * x[i__2].r;
#line 308 "zsymv.f"
		temp1.r = z__1.r, temp1.i = z__1.i;
#line 309 "zsymv.f"
		temp2.r = 0., temp2.i = 0.;
#line 310 "zsymv.f"
		i__2 = j;
#line 310 "zsymv.f"
		i__3 = j;
#line 310 "zsymv.f"
		i__4 = j + j * a_dim1;
#line 310 "zsymv.f"
		z__2.r = temp1.r * a[i__4].r - temp1.i * a[i__4].i, z__2.i = 
			temp1.r * a[i__4].i + temp1.i * a[i__4].r;
#line 310 "zsymv.f"
		z__1.r = y[i__3].r + z__2.r, z__1.i = y[i__3].i + z__2.i;
#line 310 "zsymv.f"
		y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 311 "zsymv.f"
		i__2 = *n;
#line 311 "zsymv.f"
		for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 312 "zsymv.f"
		    i__3 = i__;
#line 312 "zsymv.f"
		    i__4 = i__;
#line 312 "zsymv.f"
		    i__5 = i__ + j * a_dim1;
#line 312 "zsymv.f"
		    z__2.r = temp1.r * a[i__5].r - temp1.i * a[i__5].i, 
			    z__2.i = temp1.r * a[i__5].i + temp1.i * a[i__5]
			    .r;
#line 312 "zsymv.f"
		    z__1.r = y[i__4].r + z__2.r, z__1.i = y[i__4].i + z__2.i;
#line 312 "zsymv.f"
		    y[i__3].r = z__1.r, y[i__3].i = z__1.i;
#line 313 "zsymv.f"
		    i__3 = i__ + j * a_dim1;
#line 313 "zsymv.f"
		    i__4 = i__;
#line 313 "zsymv.f"
		    z__2.r = a[i__3].r * x[i__4].r - a[i__3].i * x[i__4].i, 
			    z__2.i = a[i__3].r * x[i__4].i + a[i__3].i * x[
			    i__4].r;
#line 313 "zsymv.f"
		    z__1.r = temp2.r + z__2.r, z__1.i = temp2.i + z__2.i;
#line 313 "zsymv.f"
		    temp2.r = z__1.r, temp2.i = z__1.i;
#line 314 "zsymv.f"
/* L90: */
#line 314 "zsymv.f"
		}
#line 315 "zsymv.f"
		i__2 = j;
#line 315 "zsymv.f"
		i__3 = j;
#line 315 "zsymv.f"
		z__2.r = alpha->r * temp2.r - alpha->i * temp2.i, z__2.i = 
			alpha->r * temp2.i + alpha->i * temp2.r;
#line 315 "zsymv.f"
		z__1.r = y[i__3].r + z__2.r, z__1.i = y[i__3].i + z__2.i;
#line 315 "zsymv.f"
		y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 316 "zsymv.f"
/* L100: */
#line 316 "zsymv.f"
	    }
#line 317 "zsymv.f"
	} else {
#line 318 "zsymv.f"
	    jx = kx;
#line 319 "zsymv.f"
	    jy = ky;
#line 320 "zsymv.f"
	    i__1 = *n;
#line 320 "zsymv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 321 "zsymv.f"
		i__2 = jx;
#line 321 "zsymv.f"
		z__1.r = alpha->r * x[i__2].r - alpha->i * x[i__2].i, z__1.i =
			 alpha->r * x[i__2].i + alpha->i * x[i__2].r;
#line 321 "zsymv.f"
		temp1.r = z__1.r, temp1.i = z__1.i;
#line 322 "zsymv.f"
		temp2.r = 0., temp2.i = 0.;
#line 323 "zsymv.f"
		i__2 = jy;
#line 323 "zsymv.f"
		i__3 = jy;
#line 323 "zsymv.f"
		i__4 = j + j * a_dim1;
#line 323 "zsymv.f"
		z__2.r = temp1.r * a[i__4].r - temp1.i * a[i__4].i, z__2.i = 
			temp1.r * a[i__4].i + temp1.i * a[i__4].r;
#line 323 "zsymv.f"
		z__1.r = y[i__3].r + z__2.r, z__1.i = y[i__3].i + z__2.i;
#line 323 "zsymv.f"
		y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 324 "zsymv.f"
		ix = jx;
#line 325 "zsymv.f"
		iy = jy;
#line 326 "zsymv.f"
		i__2 = *n;
#line 326 "zsymv.f"
		for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 327 "zsymv.f"
		    ix += *incx;
#line 328 "zsymv.f"
		    iy += *incy;
#line 329 "zsymv.f"
		    i__3 = iy;
#line 329 "zsymv.f"
		    i__4 = iy;
#line 329 "zsymv.f"
		    i__5 = i__ + j * a_dim1;
#line 329 "zsymv.f"
		    z__2.r = temp1.r * a[i__5].r - temp1.i * a[i__5].i, 
			    z__2.i = temp1.r * a[i__5].i + temp1.i * a[i__5]
			    .r;
#line 329 "zsymv.f"
		    z__1.r = y[i__4].r + z__2.r, z__1.i = y[i__4].i + z__2.i;
#line 329 "zsymv.f"
		    y[i__3].r = z__1.r, y[i__3].i = z__1.i;
#line 330 "zsymv.f"
		    i__3 = i__ + j * a_dim1;
#line 330 "zsymv.f"
		    i__4 = ix;
#line 330 "zsymv.f"
		    z__2.r = a[i__3].r * x[i__4].r - a[i__3].i * x[i__4].i, 
			    z__2.i = a[i__3].r * x[i__4].i + a[i__3].i * x[
			    i__4].r;
#line 330 "zsymv.f"
		    z__1.r = temp2.r + z__2.r, z__1.i = temp2.i + z__2.i;
#line 330 "zsymv.f"
		    temp2.r = z__1.r, temp2.i = z__1.i;
#line 331 "zsymv.f"
/* L110: */
#line 331 "zsymv.f"
		}
#line 332 "zsymv.f"
		i__2 = jy;
#line 332 "zsymv.f"
		i__3 = jy;
#line 332 "zsymv.f"
		z__2.r = alpha->r * temp2.r - alpha->i * temp2.i, z__2.i = 
			alpha->r * temp2.i + alpha->i * temp2.r;
#line 332 "zsymv.f"
		z__1.r = y[i__3].r + z__2.r, z__1.i = y[i__3].i + z__2.i;
#line 332 "zsymv.f"
		y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 333 "zsymv.f"
		jx += *incx;
#line 334 "zsymv.f"
		jy += *incy;
#line 335 "zsymv.f"
/* L120: */
#line 335 "zsymv.f"
	    }
#line 336 "zsymv.f"
	}
#line 337 "zsymv.f"
    }

#line 339 "zsymv.f"
    return 0;

/*     End of ZSYMV */

} /* zsymv_ */

