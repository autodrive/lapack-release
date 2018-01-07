#line 1 "zhpmv.f"
/* zhpmv.f -- translated by f2c (version 20100827).
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

#line 1 "zhpmv.f"
/* > \brief \b ZHPMV */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZHPMV(UPLO,N,ALPHA,AP,X,INCX,BETA,Y,INCY) */

/*       .. Scalar Arguments .. */
/*       COMPLEX*16 ALPHA,BETA */
/*       INTEGER INCX,INCY,N */
/*       CHARACTER UPLO */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16 AP(*),X(*),Y(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZHPMV  performs the matrix-vector operation */
/* > */
/* >    y := alpha*A*x + beta*y, */
/* > */
/* > where alpha and beta are scalars, x and y are n element vectors and */
/* > A is an n by n hermitian matrix, supplied in packed form. */
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
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >           On entry, N specifies the order of the matrix A. */
/* >           N must be at least zero. */
/* > \endverbatim */
/* > */
/* > \param[in] ALPHA */
/* > \verbatim */
/* >          ALPHA is COMPLEX*16 */
/* >           On entry, ALPHA specifies the scalar alpha. */
/* > \endverbatim */
/* > */
/* > \param[in] AP */
/* > \verbatim */
/* >          AP is COMPLEX*16 array of DIMENSION at least */
/* >           ( ( n*( n + 1 ) )/2 ). */
/* >           Before entry with UPLO = 'U' or 'u', the array AP must */
/* >           contain the upper triangular part of the hermitian matrix */
/* >           packed sequentially, column by column, so that AP( 1 ) */
/* >           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 ) */
/* >           and a( 2, 2 ) respectively, and so on. */
/* >           Before entry with UPLO = 'L' or 'l', the array AP must */
/* >           contain the lower triangular part of the hermitian matrix */
/* >           packed sequentially, column by column, so that AP( 1 ) */
/* >           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 ) */
/* >           and a( 3, 1 ) respectively, and so on. */
/* >           Note that the imaginary parts of the diagonal elements need */
/* >           not be set and are assumed to be zero. */
/* > \endverbatim */
/* > */
/* > \param[in] X */
/* > \verbatim */
/* >          X is COMPLEX*16 array of dimension at least */
/* >           ( 1 + ( n - 1 )*abs( INCX ) ). */
/* >           Before entry, the incremented array X must contain the n */
/* >           element vector x. */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* >          INCX is INTEGER */
/* >           On entry, INCX specifies the increment for the elements of */
/* >           X. INCX must not be zero. */
/* > \endverbatim */
/* > */
/* > \param[in] BETA */
/* > \verbatim */
/* >          BETA is COMPLEX*16 */
/* >           On entry, BETA specifies the scalar beta. When BETA is */
/* >           supplied as zero then Y need not be set on input. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Y */
/* > \verbatim */
/* >          Y is COMPLEX*16 array of dimension at least */
/* >           ( 1 + ( n - 1 )*abs( INCY ) ). */
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
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup complex16_blas_level2 */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  Level 2 Blas routine. */
/* >  The vector and matrix arguments are not referenced when N = 0, or M = 0 */
/* > */
/* >  -- Written on 22-October-1986. */
/* >     Jack Dongarra, Argonne National Lab. */
/* >     Jeremy Du Croz, Nag Central Office. */
/* >     Sven Hammarling, Nag Central Office. */
/* >     Richard Hanson, Sandia National Labs. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int zhpmv_(char *uplo, integer *n, doublecomplex *alpha, 
	doublecomplex *ap, doublecomplex *x, integer *incx, doublecomplex *
	beta, doublecomplex *y, integer *incy, ftnlen uplo_len)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5;
    doublereal d__1;
    doublecomplex z__1, z__2, z__3, z__4;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, j, k, kk, ix, iy, jx, jy, kx, ky, info;
    static doublecomplex temp1, temp2;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);


/*  -- Reference BLAS level2 routine (version 3.7.0) -- */
/*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    -- */
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
/*     .. Intrinsic Functions .. */
/*     .. */

/*     Test the input parameters. */

#line 191 "zhpmv.f"
    /* Parameter adjustments */
#line 191 "zhpmv.f"
    --y;
#line 191 "zhpmv.f"
    --x;
#line 191 "zhpmv.f"
    --ap;
#line 191 "zhpmv.f"

#line 191 "zhpmv.f"
    /* Function Body */
#line 191 "zhpmv.f"
    info = 0;
#line 192 "zhpmv.f"
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
#line 193 "zhpmv.f"
	info = 1;
#line 194 "zhpmv.f"
    } else if (*n < 0) {
#line 195 "zhpmv.f"
	info = 2;
#line 196 "zhpmv.f"
    } else if (*incx == 0) {
#line 197 "zhpmv.f"
	info = 6;
#line 198 "zhpmv.f"
    } else if (*incy == 0) {
#line 199 "zhpmv.f"
	info = 9;
#line 200 "zhpmv.f"
    }
#line 201 "zhpmv.f"
    if (info != 0) {
#line 202 "zhpmv.f"
	xerbla_("ZHPMV ", &info, (ftnlen)6);
#line 203 "zhpmv.f"
	return 0;
#line 204 "zhpmv.f"
    }

/*     Quick return if possible. */

#line 208 "zhpmv.f"
    if (*n == 0 || alpha->r == 0. && alpha->i == 0. && (beta->r == 1. && 
	    beta->i == 0.)) {
#line 208 "zhpmv.f"
	return 0;
#line 208 "zhpmv.f"
    }

/*     Set up the start points in  X  and  Y. */

#line 212 "zhpmv.f"
    if (*incx > 0) {
#line 213 "zhpmv.f"
	kx = 1;
#line 214 "zhpmv.f"
    } else {
#line 215 "zhpmv.f"
	kx = 1 - (*n - 1) * *incx;
#line 216 "zhpmv.f"
    }
#line 217 "zhpmv.f"
    if (*incy > 0) {
#line 218 "zhpmv.f"
	ky = 1;
#line 219 "zhpmv.f"
    } else {
#line 220 "zhpmv.f"
	ky = 1 - (*n - 1) * *incy;
#line 221 "zhpmv.f"
    }

/*     Start the operations. In this version the elements of the array AP */
/*     are accessed sequentially with one pass through AP. */

/*     First form  y := beta*y. */

#line 228 "zhpmv.f"
    if (beta->r != 1. || beta->i != 0.) {
#line 229 "zhpmv.f"
	if (*incy == 1) {
#line 230 "zhpmv.f"
	    if (beta->r == 0. && beta->i == 0.) {
#line 231 "zhpmv.f"
		i__1 = *n;
#line 231 "zhpmv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 232 "zhpmv.f"
		    i__2 = i__;
#line 232 "zhpmv.f"
		    y[i__2].r = 0., y[i__2].i = 0.;
#line 233 "zhpmv.f"
/* L10: */
#line 233 "zhpmv.f"
		}
#line 234 "zhpmv.f"
	    } else {
#line 235 "zhpmv.f"
		i__1 = *n;
#line 235 "zhpmv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 236 "zhpmv.f"
		    i__2 = i__;
#line 236 "zhpmv.f"
		    i__3 = i__;
#line 236 "zhpmv.f"
		    z__1.r = beta->r * y[i__3].r - beta->i * y[i__3].i, 
			    z__1.i = beta->r * y[i__3].i + beta->i * y[i__3]
			    .r;
#line 236 "zhpmv.f"
		    y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 237 "zhpmv.f"
/* L20: */
#line 237 "zhpmv.f"
		}
#line 238 "zhpmv.f"
	    }
#line 239 "zhpmv.f"
	} else {
#line 240 "zhpmv.f"
	    iy = ky;
#line 241 "zhpmv.f"
	    if (beta->r == 0. && beta->i == 0.) {
#line 242 "zhpmv.f"
		i__1 = *n;
#line 242 "zhpmv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 243 "zhpmv.f"
		    i__2 = iy;
#line 243 "zhpmv.f"
		    y[i__2].r = 0., y[i__2].i = 0.;
#line 244 "zhpmv.f"
		    iy += *incy;
#line 245 "zhpmv.f"
/* L30: */
#line 245 "zhpmv.f"
		}
#line 246 "zhpmv.f"
	    } else {
#line 247 "zhpmv.f"
		i__1 = *n;
#line 247 "zhpmv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 248 "zhpmv.f"
		    i__2 = iy;
#line 248 "zhpmv.f"
		    i__3 = iy;
#line 248 "zhpmv.f"
		    z__1.r = beta->r * y[i__3].r - beta->i * y[i__3].i, 
			    z__1.i = beta->r * y[i__3].i + beta->i * y[i__3]
			    .r;
#line 248 "zhpmv.f"
		    y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 249 "zhpmv.f"
		    iy += *incy;
#line 250 "zhpmv.f"
/* L40: */
#line 250 "zhpmv.f"
		}
#line 251 "zhpmv.f"
	    }
#line 252 "zhpmv.f"
	}
#line 253 "zhpmv.f"
    }
#line 254 "zhpmv.f"
    if (alpha->r == 0. && alpha->i == 0.) {
#line 254 "zhpmv.f"
	return 0;
#line 254 "zhpmv.f"
    }
#line 255 "zhpmv.f"
    kk = 1;
#line 256 "zhpmv.f"
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {

/*        Form  y  when AP contains the upper triangle. */

#line 260 "zhpmv.f"
	if (*incx == 1 && *incy == 1) {
#line 261 "zhpmv.f"
	    i__1 = *n;
#line 261 "zhpmv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 262 "zhpmv.f"
		i__2 = j;
#line 262 "zhpmv.f"
		z__1.r = alpha->r * x[i__2].r - alpha->i * x[i__2].i, z__1.i =
			 alpha->r * x[i__2].i + alpha->i * x[i__2].r;
#line 262 "zhpmv.f"
		temp1.r = z__1.r, temp1.i = z__1.i;
#line 263 "zhpmv.f"
		temp2.r = 0., temp2.i = 0.;
#line 264 "zhpmv.f"
		k = kk;
#line 265 "zhpmv.f"
		i__2 = j - 1;
#line 265 "zhpmv.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 266 "zhpmv.f"
		    i__3 = i__;
#line 266 "zhpmv.f"
		    i__4 = i__;
#line 266 "zhpmv.f"
		    i__5 = k;
#line 266 "zhpmv.f"
		    z__2.r = temp1.r * ap[i__5].r - temp1.i * ap[i__5].i, 
			    z__2.i = temp1.r * ap[i__5].i + temp1.i * ap[i__5]
			    .r;
#line 266 "zhpmv.f"
		    z__1.r = y[i__4].r + z__2.r, z__1.i = y[i__4].i + z__2.i;
#line 266 "zhpmv.f"
		    y[i__3].r = z__1.r, y[i__3].i = z__1.i;
#line 267 "zhpmv.f"
		    d_cnjg(&z__3, &ap[k]);
#line 267 "zhpmv.f"
		    i__3 = i__;
#line 267 "zhpmv.f"
		    z__2.r = z__3.r * x[i__3].r - z__3.i * x[i__3].i, z__2.i =
			     z__3.r * x[i__3].i + z__3.i * x[i__3].r;
#line 267 "zhpmv.f"
		    z__1.r = temp2.r + z__2.r, z__1.i = temp2.i + z__2.i;
#line 267 "zhpmv.f"
		    temp2.r = z__1.r, temp2.i = z__1.i;
#line 268 "zhpmv.f"
		    ++k;
#line 269 "zhpmv.f"
/* L50: */
#line 269 "zhpmv.f"
		}
#line 270 "zhpmv.f"
		i__2 = j;
#line 270 "zhpmv.f"
		i__3 = j;
#line 270 "zhpmv.f"
		i__4 = kk + j - 1;
#line 270 "zhpmv.f"
		d__1 = ap[i__4].r;
#line 270 "zhpmv.f"
		z__3.r = d__1 * temp1.r, z__3.i = d__1 * temp1.i;
#line 270 "zhpmv.f"
		z__2.r = y[i__3].r + z__3.r, z__2.i = y[i__3].i + z__3.i;
#line 270 "zhpmv.f"
		z__4.r = alpha->r * temp2.r - alpha->i * temp2.i, z__4.i = 
			alpha->r * temp2.i + alpha->i * temp2.r;
#line 270 "zhpmv.f"
		z__1.r = z__2.r + z__4.r, z__1.i = z__2.i + z__4.i;
#line 270 "zhpmv.f"
		y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 271 "zhpmv.f"
		kk += j;
#line 272 "zhpmv.f"
/* L60: */
#line 272 "zhpmv.f"
	    }
#line 273 "zhpmv.f"
	} else {
#line 274 "zhpmv.f"
	    jx = kx;
#line 275 "zhpmv.f"
	    jy = ky;
#line 276 "zhpmv.f"
	    i__1 = *n;
#line 276 "zhpmv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 277 "zhpmv.f"
		i__2 = jx;
#line 277 "zhpmv.f"
		z__1.r = alpha->r * x[i__2].r - alpha->i * x[i__2].i, z__1.i =
			 alpha->r * x[i__2].i + alpha->i * x[i__2].r;
#line 277 "zhpmv.f"
		temp1.r = z__1.r, temp1.i = z__1.i;
#line 278 "zhpmv.f"
		temp2.r = 0., temp2.i = 0.;
#line 279 "zhpmv.f"
		ix = kx;
#line 280 "zhpmv.f"
		iy = ky;
#line 281 "zhpmv.f"
		i__2 = kk + j - 2;
#line 281 "zhpmv.f"
		for (k = kk; k <= i__2; ++k) {
#line 282 "zhpmv.f"
		    i__3 = iy;
#line 282 "zhpmv.f"
		    i__4 = iy;
#line 282 "zhpmv.f"
		    i__5 = k;
#line 282 "zhpmv.f"
		    z__2.r = temp1.r * ap[i__5].r - temp1.i * ap[i__5].i, 
			    z__2.i = temp1.r * ap[i__5].i + temp1.i * ap[i__5]
			    .r;
#line 282 "zhpmv.f"
		    z__1.r = y[i__4].r + z__2.r, z__1.i = y[i__4].i + z__2.i;
#line 282 "zhpmv.f"
		    y[i__3].r = z__1.r, y[i__3].i = z__1.i;
#line 283 "zhpmv.f"
		    d_cnjg(&z__3, &ap[k]);
#line 283 "zhpmv.f"
		    i__3 = ix;
#line 283 "zhpmv.f"
		    z__2.r = z__3.r * x[i__3].r - z__3.i * x[i__3].i, z__2.i =
			     z__3.r * x[i__3].i + z__3.i * x[i__3].r;
#line 283 "zhpmv.f"
		    z__1.r = temp2.r + z__2.r, z__1.i = temp2.i + z__2.i;
#line 283 "zhpmv.f"
		    temp2.r = z__1.r, temp2.i = z__1.i;
#line 284 "zhpmv.f"
		    ix += *incx;
#line 285 "zhpmv.f"
		    iy += *incy;
#line 286 "zhpmv.f"
/* L70: */
#line 286 "zhpmv.f"
		}
#line 287 "zhpmv.f"
		i__2 = jy;
#line 287 "zhpmv.f"
		i__3 = jy;
#line 287 "zhpmv.f"
		i__4 = kk + j - 1;
#line 287 "zhpmv.f"
		d__1 = ap[i__4].r;
#line 287 "zhpmv.f"
		z__3.r = d__1 * temp1.r, z__3.i = d__1 * temp1.i;
#line 287 "zhpmv.f"
		z__2.r = y[i__3].r + z__3.r, z__2.i = y[i__3].i + z__3.i;
#line 287 "zhpmv.f"
		z__4.r = alpha->r * temp2.r - alpha->i * temp2.i, z__4.i = 
			alpha->r * temp2.i + alpha->i * temp2.r;
#line 287 "zhpmv.f"
		z__1.r = z__2.r + z__4.r, z__1.i = z__2.i + z__4.i;
#line 287 "zhpmv.f"
		y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 288 "zhpmv.f"
		jx += *incx;
#line 289 "zhpmv.f"
		jy += *incy;
#line 290 "zhpmv.f"
		kk += j;
#line 291 "zhpmv.f"
/* L80: */
#line 291 "zhpmv.f"
	    }
#line 292 "zhpmv.f"
	}
#line 293 "zhpmv.f"
    } else {

/*        Form  y  when AP contains the lower triangle. */

#line 297 "zhpmv.f"
	if (*incx == 1 && *incy == 1) {
#line 298 "zhpmv.f"
	    i__1 = *n;
#line 298 "zhpmv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 299 "zhpmv.f"
		i__2 = j;
#line 299 "zhpmv.f"
		z__1.r = alpha->r * x[i__2].r - alpha->i * x[i__2].i, z__1.i =
			 alpha->r * x[i__2].i + alpha->i * x[i__2].r;
#line 299 "zhpmv.f"
		temp1.r = z__1.r, temp1.i = z__1.i;
#line 300 "zhpmv.f"
		temp2.r = 0., temp2.i = 0.;
#line 301 "zhpmv.f"
		i__2 = j;
#line 301 "zhpmv.f"
		i__3 = j;
#line 301 "zhpmv.f"
		i__4 = kk;
#line 301 "zhpmv.f"
		d__1 = ap[i__4].r;
#line 301 "zhpmv.f"
		z__2.r = d__1 * temp1.r, z__2.i = d__1 * temp1.i;
#line 301 "zhpmv.f"
		z__1.r = y[i__3].r + z__2.r, z__1.i = y[i__3].i + z__2.i;
#line 301 "zhpmv.f"
		y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 302 "zhpmv.f"
		k = kk + 1;
#line 303 "zhpmv.f"
		i__2 = *n;
#line 303 "zhpmv.f"
		for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 304 "zhpmv.f"
		    i__3 = i__;
#line 304 "zhpmv.f"
		    i__4 = i__;
#line 304 "zhpmv.f"
		    i__5 = k;
#line 304 "zhpmv.f"
		    z__2.r = temp1.r * ap[i__5].r - temp1.i * ap[i__5].i, 
			    z__2.i = temp1.r * ap[i__5].i + temp1.i * ap[i__5]
			    .r;
#line 304 "zhpmv.f"
		    z__1.r = y[i__4].r + z__2.r, z__1.i = y[i__4].i + z__2.i;
#line 304 "zhpmv.f"
		    y[i__3].r = z__1.r, y[i__3].i = z__1.i;
#line 305 "zhpmv.f"
		    d_cnjg(&z__3, &ap[k]);
#line 305 "zhpmv.f"
		    i__3 = i__;
#line 305 "zhpmv.f"
		    z__2.r = z__3.r * x[i__3].r - z__3.i * x[i__3].i, z__2.i =
			     z__3.r * x[i__3].i + z__3.i * x[i__3].r;
#line 305 "zhpmv.f"
		    z__1.r = temp2.r + z__2.r, z__1.i = temp2.i + z__2.i;
#line 305 "zhpmv.f"
		    temp2.r = z__1.r, temp2.i = z__1.i;
#line 306 "zhpmv.f"
		    ++k;
#line 307 "zhpmv.f"
/* L90: */
#line 307 "zhpmv.f"
		}
#line 308 "zhpmv.f"
		i__2 = j;
#line 308 "zhpmv.f"
		i__3 = j;
#line 308 "zhpmv.f"
		z__2.r = alpha->r * temp2.r - alpha->i * temp2.i, z__2.i = 
			alpha->r * temp2.i + alpha->i * temp2.r;
#line 308 "zhpmv.f"
		z__1.r = y[i__3].r + z__2.r, z__1.i = y[i__3].i + z__2.i;
#line 308 "zhpmv.f"
		y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 309 "zhpmv.f"
		kk += *n - j + 1;
#line 310 "zhpmv.f"
/* L100: */
#line 310 "zhpmv.f"
	    }
#line 311 "zhpmv.f"
	} else {
#line 312 "zhpmv.f"
	    jx = kx;
#line 313 "zhpmv.f"
	    jy = ky;
#line 314 "zhpmv.f"
	    i__1 = *n;
#line 314 "zhpmv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 315 "zhpmv.f"
		i__2 = jx;
#line 315 "zhpmv.f"
		z__1.r = alpha->r * x[i__2].r - alpha->i * x[i__2].i, z__1.i =
			 alpha->r * x[i__2].i + alpha->i * x[i__2].r;
#line 315 "zhpmv.f"
		temp1.r = z__1.r, temp1.i = z__1.i;
#line 316 "zhpmv.f"
		temp2.r = 0., temp2.i = 0.;
#line 317 "zhpmv.f"
		i__2 = jy;
#line 317 "zhpmv.f"
		i__3 = jy;
#line 317 "zhpmv.f"
		i__4 = kk;
#line 317 "zhpmv.f"
		d__1 = ap[i__4].r;
#line 317 "zhpmv.f"
		z__2.r = d__1 * temp1.r, z__2.i = d__1 * temp1.i;
#line 317 "zhpmv.f"
		z__1.r = y[i__3].r + z__2.r, z__1.i = y[i__3].i + z__2.i;
#line 317 "zhpmv.f"
		y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 318 "zhpmv.f"
		ix = jx;
#line 319 "zhpmv.f"
		iy = jy;
#line 320 "zhpmv.f"
		i__2 = kk + *n - j;
#line 320 "zhpmv.f"
		for (k = kk + 1; k <= i__2; ++k) {
#line 321 "zhpmv.f"
		    ix += *incx;
#line 322 "zhpmv.f"
		    iy += *incy;
#line 323 "zhpmv.f"
		    i__3 = iy;
#line 323 "zhpmv.f"
		    i__4 = iy;
#line 323 "zhpmv.f"
		    i__5 = k;
#line 323 "zhpmv.f"
		    z__2.r = temp1.r * ap[i__5].r - temp1.i * ap[i__5].i, 
			    z__2.i = temp1.r * ap[i__5].i + temp1.i * ap[i__5]
			    .r;
#line 323 "zhpmv.f"
		    z__1.r = y[i__4].r + z__2.r, z__1.i = y[i__4].i + z__2.i;
#line 323 "zhpmv.f"
		    y[i__3].r = z__1.r, y[i__3].i = z__1.i;
#line 324 "zhpmv.f"
		    d_cnjg(&z__3, &ap[k]);
#line 324 "zhpmv.f"
		    i__3 = ix;
#line 324 "zhpmv.f"
		    z__2.r = z__3.r * x[i__3].r - z__3.i * x[i__3].i, z__2.i =
			     z__3.r * x[i__3].i + z__3.i * x[i__3].r;
#line 324 "zhpmv.f"
		    z__1.r = temp2.r + z__2.r, z__1.i = temp2.i + z__2.i;
#line 324 "zhpmv.f"
		    temp2.r = z__1.r, temp2.i = z__1.i;
#line 325 "zhpmv.f"
/* L110: */
#line 325 "zhpmv.f"
		}
#line 326 "zhpmv.f"
		i__2 = jy;
#line 326 "zhpmv.f"
		i__3 = jy;
#line 326 "zhpmv.f"
		z__2.r = alpha->r * temp2.r - alpha->i * temp2.i, z__2.i = 
			alpha->r * temp2.i + alpha->i * temp2.r;
#line 326 "zhpmv.f"
		z__1.r = y[i__3].r + z__2.r, z__1.i = y[i__3].i + z__2.i;
#line 326 "zhpmv.f"
		y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 327 "zhpmv.f"
		jx += *incx;
#line 328 "zhpmv.f"
		jy += *incy;
#line 329 "zhpmv.f"
		kk += *n - j + 1;
#line 330 "zhpmv.f"
/* L120: */
#line 330 "zhpmv.f"
	    }
#line 331 "zhpmv.f"
	}
#line 332 "zhpmv.f"
    }

#line 334 "zhpmv.f"
    return 0;

/*     End of ZHPMV . */

} /* zhpmv_ */

