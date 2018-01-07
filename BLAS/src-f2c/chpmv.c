#line 1 "chpmv.f"
/* chpmv.f -- translated by f2c (version 20100827).
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

#line 1 "chpmv.f"
/* > \brief \b CHPMV */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CHPMV(UPLO,N,ALPHA,AP,X,INCX,BETA,Y,INCY) */

/*       .. Scalar Arguments .. */
/*       COMPLEX ALPHA,BETA */
/*       INTEGER INCX,INCY,N */
/*       CHARACTER UPLO */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX AP(*),X(*),Y(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CHPMV  performs the matrix-vector operation */
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
/* >          ALPHA is COMPLEX */
/* >           On entry, ALPHA specifies the scalar alpha. */
/* > \endverbatim */
/* > */
/* > \param[in] AP */
/* > \verbatim */
/* >          AP is COMPLEX array of DIMENSION at least */
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
/* >          X is COMPLEX array of dimension at least */
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
/* >          BETA is COMPLEX */
/* >           On entry, BETA specifies the scalar beta. When BETA is */
/* >           supplied as zero then Y need not be set on input. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Y */
/* > \verbatim */
/* >          Y is COMPLEX array of dimension at least */
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

/* > \ingroup complex_blas_level2 */

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
/* Subroutine */ int chpmv_(char *uplo, integer *n, doublecomplex *alpha, 
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

#line 191 "chpmv.f"
    /* Parameter adjustments */
#line 191 "chpmv.f"
    --y;
#line 191 "chpmv.f"
    --x;
#line 191 "chpmv.f"
    --ap;
#line 191 "chpmv.f"

#line 191 "chpmv.f"
    /* Function Body */
#line 191 "chpmv.f"
    info = 0;
#line 192 "chpmv.f"
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
#line 193 "chpmv.f"
	info = 1;
#line 194 "chpmv.f"
    } else if (*n < 0) {
#line 195 "chpmv.f"
	info = 2;
#line 196 "chpmv.f"
    } else if (*incx == 0) {
#line 197 "chpmv.f"
	info = 6;
#line 198 "chpmv.f"
    } else if (*incy == 0) {
#line 199 "chpmv.f"
	info = 9;
#line 200 "chpmv.f"
    }
#line 201 "chpmv.f"
    if (info != 0) {
#line 202 "chpmv.f"
	xerbla_("CHPMV ", &info, (ftnlen)6);
#line 203 "chpmv.f"
	return 0;
#line 204 "chpmv.f"
    }

/*     Quick return if possible. */

#line 208 "chpmv.f"
    if (*n == 0 || alpha->r == 0. && alpha->i == 0. && (beta->r == 1. && 
	    beta->i == 0.)) {
#line 208 "chpmv.f"
	return 0;
#line 208 "chpmv.f"
    }

/*     Set up the start points in  X  and  Y. */

#line 212 "chpmv.f"
    if (*incx > 0) {
#line 213 "chpmv.f"
	kx = 1;
#line 214 "chpmv.f"
    } else {
#line 215 "chpmv.f"
	kx = 1 - (*n - 1) * *incx;
#line 216 "chpmv.f"
    }
#line 217 "chpmv.f"
    if (*incy > 0) {
#line 218 "chpmv.f"
	ky = 1;
#line 219 "chpmv.f"
    } else {
#line 220 "chpmv.f"
	ky = 1 - (*n - 1) * *incy;
#line 221 "chpmv.f"
    }

/*     Start the operations. In this version the elements of the array AP */
/*     are accessed sequentially with one pass through AP. */

/*     First form  y := beta*y. */

#line 228 "chpmv.f"
    if (beta->r != 1. || beta->i != 0.) {
#line 229 "chpmv.f"
	if (*incy == 1) {
#line 230 "chpmv.f"
	    if (beta->r == 0. && beta->i == 0.) {
#line 231 "chpmv.f"
		i__1 = *n;
#line 231 "chpmv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 232 "chpmv.f"
		    i__2 = i__;
#line 232 "chpmv.f"
		    y[i__2].r = 0., y[i__2].i = 0.;
#line 233 "chpmv.f"
/* L10: */
#line 233 "chpmv.f"
		}
#line 234 "chpmv.f"
	    } else {
#line 235 "chpmv.f"
		i__1 = *n;
#line 235 "chpmv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 236 "chpmv.f"
		    i__2 = i__;
#line 236 "chpmv.f"
		    i__3 = i__;
#line 236 "chpmv.f"
		    z__1.r = beta->r * y[i__3].r - beta->i * y[i__3].i, 
			    z__1.i = beta->r * y[i__3].i + beta->i * y[i__3]
			    .r;
#line 236 "chpmv.f"
		    y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 237 "chpmv.f"
/* L20: */
#line 237 "chpmv.f"
		}
#line 238 "chpmv.f"
	    }
#line 239 "chpmv.f"
	} else {
#line 240 "chpmv.f"
	    iy = ky;
#line 241 "chpmv.f"
	    if (beta->r == 0. && beta->i == 0.) {
#line 242 "chpmv.f"
		i__1 = *n;
#line 242 "chpmv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 243 "chpmv.f"
		    i__2 = iy;
#line 243 "chpmv.f"
		    y[i__2].r = 0., y[i__2].i = 0.;
#line 244 "chpmv.f"
		    iy += *incy;
#line 245 "chpmv.f"
/* L30: */
#line 245 "chpmv.f"
		}
#line 246 "chpmv.f"
	    } else {
#line 247 "chpmv.f"
		i__1 = *n;
#line 247 "chpmv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 248 "chpmv.f"
		    i__2 = iy;
#line 248 "chpmv.f"
		    i__3 = iy;
#line 248 "chpmv.f"
		    z__1.r = beta->r * y[i__3].r - beta->i * y[i__3].i, 
			    z__1.i = beta->r * y[i__3].i + beta->i * y[i__3]
			    .r;
#line 248 "chpmv.f"
		    y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 249 "chpmv.f"
		    iy += *incy;
#line 250 "chpmv.f"
/* L40: */
#line 250 "chpmv.f"
		}
#line 251 "chpmv.f"
	    }
#line 252 "chpmv.f"
	}
#line 253 "chpmv.f"
    }
#line 254 "chpmv.f"
    if (alpha->r == 0. && alpha->i == 0.) {
#line 254 "chpmv.f"
	return 0;
#line 254 "chpmv.f"
    }
#line 255 "chpmv.f"
    kk = 1;
#line 256 "chpmv.f"
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {

/*        Form  y  when AP contains the upper triangle. */

#line 260 "chpmv.f"
	if (*incx == 1 && *incy == 1) {
#line 261 "chpmv.f"
	    i__1 = *n;
#line 261 "chpmv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 262 "chpmv.f"
		i__2 = j;
#line 262 "chpmv.f"
		z__1.r = alpha->r * x[i__2].r - alpha->i * x[i__2].i, z__1.i =
			 alpha->r * x[i__2].i + alpha->i * x[i__2].r;
#line 262 "chpmv.f"
		temp1.r = z__1.r, temp1.i = z__1.i;
#line 263 "chpmv.f"
		temp2.r = 0., temp2.i = 0.;
#line 264 "chpmv.f"
		k = kk;
#line 265 "chpmv.f"
		i__2 = j - 1;
#line 265 "chpmv.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 266 "chpmv.f"
		    i__3 = i__;
#line 266 "chpmv.f"
		    i__4 = i__;
#line 266 "chpmv.f"
		    i__5 = k;
#line 266 "chpmv.f"
		    z__2.r = temp1.r * ap[i__5].r - temp1.i * ap[i__5].i, 
			    z__2.i = temp1.r * ap[i__5].i + temp1.i * ap[i__5]
			    .r;
#line 266 "chpmv.f"
		    z__1.r = y[i__4].r + z__2.r, z__1.i = y[i__4].i + z__2.i;
#line 266 "chpmv.f"
		    y[i__3].r = z__1.r, y[i__3].i = z__1.i;
#line 267 "chpmv.f"
		    d_cnjg(&z__3, &ap[k]);
#line 267 "chpmv.f"
		    i__3 = i__;
#line 267 "chpmv.f"
		    z__2.r = z__3.r * x[i__3].r - z__3.i * x[i__3].i, z__2.i =
			     z__3.r * x[i__3].i + z__3.i * x[i__3].r;
#line 267 "chpmv.f"
		    z__1.r = temp2.r + z__2.r, z__1.i = temp2.i + z__2.i;
#line 267 "chpmv.f"
		    temp2.r = z__1.r, temp2.i = z__1.i;
#line 268 "chpmv.f"
		    ++k;
#line 269 "chpmv.f"
/* L50: */
#line 269 "chpmv.f"
		}
#line 270 "chpmv.f"
		i__2 = j;
#line 270 "chpmv.f"
		i__3 = j;
#line 270 "chpmv.f"
		i__4 = kk + j - 1;
#line 270 "chpmv.f"
		d__1 = ap[i__4].r;
#line 270 "chpmv.f"
		z__3.r = d__1 * temp1.r, z__3.i = d__1 * temp1.i;
#line 270 "chpmv.f"
		z__2.r = y[i__3].r + z__3.r, z__2.i = y[i__3].i + z__3.i;
#line 270 "chpmv.f"
		z__4.r = alpha->r * temp2.r - alpha->i * temp2.i, z__4.i = 
			alpha->r * temp2.i + alpha->i * temp2.r;
#line 270 "chpmv.f"
		z__1.r = z__2.r + z__4.r, z__1.i = z__2.i + z__4.i;
#line 270 "chpmv.f"
		y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 271 "chpmv.f"
		kk += j;
#line 272 "chpmv.f"
/* L60: */
#line 272 "chpmv.f"
	    }
#line 273 "chpmv.f"
	} else {
#line 274 "chpmv.f"
	    jx = kx;
#line 275 "chpmv.f"
	    jy = ky;
#line 276 "chpmv.f"
	    i__1 = *n;
#line 276 "chpmv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 277 "chpmv.f"
		i__2 = jx;
#line 277 "chpmv.f"
		z__1.r = alpha->r * x[i__2].r - alpha->i * x[i__2].i, z__1.i =
			 alpha->r * x[i__2].i + alpha->i * x[i__2].r;
#line 277 "chpmv.f"
		temp1.r = z__1.r, temp1.i = z__1.i;
#line 278 "chpmv.f"
		temp2.r = 0., temp2.i = 0.;
#line 279 "chpmv.f"
		ix = kx;
#line 280 "chpmv.f"
		iy = ky;
#line 281 "chpmv.f"
		i__2 = kk + j - 2;
#line 281 "chpmv.f"
		for (k = kk; k <= i__2; ++k) {
#line 282 "chpmv.f"
		    i__3 = iy;
#line 282 "chpmv.f"
		    i__4 = iy;
#line 282 "chpmv.f"
		    i__5 = k;
#line 282 "chpmv.f"
		    z__2.r = temp1.r * ap[i__5].r - temp1.i * ap[i__5].i, 
			    z__2.i = temp1.r * ap[i__5].i + temp1.i * ap[i__5]
			    .r;
#line 282 "chpmv.f"
		    z__1.r = y[i__4].r + z__2.r, z__1.i = y[i__4].i + z__2.i;
#line 282 "chpmv.f"
		    y[i__3].r = z__1.r, y[i__3].i = z__1.i;
#line 283 "chpmv.f"
		    d_cnjg(&z__3, &ap[k]);
#line 283 "chpmv.f"
		    i__3 = ix;
#line 283 "chpmv.f"
		    z__2.r = z__3.r * x[i__3].r - z__3.i * x[i__3].i, z__2.i =
			     z__3.r * x[i__3].i + z__3.i * x[i__3].r;
#line 283 "chpmv.f"
		    z__1.r = temp2.r + z__2.r, z__1.i = temp2.i + z__2.i;
#line 283 "chpmv.f"
		    temp2.r = z__1.r, temp2.i = z__1.i;
#line 284 "chpmv.f"
		    ix += *incx;
#line 285 "chpmv.f"
		    iy += *incy;
#line 286 "chpmv.f"
/* L70: */
#line 286 "chpmv.f"
		}
#line 287 "chpmv.f"
		i__2 = jy;
#line 287 "chpmv.f"
		i__3 = jy;
#line 287 "chpmv.f"
		i__4 = kk + j - 1;
#line 287 "chpmv.f"
		d__1 = ap[i__4].r;
#line 287 "chpmv.f"
		z__3.r = d__1 * temp1.r, z__3.i = d__1 * temp1.i;
#line 287 "chpmv.f"
		z__2.r = y[i__3].r + z__3.r, z__2.i = y[i__3].i + z__3.i;
#line 287 "chpmv.f"
		z__4.r = alpha->r * temp2.r - alpha->i * temp2.i, z__4.i = 
			alpha->r * temp2.i + alpha->i * temp2.r;
#line 287 "chpmv.f"
		z__1.r = z__2.r + z__4.r, z__1.i = z__2.i + z__4.i;
#line 287 "chpmv.f"
		y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 288 "chpmv.f"
		jx += *incx;
#line 289 "chpmv.f"
		jy += *incy;
#line 290 "chpmv.f"
		kk += j;
#line 291 "chpmv.f"
/* L80: */
#line 291 "chpmv.f"
	    }
#line 292 "chpmv.f"
	}
#line 293 "chpmv.f"
    } else {

/*        Form  y  when AP contains the lower triangle. */

#line 297 "chpmv.f"
	if (*incx == 1 && *incy == 1) {
#line 298 "chpmv.f"
	    i__1 = *n;
#line 298 "chpmv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 299 "chpmv.f"
		i__2 = j;
#line 299 "chpmv.f"
		z__1.r = alpha->r * x[i__2].r - alpha->i * x[i__2].i, z__1.i =
			 alpha->r * x[i__2].i + alpha->i * x[i__2].r;
#line 299 "chpmv.f"
		temp1.r = z__1.r, temp1.i = z__1.i;
#line 300 "chpmv.f"
		temp2.r = 0., temp2.i = 0.;
#line 301 "chpmv.f"
		i__2 = j;
#line 301 "chpmv.f"
		i__3 = j;
#line 301 "chpmv.f"
		i__4 = kk;
#line 301 "chpmv.f"
		d__1 = ap[i__4].r;
#line 301 "chpmv.f"
		z__2.r = d__1 * temp1.r, z__2.i = d__1 * temp1.i;
#line 301 "chpmv.f"
		z__1.r = y[i__3].r + z__2.r, z__1.i = y[i__3].i + z__2.i;
#line 301 "chpmv.f"
		y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 302 "chpmv.f"
		k = kk + 1;
#line 303 "chpmv.f"
		i__2 = *n;
#line 303 "chpmv.f"
		for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 304 "chpmv.f"
		    i__3 = i__;
#line 304 "chpmv.f"
		    i__4 = i__;
#line 304 "chpmv.f"
		    i__5 = k;
#line 304 "chpmv.f"
		    z__2.r = temp1.r * ap[i__5].r - temp1.i * ap[i__5].i, 
			    z__2.i = temp1.r * ap[i__5].i + temp1.i * ap[i__5]
			    .r;
#line 304 "chpmv.f"
		    z__1.r = y[i__4].r + z__2.r, z__1.i = y[i__4].i + z__2.i;
#line 304 "chpmv.f"
		    y[i__3].r = z__1.r, y[i__3].i = z__1.i;
#line 305 "chpmv.f"
		    d_cnjg(&z__3, &ap[k]);
#line 305 "chpmv.f"
		    i__3 = i__;
#line 305 "chpmv.f"
		    z__2.r = z__3.r * x[i__3].r - z__3.i * x[i__3].i, z__2.i =
			     z__3.r * x[i__3].i + z__3.i * x[i__3].r;
#line 305 "chpmv.f"
		    z__1.r = temp2.r + z__2.r, z__1.i = temp2.i + z__2.i;
#line 305 "chpmv.f"
		    temp2.r = z__1.r, temp2.i = z__1.i;
#line 306 "chpmv.f"
		    ++k;
#line 307 "chpmv.f"
/* L90: */
#line 307 "chpmv.f"
		}
#line 308 "chpmv.f"
		i__2 = j;
#line 308 "chpmv.f"
		i__3 = j;
#line 308 "chpmv.f"
		z__2.r = alpha->r * temp2.r - alpha->i * temp2.i, z__2.i = 
			alpha->r * temp2.i + alpha->i * temp2.r;
#line 308 "chpmv.f"
		z__1.r = y[i__3].r + z__2.r, z__1.i = y[i__3].i + z__2.i;
#line 308 "chpmv.f"
		y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 309 "chpmv.f"
		kk += *n - j + 1;
#line 310 "chpmv.f"
/* L100: */
#line 310 "chpmv.f"
	    }
#line 311 "chpmv.f"
	} else {
#line 312 "chpmv.f"
	    jx = kx;
#line 313 "chpmv.f"
	    jy = ky;
#line 314 "chpmv.f"
	    i__1 = *n;
#line 314 "chpmv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 315 "chpmv.f"
		i__2 = jx;
#line 315 "chpmv.f"
		z__1.r = alpha->r * x[i__2].r - alpha->i * x[i__2].i, z__1.i =
			 alpha->r * x[i__2].i + alpha->i * x[i__2].r;
#line 315 "chpmv.f"
		temp1.r = z__1.r, temp1.i = z__1.i;
#line 316 "chpmv.f"
		temp2.r = 0., temp2.i = 0.;
#line 317 "chpmv.f"
		i__2 = jy;
#line 317 "chpmv.f"
		i__3 = jy;
#line 317 "chpmv.f"
		i__4 = kk;
#line 317 "chpmv.f"
		d__1 = ap[i__4].r;
#line 317 "chpmv.f"
		z__2.r = d__1 * temp1.r, z__2.i = d__1 * temp1.i;
#line 317 "chpmv.f"
		z__1.r = y[i__3].r + z__2.r, z__1.i = y[i__3].i + z__2.i;
#line 317 "chpmv.f"
		y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 318 "chpmv.f"
		ix = jx;
#line 319 "chpmv.f"
		iy = jy;
#line 320 "chpmv.f"
		i__2 = kk + *n - j;
#line 320 "chpmv.f"
		for (k = kk + 1; k <= i__2; ++k) {
#line 321 "chpmv.f"
		    ix += *incx;
#line 322 "chpmv.f"
		    iy += *incy;
#line 323 "chpmv.f"
		    i__3 = iy;
#line 323 "chpmv.f"
		    i__4 = iy;
#line 323 "chpmv.f"
		    i__5 = k;
#line 323 "chpmv.f"
		    z__2.r = temp1.r * ap[i__5].r - temp1.i * ap[i__5].i, 
			    z__2.i = temp1.r * ap[i__5].i + temp1.i * ap[i__5]
			    .r;
#line 323 "chpmv.f"
		    z__1.r = y[i__4].r + z__2.r, z__1.i = y[i__4].i + z__2.i;
#line 323 "chpmv.f"
		    y[i__3].r = z__1.r, y[i__3].i = z__1.i;
#line 324 "chpmv.f"
		    d_cnjg(&z__3, &ap[k]);
#line 324 "chpmv.f"
		    i__3 = ix;
#line 324 "chpmv.f"
		    z__2.r = z__3.r * x[i__3].r - z__3.i * x[i__3].i, z__2.i =
			     z__3.r * x[i__3].i + z__3.i * x[i__3].r;
#line 324 "chpmv.f"
		    z__1.r = temp2.r + z__2.r, z__1.i = temp2.i + z__2.i;
#line 324 "chpmv.f"
		    temp2.r = z__1.r, temp2.i = z__1.i;
#line 325 "chpmv.f"
/* L110: */
#line 325 "chpmv.f"
		}
#line 326 "chpmv.f"
		i__2 = jy;
#line 326 "chpmv.f"
		i__3 = jy;
#line 326 "chpmv.f"
		z__2.r = alpha->r * temp2.r - alpha->i * temp2.i, z__2.i = 
			alpha->r * temp2.i + alpha->i * temp2.r;
#line 326 "chpmv.f"
		z__1.r = y[i__3].r + z__2.r, z__1.i = y[i__3].i + z__2.i;
#line 326 "chpmv.f"
		y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 327 "chpmv.f"
		jx += *incx;
#line 328 "chpmv.f"
		jy += *incy;
#line 329 "chpmv.f"
		kk += *n - j + 1;
#line 330 "chpmv.f"
/* L120: */
#line 330 "chpmv.f"
	    }
#line 331 "chpmv.f"
	}
#line 332 "chpmv.f"
    }

#line 334 "chpmv.f"
    return 0;

/*     End of CHPMV . */

} /* chpmv_ */

