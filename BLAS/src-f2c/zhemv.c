#line 1 "zhemv.f"
/* zhemv.f -- translated by f2c (version 20100827).
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

#line 1 "zhemv.f"
/* > \brief \b ZHEMV */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZHEMV(UPLO,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY) */

/*       .. Scalar Arguments .. */
/*       COMPLEX*16 ALPHA,BETA */
/*       INTEGER INCX,INCY,LDA,N */
/*       CHARACTER UPLO */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16 A(LDA,*),X(*),Y(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZHEMV  performs the matrix-vector  operation */
/* > */
/* >    y := alpha*A*x + beta*y, */
/* > */
/* > where alpha and beta are scalars, x and y are n element vectors and */
/* > A is an n by n hermitian matrix. */
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
/* > \param[in] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array of DIMENSION ( LDA, n ). */
/* >           Before entry with  UPLO = 'U' or 'u', the leading n by n */
/* >           upper triangular part of the array A must contain the upper */
/* >           triangular part of the hermitian matrix and the strictly */
/* >           lower triangular part of A is not referenced. */
/* >           Before entry with UPLO = 'L' or 'l', the leading n by n */
/* >           lower triangular part of the array A must contain the lower */
/* >           triangular part of the hermitian matrix and the strictly */
/* >           upper triangular part of A is not referenced. */
/* >           Note that the imaginary parts of the diagonal elements need */
/* >           not be set and are assumed to be zero. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >           On entry, LDA specifies the first dimension of A as declared */
/* >           in the calling (sub) program. LDA must be at least */
/* >           max( 1, n ). */
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

/* > \date November 2011 */

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
/* Subroutine */ int zhemv_(char *uplo, integer *n, doublecomplex *alpha, 
	doublecomplex *a, integer *lda, doublecomplex *x, integer *incx, 
	doublecomplex *beta, doublecomplex *y, integer *incy, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1;
    doublecomplex z__1, z__2, z__3, z__4;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, j, ix, iy, jx, jy, kx, ky, info;
    static doublecomplex temp1, temp2;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);


/*  -- Reference BLAS level2 routine (version 3.4.0) -- */
/*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */

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

#line 196 "zhemv.f"
    /* Parameter adjustments */
#line 196 "zhemv.f"
    a_dim1 = *lda;
#line 196 "zhemv.f"
    a_offset = 1 + a_dim1;
#line 196 "zhemv.f"
    a -= a_offset;
#line 196 "zhemv.f"
    --x;
#line 196 "zhemv.f"
    --y;
#line 196 "zhemv.f"

#line 196 "zhemv.f"
    /* Function Body */
#line 196 "zhemv.f"
    info = 0;
#line 197 "zhemv.f"
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
#line 198 "zhemv.f"
	info = 1;
#line 199 "zhemv.f"
    } else if (*n < 0) {
#line 200 "zhemv.f"
	info = 2;
#line 201 "zhemv.f"
    } else if (*lda < max(1,*n)) {
#line 202 "zhemv.f"
	info = 5;
#line 203 "zhemv.f"
    } else if (*incx == 0) {
#line 204 "zhemv.f"
	info = 7;
#line 205 "zhemv.f"
    } else if (*incy == 0) {
#line 206 "zhemv.f"
	info = 10;
#line 207 "zhemv.f"
    }
#line 208 "zhemv.f"
    if (info != 0) {
#line 209 "zhemv.f"
	xerbla_("ZHEMV ", &info, (ftnlen)6);
#line 210 "zhemv.f"
	return 0;
#line 211 "zhemv.f"
    }

/*     Quick return if possible. */

#line 215 "zhemv.f"
    if (*n == 0 || alpha->r == 0. && alpha->i == 0. && (beta->r == 1. && 
	    beta->i == 0.)) {
#line 215 "zhemv.f"
	return 0;
#line 215 "zhemv.f"
    }

/*     Set up the start points in  X  and  Y. */

#line 219 "zhemv.f"
    if (*incx > 0) {
#line 220 "zhemv.f"
	kx = 1;
#line 221 "zhemv.f"
    } else {
#line 222 "zhemv.f"
	kx = 1 - (*n - 1) * *incx;
#line 223 "zhemv.f"
    }
#line 224 "zhemv.f"
    if (*incy > 0) {
#line 225 "zhemv.f"
	ky = 1;
#line 226 "zhemv.f"
    } else {
#line 227 "zhemv.f"
	ky = 1 - (*n - 1) * *incy;
#line 228 "zhemv.f"
    }

/*     Start the operations. In this version the elements of A are */
/*     accessed sequentially with one pass through the triangular part */
/*     of A. */

/*     First form  y := beta*y. */

#line 236 "zhemv.f"
    if (beta->r != 1. || beta->i != 0.) {
#line 237 "zhemv.f"
	if (*incy == 1) {
#line 238 "zhemv.f"
	    if (beta->r == 0. && beta->i == 0.) {
#line 239 "zhemv.f"
		i__1 = *n;
#line 239 "zhemv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 240 "zhemv.f"
		    i__2 = i__;
#line 240 "zhemv.f"
		    y[i__2].r = 0., y[i__2].i = 0.;
#line 241 "zhemv.f"
/* L10: */
#line 241 "zhemv.f"
		}
#line 242 "zhemv.f"
	    } else {
#line 243 "zhemv.f"
		i__1 = *n;
#line 243 "zhemv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 244 "zhemv.f"
		    i__2 = i__;
#line 244 "zhemv.f"
		    i__3 = i__;
#line 244 "zhemv.f"
		    z__1.r = beta->r * y[i__3].r - beta->i * y[i__3].i, 
			    z__1.i = beta->r * y[i__3].i + beta->i * y[i__3]
			    .r;
#line 244 "zhemv.f"
		    y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 245 "zhemv.f"
/* L20: */
#line 245 "zhemv.f"
		}
#line 246 "zhemv.f"
	    }
#line 247 "zhemv.f"
	} else {
#line 248 "zhemv.f"
	    iy = ky;
#line 249 "zhemv.f"
	    if (beta->r == 0. && beta->i == 0.) {
#line 250 "zhemv.f"
		i__1 = *n;
#line 250 "zhemv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 251 "zhemv.f"
		    i__2 = iy;
#line 251 "zhemv.f"
		    y[i__2].r = 0., y[i__2].i = 0.;
#line 252 "zhemv.f"
		    iy += *incy;
#line 253 "zhemv.f"
/* L30: */
#line 253 "zhemv.f"
		}
#line 254 "zhemv.f"
	    } else {
#line 255 "zhemv.f"
		i__1 = *n;
#line 255 "zhemv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 256 "zhemv.f"
		    i__2 = iy;
#line 256 "zhemv.f"
		    i__3 = iy;
#line 256 "zhemv.f"
		    z__1.r = beta->r * y[i__3].r - beta->i * y[i__3].i, 
			    z__1.i = beta->r * y[i__3].i + beta->i * y[i__3]
			    .r;
#line 256 "zhemv.f"
		    y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 257 "zhemv.f"
		    iy += *incy;
#line 258 "zhemv.f"
/* L40: */
#line 258 "zhemv.f"
		}
#line 259 "zhemv.f"
	    }
#line 260 "zhemv.f"
	}
#line 261 "zhemv.f"
    }
#line 262 "zhemv.f"
    if (alpha->r == 0. && alpha->i == 0.) {
#line 262 "zhemv.f"
	return 0;
#line 262 "zhemv.f"
    }
#line 263 "zhemv.f"
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {

/*        Form  y  when A is stored in upper triangle. */

#line 267 "zhemv.f"
	if (*incx == 1 && *incy == 1) {
#line 268 "zhemv.f"
	    i__1 = *n;
#line 268 "zhemv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 269 "zhemv.f"
		i__2 = j;
#line 269 "zhemv.f"
		z__1.r = alpha->r * x[i__2].r - alpha->i * x[i__2].i, z__1.i =
			 alpha->r * x[i__2].i + alpha->i * x[i__2].r;
#line 269 "zhemv.f"
		temp1.r = z__1.r, temp1.i = z__1.i;
#line 270 "zhemv.f"
		temp2.r = 0., temp2.i = 0.;
#line 271 "zhemv.f"
		i__2 = j - 1;
#line 271 "zhemv.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 272 "zhemv.f"
		    i__3 = i__;
#line 272 "zhemv.f"
		    i__4 = i__;
#line 272 "zhemv.f"
		    i__5 = i__ + j * a_dim1;
#line 272 "zhemv.f"
		    z__2.r = temp1.r * a[i__5].r - temp1.i * a[i__5].i, 
			    z__2.i = temp1.r * a[i__5].i + temp1.i * a[i__5]
			    .r;
#line 272 "zhemv.f"
		    z__1.r = y[i__4].r + z__2.r, z__1.i = y[i__4].i + z__2.i;
#line 272 "zhemv.f"
		    y[i__3].r = z__1.r, y[i__3].i = z__1.i;
#line 273 "zhemv.f"
		    d_cnjg(&z__3, &a[i__ + j * a_dim1]);
#line 273 "zhemv.f"
		    i__3 = i__;
#line 273 "zhemv.f"
		    z__2.r = z__3.r * x[i__3].r - z__3.i * x[i__3].i, z__2.i =
			     z__3.r * x[i__3].i + z__3.i * x[i__3].r;
#line 273 "zhemv.f"
		    z__1.r = temp2.r + z__2.r, z__1.i = temp2.i + z__2.i;
#line 273 "zhemv.f"
		    temp2.r = z__1.r, temp2.i = z__1.i;
#line 274 "zhemv.f"
/* L50: */
#line 274 "zhemv.f"
		}
#line 275 "zhemv.f"
		i__2 = j;
#line 275 "zhemv.f"
		i__3 = j;
#line 275 "zhemv.f"
		i__4 = j + j * a_dim1;
#line 275 "zhemv.f"
		d__1 = a[i__4].r;
#line 275 "zhemv.f"
		z__3.r = d__1 * temp1.r, z__3.i = d__1 * temp1.i;
#line 275 "zhemv.f"
		z__2.r = y[i__3].r + z__3.r, z__2.i = y[i__3].i + z__3.i;
#line 275 "zhemv.f"
		z__4.r = alpha->r * temp2.r - alpha->i * temp2.i, z__4.i = 
			alpha->r * temp2.i + alpha->i * temp2.r;
#line 275 "zhemv.f"
		z__1.r = z__2.r + z__4.r, z__1.i = z__2.i + z__4.i;
#line 275 "zhemv.f"
		y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 276 "zhemv.f"
/* L60: */
#line 276 "zhemv.f"
	    }
#line 277 "zhemv.f"
	} else {
#line 278 "zhemv.f"
	    jx = kx;
#line 279 "zhemv.f"
	    jy = ky;
#line 280 "zhemv.f"
	    i__1 = *n;
#line 280 "zhemv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 281 "zhemv.f"
		i__2 = jx;
#line 281 "zhemv.f"
		z__1.r = alpha->r * x[i__2].r - alpha->i * x[i__2].i, z__1.i =
			 alpha->r * x[i__2].i + alpha->i * x[i__2].r;
#line 281 "zhemv.f"
		temp1.r = z__1.r, temp1.i = z__1.i;
#line 282 "zhemv.f"
		temp2.r = 0., temp2.i = 0.;
#line 283 "zhemv.f"
		ix = kx;
#line 284 "zhemv.f"
		iy = ky;
#line 285 "zhemv.f"
		i__2 = j - 1;
#line 285 "zhemv.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 286 "zhemv.f"
		    i__3 = iy;
#line 286 "zhemv.f"
		    i__4 = iy;
#line 286 "zhemv.f"
		    i__5 = i__ + j * a_dim1;
#line 286 "zhemv.f"
		    z__2.r = temp1.r * a[i__5].r - temp1.i * a[i__5].i, 
			    z__2.i = temp1.r * a[i__5].i + temp1.i * a[i__5]
			    .r;
#line 286 "zhemv.f"
		    z__1.r = y[i__4].r + z__2.r, z__1.i = y[i__4].i + z__2.i;
#line 286 "zhemv.f"
		    y[i__3].r = z__1.r, y[i__3].i = z__1.i;
#line 287 "zhemv.f"
		    d_cnjg(&z__3, &a[i__ + j * a_dim1]);
#line 287 "zhemv.f"
		    i__3 = ix;
#line 287 "zhemv.f"
		    z__2.r = z__3.r * x[i__3].r - z__3.i * x[i__3].i, z__2.i =
			     z__3.r * x[i__3].i + z__3.i * x[i__3].r;
#line 287 "zhemv.f"
		    z__1.r = temp2.r + z__2.r, z__1.i = temp2.i + z__2.i;
#line 287 "zhemv.f"
		    temp2.r = z__1.r, temp2.i = z__1.i;
#line 288 "zhemv.f"
		    ix += *incx;
#line 289 "zhemv.f"
		    iy += *incy;
#line 290 "zhemv.f"
/* L70: */
#line 290 "zhemv.f"
		}
#line 291 "zhemv.f"
		i__2 = jy;
#line 291 "zhemv.f"
		i__3 = jy;
#line 291 "zhemv.f"
		i__4 = j + j * a_dim1;
#line 291 "zhemv.f"
		d__1 = a[i__4].r;
#line 291 "zhemv.f"
		z__3.r = d__1 * temp1.r, z__3.i = d__1 * temp1.i;
#line 291 "zhemv.f"
		z__2.r = y[i__3].r + z__3.r, z__2.i = y[i__3].i + z__3.i;
#line 291 "zhemv.f"
		z__4.r = alpha->r * temp2.r - alpha->i * temp2.i, z__4.i = 
			alpha->r * temp2.i + alpha->i * temp2.r;
#line 291 "zhemv.f"
		z__1.r = z__2.r + z__4.r, z__1.i = z__2.i + z__4.i;
#line 291 "zhemv.f"
		y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 292 "zhemv.f"
		jx += *incx;
#line 293 "zhemv.f"
		jy += *incy;
#line 294 "zhemv.f"
/* L80: */
#line 294 "zhemv.f"
	    }
#line 295 "zhemv.f"
	}
#line 296 "zhemv.f"
    } else {

/*        Form  y  when A is stored in lower triangle. */

#line 300 "zhemv.f"
	if (*incx == 1 && *incy == 1) {
#line 301 "zhemv.f"
	    i__1 = *n;
#line 301 "zhemv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 302 "zhemv.f"
		i__2 = j;
#line 302 "zhemv.f"
		z__1.r = alpha->r * x[i__2].r - alpha->i * x[i__2].i, z__1.i =
			 alpha->r * x[i__2].i + alpha->i * x[i__2].r;
#line 302 "zhemv.f"
		temp1.r = z__1.r, temp1.i = z__1.i;
#line 303 "zhemv.f"
		temp2.r = 0., temp2.i = 0.;
#line 304 "zhemv.f"
		i__2 = j;
#line 304 "zhemv.f"
		i__3 = j;
#line 304 "zhemv.f"
		i__4 = j + j * a_dim1;
#line 304 "zhemv.f"
		d__1 = a[i__4].r;
#line 304 "zhemv.f"
		z__2.r = d__1 * temp1.r, z__2.i = d__1 * temp1.i;
#line 304 "zhemv.f"
		z__1.r = y[i__3].r + z__2.r, z__1.i = y[i__3].i + z__2.i;
#line 304 "zhemv.f"
		y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 305 "zhemv.f"
		i__2 = *n;
#line 305 "zhemv.f"
		for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 306 "zhemv.f"
		    i__3 = i__;
#line 306 "zhemv.f"
		    i__4 = i__;
#line 306 "zhemv.f"
		    i__5 = i__ + j * a_dim1;
#line 306 "zhemv.f"
		    z__2.r = temp1.r * a[i__5].r - temp1.i * a[i__5].i, 
			    z__2.i = temp1.r * a[i__5].i + temp1.i * a[i__5]
			    .r;
#line 306 "zhemv.f"
		    z__1.r = y[i__4].r + z__2.r, z__1.i = y[i__4].i + z__2.i;
#line 306 "zhemv.f"
		    y[i__3].r = z__1.r, y[i__3].i = z__1.i;
#line 307 "zhemv.f"
		    d_cnjg(&z__3, &a[i__ + j * a_dim1]);
#line 307 "zhemv.f"
		    i__3 = i__;
#line 307 "zhemv.f"
		    z__2.r = z__3.r * x[i__3].r - z__3.i * x[i__3].i, z__2.i =
			     z__3.r * x[i__3].i + z__3.i * x[i__3].r;
#line 307 "zhemv.f"
		    z__1.r = temp2.r + z__2.r, z__1.i = temp2.i + z__2.i;
#line 307 "zhemv.f"
		    temp2.r = z__1.r, temp2.i = z__1.i;
#line 308 "zhemv.f"
/* L90: */
#line 308 "zhemv.f"
		}
#line 309 "zhemv.f"
		i__2 = j;
#line 309 "zhemv.f"
		i__3 = j;
#line 309 "zhemv.f"
		z__2.r = alpha->r * temp2.r - alpha->i * temp2.i, z__2.i = 
			alpha->r * temp2.i + alpha->i * temp2.r;
#line 309 "zhemv.f"
		z__1.r = y[i__3].r + z__2.r, z__1.i = y[i__3].i + z__2.i;
#line 309 "zhemv.f"
		y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 310 "zhemv.f"
/* L100: */
#line 310 "zhemv.f"
	    }
#line 311 "zhemv.f"
	} else {
#line 312 "zhemv.f"
	    jx = kx;
#line 313 "zhemv.f"
	    jy = ky;
#line 314 "zhemv.f"
	    i__1 = *n;
#line 314 "zhemv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 315 "zhemv.f"
		i__2 = jx;
#line 315 "zhemv.f"
		z__1.r = alpha->r * x[i__2].r - alpha->i * x[i__2].i, z__1.i =
			 alpha->r * x[i__2].i + alpha->i * x[i__2].r;
#line 315 "zhemv.f"
		temp1.r = z__1.r, temp1.i = z__1.i;
#line 316 "zhemv.f"
		temp2.r = 0., temp2.i = 0.;
#line 317 "zhemv.f"
		i__2 = jy;
#line 317 "zhemv.f"
		i__3 = jy;
#line 317 "zhemv.f"
		i__4 = j + j * a_dim1;
#line 317 "zhemv.f"
		d__1 = a[i__4].r;
#line 317 "zhemv.f"
		z__2.r = d__1 * temp1.r, z__2.i = d__1 * temp1.i;
#line 317 "zhemv.f"
		z__1.r = y[i__3].r + z__2.r, z__1.i = y[i__3].i + z__2.i;
#line 317 "zhemv.f"
		y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 318 "zhemv.f"
		ix = jx;
#line 319 "zhemv.f"
		iy = jy;
#line 320 "zhemv.f"
		i__2 = *n;
#line 320 "zhemv.f"
		for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 321 "zhemv.f"
		    ix += *incx;
#line 322 "zhemv.f"
		    iy += *incy;
#line 323 "zhemv.f"
		    i__3 = iy;
#line 323 "zhemv.f"
		    i__4 = iy;
#line 323 "zhemv.f"
		    i__5 = i__ + j * a_dim1;
#line 323 "zhemv.f"
		    z__2.r = temp1.r * a[i__5].r - temp1.i * a[i__5].i, 
			    z__2.i = temp1.r * a[i__5].i + temp1.i * a[i__5]
			    .r;
#line 323 "zhemv.f"
		    z__1.r = y[i__4].r + z__2.r, z__1.i = y[i__4].i + z__2.i;
#line 323 "zhemv.f"
		    y[i__3].r = z__1.r, y[i__3].i = z__1.i;
#line 324 "zhemv.f"
		    d_cnjg(&z__3, &a[i__ + j * a_dim1]);
#line 324 "zhemv.f"
		    i__3 = ix;
#line 324 "zhemv.f"
		    z__2.r = z__3.r * x[i__3].r - z__3.i * x[i__3].i, z__2.i =
			     z__3.r * x[i__3].i + z__3.i * x[i__3].r;
#line 324 "zhemv.f"
		    z__1.r = temp2.r + z__2.r, z__1.i = temp2.i + z__2.i;
#line 324 "zhemv.f"
		    temp2.r = z__1.r, temp2.i = z__1.i;
#line 325 "zhemv.f"
/* L110: */
#line 325 "zhemv.f"
		}
#line 326 "zhemv.f"
		i__2 = jy;
#line 326 "zhemv.f"
		i__3 = jy;
#line 326 "zhemv.f"
		z__2.r = alpha->r * temp2.r - alpha->i * temp2.i, z__2.i = 
			alpha->r * temp2.i + alpha->i * temp2.r;
#line 326 "zhemv.f"
		z__1.r = y[i__3].r + z__2.r, z__1.i = y[i__3].i + z__2.i;
#line 326 "zhemv.f"
		y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 327 "zhemv.f"
		jx += *incx;
#line 328 "zhemv.f"
		jy += *incy;
#line 329 "zhemv.f"
/* L120: */
#line 329 "zhemv.f"
	    }
#line 330 "zhemv.f"
	}
#line 331 "zhemv.f"
    }

#line 333 "zhemv.f"
    return 0;

/*     End of ZHEMV . */

} /* zhemv_ */

