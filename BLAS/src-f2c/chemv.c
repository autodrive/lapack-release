#line 1 "chemv.f"
/* chemv.f -- translated by f2c (version 20100827).
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

#line 1 "chemv.f"
/* > \brief \b CHEMV */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CHEMV(UPLO,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY) */

/*       .. Scalar Arguments .. */
/*       COMPLEX ALPHA,BETA */
/*       INTEGER INCX,INCY,LDA,N */
/*       CHARACTER UPLO */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX A(LDA,*),X(*),Y(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CHEMV  performs the matrix-vector  operation */
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
/* >          ALPHA is COMPLEX */
/* >           On entry, ALPHA specifies the scalar alpha. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is COMPLEX array of DIMENSION ( LDA, n ). */
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
/* Subroutine */ int chemv_(char *uplo, integer *n, doublecomplex *alpha, 
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

#line 196 "chemv.f"
    /* Parameter adjustments */
#line 196 "chemv.f"
    a_dim1 = *lda;
#line 196 "chemv.f"
    a_offset = 1 + a_dim1;
#line 196 "chemv.f"
    a -= a_offset;
#line 196 "chemv.f"
    --x;
#line 196 "chemv.f"
    --y;
#line 196 "chemv.f"

#line 196 "chemv.f"
    /* Function Body */
#line 196 "chemv.f"
    info = 0;
#line 197 "chemv.f"
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
#line 198 "chemv.f"
	info = 1;
#line 199 "chemv.f"
    } else if (*n < 0) {
#line 200 "chemv.f"
	info = 2;
#line 201 "chemv.f"
    } else if (*lda < max(1,*n)) {
#line 202 "chemv.f"
	info = 5;
#line 203 "chemv.f"
    } else if (*incx == 0) {
#line 204 "chemv.f"
	info = 7;
#line 205 "chemv.f"
    } else if (*incy == 0) {
#line 206 "chemv.f"
	info = 10;
#line 207 "chemv.f"
    }
#line 208 "chemv.f"
    if (info != 0) {
#line 209 "chemv.f"
	xerbla_("CHEMV ", &info, (ftnlen)6);
#line 210 "chemv.f"
	return 0;
#line 211 "chemv.f"
    }

/*     Quick return if possible. */

#line 215 "chemv.f"
    if (*n == 0 || alpha->r == 0. && alpha->i == 0. && (beta->r == 1. && 
	    beta->i == 0.)) {
#line 215 "chemv.f"
	return 0;
#line 215 "chemv.f"
    }

/*     Set up the start points in  X  and  Y. */

#line 219 "chemv.f"
    if (*incx > 0) {
#line 220 "chemv.f"
	kx = 1;
#line 221 "chemv.f"
    } else {
#line 222 "chemv.f"
	kx = 1 - (*n - 1) * *incx;
#line 223 "chemv.f"
    }
#line 224 "chemv.f"
    if (*incy > 0) {
#line 225 "chemv.f"
	ky = 1;
#line 226 "chemv.f"
    } else {
#line 227 "chemv.f"
	ky = 1 - (*n - 1) * *incy;
#line 228 "chemv.f"
    }

/*     Start the operations. In this version the elements of A are */
/*     accessed sequentially with one pass through the triangular part */
/*     of A. */

/*     First form  y := beta*y. */

#line 236 "chemv.f"
    if (beta->r != 1. || beta->i != 0.) {
#line 237 "chemv.f"
	if (*incy == 1) {
#line 238 "chemv.f"
	    if (beta->r == 0. && beta->i == 0.) {
#line 239 "chemv.f"
		i__1 = *n;
#line 239 "chemv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 240 "chemv.f"
		    i__2 = i__;
#line 240 "chemv.f"
		    y[i__2].r = 0., y[i__2].i = 0.;
#line 241 "chemv.f"
/* L10: */
#line 241 "chemv.f"
		}
#line 242 "chemv.f"
	    } else {
#line 243 "chemv.f"
		i__1 = *n;
#line 243 "chemv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 244 "chemv.f"
		    i__2 = i__;
#line 244 "chemv.f"
		    i__3 = i__;
#line 244 "chemv.f"
		    z__1.r = beta->r * y[i__3].r - beta->i * y[i__3].i, 
			    z__1.i = beta->r * y[i__3].i + beta->i * y[i__3]
			    .r;
#line 244 "chemv.f"
		    y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 245 "chemv.f"
/* L20: */
#line 245 "chemv.f"
		}
#line 246 "chemv.f"
	    }
#line 247 "chemv.f"
	} else {
#line 248 "chemv.f"
	    iy = ky;
#line 249 "chemv.f"
	    if (beta->r == 0. && beta->i == 0.) {
#line 250 "chemv.f"
		i__1 = *n;
#line 250 "chemv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 251 "chemv.f"
		    i__2 = iy;
#line 251 "chemv.f"
		    y[i__2].r = 0., y[i__2].i = 0.;
#line 252 "chemv.f"
		    iy += *incy;
#line 253 "chemv.f"
/* L30: */
#line 253 "chemv.f"
		}
#line 254 "chemv.f"
	    } else {
#line 255 "chemv.f"
		i__1 = *n;
#line 255 "chemv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 256 "chemv.f"
		    i__2 = iy;
#line 256 "chemv.f"
		    i__3 = iy;
#line 256 "chemv.f"
		    z__1.r = beta->r * y[i__3].r - beta->i * y[i__3].i, 
			    z__1.i = beta->r * y[i__3].i + beta->i * y[i__3]
			    .r;
#line 256 "chemv.f"
		    y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 257 "chemv.f"
		    iy += *incy;
#line 258 "chemv.f"
/* L40: */
#line 258 "chemv.f"
		}
#line 259 "chemv.f"
	    }
#line 260 "chemv.f"
	}
#line 261 "chemv.f"
    }
#line 262 "chemv.f"
    if (alpha->r == 0. && alpha->i == 0.) {
#line 262 "chemv.f"
	return 0;
#line 262 "chemv.f"
    }
#line 263 "chemv.f"
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {

/*        Form  y  when A is stored in upper triangle. */

#line 267 "chemv.f"
	if (*incx == 1 && *incy == 1) {
#line 268 "chemv.f"
	    i__1 = *n;
#line 268 "chemv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 269 "chemv.f"
		i__2 = j;
#line 269 "chemv.f"
		z__1.r = alpha->r * x[i__2].r - alpha->i * x[i__2].i, z__1.i =
			 alpha->r * x[i__2].i + alpha->i * x[i__2].r;
#line 269 "chemv.f"
		temp1.r = z__1.r, temp1.i = z__1.i;
#line 270 "chemv.f"
		temp2.r = 0., temp2.i = 0.;
#line 271 "chemv.f"
		i__2 = j - 1;
#line 271 "chemv.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 272 "chemv.f"
		    i__3 = i__;
#line 272 "chemv.f"
		    i__4 = i__;
#line 272 "chemv.f"
		    i__5 = i__ + j * a_dim1;
#line 272 "chemv.f"
		    z__2.r = temp1.r * a[i__5].r - temp1.i * a[i__5].i, 
			    z__2.i = temp1.r * a[i__5].i + temp1.i * a[i__5]
			    .r;
#line 272 "chemv.f"
		    z__1.r = y[i__4].r + z__2.r, z__1.i = y[i__4].i + z__2.i;
#line 272 "chemv.f"
		    y[i__3].r = z__1.r, y[i__3].i = z__1.i;
#line 273 "chemv.f"
		    d_cnjg(&z__3, &a[i__ + j * a_dim1]);
#line 273 "chemv.f"
		    i__3 = i__;
#line 273 "chemv.f"
		    z__2.r = z__3.r * x[i__3].r - z__3.i * x[i__3].i, z__2.i =
			     z__3.r * x[i__3].i + z__3.i * x[i__3].r;
#line 273 "chemv.f"
		    z__1.r = temp2.r + z__2.r, z__1.i = temp2.i + z__2.i;
#line 273 "chemv.f"
		    temp2.r = z__1.r, temp2.i = z__1.i;
#line 274 "chemv.f"
/* L50: */
#line 274 "chemv.f"
		}
#line 275 "chemv.f"
		i__2 = j;
#line 275 "chemv.f"
		i__3 = j;
#line 275 "chemv.f"
		i__4 = j + j * a_dim1;
#line 275 "chemv.f"
		d__1 = a[i__4].r;
#line 275 "chemv.f"
		z__3.r = d__1 * temp1.r, z__3.i = d__1 * temp1.i;
#line 275 "chemv.f"
		z__2.r = y[i__3].r + z__3.r, z__2.i = y[i__3].i + z__3.i;
#line 275 "chemv.f"
		z__4.r = alpha->r * temp2.r - alpha->i * temp2.i, z__4.i = 
			alpha->r * temp2.i + alpha->i * temp2.r;
#line 275 "chemv.f"
		z__1.r = z__2.r + z__4.r, z__1.i = z__2.i + z__4.i;
#line 275 "chemv.f"
		y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 276 "chemv.f"
/* L60: */
#line 276 "chemv.f"
	    }
#line 277 "chemv.f"
	} else {
#line 278 "chemv.f"
	    jx = kx;
#line 279 "chemv.f"
	    jy = ky;
#line 280 "chemv.f"
	    i__1 = *n;
#line 280 "chemv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 281 "chemv.f"
		i__2 = jx;
#line 281 "chemv.f"
		z__1.r = alpha->r * x[i__2].r - alpha->i * x[i__2].i, z__1.i =
			 alpha->r * x[i__2].i + alpha->i * x[i__2].r;
#line 281 "chemv.f"
		temp1.r = z__1.r, temp1.i = z__1.i;
#line 282 "chemv.f"
		temp2.r = 0., temp2.i = 0.;
#line 283 "chemv.f"
		ix = kx;
#line 284 "chemv.f"
		iy = ky;
#line 285 "chemv.f"
		i__2 = j - 1;
#line 285 "chemv.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 286 "chemv.f"
		    i__3 = iy;
#line 286 "chemv.f"
		    i__4 = iy;
#line 286 "chemv.f"
		    i__5 = i__ + j * a_dim1;
#line 286 "chemv.f"
		    z__2.r = temp1.r * a[i__5].r - temp1.i * a[i__5].i, 
			    z__2.i = temp1.r * a[i__5].i + temp1.i * a[i__5]
			    .r;
#line 286 "chemv.f"
		    z__1.r = y[i__4].r + z__2.r, z__1.i = y[i__4].i + z__2.i;
#line 286 "chemv.f"
		    y[i__3].r = z__1.r, y[i__3].i = z__1.i;
#line 287 "chemv.f"
		    d_cnjg(&z__3, &a[i__ + j * a_dim1]);
#line 287 "chemv.f"
		    i__3 = ix;
#line 287 "chemv.f"
		    z__2.r = z__3.r * x[i__3].r - z__3.i * x[i__3].i, z__2.i =
			     z__3.r * x[i__3].i + z__3.i * x[i__3].r;
#line 287 "chemv.f"
		    z__1.r = temp2.r + z__2.r, z__1.i = temp2.i + z__2.i;
#line 287 "chemv.f"
		    temp2.r = z__1.r, temp2.i = z__1.i;
#line 288 "chemv.f"
		    ix += *incx;
#line 289 "chemv.f"
		    iy += *incy;
#line 290 "chemv.f"
/* L70: */
#line 290 "chemv.f"
		}
#line 291 "chemv.f"
		i__2 = jy;
#line 291 "chemv.f"
		i__3 = jy;
#line 291 "chemv.f"
		i__4 = j + j * a_dim1;
#line 291 "chemv.f"
		d__1 = a[i__4].r;
#line 291 "chemv.f"
		z__3.r = d__1 * temp1.r, z__3.i = d__1 * temp1.i;
#line 291 "chemv.f"
		z__2.r = y[i__3].r + z__3.r, z__2.i = y[i__3].i + z__3.i;
#line 291 "chemv.f"
		z__4.r = alpha->r * temp2.r - alpha->i * temp2.i, z__4.i = 
			alpha->r * temp2.i + alpha->i * temp2.r;
#line 291 "chemv.f"
		z__1.r = z__2.r + z__4.r, z__1.i = z__2.i + z__4.i;
#line 291 "chemv.f"
		y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 292 "chemv.f"
		jx += *incx;
#line 293 "chemv.f"
		jy += *incy;
#line 294 "chemv.f"
/* L80: */
#line 294 "chemv.f"
	    }
#line 295 "chemv.f"
	}
#line 296 "chemv.f"
    } else {

/*        Form  y  when A is stored in lower triangle. */

#line 300 "chemv.f"
	if (*incx == 1 && *incy == 1) {
#line 301 "chemv.f"
	    i__1 = *n;
#line 301 "chemv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 302 "chemv.f"
		i__2 = j;
#line 302 "chemv.f"
		z__1.r = alpha->r * x[i__2].r - alpha->i * x[i__2].i, z__1.i =
			 alpha->r * x[i__2].i + alpha->i * x[i__2].r;
#line 302 "chemv.f"
		temp1.r = z__1.r, temp1.i = z__1.i;
#line 303 "chemv.f"
		temp2.r = 0., temp2.i = 0.;
#line 304 "chemv.f"
		i__2 = j;
#line 304 "chemv.f"
		i__3 = j;
#line 304 "chemv.f"
		i__4 = j + j * a_dim1;
#line 304 "chemv.f"
		d__1 = a[i__4].r;
#line 304 "chemv.f"
		z__2.r = d__1 * temp1.r, z__2.i = d__1 * temp1.i;
#line 304 "chemv.f"
		z__1.r = y[i__3].r + z__2.r, z__1.i = y[i__3].i + z__2.i;
#line 304 "chemv.f"
		y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 305 "chemv.f"
		i__2 = *n;
#line 305 "chemv.f"
		for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 306 "chemv.f"
		    i__3 = i__;
#line 306 "chemv.f"
		    i__4 = i__;
#line 306 "chemv.f"
		    i__5 = i__ + j * a_dim1;
#line 306 "chemv.f"
		    z__2.r = temp1.r * a[i__5].r - temp1.i * a[i__5].i, 
			    z__2.i = temp1.r * a[i__5].i + temp1.i * a[i__5]
			    .r;
#line 306 "chemv.f"
		    z__1.r = y[i__4].r + z__2.r, z__1.i = y[i__4].i + z__2.i;
#line 306 "chemv.f"
		    y[i__3].r = z__1.r, y[i__3].i = z__1.i;
#line 307 "chemv.f"
		    d_cnjg(&z__3, &a[i__ + j * a_dim1]);
#line 307 "chemv.f"
		    i__3 = i__;
#line 307 "chemv.f"
		    z__2.r = z__3.r * x[i__3].r - z__3.i * x[i__3].i, z__2.i =
			     z__3.r * x[i__3].i + z__3.i * x[i__3].r;
#line 307 "chemv.f"
		    z__1.r = temp2.r + z__2.r, z__1.i = temp2.i + z__2.i;
#line 307 "chemv.f"
		    temp2.r = z__1.r, temp2.i = z__1.i;
#line 308 "chemv.f"
/* L90: */
#line 308 "chemv.f"
		}
#line 309 "chemv.f"
		i__2 = j;
#line 309 "chemv.f"
		i__3 = j;
#line 309 "chemv.f"
		z__2.r = alpha->r * temp2.r - alpha->i * temp2.i, z__2.i = 
			alpha->r * temp2.i + alpha->i * temp2.r;
#line 309 "chemv.f"
		z__1.r = y[i__3].r + z__2.r, z__1.i = y[i__3].i + z__2.i;
#line 309 "chemv.f"
		y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 310 "chemv.f"
/* L100: */
#line 310 "chemv.f"
	    }
#line 311 "chemv.f"
	} else {
#line 312 "chemv.f"
	    jx = kx;
#line 313 "chemv.f"
	    jy = ky;
#line 314 "chemv.f"
	    i__1 = *n;
#line 314 "chemv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 315 "chemv.f"
		i__2 = jx;
#line 315 "chemv.f"
		z__1.r = alpha->r * x[i__2].r - alpha->i * x[i__2].i, z__1.i =
			 alpha->r * x[i__2].i + alpha->i * x[i__2].r;
#line 315 "chemv.f"
		temp1.r = z__1.r, temp1.i = z__1.i;
#line 316 "chemv.f"
		temp2.r = 0., temp2.i = 0.;
#line 317 "chemv.f"
		i__2 = jy;
#line 317 "chemv.f"
		i__3 = jy;
#line 317 "chemv.f"
		i__4 = j + j * a_dim1;
#line 317 "chemv.f"
		d__1 = a[i__4].r;
#line 317 "chemv.f"
		z__2.r = d__1 * temp1.r, z__2.i = d__1 * temp1.i;
#line 317 "chemv.f"
		z__1.r = y[i__3].r + z__2.r, z__1.i = y[i__3].i + z__2.i;
#line 317 "chemv.f"
		y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 318 "chemv.f"
		ix = jx;
#line 319 "chemv.f"
		iy = jy;
#line 320 "chemv.f"
		i__2 = *n;
#line 320 "chemv.f"
		for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 321 "chemv.f"
		    ix += *incx;
#line 322 "chemv.f"
		    iy += *incy;
#line 323 "chemv.f"
		    i__3 = iy;
#line 323 "chemv.f"
		    i__4 = iy;
#line 323 "chemv.f"
		    i__5 = i__ + j * a_dim1;
#line 323 "chemv.f"
		    z__2.r = temp1.r * a[i__5].r - temp1.i * a[i__5].i, 
			    z__2.i = temp1.r * a[i__5].i + temp1.i * a[i__5]
			    .r;
#line 323 "chemv.f"
		    z__1.r = y[i__4].r + z__2.r, z__1.i = y[i__4].i + z__2.i;
#line 323 "chemv.f"
		    y[i__3].r = z__1.r, y[i__3].i = z__1.i;
#line 324 "chemv.f"
		    d_cnjg(&z__3, &a[i__ + j * a_dim1]);
#line 324 "chemv.f"
		    i__3 = ix;
#line 324 "chemv.f"
		    z__2.r = z__3.r * x[i__3].r - z__3.i * x[i__3].i, z__2.i =
			     z__3.r * x[i__3].i + z__3.i * x[i__3].r;
#line 324 "chemv.f"
		    z__1.r = temp2.r + z__2.r, z__1.i = temp2.i + z__2.i;
#line 324 "chemv.f"
		    temp2.r = z__1.r, temp2.i = z__1.i;
#line 325 "chemv.f"
/* L110: */
#line 325 "chemv.f"
		}
#line 326 "chemv.f"
		i__2 = jy;
#line 326 "chemv.f"
		i__3 = jy;
#line 326 "chemv.f"
		z__2.r = alpha->r * temp2.r - alpha->i * temp2.i, z__2.i = 
			alpha->r * temp2.i + alpha->i * temp2.r;
#line 326 "chemv.f"
		z__1.r = y[i__3].r + z__2.r, z__1.i = y[i__3].i + z__2.i;
#line 326 "chemv.f"
		y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 327 "chemv.f"
		jx += *incx;
#line 328 "chemv.f"
		jy += *incy;
#line 329 "chemv.f"
/* L120: */
#line 329 "chemv.f"
	    }
#line 330 "chemv.f"
	}
#line 331 "chemv.f"
    }

#line 333 "chemv.f"
    return 0;

/*     End of CHEMV . */

} /* chemv_ */

