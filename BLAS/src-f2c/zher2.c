#line 1 "zher2.f"
/* zher2.f -- translated by f2c (version 20100827).
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

#line 1 "zher2.f"
/* > \brief \b ZHER2 */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZHER2(UPLO,N,ALPHA,X,INCX,Y,INCY,A,LDA) */

/*       .. Scalar Arguments .. */
/*       COMPLEX*16 ALPHA */
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
/* > ZHER2  performs the hermitian rank 2 operation */
/* > */
/* >    A := alpha*x*y**H + conjg( alpha )*y*x**H + A, */
/* > */
/* > where alpha is a scalar, x and y are n element vectors and A is an n */
/* > by n hermitian matrix. */
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
/* > \param[in] Y */
/* > \verbatim */
/* >          Y is COMPLEX*16 array of dimension at least */
/* >           ( 1 + ( n - 1 )*abs( INCY ) ). */
/* >           Before entry, the incremented array Y must contain the n */
/* >           element vector y. */
/* > \endverbatim */
/* > */
/* > \param[in] INCY */
/* > \verbatim */
/* >          INCY is INTEGER */
/* >           On entry, INCY specifies the increment for the elements of */
/* >           Y. INCY must not be zero. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array of DIMENSION ( LDA, n ). */
/* >           Before entry with  UPLO = 'U' or 'u', the leading n by n */
/* >           upper triangular part of the array A must contain the upper */
/* >           triangular part of the hermitian matrix and the strictly */
/* >           lower triangular part of A is not referenced. On exit, the */
/* >           upper triangular part of the array A is overwritten by the */
/* >           upper triangular part of the updated matrix. */
/* >           Before entry with UPLO = 'L' or 'l', the leading n by n */
/* >           lower triangular part of the array A must contain the lower */
/* >           triangular part of the hermitian matrix and the strictly */
/* >           upper triangular part of A is not referenced. On exit, the */
/* >           lower triangular part of the array A is overwritten by the */
/* >           lower triangular part of the updated matrix. */
/* >           Note that the imaginary parts of the diagonal elements need */
/* >           not be set, they are assumed to be zero, and on exit they */
/* >           are set to zero. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >           On entry, LDA specifies the first dimension of A as declared */
/* >           in the calling (sub) program. LDA must be at least */
/* >           max( 1, n ). */
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
/* > */
/* >  -- Written on 22-October-1986. */
/* >     Jack Dongarra, Argonne National Lab. */
/* >     Jeremy Du Croz, Nag Central Office. */
/* >     Sven Hammarling, Nag Central Office. */
/* >     Richard Hanson, Sandia National Labs. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int zher2_(char *uplo, integer *n, doublecomplex *alpha, 
	doublecomplex *x, integer *incx, doublecomplex *y, integer *incy, 
	doublecomplex *a, integer *lda, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5, i__6;
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

#line 190 "zher2.f"
    /* Parameter adjustments */
#line 190 "zher2.f"
    --x;
#line 190 "zher2.f"
    --y;
#line 190 "zher2.f"
    a_dim1 = *lda;
#line 190 "zher2.f"
    a_offset = 1 + a_dim1;
#line 190 "zher2.f"
    a -= a_offset;
#line 190 "zher2.f"

#line 190 "zher2.f"
    /* Function Body */
#line 190 "zher2.f"
    info = 0;
#line 191 "zher2.f"
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
#line 192 "zher2.f"
	info = 1;
#line 193 "zher2.f"
    } else if (*n < 0) {
#line 194 "zher2.f"
	info = 2;
#line 195 "zher2.f"
    } else if (*incx == 0) {
#line 196 "zher2.f"
	info = 5;
#line 197 "zher2.f"
    } else if (*incy == 0) {
#line 198 "zher2.f"
	info = 7;
#line 199 "zher2.f"
    } else if (*lda < max(1,*n)) {
#line 200 "zher2.f"
	info = 9;
#line 201 "zher2.f"
    }
#line 202 "zher2.f"
    if (info != 0) {
#line 203 "zher2.f"
	xerbla_("ZHER2 ", &info, (ftnlen)6);
#line 204 "zher2.f"
	return 0;
#line 205 "zher2.f"
    }

/*     Quick return if possible. */

#line 209 "zher2.f"
    if (*n == 0 || alpha->r == 0. && alpha->i == 0.) {
#line 209 "zher2.f"
	return 0;
#line 209 "zher2.f"
    }

/*     Set up the start points in X and Y if the increments are not both */
/*     unity. */

#line 214 "zher2.f"
    if (*incx != 1 || *incy != 1) {
#line 215 "zher2.f"
	if (*incx > 0) {
#line 216 "zher2.f"
	    kx = 1;
#line 217 "zher2.f"
	} else {
#line 218 "zher2.f"
	    kx = 1 - (*n - 1) * *incx;
#line 219 "zher2.f"
	}
#line 220 "zher2.f"
	if (*incy > 0) {
#line 221 "zher2.f"
	    ky = 1;
#line 222 "zher2.f"
	} else {
#line 223 "zher2.f"
	    ky = 1 - (*n - 1) * *incy;
#line 224 "zher2.f"
	}
#line 225 "zher2.f"
	jx = kx;
#line 226 "zher2.f"
	jy = ky;
#line 227 "zher2.f"
    }

/*     Start the operations. In this version the elements of A are */
/*     accessed sequentially with one pass through the triangular part */
/*     of A. */

#line 233 "zher2.f"
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {

/*        Form  A  when A is stored in the upper triangle. */

#line 237 "zher2.f"
	if (*incx == 1 && *incy == 1) {
#line 238 "zher2.f"
	    i__1 = *n;
#line 238 "zher2.f"
	    for (j = 1; j <= i__1; ++j) {
#line 239 "zher2.f"
		i__2 = j;
#line 239 "zher2.f"
		i__3 = j;
#line 239 "zher2.f"
		if (x[i__2].r != 0. || x[i__2].i != 0. || (y[i__3].r != 0. || 
			y[i__3].i != 0.)) {
#line 240 "zher2.f"
		    d_cnjg(&z__2, &y[j]);
#line 240 "zher2.f"
		    z__1.r = alpha->r * z__2.r - alpha->i * z__2.i, z__1.i = 
			    alpha->r * z__2.i + alpha->i * z__2.r;
#line 240 "zher2.f"
		    temp1.r = z__1.r, temp1.i = z__1.i;
#line 241 "zher2.f"
		    i__2 = j;
#line 241 "zher2.f"
		    z__2.r = alpha->r * x[i__2].r - alpha->i * x[i__2].i, 
			    z__2.i = alpha->r * x[i__2].i + alpha->i * x[i__2]
			    .r;
#line 241 "zher2.f"
		    d_cnjg(&z__1, &z__2);
#line 241 "zher2.f"
		    temp2.r = z__1.r, temp2.i = z__1.i;
#line 242 "zher2.f"
		    i__2 = j - 1;
#line 242 "zher2.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 243 "zher2.f"
			i__3 = i__ + j * a_dim1;
#line 243 "zher2.f"
			i__4 = i__ + j * a_dim1;
#line 243 "zher2.f"
			i__5 = i__;
#line 243 "zher2.f"
			z__3.r = x[i__5].r * temp1.r - x[i__5].i * temp1.i, 
				z__3.i = x[i__5].r * temp1.i + x[i__5].i * 
				temp1.r;
#line 243 "zher2.f"
			z__2.r = a[i__4].r + z__3.r, z__2.i = a[i__4].i + 
				z__3.i;
#line 243 "zher2.f"
			i__6 = i__;
#line 243 "zher2.f"
			z__4.r = y[i__6].r * temp2.r - y[i__6].i * temp2.i, 
				z__4.i = y[i__6].r * temp2.i + y[i__6].i * 
				temp2.r;
#line 243 "zher2.f"
			z__1.r = z__2.r + z__4.r, z__1.i = z__2.i + z__4.i;
#line 243 "zher2.f"
			a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 244 "zher2.f"
/* L10: */
#line 244 "zher2.f"
		    }
#line 245 "zher2.f"
		    i__2 = j + j * a_dim1;
#line 245 "zher2.f"
		    i__3 = j + j * a_dim1;
#line 245 "zher2.f"
		    i__4 = j;
#line 245 "zher2.f"
		    z__2.r = x[i__4].r * temp1.r - x[i__4].i * temp1.i, 
			    z__2.i = x[i__4].r * temp1.i + x[i__4].i * 
			    temp1.r;
#line 245 "zher2.f"
		    i__5 = j;
#line 245 "zher2.f"
		    z__3.r = y[i__5].r * temp2.r - y[i__5].i * temp2.i, 
			    z__3.i = y[i__5].r * temp2.i + y[i__5].i * 
			    temp2.r;
#line 245 "zher2.f"
		    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 245 "zher2.f"
		    d__1 = a[i__3].r + z__1.r;
#line 245 "zher2.f"
		    a[i__2].r = d__1, a[i__2].i = 0.;
#line 247 "zher2.f"
		} else {
#line 248 "zher2.f"
		    i__2 = j + j * a_dim1;
#line 248 "zher2.f"
		    i__3 = j + j * a_dim1;
#line 248 "zher2.f"
		    d__1 = a[i__3].r;
#line 248 "zher2.f"
		    a[i__2].r = d__1, a[i__2].i = 0.;
#line 249 "zher2.f"
		}
#line 250 "zher2.f"
/* L20: */
#line 250 "zher2.f"
	    }
#line 251 "zher2.f"
	} else {
#line 252 "zher2.f"
	    i__1 = *n;
#line 252 "zher2.f"
	    for (j = 1; j <= i__1; ++j) {
#line 253 "zher2.f"
		i__2 = jx;
#line 253 "zher2.f"
		i__3 = jy;
#line 253 "zher2.f"
		if (x[i__2].r != 0. || x[i__2].i != 0. || (y[i__3].r != 0. || 
			y[i__3].i != 0.)) {
#line 254 "zher2.f"
		    d_cnjg(&z__2, &y[jy]);
#line 254 "zher2.f"
		    z__1.r = alpha->r * z__2.r - alpha->i * z__2.i, z__1.i = 
			    alpha->r * z__2.i + alpha->i * z__2.r;
#line 254 "zher2.f"
		    temp1.r = z__1.r, temp1.i = z__1.i;
#line 255 "zher2.f"
		    i__2 = jx;
#line 255 "zher2.f"
		    z__2.r = alpha->r * x[i__2].r - alpha->i * x[i__2].i, 
			    z__2.i = alpha->r * x[i__2].i + alpha->i * x[i__2]
			    .r;
#line 255 "zher2.f"
		    d_cnjg(&z__1, &z__2);
#line 255 "zher2.f"
		    temp2.r = z__1.r, temp2.i = z__1.i;
#line 256 "zher2.f"
		    ix = kx;
#line 257 "zher2.f"
		    iy = ky;
#line 258 "zher2.f"
		    i__2 = j - 1;
#line 258 "zher2.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 259 "zher2.f"
			i__3 = i__ + j * a_dim1;
#line 259 "zher2.f"
			i__4 = i__ + j * a_dim1;
#line 259 "zher2.f"
			i__5 = ix;
#line 259 "zher2.f"
			z__3.r = x[i__5].r * temp1.r - x[i__5].i * temp1.i, 
				z__3.i = x[i__5].r * temp1.i + x[i__5].i * 
				temp1.r;
#line 259 "zher2.f"
			z__2.r = a[i__4].r + z__3.r, z__2.i = a[i__4].i + 
				z__3.i;
#line 259 "zher2.f"
			i__6 = iy;
#line 259 "zher2.f"
			z__4.r = y[i__6].r * temp2.r - y[i__6].i * temp2.i, 
				z__4.i = y[i__6].r * temp2.i + y[i__6].i * 
				temp2.r;
#line 259 "zher2.f"
			z__1.r = z__2.r + z__4.r, z__1.i = z__2.i + z__4.i;
#line 259 "zher2.f"
			a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 260 "zher2.f"
			ix += *incx;
#line 261 "zher2.f"
			iy += *incy;
#line 262 "zher2.f"
/* L30: */
#line 262 "zher2.f"
		    }
#line 263 "zher2.f"
		    i__2 = j + j * a_dim1;
#line 263 "zher2.f"
		    i__3 = j + j * a_dim1;
#line 263 "zher2.f"
		    i__4 = jx;
#line 263 "zher2.f"
		    z__2.r = x[i__4].r * temp1.r - x[i__4].i * temp1.i, 
			    z__2.i = x[i__4].r * temp1.i + x[i__4].i * 
			    temp1.r;
#line 263 "zher2.f"
		    i__5 = jy;
#line 263 "zher2.f"
		    z__3.r = y[i__5].r * temp2.r - y[i__5].i * temp2.i, 
			    z__3.i = y[i__5].r * temp2.i + y[i__5].i * 
			    temp2.r;
#line 263 "zher2.f"
		    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 263 "zher2.f"
		    d__1 = a[i__3].r + z__1.r;
#line 263 "zher2.f"
		    a[i__2].r = d__1, a[i__2].i = 0.;
#line 265 "zher2.f"
		} else {
#line 266 "zher2.f"
		    i__2 = j + j * a_dim1;
#line 266 "zher2.f"
		    i__3 = j + j * a_dim1;
#line 266 "zher2.f"
		    d__1 = a[i__3].r;
#line 266 "zher2.f"
		    a[i__2].r = d__1, a[i__2].i = 0.;
#line 267 "zher2.f"
		}
#line 268 "zher2.f"
		jx += *incx;
#line 269 "zher2.f"
		jy += *incy;
#line 270 "zher2.f"
/* L40: */
#line 270 "zher2.f"
	    }
#line 271 "zher2.f"
	}
#line 272 "zher2.f"
    } else {

/*        Form  A  when A is stored in the lower triangle. */

#line 276 "zher2.f"
	if (*incx == 1 && *incy == 1) {
#line 277 "zher2.f"
	    i__1 = *n;
#line 277 "zher2.f"
	    for (j = 1; j <= i__1; ++j) {
#line 278 "zher2.f"
		i__2 = j;
#line 278 "zher2.f"
		i__3 = j;
#line 278 "zher2.f"
		if (x[i__2].r != 0. || x[i__2].i != 0. || (y[i__3].r != 0. || 
			y[i__3].i != 0.)) {
#line 279 "zher2.f"
		    d_cnjg(&z__2, &y[j]);
#line 279 "zher2.f"
		    z__1.r = alpha->r * z__2.r - alpha->i * z__2.i, z__1.i = 
			    alpha->r * z__2.i + alpha->i * z__2.r;
#line 279 "zher2.f"
		    temp1.r = z__1.r, temp1.i = z__1.i;
#line 280 "zher2.f"
		    i__2 = j;
#line 280 "zher2.f"
		    z__2.r = alpha->r * x[i__2].r - alpha->i * x[i__2].i, 
			    z__2.i = alpha->r * x[i__2].i + alpha->i * x[i__2]
			    .r;
#line 280 "zher2.f"
		    d_cnjg(&z__1, &z__2);
#line 280 "zher2.f"
		    temp2.r = z__1.r, temp2.i = z__1.i;
#line 281 "zher2.f"
		    i__2 = j + j * a_dim1;
#line 281 "zher2.f"
		    i__3 = j + j * a_dim1;
#line 281 "zher2.f"
		    i__4 = j;
#line 281 "zher2.f"
		    z__2.r = x[i__4].r * temp1.r - x[i__4].i * temp1.i, 
			    z__2.i = x[i__4].r * temp1.i + x[i__4].i * 
			    temp1.r;
#line 281 "zher2.f"
		    i__5 = j;
#line 281 "zher2.f"
		    z__3.r = y[i__5].r * temp2.r - y[i__5].i * temp2.i, 
			    z__3.i = y[i__5].r * temp2.i + y[i__5].i * 
			    temp2.r;
#line 281 "zher2.f"
		    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 281 "zher2.f"
		    d__1 = a[i__3].r + z__1.r;
#line 281 "zher2.f"
		    a[i__2].r = d__1, a[i__2].i = 0.;
#line 283 "zher2.f"
		    i__2 = *n;
#line 283 "zher2.f"
		    for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 284 "zher2.f"
			i__3 = i__ + j * a_dim1;
#line 284 "zher2.f"
			i__4 = i__ + j * a_dim1;
#line 284 "zher2.f"
			i__5 = i__;
#line 284 "zher2.f"
			z__3.r = x[i__5].r * temp1.r - x[i__5].i * temp1.i, 
				z__3.i = x[i__5].r * temp1.i + x[i__5].i * 
				temp1.r;
#line 284 "zher2.f"
			z__2.r = a[i__4].r + z__3.r, z__2.i = a[i__4].i + 
				z__3.i;
#line 284 "zher2.f"
			i__6 = i__;
#line 284 "zher2.f"
			z__4.r = y[i__6].r * temp2.r - y[i__6].i * temp2.i, 
				z__4.i = y[i__6].r * temp2.i + y[i__6].i * 
				temp2.r;
#line 284 "zher2.f"
			z__1.r = z__2.r + z__4.r, z__1.i = z__2.i + z__4.i;
#line 284 "zher2.f"
			a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 285 "zher2.f"
/* L50: */
#line 285 "zher2.f"
		    }
#line 286 "zher2.f"
		} else {
#line 287 "zher2.f"
		    i__2 = j + j * a_dim1;
#line 287 "zher2.f"
		    i__3 = j + j * a_dim1;
#line 287 "zher2.f"
		    d__1 = a[i__3].r;
#line 287 "zher2.f"
		    a[i__2].r = d__1, a[i__2].i = 0.;
#line 288 "zher2.f"
		}
#line 289 "zher2.f"
/* L60: */
#line 289 "zher2.f"
	    }
#line 290 "zher2.f"
	} else {
#line 291 "zher2.f"
	    i__1 = *n;
#line 291 "zher2.f"
	    for (j = 1; j <= i__1; ++j) {
#line 292 "zher2.f"
		i__2 = jx;
#line 292 "zher2.f"
		i__3 = jy;
#line 292 "zher2.f"
		if (x[i__2].r != 0. || x[i__2].i != 0. || (y[i__3].r != 0. || 
			y[i__3].i != 0.)) {
#line 293 "zher2.f"
		    d_cnjg(&z__2, &y[jy]);
#line 293 "zher2.f"
		    z__1.r = alpha->r * z__2.r - alpha->i * z__2.i, z__1.i = 
			    alpha->r * z__2.i + alpha->i * z__2.r;
#line 293 "zher2.f"
		    temp1.r = z__1.r, temp1.i = z__1.i;
#line 294 "zher2.f"
		    i__2 = jx;
#line 294 "zher2.f"
		    z__2.r = alpha->r * x[i__2].r - alpha->i * x[i__2].i, 
			    z__2.i = alpha->r * x[i__2].i + alpha->i * x[i__2]
			    .r;
#line 294 "zher2.f"
		    d_cnjg(&z__1, &z__2);
#line 294 "zher2.f"
		    temp2.r = z__1.r, temp2.i = z__1.i;
#line 295 "zher2.f"
		    i__2 = j + j * a_dim1;
#line 295 "zher2.f"
		    i__3 = j + j * a_dim1;
#line 295 "zher2.f"
		    i__4 = jx;
#line 295 "zher2.f"
		    z__2.r = x[i__4].r * temp1.r - x[i__4].i * temp1.i, 
			    z__2.i = x[i__4].r * temp1.i + x[i__4].i * 
			    temp1.r;
#line 295 "zher2.f"
		    i__5 = jy;
#line 295 "zher2.f"
		    z__3.r = y[i__5].r * temp2.r - y[i__5].i * temp2.i, 
			    z__3.i = y[i__5].r * temp2.i + y[i__5].i * 
			    temp2.r;
#line 295 "zher2.f"
		    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 295 "zher2.f"
		    d__1 = a[i__3].r + z__1.r;
#line 295 "zher2.f"
		    a[i__2].r = d__1, a[i__2].i = 0.;
#line 297 "zher2.f"
		    ix = jx;
#line 298 "zher2.f"
		    iy = jy;
#line 299 "zher2.f"
		    i__2 = *n;
#line 299 "zher2.f"
		    for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 300 "zher2.f"
			ix += *incx;
#line 301 "zher2.f"
			iy += *incy;
#line 302 "zher2.f"
			i__3 = i__ + j * a_dim1;
#line 302 "zher2.f"
			i__4 = i__ + j * a_dim1;
#line 302 "zher2.f"
			i__5 = ix;
#line 302 "zher2.f"
			z__3.r = x[i__5].r * temp1.r - x[i__5].i * temp1.i, 
				z__3.i = x[i__5].r * temp1.i + x[i__5].i * 
				temp1.r;
#line 302 "zher2.f"
			z__2.r = a[i__4].r + z__3.r, z__2.i = a[i__4].i + 
				z__3.i;
#line 302 "zher2.f"
			i__6 = iy;
#line 302 "zher2.f"
			z__4.r = y[i__6].r * temp2.r - y[i__6].i * temp2.i, 
				z__4.i = y[i__6].r * temp2.i + y[i__6].i * 
				temp2.r;
#line 302 "zher2.f"
			z__1.r = z__2.r + z__4.r, z__1.i = z__2.i + z__4.i;
#line 302 "zher2.f"
			a[i__3].r = z__1.r, a[i__3].i = z__1.i;
#line 303 "zher2.f"
/* L70: */
#line 303 "zher2.f"
		    }
#line 304 "zher2.f"
		} else {
#line 305 "zher2.f"
		    i__2 = j + j * a_dim1;
#line 305 "zher2.f"
		    i__3 = j + j * a_dim1;
#line 305 "zher2.f"
		    d__1 = a[i__3].r;
#line 305 "zher2.f"
		    a[i__2].r = d__1, a[i__2].i = 0.;
#line 306 "zher2.f"
		}
#line 307 "zher2.f"
		jx += *incx;
#line 308 "zher2.f"
		jy += *incy;
#line 309 "zher2.f"
/* L80: */
#line 309 "zher2.f"
	    }
#line 310 "zher2.f"
	}
#line 311 "zher2.f"
    }

#line 313 "zher2.f"
    return 0;

/*     End of ZHER2 . */

} /* zher2_ */

