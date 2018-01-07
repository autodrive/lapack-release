#line 1 "chpr2.f"
/* chpr2.f -- translated by f2c (version 20100827).
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

#line 1 "chpr2.f"
/* > \brief \b CHPR2 */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CHPR2(UPLO,N,ALPHA,X,INCX,Y,INCY,AP) */

/*       .. Scalar Arguments .. */
/*       COMPLEX ALPHA */
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
/* > CHPR2  performs the hermitian rank 2 operation */
/* > */
/* >    A := alpha*x*y**H + conjg( alpha )*y*x**H + A, */
/* > */
/* > where alpha is a scalar, x and y are n element vectors and A is an */
/* > n by n hermitian matrix, supplied in packed form. */
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
/* > \param[in] X */
/* > \verbatim */
/* >          X is COMPLEX array, dimension at least */
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
/* >          Y is COMPLEX array, dimension at least */
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
/* > \param[in,out] AP */
/* > \verbatim */
/* >          AP is COMPLEX array, dimension at least */
/* >           ( ( n*( n + 1 ) )/2 ). */
/* >           Before entry with  UPLO = 'U' or 'u', the array AP must */
/* >           contain the upper triangular part of the hermitian matrix */
/* >           packed sequentially, column by column, so that AP( 1 ) */
/* >           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 ) */
/* >           and a( 2, 2 ) respectively, and so on. On exit, the array */
/* >           AP is overwritten by the upper triangular part of the */
/* >           updated matrix. */
/* >           Before entry with UPLO = 'L' or 'l', the array AP must */
/* >           contain the lower triangular part of the hermitian matrix */
/* >           packed sequentially, column by column, so that AP( 1 ) */
/* >           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 ) */
/* >           and a( 3, 1 ) respectively, and so on. On exit, the array */
/* >           AP is overwritten by the lower triangular part of the */
/* >           updated matrix. */
/* >           Note that the imaginary parts of the diagonal elements need */
/* >           not be set, they are assumed to be zero, and on exit they */
/* >           are set to zero. */
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
/* > */
/* >  -- Written on 22-October-1986. */
/* >     Jack Dongarra, Argonne National Lab. */
/* >     Jeremy Du Croz, Nag Central Office. */
/* >     Sven Hammarling, Nag Central Office. */
/* >     Richard Hanson, Sandia National Labs. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int chpr2_(char *uplo, integer *n, doublecomplex *alpha, 
	doublecomplex *x, integer *incx, doublecomplex *y, integer *incy, 
	doublecomplex *ap, ftnlen uplo_len)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5, i__6;
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

#line 185 "chpr2.f"
    /* Parameter adjustments */
#line 185 "chpr2.f"
    --ap;
#line 185 "chpr2.f"
    --y;
#line 185 "chpr2.f"
    --x;
#line 185 "chpr2.f"

#line 185 "chpr2.f"
    /* Function Body */
#line 185 "chpr2.f"
    info = 0;
#line 186 "chpr2.f"
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
#line 187 "chpr2.f"
	info = 1;
#line 188 "chpr2.f"
    } else if (*n < 0) {
#line 189 "chpr2.f"
	info = 2;
#line 190 "chpr2.f"
    } else if (*incx == 0) {
#line 191 "chpr2.f"
	info = 5;
#line 192 "chpr2.f"
    } else if (*incy == 0) {
#line 193 "chpr2.f"
	info = 7;
#line 194 "chpr2.f"
    }
#line 195 "chpr2.f"
    if (info != 0) {
#line 196 "chpr2.f"
	xerbla_("CHPR2 ", &info, (ftnlen)6);
#line 197 "chpr2.f"
	return 0;
#line 198 "chpr2.f"
    }

/*     Quick return if possible. */

#line 202 "chpr2.f"
    if (*n == 0 || alpha->r == 0. && alpha->i == 0.) {
#line 202 "chpr2.f"
	return 0;
#line 202 "chpr2.f"
    }

/*     Set up the start points in X and Y if the increments are not both */
/*     unity. */

#line 207 "chpr2.f"
    if (*incx != 1 || *incy != 1) {
#line 208 "chpr2.f"
	if (*incx > 0) {
#line 209 "chpr2.f"
	    kx = 1;
#line 210 "chpr2.f"
	} else {
#line 211 "chpr2.f"
	    kx = 1 - (*n - 1) * *incx;
#line 212 "chpr2.f"
	}
#line 213 "chpr2.f"
	if (*incy > 0) {
#line 214 "chpr2.f"
	    ky = 1;
#line 215 "chpr2.f"
	} else {
#line 216 "chpr2.f"
	    ky = 1 - (*n - 1) * *incy;
#line 217 "chpr2.f"
	}
#line 218 "chpr2.f"
	jx = kx;
#line 219 "chpr2.f"
	jy = ky;
#line 220 "chpr2.f"
    }

/*     Start the operations. In this version the elements of the array AP */
/*     are accessed sequentially with one pass through AP. */

#line 225 "chpr2.f"
    kk = 1;
#line 226 "chpr2.f"
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {

/*        Form  A  when upper triangle is stored in AP. */

#line 230 "chpr2.f"
	if (*incx == 1 && *incy == 1) {
#line 231 "chpr2.f"
	    i__1 = *n;
#line 231 "chpr2.f"
	    for (j = 1; j <= i__1; ++j) {
#line 232 "chpr2.f"
		i__2 = j;
#line 232 "chpr2.f"
		i__3 = j;
#line 232 "chpr2.f"
		if (x[i__2].r != 0. || x[i__2].i != 0. || (y[i__3].r != 0. || 
			y[i__3].i != 0.)) {
#line 233 "chpr2.f"
		    d_cnjg(&z__2, &y[j]);
#line 233 "chpr2.f"
		    z__1.r = alpha->r * z__2.r - alpha->i * z__2.i, z__1.i = 
			    alpha->r * z__2.i + alpha->i * z__2.r;
#line 233 "chpr2.f"
		    temp1.r = z__1.r, temp1.i = z__1.i;
#line 234 "chpr2.f"
		    i__2 = j;
#line 234 "chpr2.f"
		    z__2.r = alpha->r * x[i__2].r - alpha->i * x[i__2].i, 
			    z__2.i = alpha->r * x[i__2].i + alpha->i * x[i__2]
			    .r;
#line 234 "chpr2.f"
		    d_cnjg(&z__1, &z__2);
#line 234 "chpr2.f"
		    temp2.r = z__1.r, temp2.i = z__1.i;
#line 235 "chpr2.f"
		    k = kk;
#line 236 "chpr2.f"
		    i__2 = j - 1;
#line 236 "chpr2.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 237 "chpr2.f"
			i__3 = k;
#line 237 "chpr2.f"
			i__4 = k;
#line 237 "chpr2.f"
			i__5 = i__;
#line 237 "chpr2.f"
			z__3.r = x[i__5].r * temp1.r - x[i__5].i * temp1.i, 
				z__3.i = x[i__5].r * temp1.i + x[i__5].i * 
				temp1.r;
#line 237 "chpr2.f"
			z__2.r = ap[i__4].r + z__3.r, z__2.i = ap[i__4].i + 
				z__3.i;
#line 237 "chpr2.f"
			i__6 = i__;
#line 237 "chpr2.f"
			z__4.r = y[i__6].r * temp2.r - y[i__6].i * temp2.i, 
				z__4.i = y[i__6].r * temp2.i + y[i__6].i * 
				temp2.r;
#line 237 "chpr2.f"
			z__1.r = z__2.r + z__4.r, z__1.i = z__2.i + z__4.i;
#line 237 "chpr2.f"
			ap[i__3].r = z__1.r, ap[i__3].i = z__1.i;
#line 238 "chpr2.f"
			++k;
#line 239 "chpr2.f"
/* L10: */
#line 239 "chpr2.f"
		    }
#line 240 "chpr2.f"
		    i__2 = kk + j - 1;
#line 240 "chpr2.f"
		    i__3 = kk + j - 1;
#line 240 "chpr2.f"
		    i__4 = j;
#line 240 "chpr2.f"
		    z__2.r = x[i__4].r * temp1.r - x[i__4].i * temp1.i, 
			    z__2.i = x[i__4].r * temp1.i + x[i__4].i * 
			    temp1.r;
#line 240 "chpr2.f"
		    i__5 = j;
#line 240 "chpr2.f"
		    z__3.r = y[i__5].r * temp2.r - y[i__5].i * temp2.i, 
			    z__3.i = y[i__5].r * temp2.i + y[i__5].i * 
			    temp2.r;
#line 240 "chpr2.f"
		    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 240 "chpr2.f"
		    d__1 = ap[i__3].r + z__1.r;
#line 240 "chpr2.f"
		    ap[i__2].r = d__1, ap[i__2].i = 0.;
#line 242 "chpr2.f"
		} else {
#line 243 "chpr2.f"
		    i__2 = kk + j - 1;
#line 243 "chpr2.f"
		    i__3 = kk + j - 1;
#line 243 "chpr2.f"
		    d__1 = ap[i__3].r;
#line 243 "chpr2.f"
		    ap[i__2].r = d__1, ap[i__2].i = 0.;
#line 244 "chpr2.f"
		}
#line 245 "chpr2.f"
		kk += j;
#line 246 "chpr2.f"
/* L20: */
#line 246 "chpr2.f"
	    }
#line 247 "chpr2.f"
	} else {
#line 248 "chpr2.f"
	    i__1 = *n;
#line 248 "chpr2.f"
	    for (j = 1; j <= i__1; ++j) {
#line 249 "chpr2.f"
		i__2 = jx;
#line 249 "chpr2.f"
		i__3 = jy;
#line 249 "chpr2.f"
		if (x[i__2].r != 0. || x[i__2].i != 0. || (y[i__3].r != 0. || 
			y[i__3].i != 0.)) {
#line 250 "chpr2.f"
		    d_cnjg(&z__2, &y[jy]);
#line 250 "chpr2.f"
		    z__1.r = alpha->r * z__2.r - alpha->i * z__2.i, z__1.i = 
			    alpha->r * z__2.i + alpha->i * z__2.r;
#line 250 "chpr2.f"
		    temp1.r = z__1.r, temp1.i = z__1.i;
#line 251 "chpr2.f"
		    i__2 = jx;
#line 251 "chpr2.f"
		    z__2.r = alpha->r * x[i__2].r - alpha->i * x[i__2].i, 
			    z__2.i = alpha->r * x[i__2].i + alpha->i * x[i__2]
			    .r;
#line 251 "chpr2.f"
		    d_cnjg(&z__1, &z__2);
#line 251 "chpr2.f"
		    temp2.r = z__1.r, temp2.i = z__1.i;
#line 252 "chpr2.f"
		    ix = kx;
#line 253 "chpr2.f"
		    iy = ky;
#line 254 "chpr2.f"
		    i__2 = kk + j - 2;
#line 254 "chpr2.f"
		    for (k = kk; k <= i__2; ++k) {
#line 255 "chpr2.f"
			i__3 = k;
#line 255 "chpr2.f"
			i__4 = k;
#line 255 "chpr2.f"
			i__5 = ix;
#line 255 "chpr2.f"
			z__3.r = x[i__5].r * temp1.r - x[i__5].i * temp1.i, 
				z__3.i = x[i__5].r * temp1.i + x[i__5].i * 
				temp1.r;
#line 255 "chpr2.f"
			z__2.r = ap[i__4].r + z__3.r, z__2.i = ap[i__4].i + 
				z__3.i;
#line 255 "chpr2.f"
			i__6 = iy;
#line 255 "chpr2.f"
			z__4.r = y[i__6].r * temp2.r - y[i__6].i * temp2.i, 
				z__4.i = y[i__6].r * temp2.i + y[i__6].i * 
				temp2.r;
#line 255 "chpr2.f"
			z__1.r = z__2.r + z__4.r, z__1.i = z__2.i + z__4.i;
#line 255 "chpr2.f"
			ap[i__3].r = z__1.r, ap[i__3].i = z__1.i;
#line 256 "chpr2.f"
			ix += *incx;
#line 257 "chpr2.f"
			iy += *incy;
#line 258 "chpr2.f"
/* L30: */
#line 258 "chpr2.f"
		    }
#line 259 "chpr2.f"
		    i__2 = kk + j - 1;
#line 259 "chpr2.f"
		    i__3 = kk + j - 1;
#line 259 "chpr2.f"
		    i__4 = jx;
#line 259 "chpr2.f"
		    z__2.r = x[i__4].r * temp1.r - x[i__4].i * temp1.i, 
			    z__2.i = x[i__4].r * temp1.i + x[i__4].i * 
			    temp1.r;
#line 259 "chpr2.f"
		    i__5 = jy;
#line 259 "chpr2.f"
		    z__3.r = y[i__5].r * temp2.r - y[i__5].i * temp2.i, 
			    z__3.i = y[i__5].r * temp2.i + y[i__5].i * 
			    temp2.r;
#line 259 "chpr2.f"
		    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 259 "chpr2.f"
		    d__1 = ap[i__3].r + z__1.r;
#line 259 "chpr2.f"
		    ap[i__2].r = d__1, ap[i__2].i = 0.;
#line 261 "chpr2.f"
		} else {
#line 262 "chpr2.f"
		    i__2 = kk + j - 1;
#line 262 "chpr2.f"
		    i__3 = kk + j - 1;
#line 262 "chpr2.f"
		    d__1 = ap[i__3].r;
#line 262 "chpr2.f"
		    ap[i__2].r = d__1, ap[i__2].i = 0.;
#line 263 "chpr2.f"
		}
#line 264 "chpr2.f"
		jx += *incx;
#line 265 "chpr2.f"
		jy += *incy;
#line 266 "chpr2.f"
		kk += j;
#line 267 "chpr2.f"
/* L40: */
#line 267 "chpr2.f"
	    }
#line 268 "chpr2.f"
	}
#line 269 "chpr2.f"
    } else {

/*        Form  A  when lower triangle is stored in AP. */

#line 273 "chpr2.f"
	if (*incx == 1 && *incy == 1) {
#line 274 "chpr2.f"
	    i__1 = *n;
#line 274 "chpr2.f"
	    for (j = 1; j <= i__1; ++j) {
#line 275 "chpr2.f"
		i__2 = j;
#line 275 "chpr2.f"
		i__3 = j;
#line 275 "chpr2.f"
		if (x[i__2].r != 0. || x[i__2].i != 0. || (y[i__3].r != 0. || 
			y[i__3].i != 0.)) {
#line 276 "chpr2.f"
		    d_cnjg(&z__2, &y[j]);
#line 276 "chpr2.f"
		    z__1.r = alpha->r * z__2.r - alpha->i * z__2.i, z__1.i = 
			    alpha->r * z__2.i + alpha->i * z__2.r;
#line 276 "chpr2.f"
		    temp1.r = z__1.r, temp1.i = z__1.i;
#line 277 "chpr2.f"
		    i__2 = j;
#line 277 "chpr2.f"
		    z__2.r = alpha->r * x[i__2].r - alpha->i * x[i__2].i, 
			    z__2.i = alpha->r * x[i__2].i + alpha->i * x[i__2]
			    .r;
#line 277 "chpr2.f"
		    d_cnjg(&z__1, &z__2);
#line 277 "chpr2.f"
		    temp2.r = z__1.r, temp2.i = z__1.i;
#line 278 "chpr2.f"
		    i__2 = kk;
#line 278 "chpr2.f"
		    i__3 = kk;
#line 278 "chpr2.f"
		    i__4 = j;
#line 278 "chpr2.f"
		    z__2.r = x[i__4].r * temp1.r - x[i__4].i * temp1.i, 
			    z__2.i = x[i__4].r * temp1.i + x[i__4].i * 
			    temp1.r;
#line 278 "chpr2.f"
		    i__5 = j;
#line 278 "chpr2.f"
		    z__3.r = y[i__5].r * temp2.r - y[i__5].i * temp2.i, 
			    z__3.i = y[i__5].r * temp2.i + y[i__5].i * 
			    temp2.r;
#line 278 "chpr2.f"
		    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 278 "chpr2.f"
		    d__1 = ap[i__3].r + z__1.r;
#line 278 "chpr2.f"
		    ap[i__2].r = d__1, ap[i__2].i = 0.;
#line 280 "chpr2.f"
		    k = kk + 1;
#line 281 "chpr2.f"
		    i__2 = *n;
#line 281 "chpr2.f"
		    for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 282 "chpr2.f"
			i__3 = k;
#line 282 "chpr2.f"
			i__4 = k;
#line 282 "chpr2.f"
			i__5 = i__;
#line 282 "chpr2.f"
			z__3.r = x[i__5].r * temp1.r - x[i__5].i * temp1.i, 
				z__3.i = x[i__5].r * temp1.i + x[i__5].i * 
				temp1.r;
#line 282 "chpr2.f"
			z__2.r = ap[i__4].r + z__3.r, z__2.i = ap[i__4].i + 
				z__3.i;
#line 282 "chpr2.f"
			i__6 = i__;
#line 282 "chpr2.f"
			z__4.r = y[i__6].r * temp2.r - y[i__6].i * temp2.i, 
				z__4.i = y[i__6].r * temp2.i + y[i__6].i * 
				temp2.r;
#line 282 "chpr2.f"
			z__1.r = z__2.r + z__4.r, z__1.i = z__2.i + z__4.i;
#line 282 "chpr2.f"
			ap[i__3].r = z__1.r, ap[i__3].i = z__1.i;
#line 283 "chpr2.f"
			++k;
#line 284 "chpr2.f"
/* L50: */
#line 284 "chpr2.f"
		    }
#line 285 "chpr2.f"
		} else {
#line 286 "chpr2.f"
		    i__2 = kk;
#line 286 "chpr2.f"
		    i__3 = kk;
#line 286 "chpr2.f"
		    d__1 = ap[i__3].r;
#line 286 "chpr2.f"
		    ap[i__2].r = d__1, ap[i__2].i = 0.;
#line 287 "chpr2.f"
		}
#line 288 "chpr2.f"
		kk = kk + *n - j + 1;
#line 289 "chpr2.f"
/* L60: */
#line 289 "chpr2.f"
	    }
#line 290 "chpr2.f"
	} else {
#line 291 "chpr2.f"
	    i__1 = *n;
#line 291 "chpr2.f"
	    for (j = 1; j <= i__1; ++j) {
#line 292 "chpr2.f"
		i__2 = jx;
#line 292 "chpr2.f"
		i__3 = jy;
#line 292 "chpr2.f"
		if (x[i__2].r != 0. || x[i__2].i != 0. || (y[i__3].r != 0. || 
			y[i__3].i != 0.)) {
#line 293 "chpr2.f"
		    d_cnjg(&z__2, &y[jy]);
#line 293 "chpr2.f"
		    z__1.r = alpha->r * z__2.r - alpha->i * z__2.i, z__1.i = 
			    alpha->r * z__2.i + alpha->i * z__2.r;
#line 293 "chpr2.f"
		    temp1.r = z__1.r, temp1.i = z__1.i;
#line 294 "chpr2.f"
		    i__2 = jx;
#line 294 "chpr2.f"
		    z__2.r = alpha->r * x[i__2].r - alpha->i * x[i__2].i, 
			    z__2.i = alpha->r * x[i__2].i + alpha->i * x[i__2]
			    .r;
#line 294 "chpr2.f"
		    d_cnjg(&z__1, &z__2);
#line 294 "chpr2.f"
		    temp2.r = z__1.r, temp2.i = z__1.i;
#line 295 "chpr2.f"
		    i__2 = kk;
#line 295 "chpr2.f"
		    i__3 = kk;
#line 295 "chpr2.f"
		    i__4 = jx;
#line 295 "chpr2.f"
		    z__2.r = x[i__4].r * temp1.r - x[i__4].i * temp1.i, 
			    z__2.i = x[i__4].r * temp1.i + x[i__4].i * 
			    temp1.r;
#line 295 "chpr2.f"
		    i__5 = jy;
#line 295 "chpr2.f"
		    z__3.r = y[i__5].r * temp2.r - y[i__5].i * temp2.i, 
			    z__3.i = y[i__5].r * temp2.i + y[i__5].i * 
			    temp2.r;
#line 295 "chpr2.f"
		    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 295 "chpr2.f"
		    d__1 = ap[i__3].r + z__1.r;
#line 295 "chpr2.f"
		    ap[i__2].r = d__1, ap[i__2].i = 0.;
#line 297 "chpr2.f"
		    ix = jx;
#line 298 "chpr2.f"
		    iy = jy;
#line 299 "chpr2.f"
		    i__2 = kk + *n - j;
#line 299 "chpr2.f"
		    for (k = kk + 1; k <= i__2; ++k) {
#line 300 "chpr2.f"
			ix += *incx;
#line 301 "chpr2.f"
			iy += *incy;
#line 302 "chpr2.f"
			i__3 = k;
#line 302 "chpr2.f"
			i__4 = k;
#line 302 "chpr2.f"
			i__5 = ix;
#line 302 "chpr2.f"
			z__3.r = x[i__5].r * temp1.r - x[i__5].i * temp1.i, 
				z__3.i = x[i__5].r * temp1.i + x[i__5].i * 
				temp1.r;
#line 302 "chpr2.f"
			z__2.r = ap[i__4].r + z__3.r, z__2.i = ap[i__4].i + 
				z__3.i;
#line 302 "chpr2.f"
			i__6 = iy;
#line 302 "chpr2.f"
			z__4.r = y[i__6].r * temp2.r - y[i__6].i * temp2.i, 
				z__4.i = y[i__6].r * temp2.i + y[i__6].i * 
				temp2.r;
#line 302 "chpr2.f"
			z__1.r = z__2.r + z__4.r, z__1.i = z__2.i + z__4.i;
#line 302 "chpr2.f"
			ap[i__3].r = z__1.r, ap[i__3].i = z__1.i;
#line 303 "chpr2.f"
/* L70: */
#line 303 "chpr2.f"
		    }
#line 304 "chpr2.f"
		} else {
#line 305 "chpr2.f"
		    i__2 = kk;
#line 305 "chpr2.f"
		    i__3 = kk;
#line 305 "chpr2.f"
		    d__1 = ap[i__3].r;
#line 305 "chpr2.f"
		    ap[i__2].r = d__1, ap[i__2].i = 0.;
#line 306 "chpr2.f"
		}
#line 307 "chpr2.f"
		jx += *incx;
#line 308 "chpr2.f"
		jy += *incy;
#line 309 "chpr2.f"
		kk = kk + *n - j + 1;
#line 310 "chpr2.f"
/* L80: */
#line 310 "chpr2.f"
	    }
#line 311 "chpr2.f"
	}
#line 312 "chpr2.f"
    }

#line 314 "chpr2.f"
    return 0;

/*     End of CHPR2 . */

} /* chpr2_ */

