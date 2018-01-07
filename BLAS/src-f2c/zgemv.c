#line 1 "zgemv.f"
/* zgemv.f -- translated by f2c (version 20100827).
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

#line 1 "zgemv.f"
/* > \brief \b ZGEMV */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY) */

/*       .. Scalar Arguments .. */
/*       COMPLEX*16 ALPHA,BETA */
/*       INTEGER INCX,INCY,LDA,M,N */
/*       CHARACTER TRANS */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16 A(LDA,*),X(*),Y(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZGEMV  performs one of the matrix-vector operations */
/* > */
/* >    y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,   or */
/* > */
/* >    y := alpha*A**H*x + beta*y, */
/* > */
/* > where alpha and beta are scalars, x and y are vectors and A is an */
/* > m by n matrix. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] TRANS */
/* > \verbatim */
/* >          TRANS is CHARACTER*1 */
/* >           On entry, TRANS specifies the operation to be performed as */
/* >           follows: */
/* > */
/* >              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y. */
/* > */
/* >              TRANS = 'T' or 't'   y := alpha*A**T*x + beta*y. */
/* > */
/* >              TRANS = 'C' or 'c'   y := alpha*A**H*x + beta*y. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >           On entry, M specifies the number of rows of the matrix A. */
/* >           M must be at least zero. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >           On entry, N specifies the number of columns of the matrix A. */
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
/* >           Before entry, the leading m by n part of the array A must */
/* >           contain the matrix of coefficients. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >           On entry, LDA specifies the first dimension of A as declared */
/* >           in the calling (sub) program. LDA must be at least */
/* >           max( 1, m ). */
/* > \endverbatim */
/* > */
/* > \param[in] X */
/* > \verbatim */
/* >          X is COMPLEX*16 array of DIMENSION at least */
/* >           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n' */
/* >           and at least */
/* >           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise. */
/* >           Before entry, the incremented array X must contain the */
/* >           vector x. */
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
/* >          Y is COMPLEX*16 array of DIMENSION at least */
/* >           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n' */
/* >           and at least */
/* >           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise. */
/* >           Before entry with BETA non-zero, the incremented array Y */
/* >           must contain the vector y. On exit, Y is overwritten by the */
/* >           updated vector y. */
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

/* > \date November 2015 */

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
/* Subroutine */ int zgemv_(char *trans, integer *m, integer *n, 
	doublecomplex *alpha, doublecomplex *a, integer *lda, doublecomplex *
	x, integer *incx, doublecomplex *beta, doublecomplex *y, integer *
	incy, ftnlen trans_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, j, ix, iy, jx, jy, kx, ky, info;
    static doublecomplex temp;
    static integer lenx, leny;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static logical noconj;


/*  -- Reference BLAS level2 routine (version 3.6.0) -- */
/*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2015 */

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

#line 201 "zgemv.f"
    /* Parameter adjustments */
#line 201 "zgemv.f"
    a_dim1 = *lda;
#line 201 "zgemv.f"
    a_offset = 1 + a_dim1;
#line 201 "zgemv.f"
    a -= a_offset;
#line 201 "zgemv.f"
    --x;
#line 201 "zgemv.f"
    --y;
#line 201 "zgemv.f"

#line 201 "zgemv.f"
    /* Function Body */
#line 201 "zgemv.f"
    info = 0;
#line 202 "zgemv.f"
    if (! lsame_(trans, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, "T", (
	    ftnlen)1, (ftnlen)1) && ! lsame_(trans, "C", (ftnlen)1, (ftnlen)1)
	    ) {
#line 204 "zgemv.f"
	info = 1;
#line 205 "zgemv.f"
    } else if (*m < 0) {
#line 206 "zgemv.f"
	info = 2;
#line 207 "zgemv.f"
    } else if (*n < 0) {
#line 208 "zgemv.f"
	info = 3;
#line 209 "zgemv.f"
    } else if (*lda < max(1,*m)) {
#line 210 "zgemv.f"
	info = 6;
#line 211 "zgemv.f"
    } else if (*incx == 0) {
#line 212 "zgemv.f"
	info = 8;
#line 213 "zgemv.f"
    } else if (*incy == 0) {
#line 214 "zgemv.f"
	info = 11;
#line 215 "zgemv.f"
    }
#line 216 "zgemv.f"
    if (info != 0) {
#line 217 "zgemv.f"
	xerbla_("ZGEMV ", &info, (ftnlen)6);
#line 218 "zgemv.f"
	return 0;
#line 219 "zgemv.f"
    }

/*     Quick return if possible. */

#line 223 "zgemv.f"
    if (*m == 0 || *n == 0 || alpha->r == 0. && alpha->i == 0. && (beta->r == 
	    1. && beta->i == 0.)) {
#line 223 "zgemv.f"
	return 0;
#line 223 "zgemv.f"
    }

#line 226 "zgemv.f"
    noconj = lsame_(trans, "T", (ftnlen)1, (ftnlen)1);

/*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set */
/*     up the start points in  X  and  Y. */

#line 231 "zgemv.f"
    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {
#line 232 "zgemv.f"
	lenx = *n;
#line 233 "zgemv.f"
	leny = *m;
#line 234 "zgemv.f"
    } else {
#line 235 "zgemv.f"
	lenx = *m;
#line 236 "zgemv.f"
	leny = *n;
#line 237 "zgemv.f"
    }
#line 238 "zgemv.f"
    if (*incx > 0) {
#line 239 "zgemv.f"
	kx = 1;
#line 240 "zgemv.f"
    } else {
#line 241 "zgemv.f"
	kx = 1 - (lenx - 1) * *incx;
#line 242 "zgemv.f"
    }
#line 243 "zgemv.f"
    if (*incy > 0) {
#line 244 "zgemv.f"
	ky = 1;
#line 245 "zgemv.f"
    } else {
#line 246 "zgemv.f"
	ky = 1 - (leny - 1) * *incy;
#line 247 "zgemv.f"
    }

/*     Start the operations. In this version the elements of A are */
/*     accessed sequentially with one pass through A. */

/*     First form  y := beta*y. */

#line 254 "zgemv.f"
    if (beta->r != 1. || beta->i != 0.) {
#line 255 "zgemv.f"
	if (*incy == 1) {
#line 256 "zgemv.f"
	    if (beta->r == 0. && beta->i == 0.) {
#line 257 "zgemv.f"
		i__1 = leny;
#line 257 "zgemv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 258 "zgemv.f"
		    i__2 = i__;
#line 258 "zgemv.f"
		    y[i__2].r = 0., y[i__2].i = 0.;
#line 259 "zgemv.f"
/* L10: */
#line 259 "zgemv.f"
		}
#line 260 "zgemv.f"
	    } else {
#line 261 "zgemv.f"
		i__1 = leny;
#line 261 "zgemv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 262 "zgemv.f"
		    i__2 = i__;
#line 262 "zgemv.f"
		    i__3 = i__;
#line 262 "zgemv.f"
		    z__1.r = beta->r * y[i__3].r - beta->i * y[i__3].i, 
			    z__1.i = beta->r * y[i__3].i + beta->i * y[i__3]
			    .r;
#line 262 "zgemv.f"
		    y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 263 "zgemv.f"
/* L20: */
#line 263 "zgemv.f"
		}
#line 264 "zgemv.f"
	    }
#line 265 "zgemv.f"
	} else {
#line 266 "zgemv.f"
	    iy = ky;
#line 267 "zgemv.f"
	    if (beta->r == 0. && beta->i == 0.) {
#line 268 "zgemv.f"
		i__1 = leny;
#line 268 "zgemv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 269 "zgemv.f"
		    i__2 = iy;
#line 269 "zgemv.f"
		    y[i__2].r = 0., y[i__2].i = 0.;
#line 270 "zgemv.f"
		    iy += *incy;
#line 271 "zgemv.f"
/* L30: */
#line 271 "zgemv.f"
		}
#line 272 "zgemv.f"
	    } else {
#line 273 "zgemv.f"
		i__1 = leny;
#line 273 "zgemv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 274 "zgemv.f"
		    i__2 = iy;
#line 274 "zgemv.f"
		    i__3 = iy;
#line 274 "zgemv.f"
		    z__1.r = beta->r * y[i__3].r - beta->i * y[i__3].i, 
			    z__1.i = beta->r * y[i__3].i + beta->i * y[i__3]
			    .r;
#line 274 "zgemv.f"
		    y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 275 "zgemv.f"
		    iy += *incy;
#line 276 "zgemv.f"
/* L40: */
#line 276 "zgemv.f"
		}
#line 277 "zgemv.f"
	    }
#line 278 "zgemv.f"
	}
#line 279 "zgemv.f"
    }
#line 280 "zgemv.f"
    if (alpha->r == 0. && alpha->i == 0.) {
#line 280 "zgemv.f"
	return 0;
#line 280 "zgemv.f"
    }
#line 281 "zgemv.f"
    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {

/*        Form  y := alpha*A*x + y. */

#line 285 "zgemv.f"
	jx = kx;
#line 286 "zgemv.f"
	if (*incy == 1) {
#line 287 "zgemv.f"
	    i__1 = *n;
#line 287 "zgemv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 288 "zgemv.f"
		i__2 = jx;
#line 288 "zgemv.f"
		z__1.r = alpha->r * x[i__2].r - alpha->i * x[i__2].i, z__1.i =
			 alpha->r * x[i__2].i + alpha->i * x[i__2].r;
#line 288 "zgemv.f"
		temp.r = z__1.r, temp.i = z__1.i;
#line 289 "zgemv.f"
		i__2 = *m;
#line 289 "zgemv.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 290 "zgemv.f"
		    i__3 = i__;
#line 290 "zgemv.f"
		    i__4 = i__;
#line 290 "zgemv.f"
		    i__5 = i__ + j * a_dim1;
#line 290 "zgemv.f"
		    z__2.r = temp.r * a[i__5].r - temp.i * a[i__5].i, z__2.i =
			     temp.r * a[i__5].i + temp.i * a[i__5].r;
#line 290 "zgemv.f"
		    z__1.r = y[i__4].r + z__2.r, z__1.i = y[i__4].i + z__2.i;
#line 290 "zgemv.f"
		    y[i__3].r = z__1.r, y[i__3].i = z__1.i;
#line 291 "zgemv.f"
/* L50: */
#line 291 "zgemv.f"
		}
#line 292 "zgemv.f"
		jx += *incx;
#line 293 "zgemv.f"
/* L60: */
#line 293 "zgemv.f"
	    }
#line 294 "zgemv.f"
	} else {
#line 295 "zgemv.f"
	    i__1 = *n;
#line 295 "zgemv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 296 "zgemv.f"
		i__2 = jx;
#line 296 "zgemv.f"
		z__1.r = alpha->r * x[i__2].r - alpha->i * x[i__2].i, z__1.i =
			 alpha->r * x[i__2].i + alpha->i * x[i__2].r;
#line 296 "zgemv.f"
		temp.r = z__1.r, temp.i = z__1.i;
#line 297 "zgemv.f"
		iy = ky;
#line 298 "zgemv.f"
		i__2 = *m;
#line 298 "zgemv.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 299 "zgemv.f"
		    i__3 = iy;
#line 299 "zgemv.f"
		    i__4 = iy;
#line 299 "zgemv.f"
		    i__5 = i__ + j * a_dim1;
#line 299 "zgemv.f"
		    z__2.r = temp.r * a[i__5].r - temp.i * a[i__5].i, z__2.i =
			     temp.r * a[i__5].i + temp.i * a[i__5].r;
#line 299 "zgemv.f"
		    z__1.r = y[i__4].r + z__2.r, z__1.i = y[i__4].i + z__2.i;
#line 299 "zgemv.f"
		    y[i__3].r = z__1.r, y[i__3].i = z__1.i;
#line 300 "zgemv.f"
		    iy += *incy;
#line 301 "zgemv.f"
/* L70: */
#line 301 "zgemv.f"
		}
#line 302 "zgemv.f"
		jx += *incx;
#line 303 "zgemv.f"
/* L80: */
#line 303 "zgemv.f"
	    }
#line 304 "zgemv.f"
	}
#line 305 "zgemv.f"
    } else {

/*        Form  y := alpha*A**T*x + y  or  y := alpha*A**H*x + y. */

#line 309 "zgemv.f"
	jy = ky;
#line 310 "zgemv.f"
	if (*incx == 1) {
#line 311 "zgemv.f"
	    i__1 = *n;
#line 311 "zgemv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 312 "zgemv.f"
		temp.r = 0., temp.i = 0.;
#line 313 "zgemv.f"
		if (noconj) {
#line 314 "zgemv.f"
		    i__2 = *m;
#line 314 "zgemv.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 315 "zgemv.f"
			i__3 = i__ + j * a_dim1;
#line 315 "zgemv.f"
			i__4 = i__;
#line 315 "zgemv.f"
			z__2.r = a[i__3].r * x[i__4].r - a[i__3].i * x[i__4]
				.i, z__2.i = a[i__3].r * x[i__4].i + a[i__3]
				.i * x[i__4].r;
#line 315 "zgemv.f"
			z__1.r = temp.r + z__2.r, z__1.i = temp.i + z__2.i;
#line 315 "zgemv.f"
			temp.r = z__1.r, temp.i = z__1.i;
#line 316 "zgemv.f"
/* L90: */
#line 316 "zgemv.f"
		    }
#line 317 "zgemv.f"
		} else {
#line 318 "zgemv.f"
		    i__2 = *m;
#line 318 "zgemv.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 319 "zgemv.f"
			d_cnjg(&z__3, &a[i__ + j * a_dim1]);
#line 319 "zgemv.f"
			i__3 = i__;
#line 319 "zgemv.f"
			z__2.r = z__3.r * x[i__3].r - z__3.i * x[i__3].i, 
				z__2.i = z__3.r * x[i__3].i + z__3.i * x[i__3]
				.r;
#line 319 "zgemv.f"
			z__1.r = temp.r + z__2.r, z__1.i = temp.i + z__2.i;
#line 319 "zgemv.f"
			temp.r = z__1.r, temp.i = z__1.i;
#line 320 "zgemv.f"
/* L100: */
#line 320 "zgemv.f"
		    }
#line 321 "zgemv.f"
		}
#line 322 "zgemv.f"
		i__2 = jy;
#line 322 "zgemv.f"
		i__3 = jy;
#line 322 "zgemv.f"
		z__2.r = alpha->r * temp.r - alpha->i * temp.i, z__2.i = 
			alpha->r * temp.i + alpha->i * temp.r;
#line 322 "zgemv.f"
		z__1.r = y[i__3].r + z__2.r, z__1.i = y[i__3].i + z__2.i;
#line 322 "zgemv.f"
		y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 323 "zgemv.f"
		jy += *incy;
#line 324 "zgemv.f"
/* L110: */
#line 324 "zgemv.f"
	    }
#line 325 "zgemv.f"
	} else {
#line 326 "zgemv.f"
	    i__1 = *n;
#line 326 "zgemv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 327 "zgemv.f"
		temp.r = 0., temp.i = 0.;
#line 328 "zgemv.f"
		ix = kx;
#line 329 "zgemv.f"
		if (noconj) {
#line 330 "zgemv.f"
		    i__2 = *m;
#line 330 "zgemv.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 331 "zgemv.f"
			i__3 = i__ + j * a_dim1;
#line 331 "zgemv.f"
			i__4 = ix;
#line 331 "zgemv.f"
			z__2.r = a[i__3].r * x[i__4].r - a[i__3].i * x[i__4]
				.i, z__2.i = a[i__3].r * x[i__4].i + a[i__3]
				.i * x[i__4].r;
#line 331 "zgemv.f"
			z__1.r = temp.r + z__2.r, z__1.i = temp.i + z__2.i;
#line 331 "zgemv.f"
			temp.r = z__1.r, temp.i = z__1.i;
#line 332 "zgemv.f"
			ix += *incx;
#line 333 "zgemv.f"
/* L120: */
#line 333 "zgemv.f"
		    }
#line 334 "zgemv.f"
		} else {
#line 335 "zgemv.f"
		    i__2 = *m;
#line 335 "zgemv.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 336 "zgemv.f"
			d_cnjg(&z__3, &a[i__ + j * a_dim1]);
#line 336 "zgemv.f"
			i__3 = ix;
#line 336 "zgemv.f"
			z__2.r = z__3.r * x[i__3].r - z__3.i * x[i__3].i, 
				z__2.i = z__3.r * x[i__3].i + z__3.i * x[i__3]
				.r;
#line 336 "zgemv.f"
			z__1.r = temp.r + z__2.r, z__1.i = temp.i + z__2.i;
#line 336 "zgemv.f"
			temp.r = z__1.r, temp.i = z__1.i;
#line 337 "zgemv.f"
			ix += *incx;
#line 338 "zgemv.f"
/* L130: */
#line 338 "zgemv.f"
		    }
#line 339 "zgemv.f"
		}
#line 340 "zgemv.f"
		i__2 = jy;
#line 340 "zgemv.f"
		i__3 = jy;
#line 340 "zgemv.f"
		z__2.r = alpha->r * temp.r - alpha->i * temp.i, z__2.i = 
			alpha->r * temp.i + alpha->i * temp.r;
#line 340 "zgemv.f"
		z__1.r = y[i__3].r + z__2.r, z__1.i = y[i__3].i + z__2.i;
#line 340 "zgemv.f"
		y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 341 "zgemv.f"
		jy += *incy;
#line 342 "zgemv.f"
/* L140: */
#line 342 "zgemv.f"
	    }
#line 343 "zgemv.f"
	}
#line 344 "zgemv.f"
    }

#line 346 "zgemv.f"
    return 0;

/*     End of ZGEMV . */

} /* zgemv_ */

