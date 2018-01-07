#line 1 "cgemv.f"
/* cgemv.f -- translated by f2c (version 20100827).
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

#line 1 "cgemv.f"
/* > \brief \b CGEMV */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY) */

/*       .. Scalar Arguments .. */
/*       COMPLEX ALPHA,BETA */
/*       INTEGER INCX,INCY,LDA,M,N */
/*       CHARACTER TRANS */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX A(LDA,*),X(*),Y(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CGEMV performs one of the matrix-vector operations */
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
/* >          ALPHA is COMPLEX */
/* >           On entry, ALPHA specifies the scalar alpha. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is COMPLEX array of DIMENSION ( LDA, n ). */
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
/* >          X is COMPLEX array of DIMENSION at least */
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
/* >          BETA is COMPLEX */
/* >           On entry, BETA specifies the scalar beta. When BETA is */
/* >           supplied as zero then Y need not be set on input. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Y */
/* > \verbatim */
/* >          Y is COMPLEX array of DIMENSION at least */
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

/* > \date November 2011 */

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
/* Subroutine */ int cgemv_(char *trans, integer *m, integer *n, 
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

#line 201 "cgemv.f"
    /* Parameter adjustments */
#line 201 "cgemv.f"
    a_dim1 = *lda;
#line 201 "cgemv.f"
    a_offset = 1 + a_dim1;
#line 201 "cgemv.f"
    a -= a_offset;
#line 201 "cgemv.f"
    --x;
#line 201 "cgemv.f"
    --y;
#line 201 "cgemv.f"

#line 201 "cgemv.f"
    /* Function Body */
#line 201 "cgemv.f"
    info = 0;
#line 202 "cgemv.f"
    if (! lsame_(trans, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, "T", (
	    ftnlen)1, (ftnlen)1) && ! lsame_(trans, "C", (ftnlen)1, (ftnlen)1)
	    ) {
#line 204 "cgemv.f"
	info = 1;
#line 205 "cgemv.f"
    } else if (*m < 0) {
#line 206 "cgemv.f"
	info = 2;
#line 207 "cgemv.f"
    } else if (*n < 0) {
#line 208 "cgemv.f"
	info = 3;
#line 209 "cgemv.f"
    } else if (*lda < max(1,*m)) {
#line 210 "cgemv.f"
	info = 6;
#line 211 "cgemv.f"
    } else if (*incx == 0) {
#line 212 "cgemv.f"
	info = 8;
#line 213 "cgemv.f"
    } else if (*incy == 0) {
#line 214 "cgemv.f"
	info = 11;
#line 215 "cgemv.f"
    }
#line 216 "cgemv.f"
    if (info != 0) {
#line 217 "cgemv.f"
	xerbla_("CGEMV ", &info, (ftnlen)6);
#line 218 "cgemv.f"
	return 0;
#line 219 "cgemv.f"
    }

/*     Quick return if possible. */

#line 223 "cgemv.f"
    if (*m == 0 || *n == 0 || alpha->r == 0. && alpha->i == 0. && (beta->r == 
	    1. && beta->i == 0.)) {
#line 223 "cgemv.f"
	return 0;
#line 223 "cgemv.f"
    }

#line 226 "cgemv.f"
    noconj = lsame_(trans, "T", (ftnlen)1, (ftnlen)1);

/*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set */
/*     up the start points in  X  and  Y. */

#line 231 "cgemv.f"
    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {
#line 232 "cgemv.f"
	lenx = *n;
#line 233 "cgemv.f"
	leny = *m;
#line 234 "cgemv.f"
    } else {
#line 235 "cgemv.f"
	lenx = *m;
#line 236 "cgemv.f"
	leny = *n;
#line 237 "cgemv.f"
    }
#line 238 "cgemv.f"
    if (*incx > 0) {
#line 239 "cgemv.f"
	kx = 1;
#line 240 "cgemv.f"
    } else {
#line 241 "cgemv.f"
	kx = 1 - (lenx - 1) * *incx;
#line 242 "cgemv.f"
    }
#line 243 "cgemv.f"
    if (*incy > 0) {
#line 244 "cgemv.f"
	ky = 1;
#line 245 "cgemv.f"
    } else {
#line 246 "cgemv.f"
	ky = 1 - (leny - 1) * *incy;
#line 247 "cgemv.f"
    }

/*     Start the operations. In this version the elements of A are */
/*     accessed sequentially with one pass through A. */

/*     First form  y := beta*y. */

#line 254 "cgemv.f"
    if (beta->r != 1. || beta->i != 0.) {
#line 255 "cgemv.f"
	if (*incy == 1) {
#line 256 "cgemv.f"
	    if (beta->r == 0. && beta->i == 0.) {
#line 257 "cgemv.f"
		i__1 = leny;
#line 257 "cgemv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 258 "cgemv.f"
		    i__2 = i__;
#line 258 "cgemv.f"
		    y[i__2].r = 0., y[i__2].i = 0.;
#line 259 "cgemv.f"
/* L10: */
#line 259 "cgemv.f"
		}
#line 260 "cgemv.f"
	    } else {
#line 261 "cgemv.f"
		i__1 = leny;
#line 261 "cgemv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 262 "cgemv.f"
		    i__2 = i__;
#line 262 "cgemv.f"
		    i__3 = i__;
#line 262 "cgemv.f"
		    z__1.r = beta->r * y[i__3].r - beta->i * y[i__3].i, 
			    z__1.i = beta->r * y[i__3].i + beta->i * y[i__3]
			    .r;
#line 262 "cgemv.f"
		    y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 263 "cgemv.f"
/* L20: */
#line 263 "cgemv.f"
		}
#line 264 "cgemv.f"
	    }
#line 265 "cgemv.f"
	} else {
#line 266 "cgemv.f"
	    iy = ky;
#line 267 "cgemv.f"
	    if (beta->r == 0. && beta->i == 0.) {
#line 268 "cgemv.f"
		i__1 = leny;
#line 268 "cgemv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 269 "cgemv.f"
		    i__2 = iy;
#line 269 "cgemv.f"
		    y[i__2].r = 0., y[i__2].i = 0.;
#line 270 "cgemv.f"
		    iy += *incy;
#line 271 "cgemv.f"
/* L30: */
#line 271 "cgemv.f"
		}
#line 272 "cgemv.f"
	    } else {
#line 273 "cgemv.f"
		i__1 = leny;
#line 273 "cgemv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 274 "cgemv.f"
		    i__2 = iy;
#line 274 "cgemv.f"
		    i__3 = iy;
#line 274 "cgemv.f"
		    z__1.r = beta->r * y[i__3].r - beta->i * y[i__3].i, 
			    z__1.i = beta->r * y[i__3].i + beta->i * y[i__3]
			    .r;
#line 274 "cgemv.f"
		    y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 275 "cgemv.f"
		    iy += *incy;
#line 276 "cgemv.f"
/* L40: */
#line 276 "cgemv.f"
		}
#line 277 "cgemv.f"
	    }
#line 278 "cgemv.f"
	}
#line 279 "cgemv.f"
    }
#line 280 "cgemv.f"
    if (alpha->r == 0. && alpha->i == 0.) {
#line 280 "cgemv.f"
	return 0;
#line 280 "cgemv.f"
    }
#line 281 "cgemv.f"
    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {

/*        Form  y := alpha*A*x + y. */

#line 285 "cgemv.f"
	jx = kx;
#line 286 "cgemv.f"
	if (*incy == 1) {
#line 287 "cgemv.f"
	    i__1 = *n;
#line 287 "cgemv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 288 "cgemv.f"
		i__2 = jx;
#line 288 "cgemv.f"
		if (x[i__2].r != 0. || x[i__2].i != 0.) {
#line 289 "cgemv.f"
		    i__2 = jx;
#line 289 "cgemv.f"
		    z__1.r = alpha->r * x[i__2].r - alpha->i * x[i__2].i, 
			    z__1.i = alpha->r * x[i__2].i + alpha->i * x[i__2]
			    .r;
#line 289 "cgemv.f"
		    temp.r = z__1.r, temp.i = z__1.i;
#line 290 "cgemv.f"
		    i__2 = *m;
#line 290 "cgemv.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 291 "cgemv.f"
			i__3 = i__;
#line 291 "cgemv.f"
			i__4 = i__;
#line 291 "cgemv.f"
			i__5 = i__ + j * a_dim1;
#line 291 "cgemv.f"
			z__2.r = temp.r * a[i__5].r - temp.i * a[i__5].i, 
				z__2.i = temp.r * a[i__5].i + temp.i * a[i__5]
				.r;
#line 291 "cgemv.f"
			z__1.r = y[i__4].r + z__2.r, z__1.i = y[i__4].i + 
				z__2.i;
#line 291 "cgemv.f"
			y[i__3].r = z__1.r, y[i__3].i = z__1.i;
#line 292 "cgemv.f"
/* L50: */
#line 292 "cgemv.f"
		    }
#line 293 "cgemv.f"
		}
#line 294 "cgemv.f"
		jx += *incx;
#line 295 "cgemv.f"
/* L60: */
#line 295 "cgemv.f"
	    }
#line 296 "cgemv.f"
	} else {
#line 297 "cgemv.f"
	    i__1 = *n;
#line 297 "cgemv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 298 "cgemv.f"
		i__2 = jx;
#line 298 "cgemv.f"
		if (x[i__2].r != 0. || x[i__2].i != 0.) {
#line 299 "cgemv.f"
		    i__2 = jx;
#line 299 "cgemv.f"
		    z__1.r = alpha->r * x[i__2].r - alpha->i * x[i__2].i, 
			    z__1.i = alpha->r * x[i__2].i + alpha->i * x[i__2]
			    .r;
#line 299 "cgemv.f"
		    temp.r = z__1.r, temp.i = z__1.i;
#line 300 "cgemv.f"
		    iy = ky;
#line 301 "cgemv.f"
		    i__2 = *m;
#line 301 "cgemv.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 302 "cgemv.f"
			i__3 = iy;
#line 302 "cgemv.f"
			i__4 = iy;
#line 302 "cgemv.f"
			i__5 = i__ + j * a_dim1;
#line 302 "cgemv.f"
			z__2.r = temp.r * a[i__5].r - temp.i * a[i__5].i, 
				z__2.i = temp.r * a[i__5].i + temp.i * a[i__5]
				.r;
#line 302 "cgemv.f"
			z__1.r = y[i__4].r + z__2.r, z__1.i = y[i__4].i + 
				z__2.i;
#line 302 "cgemv.f"
			y[i__3].r = z__1.r, y[i__3].i = z__1.i;
#line 303 "cgemv.f"
			iy += *incy;
#line 304 "cgemv.f"
/* L70: */
#line 304 "cgemv.f"
		    }
#line 305 "cgemv.f"
		}
#line 306 "cgemv.f"
		jx += *incx;
#line 307 "cgemv.f"
/* L80: */
#line 307 "cgemv.f"
	    }
#line 308 "cgemv.f"
	}
#line 309 "cgemv.f"
    } else {

/*        Form  y := alpha*A**T*x + y  or  y := alpha*A**H*x + y. */

#line 313 "cgemv.f"
	jy = ky;
#line 314 "cgemv.f"
	if (*incx == 1) {
#line 315 "cgemv.f"
	    i__1 = *n;
#line 315 "cgemv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 316 "cgemv.f"
		temp.r = 0., temp.i = 0.;
#line 317 "cgemv.f"
		if (noconj) {
#line 318 "cgemv.f"
		    i__2 = *m;
#line 318 "cgemv.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 319 "cgemv.f"
			i__3 = i__ + j * a_dim1;
#line 319 "cgemv.f"
			i__4 = i__;
#line 319 "cgemv.f"
			z__2.r = a[i__3].r * x[i__4].r - a[i__3].i * x[i__4]
				.i, z__2.i = a[i__3].r * x[i__4].i + a[i__3]
				.i * x[i__4].r;
#line 319 "cgemv.f"
			z__1.r = temp.r + z__2.r, z__1.i = temp.i + z__2.i;
#line 319 "cgemv.f"
			temp.r = z__1.r, temp.i = z__1.i;
#line 320 "cgemv.f"
/* L90: */
#line 320 "cgemv.f"
		    }
#line 321 "cgemv.f"
		} else {
#line 322 "cgemv.f"
		    i__2 = *m;
#line 322 "cgemv.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 323 "cgemv.f"
			d_cnjg(&z__3, &a[i__ + j * a_dim1]);
#line 323 "cgemv.f"
			i__3 = i__;
#line 323 "cgemv.f"
			z__2.r = z__3.r * x[i__3].r - z__3.i * x[i__3].i, 
				z__2.i = z__3.r * x[i__3].i + z__3.i * x[i__3]
				.r;
#line 323 "cgemv.f"
			z__1.r = temp.r + z__2.r, z__1.i = temp.i + z__2.i;
#line 323 "cgemv.f"
			temp.r = z__1.r, temp.i = z__1.i;
#line 324 "cgemv.f"
/* L100: */
#line 324 "cgemv.f"
		    }
#line 325 "cgemv.f"
		}
#line 326 "cgemv.f"
		i__2 = jy;
#line 326 "cgemv.f"
		i__3 = jy;
#line 326 "cgemv.f"
		z__2.r = alpha->r * temp.r - alpha->i * temp.i, z__2.i = 
			alpha->r * temp.i + alpha->i * temp.r;
#line 326 "cgemv.f"
		z__1.r = y[i__3].r + z__2.r, z__1.i = y[i__3].i + z__2.i;
#line 326 "cgemv.f"
		y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 327 "cgemv.f"
		jy += *incy;
#line 328 "cgemv.f"
/* L110: */
#line 328 "cgemv.f"
	    }
#line 329 "cgemv.f"
	} else {
#line 330 "cgemv.f"
	    i__1 = *n;
#line 330 "cgemv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 331 "cgemv.f"
		temp.r = 0., temp.i = 0.;
#line 332 "cgemv.f"
		ix = kx;
#line 333 "cgemv.f"
		if (noconj) {
#line 334 "cgemv.f"
		    i__2 = *m;
#line 334 "cgemv.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 335 "cgemv.f"
			i__3 = i__ + j * a_dim1;
#line 335 "cgemv.f"
			i__4 = ix;
#line 335 "cgemv.f"
			z__2.r = a[i__3].r * x[i__4].r - a[i__3].i * x[i__4]
				.i, z__2.i = a[i__3].r * x[i__4].i + a[i__3]
				.i * x[i__4].r;
#line 335 "cgemv.f"
			z__1.r = temp.r + z__2.r, z__1.i = temp.i + z__2.i;
#line 335 "cgemv.f"
			temp.r = z__1.r, temp.i = z__1.i;
#line 336 "cgemv.f"
			ix += *incx;
#line 337 "cgemv.f"
/* L120: */
#line 337 "cgemv.f"
		    }
#line 338 "cgemv.f"
		} else {
#line 339 "cgemv.f"
		    i__2 = *m;
#line 339 "cgemv.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 340 "cgemv.f"
			d_cnjg(&z__3, &a[i__ + j * a_dim1]);
#line 340 "cgemv.f"
			i__3 = ix;
#line 340 "cgemv.f"
			z__2.r = z__3.r * x[i__3].r - z__3.i * x[i__3].i, 
				z__2.i = z__3.r * x[i__3].i + z__3.i * x[i__3]
				.r;
#line 340 "cgemv.f"
			z__1.r = temp.r + z__2.r, z__1.i = temp.i + z__2.i;
#line 340 "cgemv.f"
			temp.r = z__1.r, temp.i = z__1.i;
#line 341 "cgemv.f"
			ix += *incx;
#line 342 "cgemv.f"
/* L130: */
#line 342 "cgemv.f"
		    }
#line 343 "cgemv.f"
		}
#line 344 "cgemv.f"
		i__2 = jy;
#line 344 "cgemv.f"
		i__3 = jy;
#line 344 "cgemv.f"
		z__2.r = alpha->r * temp.r - alpha->i * temp.i, z__2.i = 
			alpha->r * temp.i + alpha->i * temp.r;
#line 344 "cgemv.f"
		z__1.r = y[i__3].r + z__2.r, z__1.i = y[i__3].i + z__2.i;
#line 344 "cgemv.f"
		y[i__2].r = z__1.r, y[i__2].i = z__1.i;
#line 345 "cgemv.f"
		jy += *incy;
#line 346 "cgemv.f"
/* L140: */
#line 346 "cgemv.f"
	    }
#line 347 "cgemv.f"
	}
#line 348 "cgemv.f"
    }

#line 350 "cgemv.f"
    return 0;

/*     End of CGEMV . */

} /* cgemv_ */

