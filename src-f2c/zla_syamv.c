#line 1 "zla_syamv.f"
/* zla_syamv.f -- translated by f2c (version 20100827).
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

#line 1 "zla_syamv.f"
/* > \brief \b ZLA_SYAMV computes a matrix-vector product using a symmetric indefinite matrix to calculate err
or bounds. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZLA_SYAMV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zla_sya
mv.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zla_sya
mv.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zla_sya
mv.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZLA_SYAMV( UPLO, N, ALPHA, A, LDA, X, INCX, BETA, Y, */
/*                             INCY ) */

/*       .. Scalar Arguments .. */
/*       DOUBLE PRECISION   ALPHA, BETA */
/*       INTEGER            INCX, INCY, LDA, N */
/*       INTEGER            UPLO */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         A( LDA, * ), X( * ) */
/*       DOUBLE PRECISION   Y( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZLA_SYAMV  performs the matrix-vector operation */
/* > */
/* >         y := alpha*abs(A)*abs(x) + beta*abs(y), */
/* > */
/* > where alpha and beta are scalars, x and y are vectors and A is an */
/* > n by n symmetric matrix. */
/* > */
/* > This function is primarily used in calculating error bounds. */
/* > To protect against underflow during evaluation, components in */
/* > the resulting vector are perturbed away from zero by (N+1) */
/* > times the underflow threshold.  To prevent unnecessarily large */
/* > errors for block-structure embedded in general matrices, */
/* > "symbolically" zero components are not perturbed.  A zero */
/* > entry is considered "symbolic" if all multiplications involved */
/* > in computing that entry have at least one zero multiplicand. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is INTEGER */
/* >           On entry, UPLO specifies whether the upper or lower */
/* >           triangular part of the array A is to be referenced as */
/* >           follows: */
/* > */
/* >              UPLO = BLAS_UPPER   Only the upper triangular part of A */
/* >                                  is to be referenced. */
/* > */
/* >              UPLO = BLAS_LOWER   Only the lower triangular part of A */
/* >                                  is to be referenced. */
/* > */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >           On entry, N specifies the number of columns of the matrix A. */
/* >           N must be at least zero. */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] ALPHA */
/* > \verbatim */
/* >          ALPHA is DOUBLE PRECISION . */
/* >           On entry, ALPHA specifies the scalar alpha. */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array, DIMENSION ( LDA, n ). */
/* >           Before entry, the leading m by n part of the array A must */
/* >           contain the matrix of coefficients. */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >           On entry, LDA specifies the first dimension of A as declared */
/* >           in the calling (sub) program. LDA must be at least */
/* >           max( 1, n ). */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] X */
/* > \verbatim */
/* >          X is COMPLEX*16 array, DIMENSION at least */
/* >           ( 1 + ( n - 1 )*abs( INCX ) ) */
/* >           Before entry, the incremented array X must contain the */
/* >           vector x. */
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
/* >          BETA is DOUBLE PRECISION . */
/* >           On entry, BETA specifies the scalar beta. When BETA is */
/* >           supplied as zero then Y need not be set on input. */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Y */
/* > \verbatim */
/* >          Y is DOUBLE PRECISION array, dimension */
/* >           ( 1 + ( n - 1 )*abs( INCY ) ) */
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
/* >           Unchanged on exit. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date September 2012 */

/* > \ingroup complex16SYcomputational */

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
/* >  -- Modified for the absolute-value product, April 2006 */
/* >     Jason Riedy, UC Berkeley */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int zla_syamv__(integer *uplo, integer *n, doublereal *alpha,
	 doublecomplex *a, integer *lda, doublecomplex *x, integer *incx, 
	doublereal *beta, doublereal *y, integer *incy)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Builtin functions */
    double d_imag(doublecomplex *), d_sign(doublereal *, doublereal *);

    /* Local variables */
    static integer i__, j;
    static logical symb_zero__;
    static integer iy, jx, kx, ky, info;
    static doublereal temp, safe1;
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern integer ilauplo_(char *, ftnlen);


/*  -- LAPACK computational routine (version 3.4.2) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     September 2012 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Statement Functions .. */
/*     .. */
/*     .. Statement Function Definitions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 230 "zla_syamv.f"
    /* Parameter adjustments */
#line 230 "zla_syamv.f"
    a_dim1 = *lda;
#line 230 "zla_syamv.f"
    a_offset = 1 + a_dim1;
#line 230 "zla_syamv.f"
    a -= a_offset;
#line 230 "zla_syamv.f"
    --x;
#line 230 "zla_syamv.f"
    --y;
#line 230 "zla_syamv.f"

#line 230 "zla_syamv.f"
    /* Function Body */
#line 230 "zla_syamv.f"
    info = 0;
#line 231 "zla_syamv.f"
    if (*uplo != ilauplo_("U", (ftnlen)1) && *uplo != ilauplo_("L", (ftnlen)1)
	    ) {
#line 233 "zla_syamv.f"
	info = 1;
#line 234 "zla_syamv.f"
    } else if (*n < 0) {
#line 235 "zla_syamv.f"
	info = 2;
#line 236 "zla_syamv.f"
    } else if (*lda < max(1,*n)) {
#line 237 "zla_syamv.f"
	info = 5;
#line 238 "zla_syamv.f"
    } else if (*incx == 0) {
#line 239 "zla_syamv.f"
	info = 7;
#line 240 "zla_syamv.f"
    } else if (*incy == 0) {
#line 241 "zla_syamv.f"
	info = 10;
#line 242 "zla_syamv.f"
    }
#line 243 "zla_syamv.f"
    if (info != 0) {
#line 244 "zla_syamv.f"
	xerbla_("DSYMV ", &info, (ftnlen)6);
#line 245 "zla_syamv.f"
	return 0;
#line 246 "zla_syamv.f"
    }

/*     Quick return if possible. */

#line 250 "zla_syamv.f"
    if (*n == 0 || *alpha == 0. && *beta == 1.) {
#line 250 "zla_syamv.f"
	return 0;
#line 250 "zla_syamv.f"
    }

/*     Set up the start points in  X  and  Y. */

#line 255 "zla_syamv.f"
    if (*incx > 0) {
#line 256 "zla_syamv.f"
	kx = 1;
#line 257 "zla_syamv.f"
    } else {
#line 258 "zla_syamv.f"
	kx = 1 - (*n - 1) * *incx;
#line 259 "zla_syamv.f"
    }
#line 260 "zla_syamv.f"
    if (*incy > 0) {
#line 261 "zla_syamv.f"
	ky = 1;
#line 262 "zla_syamv.f"
    } else {
#line 263 "zla_syamv.f"
	ky = 1 - (*n - 1) * *incy;
#line 264 "zla_syamv.f"
    }

/*     Set SAFE1 essentially to be the underflow threshold times the */
/*     number of additions in each row. */

#line 269 "zla_syamv.f"
    safe1 = dlamch_("Safe minimum", (ftnlen)12);
#line 270 "zla_syamv.f"
    safe1 = (*n + 1) * safe1;

/*     Form  y := alpha*abs(A)*abs(x) + beta*abs(y). */

/*     The O(N^2) SYMB_ZERO tests could be replaced by O(N) queries to */
/*     the inexact flag.  Still doesn't help change the iteration order */
/*     to per-column. */

#line 278 "zla_syamv.f"
    iy = ky;
#line 279 "zla_syamv.f"
    if (*incx == 1) {
#line 280 "zla_syamv.f"
	if (*uplo == ilauplo_("U", (ftnlen)1)) {
#line 281 "zla_syamv.f"
	    i__1 = *n;
#line 281 "zla_syamv.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 282 "zla_syamv.f"
		if (*beta == 0.) {
#line 283 "zla_syamv.f"
		    symb_zero__ = TRUE_;
#line 284 "zla_syamv.f"
		    y[iy] = 0.;
#line 285 "zla_syamv.f"
		} else if (y[iy] == 0.) {
#line 286 "zla_syamv.f"
		    symb_zero__ = TRUE_;
#line 287 "zla_syamv.f"
		} else {
#line 288 "zla_syamv.f"
		    symb_zero__ = FALSE_;
#line 289 "zla_syamv.f"
		    y[iy] = *beta * (d__1 = y[iy], abs(d__1));
#line 290 "zla_syamv.f"
		}
#line 291 "zla_syamv.f"
		if (*alpha != 0.) {
#line 292 "zla_syamv.f"
		    i__2 = i__;
#line 292 "zla_syamv.f"
		    for (j = 1; j <= i__2; ++j) {
#line 293 "zla_syamv.f"
			i__3 = j + i__ * a_dim1;
#line 293 "zla_syamv.f"
			temp = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(
				&a[j + i__ * a_dim1]), abs(d__2));
#line 294 "zla_syamv.f"
			i__3 = j;
#line 294 "zla_syamv.f"
			symb_zero__ = symb_zero__ && (x[i__3].r == 0. && x[
				i__3].i == 0. || temp == 0.);
#line 297 "zla_syamv.f"
			i__3 = j;
#line 297 "zla_syamv.f"
			y[iy] += *alpha * ((d__1 = x[i__3].r, abs(d__1)) + (
				d__2 = d_imag(&x[j]), abs(d__2))) * temp;
#line 298 "zla_syamv.f"
		    }
#line 299 "zla_syamv.f"
		    i__2 = *n;
#line 299 "zla_syamv.f"
		    for (j = i__ + 1; j <= i__2; ++j) {
#line 300 "zla_syamv.f"
			i__3 = i__ + j * a_dim1;
#line 300 "zla_syamv.f"
			temp = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(
				&a[i__ + j * a_dim1]), abs(d__2));
#line 301 "zla_syamv.f"
			i__3 = j;
#line 301 "zla_syamv.f"
			symb_zero__ = symb_zero__ && (x[i__3].r == 0. && x[
				i__3].i == 0. || temp == 0.);
#line 304 "zla_syamv.f"
			i__3 = j;
#line 304 "zla_syamv.f"
			y[iy] += *alpha * ((d__1 = x[i__3].r, abs(d__1)) + (
				d__2 = d_imag(&x[j]), abs(d__2))) * temp;
#line 305 "zla_syamv.f"
		    }
#line 306 "zla_syamv.f"
		}
#line 308 "zla_syamv.f"
		if (! symb_zero__) {
#line 308 "zla_syamv.f"
		    y[iy] += d_sign(&safe1, &y[iy]);
#line 308 "zla_syamv.f"
		}
#line 311 "zla_syamv.f"
		iy += *incy;
#line 312 "zla_syamv.f"
	    }
#line 313 "zla_syamv.f"
	} else {
#line 314 "zla_syamv.f"
	    i__1 = *n;
#line 314 "zla_syamv.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 315 "zla_syamv.f"
		if (*beta == 0.) {
#line 316 "zla_syamv.f"
		    symb_zero__ = TRUE_;
#line 317 "zla_syamv.f"
		    y[iy] = 0.;
#line 318 "zla_syamv.f"
		} else if (y[iy] == 0.) {
#line 319 "zla_syamv.f"
		    symb_zero__ = TRUE_;
#line 320 "zla_syamv.f"
		} else {
#line 321 "zla_syamv.f"
		    symb_zero__ = FALSE_;
#line 322 "zla_syamv.f"
		    y[iy] = *beta * (d__1 = y[iy], abs(d__1));
#line 323 "zla_syamv.f"
		}
#line 324 "zla_syamv.f"
		if (*alpha != 0.) {
#line 325 "zla_syamv.f"
		    i__2 = i__;
#line 325 "zla_syamv.f"
		    for (j = 1; j <= i__2; ++j) {
#line 326 "zla_syamv.f"
			i__3 = i__ + j * a_dim1;
#line 326 "zla_syamv.f"
			temp = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(
				&a[i__ + j * a_dim1]), abs(d__2));
#line 327 "zla_syamv.f"
			i__3 = j;
#line 327 "zla_syamv.f"
			symb_zero__ = symb_zero__ && (x[i__3].r == 0. && x[
				i__3].i == 0. || temp == 0.);
#line 330 "zla_syamv.f"
			i__3 = j;
#line 330 "zla_syamv.f"
			y[iy] += *alpha * ((d__1 = x[i__3].r, abs(d__1)) + (
				d__2 = d_imag(&x[j]), abs(d__2))) * temp;
#line 331 "zla_syamv.f"
		    }
#line 332 "zla_syamv.f"
		    i__2 = *n;
#line 332 "zla_syamv.f"
		    for (j = i__ + 1; j <= i__2; ++j) {
#line 333 "zla_syamv.f"
			i__3 = j + i__ * a_dim1;
#line 333 "zla_syamv.f"
			temp = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(
				&a[j + i__ * a_dim1]), abs(d__2));
#line 334 "zla_syamv.f"
			i__3 = j;
#line 334 "zla_syamv.f"
			symb_zero__ = symb_zero__ && (x[i__3].r == 0. && x[
				i__3].i == 0. || temp == 0.);
#line 337 "zla_syamv.f"
			i__3 = j;
#line 337 "zla_syamv.f"
			y[iy] += *alpha * ((d__1 = x[i__3].r, abs(d__1)) + (
				d__2 = d_imag(&x[j]), abs(d__2))) * temp;
#line 338 "zla_syamv.f"
		    }
#line 339 "zla_syamv.f"
		}
#line 341 "zla_syamv.f"
		if (! symb_zero__) {
#line 341 "zla_syamv.f"
		    y[iy] += d_sign(&safe1, &y[iy]);
#line 341 "zla_syamv.f"
		}
#line 344 "zla_syamv.f"
		iy += *incy;
#line 345 "zla_syamv.f"
	    }
#line 346 "zla_syamv.f"
	}
#line 347 "zla_syamv.f"
    } else {
#line 348 "zla_syamv.f"
	if (*uplo == ilauplo_("U", (ftnlen)1)) {
#line 349 "zla_syamv.f"
	    i__1 = *n;
#line 349 "zla_syamv.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 350 "zla_syamv.f"
		if (*beta == 0.) {
#line 351 "zla_syamv.f"
		    symb_zero__ = TRUE_;
#line 352 "zla_syamv.f"
		    y[iy] = 0.;
#line 353 "zla_syamv.f"
		} else if (y[iy] == 0.) {
#line 354 "zla_syamv.f"
		    symb_zero__ = TRUE_;
#line 355 "zla_syamv.f"
		} else {
#line 356 "zla_syamv.f"
		    symb_zero__ = FALSE_;
#line 357 "zla_syamv.f"
		    y[iy] = *beta * (d__1 = y[iy], abs(d__1));
#line 358 "zla_syamv.f"
		}
#line 359 "zla_syamv.f"
		jx = kx;
#line 360 "zla_syamv.f"
		if (*alpha != 0.) {
#line 361 "zla_syamv.f"
		    i__2 = i__;
#line 361 "zla_syamv.f"
		    for (j = 1; j <= i__2; ++j) {
#line 362 "zla_syamv.f"
			i__3 = j + i__ * a_dim1;
#line 362 "zla_syamv.f"
			temp = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(
				&a[j + i__ * a_dim1]), abs(d__2));
#line 363 "zla_syamv.f"
			i__3 = j;
#line 363 "zla_syamv.f"
			symb_zero__ = symb_zero__ && (x[i__3].r == 0. && x[
				i__3].i == 0. || temp == 0.);
#line 366 "zla_syamv.f"
			i__3 = jx;
#line 366 "zla_syamv.f"
			y[iy] += *alpha * ((d__1 = x[i__3].r, abs(d__1)) + (
				d__2 = d_imag(&x[jx]), abs(d__2))) * temp;
#line 367 "zla_syamv.f"
			jx += *incx;
#line 368 "zla_syamv.f"
		    }
#line 369 "zla_syamv.f"
		    i__2 = *n;
#line 369 "zla_syamv.f"
		    for (j = i__ + 1; j <= i__2; ++j) {
#line 370 "zla_syamv.f"
			i__3 = i__ + j * a_dim1;
#line 370 "zla_syamv.f"
			temp = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(
				&a[i__ + j * a_dim1]), abs(d__2));
#line 371 "zla_syamv.f"
			i__3 = j;
#line 371 "zla_syamv.f"
			symb_zero__ = symb_zero__ && (x[i__3].r == 0. && x[
				i__3].i == 0. || temp == 0.);
#line 374 "zla_syamv.f"
			i__3 = jx;
#line 374 "zla_syamv.f"
			y[iy] += *alpha * ((d__1 = x[i__3].r, abs(d__1)) + (
				d__2 = d_imag(&x[jx]), abs(d__2))) * temp;
#line 375 "zla_syamv.f"
			jx += *incx;
#line 376 "zla_syamv.f"
		    }
#line 377 "zla_syamv.f"
		}
#line 379 "zla_syamv.f"
		if (! symb_zero__) {
#line 379 "zla_syamv.f"
		    y[iy] += d_sign(&safe1, &y[iy]);
#line 379 "zla_syamv.f"
		}
#line 382 "zla_syamv.f"
		iy += *incy;
#line 383 "zla_syamv.f"
	    }
#line 384 "zla_syamv.f"
	} else {
#line 385 "zla_syamv.f"
	    i__1 = *n;
#line 385 "zla_syamv.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 386 "zla_syamv.f"
		if (*beta == 0.) {
#line 387 "zla_syamv.f"
		    symb_zero__ = TRUE_;
#line 388 "zla_syamv.f"
		    y[iy] = 0.;
#line 389 "zla_syamv.f"
		} else if (y[iy] == 0.) {
#line 390 "zla_syamv.f"
		    symb_zero__ = TRUE_;
#line 391 "zla_syamv.f"
		} else {
#line 392 "zla_syamv.f"
		    symb_zero__ = FALSE_;
#line 393 "zla_syamv.f"
		    y[iy] = *beta * (d__1 = y[iy], abs(d__1));
#line 394 "zla_syamv.f"
		}
#line 395 "zla_syamv.f"
		jx = kx;
#line 396 "zla_syamv.f"
		if (*alpha != 0.) {
#line 397 "zla_syamv.f"
		    i__2 = i__;
#line 397 "zla_syamv.f"
		    for (j = 1; j <= i__2; ++j) {
#line 398 "zla_syamv.f"
			i__3 = i__ + j * a_dim1;
#line 398 "zla_syamv.f"
			temp = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(
				&a[i__ + j * a_dim1]), abs(d__2));
#line 399 "zla_syamv.f"
			i__3 = j;
#line 399 "zla_syamv.f"
			symb_zero__ = symb_zero__ && (x[i__3].r == 0. && x[
				i__3].i == 0. || temp == 0.);
#line 402 "zla_syamv.f"
			i__3 = jx;
#line 402 "zla_syamv.f"
			y[iy] += *alpha * ((d__1 = x[i__3].r, abs(d__1)) + (
				d__2 = d_imag(&x[jx]), abs(d__2))) * temp;
#line 403 "zla_syamv.f"
			jx += *incx;
#line 404 "zla_syamv.f"
		    }
#line 405 "zla_syamv.f"
		    i__2 = *n;
#line 405 "zla_syamv.f"
		    for (j = i__ + 1; j <= i__2; ++j) {
#line 406 "zla_syamv.f"
			i__3 = j + i__ * a_dim1;
#line 406 "zla_syamv.f"
			temp = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(
				&a[j + i__ * a_dim1]), abs(d__2));
#line 407 "zla_syamv.f"
			i__3 = j;
#line 407 "zla_syamv.f"
			symb_zero__ = symb_zero__ && (x[i__3].r == 0. && x[
				i__3].i == 0. || temp == 0.);
#line 410 "zla_syamv.f"
			i__3 = jx;
#line 410 "zla_syamv.f"
			y[iy] += *alpha * ((d__1 = x[i__3].r, abs(d__1)) + (
				d__2 = d_imag(&x[jx]), abs(d__2))) * temp;
#line 411 "zla_syamv.f"
			jx += *incx;
#line 412 "zla_syamv.f"
		    }
#line 413 "zla_syamv.f"
		}
#line 415 "zla_syamv.f"
		if (! symb_zero__) {
#line 415 "zla_syamv.f"
		    y[iy] += d_sign(&safe1, &y[iy]);
#line 415 "zla_syamv.f"
		}
#line 418 "zla_syamv.f"
		iy += *incy;
#line 419 "zla_syamv.f"
	    }
#line 420 "zla_syamv.f"
	}
#line 422 "zla_syamv.f"
    }

#line 424 "zla_syamv.f"
    return 0;

/*     End of ZLA_SYAMV */

} /* zla_syamv__ */

