#line 1 "cla_syamv.f"
/* cla_syamv.f -- translated by f2c (version 20100827).
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

#line 1 "cla_syamv.f"
/* > \brief \b CLA_SYAMV computes a matrix-vector product using a symmetric indefinite matrix to calculate err
or bounds. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLA_SYAMV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cla_sya
mv.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cla_sya
mv.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cla_sya
mv.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CLA_SYAMV( UPLO, N, ALPHA, A, LDA, X, INCX, BETA, Y, */
/*                             INCY ) */

/*       .. Scalar Arguments .. */
/*       REAL               ALPHA, BETA */
/*       INTEGER            INCX, INCY, LDA, N */
/*       INTEGER            UPLO */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX            A( LDA, * ), X( * ) */
/*       REAL               Y( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLA_SYAMV  performs the matrix-vector operation */
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
/* >          ALPHA is REAL . */
/* >           On entry, ALPHA specifies the scalar alpha. */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is COMPLEX array of DIMENSION ( LDA, n ). */
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
/* >          X is COMPLEX array, dimension */
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
/* >          BETA is REAL . */
/* >           On entry, BETA specifies the scalar beta. When BETA is */
/* >           supplied as zero then Y need not be set on input. */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Y */
/* > \verbatim */
/* >          Y is REAL array, dimension */
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

/* > \ingroup complexSYcomputational */

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
/* Subroutine */ int cla_syamv__(integer *uplo, integer *n, doublereal *alpha,
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
    extern doublereal slamch_(char *, ftnlen);
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

#line 230 "cla_syamv.f"
    /* Parameter adjustments */
#line 230 "cla_syamv.f"
    a_dim1 = *lda;
#line 230 "cla_syamv.f"
    a_offset = 1 + a_dim1;
#line 230 "cla_syamv.f"
    a -= a_offset;
#line 230 "cla_syamv.f"
    --x;
#line 230 "cla_syamv.f"
    --y;
#line 230 "cla_syamv.f"

#line 230 "cla_syamv.f"
    /* Function Body */
#line 230 "cla_syamv.f"
    info = 0;
#line 231 "cla_syamv.f"
    if (*uplo != ilauplo_("U", (ftnlen)1) && *uplo != ilauplo_("L", (ftnlen)1)
	    ) {
#line 233 "cla_syamv.f"
	info = 1;
#line 234 "cla_syamv.f"
    } else if (*n < 0) {
#line 235 "cla_syamv.f"
	info = 2;
#line 236 "cla_syamv.f"
    } else if (*lda < max(1,*n)) {
#line 237 "cla_syamv.f"
	info = 5;
#line 238 "cla_syamv.f"
    } else if (*incx == 0) {
#line 239 "cla_syamv.f"
	info = 7;
#line 240 "cla_syamv.f"
    } else if (*incy == 0) {
#line 241 "cla_syamv.f"
	info = 10;
#line 242 "cla_syamv.f"
    }
#line 243 "cla_syamv.f"
    if (info != 0) {
#line 244 "cla_syamv.f"
	xerbla_("SSYMV ", &info, (ftnlen)6);
#line 245 "cla_syamv.f"
	return 0;
#line 246 "cla_syamv.f"
    }

/*     Quick return if possible. */

#line 250 "cla_syamv.f"
    if (*n == 0 || *alpha == 0. && *beta == 1.) {
#line 250 "cla_syamv.f"
	return 0;
#line 250 "cla_syamv.f"
    }

/*     Set up the start points in  X  and  Y. */

#line 255 "cla_syamv.f"
    if (*incx > 0) {
#line 256 "cla_syamv.f"
	kx = 1;
#line 257 "cla_syamv.f"
    } else {
#line 258 "cla_syamv.f"
	kx = 1 - (*n - 1) * *incx;
#line 259 "cla_syamv.f"
    }
#line 260 "cla_syamv.f"
    if (*incy > 0) {
#line 261 "cla_syamv.f"
	ky = 1;
#line 262 "cla_syamv.f"
    } else {
#line 263 "cla_syamv.f"
	ky = 1 - (*n - 1) * *incy;
#line 264 "cla_syamv.f"
    }

/*     Set SAFE1 essentially to be the underflow threshold times the */
/*     number of additions in each row. */

#line 269 "cla_syamv.f"
    safe1 = slamch_("Safe minimum", (ftnlen)12);
#line 270 "cla_syamv.f"
    safe1 = (*n + 1) * safe1;

/*     Form  y := alpha*abs(A)*abs(x) + beta*abs(y). */

/*     The O(N^2) SYMB_ZERO tests could be replaced by O(N) queries to */
/*     the inexact flag.  Still doesn't help change the iteration order */
/*     to per-column. */

#line 278 "cla_syamv.f"
    iy = ky;
#line 279 "cla_syamv.f"
    if (*incx == 1) {
#line 280 "cla_syamv.f"
	if (*uplo == ilauplo_("U", (ftnlen)1)) {
#line 281 "cla_syamv.f"
	    i__1 = *n;
#line 281 "cla_syamv.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 282 "cla_syamv.f"
		if (*beta == 0.) {
#line 283 "cla_syamv.f"
		    symb_zero__ = TRUE_;
#line 284 "cla_syamv.f"
		    y[iy] = 0.;
#line 285 "cla_syamv.f"
		} else if (y[iy] == 0.) {
#line 286 "cla_syamv.f"
		    symb_zero__ = TRUE_;
#line 287 "cla_syamv.f"
		} else {
#line 288 "cla_syamv.f"
		    symb_zero__ = FALSE_;
#line 289 "cla_syamv.f"
		    y[iy] = *beta * (d__1 = y[iy], abs(d__1));
#line 290 "cla_syamv.f"
		}
#line 291 "cla_syamv.f"
		if (*alpha != 0.) {
#line 292 "cla_syamv.f"
		    i__2 = i__;
#line 292 "cla_syamv.f"
		    for (j = 1; j <= i__2; ++j) {
#line 293 "cla_syamv.f"
			i__3 = j + i__ * a_dim1;
#line 293 "cla_syamv.f"
			temp = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(
				&a[j + i__ * a_dim1]), abs(d__2));
#line 294 "cla_syamv.f"
			i__3 = j;
#line 294 "cla_syamv.f"
			symb_zero__ = symb_zero__ && (x[i__3].r == 0. && x[
				i__3].i == 0. || temp == 0.);
#line 297 "cla_syamv.f"
			i__3 = j;
#line 297 "cla_syamv.f"
			y[iy] += *alpha * ((d__1 = x[i__3].r, abs(d__1)) + (
				d__2 = d_imag(&x[j]), abs(d__2))) * temp;
#line 298 "cla_syamv.f"
		    }
#line 299 "cla_syamv.f"
		    i__2 = *n;
#line 299 "cla_syamv.f"
		    for (j = i__ + 1; j <= i__2; ++j) {
#line 300 "cla_syamv.f"
			i__3 = i__ + j * a_dim1;
#line 300 "cla_syamv.f"
			temp = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(
				&a[i__ + j * a_dim1]), abs(d__2));
#line 301 "cla_syamv.f"
			i__3 = j;
#line 301 "cla_syamv.f"
			symb_zero__ = symb_zero__ && (x[i__3].r == 0. && x[
				i__3].i == 0. || temp == 0.);
#line 304 "cla_syamv.f"
			i__3 = j;
#line 304 "cla_syamv.f"
			y[iy] += *alpha * ((d__1 = x[i__3].r, abs(d__1)) + (
				d__2 = d_imag(&x[j]), abs(d__2))) * temp;
#line 305 "cla_syamv.f"
		    }
#line 306 "cla_syamv.f"
		}
#line 308 "cla_syamv.f"
		if (! symb_zero__) {
#line 308 "cla_syamv.f"
		    y[iy] += d_sign(&safe1, &y[iy]);
#line 308 "cla_syamv.f"
		}
#line 311 "cla_syamv.f"
		iy += *incy;
#line 312 "cla_syamv.f"
	    }
#line 313 "cla_syamv.f"
	} else {
#line 314 "cla_syamv.f"
	    i__1 = *n;
#line 314 "cla_syamv.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 315 "cla_syamv.f"
		if (*beta == 0.) {
#line 316 "cla_syamv.f"
		    symb_zero__ = TRUE_;
#line 317 "cla_syamv.f"
		    y[iy] = 0.;
#line 318 "cla_syamv.f"
		} else if (y[iy] == 0.) {
#line 319 "cla_syamv.f"
		    symb_zero__ = TRUE_;
#line 320 "cla_syamv.f"
		} else {
#line 321 "cla_syamv.f"
		    symb_zero__ = FALSE_;
#line 322 "cla_syamv.f"
		    y[iy] = *beta * (d__1 = y[iy], abs(d__1));
#line 323 "cla_syamv.f"
		}
#line 324 "cla_syamv.f"
		if (*alpha != 0.) {
#line 325 "cla_syamv.f"
		    i__2 = i__;
#line 325 "cla_syamv.f"
		    for (j = 1; j <= i__2; ++j) {
#line 326 "cla_syamv.f"
			i__3 = i__ + j * a_dim1;
#line 326 "cla_syamv.f"
			temp = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(
				&a[i__ + j * a_dim1]), abs(d__2));
#line 327 "cla_syamv.f"
			i__3 = j;
#line 327 "cla_syamv.f"
			symb_zero__ = symb_zero__ && (x[i__3].r == 0. && x[
				i__3].i == 0. || temp == 0.);
#line 330 "cla_syamv.f"
			i__3 = j;
#line 330 "cla_syamv.f"
			y[iy] += *alpha * ((d__1 = x[i__3].r, abs(d__1)) + (
				d__2 = d_imag(&x[j]), abs(d__2))) * temp;
#line 331 "cla_syamv.f"
		    }
#line 332 "cla_syamv.f"
		    i__2 = *n;
#line 332 "cla_syamv.f"
		    for (j = i__ + 1; j <= i__2; ++j) {
#line 333 "cla_syamv.f"
			i__3 = j + i__ * a_dim1;
#line 333 "cla_syamv.f"
			temp = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(
				&a[j + i__ * a_dim1]), abs(d__2));
#line 334 "cla_syamv.f"
			i__3 = j;
#line 334 "cla_syamv.f"
			symb_zero__ = symb_zero__ && (x[i__3].r == 0. && x[
				i__3].i == 0. || temp == 0.);
#line 337 "cla_syamv.f"
			i__3 = j;
#line 337 "cla_syamv.f"
			y[iy] += *alpha * ((d__1 = x[i__3].r, abs(d__1)) + (
				d__2 = d_imag(&x[j]), abs(d__2))) * temp;
#line 338 "cla_syamv.f"
		    }
#line 339 "cla_syamv.f"
		}
#line 341 "cla_syamv.f"
		if (! symb_zero__) {
#line 341 "cla_syamv.f"
		    y[iy] += d_sign(&safe1, &y[iy]);
#line 341 "cla_syamv.f"
		}
#line 344 "cla_syamv.f"
		iy += *incy;
#line 345 "cla_syamv.f"
	    }
#line 346 "cla_syamv.f"
	}
#line 347 "cla_syamv.f"
    } else {
#line 348 "cla_syamv.f"
	if (*uplo == ilauplo_("U", (ftnlen)1)) {
#line 349 "cla_syamv.f"
	    i__1 = *n;
#line 349 "cla_syamv.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 350 "cla_syamv.f"
		if (*beta == 0.) {
#line 351 "cla_syamv.f"
		    symb_zero__ = TRUE_;
#line 352 "cla_syamv.f"
		    y[iy] = 0.;
#line 353 "cla_syamv.f"
		} else if (y[iy] == 0.) {
#line 354 "cla_syamv.f"
		    symb_zero__ = TRUE_;
#line 355 "cla_syamv.f"
		} else {
#line 356 "cla_syamv.f"
		    symb_zero__ = FALSE_;
#line 357 "cla_syamv.f"
		    y[iy] = *beta * (d__1 = y[iy], abs(d__1));
#line 358 "cla_syamv.f"
		}
#line 359 "cla_syamv.f"
		jx = kx;
#line 360 "cla_syamv.f"
		if (*alpha != 0.) {
#line 361 "cla_syamv.f"
		    i__2 = i__;
#line 361 "cla_syamv.f"
		    for (j = 1; j <= i__2; ++j) {
#line 362 "cla_syamv.f"
			i__3 = j + i__ * a_dim1;
#line 362 "cla_syamv.f"
			temp = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(
				&a[j + i__ * a_dim1]), abs(d__2));
#line 363 "cla_syamv.f"
			i__3 = j;
#line 363 "cla_syamv.f"
			symb_zero__ = symb_zero__ && (x[i__3].r == 0. && x[
				i__3].i == 0. || temp == 0.);
#line 366 "cla_syamv.f"
			i__3 = jx;
#line 366 "cla_syamv.f"
			y[iy] += *alpha * ((d__1 = x[i__3].r, abs(d__1)) + (
				d__2 = d_imag(&x[jx]), abs(d__2))) * temp;
#line 367 "cla_syamv.f"
			jx += *incx;
#line 368 "cla_syamv.f"
		    }
#line 369 "cla_syamv.f"
		    i__2 = *n;
#line 369 "cla_syamv.f"
		    for (j = i__ + 1; j <= i__2; ++j) {
#line 370 "cla_syamv.f"
			i__3 = i__ + j * a_dim1;
#line 370 "cla_syamv.f"
			temp = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(
				&a[i__ + j * a_dim1]), abs(d__2));
#line 371 "cla_syamv.f"
			i__3 = j;
#line 371 "cla_syamv.f"
			symb_zero__ = symb_zero__ && (x[i__3].r == 0. && x[
				i__3].i == 0. || temp == 0.);
#line 374 "cla_syamv.f"
			i__3 = jx;
#line 374 "cla_syamv.f"
			y[iy] += *alpha * ((d__1 = x[i__3].r, abs(d__1)) + (
				d__2 = d_imag(&x[jx]), abs(d__2))) * temp;
#line 375 "cla_syamv.f"
			jx += *incx;
#line 376 "cla_syamv.f"
		    }
#line 377 "cla_syamv.f"
		}
#line 379 "cla_syamv.f"
		if (! symb_zero__) {
#line 379 "cla_syamv.f"
		    y[iy] += d_sign(&safe1, &y[iy]);
#line 379 "cla_syamv.f"
		}
#line 382 "cla_syamv.f"
		iy += *incy;
#line 383 "cla_syamv.f"
	    }
#line 384 "cla_syamv.f"
	} else {
#line 385 "cla_syamv.f"
	    i__1 = *n;
#line 385 "cla_syamv.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 386 "cla_syamv.f"
		if (*beta == 0.) {
#line 387 "cla_syamv.f"
		    symb_zero__ = TRUE_;
#line 388 "cla_syamv.f"
		    y[iy] = 0.;
#line 389 "cla_syamv.f"
		} else if (y[iy] == 0.) {
#line 390 "cla_syamv.f"
		    symb_zero__ = TRUE_;
#line 391 "cla_syamv.f"
		} else {
#line 392 "cla_syamv.f"
		    symb_zero__ = FALSE_;
#line 393 "cla_syamv.f"
		    y[iy] = *beta * (d__1 = y[iy], abs(d__1));
#line 394 "cla_syamv.f"
		}
#line 395 "cla_syamv.f"
		jx = kx;
#line 396 "cla_syamv.f"
		if (*alpha != 0.) {
#line 397 "cla_syamv.f"
		    i__2 = i__;
#line 397 "cla_syamv.f"
		    for (j = 1; j <= i__2; ++j) {
#line 398 "cla_syamv.f"
			i__3 = i__ + j * a_dim1;
#line 398 "cla_syamv.f"
			temp = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(
				&a[i__ + j * a_dim1]), abs(d__2));
#line 399 "cla_syamv.f"
			i__3 = j;
#line 399 "cla_syamv.f"
			symb_zero__ = symb_zero__ && (x[i__3].r == 0. && x[
				i__3].i == 0. || temp == 0.);
#line 402 "cla_syamv.f"
			i__3 = jx;
#line 402 "cla_syamv.f"
			y[iy] += *alpha * ((d__1 = x[i__3].r, abs(d__1)) + (
				d__2 = d_imag(&x[jx]), abs(d__2))) * temp;
#line 403 "cla_syamv.f"
			jx += *incx;
#line 404 "cla_syamv.f"
		    }
#line 405 "cla_syamv.f"
		    i__2 = *n;
#line 405 "cla_syamv.f"
		    for (j = i__ + 1; j <= i__2; ++j) {
#line 406 "cla_syamv.f"
			i__3 = j + i__ * a_dim1;
#line 406 "cla_syamv.f"
			temp = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(
				&a[j + i__ * a_dim1]), abs(d__2));
#line 407 "cla_syamv.f"
			i__3 = j;
#line 407 "cla_syamv.f"
			symb_zero__ = symb_zero__ && (x[i__3].r == 0. && x[
				i__3].i == 0. || temp == 0.);
#line 410 "cla_syamv.f"
			i__3 = jx;
#line 410 "cla_syamv.f"
			y[iy] += *alpha * ((d__1 = x[i__3].r, abs(d__1)) + (
				d__2 = d_imag(&x[jx]), abs(d__2))) * temp;
#line 411 "cla_syamv.f"
			jx += *incx;
#line 412 "cla_syamv.f"
		    }
#line 413 "cla_syamv.f"
		}
#line 415 "cla_syamv.f"
		if (! symb_zero__) {
#line 415 "cla_syamv.f"
		    y[iy] += d_sign(&safe1, &y[iy]);
#line 415 "cla_syamv.f"
		}
#line 418 "cla_syamv.f"
		iy += *incy;
#line 419 "cla_syamv.f"
	    }
#line 420 "cla_syamv.f"
	}
#line 422 "cla_syamv.f"
    }

#line 424 "cla_syamv.f"
    return 0;

/*     End of CLA_SYAMV */

} /* cla_syamv__ */

