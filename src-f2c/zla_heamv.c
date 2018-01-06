#line 1 "zla_heamv.f"
/* zla_heamv.f -- translated by f2c (version 20100827).
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

#line 1 "zla_heamv.f"
/* > \brief \b ZLA_HEAMV computes a matrix-vector product using a Hermitian indefinite matrix to calculate err
or bounds. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZLA_HEAMV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zla_hea
mv.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zla_hea
mv.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zla_hea
mv.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZLA_HEAMV( UPLO, N, ALPHA, A, LDA, X, INCX, BETA, Y, */
/*                             INCY ) */

/*       .. Scalar Arguments .. */
/*       DOUBLE PRECISION   ALPHA, BETA */
/*       INTEGER            INCX, INCY, LDA, N, UPLO */
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

/* > \date December 2016 */

/* > \ingroup complex16HEcomputational */

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
/* Subroutine */ int zla_heamv__(integer *uplo, integer *n, doublereal *alpha,
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


/*  -- LAPACK computational routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
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

#line 228 "zla_heamv.f"
    /* Parameter adjustments */
#line 228 "zla_heamv.f"
    a_dim1 = *lda;
#line 228 "zla_heamv.f"
    a_offset = 1 + a_dim1;
#line 228 "zla_heamv.f"
    a -= a_offset;
#line 228 "zla_heamv.f"
    --x;
#line 228 "zla_heamv.f"
    --y;
#line 228 "zla_heamv.f"

#line 228 "zla_heamv.f"
    /* Function Body */
#line 228 "zla_heamv.f"
    info = 0;
#line 229 "zla_heamv.f"
    if (*uplo != ilauplo_("U", (ftnlen)1) && *uplo != ilauplo_("L", (ftnlen)1)
	    ) {
#line 231 "zla_heamv.f"
	info = 1;
#line 232 "zla_heamv.f"
    } else if (*n < 0) {
#line 233 "zla_heamv.f"
	info = 2;
#line 234 "zla_heamv.f"
    } else if (*lda < max(1,*n)) {
#line 235 "zla_heamv.f"
	info = 5;
#line 236 "zla_heamv.f"
    } else if (*incx == 0) {
#line 237 "zla_heamv.f"
	info = 7;
#line 238 "zla_heamv.f"
    } else if (*incy == 0) {
#line 239 "zla_heamv.f"
	info = 10;
#line 240 "zla_heamv.f"
    }
#line 241 "zla_heamv.f"
    if (info != 0) {
#line 242 "zla_heamv.f"
	xerbla_("ZHEMV ", &info, (ftnlen)6);
#line 243 "zla_heamv.f"
	return 0;
#line 244 "zla_heamv.f"
    }

/*     Quick return if possible. */

#line 248 "zla_heamv.f"
    if (*n == 0 || *alpha == 0. && *beta == 1.) {
#line 248 "zla_heamv.f"
	return 0;
#line 248 "zla_heamv.f"
    }

/*     Set up the start points in  X  and  Y. */

#line 253 "zla_heamv.f"
    if (*incx > 0) {
#line 254 "zla_heamv.f"
	kx = 1;
#line 255 "zla_heamv.f"
    } else {
#line 256 "zla_heamv.f"
	kx = 1 - (*n - 1) * *incx;
#line 257 "zla_heamv.f"
    }
#line 258 "zla_heamv.f"
    if (*incy > 0) {
#line 259 "zla_heamv.f"
	ky = 1;
#line 260 "zla_heamv.f"
    } else {
#line 261 "zla_heamv.f"
	ky = 1 - (*n - 1) * *incy;
#line 262 "zla_heamv.f"
    }

/*     Set SAFE1 essentially to be the underflow threshold times the */
/*     number of additions in each row. */

#line 267 "zla_heamv.f"
    safe1 = dlamch_("Safe minimum", (ftnlen)12);
#line 268 "zla_heamv.f"
    safe1 = (*n + 1) * safe1;

/*     Form  y := alpha*abs(A)*abs(x) + beta*abs(y). */

/*     The O(N^2) SYMB_ZERO tests could be replaced by O(N) queries to */
/*     the inexact flag.  Still doesn't help change the iteration order */
/*     to per-column. */

#line 276 "zla_heamv.f"
    iy = ky;
#line 277 "zla_heamv.f"
    if (*incx == 1) {
#line 278 "zla_heamv.f"
	if (*uplo == ilauplo_("U", (ftnlen)1)) {
#line 279 "zla_heamv.f"
	    i__1 = *n;
#line 279 "zla_heamv.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 280 "zla_heamv.f"
		if (*beta == 0.) {
#line 281 "zla_heamv.f"
		    symb_zero__ = TRUE_;
#line 282 "zla_heamv.f"
		    y[iy] = 0.;
#line 283 "zla_heamv.f"
		} else if (y[iy] == 0.) {
#line 284 "zla_heamv.f"
		    symb_zero__ = TRUE_;
#line 285 "zla_heamv.f"
		} else {
#line 286 "zla_heamv.f"
		    symb_zero__ = FALSE_;
#line 287 "zla_heamv.f"
		    y[iy] = *beta * (d__1 = y[iy], abs(d__1));
#line 288 "zla_heamv.f"
		}
#line 289 "zla_heamv.f"
		if (*alpha != 0.) {
#line 290 "zla_heamv.f"
		    i__2 = i__;
#line 290 "zla_heamv.f"
		    for (j = 1; j <= i__2; ++j) {
#line 291 "zla_heamv.f"
			i__3 = j + i__ * a_dim1;
#line 291 "zla_heamv.f"
			temp = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(
				&a[j + i__ * a_dim1]), abs(d__2));
#line 292 "zla_heamv.f"
			i__3 = j;
#line 292 "zla_heamv.f"
			symb_zero__ = symb_zero__ && (x[i__3].r == 0. && x[
				i__3].i == 0. || temp == 0.);
#line 295 "zla_heamv.f"
			i__3 = j;
#line 295 "zla_heamv.f"
			y[iy] += *alpha * ((d__1 = x[i__3].r, abs(d__1)) + (
				d__2 = d_imag(&x[j]), abs(d__2))) * temp;
#line 296 "zla_heamv.f"
		    }
#line 297 "zla_heamv.f"
		    i__2 = *n;
#line 297 "zla_heamv.f"
		    for (j = i__ + 1; j <= i__2; ++j) {
#line 298 "zla_heamv.f"
			i__3 = i__ + j * a_dim1;
#line 298 "zla_heamv.f"
			temp = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(
				&a[i__ + j * a_dim1]), abs(d__2));
#line 299 "zla_heamv.f"
			i__3 = j;
#line 299 "zla_heamv.f"
			symb_zero__ = symb_zero__ && (x[i__3].r == 0. && x[
				i__3].i == 0. || temp == 0.);
#line 302 "zla_heamv.f"
			i__3 = j;
#line 302 "zla_heamv.f"
			y[iy] += *alpha * ((d__1 = x[i__3].r, abs(d__1)) + (
				d__2 = d_imag(&x[j]), abs(d__2))) * temp;
#line 303 "zla_heamv.f"
		    }
#line 304 "zla_heamv.f"
		}
#line 306 "zla_heamv.f"
		if (! symb_zero__) {
#line 306 "zla_heamv.f"
		    y[iy] += d_sign(&safe1, &y[iy]);
#line 306 "zla_heamv.f"
		}
#line 309 "zla_heamv.f"
		iy += *incy;
#line 310 "zla_heamv.f"
	    }
#line 311 "zla_heamv.f"
	} else {
#line 312 "zla_heamv.f"
	    i__1 = *n;
#line 312 "zla_heamv.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 313 "zla_heamv.f"
		if (*beta == 0.) {
#line 314 "zla_heamv.f"
		    symb_zero__ = TRUE_;
#line 315 "zla_heamv.f"
		    y[iy] = 0.;
#line 316 "zla_heamv.f"
		} else if (y[iy] == 0.) {
#line 317 "zla_heamv.f"
		    symb_zero__ = TRUE_;
#line 318 "zla_heamv.f"
		} else {
#line 319 "zla_heamv.f"
		    symb_zero__ = FALSE_;
#line 320 "zla_heamv.f"
		    y[iy] = *beta * (d__1 = y[iy], abs(d__1));
#line 321 "zla_heamv.f"
		}
#line 322 "zla_heamv.f"
		if (*alpha != 0.) {
#line 323 "zla_heamv.f"
		    i__2 = i__;
#line 323 "zla_heamv.f"
		    for (j = 1; j <= i__2; ++j) {
#line 324 "zla_heamv.f"
			i__3 = i__ + j * a_dim1;
#line 324 "zla_heamv.f"
			temp = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(
				&a[i__ + j * a_dim1]), abs(d__2));
#line 325 "zla_heamv.f"
			i__3 = j;
#line 325 "zla_heamv.f"
			symb_zero__ = symb_zero__ && (x[i__3].r == 0. && x[
				i__3].i == 0. || temp == 0.);
#line 328 "zla_heamv.f"
			i__3 = j;
#line 328 "zla_heamv.f"
			y[iy] += *alpha * ((d__1 = x[i__3].r, abs(d__1)) + (
				d__2 = d_imag(&x[j]), abs(d__2))) * temp;
#line 329 "zla_heamv.f"
		    }
#line 330 "zla_heamv.f"
		    i__2 = *n;
#line 330 "zla_heamv.f"
		    for (j = i__ + 1; j <= i__2; ++j) {
#line 331 "zla_heamv.f"
			i__3 = j + i__ * a_dim1;
#line 331 "zla_heamv.f"
			temp = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(
				&a[j + i__ * a_dim1]), abs(d__2));
#line 332 "zla_heamv.f"
			i__3 = j;
#line 332 "zla_heamv.f"
			symb_zero__ = symb_zero__ && (x[i__3].r == 0. && x[
				i__3].i == 0. || temp == 0.);
#line 335 "zla_heamv.f"
			i__3 = j;
#line 335 "zla_heamv.f"
			y[iy] += *alpha * ((d__1 = x[i__3].r, abs(d__1)) + (
				d__2 = d_imag(&x[j]), abs(d__2))) * temp;
#line 336 "zla_heamv.f"
		    }
#line 337 "zla_heamv.f"
		}
#line 339 "zla_heamv.f"
		if (! symb_zero__) {
#line 339 "zla_heamv.f"
		    y[iy] += d_sign(&safe1, &y[iy]);
#line 339 "zla_heamv.f"
		}
#line 342 "zla_heamv.f"
		iy += *incy;
#line 343 "zla_heamv.f"
	    }
#line 344 "zla_heamv.f"
	}
#line 345 "zla_heamv.f"
    } else {
#line 346 "zla_heamv.f"
	if (*uplo == ilauplo_("U", (ftnlen)1)) {
#line 347 "zla_heamv.f"
	    i__1 = *n;
#line 347 "zla_heamv.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 348 "zla_heamv.f"
		if (*beta == 0.) {
#line 349 "zla_heamv.f"
		    symb_zero__ = TRUE_;
#line 350 "zla_heamv.f"
		    y[iy] = 0.;
#line 351 "zla_heamv.f"
		} else if (y[iy] == 0.) {
#line 352 "zla_heamv.f"
		    symb_zero__ = TRUE_;
#line 353 "zla_heamv.f"
		} else {
#line 354 "zla_heamv.f"
		    symb_zero__ = FALSE_;
#line 355 "zla_heamv.f"
		    y[iy] = *beta * (d__1 = y[iy], abs(d__1));
#line 356 "zla_heamv.f"
		}
#line 357 "zla_heamv.f"
		jx = kx;
#line 358 "zla_heamv.f"
		if (*alpha != 0.) {
#line 359 "zla_heamv.f"
		    i__2 = i__;
#line 359 "zla_heamv.f"
		    for (j = 1; j <= i__2; ++j) {
#line 360 "zla_heamv.f"
			i__3 = j + i__ * a_dim1;
#line 360 "zla_heamv.f"
			temp = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(
				&a[j + i__ * a_dim1]), abs(d__2));
#line 361 "zla_heamv.f"
			i__3 = j;
#line 361 "zla_heamv.f"
			symb_zero__ = symb_zero__ && (x[i__3].r == 0. && x[
				i__3].i == 0. || temp == 0.);
#line 364 "zla_heamv.f"
			i__3 = jx;
#line 364 "zla_heamv.f"
			y[iy] += *alpha * ((d__1 = x[i__3].r, abs(d__1)) + (
				d__2 = d_imag(&x[jx]), abs(d__2))) * temp;
#line 365 "zla_heamv.f"
			jx += *incx;
#line 366 "zla_heamv.f"
		    }
#line 367 "zla_heamv.f"
		    i__2 = *n;
#line 367 "zla_heamv.f"
		    for (j = i__ + 1; j <= i__2; ++j) {
#line 368 "zla_heamv.f"
			i__3 = i__ + j * a_dim1;
#line 368 "zla_heamv.f"
			temp = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(
				&a[i__ + j * a_dim1]), abs(d__2));
#line 369 "zla_heamv.f"
			i__3 = j;
#line 369 "zla_heamv.f"
			symb_zero__ = symb_zero__ && (x[i__3].r == 0. && x[
				i__3].i == 0. || temp == 0.);
#line 372 "zla_heamv.f"
			i__3 = jx;
#line 372 "zla_heamv.f"
			y[iy] += *alpha * ((d__1 = x[i__3].r, abs(d__1)) + (
				d__2 = d_imag(&x[jx]), abs(d__2))) * temp;
#line 373 "zla_heamv.f"
			jx += *incx;
#line 374 "zla_heamv.f"
		    }
#line 375 "zla_heamv.f"
		}
#line 377 "zla_heamv.f"
		if (! symb_zero__) {
#line 377 "zla_heamv.f"
		    y[iy] += d_sign(&safe1, &y[iy]);
#line 377 "zla_heamv.f"
		}
#line 380 "zla_heamv.f"
		iy += *incy;
#line 381 "zla_heamv.f"
	    }
#line 382 "zla_heamv.f"
	} else {
#line 383 "zla_heamv.f"
	    i__1 = *n;
#line 383 "zla_heamv.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 384 "zla_heamv.f"
		if (*beta == 0.) {
#line 385 "zla_heamv.f"
		    symb_zero__ = TRUE_;
#line 386 "zla_heamv.f"
		    y[iy] = 0.;
#line 387 "zla_heamv.f"
		} else if (y[iy] == 0.) {
#line 388 "zla_heamv.f"
		    symb_zero__ = TRUE_;
#line 389 "zla_heamv.f"
		} else {
#line 390 "zla_heamv.f"
		    symb_zero__ = FALSE_;
#line 391 "zla_heamv.f"
		    y[iy] = *beta * (d__1 = y[iy], abs(d__1));
#line 392 "zla_heamv.f"
		}
#line 393 "zla_heamv.f"
		jx = kx;
#line 394 "zla_heamv.f"
		if (*alpha != 0.) {
#line 395 "zla_heamv.f"
		    i__2 = i__;
#line 395 "zla_heamv.f"
		    for (j = 1; j <= i__2; ++j) {
#line 396 "zla_heamv.f"
			i__3 = i__ + j * a_dim1;
#line 396 "zla_heamv.f"
			temp = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(
				&a[i__ + j * a_dim1]), abs(d__2));
#line 397 "zla_heamv.f"
			i__3 = j;
#line 397 "zla_heamv.f"
			symb_zero__ = symb_zero__ && (x[i__3].r == 0. && x[
				i__3].i == 0. || temp == 0.);
#line 400 "zla_heamv.f"
			i__3 = jx;
#line 400 "zla_heamv.f"
			y[iy] += *alpha * ((d__1 = x[i__3].r, abs(d__1)) + (
				d__2 = d_imag(&x[jx]), abs(d__2))) * temp;
#line 401 "zla_heamv.f"
			jx += *incx;
#line 402 "zla_heamv.f"
		    }
#line 403 "zla_heamv.f"
		    i__2 = *n;
#line 403 "zla_heamv.f"
		    for (j = i__ + 1; j <= i__2; ++j) {
#line 404 "zla_heamv.f"
			i__3 = j + i__ * a_dim1;
#line 404 "zla_heamv.f"
			temp = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(
				&a[j + i__ * a_dim1]), abs(d__2));
#line 405 "zla_heamv.f"
			i__3 = j;
#line 405 "zla_heamv.f"
			symb_zero__ = symb_zero__ && (x[i__3].r == 0. && x[
				i__3].i == 0. || temp == 0.);
#line 408 "zla_heamv.f"
			i__3 = jx;
#line 408 "zla_heamv.f"
			y[iy] += *alpha * ((d__1 = x[i__3].r, abs(d__1)) + (
				d__2 = d_imag(&x[jx]), abs(d__2))) * temp;
#line 409 "zla_heamv.f"
			jx += *incx;
#line 410 "zla_heamv.f"
		    }
#line 411 "zla_heamv.f"
		}
#line 413 "zla_heamv.f"
		if (! symb_zero__) {
#line 413 "zla_heamv.f"
		    y[iy] += d_sign(&safe1, &y[iy]);
#line 413 "zla_heamv.f"
		}
#line 416 "zla_heamv.f"
		iy += *incy;
#line 417 "zla_heamv.f"
	    }
#line 418 "zla_heamv.f"
	}
#line 420 "zla_heamv.f"
    }

#line 422 "zla_heamv.f"
    return 0;

/*     End of ZLA_HEAMV */

} /* zla_heamv__ */

