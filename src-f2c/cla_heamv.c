#line 1 "cla_heamv.f"
/* cla_heamv.f -- translated by f2c (version 20100827).
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

#line 1 "cla_heamv.f"
/* > \brief \b CLA_HEAMV computes a matrix-vector product using a Hermitian indefinite matrix to calculate err
or bounds. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLA_HEAMV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cla_hea
mv.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cla_hea
mv.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cla_hea
mv.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CLA_HEAMV( UPLO, N, ALPHA, A, LDA, X, INCX, BETA, Y, */
/*                             INCY ) */

/*       .. Scalar Arguments .. */
/*       REAL               ALPHA, BETA */
/*       INTEGER            INCX, INCY, LDA, N, UPLO */
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
/* >          A is COMPLEX array, dimension ( LDA, n ). */
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

/* > \date June 2017 */

/* > \ingroup complexHEcomputational */

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
/* Subroutine */ int cla_heamv__(integer *uplo, integer *n, doublereal *alpha,
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


/*  -- LAPACK computational routine (version 3.7.1) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     June 2017 */

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

#line 228 "cla_heamv.f"
    /* Parameter adjustments */
#line 228 "cla_heamv.f"
    a_dim1 = *lda;
#line 228 "cla_heamv.f"
    a_offset = 1 + a_dim1;
#line 228 "cla_heamv.f"
    a -= a_offset;
#line 228 "cla_heamv.f"
    --x;
#line 228 "cla_heamv.f"
    --y;
#line 228 "cla_heamv.f"

#line 228 "cla_heamv.f"
    /* Function Body */
#line 228 "cla_heamv.f"
    info = 0;
#line 229 "cla_heamv.f"
    if (*uplo != ilauplo_("U", (ftnlen)1) && *uplo != ilauplo_("L", (ftnlen)1)
	    ) {
#line 231 "cla_heamv.f"
	info = 1;
#line 232 "cla_heamv.f"
    } else if (*n < 0) {
#line 233 "cla_heamv.f"
	info = 2;
#line 234 "cla_heamv.f"
    } else if (*lda < max(1,*n)) {
#line 235 "cla_heamv.f"
	info = 5;
#line 236 "cla_heamv.f"
    } else if (*incx == 0) {
#line 237 "cla_heamv.f"
	info = 7;
#line 238 "cla_heamv.f"
    } else if (*incy == 0) {
#line 239 "cla_heamv.f"
	info = 10;
#line 240 "cla_heamv.f"
    }
#line 241 "cla_heamv.f"
    if (info != 0) {
#line 242 "cla_heamv.f"
	xerbla_("CHEMV ", &info, (ftnlen)6);
#line 243 "cla_heamv.f"
	return 0;
#line 244 "cla_heamv.f"
    }

/*     Quick return if possible. */

#line 248 "cla_heamv.f"
    if (*n == 0 || *alpha == 0. && *beta == 1.) {
#line 248 "cla_heamv.f"
	return 0;
#line 248 "cla_heamv.f"
    }

/*     Set up the start points in  X  and  Y. */

#line 253 "cla_heamv.f"
    if (*incx > 0) {
#line 254 "cla_heamv.f"
	kx = 1;
#line 255 "cla_heamv.f"
    } else {
#line 256 "cla_heamv.f"
	kx = 1 - (*n - 1) * *incx;
#line 257 "cla_heamv.f"
    }
#line 258 "cla_heamv.f"
    if (*incy > 0) {
#line 259 "cla_heamv.f"
	ky = 1;
#line 260 "cla_heamv.f"
    } else {
#line 261 "cla_heamv.f"
	ky = 1 - (*n - 1) * *incy;
#line 262 "cla_heamv.f"
    }

/*     Set SAFE1 essentially to be the underflow threshold times the */
/*     number of additions in each row. */

#line 267 "cla_heamv.f"
    safe1 = slamch_("Safe minimum", (ftnlen)12);
#line 268 "cla_heamv.f"
    safe1 = (*n + 1) * safe1;

/*     Form  y := alpha*abs(A)*abs(x) + beta*abs(y). */

/*     The O(N^2) SYMB_ZERO tests could be replaced by O(N) queries to */
/*     the inexact flag.  Still doesn't help change the iteration order */
/*     to per-column. */

#line 276 "cla_heamv.f"
    iy = ky;
#line 277 "cla_heamv.f"
    if (*incx == 1) {
#line 278 "cla_heamv.f"
	if (*uplo == ilauplo_("U", (ftnlen)1)) {
#line 279 "cla_heamv.f"
	    i__1 = *n;
#line 279 "cla_heamv.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 280 "cla_heamv.f"
		if (*beta == 0.) {
#line 281 "cla_heamv.f"
		    symb_zero__ = TRUE_;
#line 282 "cla_heamv.f"
		    y[iy] = 0.;
#line 283 "cla_heamv.f"
		} else if (y[iy] == 0.) {
#line 284 "cla_heamv.f"
		    symb_zero__ = TRUE_;
#line 285 "cla_heamv.f"
		} else {
#line 286 "cla_heamv.f"
		    symb_zero__ = FALSE_;
#line 287 "cla_heamv.f"
		    y[iy] = *beta * (d__1 = y[iy], abs(d__1));
#line 288 "cla_heamv.f"
		}
#line 289 "cla_heamv.f"
		if (*alpha != 0.) {
#line 290 "cla_heamv.f"
		    i__2 = i__;
#line 290 "cla_heamv.f"
		    for (j = 1; j <= i__2; ++j) {
#line 291 "cla_heamv.f"
			i__3 = j + i__ * a_dim1;
#line 291 "cla_heamv.f"
			temp = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(
				&a[j + i__ * a_dim1]), abs(d__2));
#line 292 "cla_heamv.f"
			i__3 = j;
#line 292 "cla_heamv.f"
			symb_zero__ = symb_zero__ && (x[i__3].r == 0. && x[
				i__3].i == 0. || temp == 0.);
#line 295 "cla_heamv.f"
			i__3 = j;
#line 295 "cla_heamv.f"
			y[iy] += *alpha * ((d__1 = x[i__3].r, abs(d__1)) + (
				d__2 = d_imag(&x[j]), abs(d__2))) * temp;
#line 296 "cla_heamv.f"
		    }
#line 297 "cla_heamv.f"
		    i__2 = *n;
#line 297 "cla_heamv.f"
		    for (j = i__ + 1; j <= i__2; ++j) {
#line 298 "cla_heamv.f"
			i__3 = i__ + j * a_dim1;
#line 298 "cla_heamv.f"
			temp = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(
				&a[i__ + j * a_dim1]), abs(d__2));
#line 299 "cla_heamv.f"
			i__3 = j;
#line 299 "cla_heamv.f"
			symb_zero__ = symb_zero__ && (x[i__3].r == 0. && x[
				i__3].i == 0. || temp == 0.);
#line 302 "cla_heamv.f"
			i__3 = j;
#line 302 "cla_heamv.f"
			y[iy] += *alpha * ((d__1 = x[i__3].r, abs(d__1)) + (
				d__2 = d_imag(&x[j]), abs(d__2))) * temp;
#line 303 "cla_heamv.f"
		    }
#line 304 "cla_heamv.f"
		}
#line 306 "cla_heamv.f"
		if (! symb_zero__) {
#line 306 "cla_heamv.f"
		    y[iy] += d_sign(&safe1, &y[iy]);
#line 306 "cla_heamv.f"
		}
#line 309 "cla_heamv.f"
		iy += *incy;
#line 310 "cla_heamv.f"
	    }
#line 311 "cla_heamv.f"
	} else {
#line 312 "cla_heamv.f"
	    i__1 = *n;
#line 312 "cla_heamv.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 313 "cla_heamv.f"
		if (*beta == 0.) {
#line 314 "cla_heamv.f"
		    symb_zero__ = TRUE_;
#line 315 "cla_heamv.f"
		    y[iy] = 0.;
#line 316 "cla_heamv.f"
		} else if (y[iy] == 0.) {
#line 317 "cla_heamv.f"
		    symb_zero__ = TRUE_;
#line 318 "cla_heamv.f"
		} else {
#line 319 "cla_heamv.f"
		    symb_zero__ = FALSE_;
#line 320 "cla_heamv.f"
		    y[iy] = *beta * (d__1 = y[iy], abs(d__1));
#line 321 "cla_heamv.f"
		}
#line 322 "cla_heamv.f"
		if (*alpha != 0.) {
#line 323 "cla_heamv.f"
		    i__2 = i__;
#line 323 "cla_heamv.f"
		    for (j = 1; j <= i__2; ++j) {
#line 324 "cla_heamv.f"
			i__3 = i__ + j * a_dim1;
#line 324 "cla_heamv.f"
			temp = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(
				&a[i__ + j * a_dim1]), abs(d__2));
#line 325 "cla_heamv.f"
			i__3 = j;
#line 325 "cla_heamv.f"
			symb_zero__ = symb_zero__ && (x[i__3].r == 0. && x[
				i__3].i == 0. || temp == 0.);
#line 328 "cla_heamv.f"
			i__3 = j;
#line 328 "cla_heamv.f"
			y[iy] += *alpha * ((d__1 = x[i__3].r, abs(d__1)) + (
				d__2 = d_imag(&x[j]), abs(d__2))) * temp;
#line 329 "cla_heamv.f"
		    }
#line 330 "cla_heamv.f"
		    i__2 = *n;
#line 330 "cla_heamv.f"
		    for (j = i__ + 1; j <= i__2; ++j) {
#line 331 "cla_heamv.f"
			i__3 = j + i__ * a_dim1;
#line 331 "cla_heamv.f"
			temp = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(
				&a[j + i__ * a_dim1]), abs(d__2));
#line 332 "cla_heamv.f"
			i__3 = j;
#line 332 "cla_heamv.f"
			symb_zero__ = symb_zero__ && (x[i__3].r == 0. && x[
				i__3].i == 0. || temp == 0.);
#line 335 "cla_heamv.f"
			i__3 = j;
#line 335 "cla_heamv.f"
			y[iy] += *alpha * ((d__1 = x[i__3].r, abs(d__1)) + (
				d__2 = d_imag(&x[j]), abs(d__2))) * temp;
#line 336 "cla_heamv.f"
		    }
#line 337 "cla_heamv.f"
		}
#line 339 "cla_heamv.f"
		if (! symb_zero__) {
#line 339 "cla_heamv.f"
		    y[iy] += d_sign(&safe1, &y[iy]);
#line 339 "cla_heamv.f"
		}
#line 342 "cla_heamv.f"
		iy += *incy;
#line 343 "cla_heamv.f"
	    }
#line 344 "cla_heamv.f"
	}
#line 345 "cla_heamv.f"
    } else {
#line 346 "cla_heamv.f"
	if (*uplo == ilauplo_("U", (ftnlen)1)) {
#line 347 "cla_heamv.f"
	    i__1 = *n;
#line 347 "cla_heamv.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 348 "cla_heamv.f"
		if (*beta == 0.) {
#line 349 "cla_heamv.f"
		    symb_zero__ = TRUE_;
#line 350 "cla_heamv.f"
		    y[iy] = 0.;
#line 351 "cla_heamv.f"
		} else if (y[iy] == 0.) {
#line 352 "cla_heamv.f"
		    symb_zero__ = TRUE_;
#line 353 "cla_heamv.f"
		} else {
#line 354 "cla_heamv.f"
		    symb_zero__ = FALSE_;
#line 355 "cla_heamv.f"
		    y[iy] = *beta * (d__1 = y[iy], abs(d__1));
#line 356 "cla_heamv.f"
		}
#line 357 "cla_heamv.f"
		jx = kx;
#line 358 "cla_heamv.f"
		if (*alpha != 0.) {
#line 359 "cla_heamv.f"
		    i__2 = i__;
#line 359 "cla_heamv.f"
		    for (j = 1; j <= i__2; ++j) {
#line 360 "cla_heamv.f"
			i__3 = j + i__ * a_dim1;
#line 360 "cla_heamv.f"
			temp = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(
				&a[j + i__ * a_dim1]), abs(d__2));
#line 361 "cla_heamv.f"
			i__3 = j;
#line 361 "cla_heamv.f"
			symb_zero__ = symb_zero__ && (x[i__3].r == 0. && x[
				i__3].i == 0. || temp == 0.);
#line 364 "cla_heamv.f"
			i__3 = jx;
#line 364 "cla_heamv.f"
			y[iy] += *alpha * ((d__1 = x[i__3].r, abs(d__1)) + (
				d__2 = d_imag(&x[jx]), abs(d__2))) * temp;
#line 365 "cla_heamv.f"
			jx += *incx;
#line 366 "cla_heamv.f"
		    }
#line 367 "cla_heamv.f"
		    i__2 = *n;
#line 367 "cla_heamv.f"
		    for (j = i__ + 1; j <= i__2; ++j) {
#line 368 "cla_heamv.f"
			i__3 = i__ + j * a_dim1;
#line 368 "cla_heamv.f"
			temp = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(
				&a[i__ + j * a_dim1]), abs(d__2));
#line 369 "cla_heamv.f"
			i__3 = j;
#line 369 "cla_heamv.f"
			symb_zero__ = symb_zero__ && (x[i__3].r == 0. && x[
				i__3].i == 0. || temp == 0.);
#line 372 "cla_heamv.f"
			i__3 = jx;
#line 372 "cla_heamv.f"
			y[iy] += *alpha * ((d__1 = x[i__3].r, abs(d__1)) + (
				d__2 = d_imag(&x[jx]), abs(d__2))) * temp;
#line 373 "cla_heamv.f"
			jx += *incx;
#line 374 "cla_heamv.f"
		    }
#line 375 "cla_heamv.f"
		}
#line 377 "cla_heamv.f"
		if (! symb_zero__) {
#line 377 "cla_heamv.f"
		    y[iy] += d_sign(&safe1, &y[iy]);
#line 377 "cla_heamv.f"
		}
#line 380 "cla_heamv.f"
		iy += *incy;
#line 381 "cla_heamv.f"
	    }
#line 382 "cla_heamv.f"
	} else {
#line 383 "cla_heamv.f"
	    i__1 = *n;
#line 383 "cla_heamv.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 384 "cla_heamv.f"
		if (*beta == 0.) {
#line 385 "cla_heamv.f"
		    symb_zero__ = TRUE_;
#line 386 "cla_heamv.f"
		    y[iy] = 0.;
#line 387 "cla_heamv.f"
		} else if (y[iy] == 0.) {
#line 388 "cla_heamv.f"
		    symb_zero__ = TRUE_;
#line 389 "cla_heamv.f"
		} else {
#line 390 "cla_heamv.f"
		    symb_zero__ = FALSE_;
#line 391 "cla_heamv.f"
		    y[iy] = *beta * (d__1 = y[iy], abs(d__1));
#line 392 "cla_heamv.f"
		}
#line 393 "cla_heamv.f"
		jx = kx;
#line 394 "cla_heamv.f"
		if (*alpha != 0.) {
#line 395 "cla_heamv.f"
		    i__2 = i__;
#line 395 "cla_heamv.f"
		    for (j = 1; j <= i__2; ++j) {
#line 396 "cla_heamv.f"
			i__3 = i__ + j * a_dim1;
#line 396 "cla_heamv.f"
			temp = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(
				&a[i__ + j * a_dim1]), abs(d__2));
#line 397 "cla_heamv.f"
			i__3 = j;
#line 397 "cla_heamv.f"
			symb_zero__ = symb_zero__ && (x[i__3].r == 0. && x[
				i__3].i == 0. || temp == 0.);
#line 400 "cla_heamv.f"
			i__3 = jx;
#line 400 "cla_heamv.f"
			y[iy] += *alpha * ((d__1 = x[i__3].r, abs(d__1)) + (
				d__2 = d_imag(&x[jx]), abs(d__2))) * temp;
#line 401 "cla_heamv.f"
			jx += *incx;
#line 402 "cla_heamv.f"
		    }
#line 403 "cla_heamv.f"
		    i__2 = *n;
#line 403 "cla_heamv.f"
		    for (j = i__ + 1; j <= i__2; ++j) {
#line 404 "cla_heamv.f"
			i__3 = j + i__ * a_dim1;
#line 404 "cla_heamv.f"
			temp = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(
				&a[j + i__ * a_dim1]), abs(d__2));
#line 405 "cla_heamv.f"
			i__3 = j;
#line 405 "cla_heamv.f"
			symb_zero__ = symb_zero__ && (x[i__3].r == 0. && x[
				i__3].i == 0. || temp == 0.);
#line 408 "cla_heamv.f"
			i__3 = jx;
#line 408 "cla_heamv.f"
			y[iy] += *alpha * ((d__1 = x[i__3].r, abs(d__1)) + (
				d__2 = d_imag(&x[jx]), abs(d__2))) * temp;
#line 409 "cla_heamv.f"
			jx += *incx;
#line 410 "cla_heamv.f"
		    }
#line 411 "cla_heamv.f"
		}
#line 413 "cla_heamv.f"
		if (! symb_zero__) {
#line 413 "cla_heamv.f"
		    y[iy] += d_sign(&safe1, &y[iy]);
#line 413 "cla_heamv.f"
		}
#line 416 "cla_heamv.f"
		iy += *incy;
#line 417 "cla_heamv.f"
	    }
#line 418 "cla_heamv.f"
	}
#line 420 "cla_heamv.f"
    }

#line 422 "cla_heamv.f"
    return 0;

/*     End of CLA_HEAMV */

} /* cla_heamv__ */

