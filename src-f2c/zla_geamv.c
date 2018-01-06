#line 1 "zla_geamv.f"
/* zla_geamv.f -- translated by f2c (version 20100827).
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

#line 1 "zla_geamv.f"
/* > \brief \b ZLA_GEAMV computes a matrix-vector product using a general matrix to calculate error bounds. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZLA_GEAMV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zla_gea
mv.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zla_gea
mv.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zla_gea
mv.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZLA_GEAMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, */
/*                              Y, INCY ) */

/*       .. Scalar Arguments .. */
/*       DOUBLE PRECISION   ALPHA, BETA */
/*       INTEGER            INCX, INCY, LDA, M, N */
/*       INTEGER            TRANS */
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
/* > ZLA_GEAMV  performs one of the matrix-vector operations */
/* > */
/* >         y := alpha*abs(A)*abs(x) + beta*abs(y), */
/* >    or   y := alpha*abs(A)**T*abs(x) + beta*abs(y), */
/* > */
/* > where alpha and beta are scalars, x and y are vectors and A is an */
/* > m by n matrix. */
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

/* > \param[in] TRANS */
/* > \verbatim */
/* >          TRANS is INTEGER */
/* >           On entry, TRANS specifies the operation to be performed as */
/* >           follows: */
/* > */
/* >             BLAS_NO_TRANS      y := alpha*abs(A)*abs(x) + beta*abs(y) */
/* >             BLAS_TRANS         y := alpha*abs(A**T)*abs(x) + beta*abs(y) */
/* >             BLAS_CONJ_TRANS    y := alpha*abs(A**T)*abs(x) + beta*abs(y) */
/* > */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >           On entry, M specifies the number of rows of the matrix A. */
/* >           M must be at least zero. */
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
/* >          ALPHA is DOUBLE PRECISION */
/* >           On entry, ALPHA specifies the scalar alpha. */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array of DIMENSION ( LDA, n ) */
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
/* >           max( 1, m ). */
/* >           Unchanged on exit. */
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
/* >          BETA is DOUBLE PRECISION */
/* >           On entry, BETA specifies the scalar beta. When BETA is */
/* >           supplied as zero then Y need not be set on input. */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Y */
/* > \verbatim */
/* >          Y is DOUBLE PRECISION array, dimension */
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
/* >           Unchanged on exit. */
/* > */
/* >  Level 2 Blas routine. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup complex16GEcomputational */

/*  ===================================================================== */
/* Subroutine */ int zla_geamv__(integer *trans, integer *m, integer *n, 
	doublereal *alpha, doublecomplex *a, integer *lda, doublecomplex *x, 
	integer *incx, doublereal *beta, doublereal *y, integer *incy)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Builtin functions */
    double d_imag(doublecomplex *), d_sign(doublereal *, doublereal *);

    /* Local variables */
    extern integer ilatrans_(char *, ftnlen);
    static integer i__, j;
    static logical symb_zero__;
    static integer iy, jx, kx, ky, info;
    static doublereal temp;
    static integer lenx, leny;
    static doublereal safe1;
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);


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

#line 226 "zla_geamv.f"
    /* Parameter adjustments */
#line 226 "zla_geamv.f"
    a_dim1 = *lda;
#line 226 "zla_geamv.f"
    a_offset = 1 + a_dim1;
#line 226 "zla_geamv.f"
    a -= a_offset;
#line 226 "zla_geamv.f"
    --x;
#line 226 "zla_geamv.f"
    --y;
#line 226 "zla_geamv.f"

#line 226 "zla_geamv.f"
    /* Function Body */
#line 226 "zla_geamv.f"
    info = 0;
#line 227 "zla_geamv.f"
    if (! (*trans == ilatrans_("N", (ftnlen)1) || *trans == ilatrans_("T", (
	    ftnlen)1) || *trans == ilatrans_("C", (ftnlen)1))) {
#line 230 "zla_geamv.f"
	info = 1;
#line 231 "zla_geamv.f"
    } else if (*m < 0) {
#line 232 "zla_geamv.f"
	info = 2;
#line 233 "zla_geamv.f"
    } else if (*n < 0) {
#line 234 "zla_geamv.f"
	info = 3;
#line 235 "zla_geamv.f"
    } else if (*lda < max(1,*m)) {
#line 236 "zla_geamv.f"
	info = 6;
#line 237 "zla_geamv.f"
    } else if (*incx == 0) {
#line 238 "zla_geamv.f"
	info = 8;
#line 239 "zla_geamv.f"
    } else if (*incy == 0) {
#line 240 "zla_geamv.f"
	info = 11;
#line 241 "zla_geamv.f"
    }
#line 242 "zla_geamv.f"
    if (info != 0) {
#line 243 "zla_geamv.f"
	xerbla_("ZLA_GEAMV ", &info, (ftnlen)10);
#line 244 "zla_geamv.f"
	return 0;
#line 245 "zla_geamv.f"
    }

/*     Quick return if possible. */

#line 249 "zla_geamv.f"
    if (*m == 0 || *n == 0 || *alpha == 0. && *beta == 1.) {
#line 249 "zla_geamv.f"
	return 0;
#line 249 "zla_geamv.f"
    }

/*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set */
/*     up the start points in  X  and  Y. */

#line 256 "zla_geamv.f"
    if (*trans == ilatrans_("N", (ftnlen)1)) {
#line 257 "zla_geamv.f"
	lenx = *n;
#line 258 "zla_geamv.f"
	leny = *m;
#line 259 "zla_geamv.f"
    } else {
#line 260 "zla_geamv.f"
	lenx = *m;
#line 261 "zla_geamv.f"
	leny = *n;
#line 262 "zla_geamv.f"
    }
#line 263 "zla_geamv.f"
    if (*incx > 0) {
#line 264 "zla_geamv.f"
	kx = 1;
#line 265 "zla_geamv.f"
    } else {
#line 266 "zla_geamv.f"
	kx = 1 - (lenx - 1) * *incx;
#line 267 "zla_geamv.f"
    }
#line 268 "zla_geamv.f"
    if (*incy > 0) {
#line 269 "zla_geamv.f"
	ky = 1;
#line 270 "zla_geamv.f"
    } else {
#line 271 "zla_geamv.f"
	ky = 1 - (leny - 1) * *incy;
#line 272 "zla_geamv.f"
    }

/*     Set SAFE1 essentially to be the underflow threshold times the */
/*     number of additions in each row. */

#line 277 "zla_geamv.f"
    safe1 = dlamch_("Safe minimum", (ftnlen)12);
#line 278 "zla_geamv.f"
    safe1 = (*n + 1) * safe1;

/*     Form  y := alpha*abs(A)*abs(x) + beta*abs(y). */

/*     The O(M*N) SYMB_ZERO tests could be replaced by O(N) queries to */
/*     the inexact flag.  Still doesn't help change the iteration order */
/*     to per-column. */

#line 286 "zla_geamv.f"
    iy = ky;
#line 287 "zla_geamv.f"
    if (*incx == 1) {
#line 288 "zla_geamv.f"
	if (*trans == ilatrans_("N", (ftnlen)1)) {
#line 289 "zla_geamv.f"
	    i__1 = leny;
#line 289 "zla_geamv.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 290 "zla_geamv.f"
		if (*beta == 0.) {
#line 291 "zla_geamv.f"
		    symb_zero__ = TRUE_;
#line 292 "zla_geamv.f"
		    y[iy] = 0.;
#line 293 "zla_geamv.f"
		} else if (y[iy] == 0.) {
#line 294 "zla_geamv.f"
		    symb_zero__ = TRUE_;
#line 295 "zla_geamv.f"
		} else {
#line 296 "zla_geamv.f"
		    symb_zero__ = FALSE_;
#line 297 "zla_geamv.f"
		    y[iy] = *beta * (d__1 = y[iy], abs(d__1));
#line 298 "zla_geamv.f"
		}
#line 299 "zla_geamv.f"
		if (*alpha != 0.) {
#line 300 "zla_geamv.f"
		    i__2 = lenx;
#line 300 "zla_geamv.f"
		    for (j = 1; j <= i__2; ++j) {
#line 301 "zla_geamv.f"
			i__3 = i__ + j * a_dim1;
#line 301 "zla_geamv.f"
			temp = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(
				&a[i__ + j * a_dim1]), abs(d__2));
#line 302 "zla_geamv.f"
			i__3 = j;
#line 302 "zla_geamv.f"
			symb_zero__ = symb_zero__ && (x[i__3].r == 0. && x[
				i__3].i == 0. || temp == 0.);
#line 305 "zla_geamv.f"
			i__3 = j;
#line 305 "zla_geamv.f"
			y[iy] += *alpha * ((d__1 = x[i__3].r, abs(d__1)) + (
				d__2 = d_imag(&x[j]), abs(d__2))) * temp;
#line 306 "zla_geamv.f"
		    }
#line 307 "zla_geamv.f"
		}
#line 309 "zla_geamv.f"
		if (! symb_zero__) {
#line 309 "zla_geamv.f"
		    y[iy] += d_sign(&safe1, &y[iy]);
#line 309 "zla_geamv.f"
		}
#line 312 "zla_geamv.f"
		iy += *incy;
#line 313 "zla_geamv.f"
	    }
#line 314 "zla_geamv.f"
	} else {
#line 315 "zla_geamv.f"
	    i__1 = leny;
#line 315 "zla_geamv.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 316 "zla_geamv.f"
		if (*beta == 0.) {
#line 317 "zla_geamv.f"
		    symb_zero__ = TRUE_;
#line 318 "zla_geamv.f"
		    y[iy] = 0.;
#line 319 "zla_geamv.f"
		} else if (y[iy] == 0.) {
#line 320 "zla_geamv.f"
		    symb_zero__ = TRUE_;
#line 321 "zla_geamv.f"
		} else {
#line 322 "zla_geamv.f"
		    symb_zero__ = FALSE_;
#line 323 "zla_geamv.f"
		    y[iy] = *beta * (d__1 = y[iy], abs(d__1));
#line 324 "zla_geamv.f"
		}
#line 325 "zla_geamv.f"
		if (*alpha != 0.) {
#line 326 "zla_geamv.f"
		    i__2 = lenx;
#line 326 "zla_geamv.f"
		    for (j = 1; j <= i__2; ++j) {
#line 327 "zla_geamv.f"
			i__3 = j + i__ * a_dim1;
#line 327 "zla_geamv.f"
			temp = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(
				&a[j + i__ * a_dim1]), abs(d__2));
#line 328 "zla_geamv.f"
			i__3 = j;
#line 328 "zla_geamv.f"
			symb_zero__ = symb_zero__ && (x[i__3].r == 0. && x[
				i__3].i == 0. || temp == 0.);
#line 331 "zla_geamv.f"
			i__3 = j;
#line 331 "zla_geamv.f"
			y[iy] += *alpha * ((d__1 = x[i__3].r, abs(d__1)) + (
				d__2 = d_imag(&x[j]), abs(d__2))) * temp;
#line 332 "zla_geamv.f"
		    }
#line 333 "zla_geamv.f"
		}
#line 335 "zla_geamv.f"
		if (! symb_zero__) {
#line 335 "zla_geamv.f"
		    y[iy] += d_sign(&safe1, &y[iy]);
#line 335 "zla_geamv.f"
		}
#line 338 "zla_geamv.f"
		iy += *incy;
#line 339 "zla_geamv.f"
	    }
#line 340 "zla_geamv.f"
	}
#line 341 "zla_geamv.f"
    } else {
#line 342 "zla_geamv.f"
	if (*trans == ilatrans_("N", (ftnlen)1)) {
#line 343 "zla_geamv.f"
	    i__1 = leny;
#line 343 "zla_geamv.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 344 "zla_geamv.f"
		if (*beta == 0.) {
#line 345 "zla_geamv.f"
		    symb_zero__ = TRUE_;
#line 346 "zla_geamv.f"
		    y[iy] = 0.;
#line 347 "zla_geamv.f"
		} else if (y[iy] == 0.) {
#line 348 "zla_geamv.f"
		    symb_zero__ = TRUE_;
#line 349 "zla_geamv.f"
		} else {
#line 350 "zla_geamv.f"
		    symb_zero__ = FALSE_;
#line 351 "zla_geamv.f"
		    y[iy] = *beta * (d__1 = y[iy], abs(d__1));
#line 352 "zla_geamv.f"
		}
#line 353 "zla_geamv.f"
		if (*alpha != 0.) {
#line 354 "zla_geamv.f"
		    jx = kx;
#line 355 "zla_geamv.f"
		    i__2 = lenx;
#line 355 "zla_geamv.f"
		    for (j = 1; j <= i__2; ++j) {
#line 356 "zla_geamv.f"
			i__3 = i__ + j * a_dim1;
#line 356 "zla_geamv.f"
			temp = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(
				&a[i__ + j * a_dim1]), abs(d__2));
#line 357 "zla_geamv.f"
			i__3 = jx;
#line 357 "zla_geamv.f"
			symb_zero__ = symb_zero__ && (x[i__3].r == 0. && x[
				i__3].i == 0. || temp == 0.);
#line 360 "zla_geamv.f"
			i__3 = jx;
#line 360 "zla_geamv.f"
			y[iy] += *alpha * ((d__1 = x[i__3].r, abs(d__1)) + (
				d__2 = d_imag(&x[jx]), abs(d__2))) * temp;
#line 361 "zla_geamv.f"
			jx += *incx;
#line 362 "zla_geamv.f"
		    }
#line 363 "zla_geamv.f"
		}
#line 365 "zla_geamv.f"
		if (! symb_zero__) {
#line 365 "zla_geamv.f"
		    y[iy] += d_sign(&safe1, &y[iy]);
#line 365 "zla_geamv.f"
		}
#line 368 "zla_geamv.f"
		iy += *incy;
#line 369 "zla_geamv.f"
	    }
#line 370 "zla_geamv.f"
	} else {
#line 371 "zla_geamv.f"
	    i__1 = leny;
#line 371 "zla_geamv.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 372 "zla_geamv.f"
		if (*beta == 0.) {
#line 373 "zla_geamv.f"
		    symb_zero__ = TRUE_;
#line 374 "zla_geamv.f"
		    y[iy] = 0.;
#line 375 "zla_geamv.f"
		} else if (y[iy] == 0.) {
#line 376 "zla_geamv.f"
		    symb_zero__ = TRUE_;
#line 377 "zla_geamv.f"
		} else {
#line 378 "zla_geamv.f"
		    symb_zero__ = FALSE_;
#line 379 "zla_geamv.f"
		    y[iy] = *beta * (d__1 = y[iy], abs(d__1));
#line 380 "zla_geamv.f"
		}
#line 381 "zla_geamv.f"
		if (*alpha != 0.) {
#line 382 "zla_geamv.f"
		    jx = kx;
#line 383 "zla_geamv.f"
		    i__2 = lenx;
#line 383 "zla_geamv.f"
		    for (j = 1; j <= i__2; ++j) {
#line 384 "zla_geamv.f"
			i__3 = j + i__ * a_dim1;
#line 384 "zla_geamv.f"
			temp = (d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(
				&a[j + i__ * a_dim1]), abs(d__2));
#line 385 "zla_geamv.f"
			i__3 = jx;
#line 385 "zla_geamv.f"
			symb_zero__ = symb_zero__ && (x[i__3].r == 0. && x[
				i__3].i == 0. || temp == 0.);
#line 388 "zla_geamv.f"
			i__3 = jx;
#line 388 "zla_geamv.f"
			y[iy] += *alpha * ((d__1 = x[i__3].r, abs(d__1)) + (
				d__2 = d_imag(&x[jx]), abs(d__2))) * temp;
#line 389 "zla_geamv.f"
			jx += *incx;
#line 390 "zla_geamv.f"
		    }
#line 391 "zla_geamv.f"
		}
#line 393 "zla_geamv.f"
		if (! symb_zero__) {
#line 393 "zla_geamv.f"
		    y[iy] += d_sign(&safe1, &y[iy]);
#line 393 "zla_geamv.f"
		}
#line 396 "zla_geamv.f"
		iy += *incy;
#line 397 "zla_geamv.f"
	    }
#line 398 "zla_geamv.f"
	}
#line 400 "zla_geamv.f"
    }

#line 402 "zla_geamv.f"
    return 0;

/*     End of ZLA_GEAMV */

} /* zla_geamv__ */

