#line 1 "cla_gbamv.f"
/* cla_gbamv.f -- translated by f2c (version 20100827).
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

#line 1 "cla_gbamv.f"
/* > \brief \b CLA_GBAMV performs a matrix-vector operation to calculate error bounds. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CLA_GBAMV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cla_gba
mv.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cla_gba
mv.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cla_gba
mv.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CLA_GBAMV( TRANS, M, N, KL, KU, ALPHA, AB, LDAB, X, */
/*                             INCX, BETA, Y, INCY ) */

/*       .. Scalar Arguments .. */
/*       REAL               ALPHA, BETA */
/*       INTEGER            INCX, INCY, LDAB, M, N, KL, KU, TRANS */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX            AB( LDAB, * ), X( * ) */
/*       REAL               Y( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CLA_GBAMV  performs one of the matrix-vector operations */
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
/* > \param[in] KL */
/* > \verbatim */
/* >          KL is INTEGER */
/* >           The number of subdiagonals within the band of A.  KL >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] KU */
/* > \verbatim */
/* >          KU is INTEGER */
/* >           The number of superdiagonals within the band of A.  KU >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] ALPHA */
/* > \verbatim */
/* >          ALPHA is REAL */
/* >           On entry, ALPHA specifies the scalar alpha. */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] AB */
/* > \verbatim */
/* >          AB is REAL array, dimension (LDAB,n) */
/* >           Before entry, the leading m by n part of the array AB must */
/* >           contain the matrix of coefficients. */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] LDAB */
/* > \verbatim */
/* >          LDAB is INTEGER */
/* >           On entry, LDAB specifies the first dimension of AB as declared */
/* >           in the calling (sub) program. LDAB must be at least */
/* >           max( 1, m ). */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] X */
/* > \verbatim */
/* >          X is REAL array, dimension */
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
/* >          BETA is REAL */
/* >           On entry, BETA specifies the scalar beta. When BETA is */
/* >           supplied as zero then Y need not be set on input. */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Y */
/* > \verbatim */
/* >          Y is REAL array, dimension */
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

/* > \date September 2012 */

/* > \ingroup complexGBcomputational */

/*  ===================================================================== */
/* Subroutine */ int cla_gbamv__(integer *trans, integer *m, integer *n, 
	integer *kl, integer *ku, doublereal *alpha, doublecomplex *ab, 
	integer *ldab, doublecomplex *x, integer *incx, doublereal *beta, 
	doublereal *y, integer *incy)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, i__1, i__2, i__3, i__4;
    doublereal d__1, d__2;

    /* Builtin functions */
    double d_imag(doublecomplex *), d_sign(doublereal *, doublereal *);

    /* Local variables */
    extern integer ilatrans_(char *, ftnlen);
    static integer i__, j;
    static logical symb_zero__;
    static integer kd, ke, iy, jx, kx, ky, info;
    static doublereal temp;
    static integer lenx, leny;
    static doublereal safe1;
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);


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
/*     .. Statement Functions */
/*     .. */
/*     .. Statement Function Definitions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 236 "cla_gbamv.f"
    /* Parameter adjustments */
#line 236 "cla_gbamv.f"
    ab_dim1 = *ldab;
#line 236 "cla_gbamv.f"
    ab_offset = 1 + ab_dim1;
#line 236 "cla_gbamv.f"
    ab -= ab_offset;
#line 236 "cla_gbamv.f"
    --x;
#line 236 "cla_gbamv.f"
    --y;
#line 236 "cla_gbamv.f"

#line 236 "cla_gbamv.f"
    /* Function Body */
#line 236 "cla_gbamv.f"
    info = 0;
#line 237 "cla_gbamv.f"
    if (! (*trans == ilatrans_("N", (ftnlen)1) || *trans == ilatrans_("T", (
	    ftnlen)1) || *trans == ilatrans_("C", (ftnlen)1))) {
#line 240 "cla_gbamv.f"
	info = 1;
#line 241 "cla_gbamv.f"
    } else if (*m < 0) {
#line 242 "cla_gbamv.f"
	info = 2;
#line 243 "cla_gbamv.f"
    } else if (*n < 0) {
#line 244 "cla_gbamv.f"
	info = 3;
#line 245 "cla_gbamv.f"
    } else if (*kl < 0 || *kl > *m - 1) {
#line 246 "cla_gbamv.f"
	info = 4;
#line 247 "cla_gbamv.f"
    } else if (*ku < 0 || *ku > *n - 1) {
#line 248 "cla_gbamv.f"
	info = 5;
#line 249 "cla_gbamv.f"
    } else if (*ldab < *kl + *ku + 1) {
#line 250 "cla_gbamv.f"
	info = 6;
#line 251 "cla_gbamv.f"
    } else if (*incx == 0) {
#line 252 "cla_gbamv.f"
	info = 8;
#line 253 "cla_gbamv.f"
    } else if (*incy == 0) {
#line 254 "cla_gbamv.f"
	info = 11;
#line 255 "cla_gbamv.f"
    }
#line 256 "cla_gbamv.f"
    if (info != 0) {
#line 257 "cla_gbamv.f"
	xerbla_("CLA_GBAMV ", &info, (ftnlen)10);
#line 258 "cla_gbamv.f"
	return 0;
#line 259 "cla_gbamv.f"
    }

/*     Quick return if possible. */

#line 263 "cla_gbamv.f"
    if (*m == 0 || *n == 0 || *alpha == 0. && *beta == 1.) {
#line 263 "cla_gbamv.f"
	return 0;
#line 263 "cla_gbamv.f"
    }

/*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set */
/*     up the start points in  X  and  Y. */

#line 270 "cla_gbamv.f"
    if (*trans == ilatrans_("N", (ftnlen)1)) {
#line 271 "cla_gbamv.f"
	lenx = *n;
#line 272 "cla_gbamv.f"
	leny = *m;
#line 273 "cla_gbamv.f"
    } else {
#line 274 "cla_gbamv.f"
	lenx = *m;
#line 275 "cla_gbamv.f"
	leny = *n;
#line 276 "cla_gbamv.f"
    }
#line 277 "cla_gbamv.f"
    if (*incx > 0) {
#line 278 "cla_gbamv.f"
	kx = 1;
#line 279 "cla_gbamv.f"
    } else {
#line 280 "cla_gbamv.f"
	kx = 1 - (lenx - 1) * *incx;
#line 281 "cla_gbamv.f"
    }
#line 282 "cla_gbamv.f"
    if (*incy > 0) {
#line 283 "cla_gbamv.f"
	ky = 1;
#line 284 "cla_gbamv.f"
    } else {
#line 285 "cla_gbamv.f"
	ky = 1 - (leny - 1) * *incy;
#line 286 "cla_gbamv.f"
    }

/*     Set SAFE1 essentially to be the underflow threshold times the */
/*     number of additions in each row. */

#line 291 "cla_gbamv.f"
    safe1 = slamch_("Safe minimum", (ftnlen)12);
#line 292 "cla_gbamv.f"
    safe1 = (*n + 1) * safe1;

/*     Form  y := alpha*abs(A)*abs(x) + beta*abs(y). */

/*     The O(M*N) SYMB_ZERO tests could be replaced by O(N) queries to */
/*     the inexact flag.  Still doesn't help change the iteration order */
/*     to per-column. */

#line 300 "cla_gbamv.f"
    kd = *ku + 1;
#line 301 "cla_gbamv.f"
    ke = *kl + 1;
#line 302 "cla_gbamv.f"
    iy = ky;
#line 303 "cla_gbamv.f"
    if (*incx == 1) {
#line 304 "cla_gbamv.f"
	if (*trans == ilatrans_("N", (ftnlen)1)) {
#line 305 "cla_gbamv.f"
	    i__1 = leny;
#line 305 "cla_gbamv.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 306 "cla_gbamv.f"
		if (*beta == 0.) {
#line 307 "cla_gbamv.f"
		    symb_zero__ = TRUE_;
#line 308 "cla_gbamv.f"
		    y[iy] = 0.;
#line 309 "cla_gbamv.f"
		} else if (y[iy] == 0.) {
#line 310 "cla_gbamv.f"
		    symb_zero__ = TRUE_;
#line 311 "cla_gbamv.f"
		} else {
#line 312 "cla_gbamv.f"
		    symb_zero__ = FALSE_;
#line 313 "cla_gbamv.f"
		    y[iy] = *beta * (d__1 = y[iy], abs(d__1));
#line 314 "cla_gbamv.f"
		}
#line 315 "cla_gbamv.f"
		if (*alpha != 0.) {
/* Computing MAX */
#line 316 "cla_gbamv.f"
		    i__2 = i__ - *kl;
/* Computing MIN */
#line 316 "cla_gbamv.f"
		    i__4 = i__ + *ku;
#line 316 "cla_gbamv.f"
		    i__3 = min(i__4,lenx);
#line 316 "cla_gbamv.f"
		    for (j = max(i__2,1); j <= i__3; ++j) {
#line 317 "cla_gbamv.f"
			i__2 = kd + i__ - j + j * ab_dim1;
#line 317 "cla_gbamv.f"
			temp = (d__1 = ab[i__2].r, abs(d__1)) + (d__2 = 
				d_imag(&ab[kd + i__ - j + j * ab_dim1]), abs(
				d__2));
#line 318 "cla_gbamv.f"
			i__2 = j;
#line 318 "cla_gbamv.f"
			symb_zero__ = symb_zero__ && (x[i__2].r == 0. && x[
				i__2].i == 0. || temp == 0.);
#line 321 "cla_gbamv.f"
			i__2 = j;
#line 321 "cla_gbamv.f"
			y[iy] += *alpha * ((d__1 = x[i__2].r, abs(d__1)) + (
				d__2 = d_imag(&x[j]), abs(d__2))) * temp;
#line 322 "cla_gbamv.f"
		    }
#line 323 "cla_gbamv.f"
		}
#line 325 "cla_gbamv.f"
		if (! symb_zero__) {
#line 325 "cla_gbamv.f"
		    y[iy] += d_sign(&safe1, &y[iy]);
#line 325 "cla_gbamv.f"
		}
#line 328 "cla_gbamv.f"
		iy += *incy;
#line 329 "cla_gbamv.f"
	    }
#line 330 "cla_gbamv.f"
	} else {
#line 331 "cla_gbamv.f"
	    i__1 = leny;
#line 331 "cla_gbamv.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 332 "cla_gbamv.f"
		if (*beta == 0.) {
#line 333 "cla_gbamv.f"
		    symb_zero__ = TRUE_;
#line 334 "cla_gbamv.f"
		    y[iy] = 0.;
#line 335 "cla_gbamv.f"
		} else if (y[iy] == 0.) {
#line 336 "cla_gbamv.f"
		    symb_zero__ = TRUE_;
#line 337 "cla_gbamv.f"
		} else {
#line 338 "cla_gbamv.f"
		    symb_zero__ = FALSE_;
#line 339 "cla_gbamv.f"
		    y[iy] = *beta * (d__1 = y[iy], abs(d__1));
#line 340 "cla_gbamv.f"
		}
#line 341 "cla_gbamv.f"
		if (*alpha != 0.) {
/* Computing MAX */
#line 342 "cla_gbamv.f"
		    i__3 = i__ - *kl;
/* Computing MIN */
#line 342 "cla_gbamv.f"
		    i__4 = i__ + *ku;
#line 342 "cla_gbamv.f"
		    i__2 = min(i__4,lenx);
#line 342 "cla_gbamv.f"
		    for (j = max(i__3,1); j <= i__2; ++j) {
#line 343 "cla_gbamv.f"
			i__3 = ke - i__ + j + i__ * ab_dim1;
#line 343 "cla_gbamv.f"
			temp = (d__1 = ab[i__3].r, abs(d__1)) + (d__2 = 
				d_imag(&ab[ke - i__ + j + i__ * ab_dim1]), 
				abs(d__2));
#line 344 "cla_gbamv.f"
			i__3 = j;
#line 344 "cla_gbamv.f"
			symb_zero__ = symb_zero__ && (x[i__3].r == 0. && x[
				i__3].i == 0. || temp == 0.);
#line 347 "cla_gbamv.f"
			i__3 = j;
#line 347 "cla_gbamv.f"
			y[iy] += *alpha * ((d__1 = x[i__3].r, abs(d__1)) + (
				d__2 = d_imag(&x[j]), abs(d__2))) * temp;
#line 348 "cla_gbamv.f"
		    }
#line 349 "cla_gbamv.f"
		}
#line 351 "cla_gbamv.f"
		if (! symb_zero__) {
#line 351 "cla_gbamv.f"
		    y[iy] += d_sign(&safe1, &y[iy]);
#line 351 "cla_gbamv.f"
		}
#line 354 "cla_gbamv.f"
		iy += *incy;
#line 355 "cla_gbamv.f"
	    }
#line 356 "cla_gbamv.f"
	}
#line 357 "cla_gbamv.f"
    } else {
#line 358 "cla_gbamv.f"
	if (*trans == ilatrans_("N", (ftnlen)1)) {
#line 359 "cla_gbamv.f"
	    i__1 = leny;
#line 359 "cla_gbamv.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 360 "cla_gbamv.f"
		if (*beta == 0.) {
#line 361 "cla_gbamv.f"
		    symb_zero__ = TRUE_;
#line 362 "cla_gbamv.f"
		    y[iy] = 0.;
#line 363 "cla_gbamv.f"
		} else if (y[iy] == 0.) {
#line 364 "cla_gbamv.f"
		    symb_zero__ = TRUE_;
#line 365 "cla_gbamv.f"
		} else {
#line 366 "cla_gbamv.f"
		    symb_zero__ = FALSE_;
#line 367 "cla_gbamv.f"
		    y[iy] = *beta * (d__1 = y[iy], abs(d__1));
#line 368 "cla_gbamv.f"
		}
#line 369 "cla_gbamv.f"
		if (*alpha != 0.) {
#line 370 "cla_gbamv.f"
		    jx = kx;
/* Computing MAX */
#line 371 "cla_gbamv.f"
		    i__2 = i__ - *kl;
/* Computing MIN */
#line 371 "cla_gbamv.f"
		    i__4 = i__ + *ku;
#line 371 "cla_gbamv.f"
		    i__3 = min(i__4,lenx);
#line 371 "cla_gbamv.f"
		    for (j = max(i__2,1); j <= i__3; ++j) {
#line 372 "cla_gbamv.f"
			i__2 = kd + i__ - j + j * ab_dim1;
#line 372 "cla_gbamv.f"
			temp = (d__1 = ab[i__2].r, abs(d__1)) + (d__2 = 
				d_imag(&ab[kd + i__ - j + j * ab_dim1]), abs(
				d__2));
#line 373 "cla_gbamv.f"
			i__2 = jx;
#line 373 "cla_gbamv.f"
			symb_zero__ = symb_zero__ && (x[i__2].r == 0. && x[
				i__2].i == 0. || temp == 0.);
#line 376 "cla_gbamv.f"
			i__2 = jx;
#line 376 "cla_gbamv.f"
			y[iy] += *alpha * ((d__1 = x[i__2].r, abs(d__1)) + (
				d__2 = d_imag(&x[jx]), abs(d__2))) * temp;
#line 377 "cla_gbamv.f"
			jx += *incx;
#line 378 "cla_gbamv.f"
		    }
#line 379 "cla_gbamv.f"
		}
#line 381 "cla_gbamv.f"
		if (! symb_zero__) {
#line 381 "cla_gbamv.f"
		    y[iy] += d_sign(&safe1, &y[iy]);
#line 381 "cla_gbamv.f"
		}
#line 384 "cla_gbamv.f"
		iy += *incy;
#line 385 "cla_gbamv.f"
	    }
#line 386 "cla_gbamv.f"
	} else {
#line 387 "cla_gbamv.f"
	    i__1 = leny;
#line 387 "cla_gbamv.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 388 "cla_gbamv.f"
		if (*beta == 0.) {
#line 389 "cla_gbamv.f"
		    symb_zero__ = TRUE_;
#line 390 "cla_gbamv.f"
		    y[iy] = 0.;
#line 391 "cla_gbamv.f"
		} else if (y[iy] == 0.) {
#line 392 "cla_gbamv.f"
		    symb_zero__ = TRUE_;
#line 393 "cla_gbamv.f"
		} else {
#line 394 "cla_gbamv.f"
		    symb_zero__ = FALSE_;
#line 395 "cla_gbamv.f"
		    y[iy] = *beta * (d__1 = y[iy], abs(d__1));
#line 396 "cla_gbamv.f"
		}
#line 397 "cla_gbamv.f"
		if (*alpha != 0.) {
#line 398 "cla_gbamv.f"
		    jx = kx;
/* Computing MAX */
#line 399 "cla_gbamv.f"
		    i__3 = i__ - *kl;
/* Computing MIN */
#line 399 "cla_gbamv.f"
		    i__4 = i__ + *ku;
#line 399 "cla_gbamv.f"
		    i__2 = min(i__4,lenx);
#line 399 "cla_gbamv.f"
		    for (j = max(i__3,1); j <= i__2; ++j) {
#line 400 "cla_gbamv.f"
			i__3 = ke - i__ + j + i__ * ab_dim1;
#line 400 "cla_gbamv.f"
			temp = (d__1 = ab[i__3].r, abs(d__1)) + (d__2 = 
				d_imag(&ab[ke - i__ + j + i__ * ab_dim1]), 
				abs(d__2));
#line 401 "cla_gbamv.f"
			i__3 = jx;
#line 401 "cla_gbamv.f"
			symb_zero__ = symb_zero__ && (x[i__3].r == 0. && x[
				i__3].i == 0. || temp == 0.);
#line 404 "cla_gbamv.f"
			i__3 = jx;
#line 404 "cla_gbamv.f"
			y[iy] += *alpha * ((d__1 = x[i__3].r, abs(d__1)) + (
				d__2 = d_imag(&x[jx]), abs(d__2))) * temp;
#line 405 "cla_gbamv.f"
			jx += *incx;
#line 406 "cla_gbamv.f"
		    }
#line 407 "cla_gbamv.f"
		}
#line 409 "cla_gbamv.f"
		if (! symb_zero__) {
#line 409 "cla_gbamv.f"
		    y[iy] += d_sign(&safe1, &y[iy]);
#line 409 "cla_gbamv.f"
		}
#line 412 "cla_gbamv.f"
		iy += *incy;
#line 413 "cla_gbamv.f"
	    }
#line 414 "cla_gbamv.f"
	}
#line 416 "cla_gbamv.f"
    }

#line 418 "cla_gbamv.f"
    return 0;

/*     End of CLA_GBAMV */

} /* cla_gbamv__ */

