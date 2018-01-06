#line 1 "dla_gbamv.f"
/* dla_gbamv.f -- translated by f2c (version 20100827).
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

#line 1 "dla_gbamv.f"
/* > \brief \b DLA_GBAMV performs a matrix-vector operation to calculate error bounds. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLA_GBAMV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dla_gba
mv.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dla_gba
mv.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dla_gba
mv.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLA_GBAMV( TRANS, M, N, KL, KU, ALPHA, AB, LDAB, X, */
/*                             INCX, BETA, Y, INCY ) */

/*       .. Scalar Arguments .. */
/*       DOUBLE PRECISION   ALPHA, BETA */
/*       INTEGER            INCX, INCY, LDAB, M, N, KL, KU, TRANS */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   AB( LDAB, * ), X( * ), Y( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLA_GBAMV  performs one of the matrix-vector operations */
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
/* >          ALPHA is DOUBLE PRECISION */
/* >           On entry, ALPHA specifies the scalar alpha. */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] AB */
/* > \verbatim */
/* >          AB is DOUBLE PRECISION array, dimension ( LDAB, n ) */
/* >           Before entry, the leading m by n part of the array AB must */
/* >           contain the matrix of coefficients. */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] LDAB */
/* > \verbatim */
/* >          LDAB is INTEGER */
/* >           On entry, LDA specifies the first dimension of AB as declared */
/* >           in the calling (sub) program. LDAB must be at least */
/* >           max( 1, m ). */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] X */
/* > \verbatim */
/* >          X is DOUBLE PRECISION array, dimension */
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

/* > \date June 2017 */

/* > \ingroup doubleGBcomputational */

/*  ===================================================================== */
/* Subroutine */ int dla_gbamv__(integer *trans, integer *m, integer *n, 
	integer *kl, integer *ku, doublereal *alpha, doublereal *ab, integer *
	ldab, doublereal *x, integer *incx, doublereal *beta, doublereal *y, 
	integer *incy)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, i__1, i__2, i__3, i__4;
    doublereal d__1;

    /* Builtin functions */
    double d_sign(doublereal *, doublereal *);

    /* Local variables */
    extern integer ilatrans_(char *, ftnlen);
    static integer i__, j;
    static logical symb_zero__;
    static integer kd, ke, iy, jx, kx, ky, info;
    static doublereal temp;
    static integer lenx, leny;
    static doublereal safe1;
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);


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
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 227 "dla_gbamv.f"
    /* Parameter adjustments */
#line 227 "dla_gbamv.f"
    ab_dim1 = *ldab;
#line 227 "dla_gbamv.f"
    ab_offset = 1 + ab_dim1;
#line 227 "dla_gbamv.f"
    ab -= ab_offset;
#line 227 "dla_gbamv.f"
    --x;
#line 227 "dla_gbamv.f"
    --y;
#line 227 "dla_gbamv.f"

#line 227 "dla_gbamv.f"
    /* Function Body */
#line 227 "dla_gbamv.f"
    info = 0;
#line 228 "dla_gbamv.f"
    if (! (*trans == ilatrans_("N", (ftnlen)1) || *trans == ilatrans_("T", (
	    ftnlen)1) || *trans == ilatrans_("C", (ftnlen)1))) {
#line 231 "dla_gbamv.f"
	info = 1;
#line 232 "dla_gbamv.f"
    } else if (*m < 0) {
#line 233 "dla_gbamv.f"
	info = 2;
#line 234 "dla_gbamv.f"
    } else if (*n < 0) {
#line 235 "dla_gbamv.f"
	info = 3;
#line 236 "dla_gbamv.f"
    } else if (*kl < 0 || *kl > *m - 1) {
#line 237 "dla_gbamv.f"
	info = 4;
#line 238 "dla_gbamv.f"
    } else if (*ku < 0 || *ku > *n - 1) {
#line 239 "dla_gbamv.f"
	info = 5;
#line 240 "dla_gbamv.f"
    } else if (*ldab < *kl + *ku + 1) {
#line 241 "dla_gbamv.f"
	info = 6;
#line 242 "dla_gbamv.f"
    } else if (*incx == 0) {
#line 243 "dla_gbamv.f"
	info = 8;
#line 244 "dla_gbamv.f"
    } else if (*incy == 0) {
#line 245 "dla_gbamv.f"
	info = 11;
#line 246 "dla_gbamv.f"
    }
#line 247 "dla_gbamv.f"
    if (info != 0) {
#line 248 "dla_gbamv.f"
	xerbla_("DLA_GBAMV ", &info, (ftnlen)10);
#line 249 "dla_gbamv.f"
	return 0;
#line 250 "dla_gbamv.f"
    }

/*     Quick return if possible. */

#line 254 "dla_gbamv.f"
    if (*m == 0 || *n == 0 || *alpha == 0. && *beta == 1.) {
#line 254 "dla_gbamv.f"
	return 0;
#line 254 "dla_gbamv.f"
    }

/*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set */
/*     up the start points in  X  and  Y. */

#line 261 "dla_gbamv.f"
    if (*trans == ilatrans_("N", (ftnlen)1)) {
#line 262 "dla_gbamv.f"
	lenx = *n;
#line 263 "dla_gbamv.f"
	leny = *m;
#line 264 "dla_gbamv.f"
    } else {
#line 265 "dla_gbamv.f"
	lenx = *m;
#line 266 "dla_gbamv.f"
	leny = *n;
#line 267 "dla_gbamv.f"
    }
#line 268 "dla_gbamv.f"
    if (*incx > 0) {
#line 269 "dla_gbamv.f"
	kx = 1;
#line 270 "dla_gbamv.f"
    } else {
#line 271 "dla_gbamv.f"
	kx = 1 - (lenx - 1) * *incx;
#line 272 "dla_gbamv.f"
    }
#line 273 "dla_gbamv.f"
    if (*incy > 0) {
#line 274 "dla_gbamv.f"
	ky = 1;
#line 275 "dla_gbamv.f"
    } else {
#line 276 "dla_gbamv.f"
	ky = 1 - (leny - 1) * *incy;
#line 277 "dla_gbamv.f"
    }

/*     Set SAFE1 essentially to be the underflow threshold times the */
/*     number of additions in each row. */

#line 282 "dla_gbamv.f"
    safe1 = dlamch_("Safe minimum", (ftnlen)12);
#line 283 "dla_gbamv.f"
    safe1 = (*n + 1) * safe1;

/*     Form  y := alpha*abs(A)*abs(x) + beta*abs(y). */

/*     The O(M*N) SYMB_ZERO tests could be replaced by O(N) queries to */
/*     the inexact flag.  Still doesn't help change the iteration order */
/*     to per-column. */

#line 291 "dla_gbamv.f"
    kd = *ku + 1;
#line 292 "dla_gbamv.f"
    ke = *kl + 1;
#line 293 "dla_gbamv.f"
    iy = ky;
#line 294 "dla_gbamv.f"
    if (*incx == 1) {
#line 295 "dla_gbamv.f"
	if (*trans == ilatrans_("N", (ftnlen)1)) {
#line 296 "dla_gbamv.f"
	    i__1 = leny;
#line 296 "dla_gbamv.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 297 "dla_gbamv.f"
		if (*beta == 0.) {
#line 298 "dla_gbamv.f"
		    symb_zero__ = TRUE_;
#line 299 "dla_gbamv.f"
		    y[iy] = 0.;
#line 300 "dla_gbamv.f"
		} else if (y[iy] == 0.) {
#line 301 "dla_gbamv.f"
		    symb_zero__ = TRUE_;
#line 302 "dla_gbamv.f"
		} else {
#line 303 "dla_gbamv.f"
		    symb_zero__ = FALSE_;
#line 304 "dla_gbamv.f"
		    y[iy] = *beta * (d__1 = y[iy], abs(d__1));
#line 305 "dla_gbamv.f"
		}
#line 306 "dla_gbamv.f"
		if (*alpha != 0.) {
/* Computing MAX */
#line 307 "dla_gbamv.f"
		    i__2 = i__ - *kl;
/* Computing MIN */
#line 307 "dla_gbamv.f"
		    i__4 = i__ + *ku;
#line 307 "dla_gbamv.f"
		    i__3 = min(i__4,lenx);
#line 307 "dla_gbamv.f"
		    for (j = max(i__2,1); j <= i__3; ++j) {
#line 308 "dla_gbamv.f"
			temp = (d__1 = ab[kd + i__ - j + j * ab_dim1], abs(
				d__1));
#line 309 "dla_gbamv.f"
			symb_zero__ = symb_zero__ && (x[j] == 0. || temp == 
				0.);
#line 312 "dla_gbamv.f"
			y[iy] += *alpha * (d__1 = x[j], abs(d__1)) * temp;
#line 313 "dla_gbamv.f"
		    }
#line 314 "dla_gbamv.f"
		}
#line 316 "dla_gbamv.f"
		if (! symb_zero__) {
#line 316 "dla_gbamv.f"
		    y[iy] += d_sign(&safe1, &y[iy]);
#line 316 "dla_gbamv.f"
		}
#line 318 "dla_gbamv.f"
		iy += *incy;
#line 319 "dla_gbamv.f"
	    }
#line 320 "dla_gbamv.f"
	} else {
#line 321 "dla_gbamv.f"
	    i__1 = leny;
#line 321 "dla_gbamv.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 322 "dla_gbamv.f"
		if (*beta == 0.) {
#line 323 "dla_gbamv.f"
		    symb_zero__ = TRUE_;
#line 324 "dla_gbamv.f"
		    y[iy] = 0.;
#line 325 "dla_gbamv.f"
		} else if (y[iy] == 0.) {
#line 326 "dla_gbamv.f"
		    symb_zero__ = TRUE_;
#line 327 "dla_gbamv.f"
		} else {
#line 328 "dla_gbamv.f"
		    symb_zero__ = FALSE_;
#line 329 "dla_gbamv.f"
		    y[iy] = *beta * (d__1 = y[iy], abs(d__1));
#line 330 "dla_gbamv.f"
		}
#line 331 "dla_gbamv.f"
		if (*alpha != 0.) {
/* Computing MAX */
#line 332 "dla_gbamv.f"
		    i__3 = i__ - *kl;
/* Computing MIN */
#line 332 "dla_gbamv.f"
		    i__4 = i__ + *ku;
#line 332 "dla_gbamv.f"
		    i__2 = min(i__4,lenx);
#line 332 "dla_gbamv.f"
		    for (j = max(i__3,1); j <= i__2; ++j) {
#line 333 "dla_gbamv.f"
			temp = (d__1 = ab[ke - i__ + j + i__ * ab_dim1], abs(
				d__1));
#line 334 "dla_gbamv.f"
			symb_zero__ = symb_zero__ && (x[j] == 0. || temp == 
				0.);
#line 337 "dla_gbamv.f"
			y[iy] += *alpha * (d__1 = x[j], abs(d__1)) * temp;
#line 338 "dla_gbamv.f"
		    }
#line 339 "dla_gbamv.f"
		}
#line 341 "dla_gbamv.f"
		if (! symb_zero__) {
#line 341 "dla_gbamv.f"
		    y[iy] += d_sign(&safe1, &y[iy]);
#line 341 "dla_gbamv.f"
		}
#line 343 "dla_gbamv.f"
		iy += *incy;
#line 344 "dla_gbamv.f"
	    }
#line 345 "dla_gbamv.f"
	}
#line 346 "dla_gbamv.f"
    } else {
#line 347 "dla_gbamv.f"
	if (*trans == ilatrans_("N", (ftnlen)1)) {
#line 348 "dla_gbamv.f"
	    i__1 = leny;
#line 348 "dla_gbamv.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 349 "dla_gbamv.f"
		if (*beta == 0.) {
#line 350 "dla_gbamv.f"
		    symb_zero__ = TRUE_;
#line 351 "dla_gbamv.f"
		    y[iy] = 0.;
#line 352 "dla_gbamv.f"
		} else if (y[iy] == 0.) {
#line 353 "dla_gbamv.f"
		    symb_zero__ = TRUE_;
#line 354 "dla_gbamv.f"
		} else {
#line 355 "dla_gbamv.f"
		    symb_zero__ = FALSE_;
#line 356 "dla_gbamv.f"
		    y[iy] = *beta * (d__1 = y[iy], abs(d__1));
#line 357 "dla_gbamv.f"
		}
#line 358 "dla_gbamv.f"
		if (*alpha != 0.) {
#line 359 "dla_gbamv.f"
		    jx = kx;
/* Computing MAX */
#line 360 "dla_gbamv.f"
		    i__2 = i__ - *kl;
/* Computing MIN */
#line 360 "dla_gbamv.f"
		    i__4 = i__ + *ku;
#line 360 "dla_gbamv.f"
		    i__3 = min(i__4,lenx);
#line 360 "dla_gbamv.f"
		    for (j = max(i__2,1); j <= i__3; ++j) {
#line 361 "dla_gbamv.f"
			temp = (d__1 = ab[kd + i__ - j + j * ab_dim1], abs(
				d__1));
#line 362 "dla_gbamv.f"
			symb_zero__ = symb_zero__ && (x[jx] == 0. || temp == 
				0.);
#line 365 "dla_gbamv.f"
			y[iy] += *alpha * (d__1 = x[jx], abs(d__1)) * temp;
#line 366 "dla_gbamv.f"
			jx += *incx;
#line 367 "dla_gbamv.f"
		    }
#line 368 "dla_gbamv.f"
		}
#line 370 "dla_gbamv.f"
		if (! symb_zero__) {
#line 370 "dla_gbamv.f"
		    y[iy] += d_sign(&safe1, &y[iy]);
#line 370 "dla_gbamv.f"
		}
#line 373 "dla_gbamv.f"
		iy += *incy;
#line 374 "dla_gbamv.f"
	    }
#line 375 "dla_gbamv.f"
	} else {
#line 376 "dla_gbamv.f"
	    i__1 = leny;
#line 376 "dla_gbamv.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 377 "dla_gbamv.f"
		if (*beta == 0.) {
#line 378 "dla_gbamv.f"
		    symb_zero__ = TRUE_;
#line 379 "dla_gbamv.f"
		    y[iy] = 0.;
#line 380 "dla_gbamv.f"
		} else if (y[iy] == 0.) {
#line 381 "dla_gbamv.f"
		    symb_zero__ = TRUE_;
#line 382 "dla_gbamv.f"
		} else {
#line 383 "dla_gbamv.f"
		    symb_zero__ = FALSE_;
#line 384 "dla_gbamv.f"
		    y[iy] = *beta * (d__1 = y[iy], abs(d__1));
#line 385 "dla_gbamv.f"
		}
#line 386 "dla_gbamv.f"
		if (*alpha != 0.) {
#line 387 "dla_gbamv.f"
		    jx = kx;
/* Computing MAX */
#line 388 "dla_gbamv.f"
		    i__3 = i__ - *kl;
/* Computing MIN */
#line 388 "dla_gbamv.f"
		    i__4 = i__ + *ku;
#line 388 "dla_gbamv.f"
		    i__2 = min(i__4,lenx);
#line 388 "dla_gbamv.f"
		    for (j = max(i__3,1); j <= i__2; ++j) {
#line 389 "dla_gbamv.f"
			temp = (d__1 = ab[ke - i__ + j + i__ * ab_dim1], abs(
				d__1));
#line 390 "dla_gbamv.f"
			symb_zero__ = symb_zero__ && (x[jx] == 0. || temp == 
				0.);
#line 393 "dla_gbamv.f"
			y[iy] += *alpha * (d__1 = x[jx], abs(d__1)) * temp;
#line 394 "dla_gbamv.f"
			jx += *incx;
#line 395 "dla_gbamv.f"
		    }
#line 396 "dla_gbamv.f"
		}
#line 398 "dla_gbamv.f"
		if (! symb_zero__) {
#line 398 "dla_gbamv.f"
		    y[iy] += d_sign(&safe1, &y[iy]);
#line 398 "dla_gbamv.f"
		}
#line 401 "dla_gbamv.f"
		iy += *incy;
#line 402 "dla_gbamv.f"
	    }
#line 403 "dla_gbamv.f"
	}
#line 405 "dla_gbamv.f"
    }

#line 407 "dla_gbamv.f"
    return 0;

/*     End of DLA_GBAMV */

} /* dla_gbamv__ */

