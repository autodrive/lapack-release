#line 1 "sla_gbamv.f"
/* sla_gbamv.f -- translated by f2c (version 20100827).
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

#line 1 "sla_gbamv.f"
/* > \brief \b SLA_GBAMV performs a matrix-vector operation to calculate error bounds. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLA_GBAMV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sla_gba
mv.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sla_gba
mv.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sla_gba
mv.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLA_GBAMV( TRANS, M, N, KL, KU, ALPHA, AB, LDAB, X, */
/*                             INCX, BETA, Y, INCY ) */

/*       .. Scalar Arguments .. */
/*       REAL               ALPHA, BETA */
/*       INTEGER            INCX, INCY, LDAB, M, N, KL, KU, TRANS */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               AB( LDAB, * ), X( * ), Y( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLA_GBAMV  performs one of the matrix-vector operations */
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
/* >          AB is REAL array, dimension ( LDAB, n ) */
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

/* > \date June 2017 */

/* > \ingroup realGBcomputational */

/*  ===================================================================== */
/* Subroutine */ int sla_gbamv__(integer *trans, integer *m, integer *n, 
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
    extern doublereal slamch_(char *, ftnlen);
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

#line 226 "sla_gbamv.f"
    /* Parameter adjustments */
#line 226 "sla_gbamv.f"
    ab_dim1 = *ldab;
#line 226 "sla_gbamv.f"
    ab_offset = 1 + ab_dim1;
#line 226 "sla_gbamv.f"
    ab -= ab_offset;
#line 226 "sla_gbamv.f"
    --x;
#line 226 "sla_gbamv.f"
    --y;
#line 226 "sla_gbamv.f"

#line 226 "sla_gbamv.f"
    /* Function Body */
#line 226 "sla_gbamv.f"
    info = 0;
#line 227 "sla_gbamv.f"
    if (! (*trans == ilatrans_("N", (ftnlen)1) || *trans == ilatrans_("T", (
	    ftnlen)1) || *trans == ilatrans_("C", (ftnlen)1))) {
#line 230 "sla_gbamv.f"
	info = 1;
#line 231 "sla_gbamv.f"
    } else if (*m < 0) {
#line 232 "sla_gbamv.f"
	info = 2;
#line 233 "sla_gbamv.f"
    } else if (*n < 0) {
#line 234 "sla_gbamv.f"
	info = 3;
#line 235 "sla_gbamv.f"
    } else if (*kl < 0 || *kl > *m - 1) {
#line 236 "sla_gbamv.f"
	info = 4;
#line 237 "sla_gbamv.f"
    } else if (*ku < 0 || *ku > *n - 1) {
#line 238 "sla_gbamv.f"
	info = 5;
#line 239 "sla_gbamv.f"
    } else if (*ldab < *kl + *ku + 1) {
#line 240 "sla_gbamv.f"
	info = 6;
#line 241 "sla_gbamv.f"
    } else if (*incx == 0) {
#line 242 "sla_gbamv.f"
	info = 8;
#line 243 "sla_gbamv.f"
    } else if (*incy == 0) {
#line 244 "sla_gbamv.f"
	info = 11;
#line 245 "sla_gbamv.f"
    }
#line 246 "sla_gbamv.f"
    if (info != 0) {
#line 247 "sla_gbamv.f"
	xerbla_("SLA_GBAMV ", &info, (ftnlen)10);
#line 248 "sla_gbamv.f"
	return 0;
#line 249 "sla_gbamv.f"
    }

/*     Quick return if possible. */

#line 253 "sla_gbamv.f"
    if (*m == 0 || *n == 0 || *alpha == 0. && *beta == 1.) {
#line 253 "sla_gbamv.f"
	return 0;
#line 253 "sla_gbamv.f"
    }

/*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set */
/*     up the start points in  X  and  Y. */

#line 260 "sla_gbamv.f"
    if (*trans == ilatrans_("N", (ftnlen)1)) {
#line 261 "sla_gbamv.f"
	lenx = *n;
#line 262 "sla_gbamv.f"
	leny = *m;
#line 263 "sla_gbamv.f"
    } else {
#line 264 "sla_gbamv.f"
	lenx = *m;
#line 265 "sla_gbamv.f"
	leny = *n;
#line 266 "sla_gbamv.f"
    }
#line 267 "sla_gbamv.f"
    if (*incx > 0) {
#line 268 "sla_gbamv.f"
	kx = 1;
#line 269 "sla_gbamv.f"
    } else {
#line 270 "sla_gbamv.f"
	kx = 1 - (lenx - 1) * *incx;
#line 271 "sla_gbamv.f"
    }
#line 272 "sla_gbamv.f"
    if (*incy > 0) {
#line 273 "sla_gbamv.f"
	ky = 1;
#line 274 "sla_gbamv.f"
    } else {
#line 275 "sla_gbamv.f"
	ky = 1 - (leny - 1) * *incy;
#line 276 "sla_gbamv.f"
    }

/*     Set SAFE1 essentially to be the underflow threshold times the */
/*     number of additions in each row. */

#line 281 "sla_gbamv.f"
    safe1 = slamch_("Safe minimum", (ftnlen)12);
#line 282 "sla_gbamv.f"
    safe1 = (*n + 1) * safe1;

/*     Form  y := alpha*abs(A)*abs(x) + beta*abs(y). */

/*     The O(M*N) SYMB_ZERO tests could be replaced by O(N) queries to */
/*     the inexact flag.  Still doesn't help change the iteration order */
/*     to per-column. */

#line 290 "sla_gbamv.f"
    kd = *ku + 1;
#line 291 "sla_gbamv.f"
    ke = *kl + 1;
#line 292 "sla_gbamv.f"
    iy = ky;
#line 293 "sla_gbamv.f"
    if (*incx == 1) {
#line 294 "sla_gbamv.f"
	if (*trans == ilatrans_("N", (ftnlen)1)) {
#line 295 "sla_gbamv.f"
	    i__1 = leny;
#line 295 "sla_gbamv.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 296 "sla_gbamv.f"
		if (*beta == 0.) {
#line 297 "sla_gbamv.f"
		    symb_zero__ = TRUE_;
#line 298 "sla_gbamv.f"
		    y[iy] = 0.;
#line 299 "sla_gbamv.f"
		} else if (y[iy] == 0.) {
#line 300 "sla_gbamv.f"
		    symb_zero__ = TRUE_;
#line 301 "sla_gbamv.f"
		} else {
#line 302 "sla_gbamv.f"
		    symb_zero__ = FALSE_;
#line 303 "sla_gbamv.f"
		    y[iy] = *beta * (d__1 = y[iy], abs(d__1));
#line 304 "sla_gbamv.f"
		}
#line 305 "sla_gbamv.f"
		if (*alpha != 0.) {
/* Computing MAX */
#line 306 "sla_gbamv.f"
		    i__2 = i__ - *kl;
/* Computing MIN */
#line 306 "sla_gbamv.f"
		    i__4 = i__ + *ku;
#line 306 "sla_gbamv.f"
		    i__3 = min(i__4,lenx);
#line 306 "sla_gbamv.f"
		    for (j = max(i__2,1); j <= i__3; ++j) {
#line 307 "sla_gbamv.f"
			temp = (d__1 = ab[kd + i__ - j + j * ab_dim1], abs(
				d__1));
#line 308 "sla_gbamv.f"
			symb_zero__ = symb_zero__ && (x[j] == 0. || temp == 
				0.);
#line 311 "sla_gbamv.f"
			y[iy] += *alpha * (d__1 = x[j], abs(d__1)) * temp;
#line 312 "sla_gbamv.f"
		    }
#line 313 "sla_gbamv.f"
		}
#line 315 "sla_gbamv.f"
		if (! symb_zero__) {
#line 315 "sla_gbamv.f"
		    y[iy] += d_sign(&safe1, &y[iy]);
#line 315 "sla_gbamv.f"
		}
#line 317 "sla_gbamv.f"
		iy += *incy;
#line 318 "sla_gbamv.f"
	    }
#line 319 "sla_gbamv.f"
	} else {
#line 320 "sla_gbamv.f"
	    i__1 = leny;
#line 320 "sla_gbamv.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 321 "sla_gbamv.f"
		if (*beta == 0.) {
#line 322 "sla_gbamv.f"
		    symb_zero__ = TRUE_;
#line 323 "sla_gbamv.f"
		    y[iy] = 0.;
#line 324 "sla_gbamv.f"
		} else if (y[iy] == 0.) {
#line 325 "sla_gbamv.f"
		    symb_zero__ = TRUE_;
#line 326 "sla_gbamv.f"
		} else {
#line 327 "sla_gbamv.f"
		    symb_zero__ = FALSE_;
#line 328 "sla_gbamv.f"
		    y[iy] = *beta * (d__1 = y[iy], abs(d__1));
#line 329 "sla_gbamv.f"
		}
#line 330 "sla_gbamv.f"
		if (*alpha != 0.) {
/* Computing MAX */
#line 331 "sla_gbamv.f"
		    i__3 = i__ - *kl;
/* Computing MIN */
#line 331 "sla_gbamv.f"
		    i__4 = i__ + *ku;
#line 331 "sla_gbamv.f"
		    i__2 = min(i__4,lenx);
#line 331 "sla_gbamv.f"
		    for (j = max(i__3,1); j <= i__2; ++j) {
#line 332 "sla_gbamv.f"
			temp = (d__1 = ab[ke - i__ + j + i__ * ab_dim1], abs(
				d__1));
#line 333 "sla_gbamv.f"
			symb_zero__ = symb_zero__ && (x[j] == 0. || temp == 
				0.);
#line 336 "sla_gbamv.f"
			y[iy] += *alpha * (d__1 = x[j], abs(d__1)) * temp;
#line 337 "sla_gbamv.f"
		    }
#line 338 "sla_gbamv.f"
		}
#line 340 "sla_gbamv.f"
		if (! symb_zero__) {
#line 340 "sla_gbamv.f"
		    y[iy] += d_sign(&safe1, &y[iy]);
#line 340 "sla_gbamv.f"
		}
#line 342 "sla_gbamv.f"
		iy += *incy;
#line 343 "sla_gbamv.f"
	    }
#line 344 "sla_gbamv.f"
	}
#line 345 "sla_gbamv.f"
    } else {
#line 346 "sla_gbamv.f"
	if (*trans == ilatrans_("N", (ftnlen)1)) {
#line 347 "sla_gbamv.f"
	    i__1 = leny;
#line 347 "sla_gbamv.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 348 "sla_gbamv.f"
		if (*beta == 0.) {
#line 349 "sla_gbamv.f"
		    symb_zero__ = TRUE_;
#line 350 "sla_gbamv.f"
		    y[iy] = 0.;
#line 351 "sla_gbamv.f"
		} else if (y[iy] == 0.) {
#line 352 "sla_gbamv.f"
		    symb_zero__ = TRUE_;
#line 353 "sla_gbamv.f"
		} else {
#line 354 "sla_gbamv.f"
		    symb_zero__ = FALSE_;
#line 355 "sla_gbamv.f"
		    y[iy] = *beta * (d__1 = y[iy], abs(d__1));
#line 356 "sla_gbamv.f"
		}
#line 357 "sla_gbamv.f"
		if (*alpha != 0.) {
#line 358 "sla_gbamv.f"
		    jx = kx;
/* Computing MAX */
#line 359 "sla_gbamv.f"
		    i__2 = i__ - *kl;
/* Computing MIN */
#line 359 "sla_gbamv.f"
		    i__4 = i__ + *ku;
#line 359 "sla_gbamv.f"
		    i__3 = min(i__4,lenx);
#line 359 "sla_gbamv.f"
		    for (j = max(i__2,1); j <= i__3; ++j) {
#line 360 "sla_gbamv.f"
			temp = (d__1 = ab[kd + i__ - j + j * ab_dim1], abs(
				d__1));
#line 361 "sla_gbamv.f"
			symb_zero__ = symb_zero__ && (x[jx] == 0. || temp == 
				0.);
#line 364 "sla_gbamv.f"
			y[iy] += *alpha * (d__1 = x[jx], abs(d__1)) * temp;
#line 365 "sla_gbamv.f"
			jx += *incx;
#line 366 "sla_gbamv.f"
		    }
#line 367 "sla_gbamv.f"
		}
#line 369 "sla_gbamv.f"
		if (! symb_zero__) {
#line 369 "sla_gbamv.f"
		    y[iy] += d_sign(&safe1, &y[iy]);
#line 369 "sla_gbamv.f"
		}
#line 372 "sla_gbamv.f"
		iy += *incy;
#line 373 "sla_gbamv.f"
	    }
#line 374 "sla_gbamv.f"
	} else {
#line 375 "sla_gbamv.f"
	    i__1 = leny;
#line 375 "sla_gbamv.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 376 "sla_gbamv.f"
		if (*beta == 0.) {
#line 377 "sla_gbamv.f"
		    symb_zero__ = TRUE_;
#line 378 "sla_gbamv.f"
		    y[iy] = 0.;
#line 379 "sla_gbamv.f"
		} else if (y[iy] == 0.) {
#line 380 "sla_gbamv.f"
		    symb_zero__ = TRUE_;
#line 381 "sla_gbamv.f"
		} else {
#line 382 "sla_gbamv.f"
		    symb_zero__ = FALSE_;
#line 383 "sla_gbamv.f"
		    y[iy] = *beta * (d__1 = y[iy], abs(d__1));
#line 384 "sla_gbamv.f"
		}
#line 385 "sla_gbamv.f"
		if (*alpha != 0.) {
#line 386 "sla_gbamv.f"
		    jx = kx;
/* Computing MAX */
#line 387 "sla_gbamv.f"
		    i__3 = i__ - *kl;
/* Computing MIN */
#line 387 "sla_gbamv.f"
		    i__4 = i__ + *ku;
#line 387 "sla_gbamv.f"
		    i__2 = min(i__4,lenx);
#line 387 "sla_gbamv.f"
		    for (j = max(i__3,1); j <= i__2; ++j) {
#line 388 "sla_gbamv.f"
			temp = (d__1 = ab[ke - i__ + j + i__ * ab_dim1], abs(
				d__1));
#line 389 "sla_gbamv.f"
			symb_zero__ = symb_zero__ && (x[jx] == 0. || temp == 
				0.);
#line 392 "sla_gbamv.f"
			y[iy] += *alpha * (d__1 = x[jx], abs(d__1)) * temp;
#line 393 "sla_gbamv.f"
			jx += *incx;
#line 394 "sla_gbamv.f"
		    }
#line 395 "sla_gbamv.f"
		}
#line 397 "sla_gbamv.f"
		if (! symb_zero__) {
#line 397 "sla_gbamv.f"
		    y[iy] += d_sign(&safe1, &y[iy]);
#line 397 "sla_gbamv.f"
		}
#line 400 "sla_gbamv.f"
		iy += *incy;
#line 401 "sla_gbamv.f"
	    }
#line 402 "sla_gbamv.f"
	}
#line 404 "sla_gbamv.f"
    }

#line 406 "sla_gbamv.f"
    return 0;

/*     End of SLA_GBAMV */

} /* sla_gbamv__ */

