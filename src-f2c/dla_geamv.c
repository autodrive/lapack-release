#line 1 "dla_geamv.f"
/* dla_geamv.f -- translated by f2c (version 20100827).
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

#line 1 "dla_geamv.f"
/* > \brief \b DLA_GEAMV computes a matrix-vector product using a general matrix to calculate error bounds. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLA_GEAMV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dla_gea
mv.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dla_gea
mv.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dla_gea
mv.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLA_GEAMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, */
/*                              Y, INCY ) */

/*       .. Scalar Arguments .. */
/*       DOUBLE PRECISION   ALPHA, BETA */
/*       INTEGER            INCX, INCY, LDA, M, N, TRANS */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   A( LDA, * ), X( * ), Y( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLA_GEAMV  performs one of the matrix-vector operations */
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
/* >          A is DOUBLE PRECISION array of DIMENSION ( LDA, n ) */
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
/* >          Y is DOUBLE PRECISION */
/* >           Array of DIMENSION at least */
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

/* > \ingroup doubleGEcomputational */

/*  ===================================================================== */
/* Subroutine */ int dla_geamv__(integer *trans, integer *m, integer *n, 
	doublereal *alpha, doublereal *a, integer *lda, doublereal *x, 
	integer *incx, doublereal *beta, doublereal *y, integer *incy)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double d_sign(doublereal *, doublereal *);

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
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 216 "dla_geamv.f"
    /* Parameter adjustments */
#line 216 "dla_geamv.f"
    a_dim1 = *lda;
#line 216 "dla_geamv.f"
    a_offset = 1 + a_dim1;
#line 216 "dla_geamv.f"
    a -= a_offset;
#line 216 "dla_geamv.f"
    --x;
#line 216 "dla_geamv.f"
    --y;
#line 216 "dla_geamv.f"

#line 216 "dla_geamv.f"
    /* Function Body */
#line 216 "dla_geamv.f"
    info = 0;
#line 217 "dla_geamv.f"
    if (! (*trans == ilatrans_("N", (ftnlen)1) || *trans == ilatrans_("T", (
	    ftnlen)1) || *trans == ilatrans_("C", (ftnlen)1))) {
#line 220 "dla_geamv.f"
	info = 1;
#line 221 "dla_geamv.f"
    } else if (*m < 0) {
#line 222 "dla_geamv.f"
	info = 2;
#line 223 "dla_geamv.f"
    } else if (*n < 0) {
#line 224 "dla_geamv.f"
	info = 3;
#line 225 "dla_geamv.f"
    } else if (*lda < max(1,*m)) {
#line 226 "dla_geamv.f"
	info = 6;
#line 227 "dla_geamv.f"
    } else if (*incx == 0) {
#line 228 "dla_geamv.f"
	info = 8;
#line 229 "dla_geamv.f"
    } else if (*incy == 0) {
#line 230 "dla_geamv.f"
	info = 11;
#line 231 "dla_geamv.f"
    }
#line 232 "dla_geamv.f"
    if (info != 0) {
#line 233 "dla_geamv.f"
	xerbla_("DLA_GEAMV ", &info, (ftnlen)10);
#line 234 "dla_geamv.f"
	return 0;
#line 235 "dla_geamv.f"
    }

/*     Quick return if possible. */

#line 239 "dla_geamv.f"
    if (*m == 0 || *n == 0 || *alpha == 0. && *beta == 1.) {
#line 239 "dla_geamv.f"
	return 0;
#line 239 "dla_geamv.f"
    }

/*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set */
/*     up the start points in  X  and  Y. */

#line 246 "dla_geamv.f"
    if (*trans == ilatrans_("N", (ftnlen)1)) {
#line 247 "dla_geamv.f"
	lenx = *n;
#line 248 "dla_geamv.f"
	leny = *m;
#line 249 "dla_geamv.f"
    } else {
#line 250 "dla_geamv.f"
	lenx = *m;
#line 251 "dla_geamv.f"
	leny = *n;
#line 252 "dla_geamv.f"
    }
#line 253 "dla_geamv.f"
    if (*incx > 0) {
#line 254 "dla_geamv.f"
	kx = 1;
#line 255 "dla_geamv.f"
    } else {
#line 256 "dla_geamv.f"
	kx = 1 - (lenx - 1) * *incx;
#line 257 "dla_geamv.f"
    }
#line 258 "dla_geamv.f"
    if (*incy > 0) {
#line 259 "dla_geamv.f"
	ky = 1;
#line 260 "dla_geamv.f"
    } else {
#line 261 "dla_geamv.f"
	ky = 1 - (leny - 1) * *incy;
#line 262 "dla_geamv.f"
    }

/*     Set SAFE1 essentially to be the underflow threshold times the */
/*     number of additions in each row. */

#line 267 "dla_geamv.f"
    safe1 = dlamch_("Safe minimum", (ftnlen)12);
#line 268 "dla_geamv.f"
    safe1 = (*n + 1) * safe1;

/*     Form  y := alpha*abs(A)*abs(x) + beta*abs(y). */

/*     The O(M*N) SYMB_ZERO tests could be replaced by O(N) queries to */
/*     the inexact flag.  Still doesn't help change the iteration order */
/*     to per-column. */

#line 276 "dla_geamv.f"
    iy = ky;
#line 277 "dla_geamv.f"
    if (*incx == 1) {
#line 278 "dla_geamv.f"
	if (*trans == ilatrans_("N", (ftnlen)1)) {
#line 279 "dla_geamv.f"
	    i__1 = leny;
#line 279 "dla_geamv.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 280 "dla_geamv.f"
		if (*beta == 0.) {
#line 281 "dla_geamv.f"
		    symb_zero__ = TRUE_;
#line 282 "dla_geamv.f"
		    y[iy] = 0.;
#line 283 "dla_geamv.f"
		} else if (y[iy] == 0.) {
#line 284 "dla_geamv.f"
		    symb_zero__ = TRUE_;
#line 285 "dla_geamv.f"
		} else {
#line 286 "dla_geamv.f"
		    symb_zero__ = FALSE_;
#line 287 "dla_geamv.f"
		    y[iy] = *beta * (d__1 = y[iy], abs(d__1));
#line 288 "dla_geamv.f"
		}
#line 289 "dla_geamv.f"
		if (*alpha != 0.) {
#line 290 "dla_geamv.f"
		    i__2 = lenx;
#line 290 "dla_geamv.f"
		    for (j = 1; j <= i__2; ++j) {
#line 291 "dla_geamv.f"
			temp = (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 292 "dla_geamv.f"
			symb_zero__ = symb_zero__ && (x[j] == 0. || temp == 
				0.);
#line 295 "dla_geamv.f"
			y[iy] += *alpha * (d__1 = x[j], abs(d__1)) * temp;
#line 296 "dla_geamv.f"
		    }
#line 297 "dla_geamv.f"
		}
#line 299 "dla_geamv.f"
		if (! symb_zero__) {
#line 299 "dla_geamv.f"
		    y[iy] += d_sign(&safe1, &y[iy]);
#line 299 "dla_geamv.f"
		}
#line 302 "dla_geamv.f"
		iy += *incy;
#line 303 "dla_geamv.f"
	    }
#line 304 "dla_geamv.f"
	} else {
#line 305 "dla_geamv.f"
	    i__1 = leny;
#line 305 "dla_geamv.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 306 "dla_geamv.f"
		if (*beta == 0.) {
#line 307 "dla_geamv.f"
		    symb_zero__ = TRUE_;
#line 308 "dla_geamv.f"
		    y[iy] = 0.;
#line 309 "dla_geamv.f"
		} else if (y[iy] == 0.) {
#line 310 "dla_geamv.f"
		    symb_zero__ = TRUE_;
#line 311 "dla_geamv.f"
		} else {
#line 312 "dla_geamv.f"
		    symb_zero__ = FALSE_;
#line 313 "dla_geamv.f"
		    y[iy] = *beta * (d__1 = y[iy], abs(d__1));
#line 314 "dla_geamv.f"
		}
#line 315 "dla_geamv.f"
		if (*alpha != 0.) {
#line 316 "dla_geamv.f"
		    i__2 = lenx;
#line 316 "dla_geamv.f"
		    for (j = 1; j <= i__2; ++j) {
#line 317 "dla_geamv.f"
			temp = (d__1 = a[j + i__ * a_dim1], abs(d__1));
#line 318 "dla_geamv.f"
			symb_zero__ = symb_zero__ && (x[j] == 0. || temp == 
				0.);
#line 321 "dla_geamv.f"
			y[iy] += *alpha * (d__1 = x[j], abs(d__1)) * temp;
#line 322 "dla_geamv.f"
		    }
#line 323 "dla_geamv.f"
		}
#line 325 "dla_geamv.f"
		if (! symb_zero__) {
#line 325 "dla_geamv.f"
		    y[iy] += d_sign(&safe1, &y[iy]);
#line 325 "dla_geamv.f"
		}
#line 328 "dla_geamv.f"
		iy += *incy;
#line 329 "dla_geamv.f"
	    }
#line 330 "dla_geamv.f"
	}
#line 331 "dla_geamv.f"
    } else {
#line 332 "dla_geamv.f"
	if (*trans == ilatrans_("N", (ftnlen)1)) {
#line 333 "dla_geamv.f"
	    i__1 = leny;
#line 333 "dla_geamv.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 334 "dla_geamv.f"
		if (*beta == 0.) {
#line 335 "dla_geamv.f"
		    symb_zero__ = TRUE_;
#line 336 "dla_geamv.f"
		    y[iy] = 0.;
#line 337 "dla_geamv.f"
		} else if (y[iy] == 0.) {
#line 338 "dla_geamv.f"
		    symb_zero__ = TRUE_;
#line 339 "dla_geamv.f"
		} else {
#line 340 "dla_geamv.f"
		    symb_zero__ = FALSE_;
#line 341 "dla_geamv.f"
		    y[iy] = *beta * (d__1 = y[iy], abs(d__1));
#line 342 "dla_geamv.f"
		}
#line 343 "dla_geamv.f"
		if (*alpha != 0.) {
#line 344 "dla_geamv.f"
		    jx = kx;
#line 345 "dla_geamv.f"
		    i__2 = lenx;
#line 345 "dla_geamv.f"
		    for (j = 1; j <= i__2; ++j) {
#line 346 "dla_geamv.f"
			temp = (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 347 "dla_geamv.f"
			symb_zero__ = symb_zero__ && (x[jx] == 0. || temp == 
				0.);
#line 350 "dla_geamv.f"
			y[iy] += *alpha * (d__1 = x[jx], abs(d__1)) * temp;
#line 351 "dla_geamv.f"
			jx += *incx;
#line 352 "dla_geamv.f"
		    }
#line 353 "dla_geamv.f"
		}
#line 355 "dla_geamv.f"
		if (! symb_zero__) {
#line 355 "dla_geamv.f"
		    y[iy] += d_sign(&safe1, &y[iy]);
#line 355 "dla_geamv.f"
		}
#line 358 "dla_geamv.f"
		iy += *incy;
#line 359 "dla_geamv.f"
	    }
#line 360 "dla_geamv.f"
	} else {
#line 361 "dla_geamv.f"
	    i__1 = leny;
#line 361 "dla_geamv.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 362 "dla_geamv.f"
		if (*beta == 0.) {
#line 363 "dla_geamv.f"
		    symb_zero__ = TRUE_;
#line 364 "dla_geamv.f"
		    y[iy] = 0.;
#line 365 "dla_geamv.f"
		} else if (y[iy] == 0.) {
#line 366 "dla_geamv.f"
		    symb_zero__ = TRUE_;
#line 367 "dla_geamv.f"
		} else {
#line 368 "dla_geamv.f"
		    symb_zero__ = FALSE_;
#line 369 "dla_geamv.f"
		    y[iy] = *beta * (d__1 = y[iy], abs(d__1));
#line 370 "dla_geamv.f"
		}
#line 371 "dla_geamv.f"
		if (*alpha != 0.) {
#line 372 "dla_geamv.f"
		    jx = kx;
#line 373 "dla_geamv.f"
		    i__2 = lenx;
#line 373 "dla_geamv.f"
		    for (j = 1; j <= i__2; ++j) {
#line 374 "dla_geamv.f"
			temp = (d__1 = a[j + i__ * a_dim1], abs(d__1));
#line 375 "dla_geamv.f"
			symb_zero__ = symb_zero__ && (x[jx] == 0. || temp == 
				0.);
#line 378 "dla_geamv.f"
			y[iy] += *alpha * (d__1 = x[jx], abs(d__1)) * temp;
#line 379 "dla_geamv.f"
			jx += *incx;
#line 380 "dla_geamv.f"
		    }
#line 381 "dla_geamv.f"
		}
#line 383 "dla_geamv.f"
		if (! symb_zero__) {
#line 383 "dla_geamv.f"
		    y[iy] += d_sign(&safe1, &y[iy]);
#line 383 "dla_geamv.f"
		}
#line 386 "dla_geamv.f"
		iy += *incy;
#line 387 "dla_geamv.f"
	    }
#line 388 "dla_geamv.f"
	}
#line 390 "dla_geamv.f"
    }

#line 392 "dla_geamv.f"
    return 0;

/*     End of DLA_GEAMV */

} /* dla_geamv__ */

