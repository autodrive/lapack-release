#line 1 "dla_syamv.f"
/* dla_syamv.f -- translated by f2c (version 20100827).
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

#line 1 "dla_syamv.f"
/* > \brief \b DLA_SYAMV computes a matrix-vector product using a symmetric indefinite matrix to calculate err
or bounds. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLA_SYAMV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dla_sya
mv.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dla_sya
mv.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dla_sya
mv.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLA_SYAMV( UPLO, N, ALPHA, A, LDA, X, INCX, BETA, Y, */
/*                             INCY ) */

/*       .. Scalar Arguments .. */
/*       DOUBLE PRECISION   ALPHA, BETA */
/*       INTEGER            INCX, INCY, LDA, N, UPLO */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   A( LDA, * ), X( * ), Y( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLA_SYAMV  performs the matrix-vector operation */
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
/* >          A is DOUBLE PRECISION array, dimension ( LDA, n ). */
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
/* >          X is DOUBLE PRECISION array, dimension */
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

/* > \date June 2017 */

/* > \ingroup doubleSYcomputational */

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
/* Subroutine */ int dla_syamv__(integer *uplo, integer *n, doublereal *alpha,
	 doublereal *a, integer *lda, doublereal *x, integer *incx, 
	doublereal *beta, doublereal *y, integer *incy)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double d_sign(doublereal *, doublereal *);

    /* Local variables */
    static integer i__, j;
    static logical symb_zero__;
    static integer iy, jx, kx, ky, info;
    static doublereal temp, safe1;
    extern doublereal dlamch_(char *, ftnlen);
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
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 219 "dla_syamv.f"
    /* Parameter adjustments */
#line 219 "dla_syamv.f"
    a_dim1 = *lda;
#line 219 "dla_syamv.f"
    a_offset = 1 + a_dim1;
#line 219 "dla_syamv.f"
    a -= a_offset;
#line 219 "dla_syamv.f"
    --x;
#line 219 "dla_syamv.f"
    --y;
#line 219 "dla_syamv.f"

#line 219 "dla_syamv.f"
    /* Function Body */
#line 219 "dla_syamv.f"
    info = 0;
#line 220 "dla_syamv.f"
    if (*uplo != ilauplo_("U", (ftnlen)1) && *uplo != ilauplo_("L", (ftnlen)1)
	    ) {
#line 222 "dla_syamv.f"
	info = 1;
#line 223 "dla_syamv.f"
    } else if (*n < 0) {
#line 224 "dla_syamv.f"
	info = 2;
#line 225 "dla_syamv.f"
    } else if (*lda < max(1,*n)) {
#line 226 "dla_syamv.f"
	info = 5;
#line 227 "dla_syamv.f"
    } else if (*incx == 0) {
#line 228 "dla_syamv.f"
	info = 7;
#line 229 "dla_syamv.f"
    } else if (*incy == 0) {
#line 230 "dla_syamv.f"
	info = 10;
#line 231 "dla_syamv.f"
    }
#line 232 "dla_syamv.f"
    if (info != 0) {
#line 233 "dla_syamv.f"
	xerbla_("DSYMV ", &info, (ftnlen)6);
#line 234 "dla_syamv.f"
	return 0;
#line 235 "dla_syamv.f"
    }

/*     Quick return if possible. */

#line 239 "dla_syamv.f"
    if (*n == 0 || *alpha == 0. && *beta == 1.) {
#line 239 "dla_syamv.f"
	return 0;
#line 239 "dla_syamv.f"
    }

/*     Set up the start points in  X  and  Y. */

#line 244 "dla_syamv.f"
    if (*incx > 0) {
#line 245 "dla_syamv.f"
	kx = 1;
#line 246 "dla_syamv.f"
    } else {
#line 247 "dla_syamv.f"
	kx = 1 - (*n - 1) * *incx;
#line 248 "dla_syamv.f"
    }
#line 249 "dla_syamv.f"
    if (*incy > 0) {
#line 250 "dla_syamv.f"
	ky = 1;
#line 251 "dla_syamv.f"
    } else {
#line 252 "dla_syamv.f"
	ky = 1 - (*n - 1) * *incy;
#line 253 "dla_syamv.f"
    }

/*     Set SAFE1 essentially to be the underflow threshold times the */
/*     number of additions in each row. */

#line 258 "dla_syamv.f"
    safe1 = dlamch_("Safe minimum", (ftnlen)12);
#line 259 "dla_syamv.f"
    safe1 = (*n + 1) * safe1;

/*     Form  y := alpha*abs(A)*abs(x) + beta*abs(y). */

/*     The O(N^2) SYMB_ZERO tests could be replaced by O(N) queries to */
/*     the inexact flag.  Still doesn't help change the iteration order */
/*     to per-column. */

#line 267 "dla_syamv.f"
    iy = ky;
#line 268 "dla_syamv.f"
    if (*incx == 1) {
#line 269 "dla_syamv.f"
	if (*uplo == ilauplo_("U", (ftnlen)1)) {
#line 270 "dla_syamv.f"
	    i__1 = *n;
#line 270 "dla_syamv.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 271 "dla_syamv.f"
		if (*beta == 0.) {
#line 272 "dla_syamv.f"
		    symb_zero__ = TRUE_;
#line 273 "dla_syamv.f"
		    y[iy] = 0.;
#line 274 "dla_syamv.f"
		} else if (y[iy] == 0.) {
#line 275 "dla_syamv.f"
		    symb_zero__ = TRUE_;
#line 276 "dla_syamv.f"
		} else {
#line 277 "dla_syamv.f"
		    symb_zero__ = FALSE_;
#line 278 "dla_syamv.f"
		    y[iy] = *beta * (d__1 = y[iy], abs(d__1));
#line 279 "dla_syamv.f"
		}
#line 280 "dla_syamv.f"
		if (*alpha != 0.) {
#line 281 "dla_syamv.f"
		    i__2 = i__;
#line 281 "dla_syamv.f"
		    for (j = 1; j <= i__2; ++j) {
#line 282 "dla_syamv.f"
			temp = (d__1 = a[j + i__ * a_dim1], abs(d__1));
#line 283 "dla_syamv.f"
			symb_zero__ = symb_zero__ && (x[j] == 0. || temp == 
				0.);
#line 286 "dla_syamv.f"
			y[iy] += *alpha * (d__1 = x[j], abs(d__1)) * temp;
#line 287 "dla_syamv.f"
		    }
#line 288 "dla_syamv.f"
		    i__2 = *n;
#line 288 "dla_syamv.f"
		    for (j = i__ + 1; j <= i__2; ++j) {
#line 289 "dla_syamv.f"
			temp = (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 290 "dla_syamv.f"
			symb_zero__ = symb_zero__ && (x[j] == 0. || temp == 
				0.);
#line 293 "dla_syamv.f"
			y[iy] += *alpha * (d__1 = x[j], abs(d__1)) * temp;
#line 294 "dla_syamv.f"
		    }
#line 295 "dla_syamv.f"
		}
#line 297 "dla_syamv.f"
		if (! symb_zero__) {
#line 297 "dla_syamv.f"
		    y[iy] += d_sign(&safe1, &y[iy]);
#line 297 "dla_syamv.f"
		}
#line 300 "dla_syamv.f"
		iy += *incy;
#line 301 "dla_syamv.f"
	    }
#line 302 "dla_syamv.f"
	} else {
#line 303 "dla_syamv.f"
	    i__1 = *n;
#line 303 "dla_syamv.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 304 "dla_syamv.f"
		if (*beta == 0.) {
#line 305 "dla_syamv.f"
		    symb_zero__ = TRUE_;
#line 306 "dla_syamv.f"
		    y[iy] = 0.;
#line 307 "dla_syamv.f"
		} else if (y[iy] == 0.) {
#line 308 "dla_syamv.f"
		    symb_zero__ = TRUE_;
#line 309 "dla_syamv.f"
		} else {
#line 310 "dla_syamv.f"
		    symb_zero__ = FALSE_;
#line 311 "dla_syamv.f"
		    y[iy] = *beta * (d__1 = y[iy], abs(d__1));
#line 312 "dla_syamv.f"
		}
#line 313 "dla_syamv.f"
		if (*alpha != 0.) {
#line 314 "dla_syamv.f"
		    i__2 = i__;
#line 314 "dla_syamv.f"
		    for (j = 1; j <= i__2; ++j) {
#line 315 "dla_syamv.f"
			temp = (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 316 "dla_syamv.f"
			symb_zero__ = symb_zero__ && (x[j] == 0. || temp == 
				0.);
#line 319 "dla_syamv.f"
			y[iy] += *alpha * (d__1 = x[j], abs(d__1)) * temp;
#line 320 "dla_syamv.f"
		    }
#line 321 "dla_syamv.f"
		    i__2 = *n;
#line 321 "dla_syamv.f"
		    for (j = i__ + 1; j <= i__2; ++j) {
#line 322 "dla_syamv.f"
			temp = (d__1 = a[j + i__ * a_dim1], abs(d__1));
#line 323 "dla_syamv.f"
			symb_zero__ = symb_zero__ && (x[j] == 0. || temp == 
				0.);
#line 326 "dla_syamv.f"
			y[iy] += *alpha * (d__1 = x[j], abs(d__1)) * temp;
#line 327 "dla_syamv.f"
		    }
#line 328 "dla_syamv.f"
		}
#line 330 "dla_syamv.f"
		if (! symb_zero__) {
#line 330 "dla_syamv.f"
		    y[iy] += d_sign(&safe1, &y[iy]);
#line 330 "dla_syamv.f"
		}
#line 333 "dla_syamv.f"
		iy += *incy;
#line 334 "dla_syamv.f"
	    }
#line 335 "dla_syamv.f"
	}
#line 336 "dla_syamv.f"
    } else {
#line 337 "dla_syamv.f"
	if (*uplo == ilauplo_("U", (ftnlen)1)) {
#line 338 "dla_syamv.f"
	    i__1 = *n;
#line 338 "dla_syamv.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 339 "dla_syamv.f"
		if (*beta == 0.) {
#line 340 "dla_syamv.f"
		    symb_zero__ = TRUE_;
#line 341 "dla_syamv.f"
		    y[iy] = 0.;
#line 342 "dla_syamv.f"
		} else if (y[iy] == 0.) {
#line 343 "dla_syamv.f"
		    symb_zero__ = TRUE_;
#line 344 "dla_syamv.f"
		} else {
#line 345 "dla_syamv.f"
		    symb_zero__ = FALSE_;
#line 346 "dla_syamv.f"
		    y[iy] = *beta * (d__1 = y[iy], abs(d__1));
#line 347 "dla_syamv.f"
		}
#line 348 "dla_syamv.f"
		jx = kx;
#line 349 "dla_syamv.f"
		if (*alpha != 0.) {
#line 350 "dla_syamv.f"
		    i__2 = i__;
#line 350 "dla_syamv.f"
		    for (j = 1; j <= i__2; ++j) {
#line 351 "dla_syamv.f"
			temp = (d__1 = a[j + i__ * a_dim1], abs(d__1));
#line 352 "dla_syamv.f"
			symb_zero__ = symb_zero__ && (x[j] == 0. || temp == 
				0.);
#line 355 "dla_syamv.f"
			y[iy] += *alpha * (d__1 = x[jx], abs(d__1)) * temp;
#line 356 "dla_syamv.f"
			jx += *incx;
#line 357 "dla_syamv.f"
		    }
#line 358 "dla_syamv.f"
		    i__2 = *n;
#line 358 "dla_syamv.f"
		    for (j = i__ + 1; j <= i__2; ++j) {
#line 359 "dla_syamv.f"
			temp = (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 360 "dla_syamv.f"
			symb_zero__ = symb_zero__ && (x[j] == 0. || temp == 
				0.);
#line 363 "dla_syamv.f"
			y[iy] += *alpha * (d__1 = x[jx], abs(d__1)) * temp;
#line 364 "dla_syamv.f"
			jx += *incx;
#line 365 "dla_syamv.f"
		    }
#line 366 "dla_syamv.f"
		}
#line 368 "dla_syamv.f"
		if (! symb_zero__) {
#line 368 "dla_syamv.f"
		    y[iy] += d_sign(&safe1, &y[iy]);
#line 368 "dla_syamv.f"
		}
#line 371 "dla_syamv.f"
		iy += *incy;
#line 372 "dla_syamv.f"
	    }
#line 373 "dla_syamv.f"
	} else {
#line 374 "dla_syamv.f"
	    i__1 = *n;
#line 374 "dla_syamv.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 375 "dla_syamv.f"
		if (*beta == 0.) {
#line 376 "dla_syamv.f"
		    symb_zero__ = TRUE_;
#line 377 "dla_syamv.f"
		    y[iy] = 0.;
#line 378 "dla_syamv.f"
		} else if (y[iy] == 0.) {
#line 379 "dla_syamv.f"
		    symb_zero__ = TRUE_;
#line 380 "dla_syamv.f"
		} else {
#line 381 "dla_syamv.f"
		    symb_zero__ = FALSE_;
#line 382 "dla_syamv.f"
		    y[iy] = *beta * (d__1 = y[iy], abs(d__1));
#line 383 "dla_syamv.f"
		}
#line 384 "dla_syamv.f"
		jx = kx;
#line 385 "dla_syamv.f"
		if (*alpha != 0.) {
#line 386 "dla_syamv.f"
		    i__2 = i__;
#line 386 "dla_syamv.f"
		    for (j = 1; j <= i__2; ++j) {
#line 387 "dla_syamv.f"
			temp = (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 388 "dla_syamv.f"
			symb_zero__ = symb_zero__ && (x[j] == 0. || temp == 
				0.);
#line 391 "dla_syamv.f"
			y[iy] += *alpha * (d__1 = x[jx], abs(d__1)) * temp;
#line 392 "dla_syamv.f"
			jx += *incx;
#line 393 "dla_syamv.f"
		    }
#line 394 "dla_syamv.f"
		    i__2 = *n;
#line 394 "dla_syamv.f"
		    for (j = i__ + 1; j <= i__2; ++j) {
#line 395 "dla_syamv.f"
			temp = (d__1 = a[j + i__ * a_dim1], abs(d__1));
#line 396 "dla_syamv.f"
			symb_zero__ = symb_zero__ && (x[j] == 0. || temp == 
				0.);
#line 399 "dla_syamv.f"
			y[iy] += *alpha * (d__1 = x[jx], abs(d__1)) * temp;
#line 400 "dla_syamv.f"
			jx += *incx;
#line 401 "dla_syamv.f"
		    }
#line 402 "dla_syamv.f"
		}
#line 404 "dla_syamv.f"
		if (! symb_zero__) {
#line 404 "dla_syamv.f"
		    y[iy] += d_sign(&safe1, &y[iy]);
#line 404 "dla_syamv.f"
		}
#line 407 "dla_syamv.f"
		iy += *incy;
#line 408 "dla_syamv.f"
	    }
#line 409 "dla_syamv.f"
	}
#line 411 "dla_syamv.f"
    }

#line 413 "dla_syamv.f"
    return 0;

/*     End of DLA_SYAMV */

} /* dla_syamv__ */

