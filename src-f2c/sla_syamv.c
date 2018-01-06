#line 1 "sla_syamv.f"
/* sla_syamv.f -- translated by f2c (version 20100827).
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

#line 1 "sla_syamv.f"
/* > \brief \b SLA_SYAMV computes a matrix-vector product using a symmetric indefinite matrix to calculate err
or bounds. */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SLA_SYAMV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sla_sya
mv.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sla_sya
mv.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sla_sya
mv.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SLA_SYAMV( UPLO, N, ALPHA, A, LDA, X, INCX, BETA, Y, */
/*                             INCY ) */

/*       .. Scalar Arguments .. */
/*       REAL               ALPHA, BETA */
/*       INTEGER            INCX, INCY, LDA, N, UPLO */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               A( LDA, * ), X( * ), Y( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SLA_SYAMV  performs the matrix-vector operation */
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
/* >          A is REAL array of DIMENSION ( LDA, n ). */
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
/* >          X is REAL array, dimension */
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

/* > \date December 2016 */

/* > \ingroup realSYcomputational */

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
/* Subroutine */ int sla_syamv__(integer *uplo, integer *n, doublereal *alpha,
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
    extern doublereal slamch_(char *, ftnlen);
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
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 219 "sla_syamv.f"
    /* Parameter adjustments */
#line 219 "sla_syamv.f"
    a_dim1 = *lda;
#line 219 "sla_syamv.f"
    a_offset = 1 + a_dim1;
#line 219 "sla_syamv.f"
    a -= a_offset;
#line 219 "sla_syamv.f"
    --x;
#line 219 "sla_syamv.f"
    --y;
#line 219 "sla_syamv.f"

#line 219 "sla_syamv.f"
    /* Function Body */
#line 219 "sla_syamv.f"
    info = 0;
#line 220 "sla_syamv.f"
    if (*uplo != ilauplo_("U", (ftnlen)1) && *uplo != ilauplo_("L", (ftnlen)1)
	    ) {
#line 222 "sla_syamv.f"
	info = 1;
#line 223 "sla_syamv.f"
    } else if (*n < 0) {
#line 224 "sla_syamv.f"
	info = 2;
#line 225 "sla_syamv.f"
    } else if (*lda < max(1,*n)) {
#line 226 "sla_syamv.f"
	info = 5;
#line 227 "sla_syamv.f"
    } else if (*incx == 0) {
#line 228 "sla_syamv.f"
	info = 7;
#line 229 "sla_syamv.f"
    } else if (*incy == 0) {
#line 230 "sla_syamv.f"
	info = 10;
#line 231 "sla_syamv.f"
    }
#line 232 "sla_syamv.f"
    if (info != 0) {
#line 233 "sla_syamv.f"
	xerbla_("SSYMV ", &info, (ftnlen)6);
#line 234 "sla_syamv.f"
	return 0;
#line 235 "sla_syamv.f"
    }

/*     Quick return if possible. */

#line 239 "sla_syamv.f"
    if (*n == 0 || *alpha == 0. && *beta == 1.) {
#line 239 "sla_syamv.f"
	return 0;
#line 239 "sla_syamv.f"
    }

/*     Set up the start points in  X  and  Y. */

#line 244 "sla_syamv.f"
    if (*incx > 0) {
#line 245 "sla_syamv.f"
	kx = 1;
#line 246 "sla_syamv.f"
    } else {
#line 247 "sla_syamv.f"
	kx = 1 - (*n - 1) * *incx;
#line 248 "sla_syamv.f"
    }
#line 249 "sla_syamv.f"
    if (*incy > 0) {
#line 250 "sla_syamv.f"
	ky = 1;
#line 251 "sla_syamv.f"
    } else {
#line 252 "sla_syamv.f"
	ky = 1 - (*n - 1) * *incy;
#line 253 "sla_syamv.f"
    }

/*     Set SAFE1 essentially to be the underflow threshold times the */
/*     number of additions in each row. */

#line 258 "sla_syamv.f"
    safe1 = slamch_("Safe minimum", (ftnlen)12);
#line 259 "sla_syamv.f"
    safe1 = (*n + 1) * safe1;

/*     Form  y := alpha*abs(A)*abs(x) + beta*abs(y). */

/*     The O(N^2) SYMB_ZERO tests could be replaced by O(N) queries to */
/*     the inexact flag.  Still doesn't help change the iteration order */
/*     to per-column. */

#line 267 "sla_syamv.f"
    iy = ky;
#line 268 "sla_syamv.f"
    if (*incx == 1) {
#line 269 "sla_syamv.f"
	if (*uplo == ilauplo_("U", (ftnlen)1)) {
#line 270 "sla_syamv.f"
	    i__1 = *n;
#line 270 "sla_syamv.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 271 "sla_syamv.f"
		if (*beta == 0.) {
#line 272 "sla_syamv.f"
		    symb_zero__ = TRUE_;
#line 273 "sla_syamv.f"
		    y[iy] = 0.;
#line 274 "sla_syamv.f"
		} else if (y[iy] == 0.) {
#line 275 "sla_syamv.f"
		    symb_zero__ = TRUE_;
#line 276 "sla_syamv.f"
		} else {
#line 277 "sla_syamv.f"
		    symb_zero__ = FALSE_;
#line 278 "sla_syamv.f"
		    y[iy] = *beta * (d__1 = y[iy], abs(d__1));
#line 279 "sla_syamv.f"
		}
#line 280 "sla_syamv.f"
		if (*alpha != 0.) {
#line 281 "sla_syamv.f"
		    i__2 = i__;
#line 281 "sla_syamv.f"
		    for (j = 1; j <= i__2; ++j) {
#line 282 "sla_syamv.f"
			temp = (d__1 = a[j + i__ * a_dim1], abs(d__1));
#line 283 "sla_syamv.f"
			symb_zero__ = symb_zero__ && (x[j] == 0. || temp == 
				0.);
#line 286 "sla_syamv.f"
			y[iy] += *alpha * (d__1 = x[j], abs(d__1)) * temp;
#line 287 "sla_syamv.f"
		    }
#line 288 "sla_syamv.f"
		    i__2 = *n;
#line 288 "sla_syamv.f"
		    for (j = i__ + 1; j <= i__2; ++j) {
#line 289 "sla_syamv.f"
			temp = (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 290 "sla_syamv.f"
			symb_zero__ = symb_zero__ && (x[j] == 0. || temp == 
				0.);
#line 293 "sla_syamv.f"
			y[iy] += *alpha * (d__1 = x[j], abs(d__1)) * temp;
#line 294 "sla_syamv.f"
		    }
#line 295 "sla_syamv.f"
		}
#line 297 "sla_syamv.f"
		if (! symb_zero__) {
#line 297 "sla_syamv.f"
		    y[iy] += d_sign(&safe1, &y[iy]);
#line 297 "sla_syamv.f"
		}
#line 300 "sla_syamv.f"
		iy += *incy;
#line 301 "sla_syamv.f"
	    }
#line 302 "sla_syamv.f"
	} else {
#line 303 "sla_syamv.f"
	    i__1 = *n;
#line 303 "sla_syamv.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 304 "sla_syamv.f"
		if (*beta == 0.) {
#line 305 "sla_syamv.f"
		    symb_zero__ = TRUE_;
#line 306 "sla_syamv.f"
		    y[iy] = 0.;
#line 307 "sla_syamv.f"
		} else if (y[iy] == 0.) {
#line 308 "sla_syamv.f"
		    symb_zero__ = TRUE_;
#line 309 "sla_syamv.f"
		} else {
#line 310 "sla_syamv.f"
		    symb_zero__ = FALSE_;
#line 311 "sla_syamv.f"
		    y[iy] = *beta * (d__1 = y[iy], abs(d__1));
#line 312 "sla_syamv.f"
		}
#line 313 "sla_syamv.f"
		if (*alpha != 0.) {
#line 314 "sla_syamv.f"
		    i__2 = i__;
#line 314 "sla_syamv.f"
		    for (j = 1; j <= i__2; ++j) {
#line 315 "sla_syamv.f"
			temp = (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 316 "sla_syamv.f"
			symb_zero__ = symb_zero__ && (x[j] == 0. || temp == 
				0.);
#line 319 "sla_syamv.f"
			y[iy] += *alpha * (d__1 = x[j], abs(d__1)) * temp;
#line 320 "sla_syamv.f"
		    }
#line 321 "sla_syamv.f"
		    i__2 = *n;
#line 321 "sla_syamv.f"
		    for (j = i__ + 1; j <= i__2; ++j) {
#line 322 "sla_syamv.f"
			temp = (d__1 = a[j + i__ * a_dim1], abs(d__1));
#line 323 "sla_syamv.f"
			symb_zero__ = symb_zero__ && (x[j] == 0. || temp == 
				0.);
#line 326 "sla_syamv.f"
			y[iy] += *alpha * (d__1 = x[j], abs(d__1)) * temp;
#line 327 "sla_syamv.f"
		    }
#line 328 "sla_syamv.f"
		}
#line 330 "sla_syamv.f"
		if (! symb_zero__) {
#line 330 "sla_syamv.f"
		    y[iy] += d_sign(&safe1, &y[iy]);
#line 330 "sla_syamv.f"
		}
#line 333 "sla_syamv.f"
		iy += *incy;
#line 334 "sla_syamv.f"
	    }
#line 335 "sla_syamv.f"
	}
#line 336 "sla_syamv.f"
    } else {
#line 337 "sla_syamv.f"
	if (*uplo == ilauplo_("U", (ftnlen)1)) {
#line 338 "sla_syamv.f"
	    i__1 = *n;
#line 338 "sla_syamv.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 339 "sla_syamv.f"
		if (*beta == 0.) {
#line 340 "sla_syamv.f"
		    symb_zero__ = TRUE_;
#line 341 "sla_syamv.f"
		    y[iy] = 0.;
#line 342 "sla_syamv.f"
		} else if (y[iy] == 0.) {
#line 343 "sla_syamv.f"
		    symb_zero__ = TRUE_;
#line 344 "sla_syamv.f"
		} else {
#line 345 "sla_syamv.f"
		    symb_zero__ = FALSE_;
#line 346 "sla_syamv.f"
		    y[iy] = *beta * (d__1 = y[iy], abs(d__1));
#line 347 "sla_syamv.f"
		}
#line 348 "sla_syamv.f"
		jx = kx;
#line 349 "sla_syamv.f"
		if (*alpha != 0.) {
#line 350 "sla_syamv.f"
		    i__2 = i__;
#line 350 "sla_syamv.f"
		    for (j = 1; j <= i__2; ++j) {
#line 351 "sla_syamv.f"
			temp = (d__1 = a[j + i__ * a_dim1], abs(d__1));
#line 352 "sla_syamv.f"
			symb_zero__ = symb_zero__ && (x[j] == 0. || temp == 
				0.);
#line 355 "sla_syamv.f"
			y[iy] += *alpha * (d__1 = x[jx], abs(d__1)) * temp;
#line 356 "sla_syamv.f"
			jx += *incx;
#line 357 "sla_syamv.f"
		    }
#line 358 "sla_syamv.f"
		    i__2 = *n;
#line 358 "sla_syamv.f"
		    for (j = i__ + 1; j <= i__2; ++j) {
#line 359 "sla_syamv.f"
			temp = (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 360 "sla_syamv.f"
			symb_zero__ = symb_zero__ && (x[j] == 0. || temp == 
				0.);
#line 363 "sla_syamv.f"
			y[iy] += *alpha * (d__1 = x[jx], abs(d__1)) * temp;
#line 364 "sla_syamv.f"
			jx += *incx;
#line 365 "sla_syamv.f"
		    }
#line 366 "sla_syamv.f"
		}
#line 368 "sla_syamv.f"
		if (! symb_zero__) {
#line 368 "sla_syamv.f"
		    y[iy] += d_sign(&safe1, &y[iy]);
#line 368 "sla_syamv.f"
		}
#line 371 "sla_syamv.f"
		iy += *incy;
#line 372 "sla_syamv.f"
	    }
#line 373 "sla_syamv.f"
	} else {
#line 374 "sla_syamv.f"
	    i__1 = *n;
#line 374 "sla_syamv.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 375 "sla_syamv.f"
		if (*beta == 0.) {
#line 376 "sla_syamv.f"
		    symb_zero__ = TRUE_;
#line 377 "sla_syamv.f"
		    y[iy] = 0.;
#line 378 "sla_syamv.f"
		} else if (y[iy] == 0.) {
#line 379 "sla_syamv.f"
		    symb_zero__ = TRUE_;
#line 380 "sla_syamv.f"
		} else {
#line 381 "sla_syamv.f"
		    symb_zero__ = FALSE_;
#line 382 "sla_syamv.f"
		    y[iy] = *beta * (d__1 = y[iy], abs(d__1));
#line 383 "sla_syamv.f"
		}
#line 384 "sla_syamv.f"
		jx = kx;
#line 385 "sla_syamv.f"
		if (*alpha != 0.) {
#line 386 "sla_syamv.f"
		    i__2 = i__;
#line 386 "sla_syamv.f"
		    for (j = 1; j <= i__2; ++j) {
#line 387 "sla_syamv.f"
			temp = (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 388 "sla_syamv.f"
			symb_zero__ = symb_zero__ && (x[j] == 0. || temp == 
				0.);
#line 391 "sla_syamv.f"
			y[iy] += *alpha * (d__1 = x[jx], abs(d__1)) * temp;
#line 392 "sla_syamv.f"
			jx += *incx;
#line 393 "sla_syamv.f"
		    }
#line 394 "sla_syamv.f"
		    i__2 = *n;
#line 394 "sla_syamv.f"
		    for (j = i__ + 1; j <= i__2; ++j) {
#line 395 "sla_syamv.f"
			temp = (d__1 = a[j + i__ * a_dim1], abs(d__1));
#line 396 "sla_syamv.f"
			symb_zero__ = symb_zero__ && (x[j] == 0. || temp == 
				0.);
#line 399 "sla_syamv.f"
			y[iy] += *alpha * (d__1 = x[jx], abs(d__1)) * temp;
#line 400 "sla_syamv.f"
			jx += *incx;
#line 401 "sla_syamv.f"
		    }
#line 402 "sla_syamv.f"
		}
#line 404 "sla_syamv.f"
		if (! symb_zero__) {
#line 404 "sla_syamv.f"
		    y[iy] += d_sign(&safe1, &y[iy]);
#line 404 "sla_syamv.f"
		}
#line 407 "sla_syamv.f"
		iy += *incy;
#line 408 "sla_syamv.f"
	    }
#line 409 "sla_syamv.f"
	}
#line 411 "sla_syamv.f"
    }

#line 413 "sla_syamv.f"
    return 0;

/*     End of SLA_SYAMV */

} /* sla_syamv__ */

