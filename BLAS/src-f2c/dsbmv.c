#line 1 "dsbmv.f"
/* dsbmv.f -- translated by f2c (version 20100827).
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

#line 1 "dsbmv.f"
/* > \brief \b DSBMV */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSBMV(UPLO,N,K,ALPHA,A,LDA,X,INCX,BETA,Y,INCY) */

/*       .. Scalar Arguments .. */
/*       DOUBLE PRECISION ALPHA,BETA */
/*       INTEGER INCX,INCY,K,LDA,N */
/*       CHARACTER UPLO */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION A(LDA,*),X(*),Y(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DSBMV  performs the matrix-vector  operation */
/* > */
/* >    y := alpha*A*x + beta*y, */
/* > */
/* > where alpha and beta are scalars, x and y are n element vectors and */
/* > A is an n by n symmetric band matrix, with k super-diagonals. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >           On entry, UPLO specifies whether the upper or lower */
/* >           triangular part of the band matrix A is being supplied as */
/* >           follows: */
/* > */
/* >              UPLO = 'U' or 'u'   The upper triangular part of A is */
/* >                                  being supplied. */
/* > */
/* >              UPLO = 'L' or 'l'   The lower triangular part of A is */
/* >                                  being supplied. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >           On entry, N specifies the order of the matrix A. */
/* >           N must be at least zero. */
/* > \endverbatim */
/* > */
/* > \param[in] K */
/* > \verbatim */
/* >          K is INTEGER */
/* >           On entry, K specifies the number of super-diagonals of the */
/* >           matrix A. K must satisfy  0 .le. K. */
/* > \endverbatim */
/* > */
/* > \param[in] ALPHA */
/* > \verbatim */
/* >          ALPHA is DOUBLE PRECISION. */
/* >           On entry, ALPHA specifies the scalar alpha. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension ( LDA, N ) */
/* >           Before entry with UPLO = 'U' or 'u', the leading ( k + 1 ) */
/* >           by n part of the array A must contain the upper triangular */
/* >           band part of the symmetric matrix, supplied column by */
/* >           column, with the leading diagonal of the matrix in row */
/* >           ( k + 1 ) of the array, the first super-diagonal starting at */
/* >           position 2 in row k, and so on. The top left k by k triangle */
/* >           of the array A is not referenced. */
/* >           The following program segment will transfer the upper */
/* >           triangular part of a symmetric band matrix from conventional */
/* >           full matrix storage to band storage: */
/* > */
/* >                 DO 20, J = 1, N */
/* >                    M = K + 1 - J */
/* >                    DO 10, I = MAX( 1, J - K ), J */
/* >                       A( M + I, J ) = matrix( I, J ) */
/* >              10    CONTINUE */
/* >              20 CONTINUE */
/* > */
/* >           Before entry with UPLO = 'L' or 'l', the leading ( k + 1 ) */
/* >           by n part of the array A must contain the lower triangular */
/* >           band part of the symmetric matrix, supplied column by */
/* >           column, with the leading diagonal of the matrix in row 1 of */
/* >           the array, the first sub-diagonal starting at position 1 in */
/* >           row 2, and so on. The bottom right k by k triangle of the */
/* >           array A is not referenced. */
/* >           The following program segment will transfer the lower */
/* >           triangular part of a symmetric band matrix from conventional */
/* >           full matrix storage to band storage: */
/* > */
/* >                 DO 20, J = 1, N */
/* >                    M = 1 - J */
/* >                    DO 10, I = J, MIN( N, J + K ) */
/* >                       A( M + I, J ) = matrix( I, J ) */
/* >              10    CONTINUE */
/* >              20 CONTINUE */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >           On entry, LDA specifies the first dimension of A as declared */
/* >           in the calling (sub) program. LDA must be at least */
/* >           ( k + 1 ). */
/* > \endverbatim */
/* > */
/* > \param[in] X */
/* > \verbatim */
/* >          X is DOUBLE PRECISION array, dimension at least */
/* >           ( 1 + ( n - 1 )*abs( INCX ) ). */
/* >           Before entry, the incremented array X must contain the */
/* >           vector x. */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* >          INCX is INTEGER */
/* >           On entry, INCX specifies the increment for the elements of */
/* >           X. INCX must not be zero. */
/* > \endverbatim */
/* > */
/* > \param[in] BETA */
/* > \verbatim */
/* >          BETA is DOUBLE PRECISION. */
/* >           On entry, BETA specifies the scalar beta. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Y */
/* > \verbatim */
/* >          Y is DOUBLE PRECISION array, dimension at least */
/* >           ( 1 + ( n - 1 )*abs( INCY ) ). */
/* >           Before entry, the incremented array Y must contain the */
/* >           vector y. On exit, Y is overwritten by the updated vector y. */
/* > \endverbatim */
/* > */
/* > \param[in] INCY */
/* > \verbatim */
/* >          INCY is INTEGER */
/* >           On entry, INCY specifies the increment for the elements of */
/* >           Y. INCY must not be zero. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup double_blas_level2 */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  Level 2 Blas routine. */
/* >  The vector and matrix arguments are not referenced when N = 0, or M = 0 */
/* > */
/* >  -- Written on 22-October-1986. */
/* >     Jack Dongarra, Argonne National Lab. */
/* >     Jeremy Du Croz, Nag Central Office. */
/* >     Sven Hammarling, Nag Central Office. */
/* >     Richard Hanson, Sandia National Labs. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dsbmv_(char *uplo, integer *n, integer *k, doublereal *
	alpha, doublereal *a, integer *lda, doublereal *x, integer *incx, 
	doublereal *beta, doublereal *y, integer *incy, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer i__, j, l, ix, iy, jx, jy, kx, ky, info;
    static doublereal temp1, temp2;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer kplus1;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);


/*  -- Reference BLAS level2 routine (version 3.7.0) -- */
/*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    -- */
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
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */

/*     Test the input parameters. */

#line 224 "dsbmv.f"
    /* Parameter adjustments */
#line 224 "dsbmv.f"
    a_dim1 = *lda;
#line 224 "dsbmv.f"
    a_offset = 1 + a_dim1;
#line 224 "dsbmv.f"
    a -= a_offset;
#line 224 "dsbmv.f"
    --x;
#line 224 "dsbmv.f"
    --y;
#line 224 "dsbmv.f"

#line 224 "dsbmv.f"
    /* Function Body */
#line 224 "dsbmv.f"
    info = 0;
#line 225 "dsbmv.f"
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
#line 226 "dsbmv.f"
	info = 1;
#line 227 "dsbmv.f"
    } else if (*n < 0) {
#line 228 "dsbmv.f"
	info = 2;
#line 229 "dsbmv.f"
    } else if (*k < 0) {
#line 230 "dsbmv.f"
	info = 3;
#line 231 "dsbmv.f"
    } else if (*lda < *k + 1) {
#line 232 "dsbmv.f"
	info = 6;
#line 233 "dsbmv.f"
    } else if (*incx == 0) {
#line 234 "dsbmv.f"
	info = 8;
#line 235 "dsbmv.f"
    } else if (*incy == 0) {
#line 236 "dsbmv.f"
	info = 11;
#line 237 "dsbmv.f"
    }
#line 238 "dsbmv.f"
    if (info != 0) {
#line 239 "dsbmv.f"
	xerbla_("DSBMV ", &info, (ftnlen)6);
#line 240 "dsbmv.f"
	return 0;
#line 241 "dsbmv.f"
    }

/*     Quick return if possible. */

#line 245 "dsbmv.f"
    if (*n == 0 || *alpha == 0. && *beta == 1.) {
#line 245 "dsbmv.f"
	return 0;
#line 245 "dsbmv.f"
    }

/*     Set up the start points in  X  and  Y. */

#line 249 "dsbmv.f"
    if (*incx > 0) {
#line 250 "dsbmv.f"
	kx = 1;
#line 251 "dsbmv.f"
    } else {
#line 252 "dsbmv.f"
	kx = 1 - (*n - 1) * *incx;
#line 253 "dsbmv.f"
    }
#line 254 "dsbmv.f"
    if (*incy > 0) {
#line 255 "dsbmv.f"
	ky = 1;
#line 256 "dsbmv.f"
    } else {
#line 257 "dsbmv.f"
	ky = 1 - (*n - 1) * *incy;
#line 258 "dsbmv.f"
    }

/*     Start the operations. In this version the elements of the array A */
/*     are accessed sequentially with one pass through A. */

/*     First form  y := beta*y. */

#line 265 "dsbmv.f"
    if (*beta != 1.) {
#line 266 "dsbmv.f"
	if (*incy == 1) {
#line 267 "dsbmv.f"
	    if (*beta == 0.) {
#line 268 "dsbmv.f"
		i__1 = *n;
#line 268 "dsbmv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 269 "dsbmv.f"
		    y[i__] = 0.;
#line 270 "dsbmv.f"
/* L10: */
#line 270 "dsbmv.f"
		}
#line 271 "dsbmv.f"
	    } else {
#line 272 "dsbmv.f"
		i__1 = *n;
#line 272 "dsbmv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 273 "dsbmv.f"
		    y[i__] = *beta * y[i__];
#line 274 "dsbmv.f"
/* L20: */
#line 274 "dsbmv.f"
		}
#line 275 "dsbmv.f"
	    }
#line 276 "dsbmv.f"
	} else {
#line 277 "dsbmv.f"
	    iy = ky;
#line 278 "dsbmv.f"
	    if (*beta == 0.) {
#line 279 "dsbmv.f"
		i__1 = *n;
#line 279 "dsbmv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 280 "dsbmv.f"
		    y[iy] = 0.;
#line 281 "dsbmv.f"
		    iy += *incy;
#line 282 "dsbmv.f"
/* L30: */
#line 282 "dsbmv.f"
		}
#line 283 "dsbmv.f"
	    } else {
#line 284 "dsbmv.f"
		i__1 = *n;
#line 284 "dsbmv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 285 "dsbmv.f"
		    y[iy] = *beta * y[iy];
#line 286 "dsbmv.f"
		    iy += *incy;
#line 287 "dsbmv.f"
/* L40: */
#line 287 "dsbmv.f"
		}
#line 288 "dsbmv.f"
	    }
#line 289 "dsbmv.f"
	}
#line 290 "dsbmv.f"
    }
#line 291 "dsbmv.f"
    if (*alpha == 0.) {
#line 291 "dsbmv.f"
	return 0;
#line 291 "dsbmv.f"
    }
#line 292 "dsbmv.f"
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {

/*        Form  y  when upper triangle of A is stored. */

#line 296 "dsbmv.f"
	kplus1 = *k + 1;
#line 297 "dsbmv.f"
	if (*incx == 1 && *incy == 1) {
#line 298 "dsbmv.f"
	    i__1 = *n;
#line 298 "dsbmv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 299 "dsbmv.f"
		temp1 = *alpha * x[j];
#line 300 "dsbmv.f"
		temp2 = 0.;
#line 301 "dsbmv.f"
		l = kplus1 - j;
/* Computing MAX */
#line 302 "dsbmv.f"
		i__2 = 1, i__3 = j - *k;
#line 302 "dsbmv.f"
		i__4 = j - 1;
#line 302 "dsbmv.f"
		for (i__ = max(i__2,i__3); i__ <= i__4; ++i__) {
#line 303 "dsbmv.f"
		    y[i__] += temp1 * a[l + i__ + j * a_dim1];
#line 304 "dsbmv.f"
		    temp2 += a[l + i__ + j * a_dim1] * x[i__];
#line 305 "dsbmv.f"
/* L50: */
#line 305 "dsbmv.f"
		}
#line 306 "dsbmv.f"
		y[j] = y[j] + temp1 * a[kplus1 + j * a_dim1] + *alpha * temp2;
#line 307 "dsbmv.f"
/* L60: */
#line 307 "dsbmv.f"
	    }
#line 308 "dsbmv.f"
	} else {
#line 309 "dsbmv.f"
	    jx = kx;
#line 310 "dsbmv.f"
	    jy = ky;
#line 311 "dsbmv.f"
	    i__1 = *n;
#line 311 "dsbmv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 312 "dsbmv.f"
		temp1 = *alpha * x[jx];
#line 313 "dsbmv.f"
		temp2 = 0.;
#line 314 "dsbmv.f"
		ix = kx;
#line 315 "dsbmv.f"
		iy = ky;
#line 316 "dsbmv.f"
		l = kplus1 - j;
/* Computing MAX */
#line 317 "dsbmv.f"
		i__4 = 1, i__2 = j - *k;
#line 317 "dsbmv.f"
		i__3 = j - 1;
#line 317 "dsbmv.f"
		for (i__ = max(i__4,i__2); i__ <= i__3; ++i__) {
#line 318 "dsbmv.f"
		    y[iy] += temp1 * a[l + i__ + j * a_dim1];
#line 319 "dsbmv.f"
		    temp2 += a[l + i__ + j * a_dim1] * x[ix];
#line 320 "dsbmv.f"
		    ix += *incx;
#line 321 "dsbmv.f"
		    iy += *incy;
#line 322 "dsbmv.f"
/* L70: */
#line 322 "dsbmv.f"
		}
#line 323 "dsbmv.f"
		y[jy] = y[jy] + temp1 * a[kplus1 + j * a_dim1] + *alpha * 
			temp2;
#line 324 "dsbmv.f"
		jx += *incx;
#line 325 "dsbmv.f"
		jy += *incy;
#line 326 "dsbmv.f"
		if (j > *k) {
#line 327 "dsbmv.f"
		    kx += *incx;
#line 328 "dsbmv.f"
		    ky += *incy;
#line 329 "dsbmv.f"
		}
#line 330 "dsbmv.f"
/* L80: */
#line 330 "dsbmv.f"
	    }
#line 331 "dsbmv.f"
	}
#line 332 "dsbmv.f"
    } else {

/*        Form  y  when lower triangle of A is stored. */

#line 336 "dsbmv.f"
	if (*incx == 1 && *incy == 1) {
#line 337 "dsbmv.f"
	    i__1 = *n;
#line 337 "dsbmv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 338 "dsbmv.f"
		temp1 = *alpha * x[j];
#line 339 "dsbmv.f"
		temp2 = 0.;
#line 340 "dsbmv.f"
		y[j] += temp1 * a[j * a_dim1 + 1];
#line 341 "dsbmv.f"
		l = 1 - j;
/* Computing MIN */
#line 342 "dsbmv.f"
		i__4 = *n, i__2 = j + *k;
#line 342 "dsbmv.f"
		i__3 = min(i__4,i__2);
#line 342 "dsbmv.f"
		for (i__ = j + 1; i__ <= i__3; ++i__) {
#line 343 "dsbmv.f"
		    y[i__] += temp1 * a[l + i__ + j * a_dim1];
#line 344 "dsbmv.f"
		    temp2 += a[l + i__ + j * a_dim1] * x[i__];
#line 345 "dsbmv.f"
/* L90: */
#line 345 "dsbmv.f"
		}
#line 346 "dsbmv.f"
		y[j] += *alpha * temp2;
#line 347 "dsbmv.f"
/* L100: */
#line 347 "dsbmv.f"
	    }
#line 348 "dsbmv.f"
	} else {
#line 349 "dsbmv.f"
	    jx = kx;
#line 350 "dsbmv.f"
	    jy = ky;
#line 351 "dsbmv.f"
	    i__1 = *n;
#line 351 "dsbmv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 352 "dsbmv.f"
		temp1 = *alpha * x[jx];
#line 353 "dsbmv.f"
		temp2 = 0.;
#line 354 "dsbmv.f"
		y[jy] += temp1 * a[j * a_dim1 + 1];
#line 355 "dsbmv.f"
		l = 1 - j;
#line 356 "dsbmv.f"
		ix = jx;
#line 357 "dsbmv.f"
		iy = jy;
/* Computing MIN */
#line 358 "dsbmv.f"
		i__4 = *n, i__2 = j + *k;
#line 358 "dsbmv.f"
		i__3 = min(i__4,i__2);
#line 358 "dsbmv.f"
		for (i__ = j + 1; i__ <= i__3; ++i__) {
#line 359 "dsbmv.f"
		    ix += *incx;
#line 360 "dsbmv.f"
		    iy += *incy;
#line 361 "dsbmv.f"
		    y[iy] += temp1 * a[l + i__ + j * a_dim1];
#line 362 "dsbmv.f"
		    temp2 += a[l + i__ + j * a_dim1] * x[ix];
#line 363 "dsbmv.f"
/* L110: */
#line 363 "dsbmv.f"
		}
#line 364 "dsbmv.f"
		y[jy] += *alpha * temp2;
#line 365 "dsbmv.f"
		jx += *incx;
#line 366 "dsbmv.f"
		jy += *incy;
#line 367 "dsbmv.f"
/* L120: */
#line 367 "dsbmv.f"
	    }
#line 368 "dsbmv.f"
	}
#line 369 "dsbmv.f"
    }

#line 371 "dsbmv.f"
    return 0;

/*     End of DSBMV . */

} /* dsbmv_ */

