#line 1 "ssbmv.f"
/* ssbmv.f -- translated by f2c (version 20100827).
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

#line 1 "ssbmv.f"
/* > \brief \b SSBMV */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SSBMV(UPLO,N,K,ALPHA,A,LDA,X,INCX,BETA,Y,INCY) */

/*       .. Scalar Arguments .. */
/*       REAL ALPHA,BETA */
/*       INTEGER INCX,INCY,K,LDA,N */
/*       CHARACTER UPLO */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL A(LDA,*),X(*),Y(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SSBMV  performs the matrix-vector  operation */
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
/* >          ALPHA is REAL */
/* >           On entry, ALPHA specifies the scalar alpha. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is REAL array of DIMENSION ( LDA, n ). */
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
/* >          X is REAL array of DIMENSION at least */
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
/* >          BETA is REAL */
/* >           On entry, BETA specifies the scalar beta. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Y */
/* > \verbatim */
/* >          Y is REAL array of DIMENSION at least */
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

/* > \date November 2011 */

/* > \ingroup single_blas_level2 */

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
/* Subroutine */ int ssbmv_(char *uplo, integer *n, integer *k, doublereal *
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


/*  -- Reference BLAS level2 routine (version 3.4.0) -- */
/*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */

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

#line 224 "ssbmv.f"
    /* Parameter adjustments */
#line 224 "ssbmv.f"
    a_dim1 = *lda;
#line 224 "ssbmv.f"
    a_offset = 1 + a_dim1;
#line 224 "ssbmv.f"
    a -= a_offset;
#line 224 "ssbmv.f"
    --x;
#line 224 "ssbmv.f"
    --y;
#line 224 "ssbmv.f"

#line 224 "ssbmv.f"
    /* Function Body */
#line 224 "ssbmv.f"
    info = 0;
#line 225 "ssbmv.f"
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
#line 226 "ssbmv.f"
	info = 1;
#line 227 "ssbmv.f"
    } else if (*n < 0) {
#line 228 "ssbmv.f"
	info = 2;
#line 229 "ssbmv.f"
    } else if (*k < 0) {
#line 230 "ssbmv.f"
	info = 3;
#line 231 "ssbmv.f"
    } else if (*lda < *k + 1) {
#line 232 "ssbmv.f"
	info = 6;
#line 233 "ssbmv.f"
    } else if (*incx == 0) {
#line 234 "ssbmv.f"
	info = 8;
#line 235 "ssbmv.f"
    } else if (*incy == 0) {
#line 236 "ssbmv.f"
	info = 11;
#line 237 "ssbmv.f"
    }
#line 238 "ssbmv.f"
    if (info != 0) {
#line 239 "ssbmv.f"
	xerbla_("SSBMV ", &info, (ftnlen)6);
#line 240 "ssbmv.f"
	return 0;
#line 241 "ssbmv.f"
    }

/*     Quick return if possible. */

#line 245 "ssbmv.f"
    if (*n == 0 || *alpha == 0. && *beta == 1.) {
#line 245 "ssbmv.f"
	return 0;
#line 245 "ssbmv.f"
    }

/*     Set up the start points in  X  and  Y. */

#line 249 "ssbmv.f"
    if (*incx > 0) {
#line 250 "ssbmv.f"
	kx = 1;
#line 251 "ssbmv.f"
    } else {
#line 252 "ssbmv.f"
	kx = 1 - (*n - 1) * *incx;
#line 253 "ssbmv.f"
    }
#line 254 "ssbmv.f"
    if (*incy > 0) {
#line 255 "ssbmv.f"
	ky = 1;
#line 256 "ssbmv.f"
    } else {
#line 257 "ssbmv.f"
	ky = 1 - (*n - 1) * *incy;
#line 258 "ssbmv.f"
    }

/*     Start the operations. In this version the elements of the array A */
/*     are accessed sequentially with one pass through A. */

/*     First form  y := beta*y. */

#line 265 "ssbmv.f"
    if (*beta != 1.) {
#line 266 "ssbmv.f"
	if (*incy == 1) {
#line 267 "ssbmv.f"
	    if (*beta == 0.) {
#line 268 "ssbmv.f"
		i__1 = *n;
#line 268 "ssbmv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 269 "ssbmv.f"
		    y[i__] = 0.;
#line 270 "ssbmv.f"
/* L10: */
#line 270 "ssbmv.f"
		}
#line 271 "ssbmv.f"
	    } else {
#line 272 "ssbmv.f"
		i__1 = *n;
#line 272 "ssbmv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 273 "ssbmv.f"
		    y[i__] = *beta * y[i__];
#line 274 "ssbmv.f"
/* L20: */
#line 274 "ssbmv.f"
		}
#line 275 "ssbmv.f"
	    }
#line 276 "ssbmv.f"
	} else {
#line 277 "ssbmv.f"
	    iy = ky;
#line 278 "ssbmv.f"
	    if (*beta == 0.) {
#line 279 "ssbmv.f"
		i__1 = *n;
#line 279 "ssbmv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 280 "ssbmv.f"
		    y[iy] = 0.;
#line 281 "ssbmv.f"
		    iy += *incy;
#line 282 "ssbmv.f"
/* L30: */
#line 282 "ssbmv.f"
		}
#line 283 "ssbmv.f"
	    } else {
#line 284 "ssbmv.f"
		i__1 = *n;
#line 284 "ssbmv.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 285 "ssbmv.f"
		    y[iy] = *beta * y[iy];
#line 286 "ssbmv.f"
		    iy += *incy;
#line 287 "ssbmv.f"
/* L40: */
#line 287 "ssbmv.f"
		}
#line 288 "ssbmv.f"
	    }
#line 289 "ssbmv.f"
	}
#line 290 "ssbmv.f"
    }
#line 291 "ssbmv.f"
    if (*alpha == 0.) {
#line 291 "ssbmv.f"
	return 0;
#line 291 "ssbmv.f"
    }
#line 292 "ssbmv.f"
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {

/*        Form  y  when upper triangle of A is stored. */

#line 296 "ssbmv.f"
	kplus1 = *k + 1;
#line 297 "ssbmv.f"
	if (*incx == 1 && *incy == 1) {
#line 298 "ssbmv.f"
	    i__1 = *n;
#line 298 "ssbmv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 299 "ssbmv.f"
		temp1 = *alpha * x[j];
#line 300 "ssbmv.f"
		temp2 = 0.;
#line 301 "ssbmv.f"
		l = kplus1 - j;
/* Computing MAX */
#line 302 "ssbmv.f"
		i__2 = 1, i__3 = j - *k;
#line 302 "ssbmv.f"
		i__4 = j - 1;
#line 302 "ssbmv.f"
		for (i__ = max(i__2,i__3); i__ <= i__4; ++i__) {
#line 303 "ssbmv.f"
		    y[i__] += temp1 * a[l + i__ + j * a_dim1];
#line 304 "ssbmv.f"
		    temp2 += a[l + i__ + j * a_dim1] * x[i__];
#line 305 "ssbmv.f"
/* L50: */
#line 305 "ssbmv.f"
		}
#line 306 "ssbmv.f"
		y[j] = y[j] + temp1 * a[kplus1 + j * a_dim1] + *alpha * temp2;
#line 307 "ssbmv.f"
/* L60: */
#line 307 "ssbmv.f"
	    }
#line 308 "ssbmv.f"
	} else {
#line 309 "ssbmv.f"
	    jx = kx;
#line 310 "ssbmv.f"
	    jy = ky;
#line 311 "ssbmv.f"
	    i__1 = *n;
#line 311 "ssbmv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 312 "ssbmv.f"
		temp1 = *alpha * x[jx];
#line 313 "ssbmv.f"
		temp2 = 0.;
#line 314 "ssbmv.f"
		ix = kx;
#line 315 "ssbmv.f"
		iy = ky;
#line 316 "ssbmv.f"
		l = kplus1 - j;
/* Computing MAX */
#line 317 "ssbmv.f"
		i__4 = 1, i__2 = j - *k;
#line 317 "ssbmv.f"
		i__3 = j - 1;
#line 317 "ssbmv.f"
		for (i__ = max(i__4,i__2); i__ <= i__3; ++i__) {
#line 318 "ssbmv.f"
		    y[iy] += temp1 * a[l + i__ + j * a_dim1];
#line 319 "ssbmv.f"
		    temp2 += a[l + i__ + j * a_dim1] * x[ix];
#line 320 "ssbmv.f"
		    ix += *incx;
#line 321 "ssbmv.f"
		    iy += *incy;
#line 322 "ssbmv.f"
/* L70: */
#line 322 "ssbmv.f"
		}
#line 323 "ssbmv.f"
		y[jy] = y[jy] + temp1 * a[kplus1 + j * a_dim1] + *alpha * 
			temp2;
#line 324 "ssbmv.f"
		jx += *incx;
#line 325 "ssbmv.f"
		jy += *incy;
#line 326 "ssbmv.f"
		if (j > *k) {
#line 327 "ssbmv.f"
		    kx += *incx;
#line 328 "ssbmv.f"
		    ky += *incy;
#line 329 "ssbmv.f"
		}
#line 330 "ssbmv.f"
/* L80: */
#line 330 "ssbmv.f"
	    }
#line 331 "ssbmv.f"
	}
#line 332 "ssbmv.f"
    } else {

/*        Form  y  when lower triangle of A is stored. */

#line 336 "ssbmv.f"
	if (*incx == 1 && *incy == 1) {
#line 337 "ssbmv.f"
	    i__1 = *n;
#line 337 "ssbmv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 338 "ssbmv.f"
		temp1 = *alpha * x[j];
#line 339 "ssbmv.f"
		temp2 = 0.;
#line 340 "ssbmv.f"
		y[j] += temp1 * a[j * a_dim1 + 1];
#line 341 "ssbmv.f"
		l = 1 - j;
/* Computing MIN */
#line 342 "ssbmv.f"
		i__4 = *n, i__2 = j + *k;
#line 342 "ssbmv.f"
		i__3 = min(i__4,i__2);
#line 342 "ssbmv.f"
		for (i__ = j + 1; i__ <= i__3; ++i__) {
#line 343 "ssbmv.f"
		    y[i__] += temp1 * a[l + i__ + j * a_dim1];
#line 344 "ssbmv.f"
		    temp2 += a[l + i__ + j * a_dim1] * x[i__];
#line 345 "ssbmv.f"
/* L90: */
#line 345 "ssbmv.f"
		}
#line 346 "ssbmv.f"
		y[j] += *alpha * temp2;
#line 347 "ssbmv.f"
/* L100: */
#line 347 "ssbmv.f"
	    }
#line 348 "ssbmv.f"
	} else {
#line 349 "ssbmv.f"
	    jx = kx;
#line 350 "ssbmv.f"
	    jy = ky;
#line 351 "ssbmv.f"
	    i__1 = *n;
#line 351 "ssbmv.f"
	    for (j = 1; j <= i__1; ++j) {
#line 352 "ssbmv.f"
		temp1 = *alpha * x[jx];
#line 353 "ssbmv.f"
		temp2 = 0.;
#line 354 "ssbmv.f"
		y[jy] += temp1 * a[j * a_dim1 + 1];
#line 355 "ssbmv.f"
		l = 1 - j;
#line 356 "ssbmv.f"
		ix = jx;
#line 357 "ssbmv.f"
		iy = jy;
/* Computing MIN */
#line 358 "ssbmv.f"
		i__4 = *n, i__2 = j + *k;
#line 358 "ssbmv.f"
		i__3 = min(i__4,i__2);
#line 358 "ssbmv.f"
		for (i__ = j + 1; i__ <= i__3; ++i__) {
#line 359 "ssbmv.f"
		    ix += *incx;
#line 360 "ssbmv.f"
		    iy += *incy;
#line 361 "ssbmv.f"
		    y[iy] += temp1 * a[l + i__ + j * a_dim1];
#line 362 "ssbmv.f"
		    temp2 += a[l + i__ + j * a_dim1] * x[ix];
#line 363 "ssbmv.f"
/* L110: */
#line 363 "ssbmv.f"
		}
#line 364 "ssbmv.f"
		y[jy] += *alpha * temp2;
#line 365 "ssbmv.f"
		jx += *incx;
#line 366 "ssbmv.f"
		jy += *incy;
#line 367 "ssbmv.f"
/* L120: */
#line 367 "ssbmv.f"
	    }
#line 368 "ssbmv.f"
	}
#line 369 "ssbmv.f"
    }

#line 371 "ssbmv.f"
    return 0;

/*     End of SSBMV . */

} /* ssbmv_ */

