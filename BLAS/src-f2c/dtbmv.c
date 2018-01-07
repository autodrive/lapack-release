#line 1 "dtbmv.f"
/* dtbmv.f -- translated by f2c (version 20100827).
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

#line 1 "dtbmv.f"
/* > \brief \b DTBMV */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DTBMV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX) */

/*       .. Scalar Arguments .. */
/*       INTEGER INCX,K,LDA,N */
/*       CHARACTER DIAG,TRANS,UPLO */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION A(LDA,*),X(*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DTBMV  performs one of the matrix-vector operations */
/* > */
/* >    x := A*x,   or   x := A**T*x, */
/* > */
/* > where x is an n element vector and  A is an n by n unit, or non-unit, */
/* > upper or lower triangular band matrix, with ( k + 1 ) diagonals. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >           On entry, UPLO specifies whether the matrix is an upper or */
/* >           lower triangular matrix as follows: */
/* > */
/* >              UPLO = 'U' or 'u'   A is an upper triangular matrix. */
/* > */
/* >              UPLO = 'L' or 'l'   A is a lower triangular matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] TRANS */
/* > \verbatim */
/* >          TRANS is CHARACTER*1 */
/* >           On entry, TRANS specifies the operation to be performed as */
/* >           follows: */
/* > */
/* >              TRANS = 'N' or 'n'   x := A*x. */
/* > */
/* >              TRANS = 'T' or 't'   x := A**T*x. */
/* > */
/* >              TRANS = 'C' or 'c'   x := A**T*x. */
/* > \endverbatim */
/* > */
/* > \param[in] DIAG */
/* > \verbatim */
/* >          DIAG is CHARACTER*1 */
/* >           On entry, DIAG specifies whether or not A is unit */
/* >           triangular as follows: */
/* > */
/* >              DIAG = 'U' or 'u'   A is assumed to be unit triangular. */
/* > */
/* >              DIAG = 'N' or 'n'   A is not assumed to be unit */
/* >                                  triangular. */
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
/* >           On entry with UPLO = 'U' or 'u', K specifies the number of */
/* >           super-diagonals of the matrix A. */
/* >           On entry with UPLO = 'L' or 'l', K specifies the number of */
/* >           sub-diagonals of the matrix A. */
/* >           K must satisfy  0 .le. K. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array of DIMENSION ( LDA, n ). */
/* >           Before entry with UPLO = 'U' or 'u', the leading ( k + 1 ) */
/* >           by n part of the array A must contain the upper triangular */
/* >           band part of the matrix of coefficients, supplied column by */
/* >           column, with the leading diagonal of the matrix in row */
/* >           ( k + 1 ) of the array, the first super-diagonal starting at */
/* >           position 2 in row k, and so on. The top left k by k triangle */
/* >           of the array A is not referenced. */
/* >           The following program segment will transfer an upper */
/* >           triangular band matrix from conventional full matrix storage */
/* >           to band storage: */
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
/* >           band part of the matrix of coefficients, supplied column by */
/* >           column, with the leading diagonal of the matrix in row 1 of */
/* >           the array, the first sub-diagonal starting at position 1 in */
/* >           row 2, and so on. The bottom right k by k triangle of the */
/* >           array A is not referenced. */
/* >           The following program segment will transfer a lower */
/* >           triangular band matrix from conventional full matrix storage */
/* >           to band storage: */
/* > */
/* >                 DO 20, J = 1, N */
/* >                    M = 1 - J */
/* >                    DO 10, I = J, MIN( N, J + K ) */
/* >                       A( M + I, J ) = matrix( I, J ) */
/* >              10    CONTINUE */
/* >              20 CONTINUE */
/* > */
/* >           Note that when DIAG = 'U' or 'u' the elements of the array A */
/* >           corresponding to the diagonal elements of the matrix are not */
/* >           referenced, but are assumed to be unity. */
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
/* > \param[in,out] X */
/* > \verbatim */
/* >          X is DOUBLE PRECISION array of dimension at least */
/* >           ( 1 + ( n - 1 )*abs( INCX ) ). */
/* >           Before entry, the incremented array X must contain the n */
/* >           element vector x. On exit, X is overwritten with the */
/* >           transformed vector x. */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* >          INCX is INTEGER */
/* >           On entry, INCX specifies the increment for the elements of */
/* >           X. INCX must not be zero. */
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
/* Subroutine */ int dtbmv_(char *uplo, char *trans, char *diag, integer *n, 
	integer *k, doublereal *a, integer *lda, doublereal *x, integer *incx,
	 ftnlen uplo_len, ftnlen trans_len, ftnlen diag_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer i__, j, l, ix, jx, kx, info;
    static doublereal temp;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer kplus1;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static logical nounit;


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

#line 226 "dtbmv.f"
    /* Parameter adjustments */
#line 226 "dtbmv.f"
    a_dim1 = *lda;
#line 226 "dtbmv.f"
    a_offset = 1 + a_dim1;
#line 226 "dtbmv.f"
    a -= a_offset;
#line 226 "dtbmv.f"
    --x;
#line 226 "dtbmv.f"

#line 226 "dtbmv.f"
    /* Function Body */
#line 226 "dtbmv.f"
    info = 0;
#line 227 "dtbmv.f"
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
#line 228 "dtbmv.f"
	info = 1;
#line 229 "dtbmv.f"
    } else if (! lsame_(trans, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, 
	    "T", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, "C", (ftnlen)1, (
	    ftnlen)1)) {
#line 231 "dtbmv.f"
	info = 2;
#line 232 "dtbmv.f"
    } else if (! lsame_(diag, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(diag, 
	    "N", (ftnlen)1, (ftnlen)1)) {
#line 233 "dtbmv.f"
	info = 3;
#line 234 "dtbmv.f"
    } else if (*n < 0) {
#line 235 "dtbmv.f"
	info = 4;
#line 236 "dtbmv.f"
    } else if (*k < 0) {
#line 237 "dtbmv.f"
	info = 5;
#line 238 "dtbmv.f"
    } else if (*lda < *k + 1) {
#line 239 "dtbmv.f"
	info = 7;
#line 240 "dtbmv.f"
    } else if (*incx == 0) {
#line 241 "dtbmv.f"
	info = 9;
#line 242 "dtbmv.f"
    }
#line 243 "dtbmv.f"
    if (info != 0) {
#line 244 "dtbmv.f"
	xerbla_("DTBMV ", &info, (ftnlen)6);
#line 245 "dtbmv.f"
	return 0;
#line 246 "dtbmv.f"
    }

/*     Quick return if possible. */

#line 250 "dtbmv.f"
    if (*n == 0) {
#line 250 "dtbmv.f"
	return 0;
#line 250 "dtbmv.f"
    }

#line 252 "dtbmv.f"
    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);

/*     Set up the start point in X if the increment is not unity. This */
/*     will be  ( N - 1 )*INCX   too small for descending loops. */

#line 257 "dtbmv.f"
    if (*incx <= 0) {
#line 258 "dtbmv.f"
	kx = 1 - (*n - 1) * *incx;
#line 259 "dtbmv.f"
    } else if (*incx != 1) {
#line 260 "dtbmv.f"
	kx = 1;
#line 261 "dtbmv.f"
    }

/*     Start the operations. In this version the elements of A are */
/*     accessed sequentially with one pass through A. */

#line 266 "dtbmv.f"
    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {

/*         Form  x := A*x. */

#line 270 "dtbmv.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 271 "dtbmv.f"
	    kplus1 = *k + 1;
#line 272 "dtbmv.f"
	    if (*incx == 1) {
#line 273 "dtbmv.f"
		i__1 = *n;
#line 273 "dtbmv.f"
		for (j = 1; j <= i__1; ++j) {
#line 274 "dtbmv.f"
		    if (x[j] != 0.) {
#line 275 "dtbmv.f"
			temp = x[j];
#line 276 "dtbmv.f"
			l = kplus1 - j;
/* Computing MAX */
#line 277 "dtbmv.f"
			i__2 = 1, i__3 = j - *k;
#line 277 "dtbmv.f"
			i__4 = j - 1;
#line 277 "dtbmv.f"
			for (i__ = max(i__2,i__3); i__ <= i__4; ++i__) {
#line 278 "dtbmv.f"
			    x[i__] += temp * a[l + i__ + j * a_dim1];
#line 279 "dtbmv.f"
/* L10: */
#line 279 "dtbmv.f"
			}
#line 280 "dtbmv.f"
			if (nounit) {
#line 280 "dtbmv.f"
			    x[j] *= a[kplus1 + j * a_dim1];
#line 280 "dtbmv.f"
			}
#line 281 "dtbmv.f"
		    }
#line 282 "dtbmv.f"
/* L20: */
#line 282 "dtbmv.f"
		}
#line 283 "dtbmv.f"
	    } else {
#line 284 "dtbmv.f"
		jx = kx;
#line 285 "dtbmv.f"
		i__1 = *n;
#line 285 "dtbmv.f"
		for (j = 1; j <= i__1; ++j) {
#line 286 "dtbmv.f"
		    if (x[jx] != 0.) {
#line 287 "dtbmv.f"
			temp = x[jx];
#line 288 "dtbmv.f"
			ix = kx;
#line 289 "dtbmv.f"
			l = kplus1 - j;
/* Computing MAX */
#line 290 "dtbmv.f"
			i__4 = 1, i__2 = j - *k;
#line 290 "dtbmv.f"
			i__3 = j - 1;
#line 290 "dtbmv.f"
			for (i__ = max(i__4,i__2); i__ <= i__3; ++i__) {
#line 291 "dtbmv.f"
			    x[ix] += temp * a[l + i__ + j * a_dim1];
#line 292 "dtbmv.f"
			    ix += *incx;
#line 293 "dtbmv.f"
/* L30: */
#line 293 "dtbmv.f"
			}
#line 294 "dtbmv.f"
			if (nounit) {
#line 294 "dtbmv.f"
			    x[jx] *= a[kplus1 + j * a_dim1];
#line 294 "dtbmv.f"
			}
#line 295 "dtbmv.f"
		    }
#line 296 "dtbmv.f"
		    jx += *incx;
#line 297 "dtbmv.f"
		    if (j > *k) {
#line 297 "dtbmv.f"
			kx += *incx;
#line 297 "dtbmv.f"
		    }
#line 298 "dtbmv.f"
/* L40: */
#line 298 "dtbmv.f"
		}
#line 299 "dtbmv.f"
	    }
#line 300 "dtbmv.f"
	} else {
#line 301 "dtbmv.f"
	    if (*incx == 1) {
#line 302 "dtbmv.f"
		for (j = *n; j >= 1; --j) {
#line 303 "dtbmv.f"
		    if (x[j] != 0.) {
#line 304 "dtbmv.f"
			temp = x[j];
#line 305 "dtbmv.f"
			l = 1 - j;
/* Computing MIN */
#line 306 "dtbmv.f"
			i__1 = *n, i__3 = j + *k;
#line 306 "dtbmv.f"
			i__4 = j + 1;
#line 306 "dtbmv.f"
			for (i__ = min(i__1,i__3); i__ >= i__4; --i__) {
#line 307 "dtbmv.f"
			    x[i__] += temp * a[l + i__ + j * a_dim1];
#line 308 "dtbmv.f"
/* L50: */
#line 308 "dtbmv.f"
			}
#line 309 "dtbmv.f"
			if (nounit) {
#line 309 "dtbmv.f"
			    x[j] *= a[j * a_dim1 + 1];
#line 309 "dtbmv.f"
			}
#line 310 "dtbmv.f"
		    }
#line 311 "dtbmv.f"
/* L60: */
#line 311 "dtbmv.f"
		}
#line 312 "dtbmv.f"
	    } else {
#line 313 "dtbmv.f"
		kx += (*n - 1) * *incx;
#line 314 "dtbmv.f"
		jx = kx;
#line 315 "dtbmv.f"
		for (j = *n; j >= 1; --j) {
#line 316 "dtbmv.f"
		    if (x[jx] != 0.) {
#line 317 "dtbmv.f"
			temp = x[jx];
#line 318 "dtbmv.f"
			ix = kx;
#line 319 "dtbmv.f"
			l = 1 - j;
/* Computing MIN */
#line 320 "dtbmv.f"
			i__4 = *n, i__1 = j + *k;
#line 320 "dtbmv.f"
			i__3 = j + 1;
#line 320 "dtbmv.f"
			for (i__ = min(i__4,i__1); i__ >= i__3; --i__) {
#line 321 "dtbmv.f"
			    x[ix] += temp * a[l + i__ + j * a_dim1];
#line 322 "dtbmv.f"
			    ix -= *incx;
#line 323 "dtbmv.f"
/* L70: */
#line 323 "dtbmv.f"
			}
#line 324 "dtbmv.f"
			if (nounit) {
#line 324 "dtbmv.f"
			    x[jx] *= a[j * a_dim1 + 1];
#line 324 "dtbmv.f"
			}
#line 325 "dtbmv.f"
		    }
#line 326 "dtbmv.f"
		    jx -= *incx;
#line 327 "dtbmv.f"
		    if (*n - j >= *k) {
#line 327 "dtbmv.f"
			kx -= *incx;
#line 327 "dtbmv.f"
		    }
#line 328 "dtbmv.f"
/* L80: */
#line 328 "dtbmv.f"
		}
#line 329 "dtbmv.f"
	    }
#line 330 "dtbmv.f"
	}
#line 331 "dtbmv.f"
    } else {

/*        Form  x := A**T*x. */

#line 335 "dtbmv.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 336 "dtbmv.f"
	    kplus1 = *k + 1;
#line 337 "dtbmv.f"
	    if (*incx == 1) {
#line 338 "dtbmv.f"
		for (j = *n; j >= 1; --j) {
#line 339 "dtbmv.f"
		    temp = x[j];
#line 340 "dtbmv.f"
		    l = kplus1 - j;
#line 341 "dtbmv.f"
		    if (nounit) {
#line 341 "dtbmv.f"
			temp *= a[kplus1 + j * a_dim1];
#line 341 "dtbmv.f"
		    }
/* Computing MAX */
#line 342 "dtbmv.f"
		    i__4 = 1, i__1 = j - *k;
#line 342 "dtbmv.f"
		    i__3 = max(i__4,i__1);
#line 342 "dtbmv.f"
		    for (i__ = j - 1; i__ >= i__3; --i__) {
#line 343 "dtbmv.f"
			temp += a[l + i__ + j * a_dim1] * x[i__];
#line 344 "dtbmv.f"
/* L90: */
#line 344 "dtbmv.f"
		    }
#line 345 "dtbmv.f"
		    x[j] = temp;
#line 346 "dtbmv.f"
/* L100: */
#line 346 "dtbmv.f"
		}
#line 347 "dtbmv.f"
	    } else {
#line 348 "dtbmv.f"
		kx += (*n - 1) * *incx;
#line 349 "dtbmv.f"
		jx = kx;
#line 350 "dtbmv.f"
		for (j = *n; j >= 1; --j) {
#line 351 "dtbmv.f"
		    temp = x[jx];
#line 352 "dtbmv.f"
		    kx -= *incx;
#line 353 "dtbmv.f"
		    ix = kx;
#line 354 "dtbmv.f"
		    l = kplus1 - j;
#line 355 "dtbmv.f"
		    if (nounit) {
#line 355 "dtbmv.f"
			temp *= a[kplus1 + j * a_dim1];
#line 355 "dtbmv.f"
		    }
/* Computing MAX */
#line 356 "dtbmv.f"
		    i__4 = 1, i__1 = j - *k;
#line 356 "dtbmv.f"
		    i__3 = max(i__4,i__1);
#line 356 "dtbmv.f"
		    for (i__ = j - 1; i__ >= i__3; --i__) {
#line 357 "dtbmv.f"
			temp += a[l + i__ + j * a_dim1] * x[ix];
#line 358 "dtbmv.f"
			ix -= *incx;
#line 359 "dtbmv.f"
/* L110: */
#line 359 "dtbmv.f"
		    }
#line 360 "dtbmv.f"
		    x[jx] = temp;
#line 361 "dtbmv.f"
		    jx -= *incx;
#line 362 "dtbmv.f"
/* L120: */
#line 362 "dtbmv.f"
		}
#line 363 "dtbmv.f"
	    }
#line 364 "dtbmv.f"
	} else {
#line 365 "dtbmv.f"
	    if (*incx == 1) {
#line 366 "dtbmv.f"
		i__3 = *n;
#line 366 "dtbmv.f"
		for (j = 1; j <= i__3; ++j) {
#line 367 "dtbmv.f"
		    temp = x[j];
#line 368 "dtbmv.f"
		    l = 1 - j;
#line 369 "dtbmv.f"
		    if (nounit) {
#line 369 "dtbmv.f"
			temp *= a[j * a_dim1 + 1];
#line 369 "dtbmv.f"
		    }
/* Computing MIN */
#line 370 "dtbmv.f"
		    i__1 = *n, i__2 = j + *k;
#line 370 "dtbmv.f"
		    i__4 = min(i__1,i__2);
#line 370 "dtbmv.f"
		    for (i__ = j + 1; i__ <= i__4; ++i__) {
#line 371 "dtbmv.f"
			temp += a[l + i__ + j * a_dim1] * x[i__];
#line 372 "dtbmv.f"
/* L130: */
#line 372 "dtbmv.f"
		    }
#line 373 "dtbmv.f"
		    x[j] = temp;
#line 374 "dtbmv.f"
/* L140: */
#line 374 "dtbmv.f"
		}
#line 375 "dtbmv.f"
	    } else {
#line 376 "dtbmv.f"
		jx = kx;
#line 377 "dtbmv.f"
		i__3 = *n;
#line 377 "dtbmv.f"
		for (j = 1; j <= i__3; ++j) {
#line 378 "dtbmv.f"
		    temp = x[jx];
#line 379 "dtbmv.f"
		    kx += *incx;
#line 380 "dtbmv.f"
		    ix = kx;
#line 381 "dtbmv.f"
		    l = 1 - j;
#line 382 "dtbmv.f"
		    if (nounit) {
#line 382 "dtbmv.f"
			temp *= a[j * a_dim1 + 1];
#line 382 "dtbmv.f"
		    }
/* Computing MIN */
#line 383 "dtbmv.f"
		    i__1 = *n, i__2 = j + *k;
#line 383 "dtbmv.f"
		    i__4 = min(i__1,i__2);
#line 383 "dtbmv.f"
		    for (i__ = j + 1; i__ <= i__4; ++i__) {
#line 384 "dtbmv.f"
			temp += a[l + i__ + j * a_dim1] * x[ix];
#line 385 "dtbmv.f"
			ix += *incx;
#line 386 "dtbmv.f"
/* L150: */
#line 386 "dtbmv.f"
		    }
#line 387 "dtbmv.f"
		    x[jx] = temp;
#line 388 "dtbmv.f"
		    jx += *incx;
#line 389 "dtbmv.f"
/* L160: */
#line 389 "dtbmv.f"
		}
#line 390 "dtbmv.f"
	    }
#line 391 "dtbmv.f"
	}
#line 392 "dtbmv.f"
    }

#line 394 "dtbmv.f"
    return 0;

/*     End of DTBMV . */

} /* dtbmv_ */

