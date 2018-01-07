#line 1 "dtbsv.f"
/* dtbsv.f -- translated by f2c (version 20100827).
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

#line 1 "dtbsv.f"
/* > \brief \b DTBSV */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DTBSV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX) */

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
/* > DTBSV  solves one of the systems of equations */
/* > */
/* >    A*x = b,   or   A**T*x = b, */
/* > */
/* > where b and x are n element vectors and A is an n by n unit, or */
/* > non-unit, upper or lower triangular band matrix, with ( k + 1 ) */
/* > diagonals. */
/* > */
/* > No test for singularity or near-singularity is included in this */
/* > routine. Such tests must be performed before calling this routine. */
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
/* >           On entry, TRANS specifies the equations to be solved as */
/* >           follows: */
/* > */
/* >              TRANS = 'N' or 'n'   A*x = b. */
/* > */
/* >              TRANS = 'T' or 't'   A**T*x = b. */
/* > */
/* >              TRANS = 'C' or 'c'   A**T*x = b. */
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
/* >           element right-hand side vector b. On exit, X is overwritten */
/* >           with the solution vector x. */
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

/* > \date November 2011 */

/* > \ingroup double_blas_level2 */

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
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dtbsv_(char *uplo, char *trans, char *diag, integer *n, 
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

#line 229 "dtbsv.f"
    /* Parameter adjustments */
#line 229 "dtbsv.f"
    a_dim1 = *lda;
#line 229 "dtbsv.f"
    a_offset = 1 + a_dim1;
#line 229 "dtbsv.f"
    a -= a_offset;
#line 229 "dtbsv.f"
    --x;
#line 229 "dtbsv.f"

#line 229 "dtbsv.f"
    /* Function Body */
#line 229 "dtbsv.f"
    info = 0;
#line 230 "dtbsv.f"
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
#line 231 "dtbsv.f"
	info = 1;
#line 232 "dtbsv.f"
    } else if (! lsame_(trans, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, 
	    "T", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, "C", (ftnlen)1, (
	    ftnlen)1)) {
#line 234 "dtbsv.f"
	info = 2;
#line 235 "dtbsv.f"
    } else if (! lsame_(diag, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(diag, 
	    "N", (ftnlen)1, (ftnlen)1)) {
#line 236 "dtbsv.f"
	info = 3;
#line 237 "dtbsv.f"
    } else if (*n < 0) {
#line 238 "dtbsv.f"
	info = 4;
#line 239 "dtbsv.f"
    } else if (*k < 0) {
#line 240 "dtbsv.f"
	info = 5;
#line 241 "dtbsv.f"
    } else if (*lda < *k + 1) {
#line 242 "dtbsv.f"
	info = 7;
#line 243 "dtbsv.f"
    } else if (*incx == 0) {
#line 244 "dtbsv.f"
	info = 9;
#line 245 "dtbsv.f"
    }
#line 246 "dtbsv.f"
    if (info != 0) {
#line 247 "dtbsv.f"
	xerbla_("DTBSV ", &info, (ftnlen)6);
#line 248 "dtbsv.f"
	return 0;
#line 249 "dtbsv.f"
    }

/*     Quick return if possible. */

#line 253 "dtbsv.f"
    if (*n == 0) {
#line 253 "dtbsv.f"
	return 0;
#line 253 "dtbsv.f"
    }

#line 255 "dtbsv.f"
    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);

/*     Set up the start point in X if the increment is not unity. This */
/*     will be  ( N - 1 )*INCX  too small for descending loops. */

#line 260 "dtbsv.f"
    if (*incx <= 0) {
#line 261 "dtbsv.f"
	kx = 1 - (*n - 1) * *incx;
#line 262 "dtbsv.f"
    } else if (*incx != 1) {
#line 263 "dtbsv.f"
	kx = 1;
#line 264 "dtbsv.f"
    }

/*     Start the operations. In this version the elements of A are */
/*     accessed by sequentially with one pass through A. */

#line 269 "dtbsv.f"
    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {

/*        Form  x := inv( A )*x. */

#line 273 "dtbsv.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 274 "dtbsv.f"
	    kplus1 = *k + 1;
#line 275 "dtbsv.f"
	    if (*incx == 1) {
#line 276 "dtbsv.f"
		for (j = *n; j >= 1; --j) {
#line 277 "dtbsv.f"
		    if (x[j] != 0.) {
#line 278 "dtbsv.f"
			l = kplus1 - j;
#line 279 "dtbsv.f"
			if (nounit) {
#line 279 "dtbsv.f"
			    x[j] /= a[kplus1 + j * a_dim1];
#line 279 "dtbsv.f"
			}
#line 280 "dtbsv.f"
			temp = x[j];
/* Computing MAX */
#line 281 "dtbsv.f"
			i__2 = 1, i__3 = j - *k;
#line 281 "dtbsv.f"
			i__1 = max(i__2,i__3);
#line 281 "dtbsv.f"
			for (i__ = j - 1; i__ >= i__1; --i__) {
#line 282 "dtbsv.f"
			    x[i__] -= temp * a[l + i__ + j * a_dim1];
#line 283 "dtbsv.f"
/* L10: */
#line 283 "dtbsv.f"
			}
#line 284 "dtbsv.f"
		    }
#line 285 "dtbsv.f"
/* L20: */
#line 285 "dtbsv.f"
		}
#line 286 "dtbsv.f"
	    } else {
#line 287 "dtbsv.f"
		kx += (*n - 1) * *incx;
#line 288 "dtbsv.f"
		jx = kx;
#line 289 "dtbsv.f"
		for (j = *n; j >= 1; --j) {
#line 290 "dtbsv.f"
		    kx -= *incx;
#line 291 "dtbsv.f"
		    if (x[jx] != 0.) {
#line 292 "dtbsv.f"
			ix = kx;
#line 293 "dtbsv.f"
			l = kplus1 - j;
#line 294 "dtbsv.f"
			if (nounit) {
#line 294 "dtbsv.f"
			    x[jx] /= a[kplus1 + j * a_dim1];
#line 294 "dtbsv.f"
			}
#line 295 "dtbsv.f"
			temp = x[jx];
/* Computing MAX */
#line 296 "dtbsv.f"
			i__2 = 1, i__3 = j - *k;
#line 296 "dtbsv.f"
			i__1 = max(i__2,i__3);
#line 296 "dtbsv.f"
			for (i__ = j - 1; i__ >= i__1; --i__) {
#line 297 "dtbsv.f"
			    x[ix] -= temp * a[l + i__ + j * a_dim1];
#line 298 "dtbsv.f"
			    ix -= *incx;
#line 299 "dtbsv.f"
/* L30: */
#line 299 "dtbsv.f"
			}
#line 300 "dtbsv.f"
		    }
#line 301 "dtbsv.f"
		    jx -= *incx;
#line 302 "dtbsv.f"
/* L40: */
#line 302 "dtbsv.f"
		}
#line 303 "dtbsv.f"
	    }
#line 304 "dtbsv.f"
	} else {
#line 305 "dtbsv.f"
	    if (*incx == 1) {
#line 306 "dtbsv.f"
		i__1 = *n;
#line 306 "dtbsv.f"
		for (j = 1; j <= i__1; ++j) {
#line 307 "dtbsv.f"
		    if (x[j] != 0.) {
#line 308 "dtbsv.f"
			l = 1 - j;
#line 309 "dtbsv.f"
			if (nounit) {
#line 309 "dtbsv.f"
			    x[j] /= a[j * a_dim1 + 1];
#line 309 "dtbsv.f"
			}
#line 310 "dtbsv.f"
			temp = x[j];
/* Computing MIN */
#line 311 "dtbsv.f"
			i__3 = *n, i__4 = j + *k;
#line 311 "dtbsv.f"
			i__2 = min(i__3,i__4);
#line 311 "dtbsv.f"
			for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 312 "dtbsv.f"
			    x[i__] -= temp * a[l + i__ + j * a_dim1];
#line 313 "dtbsv.f"
/* L50: */
#line 313 "dtbsv.f"
			}
#line 314 "dtbsv.f"
		    }
#line 315 "dtbsv.f"
/* L60: */
#line 315 "dtbsv.f"
		}
#line 316 "dtbsv.f"
	    } else {
#line 317 "dtbsv.f"
		jx = kx;
#line 318 "dtbsv.f"
		i__1 = *n;
#line 318 "dtbsv.f"
		for (j = 1; j <= i__1; ++j) {
#line 319 "dtbsv.f"
		    kx += *incx;
#line 320 "dtbsv.f"
		    if (x[jx] != 0.) {
#line 321 "dtbsv.f"
			ix = kx;
#line 322 "dtbsv.f"
			l = 1 - j;
#line 323 "dtbsv.f"
			if (nounit) {
#line 323 "dtbsv.f"
			    x[jx] /= a[j * a_dim1 + 1];
#line 323 "dtbsv.f"
			}
#line 324 "dtbsv.f"
			temp = x[jx];
/* Computing MIN */
#line 325 "dtbsv.f"
			i__3 = *n, i__4 = j + *k;
#line 325 "dtbsv.f"
			i__2 = min(i__3,i__4);
#line 325 "dtbsv.f"
			for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 326 "dtbsv.f"
			    x[ix] -= temp * a[l + i__ + j * a_dim1];
#line 327 "dtbsv.f"
			    ix += *incx;
#line 328 "dtbsv.f"
/* L70: */
#line 328 "dtbsv.f"
			}
#line 329 "dtbsv.f"
		    }
#line 330 "dtbsv.f"
		    jx += *incx;
#line 331 "dtbsv.f"
/* L80: */
#line 331 "dtbsv.f"
		}
#line 332 "dtbsv.f"
	    }
#line 333 "dtbsv.f"
	}
#line 334 "dtbsv.f"
    } else {

/*        Form  x := inv( A**T)*x. */

#line 338 "dtbsv.f"
	if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
#line 339 "dtbsv.f"
	    kplus1 = *k + 1;
#line 340 "dtbsv.f"
	    if (*incx == 1) {
#line 341 "dtbsv.f"
		i__1 = *n;
#line 341 "dtbsv.f"
		for (j = 1; j <= i__1; ++j) {
#line 342 "dtbsv.f"
		    temp = x[j];
#line 343 "dtbsv.f"
		    l = kplus1 - j;
/* Computing MAX */
#line 344 "dtbsv.f"
		    i__2 = 1, i__3 = j - *k;
#line 344 "dtbsv.f"
		    i__4 = j - 1;
#line 344 "dtbsv.f"
		    for (i__ = max(i__2,i__3); i__ <= i__4; ++i__) {
#line 345 "dtbsv.f"
			temp -= a[l + i__ + j * a_dim1] * x[i__];
#line 346 "dtbsv.f"
/* L90: */
#line 346 "dtbsv.f"
		    }
#line 347 "dtbsv.f"
		    if (nounit) {
#line 347 "dtbsv.f"
			temp /= a[kplus1 + j * a_dim1];
#line 347 "dtbsv.f"
		    }
#line 348 "dtbsv.f"
		    x[j] = temp;
#line 349 "dtbsv.f"
/* L100: */
#line 349 "dtbsv.f"
		}
#line 350 "dtbsv.f"
	    } else {
#line 351 "dtbsv.f"
		jx = kx;
#line 352 "dtbsv.f"
		i__1 = *n;
#line 352 "dtbsv.f"
		for (j = 1; j <= i__1; ++j) {
#line 353 "dtbsv.f"
		    temp = x[jx];
#line 354 "dtbsv.f"
		    ix = kx;
#line 355 "dtbsv.f"
		    l = kplus1 - j;
/* Computing MAX */
#line 356 "dtbsv.f"
		    i__4 = 1, i__2 = j - *k;
#line 356 "dtbsv.f"
		    i__3 = j - 1;
#line 356 "dtbsv.f"
		    for (i__ = max(i__4,i__2); i__ <= i__3; ++i__) {
#line 357 "dtbsv.f"
			temp -= a[l + i__ + j * a_dim1] * x[ix];
#line 358 "dtbsv.f"
			ix += *incx;
#line 359 "dtbsv.f"
/* L110: */
#line 359 "dtbsv.f"
		    }
#line 360 "dtbsv.f"
		    if (nounit) {
#line 360 "dtbsv.f"
			temp /= a[kplus1 + j * a_dim1];
#line 360 "dtbsv.f"
		    }
#line 361 "dtbsv.f"
		    x[jx] = temp;
#line 362 "dtbsv.f"
		    jx += *incx;
#line 363 "dtbsv.f"
		    if (j > *k) {
#line 363 "dtbsv.f"
			kx += *incx;
#line 363 "dtbsv.f"
		    }
#line 364 "dtbsv.f"
/* L120: */
#line 364 "dtbsv.f"
		}
#line 365 "dtbsv.f"
	    }
#line 366 "dtbsv.f"
	} else {
#line 367 "dtbsv.f"
	    if (*incx == 1) {
#line 368 "dtbsv.f"
		for (j = *n; j >= 1; --j) {
#line 369 "dtbsv.f"
		    temp = x[j];
#line 370 "dtbsv.f"
		    l = 1 - j;
/* Computing MIN */
#line 371 "dtbsv.f"
		    i__1 = *n, i__3 = j + *k;
#line 371 "dtbsv.f"
		    i__4 = j + 1;
#line 371 "dtbsv.f"
		    for (i__ = min(i__1,i__3); i__ >= i__4; --i__) {
#line 372 "dtbsv.f"
			temp -= a[l + i__ + j * a_dim1] * x[i__];
#line 373 "dtbsv.f"
/* L130: */
#line 373 "dtbsv.f"
		    }
#line 374 "dtbsv.f"
		    if (nounit) {
#line 374 "dtbsv.f"
			temp /= a[j * a_dim1 + 1];
#line 374 "dtbsv.f"
		    }
#line 375 "dtbsv.f"
		    x[j] = temp;
#line 376 "dtbsv.f"
/* L140: */
#line 376 "dtbsv.f"
		}
#line 377 "dtbsv.f"
	    } else {
#line 378 "dtbsv.f"
		kx += (*n - 1) * *incx;
#line 379 "dtbsv.f"
		jx = kx;
#line 380 "dtbsv.f"
		for (j = *n; j >= 1; --j) {
#line 381 "dtbsv.f"
		    temp = x[jx];
#line 382 "dtbsv.f"
		    ix = kx;
#line 383 "dtbsv.f"
		    l = 1 - j;
/* Computing MIN */
#line 384 "dtbsv.f"
		    i__4 = *n, i__1 = j + *k;
#line 384 "dtbsv.f"
		    i__3 = j + 1;
#line 384 "dtbsv.f"
		    for (i__ = min(i__4,i__1); i__ >= i__3; --i__) {
#line 385 "dtbsv.f"
			temp -= a[l + i__ + j * a_dim1] * x[ix];
#line 386 "dtbsv.f"
			ix -= *incx;
#line 387 "dtbsv.f"
/* L150: */
#line 387 "dtbsv.f"
		    }
#line 388 "dtbsv.f"
		    if (nounit) {
#line 388 "dtbsv.f"
			temp /= a[j * a_dim1 + 1];
#line 388 "dtbsv.f"
		    }
#line 389 "dtbsv.f"
		    x[jx] = temp;
#line 390 "dtbsv.f"
		    jx -= *incx;
#line 391 "dtbsv.f"
		    if (*n - j >= *k) {
#line 391 "dtbsv.f"
			kx -= *incx;
#line 391 "dtbsv.f"
		    }
#line 392 "dtbsv.f"
/* L160: */
#line 392 "dtbsv.f"
		}
#line 393 "dtbsv.f"
	    }
#line 394 "dtbsv.f"
	}
#line 395 "dtbsv.f"
    }

#line 397 "dtbsv.f"
    return 0;

/*     End of DTBSV . */

} /* dtbsv_ */

