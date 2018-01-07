#line 1 "dtrsm.f"
/* dtrsm.f -- translated by f2c (version 20100827).
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

#line 1 "dtrsm.f"
/* > \brief \b DTRSM */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB) */

/*       .. Scalar Arguments .. */
/*       DOUBLE PRECISION ALPHA */
/*       INTEGER LDA,LDB,M,N */
/*       CHARACTER DIAG,SIDE,TRANSA,UPLO */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION A(LDA,*),B(LDB,*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DTRSM  solves one of the matrix equations */
/* > */
/* >    op( A )*X = alpha*B,   or   X*op( A ) = alpha*B, */
/* > */
/* > where alpha is a scalar, X and B are m by n matrices, A is a unit, or */
/* > non-unit,  upper or lower triangular matrix  and  op( A )  is one  of */
/* > */
/* >    op( A ) = A   or   op( A ) = A**T. */
/* > */
/* > The matrix X is overwritten on B. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] SIDE */
/* > \verbatim */
/* >          SIDE is CHARACTER*1 */
/* >           On entry, SIDE specifies whether op( A ) appears on the left */
/* >           or right of X as follows: */
/* > */
/* >              SIDE = 'L' or 'l'   op( A )*X = alpha*B. */
/* > */
/* >              SIDE = 'R' or 'r'   X*op( A ) = alpha*B. */
/* > \endverbatim */
/* > */
/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >           On entry, UPLO specifies whether the matrix A is an upper or */
/* >           lower triangular matrix as follows: */
/* > */
/* >              UPLO = 'U' or 'u'   A is an upper triangular matrix. */
/* > */
/* >              UPLO = 'L' or 'l'   A is a lower triangular matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] TRANSA */
/* > \verbatim */
/* >          TRANSA is CHARACTER*1 */
/* >           On entry, TRANSA specifies the form of op( A ) to be used in */
/* >           the matrix multiplication as follows: */
/* > */
/* >              TRANSA = 'N' or 'n'   op( A ) = A. */
/* > */
/* >              TRANSA = 'T' or 't'   op( A ) = A**T. */
/* > */
/* >              TRANSA = 'C' or 'c'   op( A ) = A**T. */
/* > \endverbatim */
/* > */
/* > \param[in] DIAG */
/* > \verbatim */
/* >          DIAG is CHARACTER*1 */
/* >           On entry, DIAG specifies whether or not A is unit triangular */
/* >           as follows: */
/* > */
/* >              DIAG = 'U' or 'u'   A is assumed to be unit triangular. */
/* > */
/* >              DIAG = 'N' or 'n'   A is not assumed to be unit */
/* >                                  triangular. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >           On entry, M specifies the number of rows of B. M must be at */
/* >           least zero. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >           On entry, N specifies the number of columns of B.  N must be */
/* >           at least zero. */
/* > \endverbatim */
/* > */
/* > \param[in] ALPHA */
/* > \verbatim */
/* >          ALPHA is DOUBLE PRECISION. */
/* >           On entry,  ALPHA specifies the scalar  alpha. When  alpha is */
/* >           zero then  A is not referenced and  B need not be set before */
/* >           entry. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array of DIMENSION ( LDA, k ), */
/* >           where k is m when SIDE = 'L' or 'l' */
/* >             and k is n when SIDE = 'R' or 'r'. */
/* >           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k */
/* >           upper triangular part of the array  A must contain the upper */
/* >           triangular matrix  and the strictly lower triangular part of */
/* >           A is not referenced. */
/* >           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k */
/* >           lower triangular part of the array  A must contain the lower */
/* >           triangular matrix  and the strictly upper triangular part of */
/* >           A is not referenced. */
/* >           Note that when  DIAG = 'U' or 'u',  the diagonal elements of */
/* >           A  are not referenced either,  but are assumed to be  unity. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >           On entry, LDA specifies the first dimension of A as declared */
/* >           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then */
/* >           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r' */
/* >           then LDA must be at least max( 1, n ). */
/* > \endverbatim */
/* > */
/* > \param[in,out] B */
/* > \verbatim */
/* >          B is DOUBLE PRECISION array of DIMENSION ( LDB, n ). */
/* >           Before entry,  the leading  m by n part of the array  B must */
/* >           contain  the  right-hand  side  matrix  B,  and  on exit  is */
/* >           overwritten by the solution matrix  X. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >           On entry, LDB specifies the first dimension of B as declared */
/* >           in  the  calling  (sub)  program.   LDB  must  be  at  least */
/* >           max( 1, m ). */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup double_blas_level3 */

/* > \par Further Details: */
/*  ===================== */
/* > */
/* > \verbatim */
/* > */
/* >  Level 3 Blas routine. */
/* > */
/* > */
/* >  -- Written on 8-February-1989. */
/* >     Jack Dongarra, Argonne National Laboratory. */
/* >     Iain Duff, AERE Harwell. */
/* >     Jeremy Du Croz, Numerical Algorithms Group Ltd. */
/* >     Sven Hammarling, Numerical Algorithms Group Ltd. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int dtrsm_(char *side, char *uplo, char *transa, char *diag, 
	integer *m, integer *n, doublereal *alpha, doublereal *a, integer *
	lda, doublereal *b, integer *ldb, ftnlen side_len, ftnlen uplo_len, 
	ftnlen transa_len, ftnlen diag_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, k, info;
    static doublereal temp;
    static logical lside;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer nrowa;
    static logical upper;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static logical nounit;


/*  -- Reference BLAS level3 routine (version 3.4.0) -- */
/*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Parameters .. */
/*     .. */

/*     Test the input parameters. */

#line 222 "dtrsm.f"
    /* Parameter adjustments */
#line 222 "dtrsm.f"
    a_dim1 = *lda;
#line 222 "dtrsm.f"
    a_offset = 1 + a_dim1;
#line 222 "dtrsm.f"
    a -= a_offset;
#line 222 "dtrsm.f"
    b_dim1 = *ldb;
#line 222 "dtrsm.f"
    b_offset = 1 + b_dim1;
#line 222 "dtrsm.f"
    b -= b_offset;
#line 222 "dtrsm.f"

#line 222 "dtrsm.f"
    /* Function Body */
#line 222 "dtrsm.f"
    lside = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 223 "dtrsm.f"
    if (lside) {
#line 224 "dtrsm.f"
	nrowa = *m;
#line 225 "dtrsm.f"
    } else {
#line 226 "dtrsm.f"
	nrowa = *n;
#line 227 "dtrsm.f"
    }
#line 228 "dtrsm.f"
    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);
#line 229 "dtrsm.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);

#line 231 "dtrsm.f"
    info = 0;
#line 232 "dtrsm.f"
    if (! lside && ! lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {
#line 233 "dtrsm.f"
	info = 1;
#line 234 "dtrsm.f"
    } else if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 235 "dtrsm.f"
	info = 2;
#line 236 "dtrsm.f"
    } else if (! lsame_(transa, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(transa,
	     "T", (ftnlen)1, (ftnlen)1) && ! lsame_(transa, "C", (ftnlen)1, (
	    ftnlen)1)) {
#line 239 "dtrsm.f"
	info = 3;
#line 240 "dtrsm.f"
    } else if (! lsame_(diag, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(diag, 
	    "N", (ftnlen)1, (ftnlen)1)) {
#line 241 "dtrsm.f"
	info = 4;
#line 242 "dtrsm.f"
    } else if (*m < 0) {
#line 243 "dtrsm.f"
	info = 5;
#line 244 "dtrsm.f"
    } else if (*n < 0) {
#line 245 "dtrsm.f"
	info = 6;
#line 246 "dtrsm.f"
    } else if (*lda < max(1,nrowa)) {
#line 247 "dtrsm.f"
	info = 9;
#line 248 "dtrsm.f"
    } else if (*ldb < max(1,*m)) {
#line 249 "dtrsm.f"
	info = 11;
#line 250 "dtrsm.f"
    }
#line 251 "dtrsm.f"
    if (info != 0) {
#line 252 "dtrsm.f"
	xerbla_("DTRSM ", &info, (ftnlen)6);
#line 253 "dtrsm.f"
	return 0;
#line 254 "dtrsm.f"
    }

/*     Quick return if possible. */

#line 258 "dtrsm.f"
    if (*m == 0 || *n == 0) {
#line 258 "dtrsm.f"
	return 0;
#line 258 "dtrsm.f"
    }

/*     And when  alpha.eq.zero. */

#line 262 "dtrsm.f"
    if (*alpha == 0.) {
#line 263 "dtrsm.f"
	i__1 = *n;
#line 263 "dtrsm.f"
	for (j = 1; j <= i__1; ++j) {
#line 264 "dtrsm.f"
	    i__2 = *m;
#line 264 "dtrsm.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 265 "dtrsm.f"
		b[i__ + j * b_dim1] = 0.;
#line 266 "dtrsm.f"
/* L10: */
#line 266 "dtrsm.f"
	    }
#line 267 "dtrsm.f"
/* L20: */
#line 267 "dtrsm.f"
	}
#line 268 "dtrsm.f"
	return 0;
#line 269 "dtrsm.f"
    }

/*     Start the operations. */

#line 273 "dtrsm.f"
    if (lside) {
#line 274 "dtrsm.f"
	if (lsame_(transa, "N", (ftnlen)1, (ftnlen)1)) {

/*           Form  B := alpha*inv( A )*B. */

#line 278 "dtrsm.f"
	    if (upper) {
#line 279 "dtrsm.f"
		i__1 = *n;
#line 279 "dtrsm.f"
		for (j = 1; j <= i__1; ++j) {
#line 280 "dtrsm.f"
		    if (*alpha != 1.) {
#line 281 "dtrsm.f"
			i__2 = *m;
#line 281 "dtrsm.f"
			for (i__ = 1; i__ <= i__2; ++i__) {
#line 282 "dtrsm.f"
			    b[i__ + j * b_dim1] = *alpha * b[i__ + j * b_dim1]
				    ;
#line 283 "dtrsm.f"
/* L30: */
#line 283 "dtrsm.f"
			}
#line 284 "dtrsm.f"
		    }
#line 285 "dtrsm.f"
		    for (k = *m; k >= 1; --k) {
#line 286 "dtrsm.f"
			if (b[k + j * b_dim1] != 0.) {
#line 287 "dtrsm.f"
			    if (nounit) {
#line 287 "dtrsm.f"
				b[k + j * b_dim1] /= a[k + k * a_dim1];
#line 287 "dtrsm.f"
			    }
#line 288 "dtrsm.f"
			    i__2 = k - 1;
#line 288 "dtrsm.f"
			    for (i__ = 1; i__ <= i__2; ++i__) {
#line 289 "dtrsm.f"
				b[i__ + j * b_dim1] -= b[k + j * b_dim1] * a[
					i__ + k * a_dim1];
#line 290 "dtrsm.f"
/* L40: */
#line 290 "dtrsm.f"
			    }
#line 291 "dtrsm.f"
			}
#line 292 "dtrsm.f"
/* L50: */
#line 292 "dtrsm.f"
		    }
#line 293 "dtrsm.f"
/* L60: */
#line 293 "dtrsm.f"
		}
#line 294 "dtrsm.f"
	    } else {
#line 295 "dtrsm.f"
		i__1 = *n;
#line 295 "dtrsm.f"
		for (j = 1; j <= i__1; ++j) {
#line 296 "dtrsm.f"
		    if (*alpha != 1.) {
#line 297 "dtrsm.f"
			i__2 = *m;
#line 297 "dtrsm.f"
			for (i__ = 1; i__ <= i__2; ++i__) {
#line 298 "dtrsm.f"
			    b[i__ + j * b_dim1] = *alpha * b[i__ + j * b_dim1]
				    ;
#line 299 "dtrsm.f"
/* L70: */
#line 299 "dtrsm.f"
			}
#line 300 "dtrsm.f"
		    }
#line 301 "dtrsm.f"
		    i__2 = *m;
#line 301 "dtrsm.f"
		    for (k = 1; k <= i__2; ++k) {
#line 302 "dtrsm.f"
			if (b[k + j * b_dim1] != 0.) {
#line 303 "dtrsm.f"
			    if (nounit) {
#line 303 "dtrsm.f"
				b[k + j * b_dim1] /= a[k + k * a_dim1];
#line 303 "dtrsm.f"
			    }
#line 304 "dtrsm.f"
			    i__3 = *m;
#line 304 "dtrsm.f"
			    for (i__ = k + 1; i__ <= i__3; ++i__) {
#line 305 "dtrsm.f"
				b[i__ + j * b_dim1] -= b[k + j * b_dim1] * a[
					i__ + k * a_dim1];
#line 306 "dtrsm.f"
/* L80: */
#line 306 "dtrsm.f"
			    }
#line 307 "dtrsm.f"
			}
#line 308 "dtrsm.f"
/* L90: */
#line 308 "dtrsm.f"
		    }
#line 309 "dtrsm.f"
/* L100: */
#line 309 "dtrsm.f"
		}
#line 310 "dtrsm.f"
	    }
#line 311 "dtrsm.f"
	} else {

/*           Form  B := alpha*inv( A**T )*B. */

#line 315 "dtrsm.f"
	    if (upper) {
#line 316 "dtrsm.f"
		i__1 = *n;
#line 316 "dtrsm.f"
		for (j = 1; j <= i__1; ++j) {
#line 317 "dtrsm.f"
		    i__2 = *m;
#line 317 "dtrsm.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 318 "dtrsm.f"
			temp = *alpha * b[i__ + j * b_dim1];
#line 319 "dtrsm.f"
			i__3 = i__ - 1;
#line 319 "dtrsm.f"
			for (k = 1; k <= i__3; ++k) {
#line 320 "dtrsm.f"
			    temp -= a[k + i__ * a_dim1] * b[k + j * b_dim1];
#line 321 "dtrsm.f"
/* L110: */
#line 321 "dtrsm.f"
			}
#line 322 "dtrsm.f"
			if (nounit) {
#line 322 "dtrsm.f"
			    temp /= a[i__ + i__ * a_dim1];
#line 322 "dtrsm.f"
			}
#line 323 "dtrsm.f"
			b[i__ + j * b_dim1] = temp;
#line 324 "dtrsm.f"
/* L120: */
#line 324 "dtrsm.f"
		    }
#line 325 "dtrsm.f"
/* L130: */
#line 325 "dtrsm.f"
		}
#line 326 "dtrsm.f"
	    } else {
#line 327 "dtrsm.f"
		i__1 = *n;
#line 327 "dtrsm.f"
		for (j = 1; j <= i__1; ++j) {
#line 328 "dtrsm.f"
		    for (i__ = *m; i__ >= 1; --i__) {
#line 329 "dtrsm.f"
			temp = *alpha * b[i__ + j * b_dim1];
#line 330 "dtrsm.f"
			i__2 = *m;
#line 330 "dtrsm.f"
			for (k = i__ + 1; k <= i__2; ++k) {
#line 331 "dtrsm.f"
			    temp -= a[k + i__ * a_dim1] * b[k + j * b_dim1];
#line 332 "dtrsm.f"
/* L140: */
#line 332 "dtrsm.f"
			}
#line 333 "dtrsm.f"
			if (nounit) {
#line 333 "dtrsm.f"
			    temp /= a[i__ + i__ * a_dim1];
#line 333 "dtrsm.f"
			}
#line 334 "dtrsm.f"
			b[i__ + j * b_dim1] = temp;
#line 335 "dtrsm.f"
/* L150: */
#line 335 "dtrsm.f"
		    }
#line 336 "dtrsm.f"
/* L160: */
#line 336 "dtrsm.f"
		}
#line 337 "dtrsm.f"
	    }
#line 338 "dtrsm.f"
	}
#line 339 "dtrsm.f"
    } else {
#line 340 "dtrsm.f"
	if (lsame_(transa, "N", (ftnlen)1, (ftnlen)1)) {

/*           Form  B := alpha*B*inv( A ). */

#line 344 "dtrsm.f"
	    if (upper) {
#line 345 "dtrsm.f"
		i__1 = *n;
#line 345 "dtrsm.f"
		for (j = 1; j <= i__1; ++j) {
#line 346 "dtrsm.f"
		    if (*alpha != 1.) {
#line 347 "dtrsm.f"
			i__2 = *m;
#line 347 "dtrsm.f"
			for (i__ = 1; i__ <= i__2; ++i__) {
#line 348 "dtrsm.f"
			    b[i__ + j * b_dim1] = *alpha * b[i__ + j * b_dim1]
				    ;
#line 349 "dtrsm.f"
/* L170: */
#line 349 "dtrsm.f"
			}
#line 350 "dtrsm.f"
		    }
#line 351 "dtrsm.f"
		    i__2 = j - 1;
#line 351 "dtrsm.f"
		    for (k = 1; k <= i__2; ++k) {
#line 352 "dtrsm.f"
			if (a[k + j * a_dim1] != 0.) {
#line 353 "dtrsm.f"
			    i__3 = *m;
#line 353 "dtrsm.f"
			    for (i__ = 1; i__ <= i__3; ++i__) {
#line 354 "dtrsm.f"
				b[i__ + j * b_dim1] -= a[k + j * a_dim1] * b[
					i__ + k * b_dim1];
#line 355 "dtrsm.f"
/* L180: */
#line 355 "dtrsm.f"
			    }
#line 356 "dtrsm.f"
			}
#line 357 "dtrsm.f"
/* L190: */
#line 357 "dtrsm.f"
		    }
#line 358 "dtrsm.f"
		    if (nounit) {
#line 359 "dtrsm.f"
			temp = 1. / a[j + j * a_dim1];
#line 360 "dtrsm.f"
			i__2 = *m;
#line 360 "dtrsm.f"
			for (i__ = 1; i__ <= i__2; ++i__) {
#line 361 "dtrsm.f"
			    b[i__ + j * b_dim1] = temp * b[i__ + j * b_dim1];
#line 362 "dtrsm.f"
/* L200: */
#line 362 "dtrsm.f"
			}
#line 363 "dtrsm.f"
		    }
#line 364 "dtrsm.f"
/* L210: */
#line 364 "dtrsm.f"
		}
#line 365 "dtrsm.f"
	    } else {
#line 366 "dtrsm.f"
		for (j = *n; j >= 1; --j) {
#line 367 "dtrsm.f"
		    if (*alpha != 1.) {
#line 368 "dtrsm.f"
			i__1 = *m;
#line 368 "dtrsm.f"
			for (i__ = 1; i__ <= i__1; ++i__) {
#line 369 "dtrsm.f"
			    b[i__ + j * b_dim1] = *alpha * b[i__ + j * b_dim1]
				    ;
#line 370 "dtrsm.f"
/* L220: */
#line 370 "dtrsm.f"
			}
#line 371 "dtrsm.f"
		    }
#line 372 "dtrsm.f"
		    i__1 = *n;
#line 372 "dtrsm.f"
		    for (k = j + 1; k <= i__1; ++k) {
#line 373 "dtrsm.f"
			if (a[k + j * a_dim1] != 0.) {
#line 374 "dtrsm.f"
			    i__2 = *m;
#line 374 "dtrsm.f"
			    for (i__ = 1; i__ <= i__2; ++i__) {
#line 375 "dtrsm.f"
				b[i__ + j * b_dim1] -= a[k + j * a_dim1] * b[
					i__ + k * b_dim1];
#line 376 "dtrsm.f"
/* L230: */
#line 376 "dtrsm.f"
			    }
#line 377 "dtrsm.f"
			}
#line 378 "dtrsm.f"
/* L240: */
#line 378 "dtrsm.f"
		    }
#line 379 "dtrsm.f"
		    if (nounit) {
#line 380 "dtrsm.f"
			temp = 1. / a[j + j * a_dim1];
#line 381 "dtrsm.f"
			i__1 = *m;
#line 381 "dtrsm.f"
			for (i__ = 1; i__ <= i__1; ++i__) {
#line 382 "dtrsm.f"
			    b[i__ + j * b_dim1] = temp * b[i__ + j * b_dim1];
#line 383 "dtrsm.f"
/* L250: */
#line 383 "dtrsm.f"
			}
#line 384 "dtrsm.f"
		    }
#line 385 "dtrsm.f"
/* L260: */
#line 385 "dtrsm.f"
		}
#line 386 "dtrsm.f"
	    }
#line 387 "dtrsm.f"
	} else {

/*           Form  B := alpha*B*inv( A**T ). */

#line 391 "dtrsm.f"
	    if (upper) {
#line 392 "dtrsm.f"
		for (k = *n; k >= 1; --k) {
#line 393 "dtrsm.f"
		    if (nounit) {
#line 394 "dtrsm.f"
			temp = 1. / a[k + k * a_dim1];
#line 395 "dtrsm.f"
			i__1 = *m;
#line 395 "dtrsm.f"
			for (i__ = 1; i__ <= i__1; ++i__) {
#line 396 "dtrsm.f"
			    b[i__ + k * b_dim1] = temp * b[i__ + k * b_dim1];
#line 397 "dtrsm.f"
/* L270: */
#line 397 "dtrsm.f"
			}
#line 398 "dtrsm.f"
		    }
#line 399 "dtrsm.f"
		    i__1 = k - 1;
#line 399 "dtrsm.f"
		    for (j = 1; j <= i__1; ++j) {
#line 400 "dtrsm.f"
			if (a[j + k * a_dim1] != 0.) {
#line 401 "dtrsm.f"
			    temp = a[j + k * a_dim1];
#line 402 "dtrsm.f"
			    i__2 = *m;
#line 402 "dtrsm.f"
			    for (i__ = 1; i__ <= i__2; ++i__) {
#line 403 "dtrsm.f"
				b[i__ + j * b_dim1] -= temp * b[i__ + k * 
					b_dim1];
#line 404 "dtrsm.f"
/* L280: */
#line 404 "dtrsm.f"
			    }
#line 405 "dtrsm.f"
			}
#line 406 "dtrsm.f"
/* L290: */
#line 406 "dtrsm.f"
		    }
#line 407 "dtrsm.f"
		    if (*alpha != 1.) {
#line 408 "dtrsm.f"
			i__1 = *m;
#line 408 "dtrsm.f"
			for (i__ = 1; i__ <= i__1; ++i__) {
#line 409 "dtrsm.f"
			    b[i__ + k * b_dim1] = *alpha * b[i__ + k * b_dim1]
				    ;
#line 410 "dtrsm.f"
/* L300: */
#line 410 "dtrsm.f"
			}
#line 411 "dtrsm.f"
		    }
#line 412 "dtrsm.f"
/* L310: */
#line 412 "dtrsm.f"
		}
#line 413 "dtrsm.f"
	    } else {
#line 414 "dtrsm.f"
		i__1 = *n;
#line 414 "dtrsm.f"
		for (k = 1; k <= i__1; ++k) {
#line 415 "dtrsm.f"
		    if (nounit) {
#line 416 "dtrsm.f"
			temp = 1. / a[k + k * a_dim1];
#line 417 "dtrsm.f"
			i__2 = *m;
#line 417 "dtrsm.f"
			for (i__ = 1; i__ <= i__2; ++i__) {
#line 418 "dtrsm.f"
			    b[i__ + k * b_dim1] = temp * b[i__ + k * b_dim1];
#line 419 "dtrsm.f"
/* L320: */
#line 419 "dtrsm.f"
			}
#line 420 "dtrsm.f"
		    }
#line 421 "dtrsm.f"
		    i__2 = *n;
#line 421 "dtrsm.f"
		    for (j = k + 1; j <= i__2; ++j) {
#line 422 "dtrsm.f"
			if (a[j + k * a_dim1] != 0.) {
#line 423 "dtrsm.f"
			    temp = a[j + k * a_dim1];
#line 424 "dtrsm.f"
			    i__3 = *m;
#line 424 "dtrsm.f"
			    for (i__ = 1; i__ <= i__3; ++i__) {
#line 425 "dtrsm.f"
				b[i__ + j * b_dim1] -= temp * b[i__ + k * 
					b_dim1];
#line 426 "dtrsm.f"
/* L330: */
#line 426 "dtrsm.f"
			    }
#line 427 "dtrsm.f"
			}
#line 428 "dtrsm.f"
/* L340: */
#line 428 "dtrsm.f"
		    }
#line 429 "dtrsm.f"
		    if (*alpha != 1.) {
#line 430 "dtrsm.f"
			i__2 = *m;
#line 430 "dtrsm.f"
			for (i__ = 1; i__ <= i__2; ++i__) {
#line 431 "dtrsm.f"
			    b[i__ + k * b_dim1] = *alpha * b[i__ + k * b_dim1]
				    ;
#line 432 "dtrsm.f"
/* L350: */
#line 432 "dtrsm.f"
			}
#line 433 "dtrsm.f"
		    }
#line 434 "dtrsm.f"
/* L360: */
#line 434 "dtrsm.f"
		}
#line 435 "dtrsm.f"
	    }
#line 436 "dtrsm.f"
	}
#line 437 "dtrsm.f"
    }

#line 439 "dtrsm.f"
    return 0;

/*     End of DTRSM . */

} /* dtrsm_ */

