#line 1 "strsm.f"
/* strsm.f -- translated by f2c (version 20100827).
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

#line 1 "strsm.f"
/* > \brief \b STRSM */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE STRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB) */

/*       .. Scalar Arguments .. */
/*       REAL ALPHA */
/*       INTEGER LDA,LDB,M,N */
/*       CHARACTER DIAG,SIDE,TRANSA,UPLO */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL A(LDA,*),B(LDB,*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > STRSM  solves one of the matrix equations */
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
/* >          ALPHA is REAL */
/* >           On entry,  ALPHA specifies the scalar  alpha. When  alpha is */
/* >           zero then  A is not referenced and  B need not be set before */
/* >           entry. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is REAL array of DIMENSION ( LDA, k ), */
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
/* >          B is REAL array of DIMENSION ( LDB, n ). */
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

/* > \date December 2016 */

/* > \ingroup single_blas_level3 */

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
/* Subroutine */ int strsm_(char *side, char *uplo, char *transa, char *diag, 
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


/*  -- Reference BLAS level3 routine (version 3.7.0) -- */
/*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

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

#line 222 "strsm.f"
    /* Parameter adjustments */
#line 222 "strsm.f"
    a_dim1 = *lda;
#line 222 "strsm.f"
    a_offset = 1 + a_dim1;
#line 222 "strsm.f"
    a -= a_offset;
#line 222 "strsm.f"
    b_dim1 = *ldb;
#line 222 "strsm.f"
    b_offset = 1 + b_dim1;
#line 222 "strsm.f"
    b -= b_offset;
#line 222 "strsm.f"

#line 222 "strsm.f"
    /* Function Body */
#line 222 "strsm.f"
    lside = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 223 "strsm.f"
    if (lside) {
#line 224 "strsm.f"
	nrowa = *m;
#line 225 "strsm.f"
    } else {
#line 226 "strsm.f"
	nrowa = *n;
#line 227 "strsm.f"
    }
#line 228 "strsm.f"
    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);
#line 229 "strsm.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);

#line 231 "strsm.f"
    info = 0;
#line 232 "strsm.f"
    if (! lside && ! lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {
#line 233 "strsm.f"
	info = 1;
#line 234 "strsm.f"
    } else if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 235 "strsm.f"
	info = 2;
#line 236 "strsm.f"
    } else if (! lsame_(transa, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(transa,
	     "T", (ftnlen)1, (ftnlen)1) && ! lsame_(transa, "C", (ftnlen)1, (
	    ftnlen)1)) {
#line 239 "strsm.f"
	info = 3;
#line 240 "strsm.f"
    } else if (! lsame_(diag, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(diag, 
	    "N", (ftnlen)1, (ftnlen)1)) {
#line 241 "strsm.f"
	info = 4;
#line 242 "strsm.f"
    } else if (*m < 0) {
#line 243 "strsm.f"
	info = 5;
#line 244 "strsm.f"
    } else if (*n < 0) {
#line 245 "strsm.f"
	info = 6;
#line 246 "strsm.f"
    } else if (*lda < max(1,nrowa)) {
#line 247 "strsm.f"
	info = 9;
#line 248 "strsm.f"
    } else if (*ldb < max(1,*m)) {
#line 249 "strsm.f"
	info = 11;
#line 250 "strsm.f"
    }
#line 251 "strsm.f"
    if (info != 0) {
#line 252 "strsm.f"
	xerbla_("STRSM ", &info, (ftnlen)6);
#line 253 "strsm.f"
	return 0;
#line 254 "strsm.f"
    }

/*     Quick return if possible. */

#line 258 "strsm.f"
    if (*m == 0 || *n == 0) {
#line 258 "strsm.f"
	return 0;
#line 258 "strsm.f"
    }

/*     And when  alpha.eq.zero. */

#line 262 "strsm.f"
    if (*alpha == 0.) {
#line 263 "strsm.f"
	i__1 = *n;
#line 263 "strsm.f"
	for (j = 1; j <= i__1; ++j) {
#line 264 "strsm.f"
	    i__2 = *m;
#line 264 "strsm.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 265 "strsm.f"
		b[i__ + j * b_dim1] = 0.;
#line 266 "strsm.f"
/* L10: */
#line 266 "strsm.f"
	    }
#line 267 "strsm.f"
/* L20: */
#line 267 "strsm.f"
	}
#line 268 "strsm.f"
	return 0;
#line 269 "strsm.f"
    }

/*     Start the operations. */

#line 273 "strsm.f"
    if (lside) {
#line 274 "strsm.f"
	if (lsame_(transa, "N", (ftnlen)1, (ftnlen)1)) {

/*           Form  B := alpha*inv( A )*B. */

#line 278 "strsm.f"
	    if (upper) {
#line 279 "strsm.f"
		i__1 = *n;
#line 279 "strsm.f"
		for (j = 1; j <= i__1; ++j) {
#line 280 "strsm.f"
		    if (*alpha != 1.) {
#line 281 "strsm.f"
			i__2 = *m;
#line 281 "strsm.f"
			for (i__ = 1; i__ <= i__2; ++i__) {
#line 282 "strsm.f"
			    b[i__ + j * b_dim1] = *alpha * b[i__ + j * b_dim1]
				    ;
#line 283 "strsm.f"
/* L30: */
#line 283 "strsm.f"
			}
#line 284 "strsm.f"
		    }
#line 285 "strsm.f"
		    for (k = *m; k >= 1; --k) {
#line 286 "strsm.f"
			if (b[k + j * b_dim1] != 0.) {
#line 287 "strsm.f"
			    if (nounit) {
#line 287 "strsm.f"
				b[k + j * b_dim1] /= a[k + k * a_dim1];
#line 287 "strsm.f"
			    }
#line 288 "strsm.f"
			    i__2 = k - 1;
#line 288 "strsm.f"
			    for (i__ = 1; i__ <= i__2; ++i__) {
#line 289 "strsm.f"
				b[i__ + j * b_dim1] -= b[k + j * b_dim1] * a[
					i__ + k * a_dim1];
#line 290 "strsm.f"
/* L40: */
#line 290 "strsm.f"
			    }
#line 291 "strsm.f"
			}
#line 292 "strsm.f"
/* L50: */
#line 292 "strsm.f"
		    }
#line 293 "strsm.f"
/* L60: */
#line 293 "strsm.f"
		}
#line 294 "strsm.f"
	    } else {
#line 295 "strsm.f"
		i__1 = *n;
#line 295 "strsm.f"
		for (j = 1; j <= i__1; ++j) {
#line 296 "strsm.f"
		    if (*alpha != 1.) {
#line 297 "strsm.f"
			i__2 = *m;
#line 297 "strsm.f"
			for (i__ = 1; i__ <= i__2; ++i__) {
#line 298 "strsm.f"
			    b[i__ + j * b_dim1] = *alpha * b[i__ + j * b_dim1]
				    ;
#line 299 "strsm.f"
/* L70: */
#line 299 "strsm.f"
			}
#line 300 "strsm.f"
		    }
#line 301 "strsm.f"
		    i__2 = *m;
#line 301 "strsm.f"
		    for (k = 1; k <= i__2; ++k) {
#line 302 "strsm.f"
			if (b[k + j * b_dim1] != 0.) {
#line 303 "strsm.f"
			    if (nounit) {
#line 303 "strsm.f"
				b[k + j * b_dim1] /= a[k + k * a_dim1];
#line 303 "strsm.f"
			    }
#line 304 "strsm.f"
			    i__3 = *m;
#line 304 "strsm.f"
			    for (i__ = k + 1; i__ <= i__3; ++i__) {
#line 305 "strsm.f"
				b[i__ + j * b_dim1] -= b[k + j * b_dim1] * a[
					i__ + k * a_dim1];
#line 306 "strsm.f"
/* L80: */
#line 306 "strsm.f"
			    }
#line 307 "strsm.f"
			}
#line 308 "strsm.f"
/* L90: */
#line 308 "strsm.f"
		    }
#line 309 "strsm.f"
/* L100: */
#line 309 "strsm.f"
		}
#line 310 "strsm.f"
	    }
#line 311 "strsm.f"
	} else {

/*           Form  B := alpha*inv( A**T )*B. */

#line 315 "strsm.f"
	    if (upper) {
#line 316 "strsm.f"
		i__1 = *n;
#line 316 "strsm.f"
		for (j = 1; j <= i__1; ++j) {
#line 317 "strsm.f"
		    i__2 = *m;
#line 317 "strsm.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 318 "strsm.f"
			temp = *alpha * b[i__ + j * b_dim1];
#line 319 "strsm.f"
			i__3 = i__ - 1;
#line 319 "strsm.f"
			for (k = 1; k <= i__3; ++k) {
#line 320 "strsm.f"
			    temp -= a[k + i__ * a_dim1] * b[k + j * b_dim1];
#line 321 "strsm.f"
/* L110: */
#line 321 "strsm.f"
			}
#line 322 "strsm.f"
			if (nounit) {
#line 322 "strsm.f"
			    temp /= a[i__ + i__ * a_dim1];
#line 322 "strsm.f"
			}
#line 323 "strsm.f"
			b[i__ + j * b_dim1] = temp;
#line 324 "strsm.f"
/* L120: */
#line 324 "strsm.f"
		    }
#line 325 "strsm.f"
/* L130: */
#line 325 "strsm.f"
		}
#line 326 "strsm.f"
	    } else {
#line 327 "strsm.f"
		i__1 = *n;
#line 327 "strsm.f"
		for (j = 1; j <= i__1; ++j) {
#line 328 "strsm.f"
		    for (i__ = *m; i__ >= 1; --i__) {
#line 329 "strsm.f"
			temp = *alpha * b[i__ + j * b_dim1];
#line 330 "strsm.f"
			i__2 = *m;
#line 330 "strsm.f"
			for (k = i__ + 1; k <= i__2; ++k) {
#line 331 "strsm.f"
			    temp -= a[k + i__ * a_dim1] * b[k + j * b_dim1];
#line 332 "strsm.f"
/* L140: */
#line 332 "strsm.f"
			}
#line 333 "strsm.f"
			if (nounit) {
#line 333 "strsm.f"
			    temp /= a[i__ + i__ * a_dim1];
#line 333 "strsm.f"
			}
#line 334 "strsm.f"
			b[i__ + j * b_dim1] = temp;
#line 335 "strsm.f"
/* L150: */
#line 335 "strsm.f"
		    }
#line 336 "strsm.f"
/* L160: */
#line 336 "strsm.f"
		}
#line 337 "strsm.f"
	    }
#line 338 "strsm.f"
	}
#line 339 "strsm.f"
    } else {
#line 340 "strsm.f"
	if (lsame_(transa, "N", (ftnlen)1, (ftnlen)1)) {

/*           Form  B := alpha*B*inv( A ). */

#line 344 "strsm.f"
	    if (upper) {
#line 345 "strsm.f"
		i__1 = *n;
#line 345 "strsm.f"
		for (j = 1; j <= i__1; ++j) {
#line 346 "strsm.f"
		    if (*alpha != 1.) {
#line 347 "strsm.f"
			i__2 = *m;
#line 347 "strsm.f"
			for (i__ = 1; i__ <= i__2; ++i__) {
#line 348 "strsm.f"
			    b[i__ + j * b_dim1] = *alpha * b[i__ + j * b_dim1]
				    ;
#line 349 "strsm.f"
/* L170: */
#line 349 "strsm.f"
			}
#line 350 "strsm.f"
		    }
#line 351 "strsm.f"
		    i__2 = j - 1;
#line 351 "strsm.f"
		    for (k = 1; k <= i__2; ++k) {
#line 352 "strsm.f"
			if (a[k + j * a_dim1] != 0.) {
#line 353 "strsm.f"
			    i__3 = *m;
#line 353 "strsm.f"
			    for (i__ = 1; i__ <= i__3; ++i__) {
#line 354 "strsm.f"
				b[i__ + j * b_dim1] -= a[k + j * a_dim1] * b[
					i__ + k * b_dim1];
#line 355 "strsm.f"
/* L180: */
#line 355 "strsm.f"
			    }
#line 356 "strsm.f"
			}
#line 357 "strsm.f"
/* L190: */
#line 357 "strsm.f"
		    }
#line 358 "strsm.f"
		    if (nounit) {
#line 359 "strsm.f"
			temp = 1. / a[j + j * a_dim1];
#line 360 "strsm.f"
			i__2 = *m;
#line 360 "strsm.f"
			for (i__ = 1; i__ <= i__2; ++i__) {
#line 361 "strsm.f"
			    b[i__ + j * b_dim1] = temp * b[i__ + j * b_dim1];
#line 362 "strsm.f"
/* L200: */
#line 362 "strsm.f"
			}
#line 363 "strsm.f"
		    }
#line 364 "strsm.f"
/* L210: */
#line 364 "strsm.f"
		}
#line 365 "strsm.f"
	    } else {
#line 366 "strsm.f"
		for (j = *n; j >= 1; --j) {
#line 367 "strsm.f"
		    if (*alpha != 1.) {
#line 368 "strsm.f"
			i__1 = *m;
#line 368 "strsm.f"
			for (i__ = 1; i__ <= i__1; ++i__) {
#line 369 "strsm.f"
			    b[i__ + j * b_dim1] = *alpha * b[i__ + j * b_dim1]
				    ;
#line 370 "strsm.f"
/* L220: */
#line 370 "strsm.f"
			}
#line 371 "strsm.f"
		    }
#line 372 "strsm.f"
		    i__1 = *n;
#line 372 "strsm.f"
		    for (k = j + 1; k <= i__1; ++k) {
#line 373 "strsm.f"
			if (a[k + j * a_dim1] != 0.) {
#line 374 "strsm.f"
			    i__2 = *m;
#line 374 "strsm.f"
			    for (i__ = 1; i__ <= i__2; ++i__) {
#line 375 "strsm.f"
				b[i__ + j * b_dim1] -= a[k + j * a_dim1] * b[
					i__ + k * b_dim1];
#line 376 "strsm.f"
/* L230: */
#line 376 "strsm.f"
			    }
#line 377 "strsm.f"
			}
#line 378 "strsm.f"
/* L240: */
#line 378 "strsm.f"
		    }
#line 379 "strsm.f"
		    if (nounit) {
#line 380 "strsm.f"
			temp = 1. / a[j + j * a_dim1];
#line 381 "strsm.f"
			i__1 = *m;
#line 381 "strsm.f"
			for (i__ = 1; i__ <= i__1; ++i__) {
#line 382 "strsm.f"
			    b[i__ + j * b_dim1] = temp * b[i__ + j * b_dim1];
#line 383 "strsm.f"
/* L250: */
#line 383 "strsm.f"
			}
#line 384 "strsm.f"
		    }
#line 385 "strsm.f"
/* L260: */
#line 385 "strsm.f"
		}
#line 386 "strsm.f"
	    }
#line 387 "strsm.f"
	} else {

/*           Form  B := alpha*B*inv( A**T ). */

#line 391 "strsm.f"
	    if (upper) {
#line 392 "strsm.f"
		for (k = *n; k >= 1; --k) {
#line 393 "strsm.f"
		    if (nounit) {
#line 394 "strsm.f"
			temp = 1. / a[k + k * a_dim1];
#line 395 "strsm.f"
			i__1 = *m;
#line 395 "strsm.f"
			for (i__ = 1; i__ <= i__1; ++i__) {
#line 396 "strsm.f"
			    b[i__ + k * b_dim1] = temp * b[i__ + k * b_dim1];
#line 397 "strsm.f"
/* L270: */
#line 397 "strsm.f"
			}
#line 398 "strsm.f"
		    }
#line 399 "strsm.f"
		    i__1 = k - 1;
#line 399 "strsm.f"
		    for (j = 1; j <= i__1; ++j) {
#line 400 "strsm.f"
			if (a[j + k * a_dim1] != 0.) {
#line 401 "strsm.f"
			    temp = a[j + k * a_dim1];
#line 402 "strsm.f"
			    i__2 = *m;
#line 402 "strsm.f"
			    for (i__ = 1; i__ <= i__2; ++i__) {
#line 403 "strsm.f"
				b[i__ + j * b_dim1] -= temp * b[i__ + k * 
					b_dim1];
#line 404 "strsm.f"
/* L280: */
#line 404 "strsm.f"
			    }
#line 405 "strsm.f"
			}
#line 406 "strsm.f"
/* L290: */
#line 406 "strsm.f"
		    }
#line 407 "strsm.f"
		    if (*alpha != 1.) {
#line 408 "strsm.f"
			i__1 = *m;
#line 408 "strsm.f"
			for (i__ = 1; i__ <= i__1; ++i__) {
#line 409 "strsm.f"
			    b[i__ + k * b_dim1] = *alpha * b[i__ + k * b_dim1]
				    ;
#line 410 "strsm.f"
/* L300: */
#line 410 "strsm.f"
			}
#line 411 "strsm.f"
		    }
#line 412 "strsm.f"
/* L310: */
#line 412 "strsm.f"
		}
#line 413 "strsm.f"
	    } else {
#line 414 "strsm.f"
		i__1 = *n;
#line 414 "strsm.f"
		for (k = 1; k <= i__1; ++k) {
#line 415 "strsm.f"
		    if (nounit) {
#line 416 "strsm.f"
			temp = 1. / a[k + k * a_dim1];
#line 417 "strsm.f"
			i__2 = *m;
#line 417 "strsm.f"
			for (i__ = 1; i__ <= i__2; ++i__) {
#line 418 "strsm.f"
			    b[i__ + k * b_dim1] = temp * b[i__ + k * b_dim1];
#line 419 "strsm.f"
/* L320: */
#line 419 "strsm.f"
			}
#line 420 "strsm.f"
		    }
#line 421 "strsm.f"
		    i__2 = *n;
#line 421 "strsm.f"
		    for (j = k + 1; j <= i__2; ++j) {
#line 422 "strsm.f"
			if (a[j + k * a_dim1] != 0.) {
#line 423 "strsm.f"
			    temp = a[j + k * a_dim1];
#line 424 "strsm.f"
			    i__3 = *m;
#line 424 "strsm.f"
			    for (i__ = 1; i__ <= i__3; ++i__) {
#line 425 "strsm.f"
				b[i__ + j * b_dim1] -= temp * b[i__ + k * 
					b_dim1];
#line 426 "strsm.f"
/* L330: */
#line 426 "strsm.f"
			    }
#line 427 "strsm.f"
			}
#line 428 "strsm.f"
/* L340: */
#line 428 "strsm.f"
		    }
#line 429 "strsm.f"
		    if (*alpha != 1.) {
#line 430 "strsm.f"
			i__2 = *m;
#line 430 "strsm.f"
			for (i__ = 1; i__ <= i__2; ++i__) {
#line 431 "strsm.f"
			    b[i__ + k * b_dim1] = *alpha * b[i__ + k * b_dim1]
				    ;
#line 432 "strsm.f"
/* L350: */
#line 432 "strsm.f"
			}
#line 433 "strsm.f"
		    }
#line 434 "strsm.f"
/* L360: */
#line 434 "strsm.f"
		}
#line 435 "strsm.f"
	    }
#line 436 "strsm.f"
	}
#line 437 "strsm.f"
    }

#line 439 "strsm.f"
    return 0;

/*     End of STRSM . */

} /* strsm_ */

