#line 1 "strmm.f"
/* strmm.f -- translated by f2c (version 20100827).
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

#line 1 "strmm.f"
/* > \brief \b STRMM */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE STRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB) */

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
/* > STRMM  performs one of the matrix-matrix operations */
/* > */
/* >    B := alpha*op( A )*B,   or   B := alpha*B*op( A ), */
/* > */
/* > where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or */
/* > non-unit,  upper or lower triangular matrix  and  op( A )  is one  of */
/* > */
/* >    op( A ) = A   or   op( A ) = A**T. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] SIDE */
/* > \verbatim */
/* >          SIDE is CHARACTER*1 */
/* >           On entry,  SIDE specifies whether  op( A ) multiplies B from */
/* >           the left or right as follows: */
/* > */
/* >              SIDE = 'L' or 'l'   B := alpha*op( A )*B. */
/* > */
/* >              SIDE = 'R' or 'r'   B := alpha*B*op( A ). */
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
/* >          A is REAL array, dimension ( LDA, k ), where k is m */
/* >           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'. */
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
/* >          B is REAL array, dimension ( LDB, N ) */
/* >           Before entry,  the leading  m by n part of the array  B must */
/* >           contain the matrix  B,  and  on exit  is overwritten  by the */
/* >           transformed matrix. */
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
/* >  -- Written on 8-February-1989. */
/* >     Jack Dongarra, Argonne National Laboratory. */
/* >     Iain Duff, AERE Harwell. */
/* >     Jeremy Du Croz, Numerical Algorithms Group Ltd. */
/* >     Sven Hammarling, Numerical Algorithms Group Ltd. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int strmm_(char *side, char *uplo, char *transa, char *diag, 
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

#line 218 "strmm.f"
    /* Parameter adjustments */
#line 218 "strmm.f"
    a_dim1 = *lda;
#line 218 "strmm.f"
    a_offset = 1 + a_dim1;
#line 218 "strmm.f"
    a -= a_offset;
#line 218 "strmm.f"
    b_dim1 = *ldb;
#line 218 "strmm.f"
    b_offset = 1 + b_dim1;
#line 218 "strmm.f"
    b -= b_offset;
#line 218 "strmm.f"

#line 218 "strmm.f"
    /* Function Body */
#line 218 "strmm.f"
    lside = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 219 "strmm.f"
    if (lside) {
#line 220 "strmm.f"
	nrowa = *m;
#line 221 "strmm.f"
    } else {
#line 222 "strmm.f"
	nrowa = *n;
#line 223 "strmm.f"
    }
#line 224 "strmm.f"
    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);
#line 225 "strmm.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);

#line 227 "strmm.f"
    info = 0;
#line 228 "strmm.f"
    if (! lside && ! lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {
#line 229 "strmm.f"
	info = 1;
#line 230 "strmm.f"
    } else if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 231 "strmm.f"
	info = 2;
#line 232 "strmm.f"
    } else if (! lsame_(transa, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(transa,
	     "T", (ftnlen)1, (ftnlen)1) && ! lsame_(transa, "C", (ftnlen)1, (
	    ftnlen)1)) {
#line 235 "strmm.f"
	info = 3;
#line 236 "strmm.f"
    } else if (! lsame_(diag, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(diag, 
	    "N", (ftnlen)1, (ftnlen)1)) {
#line 237 "strmm.f"
	info = 4;
#line 238 "strmm.f"
    } else if (*m < 0) {
#line 239 "strmm.f"
	info = 5;
#line 240 "strmm.f"
    } else if (*n < 0) {
#line 241 "strmm.f"
	info = 6;
#line 242 "strmm.f"
    } else if (*lda < max(1,nrowa)) {
#line 243 "strmm.f"
	info = 9;
#line 244 "strmm.f"
    } else if (*ldb < max(1,*m)) {
#line 245 "strmm.f"
	info = 11;
#line 246 "strmm.f"
    }
#line 247 "strmm.f"
    if (info != 0) {
#line 248 "strmm.f"
	xerbla_("STRMM ", &info, (ftnlen)6);
#line 249 "strmm.f"
	return 0;
#line 250 "strmm.f"
    }

/*     Quick return if possible. */

#line 254 "strmm.f"
    if (*m == 0 || *n == 0) {
#line 254 "strmm.f"
	return 0;
#line 254 "strmm.f"
    }

/*     And when  alpha.eq.zero. */

#line 258 "strmm.f"
    if (*alpha == 0.) {
#line 259 "strmm.f"
	i__1 = *n;
#line 259 "strmm.f"
	for (j = 1; j <= i__1; ++j) {
#line 260 "strmm.f"
	    i__2 = *m;
#line 260 "strmm.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 261 "strmm.f"
		b[i__ + j * b_dim1] = 0.;
#line 262 "strmm.f"
/* L10: */
#line 262 "strmm.f"
	    }
#line 263 "strmm.f"
/* L20: */
#line 263 "strmm.f"
	}
#line 264 "strmm.f"
	return 0;
#line 265 "strmm.f"
    }

/*     Start the operations. */

#line 269 "strmm.f"
    if (lside) {
#line 270 "strmm.f"
	if (lsame_(transa, "N", (ftnlen)1, (ftnlen)1)) {

/*           Form  B := alpha*A*B. */

#line 274 "strmm.f"
	    if (upper) {
#line 275 "strmm.f"
		i__1 = *n;
#line 275 "strmm.f"
		for (j = 1; j <= i__1; ++j) {
#line 276 "strmm.f"
		    i__2 = *m;
#line 276 "strmm.f"
		    for (k = 1; k <= i__2; ++k) {
#line 277 "strmm.f"
			if (b[k + j * b_dim1] != 0.) {
#line 278 "strmm.f"
			    temp = *alpha * b[k + j * b_dim1];
#line 279 "strmm.f"
			    i__3 = k - 1;
#line 279 "strmm.f"
			    for (i__ = 1; i__ <= i__3; ++i__) {
#line 280 "strmm.f"
				b[i__ + j * b_dim1] += temp * a[i__ + k * 
					a_dim1];
#line 281 "strmm.f"
/* L30: */
#line 281 "strmm.f"
			    }
#line 282 "strmm.f"
			    if (nounit) {
#line 282 "strmm.f"
				temp *= a[k + k * a_dim1];
#line 282 "strmm.f"
			    }
#line 283 "strmm.f"
			    b[k + j * b_dim1] = temp;
#line 284 "strmm.f"
			}
#line 285 "strmm.f"
/* L40: */
#line 285 "strmm.f"
		    }
#line 286 "strmm.f"
/* L50: */
#line 286 "strmm.f"
		}
#line 287 "strmm.f"
	    } else {
#line 288 "strmm.f"
		i__1 = *n;
#line 288 "strmm.f"
		for (j = 1; j <= i__1; ++j) {
#line 289 "strmm.f"
		    for (k = *m; k >= 1; --k) {
#line 290 "strmm.f"
			if (b[k + j * b_dim1] != 0.) {
#line 291 "strmm.f"
			    temp = *alpha * b[k + j * b_dim1];
#line 292 "strmm.f"
			    b[k + j * b_dim1] = temp;
#line 293 "strmm.f"
			    if (nounit) {
#line 293 "strmm.f"
				b[k + j * b_dim1] *= a[k + k * a_dim1];
#line 293 "strmm.f"
			    }
#line 294 "strmm.f"
			    i__2 = *m;
#line 294 "strmm.f"
			    for (i__ = k + 1; i__ <= i__2; ++i__) {
#line 295 "strmm.f"
				b[i__ + j * b_dim1] += temp * a[i__ + k * 
					a_dim1];
#line 296 "strmm.f"
/* L60: */
#line 296 "strmm.f"
			    }
#line 297 "strmm.f"
			}
#line 298 "strmm.f"
/* L70: */
#line 298 "strmm.f"
		    }
#line 299 "strmm.f"
/* L80: */
#line 299 "strmm.f"
		}
#line 300 "strmm.f"
	    }
#line 301 "strmm.f"
	} else {

/*           Form  B := alpha*A**T*B. */

#line 305 "strmm.f"
	    if (upper) {
#line 306 "strmm.f"
		i__1 = *n;
#line 306 "strmm.f"
		for (j = 1; j <= i__1; ++j) {
#line 307 "strmm.f"
		    for (i__ = *m; i__ >= 1; --i__) {
#line 308 "strmm.f"
			temp = b[i__ + j * b_dim1];
#line 309 "strmm.f"
			if (nounit) {
#line 309 "strmm.f"
			    temp *= a[i__ + i__ * a_dim1];
#line 309 "strmm.f"
			}
#line 310 "strmm.f"
			i__2 = i__ - 1;
#line 310 "strmm.f"
			for (k = 1; k <= i__2; ++k) {
#line 311 "strmm.f"
			    temp += a[k + i__ * a_dim1] * b[k + j * b_dim1];
#line 312 "strmm.f"
/* L90: */
#line 312 "strmm.f"
			}
#line 313 "strmm.f"
			b[i__ + j * b_dim1] = *alpha * temp;
#line 314 "strmm.f"
/* L100: */
#line 314 "strmm.f"
		    }
#line 315 "strmm.f"
/* L110: */
#line 315 "strmm.f"
		}
#line 316 "strmm.f"
	    } else {
#line 317 "strmm.f"
		i__1 = *n;
#line 317 "strmm.f"
		for (j = 1; j <= i__1; ++j) {
#line 318 "strmm.f"
		    i__2 = *m;
#line 318 "strmm.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 319 "strmm.f"
			temp = b[i__ + j * b_dim1];
#line 320 "strmm.f"
			if (nounit) {
#line 320 "strmm.f"
			    temp *= a[i__ + i__ * a_dim1];
#line 320 "strmm.f"
			}
#line 321 "strmm.f"
			i__3 = *m;
#line 321 "strmm.f"
			for (k = i__ + 1; k <= i__3; ++k) {
#line 322 "strmm.f"
			    temp += a[k + i__ * a_dim1] * b[k + j * b_dim1];
#line 323 "strmm.f"
/* L120: */
#line 323 "strmm.f"
			}
#line 324 "strmm.f"
			b[i__ + j * b_dim1] = *alpha * temp;
#line 325 "strmm.f"
/* L130: */
#line 325 "strmm.f"
		    }
#line 326 "strmm.f"
/* L140: */
#line 326 "strmm.f"
		}
#line 327 "strmm.f"
	    }
#line 328 "strmm.f"
	}
#line 329 "strmm.f"
    } else {
#line 330 "strmm.f"
	if (lsame_(transa, "N", (ftnlen)1, (ftnlen)1)) {

/*           Form  B := alpha*B*A. */

#line 334 "strmm.f"
	    if (upper) {
#line 335 "strmm.f"
		for (j = *n; j >= 1; --j) {
#line 336 "strmm.f"
		    temp = *alpha;
#line 337 "strmm.f"
		    if (nounit) {
#line 337 "strmm.f"
			temp *= a[j + j * a_dim1];
#line 337 "strmm.f"
		    }
#line 338 "strmm.f"
		    i__1 = *m;
#line 338 "strmm.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 339 "strmm.f"
			b[i__ + j * b_dim1] = temp * b[i__ + j * b_dim1];
#line 340 "strmm.f"
/* L150: */
#line 340 "strmm.f"
		    }
#line 341 "strmm.f"
		    i__1 = j - 1;
#line 341 "strmm.f"
		    for (k = 1; k <= i__1; ++k) {
#line 342 "strmm.f"
			if (a[k + j * a_dim1] != 0.) {
#line 343 "strmm.f"
			    temp = *alpha * a[k + j * a_dim1];
#line 344 "strmm.f"
			    i__2 = *m;
#line 344 "strmm.f"
			    for (i__ = 1; i__ <= i__2; ++i__) {
#line 345 "strmm.f"
				b[i__ + j * b_dim1] += temp * b[i__ + k * 
					b_dim1];
#line 346 "strmm.f"
/* L160: */
#line 346 "strmm.f"
			    }
#line 347 "strmm.f"
			}
#line 348 "strmm.f"
/* L170: */
#line 348 "strmm.f"
		    }
#line 349 "strmm.f"
/* L180: */
#line 349 "strmm.f"
		}
#line 350 "strmm.f"
	    } else {
#line 351 "strmm.f"
		i__1 = *n;
#line 351 "strmm.f"
		for (j = 1; j <= i__1; ++j) {
#line 352 "strmm.f"
		    temp = *alpha;
#line 353 "strmm.f"
		    if (nounit) {
#line 353 "strmm.f"
			temp *= a[j + j * a_dim1];
#line 353 "strmm.f"
		    }
#line 354 "strmm.f"
		    i__2 = *m;
#line 354 "strmm.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 355 "strmm.f"
			b[i__ + j * b_dim1] = temp * b[i__ + j * b_dim1];
#line 356 "strmm.f"
/* L190: */
#line 356 "strmm.f"
		    }
#line 357 "strmm.f"
		    i__2 = *n;
#line 357 "strmm.f"
		    for (k = j + 1; k <= i__2; ++k) {
#line 358 "strmm.f"
			if (a[k + j * a_dim1] != 0.) {
#line 359 "strmm.f"
			    temp = *alpha * a[k + j * a_dim1];
#line 360 "strmm.f"
			    i__3 = *m;
#line 360 "strmm.f"
			    for (i__ = 1; i__ <= i__3; ++i__) {
#line 361 "strmm.f"
				b[i__ + j * b_dim1] += temp * b[i__ + k * 
					b_dim1];
#line 362 "strmm.f"
/* L200: */
#line 362 "strmm.f"
			    }
#line 363 "strmm.f"
			}
#line 364 "strmm.f"
/* L210: */
#line 364 "strmm.f"
		    }
#line 365 "strmm.f"
/* L220: */
#line 365 "strmm.f"
		}
#line 366 "strmm.f"
	    }
#line 367 "strmm.f"
	} else {

/*           Form  B := alpha*B*A**T. */

#line 371 "strmm.f"
	    if (upper) {
#line 372 "strmm.f"
		i__1 = *n;
#line 372 "strmm.f"
		for (k = 1; k <= i__1; ++k) {
#line 373 "strmm.f"
		    i__2 = k - 1;
#line 373 "strmm.f"
		    for (j = 1; j <= i__2; ++j) {
#line 374 "strmm.f"
			if (a[j + k * a_dim1] != 0.) {
#line 375 "strmm.f"
			    temp = *alpha * a[j + k * a_dim1];
#line 376 "strmm.f"
			    i__3 = *m;
#line 376 "strmm.f"
			    for (i__ = 1; i__ <= i__3; ++i__) {
#line 377 "strmm.f"
				b[i__ + j * b_dim1] += temp * b[i__ + k * 
					b_dim1];
#line 378 "strmm.f"
/* L230: */
#line 378 "strmm.f"
			    }
#line 379 "strmm.f"
			}
#line 380 "strmm.f"
/* L240: */
#line 380 "strmm.f"
		    }
#line 381 "strmm.f"
		    temp = *alpha;
#line 382 "strmm.f"
		    if (nounit) {
#line 382 "strmm.f"
			temp *= a[k + k * a_dim1];
#line 382 "strmm.f"
		    }
#line 383 "strmm.f"
		    if (temp != 1.) {
#line 384 "strmm.f"
			i__2 = *m;
#line 384 "strmm.f"
			for (i__ = 1; i__ <= i__2; ++i__) {
#line 385 "strmm.f"
			    b[i__ + k * b_dim1] = temp * b[i__ + k * b_dim1];
#line 386 "strmm.f"
/* L250: */
#line 386 "strmm.f"
			}
#line 387 "strmm.f"
		    }
#line 388 "strmm.f"
/* L260: */
#line 388 "strmm.f"
		}
#line 389 "strmm.f"
	    } else {
#line 390 "strmm.f"
		for (k = *n; k >= 1; --k) {
#line 391 "strmm.f"
		    i__1 = *n;
#line 391 "strmm.f"
		    for (j = k + 1; j <= i__1; ++j) {
#line 392 "strmm.f"
			if (a[j + k * a_dim1] != 0.) {
#line 393 "strmm.f"
			    temp = *alpha * a[j + k * a_dim1];
#line 394 "strmm.f"
			    i__2 = *m;
#line 394 "strmm.f"
			    for (i__ = 1; i__ <= i__2; ++i__) {
#line 395 "strmm.f"
				b[i__ + j * b_dim1] += temp * b[i__ + k * 
					b_dim1];
#line 396 "strmm.f"
/* L270: */
#line 396 "strmm.f"
			    }
#line 397 "strmm.f"
			}
#line 398 "strmm.f"
/* L280: */
#line 398 "strmm.f"
		    }
#line 399 "strmm.f"
		    temp = *alpha;
#line 400 "strmm.f"
		    if (nounit) {
#line 400 "strmm.f"
			temp *= a[k + k * a_dim1];
#line 400 "strmm.f"
		    }
#line 401 "strmm.f"
		    if (temp != 1.) {
#line 402 "strmm.f"
			i__1 = *m;
#line 402 "strmm.f"
			for (i__ = 1; i__ <= i__1; ++i__) {
#line 403 "strmm.f"
			    b[i__ + k * b_dim1] = temp * b[i__ + k * b_dim1];
#line 404 "strmm.f"
/* L290: */
#line 404 "strmm.f"
			}
#line 405 "strmm.f"
		    }
#line 406 "strmm.f"
/* L300: */
#line 406 "strmm.f"
		}
#line 407 "strmm.f"
	    }
#line 408 "strmm.f"
	}
#line 409 "strmm.f"
    }

#line 411 "strmm.f"
    return 0;

/*     End of STRMM . */

} /* strmm_ */

