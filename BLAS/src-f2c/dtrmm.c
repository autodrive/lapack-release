#line 1 "dtrmm.f"
/* dtrmm.f -- translated by f2c (version 20100827).
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

#line 1 "dtrmm.f"
/* > \brief \b DTRMM */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DTRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB) */

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
/* > DTRMM  performs one of the matrix-matrix operations */
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
/* >          ALPHA is DOUBLE PRECISION. */
/* >           On entry,  ALPHA specifies the scalar  alpha. When  alpha is */
/* >           zero then  A is not referenced and  B need not be set before */
/* >           entry. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >           A is DOUBLE PRECISION array of DIMENSION ( LDA, k ), where k is m */
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
/* >          B is DOUBLE PRECISION array of DIMENSION ( LDB, n ). */
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

/* > \ingroup double_blas_level3 */

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
/* Subroutine */ int dtrmm_(char *side, char *uplo, char *transa, char *diag, 
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

#line 218 "dtrmm.f"
    /* Parameter adjustments */
#line 218 "dtrmm.f"
    a_dim1 = *lda;
#line 218 "dtrmm.f"
    a_offset = 1 + a_dim1;
#line 218 "dtrmm.f"
    a -= a_offset;
#line 218 "dtrmm.f"
    b_dim1 = *ldb;
#line 218 "dtrmm.f"
    b_offset = 1 + b_dim1;
#line 218 "dtrmm.f"
    b -= b_offset;
#line 218 "dtrmm.f"

#line 218 "dtrmm.f"
    /* Function Body */
#line 218 "dtrmm.f"
    lside = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 219 "dtrmm.f"
    if (lside) {
#line 220 "dtrmm.f"
	nrowa = *m;
#line 221 "dtrmm.f"
    } else {
#line 222 "dtrmm.f"
	nrowa = *n;
#line 223 "dtrmm.f"
    }
#line 224 "dtrmm.f"
    nounit = lsame_(diag, "N", (ftnlen)1, (ftnlen)1);
#line 225 "dtrmm.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);

#line 227 "dtrmm.f"
    info = 0;
#line 228 "dtrmm.f"
    if (! lside && ! lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {
#line 229 "dtrmm.f"
	info = 1;
#line 230 "dtrmm.f"
    } else if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 231 "dtrmm.f"
	info = 2;
#line 232 "dtrmm.f"
    } else if (! lsame_(transa, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(transa,
	     "T", (ftnlen)1, (ftnlen)1) && ! lsame_(transa, "C", (ftnlen)1, (
	    ftnlen)1)) {
#line 235 "dtrmm.f"
	info = 3;
#line 236 "dtrmm.f"
    } else if (! lsame_(diag, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(diag, 
	    "N", (ftnlen)1, (ftnlen)1)) {
#line 237 "dtrmm.f"
	info = 4;
#line 238 "dtrmm.f"
    } else if (*m < 0) {
#line 239 "dtrmm.f"
	info = 5;
#line 240 "dtrmm.f"
    } else if (*n < 0) {
#line 241 "dtrmm.f"
	info = 6;
#line 242 "dtrmm.f"
    } else if (*lda < max(1,nrowa)) {
#line 243 "dtrmm.f"
	info = 9;
#line 244 "dtrmm.f"
    } else if (*ldb < max(1,*m)) {
#line 245 "dtrmm.f"
	info = 11;
#line 246 "dtrmm.f"
    }
#line 247 "dtrmm.f"
    if (info != 0) {
#line 248 "dtrmm.f"
	xerbla_("DTRMM ", &info, (ftnlen)6);
#line 249 "dtrmm.f"
	return 0;
#line 250 "dtrmm.f"
    }

/*     Quick return if possible. */

#line 254 "dtrmm.f"
    if (*m == 0 || *n == 0) {
#line 254 "dtrmm.f"
	return 0;
#line 254 "dtrmm.f"
    }

/*     And when  alpha.eq.zero. */

#line 258 "dtrmm.f"
    if (*alpha == 0.) {
#line 259 "dtrmm.f"
	i__1 = *n;
#line 259 "dtrmm.f"
	for (j = 1; j <= i__1; ++j) {
#line 260 "dtrmm.f"
	    i__2 = *m;
#line 260 "dtrmm.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 261 "dtrmm.f"
		b[i__ + j * b_dim1] = 0.;
#line 262 "dtrmm.f"
/* L10: */
#line 262 "dtrmm.f"
	    }
#line 263 "dtrmm.f"
/* L20: */
#line 263 "dtrmm.f"
	}
#line 264 "dtrmm.f"
	return 0;
#line 265 "dtrmm.f"
    }

/*     Start the operations. */

#line 269 "dtrmm.f"
    if (lside) {
#line 270 "dtrmm.f"
	if (lsame_(transa, "N", (ftnlen)1, (ftnlen)1)) {

/*           Form  B := alpha*A*B. */

#line 274 "dtrmm.f"
	    if (upper) {
#line 275 "dtrmm.f"
		i__1 = *n;
#line 275 "dtrmm.f"
		for (j = 1; j <= i__1; ++j) {
#line 276 "dtrmm.f"
		    i__2 = *m;
#line 276 "dtrmm.f"
		    for (k = 1; k <= i__2; ++k) {
#line 277 "dtrmm.f"
			if (b[k + j * b_dim1] != 0.) {
#line 278 "dtrmm.f"
			    temp = *alpha * b[k + j * b_dim1];
#line 279 "dtrmm.f"
			    i__3 = k - 1;
#line 279 "dtrmm.f"
			    for (i__ = 1; i__ <= i__3; ++i__) {
#line 280 "dtrmm.f"
				b[i__ + j * b_dim1] += temp * a[i__ + k * 
					a_dim1];
#line 281 "dtrmm.f"
/* L30: */
#line 281 "dtrmm.f"
			    }
#line 282 "dtrmm.f"
			    if (nounit) {
#line 282 "dtrmm.f"
				temp *= a[k + k * a_dim1];
#line 282 "dtrmm.f"
			    }
#line 283 "dtrmm.f"
			    b[k + j * b_dim1] = temp;
#line 284 "dtrmm.f"
			}
#line 285 "dtrmm.f"
/* L40: */
#line 285 "dtrmm.f"
		    }
#line 286 "dtrmm.f"
/* L50: */
#line 286 "dtrmm.f"
		}
#line 287 "dtrmm.f"
	    } else {
#line 288 "dtrmm.f"
		i__1 = *n;
#line 288 "dtrmm.f"
		for (j = 1; j <= i__1; ++j) {
#line 289 "dtrmm.f"
		    for (k = *m; k >= 1; --k) {
#line 290 "dtrmm.f"
			if (b[k + j * b_dim1] != 0.) {
#line 291 "dtrmm.f"
			    temp = *alpha * b[k + j * b_dim1];
#line 292 "dtrmm.f"
			    b[k + j * b_dim1] = temp;
#line 293 "dtrmm.f"
			    if (nounit) {
#line 293 "dtrmm.f"
				b[k + j * b_dim1] *= a[k + k * a_dim1];
#line 293 "dtrmm.f"
			    }
#line 294 "dtrmm.f"
			    i__2 = *m;
#line 294 "dtrmm.f"
			    for (i__ = k + 1; i__ <= i__2; ++i__) {
#line 295 "dtrmm.f"
				b[i__ + j * b_dim1] += temp * a[i__ + k * 
					a_dim1];
#line 296 "dtrmm.f"
/* L60: */
#line 296 "dtrmm.f"
			    }
#line 297 "dtrmm.f"
			}
#line 298 "dtrmm.f"
/* L70: */
#line 298 "dtrmm.f"
		    }
#line 299 "dtrmm.f"
/* L80: */
#line 299 "dtrmm.f"
		}
#line 300 "dtrmm.f"
	    }
#line 301 "dtrmm.f"
	} else {

/*           Form  B := alpha*A**T*B. */

#line 305 "dtrmm.f"
	    if (upper) {
#line 306 "dtrmm.f"
		i__1 = *n;
#line 306 "dtrmm.f"
		for (j = 1; j <= i__1; ++j) {
#line 307 "dtrmm.f"
		    for (i__ = *m; i__ >= 1; --i__) {
#line 308 "dtrmm.f"
			temp = b[i__ + j * b_dim1];
#line 309 "dtrmm.f"
			if (nounit) {
#line 309 "dtrmm.f"
			    temp *= a[i__ + i__ * a_dim1];
#line 309 "dtrmm.f"
			}
#line 310 "dtrmm.f"
			i__2 = i__ - 1;
#line 310 "dtrmm.f"
			for (k = 1; k <= i__2; ++k) {
#line 311 "dtrmm.f"
			    temp += a[k + i__ * a_dim1] * b[k + j * b_dim1];
#line 312 "dtrmm.f"
/* L90: */
#line 312 "dtrmm.f"
			}
#line 313 "dtrmm.f"
			b[i__ + j * b_dim1] = *alpha * temp;
#line 314 "dtrmm.f"
/* L100: */
#line 314 "dtrmm.f"
		    }
#line 315 "dtrmm.f"
/* L110: */
#line 315 "dtrmm.f"
		}
#line 316 "dtrmm.f"
	    } else {
#line 317 "dtrmm.f"
		i__1 = *n;
#line 317 "dtrmm.f"
		for (j = 1; j <= i__1; ++j) {
#line 318 "dtrmm.f"
		    i__2 = *m;
#line 318 "dtrmm.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 319 "dtrmm.f"
			temp = b[i__ + j * b_dim1];
#line 320 "dtrmm.f"
			if (nounit) {
#line 320 "dtrmm.f"
			    temp *= a[i__ + i__ * a_dim1];
#line 320 "dtrmm.f"
			}
#line 321 "dtrmm.f"
			i__3 = *m;
#line 321 "dtrmm.f"
			for (k = i__ + 1; k <= i__3; ++k) {
#line 322 "dtrmm.f"
			    temp += a[k + i__ * a_dim1] * b[k + j * b_dim1];
#line 323 "dtrmm.f"
/* L120: */
#line 323 "dtrmm.f"
			}
#line 324 "dtrmm.f"
			b[i__ + j * b_dim1] = *alpha * temp;
#line 325 "dtrmm.f"
/* L130: */
#line 325 "dtrmm.f"
		    }
#line 326 "dtrmm.f"
/* L140: */
#line 326 "dtrmm.f"
		}
#line 327 "dtrmm.f"
	    }
#line 328 "dtrmm.f"
	}
#line 329 "dtrmm.f"
    } else {
#line 330 "dtrmm.f"
	if (lsame_(transa, "N", (ftnlen)1, (ftnlen)1)) {

/*           Form  B := alpha*B*A. */

#line 334 "dtrmm.f"
	    if (upper) {
#line 335 "dtrmm.f"
		for (j = *n; j >= 1; --j) {
#line 336 "dtrmm.f"
		    temp = *alpha;
#line 337 "dtrmm.f"
		    if (nounit) {
#line 337 "dtrmm.f"
			temp *= a[j + j * a_dim1];
#line 337 "dtrmm.f"
		    }
#line 338 "dtrmm.f"
		    i__1 = *m;
#line 338 "dtrmm.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 339 "dtrmm.f"
			b[i__ + j * b_dim1] = temp * b[i__ + j * b_dim1];
#line 340 "dtrmm.f"
/* L150: */
#line 340 "dtrmm.f"
		    }
#line 341 "dtrmm.f"
		    i__1 = j - 1;
#line 341 "dtrmm.f"
		    for (k = 1; k <= i__1; ++k) {
#line 342 "dtrmm.f"
			if (a[k + j * a_dim1] != 0.) {
#line 343 "dtrmm.f"
			    temp = *alpha * a[k + j * a_dim1];
#line 344 "dtrmm.f"
			    i__2 = *m;
#line 344 "dtrmm.f"
			    for (i__ = 1; i__ <= i__2; ++i__) {
#line 345 "dtrmm.f"
				b[i__ + j * b_dim1] += temp * b[i__ + k * 
					b_dim1];
#line 346 "dtrmm.f"
/* L160: */
#line 346 "dtrmm.f"
			    }
#line 347 "dtrmm.f"
			}
#line 348 "dtrmm.f"
/* L170: */
#line 348 "dtrmm.f"
		    }
#line 349 "dtrmm.f"
/* L180: */
#line 349 "dtrmm.f"
		}
#line 350 "dtrmm.f"
	    } else {
#line 351 "dtrmm.f"
		i__1 = *n;
#line 351 "dtrmm.f"
		for (j = 1; j <= i__1; ++j) {
#line 352 "dtrmm.f"
		    temp = *alpha;
#line 353 "dtrmm.f"
		    if (nounit) {
#line 353 "dtrmm.f"
			temp *= a[j + j * a_dim1];
#line 353 "dtrmm.f"
		    }
#line 354 "dtrmm.f"
		    i__2 = *m;
#line 354 "dtrmm.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 355 "dtrmm.f"
			b[i__ + j * b_dim1] = temp * b[i__ + j * b_dim1];
#line 356 "dtrmm.f"
/* L190: */
#line 356 "dtrmm.f"
		    }
#line 357 "dtrmm.f"
		    i__2 = *n;
#line 357 "dtrmm.f"
		    for (k = j + 1; k <= i__2; ++k) {
#line 358 "dtrmm.f"
			if (a[k + j * a_dim1] != 0.) {
#line 359 "dtrmm.f"
			    temp = *alpha * a[k + j * a_dim1];
#line 360 "dtrmm.f"
			    i__3 = *m;
#line 360 "dtrmm.f"
			    for (i__ = 1; i__ <= i__3; ++i__) {
#line 361 "dtrmm.f"
				b[i__ + j * b_dim1] += temp * b[i__ + k * 
					b_dim1];
#line 362 "dtrmm.f"
/* L200: */
#line 362 "dtrmm.f"
			    }
#line 363 "dtrmm.f"
			}
#line 364 "dtrmm.f"
/* L210: */
#line 364 "dtrmm.f"
		    }
#line 365 "dtrmm.f"
/* L220: */
#line 365 "dtrmm.f"
		}
#line 366 "dtrmm.f"
	    }
#line 367 "dtrmm.f"
	} else {

/*           Form  B := alpha*B*A**T. */

#line 371 "dtrmm.f"
	    if (upper) {
#line 372 "dtrmm.f"
		i__1 = *n;
#line 372 "dtrmm.f"
		for (k = 1; k <= i__1; ++k) {
#line 373 "dtrmm.f"
		    i__2 = k - 1;
#line 373 "dtrmm.f"
		    for (j = 1; j <= i__2; ++j) {
#line 374 "dtrmm.f"
			if (a[j + k * a_dim1] != 0.) {
#line 375 "dtrmm.f"
			    temp = *alpha * a[j + k * a_dim1];
#line 376 "dtrmm.f"
			    i__3 = *m;
#line 376 "dtrmm.f"
			    for (i__ = 1; i__ <= i__3; ++i__) {
#line 377 "dtrmm.f"
				b[i__ + j * b_dim1] += temp * b[i__ + k * 
					b_dim1];
#line 378 "dtrmm.f"
/* L230: */
#line 378 "dtrmm.f"
			    }
#line 379 "dtrmm.f"
			}
#line 380 "dtrmm.f"
/* L240: */
#line 380 "dtrmm.f"
		    }
#line 381 "dtrmm.f"
		    temp = *alpha;
#line 382 "dtrmm.f"
		    if (nounit) {
#line 382 "dtrmm.f"
			temp *= a[k + k * a_dim1];
#line 382 "dtrmm.f"
		    }
#line 383 "dtrmm.f"
		    if (temp != 1.) {
#line 384 "dtrmm.f"
			i__2 = *m;
#line 384 "dtrmm.f"
			for (i__ = 1; i__ <= i__2; ++i__) {
#line 385 "dtrmm.f"
			    b[i__ + k * b_dim1] = temp * b[i__ + k * b_dim1];
#line 386 "dtrmm.f"
/* L250: */
#line 386 "dtrmm.f"
			}
#line 387 "dtrmm.f"
		    }
#line 388 "dtrmm.f"
/* L260: */
#line 388 "dtrmm.f"
		}
#line 389 "dtrmm.f"
	    } else {
#line 390 "dtrmm.f"
		for (k = *n; k >= 1; --k) {
#line 391 "dtrmm.f"
		    i__1 = *n;
#line 391 "dtrmm.f"
		    for (j = k + 1; j <= i__1; ++j) {
#line 392 "dtrmm.f"
			if (a[j + k * a_dim1] != 0.) {
#line 393 "dtrmm.f"
			    temp = *alpha * a[j + k * a_dim1];
#line 394 "dtrmm.f"
			    i__2 = *m;
#line 394 "dtrmm.f"
			    for (i__ = 1; i__ <= i__2; ++i__) {
#line 395 "dtrmm.f"
				b[i__ + j * b_dim1] += temp * b[i__ + k * 
					b_dim1];
#line 396 "dtrmm.f"
/* L270: */
#line 396 "dtrmm.f"
			    }
#line 397 "dtrmm.f"
			}
#line 398 "dtrmm.f"
/* L280: */
#line 398 "dtrmm.f"
		    }
#line 399 "dtrmm.f"
		    temp = *alpha;
#line 400 "dtrmm.f"
		    if (nounit) {
#line 400 "dtrmm.f"
			temp *= a[k + k * a_dim1];
#line 400 "dtrmm.f"
		    }
#line 401 "dtrmm.f"
		    if (temp != 1.) {
#line 402 "dtrmm.f"
			i__1 = *m;
#line 402 "dtrmm.f"
			for (i__ = 1; i__ <= i__1; ++i__) {
#line 403 "dtrmm.f"
			    b[i__ + k * b_dim1] = temp * b[i__ + k * b_dim1];
#line 404 "dtrmm.f"
/* L290: */
#line 404 "dtrmm.f"
			}
#line 405 "dtrmm.f"
		    }
#line 406 "dtrmm.f"
/* L300: */
#line 406 "dtrmm.f"
		}
#line 407 "dtrmm.f"
	    }
#line 408 "dtrmm.f"
	}
#line 409 "dtrmm.f"
    }

#line 411 "dtrmm.f"
    return 0;

/*     End of DTRMM . */

} /* dtrmm_ */

