#line 1 "ssyr2k.f"
/* ssyr2k.f -- translated by f2c (version 20100827).
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

#line 1 "ssyr2k.f"
/* > \brief \b SSYR2K */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SSYR2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC) */

/*       .. Scalar Arguments .. */
/*       REAL ALPHA,BETA */
/*       INTEGER K,LDA,LDB,LDC,N */
/*       CHARACTER TRANS,UPLO */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL A(LDA,*),B(LDB,*),C(LDC,*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SSYR2K  performs one of the symmetric rank 2k operations */
/* > */
/* >    C := alpha*A*B**T + alpha*B*A**T + beta*C, */
/* > */
/* > or */
/* > */
/* >    C := alpha*A**T*B + alpha*B**T*A + beta*C, */
/* > */
/* > where  alpha and beta  are scalars, C is an  n by n  symmetric matrix */
/* > and  A and B  are  n by k  matrices  in the  first  case  and  k by n */
/* > matrices in the second case. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >           On  entry,   UPLO  specifies  whether  the  upper  or  lower */
/* >           triangular  part  of the  array  C  is to be  referenced  as */
/* >           follows: */
/* > */
/* >              UPLO = 'U' or 'u'   Only the  upper triangular part of  C */
/* >                                  is to be referenced. */
/* > */
/* >              UPLO = 'L' or 'l'   Only the  lower triangular part of  C */
/* >                                  is to be referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] TRANS */
/* > \verbatim */
/* >          TRANS is CHARACTER*1 */
/* >           On entry,  TRANS  specifies the operation to be performed as */
/* >           follows: */
/* > */
/* >              TRANS = 'N' or 'n'   C := alpha*A*B**T + alpha*B*A**T + */
/* >                                        beta*C. */
/* > */
/* >              TRANS = 'T' or 't'   C := alpha*A**T*B + alpha*B**T*A + */
/* >                                        beta*C. */
/* > */
/* >              TRANS = 'C' or 'c'   C := alpha*A**T*B + alpha*B**T*A + */
/* >                                        beta*C. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >           On entry,  N specifies the order of the matrix C.  N must be */
/* >           at least zero. */
/* > \endverbatim */
/* > */
/* > \param[in] K */
/* > \verbatim */
/* >          K is INTEGER */
/* >           On entry with  TRANS = 'N' or 'n',  K  specifies  the number */
/* >           of  columns  of the  matrices  A and B,  and on  entry  with */
/* >           TRANS = 'T' or 't' or 'C' or 'c',  K  specifies  the  number */
/* >           of rows of the matrices  A and B.  K must be at least  zero. */
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
/* >          A is REAL array of DIMENSION ( LDA, ka ), where ka is */
/* >           k  when  TRANS = 'N' or 'n',  and is  n  otherwise. */
/* >           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k */
/* >           part of the array  A  must contain the matrix  A,  otherwise */
/* >           the leading  k by n  part of the array  A  must contain  the */
/* >           matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >           On entry, LDA specifies the first dimension of A as declared */
/* >           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n' */
/* >           then  LDA must be at least  max( 1, n ), otherwise  LDA must */
/* >           be at least  max( 1, k ). */
/* > \endverbatim */
/* > */
/* > \param[in] B */
/* > \verbatim */
/* >          B is REAL array of DIMENSION ( LDB, kb ), where kb is */
/* >           k  when  TRANS = 'N' or 'n',  and is  n  otherwise. */
/* >           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k */
/* >           part of the array  B  must contain the matrix  B,  otherwise */
/* >           the leading  k by n  part of the array  B  must contain  the */
/* >           matrix B. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >           On entry, LDB specifies the first dimension of B as declared */
/* >           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n' */
/* >           then  LDB must be at least  max( 1, n ), otherwise  LDB must */
/* >           be at least  max( 1, k ). */
/* > \endverbatim */
/* > */
/* > \param[in] BETA */
/* > \verbatim */
/* >          BETA is REAL */
/* >           On entry, BETA specifies the scalar beta. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* >          C is REAL array of DIMENSION ( LDC, n ). */
/* >           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n */
/* >           upper triangular part of the array C must contain the upper */
/* >           triangular part  of the  symmetric matrix  and the strictly */
/* >           lower triangular part of C is not referenced.  On exit, the */
/* >           upper triangular part of the array  C is overwritten by the */
/* >           upper triangular part of the updated matrix. */
/* >           Before entry  with  UPLO = 'L' or 'l',  the leading  n by n */
/* >           lower triangular part of the array C must contain the lower */
/* >           triangular part  of the  symmetric matrix  and the strictly */
/* >           upper triangular part of C is not referenced.  On exit, the */
/* >           lower triangular part of the array  C is overwritten by the */
/* >           lower triangular part of the updated matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] LDC */
/* > \verbatim */
/* >          LDC is INTEGER */
/* >           On entry, LDC specifies the first dimension of C as declared */
/* >           in  the  calling  (sub)  program.   LDC  must  be  at  least */
/* >           max( 1, n ). */
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
/* Subroutine */ int ssyr2k_(char *uplo, char *trans, integer *n, integer *k, 
	doublereal *alpha, doublereal *a, integer *lda, doublereal *b, 
	integer *ldb, doublereal *beta, doublereal *c__, integer *ldc, ftnlen 
	uplo_len, ftnlen trans_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2, 
	    i__3;

    /* Local variables */
    static integer i__, j, l, info;
    static doublereal temp1, temp2;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer nrowa;
    static logical upper;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);


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

#line 233 "ssyr2k.f"
    /* Parameter adjustments */
#line 233 "ssyr2k.f"
    a_dim1 = *lda;
#line 233 "ssyr2k.f"
    a_offset = 1 + a_dim1;
#line 233 "ssyr2k.f"
    a -= a_offset;
#line 233 "ssyr2k.f"
    b_dim1 = *ldb;
#line 233 "ssyr2k.f"
    b_offset = 1 + b_dim1;
#line 233 "ssyr2k.f"
    b -= b_offset;
#line 233 "ssyr2k.f"
    c_dim1 = *ldc;
#line 233 "ssyr2k.f"
    c_offset = 1 + c_dim1;
#line 233 "ssyr2k.f"
    c__ -= c_offset;
#line 233 "ssyr2k.f"

#line 233 "ssyr2k.f"
    /* Function Body */
#line 233 "ssyr2k.f"
    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {
#line 234 "ssyr2k.f"
	nrowa = *n;
#line 235 "ssyr2k.f"
    } else {
#line 236 "ssyr2k.f"
	nrowa = *k;
#line 237 "ssyr2k.f"
    }
#line 238 "ssyr2k.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);

#line 240 "ssyr2k.f"
    info = 0;
#line 241 "ssyr2k.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 242 "ssyr2k.f"
	info = 1;
#line 243 "ssyr2k.f"
    } else if (! lsame_(trans, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, 
	    "T", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, "C", (ftnlen)1, (
	    ftnlen)1)) {
#line 246 "ssyr2k.f"
	info = 2;
#line 247 "ssyr2k.f"
    } else if (*n < 0) {
#line 248 "ssyr2k.f"
	info = 3;
#line 249 "ssyr2k.f"
    } else if (*k < 0) {
#line 250 "ssyr2k.f"
	info = 4;
#line 251 "ssyr2k.f"
    } else if (*lda < max(1,nrowa)) {
#line 252 "ssyr2k.f"
	info = 7;
#line 253 "ssyr2k.f"
    } else if (*ldb < max(1,nrowa)) {
#line 254 "ssyr2k.f"
	info = 9;
#line 255 "ssyr2k.f"
    } else if (*ldc < max(1,*n)) {
#line 256 "ssyr2k.f"
	info = 12;
#line 257 "ssyr2k.f"
    }
#line 258 "ssyr2k.f"
    if (info != 0) {
#line 259 "ssyr2k.f"
	xerbla_("SSYR2K", &info, (ftnlen)6);
#line 260 "ssyr2k.f"
	return 0;
#line 261 "ssyr2k.f"
    }

/*     Quick return if possible. */

#line 265 "ssyr2k.f"
    if (*n == 0 || (*alpha == 0. || *k == 0) && *beta == 1.) {
#line 265 "ssyr2k.f"
	return 0;
#line 265 "ssyr2k.f"
    }

/*     And when  alpha.eq.zero. */

#line 270 "ssyr2k.f"
    if (*alpha == 0.) {
#line 271 "ssyr2k.f"
	if (upper) {
#line 272 "ssyr2k.f"
	    if (*beta == 0.) {
#line 273 "ssyr2k.f"
		i__1 = *n;
#line 273 "ssyr2k.f"
		for (j = 1; j <= i__1; ++j) {
#line 274 "ssyr2k.f"
		    i__2 = j;
#line 274 "ssyr2k.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 275 "ssyr2k.f"
			c__[i__ + j * c_dim1] = 0.;
#line 276 "ssyr2k.f"
/* L10: */
#line 276 "ssyr2k.f"
		    }
#line 277 "ssyr2k.f"
/* L20: */
#line 277 "ssyr2k.f"
		}
#line 278 "ssyr2k.f"
	    } else {
#line 279 "ssyr2k.f"
		i__1 = *n;
#line 279 "ssyr2k.f"
		for (j = 1; j <= i__1; ++j) {
#line 280 "ssyr2k.f"
		    i__2 = j;
#line 280 "ssyr2k.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 281 "ssyr2k.f"
			c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
#line 282 "ssyr2k.f"
/* L30: */
#line 282 "ssyr2k.f"
		    }
#line 283 "ssyr2k.f"
/* L40: */
#line 283 "ssyr2k.f"
		}
#line 284 "ssyr2k.f"
	    }
#line 285 "ssyr2k.f"
	} else {
#line 286 "ssyr2k.f"
	    if (*beta == 0.) {
#line 287 "ssyr2k.f"
		i__1 = *n;
#line 287 "ssyr2k.f"
		for (j = 1; j <= i__1; ++j) {
#line 288 "ssyr2k.f"
		    i__2 = *n;
#line 288 "ssyr2k.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 289 "ssyr2k.f"
			c__[i__ + j * c_dim1] = 0.;
#line 290 "ssyr2k.f"
/* L50: */
#line 290 "ssyr2k.f"
		    }
#line 291 "ssyr2k.f"
/* L60: */
#line 291 "ssyr2k.f"
		}
#line 292 "ssyr2k.f"
	    } else {
#line 293 "ssyr2k.f"
		i__1 = *n;
#line 293 "ssyr2k.f"
		for (j = 1; j <= i__1; ++j) {
#line 294 "ssyr2k.f"
		    i__2 = *n;
#line 294 "ssyr2k.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 295 "ssyr2k.f"
			c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
#line 296 "ssyr2k.f"
/* L70: */
#line 296 "ssyr2k.f"
		    }
#line 297 "ssyr2k.f"
/* L80: */
#line 297 "ssyr2k.f"
		}
#line 298 "ssyr2k.f"
	    }
#line 299 "ssyr2k.f"
	}
#line 300 "ssyr2k.f"
	return 0;
#line 301 "ssyr2k.f"
    }

/*     Start the operations. */

#line 305 "ssyr2k.f"
    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {

/*        Form  C := alpha*A*B**T + alpha*B*A**T + C. */

#line 309 "ssyr2k.f"
	if (upper) {
#line 310 "ssyr2k.f"
	    i__1 = *n;
#line 310 "ssyr2k.f"
	    for (j = 1; j <= i__1; ++j) {
#line 311 "ssyr2k.f"
		if (*beta == 0.) {
#line 312 "ssyr2k.f"
		    i__2 = j;
#line 312 "ssyr2k.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 313 "ssyr2k.f"
			c__[i__ + j * c_dim1] = 0.;
#line 314 "ssyr2k.f"
/* L90: */
#line 314 "ssyr2k.f"
		    }
#line 315 "ssyr2k.f"
		} else if (*beta != 1.) {
#line 316 "ssyr2k.f"
		    i__2 = j;
#line 316 "ssyr2k.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 317 "ssyr2k.f"
			c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
#line 318 "ssyr2k.f"
/* L100: */
#line 318 "ssyr2k.f"
		    }
#line 319 "ssyr2k.f"
		}
#line 320 "ssyr2k.f"
		i__2 = *k;
#line 320 "ssyr2k.f"
		for (l = 1; l <= i__2; ++l) {
#line 321 "ssyr2k.f"
		    if (a[j + l * a_dim1] != 0. || b[j + l * b_dim1] != 0.) {
#line 322 "ssyr2k.f"
			temp1 = *alpha * b[j + l * b_dim1];
#line 323 "ssyr2k.f"
			temp2 = *alpha * a[j + l * a_dim1];
#line 324 "ssyr2k.f"
			i__3 = j;
#line 324 "ssyr2k.f"
			for (i__ = 1; i__ <= i__3; ++i__) {
#line 325 "ssyr2k.f"
			    c__[i__ + j * c_dim1] = c__[i__ + j * c_dim1] + a[
				    i__ + l * a_dim1] * temp1 + b[i__ + l * 
				    b_dim1] * temp2;
#line 327 "ssyr2k.f"
/* L110: */
#line 327 "ssyr2k.f"
			}
#line 328 "ssyr2k.f"
		    }
#line 329 "ssyr2k.f"
/* L120: */
#line 329 "ssyr2k.f"
		}
#line 330 "ssyr2k.f"
/* L130: */
#line 330 "ssyr2k.f"
	    }
#line 331 "ssyr2k.f"
	} else {
#line 332 "ssyr2k.f"
	    i__1 = *n;
#line 332 "ssyr2k.f"
	    for (j = 1; j <= i__1; ++j) {
#line 333 "ssyr2k.f"
		if (*beta == 0.) {
#line 334 "ssyr2k.f"
		    i__2 = *n;
#line 334 "ssyr2k.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 335 "ssyr2k.f"
			c__[i__ + j * c_dim1] = 0.;
#line 336 "ssyr2k.f"
/* L140: */
#line 336 "ssyr2k.f"
		    }
#line 337 "ssyr2k.f"
		} else if (*beta != 1.) {
#line 338 "ssyr2k.f"
		    i__2 = *n;
#line 338 "ssyr2k.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 339 "ssyr2k.f"
			c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
#line 340 "ssyr2k.f"
/* L150: */
#line 340 "ssyr2k.f"
		    }
#line 341 "ssyr2k.f"
		}
#line 342 "ssyr2k.f"
		i__2 = *k;
#line 342 "ssyr2k.f"
		for (l = 1; l <= i__2; ++l) {
#line 343 "ssyr2k.f"
		    if (a[j + l * a_dim1] != 0. || b[j + l * b_dim1] != 0.) {
#line 344 "ssyr2k.f"
			temp1 = *alpha * b[j + l * b_dim1];
#line 345 "ssyr2k.f"
			temp2 = *alpha * a[j + l * a_dim1];
#line 346 "ssyr2k.f"
			i__3 = *n;
#line 346 "ssyr2k.f"
			for (i__ = j; i__ <= i__3; ++i__) {
#line 347 "ssyr2k.f"
			    c__[i__ + j * c_dim1] = c__[i__ + j * c_dim1] + a[
				    i__ + l * a_dim1] * temp1 + b[i__ + l * 
				    b_dim1] * temp2;
#line 349 "ssyr2k.f"
/* L160: */
#line 349 "ssyr2k.f"
			}
#line 350 "ssyr2k.f"
		    }
#line 351 "ssyr2k.f"
/* L170: */
#line 351 "ssyr2k.f"
		}
#line 352 "ssyr2k.f"
/* L180: */
#line 352 "ssyr2k.f"
	    }
#line 353 "ssyr2k.f"
	}
#line 354 "ssyr2k.f"
    } else {

/*        Form  C := alpha*A**T*B + alpha*B**T*A + C. */

#line 358 "ssyr2k.f"
	if (upper) {
#line 359 "ssyr2k.f"
	    i__1 = *n;
#line 359 "ssyr2k.f"
	    for (j = 1; j <= i__1; ++j) {
#line 360 "ssyr2k.f"
		i__2 = j;
#line 360 "ssyr2k.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 361 "ssyr2k.f"
		    temp1 = 0.;
#line 362 "ssyr2k.f"
		    temp2 = 0.;
#line 363 "ssyr2k.f"
		    i__3 = *k;
#line 363 "ssyr2k.f"
		    for (l = 1; l <= i__3; ++l) {
#line 364 "ssyr2k.f"
			temp1 += a[l + i__ * a_dim1] * b[l + j * b_dim1];
#line 365 "ssyr2k.f"
			temp2 += b[l + i__ * b_dim1] * a[l + j * a_dim1];
#line 366 "ssyr2k.f"
/* L190: */
#line 366 "ssyr2k.f"
		    }
#line 367 "ssyr2k.f"
		    if (*beta == 0.) {
#line 368 "ssyr2k.f"
			c__[i__ + j * c_dim1] = *alpha * temp1 + *alpha * 
				temp2;
#line 369 "ssyr2k.f"
		    } else {
#line 370 "ssyr2k.f"
			c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1] 
				+ *alpha * temp1 + *alpha * temp2;
#line 372 "ssyr2k.f"
		    }
#line 373 "ssyr2k.f"
/* L200: */
#line 373 "ssyr2k.f"
		}
#line 374 "ssyr2k.f"
/* L210: */
#line 374 "ssyr2k.f"
	    }
#line 375 "ssyr2k.f"
	} else {
#line 376 "ssyr2k.f"
	    i__1 = *n;
#line 376 "ssyr2k.f"
	    for (j = 1; j <= i__1; ++j) {
#line 377 "ssyr2k.f"
		i__2 = *n;
#line 377 "ssyr2k.f"
		for (i__ = j; i__ <= i__2; ++i__) {
#line 378 "ssyr2k.f"
		    temp1 = 0.;
#line 379 "ssyr2k.f"
		    temp2 = 0.;
#line 380 "ssyr2k.f"
		    i__3 = *k;
#line 380 "ssyr2k.f"
		    for (l = 1; l <= i__3; ++l) {
#line 381 "ssyr2k.f"
			temp1 += a[l + i__ * a_dim1] * b[l + j * b_dim1];
#line 382 "ssyr2k.f"
			temp2 += b[l + i__ * b_dim1] * a[l + j * a_dim1];
#line 383 "ssyr2k.f"
/* L220: */
#line 383 "ssyr2k.f"
		    }
#line 384 "ssyr2k.f"
		    if (*beta == 0.) {
#line 385 "ssyr2k.f"
			c__[i__ + j * c_dim1] = *alpha * temp1 + *alpha * 
				temp2;
#line 386 "ssyr2k.f"
		    } else {
#line 387 "ssyr2k.f"
			c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1] 
				+ *alpha * temp1 + *alpha * temp2;
#line 389 "ssyr2k.f"
		    }
#line 390 "ssyr2k.f"
/* L230: */
#line 390 "ssyr2k.f"
		}
#line 391 "ssyr2k.f"
/* L240: */
#line 391 "ssyr2k.f"
	    }
#line 392 "ssyr2k.f"
	}
#line 393 "ssyr2k.f"
    }

#line 395 "ssyr2k.f"
    return 0;

/*     End of SSYR2K. */

} /* ssyr2k_ */

