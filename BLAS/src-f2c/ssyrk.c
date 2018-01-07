#line 1 "ssyrk.f"
/* ssyrk.f -- translated by f2c (version 20100827).
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

#line 1 "ssyrk.f"
/* > \brief \b SSYRK */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC) */

/*       .. Scalar Arguments .. */
/*       REAL ALPHA,BETA */
/*       INTEGER K,LDA,LDC,N */
/*       CHARACTER TRANS,UPLO */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL A(LDA,*),C(LDC,*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SSYRK  performs one of the symmetric rank k operations */
/* > */
/* >    C := alpha*A*A**T + beta*C, */
/* > */
/* > or */
/* > */
/* >    C := alpha*A**T*A + beta*C, */
/* > */
/* > where  alpha and beta  are scalars, C is an  n by n  symmetric matrix */
/* > and  A  is an  n by k  matrix in the first case and a  k by n  matrix */
/* > in the second case. */
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
/* >              TRANS = 'N' or 'n'   C := alpha*A*A**T + beta*C. */
/* > */
/* >              TRANS = 'T' or 't'   C := alpha*A**T*A + beta*C. */
/* > */
/* >              TRANS = 'C' or 'c'   C := alpha*A**T*A + beta*C. */
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
/* >           of  columns   of  the   matrix   A,   and  on   entry   with */
/* >           TRANS = 'T' or 't' or 'C' or 'c',  K  specifies  the  number */
/* >           of rows of the matrix  A.  K must be at least zero. */
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
/* >          A is REAL array, dimension ( LDA, ka ), where ka is */
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
/* > \param[in] BETA */
/* > \verbatim */
/* >          BETA is REAL */
/* >           On entry, BETA specifies the scalar beta. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* >          C is REAL array, dimension ( LDC, N ) */
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
/* >  -- Written on 8-February-1989. */
/* >     Jack Dongarra, Argonne National Laboratory. */
/* >     Iain Duff, AERE Harwell. */
/* >     Jeremy Du Croz, Numerical Algorithms Group Ltd. */
/* >     Sven Hammarling, Numerical Algorithms Group Ltd. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int ssyrk_(char *uplo, char *trans, integer *n, integer *k, 
	doublereal *alpha, doublereal *a, integer *lda, doublereal *beta, 
	doublereal *c__, integer *ldc, ftnlen uplo_len, ftnlen trans_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, c_dim1, c_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, l, info;
    static doublereal temp;
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

#line 210 "ssyrk.f"
    /* Parameter adjustments */
#line 210 "ssyrk.f"
    a_dim1 = *lda;
#line 210 "ssyrk.f"
    a_offset = 1 + a_dim1;
#line 210 "ssyrk.f"
    a -= a_offset;
#line 210 "ssyrk.f"
    c_dim1 = *ldc;
#line 210 "ssyrk.f"
    c_offset = 1 + c_dim1;
#line 210 "ssyrk.f"
    c__ -= c_offset;
#line 210 "ssyrk.f"

#line 210 "ssyrk.f"
    /* Function Body */
#line 210 "ssyrk.f"
    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {
#line 211 "ssyrk.f"
	nrowa = *n;
#line 212 "ssyrk.f"
    } else {
#line 213 "ssyrk.f"
	nrowa = *k;
#line 214 "ssyrk.f"
    }
#line 215 "ssyrk.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);

#line 217 "ssyrk.f"
    info = 0;
#line 218 "ssyrk.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 219 "ssyrk.f"
	info = 1;
#line 220 "ssyrk.f"
    } else if (! lsame_(trans, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, 
	    "T", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, "C", (ftnlen)1, (
	    ftnlen)1)) {
#line 223 "ssyrk.f"
	info = 2;
#line 224 "ssyrk.f"
    } else if (*n < 0) {
#line 225 "ssyrk.f"
	info = 3;
#line 226 "ssyrk.f"
    } else if (*k < 0) {
#line 227 "ssyrk.f"
	info = 4;
#line 228 "ssyrk.f"
    } else if (*lda < max(1,nrowa)) {
#line 229 "ssyrk.f"
	info = 7;
#line 230 "ssyrk.f"
    } else if (*ldc < max(1,*n)) {
#line 231 "ssyrk.f"
	info = 10;
#line 232 "ssyrk.f"
    }
#line 233 "ssyrk.f"
    if (info != 0) {
#line 234 "ssyrk.f"
	xerbla_("SSYRK ", &info, (ftnlen)6);
#line 235 "ssyrk.f"
	return 0;
#line 236 "ssyrk.f"
    }

/*     Quick return if possible. */

#line 240 "ssyrk.f"
    if (*n == 0 || (*alpha == 0. || *k == 0) && *beta == 1.) {
#line 240 "ssyrk.f"
	return 0;
#line 240 "ssyrk.f"
    }

/*     And when  alpha.eq.zero. */

#line 245 "ssyrk.f"
    if (*alpha == 0.) {
#line 246 "ssyrk.f"
	if (upper) {
#line 247 "ssyrk.f"
	    if (*beta == 0.) {
#line 248 "ssyrk.f"
		i__1 = *n;
#line 248 "ssyrk.f"
		for (j = 1; j <= i__1; ++j) {
#line 249 "ssyrk.f"
		    i__2 = j;
#line 249 "ssyrk.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 250 "ssyrk.f"
			c__[i__ + j * c_dim1] = 0.;
#line 251 "ssyrk.f"
/* L10: */
#line 251 "ssyrk.f"
		    }
#line 252 "ssyrk.f"
/* L20: */
#line 252 "ssyrk.f"
		}
#line 253 "ssyrk.f"
	    } else {
#line 254 "ssyrk.f"
		i__1 = *n;
#line 254 "ssyrk.f"
		for (j = 1; j <= i__1; ++j) {
#line 255 "ssyrk.f"
		    i__2 = j;
#line 255 "ssyrk.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 256 "ssyrk.f"
			c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
#line 257 "ssyrk.f"
/* L30: */
#line 257 "ssyrk.f"
		    }
#line 258 "ssyrk.f"
/* L40: */
#line 258 "ssyrk.f"
		}
#line 259 "ssyrk.f"
	    }
#line 260 "ssyrk.f"
	} else {
#line 261 "ssyrk.f"
	    if (*beta == 0.) {
#line 262 "ssyrk.f"
		i__1 = *n;
#line 262 "ssyrk.f"
		for (j = 1; j <= i__1; ++j) {
#line 263 "ssyrk.f"
		    i__2 = *n;
#line 263 "ssyrk.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 264 "ssyrk.f"
			c__[i__ + j * c_dim1] = 0.;
#line 265 "ssyrk.f"
/* L50: */
#line 265 "ssyrk.f"
		    }
#line 266 "ssyrk.f"
/* L60: */
#line 266 "ssyrk.f"
		}
#line 267 "ssyrk.f"
	    } else {
#line 268 "ssyrk.f"
		i__1 = *n;
#line 268 "ssyrk.f"
		for (j = 1; j <= i__1; ++j) {
#line 269 "ssyrk.f"
		    i__2 = *n;
#line 269 "ssyrk.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 270 "ssyrk.f"
			c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
#line 271 "ssyrk.f"
/* L70: */
#line 271 "ssyrk.f"
		    }
#line 272 "ssyrk.f"
/* L80: */
#line 272 "ssyrk.f"
		}
#line 273 "ssyrk.f"
	    }
#line 274 "ssyrk.f"
	}
#line 275 "ssyrk.f"
	return 0;
#line 276 "ssyrk.f"
    }

/*     Start the operations. */

#line 280 "ssyrk.f"
    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {

/*        Form  C := alpha*A*A**T + beta*C. */

#line 284 "ssyrk.f"
	if (upper) {
#line 285 "ssyrk.f"
	    i__1 = *n;
#line 285 "ssyrk.f"
	    for (j = 1; j <= i__1; ++j) {
#line 286 "ssyrk.f"
		if (*beta == 0.) {
#line 287 "ssyrk.f"
		    i__2 = j;
#line 287 "ssyrk.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 288 "ssyrk.f"
			c__[i__ + j * c_dim1] = 0.;
#line 289 "ssyrk.f"
/* L90: */
#line 289 "ssyrk.f"
		    }
#line 290 "ssyrk.f"
		} else if (*beta != 1.) {
#line 291 "ssyrk.f"
		    i__2 = j;
#line 291 "ssyrk.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 292 "ssyrk.f"
			c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
#line 293 "ssyrk.f"
/* L100: */
#line 293 "ssyrk.f"
		    }
#line 294 "ssyrk.f"
		}
#line 295 "ssyrk.f"
		i__2 = *k;
#line 295 "ssyrk.f"
		for (l = 1; l <= i__2; ++l) {
#line 296 "ssyrk.f"
		    if (a[j + l * a_dim1] != 0.) {
#line 297 "ssyrk.f"
			temp = *alpha * a[j + l * a_dim1];
#line 298 "ssyrk.f"
			i__3 = j;
#line 298 "ssyrk.f"
			for (i__ = 1; i__ <= i__3; ++i__) {
#line 299 "ssyrk.f"
			    c__[i__ + j * c_dim1] += temp * a[i__ + l * 
				    a_dim1];
#line 300 "ssyrk.f"
/* L110: */
#line 300 "ssyrk.f"
			}
#line 301 "ssyrk.f"
		    }
#line 302 "ssyrk.f"
/* L120: */
#line 302 "ssyrk.f"
		}
#line 303 "ssyrk.f"
/* L130: */
#line 303 "ssyrk.f"
	    }
#line 304 "ssyrk.f"
	} else {
#line 305 "ssyrk.f"
	    i__1 = *n;
#line 305 "ssyrk.f"
	    for (j = 1; j <= i__1; ++j) {
#line 306 "ssyrk.f"
		if (*beta == 0.) {
#line 307 "ssyrk.f"
		    i__2 = *n;
#line 307 "ssyrk.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 308 "ssyrk.f"
			c__[i__ + j * c_dim1] = 0.;
#line 309 "ssyrk.f"
/* L140: */
#line 309 "ssyrk.f"
		    }
#line 310 "ssyrk.f"
		} else if (*beta != 1.) {
#line 311 "ssyrk.f"
		    i__2 = *n;
#line 311 "ssyrk.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 312 "ssyrk.f"
			c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
#line 313 "ssyrk.f"
/* L150: */
#line 313 "ssyrk.f"
		    }
#line 314 "ssyrk.f"
		}
#line 315 "ssyrk.f"
		i__2 = *k;
#line 315 "ssyrk.f"
		for (l = 1; l <= i__2; ++l) {
#line 316 "ssyrk.f"
		    if (a[j + l * a_dim1] != 0.) {
#line 317 "ssyrk.f"
			temp = *alpha * a[j + l * a_dim1];
#line 318 "ssyrk.f"
			i__3 = *n;
#line 318 "ssyrk.f"
			for (i__ = j; i__ <= i__3; ++i__) {
#line 319 "ssyrk.f"
			    c__[i__ + j * c_dim1] += temp * a[i__ + l * 
				    a_dim1];
#line 320 "ssyrk.f"
/* L160: */
#line 320 "ssyrk.f"
			}
#line 321 "ssyrk.f"
		    }
#line 322 "ssyrk.f"
/* L170: */
#line 322 "ssyrk.f"
		}
#line 323 "ssyrk.f"
/* L180: */
#line 323 "ssyrk.f"
	    }
#line 324 "ssyrk.f"
	}
#line 325 "ssyrk.f"
    } else {

/*        Form  C := alpha*A**T*A + beta*C. */

#line 329 "ssyrk.f"
	if (upper) {
#line 330 "ssyrk.f"
	    i__1 = *n;
#line 330 "ssyrk.f"
	    for (j = 1; j <= i__1; ++j) {
#line 331 "ssyrk.f"
		i__2 = j;
#line 331 "ssyrk.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 332 "ssyrk.f"
		    temp = 0.;
#line 333 "ssyrk.f"
		    i__3 = *k;
#line 333 "ssyrk.f"
		    for (l = 1; l <= i__3; ++l) {
#line 334 "ssyrk.f"
			temp += a[l + i__ * a_dim1] * a[l + j * a_dim1];
#line 335 "ssyrk.f"
/* L190: */
#line 335 "ssyrk.f"
		    }
#line 336 "ssyrk.f"
		    if (*beta == 0.) {
#line 337 "ssyrk.f"
			c__[i__ + j * c_dim1] = *alpha * temp;
#line 338 "ssyrk.f"
		    } else {
#line 339 "ssyrk.f"
			c__[i__ + j * c_dim1] = *alpha * temp + *beta * c__[
				i__ + j * c_dim1];
#line 340 "ssyrk.f"
		    }
#line 341 "ssyrk.f"
/* L200: */
#line 341 "ssyrk.f"
		}
#line 342 "ssyrk.f"
/* L210: */
#line 342 "ssyrk.f"
	    }
#line 343 "ssyrk.f"
	} else {
#line 344 "ssyrk.f"
	    i__1 = *n;
#line 344 "ssyrk.f"
	    for (j = 1; j <= i__1; ++j) {
#line 345 "ssyrk.f"
		i__2 = *n;
#line 345 "ssyrk.f"
		for (i__ = j; i__ <= i__2; ++i__) {
#line 346 "ssyrk.f"
		    temp = 0.;
#line 347 "ssyrk.f"
		    i__3 = *k;
#line 347 "ssyrk.f"
		    for (l = 1; l <= i__3; ++l) {
#line 348 "ssyrk.f"
			temp += a[l + i__ * a_dim1] * a[l + j * a_dim1];
#line 349 "ssyrk.f"
/* L220: */
#line 349 "ssyrk.f"
		    }
#line 350 "ssyrk.f"
		    if (*beta == 0.) {
#line 351 "ssyrk.f"
			c__[i__ + j * c_dim1] = *alpha * temp;
#line 352 "ssyrk.f"
		    } else {
#line 353 "ssyrk.f"
			c__[i__ + j * c_dim1] = *alpha * temp + *beta * c__[
				i__ + j * c_dim1];
#line 354 "ssyrk.f"
		    }
#line 355 "ssyrk.f"
/* L230: */
#line 355 "ssyrk.f"
		}
#line 356 "ssyrk.f"
/* L240: */
#line 356 "ssyrk.f"
	    }
#line 357 "ssyrk.f"
	}
#line 358 "ssyrk.f"
    }

#line 360 "ssyrk.f"
    return 0;

/*     End of SSYRK . */

} /* ssyrk_ */

