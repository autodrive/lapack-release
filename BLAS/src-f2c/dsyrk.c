#line 1 "dsyrk.f"
/* dsyrk.f -- translated by f2c (version 20100827).
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

#line 1 "dsyrk.f"
/* > \brief \b DSYRK */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC) */

/*       .. Scalar Arguments .. */
/*       DOUBLE PRECISION ALPHA,BETA */
/*       INTEGER K,LDA,LDC,N */
/*       CHARACTER TRANS,UPLO */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION A(LDA,*),C(LDC,*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DSYRK  performs one of the symmetric rank k operations */
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
/* >          ALPHA is DOUBLE PRECISION. */
/* >           On entry, ALPHA specifies the scalar alpha. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is */
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
/* >          BETA is DOUBLE PRECISION. */
/* >           On entry, BETA specifies the scalar beta. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* >          C is DOUBLE PRECISION array of DIMENSION ( LDC, n ). */
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
/* Subroutine */ int dsyrk_(char *uplo, char *trans, integer *n, integer *k, 
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

#line 210 "dsyrk.f"
    /* Parameter adjustments */
#line 210 "dsyrk.f"
    a_dim1 = *lda;
#line 210 "dsyrk.f"
    a_offset = 1 + a_dim1;
#line 210 "dsyrk.f"
    a -= a_offset;
#line 210 "dsyrk.f"
    c_dim1 = *ldc;
#line 210 "dsyrk.f"
    c_offset = 1 + c_dim1;
#line 210 "dsyrk.f"
    c__ -= c_offset;
#line 210 "dsyrk.f"

#line 210 "dsyrk.f"
    /* Function Body */
#line 210 "dsyrk.f"
    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {
#line 211 "dsyrk.f"
	nrowa = *n;
#line 212 "dsyrk.f"
    } else {
#line 213 "dsyrk.f"
	nrowa = *k;
#line 214 "dsyrk.f"
    }
#line 215 "dsyrk.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);

#line 217 "dsyrk.f"
    info = 0;
#line 218 "dsyrk.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 219 "dsyrk.f"
	info = 1;
#line 220 "dsyrk.f"
    } else if (! lsame_(trans, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, 
	    "T", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, "C", (ftnlen)1, (
	    ftnlen)1)) {
#line 223 "dsyrk.f"
	info = 2;
#line 224 "dsyrk.f"
    } else if (*n < 0) {
#line 225 "dsyrk.f"
	info = 3;
#line 226 "dsyrk.f"
    } else if (*k < 0) {
#line 227 "dsyrk.f"
	info = 4;
#line 228 "dsyrk.f"
    } else if (*lda < max(1,nrowa)) {
#line 229 "dsyrk.f"
	info = 7;
#line 230 "dsyrk.f"
    } else if (*ldc < max(1,*n)) {
#line 231 "dsyrk.f"
	info = 10;
#line 232 "dsyrk.f"
    }
#line 233 "dsyrk.f"
    if (info != 0) {
#line 234 "dsyrk.f"
	xerbla_("DSYRK ", &info, (ftnlen)6);
#line 235 "dsyrk.f"
	return 0;
#line 236 "dsyrk.f"
    }

/*     Quick return if possible. */

#line 240 "dsyrk.f"
    if (*n == 0 || (*alpha == 0. || *k == 0) && *beta == 1.) {
#line 240 "dsyrk.f"
	return 0;
#line 240 "dsyrk.f"
    }

/*     And when  alpha.eq.zero. */

#line 245 "dsyrk.f"
    if (*alpha == 0.) {
#line 246 "dsyrk.f"
	if (upper) {
#line 247 "dsyrk.f"
	    if (*beta == 0.) {
#line 248 "dsyrk.f"
		i__1 = *n;
#line 248 "dsyrk.f"
		for (j = 1; j <= i__1; ++j) {
#line 249 "dsyrk.f"
		    i__2 = j;
#line 249 "dsyrk.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 250 "dsyrk.f"
			c__[i__ + j * c_dim1] = 0.;
#line 251 "dsyrk.f"
/* L10: */
#line 251 "dsyrk.f"
		    }
#line 252 "dsyrk.f"
/* L20: */
#line 252 "dsyrk.f"
		}
#line 253 "dsyrk.f"
	    } else {
#line 254 "dsyrk.f"
		i__1 = *n;
#line 254 "dsyrk.f"
		for (j = 1; j <= i__1; ++j) {
#line 255 "dsyrk.f"
		    i__2 = j;
#line 255 "dsyrk.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 256 "dsyrk.f"
			c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
#line 257 "dsyrk.f"
/* L30: */
#line 257 "dsyrk.f"
		    }
#line 258 "dsyrk.f"
/* L40: */
#line 258 "dsyrk.f"
		}
#line 259 "dsyrk.f"
	    }
#line 260 "dsyrk.f"
	} else {
#line 261 "dsyrk.f"
	    if (*beta == 0.) {
#line 262 "dsyrk.f"
		i__1 = *n;
#line 262 "dsyrk.f"
		for (j = 1; j <= i__1; ++j) {
#line 263 "dsyrk.f"
		    i__2 = *n;
#line 263 "dsyrk.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 264 "dsyrk.f"
			c__[i__ + j * c_dim1] = 0.;
#line 265 "dsyrk.f"
/* L50: */
#line 265 "dsyrk.f"
		    }
#line 266 "dsyrk.f"
/* L60: */
#line 266 "dsyrk.f"
		}
#line 267 "dsyrk.f"
	    } else {
#line 268 "dsyrk.f"
		i__1 = *n;
#line 268 "dsyrk.f"
		for (j = 1; j <= i__1; ++j) {
#line 269 "dsyrk.f"
		    i__2 = *n;
#line 269 "dsyrk.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 270 "dsyrk.f"
			c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
#line 271 "dsyrk.f"
/* L70: */
#line 271 "dsyrk.f"
		    }
#line 272 "dsyrk.f"
/* L80: */
#line 272 "dsyrk.f"
		}
#line 273 "dsyrk.f"
	    }
#line 274 "dsyrk.f"
	}
#line 275 "dsyrk.f"
	return 0;
#line 276 "dsyrk.f"
    }

/*     Start the operations. */

#line 280 "dsyrk.f"
    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {

/*        Form  C := alpha*A*A**T + beta*C. */

#line 284 "dsyrk.f"
	if (upper) {
#line 285 "dsyrk.f"
	    i__1 = *n;
#line 285 "dsyrk.f"
	    for (j = 1; j <= i__1; ++j) {
#line 286 "dsyrk.f"
		if (*beta == 0.) {
#line 287 "dsyrk.f"
		    i__2 = j;
#line 287 "dsyrk.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 288 "dsyrk.f"
			c__[i__ + j * c_dim1] = 0.;
#line 289 "dsyrk.f"
/* L90: */
#line 289 "dsyrk.f"
		    }
#line 290 "dsyrk.f"
		} else if (*beta != 1.) {
#line 291 "dsyrk.f"
		    i__2 = j;
#line 291 "dsyrk.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 292 "dsyrk.f"
			c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
#line 293 "dsyrk.f"
/* L100: */
#line 293 "dsyrk.f"
		    }
#line 294 "dsyrk.f"
		}
#line 295 "dsyrk.f"
		i__2 = *k;
#line 295 "dsyrk.f"
		for (l = 1; l <= i__2; ++l) {
#line 296 "dsyrk.f"
		    if (a[j + l * a_dim1] != 0.) {
#line 297 "dsyrk.f"
			temp = *alpha * a[j + l * a_dim1];
#line 298 "dsyrk.f"
			i__3 = j;
#line 298 "dsyrk.f"
			for (i__ = 1; i__ <= i__3; ++i__) {
#line 299 "dsyrk.f"
			    c__[i__ + j * c_dim1] += temp * a[i__ + l * 
				    a_dim1];
#line 300 "dsyrk.f"
/* L110: */
#line 300 "dsyrk.f"
			}
#line 301 "dsyrk.f"
		    }
#line 302 "dsyrk.f"
/* L120: */
#line 302 "dsyrk.f"
		}
#line 303 "dsyrk.f"
/* L130: */
#line 303 "dsyrk.f"
	    }
#line 304 "dsyrk.f"
	} else {
#line 305 "dsyrk.f"
	    i__1 = *n;
#line 305 "dsyrk.f"
	    for (j = 1; j <= i__1; ++j) {
#line 306 "dsyrk.f"
		if (*beta == 0.) {
#line 307 "dsyrk.f"
		    i__2 = *n;
#line 307 "dsyrk.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 308 "dsyrk.f"
			c__[i__ + j * c_dim1] = 0.;
#line 309 "dsyrk.f"
/* L140: */
#line 309 "dsyrk.f"
		    }
#line 310 "dsyrk.f"
		} else if (*beta != 1.) {
#line 311 "dsyrk.f"
		    i__2 = *n;
#line 311 "dsyrk.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 312 "dsyrk.f"
			c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
#line 313 "dsyrk.f"
/* L150: */
#line 313 "dsyrk.f"
		    }
#line 314 "dsyrk.f"
		}
#line 315 "dsyrk.f"
		i__2 = *k;
#line 315 "dsyrk.f"
		for (l = 1; l <= i__2; ++l) {
#line 316 "dsyrk.f"
		    if (a[j + l * a_dim1] != 0.) {
#line 317 "dsyrk.f"
			temp = *alpha * a[j + l * a_dim1];
#line 318 "dsyrk.f"
			i__3 = *n;
#line 318 "dsyrk.f"
			for (i__ = j; i__ <= i__3; ++i__) {
#line 319 "dsyrk.f"
			    c__[i__ + j * c_dim1] += temp * a[i__ + l * 
				    a_dim1];
#line 320 "dsyrk.f"
/* L160: */
#line 320 "dsyrk.f"
			}
#line 321 "dsyrk.f"
		    }
#line 322 "dsyrk.f"
/* L170: */
#line 322 "dsyrk.f"
		}
#line 323 "dsyrk.f"
/* L180: */
#line 323 "dsyrk.f"
	    }
#line 324 "dsyrk.f"
	}
#line 325 "dsyrk.f"
    } else {

/*        Form  C := alpha*A**T*A + beta*C. */

#line 329 "dsyrk.f"
	if (upper) {
#line 330 "dsyrk.f"
	    i__1 = *n;
#line 330 "dsyrk.f"
	    for (j = 1; j <= i__1; ++j) {
#line 331 "dsyrk.f"
		i__2 = j;
#line 331 "dsyrk.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 332 "dsyrk.f"
		    temp = 0.;
#line 333 "dsyrk.f"
		    i__3 = *k;
#line 333 "dsyrk.f"
		    for (l = 1; l <= i__3; ++l) {
#line 334 "dsyrk.f"
			temp += a[l + i__ * a_dim1] * a[l + j * a_dim1];
#line 335 "dsyrk.f"
/* L190: */
#line 335 "dsyrk.f"
		    }
#line 336 "dsyrk.f"
		    if (*beta == 0.) {
#line 337 "dsyrk.f"
			c__[i__ + j * c_dim1] = *alpha * temp;
#line 338 "dsyrk.f"
		    } else {
#line 339 "dsyrk.f"
			c__[i__ + j * c_dim1] = *alpha * temp + *beta * c__[
				i__ + j * c_dim1];
#line 340 "dsyrk.f"
		    }
#line 341 "dsyrk.f"
/* L200: */
#line 341 "dsyrk.f"
		}
#line 342 "dsyrk.f"
/* L210: */
#line 342 "dsyrk.f"
	    }
#line 343 "dsyrk.f"
	} else {
#line 344 "dsyrk.f"
	    i__1 = *n;
#line 344 "dsyrk.f"
	    for (j = 1; j <= i__1; ++j) {
#line 345 "dsyrk.f"
		i__2 = *n;
#line 345 "dsyrk.f"
		for (i__ = j; i__ <= i__2; ++i__) {
#line 346 "dsyrk.f"
		    temp = 0.;
#line 347 "dsyrk.f"
		    i__3 = *k;
#line 347 "dsyrk.f"
		    for (l = 1; l <= i__3; ++l) {
#line 348 "dsyrk.f"
			temp += a[l + i__ * a_dim1] * a[l + j * a_dim1];
#line 349 "dsyrk.f"
/* L220: */
#line 349 "dsyrk.f"
		    }
#line 350 "dsyrk.f"
		    if (*beta == 0.) {
#line 351 "dsyrk.f"
			c__[i__ + j * c_dim1] = *alpha * temp;
#line 352 "dsyrk.f"
		    } else {
#line 353 "dsyrk.f"
			c__[i__ + j * c_dim1] = *alpha * temp + *beta * c__[
				i__ + j * c_dim1];
#line 354 "dsyrk.f"
		    }
#line 355 "dsyrk.f"
/* L230: */
#line 355 "dsyrk.f"
		}
#line 356 "dsyrk.f"
/* L240: */
#line 356 "dsyrk.f"
	    }
#line 357 "dsyrk.f"
	}
#line 358 "dsyrk.f"
    }

#line 360 "dsyrk.f"
    return 0;

/*     End of DSYRK . */

} /* dsyrk_ */

