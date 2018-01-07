#line 1 "dgemm.f"
/* dgemm.f -- translated by f2c (version 20100827).
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

#line 1 "dgemm.f"
/* > \brief \b DGEMM */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC) */

/*       .. Scalar Arguments .. */
/*       DOUBLE PRECISION ALPHA,BETA */
/*       INTEGER K,LDA,LDB,LDC,M,N */
/*       CHARACTER TRANSA,TRANSB */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DGEMM  performs one of the matrix-matrix operations */
/* > */
/* >    C := alpha*op( A )*op( B ) + beta*C, */
/* > */
/* > where  op( X ) is one of */
/* > */
/* >    op( X ) = X   or   op( X ) = X**T, */
/* > */
/* > alpha and beta are scalars, and A, B and C are matrices, with op( A ) */
/* > an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] TRANSA */
/* > \verbatim */
/* >          TRANSA is CHARACTER*1 */
/* >           On entry, TRANSA specifies the form of op( A ) to be used in */
/* >           the matrix multiplication as follows: */
/* > */
/* >              TRANSA = 'N' or 'n',  op( A ) = A. */
/* > */
/* >              TRANSA = 'T' or 't',  op( A ) = A**T. */
/* > */
/* >              TRANSA = 'C' or 'c',  op( A ) = A**T. */
/* > \endverbatim */
/* > */
/* > \param[in] TRANSB */
/* > \verbatim */
/* >          TRANSB is CHARACTER*1 */
/* >           On entry, TRANSB specifies the form of op( B ) to be used in */
/* >           the matrix multiplication as follows: */
/* > */
/* >              TRANSB = 'N' or 'n',  op( B ) = B. */
/* > */
/* >              TRANSB = 'T' or 't',  op( B ) = B**T. */
/* > */
/* >              TRANSB = 'C' or 'c',  op( B ) = B**T. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >           On entry,  M  specifies  the number  of rows  of the  matrix */
/* >           op( A )  and of the  matrix  C.  M  must  be at least  zero. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >           On entry,  N  specifies the number  of columns of the matrix */
/* >           op( B ) and the number of columns of the matrix C. N must be */
/* >           at least zero. */
/* > \endverbatim */
/* > */
/* > \param[in] K */
/* > \verbatim */
/* >          K is INTEGER */
/* >           On entry,  K  specifies  the number of columns of the matrix */
/* >           op( A ) and the number of rows of the matrix op( B ). K must */
/* >           be at least  zero. */
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
/* >           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise. */
/* >           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k */
/* >           part of the array  A  must contain the matrix  A,  otherwise */
/* >           the leading  k by m  part of the array  A  must contain  the */
/* >           matrix A. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >           On entry, LDA specifies the first dimension of A as declared */
/* >           in the calling (sub) program. When  TRANSA = 'N' or 'n' then */
/* >           LDA must be at least  max( 1, m ), otherwise  LDA must be at */
/* >           least  max( 1, k ). */
/* > \endverbatim */
/* > */
/* > \param[in] B */
/* > \verbatim */
/* >          B is DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is */
/* >           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise. */
/* >           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n */
/* >           part of the array  B  must contain the matrix  B,  otherwise */
/* >           the leading  n by k  part of the array  B  must contain  the */
/* >           matrix B. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >           On entry, LDB specifies the first dimension of B as declared */
/* >           in the calling (sub) program. When  TRANSB = 'N' or 'n' then */
/* >           LDB must be at least  max( 1, k ), otherwise  LDB must be at */
/* >           least  max( 1, n ). */
/* > \endverbatim */
/* > */
/* > \param[in] BETA */
/* > \verbatim */
/* >          BETA is DOUBLE PRECISION. */
/* >           On entry,  BETA  specifies the scalar  beta.  When  BETA  is */
/* >           supplied as zero then C need not be set on input. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* >          C is DOUBLE PRECISION array of DIMENSION ( LDC, n ). */
/* >           Before entry, the leading  m by n  part of the array  C must */
/* >           contain the matrix  C,  except when  beta  is zero, in which */
/* >           case C need not be set on entry. */
/* >           On exit, the array  C  is overwritten by the  m by n  matrix */
/* >           ( alpha*op( A )*op( B ) + beta*C ). */
/* > \endverbatim */
/* > */
/* > \param[in] LDC */
/* > \verbatim */
/* >          LDC is INTEGER */
/* >           On entry, LDC specifies the first dimension of C as declared */
/* >           in  the  calling  (sub)  program.   LDC  must  be  at  least */
/* >           max( 1, m ). */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2015 */

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
/* Subroutine */ int dgemm_(char *transa, char *transb, integer *m, integer *
	n, integer *k, doublereal *alpha, doublereal *a, integer *lda, 
	doublereal *b, integer *ldb, doublereal *beta, doublereal *c__, 
	integer *ldc, ftnlen transa_len, ftnlen transb_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2, 
	    i__3;

    /* Local variables */
    static integer i__, j, l, info;
    static logical nota, notb;
    static doublereal temp;
    static integer ncola;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer nrowa, nrowb;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);


/*  -- Reference BLAS level3 routine (version 3.6.0) -- */
/*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2015 */

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

/*     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not */
/*     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows */
/*     and  columns of  A  and the  number of  rows  of  B  respectively. */

#line 230 "dgemm.f"
    /* Parameter adjustments */
#line 230 "dgemm.f"
    a_dim1 = *lda;
#line 230 "dgemm.f"
    a_offset = 1 + a_dim1;
#line 230 "dgemm.f"
    a -= a_offset;
#line 230 "dgemm.f"
    b_dim1 = *ldb;
#line 230 "dgemm.f"
    b_offset = 1 + b_dim1;
#line 230 "dgemm.f"
    b -= b_offset;
#line 230 "dgemm.f"
    c_dim1 = *ldc;
#line 230 "dgemm.f"
    c_offset = 1 + c_dim1;
#line 230 "dgemm.f"
    c__ -= c_offset;
#line 230 "dgemm.f"

#line 230 "dgemm.f"
    /* Function Body */
#line 230 "dgemm.f"
    nota = lsame_(transa, "N", (ftnlen)1, (ftnlen)1);
#line 231 "dgemm.f"
    notb = lsame_(transb, "N", (ftnlen)1, (ftnlen)1);
#line 232 "dgemm.f"
    if (nota) {
#line 233 "dgemm.f"
	nrowa = *m;
#line 234 "dgemm.f"
	ncola = *k;
#line 235 "dgemm.f"
    } else {
#line 236 "dgemm.f"
	nrowa = *k;
#line 237 "dgemm.f"
	ncola = *m;
#line 238 "dgemm.f"
    }
#line 239 "dgemm.f"
    if (notb) {
#line 240 "dgemm.f"
	nrowb = *k;
#line 241 "dgemm.f"
    } else {
#line 242 "dgemm.f"
	nrowb = *n;
#line 243 "dgemm.f"
    }

/*     Test the input parameters. */

#line 247 "dgemm.f"
    info = 0;
#line 248 "dgemm.f"
    if (! nota && ! lsame_(transa, "C", (ftnlen)1, (ftnlen)1) && ! lsame_(
	    transa, "T", (ftnlen)1, (ftnlen)1)) {
#line 250 "dgemm.f"
	info = 1;
#line 251 "dgemm.f"
    } else if (! notb && ! lsame_(transb, "C", (ftnlen)1, (ftnlen)1) && ! 
	    lsame_(transb, "T", (ftnlen)1, (ftnlen)1)) {
#line 253 "dgemm.f"
	info = 2;
#line 254 "dgemm.f"
    } else if (*m < 0) {
#line 255 "dgemm.f"
	info = 3;
#line 256 "dgemm.f"
    } else if (*n < 0) {
#line 257 "dgemm.f"
	info = 4;
#line 258 "dgemm.f"
    } else if (*k < 0) {
#line 259 "dgemm.f"
	info = 5;
#line 260 "dgemm.f"
    } else if (*lda < max(1,nrowa)) {
#line 261 "dgemm.f"
	info = 8;
#line 262 "dgemm.f"
    } else if (*ldb < max(1,nrowb)) {
#line 263 "dgemm.f"
	info = 10;
#line 264 "dgemm.f"
    } else if (*ldc < max(1,*m)) {
#line 265 "dgemm.f"
	info = 13;
#line 266 "dgemm.f"
    }
#line 267 "dgemm.f"
    if (info != 0) {
#line 268 "dgemm.f"
	xerbla_("DGEMM ", &info, (ftnlen)6);
#line 269 "dgemm.f"
	return 0;
#line 270 "dgemm.f"
    }

/*     Quick return if possible. */

#line 274 "dgemm.f"
    if (*m == 0 || *n == 0 || (*alpha == 0. || *k == 0) && *beta == 1.) {
#line 274 "dgemm.f"
	return 0;
#line 274 "dgemm.f"
    }

/*     And if  alpha.eq.zero. */

#line 279 "dgemm.f"
    if (*alpha == 0.) {
#line 280 "dgemm.f"
	if (*beta == 0.) {
#line 281 "dgemm.f"
	    i__1 = *n;
#line 281 "dgemm.f"
	    for (j = 1; j <= i__1; ++j) {
#line 282 "dgemm.f"
		i__2 = *m;
#line 282 "dgemm.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 283 "dgemm.f"
		    c__[i__ + j * c_dim1] = 0.;
#line 284 "dgemm.f"
/* L10: */
#line 284 "dgemm.f"
		}
#line 285 "dgemm.f"
/* L20: */
#line 285 "dgemm.f"
	    }
#line 286 "dgemm.f"
	} else {
#line 287 "dgemm.f"
	    i__1 = *n;
#line 287 "dgemm.f"
	    for (j = 1; j <= i__1; ++j) {
#line 288 "dgemm.f"
		i__2 = *m;
#line 288 "dgemm.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 289 "dgemm.f"
		    c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
#line 290 "dgemm.f"
/* L30: */
#line 290 "dgemm.f"
		}
#line 291 "dgemm.f"
/* L40: */
#line 291 "dgemm.f"
	    }
#line 292 "dgemm.f"
	}
#line 293 "dgemm.f"
	return 0;
#line 294 "dgemm.f"
    }

/*     Start the operations. */

#line 298 "dgemm.f"
    if (notb) {
#line 299 "dgemm.f"
	if (nota) {

/*           Form  C := alpha*A*B + beta*C. */

#line 303 "dgemm.f"
	    i__1 = *n;
#line 303 "dgemm.f"
	    for (j = 1; j <= i__1; ++j) {
#line 304 "dgemm.f"
		if (*beta == 0.) {
#line 305 "dgemm.f"
		    i__2 = *m;
#line 305 "dgemm.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 306 "dgemm.f"
			c__[i__ + j * c_dim1] = 0.;
#line 307 "dgemm.f"
/* L50: */
#line 307 "dgemm.f"
		    }
#line 308 "dgemm.f"
		} else if (*beta != 1.) {
#line 309 "dgemm.f"
		    i__2 = *m;
#line 309 "dgemm.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 310 "dgemm.f"
			c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
#line 311 "dgemm.f"
/* L60: */
#line 311 "dgemm.f"
		    }
#line 312 "dgemm.f"
		}
#line 313 "dgemm.f"
		i__2 = *k;
#line 313 "dgemm.f"
		for (l = 1; l <= i__2; ++l) {
#line 314 "dgemm.f"
		    temp = *alpha * b[l + j * b_dim1];
#line 315 "dgemm.f"
		    i__3 = *m;
#line 315 "dgemm.f"
		    for (i__ = 1; i__ <= i__3; ++i__) {
#line 316 "dgemm.f"
			c__[i__ + j * c_dim1] += temp * a[i__ + l * a_dim1];
#line 317 "dgemm.f"
/* L70: */
#line 317 "dgemm.f"
		    }
#line 318 "dgemm.f"
/* L80: */
#line 318 "dgemm.f"
		}
#line 319 "dgemm.f"
/* L90: */
#line 319 "dgemm.f"
	    }
#line 320 "dgemm.f"
	} else {

/*           Form  C := alpha*A**T*B + beta*C */

#line 324 "dgemm.f"
	    i__1 = *n;
#line 324 "dgemm.f"
	    for (j = 1; j <= i__1; ++j) {
#line 325 "dgemm.f"
		i__2 = *m;
#line 325 "dgemm.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 326 "dgemm.f"
		    temp = 0.;
#line 327 "dgemm.f"
		    i__3 = *k;
#line 327 "dgemm.f"
		    for (l = 1; l <= i__3; ++l) {
#line 328 "dgemm.f"
			temp += a[l + i__ * a_dim1] * b[l + j * b_dim1];
#line 329 "dgemm.f"
/* L100: */
#line 329 "dgemm.f"
		    }
#line 330 "dgemm.f"
		    if (*beta == 0.) {
#line 331 "dgemm.f"
			c__[i__ + j * c_dim1] = *alpha * temp;
#line 332 "dgemm.f"
		    } else {
#line 333 "dgemm.f"
			c__[i__ + j * c_dim1] = *alpha * temp + *beta * c__[
				i__ + j * c_dim1];
#line 334 "dgemm.f"
		    }
#line 335 "dgemm.f"
/* L110: */
#line 335 "dgemm.f"
		}
#line 336 "dgemm.f"
/* L120: */
#line 336 "dgemm.f"
	    }
#line 337 "dgemm.f"
	}
#line 338 "dgemm.f"
    } else {
#line 339 "dgemm.f"
	if (nota) {

/*           Form  C := alpha*A*B**T + beta*C */

#line 343 "dgemm.f"
	    i__1 = *n;
#line 343 "dgemm.f"
	    for (j = 1; j <= i__1; ++j) {
#line 344 "dgemm.f"
		if (*beta == 0.) {
#line 345 "dgemm.f"
		    i__2 = *m;
#line 345 "dgemm.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 346 "dgemm.f"
			c__[i__ + j * c_dim1] = 0.;
#line 347 "dgemm.f"
/* L130: */
#line 347 "dgemm.f"
		    }
#line 348 "dgemm.f"
		} else if (*beta != 1.) {
#line 349 "dgemm.f"
		    i__2 = *m;
#line 349 "dgemm.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 350 "dgemm.f"
			c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
#line 351 "dgemm.f"
/* L140: */
#line 351 "dgemm.f"
		    }
#line 352 "dgemm.f"
		}
#line 353 "dgemm.f"
		i__2 = *k;
#line 353 "dgemm.f"
		for (l = 1; l <= i__2; ++l) {
#line 354 "dgemm.f"
		    temp = *alpha * b[j + l * b_dim1];
#line 355 "dgemm.f"
		    i__3 = *m;
#line 355 "dgemm.f"
		    for (i__ = 1; i__ <= i__3; ++i__) {
#line 356 "dgemm.f"
			c__[i__ + j * c_dim1] += temp * a[i__ + l * a_dim1];
#line 357 "dgemm.f"
/* L150: */
#line 357 "dgemm.f"
		    }
#line 358 "dgemm.f"
/* L160: */
#line 358 "dgemm.f"
		}
#line 359 "dgemm.f"
/* L170: */
#line 359 "dgemm.f"
	    }
#line 360 "dgemm.f"
	} else {

/*           Form  C := alpha*A**T*B**T + beta*C */

#line 364 "dgemm.f"
	    i__1 = *n;
#line 364 "dgemm.f"
	    for (j = 1; j <= i__1; ++j) {
#line 365 "dgemm.f"
		i__2 = *m;
#line 365 "dgemm.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 366 "dgemm.f"
		    temp = 0.;
#line 367 "dgemm.f"
		    i__3 = *k;
#line 367 "dgemm.f"
		    for (l = 1; l <= i__3; ++l) {
#line 368 "dgemm.f"
			temp += a[l + i__ * a_dim1] * b[j + l * b_dim1];
#line 369 "dgemm.f"
/* L180: */
#line 369 "dgemm.f"
		    }
#line 370 "dgemm.f"
		    if (*beta == 0.) {
#line 371 "dgemm.f"
			c__[i__ + j * c_dim1] = *alpha * temp;
#line 372 "dgemm.f"
		    } else {
#line 373 "dgemm.f"
			c__[i__ + j * c_dim1] = *alpha * temp + *beta * c__[
				i__ + j * c_dim1];
#line 374 "dgemm.f"
		    }
#line 375 "dgemm.f"
/* L190: */
#line 375 "dgemm.f"
		}
#line 376 "dgemm.f"
/* L200: */
#line 376 "dgemm.f"
	    }
#line 377 "dgemm.f"
	}
#line 378 "dgemm.f"
    }

#line 380 "dgemm.f"
    return 0;

/*     End of DGEMM . */

} /* dgemm_ */

