#line 1 "sgemm.f"
/* sgemm.f -- translated by f2c (version 20100827).
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

#line 1 "sgemm.f"
/* > \brief \b SGEMM */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC) */

/*       .. Scalar Arguments .. */
/*       REAL ALPHA,BETA */
/*       INTEGER K,LDA,LDB,LDC,M,N */
/*       CHARACTER TRANSA,TRANSB */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL A(LDA,*),B(LDB,*),C(LDC,*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SGEMM  performs one of the matrix-matrix operations */
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
/* >          ALPHA is REAL */
/* >           On entry, ALPHA specifies the scalar alpha. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is REAL array of DIMENSION ( LDA, ka ), where ka is */
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
/* >          B is REAL array of DIMENSION ( LDB, kb ), where kb is */
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
/* >          BETA is REAL */
/* >           On entry,  BETA  specifies the scalar  beta.  When  BETA  is */
/* >           supplied as zero then C need not be set on input. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* >          C is REAL array of DIMENSION ( LDC, n ). */
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
/* Subroutine */ int sgemm_(char *transa, char *transb, integer *m, integer *
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

/*     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not */
/*     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows */
/*     and  columns of  A  and the  number of  rows  of  B  respectively. */

#line 230 "sgemm.f"
    /* Parameter adjustments */
#line 230 "sgemm.f"
    a_dim1 = *lda;
#line 230 "sgemm.f"
    a_offset = 1 + a_dim1;
#line 230 "sgemm.f"
    a -= a_offset;
#line 230 "sgemm.f"
    b_dim1 = *ldb;
#line 230 "sgemm.f"
    b_offset = 1 + b_dim1;
#line 230 "sgemm.f"
    b -= b_offset;
#line 230 "sgemm.f"
    c_dim1 = *ldc;
#line 230 "sgemm.f"
    c_offset = 1 + c_dim1;
#line 230 "sgemm.f"
    c__ -= c_offset;
#line 230 "sgemm.f"

#line 230 "sgemm.f"
    /* Function Body */
#line 230 "sgemm.f"
    nota = lsame_(transa, "N", (ftnlen)1, (ftnlen)1);
#line 231 "sgemm.f"
    notb = lsame_(transb, "N", (ftnlen)1, (ftnlen)1);
#line 232 "sgemm.f"
    if (nota) {
#line 233 "sgemm.f"
	nrowa = *m;
#line 234 "sgemm.f"
	ncola = *k;
#line 235 "sgemm.f"
    } else {
#line 236 "sgemm.f"
	nrowa = *k;
#line 237 "sgemm.f"
	ncola = *m;
#line 238 "sgemm.f"
    }
#line 239 "sgemm.f"
    if (notb) {
#line 240 "sgemm.f"
	nrowb = *k;
#line 241 "sgemm.f"
    } else {
#line 242 "sgemm.f"
	nrowb = *n;
#line 243 "sgemm.f"
    }

/*     Test the input parameters. */

#line 247 "sgemm.f"
    info = 0;
#line 248 "sgemm.f"
    if (! nota && ! lsame_(transa, "C", (ftnlen)1, (ftnlen)1) && ! lsame_(
	    transa, "T", (ftnlen)1, (ftnlen)1)) {
#line 250 "sgemm.f"
	info = 1;
#line 251 "sgemm.f"
    } else if (! notb && ! lsame_(transb, "C", (ftnlen)1, (ftnlen)1) && ! 
	    lsame_(transb, "T", (ftnlen)1, (ftnlen)1)) {
#line 253 "sgemm.f"
	info = 2;
#line 254 "sgemm.f"
    } else if (*m < 0) {
#line 255 "sgemm.f"
	info = 3;
#line 256 "sgemm.f"
    } else if (*n < 0) {
#line 257 "sgemm.f"
	info = 4;
#line 258 "sgemm.f"
    } else if (*k < 0) {
#line 259 "sgemm.f"
	info = 5;
#line 260 "sgemm.f"
    } else if (*lda < max(1,nrowa)) {
#line 261 "sgemm.f"
	info = 8;
#line 262 "sgemm.f"
    } else if (*ldb < max(1,nrowb)) {
#line 263 "sgemm.f"
	info = 10;
#line 264 "sgemm.f"
    } else if (*ldc < max(1,*m)) {
#line 265 "sgemm.f"
	info = 13;
#line 266 "sgemm.f"
    }
#line 267 "sgemm.f"
    if (info != 0) {
#line 268 "sgemm.f"
	xerbla_("SGEMM ", &info, (ftnlen)6);
#line 269 "sgemm.f"
	return 0;
#line 270 "sgemm.f"
    }

/*     Quick return if possible. */

#line 274 "sgemm.f"
    if (*m == 0 || *n == 0 || (*alpha == 0. || *k == 0) && *beta == 1.) {
#line 274 "sgemm.f"
	return 0;
#line 274 "sgemm.f"
    }

/*     And if  alpha.eq.zero. */

#line 279 "sgemm.f"
    if (*alpha == 0.) {
#line 280 "sgemm.f"
	if (*beta == 0.) {
#line 281 "sgemm.f"
	    i__1 = *n;
#line 281 "sgemm.f"
	    for (j = 1; j <= i__1; ++j) {
#line 282 "sgemm.f"
		i__2 = *m;
#line 282 "sgemm.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 283 "sgemm.f"
		    c__[i__ + j * c_dim1] = 0.;
#line 284 "sgemm.f"
/* L10: */
#line 284 "sgemm.f"
		}
#line 285 "sgemm.f"
/* L20: */
#line 285 "sgemm.f"
	    }
#line 286 "sgemm.f"
	} else {
#line 287 "sgemm.f"
	    i__1 = *n;
#line 287 "sgemm.f"
	    for (j = 1; j <= i__1; ++j) {
#line 288 "sgemm.f"
		i__2 = *m;
#line 288 "sgemm.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 289 "sgemm.f"
		    c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
#line 290 "sgemm.f"
/* L30: */
#line 290 "sgemm.f"
		}
#line 291 "sgemm.f"
/* L40: */
#line 291 "sgemm.f"
	    }
#line 292 "sgemm.f"
	}
#line 293 "sgemm.f"
	return 0;
#line 294 "sgemm.f"
    }

/*     Start the operations. */

#line 298 "sgemm.f"
    if (notb) {
#line 299 "sgemm.f"
	if (nota) {

/*           Form  C := alpha*A*B + beta*C. */

#line 303 "sgemm.f"
	    i__1 = *n;
#line 303 "sgemm.f"
	    for (j = 1; j <= i__1; ++j) {
#line 304 "sgemm.f"
		if (*beta == 0.) {
#line 305 "sgemm.f"
		    i__2 = *m;
#line 305 "sgemm.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 306 "sgemm.f"
			c__[i__ + j * c_dim1] = 0.;
#line 307 "sgemm.f"
/* L50: */
#line 307 "sgemm.f"
		    }
#line 308 "sgemm.f"
		} else if (*beta != 1.) {
#line 309 "sgemm.f"
		    i__2 = *m;
#line 309 "sgemm.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 310 "sgemm.f"
			c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
#line 311 "sgemm.f"
/* L60: */
#line 311 "sgemm.f"
		    }
#line 312 "sgemm.f"
		}
#line 313 "sgemm.f"
		i__2 = *k;
#line 313 "sgemm.f"
		for (l = 1; l <= i__2; ++l) {
#line 314 "sgemm.f"
		    temp = *alpha * b[l + j * b_dim1];
#line 315 "sgemm.f"
		    i__3 = *m;
#line 315 "sgemm.f"
		    for (i__ = 1; i__ <= i__3; ++i__) {
#line 316 "sgemm.f"
			c__[i__ + j * c_dim1] += temp * a[i__ + l * a_dim1];
#line 317 "sgemm.f"
/* L70: */
#line 317 "sgemm.f"
		    }
#line 318 "sgemm.f"
/* L80: */
#line 318 "sgemm.f"
		}
#line 319 "sgemm.f"
/* L90: */
#line 319 "sgemm.f"
	    }
#line 320 "sgemm.f"
	} else {

/*           Form  C := alpha*A**T*B + beta*C */

#line 324 "sgemm.f"
	    i__1 = *n;
#line 324 "sgemm.f"
	    for (j = 1; j <= i__1; ++j) {
#line 325 "sgemm.f"
		i__2 = *m;
#line 325 "sgemm.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 326 "sgemm.f"
		    temp = 0.;
#line 327 "sgemm.f"
		    i__3 = *k;
#line 327 "sgemm.f"
		    for (l = 1; l <= i__3; ++l) {
#line 328 "sgemm.f"
			temp += a[l + i__ * a_dim1] * b[l + j * b_dim1];
#line 329 "sgemm.f"
/* L100: */
#line 329 "sgemm.f"
		    }
#line 330 "sgemm.f"
		    if (*beta == 0.) {
#line 331 "sgemm.f"
			c__[i__ + j * c_dim1] = *alpha * temp;
#line 332 "sgemm.f"
		    } else {
#line 333 "sgemm.f"
			c__[i__ + j * c_dim1] = *alpha * temp + *beta * c__[
				i__ + j * c_dim1];
#line 334 "sgemm.f"
		    }
#line 335 "sgemm.f"
/* L110: */
#line 335 "sgemm.f"
		}
#line 336 "sgemm.f"
/* L120: */
#line 336 "sgemm.f"
	    }
#line 337 "sgemm.f"
	}
#line 338 "sgemm.f"
    } else {
#line 339 "sgemm.f"
	if (nota) {

/*           Form  C := alpha*A*B**T + beta*C */

#line 343 "sgemm.f"
	    i__1 = *n;
#line 343 "sgemm.f"
	    for (j = 1; j <= i__1; ++j) {
#line 344 "sgemm.f"
		if (*beta == 0.) {
#line 345 "sgemm.f"
		    i__2 = *m;
#line 345 "sgemm.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 346 "sgemm.f"
			c__[i__ + j * c_dim1] = 0.;
#line 347 "sgemm.f"
/* L130: */
#line 347 "sgemm.f"
		    }
#line 348 "sgemm.f"
		} else if (*beta != 1.) {
#line 349 "sgemm.f"
		    i__2 = *m;
#line 349 "sgemm.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 350 "sgemm.f"
			c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
#line 351 "sgemm.f"
/* L140: */
#line 351 "sgemm.f"
		    }
#line 352 "sgemm.f"
		}
#line 353 "sgemm.f"
		i__2 = *k;
#line 353 "sgemm.f"
		for (l = 1; l <= i__2; ++l) {
#line 354 "sgemm.f"
		    temp = *alpha * b[j + l * b_dim1];
#line 355 "sgemm.f"
		    i__3 = *m;
#line 355 "sgemm.f"
		    for (i__ = 1; i__ <= i__3; ++i__) {
#line 356 "sgemm.f"
			c__[i__ + j * c_dim1] += temp * a[i__ + l * a_dim1];
#line 357 "sgemm.f"
/* L150: */
#line 357 "sgemm.f"
		    }
#line 358 "sgemm.f"
/* L160: */
#line 358 "sgemm.f"
		}
#line 359 "sgemm.f"
/* L170: */
#line 359 "sgemm.f"
	    }
#line 360 "sgemm.f"
	} else {

/*           Form  C := alpha*A**T*B**T + beta*C */

#line 364 "sgemm.f"
	    i__1 = *n;
#line 364 "sgemm.f"
	    for (j = 1; j <= i__1; ++j) {
#line 365 "sgemm.f"
		i__2 = *m;
#line 365 "sgemm.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 366 "sgemm.f"
		    temp = 0.;
#line 367 "sgemm.f"
		    i__3 = *k;
#line 367 "sgemm.f"
		    for (l = 1; l <= i__3; ++l) {
#line 368 "sgemm.f"
			temp += a[l + i__ * a_dim1] * b[j + l * b_dim1];
#line 369 "sgemm.f"
/* L180: */
#line 369 "sgemm.f"
		    }
#line 370 "sgemm.f"
		    if (*beta == 0.) {
#line 371 "sgemm.f"
			c__[i__ + j * c_dim1] = *alpha * temp;
#line 372 "sgemm.f"
		    } else {
#line 373 "sgemm.f"
			c__[i__ + j * c_dim1] = *alpha * temp + *beta * c__[
				i__ + j * c_dim1];
#line 374 "sgemm.f"
		    }
#line 375 "sgemm.f"
/* L190: */
#line 375 "sgemm.f"
		}
#line 376 "sgemm.f"
/* L200: */
#line 376 "sgemm.f"
	    }
#line 377 "sgemm.f"
	}
#line 378 "sgemm.f"
    }

#line 380 "sgemm.f"
    return 0;

/*     End of SGEMM . */

} /* sgemm_ */

