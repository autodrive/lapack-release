#line 1 "dsymm.f"
/* dsymm.f -- translated by f2c (version 20100827).
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

#line 1 "dsymm.f"
/* > \brief \b DSYMM */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DSYMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC) */

/*       .. Scalar Arguments .. */
/*       DOUBLE PRECISION ALPHA,BETA */
/*       INTEGER LDA,LDB,LDC,M,N */
/*       CHARACTER SIDE,UPLO */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DSYMM  performs one of the matrix-matrix operations */
/* > */
/* >    C := alpha*A*B + beta*C, */
/* > */
/* > or */
/* > */
/* >    C := alpha*B*A + beta*C, */
/* > */
/* > where alpha and beta are scalars,  A is a symmetric matrix and  B and */
/* > C are  m by n matrices. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] SIDE */
/* > \verbatim */
/* >          SIDE is CHARACTER*1 */
/* >           On entry,  SIDE  specifies whether  the  symmetric matrix  A */
/* >           appears on the  left or right  in the  operation as follows: */
/* > */
/* >              SIDE = 'L' or 'l'   C := alpha*A*B + beta*C, */
/* > */
/* >              SIDE = 'R' or 'r'   C := alpha*B*A + beta*C, */
/* > \endverbatim */
/* > */
/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >           On  entry,   UPLO  specifies  whether  the  upper  or  lower */
/* >           triangular  part  of  the  symmetric  matrix   A  is  to  be */
/* >           referenced as follows: */
/* > */
/* >              UPLO = 'U' or 'u'   Only the upper triangular part of the */
/* >                                  symmetric matrix is to be referenced. */
/* > */
/* >              UPLO = 'L' or 'l'   Only the lower triangular part of the */
/* >                                  symmetric matrix is to be referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >           On entry,  M  specifies the number of rows of the matrix  C. */
/* >           M  must be at least zero. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >           On entry, N specifies the number of columns of the matrix C. */
/* >           N  must be at least zero. */
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
/* >           m  when  SIDE = 'L' or 'l'  and is  n otherwise. */
/* >           Before entry  with  SIDE = 'L' or 'l',  the  m by m  part of */
/* >           the array  A  must contain the  symmetric matrix,  such that */
/* >           when  UPLO = 'U' or 'u', the leading m by m upper triangular */
/* >           part of the array  A  must contain the upper triangular part */
/* >           of the  symmetric matrix and the  strictly  lower triangular */
/* >           part of  A  is not referenced,  and when  UPLO = 'L' or 'l', */
/* >           the leading  m by m  lower triangular part  of the  array  A */
/* >           must  contain  the  lower triangular part  of the  symmetric */
/* >           matrix and the  strictly upper triangular part of  A  is not */
/* >           referenced. */
/* >           Before entry  with  SIDE = 'R' or 'r',  the  n by n  part of */
/* >           the array  A  must contain the  symmetric matrix,  such that */
/* >           when  UPLO = 'U' or 'u', the leading n by n upper triangular */
/* >           part of the array  A  must contain the upper triangular part */
/* >           of the  symmetric matrix and the  strictly  lower triangular */
/* >           part of  A  is not referenced,  and when  UPLO = 'L' or 'l', */
/* >           the leading  n by n  lower triangular part  of the  array  A */
/* >           must  contain  the  lower triangular part  of the  symmetric */
/* >           matrix and the  strictly upper triangular part of  A  is not */
/* >           referenced. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >           On entry, LDA specifies the first dimension of A as declared */
/* >           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then */
/* >           LDA must be at least  max( 1, m ), otherwise  LDA must be at */
/* >           least  max( 1, n ). */
/* > \endverbatim */
/* > */
/* > \param[in] B */
/* > \verbatim */
/* >          B is DOUBLE PRECISION array of DIMENSION ( LDB, n ). */
/* >           Before entry, the leading  m by n part of the array  B  must */
/* >           contain the matrix B. */
/* > \endverbatim */
/* > */
/* > \param[in] LDB */
/* > \verbatim */
/* >          LDB is INTEGER */
/* >           On entry, LDB specifies the first dimension of B as declared */
/* >           in  the  calling  (sub)  program.   LDB  must  be  at  least */
/* >           max( 1, m ). */
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
/* >           On exit, the array  C  is overwritten by the  m by n updated */
/* >           matrix. */
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

/* > \date November 2011 */

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
/* Subroutine */ int dsymm_(char *side, char *uplo, integer *m, integer *n, 
	doublereal *alpha, doublereal *a, integer *lda, doublereal *b, 
	integer *ldb, doublereal *beta, doublereal *c__, integer *ldc, ftnlen 
	side_len, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2, 
	    i__3;

    /* Local variables */
    static integer i__, j, k, info;
    static doublereal temp1, temp2;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer nrowa;
    static logical upper;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);


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

/*     Set NROWA as the number of rows of A. */

#line 230 "dsymm.f"
    /* Parameter adjustments */
#line 230 "dsymm.f"
    a_dim1 = *lda;
#line 230 "dsymm.f"
    a_offset = 1 + a_dim1;
#line 230 "dsymm.f"
    a -= a_offset;
#line 230 "dsymm.f"
    b_dim1 = *ldb;
#line 230 "dsymm.f"
    b_offset = 1 + b_dim1;
#line 230 "dsymm.f"
    b -= b_offset;
#line 230 "dsymm.f"
    c_dim1 = *ldc;
#line 230 "dsymm.f"
    c_offset = 1 + c_dim1;
#line 230 "dsymm.f"
    c__ -= c_offset;
#line 230 "dsymm.f"

#line 230 "dsymm.f"
    /* Function Body */
#line 230 "dsymm.f"
    if (lsame_(side, "L", (ftnlen)1, (ftnlen)1)) {
#line 231 "dsymm.f"
	nrowa = *m;
#line 232 "dsymm.f"
    } else {
#line 233 "dsymm.f"
	nrowa = *n;
#line 234 "dsymm.f"
    }
#line 235 "dsymm.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);

/*     Test the input parameters. */

#line 239 "dsymm.f"
    info = 0;
#line 240 "dsymm.f"
    if (! lsame_(side, "L", (ftnlen)1, (ftnlen)1) && ! lsame_(side, "R", (
	    ftnlen)1, (ftnlen)1)) {
#line 241 "dsymm.f"
	info = 1;
#line 242 "dsymm.f"
    } else if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 243 "dsymm.f"
	info = 2;
#line 244 "dsymm.f"
    } else if (*m < 0) {
#line 245 "dsymm.f"
	info = 3;
#line 246 "dsymm.f"
    } else if (*n < 0) {
#line 247 "dsymm.f"
	info = 4;
#line 248 "dsymm.f"
    } else if (*lda < max(1,nrowa)) {
#line 249 "dsymm.f"
	info = 7;
#line 250 "dsymm.f"
    } else if (*ldb < max(1,*m)) {
#line 251 "dsymm.f"
	info = 9;
#line 252 "dsymm.f"
    } else if (*ldc < max(1,*m)) {
#line 253 "dsymm.f"
	info = 12;
#line 254 "dsymm.f"
    }
#line 255 "dsymm.f"
    if (info != 0) {
#line 256 "dsymm.f"
	xerbla_("DSYMM ", &info, (ftnlen)6);
#line 257 "dsymm.f"
	return 0;
#line 258 "dsymm.f"
    }

/*     Quick return if possible. */

#line 262 "dsymm.f"
    if (*m == 0 || *n == 0 || *alpha == 0. && *beta == 1.) {
#line 262 "dsymm.f"
	return 0;
#line 262 "dsymm.f"
    }

/*     And when  alpha.eq.zero. */

#line 267 "dsymm.f"
    if (*alpha == 0.) {
#line 268 "dsymm.f"
	if (*beta == 0.) {
#line 269 "dsymm.f"
	    i__1 = *n;
#line 269 "dsymm.f"
	    for (j = 1; j <= i__1; ++j) {
#line 270 "dsymm.f"
		i__2 = *m;
#line 270 "dsymm.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 271 "dsymm.f"
		    c__[i__ + j * c_dim1] = 0.;
#line 272 "dsymm.f"
/* L10: */
#line 272 "dsymm.f"
		}
#line 273 "dsymm.f"
/* L20: */
#line 273 "dsymm.f"
	    }
#line 274 "dsymm.f"
	} else {
#line 275 "dsymm.f"
	    i__1 = *n;
#line 275 "dsymm.f"
	    for (j = 1; j <= i__1; ++j) {
#line 276 "dsymm.f"
		i__2 = *m;
#line 276 "dsymm.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 277 "dsymm.f"
		    c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
#line 278 "dsymm.f"
/* L30: */
#line 278 "dsymm.f"
		}
#line 279 "dsymm.f"
/* L40: */
#line 279 "dsymm.f"
	    }
#line 280 "dsymm.f"
	}
#line 281 "dsymm.f"
	return 0;
#line 282 "dsymm.f"
    }

/*     Start the operations. */

#line 286 "dsymm.f"
    if (lsame_(side, "L", (ftnlen)1, (ftnlen)1)) {

/*        Form  C := alpha*A*B + beta*C. */

#line 290 "dsymm.f"
	if (upper) {
#line 291 "dsymm.f"
	    i__1 = *n;
#line 291 "dsymm.f"
	    for (j = 1; j <= i__1; ++j) {
#line 292 "dsymm.f"
		i__2 = *m;
#line 292 "dsymm.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 293 "dsymm.f"
		    temp1 = *alpha * b[i__ + j * b_dim1];
#line 294 "dsymm.f"
		    temp2 = 0.;
#line 295 "dsymm.f"
		    i__3 = i__ - 1;
#line 295 "dsymm.f"
		    for (k = 1; k <= i__3; ++k) {
#line 296 "dsymm.f"
			c__[k + j * c_dim1] += temp1 * a[k + i__ * a_dim1];
#line 297 "dsymm.f"
			temp2 += b[k + j * b_dim1] * a[k + i__ * a_dim1];
#line 298 "dsymm.f"
/* L50: */
#line 298 "dsymm.f"
		    }
#line 299 "dsymm.f"
		    if (*beta == 0.) {
#line 300 "dsymm.f"
			c__[i__ + j * c_dim1] = temp1 * a[i__ + i__ * a_dim1] 
				+ *alpha * temp2;
#line 301 "dsymm.f"
		    } else {
#line 302 "dsymm.f"
			c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1] 
				+ temp1 * a[i__ + i__ * a_dim1] + *alpha * 
				temp2;
#line 304 "dsymm.f"
		    }
#line 305 "dsymm.f"
/* L60: */
#line 305 "dsymm.f"
		}
#line 306 "dsymm.f"
/* L70: */
#line 306 "dsymm.f"
	    }
#line 307 "dsymm.f"
	} else {
#line 308 "dsymm.f"
	    i__1 = *n;
#line 308 "dsymm.f"
	    for (j = 1; j <= i__1; ++j) {
#line 309 "dsymm.f"
		for (i__ = *m; i__ >= 1; --i__) {
#line 310 "dsymm.f"
		    temp1 = *alpha * b[i__ + j * b_dim1];
#line 311 "dsymm.f"
		    temp2 = 0.;
#line 312 "dsymm.f"
		    i__2 = *m;
#line 312 "dsymm.f"
		    for (k = i__ + 1; k <= i__2; ++k) {
#line 313 "dsymm.f"
			c__[k + j * c_dim1] += temp1 * a[k + i__ * a_dim1];
#line 314 "dsymm.f"
			temp2 += b[k + j * b_dim1] * a[k + i__ * a_dim1];
#line 315 "dsymm.f"
/* L80: */
#line 315 "dsymm.f"
		    }
#line 316 "dsymm.f"
		    if (*beta == 0.) {
#line 317 "dsymm.f"
			c__[i__ + j * c_dim1] = temp1 * a[i__ + i__ * a_dim1] 
				+ *alpha * temp2;
#line 318 "dsymm.f"
		    } else {
#line 319 "dsymm.f"
			c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1] 
				+ temp1 * a[i__ + i__ * a_dim1] + *alpha * 
				temp2;
#line 321 "dsymm.f"
		    }
#line 322 "dsymm.f"
/* L90: */
#line 322 "dsymm.f"
		}
#line 323 "dsymm.f"
/* L100: */
#line 323 "dsymm.f"
	    }
#line 324 "dsymm.f"
	}
#line 325 "dsymm.f"
    } else {

/*        Form  C := alpha*B*A + beta*C. */

#line 329 "dsymm.f"
	i__1 = *n;
#line 329 "dsymm.f"
	for (j = 1; j <= i__1; ++j) {
#line 330 "dsymm.f"
	    temp1 = *alpha * a[j + j * a_dim1];
#line 331 "dsymm.f"
	    if (*beta == 0.) {
#line 332 "dsymm.f"
		i__2 = *m;
#line 332 "dsymm.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 333 "dsymm.f"
		    c__[i__ + j * c_dim1] = temp1 * b[i__ + j * b_dim1];
#line 334 "dsymm.f"
/* L110: */
#line 334 "dsymm.f"
		}
#line 335 "dsymm.f"
	    } else {
#line 336 "dsymm.f"
		i__2 = *m;
#line 336 "dsymm.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 337 "dsymm.f"
		    c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1] + 
			    temp1 * b[i__ + j * b_dim1];
#line 338 "dsymm.f"
/* L120: */
#line 338 "dsymm.f"
		}
#line 339 "dsymm.f"
	    }
#line 340 "dsymm.f"
	    i__2 = j - 1;
#line 340 "dsymm.f"
	    for (k = 1; k <= i__2; ++k) {
#line 341 "dsymm.f"
		if (upper) {
#line 342 "dsymm.f"
		    temp1 = *alpha * a[k + j * a_dim1];
#line 343 "dsymm.f"
		} else {
#line 344 "dsymm.f"
		    temp1 = *alpha * a[j + k * a_dim1];
#line 345 "dsymm.f"
		}
#line 346 "dsymm.f"
		i__3 = *m;
#line 346 "dsymm.f"
		for (i__ = 1; i__ <= i__3; ++i__) {
#line 347 "dsymm.f"
		    c__[i__ + j * c_dim1] += temp1 * b[i__ + k * b_dim1];
#line 348 "dsymm.f"
/* L130: */
#line 348 "dsymm.f"
		}
#line 349 "dsymm.f"
/* L140: */
#line 349 "dsymm.f"
	    }
#line 350 "dsymm.f"
	    i__2 = *n;
#line 350 "dsymm.f"
	    for (k = j + 1; k <= i__2; ++k) {
#line 351 "dsymm.f"
		if (upper) {
#line 352 "dsymm.f"
		    temp1 = *alpha * a[j + k * a_dim1];
#line 353 "dsymm.f"
		} else {
#line 354 "dsymm.f"
		    temp1 = *alpha * a[k + j * a_dim1];
#line 355 "dsymm.f"
		}
#line 356 "dsymm.f"
		i__3 = *m;
#line 356 "dsymm.f"
		for (i__ = 1; i__ <= i__3; ++i__) {
#line 357 "dsymm.f"
		    c__[i__ + j * c_dim1] += temp1 * b[i__ + k * b_dim1];
#line 358 "dsymm.f"
/* L150: */
#line 358 "dsymm.f"
		}
#line 359 "dsymm.f"
/* L160: */
#line 359 "dsymm.f"
	    }
#line 360 "dsymm.f"
/* L170: */
#line 360 "dsymm.f"
	}
#line 361 "dsymm.f"
    }

#line 363 "dsymm.f"
    return 0;

/*     End of DSYMM . */

} /* dsymm_ */

