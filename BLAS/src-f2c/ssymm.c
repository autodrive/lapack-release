#line 1 "ssymm.f"
/* ssymm.f -- translated by f2c (version 20100827).
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

#line 1 "ssymm.f"
/* > \brief \b SSYMM */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SSYMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC) */

/*       .. Scalar Arguments .. */
/*       REAL ALPHA,BETA */
/*       INTEGER LDA,LDB,LDC,M,N */
/*       CHARACTER SIDE,UPLO */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL A(LDA,*),B(LDB,*),C(LDC,*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SSYMM  performs one of the matrix-matrix operations */
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
/* >          ALPHA is REAL */
/* >           On entry, ALPHA specifies the scalar alpha. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is REAL array of DIMENSION ( LDA, ka ), where ka is */
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
/* >          B is REAL array of DIMENSION ( LDB, n ). */
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
/* Subroutine */ int ssymm_(char *side, char *uplo, integer *m, integer *n, 
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

/*     Set NROWA as the number of rows of A. */

#line 230 "ssymm.f"
    /* Parameter adjustments */
#line 230 "ssymm.f"
    a_dim1 = *lda;
#line 230 "ssymm.f"
    a_offset = 1 + a_dim1;
#line 230 "ssymm.f"
    a -= a_offset;
#line 230 "ssymm.f"
    b_dim1 = *ldb;
#line 230 "ssymm.f"
    b_offset = 1 + b_dim1;
#line 230 "ssymm.f"
    b -= b_offset;
#line 230 "ssymm.f"
    c_dim1 = *ldc;
#line 230 "ssymm.f"
    c_offset = 1 + c_dim1;
#line 230 "ssymm.f"
    c__ -= c_offset;
#line 230 "ssymm.f"

#line 230 "ssymm.f"
    /* Function Body */
#line 230 "ssymm.f"
    if (lsame_(side, "L", (ftnlen)1, (ftnlen)1)) {
#line 231 "ssymm.f"
	nrowa = *m;
#line 232 "ssymm.f"
    } else {
#line 233 "ssymm.f"
	nrowa = *n;
#line 234 "ssymm.f"
    }
#line 235 "ssymm.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);

/*     Test the input parameters. */

#line 239 "ssymm.f"
    info = 0;
#line 240 "ssymm.f"
    if (! lsame_(side, "L", (ftnlen)1, (ftnlen)1) && ! lsame_(side, "R", (
	    ftnlen)1, (ftnlen)1)) {
#line 241 "ssymm.f"
	info = 1;
#line 242 "ssymm.f"
    } else if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 243 "ssymm.f"
	info = 2;
#line 244 "ssymm.f"
    } else if (*m < 0) {
#line 245 "ssymm.f"
	info = 3;
#line 246 "ssymm.f"
    } else if (*n < 0) {
#line 247 "ssymm.f"
	info = 4;
#line 248 "ssymm.f"
    } else if (*lda < max(1,nrowa)) {
#line 249 "ssymm.f"
	info = 7;
#line 250 "ssymm.f"
    } else if (*ldb < max(1,*m)) {
#line 251 "ssymm.f"
	info = 9;
#line 252 "ssymm.f"
    } else if (*ldc < max(1,*m)) {
#line 253 "ssymm.f"
	info = 12;
#line 254 "ssymm.f"
    }
#line 255 "ssymm.f"
    if (info != 0) {
#line 256 "ssymm.f"
	xerbla_("SSYMM ", &info, (ftnlen)6);
#line 257 "ssymm.f"
	return 0;
#line 258 "ssymm.f"
    }

/*     Quick return if possible. */

#line 262 "ssymm.f"
    if (*m == 0 || *n == 0 || *alpha == 0. && *beta == 1.) {
#line 262 "ssymm.f"
	return 0;
#line 262 "ssymm.f"
    }

/*     And when  alpha.eq.zero. */

#line 267 "ssymm.f"
    if (*alpha == 0.) {
#line 268 "ssymm.f"
	if (*beta == 0.) {
#line 269 "ssymm.f"
	    i__1 = *n;
#line 269 "ssymm.f"
	    for (j = 1; j <= i__1; ++j) {
#line 270 "ssymm.f"
		i__2 = *m;
#line 270 "ssymm.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 271 "ssymm.f"
		    c__[i__ + j * c_dim1] = 0.;
#line 272 "ssymm.f"
/* L10: */
#line 272 "ssymm.f"
		}
#line 273 "ssymm.f"
/* L20: */
#line 273 "ssymm.f"
	    }
#line 274 "ssymm.f"
	} else {
#line 275 "ssymm.f"
	    i__1 = *n;
#line 275 "ssymm.f"
	    for (j = 1; j <= i__1; ++j) {
#line 276 "ssymm.f"
		i__2 = *m;
#line 276 "ssymm.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 277 "ssymm.f"
		    c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1];
#line 278 "ssymm.f"
/* L30: */
#line 278 "ssymm.f"
		}
#line 279 "ssymm.f"
/* L40: */
#line 279 "ssymm.f"
	    }
#line 280 "ssymm.f"
	}
#line 281 "ssymm.f"
	return 0;
#line 282 "ssymm.f"
    }

/*     Start the operations. */

#line 286 "ssymm.f"
    if (lsame_(side, "L", (ftnlen)1, (ftnlen)1)) {

/*        Form  C := alpha*A*B + beta*C. */

#line 290 "ssymm.f"
	if (upper) {
#line 291 "ssymm.f"
	    i__1 = *n;
#line 291 "ssymm.f"
	    for (j = 1; j <= i__1; ++j) {
#line 292 "ssymm.f"
		i__2 = *m;
#line 292 "ssymm.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 293 "ssymm.f"
		    temp1 = *alpha * b[i__ + j * b_dim1];
#line 294 "ssymm.f"
		    temp2 = 0.;
#line 295 "ssymm.f"
		    i__3 = i__ - 1;
#line 295 "ssymm.f"
		    for (k = 1; k <= i__3; ++k) {
#line 296 "ssymm.f"
			c__[k + j * c_dim1] += temp1 * a[k + i__ * a_dim1];
#line 297 "ssymm.f"
			temp2 += b[k + j * b_dim1] * a[k + i__ * a_dim1];
#line 298 "ssymm.f"
/* L50: */
#line 298 "ssymm.f"
		    }
#line 299 "ssymm.f"
		    if (*beta == 0.) {
#line 300 "ssymm.f"
			c__[i__ + j * c_dim1] = temp1 * a[i__ + i__ * a_dim1] 
				+ *alpha * temp2;
#line 301 "ssymm.f"
		    } else {
#line 302 "ssymm.f"
			c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1] 
				+ temp1 * a[i__ + i__ * a_dim1] + *alpha * 
				temp2;
#line 304 "ssymm.f"
		    }
#line 305 "ssymm.f"
/* L60: */
#line 305 "ssymm.f"
		}
#line 306 "ssymm.f"
/* L70: */
#line 306 "ssymm.f"
	    }
#line 307 "ssymm.f"
	} else {
#line 308 "ssymm.f"
	    i__1 = *n;
#line 308 "ssymm.f"
	    for (j = 1; j <= i__1; ++j) {
#line 309 "ssymm.f"
		for (i__ = *m; i__ >= 1; --i__) {
#line 310 "ssymm.f"
		    temp1 = *alpha * b[i__ + j * b_dim1];
#line 311 "ssymm.f"
		    temp2 = 0.;
#line 312 "ssymm.f"
		    i__2 = *m;
#line 312 "ssymm.f"
		    for (k = i__ + 1; k <= i__2; ++k) {
#line 313 "ssymm.f"
			c__[k + j * c_dim1] += temp1 * a[k + i__ * a_dim1];
#line 314 "ssymm.f"
			temp2 += b[k + j * b_dim1] * a[k + i__ * a_dim1];
#line 315 "ssymm.f"
/* L80: */
#line 315 "ssymm.f"
		    }
#line 316 "ssymm.f"
		    if (*beta == 0.) {
#line 317 "ssymm.f"
			c__[i__ + j * c_dim1] = temp1 * a[i__ + i__ * a_dim1] 
				+ *alpha * temp2;
#line 318 "ssymm.f"
		    } else {
#line 319 "ssymm.f"
			c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1] 
				+ temp1 * a[i__ + i__ * a_dim1] + *alpha * 
				temp2;
#line 321 "ssymm.f"
		    }
#line 322 "ssymm.f"
/* L90: */
#line 322 "ssymm.f"
		}
#line 323 "ssymm.f"
/* L100: */
#line 323 "ssymm.f"
	    }
#line 324 "ssymm.f"
	}
#line 325 "ssymm.f"
    } else {

/*        Form  C := alpha*B*A + beta*C. */

#line 329 "ssymm.f"
	i__1 = *n;
#line 329 "ssymm.f"
	for (j = 1; j <= i__1; ++j) {
#line 330 "ssymm.f"
	    temp1 = *alpha * a[j + j * a_dim1];
#line 331 "ssymm.f"
	    if (*beta == 0.) {
#line 332 "ssymm.f"
		i__2 = *m;
#line 332 "ssymm.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 333 "ssymm.f"
		    c__[i__ + j * c_dim1] = temp1 * b[i__ + j * b_dim1];
#line 334 "ssymm.f"
/* L110: */
#line 334 "ssymm.f"
		}
#line 335 "ssymm.f"
	    } else {
#line 336 "ssymm.f"
		i__2 = *m;
#line 336 "ssymm.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 337 "ssymm.f"
		    c__[i__ + j * c_dim1] = *beta * c__[i__ + j * c_dim1] + 
			    temp1 * b[i__ + j * b_dim1];
#line 338 "ssymm.f"
/* L120: */
#line 338 "ssymm.f"
		}
#line 339 "ssymm.f"
	    }
#line 340 "ssymm.f"
	    i__2 = j - 1;
#line 340 "ssymm.f"
	    for (k = 1; k <= i__2; ++k) {
#line 341 "ssymm.f"
		if (upper) {
#line 342 "ssymm.f"
		    temp1 = *alpha * a[k + j * a_dim1];
#line 343 "ssymm.f"
		} else {
#line 344 "ssymm.f"
		    temp1 = *alpha * a[j + k * a_dim1];
#line 345 "ssymm.f"
		}
#line 346 "ssymm.f"
		i__3 = *m;
#line 346 "ssymm.f"
		for (i__ = 1; i__ <= i__3; ++i__) {
#line 347 "ssymm.f"
		    c__[i__ + j * c_dim1] += temp1 * b[i__ + k * b_dim1];
#line 348 "ssymm.f"
/* L130: */
#line 348 "ssymm.f"
		}
#line 349 "ssymm.f"
/* L140: */
#line 349 "ssymm.f"
	    }
#line 350 "ssymm.f"
	    i__2 = *n;
#line 350 "ssymm.f"
	    for (k = j + 1; k <= i__2; ++k) {
#line 351 "ssymm.f"
		if (upper) {
#line 352 "ssymm.f"
		    temp1 = *alpha * a[j + k * a_dim1];
#line 353 "ssymm.f"
		} else {
#line 354 "ssymm.f"
		    temp1 = *alpha * a[k + j * a_dim1];
#line 355 "ssymm.f"
		}
#line 356 "ssymm.f"
		i__3 = *m;
#line 356 "ssymm.f"
		for (i__ = 1; i__ <= i__3; ++i__) {
#line 357 "ssymm.f"
		    c__[i__ + j * c_dim1] += temp1 * b[i__ + k * b_dim1];
#line 358 "ssymm.f"
/* L150: */
#line 358 "ssymm.f"
		}
#line 359 "ssymm.f"
/* L160: */
#line 359 "ssymm.f"
	    }
#line 360 "ssymm.f"
/* L170: */
#line 360 "ssymm.f"
	}
#line 361 "ssymm.f"
    }

#line 363 "ssymm.f"
    return 0;

/*     End of SSYMM . */

} /* ssymm_ */

