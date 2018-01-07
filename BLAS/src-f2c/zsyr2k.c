#line 1 "zsyr2k.f"
/* zsyr2k.f -- translated by f2c (version 20100827).
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

#line 1 "zsyr2k.f"
/* > \brief \b ZSYR2K */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZSYR2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC) */

/*       .. Scalar Arguments .. */
/*       COMPLEX*16 ALPHA,BETA */
/*       INTEGER K,LDA,LDB,LDC,N */
/*       CHARACTER TRANS,UPLO */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16 A(LDA,*),B(LDB,*),C(LDC,*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZSYR2K  performs one of the symmetric rank 2k operations */
/* > */
/* >    C := alpha*A*B**T + alpha*B*A**T + beta*C, */
/* > */
/* > or */
/* > */
/* >    C := alpha*A**T*B + alpha*B**T*A + beta*C, */
/* > */
/* > where  alpha and beta  are scalars,  C is an  n by n symmetric matrix */
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
/* >              TRANS = 'N' or 'n'    C := alpha*A*B**T + alpha*B*A**T + */
/* >                                         beta*C. */
/* > */
/* >              TRANS = 'T' or 't'    C := alpha*A**T*B + alpha*B**T*A + */
/* >                                         beta*C. */
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
/* >           TRANS = 'T' or 't',  K  specifies  the number of rows of the */
/* >           matrices  A and B.  K must be at least zero. */
/* > \endverbatim */
/* > */
/* > \param[in] ALPHA */
/* > \verbatim */
/* >          ALPHA is COMPLEX*16 */
/* >           On entry, ALPHA specifies the scalar alpha. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is COMPLEX*16 array of DIMENSION ( LDA, ka ), where ka is */
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
/* >          B is COMPLEX*16 array of DIMENSION ( LDB, kb ), where kb is */
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
/* >          BETA is COMPLEX*16 */
/* >           On entry, BETA specifies the scalar beta. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* >          C is COMPLEX*16 array of DIMENSION ( LDC, n ). */
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

/* > \date November 2011 */

/* > \ingroup complex16_blas_level3 */

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
/* Subroutine */ int zsyr2k_(char *uplo, char *trans, integer *n, integer *k, 
	doublecomplex *alpha, doublecomplex *a, integer *lda, doublecomplex *
	b, integer *ldb, doublecomplex *beta, doublecomplex *c__, integer *
	ldc, ftnlen uplo_len, ftnlen trans_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2, 
	    i__3, i__4, i__5, i__6, i__7;
    doublecomplex z__1, z__2, z__3, z__4, z__5;

    /* Local variables */
    static integer i__, j, l, info;
    static doublecomplex temp1, temp2;
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

/*     Test the input parameters. */

#line 231 "zsyr2k.f"
    /* Parameter adjustments */
#line 231 "zsyr2k.f"
    a_dim1 = *lda;
#line 231 "zsyr2k.f"
    a_offset = 1 + a_dim1;
#line 231 "zsyr2k.f"
    a -= a_offset;
#line 231 "zsyr2k.f"
    b_dim1 = *ldb;
#line 231 "zsyr2k.f"
    b_offset = 1 + b_dim1;
#line 231 "zsyr2k.f"
    b -= b_offset;
#line 231 "zsyr2k.f"
    c_dim1 = *ldc;
#line 231 "zsyr2k.f"
    c_offset = 1 + c_dim1;
#line 231 "zsyr2k.f"
    c__ -= c_offset;
#line 231 "zsyr2k.f"

#line 231 "zsyr2k.f"
    /* Function Body */
#line 231 "zsyr2k.f"
    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {
#line 232 "zsyr2k.f"
	nrowa = *n;
#line 233 "zsyr2k.f"
    } else {
#line 234 "zsyr2k.f"
	nrowa = *k;
#line 235 "zsyr2k.f"
    }
#line 236 "zsyr2k.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);

#line 238 "zsyr2k.f"
    info = 0;
#line 239 "zsyr2k.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 240 "zsyr2k.f"
	info = 1;
#line 241 "zsyr2k.f"
    } else if (! lsame_(trans, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, 
	    "T", (ftnlen)1, (ftnlen)1)) {
#line 243 "zsyr2k.f"
	info = 2;
#line 244 "zsyr2k.f"
    } else if (*n < 0) {
#line 245 "zsyr2k.f"
	info = 3;
#line 246 "zsyr2k.f"
    } else if (*k < 0) {
#line 247 "zsyr2k.f"
	info = 4;
#line 248 "zsyr2k.f"
    } else if (*lda < max(1,nrowa)) {
#line 249 "zsyr2k.f"
	info = 7;
#line 250 "zsyr2k.f"
    } else if (*ldb < max(1,nrowa)) {
#line 251 "zsyr2k.f"
	info = 9;
#line 252 "zsyr2k.f"
    } else if (*ldc < max(1,*n)) {
#line 253 "zsyr2k.f"
	info = 12;
#line 254 "zsyr2k.f"
    }
#line 255 "zsyr2k.f"
    if (info != 0) {
#line 256 "zsyr2k.f"
	xerbla_("ZSYR2K", &info, (ftnlen)6);
#line 257 "zsyr2k.f"
	return 0;
#line 258 "zsyr2k.f"
    }

/*     Quick return if possible. */

#line 262 "zsyr2k.f"
    if (*n == 0 || (alpha->r == 0. && alpha->i == 0. || *k == 0) && (beta->r 
	    == 1. && beta->i == 0.)) {
#line 262 "zsyr2k.f"
	return 0;
#line 262 "zsyr2k.f"
    }

/*     And when  alpha.eq.zero. */

#line 267 "zsyr2k.f"
    if (alpha->r == 0. && alpha->i == 0.) {
#line 268 "zsyr2k.f"
	if (upper) {
#line 269 "zsyr2k.f"
	    if (beta->r == 0. && beta->i == 0.) {
#line 270 "zsyr2k.f"
		i__1 = *n;
#line 270 "zsyr2k.f"
		for (j = 1; j <= i__1; ++j) {
#line 271 "zsyr2k.f"
		    i__2 = j;
#line 271 "zsyr2k.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 272 "zsyr2k.f"
			i__3 = i__ + j * c_dim1;
#line 272 "zsyr2k.f"
			c__[i__3].r = 0., c__[i__3].i = 0.;
#line 273 "zsyr2k.f"
/* L10: */
#line 273 "zsyr2k.f"
		    }
#line 274 "zsyr2k.f"
/* L20: */
#line 274 "zsyr2k.f"
		}
#line 275 "zsyr2k.f"
	    } else {
#line 276 "zsyr2k.f"
		i__1 = *n;
#line 276 "zsyr2k.f"
		for (j = 1; j <= i__1; ++j) {
#line 277 "zsyr2k.f"
		    i__2 = j;
#line 277 "zsyr2k.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 278 "zsyr2k.f"
			i__3 = i__ + j * c_dim1;
#line 278 "zsyr2k.f"
			i__4 = i__ + j * c_dim1;
#line 278 "zsyr2k.f"
			z__1.r = beta->r * c__[i__4].r - beta->i * c__[i__4]
				.i, z__1.i = beta->r * c__[i__4].i + beta->i *
				 c__[i__4].r;
#line 278 "zsyr2k.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 279 "zsyr2k.f"
/* L30: */
#line 279 "zsyr2k.f"
		    }
#line 280 "zsyr2k.f"
/* L40: */
#line 280 "zsyr2k.f"
		}
#line 281 "zsyr2k.f"
	    }
#line 282 "zsyr2k.f"
	} else {
#line 283 "zsyr2k.f"
	    if (beta->r == 0. && beta->i == 0.) {
#line 284 "zsyr2k.f"
		i__1 = *n;
#line 284 "zsyr2k.f"
		for (j = 1; j <= i__1; ++j) {
#line 285 "zsyr2k.f"
		    i__2 = *n;
#line 285 "zsyr2k.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 286 "zsyr2k.f"
			i__3 = i__ + j * c_dim1;
#line 286 "zsyr2k.f"
			c__[i__3].r = 0., c__[i__3].i = 0.;
#line 287 "zsyr2k.f"
/* L50: */
#line 287 "zsyr2k.f"
		    }
#line 288 "zsyr2k.f"
/* L60: */
#line 288 "zsyr2k.f"
		}
#line 289 "zsyr2k.f"
	    } else {
#line 290 "zsyr2k.f"
		i__1 = *n;
#line 290 "zsyr2k.f"
		for (j = 1; j <= i__1; ++j) {
#line 291 "zsyr2k.f"
		    i__2 = *n;
#line 291 "zsyr2k.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 292 "zsyr2k.f"
			i__3 = i__ + j * c_dim1;
#line 292 "zsyr2k.f"
			i__4 = i__ + j * c_dim1;
#line 292 "zsyr2k.f"
			z__1.r = beta->r * c__[i__4].r - beta->i * c__[i__4]
				.i, z__1.i = beta->r * c__[i__4].i + beta->i *
				 c__[i__4].r;
#line 292 "zsyr2k.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 293 "zsyr2k.f"
/* L70: */
#line 293 "zsyr2k.f"
		    }
#line 294 "zsyr2k.f"
/* L80: */
#line 294 "zsyr2k.f"
		}
#line 295 "zsyr2k.f"
	    }
#line 296 "zsyr2k.f"
	}
#line 297 "zsyr2k.f"
	return 0;
#line 298 "zsyr2k.f"
    }

/*     Start the operations. */

#line 302 "zsyr2k.f"
    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {

/*        Form  C := alpha*A*B**T + alpha*B*A**T + C. */

#line 306 "zsyr2k.f"
	if (upper) {
#line 307 "zsyr2k.f"
	    i__1 = *n;
#line 307 "zsyr2k.f"
	    for (j = 1; j <= i__1; ++j) {
#line 308 "zsyr2k.f"
		if (beta->r == 0. && beta->i == 0.) {
#line 309 "zsyr2k.f"
		    i__2 = j;
#line 309 "zsyr2k.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 310 "zsyr2k.f"
			i__3 = i__ + j * c_dim1;
#line 310 "zsyr2k.f"
			c__[i__3].r = 0., c__[i__3].i = 0.;
#line 311 "zsyr2k.f"
/* L90: */
#line 311 "zsyr2k.f"
		    }
#line 312 "zsyr2k.f"
		} else if (beta->r != 1. || beta->i != 0.) {
#line 313 "zsyr2k.f"
		    i__2 = j;
#line 313 "zsyr2k.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 314 "zsyr2k.f"
			i__3 = i__ + j * c_dim1;
#line 314 "zsyr2k.f"
			i__4 = i__ + j * c_dim1;
#line 314 "zsyr2k.f"
			z__1.r = beta->r * c__[i__4].r - beta->i * c__[i__4]
				.i, z__1.i = beta->r * c__[i__4].i + beta->i *
				 c__[i__4].r;
#line 314 "zsyr2k.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 315 "zsyr2k.f"
/* L100: */
#line 315 "zsyr2k.f"
		    }
#line 316 "zsyr2k.f"
		}
#line 317 "zsyr2k.f"
		i__2 = *k;
#line 317 "zsyr2k.f"
		for (l = 1; l <= i__2; ++l) {
#line 318 "zsyr2k.f"
		    i__3 = j + l * a_dim1;
#line 318 "zsyr2k.f"
		    i__4 = j + l * b_dim1;
#line 318 "zsyr2k.f"
		    if (a[i__3].r != 0. || a[i__3].i != 0. || (b[i__4].r != 
			    0. || b[i__4].i != 0.)) {
#line 319 "zsyr2k.f"
			i__3 = j + l * b_dim1;
#line 319 "zsyr2k.f"
			z__1.r = alpha->r * b[i__3].r - alpha->i * b[i__3].i, 
				z__1.i = alpha->r * b[i__3].i + alpha->i * b[
				i__3].r;
#line 319 "zsyr2k.f"
			temp1.r = z__1.r, temp1.i = z__1.i;
#line 320 "zsyr2k.f"
			i__3 = j + l * a_dim1;
#line 320 "zsyr2k.f"
			z__1.r = alpha->r * a[i__3].r - alpha->i * a[i__3].i, 
				z__1.i = alpha->r * a[i__3].i + alpha->i * a[
				i__3].r;
#line 320 "zsyr2k.f"
			temp2.r = z__1.r, temp2.i = z__1.i;
#line 321 "zsyr2k.f"
			i__3 = j;
#line 321 "zsyr2k.f"
			for (i__ = 1; i__ <= i__3; ++i__) {
#line 322 "zsyr2k.f"
			    i__4 = i__ + j * c_dim1;
#line 322 "zsyr2k.f"
			    i__5 = i__ + j * c_dim1;
#line 322 "zsyr2k.f"
			    i__6 = i__ + l * a_dim1;
#line 322 "zsyr2k.f"
			    z__3.r = a[i__6].r * temp1.r - a[i__6].i * 
				    temp1.i, z__3.i = a[i__6].r * temp1.i + a[
				    i__6].i * temp1.r;
#line 322 "zsyr2k.f"
			    z__2.r = c__[i__5].r + z__3.r, z__2.i = c__[i__5]
				    .i + z__3.i;
#line 322 "zsyr2k.f"
			    i__7 = i__ + l * b_dim1;
#line 322 "zsyr2k.f"
			    z__4.r = b[i__7].r * temp2.r - b[i__7].i * 
				    temp2.i, z__4.i = b[i__7].r * temp2.i + b[
				    i__7].i * temp2.r;
#line 322 "zsyr2k.f"
			    z__1.r = z__2.r + z__4.r, z__1.i = z__2.i + 
				    z__4.i;
#line 322 "zsyr2k.f"
			    c__[i__4].r = z__1.r, c__[i__4].i = z__1.i;
#line 324 "zsyr2k.f"
/* L110: */
#line 324 "zsyr2k.f"
			}
#line 325 "zsyr2k.f"
		    }
#line 326 "zsyr2k.f"
/* L120: */
#line 326 "zsyr2k.f"
		}
#line 327 "zsyr2k.f"
/* L130: */
#line 327 "zsyr2k.f"
	    }
#line 328 "zsyr2k.f"
	} else {
#line 329 "zsyr2k.f"
	    i__1 = *n;
#line 329 "zsyr2k.f"
	    for (j = 1; j <= i__1; ++j) {
#line 330 "zsyr2k.f"
		if (beta->r == 0. && beta->i == 0.) {
#line 331 "zsyr2k.f"
		    i__2 = *n;
#line 331 "zsyr2k.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 332 "zsyr2k.f"
			i__3 = i__ + j * c_dim1;
#line 332 "zsyr2k.f"
			c__[i__3].r = 0., c__[i__3].i = 0.;
#line 333 "zsyr2k.f"
/* L140: */
#line 333 "zsyr2k.f"
		    }
#line 334 "zsyr2k.f"
		} else if (beta->r != 1. || beta->i != 0.) {
#line 335 "zsyr2k.f"
		    i__2 = *n;
#line 335 "zsyr2k.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 336 "zsyr2k.f"
			i__3 = i__ + j * c_dim1;
#line 336 "zsyr2k.f"
			i__4 = i__ + j * c_dim1;
#line 336 "zsyr2k.f"
			z__1.r = beta->r * c__[i__4].r - beta->i * c__[i__4]
				.i, z__1.i = beta->r * c__[i__4].i + beta->i *
				 c__[i__4].r;
#line 336 "zsyr2k.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 337 "zsyr2k.f"
/* L150: */
#line 337 "zsyr2k.f"
		    }
#line 338 "zsyr2k.f"
		}
#line 339 "zsyr2k.f"
		i__2 = *k;
#line 339 "zsyr2k.f"
		for (l = 1; l <= i__2; ++l) {
#line 340 "zsyr2k.f"
		    i__3 = j + l * a_dim1;
#line 340 "zsyr2k.f"
		    i__4 = j + l * b_dim1;
#line 340 "zsyr2k.f"
		    if (a[i__3].r != 0. || a[i__3].i != 0. || (b[i__4].r != 
			    0. || b[i__4].i != 0.)) {
#line 341 "zsyr2k.f"
			i__3 = j + l * b_dim1;
#line 341 "zsyr2k.f"
			z__1.r = alpha->r * b[i__3].r - alpha->i * b[i__3].i, 
				z__1.i = alpha->r * b[i__3].i + alpha->i * b[
				i__3].r;
#line 341 "zsyr2k.f"
			temp1.r = z__1.r, temp1.i = z__1.i;
#line 342 "zsyr2k.f"
			i__3 = j + l * a_dim1;
#line 342 "zsyr2k.f"
			z__1.r = alpha->r * a[i__3].r - alpha->i * a[i__3].i, 
				z__1.i = alpha->r * a[i__3].i + alpha->i * a[
				i__3].r;
#line 342 "zsyr2k.f"
			temp2.r = z__1.r, temp2.i = z__1.i;
#line 343 "zsyr2k.f"
			i__3 = *n;
#line 343 "zsyr2k.f"
			for (i__ = j; i__ <= i__3; ++i__) {
#line 344 "zsyr2k.f"
			    i__4 = i__ + j * c_dim1;
#line 344 "zsyr2k.f"
			    i__5 = i__ + j * c_dim1;
#line 344 "zsyr2k.f"
			    i__6 = i__ + l * a_dim1;
#line 344 "zsyr2k.f"
			    z__3.r = a[i__6].r * temp1.r - a[i__6].i * 
				    temp1.i, z__3.i = a[i__6].r * temp1.i + a[
				    i__6].i * temp1.r;
#line 344 "zsyr2k.f"
			    z__2.r = c__[i__5].r + z__3.r, z__2.i = c__[i__5]
				    .i + z__3.i;
#line 344 "zsyr2k.f"
			    i__7 = i__ + l * b_dim1;
#line 344 "zsyr2k.f"
			    z__4.r = b[i__7].r * temp2.r - b[i__7].i * 
				    temp2.i, z__4.i = b[i__7].r * temp2.i + b[
				    i__7].i * temp2.r;
#line 344 "zsyr2k.f"
			    z__1.r = z__2.r + z__4.r, z__1.i = z__2.i + 
				    z__4.i;
#line 344 "zsyr2k.f"
			    c__[i__4].r = z__1.r, c__[i__4].i = z__1.i;
#line 346 "zsyr2k.f"
/* L160: */
#line 346 "zsyr2k.f"
			}
#line 347 "zsyr2k.f"
		    }
#line 348 "zsyr2k.f"
/* L170: */
#line 348 "zsyr2k.f"
		}
#line 349 "zsyr2k.f"
/* L180: */
#line 349 "zsyr2k.f"
	    }
#line 350 "zsyr2k.f"
	}
#line 351 "zsyr2k.f"
    } else {

/*        Form  C := alpha*A**T*B + alpha*B**T*A + C. */

#line 355 "zsyr2k.f"
	if (upper) {
#line 356 "zsyr2k.f"
	    i__1 = *n;
#line 356 "zsyr2k.f"
	    for (j = 1; j <= i__1; ++j) {
#line 357 "zsyr2k.f"
		i__2 = j;
#line 357 "zsyr2k.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 358 "zsyr2k.f"
		    temp1.r = 0., temp1.i = 0.;
#line 359 "zsyr2k.f"
		    temp2.r = 0., temp2.i = 0.;
#line 360 "zsyr2k.f"
		    i__3 = *k;
#line 360 "zsyr2k.f"
		    for (l = 1; l <= i__3; ++l) {
#line 361 "zsyr2k.f"
			i__4 = l + i__ * a_dim1;
#line 361 "zsyr2k.f"
			i__5 = l + j * b_dim1;
#line 361 "zsyr2k.f"
			z__2.r = a[i__4].r * b[i__5].r - a[i__4].i * b[i__5]
				.i, z__2.i = a[i__4].r * b[i__5].i + a[i__4]
				.i * b[i__5].r;
#line 361 "zsyr2k.f"
			z__1.r = temp1.r + z__2.r, z__1.i = temp1.i + z__2.i;
#line 361 "zsyr2k.f"
			temp1.r = z__1.r, temp1.i = z__1.i;
#line 362 "zsyr2k.f"
			i__4 = l + i__ * b_dim1;
#line 362 "zsyr2k.f"
			i__5 = l + j * a_dim1;
#line 362 "zsyr2k.f"
			z__2.r = b[i__4].r * a[i__5].r - b[i__4].i * a[i__5]
				.i, z__2.i = b[i__4].r * a[i__5].i + b[i__4]
				.i * a[i__5].r;
#line 362 "zsyr2k.f"
			z__1.r = temp2.r + z__2.r, z__1.i = temp2.i + z__2.i;
#line 362 "zsyr2k.f"
			temp2.r = z__1.r, temp2.i = z__1.i;
#line 363 "zsyr2k.f"
/* L190: */
#line 363 "zsyr2k.f"
		    }
#line 364 "zsyr2k.f"
		    if (beta->r == 0. && beta->i == 0.) {
#line 365 "zsyr2k.f"
			i__3 = i__ + j * c_dim1;
#line 365 "zsyr2k.f"
			z__2.r = alpha->r * temp1.r - alpha->i * temp1.i, 
				z__2.i = alpha->r * temp1.i + alpha->i * 
				temp1.r;
#line 365 "zsyr2k.f"
			z__3.r = alpha->r * temp2.r - alpha->i * temp2.i, 
				z__3.i = alpha->r * temp2.i + alpha->i * 
				temp2.r;
#line 365 "zsyr2k.f"
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 365 "zsyr2k.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 366 "zsyr2k.f"
		    } else {
#line 367 "zsyr2k.f"
			i__3 = i__ + j * c_dim1;
#line 367 "zsyr2k.f"
			i__4 = i__ + j * c_dim1;
#line 367 "zsyr2k.f"
			z__3.r = beta->r * c__[i__4].r - beta->i * c__[i__4]
				.i, z__3.i = beta->r * c__[i__4].i + beta->i *
				 c__[i__4].r;
#line 367 "zsyr2k.f"
			z__4.r = alpha->r * temp1.r - alpha->i * temp1.i, 
				z__4.i = alpha->r * temp1.i + alpha->i * 
				temp1.r;
#line 367 "zsyr2k.f"
			z__2.r = z__3.r + z__4.r, z__2.i = z__3.i + z__4.i;
#line 367 "zsyr2k.f"
			z__5.r = alpha->r * temp2.r - alpha->i * temp2.i, 
				z__5.i = alpha->r * temp2.i + alpha->i * 
				temp2.r;
#line 367 "zsyr2k.f"
			z__1.r = z__2.r + z__5.r, z__1.i = z__2.i + z__5.i;
#line 367 "zsyr2k.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 369 "zsyr2k.f"
		    }
#line 370 "zsyr2k.f"
/* L200: */
#line 370 "zsyr2k.f"
		}
#line 371 "zsyr2k.f"
/* L210: */
#line 371 "zsyr2k.f"
	    }
#line 372 "zsyr2k.f"
	} else {
#line 373 "zsyr2k.f"
	    i__1 = *n;
#line 373 "zsyr2k.f"
	    for (j = 1; j <= i__1; ++j) {
#line 374 "zsyr2k.f"
		i__2 = *n;
#line 374 "zsyr2k.f"
		for (i__ = j; i__ <= i__2; ++i__) {
#line 375 "zsyr2k.f"
		    temp1.r = 0., temp1.i = 0.;
#line 376 "zsyr2k.f"
		    temp2.r = 0., temp2.i = 0.;
#line 377 "zsyr2k.f"
		    i__3 = *k;
#line 377 "zsyr2k.f"
		    for (l = 1; l <= i__3; ++l) {
#line 378 "zsyr2k.f"
			i__4 = l + i__ * a_dim1;
#line 378 "zsyr2k.f"
			i__5 = l + j * b_dim1;
#line 378 "zsyr2k.f"
			z__2.r = a[i__4].r * b[i__5].r - a[i__4].i * b[i__5]
				.i, z__2.i = a[i__4].r * b[i__5].i + a[i__4]
				.i * b[i__5].r;
#line 378 "zsyr2k.f"
			z__1.r = temp1.r + z__2.r, z__1.i = temp1.i + z__2.i;
#line 378 "zsyr2k.f"
			temp1.r = z__1.r, temp1.i = z__1.i;
#line 379 "zsyr2k.f"
			i__4 = l + i__ * b_dim1;
#line 379 "zsyr2k.f"
			i__5 = l + j * a_dim1;
#line 379 "zsyr2k.f"
			z__2.r = b[i__4].r * a[i__5].r - b[i__4].i * a[i__5]
				.i, z__2.i = b[i__4].r * a[i__5].i + b[i__4]
				.i * a[i__5].r;
#line 379 "zsyr2k.f"
			z__1.r = temp2.r + z__2.r, z__1.i = temp2.i + z__2.i;
#line 379 "zsyr2k.f"
			temp2.r = z__1.r, temp2.i = z__1.i;
#line 380 "zsyr2k.f"
/* L220: */
#line 380 "zsyr2k.f"
		    }
#line 381 "zsyr2k.f"
		    if (beta->r == 0. && beta->i == 0.) {
#line 382 "zsyr2k.f"
			i__3 = i__ + j * c_dim1;
#line 382 "zsyr2k.f"
			z__2.r = alpha->r * temp1.r - alpha->i * temp1.i, 
				z__2.i = alpha->r * temp1.i + alpha->i * 
				temp1.r;
#line 382 "zsyr2k.f"
			z__3.r = alpha->r * temp2.r - alpha->i * temp2.i, 
				z__3.i = alpha->r * temp2.i + alpha->i * 
				temp2.r;
#line 382 "zsyr2k.f"
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 382 "zsyr2k.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 383 "zsyr2k.f"
		    } else {
#line 384 "zsyr2k.f"
			i__3 = i__ + j * c_dim1;
#line 384 "zsyr2k.f"
			i__4 = i__ + j * c_dim1;
#line 384 "zsyr2k.f"
			z__3.r = beta->r * c__[i__4].r - beta->i * c__[i__4]
				.i, z__3.i = beta->r * c__[i__4].i + beta->i *
				 c__[i__4].r;
#line 384 "zsyr2k.f"
			z__4.r = alpha->r * temp1.r - alpha->i * temp1.i, 
				z__4.i = alpha->r * temp1.i + alpha->i * 
				temp1.r;
#line 384 "zsyr2k.f"
			z__2.r = z__3.r + z__4.r, z__2.i = z__3.i + z__4.i;
#line 384 "zsyr2k.f"
			z__5.r = alpha->r * temp2.r - alpha->i * temp2.i, 
				z__5.i = alpha->r * temp2.i + alpha->i * 
				temp2.r;
#line 384 "zsyr2k.f"
			z__1.r = z__2.r + z__5.r, z__1.i = z__2.i + z__5.i;
#line 384 "zsyr2k.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 386 "zsyr2k.f"
		    }
#line 387 "zsyr2k.f"
/* L230: */
#line 387 "zsyr2k.f"
		}
#line 388 "zsyr2k.f"
/* L240: */
#line 388 "zsyr2k.f"
	    }
#line 389 "zsyr2k.f"
	}
#line 390 "zsyr2k.f"
    }

#line 392 "zsyr2k.f"
    return 0;

/*     End of ZSYR2K. */

} /* zsyr2k_ */

