#line 1 "zher2k.f"
/* zher2k.f -- translated by f2c (version 20100827).
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

#line 1 "zher2k.f"
/* > \brief \b ZHER2K */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZHER2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC) */

/*       .. Scalar Arguments .. */
/*       COMPLEX*16 ALPHA */
/*       DOUBLE PRECISION BETA */
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
/* > ZHER2K  performs one of the hermitian rank 2k operations */
/* > */
/* >    C := alpha*A*B**H + conjg( alpha )*B*A**H + beta*C, */
/* > */
/* > or */
/* > */
/* >    C := alpha*A**H*B + conjg( alpha )*B**H*A + beta*C, */
/* > */
/* > where  alpha and beta  are scalars with  beta  real,  C is an  n by n */
/* > hermitian matrix and  A and B  are  n by k matrices in the first case */
/* > and  k by n  matrices in the second case. */
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
/* >              TRANS = 'N' or 'n'    C := alpha*A*B**H          + */
/* >                                         conjg( alpha )*B*A**H + */
/* >                                         beta*C. */
/* > */
/* >              TRANS = 'C' or 'c'    C := alpha*A**H*B          + */
/* >                                         conjg( alpha )*B**H*A + */
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
/* >           TRANS = 'C' or 'c',  K  specifies  the number of rows of the */
/* >           matrices  A and B.  K must be at least zero. */
/* > \endverbatim */
/* > */
/* > \param[in] ALPHA */
/* > \verbatim */
/* >          ALPHA is COMPLEX*16 . */
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
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] BETA */
/* > \verbatim */
/* >          BETA is DOUBLE PRECISION . */
/* >           On entry, BETA specifies the scalar beta. */
/* > \endverbatim */
/* > */
/* > \param[in,out] C */
/* > \verbatim */
/* >          C is COMPLEX*16 array of DIMENSION ( LDC, n ). */
/* >           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n */
/* >           upper triangular part of the array C must contain the upper */
/* >           triangular part  of the  hermitian matrix  and the strictly */
/* >           lower triangular part of C is not referenced.  On exit, the */
/* >           upper triangular part of the array  C is overwritten by the */
/* >           upper triangular part of the updated matrix. */
/* >           Before entry  with  UPLO = 'L' or 'l',  the leading  n by n */
/* >           lower triangular part of the array C must contain the lower */
/* >           triangular part  of the  hermitian matrix  and the strictly */
/* >           upper triangular part of C is not referenced.  On exit, the */
/* >           lower triangular part of the array  C is overwritten by the */
/* >           lower triangular part of the updated matrix. */
/* >           Note that the imaginary parts of the diagonal elements need */
/* >           not be set,  they are assumed to be zero,  and on exit they */
/* >           are set to zero. */
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
/* > */
/* >  -- Modified 8-Nov-93 to set C(J,J) to DBLE( C(J,J) ) when BETA = 1. */
/* >     Ed Anderson, Cray Research Inc. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int zher2k_(char *uplo, char *trans, integer *n, integer *k, 
	doublecomplex *alpha, doublecomplex *a, integer *lda, doublecomplex *
	b, integer *ldb, doublereal *beta, doublecomplex *c__, integer *ldc, 
	ftnlen uplo_len, ftnlen trans_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2, 
	    i__3, i__4, i__5, i__6, i__7;
    doublereal d__1;
    doublecomplex z__1, z__2, z__3, z__4, z__5, z__6;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, j, l, info;
    static doublecomplex temp1, temp2;
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

#line 242 "zher2k.f"
    /* Parameter adjustments */
#line 242 "zher2k.f"
    a_dim1 = *lda;
#line 242 "zher2k.f"
    a_offset = 1 + a_dim1;
#line 242 "zher2k.f"
    a -= a_offset;
#line 242 "zher2k.f"
    b_dim1 = *ldb;
#line 242 "zher2k.f"
    b_offset = 1 + b_dim1;
#line 242 "zher2k.f"
    b -= b_offset;
#line 242 "zher2k.f"
    c_dim1 = *ldc;
#line 242 "zher2k.f"
    c_offset = 1 + c_dim1;
#line 242 "zher2k.f"
    c__ -= c_offset;
#line 242 "zher2k.f"

#line 242 "zher2k.f"
    /* Function Body */
#line 242 "zher2k.f"
    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {
#line 243 "zher2k.f"
	nrowa = *n;
#line 244 "zher2k.f"
    } else {
#line 245 "zher2k.f"
	nrowa = *k;
#line 246 "zher2k.f"
    }
#line 247 "zher2k.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);

#line 249 "zher2k.f"
    info = 0;
#line 250 "zher2k.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 251 "zher2k.f"
	info = 1;
#line 252 "zher2k.f"
    } else if (! lsame_(trans, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, 
	    "C", (ftnlen)1, (ftnlen)1)) {
#line 254 "zher2k.f"
	info = 2;
#line 255 "zher2k.f"
    } else if (*n < 0) {
#line 256 "zher2k.f"
	info = 3;
#line 257 "zher2k.f"
    } else if (*k < 0) {
#line 258 "zher2k.f"
	info = 4;
#line 259 "zher2k.f"
    } else if (*lda < max(1,nrowa)) {
#line 260 "zher2k.f"
	info = 7;
#line 261 "zher2k.f"
    } else if (*ldb < max(1,nrowa)) {
#line 262 "zher2k.f"
	info = 9;
#line 263 "zher2k.f"
    } else if (*ldc < max(1,*n)) {
#line 264 "zher2k.f"
	info = 12;
#line 265 "zher2k.f"
    }
#line 266 "zher2k.f"
    if (info != 0) {
#line 267 "zher2k.f"
	xerbla_("ZHER2K", &info, (ftnlen)6);
#line 268 "zher2k.f"
	return 0;
#line 269 "zher2k.f"
    }

/*     Quick return if possible. */

#line 273 "zher2k.f"
    if (*n == 0 || (alpha->r == 0. && alpha->i == 0. || *k == 0) && *beta == 
	    1.) {
#line 273 "zher2k.f"
	return 0;
#line 273 "zher2k.f"
    }

/*     And when  alpha.eq.zero. */

#line 278 "zher2k.f"
    if (alpha->r == 0. && alpha->i == 0.) {
#line 279 "zher2k.f"
	if (upper) {
#line 280 "zher2k.f"
	    if (*beta == 0.) {
#line 281 "zher2k.f"
		i__1 = *n;
#line 281 "zher2k.f"
		for (j = 1; j <= i__1; ++j) {
#line 282 "zher2k.f"
		    i__2 = j;
#line 282 "zher2k.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 283 "zher2k.f"
			i__3 = i__ + j * c_dim1;
#line 283 "zher2k.f"
			c__[i__3].r = 0., c__[i__3].i = 0.;
#line 284 "zher2k.f"
/* L10: */
#line 284 "zher2k.f"
		    }
#line 285 "zher2k.f"
/* L20: */
#line 285 "zher2k.f"
		}
#line 286 "zher2k.f"
	    } else {
#line 287 "zher2k.f"
		i__1 = *n;
#line 287 "zher2k.f"
		for (j = 1; j <= i__1; ++j) {
#line 288 "zher2k.f"
		    i__2 = j - 1;
#line 288 "zher2k.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 289 "zher2k.f"
			i__3 = i__ + j * c_dim1;
#line 289 "zher2k.f"
			i__4 = i__ + j * c_dim1;
#line 289 "zher2k.f"
			z__1.r = *beta * c__[i__4].r, z__1.i = *beta * c__[
				i__4].i;
#line 289 "zher2k.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 290 "zher2k.f"
/* L30: */
#line 290 "zher2k.f"
		    }
#line 291 "zher2k.f"
		    i__2 = j + j * c_dim1;
#line 291 "zher2k.f"
		    i__3 = j + j * c_dim1;
#line 291 "zher2k.f"
		    d__1 = *beta * c__[i__3].r;
#line 291 "zher2k.f"
		    c__[i__2].r = d__1, c__[i__2].i = 0.;
#line 292 "zher2k.f"
/* L40: */
#line 292 "zher2k.f"
		}
#line 293 "zher2k.f"
	    }
#line 294 "zher2k.f"
	} else {
#line 295 "zher2k.f"
	    if (*beta == 0.) {
#line 296 "zher2k.f"
		i__1 = *n;
#line 296 "zher2k.f"
		for (j = 1; j <= i__1; ++j) {
#line 297 "zher2k.f"
		    i__2 = *n;
#line 297 "zher2k.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 298 "zher2k.f"
			i__3 = i__ + j * c_dim1;
#line 298 "zher2k.f"
			c__[i__3].r = 0., c__[i__3].i = 0.;
#line 299 "zher2k.f"
/* L50: */
#line 299 "zher2k.f"
		    }
#line 300 "zher2k.f"
/* L60: */
#line 300 "zher2k.f"
		}
#line 301 "zher2k.f"
	    } else {
#line 302 "zher2k.f"
		i__1 = *n;
#line 302 "zher2k.f"
		for (j = 1; j <= i__1; ++j) {
#line 303 "zher2k.f"
		    i__2 = j + j * c_dim1;
#line 303 "zher2k.f"
		    i__3 = j + j * c_dim1;
#line 303 "zher2k.f"
		    d__1 = *beta * c__[i__3].r;
#line 303 "zher2k.f"
		    c__[i__2].r = d__1, c__[i__2].i = 0.;
#line 304 "zher2k.f"
		    i__2 = *n;
#line 304 "zher2k.f"
		    for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 305 "zher2k.f"
			i__3 = i__ + j * c_dim1;
#line 305 "zher2k.f"
			i__4 = i__ + j * c_dim1;
#line 305 "zher2k.f"
			z__1.r = *beta * c__[i__4].r, z__1.i = *beta * c__[
				i__4].i;
#line 305 "zher2k.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 306 "zher2k.f"
/* L70: */
#line 306 "zher2k.f"
		    }
#line 307 "zher2k.f"
/* L80: */
#line 307 "zher2k.f"
		}
#line 308 "zher2k.f"
	    }
#line 309 "zher2k.f"
	}
#line 310 "zher2k.f"
	return 0;
#line 311 "zher2k.f"
    }

/*     Start the operations. */

#line 315 "zher2k.f"
    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {

/*        Form  C := alpha*A*B**H + conjg( alpha )*B*A**H + */
/*                   C. */

#line 320 "zher2k.f"
	if (upper) {
#line 321 "zher2k.f"
	    i__1 = *n;
#line 321 "zher2k.f"
	    for (j = 1; j <= i__1; ++j) {
#line 322 "zher2k.f"
		if (*beta == 0.) {
#line 323 "zher2k.f"
		    i__2 = j;
#line 323 "zher2k.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 324 "zher2k.f"
			i__3 = i__ + j * c_dim1;
#line 324 "zher2k.f"
			c__[i__3].r = 0., c__[i__3].i = 0.;
#line 325 "zher2k.f"
/* L90: */
#line 325 "zher2k.f"
		    }
#line 326 "zher2k.f"
		} else if (*beta != 1.) {
#line 327 "zher2k.f"
		    i__2 = j - 1;
#line 327 "zher2k.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 328 "zher2k.f"
			i__3 = i__ + j * c_dim1;
#line 328 "zher2k.f"
			i__4 = i__ + j * c_dim1;
#line 328 "zher2k.f"
			z__1.r = *beta * c__[i__4].r, z__1.i = *beta * c__[
				i__4].i;
#line 328 "zher2k.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 329 "zher2k.f"
/* L100: */
#line 329 "zher2k.f"
		    }
#line 330 "zher2k.f"
		    i__2 = j + j * c_dim1;
#line 330 "zher2k.f"
		    i__3 = j + j * c_dim1;
#line 330 "zher2k.f"
		    d__1 = *beta * c__[i__3].r;
#line 330 "zher2k.f"
		    c__[i__2].r = d__1, c__[i__2].i = 0.;
#line 331 "zher2k.f"
		} else {
#line 332 "zher2k.f"
		    i__2 = j + j * c_dim1;
#line 332 "zher2k.f"
		    i__3 = j + j * c_dim1;
#line 332 "zher2k.f"
		    d__1 = c__[i__3].r;
#line 332 "zher2k.f"
		    c__[i__2].r = d__1, c__[i__2].i = 0.;
#line 333 "zher2k.f"
		}
#line 334 "zher2k.f"
		i__2 = *k;
#line 334 "zher2k.f"
		for (l = 1; l <= i__2; ++l) {
#line 335 "zher2k.f"
		    i__3 = j + l * a_dim1;
#line 335 "zher2k.f"
		    i__4 = j + l * b_dim1;
#line 335 "zher2k.f"
		    if (a[i__3].r != 0. || a[i__3].i != 0. || (b[i__4].r != 
			    0. || b[i__4].i != 0.)) {
#line 336 "zher2k.f"
			d_cnjg(&z__2, &b[j + l * b_dim1]);
#line 336 "zher2k.f"
			z__1.r = alpha->r * z__2.r - alpha->i * z__2.i, 
				z__1.i = alpha->r * z__2.i + alpha->i * 
				z__2.r;
#line 336 "zher2k.f"
			temp1.r = z__1.r, temp1.i = z__1.i;
#line 337 "zher2k.f"
			i__3 = j + l * a_dim1;
#line 337 "zher2k.f"
			z__2.r = alpha->r * a[i__3].r - alpha->i * a[i__3].i, 
				z__2.i = alpha->r * a[i__3].i + alpha->i * a[
				i__3].r;
#line 337 "zher2k.f"
			d_cnjg(&z__1, &z__2);
#line 337 "zher2k.f"
			temp2.r = z__1.r, temp2.i = z__1.i;
#line 338 "zher2k.f"
			i__3 = j - 1;
#line 338 "zher2k.f"
			for (i__ = 1; i__ <= i__3; ++i__) {
#line 339 "zher2k.f"
			    i__4 = i__ + j * c_dim1;
#line 339 "zher2k.f"
			    i__5 = i__ + j * c_dim1;
#line 339 "zher2k.f"
			    i__6 = i__ + l * a_dim1;
#line 339 "zher2k.f"
			    z__3.r = a[i__6].r * temp1.r - a[i__6].i * 
				    temp1.i, z__3.i = a[i__6].r * temp1.i + a[
				    i__6].i * temp1.r;
#line 339 "zher2k.f"
			    z__2.r = c__[i__5].r + z__3.r, z__2.i = c__[i__5]
				    .i + z__3.i;
#line 339 "zher2k.f"
			    i__7 = i__ + l * b_dim1;
#line 339 "zher2k.f"
			    z__4.r = b[i__7].r * temp2.r - b[i__7].i * 
				    temp2.i, z__4.i = b[i__7].r * temp2.i + b[
				    i__7].i * temp2.r;
#line 339 "zher2k.f"
			    z__1.r = z__2.r + z__4.r, z__1.i = z__2.i + 
				    z__4.i;
#line 339 "zher2k.f"
			    c__[i__4].r = z__1.r, c__[i__4].i = z__1.i;
#line 341 "zher2k.f"
/* L110: */
#line 341 "zher2k.f"
			}
#line 342 "zher2k.f"
			i__3 = j + j * c_dim1;
#line 342 "zher2k.f"
			i__4 = j + j * c_dim1;
#line 342 "zher2k.f"
			i__5 = j + l * a_dim1;
#line 342 "zher2k.f"
			z__2.r = a[i__5].r * temp1.r - a[i__5].i * temp1.i, 
				z__2.i = a[i__5].r * temp1.i + a[i__5].i * 
				temp1.r;
#line 342 "zher2k.f"
			i__6 = j + l * b_dim1;
#line 342 "zher2k.f"
			z__3.r = b[i__6].r * temp2.r - b[i__6].i * temp2.i, 
				z__3.i = b[i__6].r * temp2.i + b[i__6].i * 
				temp2.r;
#line 342 "zher2k.f"
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 342 "zher2k.f"
			d__1 = c__[i__4].r + z__1.r;
#line 342 "zher2k.f"
			c__[i__3].r = d__1, c__[i__3].i = 0.;
#line 344 "zher2k.f"
		    }
#line 345 "zher2k.f"
/* L120: */
#line 345 "zher2k.f"
		}
#line 346 "zher2k.f"
/* L130: */
#line 346 "zher2k.f"
	    }
#line 347 "zher2k.f"
	} else {
#line 348 "zher2k.f"
	    i__1 = *n;
#line 348 "zher2k.f"
	    for (j = 1; j <= i__1; ++j) {
#line 349 "zher2k.f"
		if (*beta == 0.) {
#line 350 "zher2k.f"
		    i__2 = *n;
#line 350 "zher2k.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 351 "zher2k.f"
			i__3 = i__ + j * c_dim1;
#line 351 "zher2k.f"
			c__[i__3].r = 0., c__[i__3].i = 0.;
#line 352 "zher2k.f"
/* L140: */
#line 352 "zher2k.f"
		    }
#line 353 "zher2k.f"
		} else if (*beta != 1.) {
#line 354 "zher2k.f"
		    i__2 = *n;
#line 354 "zher2k.f"
		    for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 355 "zher2k.f"
			i__3 = i__ + j * c_dim1;
#line 355 "zher2k.f"
			i__4 = i__ + j * c_dim1;
#line 355 "zher2k.f"
			z__1.r = *beta * c__[i__4].r, z__1.i = *beta * c__[
				i__4].i;
#line 355 "zher2k.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 356 "zher2k.f"
/* L150: */
#line 356 "zher2k.f"
		    }
#line 357 "zher2k.f"
		    i__2 = j + j * c_dim1;
#line 357 "zher2k.f"
		    i__3 = j + j * c_dim1;
#line 357 "zher2k.f"
		    d__1 = *beta * c__[i__3].r;
#line 357 "zher2k.f"
		    c__[i__2].r = d__1, c__[i__2].i = 0.;
#line 358 "zher2k.f"
		} else {
#line 359 "zher2k.f"
		    i__2 = j + j * c_dim1;
#line 359 "zher2k.f"
		    i__3 = j + j * c_dim1;
#line 359 "zher2k.f"
		    d__1 = c__[i__3].r;
#line 359 "zher2k.f"
		    c__[i__2].r = d__1, c__[i__2].i = 0.;
#line 360 "zher2k.f"
		}
#line 361 "zher2k.f"
		i__2 = *k;
#line 361 "zher2k.f"
		for (l = 1; l <= i__2; ++l) {
#line 362 "zher2k.f"
		    i__3 = j + l * a_dim1;
#line 362 "zher2k.f"
		    i__4 = j + l * b_dim1;
#line 362 "zher2k.f"
		    if (a[i__3].r != 0. || a[i__3].i != 0. || (b[i__4].r != 
			    0. || b[i__4].i != 0.)) {
#line 363 "zher2k.f"
			d_cnjg(&z__2, &b[j + l * b_dim1]);
#line 363 "zher2k.f"
			z__1.r = alpha->r * z__2.r - alpha->i * z__2.i, 
				z__1.i = alpha->r * z__2.i + alpha->i * 
				z__2.r;
#line 363 "zher2k.f"
			temp1.r = z__1.r, temp1.i = z__1.i;
#line 364 "zher2k.f"
			i__3 = j + l * a_dim1;
#line 364 "zher2k.f"
			z__2.r = alpha->r * a[i__3].r - alpha->i * a[i__3].i, 
				z__2.i = alpha->r * a[i__3].i + alpha->i * a[
				i__3].r;
#line 364 "zher2k.f"
			d_cnjg(&z__1, &z__2);
#line 364 "zher2k.f"
			temp2.r = z__1.r, temp2.i = z__1.i;
#line 365 "zher2k.f"
			i__3 = *n;
#line 365 "zher2k.f"
			for (i__ = j + 1; i__ <= i__3; ++i__) {
#line 366 "zher2k.f"
			    i__4 = i__ + j * c_dim1;
#line 366 "zher2k.f"
			    i__5 = i__ + j * c_dim1;
#line 366 "zher2k.f"
			    i__6 = i__ + l * a_dim1;
#line 366 "zher2k.f"
			    z__3.r = a[i__6].r * temp1.r - a[i__6].i * 
				    temp1.i, z__3.i = a[i__6].r * temp1.i + a[
				    i__6].i * temp1.r;
#line 366 "zher2k.f"
			    z__2.r = c__[i__5].r + z__3.r, z__2.i = c__[i__5]
				    .i + z__3.i;
#line 366 "zher2k.f"
			    i__7 = i__ + l * b_dim1;
#line 366 "zher2k.f"
			    z__4.r = b[i__7].r * temp2.r - b[i__7].i * 
				    temp2.i, z__4.i = b[i__7].r * temp2.i + b[
				    i__7].i * temp2.r;
#line 366 "zher2k.f"
			    z__1.r = z__2.r + z__4.r, z__1.i = z__2.i + 
				    z__4.i;
#line 366 "zher2k.f"
			    c__[i__4].r = z__1.r, c__[i__4].i = z__1.i;
#line 368 "zher2k.f"
/* L160: */
#line 368 "zher2k.f"
			}
#line 369 "zher2k.f"
			i__3 = j + j * c_dim1;
#line 369 "zher2k.f"
			i__4 = j + j * c_dim1;
#line 369 "zher2k.f"
			i__5 = j + l * a_dim1;
#line 369 "zher2k.f"
			z__2.r = a[i__5].r * temp1.r - a[i__5].i * temp1.i, 
				z__2.i = a[i__5].r * temp1.i + a[i__5].i * 
				temp1.r;
#line 369 "zher2k.f"
			i__6 = j + l * b_dim1;
#line 369 "zher2k.f"
			z__3.r = b[i__6].r * temp2.r - b[i__6].i * temp2.i, 
				z__3.i = b[i__6].r * temp2.i + b[i__6].i * 
				temp2.r;
#line 369 "zher2k.f"
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 369 "zher2k.f"
			d__1 = c__[i__4].r + z__1.r;
#line 369 "zher2k.f"
			c__[i__3].r = d__1, c__[i__3].i = 0.;
#line 371 "zher2k.f"
		    }
#line 372 "zher2k.f"
/* L170: */
#line 372 "zher2k.f"
		}
#line 373 "zher2k.f"
/* L180: */
#line 373 "zher2k.f"
	    }
#line 374 "zher2k.f"
	}
#line 375 "zher2k.f"
    } else {

/*        Form  C := alpha*A**H*B + conjg( alpha )*B**H*A + */
/*                   C. */

#line 380 "zher2k.f"
	if (upper) {
#line 381 "zher2k.f"
	    i__1 = *n;
#line 381 "zher2k.f"
	    for (j = 1; j <= i__1; ++j) {
#line 382 "zher2k.f"
		i__2 = j;
#line 382 "zher2k.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 383 "zher2k.f"
		    temp1.r = 0., temp1.i = 0.;
#line 384 "zher2k.f"
		    temp2.r = 0., temp2.i = 0.;
#line 385 "zher2k.f"
		    i__3 = *k;
#line 385 "zher2k.f"
		    for (l = 1; l <= i__3; ++l) {
#line 386 "zher2k.f"
			d_cnjg(&z__3, &a[l + i__ * a_dim1]);
#line 386 "zher2k.f"
			i__4 = l + j * b_dim1;
#line 386 "zher2k.f"
			z__2.r = z__3.r * b[i__4].r - z__3.i * b[i__4].i, 
				z__2.i = z__3.r * b[i__4].i + z__3.i * b[i__4]
				.r;
#line 386 "zher2k.f"
			z__1.r = temp1.r + z__2.r, z__1.i = temp1.i + z__2.i;
#line 386 "zher2k.f"
			temp1.r = z__1.r, temp1.i = z__1.i;
#line 387 "zher2k.f"
			d_cnjg(&z__3, &b[l + i__ * b_dim1]);
#line 387 "zher2k.f"
			i__4 = l + j * a_dim1;
#line 387 "zher2k.f"
			z__2.r = z__3.r * a[i__4].r - z__3.i * a[i__4].i, 
				z__2.i = z__3.r * a[i__4].i + z__3.i * a[i__4]
				.r;
#line 387 "zher2k.f"
			z__1.r = temp2.r + z__2.r, z__1.i = temp2.i + z__2.i;
#line 387 "zher2k.f"
			temp2.r = z__1.r, temp2.i = z__1.i;
#line 388 "zher2k.f"
/* L190: */
#line 388 "zher2k.f"
		    }
#line 389 "zher2k.f"
		    if (i__ == j) {
#line 390 "zher2k.f"
			if (*beta == 0.) {
#line 391 "zher2k.f"
			    i__3 = j + j * c_dim1;
#line 391 "zher2k.f"
			    z__2.r = alpha->r * temp1.r - alpha->i * temp1.i, 
				    z__2.i = alpha->r * temp1.i + alpha->i * 
				    temp1.r;
#line 391 "zher2k.f"
			    d_cnjg(&z__4, alpha);
#line 391 "zher2k.f"
			    z__3.r = z__4.r * temp2.r - z__4.i * temp2.i, 
				    z__3.i = z__4.r * temp2.i + z__4.i * 
				    temp2.r;
#line 391 "zher2k.f"
			    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + 
				    z__3.i;
#line 391 "zher2k.f"
			    d__1 = z__1.r;
#line 391 "zher2k.f"
			    c__[i__3].r = d__1, c__[i__3].i = 0.;
#line 393 "zher2k.f"
			} else {
#line 394 "zher2k.f"
			    i__3 = j + j * c_dim1;
#line 394 "zher2k.f"
			    i__4 = j + j * c_dim1;
#line 394 "zher2k.f"
			    z__2.r = alpha->r * temp1.r - alpha->i * temp1.i, 
				    z__2.i = alpha->r * temp1.i + alpha->i * 
				    temp1.r;
#line 394 "zher2k.f"
			    d_cnjg(&z__4, alpha);
#line 394 "zher2k.f"
			    z__3.r = z__4.r * temp2.r - z__4.i * temp2.i, 
				    z__3.i = z__4.r * temp2.i + z__4.i * 
				    temp2.r;
#line 394 "zher2k.f"
			    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + 
				    z__3.i;
#line 394 "zher2k.f"
			    d__1 = *beta * c__[i__4].r + z__1.r;
#line 394 "zher2k.f"
			    c__[i__3].r = d__1, c__[i__3].i = 0.;
#line 397 "zher2k.f"
			}
#line 398 "zher2k.f"
		    } else {
#line 399 "zher2k.f"
			if (*beta == 0.) {
#line 400 "zher2k.f"
			    i__3 = i__ + j * c_dim1;
#line 400 "zher2k.f"
			    z__2.r = alpha->r * temp1.r - alpha->i * temp1.i, 
				    z__2.i = alpha->r * temp1.i + alpha->i * 
				    temp1.r;
#line 400 "zher2k.f"
			    d_cnjg(&z__4, alpha);
#line 400 "zher2k.f"
			    z__3.r = z__4.r * temp2.r - z__4.i * temp2.i, 
				    z__3.i = z__4.r * temp2.i + z__4.i * 
				    temp2.r;
#line 400 "zher2k.f"
			    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + 
				    z__3.i;
#line 400 "zher2k.f"
			    c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 401 "zher2k.f"
			} else {
#line 402 "zher2k.f"
			    i__3 = i__ + j * c_dim1;
#line 402 "zher2k.f"
			    i__4 = i__ + j * c_dim1;
#line 402 "zher2k.f"
			    z__3.r = *beta * c__[i__4].r, z__3.i = *beta * 
				    c__[i__4].i;
#line 402 "zher2k.f"
			    z__4.r = alpha->r * temp1.r - alpha->i * temp1.i, 
				    z__4.i = alpha->r * temp1.i + alpha->i * 
				    temp1.r;
#line 402 "zher2k.f"
			    z__2.r = z__3.r + z__4.r, z__2.i = z__3.i + 
				    z__4.i;
#line 402 "zher2k.f"
			    d_cnjg(&z__6, alpha);
#line 402 "zher2k.f"
			    z__5.r = z__6.r * temp2.r - z__6.i * temp2.i, 
				    z__5.i = z__6.r * temp2.i + z__6.i * 
				    temp2.r;
#line 402 "zher2k.f"
			    z__1.r = z__2.r + z__5.r, z__1.i = z__2.i + 
				    z__5.i;
#line 402 "zher2k.f"
			    c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 404 "zher2k.f"
			}
#line 405 "zher2k.f"
		    }
#line 406 "zher2k.f"
/* L200: */
#line 406 "zher2k.f"
		}
#line 407 "zher2k.f"
/* L210: */
#line 407 "zher2k.f"
	    }
#line 408 "zher2k.f"
	} else {
#line 409 "zher2k.f"
	    i__1 = *n;
#line 409 "zher2k.f"
	    for (j = 1; j <= i__1; ++j) {
#line 410 "zher2k.f"
		i__2 = *n;
#line 410 "zher2k.f"
		for (i__ = j; i__ <= i__2; ++i__) {
#line 411 "zher2k.f"
		    temp1.r = 0., temp1.i = 0.;
#line 412 "zher2k.f"
		    temp2.r = 0., temp2.i = 0.;
#line 413 "zher2k.f"
		    i__3 = *k;
#line 413 "zher2k.f"
		    for (l = 1; l <= i__3; ++l) {
#line 414 "zher2k.f"
			d_cnjg(&z__3, &a[l + i__ * a_dim1]);
#line 414 "zher2k.f"
			i__4 = l + j * b_dim1;
#line 414 "zher2k.f"
			z__2.r = z__3.r * b[i__4].r - z__3.i * b[i__4].i, 
				z__2.i = z__3.r * b[i__4].i + z__3.i * b[i__4]
				.r;
#line 414 "zher2k.f"
			z__1.r = temp1.r + z__2.r, z__1.i = temp1.i + z__2.i;
#line 414 "zher2k.f"
			temp1.r = z__1.r, temp1.i = z__1.i;
#line 415 "zher2k.f"
			d_cnjg(&z__3, &b[l + i__ * b_dim1]);
#line 415 "zher2k.f"
			i__4 = l + j * a_dim1;
#line 415 "zher2k.f"
			z__2.r = z__3.r * a[i__4].r - z__3.i * a[i__4].i, 
				z__2.i = z__3.r * a[i__4].i + z__3.i * a[i__4]
				.r;
#line 415 "zher2k.f"
			z__1.r = temp2.r + z__2.r, z__1.i = temp2.i + z__2.i;
#line 415 "zher2k.f"
			temp2.r = z__1.r, temp2.i = z__1.i;
#line 416 "zher2k.f"
/* L220: */
#line 416 "zher2k.f"
		    }
#line 417 "zher2k.f"
		    if (i__ == j) {
#line 418 "zher2k.f"
			if (*beta == 0.) {
#line 419 "zher2k.f"
			    i__3 = j + j * c_dim1;
#line 419 "zher2k.f"
			    z__2.r = alpha->r * temp1.r - alpha->i * temp1.i, 
				    z__2.i = alpha->r * temp1.i + alpha->i * 
				    temp1.r;
#line 419 "zher2k.f"
			    d_cnjg(&z__4, alpha);
#line 419 "zher2k.f"
			    z__3.r = z__4.r * temp2.r - z__4.i * temp2.i, 
				    z__3.i = z__4.r * temp2.i + z__4.i * 
				    temp2.r;
#line 419 "zher2k.f"
			    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + 
				    z__3.i;
#line 419 "zher2k.f"
			    d__1 = z__1.r;
#line 419 "zher2k.f"
			    c__[i__3].r = d__1, c__[i__3].i = 0.;
#line 421 "zher2k.f"
			} else {
#line 422 "zher2k.f"
			    i__3 = j + j * c_dim1;
#line 422 "zher2k.f"
			    i__4 = j + j * c_dim1;
#line 422 "zher2k.f"
			    z__2.r = alpha->r * temp1.r - alpha->i * temp1.i, 
				    z__2.i = alpha->r * temp1.i + alpha->i * 
				    temp1.r;
#line 422 "zher2k.f"
			    d_cnjg(&z__4, alpha);
#line 422 "zher2k.f"
			    z__3.r = z__4.r * temp2.r - z__4.i * temp2.i, 
				    z__3.i = z__4.r * temp2.i + z__4.i * 
				    temp2.r;
#line 422 "zher2k.f"
			    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + 
				    z__3.i;
#line 422 "zher2k.f"
			    d__1 = *beta * c__[i__4].r + z__1.r;
#line 422 "zher2k.f"
			    c__[i__3].r = d__1, c__[i__3].i = 0.;
#line 425 "zher2k.f"
			}
#line 426 "zher2k.f"
		    } else {
#line 427 "zher2k.f"
			if (*beta == 0.) {
#line 428 "zher2k.f"
			    i__3 = i__ + j * c_dim1;
#line 428 "zher2k.f"
			    z__2.r = alpha->r * temp1.r - alpha->i * temp1.i, 
				    z__2.i = alpha->r * temp1.i + alpha->i * 
				    temp1.r;
#line 428 "zher2k.f"
			    d_cnjg(&z__4, alpha);
#line 428 "zher2k.f"
			    z__3.r = z__4.r * temp2.r - z__4.i * temp2.i, 
				    z__3.i = z__4.r * temp2.i + z__4.i * 
				    temp2.r;
#line 428 "zher2k.f"
			    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + 
				    z__3.i;
#line 428 "zher2k.f"
			    c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 429 "zher2k.f"
			} else {
#line 430 "zher2k.f"
			    i__3 = i__ + j * c_dim1;
#line 430 "zher2k.f"
			    i__4 = i__ + j * c_dim1;
#line 430 "zher2k.f"
			    z__3.r = *beta * c__[i__4].r, z__3.i = *beta * 
				    c__[i__4].i;
#line 430 "zher2k.f"
			    z__4.r = alpha->r * temp1.r - alpha->i * temp1.i, 
				    z__4.i = alpha->r * temp1.i + alpha->i * 
				    temp1.r;
#line 430 "zher2k.f"
			    z__2.r = z__3.r + z__4.r, z__2.i = z__3.i + 
				    z__4.i;
#line 430 "zher2k.f"
			    d_cnjg(&z__6, alpha);
#line 430 "zher2k.f"
			    z__5.r = z__6.r * temp2.r - z__6.i * temp2.i, 
				    z__5.i = z__6.r * temp2.i + z__6.i * 
				    temp2.r;
#line 430 "zher2k.f"
			    z__1.r = z__2.r + z__5.r, z__1.i = z__2.i + 
				    z__5.i;
#line 430 "zher2k.f"
			    c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 432 "zher2k.f"
			}
#line 433 "zher2k.f"
		    }
#line 434 "zher2k.f"
/* L230: */
#line 434 "zher2k.f"
		}
#line 435 "zher2k.f"
/* L240: */
#line 435 "zher2k.f"
	    }
#line 436 "zher2k.f"
	}
#line 437 "zher2k.f"
    }

#line 439 "zher2k.f"
    return 0;

/*     End of ZHER2K. */

} /* zher2k_ */

