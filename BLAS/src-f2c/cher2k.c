#line 1 "cher2k.f"
/* cher2k.f -- translated by f2c (version 20100827).
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

#line 1 "cher2k.f"
/* > \brief \b CHER2K */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CHER2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC) */

/*       .. Scalar Arguments .. */
/*       COMPLEX ALPHA */
/*       REAL BETA */
/*       INTEGER K,LDA,LDB,LDC,N */
/*       CHARACTER TRANS,UPLO */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX A(LDA,*),B(LDB,*),C(LDC,*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CHER2K  performs one of the hermitian rank 2k operations */
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
/* >          ALPHA is COMPLEX */
/* >           On entry, ALPHA specifies the scalar alpha. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is COMPLEX array, dimension ( LDA, ka ), where ka is */
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
/* >          B is COMPLEX array, dimension ( LDB, kb ), where kb is */
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
/* >          C is COMPLEX array, dimension ( LDC, N ) */
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

/* > \ingroup complex_blas_level3 */

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
/* >  -- Modified 8-Nov-93 to set C(J,J) to REAL( C(J,J) ) when BETA = 1. */
/* >     Ed Anderson, Cray Research Inc. */
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int cher2k_(char *uplo, char *trans, integer *n, integer *k, 
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

#line 241 "cher2k.f"
    /* Parameter adjustments */
#line 241 "cher2k.f"
    a_dim1 = *lda;
#line 241 "cher2k.f"
    a_offset = 1 + a_dim1;
#line 241 "cher2k.f"
    a -= a_offset;
#line 241 "cher2k.f"
    b_dim1 = *ldb;
#line 241 "cher2k.f"
    b_offset = 1 + b_dim1;
#line 241 "cher2k.f"
    b -= b_offset;
#line 241 "cher2k.f"
    c_dim1 = *ldc;
#line 241 "cher2k.f"
    c_offset = 1 + c_dim1;
#line 241 "cher2k.f"
    c__ -= c_offset;
#line 241 "cher2k.f"

#line 241 "cher2k.f"
    /* Function Body */
#line 241 "cher2k.f"
    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {
#line 242 "cher2k.f"
	nrowa = *n;
#line 243 "cher2k.f"
    } else {
#line 244 "cher2k.f"
	nrowa = *k;
#line 245 "cher2k.f"
    }
#line 246 "cher2k.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);

#line 248 "cher2k.f"
    info = 0;
#line 249 "cher2k.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 250 "cher2k.f"
	info = 1;
#line 251 "cher2k.f"
    } else if (! lsame_(trans, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, 
	    "C", (ftnlen)1, (ftnlen)1)) {
#line 253 "cher2k.f"
	info = 2;
#line 254 "cher2k.f"
    } else if (*n < 0) {
#line 255 "cher2k.f"
	info = 3;
#line 256 "cher2k.f"
    } else if (*k < 0) {
#line 257 "cher2k.f"
	info = 4;
#line 258 "cher2k.f"
    } else if (*lda < max(1,nrowa)) {
#line 259 "cher2k.f"
	info = 7;
#line 260 "cher2k.f"
    } else if (*ldb < max(1,nrowa)) {
#line 261 "cher2k.f"
	info = 9;
#line 262 "cher2k.f"
    } else if (*ldc < max(1,*n)) {
#line 263 "cher2k.f"
	info = 12;
#line 264 "cher2k.f"
    }
#line 265 "cher2k.f"
    if (info != 0) {
#line 266 "cher2k.f"
	xerbla_("CHER2K", &info, (ftnlen)6);
#line 267 "cher2k.f"
	return 0;
#line 268 "cher2k.f"
    }

/*     Quick return if possible. */

#line 272 "cher2k.f"
    if (*n == 0 || (alpha->r == 0. && alpha->i == 0. || *k == 0) && *beta == 
	    1.) {
#line 272 "cher2k.f"
	return 0;
#line 272 "cher2k.f"
    }

/*     And when  alpha.eq.zero. */

#line 277 "cher2k.f"
    if (alpha->r == 0. && alpha->i == 0.) {
#line 278 "cher2k.f"
	if (upper) {
#line 279 "cher2k.f"
	    if (*beta == 0.) {
#line 280 "cher2k.f"
		i__1 = *n;
#line 280 "cher2k.f"
		for (j = 1; j <= i__1; ++j) {
#line 281 "cher2k.f"
		    i__2 = j;
#line 281 "cher2k.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 282 "cher2k.f"
			i__3 = i__ + j * c_dim1;
#line 282 "cher2k.f"
			c__[i__3].r = 0., c__[i__3].i = 0.;
#line 283 "cher2k.f"
/* L10: */
#line 283 "cher2k.f"
		    }
#line 284 "cher2k.f"
/* L20: */
#line 284 "cher2k.f"
		}
#line 285 "cher2k.f"
	    } else {
#line 286 "cher2k.f"
		i__1 = *n;
#line 286 "cher2k.f"
		for (j = 1; j <= i__1; ++j) {
#line 287 "cher2k.f"
		    i__2 = j - 1;
#line 287 "cher2k.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 288 "cher2k.f"
			i__3 = i__ + j * c_dim1;
#line 288 "cher2k.f"
			i__4 = i__ + j * c_dim1;
#line 288 "cher2k.f"
			z__1.r = *beta * c__[i__4].r, z__1.i = *beta * c__[
				i__4].i;
#line 288 "cher2k.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 289 "cher2k.f"
/* L30: */
#line 289 "cher2k.f"
		    }
#line 290 "cher2k.f"
		    i__2 = j + j * c_dim1;
#line 290 "cher2k.f"
		    i__3 = j + j * c_dim1;
#line 290 "cher2k.f"
		    d__1 = *beta * c__[i__3].r;
#line 290 "cher2k.f"
		    c__[i__2].r = d__1, c__[i__2].i = 0.;
#line 291 "cher2k.f"
/* L40: */
#line 291 "cher2k.f"
		}
#line 292 "cher2k.f"
	    }
#line 293 "cher2k.f"
	} else {
#line 294 "cher2k.f"
	    if (*beta == 0.) {
#line 295 "cher2k.f"
		i__1 = *n;
#line 295 "cher2k.f"
		for (j = 1; j <= i__1; ++j) {
#line 296 "cher2k.f"
		    i__2 = *n;
#line 296 "cher2k.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 297 "cher2k.f"
			i__3 = i__ + j * c_dim1;
#line 297 "cher2k.f"
			c__[i__3].r = 0., c__[i__3].i = 0.;
#line 298 "cher2k.f"
/* L50: */
#line 298 "cher2k.f"
		    }
#line 299 "cher2k.f"
/* L60: */
#line 299 "cher2k.f"
		}
#line 300 "cher2k.f"
	    } else {
#line 301 "cher2k.f"
		i__1 = *n;
#line 301 "cher2k.f"
		for (j = 1; j <= i__1; ++j) {
#line 302 "cher2k.f"
		    i__2 = j + j * c_dim1;
#line 302 "cher2k.f"
		    i__3 = j + j * c_dim1;
#line 302 "cher2k.f"
		    d__1 = *beta * c__[i__3].r;
#line 302 "cher2k.f"
		    c__[i__2].r = d__1, c__[i__2].i = 0.;
#line 303 "cher2k.f"
		    i__2 = *n;
#line 303 "cher2k.f"
		    for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 304 "cher2k.f"
			i__3 = i__ + j * c_dim1;
#line 304 "cher2k.f"
			i__4 = i__ + j * c_dim1;
#line 304 "cher2k.f"
			z__1.r = *beta * c__[i__4].r, z__1.i = *beta * c__[
				i__4].i;
#line 304 "cher2k.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 305 "cher2k.f"
/* L70: */
#line 305 "cher2k.f"
		    }
#line 306 "cher2k.f"
/* L80: */
#line 306 "cher2k.f"
		}
#line 307 "cher2k.f"
	    }
#line 308 "cher2k.f"
	}
#line 309 "cher2k.f"
	return 0;
#line 310 "cher2k.f"
    }

/*     Start the operations. */

#line 314 "cher2k.f"
    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {

/*        Form  C := alpha*A*B**H + conjg( alpha )*B*A**H + */
/*                   C. */

#line 319 "cher2k.f"
	if (upper) {
#line 320 "cher2k.f"
	    i__1 = *n;
#line 320 "cher2k.f"
	    for (j = 1; j <= i__1; ++j) {
#line 321 "cher2k.f"
		if (*beta == 0.) {
#line 322 "cher2k.f"
		    i__2 = j;
#line 322 "cher2k.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 323 "cher2k.f"
			i__3 = i__ + j * c_dim1;
#line 323 "cher2k.f"
			c__[i__3].r = 0., c__[i__3].i = 0.;
#line 324 "cher2k.f"
/* L90: */
#line 324 "cher2k.f"
		    }
#line 325 "cher2k.f"
		} else if (*beta != 1.) {
#line 326 "cher2k.f"
		    i__2 = j - 1;
#line 326 "cher2k.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 327 "cher2k.f"
			i__3 = i__ + j * c_dim1;
#line 327 "cher2k.f"
			i__4 = i__ + j * c_dim1;
#line 327 "cher2k.f"
			z__1.r = *beta * c__[i__4].r, z__1.i = *beta * c__[
				i__4].i;
#line 327 "cher2k.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 328 "cher2k.f"
/* L100: */
#line 328 "cher2k.f"
		    }
#line 329 "cher2k.f"
		    i__2 = j + j * c_dim1;
#line 329 "cher2k.f"
		    i__3 = j + j * c_dim1;
#line 329 "cher2k.f"
		    d__1 = *beta * c__[i__3].r;
#line 329 "cher2k.f"
		    c__[i__2].r = d__1, c__[i__2].i = 0.;
#line 330 "cher2k.f"
		} else {
#line 331 "cher2k.f"
		    i__2 = j + j * c_dim1;
#line 331 "cher2k.f"
		    i__3 = j + j * c_dim1;
#line 331 "cher2k.f"
		    d__1 = c__[i__3].r;
#line 331 "cher2k.f"
		    c__[i__2].r = d__1, c__[i__2].i = 0.;
#line 332 "cher2k.f"
		}
#line 333 "cher2k.f"
		i__2 = *k;
#line 333 "cher2k.f"
		for (l = 1; l <= i__2; ++l) {
#line 334 "cher2k.f"
		    i__3 = j + l * a_dim1;
#line 334 "cher2k.f"
		    i__4 = j + l * b_dim1;
#line 334 "cher2k.f"
		    if (a[i__3].r != 0. || a[i__3].i != 0. || (b[i__4].r != 
			    0. || b[i__4].i != 0.)) {
#line 335 "cher2k.f"
			d_cnjg(&z__2, &b[j + l * b_dim1]);
#line 335 "cher2k.f"
			z__1.r = alpha->r * z__2.r - alpha->i * z__2.i, 
				z__1.i = alpha->r * z__2.i + alpha->i * 
				z__2.r;
#line 335 "cher2k.f"
			temp1.r = z__1.r, temp1.i = z__1.i;
#line 336 "cher2k.f"
			i__3 = j + l * a_dim1;
#line 336 "cher2k.f"
			z__2.r = alpha->r * a[i__3].r - alpha->i * a[i__3].i, 
				z__2.i = alpha->r * a[i__3].i + alpha->i * a[
				i__3].r;
#line 336 "cher2k.f"
			d_cnjg(&z__1, &z__2);
#line 336 "cher2k.f"
			temp2.r = z__1.r, temp2.i = z__1.i;
#line 337 "cher2k.f"
			i__3 = j - 1;
#line 337 "cher2k.f"
			for (i__ = 1; i__ <= i__3; ++i__) {
#line 338 "cher2k.f"
			    i__4 = i__ + j * c_dim1;
#line 338 "cher2k.f"
			    i__5 = i__ + j * c_dim1;
#line 338 "cher2k.f"
			    i__6 = i__ + l * a_dim1;
#line 338 "cher2k.f"
			    z__3.r = a[i__6].r * temp1.r - a[i__6].i * 
				    temp1.i, z__3.i = a[i__6].r * temp1.i + a[
				    i__6].i * temp1.r;
#line 338 "cher2k.f"
			    z__2.r = c__[i__5].r + z__3.r, z__2.i = c__[i__5]
				    .i + z__3.i;
#line 338 "cher2k.f"
			    i__7 = i__ + l * b_dim1;
#line 338 "cher2k.f"
			    z__4.r = b[i__7].r * temp2.r - b[i__7].i * 
				    temp2.i, z__4.i = b[i__7].r * temp2.i + b[
				    i__7].i * temp2.r;
#line 338 "cher2k.f"
			    z__1.r = z__2.r + z__4.r, z__1.i = z__2.i + 
				    z__4.i;
#line 338 "cher2k.f"
			    c__[i__4].r = z__1.r, c__[i__4].i = z__1.i;
#line 340 "cher2k.f"
/* L110: */
#line 340 "cher2k.f"
			}
#line 341 "cher2k.f"
			i__3 = j + j * c_dim1;
#line 341 "cher2k.f"
			i__4 = j + j * c_dim1;
#line 341 "cher2k.f"
			i__5 = j + l * a_dim1;
#line 341 "cher2k.f"
			z__2.r = a[i__5].r * temp1.r - a[i__5].i * temp1.i, 
				z__2.i = a[i__5].r * temp1.i + a[i__5].i * 
				temp1.r;
#line 341 "cher2k.f"
			i__6 = j + l * b_dim1;
#line 341 "cher2k.f"
			z__3.r = b[i__6].r * temp2.r - b[i__6].i * temp2.i, 
				z__3.i = b[i__6].r * temp2.i + b[i__6].i * 
				temp2.r;
#line 341 "cher2k.f"
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 341 "cher2k.f"
			d__1 = c__[i__4].r + z__1.r;
#line 341 "cher2k.f"
			c__[i__3].r = d__1, c__[i__3].i = 0.;
#line 343 "cher2k.f"
		    }
#line 344 "cher2k.f"
/* L120: */
#line 344 "cher2k.f"
		}
#line 345 "cher2k.f"
/* L130: */
#line 345 "cher2k.f"
	    }
#line 346 "cher2k.f"
	} else {
#line 347 "cher2k.f"
	    i__1 = *n;
#line 347 "cher2k.f"
	    for (j = 1; j <= i__1; ++j) {
#line 348 "cher2k.f"
		if (*beta == 0.) {
#line 349 "cher2k.f"
		    i__2 = *n;
#line 349 "cher2k.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 350 "cher2k.f"
			i__3 = i__ + j * c_dim1;
#line 350 "cher2k.f"
			c__[i__3].r = 0., c__[i__3].i = 0.;
#line 351 "cher2k.f"
/* L140: */
#line 351 "cher2k.f"
		    }
#line 352 "cher2k.f"
		} else if (*beta != 1.) {
#line 353 "cher2k.f"
		    i__2 = *n;
#line 353 "cher2k.f"
		    for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 354 "cher2k.f"
			i__3 = i__ + j * c_dim1;
#line 354 "cher2k.f"
			i__4 = i__ + j * c_dim1;
#line 354 "cher2k.f"
			z__1.r = *beta * c__[i__4].r, z__1.i = *beta * c__[
				i__4].i;
#line 354 "cher2k.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 355 "cher2k.f"
/* L150: */
#line 355 "cher2k.f"
		    }
#line 356 "cher2k.f"
		    i__2 = j + j * c_dim1;
#line 356 "cher2k.f"
		    i__3 = j + j * c_dim1;
#line 356 "cher2k.f"
		    d__1 = *beta * c__[i__3].r;
#line 356 "cher2k.f"
		    c__[i__2].r = d__1, c__[i__2].i = 0.;
#line 357 "cher2k.f"
		} else {
#line 358 "cher2k.f"
		    i__2 = j + j * c_dim1;
#line 358 "cher2k.f"
		    i__3 = j + j * c_dim1;
#line 358 "cher2k.f"
		    d__1 = c__[i__3].r;
#line 358 "cher2k.f"
		    c__[i__2].r = d__1, c__[i__2].i = 0.;
#line 359 "cher2k.f"
		}
#line 360 "cher2k.f"
		i__2 = *k;
#line 360 "cher2k.f"
		for (l = 1; l <= i__2; ++l) {
#line 361 "cher2k.f"
		    i__3 = j + l * a_dim1;
#line 361 "cher2k.f"
		    i__4 = j + l * b_dim1;
#line 361 "cher2k.f"
		    if (a[i__3].r != 0. || a[i__3].i != 0. || (b[i__4].r != 
			    0. || b[i__4].i != 0.)) {
#line 362 "cher2k.f"
			d_cnjg(&z__2, &b[j + l * b_dim1]);
#line 362 "cher2k.f"
			z__1.r = alpha->r * z__2.r - alpha->i * z__2.i, 
				z__1.i = alpha->r * z__2.i + alpha->i * 
				z__2.r;
#line 362 "cher2k.f"
			temp1.r = z__1.r, temp1.i = z__1.i;
#line 363 "cher2k.f"
			i__3 = j + l * a_dim1;
#line 363 "cher2k.f"
			z__2.r = alpha->r * a[i__3].r - alpha->i * a[i__3].i, 
				z__2.i = alpha->r * a[i__3].i + alpha->i * a[
				i__3].r;
#line 363 "cher2k.f"
			d_cnjg(&z__1, &z__2);
#line 363 "cher2k.f"
			temp2.r = z__1.r, temp2.i = z__1.i;
#line 364 "cher2k.f"
			i__3 = *n;
#line 364 "cher2k.f"
			for (i__ = j + 1; i__ <= i__3; ++i__) {
#line 365 "cher2k.f"
			    i__4 = i__ + j * c_dim1;
#line 365 "cher2k.f"
			    i__5 = i__ + j * c_dim1;
#line 365 "cher2k.f"
			    i__6 = i__ + l * a_dim1;
#line 365 "cher2k.f"
			    z__3.r = a[i__6].r * temp1.r - a[i__6].i * 
				    temp1.i, z__3.i = a[i__6].r * temp1.i + a[
				    i__6].i * temp1.r;
#line 365 "cher2k.f"
			    z__2.r = c__[i__5].r + z__3.r, z__2.i = c__[i__5]
				    .i + z__3.i;
#line 365 "cher2k.f"
			    i__7 = i__ + l * b_dim1;
#line 365 "cher2k.f"
			    z__4.r = b[i__7].r * temp2.r - b[i__7].i * 
				    temp2.i, z__4.i = b[i__7].r * temp2.i + b[
				    i__7].i * temp2.r;
#line 365 "cher2k.f"
			    z__1.r = z__2.r + z__4.r, z__1.i = z__2.i + 
				    z__4.i;
#line 365 "cher2k.f"
			    c__[i__4].r = z__1.r, c__[i__4].i = z__1.i;
#line 367 "cher2k.f"
/* L160: */
#line 367 "cher2k.f"
			}
#line 368 "cher2k.f"
			i__3 = j + j * c_dim1;
#line 368 "cher2k.f"
			i__4 = j + j * c_dim1;
#line 368 "cher2k.f"
			i__5 = j + l * a_dim1;
#line 368 "cher2k.f"
			z__2.r = a[i__5].r * temp1.r - a[i__5].i * temp1.i, 
				z__2.i = a[i__5].r * temp1.i + a[i__5].i * 
				temp1.r;
#line 368 "cher2k.f"
			i__6 = j + l * b_dim1;
#line 368 "cher2k.f"
			z__3.r = b[i__6].r * temp2.r - b[i__6].i * temp2.i, 
				z__3.i = b[i__6].r * temp2.i + b[i__6].i * 
				temp2.r;
#line 368 "cher2k.f"
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 368 "cher2k.f"
			d__1 = c__[i__4].r + z__1.r;
#line 368 "cher2k.f"
			c__[i__3].r = d__1, c__[i__3].i = 0.;
#line 370 "cher2k.f"
		    }
#line 371 "cher2k.f"
/* L170: */
#line 371 "cher2k.f"
		}
#line 372 "cher2k.f"
/* L180: */
#line 372 "cher2k.f"
	    }
#line 373 "cher2k.f"
	}
#line 374 "cher2k.f"
    } else {

/*        Form  C := alpha*A**H*B + conjg( alpha )*B**H*A + */
/*                   C. */

#line 379 "cher2k.f"
	if (upper) {
#line 380 "cher2k.f"
	    i__1 = *n;
#line 380 "cher2k.f"
	    for (j = 1; j <= i__1; ++j) {
#line 381 "cher2k.f"
		i__2 = j;
#line 381 "cher2k.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 382 "cher2k.f"
		    temp1.r = 0., temp1.i = 0.;
#line 383 "cher2k.f"
		    temp2.r = 0., temp2.i = 0.;
#line 384 "cher2k.f"
		    i__3 = *k;
#line 384 "cher2k.f"
		    for (l = 1; l <= i__3; ++l) {
#line 385 "cher2k.f"
			d_cnjg(&z__3, &a[l + i__ * a_dim1]);
#line 385 "cher2k.f"
			i__4 = l + j * b_dim1;
#line 385 "cher2k.f"
			z__2.r = z__3.r * b[i__4].r - z__3.i * b[i__4].i, 
				z__2.i = z__3.r * b[i__4].i + z__3.i * b[i__4]
				.r;
#line 385 "cher2k.f"
			z__1.r = temp1.r + z__2.r, z__1.i = temp1.i + z__2.i;
#line 385 "cher2k.f"
			temp1.r = z__1.r, temp1.i = z__1.i;
#line 386 "cher2k.f"
			d_cnjg(&z__3, &b[l + i__ * b_dim1]);
#line 386 "cher2k.f"
			i__4 = l + j * a_dim1;
#line 386 "cher2k.f"
			z__2.r = z__3.r * a[i__4].r - z__3.i * a[i__4].i, 
				z__2.i = z__3.r * a[i__4].i + z__3.i * a[i__4]
				.r;
#line 386 "cher2k.f"
			z__1.r = temp2.r + z__2.r, z__1.i = temp2.i + z__2.i;
#line 386 "cher2k.f"
			temp2.r = z__1.r, temp2.i = z__1.i;
#line 387 "cher2k.f"
/* L190: */
#line 387 "cher2k.f"
		    }
#line 388 "cher2k.f"
		    if (i__ == j) {
#line 389 "cher2k.f"
			if (*beta == 0.) {
#line 390 "cher2k.f"
			    i__3 = j + j * c_dim1;
#line 390 "cher2k.f"
			    z__2.r = alpha->r * temp1.r - alpha->i * temp1.i, 
				    z__2.i = alpha->r * temp1.i + alpha->i * 
				    temp1.r;
#line 390 "cher2k.f"
			    d_cnjg(&z__4, alpha);
#line 390 "cher2k.f"
			    z__3.r = z__4.r * temp2.r - z__4.i * temp2.i, 
				    z__3.i = z__4.r * temp2.i + z__4.i * 
				    temp2.r;
#line 390 "cher2k.f"
			    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + 
				    z__3.i;
#line 390 "cher2k.f"
			    d__1 = z__1.r;
#line 390 "cher2k.f"
			    c__[i__3].r = d__1, c__[i__3].i = 0.;
#line 392 "cher2k.f"
			} else {
#line 393 "cher2k.f"
			    i__3 = j + j * c_dim1;
#line 393 "cher2k.f"
			    i__4 = j + j * c_dim1;
#line 393 "cher2k.f"
			    z__2.r = alpha->r * temp1.r - alpha->i * temp1.i, 
				    z__2.i = alpha->r * temp1.i + alpha->i * 
				    temp1.r;
#line 393 "cher2k.f"
			    d_cnjg(&z__4, alpha);
#line 393 "cher2k.f"
			    z__3.r = z__4.r * temp2.r - z__4.i * temp2.i, 
				    z__3.i = z__4.r * temp2.i + z__4.i * 
				    temp2.r;
#line 393 "cher2k.f"
			    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + 
				    z__3.i;
#line 393 "cher2k.f"
			    d__1 = *beta * c__[i__4].r + z__1.r;
#line 393 "cher2k.f"
			    c__[i__3].r = d__1, c__[i__3].i = 0.;
#line 396 "cher2k.f"
			}
#line 397 "cher2k.f"
		    } else {
#line 398 "cher2k.f"
			if (*beta == 0.) {
#line 399 "cher2k.f"
			    i__3 = i__ + j * c_dim1;
#line 399 "cher2k.f"
			    z__2.r = alpha->r * temp1.r - alpha->i * temp1.i, 
				    z__2.i = alpha->r * temp1.i + alpha->i * 
				    temp1.r;
#line 399 "cher2k.f"
			    d_cnjg(&z__4, alpha);
#line 399 "cher2k.f"
			    z__3.r = z__4.r * temp2.r - z__4.i * temp2.i, 
				    z__3.i = z__4.r * temp2.i + z__4.i * 
				    temp2.r;
#line 399 "cher2k.f"
			    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + 
				    z__3.i;
#line 399 "cher2k.f"
			    c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 400 "cher2k.f"
			} else {
#line 401 "cher2k.f"
			    i__3 = i__ + j * c_dim1;
#line 401 "cher2k.f"
			    i__4 = i__ + j * c_dim1;
#line 401 "cher2k.f"
			    z__3.r = *beta * c__[i__4].r, z__3.i = *beta * 
				    c__[i__4].i;
#line 401 "cher2k.f"
			    z__4.r = alpha->r * temp1.r - alpha->i * temp1.i, 
				    z__4.i = alpha->r * temp1.i + alpha->i * 
				    temp1.r;
#line 401 "cher2k.f"
			    z__2.r = z__3.r + z__4.r, z__2.i = z__3.i + 
				    z__4.i;
#line 401 "cher2k.f"
			    d_cnjg(&z__6, alpha);
#line 401 "cher2k.f"
			    z__5.r = z__6.r * temp2.r - z__6.i * temp2.i, 
				    z__5.i = z__6.r * temp2.i + z__6.i * 
				    temp2.r;
#line 401 "cher2k.f"
			    z__1.r = z__2.r + z__5.r, z__1.i = z__2.i + 
				    z__5.i;
#line 401 "cher2k.f"
			    c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 403 "cher2k.f"
			}
#line 404 "cher2k.f"
		    }
#line 405 "cher2k.f"
/* L200: */
#line 405 "cher2k.f"
		}
#line 406 "cher2k.f"
/* L210: */
#line 406 "cher2k.f"
	    }
#line 407 "cher2k.f"
	} else {
#line 408 "cher2k.f"
	    i__1 = *n;
#line 408 "cher2k.f"
	    for (j = 1; j <= i__1; ++j) {
#line 409 "cher2k.f"
		i__2 = *n;
#line 409 "cher2k.f"
		for (i__ = j; i__ <= i__2; ++i__) {
#line 410 "cher2k.f"
		    temp1.r = 0., temp1.i = 0.;
#line 411 "cher2k.f"
		    temp2.r = 0., temp2.i = 0.;
#line 412 "cher2k.f"
		    i__3 = *k;
#line 412 "cher2k.f"
		    for (l = 1; l <= i__3; ++l) {
#line 413 "cher2k.f"
			d_cnjg(&z__3, &a[l + i__ * a_dim1]);
#line 413 "cher2k.f"
			i__4 = l + j * b_dim1;
#line 413 "cher2k.f"
			z__2.r = z__3.r * b[i__4].r - z__3.i * b[i__4].i, 
				z__2.i = z__3.r * b[i__4].i + z__3.i * b[i__4]
				.r;
#line 413 "cher2k.f"
			z__1.r = temp1.r + z__2.r, z__1.i = temp1.i + z__2.i;
#line 413 "cher2k.f"
			temp1.r = z__1.r, temp1.i = z__1.i;
#line 414 "cher2k.f"
			d_cnjg(&z__3, &b[l + i__ * b_dim1]);
#line 414 "cher2k.f"
			i__4 = l + j * a_dim1;
#line 414 "cher2k.f"
			z__2.r = z__3.r * a[i__4].r - z__3.i * a[i__4].i, 
				z__2.i = z__3.r * a[i__4].i + z__3.i * a[i__4]
				.r;
#line 414 "cher2k.f"
			z__1.r = temp2.r + z__2.r, z__1.i = temp2.i + z__2.i;
#line 414 "cher2k.f"
			temp2.r = z__1.r, temp2.i = z__1.i;
#line 415 "cher2k.f"
/* L220: */
#line 415 "cher2k.f"
		    }
#line 416 "cher2k.f"
		    if (i__ == j) {
#line 417 "cher2k.f"
			if (*beta == 0.) {
#line 418 "cher2k.f"
			    i__3 = j + j * c_dim1;
#line 418 "cher2k.f"
			    z__2.r = alpha->r * temp1.r - alpha->i * temp1.i, 
				    z__2.i = alpha->r * temp1.i + alpha->i * 
				    temp1.r;
#line 418 "cher2k.f"
			    d_cnjg(&z__4, alpha);
#line 418 "cher2k.f"
			    z__3.r = z__4.r * temp2.r - z__4.i * temp2.i, 
				    z__3.i = z__4.r * temp2.i + z__4.i * 
				    temp2.r;
#line 418 "cher2k.f"
			    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + 
				    z__3.i;
#line 418 "cher2k.f"
			    d__1 = z__1.r;
#line 418 "cher2k.f"
			    c__[i__3].r = d__1, c__[i__3].i = 0.;
#line 420 "cher2k.f"
			} else {
#line 421 "cher2k.f"
			    i__3 = j + j * c_dim1;
#line 421 "cher2k.f"
			    i__4 = j + j * c_dim1;
#line 421 "cher2k.f"
			    z__2.r = alpha->r * temp1.r - alpha->i * temp1.i, 
				    z__2.i = alpha->r * temp1.i + alpha->i * 
				    temp1.r;
#line 421 "cher2k.f"
			    d_cnjg(&z__4, alpha);
#line 421 "cher2k.f"
			    z__3.r = z__4.r * temp2.r - z__4.i * temp2.i, 
				    z__3.i = z__4.r * temp2.i + z__4.i * 
				    temp2.r;
#line 421 "cher2k.f"
			    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + 
				    z__3.i;
#line 421 "cher2k.f"
			    d__1 = *beta * c__[i__4].r + z__1.r;
#line 421 "cher2k.f"
			    c__[i__3].r = d__1, c__[i__3].i = 0.;
#line 424 "cher2k.f"
			}
#line 425 "cher2k.f"
		    } else {
#line 426 "cher2k.f"
			if (*beta == 0.) {
#line 427 "cher2k.f"
			    i__3 = i__ + j * c_dim1;
#line 427 "cher2k.f"
			    z__2.r = alpha->r * temp1.r - alpha->i * temp1.i, 
				    z__2.i = alpha->r * temp1.i + alpha->i * 
				    temp1.r;
#line 427 "cher2k.f"
			    d_cnjg(&z__4, alpha);
#line 427 "cher2k.f"
			    z__3.r = z__4.r * temp2.r - z__4.i * temp2.i, 
				    z__3.i = z__4.r * temp2.i + z__4.i * 
				    temp2.r;
#line 427 "cher2k.f"
			    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + 
				    z__3.i;
#line 427 "cher2k.f"
			    c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 428 "cher2k.f"
			} else {
#line 429 "cher2k.f"
			    i__3 = i__ + j * c_dim1;
#line 429 "cher2k.f"
			    i__4 = i__ + j * c_dim1;
#line 429 "cher2k.f"
			    z__3.r = *beta * c__[i__4].r, z__3.i = *beta * 
				    c__[i__4].i;
#line 429 "cher2k.f"
			    z__4.r = alpha->r * temp1.r - alpha->i * temp1.i, 
				    z__4.i = alpha->r * temp1.i + alpha->i * 
				    temp1.r;
#line 429 "cher2k.f"
			    z__2.r = z__3.r + z__4.r, z__2.i = z__3.i + 
				    z__4.i;
#line 429 "cher2k.f"
			    d_cnjg(&z__6, alpha);
#line 429 "cher2k.f"
			    z__5.r = z__6.r * temp2.r - z__6.i * temp2.i, 
				    z__5.i = z__6.r * temp2.i + z__6.i * 
				    temp2.r;
#line 429 "cher2k.f"
			    z__1.r = z__2.r + z__5.r, z__1.i = z__2.i + 
				    z__5.i;
#line 429 "cher2k.f"
			    c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 431 "cher2k.f"
			}
#line 432 "cher2k.f"
		    }
#line 433 "cher2k.f"
/* L230: */
#line 433 "cher2k.f"
		}
#line 434 "cher2k.f"
/* L240: */
#line 434 "cher2k.f"
	    }
#line 435 "cher2k.f"
	}
#line 436 "cher2k.f"
    }

#line 438 "cher2k.f"
    return 0;

/*     End of CHER2K. */

} /* cher2k_ */

