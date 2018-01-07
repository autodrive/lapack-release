#line 1 "cherk.f"
/* cherk.f -- translated by f2c (version 20100827).
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

#line 1 "cherk.f"
/* > \brief \b CHERK */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CHERK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC) */

/*       .. Scalar Arguments .. */
/*       REAL ALPHA,BETA */
/*       INTEGER K,LDA,LDC,N */
/*       CHARACTER TRANS,UPLO */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX A(LDA,*),C(LDC,*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CHERK  performs one of the hermitian rank k operations */
/* > */
/* >    C := alpha*A*A**H + beta*C, */
/* > */
/* > or */
/* > */
/* >    C := alpha*A**H*A + beta*C, */
/* > */
/* > where  alpha and beta  are  real scalars,  C is an  n by n  hermitian */
/* > matrix and  A  is an  n by k  matrix in the  first case and a  k by n */
/* > matrix in the second case. */
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
/* >              TRANS = 'N' or 'n'   C := alpha*A*A**H + beta*C. */
/* > */
/* >              TRANS = 'C' or 'c'   C := alpha*A**H*A + beta*C. */
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
/* >           TRANS = 'C' or 'c',  K  specifies  the number of rows of the */
/* >           matrix A.  K must be at least zero. */
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
/* Subroutine */ int cherk_(char *uplo, char *trans, integer *n, integer *k, 
	doublereal *alpha, doublecomplex *a, integer *lda, doublereal *beta, 
	doublecomplex *c__, integer *ldc, ftnlen uplo_len, ftnlen trans_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, c_dim1, c_offset, i__1, i__2, i__3, i__4, i__5, 
	    i__6;
    doublereal d__1;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, j, l, info;
    static doublecomplex temp;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer nrowa;
    static doublereal rtemp;
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

#line 215 "cherk.f"
    /* Parameter adjustments */
#line 215 "cherk.f"
    a_dim1 = *lda;
#line 215 "cherk.f"
    a_offset = 1 + a_dim1;
#line 215 "cherk.f"
    a -= a_offset;
#line 215 "cherk.f"
    c_dim1 = *ldc;
#line 215 "cherk.f"
    c_offset = 1 + c_dim1;
#line 215 "cherk.f"
    c__ -= c_offset;
#line 215 "cherk.f"

#line 215 "cherk.f"
    /* Function Body */
#line 215 "cherk.f"
    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {
#line 216 "cherk.f"
	nrowa = *n;
#line 217 "cherk.f"
    } else {
#line 218 "cherk.f"
	nrowa = *k;
#line 219 "cherk.f"
    }
#line 220 "cherk.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);

#line 222 "cherk.f"
    info = 0;
#line 223 "cherk.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 224 "cherk.f"
	info = 1;
#line 225 "cherk.f"
    } else if (! lsame_(trans, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, 
	    "C", (ftnlen)1, (ftnlen)1)) {
#line 227 "cherk.f"
	info = 2;
#line 228 "cherk.f"
    } else if (*n < 0) {
#line 229 "cherk.f"
	info = 3;
#line 230 "cherk.f"
    } else if (*k < 0) {
#line 231 "cherk.f"
	info = 4;
#line 232 "cherk.f"
    } else if (*lda < max(1,nrowa)) {
#line 233 "cherk.f"
	info = 7;
#line 234 "cherk.f"
    } else if (*ldc < max(1,*n)) {
#line 235 "cherk.f"
	info = 10;
#line 236 "cherk.f"
    }
#line 237 "cherk.f"
    if (info != 0) {
#line 238 "cherk.f"
	xerbla_("CHERK ", &info, (ftnlen)6);
#line 239 "cherk.f"
	return 0;
#line 240 "cherk.f"
    }

/*     Quick return if possible. */

#line 244 "cherk.f"
    if (*n == 0 || (*alpha == 0. || *k == 0) && *beta == 1.) {
#line 244 "cherk.f"
	return 0;
#line 244 "cherk.f"
    }

/*     And when  alpha.eq.zero. */

#line 249 "cherk.f"
    if (*alpha == 0.) {
#line 250 "cherk.f"
	if (upper) {
#line 251 "cherk.f"
	    if (*beta == 0.) {
#line 252 "cherk.f"
		i__1 = *n;
#line 252 "cherk.f"
		for (j = 1; j <= i__1; ++j) {
#line 253 "cherk.f"
		    i__2 = j;
#line 253 "cherk.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 254 "cherk.f"
			i__3 = i__ + j * c_dim1;
#line 254 "cherk.f"
			c__[i__3].r = 0., c__[i__3].i = 0.;
#line 255 "cherk.f"
/* L10: */
#line 255 "cherk.f"
		    }
#line 256 "cherk.f"
/* L20: */
#line 256 "cherk.f"
		}
#line 257 "cherk.f"
	    } else {
#line 258 "cherk.f"
		i__1 = *n;
#line 258 "cherk.f"
		for (j = 1; j <= i__1; ++j) {
#line 259 "cherk.f"
		    i__2 = j - 1;
#line 259 "cherk.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 260 "cherk.f"
			i__3 = i__ + j * c_dim1;
#line 260 "cherk.f"
			i__4 = i__ + j * c_dim1;
#line 260 "cherk.f"
			z__1.r = *beta * c__[i__4].r, z__1.i = *beta * c__[
				i__4].i;
#line 260 "cherk.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 261 "cherk.f"
/* L30: */
#line 261 "cherk.f"
		    }
#line 262 "cherk.f"
		    i__2 = j + j * c_dim1;
#line 262 "cherk.f"
		    i__3 = j + j * c_dim1;
#line 262 "cherk.f"
		    d__1 = *beta * c__[i__3].r;
#line 262 "cherk.f"
		    c__[i__2].r = d__1, c__[i__2].i = 0.;
#line 263 "cherk.f"
/* L40: */
#line 263 "cherk.f"
		}
#line 264 "cherk.f"
	    }
#line 265 "cherk.f"
	} else {
#line 266 "cherk.f"
	    if (*beta == 0.) {
#line 267 "cherk.f"
		i__1 = *n;
#line 267 "cherk.f"
		for (j = 1; j <= i__1; ++j) {
#line 268 "cherk.f"
		    i__2 = *n;
#line 268 "cherk.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 269 "cherk.f"
			i__3 = i__ + j * c_dim1;
#line 269 "cherk.f"
			c__[i__3].r = 0., c__[i__3].i = 0.;
#line 270 "cherk.f"
/* L50: */
#line 270 "cherk.f"
		    }
#line 271 "cherk.f"
/* L60: */
#line 271 "cherk.f"
		}
#line 272 "cherk.f"
	    } else {
#line 273 "cherk.f"
		i__1 = *n;
#line 273 "cherk.f"
		for (j = 1; j <= i__1; ++j) {
#line 274 "cherk.f"
		    i__2 = j + j * c_dim1;
#line 274 "cherk.f"
		    i__3 = j + j * c_dim1;
#line 274 "cherk.f"
		    d__1 = *beta * c__[i__3].r;
#line 274 "cherk.f"
		    c__[i__2].r = d__1, c__[i__2].i = 0.;
#line 275 "cherk.f"
		    i__2 = *n;
#line 275 "cherk.f"
		    for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 276 "cherk.f"
			i__3 = i__ + j * c_dim1;
#line 276 "cherk.f"
			i__4 = i__ + j * c_dim1;
#line 276 "cherk.f"
			z__1.r = *beta * c__[i__4].r, z__1.i = *beta * c__[
				i__4].i;
#line 276 "cherk.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 277 "cherk.f"
/* L70: */
#line 277 "cherk.f"
		    }
#line 278 "cherk.f"
/* L80: */
#line 278 "cherk.f"
		}
#line 279 "cherk.f"
	    }
#line 280 "cherk.f"
	}
#line 281 "cherk.f"
	return 0;
#line 282 "cherk.f"
    }

/*     Start the operations. */

#line 286 "cherk.f"
    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {

/*        Form  C := alpha*A*A**H + beta*C. */

#line 290 "cherk.f"
	if (upper) {
#line 291 "cherk.f"
	    i__1 = *n;
#line 291 "cherk.f"
	    for (j = 1; j <= i__1; ++j) {
#line 292 "cherk.f"
		if (*beta == 0.) {
#line 293 "cherk.f"
		    i__2 = j;
#line 293 "cherk.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 294 "cherk.f"
			i__3 = i__ + j * c_dim1;
#line 294 "cherk.f"
			c__[i__3].r = 0., c__[i__3].i = 0.;
#line 295 "cherk.f"
/* L90: */
#line 295 "cherk.f"
		    }
#line 296 "cherk.f"
		} else if (*beta != 1.) {
#line 297 "cherk.f"
		    i__2 = j - 1;
#line 297 "cherk.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 298 "cherk.f"
			i__3 = i__ + j * c_dim1;
#line 298 "cherk.f"
			i__4 = i__ + j * c_dim1;
#line 298 "cherk.f"
			z__1.r = *beta * c__[i__4].r, z__1.i = *beta * c__[
				i__4].i;
#line 298 "cherk.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 299 "cherk.f"
/* L100: */
#line 299 "cherk.f"
		    }
#line 300 "cherk.f"
		    i__2 = j + j * c_dim1;
#line 300 "cherk.f"
		    i__3 = j + j * c_dim1;
#line 300 "cherk.f"
		    d__1 = *beta * c__[i__3].r;
#line 300 "cherk.f"
		    c__[i__2].r = d__1, c__[i__2].i = 0.;
#line 301 "cherk.f"
		} else {
#line 302 "cherk.f"
		    i__2 = j + j * c_dim1;
#line 302 "cherk.f"
		    i__3 = j + j * c_dim1;
#line 302 "cherk.f"
		    d__1 = c__[i__3].r;
#line 302 "cherk.f"
		    c__[i__2].r = d__1, c__[i__2].i = 0.;
#line 303 "cherk.f"
		}
#line 304 "cherk.f"
		i__2 = *k;
#line 304 "cherk.f"
		for (l = 1; l <= i__2; ++l) {
#line 305 "cherk.f"
		    i__3 = j + l * a_dim1;
#line 305 "cherk.f"
		    if (a[i__3].r != 0. || a[i__3].i != 0.) {
#line 306 "cherk.f"
			d_cnjg(&z__2, &a[j + l * a_dim1]);
#line 306 "cherk.f"
			z__1.r = *alpha * z__2.r, z__1.i = *alpha * z__2.i;
#line 306 "cherk.f"
			temp.r = z__1.r, temp.i = z__1.i;
#line 307 "cherk.f"
			i__3 = j - 1;
#line 307 "cherk.f"
			for (i__ = 1; i__ <= i__3; ++i__) {
#line 308 "cherk.f"
			    i__4 = i__ + j * c_dim1;
#line 308 "cherk.f"
			    i__5 = i__ + j * c_dim1;
#line 308 "cherk.f"
			    i__6 = i__ + l * a_dim1;
#line 308 "cherk.f"
			    z__2.r = temp.r * a[i__6].r - temp.i * a[i__6].i, 
				    z__2.i = temp.r * a[i__6].i + temp.i * a[
				    i__6].r;
#line 308 "cherk.f"
			    z__1.r = c__[i__5].r + z__2.r, z__1.i = c__[i__5]
				    .i + z__2.i;
#line 308 "cherk.f"
			    c__[i__4].r = z__1.r, c__[i__4].i = z__1.i;
#line 309 "cherk.f"
/* L110: */
#line 309 "cherk.f"
			}
#line 310 "cherk.f"
			i__3 = j + j * c_dim1;
#line 310 "cherk.f"
			i__4 = j + j * c_dim1;
#line 310 "cherk.f"
			i__5 = i__ + l * a_dim1;
#line 310 "cherk.f"
			z__1.r = temp.r * a[i__5].r - temp.i * a[i__5].i, 
				z__1.i = temp.r * a[i__5].i + temp.i * a[i__5]
				.r;
#line 310 "cherk.f"
			d__1 = c__[i__4].r + z__1.r;
#line 310 "cherk.f"
			c__[i__3].r = d__1, c__[i__3].i = 0.;
#line 311 "cherk.f"
		    }
#line 312 "cherk.f"
/* L120: */
#line 312 "cherk.f"
		}
#line 313 "cherk.f"
/* L130: */
#line 313 "cherk.f"
	    }
#line 314 "cherk.f"
	} else {
#line 315 "cherk.f"
	    i__1 = *n;
#line 315 "cherk.f"
	    for (j = 1; j <= i__1; ++j) {
#line 316 "cherk.f"
		if (*beta == 0.) {
#line 317 "cherk.f"
		    i__2 = *n;
#line 317 "cherk.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 318 "cherk.f"
			i__3 = i__ + j * c_dim1;
#line 318 "cherk.f"
			c__[i__3].r = 0., c__[i__3].i = 0.;
#line 319 "cherk.f"
/* L140: */
#line 319 "cherk.f"
		    }
#line 320 "cherk.f"
		} else if (*beta != 1.) {
#line 321 "cherk.f"
		    i__2 = j + j * c_dim1;
#line 321 "cherk.f"
		    i__3 = j + j * c_dim1;
#line 321 "cherk.f"
		    d__1 = *beta * c__[i__3].r;
#line 321 "cherk.f"
		    c__[i__2].r = d__1, c__[i__2].i = 0.;
#line 322 "cherk.f"
		    i__2 = *n;
#line 322 "cherk.f"
		    for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 323 "cherk.f"
			i__3 = i__ + j * c_dim1;
#line 323 "cherk.f"
			i__4 = i__ + j * c_dim1;
#line 323 "cherk.f"
			z__1.r = *beta * c__[i__4].r, z__1.i = *beta * c__[
				i__4].i;
#line 323 "cherk.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 324 "cherk.f"
/* L150: */
#line 324 "cherk.f"
		    }
#line 325 "cherk.f"
		} else {
#line 326 "cherk.f"
		    i__2 = j + j * c_dim1;
#line 326 "cherk.f"
		    i__3 = j + j * c_dim1;
#line 326 "cherk.f"
		    d__1 = c__[i__3].r;
#line 326 "cherk.f"
		    c__[i__2].r = d__1, c__[i__2].i = 0.;
#line 327 "cherk.f"
		}
#line 328 "cherk.f"
		i__2 = *k;
#line 328 "cherk.f"
		for (l = 1; l <= i__2; ++l) {
#line 329 "cherk.f"
		    i__3 = j + l * a_dim1;
#line 329 "cherk.f"
		    if (a[i__3].r != 0. || a[i__3].i != 0.) {
#line 330 "cherk.f"
			d_cnjg(&z__2, &a[j + l * a_dim1]);
#line 330 "cherk.f"
			z__1.r = *alpha * z__2.r, z__1.i = *alpha * z__2.i;
#line 330 "cherk.f"
			temp.r = z__1.r, temp.i = z__1.i;
#line 331 "cherk.f"
			i__3 = j + j * c_dim1;
#line 331 "cherk.f"
			i__4 = j + j * c_dim1;
#line 331 "cherk.f"
			i__5 = j + l * a_dim1;
#line 331 "cherk.f"
			z__1.r = temp.r * a[i__5].r - temp.i * a[i__5].i, 
				z__1.i = temp.r * a[i__5].i + temp.i * a[i__5]
				.r;
#line 331 "cherk.f"
			d__1 = c__[i__4].r + z__1.r;
#line 331 "cherk.f"
			c__[i__3].r = d__1, c__[i__3].i = 0.;
#line 332 "cherk.f"
			i__3 = *n;
#line 332 "cherk.f"
			for (i__ = j + 1; i__ <= i__3; ++i__) {
#line 333 "cherk.f"
			    i__4 = i__ + j * c_dim1;
#line 333 "cherk.f"
			    i__5 = i__ + j * c_dim1;
#line 333 "cherk.f"
			    i__6 = i__ + l * a_dim1;
#line 333 "cherk.f"
			    z__2.r = temp.r * a[i__6].r - temp.i * a[i__6].i, 
				    z__2.i = temp.r * a[i__6].i + temp.i * a[
				    i__6].r;
#line 333 "cherk.f"
			    z__1.r = c__[i__5].r + z__2.r, z__1.i = c__[i__5]
				    .i + z__2.i;
#line 333 "cherk.f"
			    c__[i__4].r = z__1.r, c__[i__4].i = z__1.i;
#line 334 "cherk.f"
/* L160: */
#line 334 "cherk.f"
			}
#line 335 "cherk.f"
		    }
#line 336 "cherk.f"
/* L170: */
#line 336 "cherk.f"
		}
#line 337 "cherk.f"
/* L180: */
#line 337 "cherk.f"
	    }
#line 338 "cherk.f"
	}
#line 339 "cherk.f"
    } else {

/*        Form  C := alpha*A**H*A + beta*C. */

#line 343 "cherk.f"
	if (upper) {
#line 344 "cherk.f"
	    i__1 = *n;
#line 344 "cherk.f"
	    for (j = 1; j <= i__1; ++j) {
#line 345 "cherk.f"
		i__2 = j - 1;
#line 345 "cherk.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 346 "cherk.f"
		    temp.r = 0., temp.i = 0.;
#line 347 "cherk.f"
		    i__3 = *k;
#line 347 "cherk.f"
		    for (l = 1; l <= i__3; ++l) {
#line 348 "cherk.f"
			d_cnjg(&z__3, &a[l + i__ * a_dim1]);
#line 348 "cherk.f"
			i__4 = l + j * a_dim1;
#line 348 "cherk.f"
			z__2.r = z__3.r * a[i__4].r - z__3.i * a[i__4].i, 
				z__2.i = z__3.r * a[i__4].i + z__3.i * a[i__4]
				.r;
#line 348 "cherk.f"
			z__1.r = temp.r + z__2.r, z__1.i = temp.i + z__2.i;
#line 348 "cherk.f"
			temp.r = z__1.r, temp.i = z__1.i;
#line 349 "cherk.f"
/* L190: */
#line 349 "cherk.f"
		    }
#line 350 "cherk.f"
		    if (*beta == 0.) {
#line 351 "cherk.f"
			i__3 = i__ + j * c_dim1;
#line 351 "cherk.f"
			z__1.r = *alpha * temp.r, z__1.i = *alpha * temp.i;
#line 351 "cherk.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 352 "cherk.f"
		    } else {
#line 353 "cherk.f"
			i__3 = i__ + j * c_dim1;
#line 353 "cherk.f"
			z__2.r = *alpha * temp.r, z__2.i = *alpha * temp.i;
#line 353 "cherk.f"
			i__4 = i__ + j * c_dim1;
#line 353 "cherk.f"
			z__3.r = *beta * c__[i__4].r, z__3.i = *beta * c__[
				i__4].i;
#line 353 "cherk.f"
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 353 "cherk.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 354 "cherk.f"
		    }
#line 355 "cherk.f"
/* L200: */
#line 355 "cherk.f"
		}
#line 356 "cherk.f"
		rtemp = 0.;
#line 357 "cherk.f"
		i__2 = *k;
#line 357 "cherk.f"
		for (l = 1; l <= i__2; ++l) {
#line 358 "cherk.f"
		    d_cnjg(&z__3, &a[l + j * a_dim1]);
#line 358 "cherk.f"
		    i__3 = l + j * a_dim1;
#line 358 "cherk.f"
		    z__2.r = z__3.r * a[i__3].r - z__3.i * a[i__3].i, z__2.i =
			     z__3.r * a[i__3].i + z__3.i * a[i__3].r;
#line 358 "cherk.f"
		    z__1.r = rtemp + z__2.r, z__1.i = z__2.i;
#line 358 "cherk.f"
		    rtemp = z__1.r;
#line 359 "cherk.f"
/* L210: */
#line 359 "cherk.f"
		}
#line 360 "cherk.f"
		if (*beta == 0.) {
#line 361 "cherk.f"
		    i__2 = j + j * c_dim1;
#line 361 "cherk.f"
		    d__1 = *alpha * rtemp;
#line 361 "cherk.f"
		    c__[i__2].r = d__1, c__[i__2].i = 0.;
#line 362 "cherk.f"
		} else {
#line 363 "cherk.f"
		    i__2 = j + j * c_dim1;
#line 363 "cherk.f"
		    i__3 = j + j * c_dim1;
#line 363 "cherk.f"
		    d__1 = *alpha * rtemp + *beta * c__[i__3].r;
#line 363 "cherk.f"
		    c__[i__2].r = d__1, c__[i__2].i = 0.;
#line 364 "cherk.f"
		}
#line 365 "cherk.f"
/* L220: */
#line 365 "cherk.f"
	    }
#line 366 "cherk.f"
	} else {
#line 367 "cherk.f"
	    i__1 = *n;
#line 367 "cherk.f"
	    for (j = 1; j <= i__1; ++j) {
#line 368 "cherk.f"
		rtemp = 0.;
#line 369 "cherk.f"
		i__2 = *k;
#line 369 "cherk.f"
		for (l = 1; l <= i__2; ++l) {
#line 370 "cherk.f"
		    d_cnjg(&z__3, &a[l + j * a_dim1]);
#line 370 "cherk.f"
		    i__3 = l + j * a_dim1;
#line 370 "cherk.f"
		    z__2.r = z__3.r * a[i__3].r - z__3.i * a[i__3].i, z__2.i =
			     z__3.r * a[i__3].i + z__3.i * a[i__3].r;
#line 370 "cherk.f"
		    z__1.r = rtemp + z__2.r, z__1.i = z__2.i;
#line 370 "cherk.f"
		    rtemp = z__1.r;
#line 371 "cherk.f"
/* L230: */
#line 371 "cherk.f"
		}
#line 372 "cherk.f"
		if (*beta == 0.) {
#line 373 "cherk.f"
		    i__2 = j + j * c_dim1;
#line 373 "cherk.f"
		    d__1 = *alpha * rtemp;
#line 373 "cherk.f"
		    c__[i__2].r = d__1, c__[i__2].i = 0.;
#line 374 "cherk.f"
		} else {
#line 375 "cherk.f"
		    i__2 = j + j * c_dim1;
#line 375 "cherk.f"
		    i__3 = j + j * c_dim1;
#line 375 "cherk.f"
		    d__1 = *alpha * rtemp + *beta * c__[i__3].r;
#line 375 "cherk.f"
		    c__[i__2].r = d__1, c__[i__2].i = 0.;
#line 376 "cherk.f"
		}
#line 377 "cherk.f"
		i__2 = *n;
#line 377 "cherk.f"
		for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 378 "cherk.f"
		    temp.r = 0., temp.i = 0.;
#line 379 "cherk.f"
		    i__3 = *k;
#line 379 "cherk.f"
		    for (l = 1; l <= i__3; ++l) {
#line 380 "cherk.f"
			d_cnjg(&z__3, &a[l + i__ * a_dim1]);
#line 380 "cherk.f"
			i__4 = l + j * a_dim1;
#line 380 "cherk.f"
			z__2.r = z__3.r * a[i__4].r - z__3.i * a[i__4].i, 
				z__2.i = z__3.r * a[i__4].i + z__3.i * a[i__4]
				.r;
#line 380 "cherk.f"
			z__1.r = temp.r + z__2.r, z__1.i = temp.i + z__2.i;
#line 380 "cherk.f"
			temp.r = z__1.r, temp.i = z__1.i;
#line 381 "cherk.f"
/* L240: */
#line 381 "cherk.f"
		    }
#line 382 "cherk.f"
		    if (*beta == 0.) {
#line 383 "cherk.f"
			i__3 = i__ + j * c_dim1;
#line 383 "cherk.f"
			z__1.r = *alpha * temp.r, z__1.i = *alpha * temp.i;
#line 383 "cherk.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 384 "cherk.f"
		    } else {
#line 385 "cherk.f"
			i__3 = i__ + j * c_dim1;
#line 385 "cherk.f"
			z__2.r = *alpha * temp.r, z__2.i = *alpha * temp.i;
#line 385 "cherk.f"
			i__4 = i__ + j * c_dim1;
#line 385 "cherk.f"
			z__3.r = *beta * c__[i__4].r, z__3.i = *beta * c__[
				i__4].i;
#line 385 "cherk.f"
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 385 "cherk.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 386 "cherk.f"
		    }
#line 387 "cherk.f"
/* L250: */
#line 387 "cherk.f"
		}
#line 388 "cherk.f"
/* L260: */
#line 388 "cherk.f"
	    }
#line 389 "cherk.f"
	}
#line 390 "cherk.f"
    }

#line 392 "cherk.f"
    return 0;

/*     End of CHERK . */

} /* cherk_ */

