#line 1 "zsyrk.f"
/* zsyrk.f -- translated by f2c (version 20100827).
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

#line 1 "zsyrk.f"
/* > \brief \b ZSYRK */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC) */

/*       .. Scalar Arguments .. */
/*       COMPLEX*16 ALPHA,BETA */
/*       INTEGER K,LDA,LDC,N */
/*       CHARACTER TRANS,UPLO */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16 A(LDA,*),C(LDC,*) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZSYRK  performs one of the symmetric rank k operations */
/* > */
/* >    C := alpha*A*A**T + beta*C, */
/* > */
/* > or */
/* > */
/* >    C := alpha*A**T*A + beta*C, */
/* > */
/* > where  alpha and beta  are scalars,  C is an  n by n symmetric matrix */
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
/* >           TRANS = 'T' or 't',  K  specifies  the number of rows of the */
/* >           matrix A.  K must be at least zero. */
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
/* > \endverbatim */
/* > */
/*  ===================================================================== */
/* Subroutine */ int zsyrk_(char *uplo, char *trans, integer *n, integer *k, 
	doublecomplex *alpha, doublecomplex *a, integer *lda, doublecomplex *
	beta, doublecomplex *c__, integer *ldc, ftnlen uplo_len, ftnlen 
	trans_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, c_dim1, c_offset, i__1, i__2, i__3, i__4, i__5, 
	    i__6;
    doublecomplex z__1, z__2, z__3;

    /* Local variables */
    static integer i__, j, l, info;
    static doublecomplex temp;
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

#line 210 "zsyrk.f"
    /* Parameter adjustments */
#line 210 "zsyrk.f"
    a_dim1 = *lda;
#line 210 "zsyrk.f"
    a_offset = 1 + a_dim1;
#line 210 "zsyrk.f"
    a -= a_offset;
#line 210 "zsyrk.f"
    c_dim1 = *ldc;
#line 210 "zsyrk.f"
    c_offset = 1 + c_dim1;
#line 210 "zsyrk.f"
    c__ -= c_offset;
#line 210 "zsyrk.f"

#line 210 "zsyrk.f"
    /* Function Body */
#line 210 "zsyrk.f"
    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {
#line 211 "zsyrk.f"
	nrowa = *n;
#line 212 "zsyrk.f"
    } else {
#line 213 "zsyrk.f"
	nrowa = *k;
#line 214 "zsyrk.f"
    }
#line 215 "zsyrk.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);

#line 217 "zsyrk.f"
    info = 0;
#line 218 "zsyrk.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 219 "zsyrk.f"
	info = 1;
#line 220 "zsyrk.f"
    } else if (! lsame_(trans, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, 
	    "T", (ftnlen)1, (ftnlen)1)) {
#line 222 "zsyrk.f"
	info = 2;
#line 223 "zsyrk.f"
    } else if (*n < 0) {
#line 224 "zsyrk.f"
	info = 3;
#line 225 "zsyrk.f"
    } else if (*k < 0) {
#line 226 "zsyrk.f"
	info = 4;
#line 227 "zsyrk.f"
    } else if (*lda < max(1,nrowa)) {
#line 228 "zsyrk.f"
	info = 7;
#line 229 "zsyrk.f"
    } else if (*ldc < max(1,*n)) {
#line 230 "zsyrk.f"
	info = 10;
#line 231 "zsyrk.f"
    }
#line 232 "zsyrk.f"
    if (info != 0) {
#line 233 "zsyrk.f"
	xerbla_("ZSYRK ", &info, (ftnlen)6);
#line 234 "zsyrk.f"
	return 0;
#line 235 "zsyrk.f"
    }

/*     Quick return if possible. */

#line 239 "zsyrk.f"
    if (*n == 0 || (alpha->r == 0. && alpha->i == 0. || *k == 0) && (beta->r 
	    == 1. && beta->i == 0.)) {
#line 239 "zsyrk.f"
	return 0;
#line 239 "zsyrk.f"
    }

/*     And when  alpha.eq.zero. */

#line 244 "zsyrk.f"
    if (alpha->r == 0. && alpha->i == 0.) {
#line 245 "zsyrk.f"
	if (upper) {
#line 246 "zsyrk.f"
	    if (beta->r == 0. && beta->i == 0.) {
#line 247 "zsyrk.f"
		i__1 = *n;
#line 247 "zsyrk.f"
		for (j = 1; j <= i__1; ++j) {
#line 248 "zsyrk.f"
		    i__2 = j;
#line 248 "zsyrk.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 249 "zsyrk.f"
			i__3 = i__ + j * c_dim1;
#line 249 "zsyrk.f"
			c__[i__3].r = 0., c__[i__3].i = 0.;
#line 250 "zsyrk.f"
/* L10: */
#line 250 "zsyrk.f"
		    }
#line 251 "zsyrk.f"
/* L20: */
#line 251 "zsyrk.f"
		}
#line 252 "zsyrk.f"
	    } else {
#line 253 "zsyrk.f"
		i__1 = *n;
#line 253 "zsyrk.f"
		for (j = 1; j <= i__1; ++j) {
#line 254 "zsyrk.f"
		    i__2 = j;
#line 254 "zsyrk.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 255 "zsyrk.f"
			i__3 = i__ + j * c_dim1;
#line 255 "zsyrk.f"
			i__4 = i__ + j * c_dim1;
#line 255 "zsyrk.f"
			z__1.r = beta->r * c__[i__4].r - beta->i * c__[i__4]
				.i, z__1.i = beta->r * c__[i__4].i + beta->i *
				 c__[i__4].r;
#line 255 "zsyrk.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 256 "zsyrk.f"
/* L30: */
#line 256 "zsyrk.f"
		    }
#line 257 "zsyrk.f"
/* L40: */
#line 257 "zsyrk.f"
		}
#line 258 "zsyrk.f"
	    }
#line 259 "zsyrk.f"
	} else {
#line 260 "zsyrk.f"
	    if (beta->r == 0. && beta->i == 0.) {
#line 261 "zsyrk.f"
		i__1 = *n;
#line 261 "zsyrk.f"
		for (j = 1; j <= i__1; ++j) {
#line 262 "zsyrk.f"
		    i__2 = *n;
#line 262 "zsyrk.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 263 "zsyrk.f"
			i__3 = i__ + j * c_dim1;
#line 263 "zsyrk.f"
			c__[i__3].r = 0., c__[i__3].i = 0.;
#line 264 "zsyrk.f"
/* L50: */
#line 264 "zsyrk.f"
		    }
#line 265 "zsyrk.f"
/* L60: */
#line 265 "zsyrk.f"
		}
#line 266 "zsyrk.f"
	    } else {
#line 267 "zsyrk.f"
		i__1 = *n;
#line 267 "zsyrk.f"
		for (j = 1; j <= i__1; ++j) {
#line 268 "zsyrk.f"
		    i__2 = *n;
#line 268 "zsyrk.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 269 "zsyrk.f"
			i__3 = i__ + j * c_dim1;
#line 269 "zsyrk.f"
			i__4 = i__ + j * c_dim1;
#line 269 "zsyrk.f"
			z__1.r = beta->r * c__[i__4].r - beta->i * c__[i__4]
				.i, z__1.i = beta->r * c__[i__4].i + beta->i *
				 c__[i__4].r;
#line 269 "zsyrk.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 270 "zsyrk.f"
/* L70: */
#line 270 "zsyrk.f"
		    }
#line 271 "zsyrk.f"
/* L80: */
#line 271 "zsyrk.f"
		}
#line 272 "zsyrk.f"
	    }
#line 273 "zsyrk.f"
	}
#line 274 "zsyrk.f"
	return 0;
#line 275 "zsyrk.f"
    }

/*     Start the operations. */

#line 279 "zsyrk.f"
    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {

/*        Form  C := alpha*A*A**T + beta*C. */

#line 283 "zsyrk.f"
	if (upper) {
#line 284 "zsyrk.f"
	    i__1 = *n;
#line 284 "zsyrk.f"
	    for (j = 1; j <= i__1; ++j) {
#line 285 "zsyrk.f"
		if (beta->r == 0. && beta->i == 0.) {
#line 286 "zsyrk.f"
		    i__2 = j;
#line 286 "zsyrk.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 287 "zsyrk.f"
			i__3 = i__ + j * c_dim1;
#line 287 "zsyrk.f"
			c__[i__3].r = 0., c__[i__3].i = 0.;
#line 288 "zsyrk.f"
/* L90: */
#line 288 "zsyrk.f"
		    }
#line 289 "zsyrk.f"
		} else if (beta->r != 1. || beta->i != 0.) {
#line 290 "zsyrk.f"
		    i__2 = j;
#line 290 "zsyrk.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 291 "zsyrk.f"
			i__3 = i__ + j * c_dim1;
#line 291 "zsyrk.f"
			i__4 = i__ + j * c_dim1;
#line 291 "zsyrk.f"
			z__1.r = beta->r * c__[i__4].r - beta->i * c__[i__4]
				.i, z__1.i = beta->r * c__[i__4].i + beta->i *
				 c__[i__4].r;
#line 291 "zsyrk.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 292 "zsyrk.f"
/* L100: */
#line 292 "zsyrk.f"
		    }
#line 293 "zsyrk.f"
		}
#line 294 "zsyrk.f"
		i__2 = *k;
#line 294 "zsyrk.f"
		for (l = 1; l <= i__2; ++l) {
#line 295 "zsyrk.f"
		    i__3 = j + l * a_dim1;
#line 295 "zsyrk.f"
		    if (a[i__3].r != 0. || a[i__3].i != 0.) {
#line 296 "zsyrk.f"
			i__3 = j + l * a_dim1;
#line 296 "zsyrk.f"
			z__1.r = alpha->r * a[i__3].r - alpha->i * a[i__3].i, 
				z__1.i = alpha->r * a[i__3].i + alpha->i * a[
				i__3].r;
#line 296 "zsyrk.f"
			temp.r = z__1.r, temp.i = z__1.i;
#line 297 "zsyrk.f"
			i__3 = j;
#line 297 "zsyrk.f"
			for (i__ = 1; i__ <= i__3; ++i__) {
#line 298 "zsyrk.f"
			    i__4 = i__ + j * c_dim1;
#line 298 "zsyrk.f"
			    i__5 = i__ + j * c_dim1;
#line 298 "zsyrk.f"
			    i__6 = i__ + l * a_dim1;
#line 298 "zsyrk.f"
			    z__2.r = temp.r * a[i__6].r - temp.i * a[i__6].i, 
				    z__2.i = temp.r * a[i__6].i + temp.i * a[
				    i__6].r;
#line 298 "zsyrk.f"
			    z__1.r = c__[i__5].r + z__2.r, z__1.i = c__[i__5]
				    .i + z__2.i;
#line 298 "zsyrk.f"
			    c__[i__4].r = z__1.r, c__[i__4].i = z__1.i;
#line 299 "zsyrk.f"
/* L110: */
#line 299 "zsyrk.f"
			}
#line 300 "zsyrk.f"
		    }
#line 301 "zsyrk.f"
/* L120: */
#line 301 "zsyrk.f"
		}
#line 302 "zsyrk.f"
/* L130: */
#line 302 "zsyrk.f"
	    }
#line 303 "zsyrk.f"
	} else {
#line 304 "zsyrk.f"
	    i__1 = *n;
#line 304 "zsyrk.f"
	    for (j = 1; j <= i__1; ++j) {
#line 305 "zsyrk.f"
		if (beta->r == 0. && beta->i == 0.) {
#line 306 "zsyrk.f"
		    i__2 = *n;
#line 306 "zsyrk.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 307 "zsyrk.f"
			i__3 = i__ + j * c_dim1;
#line 307 "zsyrk.f"
			c__[i__3].r = 0., c__[i__3].i = 0.;
#line 308 "zsyrk.f"
/* L140: */
#line 308 "zsyrk.f"
		    }
#line 309 "zsyrk.f"
		} else if (beta->r != 1. || beta->i != 0.) {
#line 310 "zsyrk.f"
		    i__2 = *n;
#line 310 "zsyrk.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 311 "zsyrk.f"
			i__3 = i__ + j * c_dim1;
#line 311 "zsyrk.f"
			i__4 = i__ + j * c_dim1;
#line 311 "zsyrk.f"
			z__1.r = beta->r * c__[i__4].r - beta->i * c__[i__4]
				.i, z__1.i = beta->r * c__[i__4].i + beta->i *
				 c__[i__4].r;
#line 311 "zsyrk.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 312 "zsyrk.f"
/* L150: */
#line 312 "zsyrk.f"
		    }
#line 313 "zsyrk.f"
		}
#line 314 "zsyrk.f"
		i__2 = *k;
#line 314 "zsyrk.f"
		for (l = 1; l <= i__2; ++l) {
#line 315 "zsyrk.f"
		    i__3 = j + l * a_dim1;
#line 315 "zsyrk.f"
		    if (a[i__3].r != 0. || a[i__3].i != 0.) {
#line 316 "zsyrk.f"
			i__3 = j + l * a_dim1;
#line 316 "zsyrk.f"
			z__1.r = alpha->r * a[i__3].r - alpha->i * a[i__3].i, 
				z__1.i = alpha->r * a[i__3].i + alpha->i * a[
				i__3].r;
#line 316 "zsyrk.f"
			temp.r = z__1.r, temp.i = z__1.i;
#line 317 "zsyrk.f"
			i__3 = *n;
#line 317 "zsyrk.f"
			for (i__ = j; i__ <= i__3; ++i__) {
#line 318 "zsyrk.f"
			    i__4 = i__ + j * c_dim1;
#line 318 "zsyrk.f"
			    i__5 = i__ + j * c_dim1;
#line 318 "zsyrk.f"
			    i__6 = i__ + l * a_dim1;
#line 318 "zsyrk.f"
			    z__2.r = temp.r * a[i__6].r - temp.i * a[i__6].i, 
				    z__2.i = temp.r * a[i__6].i + temp.i * a[
				    i__6].r;
#line 318 "zsyrk.f"
			    z__1.r = c__[i__5].r + z__2.r, z__1.i = c__[i__5]
				    .i + z__2.i;
#line 318 "zsyrk.f"
			    c__[i__4].r = z__1.r, c__[i__4].i = z__1.i;
#line 319 "zsyrk.f"
/* L160: */
#line 319 "zsyrk.f"
			}
#line 320 "zsyrk.f"
		    }
#line 321 "zsyrk.f"
/* L170: */
#line 321 "zsyrk.f"
		}
#line 322 "zsyrk.f"
/* L180: */
#line 322 "zsyrk.f"
	    }
#line 323 "zsyrk.f"
	}
#line 324 "zsyrk.f"
    } else {

/*        Form  C := alpha*A**T*A + beta*C. */

#line 328 "zsyrk.f"
	if (upper) {
#line 329 "zsyrk.f"
	    i__1 = *n;
#line 329 "zsyrk.f"
	    for (j = 1; j <= i__1; ++j) {
#line 330 "zsyrk.f"
		i__2 = j;
#line 330 "zsyrk.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 331 "zsyrk.f"
		    temp.r = 0., temp.i = 0.;
#line 332 "zsyrk.f"
		    i__3 = *k;
#line 332 "zsyrk.f"
		    for (l = 1; l <= i__3; ++l) {
#line 333 "zsyrk.f"
			i__4 = l + i__ * a_dim1;
#line 333 "zsyrk.f"
			i__5 = l + j * a_dim1;
#line 333 "zsyrk.f"
			z__2.r = a[i__4].r * a[i__5].r - a[i__4].i * a[i__5]
				.i, z__2.i = a[i__4].r * a[i__5].i + a[i__4]
				.i * a[i__5].r;
#line 333 "zsyrk.f"
			z__1.r = temp.r + z__2.r, z__1.i = temp.i + z__2.i;
#line 333 "zsyrk.f"
			temp.r = z__1.r, temp.i = z__1.i;
#line 334 "zsyrk.f"
/* L190: */
#line 334 "zsyrk.f"
		    }
#line 335 "zsyrk.f"
		    if (beta->r == 0. && beta->i == 0.) {
#line 336 "zsyrk.f"
			i__3 = i__ + j * c_dim1;
#line 336 "zsyrk.f"
			z__1.r = alpha->r * temp.r - alpha->i * temp.i, 
				z__1.i = alpha->r * temp.i + alpha->i * 
				temp.r;
#line 336 "zsyrk.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 337 "zsyrk.f"
		    } else {
#line 338 "zsyrk.f"
			i__3 = i__ + j * c_dim1;
#line 338 "zsyrk.f"
			z__2.r = alpha->r * temp.r - alpha->i * temp.i, 
				z__2.i = alpha->r * temp.i + alpha->i * 
				temp.r;
#line 338 "zsyrk.f"
			i__4 = i__ + j * c_dim1;
#line 338 "zsyrk.f"
			z__3.r = beta->r * c__[i__4].r - beta->i * c__[i__4]
				.i, z__3.i = beta->r * c__[i__4].i + beta->i *
				 c__[i__4].r;
#line 338 "zsyrk.f"
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 338 "zsyrk.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 339 "zsyrk.f"
		    }
#line 340 "zsyrk.f"
/* L200: */
#line 340 "zsyrk.f"
		}
#line 341 "zsyrk.f"
/* L210: */
#line 341 "zsyrk.f"
	    }
#line 342 "zsyrk.f"
	} else {
#line 343 "zsyrk.f"
	    i__1 = *n;
#line 343 "zsyrk.f"
	    for (j = 1; j <= i__1; ++j) {
#line 344 "zsyrk.f"
		i__2 = *n;
#line 344 "zsyrk.f"
		for (i__ = j; i__ <= i__2; ++i__) {
#line 345 "zsyrk.f"
		    temp.r = 0., temp.i = 0.;
#line 346 "zsyrk.f"
		    i__3 = *k;
#line 346 "zsyrk.f"
		    for (l = 1; l <= i__3; ++l) {
#line 347 "zsyrk.f"
			i__4 = l + i__ * a_dim1;
#line 347 "zsyrk.f"
			i__5 = l + j * a_dim1;
#line 347 "zsyrk.f"
			z__2.r = a[i__4].r * a[i__5].r - a[i__4].i * a[i__5]
				.i, z__2.i = a[i__4].r * a[i__5].i + a[i__4]
				.i * a[i__5].r;
#line 347 "zsyrk.f"
			z__1.r = temp.r + z__2.r, z__1.i = temp.i + z__2.i;
#line 347 "zsyrk.f"
			temp.r = z__1.r, temp.i = z__1.i;
#line 348 "zsyrk.f"
/* L220: */
#line 348 "zsyrk.f"
		    }
#line 349 "zsyrk.f"
		    if (beta->r == 0. && beta->i == 0.) {
#line 350 "zsyrk.f"
			i__3 = i__ + j * c_dim1;
#line 350 "zsyrk.f"
			z__1.r = alpha->r * temp.r - alpha->i * temp.i, 
				z__1.i = alpha->r * temp.i + alpha->i * 
				temp.r;
#line 350 "zsyrk.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 351 "zsyrk.f"
		    } else {
#line 352 "zsyrk.f"
			i__3 = i__ + j * c_dim1;
#line 352 "zsyrk.f"
			z__2.r = alpha->r * temp.r - alpha->i * temp.i, 
				z__2.i = alpha->r * temp.i + alpha->i * 
				temp.r;
#line 352 "zsyrk.f"
			i__4 = i__ + j * c_dim1;
#line 352 "zsyrk.f"
			z__3.r = beta->r * c__[i__4].r - beta->i * c__[i__4]
				.i, z__3.i = beta->r * c__[i__4].i + beta->i *
				 c__[i__4].r;
#line 352 "zsyrk.f"
			z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 352 "zsyrk.f"
			c__[i__3].r = z__1.r, c__[i__3].i = z__1.i;
#line 353 "zsyrk.f"
		    }
#line 354 "zsyrk.f"
/* L230: */
#line 354 "zsyrk.f"
		}
#line 355 "zsyrk.f"
/* L240: */
#line 355 "zsyrk.f"
	    }
#line 356 "zsyrk.f"
	}
#line 357 "zsyrk.f"
    }

#line 359 "zsyrk.f"
    return 0;

/*     End of ZSYRK . */

} /* zsyrk_ */

